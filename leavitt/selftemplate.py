#!/usr/bin/env python
#
# Use data itself to design a smooth light curve template

from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
#%matplotlib
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import os
from dl import queryClient as qc
from astropy.table import Table, vstack
import utils
from collections import Counter
import psearch_py3
from scipy.signal import find_peaks, peak_prominences
from argparse import ArgumentParser
from dlnpyutils import utils as dln,bindata
import statsmodels.api as sm


def get_data(objname, bands = ['u','g','r','i','z','Y','VR']):
    """Query the object by name, extract light curves, 
       error, filters and top N estimated periods."""
    res=qc.query(sql="""SELECT mjd,mag_auto,magerr_auto,filter,fwhm
                        FROM nsc_dr2.meas 
                        WHERE objectid='{:s}'""".format(objname),
                 fmt='table')
    res['mjd'] += np.random.randn(len(res))*10**-10
    
    selbnds = [i for i, val in enumerate(res['filter']) if val in bands]
    selfwhm = np.where(res['fwhm'] <= 4.0)[0]
    sel = [x for x in selbnds if x in selfwhm]
    res = res[sel]
    
    if len(res) <= 0:
        raise Exception('No data with low fwhm')
    
    res['fltr']   = -1
    for i in range(len(res)):
        res['fltr'][i] = bands.index(res['filter'][i])
    
    res.rename_column('mag_auto', 'mag')
    res.rename_column('magerr_auto', 'err')
    res.sort(['fltr','mjd'])
    
    return res

def get_periods(mjd,mag,err,fltr,objname='',N = 5,pmin=.2,bands=['u','g','r','i','z','Y','VR']):
    
    # The filter information here uses indices determined from the order they
    # appear in bands. To run psearch we want to reassign these indices to remove
    # any unused bands. For example, if only 'g', 'r' and 'z' are used, indices
    # should be 0,1,2 and not 1,2,4.
    cnt  = Counter(fltr)
    
    mult = np.where(np.array(list(cnt.values()))>1)[0]
    sel  = np.in1d(fltr, mult)
    
    fltinds = list(set(fltr))
    replace = {fltinds[i]:i for i in range(len(fltinds))}
    newinds = np.array([replace.get(n,n) for n in fltr],dtype=np.float64)
    fltrnms = (np.array(bands))[list(set(fltr[sel]))]
    
    dphi = 0.02
    plist, psiarray, thresh = \
            psearch_py3.psearch_py( mjd[sel], mag[sel], err[sel], 
                                   newinds[sel], fltrnms, pmin, dphi )
    
    psi = psiarray.sum(0)
    
    pkinds = find_peaks(psi,distance=len(plist)/2000)[0]
    prom   = peak_prominences(psi,pkinds)[0]
    inds0  = pkinds[np.argsort(-prom)[:10*N]]
    inds   = inds0[np.argsort(-psi[inds0])[:N]]
    
    plot_periodogram(plist,psi,inds,objname)
    
    return plist[inds]

def plot_periodogram(prds,psi,inds,objname='',outdir='results/plots'):
   
    fig, ax = plt.subplots(figsize=(10,7))
        
    ax.plot(prds,psi,lw=0.1)
    ax.scatter(prds[inds[1:]],psi[inds[1:]],c='k',s=10)
    ax.scatter(prds[inds[0]],psi[inds[0]],c='r',s=12)
    
    ax.set_xlabel('log period (days)',fontsize=18)
    ax.set_ylabel('psi',fontsize=18)
    ax.set_title('{} Periodogram'.format(objname),fontsize=20)
    ax.set_xscale('log')
    ax.text(0.7,0.9,'best period = {:.3f} days'.format(prds[inds[0]]),transform=ax.transAxes,color='r')
    
#     fig.savefig(outdir+'\\{}_periodogram.png'.format(objname))
    
    # create zoomed in copy
    ax.set_title('{} Periodogram Zoomed In'.format(objname),fontsize=20)
    minp = min(prds[inds])
    maxp = max(prds[inds])
    ax.set_xlim(minp*.67,maxp*1.33)
    fig.savefig(outdir+'\\{}_periodogram_zoomedin.png'.format(objname))
    
    plt.close(fig)
    return


def selftemplate(cat,period,maxiter=5,minrmsdiff=0.02,verbose=False):
    """
    Generate template from data itself.

    Parameters
    ----------
    cat : astropy table or numpy structured array
       Input catalog of measurement data with "mjd", "mag", "err" and "fltr" columns.
    period : float
       Period to use in days.
    maxiter : int, optional
       Maximum iterations to try.  Default is 5.
    minrmsdiff : float, optional
       Maximum RMS difference of the template between iterations before the
         iteration loop is halted.  Default is 0.02.
    verbose : bool, optional
       Verbose output.  Default is False.

    Returns
    -------
    bands : numpy array
      Array of unique band names.
    pars : numpy array
      Array of parameters: [period, t0, amplitudes, mean magnitudes]
    template : numpy structured array
      Template information with columns "phase" and "flux".
    chisq : float
       Chi-squared of template, amplitudes and mean magnitudes.

    Example
    -------
    bands,pars,template,chisq = selftemplate(cat,0.6052)


    """
    if isinstance(cat,Table):
        cat = np.array(cat)
    t = cat['mjd']
    mag = cat['mag']
    err = cat['err']
    flter = cat['fltr']
    t0 = min(t)
    ph = (t - t0) / period %1

    bands = np.unique(flter)
    nbands = len(bands)

    # Generate "self" template iteratively
    flag = True
    niter = 0
    while (flag):
        if verbose:
            print('Niter = ',niter)
            
        # Initial amplitudes
        if niter==0:
            #  loop over all bands and get median of data in 0.25 phase chunks
            meds = np.zeros((nbands,4),float)
            numbin = np.zeros((nbands,4),int)
            amp = np.ones(nbands,float)
            mnmag = np.zeros(nbands,float)
            num = np.zeros(nbands,int)
            sclmag = mag.copy()*0
            for i,b in enumerate(bands):
                ind, = np.where(flter==b)
                num[i] = len(ind)
                if len(ind)>1:
                    ybin,bin_edges,binnumber = bindata.binned_statistic(ph[ind],mag[ind],statistic='median',bins=4,range=[0.0,1.0])
                    numbin1,bin_edges2,binnumber2 = bindata.binned_statistic(ph[ind],mag[ind],statistic='count',bins=4,range=[0.0,1.0])
                    meds[i,:] = ybin
                    numbin[i,:] = numbin1
                    amp[i] = np.nanstd(ybin)
                    if amp[i]>=0:
                        amp[i] = np.nanstd(mag[ind])
                    mnmag[i] = np.nanmedian(ybin)
                    if mnmag[i]<=0 or ~np.isfinite(mnmag[i]):
                        mnmag[i] = np.nanmedian(mag[ind])
                    if amp[i]>0:
                        sclmag[ind] = (mag[ind]-mnmag[i])/amp[i]
                    else:
                        sclmag[ind] = (mag[ind]-mnmag[i])                        
            # shift so the will be positive
            sclmag -= np.nanmin(sclmag)
            sclmag /= (3*np.nanstd(sclmag))  # scale to roughly a max of 1

            
        # Use existing template to get improved amplitudes and mean mags
        else:                
            # Loop over bands and solve for best amplitude and mean mag
            f = interp1d(xtemp,ytemp,kind='cubic',bounds_error=None,fill_value="extrapolate")
            amp = np.zeros(nbands,float)
            mnmag = np.zeros(nbands,float)
            num = np.zeros(nbands,int)
            sclmag = mag.copy()*0
            resid = mag.copy()*0
            for i,b in enumerate(bands):
                ind, = np.where(flter==b)
                num[i] = len(ind)
                if len(ind)>1:
                    temp = f(ph[ind])
                    # Make sure amplitude is non-negative
                    #amp[i] = np.maximum(dln.wtslope(temp,mag[ind],err[ind],reweight=False),0.0)
                    amp[i] = np.maximum(dln.mediqrslope(temp,mag[ind]),0.0)
                    mnmag[i] = np.median(mag[ind]-amp[i]*temp)
                    if amp[i]>0.0:
                        sclmag[ind] = (mag[ind]-mnmag[i])/amp[i]
                    else:
                        sclmag[ind] = (mag[ind]-mnmag[i])                        
                    resid[ind] = mag[ind]-(temp*amp[i]+mnmag[i])
            chisq = np.sum(resid**2/err**2)
            
        if verbose:
            print('Bands = ',bands)
            print('Amps = ',amp)
            print('Mnmag = ',mnmag)
            if niter>0:
                print('chisq = ',chisq)
                    
        # Add copy of data offset by one phase to left and right
        ph2 = np.concatenate((ph-1,ph,ph+1))
        sclmag2 = np.concatenate((sclmag,sclmag,sclmag))
        flter2 = np.concatenate((flter,flter,flter))
        keep, = np.where((ph2>=-0.25) & (ph2<=1.25))
        ph2 = ph2[keep]
        sclmag2 = sclmag2[keep]
        flter2 = flter2[keep]
        
        # Use LOWESS to generate empirical template
        # it will use closest frac*N data points to a given point to estimate the smooth version
        # want at least 5 points
        frac = np.maximum(10.0/len(ph2),0.05)
        lowess = sm.nonparametric.lowess(sclmag2,ph2, frac=frac)
        gd, = np.where((lowess[:,0]>=0.0) & (lowess[:,0]<=1.0))
        # interpolate onto fine grid, leave some overhang
        xtemp = np.linspace(-0.2,1.2,141)
        #xtemp = np.linspace(0.0,1.0,101)
        ytemp = interp1d(lowess[:,0],lowess[:,1],kind='quadratic',bounds_error=None,fill_value="extrapolate")(xtemp)
        
        # Shift t0 based on template minimum
        inbounds, = np.where((xtemp>=0.0) & (xtemp<=1.0))
        minind = np.argmin(ytemp[inbounds])
        phasemin = xtemp[inbounds[minind]]
        if verbose:
            print('phase min = ',phasemin)
        if phasemin>0.5:
            phasemin -= 1.0
        if np.abs(phasemin)>0.01 and np.abs(phasemin-1.0)>0.01:
            if verbose:
                print('shifting phase minimum by %8.4f' % phasemin)
            timeoffset = phasemin*period            
            t0 += timeoffset
            ph = (t - t0) / period %1
            # Repeat xtemp/ytemp left and right to keep full coverage
            xtemp = np.concatenate((xtemp[inbounds]-1,xtemp[inbounds],xtemp[inbounds]+1))
            ytemp = np.concatenate((ytemp[inbounds],ytemp[inbounds],ytemp[inbounds]))
            xtemp -= phasemin
            si = np.argsort(xtemp)  # sort again
            xtemp = xtemp[si]
            ytemp = ytemp[si]
            temp,ui = np.unique(xtemp,return_index=True)  # make sure they are unique
            xtemp = xtemp[ui]
            ytemp = ytemp[ui]
            gd, = np.where((xtemp>=-0.2) & (xtemp<=1.2))
            xtemp = xtemp[gd]
            ytemp = ytemp[gd]
        
        # Scale so the template values are between 0 and 1
        ytemp -= np.min(ytemp)
        ytemp /= np.max(ytemp)
            
        if niter>0:
            f = interp1d(xtemp,ytemp,kind='cubic',bounds_error=None,fill_value="extrapolate")
            temp = f(xtemp_last)
            rms = np.sqrt(np.mean((ytemp_last-temp)**2))
        else:
            rms = 999999.
        if verbose:
            print('RMS = ',rms)

        if niter>0 and (niter>maxiter or rms<minrmsdiff): flag=False
                
        xtemp_last = xtemp.copy()
        ytemp_last = ytemp.copy()
        niter += 1

       # return xtemp,ytemp
        
    # Trim template phase range to 0-1
    gd, = np.where((xtemp>=0) & (xtemp<=1.0))
    xtemp = xtemp[gd]
    ytemp = ytemp[gd]
    template = np.zeros(len(xtemp),dtype=np.dtype([('phase',float),('flux',float),]))
    template['phase'] = xtemp
    template['flux'] = ytemp

    # Put together the parameters
    pars = np.concatenate((np.array([period]),np.array([t0]),amp,mnmag))

    if verbose:
        print('Final parameters = ',pars)
    
    return bands,pars,template,chisq


def model(cat,template,pars):
    """
    Return the template model for data points.

    Parameters
    ----------
    cat : astropy table or numpy structured array
       Input catalog of measurement data with "mjd", "mag", "err" and "fltr" columns.
    template : numpy structured array
       Template information with columns "phase" and "flux".
    pars : numpy array
       Array of parameters: [period, t0, amplitudes, mean magnitudes]

    Returns
    -------
    modelmag : numpy array
       Model magnitudes contructed from the template.

    Example
    -------

    modelmag = model(cat,template,pars)

    """

    period = pars[0]
    t1 = pars[1]
    bands = np.unique(cat['fltr'])
    nbands = len(bands)
    amp = pars[2:2+nbands]
    mnmag = pars[-nbands:]
    ph = (cat['mjd'] - pars[1]) / period %1
    modelmag = np.zeros(len(cat),float)
    f = interp1d(template['phase'],template['flux'])
    for b,i in enumerate(bands):
        ind, = np.where(cat['fltr']==b)
        temp = f(ph[ind])
        modelmag[ind] = temp*amp[i]+mnmag[i]

    return modelmag

    
def scaledmags(cat,template,pars):
    """
    Return the scaled observed magnitudes so they
    can be easily compared to the template.

    Parameters
    ----------
    cat : astropy table or numpy structured array
       Input catalog of measurement data with "mjd", "mag", "err" and "fltr" columns.
    template : numpy structured array
       Template information with columns "phase" and "flux".
    pars : numpy array
       Array of parameters: [period, t0, amplitudes, mean magnitudes]

    Returns
    -------
    ph : numpy array
       Array of phases.
    sclmag : numpy array
       Scaled observed magnitudes.

    Example
    -------

    sclmag = scaledmags(cat,template,pars)

    """

    period = pars[0]
    t1 = pars[1]
    bands = np.unique(cat['fltr'])
    nbands = len(bands)
    amp = pars[2:2+nbands]
    mnmag = pars[-nbands:]
    ph = (cat['mjd'] - pars[1]) / period %1
    sclmag = np.zeros(len(cat),float)
    for b,i in enumerate(bands):
        ind, = np.where(cat['fltr']==b)
        sclmag[ind] = (cat['mag'][ind]-mnmag[i])/amp[i]

    return ph,sclmag
    

class RRLfitter:
    def __init__ (self, tmps, fltnames= ['u','g','r','i','z','Y','VR'],
                  ampratio=[1.81480451,1.46104910,1.0,0.79662171,0.74671563,0.718746,1.050782]):
        # constants
        self.tmps     = tmps # Table containing templates
        self.fltnames = fltnames # list of names of usable filters
        self.Nflts    = len(fltnames) # number of usable filters
        self.ampratio = np.array(ampratio)
        # model variables
        self.fltinds  = [] # list of filter index values (0:'u', 1:'g', etc.)
        self.tmpind   = 1 # index of template currently being used 1,2,...,N
        self.period   = 1
        
    def model(self, t, *args):
        """modify the template using peak-to-peak amplitude and yoffset
        input times t should be epoch folded, phase shift to match template"""
        t0 = args[0]
        amplist = (args[1] * self.ampratio)[self.fltinds]
        yofflist = np.array(args[2:])[self.fltinds]
        
        ph = (t - t0) / self.period %1
        template = interp1d(self.tmps.columns[0],self.tmps.columns[self.tmpind])(ph)
        
        mag = template * amplist + yofflist
        
        return mag

    def tmpfit(self,mjd,mag,err,fltinds,plist,initpars=None):
        self.fltinds = fltinds
        if isinstance(plist, (int,float)):
            plist = [plist]
        
        
        if initpars is None:
            initpars = np.zeros( 2 + self.Nflts )
            initpars[0]  = min(mjd)
            initpars[2:] = np.median(mag)
            ampest = []
            for f in set(fltinds):
                ampest.append( (max(mag[fltinds==f])-min(mag[fltinds==f]))/self.ampratio[f] )
            initpars[1]  = np.mean(ampest)

        bounds = ( np.zeros(2+self.Nflts), np.zeros(2+self.Nflts))
        bounds[0][0] =  0.0
        bounds[1][0] = np.inf
        bounds[0][1] =  0.0
        bounds[1][1] = 50.0
        bounds[0][2:]=-50.0
        bounds[1][2:]= 50.0

        for i in set(range(self.Nflts))-set(self.fltinds):
            initpars[2+i]  =   0
            bounds[0][2+i] = -10**-6
            bounds[1][2+i] =  10**-6
        
        minx2    = 2**99
        bestpars = np.zeros( 2 + self.Nflts )
        besttmp  =-1
        besterr  = 0
        bestprd  = 0
        for p in plist:
            self.period = p
            
            for n in range(1,len(self.tmps.columns)):
                self.tmpind = n
                
                try:
                    pars, cov = curve_fit(self.model, mjd, mag, 
                                          bounds=bounds, sigma=err,
                                          p0=initpars, maxfev=5000)
                except RuntimeError:
                    continue
                x2 = sum((self.model(mjd,*pars)-mag)**2/err**2)
                print(p,n,pars,x2)
                if x2 < minx2:
                    minx2 = x2
                    bestpars = pars
                    besterr = np.sqrt(np.diag(cov))
                    bestprd = p
                    besttmp = n
                    
        self.period = bestprd
        self.tmpind = besttmp
        
        return bestpars, bestprd, besterr, besttmp, minx2

    def fit_plot(self,objname,N=10):
        print('getting data')
        crvdat = get_data(objname,bands=self.fltnames)
        print(str(len(crvdat))+' data points')
        
        print('Getting period')
        #plist  = get_periods(crvdat['mjd'],crvdat['mag'],crvdat['err'],crvdat['fltr'],
        #                     objname=objname,bands=self.fltnames,N=10)
        #print('periods: ',plist)
        period = 0.60527109
        
        self.selftemplate(crvdat,period)
        
        # Fit curve
        print('First fitting')
        pars,p,err,tmpind,chi2 = self.tmpfit(crvdat['mjd'],crvdat['mag'],crvdat['err'],crvdat['fltr'],plist)
        
        # Reject outliers, select inliers
        resid   = np.array(abs(crvdat['mag']-self.model(crvdat['mjd'],*pars)))
        crvdat['inlier'] = resid<utils.mad(resid)*5
        print(str(len(np.sum(~crvdat['inlier'])))+' outliers rejected')
        
        # Fit with inliers only
        print('Second fitting')
        pars,p,err,tmpind,chi2 = self.tmpfit(crvdat['mjd'][crvdat['inlier']],crvdat['mag'][crvdat['inlier']],
                                             crvdat['err'][crvdat['inlier']],crvdat['fltr'][crvdat['inlier']],plist,pars)
        
        redchi2 = chi2/(sum(crvdat['inlier'])-len(set(crvdat['fltr'][crvdat['inlier']]))-2)
        
        # get the filters with inlier data (incase it's different from all data)
        inlierflts = set(crvdat['fltr'][crvdat['inlier']])
        # Add phase to crvdat and sort
        crvdat['ph'] = ph = (crvdat['mjd'] - pars[0]) / p %1
        crvdat.sort(['fltr','ph'])
        self.fltinds = crvdat['fltr']
        
        # Plot
        colors  = ['#1f77b4','#2ca02c','#d62728','#9467bd','#8c564b','y','k']
        nf      = len(inlierflts) # Number of filters with inliers
        fig, ax = plt.subplots(nf, figsize=(12,4*(nf**.75+1)), sharex=True)
        if nf == 1:
            ax  = [ax]
        
        for i,f in enumerate(inlierflts):
            sel = crvdat['fltr'] == f
            ax[i].scatter(crvdat['ph'][sel],crvdat['mag'][sel],c=colors[f])
            ax[i].scatter(crvdat['ph'][sel]+1,crvdat['mag'][sel],c=colors[f])
            tmpmag = np.tile(self.tmps.columns[tmpind]*pars[1]*self.ampratio[f]+pars[2:][f],2)
            tmpph  = np.tile(self.tmps['PH'],2)+([0]*len(self.tmps['PH'])+[1]*len(self.tmps['PH']))
            ax[i].plot(tmpph,tmpmag,c='k')
            ax[i].invert_yaxis()
            ax[i].set_ylabel(self.fltnames[f], fontsize=20)
        
        ax[-1].set_xlabel('Phase', fontsize=20)
        ax[0].set_title("Object: {}    Period: {:.3f} d    Type: {}".format(
                                            objname,p,self.tmps.colnames[tmpind]), fontsize=22)
        fig.savefig('results/plots/{}_plot.png'.format(objname))
        plt.close(fig)
        
        # save parameters and results
        res = Table([[objname]],names=['name'])
        res['period'] = p
        res['t0']     = pars[0]
        res['r amp']  = pars[1]
        for i in range(2,len(pars)):
            f = self.fltnames[i-2]
            res['{} mag'.format(f)] = pars[i]
        res['chi2']   = chi2
        res['redchi2']= redchi2
        res['template']= self.tmps.colnames[tmpind]
        res['t0 err']     = err[0]
        res['amp err']  = err[1]
        for i in range(2,len(err)):
            f = self.fltnames[i-2]
            res['{} mag err'.format(f)] = err[i]
        res['Ndat']      = len(crvdat)
        res['N inliers'] = sum(crvdat['inlier'])
        for i in range(len(self.fltnames)):
            f = self.fltnames[i]
            res['N {}'.format(f)] = sum(crvdat['fltr'][crvdat['inlier']]==i)
        res.write('results/{}_res.fits'.format(objname),format='fits',overwrite=True)
        
        return
        

# Main command-line program
if __name__ == "__main__":
    parser = ArgumentParser(description='Run self template on star')
    parser.add_argument('objectid', type=str, nargs='+', help='Object ID')

    tmps = Table.read('templates/layden_templates.fits',format='fits')['PH','RRA1','RRA2','RRA3','RRB1','RRB2','RRB3','RRC']
    fitter  = RRLfitter(tmps,['u','g','r','i','z','Y','VR'],[1.8148,1.4610,1.0,0.7966,0.7467,0.7187,1.0507])


    objectid = '93142_19513'
    fitter.fit_plot(objectid)
