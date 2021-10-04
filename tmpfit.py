import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, peak_prominences
from scipy.interpolate import interp1d
import os
from dl import queryClient as qc
from astropy.table import Table
import utils
from collections import Counter
import psearch_py3 as psearch

def get_data(objname, bands = ['u','g','r','i','z','Y','VR']):
    """Query the object by name, extract light curves, 
       error, filters and top N estimated periods."""
    res=qc.query(sql="""SELECT mjd,mag_auto,magerr_auto,filter,fwhm
                        FROM nsc_dr2.meas 
                        WHERE objectid='{:s}'""".format(objname),
                 fmt='table')
    
    selbnds = [i for i, val in enumerate(res['filter']) if val in bands]
    selfwhm = np.where(res['fwhm'] <= 4.0)[0]
    sel = [x for x in selbnds if x in selfwhm]
    res = res[sel]
    
    res['fltr']   = -1
    for i in range(len(res)):
        res['fltr'][i] = bands.index(res['filter'][i])
    
    res.rename_column('mag_auto', 'mag')
    res.rename_column('magerr_auto', 'err')
    res.sort(['fltr','mjd'])
    
    return res

def get_periods(mjd,mag,err,fltr,objname='',outdir='results/plots',N = 10,
                pmin=.2,bands=['u','g','r','i','z','Y','VR'],verbose=False):
    
    # The filter information here uses indices determined from the order they
    # appear in bands. To run psearch we want to reassign these indices to remove
    # any unused bands. For example, if only 'g', 'r' and 'z' are used, indices
    # should be 0,1,2 and not 1,2,4.
    
    
    fltlist = np.array(fltr)
    
    sel1 = [i for i,val in enumerate(err) if val<.2]
    sel  = [i for i,val in enumerate(fltlist) if np.sum(fltlist[sel1]==val)>1]
    
    fltrnms = np.array(bands)[np.unique(fltr[sel])]
    mapping = {v:k for k,v in enumerate(set(fltr[sel]))}
    modflts = np.array([mapping[y] for y in fltr[sel]]).astype(np.float64)
    
    
    if len(fltrnms) < 1:
        raise Exception("No bands with multiple data, can't run PSearch.")
        
    
    print(Counter(modflts))
    dphi = 0.02
    plist, psiarray, thresh = \
            psearch.psearch_py( mjd[sel], mag[sel], err[sel], 
                                modflts, fltrnms, pmin, dphi, verbose=verbose )
    
    if len(fltrnms) > 1:
        psi = psiarray.sum(0)
    elif len(fltrnms) == 1:
        psi = psiarray
    
    pkinds = find_peaks(psi,distance=len(plist)**.45)[0]
    prom   = peak_prominences(psi,pkinds)[0]
    inds0  = pkinds[np.argsort(-prom)[:10*N]]
    inds   = inds0[np.argsort(-psi[inds0])[:N]]
    plot_periodogram(plist,psi,inds,objname,outdir)
    
    return plist, psi, inds

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
    
    fig.savefig(outdir+'\\{}_periodogram.png'.format(objname))
    
    # create zoomed in copy
    # ax.set_title('{} Periodogram Zoomed In'.format(objname),fontsize=20)
    # minp = min(prds[inds])
    # maxp = max(prds[inds])
    # ax.set_xlim(minp*.67,maxp*1.33)
    # fig.savefig(outdir+'\\{}_periodogram_zoomedin.png'.format(objname))
    
    plt.close(fig)
    return

class tmpfitter:
    """
    Object used to fit templates to data. Initialize with the templates you plan
    to compare against and lists of the bands and amplitude ratios for those bands
    Templates should be in an Astropy table with columns for the phase and each
    unique template. Column names are assumed to be template names.
    """
    def __init__ (self, tmps, fltnames= ['u','g','r','i','z','Y','VR']):
        # constants
        self.tmps     = tmps # Table containing templates
        self.fltnames = fltnames # list of names of usable filters
        self.Nflts    = len(fltnames) # number of usable filters
        # model variables
        self.fltinds  = [] # list of filter index values (0:'u', 1:'g', etc.)
        self.tmpind   = 1 # index of template currently being used 1,2,...,N
        self.period   = 1
        
    def model(self, t, *args):
        """modify the template using peak-to-peak amplitude and yoffset
        input times t should be epoch folded, phase shift to match template"""
        t0 = args[0]
        amplist = np.array(args[1:-self.Nflts])[self.fltinds]
        yofflist = np.array(args[-self.Nflts:])[self.fltinds]
        
        ph = (t - t0) / self.period %1
        template = interp1d(self.tmps.columns[0],self.tmps.columns[self.tmpind])(ph)
        
        mag = template * amplist + yofflist
        return mag

    def tmpfit(self,mjd,mag,err,fltinds,plist,initpars=None, verbose=False):
        self.fltinds = fltinds
        if isinstance(plist, (int,float)):
            plist = [plist]
            
        if initpars is None:
            initpars = np.zeros( 1 + 2*self.Nflts )
            initpars[0]  = min(mjd)
            
            for f in set(fltinds):
                initpars[1 + f] = max(mag[fltinds==f]) - min(mag[fltinds==f])
                initpars[1 + self.Nflts + f] = min(mag[fltinds==f])
                
            if verbose:
                print("Initial Parameters:")
                print("t0:",initpars[0])
                print("amps:",initpars[1:self.Nflts+1])
                print("y offset:",initpars[-self.Nflts:])
        
        bounds = ( np.zeros(1+2*self.Nflts), np.zeros(1+2*self.Nflts))
        bounds[0][0]  =  0.0
        bounds[1][0]  = np.inf
        bounds[1][1:] = 50.0
        bounds[1][2:] = 50.0
        bounds[0][1:-self.Nflts] =  0.0
        bounds[0][-self.Nflts:]  =-50.0

        for i in set(range(self.Nflts))-set(self.fltinds):
            initpars[1+i]  =   0
            bounds[0][1+i] = -10**-6
            bounds[1][1+i] =  10**-6
            initpars[1  + self.Nflts + i] =   0
            bounds[0][1 + self.Nflts + i] = -10**-6
            bounds[1][1 + self.Nflts + i] =  10**-6
        
        minx2    = 2**99
        bestpars = np.zeros( 1 + 2*self.Nflts )
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
                                          p0=initpars, maxfev=9000)
                except RuntimeError:
                    continue
                    
                x2 = sum((self.model(mjd,*pars)-mag)**2/err**2)
                if x2 < minx2:
                    minx2 = x2
                    bestpars = pars
                    besterr = np.sqrt(np.diag(cov))
                    bestprd = p
                    besttmp = n
                    
        self.period = bestprd
        self.tmpind = besttmp
        
        if verbose:
            print("Results: ")
            print("t0:",bestpars[0])
            print("amps:",bestpars[1:self.Nflts+1])
            print("y offset:",bestpars[-self.Nflts:])
            print("Period:",bestprd)
            print("Best Template:",self.tmps.colnames[besttmp])
            print("Chi Square:",minx2)
            
        return bestpars, bestprd, besterr, besttmp, minx2

def fit_plot(fitter,objname,plist=None,N=10,verbose=False):
    if verbose:
        print('Get data')
    crvdat = get_data(objname,bands=fitter.fltnames)
    if plist is None:
        if verbose:
            print('Get Periods')
        ps,psi,inds  = get_periods(crvdat['mjd'],crvdat['mag'],crvdat['err'],crvdat['fltr'],
                             objname=objname,bands=fitter.fltnames,N=N)
        plist = ps[inds]
    if verbose:
        print('First Fit')
    # Fit curve
    pars,p,err,tmpind,chi2 = fitter.tmpfit(crvdat['mjd'],crvdat['mag'],crvdat['err'],crvdat['fltr'],plist,verbose=verbose)
    
    if verbose:
        print('Outlier Rejection')   
    # Reject outliers, select inliers
    resid   = np.array(crvdat['mag']-fitter.model(crvdat['mjd'],*pars))
    crvdat['inlier'] = abs(resid)<utils.mad(resid)*5
    
    if verbose:
        print('Second Fit')
    # Fit with inliers only
    pars,p,err,tmpind,chi2 = fitter.tmpfit(crvdat['mjd'][crvdat['inlier']],crvdat['mag'][crvdat['inlier']],
                                         crvdat['err'][crvdat['inlier']],crvdat['fltr'][crvdat['inlier']],plist,pars,verbose=verbose)
    
    redchi2 = chi2/(sum(crvdat['inlier'])-len(set(crvdat['fltr'][crvdat['inlier']]))-2)
    
    # get the filters with inlier data (incase it's different from all data)
    inlierflts = set(crvdat['fltr'][crvdat['inlier']])
    
    # Add phase to crvdat and sort
    crvdat['ph'] = (crvdat['mjd'] - pars[0]) / p %1
    crvdat.sort(['fltr','ph'])
    fitter.fltinds = crvdat['fltr']
    
    if verbose:
        print('Start Plotting')
    # Plot
    colors  = ['#1f77b4','#2ca02c','#d62728','#9467bd','#8c564b','y','#ff7f0e']
    nf      = len(inlierflts) # Number of filters with inliers
    fig, ax = plt.subplots(nf, figsize=(12,4*(nf**.75+1)), sharex=True)
    if nf == 1:
        ax  = [ax]
        
    for i,f in enumerate(inlierflts):
        sel    = crvdat['fltr'] == f
        ax[i].scatter(crvdat['ph'][sel],crvdat['mag'][sel],c=colors[f])
        ax[i].scatter(crvdat['ph'][sel]+1,crvdat['mag'][sel],c=colors[f])
        tmpmag = np.tile(fitter.tmps.columns[tmpind]*pars[1+f]+pars[1+fitter.Nflts+f],2)
        tmpph  = np.tile(fitter.tmps['PH'],2)+([0]*len(fitter.tmps['PH'])+[1]*len(fitter.tmps['PH']))
        ax[i].plot(tmpph,tmpmag,c='k')
        xsel   = sel*(~crvdat['inlier'])
        ax[i].scatter(crvdat['ph'][xsel],crvdat['mag'][xsel],c='k',marker='x')
        ax[i].scatter(crvdat['ph'][xsel]+1,crvdat['mag'][xsel],c='k',marker='x')
        ax[i].invert_yaxis()
        ax[i].set_ylabel(fitter.fltnames[f], fontsize=20)
    
    ax[-1].set_xlabel('Phase', fontsize=20)
    ax[0].set_title("Object: {}    Period: {:.3f} d    Type: {}".format(
                                        objname,p,fitter.tmps.colnames[tmpind]), fontsize=22)
    path = 'results/plots/{}_plot.png'.format(objname)
    fig.savefig(path)
    if verbose:
        print('Saved to',path)
        print('Save parameters')
    plt.close(fig)
    
    # save parameters and results
    res = Table([[objname]],names=['name'])
    res['period'] = p
    res['t0']     = pars[0]
    res['t0 err'] = err[0]
    for i in range(fitter.Nflts):
        f = fitter.fltnames[i]
        res['{} amp'.format(f)]     = pars[1+i]
        res['{} amp err'.format(f)] = err[1+i]
        res['{} mag'.format(f)]     = pars[1+fitter.Nflts+i]
        res['{} mag err'.format(f)] = err[1+fitter.Nflts+i]
    res['chi2']       = chi2
    res['redchi2']    = redchi2
    res['template']   = fitter.tmps.colnames[tmpind]
    res['Ndat']       = len(crvdat)
    res['N outliers'] = len(crvdat) - sum(crvdat['inlier'])
    for i in range(fitter.Nflts):
        f = fitter.fltnames[i]
        res['N {}'.format(f)] = sum(crvdat['fltr'][crvdat['inlier']]==i)
    path = 'results/{}_res.fits'.format(objname)
    res.write(path,format='fits',overwrite=True)
    
    if verbose:
        print('Saved to',path)
        print('end')
    return