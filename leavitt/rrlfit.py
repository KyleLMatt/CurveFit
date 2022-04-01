import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, peak_prominences
from scipy.interpolate import interp1d
from dl import queryClient as qc
from astropy.table import Table
import utils
from collections import Counter
import psearch_py3 as psearch
import statsmodels.api as sm
import os

def get_data(nms, bands = ['u','g','r','i','z','Y','VR']):
    """Query list of objects by name, extract light curves, 
       error, filters and top N estimated periods."""
    res = qc.query(sql="""select objectid,mjd,mag_auto,magerr_auto,filter,fwhm
                          from nsc_dr2.meas 
                          where objectid in {}""".format(tuple(nms)),fmt='table')
    res['mjd'] += np.random.randn(len(res))*1e-10
    
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

def rollAve(lst, k=5):
    """
    Returns the rolling average of the list,
    averaging k terms at a time.
    """
    ret = np.cumsum(lst)
    ret[k:] = ret[k:] - ret[:-k]
    
    return ret[k - 1:] / k

def gaussAveSlow(t,y,std=.01):
    """
    Does an average in y with weights based on distances to all other points in t.
    """
    assert len(y) == len(t)
    n = len(t)
    dists = (t.reshape(n,1) - t.reshape(1,n))
    dists[dists> .5] -= 1
    dists[dists<-.5] += 1
    w = np.exp(-.5*dists**2/std**2)
    return np.sum( w*y, axis=1)/np.sum(w,axis=1)

def gaussAve(t,y,std=0.01,nei=50):
    n = len(t)
    out = np.zeros(n,float)
    for i in range(n):
        if i-nei<0:
            t1 = t[0:i+nei]
            y1 = y[0:i+nei]
            t2 = t[(i-nei):]-1.0
            y2 = y[(i-nei):]
            t1 = np.hstack((t1,t2))
            y1 = np.hstack((y1,y2))
        elif i+nei>n:
            t1 = t[i-nei:]
            y1 = y[i-nei:]
            t2 = t[:(i+nei)-n]+1.0
            y2 = y[:(i+nei)-n]
            t1 = np.hstack((t1,t2))
            y1 = np.hstack((y1,y2))	    
        else:
            lo = np.maximum(i-nei,0)
            hi = np.minimum(i+nei,n)
            t1 = t[lo:hi]
            y1 = y[lo:hi]
        deltat = t1-t[i]
        w = np.exp(-.5*deltat**2/std**2)
        out[i] = np.sum(w*y1)/np.sum(w)
    return out


def checkHarmonics(t,y,p,hlim=5,std=.01,fast=None):
    testp = p/np.union1d(np.arange(1,hlim+1),1/(np.arange(1,hlim+1)))
    return p_test(t,y,testp[testp>.1],std=std,fast=fast)

def p_refine(t,y,p,pwidth=.01,N=1000,std=.01,fast=None):
    testp = np.linspace(p-p*pwidth,p+p*pwidth,N)
    testp = np.union1d(testp,p)
    return p_test(t,y,testp[testp>.1],std=std,fast=fast)

def p_test(t,y,plst,std=.01,fast=None):
    """
    Find period from plst that produces the minimum string length.
    Gaussian weighted averages may be more accurate but will be 
    slower than the rolling average.
    """
    if fast is None:
        if 3.6e-8 * len(plst) * len(t)**2 > 75:
            fast = True
        else:
            fast = False
    mindl = 9**99
    besti = 0
    for i,prd in enumerate(plst):
        phase = t/prd %1
        mlst = y[np.argsort(phase)]
        phase.sort()
        
        if fast:
            k    = int(len(mlst)**.75*std*np.log(4)**.5+3)
            mlst = rollAve(np.append(mlst,mlst[:k-1]),k=k)
            mlst = np.roll(mlst,int(k/2))
        else:
            mlst  = gaussAve(phase,mlst,std=std)
        
        dy = np.abs( np.roll(mlst,-1) - mlst )
        dt = (np.roll(phase,-1) - phase) 
        dt[-1] += 1
        dl = (dt**2+dy**2)**.5
        if np.sum(dl) < mindl:
            mindl = np.sum(dl)
            besti = i
    return plst[besti]

class RRLfitter:
    """
    Object used to fit templates to data. Initialize with the templates you plan
    to compare against and lists of the bands and amplitude ratios for those bands
    Templates should be in an Astropy table with columns for the phase and each
    unique template. Column names are assumed to be template names.
    """
    def __init__ (self, tmps, fltnames= ['u','g','r','i','z','Y','VR'], ampratio=[1.81480451,1.46104910,1.0,0.79662171,0.74671563,0.718746,1.050782]):
        # constants
        self.tmps     = tmps # Table containing templates
        self.fltnames = fltnames # list of names of usable filters
        self.Nflts    = len(fltnames) # number of usable filters
        self.ampratio = np.array(ampratio)
        # model variables
        self.fltinds  = [] # list of filter index values (0:'u', 1:'g', etc.)
        self.tmpind   = 1  # index of template (1,2,...,N) being used.
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

    def tmpfit(self,mjd,mag,err,fltinds,plist,initpars=None, verbose=False):
        self.fltinds = fltinds
        if isinstance(plist, (int,float)):
            plist = [plist]
            
        if initpars is None:
            initpars = np.zeros( 2 + self.Nflts )
            initpars[0]  = min(mjd)
            ampest = []
            for f in set(fltinds):
                initpars[f+2] = min(mag[fltinds==f])
                fampest = (max(mag[fltinds==f])-min(mag[fltinds==f]))/self.ampratio[f]
                ampest.append(fampest)
                
            initpars[1]  = np.mean(ampest)
        
            if verbose:
                print("Initial Parameters:")
                print("t0:",initpars[0])
                print("r amp:",initpars[1])
                print("y offset:",initpars[2:])
        
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
            print("r amp:",bestpars[1])
            print("y offset:",bestpars[2:])
            print("Period:",bestprd)
            print("Best Template:",self.tmps.colnames[besttmp])
            print("Chi Square:",minx2)
            
        return bestpars, bestprd, besterr, besttmp, minx2

def fit_plot(fitter,objname,cat=None,plist=None,N=10,verbose=False,dirpath=''):

    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
    path = dirpath

    
    if cat is None:
        if verbose:
            print('Get data')
        cat = get_data(objname,bands=fitter.fltnames)
    if plist is None:
        if verbose:
            print('Get Periods')
        ps,psi,inds  = get_periods(cat['mjd'],cat['mag'],cat['err'],cat['fltr'],
                             objname=objname,bands=fitter.fltnames,N=N)
        plist = ps[inds]
    if verbose:
        print('First Fit')
    # Fit curve
    pars,p,err,tmpind,chi2 = fitter.tmpfit(cat['mjd'],cat['mag'],cat['err'],cat['fltr'],plist,verbose=verbose)
    
    if verbose:
        print(pars)
        print('Outlier Rejection')   
    # Reject outliers, select inliers
    resid   = np.array(cat['mag']-fitter.model(cat['mjd'],*pars))
    cat['inlier'] = abs(resid)<utils.mad(resid)*5
    
    if verbose:
        print(sum(~cat['inlier']),'outliers rejected')
        print('Second Fit')
    # Fit with inliers only
    pars,p,err,tmpind,chi2 = fitter.tmpfit(cat['mjd'][cat['inlier']],cat['mag'][cat['inlier']],
                                         cat['err'][cat['inlier']],cat['fltr'][cat['inlier']],plist,pars,verbose=verbose)
    
    redchi2 = chi2/(sum(cat['inlier'])-len(set(cat['fltr'][cat['inlier']]))-2)
    if verbose:
        print(pars)
    # get the filters with inlier data (incase it's different from all data)
    inlierflts = set(cat['fltr'][cat['inlier']])
    
    # Add phase to cat and sort
    cat['ph'] = (cat['mjd'] - pars[0]) / p %1
    cat.sort(['fltr','ph'])
    fitter.fltinds = cat['fltr']
    
    if verbose:
        print('Start Plotting')
    # Plot
    colors  = ['#1f77b4','#2ca02c','#d62728','#9467bd','#8c564b','y','#ff7f0e']
    nf      = len(inlierflts) # Number of filters with inliers
    fig, ax = plt.subplots(nf, figsize=(12,4*(nf**.75+1)), sharex=True)
    if nf == 1:
        ax  = [ax]
        
    for i,f in enumerate(inlierflts):
        sel    = cat['fltr'] == f
        ax[i].scatter(cat['ph'][sel],cat['mag'][sel],c=colors[f])
        ax[i].scatter(cat['ph'][sel]+1,cat['mag'][sel],c=colors[f])
        tmpmag = np.tile(fitter.tmps.columns[tmpind]*pars[1]*fitter.ampratio[f]+pars[2:][f],2)
        tmpph  = np.tile(fitter.tmps['PH'],2)+([0]*len(fitter.tmps['PH'])+[1]*len(fitter.tmps['PH']))
        ax[i].plot(tmpph,tmpmag,c='k')
        xsel   = sel*(~cat['inlier'])
        ax[i].scatter(cat['ph'][xsel],cat['mag'][xsel],c='k',marker='x')
        ax[i].scatter(cat['ph'][xsel]+1,cat['mag'][xsel],c='k',marker='x')
        ax[i].invert_yaxis()
        ax[i].set_ylabel(fitter.fltnames[f], fontsize=20)
    
    ax[-1].set_xlabel('Phase', fontsize=20)
    ax[0].set_title("Object: {}    Period: {:.3f} d    Type: {}".format(
                                        objname,p,fitter.tmps.colnames[tmpind]), fontsize=22)
    
    path = dirpath  + '{}_plot.png'.format(objname)
    fig.savefig(path)
    if verbose:
        print('Saved to',path)
        print('Save parameters')
    plt.close(fig)
    
    # save parameters and results
    res = Table([[objname]],names=['name'])
    res['period'] = p
    res['t0']     = pars[0]
    res['r amp']  = pars[1]
    for i in range(2,len(pars)):
        f = fitter.fltnames[i-2]
        res['{} mag'.format(f)] = pars[i]
    res['chi2']   = chi2
    res['redchi2']= redchi2
    res['template']= fitter.tmps.colnames[tmpind]
    res['t0 err']     = err[0]
    res['amp err']  = err[1]
    for i in range(2,len(err)):
        f = fitter.fltnames[i-2]
        res['{} mag err'.format(f)] = err[i]
    res['Ndat']      = len(cat)
    res['N inliers'] = sum(cat['inlier'])
    for i in range(len(fitter.fltnames)):
        f = fitter.fltnames[i]
        res['N {}'.format(f)] = sum(cat['fltr'][cat['inlier']]==i)
    
    path = dirpath  + '{}_res.fits'.format(objname)

    res.write(path,format='fits',overwrite=True)
    
    if verbose:
        print('Saved to',path)
        print('end')
    return
