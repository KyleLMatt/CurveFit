#!/usr/bin/env python

"""SAMPLER.PY - Variable star sampler

"""

from __future__ import print_function

__authors__ = 'David Nidever <dnidever@montana.edu>'
__version__ = '20220320'  # yyyymmdd

import time
import numpy as np
from dlnpyutils import utils as dln
from astropy.table import Table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import corner

def solveone(data,template,ampratios,bandindex,period,offset,totwtdict,totwtydict):

    ndata = len(data)
        
    # Calculate phase for each data point
    phase = (data['jd']/period + offset) % 1
            
    # Calculate template values for this set of period and phase
    tmpl = np.interp(phase,template['phase'],template['mag'])
            
    # -- Find best fitting values for linear parameters ---
    # Calculate amplitude
    # term1 = Sum of XY
    # term2 = Sum of X * Y / W 
    # term3 = Sum of X^2
    # term4 = Sum of X * X / W
    # amplitude = (term1 - term2)/(term3 - term4)
    term1 = 0.0
    term2 = 0.0
    term3 = 0.0
    term4 = 0.0
    totwtxdict = {}
    for b in bandindex.keys():
        ind = bandindex[b]
        totwtx1 = np.sum(data['wt'][ind] * tmpl[ind]*ampratios[b])
        totwtxdict[b] = totwtx1
        totwtx2 = np.sum(data['wt'][ind] * (tmpl[ind]*ampratios[b])**2)
        totwtxy = np.sum(data['wt'][ind] * tmpl[ind]*ampratios[b] * data['mag'][ind])      
        term1 += totwtxy
        term2 += totwtx1 * totwtydict[b] / totwtdict[b]
        term3 += totwtx2
        term4 += totwtx1**2 / totwtdict[b]
    amplitude = (term1-term2)/(term3-term4)
            
    # Calculate best mean magnitudes
    # mean mag = (Y - amplitude * X)/W
    meanmag = {}
    for b in bandindex.keys():
        meanmag1 = (totwtydict[b] - amplitude * totwtxdict[b])/totwtdict[b]
        meanmag[b] = meanmag1
            
    # Calculate likelihood/chisq
    model = np.zeros(ndata,float)
    resid = np.zeros(ndata,float)
    wtresid = np.zeros(ndata,float)        
    for b in bandindex.keys():
        ind = bandindex[b]          
        model1 = tmpl[ind]*ampratios[b]*amplitude+meanmag[b]
        model[ind] = model1
        resid[ind] = data['mag'][ind]-model1
        wtresid[ind] = resid[ind]**2 * data['wt'][ind]
    lnlkhood = -0.5*np.sum(wtresid + np.log(2*np.pi*data['err']**2))

    return amplitude,meanmag,model,lnlkhood


def sampler(catalog,template,pmin=0.1,pmax=None,ampratios=None,minerror=0.02,
            minsample=128,npoints=200000,plotbase='sampler'):
    """

    catalog : table
       Catalog of data points, just have mag, err, jd, band
    template : table
       Template information.
    pmin : float, optional
       Minimum period to search in days.  Default is 0.1 days.
    pmax : float, optional
       Maximum period to search in days.  Default is 2 x time baseline.
    ampratios : dict, optional
       Amplitude ratios.  Keys should be the unique band names
         and values should be the amplitue ratios.
         If this is not input, then a ratio of 1.0 is used.
    minerror : float, optional
       Minimum error to use.  Default is 0.02.
    minsample : int, optional
       Mininum number of samples to return.  Default is 128.
    npoints : int, optional
       Number of points to use per loop.  Default is 200,000.
    plotbase : str, optional
       Base name for output plots.  Default is "sampler".
    
    """

    t0 = time.time()

    # Create the sampling for Period (pmin to pmax) and phase offset (0-1)

    # Internal catalog
    data = Table(catalog).copy()
    data['wt'] = 1/np.maximum(data['err'],minerror)**2
    
    # Only keep bands with 2+ observations
    uband = np.unique(data['band'])
    badind = np.array([],int)
    for i,b in enumerate(uband):
        ind, = np.where(data['band']==b)
        if len(ind)<2:
            print('band '+str(b)+' only has '+str(len(ind))+' observations.  Not using')
            badind = np.hstack((badind,ind))
    if len(badind)>0:
        data.remove_rows(badind)
    ndata = len(data)

    print(str(ndata)+' data points')
    print('time baselines = %.2f' % (np.max(data['jd'])-np.min(data['jd'])))
    
    # Get band index
    uband = np.unique(data['band'])
    nband = len(uband)
    bandindex = {}
    for i,b in enumerate(uband):
        ind, = np.where(data['band']==b)
        bandindex[b] = ind            

    print(str(len(uband))+' bands = ',', '.join(np.char.array(uband).astype(str)))
        
    # No amplitude ratios input
    if ampratios is None:
        ampratios = {}
        for b in uband:
            ampratios[b] = 1.0

    # Period range
    if pmax is None:
        pmax = (np.max(data['jd'])-np.min(data['jd']))*2
    lgminp = np.log10(pmin)
    lgmaxp = np.log10(pmax)
    
    print('Pmin = %.3f' % pmin)
    print('Pmax = %.3f' % pmax)    
    
    # Pre-calculate some terms that are constant
    totwtdict = {}
    totwtydict = {}
    for b in uband:
        ind = bandindex[b]
        totwtdict[b] = np.sum(data['wt'][ind])
        totwtydict[b] = np.sum(data['wt'][ind] * data['mag'][ind])
    
    # Loop until we have enough samples
    nsamples = 0
    samplelist = []
    count = 0
    dtt = [('period',float),('offset',float),('amplitude',float),('lnlikelihood',float),('lnprob',float)]
    for b in uband:
        dtt += [('mag'+str(b),float)]
    trials = None
    while (nsamples<minsample):
    
        # Uniformly sample from log(pmin) to log(pmax)
        period = np.random.rand(npoints)*(lgmaxp-lgminp)+lgminp    
        period = 10**period
        # Uniformly sample from 0 to 1
        offset = np.random.rand(npoints)


        # Get phase and template points
        phase = (data['jd'].reshape(-1,1)/period.reshape(1,-1) + offset.reshape(1,-1)) % 1
        tmpl = np.interp(phase.ravel(),template['phase'],template['mag'])
        tmpl = tmpl.reshape(ndata,npoints)
            
        # -- Find best fitting values for linear parameters ---
        # Calculate amplitude
        # term1 = Sum of XY
        # term2 = Sum of X * Y / W 
        # term3 = Sum of X^2
        # term4 = Sum of X * X / W
        # amplitude = (term1 - term2)/(term3 - term4)
        term1,term2,term3,term4 = 0,0,0,0
        totwtxdict = {}
        for b in uband:
            ind = bandindex[b]
            totwtx1 = np.sum(data['wt'][ind].reshape(-1,1) * tmpl[ind,:]*ampratios[b],axis=0)
            totwtxdict[b] = totwtx1
            totwtx2 = np.sum(data['wt'][ind].reshape(-1,1) * (tmpl[ind,:]*ampratios[b])**2,axis=0)
            totwtxy = np.sum(data['wt'][ind].reshape(-1,1) * tmpl[ind,:]*ampratios[b] * data['mag'][ind].reshape(-1,1),axis=0)      
            term1 += totwtxy
            term2 += totwtx1 * totwtydict[b] / totwtdict[b]
            term3 += totwtx2
            term4 += totwtx1**2 / totwtdict[b]
        amplitude = (term1-term2)/(term3-term4)
    
        # Calculate best mean magnitudes
        # mean mag = (Y - amplitude * X)/W
        meanmag = {}
        for b in uband:
            meanmag1 = (totwtydict[b] - amplitude * totwtxdict[b])/totwtdict[b]
            meanmag[b] = meanmag1
            
        # Calculate likelihood/chisq
        model = np.zeros((ndata,npoints),float)
        resid = np.zeros((ndata,npoints),float)
        wtresid = np.zeros((ndata,npoints),float)        
        for b in uband:
            ind = bandindex[b]
            model1 = tmpl[ind,:]*ampratios[b]*amplitude+meanmag[b]
            model[ind,:] = model1
            resid[ind,:] = data['mag'][ind].reshape(-1,1)-model1
            wtresid[ind,:] = resid[ind,:]**2 * data['wt'][ind].reshape(-1,1)
        lnlikelihood = -0.5*np.sum(wtresid,axis=0)
        lnlikelihood += -0.5*np.sum(np.log(2*np.pi*data['err']**2))

        # Calculate ln probability = ln prior + ln likelihood
        # use flat prior, divide by area
        lnprior = np.ones(npoints,float) + np.log(1/(1.0*(lgmaxp-lgminp)))
        lnprob = lnprior + lnlikelihood

        # Save the information
        trials1 = np.zeros(npoints,dtype=dtt)
        trials1['period'] = period
        trials1['offset'] = offset
        trials1['amplitude'] = amplitude
        for k in meanmag.keys():
            trials1['mag'+str(k)] = meanmag[k]
        trials1['lnlikelihood'] = lnlikelihood
        trials1['lnprob'] = lnprob        
        if trials is None:
            trials = trials1
        else:
            trials = np.hstack((trials,trials1))
        
        # REJECT NEGATIVE AMPLITUDES??
        
        # Rejection sampling
        draw = np.random.rand(npoints)
        #maxprob = np.max(np.exp(lnprob))
        #ind, = np.where(draw < np.exp(lnprob)/maxprob)
        ind, = np.where(draw < np.exp(lnprob))  
        if len(ind)>0:
            for i in ind:
                samp = {'period':period[i],'offset':offset[i],'amplitude':amplitude[i]}
                for k in meanmag.keys():
                    samp[k] = meanmag[k][i]
                samp['lnlikelihood'] = lnlikelihood[i]
                samp['lnprob'] = lnprob[i]
                samplelist.append(samp)
            nsamples += len(ind)
            
        print(count+1,nsamples)
        count += 1
        
    # Convert sample list to table
    dt = [('period',float),('offset',float),('amplitude',float)]
    for k in meanmag.keys():
        dt += [('mag'+str(k),float)]
    dt += [('lnlikelihood',float),('lnprob',float)]
    samples = np.zeros(len(samplelist),dtype=dt)
    for i,samp in enumerate(samplelist):
        samples['period'][i] = samp['period']
        samples['offset'][i] = samp['offset']
        samples['amplitude'][i] = samp['amplitude']
        samples['lnlikelihood'][i] = samp['lnlikelihood']
        samples['lnprob'][i] = samp['lnprob']
        for k in meanmag.keys():
            samples['mag'+str(k)][i] = samp[k]

    # Convert to astropy tables
    samples = Table(samples)
    trials = Table(trials)
            
    best = np.argmax(trials['lnprob'])
    bestperiod = trials['period'][best]
    bestoffset = trials['offset'][best]
    bestlnprob = trials['lnprob'][best]    
    print('Best period = ',bestperiod)
    print('Best offset = ',bestoffset)    
    print('Best lnprob = ',bestlnprob)

    ntrials = npoints*count
    print('ntrials = ',ntrials)
    
    print('dt=',time.time()-t0)

    # Make plots
    matplotlib.use('Agg')
    fig,ax = plt.subplots(2,1)
    fig.set_figheight(10)
    fig.set_figwidth(10)
    # 2D density map
    im,b,c,d = stats.binned_statistic_2d(trials['offset'],np.log10(trials['period']),trials['lnprob'],statistic='mean',bins=(250,250))
    z1 = ax[0].imshow(im,aspect='auto',origin='lower',extent=(c[0],c[-1],b[0],b[-1]))
    ax[0].set_xlabel('log(Period)')
    ax[0].set_ylabel('Offset')
    plt.colorbar(z1,ax=ax[0],label='Mean ln(Prob)')
    # Period histogram
    hist,a,b = stats.binned_statistic(np.log10(trials['period']),trials['lnprob'],statistic='mean',bins=1000)
    ax[1].plot(a[0:-1],hist)
    ax[1].set_xlabel('log(Period)')
    ax[1].set_ylabel('Mean ln(Prob)')
    fig.savefig(plotbase+'_trials.png',bbox_inches='tight')
    plt.close(fig)
    print('Saving to '+plotbase+'_trials.png')

    sampdata = np.zeros((len(samples),3),float)
    sampdata[:,0] = samples['period']
    sampdata[:,1] = samples['offset']
    sampdata[:,2] = samples['amplitude']
    samplabels = ['Period','Offset','Amplitude']    
    fig = corner.corner(sampdata, labels=samplabels)
    plt.savefig(plotbase+'_corner.png',bbox_inches='tight')
    plt.close(fig)
    print('Corner plot saved to '+plotbase+'_corner.png')


    
    return samples, trials
