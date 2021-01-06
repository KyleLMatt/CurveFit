#!/usr/bin/env python

from glob import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#%matplotlib inline
import pandas as pd
from scipy.interpolate import interp1d
from scipy.signal import gaussian, convolve
from statistics import mean, median
from astropy import stats
from scipy.optimize import curve_fit, least_squares
import collections
import os

import utils
from dl import queryClient as qc
from argparse import ArgumentParser
import time
import socket
from datetime import datetime

def tempdir():
    """ Return the template directory."""
    fil = os.path.abspath(__file__)
    codedir = os.path.dirname(fil)
    datadir = codedir+'/templets/'
    return datadir

def get_data(df,objname,outdir='results'):
    order = ['u','g','r','i','z']
    best_periods = []
    crv=[]
    fltrs=[]
    for f in order:
        selfltr = (df['filter'] == f)
        selfwhm = (df['fwhm'] <= 4.0)
        sel = selfltr & selfwhm
        t = df['mjd'][sel].values
        y = df['mag_auto'][sel].values
        dy = df['magerr_auto'][sel].values
        if len(t) < 25:
            continue

        pout = get_ls_period(t,y,objname=objname+'_'+f,outdir=outdir)
        best_periods.append(pout)
        crvi = np.vstack((t,y,dy)).T
        crv.append(crvi[np.argsort(crvi[:,0])])
        fltrs.append(f)
    period = 0
    for p in best_periods:
        period += p/len(best_periods)
    return crv, period, fltrs

def get_tmps(fltrs):
    tmps=[]
    typs =[]
    names=[]
    templatedir = tempdir()
    #print('templatedir = '+templatedir)
    for fltr in fltrs:
        typ = []
        templets = glob(templatedir+'/*{}.dat'.format(fltr))
        #templets = glob('templets/*{}.dat'.format(fltr))
        tmp = np.zeros((len(templets),501,2))
        for i in range(len(templets)):
            tmp[i] = np.concatenate((np.array([[0,0]]),
                                     np.array(pd.read_csv(templets[i],sep=' ')),
                                     np.array([[1,0]])))
            #adjust if filepath to templets changes
            if len(os.path.basename(templets[i]))==8:
                typ.append('RRab')
            elif len(os.path.basename(templets[i]))==6:
                typ.append('RRc')
        typs.append(typ)
        names.append(templets)
        tmps.append(tmp)
    return tmps, names, typs

def double_tmps(tmps):
    tmps2=[]
    for f in range(len(tmps)):
        tmps2.append(np.tile(tmps[f],(2,1)))
        tmps2[f][:,int(len(tmps2[f][0])/2):,0] += 1
    return tmps2

def plot_periodogram(period,power,best_period=None,objname='',ax=None,outdir='results'):
   
    if ax is None:
        fig, ax = plt.subplots(figsize=(10,7))
        
    ax.plot(period,power,lw=0.1)
    ax.set_xlabel('period (days)')
    ax.set_ylabel('relative power')
    ax.set_title(objname)
    
    if best_period is not None:
        ax.axvline(best_period,color='r');
        ax.text(0.03,0.93,'period = {:.3f} days'.format(best_period),transform=ax.transAxes,color='r')
    fig.savefig(outdir+'/{}_periodogram.png'.format(objname))
    plt.close(fig)

def get_ls_period(t,y,min_freq=1./1.,max_freq=1./0.1,objname='_',outdir='results'):
    """Use Lomb-Scargle periodogram to get an estimate on period"""
    
    ls = stats.LombScargle(t, y)
    frequency, power = ls.autopower(minimum_frequency=min_freq,maximum_frequency=max_freq)
    period = 1./frequency # period is the inverse of frequency
    
    best_period = period[np.argmax(power)]
    
    plot_periodogram(period,power,best_period,objname=objname,outdir=outdir)
    return best_period

def get_pinit(crv,period):
    pinit = ()
    for ltcrv in crv:
        pinit += ((0.0,max(ltcrv[:,1])-min(ltcrv[:,1]),0.0),)
    pinit += (period,)
    return pinit

def update_pinit(pars,period):
    pinit = ()
    for i in range(len(pars)):
        pinit += (tuple(pars[i,:-1]),)
    pinit += (period,)
    return pinit

def RemoveOutliers(crv,tmps,pars,period):
    n = pars[:,-1].astype(int)
    crv_in = []
    for i in range(len(crv)):
        f = interp1d(tmps[i][n[i],:,0],tmps[i][n[i],:,1]*pars[i,1]+pars[i,2])
        phase = (crv[i][:,0]/period-pars[i,0]) %1
        dif = abs(crv[i][:,1]-f(phase))
        crv_in.append(crv[i][dif<utils.mad(dif)*5])
    return crv_in

def double_period(crv,pars,period):
    crv2 = []
    for i in range(len(crv)):
        crv2.append(crv[i].copy())
        crv2[i][:,1] -= pars[i,2]
        
        crv2[i][:,0] = (crv2[i][:,0]/period-pars[i,0])%1
        crv2[i] = np.tile(crv2[i].T,2).T
        crv2[i][int(len(crv2[i])/2):,0] += 1
        crv2[i] = crv2[i][crv2[i][:,0].argsort()]
        
    return crv2

class tmpfitter:
    def __init__ (self, tmps):
        self.fltr=0
        self.n=0
        self.tmps=tmps

    def model(self, t, t0, amplitude, yoffset):
        # modify the template using peak-to-peak amplitude, yoffset
        # fold input times t by period, phase shift to match template
        xtemp = self.tmps[self.fltr][self.n,:,0]
        ytemp = self.tmps[self.fltr][self.n,:,1]*amplitude + yoffset
        ph = (t - t0) %1
        #print((ph[0],period,t0%1))
        #print((period,t0,amplitude,yoffset))
        # interpolate the modified template to the phase we want
        return interp1d(xtemp,ytemp)(ph)


def tmpfit(crv,tmps,pinit,w=.1,steps=21,n=1):
    fitter = tmpfitter(tmps)
    
    lsteps = int(steps/2+.5)
    rsteps = steps - lsteps
    pl = np.linspace(pinit[-1]-w,pinit[-1],lsteps)
    pr = np.linspace(pinit[-1]+w,pinit[-1],rsteps,endpoint=False)
    plist = np.zeros(pl.size+pr.size)
    plist[0::2] = np.flip(pl)
    plist[1::2] = np.flip(pr)
    plist = plist[plist>0]
    
    pars = np.zeros((len(tmps),4))
    minsumx2 = 10**50
    minp = 0
    for p in plist:
        sumx2=0
        ppars=np.zeros((len(tmps),4))
        for f in range(len(tmps)):
            fitter.fltr = f
            phase = crv[f][:,0]/p%n #1 for one period, 2 for two periods
            minx2 = 10**50
            for i in range(len(tmps[f])):
                fitter.n = i
                try:
                    tpars, cov = curve_fit(fitter.model, phase, crv[f][:,1], 
                                          bounds = ((-.5,0,-50),(.5,10,50)),
                                          sigma=crv[f][:,2], p0=pinit[f], maxfev=500)
                except RuntimeError:
                    #print('Error: Curve_fit failed on templet={}-{}, p={:.4}'.format(f,i,p))
                    continue
                
                x2 = sum((fitter.model(phase,tpars[0],tpars[1],tpars[2])-crv[f][:,1])**2/crv[f][:,2]**2)
                if x2 < minx2:
                    ppars[f,:-1] = tpars
                    ppars[f,-1] = i
                    minx2 = x2
            
            sumx2 += minx2
            if sumx2 > minsumx2:
                break
        if sumx2 < minsumx2:
            minsumx2 = sumx2
            minp = p
            pars = ppars
    npoints=0
    for i in range(len(crv)):
        npoints += len(crv[i])
    return pars, minp, minsumx2, minsumx2/npoints

def fit_plot(objname,outdir='results'):
    # Create output directory if necessary
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # Get the measurements for this star
    star=qc.query(sql="""SELECT mag_auto,magerr_auto,measid,exposure,filter,fwhm,mjd
                     FROM nsc_dr2.meas
                     WHERE objectid='{:s}'""".format(objname),
              fmt='pandas',
              profile='db01')
    nmeas = len(star)
    print('  '+str(nmeas)+' measurements')
    #print(collections.Counter(star['filter']))
    crv,period,fltrs = get_data(star,objname,outdir=outdir)
    if len(fltrs) == 0:
        return
    tmps, tmpnames, typs = get_tmps(fltrs)
    
    pinit = get_pinit(crv,period)
    pars, p, chi2_1, rchi2_1 = tmpfit(crv,tmps,pinit,w=.1,steps=25)
    crv_in = RemoveOutliers(crv,tmps,pars,p)
    pinit = update_pinit(pars,p)
    pars_in,p_in,chi2,rchi2 = tmpfit(crv_in,tmps, pinit,w=.01,steps=25)
    
    crv2 = double_period(crv,pars_in,p_in)
    tmps2= double_tmps(tmps)
    n = pars[:,-1].astype(int)
    
    colors = []
    for f in fltrs:
        if f == 'r' or f == 'g':
            colors.append(f)
        else:
            colors.append('black')

    #Check if each filter is consistent with RR type (RRab or RRc)
    consistent = True
    for i in range(len(typs)):
        for j in range(i+1,len(typs)):
            if typs[i][n[i]] != typs[j][n[j]]:
                consistent = False
                break
        if not consistent:
            break
    if consistent:
        typ = typs[0][n[0]]
    else:
        typ = '???'
    # Make the lightcurve plot
    matplotlib.use('Agg')
    fig, ax = plt.subplots(len(fltrs), figsize=(10,7.5), sharex=True, sharey=True)
    if len(fltrs) == 1:
        ax = [ax]
    for i in range(len(fltrs)):
        crvmean = mean(crv2[i][:,1])
        ax[i].scatter(crv2[i][:,0],crv2[i][:,1]-crvmean,c=colors[i])
        ax[i].plot(tmps2[i][n[i],:,0],tmps2[i][n[i],:,1]*pars_in[i,1]-crvmean,c='black')
        ax[i].invert_yaxis()
        ax[i].set_ylabel(fltrs[i], fontsize=18)

    ax[-1].set_xlabel('Phase', fontsize=16)
    ax[0].set_title("Object: {}    Period: {:.3f} d    Type: {}".format(objname,p_in,typ), fontsize=20)
    fig.savefig(outdir+'/{}_lightcurve.png'.format(objname))
    # Write out the parameters
    ofile = open(outdir+'/'+str(objname)+"_parameters.csv",'w')
    ofile.write("{},{:d},{:.3f},{:.3f},{:.3f}\n".format(objname,nmeas,chi2,rchi2,p_in))
    print("  {},{:d},{:.3f},{:.3f},{:.3f}".format(objname,nmeas,chi2,rchi2,p_in))
    for i in range(len(fltrs)):
        ofile.write("{:s},{:.3f},{:.3f},{:.3f},{}\n".format(fltrs[i],pars_in[i][0],pars_in[i][1]/2,pars_in[i][2],os.path.basename(tmpnames[i][n[i]][9:])))
        print("  {:s},{:.3f},{:.3f},{:.3f},{}".format(fltrs[i],pars_in[i][0],pars_in[i][1]/2,pars_in[i][2],os.path.basename(tmpnames[i][n[i]][9:])))
    #ofile.write("---\n")
    plt.close(fig)


if __name__ == "__main__":
    parser = ArgumentParser(description='Fit RR Lyrae templates to NSC DR2 stars.')
    parser.add_argument('objectid', type=str, nargs='*', help='HEALPix pixel number')
    parser.add_argument('--outdir', type=str, default='', help='Output directory')

    args = parser.parse_args()

    t0 = time.time()
    hostname = socket.gethostname()
    host = hostname.split('.')[0]

    print('Running rrlfit.py on '+host)
    print(datetime.now())
    print(' ')

    # Inputs
    objectid = args.objectid
    print(str(len(objectid))+' objects')
    print(' ')

    basedir = '/net/dl2/dnidever/nsc/instcal/v3/rrl/' 
    # Loop over the objects
    for i,obj in enumerate(objectid):
        print(str(i+1)+' '+obj)
        t0 = time.time()
        hpix = int(obj.split('_')[0])
        hpixgrp = hpix//1000
        outdir = basedir+str(hpixgrp)
        print('  outdir='+outdir)
        try:
            fit_plot(obj,outdir=outdir)
        except:
            print('  problem')
        print('  dt = '+str(time.time()-t0)+' sec.')
        print(' ')
