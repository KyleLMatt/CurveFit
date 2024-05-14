import os
import time
import numpy as np
from astropy.io import fits
from astropy.table import Table
from scipy.optimize import curve_fit
# from numba import njit
import matplotlib.pyplot as plt


def filtersort(filts):
    """
    Sort the filter array.
    """
    # Get the measurements in each band
    findx = dln.create_index(filts)
    ni,nv,nb = 0,0,0
    index = np.array([],int)
    iind, = np.where(findx['value']==1)
    if len(iind)>0:
        ni = findx['num'][iind[0]]
        ind = findx['index'][findx['lo'][iind[0]]:findx['hi'][iind[0]]+1]
        index = np.concatenate((index,ind))
    vind, = np.where(findx['value']==2)
    if len(vind)>0:
        nv = findx['num'][vind[0]]
        ind = findx['index'][findx['lo'][vind[0]]:findx['hi'][vind[0]]+1]
        index = np.concatenate((index,ind))        
    bind, = np.where(findx['value']==3)    
    if len(bind)>0:
        nb = findx['num'][bind[0]]
        ind = findx['index'][findx['lo'][bind[0]]:findx['hi'][bind[0]]+1]
        index = np.concatenate((index,ind))        
    # Breaks
    break1 = ni
    break2 = ni+nv
    nmeas = [ni,nv,nb]
    return index,break1,break2,nmeas

def lightcurve(t,filts,pars,tempdict,doderiv=False):
    """
    Create Cepheid lightcurve.

    Parameters
    ----------
    t : numpy array
       Time array in days.
    filts : numpy array
       Filters IDs. 1: I-band, 2: V-band, 3: B-band
    pars : numpy array
       Parameters.
         Period (days)  
         phase shift (unitless, 0-1)
         mag0_i -- average I-band magnitude
         mag0_v -- average V-band magnitude
         mag0_b -- average B-band magnitude
    tempdict : dictionary
       Dictionary with template information.
    doderiv : bool, optional
       Return phase derivative as well. Default is False.

    Returns
    -------
    result : numpy array
       The model magnitudes.
    pars : numpy array

    Example
    -------

    result = lightcurve(t,filts,pars,tempdict)

    """

    if len(t) != len(filts):
        raise ValueError('t and filts must have the same length')
    
    # Sort the filter information
    findex,break1,break2,nmeas = filtersort(filts)
    ni,nv,nb = nmeas

    # Put together the parameters with the breaks
    # p[0] = m_0_i -- average I-band magnitude
    # p[1] = m_0_v -- average V-band magnitude
    # p[2] = m_0_b -- average B-band mag
    # p[3] = Period (days)  
    # p[4] = phase shift (days)
    # p[5] = t_break (first element that has V-band data)
    # p[6] = t_break2 (first element that has B-band data)
    pars_ord = np.zeros(7,float)
    pars_ord[0:3] = pars[2:5]      # mag0_i, mag0_v, mag0_b (magnitude)
    pars_ord[3] = pars[0]          # period (days)
    pars_ord[4] = pars[1]*pars[0]  # phase offset in days
    pars_ord[5] = break1
    pars_ord[6] = break2

    # Filter-order the times as well
    tord = t[findex]
        
    # Run the normal code with filter-ordered measurements
    result_ord,pcoeff,phderiv_ord = _lightcurve(tord, pars_ord, tempdict, doderiv)
    
    # Reorder
    result = np.zeros(len(t),float)
    result[findex] = result_ord

    # Have option to return the derivative
    # period, phase offset, mag0's
    # period - need to run a second period
    # phase - used in nb_fouriercomp()
    # mag0's  these are just 1
    if doderiv:
        phderiv = np.zeros(len(t),float)
        phderiv[findex] = phderiv_ord
        # phderiv is in mag/days, convert to unitless phase units
        # by multiplying by the period
        phderiv *= pars[0]
        return result,phderiv
        
    return result

class Cepheid():
    """
    Cepheid object to model and fit cepheid light curves.
    """
    
    def __init__(self):
        """ Initialize the Cepheid object."""
        # Load the templates
        # convert to dictionaries, faster retrieval        
        tempdir = datadir()
        short_struc = fits.getdata(tempdir+'short_period.fits')
        short_struc = dict(zip(np.char.array(short_struc.names).lower(),short_struc[0]))
        long_struc = fits.getdata(tempdir+'long_period.fits')
        long_struc = dict(zip(np.char.array(long_struc.names).lower(),long_struc[0]))
        self._data = [short_struc, long_struc]

    def gettemplate(self,period):
        """ Get the template for the right period."""
        if period < 10:
            return self._data[0]
        else:
            return self._data[1]
    
    def __call__(self,xdata=None,pars=None):
        """
        Create a template for the given input values.

        Parameters
        ----------
        xdata : list or tuple
           List/tuple containing the time and filter information.
        pars : numpy array or list
           The parameters. [period, phase, mag0_i, mag0_v, mag0_b]

        Returns
        -------
        m : numpy array
           The model magnitudes.

        Example
        -------

        mags = ceph((t,filt),pars)

        """

        if xdata is None and pars is None:
            period = 15.0
            pars = [period, 0.10, 5.0, 0.0, 0.0]
            n = 100
            t = np.linspace(0,period,n)
            filt = np.ones(n,int)
            xdata = (t,filt)
        if xdata is not None and pars is None:
            raise ValueError('Need pars')
        if xdata is None and pars is not None:
            period = pars[0]
            n = 100
            t = np.linspace(0,period,n)
            filt = np.ones(n,int)
            xdata = (t,filt)
        if len(xdata) != 2:
            raise ValueError('xdata must have t and filts')
        if len(pars) < 5:
            origpars = np.array(pars).copy()
            pars = np.zeros(5,float)
            pars[:len(origpars)] = origpars

        return self.model(xdata,*pars)
        
    def model(self,xdata,*pars,doderiv=False):
        """
        Generate a model for a given input parameters.  This is used
        with fit() and or curve_fit().

        Parameters
        ----------
        xdata : numpy array, list or tuple
           This should contain the time and filter information.
        pars : args
           The parameters.  period, phase, magnitude offsets.
        doderiv : bool, optional
           Return the derivative as well.  Default is False.

        Returns
        -------
        result : numpy array
           The model magnitudes.
        phderiv : numpy array
           The phase derivative.  Only if doderiv=True.

        Example
        -------

        result = ceph.model(xdata,*pars)

        """
        
        tempdict = self.gettemplate(pars[0])
        t,filts = xdata
        findex,break1,break2,nmeas = filtersort(filts)        
        result = lightcurve(t,filts,pars,tempdict,doderiv)
        return result

    def jac(self,xdata,*pars,retmodel=False):
        """
        Return Jacobian matrix, the derivatives for each data point.

        Parameters
        ----------
        xdata : numpy array, list or tuple
           This should contain the time and filter information.
        pars : args
           The parameters.  period, phase, magnitude offsets.
        retmodel : boolean, optional
            Return the model as well.  Default is retmodel=False.

        Returns
        -------
        jac : numpy array
           Jacobian matrix of partial derivatives [N,Npars].
        model : numpy array
           The model magnitudes (if retmodel=True).

        Example
        -------

        jac = ceph.jac(xdata,*pars)

        """

        npars = len(pars)
        t,filts = xdata
        nt = len(t)

        # Run model and get the phase derivative
        model,phderiv = self.model(xdata,*pars,doderiv=True)
        # Get second model with period offset
        dperiod = 0.0001
        pars2 = np.array(pars).copy()
        pars2[0] += dperiod
        model2 = self.model(xdata,*pars2)
        perderiv = (model2-model)/dperiod
        
        # Initialize jacobian matrix
        jac = np.zeros((nt,npars),float)
        jac[:,0] = perderiv
        jac[:,1] = phderiv
        jac[:,2] = 1  # mag0_i
        jac[:,3] = 1  # mag0_v
        jac[:,4] = 1  # mag0_b

        if retmodel:
            return jac,model
        else:
            return jac
    
    def fit(self,data):
        """
        Fit Cepheid model to the data.

        Parameters
        ----------
        data : numpy array, list or tuple
          The data to fit.  This should be time, magnitude, filter, and uncertainties.

        Returns
        -------
        pars : numpy array
           Best-fitting parameters.
        pcov : numpy array
           Covariance matrix.

        Example
        -------

        pars,pcov = ceph.fit(data)

        """

        pass