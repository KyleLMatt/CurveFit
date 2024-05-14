import numpy as np
from astropy.timeseries import LombScargle, TimeSeries, LombScargleMultiband
from astropy.time import Time
import astropy.units as u
from .utils import *

# Data Lab
from dl import authClient as ac, queryClient as qc



class Variable:
    """
    This class gives the core functionality to look for
    and analyse variability inside of the NSC. It operates 
    on the basis of stars being objects.
    
    Parameters
    ----------
    objid: str
        Unique object ID inside of the NSC catalog.
    period: float, optional
        Period of the variable star, if known.
    variclass: str, optional
        Variability class (e.g. RRLab, Classical Cepheid, etc.). No functionality implemented for it at the moment.
    timeseries: TimeSeries object, optional
        Object with data including times and magnitudes. Additional data is also possible using the Astropy TimeSeries functionality.
        If not given, the data will be retrieved automatically based on the given Object ID and data release.
    datarelease: str, optional
        Data release from which the data comes from, or where it should be taken from. Default is DR2.
    """
    
    def __init__(self, objid, period=None, variclass=None, timeseries=None, datarelease="dr2"):
        """ 
        Initialize the Star object. If data for the time series is not given, then
        get it from Datalab.
        """
        
        accepted_variclass = ['Cepheid', 'RRLyrae', 'RRLab', 'RRLc', 'Mira', 'LPV']
        
        self.objid = objid
        self.period = period
        self.variclass = variclass
        self.datarelease = datarelease
        if timeseries is None:
            self.timeseries = self.get_timeseries_data()
        else:
            self.timeseries = timeseries
        

    def get_timeseries_data(self, datarelease=None, datalab_token=None):
        """
        For a given Object ID, return time series data.
        
        Parameters
        ----------
        objid : str, optional
            ID of the star. Keep in mind it does not carry over different data releases.
        datarelease: str, optional
            Data release from which to get the data from. Default is the class default, currently "dr2".
        """
        
        if datarelease==None: datarelease = self.datarelease
        
        if datarelease=='dr2' or datarelease=='dr1' or datarelease=='DR2' or datarelease=='DR1':
            mag_name = 'mag_auto'
            magerr_name = 'magerr_auto'
        else:
            mag_name = 'mag'
            magerr_name = 'magerr'
        
        query = f"""SELECT m.mjd,m.{mag_name},m.{magerr_name},m.filter,e.exptime 
                FROM nsc_dr2.meas AS m JOIN nsc_dr2.exposure as e ON m.exposure=e.exposure
                WHERE m.objectid='{self.objid}'"""
        # f"SELECT mjd,{mag_name},{magerr_name},filter FROM nsc_{datarelease}.meas WHERE objectid='{self.objid}'"
        table_res = qc.query(query,fmt='table')
        timeseries_obj = TimeSeries(data=[table_res[mag_name],table_res[magerr_name],table_res['filter'],table_res['exptime']],time=Time(table_res['mjd'], format='mjd'))
        
        return timeseries_obj
         
    def issp(self):
        """
        Function to determine whether a variable has or could have
        a short or long period. Mostly for internal use, to establish
        period search ranges.
        
        It will use the actual period value if available, if not the variability class
        is used. If nothing is available, a short period is assumed.
        
        Returns
        -------
        shortperiod: bool
            True if the variable has a short period, False if not.
        """
        
        if self.period!=None:
            if self.period <= 10*u.day:
                return True
            else:
                return False
            
        if self.variclass in ['Cepheid', 'RRLyrae', 'RRLab', 'RRLc']:
            return True
        elif self.variclass=='Mira' or self.variclass=='LPV':
            return False
        else:
            return True
        
    def franges(self):
        '''
        Function to determine frequency ranges if they are not provided.
        
        Returns
        -------
        min_frequency: Time
            Minimum frequency, in 1/days (default).
        max_frequency: Time
            Maximum frequency, in 1/days (default).
        '''
        
        if self.issp():
            minimum_frequency = 1/(50*u.day)
            maximum_frequency = 1/(0.04*u.day)
        else:
            minimum_frequency = 1/(1000*u.day)
            maximum_frequency = 1/(50*u.day)
            
        return minimum_frequency, maximum_frequency
        
    
    def frequency_array(self, nbins=100, minimum_frequency=None, maximum_frequency=None):
        '''
        Function to define the frequency array used for calculating the 
        different periodograms. It looks for information on whether the 
        suspected variable has a long or short period to define the 
        frequency ranges.
        
        Returns:
        --------
        farray: Time array
            Array of frequencies (in 1/days, default).
        
        '''
        
        if maximum_frequency is None or minimum_frequency is None:
            min_freq, max_freq = franges()
        
        farray = np.linspace(min_freq, max_freq, nbins)
        
        return farray
        
        
    def ls_mb_periodogram(self, method='flexible', normalization='standard', minimum_frequency=None, maximum_frequency=None):
        """
        Calculates a multi-band Lomb-Scargle periodogram 
        based on the data stored in timeseries.
        
        Returns
        -------
        frequency : ndarray
            Frequencies for the periodogram.
        power : ndarray
            Power for the corresponding frequencies.
        """

        if maximum_frequency is None or minimum_frequency is None:
            minimum_frequency, maximum_frequency = self.franges()
        
        if self.datarelease=='dr1' or self.datarelease=='dr2':
            mags = self.timeseries['mag_auto']
            mags_errs = self.timeseries['magerr_auto']
        else:
            mags = self.timeseries['mag']
            mags_errs = self.timeseries['magerr']
        
        frequency, power = LombScargleMultiband(self.timeseries['time'],mags,self.timeseries['filter'],dy=mags_errs).autopower(
            method=method,
            normalization=normalization,
            minimum_frequency=minimum_frequency,
            maximum_frequency=maximum_frequency
            )
        
        return frequency, power
    
    def ls_periodogram(self, band=None, method='flexible', normalization='standard', minimum_frequency=None, maximum_frequency=None):
        """
        Calculates a Lomb-Scargle periodogram for a single
        band, based on the data stored in timeseries.
        
        Returns
        -------
        frequency : ndarray
            Frequencies for the periodogram.
        power : ndarray
            Power for the corresponding frequencies.
        """
        
        if band==None:
            band = most_frequent(self.timeseries['filter'])
        
        selection = self.timeseries[self.timeseries['filter']==band]
        
        if maximum_frequency is None or minimum_frequency is None:
            minimum_frequency, maximum_frequency = self.franges()
        
        if self.datarelease=='dr1' or self.datarelease=='dr2':
            mags = selection['mag_auto']
            mags_errs = selection['magerr_auto']
        else:
            mags = selection['mag']
            mags_errs = selection['magerr']
            
        frequency, power = LombScargle(selection['time'],mags,dy=mags_errs).autopower(
            method=method,
            normalization=normalization,
            minimum_frequency=minimum_frequency,
            maximum_frequency=maximum_frequency
            )
        
        return frequency, power

    
    def lk_periodogram(self, minimum_frequency=None, maximum_frequency=None):
        '''
        Calculates a Lafler-Kinman periodogram for a single band, 
        based on the data stored in timeseries.
        
        Returns
        -------
        frequency : ndarray
            Frequencies for the periodogram.
        power : ndarray
            Power for the corresponding frequencies.
        '''
        
        if maximum_frequency is None or minimum_frequency is None:
            minimum_frequency, maximum_frequency = self.franges()
        farray = frequency_array(minimum_frequency, maximum_frequency)
        parray = 1./farray
        
        tobs = self.timeseries['time']
        if self.datarelease=='dr1' or self.datarelease=='dr2':
            mags = self.timeseries['mag_auto']
        else:
            mags = self.timeseries['mag']
        
        t0 = np.min(tobs)
        tt = tobs - t0
        theta = np.zeros_like(parray)
        mmplus_km = np.zeros_like(mags)
        avm_km = np.sum(mags)/len(mags) 
        denom_km = np.sum( (mags-avm_km)**2 )
        m = len(parray)
        for k in range(m):
            period = parray[k]
            phi = tt / period
            #nphi = np.fix(phi)          #KJM: literal but slower
            nphi = phi.astype(np.int64)  #KJM: ~25% faster
            phi = phi - nphi
            ss = np.argsort(phi)  #KJM: BEWARE the IDL sort gotcha!
            mm  = mags[ss]
            #mmplus = np.append(mm[1:], mm[0])   #KJM: literal but slower
            #numer = np.sum( (mmplus - mm)**2 )  #KJM: uses mmplus
            mmplus_km[:-1] = mm[1:]  #KJM: Don't use np.append within loops!
            mmplus_km[-1] = mm[0]    #KJM: Don't use np.append within loops!
            #assert np.allclose(mmplus,mmplus_km)  #KJM: NUM? ME VEXO?
            numer = np.sum( (mmplus_km - mm)**2 )  #KJM: uses mmplus_km
            #KJM: using mmplus_km is ~24% faster
            theta[k] = numer/denom_km
        
        return theta, phi
    
    
    def get_folded_ts(self, period=None):
        """
        Folds time series data according to the period of the star.
        
        Returns
        -------
        phase : ndarray
            
        """
        
        if period==None: period = self.period
        
        try:
            phase = phase_fold(self.timeseries['mjd'], period)
        except:
            print('The star has no calculated period.')
            return None
            
        return phase
    
    def get_period(self, frequency, power):
        """
        Returns the most likely period for the given periodogram.
        Looks for the absolute maximum in the periodogram. It also 
        sets the period attribute for the Class.
        
        Parameters
        ----------
        frequency: array-like
            Frequencies evaluated in the periodogram.
        power: array-like
            Power at each frequency. Must have the same shape.
            
        Returns
        -------
        period: float
            Most likely period.
        error: float
            Error based on the precision in the given frequencies.
        """
        
        bestind = np.argmax(power)
        period = 1./frequency[bestind]
        
        self.period = period # The period is an attribute of Variable, set it here too.
        
        error_f = np.abs(frequency[bestind-1]-frequency[bestind+1])
        error = error_f / frequency[bestind]**2
        
        return period, error
