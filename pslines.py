# Power spectra lines plotting tasks.
# Author: S. Munshi

import numpy as np
import matplotlib.pyplot as plt
import datetime
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.cosmology import Planck15 as cosmo
import astropy.constants as const

def dec_to_str(decimal_hours):
    if type(decimal_hours) == str:
        return decimal_hours
    else:
        hours = int(decimal_hours)
        minutes = int((decimal_hours - hours) * 60)
        seconds = int(((decimal_hours - hours) * 60 - minutes) * 60)
        time_obj = datetime.timedelta(hours=hours, minutes=minutes, seconds=seconds)
        time_str = str(time_obj)
        return time_str
    
def time_to_h(time_str):
    if type(time_str) != str:
        return time_str
    else:
        hh, mm, ss = map(float, time_str.split(':'))
        total_hours = hh + (mm / 60) + (ss / 3600)
        return total_hours

def utc_to_lst(utc, date, lat, lon):
    loc = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=0*u.m)
    time = date+'T'+dec_to_str(utc)
    t = Time(time, scale='utc', location=loc)
    return t.sidereal_time('mean').rad

def angular_distance(ra1, dec1, ra2, dec2):
    dlon = ra2 - ra1
    dlat = dec2 - dec1
    a = np.sin(dlat/2)**2+np.cos(dec1)*np.cos(dec2)*np.sin(dlon/2)**2
    angle= 2 * np.arctan(np.sqrt(a)/np.sqrt(1-a))
    return angle

class get_horizon:
    """
    Make and plot horizon lines for cylindrical power spectra.

    Arguments:
    freq_cen (float): Central frequency in Hz
    lat (float): Latitude in degree
    lon (float): Longitude in degree
    ra0 (float): RA of phase center in hour
    dec0 (float): Declination of phase center in degree
    """

    def __init__(self, freq_cen, lat, lon, ra0, dec0):
        self.freq_cen = freq_cen
        self.z = (1.4204057e9/freq_cen)-1
        self.lat = lat*np.pi/180
        self.lon = lon*np.pi/180
        self.ra0 = ra0*np.pi/12
        self.dec0 = dec0*np.pi/180
        self.date = None
        self.utc = None
        self.factor_cosmo = ((cosmo.comoving_transverse_distance(self.z)/((const.c*(1+self.z)**2)/(cosmo.H(self.z)*1.4204057e9*u.Hz))).to(u.Hz).value)/freq_cen
    
    def get_factor_lst(self, lst):
        factor = (1+np.sqrt(1-(np.sin(self.lat)*np.sin(self.dec0)+np.cos(self.lat)*np.cos(self.dec0)*np.cos(lst-self.ra0))**2))/(np.sin(self.lat)*np.sin(self.dec0)+np.cos(self.lat)*np.cos(self.dec0)*np.cos(lst-self.ra0))
        if factor < 0:
            raise ValueError('Phase Center goes below the horizon!')
        else:
            return factor
        
    def get_factor_max(self, N=100):
        utc_start = time_to_h(self.utc[0])
        utc_end = time_to_h(self.utc[1])
        if isinstance(self.date, tuple):
            if len(self.date) > 2 or (dec_to_str(self.utc[1]) > dec_to_str(self.utc[0])):
                raise ValueError('If observation spans more than 24 hours, use full synthesis line.')
            else:
                utc_range = np.linspace(utc_start, 24+utc_end, N)
                factors = []
                for i in range(N):
                    if utc_range[i] < 24:
                        lst = utc_to_lst(utc_range[i], self.date[0], self.lat, self.lon)
                    else:
                        lst = utc_to_lst(utc_range[i]-24, self.date[1], self.lat, self.lon)
                    factors.append(self.get_factor_lst(lst))
                return np.amax(factors)
        else:
            if (time_to_h(self.utc[1]) <= time_to_h(self.utc[0])):
                raise ValueError('End UTC must be greater than start UTC for an observation spanning a single day.')
            else:
                utc_range = np.linspace(utc_start, utc_end, N)
                factors = []
                for i in utc_range:
                    lst = utc_to_lst(i, self.date, self.lat, self.lon)
                    factors.append(self.get_factor_lst(lst))
                return np.amax(factors)
    
    def get_horizon_lst(self, kper_list):
        if isinstance(self.utc, tuple):
            if len(self.utc) != 2:
                raise ValueError('UTC should be either a single value or a tuple of two values corresponding to the start and end of the observation.')
            else:
                factor = self.get_factor_max()
        else:
            lst = utc_to_lst(self.utc, self.date, self.lat, self.lon)
            factor = self.get_factor_lst(lst)
        return factor*self.factor_cosmo*kper_list
    
    def get_horizon_full(self, kper_list):    
        factor = (1+np.sin(abs(self.lat)+abs(self.dec0)))/(np.sin(abs(self.dec0)-(np.pi/2-abs(self.lat))))
        return factor*self.factor_cosmo*kper_list
    
    def get_horizon_flat(self, kper_list):
        return self.factor_cosmo*np.array(kper_list)
    
    def plot(self, kper_list, kind='lst', utc=None, date=None, ax=None, **kargs):
        """
        Plot horizon line.

        Arguments:
        kper_list (array): 1D Array containing the k_per values for which horizon line is plotted.
        kind (optional, string): Three options (Default = 'lst'):
            'lst': Plot horizon lines at a particular lst or for a given lst range. 
            'full': Plot horizon line assuming full synthesis of 24 hours. 
            'flat': Plot horizon line with flat sky approximation.
        ax (optional, None): Axis on which to plot horizon line. If no axis is given, a new figure is created and the horizon line is plotted on it.
        date (string or 2-tuple): Necessary only when using kind='lst'. Date of observation (string in the format 'yyyy-mm-dd'). 
            If observation spans two days, use a 2-tuple: (day1, day2) with day1 and day2 in the format 'yyyy-mm-dd'. It is then necessary to use a 2-tuple for utc: (start_utc, end_utc).
        utc (optional, float or string or 2-tuple): Necessary only when using kind='lst'. There are 3 options:
            If 2-tuple (start_utc, end_utc) is given, horizon line for that observation duration is plotted when using kind = 'lst'. start_utc and end_utc can be either strings in the format 'hh:mm:ss' or floats in decimal hours.
            If single value is given (string in the format 'hh:mm:ss' or float in decimal hours), horizon line at that utc is plotted when using kind = 'lst'. 
            If blank, only full synthesis line and flat sky line can be plotted.
        """
        if kind=='lst' and (utc==None or date==None):
            raise ValueError("Argument utc and date are necessary if using kind='lst'")
        self.utc = utc
        self.date = date
        if ax is None:
            fig,ax = plt.subplots()    
        if kind == 'lst':
            kpar_list = self.get_horizon_lst(kper_list)
            kargs.setdefault('color', 'black')
            kargs.setdefault('linestyle', 'solid')
            kargs.setdefault('linewidth', 3)    
        elif kind == 'full':
            kpar_list = self.get_horizon_full(kper_list)
            kargs.setdefault('color', 'red')
            kargs.setdefault('linestyle', 'dashed')
            kargs.setdefault('linewidth', 3)
        elif kind == 'flat':
            kpar_list = self.get_horizon_flat(kper_list)
            kargs.setdefault('color', 'black')
            kargs.setdefault('linestyle', 'dotted')
            kargs.setdefault('linewidth', 2)    
        else:
            raise ValueError('kind = "lst"/"full"/"flat"')
        if kind != 'full' or abs(self.dec0) > np.pi/2-abs(self.lat):
            ax.plot(kper_list, kpar_list, **kargs)
        else:
            print('Full synthesis line is vertical!')
            
    def get_line(self, kper_list, kind='lst', utc=None, date=None):
        """
        Return k_parallel values for the horizon line.
        
        Arguments:
        kper_list, kind, date, utc: Same description as plot function
        """
        if kind=='lst' and (utc==None or date==None):
            raise ValueError("Argument utc and date are necessary if using kind='lst'")
        self.utc = utc
        self.date = date
        if kind == 'lst':
            kpar_list = self.get_horizon_lst(kper_list)
        elif kind == 'full':
            kpar_list = self.get_horizon_full(kper_list)
        elif kind == 'flat':
            kpar_list = self.get_horizon_flat(kper_list)
        else:
            raise ValueError('kind = "lst"/"full"/"flat"')
        if kind != 'full' or abs(self.dec0) > np.pi/2-abs(self.lat):
            return kpar_list
        else:
            print('Full synthesis line is vertical!')

class get_source:
    """
    Make and plot source lines for cylindrical power spectra.

    Arguments:
    freq_cen (float): Central frequency in Hz
    lat (float): Latitude in degree
    lon (float): Longitude in degree
    ra0 (float): RA of phase center in hour
    dec0 (float): Declination of phase center in degree
    ra (float): RA of source in hour
    dec (float): Declination of source in degree
    """
    def __init__(self, freq_cen, lat, lon, ra0, dec0, ra, dec):
        self.freq_cen = freq_cen
        self.z = (1.4204057e9/freq_cen)-1
        self.lat = lat*np.pi/180
        self.lon = lon*np.pi/180
        self.ra0 = ra0*np.pi/12
        self.dec0 = dec0*np.pi/180
        self.ra = ra*np.pi/12
        self.dec = dec*np.pi/180
        self.date = None
        self.utc = None
        self.factor_cosmo = ((cosmo.comoving_transverse_distance(self.z)/((const.c*(1+self.z)**2)/(cosmo.H(self.z)*1.4204057e9*u.Hz))).to(u.Hz).value)/freq_cen
    
    def get_factor_lst(self, lst):
        return ((np.sin(lst-self.ra)*np.cos(self.dec) - np.sin(lst-self.ra0)*np.cos(self.dec0))**2 + (np.sin(self.lat)*np.cos(self.dec)*np.cos(lst-self.ra) - np.sin(self.lat)*np.cos(self.dec0)*np.cos(lst-self.ra0) - np.sin(self.dec)*np.cos(self.lat) + np.sin(self.dec0)*np.cos(self.lat))**2)/np.sqrt((np.sin(lst-self.ra)*np.cos(self.dec) - np.sin(lst-self.ra0)*np.cos(self.dec0))**2 - ((np.sin(lst-self.ra)*np.cos(self.dec) - np.sin(lst-self.ra0)*np.cos(self.dec0))*np.sin(lst-self.ra0)*np.cos(self.dec0) + (np.sin(self.lat)*np.cos(self.dec0)*np.cos(lst-self.ra0) - np.sin(self.dec0)*np.cos(self.lat))*(np.sin(self.lat)*np.cos(self.dec)*np.cos(lst-self.ra) - np.sin(self.lat)*np.cos(self.dec0)*np.cos(lst-self.ra0) - np.sin(self.dec)*np.cos(self.lat) + np.sin(self.dec0)*np.cos(self.lat)))**2 + (np.sin(self.lat)*np.cos(self.dec)*np.cos(lst-self.ra) - np.sin(self.lat)*np.cos(self.dec0)*np.cos(lst-self.ra0) - np.sin(self.dec)*np.cos(self.lat) + np.sin(self.dec0)*np.cos(self.lat))**2)
    
    def get_factor_range(self, N=100):
        utc_start = time_to_h(self.utc[0])
        utc_end = time_to_h(self.utc[1])
        if isinstance(self.date, tuple):
            if len(self.date) > 2 or (dec_to_str(self.utc[1]) > dec_to_str(self.utc[0])):
                raise ValueError('Observation should be less than 24 hours.')
            else:
                utc_range = np.linspace(utc_start, 24+utc_end, N)
                factors = []
                for i in range(N):
                    if utc_range[i] < 24:
                        lst = utc_to_lst(utc_range[i], self.date[0], self.lat, self.lon)
                    else:
                        lst = utc_to_lst(utc_range[i]-24, self.date[1], self.lat, self.lon)
                    factors.append(self.get_factor_lst(lst))
                return np.amin(factors),np.amax(factors)
        else:
            if (time_to_h(self.utc[1]) <= time_to_h(self.utc[0])):
                raise ValueError('End UTC must be greater than start UTC for an observation spanning a single day.')
            else:
                utc_range = np.linspace(utc_start, utc_end, N)
                factors = []
                for i in utc_range:
                    lst = utc_to_lst(i, self.date, self.lat, self.lon)
                    factors.append(self.get_factor_lst(lst))
                return np.amin(factors), np.amax(factors)
    
    def get_source_range(self, kper_list):
        if len(self.utc) != 2:
            raise ValueError('UTC should be either a single value or a tuple of two values corresponding to the start and end of the observation')
        else:
            factor_min, factor_max = self.get_factor_range()
        return factor_min*self.factor_cosmo*kper_list, factor_max*self.factor_cosmo*kper_list
    
    def get_source_lst(self, kper_list):    
        lst = utc_to_lst(self.utc, self.date, self.lat, self.lon)
        factor = self.get_factor_lst(lst)
        return factor*self.factor_cosmo*kper_list
    
    def get_source_flat(self, kper_list):
        angle = angular_distance(self.ra,self.dec,self.ra0,self.dec0)
        return self.factor_cosmo*np.array(kper_list)*np.sin(angle)
        
    def plot(self, kper_list, kind='lst', ax=None, utc=None, date=None, **kargs):
        """
        Plot source line.

        Arguments:
        kper_list (array): 1D Array containing the k_per values for which horizon line is plotted.
        kind (optional, string): Two options (Default = 'lst'): 
            'lst': Plot source lines at a particular lst or for a given lst range. 
            'flat': Plot horizon line with flat sky approximation.
        ax (optional, None): Axis on which to plot source line. If no axis is given, a new figure is created and the source line is plotted on it.
        date (string or 2-tuple): Necessary only when using kind='lst'. Date of observation (string in the format 'yyyy-mm-dd'). 
            If observation spans two days, use a 2-tuple: (day1, day2) with day1 and day2 in the format 'yyyy-mm-dd'. It is then necessary to use a 2-tuple for utc: (start_utc, end_utc).
        utc (optional, float or string or 2-tuple): Necessary only when using kind='lst'. 3 options:
            If 2-tuple (start_utc, end_utc) is given, two source lines (minimum and maximum) for that observation duration is plotted when using kind = 'lst'. start_utc and end_utc can be either strings in the format 'hh:mm:ss' or floats in decimal hours.
            If single value is given (string in the format 'hh:mm:ss' or float in decimal hours), source line at that utc is plotted when using kind = 'lst'.
            If blank, only flat sky approximation line can be plotted        
        """
        if kind=='lst' and (utc==None or date==None):
            raise ValueError("Argument utc and date are necessary if using kind='lst'")
        self.utc = utc
        self.date = date
        if ax is None:
            fig,ax = plt.subplots()    
        if kind == 'lst':
            if isinstance(self.utc, tuple):
                kpar_list_min, kpar_list_max = self.get_source_range(kper_list)
                kargs.setdefault('color', 'black')
                kargs.setdefault('linestyle', 'solid')
                kargs.setdefault('linewidth', 2)
                ax.plot(kper_list, kpar_list_min, **kargs)
                ax.plot(kper_list, kpar_list_max, **kargs)
            else:
                kpar_list = self.get_source_lst(kper_list)
                kargs.setdefault('color', 'black')
                kargs.setdefault('linestyle', 'solid')
                kargs.setdefault('linewidth', 2)
                ax.plot(kper_list, kpar_list, **kargs)             
        elif kind == 'flat':
            kpar_list = self.get_source_flat(kper_list)
            kargs.setdefault('color', 'black')
            kargs.setdefault('linestyle', 'dotted')
            kargs.setdefault('linewidth', 2)    
            ax.plot(kper_list, kpar_list, **kargs)             
        else:
            raise ValueError('kind = "lst"/flat"')
            
    def get_line(self, kper_list, kind='lst', utc=None, date=None):
        """
        Return k_parallel values for the source line.
        
        Arguments:
        kper_list, kind, date, utc: Same description as plot function
        """
        if kind=='lst' and (utc==None or date==None):
            raise ValueError("Argument utc and date are necessary if using kind='lst'")
        self.utc = utc
        self.date = date
        if kind == 'lst':
            if isinstance(self.utc, tuple):
                kpar_list_min, kpar_list_max = self.get_source_range(kper_list)
                return kpar_list_min, kpar_list_max
            else:
                kpar_list = self.get_source_lst(kper_list)
                return kpar_list
        elif kind == 'flat':
            kpar_list = self.get_source_flat(kper_list)
            return kpar_list
        else:
            raise ValueError('kind = "lst"/flat"')