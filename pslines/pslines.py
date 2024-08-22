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
    loc = EarthLocation(lat=lat*u.rad, lon=lon*u.rad, height=0*u.m)
    time = date+'T'+dec_to_str(utc)
    t = Time(time, scale='utc', location=loc)
    return t.sidereal_time('mean').rad

def angular_distance(ra1, dec1, ra2, dec2):
    dlon = ra2 - ra1
    dlat = dec2 - dec1
    a = np.sin(dlat/2)**2+np.cos(dec1)*np.cos(dec2)*np.sin(dlon/2)**2
    angle = 2 * np.arctan(np.sqrt(a)/np.sqrt(1-a))
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
        self.duration = None
        self.factor_cosmo = ((cosmo.comoving_transverse_distance(self.z)/((const.c*(1+self.z)**2)/(cosmo.H(self.z)*1.4204057e9*u.Hz))).to(u.Hz).value)/freq_cen
    
    def get_factor_lst(self, lst):
        factor = (1+np.sqrt(1-(np.sin(self.lat)*np.sin(self.dec0)+np.cos(self.lat)*np.cos(self.dec0)*np.cos(lst-self.ra0))**2))/(np.sin(self.lat)*np.sin(self.dec0)+np.cos(self.lat)*np.cos(self.dec0)*np.cos(lst-self.ra0))
        if factor < 0:
            print('Phase Center goes below the horizon at lst = %0.3fh!'%(lst*12/np.pi))
        return factor
        
    def get_factor_max(self, N=100):
        date_start = self.date
        utc_start = time_to_h(self.utc)
        date_end = self.date
        utc_end = time_to_h(self.utc) + self.duration
        if utc_end >= 24:
            utc_end -= 24
            date_end = (datetime.datetime.strptime(date_start, '%Y-%m-%d')+datetime.timedelta(days=1)).strftime('%Y-%m-%d')
        if utc_end > utc_start:
            utc_range = np.linspace(utc_start, utc_end, N)
            factors = []
            for i in utc_range:
                lst = utc_to_lst(i, date_start, self.lat, self.lon)
                factors.append(self.get_factor_lst(lst))
            return np.amax(factors)
        else:
            utc_range = np.linspace(utc_start, 24+utc_end, N)
            factors = []
            for i in range(N):
                if utc_range[i] < 24:
                    lst = utc_to_lst(utc_range[i], date_start, self.lat, self.lon)
                else:
                    lst = utc_to_lst(utc_range[i]-24, date_end, self.lat, self.lon)
                factors.append(self.get_factor_lst(lst))
            return np.amax(factors)
    
    def get_horizon_lst(self, kper_list):
        if self.duration == None or self.duration == 0:
            lst = utc_to_lst(self.utc, self.date, self.lat, self.lon)
            factor = self.get_factor_lst(lst)
        elif self.duration > 24:
            raise ValueError('If observation spans more than 24 hours, use full synthesis line.')
        elif self.duration < 0:
            raise ValueError('Duration needs to be positive.')
        else:
            factor = self.get_factor_max()
        return factor*self.factor_cosmo*kper_list
    
    def get_horizon_full(self, kper_list):    
        factor = (1+np.sin(abs(self.lat+self.dec0)))/(-np.cos(self.lat+self.dec0))
        return factor*self.factor_cosmo*kper_list
    
    def get_horizon_flat(self, kper_list):
        return self.factor_cosmo*np.array(kper_list)
    
    def plot(self, kper_list, kind='lst', utc=None, date=None, duration=None, ax=None, **kargs):
        """
        Plot horizon line.

        Arguments:
        kper_list (array): 1D Array containing the k_per values for which horizon line is plotted.
        kind (optional, string): Three options (Default = 'lst'):
            'lst': Plot horizon lines at a particular lst or for a given lst range. 
            'full': Plot horizon line assuming full synthesis of 24 hours. 
            'flat': Plot horizon line with flat sky approximation.
        ax (optional, None): Axis on which to plot horizon line. If no axis is given, a new figure is created and the horizon line is plotted on it.
        date (string): Necessary only when using kind='lst'. Start date of observation (string in the format 'yyyy-mm-dd').
        utc (optional, float or string): Necessary only when using kind='lst'. Start UTC of the observation (string in the format 'hh:mm:ss' or float in decimal hours).
        duration (optional, float): Duration of the observation in hours. 
            If duration = 0 or None, horizon line at the given utc is plotted.
            Otherwise horizon line for the observation duration is plotted.
        """
        if kind=='lst' and (utc==None or date==None):
            raise ValueError("Argument utc and date are necessary if using kind='lst'")
        self.utc = utc
        self.date = date
        self.duration = duration
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
            kargs.setdefault('linewidth', 3)    
        else:
            raise ValueError('kind = "lst"/"full"/"flat"')
        if kind != 'full' or abs(self.dec0+self.lat) > np.pi/2:
            ax.plot(kper_list, kpar_list, **kargs)
        else:
            print('Full synthesis line is vertical!')
            
    def get_line(self, kper_list, kind='lst', utc=None, date=None, duration=None):
        """
        Return k_parallel values for the horizon line.
        
        Arguments:
        kper_list, kind, date, utc: Same description as plot function
        """
        if kind=='lst' and (utc==None or date==None):
            raise ValueError("Argument utc and date are necessary if using kind='lst'")
        self.utc = utc
        self.date = date
        self.duration = duration
        if kind == 'lst':
            kpar_list = self.get_horizon_lst(kper_list)
        elif kind == 'full':
            kpar_list = self.get_horizon_full(kper_list)
        elif kind == 'flat':
            kpar_list = self.get_horizon_flat(kper_list)
        else:
            raise ValueError('kind = "lst"/"full"/"flat"')
        if kind != 'full' or abs(self.dec0+self.lat) > np.pi/2:
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
    def __init__(self, freq_cen, lat, lon, ra0, dec0, ra, dec, min_alt=0):
        self.freq_cen = freq_cen
        self.z = (1.4204057e9/freq_cen)-1
        self.lat = lat*np.pi/180
        self.lon = lon*np.pi/180
        self.ra0 = ra0*np.pi/12
        self.dec0 = dec0*np.pi/180
        self.ra = ra*np.pi/12
        self.dec = dec*np.pi/180
        self.min_alt = min_alt*np.pi/180
        self.date = None
        self.utc = None
        self.duration = None
        self.factor_cosmo = ((cosmo.comoving_transverse_distance(self.z)/((const.c*(1+self.z)**2)/(cosmo.H(self.z)*1.4204057e9*u.Hz))).to(u.Hz).value)/freq_cen

    def get_alt(self, lst):
        return np.arcsin(np.sin(self.lat)*np.sin(self.dec)+np.cos(self.lat)*np.cos(self.dec)*np.cos(lst-self.ra))
    
    def get_factor_lst(self, lst):
        return ((np.sin(lst-self.ra)*np.cos(self.dec) - np.sin(lst-self.ra0)*np.cos(self.dec0))**2 + (np.sin(self.lat)*np.cos(self.dec)*np.cos(lst-self.ra) - np.sin(self.lat)*np.cos(self.dec0)*np.cos(lst-self.ra0) - np.sin(self.dec)*np.cos(self.lat) + np.sin(self.dec0)*np.cos(self.lat))**2)/np.sqrt((np.sin(lst-self.ra)*np.cos(self.dec) - np.sin(lst-self.ra0)*np.cos(self.dec0))**2 - ((np.sin(lst-self.ra)*np.cos(self.dec) - np.sin(lst-self.ra0)*np.cos(self.dec0))*np.sin(lst-self.ra0)*np.cos(self.dec0) + (np.sin(self.lat)*np.cos(self.dec0)*np.cos(lst-self.ra0) - np.sin(self.dec0)*np.cos(self.lat))*(np.sin(self.lat)*np.cos(self.dec)*np.cos(lst-self.ra) - np.sin(self.lat)*np.cos(self.dec0)*np.cos(lst-self.ra0) - np.sin(self.dec)*np.cos(self.lat) + np.sin(self.dec0)*np.cos(self.lat)))**2 + (np.sin(self.lat)*np.cos(self.dec)*np.cos(lst-self.ra) - np.sin(self.lat)*np.cos(self.dec0)*np.cos(lst-self.ra0) - np.sin(self.dec)*np.cos(self.lat) + np.sin(self.dec0)*np.cos(self.lat))**2)

    def get_factor_range(self, N=100):
        date_start = self.date
        utc_start = time_to_h(self.utc)
        date_end = self.date
        utc_end = time_to_h(self.utc) + self.duration
        if utc_end >= 24:
            utc_end -= 24
            date_end = (datetime.datetime.strptime(date_start, '%Y-%m-%d')+datetime.timedelta(days=1)).strftime('%Y-%m-%d')        
        if utc_end > utc_start:
            utc_range = np.linspace(utc_start, utc_end, N)
            factors = []
            for i in utc_range:
                lst = utc_to_lst(i, date_start, self.lat, self.lon)
                alt = self.get_alt(lst)
                if alt >= self.min_alt:
                    factors.append(self.get_factor_lst(lst))
            return np.amin(factors), np.amax(factors)
        else:
            utc_range = np.linspace(utc_start, 24+utc_end, N)
            factors = []
            for i in range(N):
                if utc_range[i] < 24:
                    lst = utc_to_lst(utc_range[i], date_start, self.lat, self.lon)
                    alt = self.get_alt(lst)
                else:
                    lst = utc_to_lst(utc_range[i]-24, date_end, self.lat, self.lon)
                    alt = self.get_alt(lst)
                if alt >= self.min_alt:
                    factors.append(self.get_factor_lst(lst))
            return np.amin(factors),np.amax(factors)
    
    def get_source_range(self, kper_list):
        factor_min, factor_max = self.get_factor_range()
        return factor_min*self.factor_cosmo*kper_list, factor_max*self.factor_cosmo*kper_list
    
    def get_source_lst(self, kper_list):    
        lst = utc_to_lst(self.utc, self.date, self.lat, self.lon)
        factor = self.get_factor_lst(lst)
        return factor*self.factor_cosmo*kper_list
    
    def get_source_flat(self, kper_list):
        angle = angular_distance(self.ra,self.dec,self.ra0,self.dec0)
        return self.factor_cosmo*np.array(kper_list)*np.sin(angle)
        
    def plot(self, kper_list, kind='lst', ax=None, utc=None, date=None, duration=None, **kargs):
        """
        Plot source line.

        Arguments:
        kper_list (array): 1D Array containing the k_per values for which horizon line is plotted.
        kind (optional, string): Two options (Default = 'lst'): 
            'lst': Plot source lines at a particular lst or for a given lst range. 
            'flat': Plot horizon line with flat sky approximation.
        ax (optional, None): Axis on which to plot source line. If no axis is given, a new figure is created and the source line is plotted on it.
        date (string): Necessary only when using kind='lst'. Start date of observation (string in the format 'yyyy-mm-dd').
        utc (optional, float or string): Necessary only when using kind='lst'. Start UTC of the observation (string in the format 'hh:mm:ss' or float in decimal hours).
        duration (optional, float): Duration of the observation in hours. 
            If duration = 0 or None, source line at the given utc is plotted.
            Otherwise two source lines (minimum and maximum) for that observation duration is plotted.      
        """
        if kind=='lst' and (utc==None or date==None):
            raise ValueError("Argument utc and date are necessary if using kind='lst'")
        self.utc = utc
        self.date = date
        self.duration = duration
        if ax is None:
            fig,ax = plt.subplots()    
        if kind == 'lst':
            if duration == None or duration == 0:
                kpar_list = self.get_source_lst(kper_list)
                kargs.setdefault('color', 'black')
                kargs.setdefault('linestyle', 'solid')
                kargs.setdefault('linewidth', 3)
                ax.plot(kper_list, kpar_list, **kargs)   
            elif duration < 0 or duration > 24:
                raise ValueError('Duration needs to be between 0 and 24.')
            else:
                kpar_list_min, kpar_list_max = self.get_source_range(kper_list)
                kargs.setdefault('color', 'black')
                kargs.setdefault('linestyle', 'solid')
                kargs.setdefault('linewidth', 3)
                ax.plot(kper_list, kpar_list_min, **kargs)
                ax.plot(kper_list, kpar_list_max, **kargs)
        elif kind == 'flat':
            kpar_list = self.get_source_flat(kper_list)
            kargs.setdefault('color', 'black')
            kargs.setdefault('linestyle', 'dotted')
            kargs.setdefault('linewidth', 3)    
            ax.plot(kper_list, kpar_list, **kargs)             
        else:
            raise ValueError('kind = "lst"/flat"')
            
    def get_line(self, kper_list, kind='lst', utc=None, date=None, duration=None):
        """
        Return k_parallel values for the source line.
        
        Arguments:
        kper_list, kind, date, utc: Same description as plot function
        """
        if kind=='lst' and (utc==None or date==None):
            raise ValueError("Argument utc and date are necessary if using kind='lst'")
        self.utc = utc
        self.date = date
        self.duration = duration
        if kind == 'lst':
            if duration == None or duration == 0:
                kpar_list = self.get_source_lst(kper_list)
                return kpar_list
            elif duration < 0 or duration > 24:
                raise ValueError('Duration needs to be between 0 and 24.')
            else:
                kpar_list_min, kpar_list_max = self.get_source_range(kper_list)
                return kpar_list_min, kpar_list_max
        elif kind == 'flat':
            kpar_list = self.get_source_flat(kper_list)
            return kpar_list
        else:
            raise ValueError('kind = "lst"/flat"')