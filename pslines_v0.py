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
    hh, mm, ss = map(int, time_str.split(':'))
    total_hours = hh + (mm / 60) + (ss / 3600)
    return total_hours

def utc_to_lst(utc, date, lat, lon):
    nancay = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=0*u.m)
    time = date+'T'+dec_to_str(utc)
    t = Time(time, scale='utc', location=nancay)
    return t.sidereal_time('mean').rad

class get_horizon:
    def __init__(self, freq_cen, lat, lon, ra0, dec0, date, utc=None):
        self.freq_cen = freq_cen
        self.z = (1.4204057e9/freq_cen)-1
        self.lat = lat*np.pi/180
        self.lon = lon*np.pi/180
        self.ra0 = ra0*np.pi/180
        self.dec0 = dec0*np.pi/180
        self.date = date
        self.utc = utc
        self.factor_cosmo = ((cosmo.comoving_transverse_distance(self.z)/((const.c*(1+self.z)**2)/(cosmo.H(self.z)*1.4204057e9*u.Hz))).to(u.Hz).value)/freq_cen
    
    def get_factor_lst(self, lst):
        return (1+np.sqrt(1-(np.sin(self.lat)*np.sin(self.dec0)+np.cos(self.lat)*np.cos(self.dec0)*np.cos(lst-self.ra0))**2))/(np.sin(self.lat)*np.sin(self.dec0)+np.cos(self.lat)*np.cos(self.dec0)*np.cos(lst-self.ra0))
    
    def get_factor_max(self, N=100):
        utc_start = time_to_h(self.utc[0])
        utc_end = time_to_h(self.utc[1])
        if isinstance(self.date, tuple):
            if len(self.date) > 2:
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
                return np.max(factors)
        else:
            utc_range = np.linspace(utc_start, utc_end, N)
            factors = []
            for i in utc_range:
                lst = utc_to_lst(i, self.date, self.lat, self.lon)
                factors.append(self.get_factor_lst(lst))
            return np.max(factors)
    
    def get_horizon_lst(self, kper_list):
        if isinstance(self.utc, tuple):
            if len(self.utc) != 2:
                raise ValueError('utc should be either a single value or a tuple of two values corresponding to the start and end lst')
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
    
    def plot(self, kper_list, kind='lst', ax=None, **kargs):
        if ax is None:
            fig,ax = plt.subplots()    
        if kind == 'lst':
            kpar_list = self.get_horizon_lst(kper_list)
            kargs.setdefault('color', 'white')
            kargs.setdefault('linestyle', 'solid')
            kargs.setdefault('linewidth', 2)    
        elif kind == 'full':
            kpar_list = self.get_horizon_full(kper_list)
            kargs.setdefault('color', 'red')
            kargs.setdefault('linestyle', 'dashed')
            kargs.setdefault('linewidth', 2)
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