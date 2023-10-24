# Plotting horizon and source lines with pslines

## Import module
```
from pslines import *
```

## Horizon line

#### Step 1: Make horizon line object
```
horizon = get_horizon(freq_cen, lat, lon, ra0, dec0)
```
Arguments  
- freq_cen (float): Central frequency in Hz  
- lat (float): Latitude in degree  
- lon (float): Longitude in degree  
- ra0 (float): RA of phase center in hour  
- dec0 (float): Declination of phase center in degree  

#### Step 2: Plot horizon line
```
horizon.plot(kper_list, ax=ax, kind='lst', utc=(utc_start, utc_end), date=date)
```
Arguments:  
- kper_list (array): 1D Array containing the k_per values for which horizon line is plotted.  
- kind (optional, string): Three options (Default = 'lst'):  
    'lst': Plot horizon lines at a particular lst or for a given lst range.   
    'full': Plot horizon line assuming full synthesis of 24 hours.  
    'flat': Plot horizon line with flat sky approximation.  
- ax (optional, None): Axis on which to plot horizon line. If no axis is given, a new figure is created and the horizon line is plotted on it.  
- date (string or 2-tuple): Necessary only when using kind='lst'. Date of observation (string in the format 'yyyy-mm-dd').   
    If observation spans two days, use a 2-tuple: (day1, day2) with day1 and day2 in the format 'yyyy-mm-dd'. It is then necessary to use a 2-tuple for utc: (start_utc, end_utc).  
- utc (optional, float or string or 2-tuple): Necessary only when using kind='lst'. There are 3 options:  
    If 2-tuple (start_utc, end_utc) is given, horizon line for that observation duration is plotted when using kind = 'lst'. start_utc and end_utc can be either strings in the format 'hh:mm:ss' or floats in decimal hours.  
    If single value is given (string in the format 'hh:mm:ss' or float in decimal hours), horizon line at that utc is plotted when using kind = 'lst'.  
    If blank, only full synthesis line and flat sky line can be plotted.

## Source line

#### Step 1: Make source line object
```
source = get_source(freq_cen, lat, lon, ra0, dec0, ra, dec)
```
Arguments:  
- freq_cen (float): Central frequency in Hz  
- lat (float): Latitude in degree  
- lon (float): Longitude in degree  
- ra0 (float): RA of phase center in hour  
- dec0 (float): Declination of phase center in degree  
- ra (float): RA of source in hour  
- dec (float): Declination of source in degree  

#### Step 2: Plot source line
```
source.plot(kper_list, ax=ax, kind='lst', utc=(utc_start,utc_end), date=date)
```
Arguments:  
- kper_list (array): 1D Array containing the k_per values for which horizon line is plotted.  
- kind (optional, string): Two options (Default = 'lst'):  
    'lst': Plot source lines at a particular lst or for a given lst range.  
    'flat': Plot horizon line with flat sky approximation.  
- ax (optional, None): Axis on which to plot source line. If no axis is given, a new figure is created and the source line is plotted on it.  
- date (string or 2-tuple): Necessary only when using kind='lst'. Date of observation (string in the format 'yyyy-mm-dd').  
    If observation spans two days, use a 2-tuple: (day1, day2) with day1 and day2 in the format 'yyyy-mm-dd'. It is then necessary to use a 2-tuple for utc: (start_utc, end_utc).  
- utc (optional, float or string or 2-tuple): Necessary only when using kind='lst'. 3 options:  
    If 2-tuple (start_utc, end_utc) is given, two source lines (minimum and maximum) for that observation duration is plotted when using kind = 'lst'. start_utc and end_utc can be either strings in the format 'hh:mm:ss' or floats in decimal hours.  
    If single value is given (string in the format 'hh:mm:ss' or float in decimal hours), source line at that utc is plotted when using kind = 'lst'.  
    If blank, only flat sky approximation line can be plotted
