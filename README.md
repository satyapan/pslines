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
- date (string): Necessary only when using kind='lst'. Start date of observation (string in the format 'yyyy-mm-dd').
- utc (optional, float or string): Necessary only when using kind='lst'. Start UTC of the observation (string in the format 'hh:mm:ss' or float in decimal hours).
- duration (optional, float): Duration of the observation in hours. 
      If duration = 0 or None, horizon line at the given utc is plotted.
      Otherwise horizon line for the observation duration is plotted.

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
- date (string): Necessary only when using kind='lst'. Start date of observation (string in the format 'yyyy-mm-dd').
- utc (optional, float or string): Necessary only when using kind='lst'. Start UTC of the observation (string in the format 'hh:mm:ss' or float in decimal hours).
- duration (optional, float): Duration of the observation in hours. 
      If duration = 0 or None, source line at the given utc is plotted.
      Otherwise two source lines (minimum and maximum) for that observation duration is plotted. 
