"""
readFL
======

Tool for importing and optionally slicing flight-level data
from the NOAA P-3. Additional functionality may be added in
the future as needed.

    getP3

"""
import xarray as xr
import numpy as np
from datetime import datetime as dt
import datetime
import sys
import warnings

warnings.filterwarnings('ignore', 'invalid value encountered in less')
warnings.filterwarnings('ignore', 'invalid value encountered in greater')

def getP3(flFile,strtDT=None,endDT=None,readInsitu=False,projID=None):
    """
    Load in a specified flight-level file (from the NOAA P-3)
    
    
    Parameters
    ----------
    flFile : string
        String specifying the path where the flight level file
        is located (including the full file name).
    strtDT : string, optional
        String specifying the start date and time of the period
        desired. Must be in 'YYYYmmdd-HHMMSS' format. Defaults
        to None, in which case date from the beginning of the flight
        to endDT will be included.
    endDT : string, optional
        String specifying the end date and time of the period
        desired. Must be in 'YYYYmmdd-HHMMSS' format. Defaults
        to None, in which case date through the end of the flight
        will be included.
    readInsitu : bool, optional
        If `True` will read in basic in situ obs from the P-3, including
        temperature, dewpoint, and wind speed and direction.
        
    Returns
    -------
    flData_out : dict
        Dictionary containing FL data bounded by the start/end times (if given).

    """
    
    # Convert input start/end strings into datetimes (if not None)
    if strtDT is not None:
        strtDT = dt.strptime(strtDT,'%Y%m%d-%H%M%S')
    if endDT is not None:
        endDT = dt.strptime(endDT,'%Y%m%d-%H%M%S')
    
    # Define the names of import variables as given in the FL
    # netCDF input file
    latVar = 'LATref'
    lonVar = 'LONref'
    altVar = 'ALTref'
    headVar = 'THDGref'
    rollVar = 'ROLLref'
    tempVar = 'TA.d'
    dewVar = 'TDM.2' #20170430 - VSE
    wdVar = 'WD.d'
    wsVar = 'WS.d' # m/s
    
    # Import all needed variables
    flData = xr.open_dataset(flFile,decode_times=False)

    flTime = flData.get('Time')
    
    if projID == 'VSE':
        baseT = dt.strptime(flTime.units,'seconds since %Y-%m-%d %H:%M:%S %z') #2017 data
    elif projID == 'TORUS':
        baseT = dt.utcfromtimestamp(flData.StartTime) # 2019 data
    else:
        try:
            baseT = dt.strptime(flTime.units,'seconds since %Y-%m-%d %H:%M:%S %z')
        except RuntimeError:
            sys.exit('Flight-level base time variable does not match expected format. Exiting...')
            
    flTimeSec = flTime.to_masked_array()[10:] # Temp hack for 23 May 2022
    dtArr = np.asarray([baseT + datetime.timedelta(seconds=int(t)) for t in flTimeSec])
    
    dtArr = np.asarray([d.replace(tzinfo=None) for d in dtArr])
    
    #print('\tFlight-level begDT: {:%Y-%m-%d %H:%M:%S}'.format(dtArr[0]))
    #print('\tFlight-level endDT: {:%Y-%m-%d %H:%M:%S}'.format(dtArr[-1]))

    lat = flData.get(latVar).to_masked_array()[10:] # Temp hack for 23 May 2022
    lon = flData.get(lonVar).to_masked_array()[10:] # Temp hack for 23 May 2022
    alt = flData.get(altVar).to_masked_array()[10:] # Temp hack for 23 May 2022
    hdng = flData.get(headVar).to_masked_array()[10:] # Temp hack for 23 May 2022
    roll = flData.get(rollVar).to_masked_array()[10:] # Temp hack for 23 May 2022
    
    if readInsitu:
        tempC = flData.get(tempVar).to_masked_array()
        tdewC = flData.get(dewVar).to_masked_array()
        ws = flData.get(wsVar).to_masked_array()
        wd = flData.get(wdVar).to_masked_array()
    
    latDiff = np.append(0,np.diff(lat))
    lonDiff = np.append(0,np.diff(lon))
    badLat = np.squeeze(np.where(np.logical_or(latDiff > 0.1,latDiff < -0.1)))
    badLon = np.squeeze(np.where(np.logical_or(lonDiff > 0.1,lonDiff < -0.1)))
    lat[badLat] = np.nan
    lon[badLon] = np.nan
    
    np.ma.set_fill_value(dtArr,np.nan)
    np.ma.set_fill_value(lat,np.nan)
    np.ma.set_fill_value(lon,np.nan)
    np.ma.set_fill_value(alt,np.nan)
    np.ma.set_fill_value(hdng,np.nan)
    np.ma.set_fill_value(roll,np.nan)
    
    if readInsitu:
        np.ma.set_fill_value(tempC,np.nan)
        np.ma.set_fill_value(tdewC,np.nan)
        np.ma.set_fill_value(ws,np.nan)
        np.ma.set_fill_value(wd,np.nan)
    
    
    # If start and/or end dates and times are given, slice
    # data arrays appropriately
    if isinstance(strtDT,dt):
        strtIx = np.squeeze(np.where(dtArr == strtDT))
    elif strtDT is None:
        strtIx = None
    else:
        sys.exit('Starting date/time argument must be in datetime format, or None to include data from begin of flight')

    if isinstance(endDT,dt):
        endIx = np.squeeze(np.where(dtArr == endDT)) + 1
    elif endDT is None:
        endIx = None
    else:
        sys.exit('Ending date/time argument must be in datetime format, or None to include data to end of flight')
        
    
    dt_out = dtArr[strtIx:endIx]
    lat_out = lat[strtIx:endIx]
    lon_out = lon[strtIx:endIx]
    alt_out = alt[strtIx:endIx]
    hdng_out = hdng[strtIx:endIx]
    roll_out = roll[strtIx:endIx]
    
    if readInsitu:
        tempC_out = tempC[strtIx:endIx]
        tdewC_out = tdewC[strtIx:endIx]
        ws_out = ws[strtIx:endIx]
        wd_out = wd[strtIx:endIx]
    
        flData_out = {'flDT': dt_out, 'flLat': lat_out, 'flLon': lon_out, 'flAlt': alt_out, 'flHdng': hdng_out, 'flRoll': roll_out,
                      'flTempC': tempC_out, 'flTdewC': tdewC_out, 'flWS': ws_out, 'flWD': wd_out}
    else:
        flData_out = {'flDT': dt_out, 'flLat': lat_out, 'flLon': lon_out, 'flAlt': alt_out, 'flHdng': hdng_out, 'flRoll': roll_out}
    
    return flData_out