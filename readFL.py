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

def getP3(flFile,strtDT=None,endDT=None):
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
    
    # Import all needed variables
    flData = xr.open_dataset(flFile,decode_times=False)

    flDate = flData.FlightDate
    HH = flData.get('HH').to_masked_array()
    
    # Find indices where time variables (HH as proxy) are missing
    #   and use that to slice the other variables appropriately
    flNoMsk = np.where(HH.mask == False)[0]

    HH = HH[flNoMsk]
    [flNoMsk]
    MM = flData.get('MM').to_masked_array()[flNoMsk]
    SS = flData.get('SS').to_masked_array()[flNoMsk]
    lat = flData.get(latVar).to_masked_array()[flNoMsk]
    lon = flData.get(lonVar).to_masked_array()[flNoMsk]
    alt = flData.get(altVar).to_masked_array()[flNoMsk]
    hdng = flData.get(headVar).to_masked_array()[flNoMsk]
    roll = flData.get(rollVar).to_masked_array()[flNoMsk]
    
    latDiff = np.append(0,np.diff(lat))
    lonDiff = np.append(0,np.diff(lon))
    badLat = np.squeeze(np.where(np.logical_or(latDiff > 0.1,latDiff < -0.1)))
    badLon = np.squeeze(np.where(np.logical_or(lonDiff > 0.1,lonDiff < -0.1)))
    lat[badLat] = np.nan
    lon[badLon] = np.nan
    
    np.ma.set_fill_value(lat,np.nan)
    np.ma.set_fill_value(lon,np.nan)
    np.ma.set_fill_value(alt,np.nan)
    np.ma.set_fill_value(hdng,np.nan)
    np.ma.set_fill_value(roll,np.nan)
    
    
    # flTS = (flT - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
#     dtArr = np.asarray([dt.utcfromtimestamp(t) for t in flTS])
    
    # Combine flight date and each time variable into a datetime array
    try:
        crntDay = dt.strptime(flDate,'%m/%d/%Y')
    except:
        crntDay = dt.strptime(flDate,'%Y-%m-%d')
    if np.min(np.diff(HH)) == -23.0:
        crsMidnt = True
        midntIx = np.where(np.diff(HH) == -23.0)[0][0] + 1
        nxtDay = crntDay + datetime.timedelta(days=1)
    else:
        crsMidnt = False

    dtArr = np.empty(np.shape(HH),dtype=object)

    if crsMidnt:
        for idt in np.arange(np.size(HH[:midntIx])):
            dtArr[idt] = dt.strptime('{:%Y%m%d}{:02.0f}{:02.0f}{:02.0f}'.format(crntDay,HH[idt],MM[idt],SS[idt]),'%Y%m%d%H%M%S')
        for idt in range(np.size(HH[:midntIx]),np.size(HH)):
            dtArr[idt] = dt.strptime('{:%Y%m%d}{:02.0f}{:02.0f}{:02.0f}'.format(nxtDay,HH[idt],MM[idt],SS[idt]),'%Y%m%d%H%M%S')
    else:
        for idt in np.arange(np.size(HH)):
            dtArr[idt] = dt.strptime('{:%Y%m%d}{:02.0f}{:02.0f}{:02.0f}'.format(crntDay,HH[idt],MM[idt],SS[idt]),'%Y%m%d%H%M%S')

    # dtArrMsked = np.ma.MaskedArray(dtArr,mask=np.ma.getmaskarray(HH),fill_value=np.nan)
    
    
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
    
    flData_out = {'flDT': dt_out, 'flLat': lat_out, 'flLon': lon_out, 'flAlt': alt_out, 'flHdng': hdng_out, 'flRoll': roll_out}
    
    
    return flData_out