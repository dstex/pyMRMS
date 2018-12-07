"""
retrieveMRMS
======

Download MRMS Seamless Hybrid Scan Reflectivity data.

    fetch
    extract

"""
import urllib.request as request
from datetime import datetime as dt
import datetime
from subprocess import call
import os
from glob import glob
import numpy as np
import xarray as xr

def fetch(strtDT,endDT,saveDir):
    """
    Download and unzip MRMS files from the Iowa State data 
    archive between a given start and end date/time. These 
    files exist at 2-minute intervals.
    
    
    Parameters
    ----------
    strtDT : string
        String specifying the start date and time of the period
        desired. Must be in 'YYYYmmdd-HHMMSS' format.
    endDT : string
        String specifying the end date and time of the period
        desired. Must be in 'YYYYmmdd-HHMMSS' format.
    saveDir : string
        String specifying the path where the downloaded files
        are to be saved. If the given path does not exist, it
        will be created.

    """


    # strtDT = '20150620-040000'
    # endDT = '20150620-050000'
    sDT = dt.strptime(strtDT,'%Y%m%d-%H%M%S')
    eDT = dt.strptime(endDT,'%Y%m%d-%H%M%S')

    # saveDir = '/Users/danstechman/Desktop/MRMS/'
    rmtBase = 'http://mtarchive.geol.iastate.edu/'

    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    tStep = sDT
    while tStep <= eDT:
        locPath = '{}SeamlessHSR_00.00_{:%Y%m%d-%H%M}00.grib2.gz'.format(saveDir,tStep)
        rmtPath = '{}{:%Y/%m/%d}/mrms/ncep/SeamlessHSR/SeamlessHSR_00.00_{:%Y%m%d-%H%M}00.grib2.gz'.format(rmtBase,tStep,tStep)
        
        if os.path.isfile('{}nc'.format(locPath[:-8])):
            print('\tSkipping download and conversion as final netCDF file already exists.')
        elif os.path.isfile(locPath[:-3]):
            print('Skipping download of {:%Y%m%d - %H:%M}. File already exists locally.'.format(tStep))

            print('\tConverting existing grib file to netCDF...')
            call(['ncl_convert2nc',locPath[:-3],'-o',format(locPath[:-42])])
            call(['rm',locPath[:-3]])
        else:
            print('Currently downloading {:%Y%m%d - %H:%M}...'.format(tStep))

            request.urlretrieve(rmtPath,locPath)
            call(['gunzip',locPath])
            call(['ncl_convert2nc',locPath[:-3],'-o',format(locPath[:-42])])
            call(['rm',locPath[:-3]])
        tStep += datetime.timedelta(minutes=2)
        
def extract(mrmsDir,strtDT=None,endDT=None,llCrds=None,urCrds=None):
    """
    Create a masked array containing all MRMS data
    between a given set of times, and optionally limit
    the domain to a given region.
    
    
    Parameters
    ----------
    mrmsDir : string
        String specifying the path to the directory containing
        the MRMS netCDF files.
    strtDT : string, optional
        String specifying the start date and time of the period
        desired. Must be in 'YYYYmmdd-HHMMSS' format. Defaults
        to None, in which case all data from the first file through 
        endDT will be used.
    endDT : string, optional
        String specifying the end date and time of the period
        desired. Must be in 'YYYYmmdd-HHMMSS' format. Defaults
        to None, in which case all data from strtDT through the
        end of available files will be used.
    llCrds : tuple
        Tuple of (lat,lon) at the lower left corner of the desired
        domain. Defaults to None, in which case the entire dataset
        is returned. Both llCrds and urCrds must have valid values
        to be used.
    urCrds : tuple
        Tuple of (lat,lon) at the upper right corner of the desired
        domain. Defaults to None, in which case the entire dataset
        is returned. Both llCrds and urCrds must have valid values
        to be used.
        
    Returns
    -------
    mrms_out : dict
        Dictionary containing MRMS data bounded by the start/end times (if given)
        and within the lat/lon limits (if given).

    """

    if (strtDT is not None and len(strtDT) == 15) and (endDT is not None and len(endDT) == 15):
        sDT = dt.strptime(strtDT,'%Y%m%d-%H%M%S')
        eDT = dt.strptime(endDT,'%Y%m%d-%H%M%S')
    elif strtDT is not None and len(strtDT) == 15:
        sDT = dt.strptime(strtDT,'%Y%m%d-%H%M%S')
        eDT = None
    elif endDT is not None and len(endDT) == 15:
        sDT = None
        eDT = dt.strptime(endDT,'%Y%m%d-%H%M%S')
    else:
        sDT = None
        eDT = None
        print('Extracting all data. Ensure strtDT and/or endDT were properly formatted if passed as arguments.')
        
    
    if mrmsDir[-1] == '/':
        mrmsDir = mrmsDir[:-1]
    mrmsFiles = sorted(glob('{}/*.nc'.format(mrmsDir)))
    dtStrs = np.asarray([dt.strptime(mDT[-18:-3],'%Y%m%d-%H%M%S') for mDT in mrmsFiles])    

    
    if sDT is not None:
        sIX = np.squeeze(np.where(dtStrs == sDT))
    else:
        sIX = None
        
    if eDT is not None:
        eIX = np.squeeze(np.where(dtStrs == eDT)) + 1
    else:
        eIX = None
        

    da_DT = xr.DataArray(dtStrs[sIX:eIX],name='DateTime',dims='DateTime',coords={'DateTime': dtStrs[sIX:eIX]})
        
    
    mrms_xrs = []
    for mFile in mrmsFiles[sIX:eIX]:
        mrmsSngl = xr.open_dataset(mFile,mask_and_scale=True)
        mrms_xrs.append(mrmsSngl)
        
    mrms_All = xr.concat(mrms_xrs,dim=da_DT)
    
    
    if llCrds is not None and urCrds is not None:
        lat1 = llCrds[0]
        lat2 = urCrds[0]
        lon1 = urCrds[1] + 360
        lon2 = llCrds[1] + 360
        mrms_sel = mrms_All.sel(lon_0=slice(lon2,lon1)).sel(lat_0=slice(lat2,lat1))
    else:
        mrms_sel = mrms_All
        
    
    
    
    shsr_out = mrms_sel.SeamlessHSR_P0_L102_GLL0.to_masked_array()
    lat_out = mrms_sel.lat_0.to_masked_array()
    lon_out = mrms_sel.lon_0.to_masked_array()
    dt_out1 = mrms_sel.DateTime.to_masked_array()
    
    mTS = (dt_out1 - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
    dt_out = np.asarray([dt.utcfromtimestamp(t) for t in mTS])
    
    mrms_out = {'mSHSR': shsr_out, 'mLat': lat_out, 'mLon': lon_out, 'mDT': dt_out}
    
    return mrms_out