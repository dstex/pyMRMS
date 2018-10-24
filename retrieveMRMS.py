"""
retrieveMRMS
======

Download MRMS Seamless Hybrid Scan Reflectivity data.

    fetch

"""
import urllib.request as request
from datetime import datetime as dt
import datetime
from subprocess import call
import os

def fetch(strtDTstr,endDTstr,saveDir):
    """
    Download and unzip MRMS files from the Iowa State data 
    archive between a given start and end date/time. These 
    files exist at 2-minute intervals.
    
    
    Parameters
    ----------
    strtDTstr : string
        String specifying the start date and time of the period
        desired. Must be in 'YYYYmmdd-HHMMSS' format.
    endDTstr : string
        String specifying the end date and time of the period
        desired. Must be in 'YYYYmmdd-HHMMSS' format.
    saveDir : string
        String specifying the path where the downloaded files
        are to be saved. If the given path does not exist, it
        will be created.

    """


    # strtDTstr = '20150620-040000'
    # endDTstr = '20150620-050000'
    sDT = dt.strptime(strtDTstr,'%Y%m%d-%H%M%S')
    eDT = dt.strptime(endDTstr,'%Y%m%d-%H%M%S')

    # saveDir = '/Users/danstechman/Desktop/MRMS/'
    rmtBase = 'http://mtarchive.geol.iastate.edu/'

    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    tStep = sDT
    while tStep <= eDT:
        locPath = '{}SeamlessHSR_00.00_{:%Y%m%d-%H%M}00.grib2.gz'.format(saveDir,tStep)
        rmtPath = '{}{:%Y/%m/%d}/mrms/ncep/SeamlessHSR/SeamlessHSR_00.00_{:%Y%m%d-%H%M}00.grib2.gz'.format(rmtBase,tStep,tStep)
        
        if os.path.isfile(locPath[:-3]):
            print('Skipping download of {:%Y%m%d - %H:%M}. File already exists locally.'.format(tStep))
            if os.path.isfile('{}nc'.format(locPath[:-8])):
                print('\tAlso skipping grib2nc conversion as netCDF file already exists.')
            else:
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