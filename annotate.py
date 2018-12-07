"""
annotate
======

Tools for calculating and/or transforming coordinates for various
plot annotations.

    calcSwpLocs
    rangeRingCalc
    

"""

import numpy as np


def calcSwpLocs(lat,lon,hdng,radRange=47966.79):
    """
    Calculates the lat,lon values at the terminus of the fore and aft
    TDR beams, both left and right of the plane in the horizontal.
    This is used when creating the so-called "Dragonfly" for visualizing
    radar coverage for a P-3 mission.
    
    
    Parameters
    ----------
    lat,lon : floats
        Lat and lon of the P-3 at a given moment in time, in decimal degrees
    hdng : float
        Heading (direction aircraft is flying) at a given moment in time (in degrees)
    radRange : float, optional
        Maximum unambiguous range of the TDR in meters. Defaults to value used by the
        TDR during the 2017 VSE project.
    
    Returns
    -------
    foreSwpLocs,aftSwpLocs : dicts
        Dictionaries containing the left and right lat and lon values
        for the fore and after TDR beams.

    """
    
    # Sweep unambiguous range (for NOAA TDR - 2017 VSE, 47966.79 m) converted to radians
    # (relative to spherical Earth, of radius 6372797.6 m)
    beamDistR = radRange/6372797.6
    
    # Calculate hdngs of both right and left (relative to motion) fore and aft sweeps
    aftLeftHead = hdng - (90+20)
    aftRightHead = hdng + (90+20)
    foreLeftHead = hdng - (90-20)
    foreRightHead = hdng + (90-20)
    if aftLeftHead < 0:
        aftLeftHead += 360
    if aftRightHead < 0:
        aftRightHead += 360
    if foreLeftHead < 0:
        foreLeftHead += 360
    if foreRightHead < 0:
        foreRightHead += 360
        
    # Convert hdngs and coords to radians
    latR = np.deg2rad(lat)
    lonR = np.deg2rad(lon)
    aftLeftHeadR = np.deg2rad(aftLeftHead)
    aftRightHeadR = np.deg2rad(aftRightHead)
    foreLeftHeadR = np.deg2rad(foreLeftHead)
    foreRightHeadR = np.deg2rad(foreRightHead)
    
    # Calculate destination coordinates and convert back to degrees
    lat2aftLeftR = np.arcsin(np.sin(latR)*np.cos(beamDistR) + np.cos(latR)*np.sin(beamDistR)*np.cos(aftLeftHeadR))
    lon2aftLeftR = lonR + np.arctan2(np.sin(aftLeftHeadR)*np.sin(beamDistR)*np.cos(latR), np.cos(beamDistR) - np.sin(latR)*np.sin(lat2aftLeftR))
    lat2aftRightR = np.arcsin(np.sin(latR)*np.cos(beamDistR) + np.cos(latR)*np.sin(beamDistR)*np.cos(aftRightHeadR))
    lon2aftRightR = lonR + np.arctan2(np.sin(aftRightHeadR)*np.sin(beamDistR)*np.cos(latR), np.cos(beamDistR) - np.sin(latR)*np.sin(lat2aftRightR))
    lat2foreLeftR = np.arcsin(np.sin(latR)*np.cos(beamDistR) + np.cos(latR)*np.sin(beamDistR)*np.cos(foreLeftHeadR))
    lon2foreLeftR = lonR + np.arctan2(np.sin(foreLeftHeadR)*np.sin(beamDistR)*np.cos(latR), np.cos(beamDistR) - np.sin(latR)*np.sin(lat2foreLeftR))
    lat2foreRightR = np.arcsin(np.sin(latR)*np.cos(beamDistR) + np.cos(latR)*np.sin(beamDistR)*np.cos(foreRightHeadR))
    lon2foreRightR = lonR + np.arctan2(np.sin(foreRightHeadR)*np.sin(beamDistR)*np.cos(latR), np.cos(beamDistR) - np.sin(latR)*np.sin(lat2foreRightR))
    
    latAL = np.rad2deg(lat2aftLeftR)
    lonAL = np.rad2deg(lon2aftLeftR)
    latAR = np.rad2deg(lat2aftRightR)
    lonAR = np.rad2deg(lon2aftRightR)

    aftSwpLocs = {'latL': np.rad2deg(lat2aftLeftR),
                  'lonL': np.rad2deg(lon2aftLeftR),
                  'latR': np.rad2deg(lat2aftRightR),
                  'lonR': np.rad2deg(lon2aftRightR)}
    
    foreSwpLocs = {'latL': np.rad2deg(lat2foreLeftR),
                   'lonL': np.rad2deg(lon2foreLeftR),
                   'latR': np.rad2deg(lat2foreRightR),
                   'lonR': np.rad2deg(lon2foreRightR)}
    
    return aftSwpLocs,foreSwpLocs
    

def rangeRingCalc(lat1,lon1,maxRngKm):
    """
    Calculates lat/lon pairs at some distance from a given point.

    Parameters
    ----------
    lat1,lon1 : float
        Lon/lat coordinates (degrees) of the initial point.
    maxRngKm : float
        Max range of radar (in KM).
     

    Returns
    -------
    latCirc,lonCirc : 1D arrays of floats
        Lat/lon (in degrees) following some radar max range
    """

    lat1r = np.deg2rad(lat1)
    lon1r = np.deg2rad(lon1)
    
    rEarth = 6371
    
    
    latCirc = []
    lonCirc = []
    for az in range(0,360):
        hdngRad = np.deg2rad(az)
        
        lat2r = np.arcsin( np.sin(lat1r)*np.cos(maxRngKm/rEarth) + np.cos(lat1r)*np.sin(maxRngKm/rEarth)*np.cos(hdngRad) )
        lon2r = lon1r + np.arctan2( np.sin(hdngRad)*np.sin(maxRngKm/rEarth)*np.cos(lat1r), np.cos(maxRngKm/rEarth)-np.sin(lat1r)*np.sin(lat2r) )
        
        latCirc.append(np.rad2deg(lat2r))
        lonCirc.append(np.rad2deg(lon2r))
    
    
    return latCirc,lonCirc