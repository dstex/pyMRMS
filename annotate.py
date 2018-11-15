"""
annotate
======

Tools for calculating and/or transforming coordinates for various
plot annotations.

    rangeRingCalc

"""

import numpy as np


    

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