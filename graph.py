"""
graph
======

Tools for calculating and/or transforming coordinates for various
plot annotations.

    initMap
    calcSwpLocs
    plotSHSR
    
    

"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import os
import sys
from datetime import datetime as dt
from pyart import graph as pag
import cartopy.crs as ccrs
import cartopy.feature as cf
from pyproj import Proj as Proj
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


def initMap(llCrds,urCrds,figsize=(16,16)):
    """
    Function for initializing a figure using a cartopy map projection.
    
    
    Parameters
    ----------
    llCrds : tuple
        Tuple of (lat,lon) at the lower left corner of the desired
        domain.
    urCrds : tuple
        Tuple of (lat,lon) at the upper right corner of the desired
        domain.
    figsize : tuple, optional
        Tuple specifying the size of the output figure in inches (W,H).
        Defaults to (16,16).
    
    Returns
    -------
    fig,ax,grd,proj : handles
        Figure and axes handles used for adding to and further customizing the plot.
    """
    
    proj = ccrs.PlateCarree()

    states_provinces = cf.NaturalEarthFeature(
            category='cultural',name='admin_1_states_provinces_lines',
            scale='50m',facecolor='none')
    
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111,projection=proj)
    ax.coastlines()
    ax.add_feature(cf.OCEAN)
    ax.add_feature(cf.LAKES)
    ax.add_feature(cf.BORDERS)
    ax.add_feature(states_provinces, edgecolor='gray')
    ax.set_extent([llCrds[1],urCrds[1],llCrds[0],urCrds[0]])

    grd = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    grd.xlabels_top = False
    grd.ylabels_right = False
    grd.xformatter = LONGITUDE_FORMATTER
    grd.yformatter = LATITUDE_FORMATTER
    grd.xlabel_style = {'size': 18}
    grd.ylabel_style = {'size': 18}
    
    ax.text(-0.1, 0.5, 'Latitude', va='bottom', ha='center',
            rotation='vertical', rotation_mode='anchor', size=20,
            transform=ax.transAxes)
    ax.text(0.5, -0.07, 'Longitude', va='bottom', ha='center',
            rotation='horizontal', rotation_mode='anchor', size=20,
            transform=ax.transAxes)
            
    return fig,ax,grd,proj


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


def plotSHSR(mLon,mLat,mSHSR,mDT,plotFltTrk=False,flLon=None,flLat=None,flDT=None,
                projID='',fType='png',fig=None,ax=None,grd=None,proj=None):
    """
    Function for plotting the MRMS Seamless Hybrid Scan Reflectivity (SHSR)
    atop a cartopy map projection.
    
    
    Parameters
    ----------
    mLon,mLat : 1D float arrays
        Arrays containing the lat/lon of the MRMS data to be plotted.
    mSHSR : 3D float array
        3D array of MRMS SHSR data of shape (time,lon,lat)
    mDT : 1D datetime array
        Datetime array of the MRMS data to be plotted.
    plotFltTrk : bool, optional
        Determines whether to plot the P-3 flight track atop the SHSR plot. If True,
        the associated flight-level data must be specified. Defaults to False.
    flLon,flLat : 1D float arrays, conditionally optional
        Arrays containing the lat/lon of the flight-level data to be plotted.
    flDT : 1D datetime array, conditionally optional
        Datetime array for the flight-level data to be plotted.
    projID : string, optional
        String specifying project name. Will be included in plot titles and filenames.
    fType : string, optional
        String specifying the output figure file type. Defaults to 'png'
    fig,ax,grd,proj : handles, optional
        Figure and axes handles for an existing plot. If all are None (default),
        this function will generate these handles through a call to initMap.
    
    Returns
    -------
    fig,ax,grd,proj : handles
        Figure and axes handles used for further customizing the plot.
    """
    ### Determine map boundaries ###
    llCrds = (np.min(mLat),np.min(mLon))
    urCrds = (np.max(mLat),np.max(mLon))
    
    
    ### Set contour plot options ###
    vMin = -4 # Max and min values to plot
    vMax = 60
    cmap = pag.cm_colorblind.HomeyerRainbow
    cmap.set_under('w')
    bounds = np.linspace(vMin,vMax,(np.abs(vMin)+np.abs(vMax)+1))
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    # norm = None
    cbarStr = 'Reflectivity (dBZ)'
    
    
    ### Initialize a map plot if one was not specified ###
    if fig is None or ax is None or grd is None or proj is None:
        fig,ax,grd,proj = initMap(llCrds,urCrds)
        
    
    ### Plot the MRMS SHSR data ###
    plt.pcolormesh(mLon,mLat,mSHSR,cmap=cmap,norm=norm,vmin=vMin,vmax=vMax,transform=proj)
    
    cb = plt.colorbar(shrink=0.7, pad = 0.01, aspect=25)
    cb.set_label('Reflectivity (dBZ)',size=17)
    cb.ax.tick_params(labelsize=14)
    
    dtTtlStr = dt.strftime(dt.utcfromtimestamp(int((mDT.astype(dt))*1e-9)),'%Y/%m/%d - %H:%M')
    dtSvStr = dt.strftime(dt.utcfromtimestamp(int((mDT.astype(dt))*1e-9)),'%Y%m%d-%H%M')
    
    plt.title('MRMS Seamless Hybrid Scan Reflectivity\n{} - {}'.format(projID,dtTtlStr),size=22)
                
    
    ### Plot flight track if desired ###
    if plotFltTrk:
        if flLon is None or flLat is None or flDT is None:
            sys.exit('flLon, flLat, and flDT must all receive arguments when plotFltTrk is True')
        else:
            # Get difference in minutes so that we know how many minutes of flight track
            # to plot on each given radar composite
            compDiff = np.ones(len(radDT))
            for ix in range(1,len(radDT)):
                compDiff[ix-1] = (radDT[ix]-radDT[ix-1]).total_seconds()/60


            ## Plotting
            pastTrack = False # We'll set this to true the first time the flight track is in the domain
    
            ## Find the closest FL data index correlating to the current radar grid time
            domMatch = min(flDT, key=lambda x: abs(x - mDT))
            flDomIx = np.squeeze(np.where(flDT == domMatch))

            crntFLlat = flLat[flDomIx]
            crntFLlon = flLon[flDomIx]
            crntFLheading = flHeading[flDomIx]
            plotFLlat = flLat[0:flDomIx]
            plotFLlon = flLon[0:flDomIx]


            # Expand domain by half degree all around - this will help ensure sweep locations
            # are plotted even if plane itself is out of domain
            mapLatMin = llCrds[0]-0.5
            mapLatMax = urCrds[0]+0.5
            mapLonMin = llCrds[1]-0.5
            mapLonMax = urCrds[1]+0.5

            # We want to replot a given composite for every minute we have flight (unless plotFltStatic 
            # is True) data until the next composite
            # We'll only make a plot for every minute when the plane is actually in the 
            # air and within 0.5 deg of the domain
            if ((flDT[0]<= radDT[ix] <= flDT[-1]) & (gridLatMin <= crntFLlat <= gridLatMax) & (gridLonMin <= crntFLlon <= gridLonMax)):
                inDomain = True
                inLoop = int(compDiff[ix])
                if inLoop > 9:
                    inLoop = 1;
                    print('Difference between grids exceeded 9 minutes - skipping 1-min flt track')
                pastTrack = True
            else:
                inDomain = False
                inLoop = 1
    
    
    return fig,ax,grd,proj
    
    
def plotFT(flLon,flLat,flDT,mDT,mllCrds,murCrds,fig=None,ax=None,grd=None,proj=None):
    """
    Function for plotting the P-3 flight track, current location, and TDR
    fore/aft beam extents.
    
    
    Parameters
    ----------
    flLon,flLat : 1D float arrays
        Arrays containing the lat/lon of the flight-level data to be plotted.
    flDT : 1D datetime array
        Datetime array for the flight-level data to be plotted.
    mDT : datetime
        Datetime of the MRMS data we're plotting the flight track on.
    mllCrds,murCrds : tuples
        Tuples containing the (lat,lon) of the lower left corner (mllCrds) and
        upper right corner (murCrds) of the chosen MRMS SHSR domain. These are
        used as a reference to determine when the P-3 (or TDR beams) are within 
        the plot domain.
    fig,ax,grd,proj : handles, optional
        Figure and axes handles for an existing plot. If all are None (default),
        this function will generate these handles through a call to initMap.
    
    Returns
    -------
    fig,ax,grd,proj : handles
        Figure and axes handles used for further customizing the plot.
    """
    
    