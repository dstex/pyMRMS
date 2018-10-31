"""
graph
======

Tools for calculating and/or transforming coordinates for various
plot annotations.

    plotSHSR
    

"""

import numpy as np
from matplotlib import pyplot as plt
from pyart import graph
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
    
    ax.text(-0.07, 0.5, 'Latitude', va='bottom', ha='center',
            rotation='vertical', rotation_mode='anchor', size=20,
            transform=ax.transAxes)
    ax.text(0.5, -0.12, 'Longitude', va='bottom', ha='center',
            rotation='horizontal', rotation_mode='anchor', size=20,
            transform=ax.transAxes)
            
    return fig,ax,grd,proj


def plotSHSR(mLon,mLat,mSHSR,mDT,saveFigs=True,savePath='',fig=None,ax=None,grd=None,proj=None):
    """
    Function for plotting the MRMS Seamless Hybrid Scan Reflectivity (SHSR)
    atop a cartopy map projection.
    
    
    Parameters
    ----------
    mLon,mLat : 1D float arrays
        Arrays containing the lat/lon of the MRMS data to be plotted.
    mSHSR : 2D float array
        2D array of MRMS SHSR data. This function only plots a single period
        at a time (do not pass in the 3D version of mSHSR)
    mDT : datetime
        Datetime of the MRMS data to be plotted.
    saveFigs : bool, optional
        Determines whether figures are saved. Defaults to True.
    savePath : string, optional
        Location where figures will be saved to. If output directory doesn't exist,
        it (and all parents needed) will be created. Defaults to empty string.
    fig,ax,grd,proj : handles, optional
        Figure and axes handles for an existing plot. If all are None (default),
        this function will generate these handles through a call to initMap.
    
    Returns
    -------
    fig,ax,grd,proj : handles
        Figure and axes handles used for further customizing the plot.
    """
    
    vMin = -4 # Max and min values to plot
    vMax = 60
    cmap = graph.cm_colorblind.HomeyerRainbow
    cmap.set_under('w')
    bounds = np.linspace(vMin,vMax,(np.abs(vMin)+np.abs(vMax)+1))
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    # norm = None
    cbarStr = 'Reflectivity (dBZ)'
    
    
    llCrds = (np.min(mLat),np.min(mLon))
    urCrds = (np.max(mLat),np.max(mLon))
    
    if fig is None or ax is None or grd is None or proj is None:
        fig,ax,grd,proj = initMap(llCrds,urCrds)
        
    plt.pcolormesh(mLon,mLat,mSHSR,cmap=cmap,norm=norm,vmin=vMin,vmax=vMax,transform=proj)
    
    cb = plt.colorbar(shrink=0.7, pad = 0.01, aspect=25)
    cb.set_label('Reflectivity (dBZ)',size=17)
    cb.ax.tick_params(labelsize=14)
    
    plt.title('MRMS Seamless Hybrid Scan Reflectivity',size=22)