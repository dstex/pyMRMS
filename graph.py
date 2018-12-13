"""
graph
======

Tools for calculating and/or transforming coordinates for various
plot annotations.

    initMap
    plotSHSR
    
    

"""

import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import sys
import os
from datetime import datetime as dt
import datetime
from pyart import graph as pag
import cartopy.crs as ccrs
import cartopy.feature as cf
from pyproj import Proj as Proj
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from tqdm import tqdm
import math
from adjustText import adjust_text

from . import annotate


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
    ax.text(0.5, -0.1, 'Longitude', va='bottom', ha='center',
            rotation='horizontal', rotation_mode='anchor', size=20,
            transform=ax.transAxes)
            
    return fig,ax,grd,proj


def plotAssets(aName,aLon,aLat,aMark,aCol,proj):
    """
    Function for plotting the asset locations.
    
    
    Parameters
    ----------
    
    aName : string
        String to label asset with.
    aLon,aLat : floats
        Floats specifying the coordinates of a given asset location.
    aMark : string
        Any valid matplotlib marker specifier.
    aCol : string/float/tuple
        Any single valid matplotlib color specifier.
    proj : map projection handle
        Cartopy map projection handle for transforming asset coordinates.

    """
    
    plt.plot(aLon,aLat,marker=aMark,markerfacecolor=aCol,markersize=15,markeredgecolor='w',markeredgewidth=1,transform=proj)
    txt = plt.annotate(aName,xy=(aLon,aLat),zorder=10,fontsize=13,color='w',
                 bbox=dict(boxstyle="round", fc=aCol,alpha=0.6,pad=0.01))
    
    return txt
    

def plotRangeRing(rLon,rLat,rRange,rName,proj,linestyle='-'):
    """
    Function for plotting the asset locations.
    
    
    Parameters
    ----------
    
    rLon,rLat : floats
        Floats specifying the coordinates of a given radar.
    rRange : float
        Max unambiguous range of given radar in km.
    rCol : string/float/tuple
        Any single valid matplotlib color specifier.
    proj : map projection handle
        Cartopy map projection handle for transforming radar/range coordinates.

    """
    
    latRing,lonRing = annotate.rangeRingCalc(rLat,rLon,rRange)
    plt.plot(lonRing,latRing,linewidth=2,color='w',transform=proj)
    plt.plot(lonRing,latRing,linewidth=1,transform=proj,linestyle=linestyle)
    
    

def plotSHSR(mLon,mLat,mSHSR,mDT,plotFltTrk=False,plotSwpLocs=True,plotRng=False,
                flLon=None,flLat=None,flDT=None,flHdng=None,fltTrkMin=20,fadeFltTrk=None,fadeFac=10,
                projID='',fType='png',saveDir='',figsize=(16,16)):
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
    
    Returns
    -------
    fig,ax,grd,proj : handles
        Figure and axes handles used for further customizing the plot.
    """
    
    if plotFltTrk:
        if flLon is None or flLat is None or flDT is None:
            sys.exit('flLon, flLat, and flDT must all receive arguments when plotFltTrk is True')
        fltTrkSec = fltTrkMin*60
    
    ### Determine map boundaries ###
    llCrds = (np.min(mLat),np.min(mLon-360))
    urCrds = (np.max(mLat),np.max(mLon-360))
    
    
    ### Set contour plot options ###
    vMin = -4 # Max and min values to plot
    vMax = 60
    cmap = pag.cm_colorblind.HomeyerRainbow
    cmap.set_under('w')
    bounds = np.linspace(vMin,vMax,(np.abs(vMin)+np.abs(vMax)+1))
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    # norm = None
    cbarStr = 'Reflectivity (dBZ)'
    
    
    for t in tqdm(range(len(mDT)),dynamic_ncols=True):
        mDTt = mDT[t]
        mSHSRt = mSHSR[t]
                
    
        ### Plot flight track if desired ###
        if plotFltTrk:
            pastTrack = False # We'll set this to true the first time the flight track is in the domain

            ## Find the closest FL data index correlating to the current radar grid time
            domMatch = min(flDT, key=lambda x: abs(x - mDTt))
            flDomIx = np.squeeze(np.where(flDT == domMatch))

            crntFLlat = flLat[flDomIx]
            crntFLlon = flLon[flDomIx]

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
            if ((flDT[0]<= mDTt <= flDT[-1]) & (mapLatMin <= crntFLlat <= mapLatMax) & (mapLonMin <= crntFLlon <= mapLonMax)):
                inDomain = True
                if t < len(mDT)-1:
                    nxtDomMatch = min(flDT, key=lambda x: abs(x - mDT[t+1]))
                    nxtFLdomIX = np.squeeze(np.where(flDT == nxtDomMatch))
                    # If the plane is stationary, we don't need to plot every minute
                    if math.isclose(crntFLlat,flLat[nxtFLdomIX],abs_tol=0.00005) and math.isclose(crntFLlon,flLon[nxtFLdomIX],abs_tol=0.00001):
                        inLoop = 1
                    else:
                        inLoop = 2 # Num minutes between MRMS composites to plot flt track over - can be dynamic if desired
                else:
                    inLoop = 1
                pastTrack = True
            else:
                inDomain = False
                inLoop = 1
                
            for iz in range(inLoop):
                pltT = mDTt + datetime.timedelta(minutes=iz)
                
                ### Initialize a map plot ###
                fig,ax,grd,proj = initMap(llCrds,urCrds,figsize=figsize)
        
    
                ### Plot the MRMS SHSR data ###
                plt.pcolormesh(mLon,mLat,mSHSRt,cmap=cmap,norm=norm,vmin=vMin,vmax=vMax,transform=proj)
    
                cb = plt.colorbar(shrink=0.7, pad = 0.01, aspect=25)
                cb.set_label('Reflectivity (dBZ)',size=17)
                cb.ax.tick_params(labelsize=14)

                dtTtlStr = dt.strftime(pltT,'%Y/%m/%d - %H:%M')
                dtSvStr = dt.strftime(pltT,'%Y%m%d-%H%M')
    
                plt.title('MRMS Seamless Hybrid Scan Reflectivity\n{} - {}'.format(projID,dtTtlStr),size=22)
    
    
                # If our radar grid time is within bounds of the flight, plot the flight track
                # This will automatically be False if plot_FltTrk is False
                if pastTrack:
                    domMatch = min(flDT, key=lambda x: abs(x - pltT))
                    flDomIx = np.squeeze(np.where(flDT == domMatch))

                    crntFLlat = flLat[flDomIx]
                    crntFLlon = flLon[flDomIx]
                    crntFLheading = flHdng[flDomIx]
                    
                    if fadeFltTrk == 'fade':
                        fadeIntvl = int(fltTrkSec/fadeFac)
                        fadeDiff = flDomIx-fltTrkSec
                        
                        numIntvls = int(fadeDiff/fadeIntvl)
                        remIntvl = fadeDiff%fadeIntvl
                    
                        if fadeDiff < fltTrkSec:
                            plotFLlat = flLat[0:flDomIx+1]
                            plotFLlon = flLon[0:flDomIx+1]
                            
                            plt.plot(plotFLlon,plotFLlat,linewidth=2.25,color='k',transform=proj)
                            plt.plot(plotFLlon,plotFLlat,linewidth=0.75,color='w',transform=proj)
                            
                        elif fadeDiff >= fltTrkSec:
                            plotFLlat = flLat[fadeDiff:flDomIx+1]
                            plotFLlon = flLon[fadeDiff:flDomIx+1]
                            plt.plot(plotFLlon,plotFLlat,linewidth=2.25,color='k',transform=proj)
                            plt.plot(plotFLlon,plotFLlat,linewidth=0.75,color='w',transform=proj)
                            
                            if numIntvls > fadeFac:
                                numIntvls_orig = numIntvls
                                numIntvls = fadeFac
                                remIntvl = 0
                            for itvl in np.arange(1,numIntvls+1):
                                ixStrt = fadeDiff-(fadeIntvl*itvl)
                                ixEnd = fadeDiff-(fadeIntvl*(itvl-1))
                                
                                plotFLlat = flLat[ixStrt:ixEnd]
                                plotFLlon = flLon[ixStrt:ixEnd]
                                plt.plot(plotFLlon,plotFLlat,linewidth=2.25,color='k',alpha=(1.05-(itvl/fadeFac)),transform=proj)
                                plt.plot(plotFLlon,plotFLlat,linewidth=0.75,color='w',alpha=(1.05-(itvl/fadeFac)),transform=proj)
                                if (itvl == numIntvls) and (remIntvl != 0):
                                    ixEnd = remIntvl
                                    ixStrt = 0
                                    
                                    plotFLlat = flLat[ixStrt:ixEnd]
                                    plotFLlon = flLon[ixStrt:ixEnd]
                                    plt.plot(plotFLlon,plotFLlat,linewidth=2.25,color='k',alpha=(1.05-(itvl/fadeFac)),transform=proj)
                                    plt.plot(plotFLlon,plotFLlat,linewidth=0.75,color='w',alpha=(1.05-(itvl/fadeFac)),transform=proj)
                    
                    else:
                        if fadeFltTrk is None:
                            plotFLlat = flLat[0:flDomIx+1]
                            plotFLlon = flLon[0:flDomIx+1]
                        elif fadeFltTrk == 'cut':
                            if fadeDiff >= fltTrkSec:
                                plotFLlat = flLat[fadeDiff:flDomIx+1]
                                plotFLlon = flLon[fadeDiff:flDomIx+1]
                            else:
                                plotFLlat = flLat[0:flDomIx+1]
                                plotFLlon = flLon[0:flDomIx+1]
                        
                        plt.plot(plotFLlon,plotFLlat,linewidth=2.25,color='k',transform=proj)
                        plt.plot(plotFLlon,plotFLlat,linewidth=0.75,color='w',transform=proj)
                    
            
                    # If the plane is within 0.5 deg box of domain, plot at least the plane
                    # location, and optionally, the sweep locations
                    if inDomain:
                        if plotSwpLocs:
                            aftSwpLocs,foreSwpLocs = annotate.calcSwpLocs(crntFLlat,crntFLlon,crntFLheading)
                            

                            # Create x/y coordinate arrays with the start/end values for each sweep
                            xAL = [crntFLlon,aftSwpLocs['lonL']]
                            xAR = [crntFLlon,aftSwpLocs['lonR']]
                            xFL = [crntFLlon,foreSwpLocs['lonL']]
                            xFR = [crntFLlon,foreSwpLocs['lonR']]
                            yAL = [crntFLlat,aftSwpLocs['latL']]
                            yAR = [crntFLlat,aftSwpLocs['latR']]
                            yFL = [crntFLlat,foreSwpLocs['latL']]
                            yFR = [crntFLlat,foreSwpLocs['latR']]

                            # Plot lines indicating sweep locations
                            plt.plot(xAL,yAL,linewidth=3,color='w',zorder=15,transform=proj)
                            plt.plot(xAL,yAL,linewidth=1.75,color='MediumBlue',zorder=15,transform=proj)
                            plt.plot(xAR,yAR,linewidth=3,color='w',zorder=15,transform=proj)
                            plt.plot(xAR,yAR,linewidth=1.75,color='MediumBlue',zorder=15,transform=proj)
                            plt.plot(xFL,yFL,linewidth=3,color='w',zorder=15,transform=proj)
                            plt.plot(xFL,yFL,linewidth=1.75,color='DarkRed',zorder=15,transform=proj)
                            plt.plot(xFR,yFR,linewidth=3,color='w',zorder=15,transform=proj)
                            plt.plot(xFR,yFR,linewidth=1.75,color='DarkRed',zorder=15,transform=proj)
                    
                        # Plot marker at current location of plane
                        plt.plot(crntFLlon,crntFLlat,'wo',markersize=12,zorder=15,transform=proj)
                        plt.plot(crntFLlon,crntFLlat,'ko',markersize=8,zorder=15,transform=proj)
                        
                assets = annotate.readYaml('/Users/danstechman/PersonalGoogleDrive/Code/radar/pyMRMS/assets_template.yml')
                w88Ds = annotate.readYaml('/Users/danstechman/PersonalGoogleDrive/Code/radar/pyMRMS/w88Ds.yml')
                
                aNames = assets['assetNames']
                aLons = assets['assetLons']
                aLats = assets['assetLats']
                aSymbs = assets['assetSymbs']
                aCols = assets['assetColors']
                rNames = assets['radarNames']
                rLons = assets['radarLons']
                rLats = assets['radarLats']
                rRanges = assets['radarRanges']
                rSymbs = assets['radarSymbs']
                rCols = assets['radarColors']
                wNames = assets['wsrNames']
                wCol = assets['wsrColor']
                wSymb = assets['wsrSymb']
            
                txts = []
                for ia in range(len(aNames)):
                    txts.append(plotAssets(aNames[ia],aLons[ia],aLats[ia],aSymbs[ia],aCols[ia],proj))
                for ir in range(len(rNames)):
                    txts.append(plotAssets(rNames[ir],rLons[ir],rLats[ir],rSymbs[ir],rCols[ir],proj))
                for iw in range(len(wNames)):
                    wName = wNames[iw]
                    wLon = w88Ds[wName]['coords'][1]
                    wLat = w88Ds[wName]['coords'][0]
                    txts.append(plotAssets(wName,wLon,wLat,wSymb,wCol,proj))
                    
                if plotRng:
                    for iGR in range(len(rNames)):
                        if iGR%2 == 0:
                            plotRangeRing(rLons[iGR],rLats[iGR],rRanges[iGR],rNames[iGR],proj)
                        else:
                            plotRangeRing(rLons[iGR],rLats[iGR],rRanges[iGR],rNames[iGR],proj,linestyle='--')
                    for iGR in range(len(wNames)):
                        wName = wNames[iGR]
                        wLon = w88Ds[wName]['coords'][1]
                        wLat = w88Ds[wName]['coords'][0]
                        wRange = w88Ds[wName]['range']
                        if iGR%2 == 0:
                            plotRangeRing(wLon,wLat,wRange,wName,proj)
                        else:
                            plotRangeRing(wLon,wLat,wRange,wName,proj,linestyle='--')
                        
                    #plt.legend(fontsize=14,loc='upper right')
                
                adjust_text(txts,expand_text=(1.2, 1.2), expand_points=(1.2, 1.2),force_points=(0.5,0.5))

            
                if not os.path.exists(saveDir):
                    os.makedirs(saveDir)
                fig.savefig('{}/{}_mSHSR.{}'.format(saveDir,dtSvStr,fType),bbox_inches='tight',dpi=100)
                plt.close('all')