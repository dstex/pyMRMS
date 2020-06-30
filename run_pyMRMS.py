from pyMRMS import retrieveMRMS,readFL,annotate,graph


projID = 'VSE'
tdrRng = 47966.79

plotMRMS = True
getMRMS = False


pltFT = True

plotAsts = False
plotRng = False
pltMultRng = False

plotStrmRpts = True
rprtWindow = -300

animAll = False

fltIntv = 60 # Plotting interval for flight track (seconds)

plotDom = False
domOrig = (34.0, -87.75)#(34.15, -87.3) 
domX = 240
domY = 170



strtMRMS = '20170430-190000'
endMRMS = '20170430-192000'

mrmsSaveDir = '/Volumes/Pegasus/projData/vortexse17/MRMS/20170430/1800-2200Z'
mrmsPlotDir = '/Volumes/Pegasus/projData/vortexse17/MRMS/20170430/plots_StrmRpts'

flFile = '/Volumes/Pegasus/projData/vortexse17/flightLevel/20170430H1_AC.nc'

astsPrmsF = '/Volumes/Pegasus/projData/vortexse17/MRMS/20170430_assets_soundingLocs.yml'
strmRptsF = '/Volumes/Pegasus/projData/vortexse17/MRMS/20170430_stormReports.yml'

minLat = 33.1
minLon = -88.0
maxLat = 36.4
maxLon = -84.0



strtFL = strtMRMS
endFL = endMRMS

if getMRMS:
    retrieveMRMS.fetch(strtMRMS,endMRMS,mrmsSaveDir)



if plotMRMS:
    if pltFT:
        print('Reading flight-level data...')
        flData = readFL.getP3(flFile)#,endDT=endFL)
   
    print('Extracting MRMS data...')
    mData = retrieveMRMS.extract(mrmsSaveDir,strtDT=strtMRMS,endDT=endMRMS,llCrds=(minLat,minLon),urCrds=(maxLat,maxLon))
    
    print('Plotting...')
    graph.plotSHSR(mData['mLon'],mData['mLat'],mData['mSHSR'],mData['mDT'],radRange=tdrRng,
                   plotAsts=plotAsts,prmsFile=astsPrmsF,plotRng=plotRng,plotFltTrk=pltFT,multRings=pltMultRng,
                   flLon=flData['flLon'],flLat=flData['flLat'],flDT=flData['flDT'],
                   flHdng=flData['flHdng'],fltIntv=fltIntv,fadeFltTrk='fade',
                   projID=projID,plotDom=plotDom,domOrig=domOrig,domX=domX,domY=domY,
                   fType='png',saveDir=mrmsPlotDir,figsize=(19.5,15),animAll=animAll,
                   plotStrmRpts=plotStrmRpts,strmRptsF=strmRptsF,rprtWindow=rprtWindow)

