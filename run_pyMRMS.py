from pyMRMS import retrieveMRMS,readFL,annotate,graph


projID = 'VSE'
tdrRng = 47966.79

getMRMS = False
plotMRMS = True

plotAsts = True
plotRng = False
pltMultRng = False
pltFT = True

animAll = True

fltIntv = 120 # Plotting interval for flight track (seconds)

plotDom = False
domOrig = (34.15, -87.3)
domX = 180
domY = 155



strtMRMS = '20170430-180000'
endMRMS = '20170430-220000'

mrmsSaveDir = '/Volumes/Pegasus/projData/vortexse17/MRMS/20170430/test'
mrmsPlotDir = '/Volumes/Pegasus/projData/vortexse17/MRMS/20170430/plots_SoundLocs'
flFile = '/Volumes/Pegasus/projData/vortexse17/flightLevel/20170430H1_AC.nc'
astsPrmsF = '/Volumes/Pegasus/projData/vortexse17/MRMS/20170430_assets_soundingLocs.yml'
# minLat = 33.6
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
                   fType='png',saveDir=mrmsPlotDir,figsize=(19.5,15),animAll=animAll)

