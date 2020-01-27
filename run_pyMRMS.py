from pyMRMS import retrieveMRMS,readFL,annotate,graph


projID = 'VSE'
tdrRng = 47966.79

getMRMS = False
plotMRMS = True

plotAsts = True
plotRng = True
pltFT = True



#strtMRMS = '20170430-181600'
strtMRMS = '20170430-194000'
endMRMS = '20170430-232600'
#endMRMS = '20170430-185600'
mrmsSaveDir = '/Users/danstechman/Research/VORTEXSE/MRMS/20170430'
mrmsPlotDir = mrmsSaveDir + '/plots'
flFile = '/Users/danstechman/Research/VORTEXSE/FlightLevel/20170430H1_AC.nc'
astsPrmsF = '/Users/danstechman/Research/VORTEXSE/MRMS/20170430_assets.yml'
minLat = 33.6
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
                   plotAsts=plotAsts,prmsFile=astsPrmsF,plotRng=plotRng,
                   plotFltTrk=pltFT,flLon=flData['flLon'],flLat=flData['flLat'],flDT=flData['flDT'],
                   flHdng=flData['flHdng'],fltIntv=15,fadeFltTrk='fade',
                   projID=projID,fType='png',saveDir=mrmsPlotDir,figsize=(13,10),progDisp=True)

