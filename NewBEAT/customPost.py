import tools
from tools import constants as const
import sys
import os
import numpy as np
import glob

def customPost():
    rootDir = '/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/NewBEAT/'
    settingsDict = tools.startup(sys.argv,rootDir,rePlot=True)
    allFname = os.path.join(settingsDict['finalFolder'],"combined-BIstripped-MCMCdata.fits")
    skipBurnInStrip=True
    if os.path.exists(allFname)==False:
        allFname = os.path.join(settingsDict['finalFolder'],'combinedMCMCdata.fits')
        skipBurnInStrip=False
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=10)
    ##make list of stages to run
    stageList = []
    if settingsDict['symMode'][0]=='MC':
        stageList = ['MC']
    elif settingsDict['symMode'][0]=='SA':
        stageList = ['SA']
    elif settingsDict['symMode'][0]=='MCMC':
        stageList = ['SA','ST','MCMC']

    ##make hack list of output files
    outFiles = np.sort(glob.glob(os.path.join(settingsDict['finalFolder'],"outputDataMCMC*.fits")))
    
    ## calc and strip burn-in?
    if True:
        
        burnInStr = ''
        if skipBurnInStrip==False:
            if (len(outFiles)>1)and(settingsDict['CalcBurn'] and(settingsDict['symMode'][0]=='MCMC')):
                (burnInStr,burnInLengths) = tools.burnInCalc(outFiles,allFname)    
                if settingsDict['rmBurn'][0]:
                    strippedFnames = tools.burnInStripper(outFiles,burnInLengths)
                    outFiles = strippedFnames
                    ## combine stripped files to make final file?
                    if len(strippedFnames)>0:
                        strippedAllFname = os.path.join(os.path.dirname(strippedFnames[0]),"combined-BIstripped-MCMCdata.fits")
                        tools.combineFits(strippedFnames,strippedAllFname)
                        ## replace final combined filename with new stripped version
                        allFname = strippedAllFname
        
    ## find best fit
    if True:
        if os.path.exists(allFname):
            bestFit = tools.findBestOrbit(allFname)
    else:
        bestFit = np.array([  1.32366515e+00,   2.73842875e-01,   3.59629752e+01,
                             2.45465119e+02,   3.72512416e-01,   2.45234944e+06,
                             2.45234944e+06,   2.09427724e+01,   1.58152772e+02,
                             3.50344680e+02,   8.88140043e+00,   8.19243621e+01,
                             8.68267992e+02,   6.19627726e+03,   3.85209089e+02,
                             6.31690499e+03])
        
    #effPtsStr = tools.mcmcEffPtsCalc(allFname)
    
    if True:
        ##for reference: DIlims=[[[xMin,xMax],[yMin,yMax]],[[xCropMin,xCropMax],[yCropMin,yCropMax]]]   [[[,],[,]],[[,],[]]]
        ##               RVlims=[[yMin,yMax],[yResidMin,yResidMax],[xMin,xMax]]
        plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbitPlot-MANUAL-'+settingsDict['symMode'][0])
        tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps',DIlims=[],RVlims=[])
        #tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps',DIlims=[],RVlims=[[-9.9,9.9],[-0.7,0.8],[-0.515,0.515]])
        
    clStr=''
    if True:
        plotFilename = os.path.join(settingsDict['finalFolder'],'summaryPlot-MANUAL-'+settingsDict['symMode'][0])
        clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[],xLims=[],stage=settingsDict['symMode'][0], shadeConfLevels=True,forceRecalc=False)
        #clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[0,1,4,8],xLims=[[0.5,2.1],[0.5,1.7],[-0.005,0.125],[37,53]],bestVals=[1.0,1.0,0.048,45.0],stage=settingsDict['symMode'][0], shadeConfLevels=True,forceRecalc=False)
    if False: 
        plotFilename = os.path.join(settingsDict['finalFolder'],'summaryPlot-MANUAL-'+settingsDict['symMode'][0])
        tools.cornerPlotter(allFname, plotFilename)#,paramsToPlot=[0,1,7,8],xLims=[[0.5,2.1],[0.5,1.7],[11.5,13.1],[37,53]],bestVals=[1.0,1.0,11.9,45.0])
        
    ##calc R?
    grStr = ''
    if (len(outFiles)>1) and (settingsDict['CalcGR'] and (settingsDict['symMode'][0]=='MCMC')):
        (GRs,Ts,grStr) = tools.gelmanRubinCalc(outFiles,settingsDict['nSamples'][0])
        
    ## custom re check of the orbit fit     
    if False: 
        orbParams = [1.16679785333,0.263640913869,37.2109325013,245.327413195,0.374265026422,2452348.33995,2452348.33995,20.9645791913,158.909378023,350.3113539,8.57166434448,92.2373136221,868.357550689,6196.81653647,385.187801642,6318.18670944]
        finalFits=''
        nus = [70, 2, 66]
        epochs=[2457259.0]
        tools.predictLocation(orbParams,settingsDict,epochs)
        #tools.recheckFit3D(orbParams,settingsDict,finalFits,nus)
    
    ## following post-processing stages can take a long time, so write the current
    ## summary information to the summary file and add the rest later
    if True:
        if os.path.exists(allFname):
            tools.summaryFilePart1(settingsDict,stageList,allFname,clStr,burnInStr,bestFit,grStr)
            
    ## calc correlation length & number effective points? # This one takes a long time for long runs!!!
    effPtsStr = ''
    if ((len(outFiles)>1)and(settingsDict['symMode'][0]=='MCMC'))and (settingsDict['calcCL'] and os.path.exists(allFname)):
        effPtsStr = tools.mcmcEffPtsCalc(allFname)
        
    if os.path.exists(allFname):
        tools.summaryFilePart2(settingsDict,effPtsStr,22,22)
        
    ##clean up files (move to folders or delete them)
    tools.cleanUp(settingsDict,stageList,allFname)
    if settingsDict['CopyToDB']:
        tools.copyToDB(settingsDict)
    
def stackedPosteriorsPlotterHackStarter():
    outputDataFilenames = []
    
    #outputDataFilenames.append('/run/media/kmede/Data1/Todai_Work/Data/data_SMODT/SMODT2-SyntheticJUPITER-3D-20percent-startAtBest-lowEccTrue/combined-BIstripped-MCMCdata.fits')
    outputDataFilenames.append('/run/media/kmede/Data1/Todai_Work/Data/data_SMODT/JUPITER2-3D-10percent-startAtBest-lowEccTrue-newPriors-gaussParallax/combined-BIstripped-MCMCdata.fits')
    outputDataFilenames.append('/run/media/kmede/Data1/Todai_Work/Data/data_SMODT/JUPITER2-3D-5percent-startAtBest-lowEccTrue-newPriors-gaussParallax/combined-BIstripped-MCMCdata.fits')
    outputDataFilenames.append('/run/media/kmede/Data1/Todai_Work/Data/data_SMODT/JUPITER2-3D-1percent-startAtBest-lowEccTrue-newPriors-gaussParallax/combined-BIstripped-MCMCdata.fits')
    
    plotFilename = os.path.join(os.path.abspath('/run/media/kmede/Data1/Todai_Work/Data/data_SMODT'),'stackedPosteriors3-lowEccTrue')
    tools.stackedPosteriorsPlotter(outputDataFilenames, plotFilename,paramsToPlot=[1,4,8],xLims=[[0.55,1.8],[-0.005,0.125],[28,61]])
    #tools.stackedPosteriorsPlotter(outputDataFilenames, plotFilename,paramsToPlot=[],xLims=[])
    #print 'Final stacked plot file written to:\n'+plotFilename
    if True:
        print 'converted to PDF as well'
        os.system("epstopdf "+plotFilename+'.eps')
        
        
def paramConverterTest():
    
    ## Make Orbit cpp obj
    Orbit = tools.cppTools.Orbit()
    Orbit.loadStaticVars(0,0,True)
    Orbit.loadConstants(const.Grav,const.pi,const.KGperMsun, const.daysPerYear,const.secPerYear,const.MperAU)
    bestFit = np.array([ 1.08040127e+00,   9.60493873e-04,   5.02111885e+01,
          1.00808667e+02,   6.29123973e-02,   2.45151739e+06,
          2.45151739e+06,   1.19376399e+01,   4.75040520e+01,
          2.71592539e+02,   5.36101436e+00,   4.15438649e+01,
          8.77774317e+00,  -5.73899482e-02]) 
    ## ensure orbParams are in required format for Orbit
    params = []
    for par in bestFit:
        params.append(par)
    params=np.array(params,dtype=np.dtype('d'),order='C')
    
    es = [0.1,0.25,0.5,0.75,0.99]
    omegas=[1.,45.,90.,135.,180.,225.,270.,315.,360.]
    
    for i in range(len(omegas)):
        params[4]=0.25
        params[9]=omegas[i]
        print 'before: e= '+str(params[4])+', omega= '+str(params[9])
        Orbit.convertParsToRaw(params)
        print 'raw: par[4]= '+str(params[4])+', par[9]= '+str(params[9])
        Orbit.convertParsFromRaw(params)
        print 'after: e= '+str(params[4])+', omega= '+str(params[9])
    
if __name__ == '__main__':
    customPost()
    #stackedPosteriorsPlotterHackStarter()
    #paramConverterTest()
    