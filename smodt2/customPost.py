import tools
from tools import constants as const
import sys
import os
import numpy as np
import glob

def customPost():
    rootDir = '/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/'
    settingsDict = tools.startup(sys.argv,rootDir,rePlot=True)
    #allFname = os.path.join(settingsDict['finalFolder'],'combinedMCMCdata.fits')
    allFname = os.path.join(settingsDict['finalFolder'],"combined-BIstripped-MCMCdata.fits")
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=10)
    ##make list of stages to run
    stageList = []
    if settingsDict['symMode'][0]=='MC':
        stageList = ['MC']
    elif settingsDict['symMode'][0]=='SA':
        stageList = ['SA']
    elif settingsDict['symMode'][0]=='MCMC':
        stageList = ['SA','ST','MCMC']
#     bestFit = np.array([ 1.08040127e+00,   9.60493873e-04,   5.02111885e+01,
#          1.00808667e+02,   6.29123973e-02,   2.45151739e+06,
#          2.45151739e+06,   1.19376399e+01,   4.75040520e+01,
#          2.71592539e+02,   5.36101436e+00,   4.15438649e+01,
#          8.77774317e+00,  -5.73899482e-02])  
    ##make hack list of output files
    outFiles = np.sort(glob.glob(os.path.join(settingsDict['finalFolder'],"outputDataMCMC*.fits")))
    
    ## calc and strip burn-in?
    burnInStr = ''
    if (len(outFiles)>1)and(settingsDict['CalcBurn'] and(settingsDict['symMode'][0]=='MCMC')):
        (burnInStr,burnInLengths) = tools.burnInCalc(outFiles,allFname)    
        if settingsDict['delBurn'][0]:
            strippedFnames = tools.burnInStripper(outFiles,burnInLengths)
            outFiles = strippedFnames
            ## combine stripped files to make final file?
            if len(strippedFnames)>0:
                strippedAllFname = os.path.join(os.path.dirname(strippedFnames[0]),"combined-BIstripped-MCMCdata.fits")
                tools.combineFits(strippedFnames,strippedAllFname)
                ## replace final combined filename with new stripped version
                allFname = strippedAllFname
    
    ## find best fit
    if os.path.exists(allFname):
        bestFit = tools.findBestOrbit(allFname)
        
    #effPtsStr = tools.mcmcEffPtsCalc(allFname)
    
    ##for reference: DIlims=[[[xMin,xMax],[yMin,yMax]],[[xCropMin,xCropMax],[yCropMin,yCropMax]]]   [[[,],[,]],[[,],[]]]
    ##               RVlims=[[yMin,yMax],[yResidMin,yResidMax],[xMin,xMax]]
    plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbitPlot-MANUAL-'+settingsDict['symMode'][0])
    #tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps',DIlims=[],RVlims=[])
    tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps',DIlims=[],RVlims=[[-9.9,9.9],[-0.7,0.8],[-0.515,0.515]])
    
    clStr=''
    plotFilename = os.path.join(settingsDict['finalFolder'],'summaryPlot-MANUAL-'+settingsDict['symMode'][0])
    #clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[],xLims=[],stage=settingsDict['symMode'][0], shadeConfLevels=True,forceRecalc=False)
    clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[0,1,7,8],xLims=[[0.0,3.0],[0.5,1.7],[11.5,12.5],[30,60]],stage=settingsDict['symMode'][0], shadeConfLevels=True,forceRecalc=False)
    
    ##calc R?
    grStr = ''
    if (len(outFiles)>1) and (settingsDict['CalcGR'] and (settingsDict['symMode'][0]=='MCMC')):
        (GRs,Ts,grStr) = tools.gelmanRubinCalc(outFiles,settingsDict['nSamples'][0])
    
    ## following post-processing stages can take a long time, so write the current
    ## summary information to the summary file and add the rest later
    if os.path.exists(allFname):
        tools.summaryFilePart1(settingsDict,stageList,allFname,clStr,burnInStr,bestFit,grStr)
    
def stackedPosteriorsPlotterHackStarter():
    outputDataFilenames = []
    
    #outputDataFilenames.append('/run/media/kmede/Data1/Todai_Work/Data/data_SMODT/SMODT2-SyntheticJUPITER-3D-20percent-startAtBest-lowEccTrue/combined-BIstripped-MCMCdata.fits')
    outputDataFilenames.append('/run/media/kmede/Data1/Todai_Work/Data/data_SMODT/SMODT2-JUPITER2-3D-10percent-startAtBest-lowEccTrue-3/combined-BIstripped-MCMCdata.fits')
    outputDataFilenames.append('/run/media/kmede/Data1/Todai_Work/Data/data_SMODT/SMODT2-JUPITER2-3D-5percent-startAtBest-lowEccTrue-3/combined-BIstripped-MCMCdata.fits')
    outputDataFilenames.append('/run/media/kmede/Data1/Todai_Work/Data/data_SMODT/SMODT2-JUPITER2-3D-1percent-tight-startAtBest-lowEccTrue-3/combined-BIstripped-MCMCdata.fits')
    
    plotFilename = os.path.join(os.path.abspath('/run/media/kmede/Data1/Todai_Work/Data/data_SMODT'),'stackedPosteriors3-lowEccTrue')
    tools.stackedPosteriorsPlotter(outputDataFilenames, plotFilename,paramsToPlot=[2,4,8],xLims=[[35,65],[-0.005,0.125],[28,61]])
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
    #customPost()
    stackedPosteriorsPlotterHackStarter()
    #paramConverterTest()
    