import tools
import sys
import os
import numpy as np

def rePlot():
    rootDir = '/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/'
    bestFit = np.array([9.75380123e-01,   9.06431634e-04,   5.16244701e+01,
         1.00505643e+02,   6.32600194e-02,   2.45148405e+06,
         2.45148405e+06,   1.19393851e+01,   4.70559841e+01,
         2.68872348e+02,   5.18192825e+00,   4.15189352e+01,
         8.80382274e+00,  -6.06784028e-02 ])  
    
    settingsDict = tools.startup(sys.argv,rootDir,rePlot=True)
    allFname = os.path.join(settingsDict['finalFolder'],'combined-BIstripped-MCMCdata.fits')
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=settingsDict['logLevel'])
    #effPtsStr = tools.mcmcEffPtsCalc(allFname)
    ##DIlims=[[[xMin,xMax],[yMin,yMax]],[[xCropMin,xCropMax],[yCropMin,yCropMax]]]   [[[,],[,]],[[,],[]]]
    ##RVlims=[[yMin,yMax],[yResidMin,yResidMax],[xMin,xMax]]
    plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbitPlot-MANUAL-'+settingsDict['symMode'][0])
    tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps',DIlims=[],RVlims=[[-9.9,9.9],[-0.7,0.8],[-0.515,0.515]])
    #plotFilename = os.path.join(settingsDict['finalFolder'],'summaryPlot-MANUAL-'+settingsDict['symMode'][0])
    #clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[0,1,7,8],xLims=[[0.0,3.0],[0.5,1.7],[11.5,12.5],[30,60]],stage=settingsDict['symMode'][0], shadeConfLevels=True,forceRecalc=False)
    
if __name__ == '__main__':
    rePlot()