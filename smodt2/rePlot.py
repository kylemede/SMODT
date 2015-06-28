import tools
import sys
import os
import numpy as np

def rePlot():
    rootDir = '/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/'
    bestFit = np.array([1.12808647e+00,   1.05725558e-03,   4.81033871e+01,
         9.97488899e+01,   4.48211024e-02,   2.45694802e+06,
         2.45694802e+06,   1.19839147e+01,   4.47671348e+01,
         2.70731658e+02,   5.45288743e+00,   7.32365555e+01,
         8.94578734e+00,   7.87455080e-02 ])  
    
    settingsDict = tools.startup(sys.argv,rootDir,rePlot=True)
    allFname = os.path.join(settingsDict['finalFolder'],'combined-BIstripped-MCMCdata.fits')
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=settingsDict['logLevel'])
    #effPtsStr = tools.mcmcEffPtsCalc(allFname)
    #plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbitPlot-MANUAL-'+settingsDict['symMode'][0])
    #tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps')
    plotFilename = os.path.join(settingsDict['finalFolder'],'summaryPlot-MANUAL-'+settingsDict['symMode'][0])
    clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[],xLims=[[0.001,3.0],[0.2,0.5999],[0.78,0.82],[79.001,80.5]],stage=settingsDict['symMode'][0], shadeConfLevels=True,forceRecalc=False)
    
if __name__ == '__main__':
    rePlot()