import tools
import sys
import os
import numpy as np

def rePlot():
    rootDir = '/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/'
    bestFit = np.array([ 1.00000000e+00,   3.28506077e-01,   3.77924938e+01,
         2.43863828e+02,   3.80925182e-01,   2.45232874e+06,
         2.45232874e+06,   2.10947212e+01,   1.64122595e+02,
         3.49158912e+02,   8.39230695e+00,   2.07005517e+02,
         8.67012992e+02,   6.19531253e+03,   3.87337955e+02,
         6.32026813e+03])  
    
    settingsDict = tools.startup(sys.argv,rootDir,rePlot=True)
    allFname = os.path.join(settingsDict['finalFolder'],'combined-BIstripped-MCMCdata.fits')
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=settingsDict['logLevel'])
    #effPtsStr = tools.mcmcEffPtsCalc(allFname)
    plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbitPlot-MANUAL-'+settingsDict['symMode'][0])
    tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps')
    plotFilename = os.path.join(settingsDict['finalFolder'],'summaryPlot-MANUAL-'+settingsDict['symMode'][0])
    clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[0,1,2,3,4,5,6,7],stage=settingsDict['symMode'][0], shadeConfLevels=True,forceRecalc=False)
    
if __name__ == '__main__':
    rePlot()