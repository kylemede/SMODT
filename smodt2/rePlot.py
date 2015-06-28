import tools
import sys
import os
import numpy as np

def rePlot():
    rootDir = '/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/'
    bestFit = np.array([1.00000000e+00,   3.33000000e-01,   5.00000000e+01,
         5.03505022e+01,   1.05420763e-03,   2.45701457e+06,
         2.45701457e+06,   2.00324241e-01,   1.00390101e+01,
         1.62404414e+02,   3.76770340e-01,   6.30389795e+01,
         2.43956258e+03,   8.59658550e+00])  
    
    settingsDict = tools.startup(sys.argv,rootDir,rePlot=True)
    allFname = os.path.join(settingsDict['finalFolder'],'combined-BIstripped-MCMCdata.fits')
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=settingsDict['logLevel'])
    #effPtsStr = tools.mcmcEffPtsCalc(allFname)
    plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbitPlot-MANUAL-'+settingsDict['symMode'][0])
    tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps')
    #plotFilename = os.path.join(settingsDict['finalFolder'],'summaryPlot-MANUAL-'+settingsDict['symMode'][0])
    #clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[0,1,3,4],stage=settingsDict['symMode'][0], shadeConfLevels=True,forceRecalc=False)
    
if __name__ == '__main__':
    rePlot()