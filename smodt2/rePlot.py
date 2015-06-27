import tools
import sys
import os
import numpy as np

def rePlot():
    rootDir = '/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/'
    bestFit = np.array([1.15017844e-01,   0.00000000e+00,   5.40662305e+01,
         3.13364841e+02,   1.73993749e-01,   2.45463038e+06,
         2.45463038e+06,   1.01751015e+01,   9.59430071e+01,
         2.46721643e+02,   2.28345203e+00,   2.92536522e+01,
         0.00000000e+00])  
    
    settingsDict = tools.startup(sys.argv,rootDir,rePlot=True)
    allFname = os.path.join(settingsDict['finalFolder'],'combinedMCMCdata.fits')
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=settingsDict['logLevel'])
    effPtsStr = tools.mcmcEffPtsCalc(allFname)
    #plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbitPlot-MANUAL2-'+settingsDict['symMode'][0])
    #tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps')
    #plotFilename = os.path.join(settingsDict['finalFolder'],'summaryPlot-MANUAL-'+settingsDict['symMode'][0])
    #clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[],stage=settingsDict['symMode'][0], shadeConfLevels=True,forceRecalc = False)
    
if __name__ == '__main__':
    rePlot()