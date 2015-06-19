import tools
import sys
import os
import numpy as np

def rePlot():
    rootDir = '/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/'
    bestFit = np.array([  1.18372818e+00,   2.42381915e-01,   3.75982532e+01,
         2.46284056e+02,   3.73239929e-01,   2.45234635e+06,
         2.45234635e+06,   2.07550963e+01,   1.57255024e+02,
         3.50285017e+02,   8.50051370e+00,   1.67851673e+02,
         8.64072059e+02,   6.18912203e+03,   3.77111278e+02,
         6.32174239e+03])    #[0.055,0.054,55,313,0.17,2454623,2454623,10.23,95,245,0,0])
    
    settingsDict = tools.startup(sys.argv,rootDir,rePlot=True)
    allFname = os.path.join(settingsDict['finalFolder'],'combinedMCMCdata.fits')
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=settingsDict['logLevel'])
    plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbitPlot-MANUAL2-'+settingsDict['symMode'][0])
    tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps')
    plotFilename = os.path.join(settingsDict['finalFolder'],'summaryPlot-MANUAL-'+settingsDict['symMode'][0])
    clStr = tools.summaryPlotter(allFname, plotFilename,stage=settingsDict['symMode'][0], shadeConfLevels=True,forceRecalc = False)
    
if __name__ == '__main__':
    rePlot()