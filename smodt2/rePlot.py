import tools
import sys
import os
import numpy as np

def rePlot():
    rootDir = '/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/'
    bestFit = np.array([  1.22621459e+00,   2.38591813e-01,   3.72900583e+01,
         2.46364239e+02,   3.73625558e-01,   2.45234549e+06,
         2.45234549e+06,   2.07372326e+01,   1.56484523e+02,
         3.50238668e+02,   8.57179143e+00,   1.68032374e+02,
         8.62630643e+02,   6.18726429e+03,   3.75940155e+02,
         6.32321901e+03])    #[0.055,0.054,55,313,0.17,2454623,2454623,10.23,95,245,0,0])
    
    settingsDict = tools.startup(sys.argv,rootDir,rePlot=True)
    allFname = os.path.join(settingsDict['finalFolder'],'combinedMCMCdata.fits')
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=settingsDict['logLevel'])
    plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbitPlot-MANUAL-'+settingsDict['symMode'][0])
    tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps')
    #plotFilename = os.path.join(os.path.dirname(allFname),'summaryPlot-MANUAL'+settingsDict['symMode'][0])
    #clStr = tools.summaryPlotter(allFname, plotFilename,stage=settingsDict['symMode'][0], shadeConfLevels=True,forceRecalc = False)
    
if __name__ == '__main__':
    rePlot()