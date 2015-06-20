import tools
import sys
import os
import numpy as np

def rePlot():
    rootDir = '/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/'
    bestFit = np.array([ 1.24369022e+00,   2.38859703e-01,   3.71435905e+01,
         2.46358039e+02,   3.72328628e-01,   2.45234737e+06,
         2.45234737e+06,   2.07526388e+01,   1.56254637e+02,
         3.50362217e+02,   8.61052528e+00,   1.67949208e+02,
         8.63887801e+02,   6.18916868e+03,   3.77725887e+02,
         6.32151151e+03])  
    
    settingsDict = tools.startup(sys.argv,rootDir,rePlot=True)
    allFname = os.path.join(settingsDict['finalFolder'],'combinedMCMCdata.fits')
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=settingsDict['logLevel'])
    plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbitPlot-MANUAL2-'+settingsDict['symMode'][0])
    tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps')
    #plotFilename = os.path.join(settingsDict['finalFolder'],'summaryPlot-MANUAL-'+settingsDict['symMode'][0])
    #clStr = tools.summaryPlotter(allFname, plotFilename,stage=settingsDict['symMode'][0], shadeConfLevels=True,forceRecalc = False)
    
if __name__ == '__main__':
    rePlot()