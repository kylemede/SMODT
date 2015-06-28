import tools
import sys
import os
import numpy as np

def rePlot():
    rootDir = '/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/'
    bestFit = np.array([1.33134207e+00,   4.00500108e-01,   4.52060738e+01,
         3.00189064e+02,   8.02505359e-01,   2.45697536e+06,
         2.45697536e+06,   8.01343295e+01,   7.98487024e+01,
         2.09698594e+02,   2.23199300e+01,   3.83265956e+01,
         3.16527083e+03,   9.97265557e+00])  
    
    settingsDict = tools.startup(sys.argv,rootDir,rePlot=True)
    allFname = os.path.join(settingsDict['finalFolder'],'combined-BIstripped-MCMCdata.fits')
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=settingsDict['logLevel'])
    #effPtsStr = tools.mcmcEffPtsCalc(allFname)
    #plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbitPlot-MANUAL2-'+settingsDict['symMode'][0])
    #tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps')
    plotFilename = os.path.join(settingsDict['finalFolder'],'summaryPlot-MANUAL-'+settingsDict['symMode'][0])
    clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[0,1,3,4],stage=settingsDict['symMode'][0], shadeConfLevels=True,forceRecalc=False)
    
if __name__ == '__main__':
    rePlot()