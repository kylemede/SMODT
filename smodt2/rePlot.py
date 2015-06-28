import tools
import sys
import os
import numpy as np

def rePlot():
    rootDir = '/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/'
    bestFit = np.array([1.29317263e+00,   1.12053685e-03,   4.67045643e+01,
         9.91969983e+01,   5.63000318e-02,   2.45157942e+06,
         2.45157942e+06,   1.20863533e+01,   4.57675756e+01,
         2.78900707e+02,   5.73920929e+00,   3.73509461e+01,
         8.78771777e+00,   3.47375918e-02 ])  
    
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