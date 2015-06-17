import tools
import sys
import os
import numpy as np

def rePlot():
    rootDir = '/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/'
    bestFit = np.array([0.055,0.054,55,313,0.17,2454623,2454623,10.23,95,245,0,0])
    settingsDict = tools.startup(sys.argv,rootDir)
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=settingsDict['logLevel'])
    plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbitPlot-MANUAL-'+settingsDict['symMode'][0])
    tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps')
    
if __name__ == '__main__':
    rePlot()