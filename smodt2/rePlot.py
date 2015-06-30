import tools
import sys
import os
import numpy as np

def rePlot():
    rootDir = '/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/'
    bestFit = np.array([ 1.18210269e+00,   1.02060609e-03,   4.85779987e+01,
         1.00627225e+02,   6.15042552e-02,   2.45150854e+06,
         2.45150854e+06,   1.19477854e+01,   4.74681216e+01,
         2.71046459e+02,   5.52729440e+00,   4.14134613e+01,
         8.77599249e+00,  -6.65358481e-02 ])  
    
    settingsDict = tools.startup(sys.argv,rootDir,rePlot=True)
    allFname = os.path.join(settingsDict['finalFolder'],'combined-BIstripped-MCMCdata.fits')
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=settingsDict['logLevel'])
    #effPtsStr = tools.mcmcEffPtsCalc(allFname)
    ##DIlims=[[[xMin,xMax],[yMin,yMax]],[[xCropMin,xCropMax],[yCropMin,yCropMax]]]   [[[,],[,]],[[,],[]]]
    ##RVlims=[[yMin,yMax],[yResidMin,yResidMax],[xMin,xMax]]
    plotFnameBase = os.path.join(settingsDict['finalFolder'],'orbitPlot-MANUAL-'+settingsDict['symMode'][0])
    tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps',DIlims=[],RVlims=[[-9.9,9.9],[-0.7,0.8],[-0.515,0.515]])
    plotFilename = os.path.join(settingsDict['finalFolder'],'summaryPlot-MANUAL-'+settingsDict['symMode'][0])
    clStr = tools.summaryPlotter(allFname, plotFilename,paramsToPlot=[0,1,7,8],xLims=[[0.0,3.0],[0.5,1.7],[11.5,12.5],[30,60]],stage=settingsDict['symMode'][0], shadeConfLevels=True,forceRecalc=False)
    
def stackedPosteriorsPlotterHackStarter():
    outputDataFilenames = []
    if False:
        outputDataFilenames.append('/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-tight-20PercentRealizedError--14-Million-in_Total/outputData-ALL.dat')
        outputDataFilenames.append('/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-tight-10PercentRealizedError--14-Million-in_Total/outputData-ALL.dat')
        outputDataFilenames.append('/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-tight-5PercentRealizedError--14-Million-in_Total/outputData-ALL.dat')
        outputDataFilenames.append('/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-tight-1PercentRealizedError--14-Million-in_Total/outputData-ALL.dat')
    else:
        outputDataFilenames.append('/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-veryOpen-20PercentRealizedError2--50-Million-in_Total/outputData-ALL.dat')
        outputDataFilenames.append('/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-veryOpen-10PercentRealizedError2--50-Million-in_Total/outputData-ALL.dat')
        outputDataFilenames.append('/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-veryOpen-5PercentRealizedError2--50-Million-in_Total/outputData-ALL.dat')
        outputDataFilenames.append('/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-veryOpen-1PercentRealizedError2--50-Million-in_Total/outputData-ALL.dat')
    
    plotFilename = os.path.join(os.path.abspath('/run/media/kmede/Data1/Todai_Work/Data/data_SMODT'),'stackedPosteriorsPlot')
    tools.stackedPosteriorsPlotter(outputDataFilenames, plotFilename,ALLparams=False)
    #print 'Final stacked plot file written to:\n'+plotFilename
    if True:
        print 'converted to PDF as well'
        os.system("epstopdf "+plotFilename+'.eps')
    
if __name__ == '__main__':
    rePlot()