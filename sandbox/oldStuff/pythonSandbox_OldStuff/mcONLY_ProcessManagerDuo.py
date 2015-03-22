#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import os
import time
from multiprocessing import Process
from Toolboxes.orbitToolboxDuo import *
from Toolboxes.plottingToolbox import *
from Toolboxes.DItoolbox import *
from Toolboxes.RVtoolbox import *

class mcProcessManager(Process):
    '''
    This is the Manager object that controls the processes for an MCMC simulator run 
    initiated using the mcmcSimStarter function.
    '''
    def __init__(self, simSettingsFilename, filename, cplusplusCodeDir):
        
        Process.__init__(self)
        self.simSettingsFilename = simSettingsFilename
        self.filename = filename
        self.cplusplusCodeDir = cplusplusCodeDir        
        
    def run(self):
        self.process(simSettingsFilename=self.simSettingsFilename, 
                     filename=self.filename, cplusplusCodeDir=self.cplusplusCodeDir)
        
    
    def process(self, simSettingsFilename, filename, cplusplusCodeDir):
    
        print 'Starting to run process for file title: '+filename
        
        CplusplusCodeCALL = os.path.join(cplusplusCodeDir,'mcONLYorbSimulator')
        CplusplusCodeCALL = CplusplusCodeCALL+' '+simSettingsFilename+' '+filename
        print 'CplusplusCodeCALL : ',CplusplusCodeCALL
        print "call being made"
        os.system(CplusplusCodeCALL)
        print 'call completed'
        
def mcSimStarter(paramSettingsDict):
    """
    This is the master function to start a multiprocess Monte Carlo simulator run.
    
    """
    # TEMP!!!!!!!!!!!!!! REMOVE FOLLOWING COMMENT BLOCK ASAP!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    """
    # set simulation type
    mcONLY = paramSettingsDict['mcONLY'] # using only standard 'shotgun' monte carlo?
    
    # set simulation general settings
    numProcesses = paramSettingsDict['numProcesses']    
    silent = paramSettingsDict['silent']
    verbose = paramSettingsDict['verbose']
    data_dir = paramSettingsDict['outputData_filenameRoot']
    filenameRoot = paramSettingsDict['outputData_filenameRoot']  
    stripBurnIn = paramSettingsDict['stripBurnIn']
    calcCorrLengths = paramSettingsDict['calcCorrLengths']
    DIonly = paramSettingsDict['DIonly']
    RVonly = paramSettingsDict['RVonly']
    """
    # record the time the chain started
    startTime = time.clock()

    master = []
    
    # make call to 'make' to build the C++ code for this simulation
    # to ensure it is good to go
    pwd = os.curdir
    makeDir = paramSettingsDict['cplusplusCodeDir'] 
    print "\n"+'*'*95+"\n"+'*'*33+"  About to compile C++ code  " +'*'*33+"\n"+'*'*95+"\n" 
    print 'makeDir:'+makeDir
    os.chdir(makeDir)
    os.system('make mcONLYorbSimulator')
    os.chdir(pwd)
    print "\n"+'*'*95+"\n"+'*'*38+' C++ code compiled '+'*'*38+"\n"+'*'*95+'\n'
        
    print '\nMultiprocess: $$$$$$$$$$$$$ STARTED Multiprocess mcONLY $$$$$$$$$$$$$$$\n'
    
    simSettingsFilename = paramSettingsDict['UpdatedSettingsFile']
    
    ### start running threads ###
    numProcesses = paramSettingsDict['numProcesses']
    for processNumber in range(numProcesses):
        # create title for this thread's output data file
        [filenameRootNOext,ext] = os.path.splitext(paramSettingsDict['outputData_filenameRoot']) 
        
        if ext=='':
            ext='.dat' # just a default value, '.txt' is also fine
            
        filenameUSE = 'outputData-chain_'+str(processNumber+1)+ext
        filenameUSE_full = os.path.join(paramSettingsDict['outputData_dir'],filenameUSE)
        # start process
        master.append(mcProcessManager(simSettingsFilename=simSettingsFilename, 
                                         filename=filenameUSE_full, 
                                         cplusplusCodeDir=paramSettingsDict['cplusplusCodeDir']))
        master[processNumber].start()
        
    # wait for completion of all process $$$$ STILL DON'T UNDERSTAND HOW THIS WORKS $$$$
    for processNumber in range(numProcesses):
        master[processNumber].join()    
        
    # collect all temp filenames together and write to 'burn in included' final files in finalFolder
    dataFiles = []
    dataFinalFilename = os.path.join(paramSettingsDict['outputData_dir'],'outputData-ALL.dat') 
    
    print 'Starting to write original data to final combined file'
    for processNumber in range(numProcesses):
        dataFiles.append(master[processNumber].filename)
        
    print '\nNow combining all files into one final file\n'
    # combine the input files into one final file
    dataFileCombiner(dataFiles, dataFinalFilename)
    
    print '\nMultiprocess mcONLY: $$$$$$$$$$$$$ FINISHED Multiprocess MC Sim $$$$$$$$$$$$$$$\n' 
    ## make general parameter result summary figures
    summaryPlotFile = os.path.join(paramSettingsDict['outputData_dir'],'summaryPlot')
    weightedSummaryPlotFile = os.path.join(paramSettingsDict['outputData_dir'],'summaryPlot-weighted')
    
    sysDatafilename = os.path.join(paramSettingsDict['settings_and_InputDataDir'],'SystemData.txt')
    sysDataDict = sysDataToDict(sysDatafilename)
    
    bestOrbit = bestOrbitFinderNEW(dataFinalFilename, printToScreen=False, saveToFile=False, returnAsList=True)
    
    ## Make DI ellipse plot if DI data exists
    DIdatafilename = os.path.join(paramSettingsDict['settings_and_InputDataDir'],'DIdata.dat')
    numEpochs_DI = 0
    if os.path.exists(DIdatafilename)and (paramSettingsDict['RVonly']==False):
        DIdataDict = DIdataToDict(DIdatafilename)
        numEpochs_DI = DIdataDict['numEpochs']
        
        orbitEllipsePlotFilename = os.path.join(paramSettingsDict['outputData_dir'],'orbitEllipsePlot')
        orbitEllipsePlotterDuo(bestOrbit[0],bestOrbit[1],bestOrbit[3],bestOrbit[4],bestOrbit[5],bestOrbit[6],\
                             sysDataDict,DIdataDict,plotFilename=orbitEllipsePlotFilename,show=False)
     
    ## Make RV scatter.trend plots if RV data exists
    #print '\nWARNING: RV plotting still disabled in mcONLY_ProcessManagerDuo!'
    RVdatafilename = os.path.join(paramSettingsDict['settings_and_InputDataDir'],'RVdata.dat')
    numEpochs_RV = 0
    if os.path.exists(RVdatafilename) and (paramSettingsDict['DIonly']==False):
        RVdataDict = RVdataToDict(RVdatafilename)
        numEpochs_RV = RVdataDict['numEpochs']
            
        rvPlotFilename = os.path.join(paramSettingsDict['outputData_dir'],'orbitRVplot')

        # full orbit
        rvPlotterDuo(bestOrbit[1],bestOrbit[2],bestOrbit[3],bestOrbit[4],bestOrbit[5],bestOrbit[6], \
                  sysDataDict,RVdataDict,paramSettingsDict,K=bestOrbit[7],RVoffsets=bestOrbit[8],\
                  plotFilename=rvPlotFilename+'-FullOrbit', show=False, plotFullOrbit=True)
        if True:
                # half orbit
                rvPlotterDuo(bestOrbit[1],bestOrbit[2],bestOrbit[3],bestOrbit[4],bestOrbit[5],bestOrbit[6], \
                      sysDataDict,RVdataDict,paramSettingsDict,K=bestOrbit[7],RVoffsets=bestOrbit[8],\
                      plotFilename=rvPlotFilename, show=False, plotFullOrbit=False)
    
    
    print '**** Now starting to make a weighted and non-weighted summary plots of data in final file  ****'
    
    ## calculate nu for converting reducedChiSquareds to non-reduced for weighting
    if paramSettingsDict['RVonly']==False:
        numParams = 7
    else:
        numParams = 6
    nu = (1.0*numEpochs_RV)+(2.0*numEpochs_DI)-numParams
    
    print "\n"+"*"*50
    print "For weighted plots, using nu = "+str(nu)
    print "OR 1/nu = "+str(1.0/nu)
    print "\n"+"*"*50+"\n"
    
    plot4x1 = False
    if ((paramSettingsDict["longAN_degMAX"]==0)and(paramSettingsDict["periodMAX"]==0)and(paramSettingsDict["inclination_degMAX"]==0)):
        plot4x1 = True
        
    #dataReadAndPlotNEW3(dataFinalFilename,weightedSummaryPlotFile, weight=True, confLevels=True, nu=nu , plot4x1=plot4x1)
    summaryPlotter2(dataFinalFilename, weightedSummaryPlotFile+"-NewConfLevelStyle", weight=True, confLevels=True, nu=nu, plot4x1=plot4x1)
    if True:
        #dataReadAndPlotNEW3(dataFinalFilename,summaryPlotFile, weight=False, confLevels=True, plot4x1=plot4x1)
        summaryPlotter2(dataFinalFilename, summaryPlotFile+"-NewConfLevelStyle", weight=False, confLevels=True, nu=1, plot4x1=plot4x1)
    

    print '\n**** EVERYTHING FINISHED ****\n'
    
    # record the time the chain finished and print
    endTime = time.clock()
    totalTime = (endTime-startTime) # in seconds
    totalTimeString = timeString(totalTime)
    print 'MCONLY: Total simulation took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
    
    
#def main():
#    
#    cplusplusCodeDir = "/run/media/Kyle/Data1/Todai_Work/Data/data_Binary/data_C++/"
#    CplusplusFunc = 'mathyTEST'
#    # compile code just to ensure it is good to go
#    os.system('g++ -Wall -pedantic -lrt -o '+cplusplusCodeDir+CplusplusFunc+' '+cplusplusCodeDir+CplusplusFunc+'.cpp')
#    
#    filename = 'sillyWillyTest.txt'
#    numSamples = int(1e6)
#    numProcesses = int(4)
#    
#    mcmcSimStarter(filename, numSamples, numProcesses=numProcesses)
#                
#        
##############################################################
## end
##############################################################
#
#if __name__ == '__main__':
#    main()       
        
        