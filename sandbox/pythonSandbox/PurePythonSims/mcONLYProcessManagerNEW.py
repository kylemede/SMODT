import os
import numpy as np
import shutil
import time
from multiprocessing import Process
from mcONLY3DorbSimUniform import mcOrbSimUniform
from mcONLY3DorbSimUniformTH_I import mcOrbSimUniformTH_I
from orbitToolbox import dataReaderNEW, burnInStripper, dataFileCombiner, dataReadAndPlotNEW3, bestOrbitFinderNEW, orbitEllipsePlotter3, rvPlotter

class MCProcessManager(Process):
    '''
    This is the Manager object that controls the processes for an MC simulator run 
    initiated using the mcSimStarter function.
    '''
    def __init__(self, tempTitle, numSamples, silent=True, verbose=False,\
                        dataTempDir='/home/Kyle/DataTempForSims/', dataFinalDir='/home/Kyle/DataTempForSims/', TH_I=False):
        
        Process.__init__(self)
        self.tempTitle = tempTitle
        self.numSamples = numSamples
        self.silent = silent
        self.verbose = verbose
        self.dataTempDir = dataTempDir
        self.dataFinalDir = dataFinalDir
        self.TH_I = TH_I
        
    def run(self):
        self.process(tempTitle=self.tempTitle, numSamples=self.numSamples, silent=self.silent, verbose=self.verbose, dataTempDir=self.dataTempDir, dataFinalDir=self.dataFinalDir, TH_I=self.TH_I)
    	
    
    def process(self, tempTitle, numSamples, silent=True, verbose=False, dataTempDir='/home/Kyle/DataTempForSims/', dataFinalDir='/home/Kyle/DataTempForSims/', TH_I=False):
        
        if not silent:
            print 'Starting to run process for file title: '+tempTitle
        
        ## Run simulator for this process
        if TH_I:
        	mcOrbSimUniformTH_I(tempTitle, numSamples, silent=silent, verbose=verbose, dataTempDir=dataTempDir, dataFinalDir=dataFinalDir)
        else:
         	mcOrbSimUniform(tempTitle, numSamples, silent=silent, verbose=verbose, dataTempDir=dataTempDir, dataFinalDir=dataFinalDir)

def mcSimStarter(filename, numSamples, numProcesses=1, silent=True, verbose=False,\
                        dataTempDir='/home/Kyle/DataTempForSims/', dataFinalDir='/run/media/Kyle/Data1/Todai_Work/workspace/Binary-project/data/', TH_I=False):
    """
    This is the master function to start a multiprocess MCMC simulator run.
    NOTE: The dataFinalDir used in the processes will be equal to the dataTempDir for faster file combining after processes
          finish and will the single final combined file will be written to the dataFinalDir in the end.
    
    """
    #dataPath =  #$$$$$$ DO SOMETHING TO MAKE THIS KINDA THING WORK, GRAB FROM PLOTFILETITLE OR...?? $$$$$$$$$$$$$
    #numProcesses = 2  #$$$$$ MAKE THIS AUTOMATIC BASED ON NUMBER OF PROCESSORS/THREADS AVAILABLE
    master = []
    
    # Ensure the filename has .txt at the end
    if filename[-4:]!='.txt':
        filename = filename+".txt"
    
    # Ensure directories provided have a final '/' character
    if dataTempDir[-1:] is not '/':
         dataTempDir = dataTempDir+'/'
    if dataFinalDir[-1:] is not '/':
         dataFinalDir = dataFinalDir+'/'
    
#    # record the time the chain started
#    startTimeMAIN = time.clock() # Sadly, these start and end times don't work right in multi-processing right now ?!?!
    
    if not silent:
        print '\nMultiprocess mcONLY: $$$$$$$$$$$$$ STARTED Multiprocess MC Sim $$$$$$$$$$$$$$$\n'
    
    ### start running threads ###
    for processNumber in range(numProcesses):
        # create temp title for this thread's temp data file
        tempProcessTitle = 'tempProcessTitle-process_'+str(processNumber+1)+'.txt'
        # start process
        master.append(MCProcessManager(tempTitle=tempProcessTitle, numSamples=numSamples, \
                      silent=silent, verbose=verbose, dataTempDir=dataTempDir, dataFinalDir=dataTempDir, TH_I=TH_I))
        master[processNumber].start()
        
    # wait for completion of all process $$$$ STILL DON'T UNDERSTAND HOW THIS WORKS $$$$
    for processNumber in range(numProcesses):
        master[processNumber].join()
        
    # collect all temp filenames together and write to 'burn in included' final files in finalFolder
    tempFiles = []
    finalFolder = dataFinalDir+filename[:-4]+'/'
    dataFinalFilename = finalFolder+'outputData-ALL.txt' #filename
    
    print 'Starting to write original data and resultsDicts to final files'
    for processNumber in range(numProcesses):
        tmpFileWithDir = dataTempDir+master[processNumber].tempTitle
        tempFiles.append(tmpFileWithDir)
        
        # retrieve results dictionary file and move it to final location and name 
        resultsDictTempFilename = tmpFileWithDir[:-4]+'-resultsDict.txt'
        resultsDictFinalFilename = finalFolder+'results-chain'+str(processNumber+1)+'.txt'
        shutil.move(resultsDictTempFilename,resultsDictFinalFilename)   
        
    
    print '\nNow combining all files into one final file\n'
    # combine the input files that burned in successfully into one final file
    dataFileCombiner(tempFiles, dataFinalFilename)
        
    # kill off INS temp files
    for tempFilename in tempFiles:
        if os.path.exists(tempFilename):
            os.remove(tempFilename)
        else:
            print 'Problem: '+tempFilename+' could not be found to remove it!!!'    
    
#    # record the time the chain finished and print
#    endTimeMAIN = time.clock()
#    totalTimeMAIN = (endTimeMAIN-startTimeMAIN) # in seconds
#    if totalTimeMAIN<60:
#        print 'Multiprocess MCMC: The took '+str(totalTimeMAIN)+' seconds to complete.\n'  ##### only print in 'silent' mode to track time
#    elif totalTimeMAIN>3600:
#        totalTimeMAIN = totalTimeMAIN/3600
#        print 'Multiprocess MCMC: The took '+str(totalTimeMAIN)+' hours and '+str(60*((totalTimeMAIN/3600.0)-(totalTimeMAIN/3600)))+' minutes to complete.\n'  ##### only print in 'silent' mode to track time
#    else:
#        totalTimeMAIN = totalTimeMAIN/60
#        print 'Multiprocess MCMC: The took '+str(totalTimeMAIN)+' minutes and '+str(60*((totalTimeMAIN/60.0)-(totalTimeMAIN/60)))+' seconds to complete.\n'  ##### only print in 'silent' mode to track time
    
    if not silent:
        print '\nMultiprocess mcONLY: $$$$$$$$$$$$$ FINISHED Multiprocess MC Sim $$$$$$$$$$$$$$$\n'
        
    print '**** Now starting to make a weighted and non-weighted summary plots of data in final file  ****'
    
    dataReadAndPlotNEW3(dataFinalFilename,finalFolder+'summaryPlot', weight=False, confLevels=False)
    dataReadAndPlotNEW3(dataFinalFilename,finalFolder+'summaryPlot-weighted', weight=True, confLevels=True)
    
    bestOrbit = bestOrbitFinderNEW(dataFinalFilename, printToScreen=False, saveToFile=False, returnAsList=True)
    orbitEllipsePlotFilename = finalFolder+"orbitEllipsePlot"
    orbitEllipsePlotter3(bestOrbit[0],bestOrbit[1],bestOrbit[3],bestOrbit[4],bestOrbit[5],bestOrbit[6],\
                         plotFilename=orbitEllipsePlotFilename,show=True)
     
    rvPlotter(bestOrbit[0], bestOrbit[1],bestOrbit[2],bestOrbit[3],bestOrbit[4],bestOrbit[5],bestOrbit[6], \
              RV_origin_vel_0_proposed=bestOrbit[7], RV_origin_vel_1_proposed=bestOrbit[8],\
              plotFilename=finalFolder+"orbitRVplot", show=True)
    print '\n**** EVERYTHING FINISHED ****\n'    
    
        