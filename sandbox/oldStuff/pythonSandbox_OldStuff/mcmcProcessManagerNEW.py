import os
import numpy as np
import time
from multiprocessing import Process
from mcmcOrbSimUniform7_3 import mcmcOrbSimUniform
from orbitToolbox import bestOrbitFinderNEW
from orbitToolbox2 import dataFileCombiner

class MCMCProcessManager(Process):
    '''
    This is the Manager object that controls the processes for an MCMC simulator run 
    initiated using the mcmcSimStarter function.
    '''
    def __init__(self, tempTitle, numSamples, temperature=1.0, sigmaPercent=1.0, startParams=False, silent=True, verbose=False,\
                        dataTempDir='/home/Kyle/DataTempForSims/', dataFinalDir='/home/Kyle/DataTempForSims/'):
        
        Process.__init__(self)
        self.tempTitle = tempTitle
        self.numSamples = numSamples
        self.temperature = temperature
        self.sigmaPercent = sigmaPercent
        self.startParams = startParams
        self.silent = silent
        self.verbose = verbose
        self.dataTempDir = dataTempDir
        self.dataFinalDir = dataFinalDir
        
    def run(self):
        self.process(tempTitle=self.tempTitle, numSamples=self.numSamples, temperature=self.temperature, sigmaPercent=self.sigmaPercent, \
                     startParams=self.startParams, silent=self.silent, verbose=self.verbose, dataTempDir=self.dataTempDir, dataFinalDir=self.dataFinalDir)
    
    def process(self, tempTitle, numSamples, temperature=1.0, sigmaPercent=1.0, startParams=False, silent=True, verbose=False,\
                        dataTempDir='/home/Kyle/DataTempForSims/', dataFinalDir='/home/Kyle/DataTempForSims/'):
        
        if not silent:
            print 'Starting to run process for file title: '+tempTitle
        
        ## Run simulator for this process
        mcmcOrbSimUniform(tempTitle, numSamples, temperature=temperature, sigmaPercent=sigmaPercent, startParams=startParams, silent=silent,\
                          verbose=verbose, dataTempDir=dataTempDir, dataFinalDir=dataFinalDir)

def mcmcSimStarter(filename, numSamples, numProcesses=1, temperature=1.0, sigmaPercent=1.0, startParams=False, silent=True, verbose=False,\
                        dataTempDir='/home/Kyle/DataTempForSims/', dataFinalDir='/media/Data1/Todai_Work/workspace/Binary-project/data/'):
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
    
#    # record the time the chain started
#    startTimeMAIN = time.clock() # Sadly, these start and end times don't work right in multi-processing right now ?!?!
    
    if not silent:
        print '\nMultiprocess MCMC: $$$$$$$$$$$$$ STARTED Multiprocess MCMC Sim $$$$$$$$$$$$$$$\n'
    
    # Check number of processes matches number of startParams sets passed in
    if len(startParams)<numProcesses:
        print 'Number of processes = '+str(numProcesses)+' is greater than the number of startParam sets = '+str(len(startParams))
        print 'Thus, using all random starting locations for each process'
        startParams = False
        
    ### start running threads ###
    for processNumber in range(numProcesses):
        # create temp title for this thread's temp data file
        tempProcessTitle = 'tempProcessTitle-process_'+str(processNumber+1)+'.txt'
        # pull out startParams for this process if they exist
        if startParams is not False:
            startParamsUse = startParams[processNumber]
        else:
            startParamsUse = False
        # start process
        master.append(MCMCProcessManager(tempTitle=tempProcessTitle, numSamples=numSamples, temperature=temperature, sigmaPercent=sigmaPercent, \
                     startParams=startParamsUse, silent=silent, verbose=verbose, dataTempDir=dataTempDir, dataFinalDir=dataTempDir))
        master[processNumber].start()
        
   
    # wait for completion of all process $$$$ STILL DON'T UNDERSTAND HOW THIS WORKS $$$$
    for processNumber in range(numProcesses):
        master[processNumber].join()
        
    ### put all the outputs together into a single set of data files ###
    tempFiles = []
    filesFinalFilename = dataFinalDir+filename
    for processNumber in range(numProcesses):
        tempFiles.append(dataTempDir+master[processNumber].tempTitle)
    # combine the input files into one final file
    dataFileCombiner(tempFiles, filesFinalFilename)
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
        print '\nMultiprocess MCMC: $$$$$$$$$$$$$ FINISHED Multiprocess MCMC Sim $$$$$$$$$$$$$$$\n'
        
        
        
        
        
        
        