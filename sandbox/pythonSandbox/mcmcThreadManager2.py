import os
import numpy as np
import time
import threading
from mcmcOrbSimulatorUniform4 import mcmcUniformOrbSim
from orbitToolbox2 import dataWriter, dataFileCombiner

class MCMCThreadManager(threading.Thread):
    '''
    This is the Manager object that controls the threading for an MCMC simulator run 
    initiated using the mcmcSimStarter function.
    '''
    def __init__(self, sema, tempTitle, numSamples, silent=True):
        threading.Thread.__init__(self)
        self.semaphore = sema
        self.tempTitle = tempTitle
        self.numSamples = numSamples
        self.silent = silent
        
    def run(self):
        self.semaphore.acquire()
        self.process(tempTitle=self.tempTitle, numSamples=self.numSamples, silent=self.silent)
        self.semaphore.release()
        
    def process(self, tempTitle, numSamples, silent=True):
        print 'Starting to run thread for file title: '+tempTitle
        
        # Run simulator for this thread
        (longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds) = mcmcUniformOrbSim(numSamples, silent=silent)
                
        # store the value of the number of epochs in self object        
        self.numEpochs = np.shape(a1s2[:][:])[1]
        
        # write data to temp file    
        dataWriter(tempTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                  ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds)
        

def mcmcSimStarter(plotFileTitle, numSamples, silent=True):
    """
    This is the master function to start a multi-threaded MCMC simulator run.
    
    """
    #dataPath =  #$$$$$$ DO SOMETHING TO MAKE THIS KINDA THING WORK, GRAB FROM PLOTFILETITLE OR...?? $$$$$$$$$$$$$
    numThreads = 6  #$$$$$ MAKE THIS AUTOMATIC BASED ON NUMBER OF PROCESSORS/THREADS AVAILABLE
    master = []
    dataPath = '../data/'
    pool_sema = threading.BoundedSemaphore(numThreads)
  
    # set verbose boolean to be opposite of silent
    verbose = True
    #if silent:
    #    verbose = False
    
    # record the time the chain started
    startTime = time.clock()
    
    print '\nMulti-Threaded MCMC: $$$$$$$$$$$$$ STARTED Multi-Threaded MCMC Sim $$$$$$$$$$$$$$$\n'
    ### start running threads ###
    for threadNumber in range(numThreads):
        # create temp title for this thread's temp data file
        tempThreadTitle = 'tempThreadTitle-thread_'+str(threadNumber+1)
        # start thread
        master.append(MCMCThreadManager(sema=pool_sema, tempTitle=tempThreadTitle, numSamples=numSamples,silent=silent))
        master[threadNumber].start()
        
    # wait for completion of all threads $$$$ STILL DON'T UNDERSTAND HOW THIS WORKS $$$$
    master[numThreads-1].join()
    
    ### put all the outputs together into a single set of data files ###
    ## first put all the 'INS' together
    INS = []
    INSfinalFilename = dataPath+plotFileTitle+'_INS.txt'
    for threadNumber in range(numThreads):
        INS.append(dataPath+master[threadNumber].tempTitle+'_INS.txt')
    # combine the input files into one final file
    dataFileCombiner(INS, INSfinalFilename)
    # kill off INS temp files
    for filename in INS:
        if os.path.exists(filename):
            os.remove(filename)
        else:
            print 'Problem: '+filename+' could not be found to remove it!!!'
    
    ## now put all the outputs together for each epoch
    # get number of epochs from first object in master list
    numEpochs = master[0].numEpochs
    for epoch in range(numEpochs):
        OUTSepoch = []
        OUTSepochfinalFilename = dataPath+plotFileTitle+'_OUTSepoch'+str(epoch+1)+'.txt'
        for threadNumber in range(numThreads):
            OUTSepoch.append(dataPath+master[threadNumber].tempTitle+'_OUTSepoch'+str(epoch+1)+'.txt')
        # combine the temp output files for this epoch into one output file
        dataFileCombiner(OUTSepoch, OUTSepochfinalFilename)
        # kill off OUTSepoch temp files
        for filename in OUTSepoch:
            if os.path.exists(filename):
                os.remove(filename)
            else:
                print 'Problem: '+filename+' could not be found to remove it!!!'
    
    # record the time the chain finished and print
    endTime = time.clock()
    totalTime = (endTime-startTime)/60
    if totalTime>60:
        totalTime = totalTime/60
        print 'Multi-Threaded MCMC: The took '+str(totalTime)+' hours to complete.\n'  ##### only print in 'silent' mode to track time
    else:
        print 'Multi-Threaded MCMC: The chain took '+str(totalTime)+' minutes to complete.\n'  ##### only print in 'silent' mode to track time
    
        
    if not silent:
        print '\nMulti-Threaded MCMC: $$$$$$$$$$$$$ FINISHED Multi-Threaded MCMC Sim $$$$$$$$$$$$$$$\n'
        
        
        
        
        
        
        