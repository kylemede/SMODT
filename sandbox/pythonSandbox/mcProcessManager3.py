import os
import numpy as np
import time
from multiprocessing import Process
from mcmcOrbSimulatorUniform4 import mcmcUniformOrbSim
from mcmcOrbSimulatorTriangular4 import mcmcTriangularOrbSim
#from mcONLYOrbSimulatorUniform4 import mcONLYUniformOrbSim
from mcONLYOrbSimulatorUniform1 import mcONLYUniformOrbSim
#from orbitToolbox import dataWriter
#from orbitToolbox2 import dataFileCombiner
from orbitToolbox2 import dataWriter, dataFileCombiner

class MCMCProcessManager(Process):
    '''
    This is the Manager object that controls the processes for an MCMC simulator run 
    initiated using the mcmcSimStarter function.
    '''
    def __init__(self, tempTitle, numSamples, silent=True, lowRAM=False):
        Process.__init__(self)
        self.tempTitle = tempTitle
        self.numSamples = numSamples
        self.silent = silent
        self.lowRAM = lowRAM
        
    def run(self):
        self.process(tempTitle=self.tempTitle, numSamples=self.numSamples, silent=self.silent, lowRAM=self.lowRAM)
    
    def process(self, tempTitle, numSamples, silent=True, lowRAM=False):
        if not silent:
            print 'Starting to run process for file title: '+tempTitle
        
        ## Run simulator for this process
        #(longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
        #        ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds) = mcmcUniformOrbSim(numSamples, silent=silent)
        #(longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
        #        ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds) = mcmcTriangularOrbSim(numSamples, silent=silent)
        (longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds) = \
                                    mcONLYUniformOrbSim(numSamples, silent=silent, lowRAM=lowRAM)
        
        
        # write data to temp file    
        dataWriter(tempTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                  ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds)
        
        # store the value of the number of epochs to a temp file in same place as rest of data  
        numEpochs =  np.shape(a1s2[:][:])[1]
        #print ' number of epochs = '+str(numEpochs)
        numEpochsFilename = tempTitle[:-1]+'numEpochs.txt'
        #print 'writing numEpochs to file '+numEpochsFilename
        f = open(numEpochsFilename,'w')
        f.write(str(numEpochs))
        f.close()

def mcmcSimStarter(plotFileTitle, numSamples, numProcesses=2, silent=True, lowRAM=False):
    """
    This is the master function to start a multiprocess MCMC simulator run.
    
    """
    #dataPath =  #$$$$$$ DO SOMETHING TO MAKE THIS KINDA THING WORK, GRAB FROM PLOTFILETITLE OR...?? $$$$$$$$$$$$$
    #numProcesses = 2  #$$$$$ MAKE THIS AUTOMATIC BASED ON NUMBER OF PROCESSORS/THREADS AVAILABLE
    master = []
    dataPath = '../data/'
    # set verbose boolean to be opposite of silent
    verbose = True
    #if silent:
    #    verbose = False
    
    # record the time the chain started
    startTimeMAIN = time.clock()
    
    if not silent:
        print '\nMultiprocess MCMC: $$$$$$$$$$$$$ STARTED Multiprocess MCMC Sim $$$$$$$$$$$$$$$\n'
        
    ### start running threads ###
    for processNumber in range(numProcesses):
        # create temp title for this thread's temp data file
        tempProcessTitle = 'tempProcessTitle-process_'+str(processNumber+1)
        # start process
        master.append(MCMCProcessManager(tempTitle=tempProcessTitle, numSamples=numSamples,silent=silent,lowRAM=lowRAM))
        master[processNumber].start()
        
   
    # wait for completion of all process $$$$ STILL DON'T UNDERSTAND HOW THIS WORKS $$$$
    for processNumber in range(numProcesses):
        master[processNumber].join()
        
    ### put all the outputs together into a single set of data files ###
    ## first put all the 'INS' together
    INS = []
    INSfinalFilename = dataPath+plotFileTitle+'_INS.txt'
    for processNumber in range(numProcesses):
        INS.append(dataPath+master[processNumber].tempTitle+'_INS.txt')
    # combine the input files into one final file
    dataFileCombiner(INS, INSfinalFilename)
    # kill off INS temp files
    for filename in INS:
        if os.path.exists(filename):
            os.remove(filename)
        else:
            print 'Problem: '+filename+' could not be found to remove it!!!'
    
    ## now put all the outputs together for each epoch
    # get number of epochs from first temp numEpochs file
    numEpochsFilename = tempProcessTitle[:-1]+'numEpochs.txt'
    if os.path.exists(numEpochsFilename):
        f = open(numEpochsFilename,'r')
        numEpochs = int(f.readlines()[0])
        f.close()
        os.remove(numEpochsFilename)
        if not silent:
            print 'numEpochs value of '+str(numEpochs)+' retrieved and temp file deleted'
    else:
        print 'PROBLEM: temp numEpochs file '+numEpochsFilename+' does not exist!!!'
    # combine the files made for each epoch during all processes    
    for epoch in range(numEpochs):
        OUTSepoch = []
        OUTSepochfinalFilename = dataPath+plotFileTitle+'_OUTSepoch'+str(epoch+1)+'.txt'
        for processNumber in range(numProcesses):
            OUTSepoch.append(dataPath+master[processNumber].tempTitle+'_OUTSepoch'+str(epoch+1)+'.txt')
        # combine the temp output files for this epoch into one output file
        dataFileCombiner(OUTSepoch, OUTSepochfinalFilename)
        # kill off OUTSepoch temp files
        for filename in OUTSepoch:
            if os.path.exists(filename):
                os.remove(filename)
            else:
                print 'Problem: '+filename+' could not be found to remove it!!!'
    
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
        
        
        
        
        
        
        