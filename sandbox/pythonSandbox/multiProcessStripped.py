import os
import numpy as np
import shutil
import time
from multiprocessing import Process

class mcProcessManager(Process):
    '''
    This is the Manager object that controls the processes for an MCMC simulator run 
    initiated using the mcmcSimStarter function.
    '''
    def __init__(self, tempTitle, numSamples):
        
        Process.__init__(self)
        self.tempTitle = tempTitle
        self.numSamples = numSamples
        
        
    def run(self):
        self.process(tempTitle=self.tempTitle, numSamples=self.numSamples)
        
    
    def process(self, tempTitle, numSamples):
        CplusplusCodeDir = "/run/media/Kyle/Data1/Todai_Work/Dropbox/workspace/C++project/"
    
        print 'Starting to run process for file title: '+tempTitle
        
        CplusplusCodeCALL = CplusplusCodeDir+"mathyTEST"+' '+str(numSamples)
        print 'CplusplusCodeCALL : ',CplusplusCodeCALL
        print "call being made"
        os.system(CplusplusCodeCALL)
        print 'call completed'
        
def mcSimStarter(paramSettingsDict):
    """
    This is the master function to start a multiprocess MCMC simulator run.
    
    """
    # TEMP!!!!!!!!!!!!!! REMOVE FOLLOWING COMMENT BLOCK ASAP!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    """
    # set simulation type
    mcONLY = paramSettingsDict['mcONLY'] # using only standard 'shotgun' monte carlo?
    
    # set simulation general settings
    numProcesses = paramSettingsDict['numProcesses']    
    silent = paramSettingsDict['silent']
    verbose = paramSettingsDict['verbose']
    data_dir = paramSettingsDict['data_dir']
    filenameRoot = paramSettingsDict['filenameRoot']  
    stripBurnIn = paramSettingsDict['stripBurnIn']
    calcCorrLengths = paramSettingsDict['calcCorrLengths']
    DIonly = paramSettingsDict['DIonly']
    RVonly = paramSettingsDict['RVonly']
    """
    #dataPath =  #$$$$$$ DO SOMETHING TO MAKE THIS KINDA THING WORK, GRAB FROM PLOTFILETITLE OR...?? $$$$$$$$$$$$$
    #numProcesses = 2  #$$$$$ MAKE THIS AUTOMATIC BASED ON NUMBER OF PROCESSORS/THREADS AVAILABLE
    master = []
    
    # Ensure the filename has .txt at the end
    if filename[-4:]!='.txt':
        filename = filename+".txt"
        
        print '\nMultiprocess: $$$$$$$$$$$$$ STARTED Multiprocess $$$$$$$$$$$$$$$\n'
    
    
    ### start running threads ###
    for processNumber in range(numProcesses):
        # create temp title for this thread's temp data file
        tempProcessTitle = 'tempProcessTitle-process_'+str(processNumber+1)+'.txt'
        
        # start process
        master.append(MCMCProcessManager(tempTitle=tempProcessTitle, numSamples=numSamples))
        master[processNumber].start()
        
    # wait for completion of all process $$$$ STILL DON'T UNDERSTAND HOW THIS WORKS $$$$
    for processNumber in range(numProcesses):
        master[processNumber].join()
        
    # just examples of how to call preset values for each processes object
    for processNumber in range(numProcesses):
        master[processNumber].tempTitle
        master[processNumber].numSamples
        
        
    print '\n**** EVERYTHING FINISHED ****\n'
    
def main():
    
    CplusplusCodeDir = "/run/media/Kyle/Data1/Todai_Work/Data/data_Binary/data_C++/"
    CplusplusFunc = 'mathyTEST'
    # compile code just to ensure it is good to go
    os.system('g++ -Wall -pedantic -lrt -o '+CplusplusCodeDir+CplusplusFunc+' '+CplusplusCodeDir+CplusplusFunc+'.cpp')
    
    filename = 'sillyWillyTest.txt'
    numSamples = int(1e6)
    numProcesses = int(4)
    
    mcmcSimStarter(filename, numSamples, numProcesses=numProcesses)
        
        
        
#############################################################
# end
#############################################################

if __name__ == '__main__':
    main()       
        
        