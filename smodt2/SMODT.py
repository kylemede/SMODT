#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
#import numpy as np
import tools
import simulator
import sys
import os
import numpy as np
from multiprocessing import Process

"""
    This is the 'main' of SMODT. 
    It will start things off, call the appropriate set of 
    simulation and post-processing steps.
"""


class singleProc(Process):
    """
    This is the Manager object that controls the a single processes for a 
    SMODT2.0 simulation run.  It is called by the multiProcessStarter once for 
    each chain/process requested by the user through the simulation settings 
    file.
    
    :param str settingsDict: settings Dictionary
    :param str fNameBase: File name, including the full path, for the output 
        data files.
    :param list stageList: List of stages to run ex.['MC','SA','ST','MCMC'] 
        lives.
    :param int chainNum: number of this chain
    """
    def __init__(self, settingsDict, stageList, chainNum=1):
        
        Process.__init__(self)
        self.chainNum = chainNum
        self.log = tools.getLogger('main.singleProcess',lvl=100,addFH=False)
        self.settingsDict = settingsDict 
        self.stageList = stageList
        self.Sim = simulator.Simulator(settingsDict)
        
    def run(self):
        #$$$$$$$$$$$$$$$$$$$$$$$ TEMP $$$$$$$$$$$$$$$$$$$$$$$$$$
        ##temp params that are good enough to start an MCMC run if that is all we want to test.
        paramsST = np.array([  9.50122214e-01,   2.21670225e-01,   5.05405775e+00,
             6.07928962e+01,   4.08889781e-01,   2.45700081e+06,
             2.45700081e+06,   1.50094339e+01,   2.71264582e+01,
             1.09266595e+02,   6.41461741e+00,   3.40786762e+01,
             1.20316127e+03,   1.83558842e+00])
        sigmasST = np.array([ 0.03,  0.07,  0.13,  0.13,  0.09,  0.17,  0.01,  0.09,  0.09,
            0.21,  0.01,  0.01,  0.01,  0.05])
        bestRedChiSqr=1.0
        #$$$$$$$$$$$$$$$$$$$$$$$ TEMP $$$$$$$$$$$$$$$$$$$$$$$$$$
        
        ## run each requested stage
        self.log.info('Starting to run process #'+str(self.chainNum))
        if 'MC' in self.stageList:
            outMCFname = self.Sim.simulatorFunc('MC',self.chainNum)
            self.log.info('chain #'+str(self.chainNum)+' MC OUTFILE :\n'+outMCFname)
        if 'SA' in self.stageList:
            (paramsSA,sigmasSA,bestRedChiSqr) = self.Sim.simulatorFunc('SA',self.chainNum)
        if bestRedChiSqr<self.settingsDict['chiMAX'][0]:
            if 'ST' in self.stageList:
                (paramsST,sigmasST) = self.Sim.simulatorFunc('ST',self.chainNum,paramsSA,sigmasSA)
            if 'MCMC' in self.stageList:
                outMCMCFname = self.Sim.simulatorFunc('MCMC',self.chainNum,paramsST,sigmasST)
                self.log.info('chain #'+str(self.chainNum)+' MCMC OUTFILE :\n'+outMCMCFname)
        else:
            self.log.critical("NO ORBIT WITH A CHISQUARED < "+str(self.settingsDict['chiMAX'][0])+\
                              " WAS FOUND FOR CHAIN #"+str(self.chainNum))  
               
def smodt():
    """
    'main'
    """
    settingsDict = tools.startup(sys.argv)
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=settingsDict['logLevel'])
    log.debug("Prepend string passed in was '"+settingsDict['prepend']+"'")
             
    ##IDEA: could call the set of processes to only perform on stage at a time
    ##     Then choose the best output from all of them as the start of the next
    ##     stage.  good idea???
    
    #make list of stages to run
    stageList = []
    if settingsDict['symMode'][0]=='MC':
        stageList = ['MC']
    elif settingsDict['symMode'][0]=='SA':
        stageList = ['SA']
    elif settingsDict['symMode'][0]=='MCMC':
        stageList = ['SA','ST','MCMC']
    stageList = ['MC','SA','ST','MCMC']#$$$$$$$$$$$$$$$$$$$$$$$$$
    ## Start the number of processes/chains requested
    master = []
    log.info("Going to start "+str(settingsDict['nChains'][0])+" chains, with each running these stages: "+repr(stageList))
    for procNumber in range(settingsDict['nChains'][0]):
        master.append(singleProc(settingsDict,stageList,procNumber))
        master[procNumber].start()
    for procNumber in range(settingsDict['nChains'][0]):
        master[procNumber].join()    
    
    ## load up lists of output files
    mcmcFiles = []
    mcFiles = []
    for procNumber in range(settingsDict['nChains'][0]):
        if 'MCMC' in stageList:
            fname = os.path.join(settingsDict['finalFolder'],'outputDataMCMC'+str(procNumber)+'.fits')
            if os.path.exists(fname):
                mcmcFiles.append(fname)
        if 'MC' in stageList:
            fname = os.path.join(settingsDict['finalFolder'],'outputDataMC'+str(procNumber)+'.fits')
            if os.path.exists(fname):
                mcFiles.append(fname)
    ##################
    ## Post-processing 
    ##################
    ## combine the data files
    outMCFname = ''
    outMCMCFname = ''
    if len(mcFiles)>0:
        outMCFname = os.path.join(os.path.dirname(mcFiles[0]),"outputMC-ALL.fits")
        tools.combineFits(mcFiles,outMCFname)
    if len(mcmcFiles)>0:
        outMCMCFname = os.path.join(os.path.dirname(mcmcFiles[0]),"outputMCMC-ALL.fits")
        tools.combineFits(mcmcFiles,outMCMCFname)
    ## plot posteriors
    if settingsDict['pltDists']:
        if os.path.exists(outMCFname):
            plotFilename = os.path.join(os.path.dirname(outMCFname),'SummaryPlotMC')
            tools.summaryPlotter(outMCFname, plotFilename, shadeConfLevels=True)
        if os.path.exists(outMCMCFname):
            plotFilename = os.path.join(os.path.dirname(outMCMCFname),'SummaryPlotMCMC')
            tools.summaryPlotter(outMCMCFname, plotFilename, shadeConfLevels=True)
   
        
    log.info("End of SMODT2.0 main")

if __name__ == '__main__':
    smodt()