#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
#import numpy as np
import tools
import simulator
import sys
import os
import timeit
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
    def __init__(self, settingsDict, SimObj, stageList, chainNum=1):
        
        Process.__init__(self)
        self.chainNum = chainNum
        self.log = tools.getLogger('main.singleProcess',lvl=100,addFH=False)
        self.settingsDict = settingsDict 
        self.stageList = stageList
        self.Sim = SimObj
        
    def run(self):
        #$$$$$$$$$$$$$$$$$$$$$$$ TEMP $$$$$$$$$$$$$$$$$$$$$$$$$$
        ##temp params that are good enough to start an MCMC run if that is all we want to test.
        paramsST = np.array([9.91460200e-01,   1.98628676e-01,   5.02754554e+00,
         5.91088681e+01,   4.08044156e-01,   2.45701883e+06,
         2.45701883e+06,   1.49596676e+01,   2.99163981e+01,
         1.11421670e+02,   6.43357074e+00,   6.70551229e+01,
         1.16794139e+03,   5.20187725e+00])
        sigmasST = np.array([0.09,  0.05,  0.03,  0.01,  0.03,  0.05,  0.01,  0.05,  0.03,
        0.07,  0.01,  0.01,  0.01,  0.05])
        bestRedChiSqr=1.0
        #$$$$$$$$$$$$$$$$$$$$$$$ TEMP $$$$$$$$$$$$$$$$$$$$$$$$$$
        
        ## run each requested stage
        self.log.info('Starting to run process #'+str(self.chainNum))
        if 'MC' in self.stageList:
            outMCFname = self.Sim.simulatorFunc('MC',self.chainNum)
            self.log.info('chain #'+str(self.chainNum)+' MC OUTFILE :\n'+outMCFname)
        if 'SA' in self.stageList:
            (paramsSA,sigmasSA,bestRedChiSqr) = self.Sim.simulatorFunc('SA',self.chainNum)
        if bestRedChiSqr<self.settingsDict['cMaxMCMC'][0]:
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
    ## Call startup to get dict and load up final directories into it.
    settingsDict = tools.startup(sys.argv)
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=settingsDict['logLevel'])
    log.debug("Prepend string passed in was '"+settingsDict['prepend']+"'")
    Sim = simulator.Simulator(settingsDict)
       
    ##################################
    # Run nChains for mode requested #
    ##################################     
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ##IDEA: could call the set of processes to only perform on stage at a time
    ##     Then choose the best output from all of them as the start of the next
    ##     stage.  good idea???
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    tic=timeit.default_timer()
    if True: #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        ##make list of stages to run
        stageList = []
        if settingsDict['symMode'][0]=='MC':
            stageList = ['MC']
        elif settingsDict['symMode'][0]=='SA':
            stageList = ['SA']
        elif settingsDict['symMode'][0]=='MCMC':
            stageList = ['SA','ST','MCMC']
        #stageList=['MCMC']##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        ## Start the number of processes/chains requested
        master = []
        log.info("Going to start "+str(settingsDict['nChains'][0])+" chains, with each running these stages: "+repr(stageList))
        for procNumber in range(settingsDict['nChains'][0]):
            master.append(singleProc(settingsDict,Sim,stageList,procNumber))
            master[procNumber].start()
        for procNumber in range(settingsDict['nChains'][0]):
            master[procNumber].join()    
        toc=timeit.default_timer()
        log.info("ALL stages took a total of "+str(int(toc-tic))+' seconds')
        
        ###################
        # Post-processing # 
        ###################
        ## load up list of output files
        tic2=timeit.default_timer()
        outFiles = []
        for procNumber in range(settingsDict['nChains'][0]):
            fname = os.path.join(settingsDict['finalFolder'],'outputData'+settingsDict['symMode'][0]+str(procNumber)+'.fits')
            if os.path.exists(fname):
                outFiles.append(fname)      
        
        ## calc and strip burn-in?
        burnInStr = ''
        
        ## combine the data files
        allFname = ''
        if len(outFiles)>0:
            allFname = os.path.join(os.path.dirname(outFiles[0]),"outputData"+settingsDict['symMode'][0]+"-ALL.fits")
            tools.combineFits(outFiles,allFname)
            
        ## find best fit
        if os.path.exists(allFname):
            bestFit = tools.findBestOrbit(allFname)
            
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    if False:
        tic2=timeit.default_timer()
        allFname = '/run/media/kmede/Data1/Todai_Work/Data/data_SMODT/mockdata-SMODT2format-RV.dat'#$$$$$$$$$$$$$$$$$
        bestFit = np.array([1.07274714e+00,   1.77527538e-01,   4.96417753e+00,
             6.87573708e+01,   4.09298027e-01,   2.45701635e+06,
             2.45701635e+06,   1.49448020e+01,   3.51095143e+01,
             1.10365596e+02,   6.53591266e+00,   1.97384721e+01,
             1.16593874e+03,  -3.50058385e+00])
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ## orbit plots?
    if settingsDict['pltOrbit'] and os.path.exists(allFname):
        plotFnameBase = os.path.join(os.path.dirname(allFname),'orbitPlot'+settingsDict['symMode'][0])
        tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps')
    ## plot posteriors?
    clStr = ''
    if settingsDict['pltDists'] and os.path.exists(allFname):
        plotFilename = os.path.join(os.path.dirname(allFname),'summaryPlot'+settingsDict['symMode'][0])
        clStr = tools.summaryPlotter(allFname, plotFilename,stage=settingsDict['symMode'][0], shadeConfLevels=True)
    
    ## progress plots?
    
    ##calc R?
    grStr = ''
    if (len(outFiles)>0) and (settingsDict['CalcGR'] and (settingsDict['symMode'][0]=='MCMC')):
        (GRs,Ts,grStr) = gelmanRubinCalc(mcmcFileList)
        #grStr+="GRs = "+repr(GRs)+"\n"
        #rStr+="Ts = "+repr(Ts)+'\n'
        
    ## calc correlation length & number effective points?
    effPtsStr = ''
    if (settingsDict['symMode'][0]=='MCMC')and (settingsDict['calcCL'] and os.path.exists(allFname)):
        effPtsStr = tools.mcmcEffPtsCalc(allFname)
    
    ## make summary file
    if os.path.exists(allFname):
        summaryFname = os.path.join(os.path.dirname(allFname),'SUMMARY-'+settingsDict['symMode'][0]+".txt")
        tools.summaryFile(allFname,summaryFname,grStr,effPtsStr,clStr,burnInStr,bestFit)
    
            
    ##clean up files (move to folders or delete them)
    
    ## Final log messages and end
    toc=timeit.default_timer()
    log.info("Post-processing took a total of "+str(int(toc-tic2))+' seconds')
    log.info("\n\nEVERYTHING took a total of "+str(int(toc-tic))+' seconds\n\n')
    log.info("End of SMODT2.0 main")
    ##END MAIN 

if __name__ == '__main__':
    smodt()