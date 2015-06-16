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
        bestRedChiSqrSA=1.0
        bestRedChiSqrST=1.0
        #$$$$$$$$$$$$$$$$$$$$$$$ TEMP $$$$$$$$$$$$$$$$$$$$$$$$$$
        
        ## run each requested stage
        self.log.debug('Starting to run process #'+str(self.chainNum))
        if 'MC' in self.stageList:
            outMCFname = self.Sim.simulatorFunc('MC',self.chainNum)
            self.log.info('chain #'+str(self.chainNum)+' MC OUTFILE :\n'+outMCFname)
        if 'SA' in self.stageList:
            (paramsSA,sigmasSA,bestRedChiSqrSA) = self.Sim.simulatorFunc('SA',self.chainNum)
        if bestRedChiSqrSA<self.settingsDict['chiMaxST'][0]:
            if 'ST' in self.stageList:
                self.log.warning('chain #'+str(self.chainNum)+" made it to the ST stage, with bestRedChiSqrSA = "+str(bestRedChiSqrSA)+" :-)")
                (paramsST,sigmasST,bestRedChiSqrST) = self.Sim.simulatorFunc('ST',self.chainNum,paramsSA,sigmasSA)
            if bestRedChiSqrST<self.settingsDict['cMaxMCMC'][0]:
                if 'MCMC' in self.stageList:
                    self.log.warning('chain #'+str(self.chainNum)+" made it to the MCMC stage, with bestRedChiSqrST = "+str(bestRedChiSqrST)+" :-D")
                    outMCMCFname = self.Sim.simulatorFunc('MCMC',self.chainNum,paramsST,sigmasST)
                    self.log.info('chain #'+str(self.chainNum)+' MCMC OUTFILE :\n'+outMCMCFname)
            else:
                self.log.warning("NO ST SOLUTION WITH A REDUCED CHISQUARED < "+str(self.settingsDict['cMaxMCMC'][0])+\
                                  " WAS FOUND FOR CHAIN #"+str(self.chainNum)+'\n')  
        else:
            self.log.warning("NO SA SOLUTION WITH A REDUCED CHISQUARED < "+str(self.settingsDict['chiMaxST'][0])+\
                              " WAS FOUND FOR CHAIN #"+str(self.chainNum)+'\n')  
               
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
    ##IDEA: could call the set of processes to only perform one stage at a time
    ##     Then choose the best output from all of them as the start of the next
    ##     stage.  good idea???
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    tic=timeit.default_timer()
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
    
    ## combine the data files
    allFname = ''
    if len(outFiles)>0:
        allFname = os.path.join(os.path.dirname(outFiles[0]),"combined"+settingsDict['symMode'][0]+"data.fits")
        tools.combineFits(outFiles,allFname)
    
    ## calc and strip burn-in?
    burnInStr = ''
    strippedAllFname = ''
    if (len(outFiles)>1)and(settingsDict['CalcBurn'] and(settingsDict['symMode'][0]=='MCMC')):
        (burnInStr,burnInLengths) = tools.burnInCalc(outFiles,allFname)    
        if settingsDict['delBurn'][0]:
            strippedFnames = tools.burnInStripper(outFiles,burnInLengths)
            ## combine stripped files to make final file?
            if len(strippedFnames)>0:
                strippedAllFname = os.path.join(os.path.dirname(strippedFnames[0]),"combined-BIstripped-MCMCdata.fits")
                tools.combineFits(strippedFnames,strippedAllFname)
    
    ## find best fit
    if os.path.exists(allFname):
        bestFit = tools.findBestOrbit(allFname)
            
    ## orbit plots?
    if settingsDict['pltOrbit'] and os.path.exists(allFname):
        plotFnameBase = os.path.join(os.path.dirname(allFname),'orbitPlot'+settingsDict['symMode'][0])
        tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps')
    ## plot posteriors?
    clStr = ''
    if settingsDict['pltDists'] and os.path.exists(allFname):
        plotFilename = os.path.join(os.path.dirname(allFname),'summaryPlot'+settingsDict['symMode'][0])
        clStr = tools.summaryPlotter(allFname, plotFilename,stage=settingsDict['symMode'][0], shadeConfLevels=True)
    
    ## progress plots?  INCLUDE?? maybe kill this one.
    
    ## calc correlation length & number effective points?
    effPtsStr = ''
    if ((len(outFiles)>1)and(settingsDict['symMode'][0]=='MCMC'))and (settingsDict['calcCL'] and os.path.exists(allFname)):
        effPtsStr = tools.mcmcEffPtsCalc(allFname)
        
    ##calc R?
    grStr = ''
    if (len(outFiles)>1) and (settingsDict['CalcGR'] and (settingsDict['symMode'][0]=='MCMC')):
        (GRs,Ts,grStr) = tools.gelmanRubinCalc(outFiles,settingsDict['nSamples'][0])
        
    ## make summary file
    toc=timeit.default_timer()
    postTime = toc-tic2
    allTime = toc-tic
    if os.path.exists(allFname):
        tools.summaryFile(settingsDict,stageList,allFname,grStr,effPtsStr,clStr,burnInStr,bestFit,allTime,postTime)
        
    ##clean up files (move to folders or delete them)
    tools.cleanUp(settingsDict,stageList,allFname)
    if settingsDict['CopyToDB']:
        tools.copyToDB(settingsDict)
    ## Final log messages and end
    log.info("Post-processing took a total of "+tools.timeStrMaker(postTime))
    log.info("\n\nEVERYTHING took a total of "+tools.timeStrMaker(allTime)+'\n\n')
    log.info("End of SMODT2.0 main")
    ##END MAIN 

if __name__ == '__main__':
    smodt()