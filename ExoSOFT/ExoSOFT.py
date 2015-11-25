#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
#import numpy as np
import tools
import simulator
from exosoftpath import rootDir
import sys
import os
import timeit
import numpy as np
from multiprocessing import Process
import pickle


"""
    This is the 'main' of ExoSOFT. 
    It will start things off, call the appropriate set of 
    simulation and post-processing steps.
"""
class singleProc(Process):
    """
    This is the Manager object that controls the a single processes for a 
    ExoSOFT simulation run.  It is called by the multiProcessStarter once for 
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
        #bestRedChiSqrSA=1.0
        #bestRedChiSqrST=1.0
        #$$$$$$$$$$$$$$$$$$$$$$$ TEMP $$$$$$$$$$$$$$$$$$$$$$$$$$
        
        ## run each requested stage
        self.log.debug('Starting to run process #'+str(self.chainNum))
        if 'MC' in self.stageList:
            outMCFname = self.Sim.simulatorFunc('MC',self.chainNum)
            self.log.info('chain #'+str(self.chainNum)+' MC OUTFILE :\n'+outMCFname)
        if 'SA' in self.stageList:
            (paramsSA,sigmasSA,bestRedChiSqrSA) = self.Sim.simulatorFunc('SA',self.chainNum)
        if 'ST' in self.stageList:
            if bestRedChiSqrSA<self.settingsDict['chiMaxST'][0]:
                if 'ST' in self.stageList:
                    self.log.warning('chain #'+str(self.chainNum)+" made it to the ST stage, with bestRedChiSqrSA = "+str(bestRedChiSqrSA)+" :-)\n")
                    (paramsST,sigmasST,bestRedChiSqrST) = self.Sim.simulatorFunc('ST',self.chainNum,paramsSA,sigmasSA)
                if bestRedChiSqrST<self.settingsDict['cMaxMCMC'][0]:
                    print '\after: numMCMCchains = '+str(numMCMCchains)+", maxNumMCMCchains = "+str(maxNumMCMCchains)
                    if ('MCMC' in self.stageList):
                        self.log.warning('chain #'+str(self.chainNum)+" made it to the MCMC stage, with bestRedChiSqrST = "+str(bestRedChiSqrST)+" :-D\n")
                        outMCMCFname = self.Sim.simulatorFunc('MCMC',self.chainNum,paramsST,sigmasST)
                        self.log.info('chain #'+str(self.chainNum)+' MCMC OUTFILE :\n'+outMCMCFname)
                else:
                    self.log.info("NO ST SOLUTION WITH A REDUCED CHISQUARED < "+str(self.settingsDict['cMaxMCMC'][0])+\
                                      " WAS FOUND FOR CHAIN #"+str(self.chainNum)+'\n')  
            else:
                self.log.info("NO SA SOLUTION WITH A REDUCED CHISQUARED < "+str(self.settingsDict['chiMaxST'][0])+\
                                  " WAS FOUND FOR CHAIN #"+str(self.chainNum)+'\n')     

class singleProc2(Process):
    """
    This is the Manager object that controls the a single processes for a 
    ExoSOFT simulation run.  It is called by the multiProcessStarter once for 
    each chain/process requested by the user through the simulation settings 
    file.
    
    :param str settingsDict: settings Dictionary
    :param str fNameBase: File name, including the full path, for the output 
        data files.
    :param list stageList: List of stages to run ex.['MC','SA','ST','MCMC'] 
        lives.
    :param int chainNum: number of this chain
    """
    def __init__(self, settingsDict, SimObj, stage, chainNum, pklFilename = '', params=[],sigmas=[]):
        Process.__init__(self)
        self.chainNum = chainNum
        self.log = tools.getLogger('main.singleProcess',lvl=100,addFH=False)
        self.settingsDict = settingsDict 
        self.stage = stage
        self.params = params
        self.sigmas = sigmas
        self.Sim = SimObj
        self.pklFilename = pklFilename
        
    def run(self):        
        ## run each requested stage
        self.log.debug('Starting to run process #'+str(self.chainNum))
        (outFname,params,sigmas,bestRedChiSqr) = self.Sim.simulatorFunc(self.stage,self.chainNum,self.params,self.sigmas)
        self.log.debug('chain #'+str(self.chainNum)+self.stage+'  OUTFILE :\n'+outFname)
        pickle.dump([outFname,params,sigmas,bestRedChiSqr], open(self.pklFilename,'wb'))
        
def multiProc(settingsDict,Sim,stage,numProcs,params=[],sigmas=[]):
    log = tools.getLogger('main.multiProc',lvl=100,addFH=False)
    if stage in ['ST','MCMC']:
        if len(params)==0:
            log.error("params and sigmas parameter inputs to multiproc must have a set of values for each requested proc to start for ST and MCMC modes!")
    else:
        #load up basically empty arrays for params and sigmas for MC and SA cases
        sigmas = params = range(numProcs)
    master = []
    tic=timeit.default_timer()
    log.warning("Going to start "+str(numProcs)+" chains for the  "+stage+" stage")
    for procNumber in range(numProcs):
        pklFilename = os.path.join(settingsDict['finalFolder'],'pklTemp'+"-"+stage+'-'+str(procNumber)+".p")
        master.append(singleProc2(settingsDict,Sim,stage,procNumber,pklFilename=pklFilename,params=params[procNumber],sigmas=sigmas[procNumber]))
        master[procNumber].start()
    for procNumber in range(numProcs):
        master[procNumber].join()    
    toc=timeit.default_timer()
    tools.timeStrMaker(int(toc-tic))
    log.warning("ALL "+stage+" chains took a total of "+tools.timeStrMaker(int(toc-tic)))
    rets = [[],[],[],[]]
    for procNumber in range(numProcs):
        ret = pickle.load(open(master[procNumber].pklFilename,'rb'))
        for i in range(4):
            rets[i].append(ret[i])
    return rets

def exoSOFT():
    """
    'main'
    """
    ## Call startup to get dict and load up final directories into it.
    settingsDict = tools.startup(sys.argv,rootDir)
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=settingsDict['logLevel'])
    log.debug("Prepend string passed in was '"+settingsDict['prepend']+"'")
    Sim = simulator.Simulator(settingsDict)
       
    ##################################
    # Run nChains for mode requested #
    ##################################     
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ##IDEA: could call the set of processes to only perform one stage at a time
    ##      Then choose the best output from all of them as the start of the next
    ##      stage.  good idea???
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    tic=timeit.default_timer()
    ##make list of stages to run
    stgLstDict = {'MC':['MC'],'SA':['SA'],'MCMC':['SA','ST','MCMC']}
    stageList = stgLstDict[settingsDict['symMode'][0]]
    ## Start the number of processes/chains requested
#     master = []
#     log.info("Going to start "+str(settingsDict['nChains'][0])+" chains, with each running these stages: "+repr(stageList))
#     for procNumber in range(settingsDict['nChains'][0]):
#         master.append(singleProc(settingsDict,Sim,stageList,procNumber))
#         master[procNumber].start()
#     for procNumber in range(settingsDict['nChains'][0]):
#         master[procNumber].join()    
#     toc=timeit.default_timer()
#     log.info("ALL stages took a total of "+tools.timeStrMaker(int(toc-tic)))
    
    ###-------------------------------------------------------
    ###-------------------------------------------------------
    returnsSA = []
    returnsST = []
    tic=timeit.default_timer()
    maxNumMCMCprocs = settingsDict['nMCMCcns'][0]
    if 'MC' in stageList:
        returns = returnsMC = multiProc(settingsDict,Sim,'MC',settingsDict['nChains'][0])
    if 'SA' in stageList:
        returns = returnsSA = multiProc(settingsDict,Sim,'SA',settingsDict['nChains'][0])
    if ('ST' in stageList)and(len(returnsSA)>0):
        startParams = []
        startSigmas = []
        #print 'returnsSA = '+repr(returnsSA)
        for i in range(len(returnsSA[0])):
            #[outFname,params,sigmas,bestRedChiSqr]
            if returnsSA[3][i]<settingsDict['chiMaxST'][0]:
                startParams.append(returnsSA[1][i])
                startSigmas.append(returnsSA[2][i])
        if len(startSigmas)>0:
            returns = returnsST = multiProc(settingsDict,Sim,'ST',len(startSigmas),startParams,startSigmas)
    else:
        log.critical("No SA results available to start the ST chains with.")
    if ('MCMC' in stageList)and(len(returnsST)>0):
        startParams = []
        startSigmas = []
        #Filter inputs if more than max num MCMC proc available to use the best ones
        chisSorted = np.sort(returnsST[3])
        #print "# orig chisSorted = "+repr(len(chisSorted))
        chisSorted = chisSorted[np.where(chisSorted<settingsDict['cMaxMCMC'][0])]
        #print "# trimmed chisSorted = "+repr(len(chisSorted))
        if len(chisSorted)>maxNumMCMCprocs:
            chisSorted = chisSorted[:maxNumMCMCprocs]
        #print "# using chisSorted = "+repr(len(chisSorted))
        for i in range(len(returnsST[0])):
            if returnsST[3][i] in chisSorted:
                startParams.append(returnsST[1][i])
                startSigmas.append(returnsST[2][i])
        if len(chisSorted)>0:
            returns = returnsMCMC = multiProc(settingsDict,Sim,'MCMC',len(chisSorted),startParams,startSigmas)
    else:
        log.critical("No ST results available to start the MCMC chains with.")
    outFiles = returns[0]
    toc=timeit.default_timer()
    log.info("ALL stages took a total of "+tools.timeStrMaker(int(toc-tic)))
    ###-------------------------------------------------------
    ###-------------------------------------------------------
    
    ###################
    # Post-processing # 
    ###################
    log.warning("Starting Post-Processing")
    tic2=timeit.default_timer()
    ## load up list of output files
#     outFiles = []
#     for procNumber in range(settingsDict['nChains'][0]):
#         fname = os.path.join(settingsDict['finalFolder'],'outputData'+settingsDict['symMode'][0]+str(procNumber)+'.fits')
#         if os.path.exists(fname):
#             outFiles.append(fname)      
    
    ## combine the data files
    allFname = ''
    if len(outFiles)>0:
        allFname = os.path.join(os.path.dirname(outFiles[0]),"combined"+settingsDict['symMode'][0]+"data.fits")
        tools.combineFits(outFiles,allFname)
    
    ## calc and strip burn-in?
    burnInStr = ''
    if (len(outFiles)>1)and(settingsDict['CalcBurn'] and(settingsDict['symMode'][0]=='MCMC')):
        (burnInStr,burnInLengths) = tools.burnInCalc(outFiles,allFname)    
        if settingsDict['rmBurn'][0]:
            strippedFnames = tools.burnInStripper(outFiles,burnInLengths)
            outFiles = strippedFnames
            ## combine stripped files to make final file?
            if len(strippedFnames)>0:
                strippedAllFname = os.path.join(os.path.dirname(strippedFnames[0]),"combined-BIstripped-MCMCdata.fits")
                tools.combineFits(strippedFnames,strippedAllFname)
                ## replace final combined filename with new stripped version
                allFname = strippedAllFname
                
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
        clStr = tools.summaryPlotter(allFname,plotFilename,bestVals=bestFit,stage=settingsDict['symMode'][0],shadeConfLevels=True,plotALLpars=True)
    
    ##calc R?
    grStr = ''
    if (len(outFiles)>1) and (settingsDict['CalcGR'] and (settingsDict['symMode'][0]=='MCMC')):
        (GRs,Ts,grStr) = tools.gelmanRubinCalc(outFiles,settingsDict['nSamples'][0])
    
    ## progress plots?  INCLUDE?? maybe kill this one. Function exists, but not decided how to use it here.
    
    ## calc correlation length & number effective points? 
    effPtsStr = ''
    if ((len(outFiles)>1)and(settingsDict['symMode'][0]=='MCMC'))and (settingsDict['calcCL'] and os.path.exists(allFname)):
        effPtsStr = tools.mcmcEffPtsCalc(allFname)

    ## Make a summary file of results 
    toc=timeit.default_timer()
    postTime = toc-tic2
    allTime = toc-tic
    if os.path.exists(allFname):
        tools.summaryFile(settingsDict,stageList,allFname,clStr,burnInStr,bestFit,grStr,effPtsStr,allTime,postTime)
    
    ##clean up files (move to folders or delete them)
    tools.cleanUp(settingsDict,stageList,allFname)
    if settingsDict['CopyToDB']:
        tools.copyToDB(settingsDict)
        
    ## Final log messages and end
    log.info("Post-processing took a total of "+tools.timeStrMaker(postTime))
    log.info("\n\nEVERYTHING took a total of "+tools.timeStrMaker(allTime)+'\n\n')
    log.info("End of ExoSOFT main")
    ##END MAIN 

if __name__ == '__main__':
    exoSOFT()