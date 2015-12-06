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
    def __init__(self, settingsDict, SimObj, stage, chainNum, pklFilename = '', params=[],sigmas=[],strtTemp=1.0):
        Process.__init__(self)
        self.chainNum = chainNum
        self.log = tools.getLogger('main.singleProcess',lvl=100,addFH=False)
        self.settingsDict = settingsDict 
        self.stage = stage
        self.params = params
        self.sigmas = sigmas
        self.strtTemp = strtTemp
        self.Sim = SimObj
        self.pklFilename = pklFilename
        
    def run(self):        
        ## run the requested stage and [ickle its return values
        self.log.debug('Starting to run process #'+str(self.chainNum))
        (outFname,params,sigmas,bestRedChiSqr) = self.Sim.simulatorFunc(self.stage,self.chainNum,self.params,self.sigmas,self.strtTemp)
        self.log.debug('chain #'+str(self.chainNum)+" of "+self.stage+' stage  OUTFILE :\n'+outFname)
        pickle.dump([outFname,params,sigmas,bestRedChiSqr], open(self.pklFilename,'wb'))
        
def multiProc(settingsDict,Sim,stage,numProcs,params=[],sigmas=[],strtTemp=1.0):
    log = tools.getLogger('main.multiProc',lvl=100,addFH=False)
    if stage in ['ST','MCMC']:
        if len(params)==0:
            log.error("params and sigmas parameter inputs to multiproc must have a set of values for each requested proc to start for ST and MCMC modes!")
    else:
        #load up basically empty arrays for params and sigmas for MC and SA cases
        sigmas = params = range(numProcs)
    master = []
    tic=timeit.default_timer()
    extra = ''
    if stage=='SA':
        extra+=" with a starting temperature of "+str(strtTemp)
    log.warning("Going to start "+str(numProcs)+" chains for the "+stage+" stage"+extra)
    for procNumber in range(numProcs):
        pklFilename = os.path.join(settingsDict['finalFolder'],'pklTemp'+"-"+stage+'-'+str(procNumber)+".p")
        master.append(singleProc(settingsDict,Sim,stage,procNumber,pklFilename=pklFilename,params=params[procNumber],sigmas=sigmas[procNumber],strtTemp=strtTemp))
        master[procNumber].start()
    for procNumber in range(numProcs):
        master[procNumber].join()    
    toc=timeit.default_timer()
    s = "ALL "+str(numProcs)+" chains of the "+stage+" stage took a total of "+tools.timeStrMaker(int(toc-tic))
    retStr =s+"\n"
    log.warning(s)
    retAry = [[],[],[],[]]
    for procNumber in range(numProcs):
        ret = pickle.load(open(master[procNumber].pklFilename,'rb'))
        for i in range(4):
            retAry[i].append(ret[i])
    return (retAry,retStr)

def iterativeSA(settingsDict,Sim):
    """
    Perform SA with multiProc nSAiters times, droping the starting temperature each time by strtTemp/nSAiters.
    """
    log = tools.getLogger('main.iterativeSA',lvl=100,addFH=False)
    maxNumMCMCprocs = settingsDict['nMCMCcns'][0]
    numProcs = settingsDict['nChains'][0]
    nSAiters = settingsDict['nSAiters'][0]
    strtT = settingsDict['strtTemp'][0]
    strtPars = range(0,numProcs)
    bestRetAry = [[],[],[],[]]
    for iter in range(0,nSAiters):
        if iter>0:
            #clean up previous SA data files on disk to avoid clash
            #[outFname,params,sigmas,bestRedChiSqr]
            tools.rmFiles(retAry[0][:])
        temp = strtT-(strtT/float(nSAiters))*float(iter)
        #print 'repr(strtPars) = '+repr(strtPars)
        (retAry,retStr) = multiProc(settingsDict,Sim,'SA',numProcs,params=strtPars,sigmas=[],strtTemp=temp)
        if len(retAry)>0:
            chisSorted = [] 
            goodParams = []           
            #Filter inputs if more than max num MCMC proc available to use the best ones
            chisSorted = np.sort(retAry[3])
            chisSorted = chisSorted[np.where(chisSorted<settingsDict['cMaxMCMC'][0])]
            #print 'Full list of chiSquareds back from SA '+repr(chisSorted)
            #first updated bestRetAry
            for i in range(len(retAry[0])):
                if retAry[3][i] in chisSorted:                       
                    if len(bestRetAry[0])<numProcs:
                        for j in range(4):
                            bestRetAry[j].append(retAry[j][i])
                    else:
                        bestChis = np.sort(bestRetAry[3])
                        if retAry[3][i]<bestChis[-1]:
                            for j in range(4):
                                bestRetAry[j].append(retAry[j][i])
            #now make list of best ones to use in next round
            if len(bestRetAry[0])>numProcs:
                bestChis = np.sort(bestRetAry[3])
                bestRetAry2 = [[],[],[],[]]
                #trim best lists down to size
                for i in range(0,len(bestChis)):
                    if bestRetAry[3][i] in bestChis[:numProcs]:
                        for j in range(4):
                            bestRetAry2[j].append(bestRetAry[j][i])
                bestRetAry = bestRetAry2
            #print 'repr(np.sort(bestRetAry[3])) = '+repr(np.sort(bestRetAry[3]))
            goodParams = bestRetAry[1]
            log.info("Iteration #"+str(iter+1)+" resulted in the chiSquareds "+repr(chisSorted))
            ## Now fill out an array of starting parameter sets from the best above.
            ## first load up with one set of goodParams, then randomly from it till full.
            strtPars=[]
            if len(goodParams)>1:
                for i in range(0,len(goodParams)):
                    strtPars.append(goodParams[i])
                while len(strtPars)<numProcs:
                    rndVal = np.random.randint(0,len(goodParams))
                    #print 'rndVal = '+str(rndVal)
                    strtPars.append(goodParams[rndVal])
    return (retAry,retStr)

def exoSOFT():
    """
    'main'
    """
    ## Call startup to get dict and load up final directories into it.
    settingsDict = tools.startup(sys.argv,rootDir)
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=settingsDict['logLevel'])
    log.debug("Prepend string passed in was '"+settingsDict['prepend']+"'")
    Sim = simulator.Simulator(settingsDict)
       
    ###########################################
    # Run nChains for MC/SA/ST mode requested #
    #  Then up to nMCMCcns if MCMC requested  #
    ###########################################     
    tic=timeit.default_timer()
    stageList = settingsDict['stageList']
    maxNumMCMCprocs = settingsDict['nMCMCcns'][0]
    durationStrings = ''
    if 'MC' in stageList:
        (returns,b) = (returnsMC,durStr) = multiProc(settingsDict,Sim,'MC',settingsDict['nChains'][0])
        durationStrings+=durStr
    if 'SA' in stageList:
        (returns,b) = (returnsSA,durStr) = iterativeSA(settingsDict,Sim)
        durationStrings+=durStr
    if 'ST' in stageList:
        startParams = []
        startSigmas = []
        if settingsDict['stages'] in ['ST','STMCMC']:
            for i in range(0,maxNumMCMCprocs):
                startParams.append(settingsDict['startParams'])
                startSigmas.append(settingsDict['startSigmas'])
        elif len(returnsSA)>0:
            for i in range(len(returnsSA[0])):
                if returnsSA[3][i]<settingsDict['chiMaxST'][0]:
                    startParams.append(returnsSA[1][i])
                    startSigmas.append(returnsSA[2][i])
        else:
            log.critical("No SA results available to start the ST chains with.")
        if len(startSigmas)>0:
            (returns,b) = (returnsST,durStr) = multiProc(settingsDict,Sim,'ST',len(startSigmas),startParams,startSigmas)
            durationStrings+=durStr
        # check best results of ST and store to a file.
        # Maybe replace pars and sigs in original settings files?
        if len(returnsST[0])>0:
            bstChiSqr = np.sort(returnsST[3])[0]
            for i in range(len(returnsST[0])):
                if returnsST[3][i] == bstChiSqr:
                    bestSTpars = returnsST[1][i]
                    bestSTsigs = returnsST[2][i]
            tools.writeBestSTtoFile(settingsDict,bestSTpars,bestSTsigs,bstChiSqr)
            tools.pushIntoOrigSettFiles(settingsDict,bestSTpars,sigs=bestSTsigs)
        
    if 'MCMC' in stageList:
        startParams = []
        startSigmas = []
        if settingsDict['stages']=='MCMC':
            chisSorted = range(0,maxNumMCMCprocs)
            for i in range(0,maxNumMCMCprocs):
                startParams.append(settingsDict['startParams'])
                startSigmas.append(settingsDict['startSigmas'])
        elif len(returnsST)>0:
            chisSorted = []            
            #Filter inputs if more than max num MCMC proc available to use the best ones
            chisSorted = np.sort(returnsST[3])
            chisSorted = chisSorted[np.where(chisSorted<settingsDict['cMaxMCMC'][0])]
            if len(chisSorted)>maxNumMCMCprocs:
                chisSorted = chisSorted[:maxNumMCMCprocs]
            for i in range(len(returnsST[0])):
                if returnsST[3][i] in chisSorted:
                    startParams.append(returnsST[1][i])
                    startSigmas.append(returnsST[2][i])
        else:
            log.critical("No ST results available to start the MCMC chains with.")
        if len(chisSorted)>0:
            (returns,b) = (returnsMCMC,durStr) = multiProc(settingsDict,Sim,'MCMC',len(chisSorted),startParams,startSigmas)
            durationStrings+=durStr
            # Maybe replace pars in original settings files?
            bstChiSqr = np.sort(returnsMCMC[3])[0]
            for i in range(len(returnsMCMC[0])):
                if returnsMCMC[3][i] == bstChiSqr:
                    bestMCMCpars = returnsMCMC[1][i]
            tools.pushIntoOrigSettFiles(settingsDict,bestMCMCpars,sigs=[])
    outFiles = returns[0]
    toc=tic2=timeit.default_timer()
    s = "ALL stages took a total of "+tools.timeStrMaker(int(toc-tic))
    durationStrings+=s+'\n'
    log.info(s)
    
    ###################
    # Post-processing # 
    ###################
    log.warning("Starting Post-Processing")  
    
    #figure out
    ## combine the data files
    allFname = ''
    if len(outFiles)>0:
        allFname = os.path.join(os.path.dirname(outFiles[0]),"combined"+stageList[-1]+"data.fits")
        tools.combineFits(outFiles,allFname)
    
    ## calc and strip burn-in?
    burnInStr = ''
    if (len(outFiles)>1)and(settingsDict['CalcBurn'] and ('MCMC' in stageList)):
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
        plotFnameBase = os.path.join(os.path.dirname(allFname),'orbitPlot'+stageList[-1])
        tools.orbitPlotter(bestFit,settingsDict,plotFnameBase,format='eps')
    
    ## plot posteriors?
    clStr = ''
    if settingsDict['pltDists'] and os.path.exists(allFname):
        plotFilename = os.path.join(os.path.dirname(allFname),'summaryPlot'+stageList[-1])
        clStr = tools.summaryPlotter(allFname,plotFilename,bestVals=bestFit,stage=settingsDict['stages'],shadeConfLevels=True,plotALLpars=True)
    
    ##calc R?
    grStr = ''
    if (len(outFiles)>1) and (settingsDict['CalcGR'] and ('MCMC' in stageList)):
        (GRs,Ts,grStr) = tools.gelmanRubinCalc(outFiles,settingsDict['nSamples'][0])
    
    ## progress plots?  INCLUDE?? maybe kill this one. Function exists, but not decided how to use it here.
    
    ## calc correlation length & number effective points? 
    effPtsStr = ''
    if ((len(outFiles)>1)and('MCMC' in stageList))and (settingsDict['calcCL'] and os.path.exists(allFname)):
        effPtsStr = tools.mcmcEffPtsCalc(allFname)

    ## Make a summary file of results 
    toc=timeit.default_timer()
    postTime = toc-tic2
    allTime = toc-tic
    if os.path.exists(allFname):
        tools.summaryFile(settingsDict,stageList,allFname,clStr,burnInStr,bestFit,grStr,effPtsStr,allTime,postTime,durationStrings)
    
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