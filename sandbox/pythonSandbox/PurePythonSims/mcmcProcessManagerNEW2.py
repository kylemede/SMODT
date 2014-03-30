import os
import numpy as np
import shutil
import time
from multiprocessing import Process
from mcmc3DorbSimUniform import mcmcOrbSimUniform
from mcmc3DorbSimUniformTH_I import mcmcOrbSimUniformTH_I
#from mcmcOrbSimUniform7_3 import mcmcOrbSimUniform
from orbitToolbox import burnInCalc3, dataReaderNEW, burnInStripper, dataFileCombiner, dataReadAndPlotNEW3, CorrLengthCalc, bestOrbitFinderNEW, orbitEllipsePlotter3

class MCMCProcessManager(Process):
    '''
    This is the Manager object that controls the processes for an MCMC simulator run 
    initiated using the mcmcSimStarter function.
    '''
    def __init__(self, tempTitle, numSamples, temperature=1.0, sigmaPercent=1.0, startParams=False, silent=True, verbose=False,\
                        dataTempDir='/home/Kyle/DataTempForSims/', dataFinalDir='/home/Kyle/DataTempForSims/', TH_I=False):
        
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
        self.TH_I = TH_I
        
    def run(self):
        self.process(tempTitle=self.tempTitle, numSamples=self.numSamples, temperature=self.temperature, sigmaPercent=self.sigmaPercent, \
                     startParams=self.startParams, silent=self.silent, verbose=self.verbose, dataTempDir=self.dataTempDir, dataFinalDir=self.dataFinalDir, TH_I=self.TH_I)
        
    
    def process(self, tempTitle, numSamples, temperature=1.0, sigmaPercent=1.0, startParams=False, silent=True, verbose=False,\
                        dataTempDir='/home/Kyle/DataTempForSims/', dataFinalDir='/home/Kyle/DataTempForSims/', TH_I=False):
        
        if not silent:
            print 'Starting to run process for file title: '+tempTitle
        
        ## Run simulator for this process
        if TH_I:
            mcmcOrbSimUniformTH_I(tempTitle, numSamples, temperature=temperature, sigmaPercent=sigmaPercent, startParams=startParams, silent=silent,\
                          verbose=verbose, dataTempDir=dataTempDir, dataFinalDir=dataFinalDir)
        else:
            mcmcOrbSimUniform(tempTitle, numSamples, temperature=temperature, sigmaPercent=sigmaPercent, startParams=startParams, silent=silent,\
                          verbose=verbose, dataTempDir=dataTempDir, dataFinalDir=dataFinalDir)

def mcmcSimStarter(filename, numSamples, numProcesses=1, temperature=1.0, sigmaPercent=1.0, \
                    startParams=False, silent=True, verbose=False, stripBurnIn=False, calcCorrLengths=False,\
                    dataTempDir='/home/Kyle/DataTempForSims/', dataFinalDir='/media/Data1/Todai_Work/workspace/Binary-project/data/', TH_I=False):
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
        print '\nMultiprocess MCMC: $$$$$$$$$$$$$ STARTED Multiprocess MCMC Sim $$$$$$$$$$$$$$$\n'
    
#    # Check number of processes matches number of startParams sets passed in
#    if len(startParams)<numProcesses:
#        print 'Number of processes = '+str(numProcesses)+' is greater than the number of startParam sets = '+str(len(startParams))
#        print 'Thus, using all random starting locations for each process'
#        startParams = False
    
    ### start running threads ###
    for processNumber in range(numProcesses):
        # create temp title for this thread's temp data file
        tempProcessTitle = 'tempProcessTitle-process_'+str(processNumber+1)+'.txt'
        # pull out startParams for this process if they exist
        if startParams is not False:
            if len(startParams)>=numProcesses:
                startParamsUse = startParams[processNumber]
            elif (len(startParams)>=1.0) and (len(startParams)<numProcesses):
                print 'Only using first set of start parameters as there are less of them than the number of chains requested'
                startParamsUse = startParams[0]
            else:
                print 'No start parameters provided, so starting at random locations'
                startParamsUse = False
        else:
            startParamsUse = False
        # start process
        master.append(MCMCProcessManager(tempTitle=tempProcessTitle, numSamples=numSamples, temperature=temperature, sigmaPercent=sigmaPercent, \
                     startParams=startParamsUse, silent=silent, verbose=verbose, dataTempDir=dataTempDir, dataFinalDir=dataTempDir, TH_I=TH_I))
        master[processNumber].start()
        
    # wait for completion of all process $$$$ STILL DON'T UNDERSTAND HOW THIS WORKS $$$$
    for processNumber in range(numProcesses):
        master[processNumber].join()
        
    # collect all temp filenames together and write to 'burn in included' final files in finalFolder
    tempFiles = []
    burnedInFiles = []
#    burnIns = []
#    burnInChainNumbers = []
    finalFolder = dataFinalDir+filename[:-4]+'/'
    dataFinalFilename = finalFolder+'outputData-WITHOUT-Burn-in-ALL.txt'#filename[:-4]+'-WITHOUT-Burn-in-ALL.txt'
    
    print '*Starting to write original data and resultsDicts to final files AND possibly* \n* calculate the burn-in and correlation length of all the data in the chains *\n'
        
    if stripBurnIn:
        ## Calculate the burn in length and strip from files
        # FIRST calculate the median burn in from each chain then 
        # median for all together
        
        # Create a final name for the combination of all the stripped chains
        dataFinalFilename = finalFolder+'outputData-WITHOUT-Burn-in-ALL.txt'#filename[:-4]+'-WITHOUT-Burn-in-ALL.txt'
        
        mediansPerChain = []
        for processNumber in range(numProcesses):
            tmpFileWithDir = dataTempDir+master[processNumber].tempTitle
            chiSquareds = dataReaderNEW(tmpFileWithDir,7)
            median = np.median(chiSquareds)
            del chiSquareds
            mediansPerChain.append(median)
            
        medianALLchains = np.median(mediansPerChain)
        
        for processNumber in range(numProcesses):
            burnInSuccessful = False
            tmpFileWithDir = dataTempDir+master[processNumber].tempTitle
            tempFiles.append(tmpFileWithDir)
            filenameWithBurnIN = finalFolder+'outputData-WITH-Burn-in-chain'+str(processNumber+1)+'.txt'#filename[:-4]+'-WITH-Burn-in-chain'+str(processNumber+1)+'.txt'
            shutil.copy(tmpFileWithDir,filenameWithBurnIN)
            print "Original version of this chain's data written to "+filenameWithBurnIN
            # retrieve results dictionary file and move it to final location and name 
            resultsDictTempFilename = tmpFileWithDir[:-4]+'-resultsDict.txt'
            resultsDictFinalFilename = finalFolder+'results-chain'+str(processNumber+1)+'.txt'
            shutil.move(resultsDictTempFilename,resultsDictFinalFilename)
            # open final results dict file and update with burn in info for each param in loop below
            resultsDictFile = open(resultsDictFinalFilename,'a')
            
            # Calculate the burn-in length for this chain
            resultsDictFile.write('/n/n***** Burn-in calculation information for this chain *****')
            resultsDictFile.write('\nChain '+str(processNumber+1)+' burn in length:')
            print '\nChain '+str(processNumber+1)+' burn in length:'
            chiSquareds = dataReaderNEW(tmpFileWithDir,7)
            burnInLength = burnInCalc3(chiSquareds,medianALLchains)
    
            resultsDictFile.write('\n burn in length = '+str(burnInLength))
            print 'burn in length = '+str(burnInLength)
               
            if burnInMAX<paramLength:
                burnInSuccessful = True
            
            resultsDictFile.write('Total param lengths = '+str(paramLength))
            if burnInSuccessful:
                burnInStrippedFilename = finalFolder+'outputData-WITHOUT-Burn-in-chain'+str(processNumber+1)+'.txt'#filename[:-4]+'-WITHOUT-Burn-in-chain'+str(processNumber+1)+'.txt'
                burnInStripper(tmpFileWithDir, burnInMAX, burnInStrippedFilename)
                burnedInFiles.append(burnInStrippedFilename)
                resultsDictFile.write('\n\nBurn in stripped data written to '+burnInStrippedFilename)
            
            print 'Results of chain written to '+resultsDictFinalFilename
            resultsDictFile.close()    
    
    else:
        # Create a final name for the combination of all the stripped chains
        dataFinalFilename = finalFolder+'outputData-ALL.txt'#filename[:-4]+'-ALL.txt'
        
        # Don't strip the burn in, thus just use the original data.
        # So, just copy temp data files to final for each chain.
        for processNumber in range(numProcesses):
            tmpFileWithDir = dataTempDir+master[processNumber].tempTitle
            tempFiles.append(tmpFileWithDir)
            finalChainFilename = finalFolder+'outputData-chain'+str(processNumber+1)+'.txt'#filename[:-4]+'-chain'+str(processNumber+1)+'.txt'
            shutil.copy(tmpFileWithDir,finalChainFilename)
            print "Original version of this chain's data written to "+finalChainFilename
            # retrieve results dictionary file and move it to final location and name 
            resultsDictTempFilename = tmpFileWithDir[:-4]+'-resultsDict.txt'
            resultsDictFinalFilename = finalFolder+'results-chain'+str(processNumber+1)+'.txt'
            shutil.move(resultsDictTempFilename,resultsDictFinalFilename)

    ## Calculate the correlation lengths of each parameter in each chain 
    ## and save info to the results dictionaries for later reference
    if calcCorrLengths:
        if stripBurnIn:
            dataFiles = burnedInFiles
        else:
            dataFiles = tempFiles
        
        for processNumber in range(numProcesses):
            dataFile = dataFiles[processNumber]
            resultsDictFinalFilename = finalFolder+'results-chain'+str(processNumber+1)+'.txt'
            # open final results dict file and update with Correlation length info for each param in loop below
            resultsDictFile = open(resultsDictFinalFilename,'a')
            # Calculate the correlation length for each param in this chain
            corrLengthMAX = 0
            corrLengthMAXparamNumber = 0
            resultsDictFile.write('/n/n***** Correlation Length calculation information for this chain *****')
            resultsDictFile.write('\nChain '+str(processNumber+1)+' Correlation Lengths:')
            print '\nChain '+str(processNumber+1)+' Correlation Lengths:'
            for param in range(0,7):
                paramData = dataReaderNEW(tmpFileWithDir,param)
                paramLength = len(paramData)
                corrLength = CorrLengthCalc(paramData)
                del paramData
                resultsDictFile.write('\nParam '+str(param)+' burn in length = '+str(corrLength))
                print 'Param '+str(param)+' correlation length = '+str(corrLength)
                if corrLength > corrLengthMAX:
                    corrLengthMAX = corrLength
                    corrLengthMAXparamNumber = param
                    
            maxString = '\nMax Correlation Length found to be = '+str(corrLengthMAX)+\
                        ' for zero-indexed param number, '+str(corrLengthMAXparamNumber)
            print maxString
            resultsDictFile.write('\n'+maxString) 
            resultsDictFile.write('Total param lengths = '+str(paramLength))
            # Done this chain, so close dict 
            resultsDictFile.close()

    ## Combine the chain files that burned in successfully (OR just ALL the chain files) into one final file
    if stripBurnIn:        
        print '\nNow combining all stripped files into one final file'
        dataFileCombiner(burnedInFiles, dataFinalFilename)
    else:
        print '\nNow combining all files into one final file'
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
        print '\nMultiprocess MCMC: $$$$$$$$$$$$$ FINISHED Multiprocess MCMC Sim $$$$$$$$$$$$$$$\n'
    
    print '**** Now starting to make a weighted and non-weighted summary plots of data in final file  ****'
    
    dataReadAndPlotNEW3(dataFinalFilename,finalFolder+'summaryPlot', weight=False, confLevels=False)
    dataReadAndPlotNEW3(dataFinalFilename,finalFolder+'summaryPlot-weighted', weight=True, confLevels=True)
    
    bestOrbit = bestOrbitFinderNEW(dataFinalFilename, printToScreen=False, saveToFile=False, returnAsList=True)
    orbitEllipsePlotFilename = finalFolder+"orbitEllipsePlot"
    orbitEllipsePlotter3(bestOrbit[0],bestOrbit[1],bestOrbit[3],bestOrbit[4],bestOrbit[5],bestOrbit[6],plotFilename=orbitEllipsePlotFilename,show=True)
    
    print '\n**** EVERYTHING FINISHED ****\n'
        
        
        
        
        
        