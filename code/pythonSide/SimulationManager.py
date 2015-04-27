#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import os
import shutil
import timeit
from multiprocessing import Process
import toolboxes as tools

class singleProcessStarter(Process):
    """
    This is the Manager object that controls the a single processes for a 
    SMODT simulation run.  It is called by the multiProcessStarter once for 
    each chain/process requested by the user through the simulation settings 
    file.
    
    :param str simSettingsFilename: File name, including the full path, to the 
        input settings file to use.
    :param str filename: File name, including the full path, for the output 
        data files.
    :param str cppCodeDir: Full path the the directory where there C++ code 
        lives.
    :param bool simAnneal: Perform Simulated Annealing instead of MCMC?
    :param bool loopedMCMC: Perform 'looped' MCMC instead of standard MCMC?
    """
    def __init__(self, simSettingsFilename, filename, cppCodeDir, simAnneal, mcONLY, loopedMCMC):
        
        Process.__init__(self)
        self.simSettingsFilename = simSettingsFilename
        self.filename = filename
        self.cppCodeDir = cppCodeDir    
        self.simAnneal = simAnneal  
        self.mcONLY = mcONLY
        self.loopedMCMC = loopedMCMC  
        
    def run(self):
        self.process(simSettingsFilename=self.simSettingsFilename, 
                     filename=self.filename, cppCodeDir=self.cppCodeDir, 
                     simAnneal=self.simAnneal, mcONLY=self.mcONLY ,loopedMCMC=self.loopedMCMC)
    
    def process(self, simSettingsFilename, filename, cppCodeDir, simAnneal, mcONLY, loopedMCMC):
        verbose = False
        if False:
            print 'Starting to run process for file title: '+filename
        
        #if simAnneal:
        #    CplusplusCodeCALL = os.path.join(cppCodeDir,'simAnnealOrbSimulator')
        #if loopedMCMC:
        #    CplusplusCodeCALL = os.path.join(cppCodeDir,'looped_MCMCorbSimulator')
        if mcONLY:
            CplusplusCodeCALL = os.path.join(cppCodeDir,'mcONLYorbSimulator')
        else:
            CplusplusCodeCALL = os.path.join(cppCodeDir,'MCMCorbSimulator')
        CplusplusCodeCALL = CplusplusCodeCALL+' '+simSettingsFilename+' '+filename
        if verbose:
            print 'CplusplusCodeCALL : ',CplusplusCodeCALL
            print "call being made"
        os.system(CplusplusCodeCALL)
        if verbose:
            print 'call completed'
                   
def multiProcessStarter(paramSettingsDict):
    """
    This will call singleProcessStarter to handle each individual chain.   
    The output dataFile for each chain will be returned as a list of filenames.
        
    :param dict paramSettingsDict: A dictionary created by using the 
        cFileToSimSettingsDict function that converts the SimSettings.txt into
        a Python dictionary.
    """
    # Get copied settings filename with its directory
    outSimSettingsFilename = os.path.join(paramSettingsDict['outputData_dir'],'code-used/'+os.path.basename(paramSettingsDict['UpdatedSettingsFile']))
    numChains = paramSettingsDict['numChains']
    master = []
    for processNumber in range(numChains):
        # create title for this thread's output data file
        ext = os.path.splitext(paramSettingsDict['outputData_filenameRoot'])[1]
        
        if ext=='':
            ext='.dat' # just a default value, '.txt' is also fine
        filenameUSE = 'outputData-chain_'+str(processNumber+1)+ext
        filenameUSE_full = os.path.join(paramSettingsDict['outputData_dir'],filenameUSE)
        # start process
        master.append(singleProcessStarter(simSettingsFilename=outSimSettingsFilename, 
                                         filename=filenameUSE_full, 
                                         cppCodeDir=paramSettingsDict['cppCodeDir'],
                                         simAnneal=paramSettingsDict['simAnneal'],
                                         mcONLY=paramSettingsDict['mcONLY'],
                                         loopedMCMC=paramSettingsDict['loopedMCMC']))
        master[processNumber].start()
    # wait for completion of all process $$$$ STILL DON'T UNDERSTAND HOW THIS WORKS $$$$
    for processNumber in range(numChains):
        master[processNumber].join()    
    
    dataFiles = []
    for processNumber in range(numChains):
        dataFiles.append(master[processNumber].filename)
   
    return dataFiles
    
def simulator(paramSettingsDict):
    """
    This is the master function to start a multiprocess SMODT simulator run.
    The singleProcessStarter will be called by multiProcessStarter to handle each individual chain
    and the toolboxes will then be utilized to perform post-completion 
    statistical calculations and plotting to summarize the results.  The user 
    can choose which calculations and plots to produce with the simulation 
    settings.    
    
        
    :param dict paramSettingsDict: A dictionary created by using the 
        cFileToSimSettingsDict function that converts the SimSettings.txt into
        a Python dictionary.
    """
    
    ##set up Python log file name and open file
    logFilename = os.path.join(paramSettingsDict['outputData_dir'],"processManagerLogFile.txt")
    PMlogFile = open(logFilename,'a')
    
    ##Start a process that tracks the memory usage in the background
    RT = tools.RAMtracker(paramSettingsDict)
    
    ###############################################################
    # make call to 'make' to build the C++ code for this simulation
    # to ensure it is good to go
    ###############################################################
    pwd = os.curdir
    makeDir = paramSettingsDict['cppCodeDir']
    s =  '*'*23+"  About to compile C++ code  " +'*'*23+"\n"+'*'*75+"\n"
    s = s+ 'makeDir:'+makeDir
    print s
    PMlogFile.write(s)
    os.chdir(makeDir)
    #if paramSettingsDict['simAnneal']:
    #    os.system('make simAnnealOrbSimulator')
    #elif paramSettingsDict['loopedMCMC']:
    #    os.system('make looped_MCMCorbSimulator')
    if paramSettingsDict['mcONLY']:
        os.system('make mcONLYorbSimulator')
    else:
        os.system('make MCMCorbSimulator')
    os.chdir(pwd)
    s =  '*'*75+"\n"+'*'*28+' C++ code compiled '+'*'*28+"\n"
    s = s+ '\nsimulator: *************  STARTED Multiprocess Simulation *************'
    print s
    PMlogFile.write(s)
    
    ############################################################
    ### start running threads ###
    ############################################################
    # record the time the Total Simulation started
    tic=timeit.default_timer()
    ## Use multiProcessStarter func to handle starting and running threads  
    s =  'simulator: *************        Entering C++ Side         *************' 
    print s
    PMlogFile.write(s)
    dataFiles = multiProcessStarter(paramSettingsDict)
    s =  'simulator: *************        Back From C++ Side        *************' 
    s+= '\nsimulator: ************* FINISHED Multiprocess Simulation *************' 
    # write total elapsed time to screen and log.
    toc=timeit.default_timer()
    totalTimeString2 = tools.gen.timeString(toc - tic)
    s= s+'\n\nTotal simulation took '+totalTimeString2+' to complete.\n'
    s+='\n'+"*"*50+'\n*********    STARTING Post-Processing    *********\n'+"*"*50+'\n'
    print s
    PMlogFile.write(s)
    RT.chainsDonePrint()
    ############################################################
    # combine the input files into one final file
    ############################################################
    dataFinalFilename = os.path.join(paramSettingsDict['outputData_dir'],'outputData-ALL.dat') 
    s= '**** Starting to write original data to final combined file ****'
    if paramSettingsDict['SILENT']==False:
        print s
    PMlogFile.write(s)
    tools.gen.dataFileCombiner(dataFiles, dataFinalFilename)
    
    ############################################################
    # Prepare system data dictionary for plotting functions
    ############################################################
    sysDatafilename = os.path.join(paramSettingsDict['outputData_dir'],'code-used/'+paramSettingsDict['SystemDataFilename'])
    sysDataDict = tools.gen.sysDataToDict(sysDatafilename)
    
    ############################################################
    # determine if to plot a 3x1 or 3x2 plot
    ############################################################
    plot4x1 = False
    if ((paramSettingsDict["longAN_degMAX"]==0)and(paramSettingsDict["periodMAX"]==0)and(paramSettingsDict["inclination_degMAX"]==0)):
        plot4x1 = True
    if paramSettingsDict['SILENT']==False:
        if plot4x1:
            s='plot4x1 found to be True, so only plotting 4 key varied params.'
        else:
            s='plot4x1 found to be False, so plotting all params.'
        
    if paramSettingsDict['TcStepping']==True:
        s=s+'\nTcStepping also found to be true, so making plots of Tc instead of To'
    else:
        s=s+'\nTcStepping also found to be false, so making plots of To instead of Tc'
    if paramSettingsDict['SILENT']==False:
        print s
    PMlogFile.write(s)         
    
    ################################################################
    # Call the function to find the best orbit values for reference.
    ################################################################
    if paramSettingsDict['SILENT']==False:
        print '#'*50
    p = True
    if paramSettingsDict['SILENT']:
        p = False

    bestOrbit = tools.gen.bestOrbitFinder(dataFinalFilename, printToScreen=p, saveToFile=True, returnAsList=True)
    if False:#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        bestOrbit = [70.0, 0.5, 2457000., 2457000., 5.0, 40.0, 50.0, 3.34718746581, 0, [0], 0]#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    logFilename = os.path.join(paramSettingsDict['outputData_dir'],'log-chain_1.txt')
    [nu,nuRV,nuDI,printStr] = tools.gen.findNuFromLog(logFilename)
    nus = [nu,nuRV,nuDI]
    if paramSettingsDict['SILENT']==False:
        print printStr+'\n'+'#'*50
    PMlogFile.write(printStr+'\n'+'#'*50+'\n')
    
    ############################################################
    ## make general parameter result summary figures
    ############################################################
    summaryPlotFile = os.path.join(paramSettingsDict['outputData_dir'],'summaryPlot')
    keyPosteriorsPlotFile = os.path.join(paramSettingsDict['outputData_dir'],'KeyPosteriorsPlot')
    cleanDataFilename=''
    if True:
        if (paramSettingsDict['simAnneal']==False)and(paramSettingsDict['makePosteriorsPlot']):
            s= '**** Starting to plot the posterior distributions of the final data file ****'
            print s
            PMlogFile.write(s)
            if True:
                tools.plot.summaryPlotter(dataFinalFilename, summaryPlotFile, weight=False, confLevels=True, nu=nu, plot4x1=plot4x1, TcStepping=paramSettingsDict['TcStepping'] )         
            if False:
                if paramSettingsDict['SILENT']==False:
                    print "\n\n"+"!"*75+'\nNOTE: Making Posteriors plot with the makeCleanSummaryPlot function instead of standard summaryPlotter\n'+"!"*75+"\n\n"
                cleanDataFilename = tools.plot.makeCleanSummaryPlot(dataFinalFilename)
            s = '**** Back from making summary plot if requested ****\n'
            if paramSettingsDict['SILENT']==False:
                print s
            PMlogFile.write(s)
    if True:
        tools.plot.stackedPosteriorsPlotterFunc([dataFinalFilename],keyPosteriorsPlotFile)
    
    ############################################################
    ## Make DI ellipse plot if DI data exists
    ############################################################
    DIdatafilename = os.path.join(paramSettingsDict['outputData_dir'],'code-used/'+paramSettingsDict['DIdataFilename'])
    chiSquaredStrDI=''
    if os.path.exists(DIdatafilename)and ((paramSettingsDict['RVonly']==False)and(paramSettingsDict['makeOrbitPlots'])):
        s = '**** Starting to make a DI orbit plot ****'
        print s
        PMlogFile.write(s)
        DIdataDict = tools.di.DIdataToDict(DIdatafilename)
        orbitEllipsePlotFilename = os.path.join(paramSettingsDict['outputData_dir'],'orbitEllipsePlot')
        #update argPeri value to take offset into account
        argPeriUse = bestOrbit[6]+paramSettingsDict['argPeriOffsetDI']
        chiSquaredStrDI = tools.plot.orbitEllipsePlotter(bestOrbit[0],bestOrbit[1],bestOrbit[4],bestOrbit[5],argPeriUse,bestOrbit[7],\
                             sysDataDict,DIdataDict,plotFilename=orbitEllipsePlotFilename,show=False,To=bestOrbit[2], nuDI=nuDI)          
        s = '****   Back from making a DI orbit plot   ****\n'
        if paramSettingsDict['SILENT']==False:
            print s
        PMlogFile.write(s)
        
    ############################################################
    ## Make RV scatter.trend plots if RV data exists
    ############################################################
    RVdatafilename = os.path.join(paramSettingsDict['outputData_dir'],'code-used/'+paramSettingsDict['RVdataFilename'])
    chiSquaredStrRV=''
    if True:
        s= "RVdatafilename = "+paramSettingsDict['RVdataFilename']
        if paramSettingsDict['SILENT']==False:
            print s
        PMlogFile.write(s)
    if os.path.exists(RVdatafilename) and (paramSettingsDict['DIonly']==False)and(paramSettingsDict['makeOrbitPlots']):
        s = '**** Starting to make a RV orbit plot ****'
        print s
        PMlogFile.write(s)
        RVdataDict = tools.rv.RVdataToDict(RVdatafilename)
        rvPlotFilename = os.path.join(paramSettingsDict['outputData_dir'],'orbitRVplot')
        if paramSettingsDict['loopedMCMC']==False:
            #update argPeri value to take offset into account
            argPeriUse = bestOrbit[6]+paramSettingsDict['argPeriOffsetRV']
            # full orbit
            chiSquaredStrRV = tools.plot.rvPlotter(bestOrbit[1],bestOrbit[2],bestOrbit[3],bestOrbit[4],bestOrbit[5],argPeriUse,bestOrbit[7], \
                  sysDataDict,RVdataDict,paramSettingsDict,K=bestOrbit[8],RVoffsets=bestOrbit[9],\
                  nuRV=nuRV,plotFilename=rvPlotFilename+'-FullOrbit', show=False, plotFullOrbit=True)              
            if False:
                # data range limited orbit
                tools.plot.rvPlotter(bestOrbit[1],bestOrbit[2],bestOrbit[3],bestOrbit[4],bestOrbit[5],argPeriUse,bestOrbit[7], \
                        sysDataDict,RVdataDict,paramSettingsDict,K=bestOrbit[8],RVoffsets=bestOrbit[9],\
                        nuRV=nuRV, plotFilename=rvPlotFilename, show=False, plotFullOrbit=False)                            
        s = '****   Back from making a RV orbit plot   ****'
        if paramSettingsDict['SILENT']==False:
            print s
        PMlogFile.write(s)
        
        #################################################################################################
        ##$$$$$$$$$$$$ This section is for the 'looped MCMC' method.  Not sure if we should #$$$$$$$$$$$$ 
        ##$$$$$$$$$$$$ clean it up and get it working again.  So, doesn't work for now!!!!! #$$$$$$$$$$$$
#        ## Handle 100 MOD dataset plots
#        if paramSettingsDict['loopedMCMC']:
#            
#            rvPlotter(bestOrbit[1],bestOrbit[2],bestOrbit[3],bestOrbit[4],bestOrbit[5],bestOrbit[6], \
#                  sysDataDict,RVdataDict,paramSettingsDict,K=bestOrbit[7],RVoffsets=bestOrbit[8],\
#                  plotFilename=rvPlotFilename+'-FullOrbit', show=False, plotFullOrbit=True)
#            
#            if False:
#                RVmodDatafilename = os.path.join(paramSettingsDict['settings_and_InputDataDir'],'mod1000Dataset.dat')
#                RVmodDataDict = RVdataToDict(RVmodDatafilename)
#                AllOutputsDict = outputDatafileToDict(dataFinalFilename)
#                RVmodDataDictUSE = RVmodDataDict
#                s= 'Number of datasets in modDataset file found = ',len(RVmodDataDict['RVs'])
#                s=s+ 'Number of epochs in modDataset file found = ',len(RVmodDataDict['RVs'][0])                                                     
#                s=s+ 'Making RV plots for '+str(len(AllOutputsDict['es']))+' MOD datasets'
#                print s
#                PMlogFile.write(s)
#                
#                for modDatasetNum in range(0,len(AllOutputsDict['es'])):
#                    RVmodDataDict = RVdataToDict(RVmodDatafilename)
#                    s= 'Making RV plot for MOD dataset # ',(modDatasetNum+1)
#                    print s
#                    PMlogFile.write(s)
#                    CURrvPlotFilename = os.path.join(paramSettingsDict['outputData_dir'],'RVplot_for_RVmodData_'+str(modDatasetNum+1))
#                    RVmodDataDictUSE["RV_epochs"] = RVmodDataDict['RV_epochs'][modDatasetNum]
#                    RVmodDataDictUSE["RV_errors"] = RVmodDataDict['RV_errors'][modDatasetNum]
#                    RVmodDataDictUSE["RVs"] = RVmodDataDict['RVs'][modDatasetNum]
#                    
#                    rvPlotterDuo(AllOutputsDict["es"][modDatasetNum],AllOutputsDict["Ts"][modDatasetNum],AllOutputsDict["periods"][modDatasetNum],
#                                 AllOutputsDict["inclinations"][modDatasetNum],AllOutputsDict["argPeris"][modDatasetNum],AllOutputsDict["a_totals"][modDatasetNum], \
#                                 sysDataDict,RVmodDataDictUSE,paramSettingsDict,RVoffsets=AllOutputsDict["RVoffsets"][modDatasetNum],\
#                                 plotFilename=CURrvPlotFilename+'-FullOrbit', show=False, plotFullOrbit=True)
        #################################################################################################
    
    ############################################################
    ## Make progress summary plot of parameters   
    ############################################################
    effectivePointsStr = ''
    burnInStr = ''
    # set up files and make plots for simAnneal data as well if MCMC is being ran
    if paramSettingsDict['mcONLY']==False:
        if paramSettingsDict['makeSimAnnealProgPlots'] or paramSettingsDict['makeMCMCprogPlots']:
            s= '**** Starting to make a parameter progress summary plots for each chain ****'
            print s
            PMlogFile.write(s)
        if paramSettingsDict['simAnneal']==False:
            MCMCsummaryRootPlotFilename = os.path.join(paramSettingsDict['outputData_dir'],'MCMCprogressSummary')   
            simAnnealSummaryRootPlotFilename = os.path.join(paramSettingsDict['outputData_dir'],'simAnnealProgressSummary')
            # make list of simAnneal data files for plotting
            simAnnealDataFiles = []
            for f in dataFiles:
                simAnnealFile = "SimAnneal_"+os.path.basename(f)
                simAnnealFileFull = os.path.join(paramSettingsDict['outputData_dir'],simAnnealFile)
                simAnnealDataFiles.append(simAnnealFileFull)
            if paramSettingsDict['makeSimAnnealProgPlots']:
                tools.plot.mcmcProgressPlotter(simAnnealDataFiles,simAnnealSummaryRootPlotFilename, nu=nu, plot4x1=plot4x1,TcStepping=paramSettingsDict['TcStepping'])
            if paramSettingsDict['makeMCMCprogPlots']:
                tools.plot.mcmcProgressPlotter(dataFiles,MCMCsummaryRootPlotFilename, nu=nu, plot4x1=plot4x1,TcStepping=paramSettingsDict['TcStepping'])
            #check burn-in lengths of MCMC part of simAnneal and proper MCMC chains
            if ((len(dataFiles)>1)and(paramSettingsDict["mcONLY"]==False))and(paramSettingsDict['CalcBurnIn']and(paramSettingsDict['simAnneal']==False)):
                if False:
                    if paramSettingsDict['SILENT']==False:
                        print '\nFor The Simulated Annealing files, the burn in values are:'
                    (s,burnInLengths) = tools.gen.burnInCalcMultiFile(simAnnealDataFiles,simAnneal=True)
                    burnInStr+=s 
                if paramSettingsDict['SILENT']==False:
                    print '\nFor The MCMC files, the burn in values are:'
                (s,burnInLengths) = tools.gen.burnInCalcMultiFile(dataFiles,simAnneal=False)
                burnInStr+=s 
            # calculate the correlation lengths of the MCMC part of the simAnneal and proper MCMC chains
            if paramSettingsDict['calcCorrLengths']:
                if False: 
                    print '***** Calculating the number of effective points for the MCMC part of simAnneal chains *****'
                    effectivePointsStr+= tools.gen.mcmcEffectivePointsCalc(simAnnealDataFiles,simAnneal=True)
                print '***** Calculating the number of effective points for the MCMC chains *****'
                effectivePointsStr+= tools.gen.mcmcEffectivePointsCalc(dataFiles,simAnneal=False)
            if (len(dataFiles)>1)and(paramSettingsDict['removeBurnIn'] and (paramSettingsDict["mcONLY"]==False))and(paramSettingsDict['CalcBurnIn']and (paramSettingsDict['simAnneal']==False)):
                #########################################################################
                ## make general parameter result summary figures AFTER BURN-IN STRIPPED!!
                #########################################################################
                ## strip burn ins and combine into new final file
                strippedNames = []
                for i in range(0,len(dataFiles)):
                    strippedName = dataFiles[i][:-4]+"_burnInStripped.dat"
                    tools.gen.burnInStripper(dataFiles[i], burnInLengths[i], strippedName)
                    strippedNames.append(strippedName)
                dataFinalFilename2 = os.path.join(paramSettingsDict['outputData_dir'],'outputData-ALL-burnInRemoved.dat') 
                tools.gen.dataFileCombiner(strippedNames, dataFinalFilename2)
                summaryPlotFile2 = os.path.join(paramSettingsDict['outputData_dir'],'summaryPlot-burnInRemoved')
                s= '**** Starting to plot the posterior distributions of data in final file AFTER BURN-IN REMOVED  ****'
                print s
                PMlogFile.write(s)
                cleanDataFilename=''
                if (paramSettingsDict['simAnneal']==False)and(paramSettingsDict['makePosteriorsPlot']):
                    if True:
                        tools.plot.summaryPlotter(dataFinalFilename2, summaryPlotFile2, weight=False, confLevels=True, nu=nu, plot4x1=plot4x1, TcStepping=paramSettingsDict['TcStepping'] )         
                    if False:
                        print "\n\n"+"!"*75+'\nNOTE: Making Posteriors plot with the makeCleanSummaryPlot function instead of standard summaryPlotter\n'+"!"*75+"\n\n"
                        cleanDataFilename = tools.plot.makeCleanSummaryPlot(dataFinalFilename2)
                s = '**** Back from making summary plot AFTER BURN-IN REMOVED if requested ***\n'
                if paramSettingsDict['SILENT']==False:
                    print s
                PMlogFile.write(s)
        
        else:
            if paramSettingsDict['makeSimAnnealProgPlots']:
                MCMCsummaryRootPlotFilename = os.path.join(paramSettingsDict['outputData_dir'],'simAnnealProgressSummary')
                tools.plot.mcmcProgressPlotter(dataFiles,MCMCsummaryRootPlotFilename, nu=nu, plot4x1=plot4x1,TcStepping=paramSettingsDict['TcStepping'])
        
    ############################################################
    # finish the Gelman-Rubin statistic calculations if requested
    ############################################################
    if (paramSettingsDict['CalcGelmanRubin']and paramSettingsDict['useMultiProcessing'])and((paramSettingsDict["mcONLY"]==False)and(paramSettingsDict['simAnneal']==False)):
        tools.gen.gelmanRubinStage2(dataFiles)    
    
    ####################################################################################################
    ## perform final wrap up functions, ie finish tracking RAM use and recored simulations total results
    ####################################################################################################    
    ##wrap up background process tacking RAM use and plot results
    maxRAMuse = RT.wrapUp()
    ## recordResults MUST be done before deleting all the extra files!!!!
    tools.gen.recordResults(paramSettingsDict,maxRAMuse,nus,chiSquaredStrDI,chiSquaredStrRV,effectivePointsStr,burnInStr)
    
    ############################################################
    ## Move useful results files into new results sub-folder
    ############################################################
    print '\n'+"*"*50+"\n************* Starting clean up steps ************\n"+"*"*50
    resultsFolder = os.path.join(paramSettingsDict['outputData_dir'],'RESULTS/')
    os.mkdir(resultsFolder)
    origFolder = paramSettingsDict['outputData_dir']
    ## build up list of files to copy over
    origFiles = []
    origFiles.append("bestOrbit.txt")
    # get Summary plot filenames
    if (paramSettingsDict['simAnneal']==False)and(paramSettingsDict['makePosteriorsPlot']):
        origFiles.append(os.path.basename(summaryPlotFile)+".png")
        origFiles.append(os.path.basename(summaryPlotFile)+"-ChiSquaredDist.png")
        origFiles.append(os.path.basename(keyPosteriorsPlotFile)+".png")
        if paramSettingsDict['DIonly']==False:
            origFiles.append(os.path.basename(summaryPlotFile)+"-RVoffsets.png")
        if (len(dataFiles)>1)and(paramSettingsDict['removeBurnIn'] and (paramSettingsDict["mcONLY"]==False))and(paramSettingsDict['CalcBurnIn']and (paramSettingsDict['simAnneal']==False)):
            origFiles.append(os.path.basename(summaryPlotFile)+"-burnInRemoved.png")
            origFiles.append(os.path.basename(summaryPlotFile)+"-burnInRemoved-ChiSquaredDist.png")
            if paramSettingsDict['DIonly']==False:
                origFiles.append(os.path.basename(summaryPlotFile)+"-burnInRemoved-RVoffsets.png")
    # get RV plot filenames
    if os.path.exists(RVdatafilename) and (paramSettingsDict['DIonly']==False)and(paramSettingsDict['makeOrbitPlots']):
        origFiles.append(os.path.basename(rvPlotFilename)+'-FullOrbit.png')
        origFiles.append(os.path.basename(rvPlotFilename)+'-FullOrbit-paramInfo.png')
        origFiles.append(os.path.basename(rvPlotFilename)+'-FullOrbit_TREND.png')
        origFiles.append(os.path.basename(rvPlotFilename)+'_DataDist.png')
        origFiles.append(os.path.basename(rvPlotFilename)+'-FullOrbit_DataDist.png')
        origFiles.append(os.path.basename(rvPlotFilename)+'_TREND.png')
        origFiles.append(os.path.basename(rvPlotFilename)+'.png')
        origFiles.append(os.path.basename(rvPlotFilename)+'-paramInfo.png')
    # get DI plot filenames
    if os.path.exists(DIdatafilename)and ((paramSettingsDict['RVonly']==False)and(paramSettingsDict['makeOrbitPlots'])):
        origFiles.append(os.path.basename(orbitEllipsePlotFilename)+'.png')
        origFiles.append(os.path.basename(orbitEllipsePlotFilename)+'-CROPPED.png')
        origFiles.append(os.path.basename(orbitEllipsePlotFilename)+'-paramInfo.png')
    # get progress plot filenames
    if paramSettingsDict['simAnneal']:
        if paramSettingsDict['makeSimAnnealProgPlots']:
            for chainNum in range(1,len(dataFiles)+1):
                origFiles.append(os.path.basename(MCMCsummaryRootPlotFilename)+"-chain_"+str(chainNum)+".png" )
                if paramSettingsDict['DIonly']==False:
                    origFiles.append(os.path.basename(MCMCsummaryRootPlotFilename)+"-chain_"+str(chainNum)+"-RVoffsets.png")
    elif paramSettingsDict['makeMCMCprogPlots']:
        for chainNum in range(1,len(dataFiles)+1):
            origFiles.append(os.path.basename(MCMCsummaryRootPlotFilename)+"-chain_"+str(chainNum)+".png")
            if paramSettingsDict['DIonly']==False:
                origFiles.append(os.path.basename(MCMCsummaryRootPlotFilename)+"-chain_"+str(chainNum)+"-RVoffsets.png")
        if paramSettingsDict['makeSimAnnealProgPlots']:
            for chainNum in range(1,len(dataFiles)+1):
                origFiles.append(os.path.basename(simAnnealSummaryRootPlotFilename)+"-chain_"+str(chainNum)+".png")
                if paramSettingsDict['DIonly']==False:
                    origFiles.append(os.path.basename(simAnnealSummaryRootPlotFilename)+"-chain_"+str(chainNum)+"-RVoffsets.png")
    # get GR filenames
    if (paramSettingsDict['CalcGelmanRubin']and paramSettingsDict['useMultiProcessing'])and((paramSettingsDict["mcONLY"]==False)and(paramSettingsDict['simAnneal']==False)):
        origFiles.append('GRvalues.txt')
        origFiles.append('Tvalues.txt')
    # get RESULTS.txt filenames
    origFiles.append('RESULTS.txt')
    origFiles.append('RAMusage_clean.png')
    ## move the files over 
    if paramSettingsDict['SILENT']==False:
        print 'Moving all useful results files into a new RESULTS folder:\n'+resultsFolder
    else:
        print "FINAL results moved to:\n"+resultsFolder
    for file in origFiles:
        try:
            shutil.move(os.path.join(origFolder,file),os.path.join(resultsFolder,file))
        except:
            nothing=True
            
    ## build up list of log files and move to new folder
    origFiles = []
    for chainNum in range(1,len(dataFiles)+1):
        origFiles.append("log-chain_"+str(chainNum)+".txt")
    origFiles.append('processManagerLogFile.txt')
    origFiles.append("RAMusage_clean.log")
    ## move the files over 
    logFolder = os.path.join(paramSettingsDict['outputData_dir'],'logs/')
    os.mkdir(logFolder)
    if paramSettingsDict['SILENT']==False:
        print 'Moving all log files into a new logs folder:\n '+logFolder
    for file in origFiles:
        try:
            shutil.move(os.path.join(origFolder,file),os.path.join(logFolder,file))
        except:
            nothing=True
        
    ############################################################
    ## COPY OUTPUT PLOTS TO DROPBOX FOLDER?
    ############################################################
    if paramSettingsDict['CopyToDrobox']:
        f = paramSettingsDict['outputData_dir']
        if f[-1]=='/':
            f = f[:-1]
        fs = f.split('/')
        DBdir = os.path.join('/mnt/Data1/Todai_Work/Dropbox/SMODT-outputCopies/',fs[-1])
        print 'Copying all files in the RESULTS folder over to DropBox folder:\n '+DBdir
        tools.gen.copytree(resultsFolder, DBdir)
            
    ############################################################
    ## Delete chain and/or combined data files?
    ############################################################
    ## check if the user wanted the individual chain data files deleted
    if paramSettingsDict['delChainsAfter']:
        if (paramSettingsDict['simAnneal']==False)and(paramSettingsDict["mcONLY"]==False):
                print '***** Deleting Simulated Annealing chain data files *****'
                if paramSettingsDict['SILENT']==False:
                    print "-"*40
                for filename in simAnnealDataFiles:
                    if paramSettingsDict['SILENT']==False:
                        print 'Deleting file: '+os.path.basename(filename)
                    os.remove(filename)
                print '************ Deleting MCMC chain data files *************'
                if paramSettingsDict['SILENT']==False:
                    print "-"*40
                for filename in dataFiles:
                    if paramSettingsDict['SILENT']==False:
                        print 'Deleting file: '+os.path.basename(filename)
                    os.remove(filename)     
                if (len(dataFiles)>1)and(paramSettingsDict['removeBurnIn'] and (paramSettingsDict["mcONLY"]==False))and paramSettingsDict['CalcBurnIn']:
                    print '***** Deleting Burn-In removed MCMC chain data files ****'
                    if paramSettingsDict['SILENT']==False:
                        print "-"*40
                    for filename in strippedNames:
                        if paramSettingsDict['SILENT']==False:
                            print 'Deleting file: '+os.path.basename(filename)
                        os.remove(filename)           
        else:
            print '********* Deleting final output chain data files *********'
            if paramSettingsDict['SILENT']==False:
                print "-"*40
            for filename in dataFiles:
                if paramSettingsDict['SILENT']==False:
                    print 'Deleting file: '+os.path.basename(filename)
                os.remove(filename)
            if (len(dataFiles)>1)and(paramSettingsDict['removeBurnIn'] and (paramSettingsDict["mcONLY"]==False))and(paramSettingsDict['CalcBurnIn']and (paramSettingsDict['simAnneal']==False)):
                print '***** Deleting Burn-In removed MCMC chain data files *****'
                if paramSettingsDict['SILENT']==False:
                    print "-"*40
                for filename in strippedNames:
                    if paramSettingsDict['SILENT']==False:
                        print 'Deleting file: '+os.path.basename(filename)
                    os.remove(filename)
            
    ## delete GR chain files if requested
    if (paramSettingsDict['CalcGelmanRubin'] and(paramSettingsDict['delGRchainFiles'] and (len(dataFiles)>1)))and((paramSettingsDict["mcONLY"]==False)and (paramSettingsDict['simAnneal']==False)):
        print '******* Deleting final output GR chain value files ******'
        if paramSettingsDict['SILENT']==False:
            print "-"*40
        for chainNum in range(1,len(dataFiles)+1):
            filename = os.path.join(paramSettingsDict['outputData_dir'],"gelmanRubin-chain_"+str(chainNum)+".txt")
            if paramSettingsDict['SILENT']==False:
                print 'Deleting file: '+os.path.basename(filename)
            os.remove(filename) 
    ## delete combined data files if requested
    if paramSettingsDict['delCombinedDataAfter']:
        print '************* Deleting combined data files **************'
        if paramSettingsDict['SILENT']==False:
            print "-"*40
        if paramSettingsDict['SILENT']==False:
            print 'Deleting file: '+dataFinalFilename
        os.remove(dataFinalFilename)
        if cleanDataFilename!='':
            if paramSettingsDict['SILENT']==False:
                print 'Deleting file: '+cleanDataFilename
            os.remove(cleanDataFilename)
        if ((paramSettingsDict['removeBurnIn'] and (paramSettingsDict["mcONLY"]==False))and(len(dataFiles)>1))and (paramSettingsDict['CalcBurnIn']and (paramSettingsDict['simAnneal']==False)):
            if paramSettingsDict['SILENT']==False:
                print 'Deleting file: '+dataFinalFilename2
            os.remove(dataFinalFilename2)
            
    s= '\n'+'*'*75+'\n'+25*'*'+' EVERYTHING FINISHED '+'*'*29+'\n'+'*'*75+"\n"
    print s
    PMlogFile.write(s)
    
    ############################################################
    # write total elapsed time to screen and log.
    ############################################################
    toc=timeit.default_timer()
    totalTimeString2 = tools.gen.timeString(toc - tic)
    s= 'simulator:  Total simulation + post processing took '+totalTimeString2+' to complete.\n'
    print s
    PMlogFile.write(s+'\n\nCLOSING PM LOG NOW!!\n\n')
    PMlogFile.close()
    
##############################################################
## end
##############################################################
     
        
        