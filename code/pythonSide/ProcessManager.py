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
    
        if False:
            print 'Starting to run process for file title: '+filename
        
        if simAnneal:
            CplusplusCodeCALL = os.path.join(cppCodeDir,'simAnnealOrbSimulator')
        elif loopedMCMC:
            CplusplusCodeCALL = os.path.join(cppCodeDir,'looped_MCMCorbSimulator')
        elif mcONLY:
            CplusplusCodeCALL = os.path.join(cppCodeDir,'mcONLYorbSimulator')
        else:
            CplusplusCodeCALL = os.path.join(cppCodeDir,'MCMCorbSimulator')
        CplusplusCodeCALL = CplusplusCodeCALL+' '+simSettingsFilename+' '+filename
        if False:
            print 'CplusplusCodeCALL : ',CplusplusCodeCALL
            print "call being made"
        os.system(CplusplusCodeCALL)
        if False:
            print 'call completed'
        
def multiProcessStarter(paramSettingsDict):
    """
    This is the master function to start a multiprocess SMODT simulator run.
    The singleProcessStarter will be called to handle each individual chain
    and the toolboxes will then be utilized to perform post-completion 
    statistical calculations and plotting to summarize the results.  The user 
    can choose which calculations and plots to produce with the simulation 
    settings.    
    
        
    :param dict paramSettingsDict: A dictionary created by using the 
        cFileToSimSettingsDict function that converts the SimSettings.txt into
        a Python dictionary.
    """
    #set up Python log file name and open file
    logFilename = os.path.join(paramSettingsDict['outputData_dir'],"processManagerLogFile.txt")
    PMlogFile = open(logFilename,'a')
    
    ###############################################################
    # make call to 'make' to build the C++ code for this simulation
    # to ensure it is good to go
    ###############################################################
    pwd = os.curdir
    makeDir = paramSettingsDict['cppCodeDir']
    s =  "\n"+'*'*95+"\n"+'*'*33+"  About to compile C++ code  " +'*'*33+"\n"+'*'*95+"\n"
    s = s+ 'makeDir:'+makeDir
    print s
    PMlogFile.write(s)
    os.chdir(makeDir)
    if paramSettingsDict['simAnneal']:
        os.system('make simAnnealOrbSimulator')
    elif paramSettingsDict['mcONLY']:
        os.system('make mcONLYorbSimulator')
    elif paramSettingsDict['loopedMCMC']:
        os.system('make looped_MCMCorbSimulator')
    else:
        os.system('make MCMCorbSimulator')
    os.chdir(pwd)
    s =  "\n"+'*'*95+"\n"+'*'*38+' C++ code compiled '+'*'*38+"\n"+'*'*95+'\n'
    s = s+ '\nMultiprocess: $$$$$$$$$$$$$ STARTED Multiprocess Sim $$$$$$$$$$$$$$$\n'
    print s
    PMlogFile.write(s)
    
    # Get copied settings filename with its directory
    outSimSettingsFilename = os.path.join(paramSettingsDict['outputData_dir'],'code-used/'+os.path.basename(paramSettingsDict['UpdatedSettingsFile']))
    
    ############################################################
    ### start running threads ###
    ############################################################
    # record the time the Total Simulation started
    tic=timeit.default_timer()
    numProcesses = paramSettingsDict['numProcesses']
    master = []
    for processNumber in range(numProcesses):
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
    for processNumber in range(numProcesses):
        master[processNumber].join()    
        
    s= '\nMultiprocess: $$$$$$$$$$$$$ FINISHED Multiprocess Sim $$$$$$$$$$$$$$$\n' 
    # write total elapsed time to screen and log.
    toc=timeit.default_timer()
    totalTimeString2 = tools.gen.timeString(toc - tic)
    s= s+'\n\nTotal simulation took '+totalTimeString2+' to complete.\n'
    print s
    PMlogFile.write(s)
    
    ############################################################
    # combine the input files into one final file
    ############################################################
    dataFiles = []
    dataFinalFilename = os.path.join(paramSettingsDict['outputData_dir'],'outputData-ALL.dat') 
    s= '\nStarting to write original data to final combined file'
    print s
    PMlogFile.write(s)
    for processNumber in range(numProcesses):
        dataFiles.append(master[processNumber].filename)
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
    if plot4x1:
        s='plot4x1 found to be True, so only plotting 4 key varied params.'
    else:
        s='plot4x1 found to be False, so plotting all params.'
        
    if paramSettingsDict['TcStepping']==True:
        s=s+'\nTcStepping also found to be true, so making plots of Tc instead of To'
    else:
        s=s+'\nTcStepping also found to be false, so making plots of To instead of Tc'
    print s
    PMlogFile.write(s)      

    s= '\n**** Now combining all files into one final file ****\n'
    print s
    PMlogFile.write(s)   
    
    ################################################################
    # Call the function to find the best orbit values for reference.
    ################################################################
    print '#'*50
    bestOrbit = tools.gen.bestOrbitFinder(dataFinalFilename, printToScreen=True, saveToFile=True, returnAsList=True)
    logFilename = os.path.join(paramSettingsDict['outputData_dir'],'log-chain_1.txt')
    [nu,nuRV,nuDI,printStr] = tools.gen.findNuFromLog(logFilename)
    print printStr+'\n'+'#'*50
    PMlogFile.write(printStr+'\n'+'#'*50+'\n')
    
    ############################################################
    ## make general parameter result summary figures
    ############################################################
    summaryPlotFile = os.path.join(paramSettingsDict['outputData_dir'],'summaryPlot')
    s= '\n**** Now starting to make a non-weighted, conf Levels summary plots of data in final file if requested  ****'
    print s
    PMlogFile.write(s)
    cleanDataFilename=''
    if (paramSettingsDict['simAnneal']==False)and(paramSettingsDict['makePosteriorsPlot']):
        if True:
            tools.plot.summaryPlotter(dataFinalFilename, summaryPlotFile, weight=False, confLevels=True, nu=nu, plot4x1=plot4x1, TcStepping=paramSettingsDict['TcStepping'] )         
        if False:
            print "\n\n"+"!"*75+'\nNOTE: Making Posteriors plot with the makeCleanSummaryPlot function instead of standard summaryPlotter2\n'+"!"*75+"\n\n"
            cleanDataFilename = tools.plot.makeCleanSummaryPlot(dataFinalFilename)
    s = '\n**** Back from making summary plot if requested ***\n'
    print s
    PMlogFile.write(s)
    
    ############################################################
    ## Make DI ellipse plot if DI data exists
    ############################################################
    DIdatafilename = os.path.join(paramSettingsDict['outputData_dir'],'code-used/'+paramSettingsDict['DIdataFilename'])
    if os.path.exists(DIdatafilename)and ((paramSettingsDict['RVonly']==False)and(paramSettingsDict['makeOrbitPlots'])):
        s = '\n**** Now starting to make a DI orbit plot ***\n'
        print s
        PMlogFile.write(s)
        DIdataDict = tools.di.DIdataToDict(DIdatafilename)
        orbitEllipsePlotFilename = os.path.join(paramSettingsDict['outputData_dir'],'orbitEllipsePlot')
        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        #$$$$$$$$$$$$  NOTE: We were not sure if there was a 180deg shift  $$$$$$$$$$$$
        #$$$$$$$$$$$$        in argPeri between DI and RV models...        $$$$$$$$$$$$
        if False:#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            argPeriUse = bestOrbit[6]+180.0#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        else:#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            argPeriUse = bestOrbit[6]#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        tools.plot.orbitEllipsePlotter(bestOrbit[0],bestOrbit[1],bestOrbit[4],bestOrbit[5],argPeriUse,bestOrbit[7],\
                             sysDataDict,DIdataDict,plotFilename=orbitEllipsePlotFilename,show=False,To=bestOrbit[2], nuDI=nuDI)          
        s = '\n**** Back from making a DI orbit plot ***\n'
        print s
        PMlogFile.write(s)
        
    ############################################################
    ## Make RV scatter.trend plots if RV data exists
    ############################################################
    RVdatafilename = os.path.join(paramSettingsDict['outputData_dir'],'code-used/'+paramSettingsDict['RVdataFilename'])
    if True:
        s= "RVdatafilename = "+paramSettingsDict['RVdataFilename']
        print s
        PMlogFile.write(s)
    if os.path.exists(RVdatafilename) and (paramSettingsDict['DIonly']==False)and(paramSettingsDict['makeOrbitPlots']):
        s = '\n**** Now starting to make a RV orbit plot ***\n'
        print s
        PMlogFile.write(s)
        RVdataDict = tools.rv.RVdataToDict(RVdatafilename)
        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        #$$$$$$$$$$$$  NOTE: We were not sure if there was a 180deg shift  $$$$$$$$$$$$
        #$$$$$$$$$$$$        in argPeri between DI and RV models...        $$$$$$$$$$$$
        if True:#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            argPeriUse = bestOrbit[6]+180.0#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        else:#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            argPeriUse = bestOrbit[6]#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        rvPlotFilename = os.path.join(paramSettingsDict['outputData_dir'],'orbitRVplot')
        if paramSettingsDict['loopedMCMC']==False:
            # full orbit
            tools.plot.rvPlotter(bestOrbit[1],bestOrbit[2],bestOrbit[3],bestOrbit[4],bestOrbit[5],argPeriUse,bestOrbit[7], \
                  sysDataDict,RVdataDict,paramSettingsDict,K=bestOrbit[8],RVoffsets=bestOrbit[9],\
                  nuRV=nuRV,plotFilename=rvPlotFilename+'-FullOrbit', show=False, plotFullOrbit=True)              
            if True:
                # data range limited orbit
                tools.plot.rvPlotter(bestOrbit[1],bestOrbit[2],bestOrbit[3],bestOrbit[4],bestOrbit[5],argPeriUse,bestOrbit[7], \
                        sysDataDict,RVdataDict,paramSettingsDict,K=bestOrbit[8],RVoffsets=bestOrbit[9],\
                        nuRV=nuRV, plotFilename=rvPlotFilename, show=False, plotFullOrbit=False)                            
        s = '\n**** Back from making a RV orbit plot ***\n'
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
    s= '\n**** Now starting to make a parameter progress summary plots for each chain (if requested) ****'
    print s
    PMlogFile.write(s)
    # set up files and make plots for simAnneal data as well if MCMC is being ran
    if paramSettingsDict['mcONLY']==False:
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
            if paramSettingsDict['CalcBurnIn']:
                print '\nFor The Simulated Annealing files, the burn in values are:'
                tools.gen.burnInCalcMultiFile(simAnnealDataFiles,simAnneal=True)
                print '\nFor The MCMC files, the burn in values are:'
                tools.gen.burnInCalcMultiFile(dataFiles,simAnneal=False)
            # calculate the correlation lengths of the MCMC part of the simAnneal and proper MCMC chains
            if paramSettingsDict['calcCorrLengths']and False:
                print '\nCalculating the number of effective points for the MCMC part of simAnneal chains\n'
                tools.gen.MCMCeffectivePointsCalc(simAnnealDataFiles,simAnneal=True)
                print '\nCalculating the number of effective points for the MCMC chains\n'
                tools.gen.MCMCeffectivePointsCalc(dataFiles,simAnneal=False)
        else:
            if paramSettingsDict['makeSimAnnealProgPlots']:
                MCMCsummaryRootPlotFilename = os.path.join(paramSettingsDict['outputData_dir'],'simAnnealProgressSummary')
                tools.plot.mcmcProgressPlotter(dataFiles,MCMCsummaryRootPlotFilename, nu=nu, plot4x1=plot4x1,TcStepping=paramSettingsDict['TcStepping'])
        
    ############################################################
    # finish the Gelman-Rubin statistic calculations if requested
    ############################################################
    if paramSettingsDict['CalcGelmanRubin']and(paramSettingsDict['simAnneal']==False):
        tools.gen.gelmanRubinStage2(dataFiles)    
    
    ############################################################
    ## Delete chain and/or combined data files?
    ############################################################
    ## check if the user wanted the individual chain data files deleted
    if paramSettingsDict['delChainsAfter']:
        if (paramSettingsDict['simAnneal']==False)and(paramSettingsDict["mcONLY"]==False):
                print '\n\nDeleting Simulated Annealing chain data files\n'+"-"*40
                for filename in simAnnealDataFiles:
                    print 'Deleting file: '+os.path.basename(filename)
                    os.remove(filename)
                print '\nDeleting MCMC chain data files\n'+"-"*40
                for filename in dataFiles:
                    print 'Deleting file: '+os.path.basename(filename)
                    os.remove(filename)
        else:
            print '\n\nDeleting final output chain data files\n'+"-"*40
            for filename in dataFiles:
                print 'Deleting file: '+os.path.basename(filename)
                os.remove(filename)
    ## delete combined data files if requested
    if paramSettingsDict['delCombinedDataAfter']:
        print '\nDeleting combined data files\n'+"-"*40
        print 'Deleting file: '+dataFinalFilename
        os.remove(dataFinalFilename)
        if cleanDataFilename!='':
            print 'Deleting file: '+cleanDataFilename
            os.remove(cleanDataFilename)    
    
    s= '\n**** EVERYTHING FINISHED ****\n'
    print s
    PMlogFile.write(s)
    
    ############################################################
    # write total elapsed time to screen and log.
    ############################################################
    toc=timeit.default_timer()
    totalTimeString2 = tools.gen.timeString(toc - tic)
    s= '\n\nMultiprocess:  Total simulation + post processing took '+totalTimeString2+' to complete.\n'
    print s
    PMlogFile.write(s+'\n\nCLOSING PM LOG NOW!!\n\n')
    PMlogFile.close()
    
    ############################################################
    ## COPY OUTPUT PLOTS TO DROPBOX FOLDER?
    ############################################################
    if paramSettingsDict['CopyToDrobox']:
        f = paramSettingsDict['outputData_dir']
        if f[-1]=='/':
            f = f[:-1]
        fs = f.split('/')
        DBdir = os.path.join('/mnt/Data1/Todai_Work/Dropbox/SMODT-outputCopies/',fs[-1])
        origFiles = []
        copyFiles = []
        
        # get Summary plot filenames
        if (paramSettingsDict['simAnneal']==False)and(paramSettingsDict['makePosteriorsPlot']):
            origFiles.append(summaryPlotFile+".png")
            copyFiles.append(os.path.join(DBdir,os.path.basename(summaryPlotFile)+".png"))
            origFiles.append(summaryPlotFile+"-ChiSquaredDist.png")
            copyFiles.append(os.path.join(DBdir,os.path.basename(summaryPlotFile)+"-ChiSquaredDist.png"))
        # get RV plot filenames
        if os.path.exists(RVdatafilename) and (paramSettingsDict['DIonly']==False)and(paramSettingsDict['makeOrbitPlots']):
            origFiles.append(rvPlotFilename+'-FullOrbit.png')
            copyFiles.append(os.path.join(DBdir,os.path.basename(rvPlotFilename+'-FullOrbit.png')))
            origFiles.append(rvPlotFilename+'_DataDist.png')
            copyFiles.append(os.path.join(DBdir,os.path.basename(rvPlotFilename+'_DataDist.png')))
            origFiles.append(rvPlotFilename+'_TREND.png')
            copyFiles.append(os.path.join(DBdir,os.path.basename(rvPlotFilename+'_TREND.png')))
            origFiles.append(rvPlotFilename+'.png')
            copyFiles.append(os.path.join(DBdir,os.path.basename(rvPlotFilename+'.png')))
        # get DI plot filenames
        if os.path.exists(DIdatafilename)and ((paramSettingsDict['RVonly']==False)and(paramSettingsDict['makeOrbitPlots'])):
            origFiles.append(orbitEllipsePlotFilename+'.png')
            copyFiles.append(os.path.join(DBdir,os.path.basename(orbitEllipsePlotFilename+'.png')))
        # get progress plot filenames
        if paramSettingsDict['simAnneal']:
            if paramSettingsDict['makeSimAnnealProgPlots']:
                MCMCsummaryRootPlotFilename2 = os.path.join(DBdir,'simAnnealProgressSummary')
                for chainNum in range(1,len(dataFiles)+1):
                    origFiles.append(MCMCsummaryRootPlotFilename+"-chain_"+str(chainNum)+".png" )
                    copyFiles.append(MCMCsummaryRootPlotFilename2+"-chain_"+str(chainNum)+".png" )
        if paramSettingsDict['makeMCMCprogPlots']:
            MCMCsummaryRootPlotFilename2 = os.path.join(DBdir,'MCMCprogressSummary')
            for chainNum in range(1,len(dataFiles)+1):
                origFiles.append(MCMCsummaryRootPlotFilename+"-chain_"+str(chainNum)+".png")
                copyFiles.append(MCMCsummaryRootPlotFilename2+"-chain_"+str(chainNum)+".png")
            if paramSettingsDict['makeSimAnnealProgPlots']:
                simAnnealSummaryRootPlotFilename2 = os.path.join(DBdir,'simAnnealProgressSummary')
                for chainNum in range(1,len(dataFiles)+1):
                    origFiles.append(simAnnealSummaryRootPlotFilename+"-chain_"+str(chainNum)+".png")
                    copyFiles.append(simAnnealSummaryRootPlotFilename2+"-chain_"+str(chainNum)+".png")
        # get log filenames
        chainLogRootName = os.path.join(paramSettingsDict['outputData_dir'],'log')
        chainLogRootName2 = os.path.join(DBdir,'log')
        for chainNum in range(1,len(dataFiles)+1):
            origFiles.append(chainLogRootName+"-chain_"+str(chainNum)+".txt")
            copyFiles.append(chainLogRootName2+"-chain_"+str(chainNum)+".txt")
        origFiles.append(os.path.join(paramSettingsDict['outputData_dir'],'processManagerLogFile.txt'))
        copyFiles.append(os.path.join(DBdir,'processManagerLogFile.txt'))
        # get GR filenames
        origFiles.append(os.path.join(paramSettingsDict['outputData_dir'],'GRvalues.txt'))
        copyFiles.append(os.path.join(DBdir,'GRvalues.txt'))
        
        ## copy all files in lists
        for fileNum in range(0,len(origFiles)):
            if os.path.exists(origFiles[fileNum]):
                shutil.copy(origFiles[fileNum],copyFiles[fileNum])
            
        
   
##############################################################
## end
##############################################################
     
        
        