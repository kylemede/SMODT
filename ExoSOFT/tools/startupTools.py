import sys
import shutil
import os
import copy
import numpy as np
import exoSOFTlogger
import generalTools as genTools
import warnings
warnings.simplefilter("error")


log = exoSOFTlogger.getLogger('main.startupTools',lvl=100,addFH=False) 

def startup(argv,rootDir,rePlot=False):
    """
    Perform the following vital start up steps (and return filled out and cleaned up settingsDict):
    -Figure out important directories
    -Copy settings files to temp directory to combine them into the master settings dict
    -Get master settings dict
    -Make output folder
    -Remake the SWIG tools?
    -Copy all ExoSOFT code into output dir for emergencies.    
    -Check range min and max values, including updates for 'lowecc' mode.
    -find which parameters will be varying.
    -Check if start startParams and startSigmas in dictionary make sense.
    NOTE: have rootDir handled with setup.py??
    """    
    ## Pull in settings filename prepend from command line args, if provided
    prepend = ''
    if len(argv)>1:
        try:
            prepend = argv[1]
        except:
            print '\nWarning: the settings file prepended feature is not working correctly !!\n'    
    ## Load up the required specific directory paths in dict
    settingsDict = genTools.loadSettingsDict(rootDir+'settings_and_inputData/'+prepend)
    settingsDict['ExoSOFTdir']=rootDir
    settingsDict['settingsDir']=os.path.join(settingsDict['ExoSOFTdir'],'settings_and_inputData/')
    settingsDict['prepend']=prepend
    ## Make a directory (folder) to place all the files from this simulation run
    settingsDict['finalFolder'] = os.path.join(settingsDict['outDir'],settingsDict['outRoot'])
    ##if not doing a re-post analysis with rePlot.py
    if rePlot==False:
        if os.path.exists(settingsDict['finalFolder']):
            if settingsDict['logLevel']<50: ## Handle this with a 'clob' bool in dict??? $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                print '\n'+'$'*50
                print 'WARNING!! the folder:\n"'+settingsDict['finalFolder']+'"\nALREADY EXISTS!'
                print 'You can overwrite the data in it, or exit this simulation.'
                YN = raw_input('OVERWRITE current folder (y/n):')
            else:
                YN = 'y'
            if (('y' in YN) or ('Y' in YN)):
                shutil.rmtree(settingsDict['finalFolder'])
                os.mkdir(settingsDict['finalFolder'])
                dbDir = os.path.join(settingsDict['dbFolder'],settingsDict['outRoot'])
                if os.path.exists(dbDir):
                    shutil.rmtree(dbDir)
            else: #elif (('n' in YN) or ('N' in YN)):
                sys.exit()
            if settingsDict['logLevel']<50:
                print '$'*50+'\n'
        else:
            os.mkdir(settingsDict['finalFolder'])
        if False:
            for key in settingsDict:
               print key+' = '+repr(settingsDict[key])
        ## run make for swig if requested
        if settingsDict['remake']:
            cwd = os.getcwd()
            log.debug("-"*45+" Starting to remake CPP/SWIG tools "+45*"-")
            os.chdir(os.path.join(settingsDict['ExoSOFTdir'],'tools/cppTools/'))
            os.system('make clean')
            os.system('make')
            os.chdir(cwd)
            log.debug("-"*45+" Done re-making CPP/SWIG tools "+45*"-")
        ## copy all of current code to output directory
        codeCopyDir = os.path.join(settingsDict['finalFolder'],'codeUsed')
        os.mkdir(codeCopyDir)
        log.debug('Copying all files in the RESULTS folder over to DropBox folder:\n '+codeCopyDir)
        genTools.copytree(settingsDict['ExoSOFTdir'], codeCopyDir)
    #########################################################################################
    ## Check parameter range settings make sense for data provided and mode of operation.   #
    ## Then load up a list of the parameters to vary during simulation.                     #
    ## Note: Originally this was done in simulator startup, but thought better to move here.#
    #########################################################################################
    filenameRoot = os.path.join(genTools.getSimpleDictVal(settingsDict,'settingsDir'),genTools.getSimpleDictVal(settingsDict,'prepend'))
    realData = genTools.loadRealData(filenameRoot,dataMode=genTools.getSimpleDictVal(settingsDict,'dataMode'))
    ##check there are matching number of RV datasets and provided min/max vals for offsets
    if np.min(realData[:,6])<1e6:
        numVmins=len(genTools.getSimpleDictVal(settingsDict,'vMINs'))
        if numVmins==0:
            numVmins=1
        if np.max(realData[:,7])!=(numVmins-1):
            log.error("THE NUMBER OF vMINs DOES NOT MATCH THE NUMBER OF RV DATASETS!!!\n"+\
                           "please check the vMINs/vMAXs arrays in the simple settings file\n"+\
                           "to make sure they have matching lengths to the number of RV datasets.")
    if genTools.getSimpleDictVal(settingsDict,'TMAX')==genTools.getSimpleDictVal(settingsDict,'TMIN')==-1:
        ## set T range to [earliest Epoch-max period,earliest epoch]
        settingsDict['TMAX']=np.min(realData[:,0])
        settingsDict['TMIN']=np.min(realData[:,0])-genTools.getSimpleDictVal(settingsDict,'PMAX')*const.daysPerYear
    ##load up range min,max and sigma arrayS
    rangeMaxs = [genTools.getSimpleDictVal(settingsDict,'mass1MAX'),\
               genTools.getSimpleDictVal(settingsDict,'mass2MAX'),\
               genTools.getSimpleDictVal(settingsDict,'paraMAX'),\
               genTools.getSimpleDictVal(settingsDict,'OmegaMAX'),\
               genTools.getSimpleDictVal(settingsDict,'eMAX'),\
               genTools.getSimpleDictVal(settingsDict,'TMAX'),\
               genTools.getSimpleDictVal(settingsDict,'TMAX'),\
               genTools.getSimpleDictVal(settingsDict,'PMAX'),\
               genTools.getSimpleDictVal(settingsDict,'incMAX'),\
               genTools.getSimpleDictVal(settingsDict,'omegaMAX'),\
               0,\
               0,\
               genTools.getSimpleDictVal(settingsDict,'KMAX')]
    rangeMins = [genTools.getSimpleDictVal(settingsDict,'mass1MIN'),\
               genTools.getSimpleDictVal(settingsDict,'mass2MIN'),\
               genTools.getSimpleDictVal(settingsDict,'paraMIN'),\
               genTools.getSimpleDictVal(settingsDict,'OmegaMIN'),\
               genTools.getSimpleDictVal(settingsDict,'eMIN'),\
               genTools.getSimpleDictVal(settingsDict,'TMIN'),\
               genTools.getSimpleDictVal(settingsDict,'TMIN'),\
               genTools.getSimpleDictVal(settingsDict,'PMIN'),\
               genTools.getSimpleDictVal(settingsDict,'incMIN'),\
               genTools.getSimpleDictVal(settingsDict,'omegaMIN'),\
               0,\
               0,\
               genTools.getSimpleDictVal(settingsDict,'KMIN')]
    ##start with uniform sigma values
    sigSize = genTools.getSimpleDictVal(settingsDict,'strtSig')
    sigmas = [sigSize,sigSize,sigSize,sigSize,sigSize,sigSize,sigSize,sigSize,sigSize,sigSize,0,0,sigSize]
    if len(genTools.getSimpleDictVal(settingsDict,'vMINs'))!=len(genTools.getSimpleDictVal(settingsDict,'vMAXs')):
        log.critical("THE NUMBER OF vMINs NOT EQUAL TO NUMBER OF vMAXs!!!")
    for i in range(0,len(genTools.getSimpleDictVal(settingsDict,'vMINs'))):
        sigmas.append(sigSize)
        rangeMins.append(genTools.getSimpleDictVal(settingsDict,'vMINs')[i])
        rangeMaxs.append(genTools.getSimpleDictVal(settingsDict,'vMAXs')[i])
    rangeMaxs = np.array(rangeMaxs)
    rangeMins = np.array(rangeMins)
    ##For lowEcc case, make Raw min/max vals for param drawing during MC mode
    rangeMaxsRaw = copy.deepcopy(rangeMaxs)
    rangeMinsRaw = copy.deepcopy(rangeMins)
    if genTools.getSimpleDictVal(settingsDict,'lowEcc'):
        ## run through the possible numbers for e and omega to find min/max for RAW versions
        fourMin=1e6
        fourMax=-1e6
        nineMin=1e6
        nineMax=-1e6
        for omeg in range(int(rangeMins[9]*10),int(rangeMaxs[9]*10),1):
            omega = float(omeg)/10.0
            for e in range(int(rangeMins[4]*100),int(rangeMaxs[4]*100),1):
                ecc = float(e)/100.0
                four = np.sqrt(ecc)*np.sin((np.pi/180.0)*omega)
                nine = np.sqrt(ecc)*np.cos((np.pi/180.0)*omega)
                if four>fourMax:
                    fourMax = four
                if four<fourMin:
                    fourMin = four
                if nine>nineMax:
                    nineMax = nine
                if nine<nineMin:
                    nineMin = nine
        rangeMaxsRaw[9] = nineMax
        rangeMaxsRaw[4] = fourMax
        rangeMinsRaw[9] = nineMin
        rangeMinsRaw[4] = fourMin
    ## figure out which parameters are varying in this run.
    ## Don't vary atot or chiSquared ever, and take care of TcEqualT and Kdirect cases
    paramInts = []
    for i in range(0,len(rangeMins)):
        if (i!=10)and(i!=11):
            if (i>12):
                if genTools.getSimpleDictVal(settingsDict,'dataMode')!='DI':
                    if rangeMaxs[i]!=0:
                        paramInts.append(i) 
            elif (i==8)or(i==12):
                if (genTools.getSimpleDictVal(settingsDict,'dataMode')!='RV')or(genTools.getSimpleDictVal(settingsDict,'Kdirect')==False):
                    if (rangeMaxs[8]!=0)and(i==8):
                        paramInts.append(8)
                elif genTools.getSimpleDictVal(settingsDict,'Kdirect'):
                    if (rangeMaxs[12]!=0)and(i==12):
                        paramInts.append(12)                                           
            elif (i==2)or(i==3)or(i==0)or(i==1):
                if (genTools.getSimpleDictVal(settingsDict,'dataMode')!='RV'):
                    if(rangeMaxs[i]!=0):
                        paramInts.append(i)
                elif (genTools.getSimpleDictVal(settingsDict,'Kdirect')==False)and((i!=3)and(i!=2)):
                    if(rangeMaxs[i]!=0):
                        paramInts.append(i)
            elif rangeMaxs[i]!=0:
                if (i==5)or(i==6):
                    if genTools.getSimpleDictVal(settingsDict,'TcStep')and(i!=5):
                            paramInts.append(i)
                    elif (genTools.getSimpleDictVal(settingsDict,'TcStep')==False)and(i!=6):
                            paramInts.append(i)
                else:
                    paramInts.append(i)
    
        
    ## push all these important parameter related items into the dict for later use.
    settingsDict['realData'] = realData
    settingsDict['rangeMinsRaw'] = rangeMinsRaw
    settingsDict['rangeMaxsRaw'] = rangeMaxsRaw
    settingsDict['rangeMins'] = rangeMins
    settingsDict['rangeMaxs'] = rangeMaxs
    settingsDict['paramInts'] = np.array(paramInts)
    
    settingsDict = modePrep(settingsDict)
    
    return settingsDict

def modePrep(settingsDict):
    """
    Check if start startParams and startSigmas in dictionary make sense.
    Arrays must exist, be the right length, and have non-zeros values for all varying parameters.     
    """

    startParams = genTools.getSimpleDictVal(settingsDict,'startParams')
    startSigmas = genTools.getSimpleDictVal(settingsDict,'startSigmas')
    paramInts = genTools.getSimpleDictVal(settingsDict,'paramInts')
    rangeMaxs = genTools.getSimpleDictVal(settingsDict,'rangeMaxs')
    passed = False
    if (type(startParams)==list)or(type(startParams)==np.ndarray):
        if type(startParams)==list:
            startParams = np.array(startParams)
        if len(startParams)==len(rangeMaxs):
            i=0
            passed = True
            while (i<len(startParams))and(passed==True):
                if i in paramInts:
                    if startParams[i]==0:
                        passed=False
                i+=1
        else:
            passed=False
    if passed==False:
        startParams = False
        log.info("Original startParams in settings files were not usable, so setting to False.")
    if (startParams==False)and(settingsDict['symMode'][0] in ['ST','MCMC']):
        log.critical('ST or MCMC mode requested, but not starting parameters provided.  Quiting ExoSOFT!!')
        #***************************************************************************************************
        sys.exit("MUST PROVIDE USEFUL STARTPARAMS IN SIMPLE SETTINGS DICT FOR ST or MCMC MODE, ELSE NOTHING TO START CHAINS WITH!\n\n!!EXITING ExoSOFT!!")
        #***************************************************************************************************
    passed = False
    if (type(startSigmas)==list)or(type(startSigmas)==np.ndarray):
        if len(startSigmas)==len(rangeMaxs):
            i=0
            passed = True
            while (i<len(startSigmas))and(passed==True):
                if i in paramInts:
                    if startSigmas[i]==0:
                        passed=False
                i+=1
        else:
            passed=False
    if passed==False:
        log.info("Original startSigmas in settings files were not usable, so setting to default values.")
        startSigmas = sigmas
    # clean up sigs of pars that were not varying, and that ary is a ndarray.
    for i in range(0,len(startSigmas)):
        if i not in paramInts:
            startSigmas[i] = 0
    if type(startSigmas)!=np.ndarray:
        startSigmas = np.array(startSigmas)
        
    ##make list of stages to run
    stgLstDict = {'MC':['MC'],'SA':['SA'],'SAST':['SA','ST'],'ST':['ST'],'SASTMCMC':['SA','ST','MCMC'],'MCMC':['MCMC']}
    stageList = stgLstDict[genTools.getSimpleDictVal(settingsDict,'symMode')]
        
    settingsDict['startParams'] = startParams
    settingsDict['startSigmas'] = startSigmas
    settingsDict['stageList'] = (stageList, 'List of stages to run')
    
    return settingsDict

#END OF FILE