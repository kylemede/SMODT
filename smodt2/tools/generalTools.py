#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
#import numpy as np
import smodtLogger
import cppTools
import os
import shutil
import timeit
import glob
import numpy as np
import sys
import pyfits
import warnings
warnings.simplefilter("error")

log = smodtLogger.getLogger('main.genTools',lvl=100,addFH=False)

def test():
    log.info("inside the tools test func")
    
def corrLengthCalcVar(paramIN):
    """
    This version uses np.var
    This is the most ideal way to calculate the correlation length, instead of the std based version.
    
    This function will calculate the mean correlation length and return its value.
    This is equal to the number of steps it takes for the variance to equal half of the total chain's variance.
    This is done in a loop, calculating it in an end-to-end fashion with the result being the mean of those 
    correlation lengths.  
    
    :param paramIN:     parameter array after burn in data stripped
    :type paramIN:      array (list) of doubles
    """
    verbose = False
    if verbose:
        print "Entered corrLengthCalcVar"
    try:
        varALL = np.var(paramIN)
        if verbose:
            print 'varALL = '+str(varALL)
    except:
        useless=0
    #print 'varALL = '+str(varALL)
    halfVarALL = varALL/2.0
    CorrLength = meanCorrLength = len(paramIN)
    varCur=0
    
    if paramIN[0]==paramIN[-1]:
        if verbose:
            print 'First and last parameters were the same, so returning a length of the input array.'
        CorrLength=meanCorrLength=len(paramIN)
    else:
        startLoc = 0
        corrLengths = []
        notFinished=True
        while notFinished:
            for i in range(startLoc+1,len(paramIN)+1):
                if i>=len(paramIN):
                    #hit the end, so stop
                    notFinished=False
                    break
                try:
                    varCur = np.var(paramIN[startLoc:i])
                    if verbose:
                        print 'varCur = '+str(varCur)
                except:
                    useless=1
                if varCur>halfVarALL:
                    CorrLength = i-startLoc
                    corrLengths.append(CorrLength)
                    if verbose:
                        print 'CorrLength = '+str(CorrLength)
                    startLoc = i
                    break
        if (startLoc==0)and(CorrLength == len(paramIN)):
            print "PROBLEM: Param had a correlation length equal to param length, ie. the chain never burned in"
        meanCorrLength = int(np.mean(corrLengths))
    if verbose:
        print 'mean Correlation length found to be = ',meanCorrLength
        print "Leaving corrLengthCalcVar"
    return meanCorrLength    
    
def mcmcEffPtsCalc(outputDataFilename):
    """
    Calculate correlation length and the number of effective steps for each parameter 
    that was varying during the simulation.  The results are put into the log.
    """
    log.info("Starting to calculate correlation lengths")
    (head,data) = loadFits(outputDataFilename)
    numSteps = data.shape[0]
    (paramList,paramStrs,paramFileStrs) = getParStrs(head,latex=False)
    ### Make post tools cpp obj
    #PostCTools = cppTools.PostCtools()
    #dataC = []
    #for i in range(0,numSteps):
    #    dataC.append(data[i,:])
    #dataC=np.array(dataC,dtype=np.dtype('d'),order='C')
    #PostCTools.loadParamData(dataC)
    
    completeStr = '\n'+'-'*45+'\nThe correlation lengths of all params are:\n'+'-'*45+'\nparam #, param name, correlation length'
    completeStr+= ' -> total number of steps/correlation length = number of effective points\n'
    for i in range(0,len(paramList)):
        #print 'starting to calculate corr length for '+paramStrs[paramList[i]]+' with CPP'
        #tic=timeit.default_timer()
        #meanCorrLengthC = PostCTools.corrLenCalc(paramList[i])
        #print 'meanCorrLengthC = '+str(meanCorrLengthC)
        #print 'it took: '+timeStrMaker(timeit.default_timer()-tic)
        log.debug('starting to calculate corr length for '+paramStrs[paramList[i]]+' with Python')
        tic=timeit.default_timer()
        meanCorrLength = corrLengthCalcVar(data[:,paramList[i]])
        log.debug('it took: '+timeStrMaker(timeit.default_timer()-tic))
        currParamStr = str(paramList[i])+', '+paramStrs[paramList[i]]+", "+str(meanCorrLength)
        currParamStr+=    ' -> '+str(numSteps)+'/'+str(meanCorrLength)+' = '+str(numSteps/meanCorrLength)+'\n'
        completeStr+=currParamStr
        log.debug(currParamStr)
    log.debug(completeStr)
    return completeStr

def burnInCalc(mcmcFnames,combinedFname):
    """
    NOTE: SMODT was designed to start the full MCMC chain from the last point of the 
        Sigma Tuning stage.  As this stage effectively acts as a form of burn-in period
        the burn-in value found from the pure MCMC tends to be very short.
 
    Calculate the burn in for a set of MCMC chains following the formulation of Tegmark.
     
    Burn-in is defined as the first point in a chain where the likelihood is greater than 
    the median value of all the chains.  Thus, there MUST be more than 1 chain to perform this calculation.
    """
    log.info("Starting to calculate burn-in.")
     
    chiSquaredsALL = np.array([])
    burnInLengths = []
    # calculate median of combined data ary
    (head0,data0) = loadFits(combinedFname)
    nu = float(head0['nu'])
    chiSqsALL = data0[:,11]
    if type(chiSqsALL)!=np.ndarray:
        chiSqsALL = np.array(chiSqsALL)
    likelihoodsALL = np.exp(-nu*chiSqsALL/2.0)
    medainALL = np.median(likelihoodsALL,axis=0)         
    log.debug("medainALL = "+str(medainALL))
    s = '\nmedian value for all chains = '+str(medainALL)
    ## calculate location of medianALL in each chain
    for filename in mcmcFnames:
        if os.path.exists(filename):
            (head,data) = loadFits(filename)
            chiSqsChain = data[:,11]
            likelihoodsChain = np.exp(-nu*chiSqsChain/2.0)
            #medianChain = np.median(chiSquaredsChain)
            for i in range(1,len(likelihoodsChain)):
                if likelihoodsChain[i]>medainALL:
                    #print 'chiSqsChain[i] = '+str(chiSqsChain[i])
                    burnInLength = i+1
                    break
            burnInLengths.append(burnInLength)
            s2 = "\nfor chain #"+str(head['chainNum'])
            s2 += "\nTotal number of points in the chain = "+str(len(chiSqsChain))+"\n"
            s2 += "Burn-in length = "+str(burnInLength)+"\n"
            s+=s2
            log.debug(s2)
    log.debug(s)
    return (s,burnInLengths)

def burnInStripper(mcmcFnames,burnInLengths):
    """
    Strip the initial burn-in off each chain.
    """
    newFnames=[]
    for i in range(0,len(mcmcFnames)):
        filename = mcmcFnames[i]
        burnIn = burnInLengths[i]
        if os.path.exists(filename):
            (head,data) = loadFits(filename)
            ##strip burn-in and write to new fits             
            hdu = pyfits.PrimaryHDU(data[burnIn:,:])
            hdulist = pyfits.HDUList([hdu])
            newHead = hdulist[0].header
            for key in head:
                newHead[key]=(head[key],head.comments[key])
            n = os.path.splitext(os.path.basename(filename))
            newFname = os.path.join(os.path.dirname(mcmcFnames[0]),n[0]+'_BIstripped.fits')
            hdulist.writeto(newFname)
            log.info("output file written to:below\n"+newFname)
            hdulist.close()
            newFnames.append(newFname)
            log.debug("burn-in stripped file written to:\n"+newFname)
    return newFnames

def gelmanRubinCalc(mcmcFileList,nMCMCsamp=1):
    """
    Calculate Gelman-Rubin statistic for each varying param.
    Input MUST be the list of more than one MCMC chain files.
    """
    GRs=[]
    Ts = []
    grStr = '\n'+'-'*30+"\nGelman-Rubin Results:\n"+'-'*30+'\n'
    Lcfloat = float(nMCMCsamp)
    if os.path.exists(mcmcFileList[0]):
        log.info("Starting to calculate R&T")
        ###########################################################
        ## stage 1 ->  load up values for each param in each chain.
        ## allStg1vals = [chain#, param#, (mean,variance,Lc)]
        ## stage 2 ->  Use them to compare between chains 
        ##             then calc R and T.
        ###########################################################
        (head,data) = loadFits(mcmcFileList[0])
        (paramList,paramStrs,paramFileStrs) = getParStrs(head,latex=False)
        Nc = len(mcmcFileList)
        ##start stage 1
        allStg1vals=np.zeros((Nc,len(paramList),3))
        for i in range(0,len(mcmcFileList)):
            log.debug("Starting to calc chain #"+str(i)+' GR values')
            (head,data) = loadFits(mcmcFileList[i])
            allStg1vals[i,:,2]=data.shape[0]
            for j in range(0,len(paramList)):
                log.debug("calculating stage 1 of GR for chain #"+str(i)+", param: "+paramStrs[j])
                allStg1vals[i,j,0]=np.mean(data[:,paramList[j]])
                allStg1vals[i,j,1]=np.var(data[:,paramList[j]])
        ##start stage 2         
        rHighest = 0
        tLowest = 1e9       
        rHighStr = ''
        tLowStr = ''
        for j in range(0,len(paramList)):
            log.debug("Starting stage 2 for param: "+paramStrs[paramList[j]])
            Ncfloat = float(Nc)
            ##calc R
            W = 0
            for i in range(0,Nc):
                #Lcfloat = float(allStg1vals[i,j,2])
                W+=(Lcfloat/(Lcfloat-1.0))*allStg1vals[i,j,1]
            W=W/Ncfloat
            V = np.mean(allStg1vals[:,j,1])+(Ncfloat/(Ncfloat-1.0))*np.var(allStg1vals[:,j,0])
            R=np.NaN
            if W!=0:
                R = np.sqrt(V/W)
            GRs.append(R)
            ##calc T
            #B = (Ncfloat/(Ncfloat-1.0))*np.var(allStg1vals[:,j,0]*allStg1vals[:,j,2])
            B = (Ncfloat/(Ncfloat-1.0))*np.var(allStg1vals[:,j,0])*Lcfloat
            #Uses the mean of Lc values
            T = np.NAN
            if B!=0:
                #T = np.mean(allStg1vals[:,j,2])*Ncfloat*np.min([(V/B),1.0])
                T = Lcfloat*Ncfloat*np.min([(V/B),1.0])
            Ts.append(T)       
            grStr+=paramStrs[paramList[j]]+" had R = "+str(R)+", T = "+str(T)+'\n'
            if T<tLowest:
                tLowest=T
                tLowStr="Lowest T = "+str(T)+' for '+paramStrs[paramList[j]]+'\n'
            if R>rHighest:
                rHighest=R
                rHighStr="Highest R = "+str(R)+' for '+paramStrs[paramList[j]]+'\n'
        grStr+='\nWorst R&T values were:\n'+rHighStr+tLowStr+'\n'
    else:
        log.critical("Gelman-Rubin stat can NOT be calculated as file does not exist!!:\n"+chainDataFileList[0])
    
    return (GRs,Ts,grStr)

def timeStrMaker(deltaT):
    """
    Convert a time in seconds into a nicer string.
    """
    timeStr = ''
    if deltaT>60:
        if deltaT>(60*60):
            hr = int(deltaT//(60*60))
            min = int((deltaT-hr*60*60)/60.0)
            timeStr = str(hr)+" hours and "+str(min)+" minutes"
        else:
            min = int(deltaT//(60))
            timeStr = str(min)+" minutes and "+str(int(deltaT-min*60))+" seconds"
    else:
        timeStr = str(int(deltaT))+' seconds'
    return timeStr

def getParStrs(head,latex=True):
    """
    Return matching paramList, paramStrs, paramFileStrs for provided header.
    """
    paramList = getParInts(head)    
    paramFileStrs = ['M1','M2','parallax','Omega','e','To', 'Tc','P','i','omega','a_total','chiSquared','K']
    paramStrs = ['M1 [Msun]','M2 [Msun]','Parallax [mas]','Omega [deg]','e','To [JD]', 'Tc [JD]','P [Yrs]','i [deg]','omega [deg]','a_total [AU]','chiSquared','K [m/s]']
    if latex:
        paramStrs = ['$M_1$ [$M_{\odot}$]','$M_2$ [$M_{\odot}$]','$Parallax$ [mas]','$\Omega$ [deg]','$e$','$T_o$ [JD]', '$T_c$ [JD]','$P$ [Yrs]','$i$ [deg]','$\omega$ [deg]','$a_{total}$ [AU]','$\chi^2$','$K$ [m/s]']

    if head["nRVdsets"]>0:
        for dataset in range(1,head["nRVdsets"]+1):
            paramFileStrs.append('offset_'+str(dataset))
            if latex:
                paramStrs.append("$offset_{"+str(dataset)+"}$ [m/s]")
            else:
                paramStrs.append('offset '+str(dataset)+' [m/s]')
    return (paramList,paramStrs,paramFileStrs)

def loadDIdata(filename):
    """
    Load the astrometry data into a numpy array.
    
    file format:
    title 
    column headers
    data
    .
    .
    .
    
    Data must be in the columns:
    obsDate[JD] x["] x_error["] y["] y_error["]
    """
    if filename[-4:]!='.dat':
        filename = filename+'.dat'
    file = open(filename, 'r')
    diData = []
    lines = file.readlines()
    file.close()
    for line in lines:
        #log.debug("line was:'"+line+"'")#$$$$$$$$$$$$$$$$$$$$$$$$
        if len(line.split())>2:
            if line.split()[0].replace('.','',1).isdigit() and line.split()[3].replace('.','',1).replace('-','',1).isdigit():
                diData.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4])])  
    return np.array(diData)
    
def loadRVdata(filename):
    """
    Load the radial velocity data into a numpy array.  Provided jitter values will be added in quadrature with 
    the errors.
    
    file format:
    title 
    column headers
    data
    .
    .
    Empty line between data sets
    data
    .
    .
    
    Data must be in the columns:
    obsDate[JD] RV[m/s] RV_error[m/s] jitter[m/s] datasetNumber[int]
    NOTE: datasetNumber is optional, if not provided they will be automatically set to 0,1,2... following the order of the data in the file.
          If jitter is not provided, it will be assumed zero.
    """
    if filename[-4:]!='.dat':
        filename = filename+'.dat'
    file = open(filename, 'r')
    lines = file.readlines()
    file.close()
    rvData = []
    datasetNumLast = 0
    jitterLast = 0
    lastWasDataLine=False
    thisIsDataLine = False
    for line in lines:
        #print "line = "+line
        lastWasDataLine=thisIsDataLine
        thisIsDataLine=False
        if len(line.split())>2:
            if line.split()[0].replace('.','',1).isdigit() and line.split()[1].replace('.','',1).replace('-','',1).isdigit():
                thisIsDataLine=True
                curDataAry = [float(line.split()[0]),float(line.split()[1])]
                #if jitter was provided on first line of data set
                if len(line.split())>3:
                    try:
                        jitterLast = float(line.split()[3])
                    except:
                         log.error("could not convert 4th element of split into jitter.  4th element was: "+str(line.split()[3]))
                curDataAry.append(np.sqrt(float(line.split()[2])**2+jitterLast**2))
                #if datasetNum was provided on first line of data set
                if len(line.split())>4:
                    try:
                        datasetNumLast = float(line.split()[4])
                    except:
                        log.error("could not convert 5th element of split into datasetNum.  5th element was: "+str(line.split()[4]))
                curDataAry.append(datasetNumLast)
                #print repr(curDataAry)
                rvData.append(curDataAry)
        if lastWasDataLine and (thisIsDataLine==False):
            jitterLast = 0
            datasetNumLast+=1
            #print 'incrementing datasetNum'
    return np.array(rvData)
    
def loadRealData(filenameRoot,dataMode='3D'):
    """
    Load the observed real data into a numpy array.
    This will be a combination of the RV and DI data,sorted into cronological order.
    filenameRoot would be the absolute path plus the prepend to the settings files.
    ex. '/run/..../SMODT/settings_and_inputData/FakeData_'
    """
    diEpochs = []
    rvEpochs = []
    if dataMode!='RV':
        diFilename = filenameRoot+'DIdata.dat'
        #print 'using diFilename = '+diFilename        
        if os.path.exists(diFilename):
            diData = loadDIdata(diFilename)
            diEpochs = diData[:,0]
    if dataMode!='DI':
        rvFilename = filenameRoot+'RVdata.dat'
        #print 'using rvFilename = '+rvFilename
        if os.path.exists(rvFilename):
            rvData = loadRVdata(rvFilename)
            rvEpochs = rvData[:,0]
    #print 'rvData = '+repr(rvData)
    #for i in range(0,rvData.shape[0]):
    #    print 'ORIG rv data = '+str(rvData[i,0])+', '+str(rvData[i,1])+", "+str(rvData[i,2])+", "+str(rvData[i,3])
    ##load in epochs from both sets, sort and kill double entries
    epochsTemp = np.concatenate((diEpochs,rvEpochs))
    epochsTemp.sort()
    epochs = []
    for epoch in epochsTemp:
        if epoch not in epochs:
            epochs.append(epoch)
    epochs = np.array(epochs)
    realData = np.zeros((epochs.shape[0],8))
    ##set error values to 1e6 which signals not to calculate the predicted version in orbit.cc
    realData[:,2]=realData[:,4]=realData[:,6]=1e6
    realData[:,0]=epochs[:]
    for i in range(epochs.shape[0]):
        if len(diEpochs)>0:
            if epochs[i] in diData[:,0]:
                realData[i,1:5]=diData[np.where(diData[:,0]==epochs[i])[0],1:]
        if len(rvEpochs)>0:
            if epochs[i] in rvData[:,0]:
                realData[i,5:]=rvData[np.where(rvData[:,0]==epochs[i])[0],1:]
    #print 'dataMode'+dataMode+'->realData = '+repr(realData)
    #for i in range(0,realData.shape[0]):
    #    print 'realData = '+str(realData[i,0])+', '+str(realData[i,5])+", "+str(realData[i,6])+", "+str(realData[i,7])
    return realData
            
def loadSettingsDict(filenameRoot):
    """
    Load the values from both the simple (symSettingsSimple.py) and advanced (symSettingsAdvanced.py)
    into a dictionary for use throughout the simulation and post-processing.
    Those that are deemed useful will be loaded in as a tuple with a comment for later adding to 
    the resulting simulation data file fits header.
    NOTE: the first step is to copy these files to standardized names so they can be called in to 
          use.  They will overwrite the files:
          SMODT/tools/temp/simpleSettings.py   &   advancedSettings.py 
    
    filenameRoot would be the absolute path plus the prepend to the settings files.
    ex. '/run/..../SMODT/settings_and_inputData/FakeData_'
    """
    ## A BIT HACKY FOR NOW, NEED TO FIND A CLEANER WAY TO DO THIS!?!?! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    cwd = os.getcwd()
    smodtHeadDir = filenameRoot.split("smodt2")[0]
    try:
        os.remove(os.path.join(smodtHeadDir,'smodt2/tools/temp/settingsSimple.py'))
        os.remove(os.path.join(smodtHeadDir,'smodt2/tools/temp/settingsAdvanced.py'))
        os.remove(os.path.join(smodtHeadDir,'smodt2/tools/temp/constants.py'))
    except:
        temp=True
    shutil.copy(filenameRoot+'settingsSimple.py',os.path.join(smodtHeadDir,'smodt2/tools/temp/settingsSimple.py'))
    shutil.copy(filenameRoot+'settingsAdvanced.py',os.path.join(smodtHeadDir,'smodt2/tools/temp/settingsAdvanced.py'))
    shutil.copy(os.path.join(smodtHeadDir,'smodt2/tools/constants.py'),os.path.join(smodtHeadDir,'smodt2/tools/temp/constants.py'))
    if False:
        print 'Copied simple to:\n'+os.path.join(smodtHeadDir,'smodt2/tools/temp/settingsSimple.py')
        print 'Copied advanced to:\n'+os.path.join(smodtHeadDir,'smodt2/tools/temp/settingsAdvanced.py')
        print 'Copied constants to:\n'+os.path.join(smodtHeadDir,'smodt2/tools/temp/constants.py')
    os.chdir(os.path.join(smodtHeadDir,'smodt2'))
    from tools.temp.settingsAdvanced import settingsDict
    os.chdir(cwd)
    
    #######################################################
    ## determine argPeriOffsetRV and argPeriOffsetDI values
    #######################################################
    omegaFdi = 0
    omegaFrv = 0
    #first using RV special bools
    if (settingsDict['primeRVs'][0] and settingsDict['fitPrime'][0]):
        omegaFdi=-180.0
    elif (settingsDict['primeRVs'][0] and(settingsDict['fitPrime'][0]==False)):
        omegaFrv=180.0
    #now update due to fixed argPeriPlus values
    omegaFdi+=settingsDict['omegaPdi'][0]
    omegaFrv+=settingsDict['omegaPrv'][0]
    settingsDict['omegaFdi'] = (omegaFdi,"Total fixed val added to DI omega in model")
    settingsDict['omegaFrv'] = (omegaFrv,"Total fixed val added to RV omega in model")
    log.debug("Setting fixed omega offsets to:\nomegaFdi = "+str(omegaFdi)+"\nomegaFrv = "+str(omegaFrv))
    ##In DI mode can only find Mtotal, thus push all mass into M1 and kill M2
    if settingsDict['dataMode'][0]=='DI':
        settingsDict['mass1MIN']=settingsDict['mass1MIN']+settingsDict['mass2MIN']
        settingsDict['mass1MAX']=settingsDict['mass1MAX']+settingsDict['mass2MAX']
        settingsDict['mass2MIN']=0
        settingsDict['mass2MAX']=0
        log.debug("DI dataMode, so pushed all mass range vals into M1 and set ones for M2 to zero")
        
    return settingsDict
    
def startup(argv,rootDir,rePlot=False):
    """
    -Figure out important directories
    -copy settings files to temp directory to combine them into the master settings dict
    -get master settings dict
    -make output folder
    -remake the SWIG tools?
    -copy all SMODT2.0 code into output dir for emergencies.    
    
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
    settingsDict = loadSettingsDict(rootDir+'settings_and_inputData/'+prepend)
    settingsDict['smodtdir']=rootDir
    settingsDict['settingsDir']=os.path.join(settingsDict['smodtdir'],'settings_and_inputData/')
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
            os.chdir(os.path.join(settingsDict['smodtdir'],'tools/cppTools/'))
            os.system('make clean')
            os.system('make')
            os.chdir(cwd)
            log.debug("-"*45+" Done re-making CPP/SWIG tools "+45*"-")
        ## copy all of current code to output directory
        codeCopyDir = os.path.join(settingsDict['finalFolder'],'codeUsed')
        os.mkdir(codeCopyDir)
        log.debug('Copying all files in the RESULTS folder over to DropBox folder:\n '+codeCopyDir)
        copytree(settingsDict['smodtdir'], codeCopyDir)
        
    return settingsDict
        

def cleanUp(settingsDict,stageList,allFname):
    """
    Clean up final directory after simulation completes
    """
    #os.mkdir(settingsDict['finalFolder'])
    
    delFiles = []
    ##get chain data filenames to delete
    if settingsDict["delChains"]:
        for stage in stageList:
            fnames = glob.glob(os.path.join(settingsDict['finalFolder'],"outputData"+stage+"*.fits"))
            for i in range(0,len(fnames)):
                delFiles.append(fnames[i])
    ##get combined data filename to delete
    if settingsDict["delCombined"]:
        delFiles.append(allFname)
        if settingsDict['delBurn'][0]:
            ##the burn-in was stripped in the final file, so kill the non-stripped version if it exists
            nm = os.path.join(os.path.dirname(allFname),'combined*data.fits')
            if os.path.exists(nm):
                delFiles.append(nm)
                
    ##try to delete files
    for fname in delFiles:
        try:
            log.debug('Deleting file: '+os.path.basename(fname))
            os.remove(fname) 
        except:
            log.error('Failed to delete file: '+os.path.basename(fname))

def writeFits(baseFilename,data,settingsDict):
    """
    Data will be written to a fits file with a single PrimaryHDU,
    with the .header loaded up with the tuples from the settingsDict 
    and .data = provided data.
    File will be stored in the 'finalFolder' directory from the settingsDict.
    """
    if '.fits' not in baseFilename:
        baseFilename=baseFilename+'.fits'
    outFname = os.path.join(settingsDict['finalFolder'],baseFilename)
    hdu = pyfits.PrimaryHDU(data)
    hdulist = pyfits.HDUList([hdu])
    header = hdulist[0].header
    ##load up header with tuples from settingsDict
    for key in settingsDict:
        if type(settingsDict[key])==tuple:
            header[key]=settingsDict[key]
            if len(settingsDict[key][1])>47:
                log.warning("comment too long for pyfits headers:"+settingsDict[key][1])
    hdulist.writeto(outFname)
    log.info("output file written to:below\n"+outFname)
    hdulist.close()
    ## check resulting fits file header
    if False:
        f = pyfits.open(os.path.join(settingsDict['finalFolder'],baseFilename),'readonly')
        head = f[0].header
        f.close()
        if False:
            for key in head:
                print key+' = '+repr((header[key],header.comments[key]))
                #print 'type(header[key] = '+repr(type(header[key]))
        print '\n\nEntire Header as a repr:\n'+repr(head)
    return outFname
    
def loadFits(filename):
    """
    Load in a fits file written by SMODT.
    Return (header dict, data)
    """
    if os.path.exists(filename):
        f = pyfits.open(filename,'readonly')
        head = f[0].header
        data = f[0].data
        f.close()
    else:
        log.critical("fits file does not exist!!! filename =\n"+str(filename))
        head=data=False
    return (head,data)

def combineFits(filenames,outFname):
    """
    combine the data in multiple SMODT2 fits files together.
    Used primarily for after multi-process runs.
    """
    nFiles = len(filenames)
    (head0,dataALL) = loadFits(filenames[0])
    for filename in filenames:
        (head,data) = loadFits(filename)
        dataALL = np.concatenate((dataALL,data))
    hdu = pyfits.PrimaryHDU(dataALL)
    hdulist = pyfits.HDUList([hdu])
    header = hdulist[0].header
    for key in head0:
        if key=='NSAMPLES':
            ##load in total number of samples for this combined file
            header['NSAMPLES'] = (int(head0['NSAMPLES'])*nFiles,head0.comments['NSAMPLES'])
        else:
            header[key] = (head0[key],head0.comments[key])
    hdulist.writeto(outFname)
    hdulist.close()
    log.info("output file written to:below\n"+outFname)
    
def summaryFilePart1(settingsDict,stageList,finalFits,clStr,burnInStr,bestFit,grStr):
    """
    Make a txt file that summarizes the results nicely.
    """
    summaryFname = os.path.join(settingsDict['finalFolder'],'RESULTS.txt')
    f = open(summaryFname,'w')
    (head,data) = loadFits(finalFits)
    totalSamps = head['NSAMPLES']
    (paramList,paramStrs,paramFileStrs) = getParStrs(head,latex=False)
    f.write('\n'+'-'*7+'\nBasics:\n'+'-'*7)
    f.write('\nparamList:\n'+repr(paramList))
    f.write('\nparamStrs:\n'+repr(paramStrs))
    f.write('\nparamFileStrs:\n'+repr(paramFileStrs))
    try:
        ## try to make and write the more advanced summary strings to the file
        nusStr = "\nnu values were: [total,DI,RV] = ["+str(head['NU'])+", "+str(head['NUDI'])+", "+str(head['NURV'])+"]\n"
        f.write(nusStr)
        stgNsampStrDict = {"MC":"nSamples","SA":"nSAsamp","ST":"nSTsamp","MCMC":"nSamples"}
        numFilesStr = '\nTotal # of files that finished each stage were:\n'
        chiSquaredsStr = '\nBest Reduced Chi Squareds for each stage were:\n'
        log.debug('ln680:summaryFilePart1')  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        for stage in stageList:
            fnames = np.sort(glob.glob(os.path.join(settingsDict['finalFolder'],"outputData"+stage+"*.fits")))
            if (stage=="MCMC")and settingsDict["delBurn"][0]:
                fnames = np.sort(glob.glob(os.path.join(settingsDict['finalFolder'],"outputData"+stage+"*BIstripped.fits")))
            numFilesStr+=stage+' = '+str(len(fnames))+", each with "+str(settingsDict[stgNsampStrDict[stage]][0])+" samples\n"
            if len(fnames)>0:
                log.debug('len(fnames) = '+repr(len(fnames)))#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                chiSquaredsStr+=stage+" = ["
                for fname in fnames: 
                    try:
                        bestFit2 = findBestOrbit(fname)
                        chiSquaredsStr+=str(bestFit2[11]/float(head['NU']))+', '
                    except:
                        log.error("A problem occurred while trying to find best fit of:\n"+fname)
                chiSquaredsStr = chiSquaredsStr[:-2]+']\n'
        log.debug('ln696:summaryFilePart1')#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        numFilesStr+="\n"+"*"*65+"\nThe final combined file was for a total of "+str(totalSamps)+" samples\n"+"*"*65+'\n'
        f.write(numFilesStr)
        f.write(chiSquaredsStr)
        bestStr = '\n'+'-'*21+'\nBest fit values were:\n'+'-'*21+'\n'
        log.debug('ln701:summaryFilePart1')#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        for i in range(len(bestFit)):
            if i==2:
                bestStr+=paramStrs[2]+" = "+str(bestFit[2])
                if (bestFit[2]!=0):
                    bestStr+=", OR  "+str(1.0/(bestFit[2]/1000.0))+'[PC]\n'
                    #print 'bestStr '+bestStr#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                else:
                    bestStr+='\n'
            elif i==1:
                if bestFit[i]<0.02:
                    mJupMult=(const.KGperMsun/const.KGperMjupiter)
                    bestStr+=paramStrs[i]+" = "+str(bestFit[i])+", OR "+str(bestFit[i]*mJupMult)+' in [Mjupiter]\n'
                else:
                    bestStr+=paramStrs[i]+" = "+str(bestFit[i])+'\n'
            else:
                bestStr+=paramStrs[i]+" = "+str(bestFit[i])+'\n'
        bestStr+='\n'+'*'*45+'\nBEST REDUCED CHISQUARED = '+str(bestFit[11]/float(head['nu']))+'\n'+'*'*45
        f.write(bestStr)
        log.debug('ln714:summaryFilePart1')#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    except:
        log.critical("A problem occured while trying to produce advanced summary strings.")
    f.write('\n'+clStr)
    f.write('\n'+burnInStr)
    f.write('\n'+grStr)
    f.close()
    log.info("Summary file written to:\n"+summaryFname)  

def summaryFilePart2(settingsDict,effPtsStr,allTime,postTime):
    """
    Finish writting the important summary information for the post-processing 
    steps that can take along time.
    """
    summaryFname = os.path.join(settingsDict['finalFolder'],'RESULTS.txt')
    f = open(summaryFname,'a')
    f.write('\n\nPost-Processing took: '+timeStrMaker(postTime)+'\n')
    f.write('Total simulation took: '+timeStrMaker(allTime)+'\n')
    f.write('\n'+effPtsStr)
    f.close()
    log.info("Summary file written to:\n"+summaryFname)
    
    
def copyToDB(settingsDict):
    """
    Copy vital results files to Dropbox.
    """
    fnamesALL = []
    for extension in ['pdf','txt','log']:
        fnames = glob.glob(os.path.join(settingsDict['finalFolder'],"*."+extension))
        log.debug('found files to copy to DB:\n'+repr(fnames))
        for name in fnames:
            fnamesALL.append(name)
    dbDir = os.path.join(settingsDict['dbFolder'],settingsDict['outRoot'])
    os.mkdir(dbDir)
    log.debug('DB dir is:\n'+repr(dbDir))
    for f in fnamesALL:
        try:
            log.debug('trying to copy file:\n'+repr(os.path.basename(f)))
            shutil.copy(f,os.path.join(dbDir,os.path.basename(f)))
        except:
            log.error('failed to move file:\n'+f+'\nintto DB folder:\n'+dbDir)
    log.info("vital results files copied to DB folder:\n"+dbDir)
    
def getParInts(head):
    """
    convert string version of paramInts into a list again.
    """
    s = head['parInts'] 
    ints = s.split("[")[1].split("]")[0].split(',')
    parInts = []
    for i in ints:
        parInts.append(int(i))  
    return parInts
        
def confLevelFinder(filename, colNum=False, returnData=False, returnChiSquareds=False, returnBestDataVal=False):
    """
    A function to find the 68.3 and 95.4% confidence levels in a given output data file's column.
    
    return [[68.3% minimum, 68.3% maximum],[95.5% minimum, 95.5% maximum]]
    
    columnNum must be an int.    
    """
    verboseInternal = False
    bestCentered = False
    log.debug('Inside confLevelFinder')
    outStr=''
    if os.path.exists(filename):
        (dataAry,chiSquareds,[bestDataVal,dataMedian,dataValueStart,dataValueMid,dataValueEnd]) = dataReader(filename, colNum)
    
        if len(dataAry>0) or (dataValueStart!=dataValueMid!=dataValueEnd):
            #Convert data array to a sorted numpy array
            dataAry = np.sort(dataAry)
            size = dataAry.size
        
            if bestCentered:
                mid = np.where(dataAry==bestDataVal)[0][0]
            else:
                mid=size//2
                
            minLoc68=mid-int(float(size)*0.683)//2
            if minLoc68<0:
                minLoc68 = 0
            maxLoc68 = mid+int(float(size)*0.683)//2
            if maxLoc68>(size-1):
                maxLoc68 = size
            minLoc95=mid-int(float(size)*0.958)//2
            if minLoc95<0:
                minLoc95 = 0
            maxLoc95= mid+int(float(size)*0.958)//2
            if maxLoc95>(size-1):
                maxLoc95 = size
            
            conf68Vals = [dataAry[minLoc68],dataAry[maxLoc68]]
            conf95Vals = [dataAry[minLoc95],dataAry[maxLoc95]]
            conf68ValsRough=[]
            conf95ValsRough=[]
            
            if ((len(conf68Vals)==0) or (len(conf95Vals)==0)):
                if (len(conf68Vals)==0):
                    log.error('confLevelFinder: ERROR!!! No FINE 68.3% confidence levels were found')
                    if (len(conf68ValsRough)==0):
                        log.error('confLevelFinder: ERROR!!! No ROUGH 68% confidence levels were found, so returning [0,0]')
                        conf68Vals = [0,0]
                    else:
                        conf68Vals = conf68ValsRough
                        log.error("confLevelFinder: Had to use ROUGH 68% [68,69] as no FINE 68.3% was found. So, using range "+repr(conf68Vals))                
                if (len(conf95Vals)==0):
                    log.error('confLevelFinder: ERROR!!! No FINE 95.4% confidence levels were found')
                    if (len(conf95ValsRough)==0):
                        log.error('confLevelFinder: ERROR!!! No ROUGH 95% confidence levels were found, so returning [0,0]')
                        conf95Vals = [0,0]
                    else:
                        conf95Vals = conf95ValsRough
                        log.error("confLevelFinder: Had to use ROUGH 95% [95,96] as no FINE 95.4% was found. So, using range "+repr(conf95Vals))
        else:
            ## There was no useful data, so return values indicating that
            dataAry=bestDataVal=dataMedian=dataValueStart
            conf68Vals = [dataValueStart,dataValueStart]
            conf95Vals = [dataValueStart,dataValueStart]
            chiSquareds = 0
            
        mJupMult=(const.KGperMsun/const.KGperMjupiter)
        s= "\nFinal 68% range values are: "+repr(conf68Vals)+'\n'
        s=s+"Final 95% range values are: "+repr(conf95Vals)+'\n'
        if bestCentered:
            s=s+ "\nerror is centered on best \n"
            s=s+"68.3% error level = "+str(bestDataVal-conf68Vals[0])+'\n'
            s=s+" ->   "+str(bestDataVal)+'  +/-  '+str(bestDataVal-conf68Vals[0])+'\n'
            if (colNum==1) and (bestDataVal<0.02):
                s=s+"Or in Mjup: ->   "+str(bestDataVal*mJupMult)+'  +/-  '+str(bestDataVal*mJupMult-conf68Vals[0]*mJupMult)+'\n'
            outStr+=s
        else:
            s=s+ "\nerror is centered on Median \n"
            s=s+"68.3% error level = "+str(dataMedian-conf68Vals[0])
            s=s+" ->   "+str(dataMedian)+'  +/-  '+str(dataMedian-conf68Vals[0])+'\n'
            if (colNum==1) and (dataMedian<0.02):
                s=s+"Or in Mjup: ->   "+str(dataMedian*mJupMult)+'  +/-  '+str(dataMedian*mJupMult-conf68Vals[0]*mJupMult)+'\n'
            outStr+=s
        s=s+'\n'+75*'-'+'\n Leaving confLevelFinder \n'+75*'-'+'\n'
        log.debug(s)
        
        if verboseInternal:
            print 'returnData = '+repr(returnData)+', returnChiSquareds = '+repr(returnChiSquareds)+', returnBestDataVal = '+repr(returnBestDataVal)
        ## return requested form of results
        if (returnData and returnChiSquareds and (returnBestDataVal==False)):
            if verboseInternal:
                print 'returning first 3'
            returnList =  ([conf68Vals,conf95Vals],dataAry, chiSquareds)
        elif (returnData and returnChiSquareds and returnBestDataVal):
            if verboseInternal:
                print 'returning all 4'
            returnList =   ([conf68Vals,conf95Vals],dataAry, chiSquareds, bestDataVal)
        elif (returnData and (returnChiSquareds==False)and (returnBestDataVal==False)):
            if verboseInternal:
                print 'returning data only'
            returnList =   ([conf68Vals,conf95Vals],dataAry)
        elif (returnData and (returnChiSquareds==False) and returnBestDataVal):
            if verboseInternal:
                print 'returning data and bestval'
            returnList =   ([conf68Vals,conf95Vals],dataAry, bestDataVal,outStr) ##MODIFIED
        elif ((returnData==False) and returnChiSquareds):
            if verboseInternal:
                print 'returning just chiSquareds'
            returnList =   ([conf68Vals,conf95Vals], chiSquareds)
        elif ((returnData==False) and returnChiSquareds and returnBestDataVal):
            if verboseInternal:
                print 'returning chiSquareds and bestval'
            returnList =   ([conf68Vals,conf95Vals], chiSquareds, bestDataVal)
        elif ((returnData==False)and(returnChiSquareds==False) and returnBestDataVal):
            if verboseInternal:
                print 'returning CLevels and bestval'
            returnList = ([conf68Vals,conf95Vals], bestDataVal)
        elif ((returnData==False)and(returnChiSquareds==False) and (returnBestDataVal==False)):
            if verboseInternal:
                print 'returning only CLevels'
            returnList =   [conf68Vals,conf95Vals]
        
        return returnList 
        
    else:
        log.critical( "confLevelFinder: ERROR!!!! file doesn't exist")
        
def dataReader(filename, colNum=0):
    """
    Read in the data for a single column of data.
    """
    verboseInternal = False
    ## First get ranges of param and ChiSquared values
    log.debug('\nOpening and finding ranges for data in column # '+str(colNum))
    
    ## Check if file has useful data for that column#
    (head,data) = loadFits(filename)
    if head!=False:
        TotalSamples=data.shape[0]
        dataAry = data[:,colNum]
        chiSquareds = data[:,11]
        bestDataVal = dataAry[np.where(chiSquareds==np.min(chiSquareds))][0]          
        return (dataAry,chiSquareds,[bestDataVal,np.median(dataAry),dataAry[0],dataAry[len(dataAry)//2],dataAry[-1]])                
                                         
def findBestOrbit(filename):        
    """
    Find the orbital elements for the best fit in a SMODT format fits file.
    """             
    log.debug("trying to find best orbit in file:\n"+filename)     
    (head,data) = loadFits(filename)
    chiBest = np.min(data[:,11])
    loc = np.where(data[:,11]==chiBest)
    orbBest = data[loc[0][0],:]
    log.info("Best fit found to be:\n"+repr(orbBest))
    return orbBest
                                         
                                         
def copytree(src, dst):
    """
    Recursively copy a directory and its contents to another directory.
    
    WARNING: this is not advised for higher level folders as it can also copy subfolders 
    thus leading to a very large copy command if not careful.
    
    Code taken and simplified from:
    http://stackoverflow.com/questions/1868714/how-do-i-copy-an-entire-directory-of-files-into-an-existing-directory-using-pyth
    """
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            try:
                shutil.copytree(s, d)
                log.debug("Copying:\n "+repr(s)+'\nto:\n'+repr(d))
            except:
                log.error('FAILED while copying:\n'+repr(s)+'\nto:\n'+repr(d))
        else:
            try:
                shutil.copy2(s, d)
                log.debug("Copying:\n "+repr(s)+'\nto:\n'+repr(d))
            except:
                log.error('FAILED while copying:\n'+repr(s)+'\nto:\n'+repr(d))                        
    