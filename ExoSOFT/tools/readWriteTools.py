#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import numpy as np
import exoSOFTlogger
import generalTools as genTools
import pyfits
import os
import sys
import shutil
import warnings
warnings.simplefilter("error")

log = exoSOFTlogger.getLogger('main.rwTools',lvl=100,addFH=False)  

def dataReader(filename, colNum=0):
    """
    Read in the data for a single column of data.
    """
    try:
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
    except:
        log.critical("a problem occured while trying to load data file:\n"+filename)
        return False

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
    OR if pasa key in settingsAdvanced ==True, then:
    obsDate[JD] PA[deg] PA_error[deg] SA["] SA_error["]
    """
    diData = []
    try:
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
        diData = np.array(diData)
    except:
        log.critical("a problem occured while trying to load DI data. \nlease check it is formatted correctly.")
    return diData
    
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
    rvData=[]
    try:
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
        rvData = np.array(rvData)
    except:
        log.critical("a problem occured while trying to load RV data.  \nPlease check it is formatted correctly.")
    return rvData
    
def loadRealData(diFilename='',rvFilename='',dataMode='3D'):
    """
    Load the observed real data into a numpy array.
    This will be a combination of the RV and DI data,sorted into cronological order.
    Both filenames must be the full path to the files.
    """
    realData = False
    try:
        diEpochs = []
        rvEpochs = []
        if dataMode!='RV':
            #print 'using diFilename = '+diFilename        
            if os.path.exists(diFilename):
                diData = loadDIdata(diFilename)
                diEpochs = diData[:,0]
        if dataMode!='DI':
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
    except:
        log.critical("An error occured while trying to load data!!")
    return realData
            
def loadSettingsDict(ExoSOFTdir,settFilePath):
    """
    Load the values from both the simple (symSettingsSimple.py) and advanced (symSettingsAdvanced.py)
    into a dictionary for use throughout the simulation and post-processing.
    Those that are deemed useful will be loaded in as a tuple with a comment for later adding to 
    the resulting simulation data file fits header.
    NOTE: the first step is to copy these files to standardized names so they can be called in to 
          use.  They will overwrite the files:
          ExoSOFT/tools/temp/simpleSettings.py   &   advancedSettings.py 
    
    filenameRoot would be the absolute path plus the prepend to the settings files.
    ex. '/run/..../ExoSOFT/settings_and_inputData/FakeData_'
    """
    ## A BIT HACKY FOR NOW, NEED TO FIND A CLEANER WAY TO DO THIS!?!?! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ## at same time, maybe consider entirely new way to handle the settings and funcs.  
    ## Looked at YAML, but it is just convoluted python code with no advantages.
    ## considered default files for the priors and ranges that could be overriden by those in the advanced dict, but gave up.  Also seems too convoluted.
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ## idea: user can pass in a filename root, OR a full path to the settings file.
    ## If just root, assume it is in standard settings dir.
    ## copy it to the temp dir with the constants file, import and then delete temp versions.
    ##
    cwd = os.getenv('PWD')
    if False:
        try:
            os.remove(os.path.join(ExoSOFTdir,'tools/temp/settingsSimple.py'))
            os.remove(os.path.join(ExoSOFTdir,'tools/temp/settingsAdvanced.py'))
            os.remove(os.path.join(ExoSOFTdir,'tools/temp/constants.py'))
        except:
            temp=True
        shutil.copy(filenameRoot+'settingsSimple.py',os.path.join(ExoSOFTdir,'tools/temp/settingsSimple.py'))
        shutil.copy(filenameRoot+'settingsAdvanced.py',os.path.join(ExoSOFTdir,'tools/temp/settingsAdvanced.py'))
        shutil.copy(os.path.join(ExoSOFTdir,'tools/constants.py'),os.path.join(ExoSOFTdir,'tools/temp/constants.py'))
        if False:
            print 'Copied simple to:\n'+os.path.join(ExoSOFTdir,'tools/temp/settingsSimple.py')
            print 'Copied advanced to:\n'+os.path.join(ExoSOFTdir,'tools/temp/settingsAdvanced.py')
            print 'Copied constants to:\n'+os.path.join(ExoSOFTdir,'tools/temp/constants.py')
        os.chdir(ExoSOFTdir)
        from tools.temp.settingsAdvanced import settingsDict
        os.chdir(cwd)
    else:
        try:
            os.remove(os.path.join(ExoSOFTdir,'tools/temp/settings.py'))
            os.remove(os.path.join(ExoSOFTdir,'tools/temp/constants.py'))
            os.remove(os.path.join(ExoSOFTdir,'tools/temp/priors.py'))
        except:
            pass
        shutil.copy(settFilePath,os.path.join(ExoSOFTdir,'tools/temp/settings.py'))
        shutil.copy(os.path.join(os.path.dirname(settFilePath),'priors.py'),os.path.join(ExoSOFTdir,'tools/temp/priors.py'))
        shutil.copy(os.path.join(ExoSOFTdir,'tools/constants.py'),os.path.join(ExoSOFTdir,'tools/temp/constants.py'))
        os.chdir(ExoSOFTdir)
        from tools.temp.settings import settingsDict
        ## load in functions for priors
        #from tools.temp.priors import ePriorRatio,pPriorRatio,incPriorRatio,mass1PriorRatio,mass2PriorRatio,paraPriorRatio,chabrierPrior
        try:
            os.remove(os.path.join(ExoSOFTdir,'tools/temp/settings.py'))
            os.remove(os.path.join(ExoSOFTdir,'tools/temp/constants.py'))
            os.remove(os.path.join(ExoSOFTdir,'tools/temp/priors.py'))
        except:
            pass
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
    #for key in settingsDict:
    #    print key+' = '+repr(settingsDict[key])
    #sys.exit('shirt')
    return settingsDict

def loadFits(filename):
    """
    Load in a fits file written by ExoSOFT.
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

def writeFits(baseFilename,data,settingsDict):
    """
    Data will be written to a fits file with a single PrimaryHDU,
    with the .header loaded up with the tuples from the settingsDict 
    and .data = provided data.
    File will be stored in the 'finalFolder' directory from the settingsDict.
    If data variable is a string, this function will assume it is a filename 
    of where the data is stored in a .npy file, and load it in.
    """
    outFname=''
    try:
        ##check if data is a .npy filename
        if type(data)==str:
            if os.path.exists(data):
                dataFname = data
                data = np.load(dataFname)
                os.remove(dataFname)
                log.debug("just removed data file from disk:\n"+dataFname)
        if len(data)>0:
            if '.fits' not in baseFilename:
                baseFilename=baseFilename+'.fits'
            outFname = os.path.join(settingsDict['finalFolder'],baseFilename)
            hdu = pyfits.PrimaryHDU(data)
            hdulist = pyfits.HDUList([hdu])
            header = hdulist[0].header
            ##load up header with tuples from settingsDict
            for key in settingsDict:
                if type(settingsDict[key])==tuple:
                    header[key]=settingsDict[key][0]
                    if len(settingsDict[key][1])>47:
                        log.warning("comment too long for pyfits headers:"+settingsDict[key][1])
                    else:
                        header.comments[key] = settingsDict[key][1]
                        #print key+' = '+repr((header[key],header.comments[key]))
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
        else:
            log.error("No data to write to file:\n"+baseFilename)
    except:
        log.error("could not write file to disk for some reason")
    return outFname

def periodicDataDump(filename,d):
    """
    dump a ndarray to disk.  If first time, just dump it.
    Else, load current ary and cat d to it before dumping.
    """
    if len(d)!=0:
        if os.path.exists(filename):
            d0 = np.load(filename)
            np.save(filename,np.concatenate((d0,d)))
        else:
            np.save(filename,d)

def combineFits(filenames,outFname):
    """
    combine the data in multiple ExoSOFT fits files together.
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
        
def renameFits(filenameIn,filenameOut,killInput=True,overwrite=True):
    """
    Load input into a new fits HDU, write to output name, close, rm input file.
    """
    goodToGo = True
    if os.path.exists(filenameOut):
        if overwrite:
            rmFiles([filenameOut])
        else:
            goodToGo=False
            log.critical("File already exists, either set overwrite==True, or there be more careful.")
    if goodToGo:
        (head,data) = loadFits(filenameIn)
        hdu = pyfits.PrimaryHDU(data)
        hdulist = pyfits.HDUList([hdu])
        header = hdulist[0].header
        for key in head:
            header[key] = (head[key],head.comments[key])
        hdulist.writeto(filenameOut)
        hdulist.close()
        log.info("output file written to:below\n"+filenameOut)
        if killInput:
            rmFiles([filenameIn])
    
def rmFiles(files):
    ##try to delete files
    for fname in files:
        try:
            if os.path.exists(fname):
                log.debug('Deleting file: '+os.path.basename(fname))
                os.remove(fname) 
        except:
            log.error('Failed to delete file: '+os.path.basename(fname))
    
    
def writeBestsFile(settingsDict,pars,sigs,bstChiSqr,stage):
    
    filename = os.path.join(settingsDict['finalFolder'],'best'+stage+'paramsAndSigs.txt')
    f = open(filename,'w')
    f.write("Best-fit between all "+stage+" chains had a reduce chi squared of "+str(bstChiSqr)+'\n')
    f.write("\nIts parameters were:\n")
    f.write(genTools.nparyTolistStr(pars)+'\n')
    f.write("\n\nIts sigmas were:\n")
    #double check clean up sigs of pars that were not varying
    paramInts = settingsDict['paramInts']
    if stage=='MC':
        sigs = np.zeros(len(pars))
    else:
        for i in range(0,len(pars)):
            if i not in paramInts:
                sigs[i] = 0
    f.write(genTools.nparyTolistStr(sigs)+'\n')
    f.close()
    log.info("Best fit params and sigmas from "+stage+" stage were written to :\n"+filename)  
    
def pushIntoOrigSettFiles(settingsDict,pars,sigs=[]):
    """
    Push the latest found parameters into startParams if better, also replace sigmas?
    Set 'pushToSettFiles'=False to avoid this.
    """
    if settingsDict["pushToSettFiles"]==True:
        origSettFileRoot = settingsDict['origSettFileRoot']
        simpFname = origSettFileRoot+"settingsSimple.py"
        advFname = origSettFileRoot+"settingsAdvanced.py" 
        f = open(simpFname,'r')
        lines = f.readlines()
        f.close()
        linesOut =[]
        replacePars = False
        #go through lines and replace the lines if new vals better than old.
        if settingsDict['startParams'][11]>float(str(pars[11])):
            replacePars = True
            log.debug(repr(settingsDict['startParams'][11])+' > '+str(pars[11]))
        if replacePars==False:
            log.info("Best pars were not better than those in the current settings file.")
        for line in lines:
            key = line.split(':')[0].replace(' ','').replace('"','').replace("'",'')
            if (key=='startParams')and replacePars:
                linesOut.append("'"+key+"'"+" : "+genTools.nparyTolistStr(pars)+',\n')
                log.info('Updated original *settingsSimple.py startParams value to latest best-fit values.')
            elif (key=='startSigmas')and(len(sigs)>1):
                linesOut.append("'"+key+"'"+" : "+genTools.nparyTolistStr(sigs)+',\n')
                log.info('Updated original *settingsSimple.py startSigmas value to latest best-fit values.')
            else:
                linesOut.append(line)
        #Write updated lines to file by the same name, ie replace old version.
        f = open(simpFname,'w')
        f.writelines(linesOut)
        f.close()
        log.info('just overwrote an updated set of values to settings dict:\n'+simpFname+'\n')
        
    
    
    
    
    
    
    
#END OF FILE
    
      