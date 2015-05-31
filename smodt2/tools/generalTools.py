#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
#import numpy as np
import smodtLogger
import os
import shutil
import numpy as np
import sys
import pyfits

log = smodtLogger.getLogger('main.genTools',lvl=100,addFH=False)

def test():
    log.info("inside the tools test func")
    
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
    for line in lines:
        if line=='\n':
            jitterLast = 0
            datasetNumLast+=1
        if len(line.split())>2:
            if line.split()[0].replace('.','',1).isdigit() and line.split()[1].replace('.','',1).replace('-','',1).isdigit():
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
                rvData.append(curDataAry)
                
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
        if False:
            print 'using diFilename = '+diFilename        
        if os.path.exists(diFilename):
            diData = loadDIdata(diFilename)
            diEpochs = diData[:,0]
    if dataMode!='DI':
        rvFilename = filenameRoot+'RVdata.dat'
        if False:
            print 'using rvFilename = '+rvFilename
        if os.path.exists(rvFilename):
            rvData = loadRVdata(rvFilename)
            rvEpochs = rvData[:,0]
            
    #load in epochs from both sets, sort and kill double entries
    epochsTemp = np.concatenate((diEpochs,rvEpochs))
    epochsTemp.sort()
    epochs = []
    for epoch in epochsTemp:
        if epoch not in epochs:
            epochs.append(epoch)
    epochs = np.array(epochs)
    realData = np.zeros((epochs.shape[0],8))
    realData[:,2]=realData[:,4]=realData[:,6]=np.inf
    realData[:,0]=epochs[:]
    diCounter = 0
    rvCounter = 0
    for i in range(epochs.shape[0]):
        if len(diEpochs)>0:
            if epochs[i]==diData[diCounter,0]:
                realData[i,1:5]=diData[diCounter,1:]
                diCounter+=1
        if len(rvEpochs)>0:
            if epochs[i]==rvData[rvCounter,0]:
                realData[i,5:]=rvData[rvCounter,1:]
                rvCounter+=1
    #print 'dataMode'+dataMode+'->realData = '+repr(realData)
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
    smodtHeadDir = filenameRoot.split("SMODT")[0]
    try:
        os.remove(os.path.join(smodtHeadDir,'SMODT/smodt2/tools/temp/settingsSimple.py'))
        os.remove(os.path.join(smodtHeadDir,'SMODT/smodt2/tools/temp/settingsAdvanced.py'))
        os.remove(os.path.join(smodtHeadDir,'SMODT/smodt2/tools/temp/constants.py'))
    except:
        temp=True
    shutil.copy(filenameRoot+'settingsSimple.py',os.path.join(smodtHeadDir,'SMODT/smodt2/tools/temp/settingsSimple.py'))
    shutil.copy(filenameRoot+'settingsAdvanced.py',os.path.join(smodtHeadDir,'SMODT/smodt2/tools/temp/settingsAdvanced.py'))
    shutil.copy(os.path.join(smodtHeadDir,'SMODT/smodt2/tools/constants.py'),os.path.join(smodtHeadDir,'SMODT/smodt2/tools/temp/constants.py'))
    if False:
        print 'Copied simple to:\n'+os.path.join(smodtHeadDir,'SMODT/smodt2/tools/temp/settingsSimple.py')
        print 'Copied advanced to:\n'+os.path.join(smodtHeadDir,'SMODT/smodt2/tools/temp/settingsAdvanced.py')
        print 'Copied constants to:\n'+os.path.join(smodtHeadDir,'SMODT/smodt2/tools/temp/constants.py')
    os.chdir(os.path.join(smodtHeadDir,'SMODT'))
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
    
    return settingsDict
    
def startup(argv):    
    ## Pull in settings filename prepend from command line args, if provided
    prepend = ''
    if len(argv)>1:
        try:
            prepend = argv[1]
        except:
            print '\nWarning: the settings file prepended feature is not working correctly !!\n'    
    ######################### FIND MORE ELEGANT WAY TO DO THIS!!!!############################
    tempRoot = '/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/'###$$$$ this will be handled with setup.py??? How to know where SMODT is on disk??
    settingsDict = loadSettingsDict(tempRoot+'settings_and_inputData/'+prepend)
    settingsDict['smodtdir']=tempRoot
    settingsDict['settingsDir']=os.path.join(settingsDict['smodtdir'],'settings_and_inputData/')
    settingsDict['prepend']=prepend
    #####################################################################################
    
    ## Make a directory (folder) to place all the files from this simulation run
    settingsDict['finalFolder'] = os.path.join(settingsDict['outDir'],settingsDict['outRoot'])
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
        #print 'moved back to:\n'+cwd
    return settingsDict
        
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
                print key+' = '+repr(header[key])
                print 'type(header[key] = '+repr(type(header[key]))
        print '\n\n'+repr(head)
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
    (head0,dataALL) = loadFits(filenames[0])
    for filename in filenames:
        (head,data) = loadFits(filename)
        dataALL = np.concatenate((dataALL,data))
    hdu = pyfits.PrimaryHDU(dataALL)
    hdulist = pyfits.HDUList([hdu])
    header = hdulist[0].header
    for key in head0:
        header[key]=head0[key]
    hdulist.writeto(outFname)
    hdulist.close()
    log.info("output file written to:below\n"+outFname)
    
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
        
def confLevelFinder(filename, colNum=False, returnData=False, returnChiSquareds=False, returnBestDataVal=False,fast=True):
    """
    A function to find the 68.3 and 95.4% confidence levels in a given output data file's column.
    
    return [[68.3% minimum, 68.3% maximum],[95.5% minimum, 95.5% maximum]]
    
    columnNum must be an int.    
    """
    verboseInternal = False
    bestCentered = False
    log.debug('Inside confLevelFinder')
    if os.path.exists(filename):
        (dataAry,chiSquareds,[bestDataVal,dataMedian,dataValueStart,dataValueMid,dataValueEnd]) = dataReader(filename, colNum)
    
        if len(dataAry>0):
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
            
        s= "Final 68% range values are: "+repr(conf68Vals)+'\n'
        s=s+"Final 95% range values are: "+repr(conf95Vals)+'\n'
        if bestCentered:
            s=s+ "\nerror is centered on best \n"
            s=s+"68.3% error level = "+str(bestDataVal-conf68Vals[0])+'\n'
            s=s+" =>   "+str(dataMedian)+'  +/-  '+str(bestDataVal-conf68Vals[0])+'\n'
        else:
            s=s+ "\nerror is centered on Median \n"
            s=s+"68.3% error level = "+str(dataMedian-conf68Vals[0])
            s=s+" =>   "+str(dataMedian)+'  +/-  '+str(dataMedian-conf68Vals[0])+'\n'
        s=s+'\n'+75*'-'+'\n Leaving confLevelFinder \n'+75*'-'+'\n'
        log.debug(s)
        
        if verboseInternal:
            print 'returnData = '+repr(returnData)+', returnChiSquareds = '+repr(returnChiSquareds)+', returnBestDataVal = '+repr(returnBestDataVal)
        
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
            returnList =   ([conf68Vals,conf95Vals],dataAry, bestDataVal)
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
    (head,data) = loadFits(filename)
    chiBest = np.min(data[:,11])
    loc = np.where(data[:,11]==chiBest)
    orbBest = data[loc[0][0],:]
    log.info("Best fit found to be:\n"+repr(orbBest))
    return orbBest
                                         
                                         
                                         
    