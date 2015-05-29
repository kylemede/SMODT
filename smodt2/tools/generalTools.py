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
    tempRoot = '/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/settings_and_inputData/'+prepend###$$$$ this will be handled with setup.py??? How to know where SMODT is on disk??
    settingsDict = loadSettingsDict(tempRoot)
    settingsDict['smodtdir']='/run/media/kmede/HOME/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/'###$$$$ this will be handled with setup.py??? How to know where SMODT is on disk??
    settingsDict['settingsDir']=os.path.join(settingsDict['smodtdir'],'settings_and_inputData/')
    settingsDict['prepend']=prepend
    #####################################################################################
    
    ## Make a directory (folder) to place all the files from this simulation run
    settingsDict['finalFolder'] = os.path.join(settingsDict['outDir'],settingsDict['outRoot'])
    if os.path.exists(settingsDict['finalFolder']):
        if settingsDict['logLevel']<100: ## Handle this with a 'clob' bool in dict??? $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            print '\n'+'$'*50
            print 'WARNING!! the folder:\n"'+settingsDict['finalFolder']+'"\nALREADY EXISTS!'
            print 'You can overwrite the data in it, or exit this simulation.'
            YN = raw_input('OVERWRITE current folder (y/n):')
        else:
            YN = 'y'
        if (('y' in YN) or ('Y' in YN)):
            shutil.rmtree(settingsDict['finalFolder'])
            os.mkdir(settingsDict['finalFolder'])
        elif (('n' in YN) or ('N' in YN)):
            sys.exit()
        if settingsDict['logLevel']<100:
            print '$'*50+'\n'
    else:
        os.mkdir(settingsDict['finalFolder'])
    if False:
        for key in settingsDict:
           print key+' = '+repr(settingsDict[key])
    ## run make for swig if requested
    if settingsDict['remake']:
        cwd = os.getcwd()
        os.chdir(os.path.join(settingsDict['smodtdir'],'tools/cppTools/'))
        os.system('make clean')
        os.system('make')
        os.chdir(cwd)
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
        
def confLevelFinder(filename, columNum=False, returnData=False, returnChiSquareds=False, returnBestDataVal=False,fast=True):
    """
    A function to find the 68.3 and 95.4% confidence levels in a given output data file's column.
    
    return [[68.3% minimum, 68.3% maximum],[95.5% minimum, 95.5% maximum]]
    
    columnNum must be an int.
    
    file format must be:
        
        line1: title
        line2: data headings
        line3: space delimited data
        line4: space delimited data
    
    """
    verboseInternal = False
    bestCentered = False
    log.debug('Inside confLevelFinder')
    if os.path.exists(filename):
        if fast:
            ignoreConstParam = False
        else:
            ignoreConstParam = True
        (log,dataAry,chiSquareds,[bestDataVal,dataMedian,dataValueStart,dataValueMid,dataValueEnd]) = dataReader(filename, columNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True, ignoreConstParam=ignoreConstParam)
    
        if len(dataAry>0):
            #Convert data array to a numpy array
            dataAry = np.sort(dataAry)
            # find range of data's values
            dataMax = dataAry[-1]
            dataMin = dataAry[0]
            size = dataAry.size
            if (size%2)==0:
                dataMedian = (dataAry[size / 2 - 1] + dataAry[size / 2]) / 2
            else:
                dataMedian = dataAry[size / 2]
        
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
        log.info(s)
        
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
        
def dataReader(filename, columNum=False, returnData=False, returnChiSquareds=False, returnBestDataVal=False, ignoreConstParam=False):
    """
    Read in the data for a single column of data.
    """
    verboseInternal = False
    ## First get ranges of param and ChiSquared values
    log.debug('\nOpening and finding ranges for data in column # '+str(columNum))
    
    ## Check if file has useful data for that column#
    (head,data) = loadFits(filename)
    if head!=False:
        TotalSamples=data.shape[0]

        
        
        # find values at start, mid and end of file
        fp = open(filename,'r')
        lastColLoc = 0
        dataValueStart = dataValueMid = dataValueEnd =0
        for i,line in enumerate(fp):
            if i==(0+2):
                splitAry = line.split()
                lastColLoc = len(splitAry)-1
                dataValueStart = float(splitAry[columNum])
                if verboseInternal:
                    print '\nstart = '+str(dataValueStart)+'\n'
            elif i==((numDataLines//2)+2):
                splitAry = line.split()
                dataValueMid = float(splitAry[columNum])
                if verboseInternal:
                    print '\nmid = '+str(dataValueMid)+'\n'
            elif i==numDataLines:
                splitAry = line.split()
                dataValueEnd = float(splitAry[columNum])
                if verboseInternal:
                    print '\nend = '+str(dataValueEnd)+'\n'
        fp.close()
        
        doesntVary = True
        dataAry = []
        chiSquareds = []
        bestOrbit = 0
        bestDataVal = 0
        totalAccepted = 0
        chiSquaredMin=1e6
        dataMax = 0
        dataMin = 1e9
        if ((dataValueStart!=dataValueMid)and(dataValueStart!=dataValueEnd)):
            if gotLog:
                log.write("Values for parameter found to be constant!!")
            if verboseInternal:
                print "Values for parameter found to be constant!!"
            doesntVary = False
            
        if ((doesntVary==True)and(ignoreConstParam==False)):
            if returnData:
                dataAry = [dataValueStart]*TotalSamples
            if returnChiSquareds:
                chiSquareds = [0]*TotalSamples
        elif ((doesntVary==False)or(ignoreConstParam==True)):#or(fast==False):  
            s=''
            fp = open(filename,'r')
            #Old string parsing directly version UPDATED
            startTime2 = timeit.default_timer()
            totalAccepted = 0
            dataAry = [None]*TotalSamples
            if returnChiSquareds:
                chiSquareds = [None]*TotalSamples
            j = 0
            lineNum=0
            numNoDataLines=0
            firstDataLine = ""
            lastDataLine = ""
            firstJ = ""
            lastJ = ""
            for i,line in enumerate(fp):
                lineNum+=1
                if line[0].isdigit():
                    s2 = "?"#$$$$$$$$$$$ DEBUGGING $$$$$$$$$$
                    if firstDataLine=="":
                        firstDataLine=line
                    try:
                        dataLineCols = line.split()
                        ## this should never happen, but it is a check for a double decimal value
                        decimalSplit = dataLineCols[columNum].split('.')
                        if len(decimalSplit)>2:
                            dataValue = float(decimalSplit[0]+'.'+decimalSplit[1])
                        else:
                            dataValue = float(dataLineCols[columNum])
                        #s = "1060"#$$$$$$$$$$$ DEBUGGING $$$$$$$$$$
                        chiSquared = float(dataLineCols[8])
                        #s = "1062"#$$$$$$$$$$$ DEBUGGING $$$$$$$$$$
                        s2 =" chiSquared = "+str(chiSquared)+", dataValue = "+str(dataValue)#$$$$$$$$$$$ DEBUGGING $$$$$$$$$$
                        #if (chiSquared==0)or(dataValue==0):
                            # if verboseInternal:
                            #     print line
                        if firstJ=="":
                            firstJ=j
                        s2+="\nIn itter loop, j="+str(j)+", totalAccepted="+str(totalAccepted)+", len(dataAry)="+str(len(dataAry))
                        if totalAccepted>len(dataAry):
                            print "\n*** totalAccepted>len(dataAry) ***"
                            print s2
                            break
                        else:
                            try:
                                dataAry[j]=dataValue
                            except:
                                print s2
                                print "\nfailed to load data into dataArray"+", chiSquared = "+str(chiSquared)+", dataValue = "+str(dataValue)
                            if returnChiSquareds:
                                try:
                                    chiSquareds[j]=chiSquared
                                except:
                                    print s2
                                    print "\nfailed to load chiSquared into chiSquareds array"+", chiSquared = "+str(chiSquared)+", dataValue = "+str(dataValue)
                            totalAccepted+=1
                            j+=1
                        #s = "1074"#$$$$$$$$$$$ DEBUGGING $$$$$$$$$$
                        if dataValue>dataMax:
                            dataMax = dataValue
                        #s = "1077"#$$$$$$$$$$$ DEBUGGING $$$$$$$$$$
                        if dataValue<dataMin:
                            dataMin = dataValue
                        #s = "1080"#$$$$$$$$$$$ DEBUGGING $$$$$$$$$$
                        if chiSquared<chiSquaredMin:
                            chiSquaredMin = chiSquared
                            bestDataVal = dataValue
                            bestOrbit=lineNum   
                        #s = "1085"#$$$$$$$$$$$ DEBUGGING $$$$$$$$$$   
                    except:
                        print "code line failed = "+s2
                        print 'Failed for line: '+line 
                else:
                    numNoDataLines+=1
            endTime2 = timeit.default_timer()
            totalTime = (endTime2-startTime2) # in seconds
            totalTimeString = timeString(totalTime)
            s=s+ '\nUPDATED Direct data loading took '+totalTimeString+' to complete.\n' 
            s=s+'The resulting arrays had '+str(totalAccepted)+' elements, with best value '+str(bestDataVal)+', and minChiSquared '+str(chiSquaredMin)
            if gotLog:
                log.write(s+'\n')
            if verboseInternal:
                print s+"\n"
            lastDataLine = line
            lastJ = j
            fp.close()
        dataAry = np.array(dataAry)
        dataMedian = np.median(dataAry)
        s=  '\nTotal number of orbits = '+str(totalAccepted)
        if verboseInternal:
            s+=", len(dataAry)="+str(len(dataAry))+", i = "+str(i)+", j = "+str(j)
            s+=", fistJ = "+str(firstJ)+", lastJ = "+str(lastJ)
            s+=", lineNum = "+str(lineNum)+", numDataLines = "+str(numDataLines)+", numNoDataLines = "+str(numNoDataLines)
            s+="\nfirstDataLine = "+firstDataLine+"\nlastDataLine = "+lastDataLine+"\n"
        s=s+'\nBest value found was '+str(bestDataVal)+", at line Number "+str(bestOrbit)+", and had a chiSquared = "+str(chiSquaredMin)
        s=s+'\nMedian value = '+str(dataMedian)
        s=s+'\n[Min,Max] values found for data were '+repr([dataMin,dataMax])
        if gotLog:
            log.write(s+'\n')
        if verboseInternal:
            print s
            print "first and last elements of dataAry are "+str(dataAry[0])+", "+str(dataAry[-1])+"\n"
        
        return (log,dataAry,chiSquareds,[bestDataVal,dataMedian,dataValueStart,dataValueMid,dataValueEnd])
    