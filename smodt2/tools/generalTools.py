#import numpy as np
import smodtLogger
import os
import shutil
import numpy as np
import sys
import pyfits
#np.set_printoptions(precision=15)

log = smodtLogger.getLogger('main.tools',lvl=100,addFH=False)

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
        #log.debug("line was:'"+line+"'")#$$$$$$$$$$$$$$$$$$$$$$$$
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
#             else:
#                 print '0th and 1st column vals not numbers, \nline:'+line
#                 print '0th = '+repr(line.split()[0].replace('.','',1).isdigit())
#                 print '1st = '+repr(line.split()[1].replace('.','',1).replace('-','',1).isdigit())
#         else:
#             print "line had less than 2 columns, \nline:"+line
                
    return np.array(rvData)
    
def loadRealData(filenameRoot):
    """
    Load the observed real data into a numpy array.
    This will be a combination of the RV and DI data,sorted into cronological order.
    filenameRoot would be the absolute path plus the prepend to the settings files.
    ex. '/run/..../SMODT/settings_and_inputData/FakeData_'
    """
    diFilename = filenameRoot+'DIdata.dat'
    rvFilename = filenameRoot+'RVdata.dat'
    if False:
        print 'using diFilename = '+diFilename
        print 'using rvFilename = '+rvFilename
    diEpochs = []
    rvEpochs = []
    if os.path.exists(diFilename):
        diData = loadDIdata(diFilename)
        #print 'diData = '+repr(diData)
        diEpochs = diData[:,0]
    if os.path.exists(rvFilename):
        rvData = loadRVdata(rvFilename)
        #print 'rvData = '+repr(rvData)
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
    #print 'epochs = '+repr(epochs)
    diCounter = 0
    rvCounter = 0
    for i in range(epochs.shape[0]):
        if len(diEpochs)>0:
            if epochs[i]==diData[diCounter,0]:
                realData[i,0:5]=diData[diCounter]
                diCounter+=1
        if len(rvEpochs)>0:
            if epochs[i]==rvData[rvCounter,0]:
                realData[i,5:]=rvData[rvCounter,1:]
                rvCounter+=1
    #print 'realData = '+repr(realData)
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
    print '\n'+filenameRoot+'settingsSimple.py\n'
    shutil.copy(filenameRoot+'settingsSimple.py',os.path.join(smodtHeadDir,'SMODT/smodt2/tools/temp/settingsSimple.py'))
    shutil.copy(filenameRoot+'settingsAdvanced.py',os.path.join(smodtHeadDir,'SMODT/smodt2/tools/temp/settingsAdvanced.py'))
    shutil.copy(os.path.dirname(filenameRoot)+'/constants.py',os.path.join(smodtHeadDir,'SMODT/smodt2/tools/temp/constants.py'))
    if False:
        print 'Copied simple to:\n'+os.path.join(smodtHeadDir,'SMODT/smodt2/tools/temp/settingsSimple.py')
        print 'Copied advanced to:\n'+os.path.join(smodtHeadDir,'SMODT/smodt2/tools/temp/settingsAdvanced.py')
        print 'Copied constants to:\n'+os.path.join(smodtHeadDir,'SMODT/smodt2/tools/temp/constants.py')
    os.chdir(os.path.join(smodtHeadDir,'SMODT'))
    from tools.temp.settingsAdvanced import settingsDict
    os.chdir(cwd)
    return settingsDict
    #print "settingsDict['pPrior'](3.0) ="+str(settingsDict['pPrior'](3.0))

    
def startup(argv):    
    ## Pull in settings filename prepend from command line args, if provided
    prepend = ''
    if len(argv)>1:
        try:
            prepend = argv[1]
        except:
            print '\nWarning: the settings file prepended feature is not working correctly !!\n'    
    ######################### FIND MORE ELEGANT WAY TO DO THIS!!!!############################
    tempRoot = '/run/media/kmede/Data1/Todai_Work/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/settings_and_inputData/'+prepend###$$$$ this will be handled with setup.py??? How to know where SMODT is on disk??
    settingsDict = loadSettingsDict(tempRoot)
    settingsDict['smodtdir']='/run/media/kmede/Data1/Todai_Work/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/'###$$$$ this will be handled with setup.py??? How to know where SMODT is on disk??
    settingsDict['settingsDir']=os.path.join(settingsDict['smodtdir'],'settings_and_inputData/')
    settingsDict['prepend']=prepend
    #####################################################################################
    
    ## Make a directory (folder) to place all the files from this simulation run
    settingsDict['finalFolder'] = os.path.join(settingsDict['outDir'],settingsDict['outRoot'])
    if os.path.exists(settingsDict['finalFolder']):
        if settingsDict['SILENT']==False:
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
        if settingsDict['SILENT']==False:
            print '$'*50+'\n'
    else:
        os.mkdir(settingsDict['finalFolder'])
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
    hdu = pyfits.PrimaryHDU(data)
    hdulist = pyfits.HDUList([hdu])
    header = hdulist[0].header
    ##load up header with tuples from settingsDict
    #print "starting to load up header with tuples from settingsDict"
    for key in settingsDict:
        if type(settingsDict[key])==tuple:
            header[key]=settingsDict[key]
            #print "adding '"+key+"' = "+repr(settingsDict[key])
    hdulist.writeto(os.path.join(settingsDict['finalFolder'],baseFilename))
    hdulist.close()
    ## check resulting fits file header
    if False:
        f = pyfits.open(os.path.join(settingsDict['finalFolder'],baseFilename))
        head = f[0].header
        f.close()
        for key in head:
            print key+' = '+repr(header[key])
        print '\n\n'+repr(head)
    
    
    
    
    