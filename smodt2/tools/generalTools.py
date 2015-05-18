#import numpy as np
import smodtLogger
import os
import numpy as np

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
            if line.split()[0].replace('.','',1).isdigit() and line.split()[2].replace('.','',1).isdigit():
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
            if line.split()[0].replace('.','',1).isdigit() and line.split()[2].replace('.','',1).isdigit():
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
    
    
def loadRealData(filenameRoot):
    """
    Load the observed real data into a numpy array.
    This will be a combination of the RV and DI data,sorted into cronological order.
    filenameRoot would be the absolute path plus the prepend to the settings files.
    ex. '/run/..../SMODT/settings_and_inputData/FakeData_'
    """
    diFilename = filenameRoot+'DIdata.dat'
    rvFilename = filenameRoot+'RVdata.dat'
    diEpochs = []
    rvEpochs = []
    if os.path.exists(diFilename):
        diData = loadDIdata(diFilename)
        diEpochs = diData[:,0]
        #print 'diData = '+repr(diData)
    if os.path.exists(rvFilename):
        rvData = loadRVdata(rvFilename)
        rvEpochs = rvData[:,0]
        #print 'rvData = '+repr(rvData)
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
    print 'realData = '+repr(realData)
    return realData
            
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    