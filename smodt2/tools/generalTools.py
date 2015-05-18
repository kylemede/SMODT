#import numpy as np
import smodtLogger

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
    # First load first file into memory
    file = open(filename, 'r')
    diData = []
    while (isnumeric(line.split()[0]) and isnumeric(line.split()[2])):
        diData.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2]),float(line.split()[3]),float(line.split()[4])])
        line = file.readline()
    file.close()    
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
    if filename[-4:]!='.dat':
        filename = filename+'.dat'
    # First load first file into memory
    file = open(filename, 'r')
    rvData = []
    datasetNumLast = 0
    jitterLast = 0
    while (isnumeric(line.split()[0]) and isnumeric(line.split()[2])):
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
        line = file.readline()
        log.debug("line was:'"+line+"'")#$$$$$$$$$$$$$$$$$$$$$$$$
        if line=='':
            jitterLast = 0
            datasetNum+=1
    file.close()    
    return RVData 
    
    
def loadRealData(filenameRoot):
    """
    Load the observed real data into a numpy array.
    This will be a combination of the RV and DI data.
    """
    