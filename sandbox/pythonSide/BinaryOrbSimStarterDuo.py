#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import os
import sys
import shutil
import toolboxes as tools
from MCMC_ProcessManagerDuo import mcmcSimStarter
#from mcONLY_ProcessManagerDuo import mcSimStarter
import dicts

# Hard coded paths
pythonCodeDir = dicts.hardcodedSettingsDict['pythonCodeDir']
cplusplusCodeDir = dicts.hardcodedSettingsDict['cplusplusCodeDir']
settings_and_InputDataDir = dicts.hardcodedSettingsDict['settings_and_InputDataDir']

###########################################################################
## VERY IMPORTANT: -Settings file must be called 'SimSettings.txt'
##                 -RV data file must be called 'RVdata.dat'
##                 -DI data file must be called 'DIdata.dat'
##                 -Binary system data file must be called 'SystemData.txt'
###########################################################################

print '\n\n\n'+'='*100+'\n'+'='*100+'\n'+'='*35+'  Starting Orbit Simulation   '+'='*35+'\n'+'='*100+'\n'+'='*100+'\n\n\n'
verboseInternal = False
if verboseInternal:
    print '*'*50
    print 'About to call tools.gen.cFileToSimSettingsDict to parse settings file and load up paramSettingsDict'
# get the input settings file name with dir
inputSettingsFile = os.path.join(settings_and_InputDataDir,'SimSettings.txt')
# make an output settings file name for the duo version output
[fName,ext] = os.path.splitext(inputSettingsFile)
outputSettingsFile = fName+'_DuoVersion'+ext
paramSettingsDict = tools.gen.cFileToSimSettingsDict(inputSettingsFile, outputSettingsFile)
paramSettingsDict["UpdatedSettingsFile"] = outputSettingsFile
if verboseInternal:
    print 'DONE loading up paramSettingsDict'
    print '*'*50
#print 'Repr of paramSettingsDict:'
#print repr(paramSettingsDict)
# set simulation type
mcONLY = paramSettingsDict['mcONLY'] # using only standard 'shotgun' monte carlo?

# set simulation general settings  
numSamples = paramSettingsDict['numSamples']  
silent = paramSettingsDict['silent']
verbose = paramSettingsDict['verbose']
data_dir = paramSettingsDict['outputData_dir']
filenameRoot = paramSettingsDict['outputData_filenameRoot']  
CalcBurnIn = paramSettingsDict['CalcBurnIn']
calcCorrLengths = paramSettingsDict['calcCorrLengths']
DIonly = paramSettingsDict['DIonly']
RVonly = paramSettingsDict['RVonly']
paramSettingsDict['pythonCodeDir'] = pythonCodeDir
paramSettingsDict['cplusplusCodeDir'] = cplusplusCodeDir
paramSettingsDict['settings_and_InputDataDir'] = settings_and_InputDataDir

if paramSettingsDict['useMultiProcessing']:
    ## find the number of available cpus
    import multiprocessing
    numCores = multiprocessing.cpu_count()
    # set to max-1
    numProcesses = numCores-1
else:
    numProcesses = 1
    
paramSettingsDict['numProcesses'] = numProcesses 

# convert number of samples into easy to read string for file naming
numSamplesTOTAL = numSamples*numProcesses
if numSamplesTOTAL>=int(1e9):
    numSamplesString = str(int(numSamplesTOTAL/int(1e9)))+'-Billion-in_Total'
elif numSamplesTOTAL>=int(1e6):
    numSamplesString = str(int(numSamplesTOTAL/int(1e6)))+'-Million-in_Total'
elif numSamplesTOTAL>=int(1e3):
    numSamplesString = str(int(numSamplesTOTAL/int(1e3)))+'-Thousand-in_Total'
else:
    numSamplesString = str(int(numSamplesTOTAL))+'-in_Total'

# add num sample string to filename
if (filenameRoot[-1]!='-')and(filenameRoot[-1]!='_'):
    filenameRootModed = filenameRoot+'--'
else:
    filenameRootModed = filenameRoot
filename = filenameRootModed+numSamplesString+'.txt'

## Make a directory (folder) to place all the files from this simulation run
# Ensure directories provided have a final '/' character
if data_dir[-1:] is not '/':
     data_dir = data_dir+'/'
finalFolder = data_dir+filename[:-4]+'/'
if os.path.exists(finalFolder):
    print '\n'+'$'*100+'\n'
    print 'WARNING!! the folder:'+finalFolder+', all ready exists!'
    print 'You can overwrite the data in it, or exit this simulation.'
    YN = raw_input('\nOVERWRITE current folder (y/n):')
    if (('y' in YN) or ('Y' in YN)):
        print '\nDELETING all contents of folder:'+finalFolder
        shutil.rmtree(finalFolder)
        print 'MAKING new empty folder:'+finalFolder
        os.mkdir(finalFolder)
    elif (('n' in YN) or ('N' in YN)):
        sys.exit()
    print '\n'+'$'*100+'\n'
else:
    os.mkdir(finalFolder)
    
paramSettingsDict['outputData_dir'] = finalFolder

# make a directory inside the data folder for the code used during this simulation
os.mkdir(finalFolder+'code-used/')

if paramSettingsDict['CopyToDrobox']:
    finalFolder2 = os.path.join('/mnt/Data1/Todai_Work/Dropbox/SMODT-outputCopies/',filename[:-4]+'/')
    if os.path.exists(finalFolder2):
        print '\n'+'$'*100+'\n'
        print 'WARNING!! the folder:'+finalFolder2+', all ready exists!'
        print 'You can overwrite the data in it, or exit this simulation.'
        YN = raw_input('\nOVERWRITE current folder (y/n):')
        if (('y' in YN) or ('Y' in YN)):
            print '\nDELETING all contents of folder:'+finalFolder2
            shutil.rmtree(finalFolder2)
            print 'MAKING new empty folder:'+finalFolder2
            os.mkdir(finalFolder2)
        elif (('n' in YN) or ('N' in YN)):
            sys.exit()
        print '\n'+'$'*100+'\n'
    else:
        os.mkdir(finalFolder2)
    # make a directory inside the data folder for the code used during this simulation
    os.mkdir(finalFolder2+'code-used/')


# figure out what files are used in this simulation run
pythonFiles = ['BinaryOrbSimStarterDuo.py','Dicts/HARDCODEDsettingsDict.py']
pythonToolboxFiles = ['orbitToolboxDuo.py','DItoolbox.py','plottingToolbox.py','RVtoolbox.py']

cplusplusFiles = ['Toolboxes/SimSettingsObj.cpp','Toolboxes/RVtools.cpp',
                  'Toolboxes/orbSimStructures.cpp','Toolboxes/OrbSimGeneralToolbox.cpp',
                  'Toolboxes/DItools.cpp','Toolboxes/DataObj.cpp',
                  'Toolboxes/DItools.h','Toolboxes/DataObj.h',
                  'Toolboxes/orbSimStructures.h','Toolboxes/orbToolboxes.h',
                  'Toolboxes/RVtools.h','Toolboxes/SimSettingsObj.h']

settingsFiles = ['SystemData.txt','SimSettings.txt']

if DIonly:
    settingsFiles.append('DIdata.dat')
if RVonly:
    settingsFiles.append('RVdata.dat')
if ((DIonly==False)and(RVonly==False)):
    settingsFiles.append('DIdata.dat')
    settingsFiles.append('RVdata.dat')

if mcONLY:
    pythonFiles.append('mcONLY_ProcessManagerDuo.py')
    cplusplusFiles.append('mcONLYorbSimulator.cpp')
else:
    pythonFiles.append('MCMC_ProcessManagerDuo.py')
    cplusplusFiles.append('MCMCorbSimulator.cpp') 
    cplusplusFiles.append('MCMCorbSimFunc.cpp')  
    cplusplusFiles.append('MCMCorbSimFunc.h') 
    cplusplusFiles.append('simAnnealOrbSimulator.cpp') 
    cplusplusFiles.append('simAnnealOrbSimFunc.cpp') 
    cplusplusFiles.append('simAnnealOrbSimFunc.h') 

# add directories and add all complete filenames to overall list
allFiles = []
for file in pythonFiles:
    allFiles.append(os.path.join(pythonCodeDir , file))
for file in pythonToolboxFiles:
    toolboxDir = os.path.join(pythonCodeDir,'Toolboxes')
    fullpathUse = os.path.join(toolboxDir , file)
    allFiles.append(fullpathUse)
for file in cplusplusFiles:
    allFiles.append(os.path.join(cplusplusCodeDir , file))
for file in settingsFiles:
    allFiles.append(os.path.join(settings_and_InputDataDir , file))
    
# move all final files to code-used directory in output data folder
for file in allFiles:
    #print '\nCopying code files to output data directory.'
    filebasename = os.path.basename(file)
    outputfilename = finalFolder+'code-used/'+filebasename 
    shutil.copy(file,outputfilename)
    if paramSettingsDict['CopyToDrobox']:
        outputfilename2 = finalFolder2+'code-used/'+filebasename  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        print 'file being copied:'+file
        print 'Being copied to:'+outputfilename
        shutil.copy(file,outputfilename2)#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
    
## finally, make print that sim is starting and then start it.
print '\n** Starting '+str(numProcesses)+' processes of '+numSamplesString+' samples each ** \n'
if mcONLY:
    mcSimStarter(paramSettingsDict)
else:
    mcmcSimStarter(paramSettingsDict)

