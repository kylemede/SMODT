from mcmcProcessManagerNEW2 import mcmcSimStarter
from mcONLYProcessManagerNEW import mcSimStarter
from paramSettingsDict import paramSettingsDict
import os
import shutil
from orbitToolbox import dictToFile

# set simulation type
mcONLY = paramSettingsDict['General']['mcONLY'] # using only standard 'shotgun' monte carlo?
TH_I = paramSettingsDict['General']['TH_I'] # using Thiele-Innes version?

# set simulation general settings
numSamples = paramSettingsDict['General']['numSamples']  
numProcesses = paramSettingsDict['General']['numProcesses']
temperature = paramSettingsDict['General']['temperature']  
sigmaPercent = paramSettingsDict['General']['sigmaPercent']  
startParams = paramSettingsDict['General']['startParams']  
silent = paramSettingsDict['General']['silent']
verbose = paramSettingsDict['General']['verbose']
dataTempDir = paramSettingsDict['General']['dataTempDir']
dataFinalDir = paramSettingsDict['General']['dataFinalDir']
stripBurnIn = paramSettingsDict['General']['stripBurnIn']
calcCorrLengths = paramSettingsDict['General']['calcCorrLengths']
codeDir = paramSettingsDict['General']['codeDir']
# convert number of samples into easy to read string for file naming
if numSamples<int(1e6):
    numSamplesString = str(numSamples/int(1e3))+'-thousand'
else:
    numSamplesString = str(numSamples/int(1e6))+'-million'
#roundString = 'Round-1'
filenameRoot = paramSettingsDict['General']['filenameRoot']  
filename = filenameRoot+numSamplesString+'.txt'#+'-'+roundString

## Make a directory (folder) to place all the files from this simulation run
# Ensure directories provided have a final '/' character
if dataFinalDir[-1:] is not '/':
     dataFinalDir = dataFinalDir+'/'
finalFolder = dataFinalDir+filename[:-4]+'/'
os.mkdir(finalFolder)
# write settings dict values to finalFolder for reference
#dictToFile(paramSettingsDict,finalFolder+'settingsDictUsed.txt')

#codeDir = '/media/Data1-500GB-NTFS/Todai_Work/workspace/Binary-project/code/'
os.mkdir(finalFolder+'code-used/')
# figure out what files are used in this simulation run
files = ['orbitToolbox.py','mcmcOrbSimStarterNEW.py','paramSettingsDict.py']
if mcONLY:
    files.append('mcONLYProcessManagerNEW.py')
    if TH_I:
        files.append('mcONLY3DorbSimUniformTH_I.py')
    else:
        files.append('mcONLY3DorbSimUniform.py')
else:
    files.append('mcmcProcessManagerNEW2.py')
    if TH_I:
        files.append('mcmc3DorbSimUniformTH_I.py')
    else:
        files.append('mcmc3DorbSimUniform.py')
        
for file in files:
    shutil.copy(codeDir+file,finalFolder+'code-used/'+file)
        
## finally, make print that sim is starting and then start it.
print '\n** Starting '+str(numProcesses)+' processes of '+numSamplesString+' samples each ** \n'
if mcONLY:
    mcSimStarter(filename, numSamples, numProcesses, \
                    silent=silent, verbose=verbose, dataTempDir=dataTempDir, dataFinalDir=dataFinalDir, TH_I=TH_I)
else:
    mcmcSimStarter(filename, numSamples, numProcesses, temperature=temperature, sigmaPercent=sigmaPercent, startParams=startParams,\
                    silent=silent, verbose=verbose, stripBurnIn=stripBurnIn, calcCorrLengths=calcCorrLengths, dataTempDir=dataTempDir, \
                    dataFinalDir=dataFinalDir, TH_I=TH_I)

