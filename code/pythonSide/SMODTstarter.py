#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import os
import sys
import shutil
import toolboxes as tools
from SimulationManager import simulator
import dicts

def main():
    """
    This is the 'main' that starts a SMODT simulation based on the settings
    in the input settings file.
    
    The simulation is started by calling:
    $ python SMODTstarter.py
    from the directory it is contained in (ie. SMODT/code/pythonSide/)
    """

    # Hard coded paths
    pythonCodeDir = dicts.hardcodedSettingsDict['pythonCodeDir']
    cppCodeDir = dicts.hardcodedSettingsDict['cppCodeDir']
    settings_and_InputDataDir = dicts.hardcodedSettingsDict['settings_and_InputDataDir']
    
    ###########################################################################
    ## VERY IMPORTANT: -Settings file must be called '*SimSettings.txt'
    ##                 -RV data file must be called '*RVdata.dat'
    ##                 -DI data file must be called '*DIdata.dat'
    ##                 -Binary system data file must be called '*SystemData.txt'
    ## Where * represents the string input as the second argument in the terminal 
    ## call, such as:
    ##     $ python SMODTstarter.py HR8799a_
    ## with the file names HR8799a_SimSettings.txt, HR8799a_RVdata.dat...
    ## The default is '' if nothing is entered as the second argument.
    ###########################################################################
    
    print '\n'+'*'*75+'\n'+'*'*22+'  Starting Orbit Simulation   '+'*'*23+'\n'+'*'*75+'\n'
    verboseInternal = False
    if verboseInternal:
        print '*'*50
        print 'About to call tools.gen.cFileToSimSettingsDict to parse settings file and load up paramSettingsDict'
        print '\n'+'*'*50+'\n'
        for i in range(len(sys.argv)):
            print 'arg '+str(i)+' = '+repr(sys.argv[i])
        print '\n'+'*'*50+'\n'
    #Pull in settings filename prepend from command line args, if provided
    settingsNamePrependStr = ''
    if len(sys.argv)>1:
        try:
            settingsNamePrependStr = sys.argv[1]
        except:
            print '\nWarning: the settings file prepended feature is not working correctly !!\n'
    if len(sys.argv)>2:
        print '\nWarning: SMODTstarter only handles a single extra string argument from the command line, and '\
                +str(len(sys.argv)-1)+' were provided !!\n'
                
    # get the input settings file name with dir
    inputSettingsFile = os.path.join(settings_and_InputDataDir,settingsNamePrependStr+'SimSettings.txt')
    # make an output settings file name for the duo version output
    [fName,ext] = os.path.splitext(inputSettingsFile)
    outputSettingsFile = fName+'_DuoVersion'+ext
    paramSettingsDict = tools.gen.cFileToSimSettingsDict(inputSettingsFile, outputSettingsFile, settingsNamePrependStr)
    paramSettingsDict["UpdatedSettingsFile"] = outputSettingsFile
    if verboseInternal:
        print 'DONE loading up paramSettingsDict'
        print '*'*50
    #print 'Repr of paramSettingsDict:'
    #print repr(paramSettingsDict)
    # set simulation type
    mcONLY = paramSettingsDict['mcONLY'] # using only standard 'shotgun' monte carlo?
    
    ## set simulation general settings  
    numSamples = paramSettingsDict['numSamples']  
    #silent = paramSettingsDict['silent']
    #verbose = paramSettingsDict['verbose']
    data_dir = paramSettingsDict['outputData_dir']
    filenameRoot = paramSettingsDict['outputData_filenameRoot']  
    #CalcBurnIn = paramSettingsDict['CalcBurnIn']
    #calcCorrLengths = paramSettingsDict['calcCorrLengths']
    DIonly = paramSettingsDict['DIonly']
    RVonly = paramSettingsDict['RVonly']
    paramSettingsDict['pythonCodeDir'] = pythonCodeDir
    paramSettingsDict['cppCodeDir'] = cppCodeDir
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
    numSamplesString = tools.gen.totalSamplesStr(numSamplesTOTAL)
    
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
        if paramSettingsDict['SILENT']==False:
            print '\n'+'$'*100+'\n'
            print 'WARNING!! the folder:'+finalFolder+', all ready exists!'
        if paramSettingsDict['SILENT']==False:
            print 'You can overwrite the data in it, or exit this simulation.'
        YN = raw_input('\nOVERWRITE current folder (y/n):')
        if (('y' in YN) or ('Y' in YN)):
            if paramSettingsDict['SILENT']==False:
                print '\nDELETING all contents of folder:'+finalFolder
            shutil.rmtree(finalFolder)
            if paramSettingsDict['SILENT']==False:
                print 'MAKING new empty folder:'+finalFolder
            os.mkdir(finalFolder)
        elif (('n' in YN) or ('N' in YN)):
            sys.exit()
        if paramSettingsDict['SILENT']==False:
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
    pythonFiles = ['SMODTstarter.py','dicts/HARDCODEDsettingsDict.py','SimulationManager.py']
    pythonToolboxFiles = ['generalToolbox.py','DItoolbox.py','plotToolbox.py','RVtoolbox.py']
    cppFiles = []
    cppToolboxFiles = ['SimSettingsObj.cpp','RVtools.cpp', 'generalTools.h',
                       'generalTools.cpp','DItools.cpp','DataObj.cpp','DItools.h',
                       'DataObj.h','orbToolboxes.h','RVtools.h','SimSettingsObj.h']
    
    settingsFiles = [paramSettingsDict['SystemDataFilename'],paramSettingsDict["UpdatedSettingsFile"]]
    
    if DIonly:
        settingsFiles.append(paramSettingsDict['DIdataFilename'])
    if RVonly:
        settingsFiles.append(paramSettingsDict['RVdataFilename'])
    if ((DIonly==False)and(RVonly==False)):
        settingsFiles.append(paramSettingsDict['DIdataFilename'])
        settingsFiles.append(paramSettingsDict['RVdataFilename'])
    
    if mcONLY:
        cppFiles.append('mcONLYorbSimulator.cpp')
    else:
        cppFiles.append('MCMCorbSimulator.cpp') 
        cppFiles.append('MCMCorbSimFunc.cpp')  
        cppFiles.append('MCMCorbSimFunc.h') 
        cppFiles.append('simAnnealOrbSimFunc.cpp') 
        cppFiles.append('simAnnealOrbSimFunc.h')    
    
    # add directories and add all complete filenames to overall list
    allFiles = []
    for f in pythonFiles:
        allFiles.append(os.path.join(pythonCodeDir , f))
    for f in pythonToolboxFiles:
        toolboxDir = os.path.join(pythonCodeDir,'toolboxes')
        fullPathUse = os.path.join(toolboxDir , f)
        allFiles.append(fullPathUse)
    for f in cppToolboxFiles:
        cppToolboxDir = os.path.join(cppCodeDir,'toolboxes')
        fullPathUse = os.path.join(cppToolboxDir,f)
        allFiles.append(fullPathUse)
    for f in cppFiles:
        allFiles.append(os.path.join(cppCodeDir , f))
    for f in settingsFiles:
        allFiles.append(os.path.join(settings_and_InputDataDir , f))
        
    # move all final files to code-used directory in output data folder
    for f in allFiles:
        #print '\nCopying code files to output data directory.'
        filebasename = os.path.basename(f)
        outputfilename = finalFolder+'code-used/'+filebasename 
        shutil.copy(f,outputfilename)
        if paramSettingsDict['CopyToDrobox']:
            outputfilename2 = finalFolder2+'code-used/'+filebasename  
            #print 'file being copied:'+f
            #print 'Being copied to:'+outputfilename
            shutil.copy(f,outputfilename2)
        
    ## finally, make print that sim is starting and then start it.
    print '\n** Starting '+str(numProcesses)+' processes of '+numSamplesString+' samples each ** \n'
#     if mcONLY:
#         mcSimStarter(paramSettingsDict)
#     else:
    simulator(paramSettingsDict)

if __name__ == '__main__':
    main()  
