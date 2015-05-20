#import numpy as np
import tools
import simulator
import sys
import os
import shutil

"""
    This is the 'main' of SMODT. 
    It will start things off, call the appropriate set of 
    simulation and post-processing steps.
"""

def smodt():
    """
    'main'
    """
    ###################################################################################
    ######################## START OF START UP STUFF ################################## $$$ Keep here or what???? $$$
    ###################################################################################
    #Pull in settings filename prepend from command line args, if provided
    prepend = ''
    if len(sys.argv)>1:
        try:
            prepend = sys.argv[1]
        except:
            print '\nWarning: the settings file prepended feature is not working correctly !!\n'    
    tempRoot = '/run/media/kmede/Data1/Todai_Work/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/settings_and_inputData/'+prepend###$$$$ this will be handled with setup.py??? How to know where SMODT is on disk??
    settingsDict = tools.loadSettingsDict(tempRoot)
    settingsDict['settingsDir']='/run/media/kmede/Data1/Todai_Work/Dropbox/EclipseWorkspaceDB/SMODT/smodt2/settings_and_inputData/'###$$$$ this will be handled with setup.py??? How to know where SMODT is on disk??
    settingsDict['prepend']=prepend
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
    ##################################################################################
    ######################## END OF START UP STUFF ###################################
    ##################################################################################
        
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=100)
    log.info("Prepend string passed in was '"+prepend+"'")
    Sim = simulator.Simulator(settingsDict)
    Sim.monteCarlo()#$$doesn't do anything yet!!
    Sim.simAnneal()#$$doesn't do anything yet!!
    Sim.mcmc()#$$doesn't do anything yet!!
   
    
    
    
    log.info("End of SMODT2.0 main")

if __name__ == '__main__':
    smodt()