#import numpy as np
import tools
import simulator
import sys

"""
    This is the 'main' of SMODT. 
    It will start things off, call the appropriate set of 
    simulation and post-processing steps.
"""

def smodt():
    """
    'main'
    """
    settingsDict = tools.startup(sys.argv)
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=100)
    log.info("Prepend string passed in was '"+settingsDict['prepend']+"'")
    Sim = simulator.Simulator(settingsDict)
    Sim.monteCarlo()#$$doesn't do anything yet!!
    Sim.simAnneal()#$$doesn't do anything yet!!
    Sim.mcmc()#$$doesn't do anything yet!!
   
    
    
    
    log.info("End of SMODT2.0 main")

if __name__ == '__main__':
    smodt()