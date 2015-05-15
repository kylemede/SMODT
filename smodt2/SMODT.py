#import numpy as np
import tools
import simulator

"""
    This is the 'main' of SMODT. 
    It will start things off, call the appropriate set of 
    simulation and post-processing steps.
"""

def smodt():
    """
    'main'
    """
    log = tools.getLogger('main',lvl=100)
    log.info("just entered smodt main")
    settingsDict = dict()#$$$$$$$$$$$$$$$ load this up or get it with a tool func
    Sim = simulator.Simulator(settingsDict)
    Sim.monteCarlo()#$$doesn't do anything yet!!
    Sim.simAnneal()#$$doesn't do anything yet!!
    Sim.mcmc()#$$doesn't do anything yet!!
    log.info("End of smodt main")
    



if __name__ == '__main__':
    smodt()