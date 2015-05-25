#import numpy as np
import tools
import simulator
import sys
import numpy as np

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
    e = 0.4
    Sys_Dist_PC = 5.0
    Mass1 = 1.0
    Mass2 = 0.2
    Omega = 60.0
    omega = 110.0
    T = 2457000.0
    T_center = 2457000.0
    P = 15.0
    inc =  30.0
    offset = 0.0
    params = np.array([Mass1,Mass2,Sys_Dist_PC,Omega,e,T,T_center,P,inc,omega,0,0,0,offset])
    ##mcONLYcall
    Sim.simulatorFunc('MC')
    if False:
        ## SA call
        (params,sigmas) = Sim.simulatorFunc('SA',params,sigmas)
        ## ST call
        (params,sigmas) = Sim.simulatorFunc('ST',params,sigmas)
        ##MCMC call
        Sim.simulatorFunc('MCMC',params,sigmas)
    #Sim.monteCarlo()#$$doesn't do anything yet!!
    #Sim.simAnneal()#$$doesn't do anything yet!!
    #Sim.mcmc()#$$doesn't do anything yet!!
   
    
    
    
    log.info("End of SMODT2.0 main")

if __name__ == '__main__':
    smodt()