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
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=20)
    log.debug("Prepend string passed in was '"+settingsDict['prepend']+"'")
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
    paramsST = np.array( [  9.37145408e-01,   1.98238112e-01,   5.02014226e+00,\
         6.15361129e+01,   4.15495712e-01,   2.45700439e+06,\
         2.45700439e+06,   1.50199284e+01,   2.93666210e+01,\
         1.09005939e+02,   6.35043912e+00,   4.04673143e+01,\
         1.18546830e+03,   1.12011498e+02])
    sigmasST = np.array([ 0.06 ,  0.07 ,  0.05 ,  0.06 ,  0.06 ,  0.055,  0.01 ,  0.07 ,\
        0.06 ,  0.065,  0.01 ,  0.01 ,  0.01 ,  0.085])
    ##mcONLYcall
    #outFname = Sim.simulatorFunc('MC')
    ## SA call
    #(paramsSA,sigmasSA) = Sim.simulatorFunc('SA')
    ## ST call
    #print '-'*50+'\n\n\n'
    #(paramsST,sigmasST) = Sim.simulatorFunc('ST',paramsSA,sigmasSA)
    ##MCMC call
    #print '-'*50+'\n\n\n'
    outFname = Sim.simulatorFunc('MCMC',paramsST,sigmasST)
    print 'FINAL OUTFILE :\n'+outFname    
    
    
    log.info("End of SMODT2.0 main")
    ## move local log file to output dir!!!!

if __name__ == '__main__':
    smodt()