#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
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
    log = tools.getLogger('main',dir=settingsDict['finalFolder'],lvl=settingsDict['logLevel'])
    log.debug("Prepend string passed in was '"+settingsDict['prepend']+"'")
    Sim = simulator.Simulator(settingsDict)
    
#     e = 0.4
#     Sys_Dist_PC = 5.0
#     Mass1 = 1.0
#     Mass2 = 0.2
#     Omega = 60.0
#     omega = 110.0
#     T = 2457000.0
#     T_center = 2457000.0
#     P = 15.0
#     inc =  30.0
#     offset = 0.0
    
#     params = np.array([Mass1,Mass2,Sys_Dist_PC,Omega,e,T,T_center,P,inc,omega,0,0,0,offset])
#     paramsST = np.array( [  9.37145408e-01,   1.98238112e-01,   5.02014226e+00,\
#          6.15361129e+01,   4.15495712e-01,   2.45700439e+06,\
#          2.45700439e+06,   1.50199284e+01,   2.93666210e+01,\
#          1.09005939e+02,   6.35043912e+00,   4.04673143e+01,\
#          1.18546830e+03,   1.12011498e+02])
#     sigmasST = np.array([ 0.06 ,  0.07 ,  0.05 ,  0.06 ,  0.06 ,  0.055,  0.01 ,  0.07 ,\
#         0.06 ,  0.065,  0.01 ,  0.01 ,  0.01 ,  0.085])
    if True:
        ##mcONLYcall
        #outFname = Sim.simulatorFunc('MC')
        ## SA call
        #print '-'*50+'\n\n\n'
        (paramsSA,sigmasSA) = Sim.simulatorFunc('SA')
        ## ST call
        #print '-'*50+'\n\n\n'
        (paramsST,sigmasST) = Sim.simulatorFunc('ST',paramsSA,sigmasSA)
        ##MCMC call
        #print '-'*50+'\n\n\n'
        outFname = Sim.simulatorFunc('MCMC',paramsST,sigmasST)
        print 'FINAL OUTFILE :\n'+outFname    
    
    ## Post-processing goes here!!
    if False:
        ## Test to kill soon after done $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        omega = 20.0
        e = 0.5
        P = 5.0
        T = 2457000.0
        Tc=0
        omegas = [0,45,90,135,180,225,270,315]
        sjCenter=True
        print 'P = '+str(P)+', e = '+str(e)+", sjCenter = "+repr(sjCenter)
        for omega in omegas:
            Tc=0
            if sjCenter:
                print '\nBefore: T = '+str(T)+", Tc = "+str(Tc)+", omega = "+str(omega)+" -> "+str(omega/(360*0.01))+"%, 90-omega -> "+str((90-omega)/(360*0.01))+"%"
                ta = np.pi/2.0 - omega*(np.pi/180.0)
            else:
                print '\nBefore: T = '+str(T)+", Tc = "+str(Tc)+", omega = "+str(omega)+" -> "+str(omega/(360*0.01))+"%, 270-omega -> "+str((270-omega)/(360*0.01))+"%"
                ta = 1.5*np.pi - omega*(np.pi/180.0)
            halfE = np.arctan2(np.sqrt(1.0-e)*np.sin(ta/2.0),np.sqrt(1.0+e)*np.cos(ta/2.0))
            mTTc = 2.0*halfE-e*np.sin(2.0*halfE)
            deltaT = (mTTc*P*365.2442)/(2.0*np.pi)
            if Tc==0:
                Tc = T+deltaT
            else:
                T = Tc-deltaT
            print 'After: T = '+str(T)+", Tc = "+str(Tc)+", deltaT = "+str(deltaT)+" -> %orbit = "+str(deltaT/(P*365.2442*0.01))
            print '-'*100
    
    
    
        
    log.info("End of SMODT2.0 main")

if __name__ == '__main__':
    smodt()