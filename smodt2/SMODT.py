#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
#import numpy as np
import tools
import simulator
import sys
import os
import numpy as np
from multiprocessing import Process

"""
    This is the 'main' of SMODT. 
    It will start things off, call the appropriate set of 
    simulation and post-processing steps.
"""


class singleProc(Process):
    """
    This is the Manager object that controls the a single processes for a 
    SMODT2.0 simulation run.  It is called by the multiProcessStarter once for 
    each chain/process requested by the user through the simulation settings 
    file.
    
    :param str settingsDict: settings Dictionary
    :param str fNameBase: File name, including the full path, for the output 
        data files.
    :param list stageList: List of stages to run ex.['MC','SA','ST','MCMC'] 
        lives.
    :param int chainNum: number of this chain
    """
    def __init__(self, settingsDict, stageList, chainNum=1):
        
        Process.__init__(self)
        self.chainNum = chainNum
        self.log = tools.getLogger('main.singleProcess',lvl=100,addFH=False)
        self.settingsDict = settingsDict 
        self.stageList = stageList
        self.Sim = simulator.Simulator(settingsDict)
        
    def run(self):
        self.settingsDict['']
        self.log.info('Starting to run process for file: '+fNameBase)
        if 'MC' in self.stageList:
            outMCFname = self.Sim.simulatorFunc('MC')
        if 'SA' in self.stageList:
            (paramsSA,sigmasSA,bestRedChiSqr) = self.Sim.simulatorFunc('SA')
        if bestRedChiSqr<self.settingsDict['chiMAX'][0]:
            if 'ST' in self.stageList:
                (paramsST,sigmasST) = self.Sim.simulatorFunc('ST',paramsSA,sigmasSA)
            if 'MCMC' in self.stageList:
                outMCMCFname = self.Sim.simulatorFunc('MCMC',paramsST,sigmasST)
                self.log.info('FINAL MCMC OUTFILE :\n'+outMCMCFname)
        else:
            self.log.critical()  
               
def multiProcStarter(settingsDict):
    """
    This will call singleProcessStarter to handle each individual chain.   
    The output dataFile for each chain will be returned as a list of filenames.
        
    :param dict settingsDict: A standard SMODT2.0 settings dictionary.
    """
    master = []
    for procNumber in range(settingsDict['numChains'][0]):
        master.append(singleProc(settingsDict,stageList,procNumber))
        master[procNumber].start()
    for procNumber in range(settingsDict['numChains'][0]):
        master[procNumber].join()    
    
    dataFiles = []
    for processNumber in range(settingsDict['numChains'][0]):
        dataFiles.append(master[processNumber].filename)
   
    return dataFiles

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
    paramsST = np.array([  9.50122214e-01,   2.21670225e-01,   5.05405775e+00,
         6.07928962e+01,   4.08889781e-01,   2.45700081e+06,
         2.45700081e+06,   1.50094339e+01,   2.71264582e+01,
         1.09266595e+02,   6.41461741e+00,   3.40786762e+01,
         1.20316127e+03,   1.83558842e+00])
    sigmasST = np.array([ 0.03,  0.07,  0.13,  0.13,  0.09,  0.17,  0.01,  0.09,  0.09,
        0.21,  0.01,  0.01,  0.01,  0.05])
    bestRedChiSqr=1.0
    outMCFname=''
    outMCMCFname=''
    if True:
        ##mcONLYcall
        #outMCFname = Sim.simulatorFunc('MC')
        ## SA call
        #(paramsSA,sigmasSA,bestRedChiSqr) = Sim.simulatorFunc('SA')
        if bestRedChiSqr<settingsDict['chiMAX'][0]:
            ## ST call
            #(paramsST,sigmasST) = Sim.simulatorFunc('ST',paramsSA,sigmasSA)
            ##MCMC call
            outMCMCFname = Sim.simulatorFunc('MCMC',paramsST,sigmasST)
            print 'FINAL MCMC OUTFILE :\n'+outMCMCFname    
        
            if True:
                ## Post-processing goes here!!
                if os.path.exists(outMCFname):
                    plotFilename = os.path.join(os.path.dirname(outMCFname),'SummaryPlotMC')
                    tools.summaryPlotter(outMCFname, plotFilename,stage="MC", shadeConfLevels=True)
                if os.path.exists(outMCMCFname):
                    plotFilename = os.path.join(os.path.dirname(outMCMCFname),'SummaryPlotMCMC')
                    tools.summaryPlotter(outMCMCFname, plotFilename,stage="MCMC", shadeConfLevels=True)
        else:
            log.critical("NO ORBIT WITH REDUCED CHISQUARED BELOW "+str(settingsDict['chiMAX'][0])+" WAS FOUND!!!")
    
    
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