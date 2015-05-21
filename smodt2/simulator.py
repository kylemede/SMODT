import numpy as np
import os
from scipy.constants.codata import precision
#np.set_printoptions(precision=15)
import tools
import timit
from settings_and_inputData import constants


class Simulator(object):
    """
    This is the Simulator parent class.  
    It contains the functions to perform basic 'shotgun' Monte Carlo, 
    Simulated Annealing, Sigma Tunning, and pure MCMC simulations.
    """
    def __init__(self,settingsDict):
        self.settingsDict = settingsDict
        self.log = tools.getLogger('main.simulator',lvl=100,addFH=False)
        tools.logSystemInfo(self.log)
        self.Orbit = tools.cppTools.Orbit()
        self.realData = tools.loadRealData(os.path.join(settingsDict['settingsDir'],settingsDict['prepend']))
        self.Orbit.loadRealData(self.realData)
        self.Orbit.loadConstants(constants.Grav,constants.pi,constants.KGperMsun, constants.daysPerYear,constants.secPerYear,constants.MperAU)
        (self.rangeMaxs,self.rangeMins,self.starterSigmas) = self.starter()
        self.seed = timeit.default_timer()
        self.log.info("random number seed = "+str(self.seed))
        np.random.seed(self.seed)
        #Load real data here? doesn't change so might as well
        #could also make the empty model data array here too
        ##Examples
        #(memTrackProc,memLogFilename) = self.starter()
        #self._memLogFilename = memLogFilename
        #self._memTrackProc = memTrackProc
           
    def starter(self):
        """
        Things needed by all simulator modes that can be internal, but 
        need code to load them up.
        """
        self.log.info("In Simulator.starter")
        ##load up range min,max and sigma arrays
        #params order =[Mass1,Mass2,Sys_Dist_PC,Omega,e,T,T_center,P,inc,omega,aTot,chiSquared,K,offset]
        rangeMaxs = [self.dictVal('mass1MAX'),\
               self.dictVal('mass2MAX'),\
               self.dictVal('distMAX'),\
               self.dictVal('OmegaMAX'),\
               self.dictVal('eMAX'),\
               self.dictVal('TMAX'),\
               self.dictVal('TMAX'),\
               self.dictVal('PMAX'),\
               self.dictVal('incMAX'),\
               self.dictVal('omegaMAX'),\
               0,\
               0,\
               self.dictVal('KMAX'),\
               self.dictVal('vMAXs')]
        rangeMins = [self.dictVal('mass1MIN'),\
               self.dictVal('mass2MIN'),\
               self.dictVal('distMIN'),\
               self.dictVal('OmegaMIN'),\
               self.dictVal('eMIN'),\
               self.dictVal('TMIN'),\
               self.dictVal('TMIN'),\
               self.dictVal('PMIN'),\
               self.dictVal('incMIN'),\
               self.dictVal('omegaMIN'),\
               0,\
               0,\
               self.dictVal('KMIN'),\
               self.dictVal('vMINs')]
        sigmas = [0.1*(self.dictVal('mass1MAX')-self.dictVal('mass1MIN')),\
               0.1*(self.dictVal('mass2MAX')-self.dictVal('mass2MIN')),\
               0.1*(self.dictVal('distMAX')-self.dictVal('distMIN')),\
               0.1*(self.dictVal('OmegaMAX')-self.dictVal('OmegaMIN')),\
               0.1*(self.dictVal('eMAX')-self.dictVal('eMIN')),\
               0.1*(self.dictVal('TMAX')-self.dictVal('TMIN')),\
               0.1*(self.dictVal('TMAX')-self.dictVal('TMIN')),\
               0.1*(self.dictVal('PMAX')-self.dictVal('PMIN')),\
               0.1*(self.dictVal('incMAX')-self.dictVal('incMIN')),\
               0.1*(self.dictVal('omegaMAX')-self.dictVal('omegaMIN')),\
               0,\
               0,\
               0.1*(self.dictVal('KMAX')-self.dictVal('KMIN')),\
               []]
        for i in range(0,len(self.dictVal('vMINs'))):
            sigmas[-1].append(0.1*(self.dictVal('vMAXs')[i]-self.dictVal('vMINs')[i]))
        
        return (rangeMaxs,rangeMins,sigmas)
    
    def dictVal(self,key):
        """
        Get the value for a key in the settings dictionary.
        This will handle the values that are tuples and not
        returning the value.
        """
        if type(self.settingsDict[key])==tuple:
            return self.settingsDict[key][0]
        else:
            return self.settingsDict[key]
    
    def increment(self,params,sigmas=[],mcOnly=False):
        """
        Randomly increment one of the parameters
        """
        #params order =[Mass1,Mass2,Sys_Dist_PC,Omega,e,T,T_center,P,inc,omega,aTot,chiSquared,K,offset]
        #figure out which parameters are varying in this run
        paramInts = []
        for i in range(0,len(params)-1):
            if self.rangeMaxs[i]!=0:
                paramInts.append(i)

            
        if mcOnly:
            for i in range(0,len(params)-1):
                params[i]=np.random.uniform(self.rangeMins[i],self.rangeMaxs[i])
            for i in range(0,len(params[-1])):
                params[-1][i]=np.random.uniform(self.rangeMins[-1][i],self.rangeMaxs[-1][i])
        else:
            """
            """
        return params
    
    def monteCarlo(self):
        """
        Performs 'shotgun' Monte Carlo.
        Returns ?? not decided yet!!!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$ 
        """
        self.log.info("In Simulator.monteCarlo")
        self.log.info('starting c++ obj test')
        modelData = np.zeros((self.realData.shape[0],3))
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
        self.Orbit.calculate(modelData,params)
        params = self.increment(params,mcOnly=True)
        #atot = params[10]
        #K = params[12]
        print '\nmodelData = \n'+repr(modelData)
        baseFilename = 'testData.fits'
        tools.writeFits(baseFilename,modelData,self.settingsDict)
        #self.log.info('\nmodelData = \n'+repr(modelData))
        
        
    def simAnneal(self):
        """
        Performs Simulated Annealing.
        Returns ?? not decided yet!!!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$ 
        """
        self.log.info("In Simulator.simAnneal")
        
        
    def mcmc(self):
        """
        Performs pure Markov Chain Monte Carlo.
        Returns ?? not decided yet!!!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$ 
        """
        self.log.info("In Simulator.mcmc")
        