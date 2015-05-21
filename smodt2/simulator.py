import numpy as np
import os
from scipy.constants.codata import precision
#np.set_printoptions(precision=15)
import tools
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
        #Load real data here? doesn't change so might as well
        #could also make the empty model data array here too
        ##Examples
        #(memTrackProc,memLogFilename) = self.starter()
        #self._memLogFilename = memLogFilename
        #self._memTrackProc = memTrackProc
           
    def starter(self):
        """
        Start things off??$$$$$$$$$$$$$$$$$$$$$$$$
        """
        self.log.info("In Simulator.starter")
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
        