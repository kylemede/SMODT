import numpy as np
import tools
from settings_and_data import constants


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
    
    def monteCarlo(self):
        """
        Performs 'shotgun' Monte Carlo.
        Returns ?? not decided yet!!!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$ 
        """
        self.log.info("In Simulator.monteCarlo")
        #tools.test()
        self.log.info('starting c++ obj test')
        #self.Orbit.testDouble = 20.0
        a = np.ones((10,7))
        #a = np.array([[1.0,2.0],[3.0,4.0]])
        self.Orbit.loadRealData(a)
        self.Orbit.loadConstants(constants.Grav,constants.pi,constants.KGperMsun, constants.daysPerYear,constants.secPerYear,constants.MperAU)
        b = np.zeros((10,3))
        params = np.ones((14))
        self.Orbit.calculate(b,params)
        self.log.info('\nb = \n'+repr(b))
        
        
    def simAnneal(self):
        """
        Performs Simulated Annealing.
        Returns ?? not decided yet!!!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$ 
        """
        self.log.info("In Simulator.simAnneal")
        tools.test()
        
    def mcmc(self):
        """
        Performs pure Markov Chain Monte Carlo.
        Returns ?? not decided yet!!!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$ 
        """
        self.log.info("In Simulator.mcmc")
        tools.test()