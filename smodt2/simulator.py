import numpy as np
import tools


class Simulator(object):
    """
    This is the Simulator parent class.  
    It contains the functions to perform basic 'shotgun' Monte Carlo, 
    Simulated Annealing, Sigma Tunning, and pure MCMC simulations.
    """
    def __init__(self,settingsDict):
        self.settingsDict = settingsDict
        self.log = tools.getLogger('main.simulator',lvl=10,addFH=False)
        tools.logSystemInfo(self.log)
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
        tools.test()
        
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