import threading
from mcmcOrbSimulatorUniform4 import mcmcUniformOrbSim

numSamples=1 
silent=True

class MCMCThreadManager(threading.Thread):

	def run(self):
		
		## importing global variables
		global numSamples
		global silent
		
		(longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds) = mcmcUniformOrbSim(numSamples, silent=True)
         
        return (longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds)
         
	def setVars(self, numSamples1, silent1=True):
		
		global numSamples
		global silent
		self.numSamples = numSamples1
		self.silent = silent1 
		print '\n input variables set'
		 
		
		
		
		
		
		
	
                
		
		
 