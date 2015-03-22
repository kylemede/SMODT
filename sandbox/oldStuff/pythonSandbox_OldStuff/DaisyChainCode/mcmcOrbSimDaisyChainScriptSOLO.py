#from mcmcProcessManager3 import mcmcSimStarter
# 
#numSamples = int(1e8)
#numProcesses = 1
#numSamplesString = str(numSamples/int(1e6))+'-million'
#roundString = 'SOLO'
#title = '2MASS-J1534-Triangular-'+numSamplesString+'-'+roundString
#print '\n** Starting Round '+roundString+' of '+str(numProcesses)+' processes of '+numSamplesString+' samples each ** \n'
#mcmcSimStarter(title, numSamples, numProcesses, silent=True)

from mcmcOrbSimulatorUniform5 import mcmcUniformOrbSim
 
numSamples = int(1e6)
numSamplesString = str(numSamples/int(1e6))+'-million'
roundString = 'SOLO'
title = '2MASS-J1534-Triangular-'+numSamplesString+'-'+roundString
print '\n** Starting Round '+roundString+' of '+str(0)+' processes of '+numSamplesString+' samples each ** \n'
mcmcUniformOrbSim(numSamples, silent=False)