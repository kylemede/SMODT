from mcProcessManager3 import mcmcSimStarter
 
numSamples = int(6.2e7)
numProcesses = 7
numSamplesString = str(numSamples/int(1e6))+'-million'
roundString = 'Round-2'
title = 'Tau-Boo-ALL-Uniform-mcONLY-trial13-'+numSamplesString+'-'+roundString
print '\n** Starting Round '+roundString+' of '+str(numProcesses)+' processes of '+numSamplesString+' samples each ** \n'
mcmcSimStarter(title, numSamples, numProcesses, silent=False)

