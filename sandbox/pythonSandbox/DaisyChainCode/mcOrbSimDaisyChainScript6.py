from mcProcessManager3 import mcmcSimStarter
 
numSamples = int(1.25e8)
numProcesses = 5
numSamplesString = str(numSamples/int(1e6))+'-million'
roundString = 'Round-6'
title = 'Tau-Boo-NEW-Uniform-mcONLY-trial11-'+numSamplesString+'-'+roundString
print '\n** Starting Round '+roundString+' of '+str(numProcesses)+' processes of '+numSamplesString+' samples each ** \n'
mcmcSimStarter(title, numSamples, numProcesses, silent=True)
