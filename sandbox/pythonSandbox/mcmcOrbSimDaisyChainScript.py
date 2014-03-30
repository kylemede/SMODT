from mcmcOrbSimulatorUniform3 import mcmcUniformOrbSim

import gc
 
print '\nStarting Round 1\n'
mcmcUniformOrbSim(int(1e5), 'Round-1-', silent=True)

gc.collect()

print 'Starting Round 2\n'
mcmcUniformOrbSim(int(1e5), 'Round-2-', silent=True)

gc.collect()

print 'Starting Round 2\n'
mcmcUniformOrbSim(int(1e5), 'Round-3-', silent=True)

gc.collect()

print 'Starting Round 2\n'
mcmcUniformOrbSim(int(1e5), 'Round-4-', silent=True)