
from mcmcThreadManager import MCMCThreadManager
from mcmcOrbPlotterAndWriter import plotAndWriter

def mcmcSimulationStarter(numThreads, numSamples1, titleMod1, silent1=True):
    
    # init MASTER input lists
    longAN_degs_MASTER = []
    es_MASTER = []
    Ts_MASTER = []
    periods_MASTER = []
    inclination_degs_MASTER = []
    argPeri_degs_MASTER = []
    
    # init MASTER output lists
    a1s2_MASTER = []
    a2s2_MASTER = []
    chiSquareds_MASTER = []
    Sep_Dists2_MASTER = []
    ns2_MASTER = []
    Ms2_MASTER = []
    Es2_MASTER = []
    thetas2_MASTER = []
    PA_deg_measured_models2_MASTER = []
    
    # create threat manager object and start thread
    threadManager = MCMCThreadManager()
    threadManager.setVars(self.numSamples1, silent1)
    
    for threadNum in range(1,numThreads+1):
        print '\n Starting thread '+str(threadNum)
        
        (longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds) = threadManager.start()
                
        # append to MASTER input lists
        longAN_degs_MASTER.append(longAN_degs)
        es_MASTER.append(es) 
        Ts_MASTER.append(Ts) 
        periods_MASTER.append(periods)
        inclination_degs_MASTER.append(inclination_degs)
        argPeri_degs_MASTER.append(argPeri_degs)
        
        # append to MASTER output lists
        a1s2_MASTER.append(a1s2)
        a2s2_MASTER.append(a2s2)
        chiSquareds_MASTER.append(chiSquareds)
        Sep_Dists2_MASTER.append(Sep_Dists2)
        ns2_MASTER.append(ns2)
        Ms2_MASTER.append(Ms2)
        Es2_MASTER.append(Es2)
        thetas2_MASTER.append(thetas2)
    
        print '\n Finished thread '+str(threadNum)
    
    # All threads done, so plot and write the data to disk
    plotAndWriter(plotFileTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
        ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds, showPlots=False, summaryOnly=False, verbose=False)

