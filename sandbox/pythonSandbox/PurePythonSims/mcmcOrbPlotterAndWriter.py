
from orbitToolbox2 import *

def plotAndWrite(plotFileTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
        ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds, showPlots=False, summaryOnly=False, verbose=False):
    """
    This is a function to accept the input and output parameters from the mcmcOrbSimulatorUniform4 and then 
    plot them and write the lists to disk.
    """
    
    ## write resultant data to files
    dataWriter(plotFileTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds)

    ## plot results
    orbElementPlotter(plotFileTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
        ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds, showPlots=showPlots, summaryOnly=False, verbose=verbose)