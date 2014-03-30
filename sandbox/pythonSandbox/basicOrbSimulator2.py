#@Author: Kyle Mede, kylemede@gmail.com, written for the University of Tokyo during UTRIP 2011.
import math
import numpy as np
import random as rand
import pylab 
from orbitToolbox2 import *


def basicOrbSim(numOrbs, Sep_Angle_arcsec_measured_REAL, SA_mean_error, PA_deg_measured_REAL, PA_mean_error, \
                                            chiSquaredMax, t, Sys_Dist_PC, Mass1=1, Mass2=1, verbose=False):
    """
    This basic orbit simulator will just create numOrbs worth of random orbits and return the median of the parameters
    of the orbits which passed the chiSquaredMax test.  This uses a UNIFORM random number generator, NOT Gaussian.
    
    Note:
    This was created as an early form of testing the orbCalculator
    to ensure it's equations were correct and is now used to determine a set of starting parameters for a 
    Markov Chain Monte Carlo simulation of binary orbits.
    """
    
    silent = True
    
    if not silent:
        print '\nbasicSim: ********  STARTING TO CREATE ORBITS  *********'
        print 'Number of orbits being created = ',numOrbs

    # init input lists
    longAN_degs = []
    es = []
    Ts = []
    periods = []
    inclination_degs = []
    argPeri_degs = []

    # init output lists
    ns = []
    Ms = []
    Es = []
    thetas = []
    Sep_Dists = []
    PA_deg_measured_models = []
    a1s = []
    a2s = []

    for i in range(1,numOrbs+1):
        ## Input variables made with uniform random numbers
        # Note: These will be stored if they satisfy the chi**2 max
        #       so only the good ones are stored to save processor time
        longAN_deg = rand.uniform(0.0,180.0)
        e = rand.uniform(0.001,0.999)
        period = rand.uniform(1.0,100.0) # [yrs]
        T = rand.uniform(t-period*365.0,t) # thus between (now-one period) and now
        inclination_deg = rand.uniform(0.0,180.0)
        argPeri_deg = rand.uniform(0.0,90.0)
        
        # call orbitCalculator to take random variables and calc orbital elements
        (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU, PA_deg_measured_model, a1, a2) = \
        orbitCalculator(t, Sys_Dist_PC, Sep_Angle_arcsec_measured_REAL, inclination_deg, longAN_deg, \
                                e, T, period, argPeri_deg, Mass1=Mass1, Mass2=Mass2, verbose=verbose)
                     
        ## calc both PA and SA kai's and check they are less than chiSquaredMax
        #SA_chi_squared = chiSquaredCalc(Sep_Angle_arcsec_measured_REAL, SA_mean_error, Sep_Angle_arcsec_measured_model)
        PA_chi_squared = chiSquaredCalc(PA_deg_measured_REAL, PA_mean_error, PA_deg_measured_model)
    
        if chiSquaredMax>=(PA_chi_squared):#+SA_chi_squared):  
            # store input values to the model
            longAN_degs.append(longAN_deg)
            es.append(e)
            Ts.append(T)            
            periods.append(period)
            inclination_degs.append(inclination_deg)
            argPeri_degs.append(argPeri_deg)
        
            # store output orbital elements in lists for later plotting
            ns.append(n)
            Ms.append(M_deg)
            Es.append(E_latest_deg)
            thetas.append(TA_deg)
            Sep_Dists.append(Sep_Dist_AU)
            PA_deg_measured_models.append(PA_deg_measured_model)
            a1s.append(a1)
            a2s.append(a2)
            
    ## for loop done, so calc mean values of model inputs that satisfied chiSquaredMax      
    longAN_degs_mean = np.mean(longAN_degs)
    es_mean = np.mean(es)
    Ts_mean = np.mean(Ts)    
    periods_mean = np.mean(periods)
    inclination_degs_mean = np.mean(inclination_degs)
    argPeri_degs_mean = np.mean(argPeri_degs)
    ## ditto for the outputs
    ns_mean = np.mean(ns)
    Ms_mean = np.mean(Ms)
    Es_mean = np.mean(Es)
    thetas_mean = np.mean(thetas)
    Sep_Dists_mean = np.mean(Sep_Dists)
    PA_deg_measured_models_mean = np.mean(PA_deg_measured_models)  
    a1s_mean = np.mean(a1s)
    a2s_mean = np.mean(a2s)  
          
    if verbose: 
        print 'basicSim: ********** RESULTS SUMMARY  ************'
        print 'Measured Position Angle in image REAL [deg] = ',PA_deg_measured_REAL
        print 'Number of orbits being created ',numOrbs
        print 'chi**2 max = ',chiSquaredMax
        print 'number of orbits with chi**2 under max = ',len(Sep_Dists)

        print '\nKnown/Measured Input Parameters:'
        print 'Epoch of observation/image [yrs] = ',t
        print 'Measured system distance from Earth [PC] = ',Sys_Dist_PC
        print 'Measured Separation Angle of stars REAL ["] = ',Sep_Angle_arcsec_measured_REAL

        print '\nRandom Input Parameters means:'
        print 'median Longitude of Ascending Node [deg]= ', longAN_degs_mean
        print 'median e = ', es_mean
        print 'median T [julian date] = ', Ts_mean
        print 'median period [yrs] = ',periods_mean
        print 'median inclination [deg] = ',inclination_degs_mean
        (n,bins,patches) = pylab.matplotlib.pyplot.hist(inclination_degs,bins=100)
        #pylab.plot(inclination_degs)
        pylab.show()
        print 'median Argument of Pariapsis [deg] = ',argPeri_degs_mean
    
        print '\nOutput Parameters means:'
        print 'median Sep_Dist [AU] = ',Sep_Dists_mean
        print 'median True Anomaly [deg] = ', thetas_mean
        print 'median Eccentric Anomaly [deg] = ',Es_mean
        print 'median Mean Anomaly [deg] = ', Ms_mean
        print 'median Mean Motion [rad/yr] = ',ns_mean
        print 'median Measured Position Angle in image MODEL [deg] = ',PA_deg_measured_models_mean
        print 'median a1 [AU] = ',a1s_mean
        print 'median a2 [AU] = ',a2s_mean
    if not silent:
        print 'basicSim: ********   FINISHED CREATING ORBITS   *********'

    return (longAN_degs_mean, es_mean, Ts_mean, periods_mean, inclination_degs_mean,\
             argPeri_degs_mean, ns_mean, Ms_mean, Es_mean, thetas_mean, Sep_Dists_mean,\
              PA_deg_measured_models_mean, a1s_mean, a2s_mean,)
       
    