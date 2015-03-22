#@Author: Kyle Mede, kylemede@gmail.com, written for the University of Tokyo during UTRIP 2011.
""" This is a Markov Chain Monte Carlo simulator using a Metropolis-Hastings algorithm
    to find the orbital parameters for either a planetary or binary star system 
    with provided separation angle, system distance and position angle.
    This version uses UNIFORM random variables as inputs."""

import math
import numpy as np
import random as rand
from orbitToolbox import *
from basicOrbSimulator import basicOrbSim

#*********** In this brick is all you need to play with *************    
numSamples = int(1e5)
verbose = False
numSamplePrints = 25.0
# Ranges for acceptable random number inputs
longAN_degMIN = 0.001 # [deg]
longAN_degMAX = 179.999 # [deg]
eMIN = 0.001
eMAX = 0.999 
a1MIN = 0.1  # for binary systems. a1=0.0 always in planet systems. [AU]
a1MAX = 20.0 # for binary systems. a1=0.0 always in planet systems. [AU]
a1Multiplier = 10.0 # for binaries, a2MAX = a1*a1Multiplier
a2MIN = 1.0  # for planets. [AU]
a2MAX = 50.0 # for planets. [AU]
periodMIN = 1.0 # [yrs]
periodMAX = 100.0 # [yrs]
inclination_degMIN = 0.001 # [deg]
inclination_degMAX = 179.999 # [deg]
argPeri_degMIN = 0.001 # [deg]
argPeri_degMAX = 89.999 # [deg]
#### MEASURED VALUES ########################
Sep_Angle_arcsec_measured_REALs = [0.210,0.326] # used to calc kai^2 
SA_mean_errors = [0.027,0.013]
PA_deg_measured_REALs = [211.49,210.64] # used to calc kai^2
PA_mean_errors = [1.9, 1.2]
epochs = [2454781.50,2455194.50]
chiSquaredMax = 80
Sys_Dist_PC = 19.3
binary = False # Calculate for a binary star system (True) or a planet (False)

plotFileTitle = 'mcmcUniform-numSamples_'+str(numSamples)+'_chiSquaredMax_'+str(chiSquaredMax)
showPlots = True
#*********************************************************************

print '\n*************  STARTING THE MARKOV CHAIN  *****************'
print 'Number of sample orbits being created = ',numSamples

# init input lists
longAN_degs = []
es = []
Ts = []
a1s = []
a2s = []
periods = []
inclination_degs = []
argPeri_degs = []

# init output lists
ns2 = []
Ms2 = []
Es2 = []
thetas2 = []
Sep_Dists2 = []
Sep_Angle_arcsec_measured_models2 = []
PA_deg_measured_models2 = []

inclination_degs_lastOnes = [] #$$$$$$$$$$$$$$$$$$$$$$$$$

## Use the basicOrbSim to get initial guess values for the model inputs for first epoch
(longAN_deg_initial, e_initial, T_initial, a1_initial, a2_initial, period_initial, \
 inclination_deg_initial, argPeri_deg_initial, ns_initial, Ms_initial, Es_initial, thetas_initial, Sep_Dists_initial,\
 Sep_Angle_arcsec_measured_models_initial, PA_deg_measured_models_initial) =  \
 basicOrbSim(10000, Sep_Angle_arcsec_measured_REALs[0], SA_mean_errors[0], PA_deg_measured_REALs[0],\
  PA_mean_errors[0], chiSquaredMax, epochs[0], Sys_Dist_PC, binary=binary, verbose=verbose)

## Make an initial chi**2 value to start chain from
chi_squared_total_last = chiSquaredMax-0.1

# variables for the success rate print block in chain loop
printTime = numSamples/numSamplePrints
printCount = 0

## Start to create samples and check if they are are accepted
## ie. Start the chain
for i in range(1,numSamples+1):
    
    # block to control printing success rate to screen
    printCount = printCount + 1
    if printCount==printTime:
        printCount = 0
        print str(len(es))+' successful samples out of '+str(i)
        
    ## Input variables made with uniform random numbers
    #  Note: These will be stored if they satisfy the chi**2 max
    #        so only the good ones are stored to save processor time
    inclination_deg = rand.uniform(inclination_degMIN,inclination_degMAX)
    longAN_deg = rand.uniform(longAN_degMIN,longAN_degMAX)
    e = rand.uniform(eMIN,eMAX)
    period = rand.uniform(periodMIN,periodMAX) # [yrs]
    T = rand.uniform(epochs[0]-period*365.0,epochs[0]) # thus between a full period ago and now
    if binary:
        # binary star system with primary being the larger star so smaller semi-major
        a1 = rand.uniform(a1MIN,a1MAX)    
        a2 = rand.uniform(a1,a1Multiplier*a1)
    elif binary is False:
        a1 = 0.0 # always the case for a planet-star system as foci commonly inside star.
        a2 = rand.uniform(a2MIN,a2MAX) 
    argPeri_deg = rand.uniform(argPeri_degMIN,argPeri_degMAX)
    
    (chi_squared_total_cur, ns, Ms, Es, thetas, Sep_Dists, Sep_Angle_arcsec_measured_models, PA_deg_measured_models) =\
     multiEpochOrbCalc(Sep_Angle_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors, chiSquaredMax,\
                        epochs, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, a1, a2, period, argPeri_deg, verbose=verbose)
    
    # Does current orbit satisfy chi**2?
    if (chi_squared_total_cur != False) and (chi_squared_total_cur<=chiSquaredMax):
        ##### In my original model, I would calc the gauss dist in each direction, BUT
        ##### as the gauss dist is symmetric, the fraction of the two will ALWAYS = 1.
        ##### Thus I am only using a ratio of the probability densities (ie. a likelihood ratio).
        probRatio = (chiSquaredMax-chi_squared_total_cur)/(chiSquaredMax-chi_squared_total_last)
        alpha = rand.uniform(0.0,1.0) #### maybe set this to 1.0 !?!?!?!?!? $$$$$$$$$$$$$
        
        if alpha <= probRatio:
            # keep the NEW value, ie. sample accepted!
            chi_squared_total_last = chi_squared_total_cur
        
            # store input values to the model
            longAN_degs.append(longAN_deg)
            es.append(e)
            Ts.append(T)
            a1s.append(a1)
            a2s.append(a2)
            periods.append(period)
            inclination_degs.append(inclination_deg)
            argPeri_degs.append(argPeri_deg)
            # store output orbital elements 
            ns2.append(ns)
            Ms2.append(Ms)
            Es2.append(Es)
            thetas2.append(thetas)
            Sep_Dists2.append(Sep_Dists)
            Sep_Angle_arcsec_measured_models2.append(Sep_Angle_arcsec_measured_models)
            PA_deg_measured_models2.append(PA_deg_measured_models)
            
            if i>(numSamples/4): #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                inclination_degs_lastOnes.append(inclination_deg)

print '*************** Markov Chain complete *******************'

#TEMP print mean orbital element values for kicks
print '\n*************** RESULTS SUMMARY  *********************'
print 'Measured Separation Angle of stars REAL ["] = ',Sep_Angle_arcsec_measured_REALs[0]
print 'Measured Position Angle in image REAL [deg] = ',PA_deg_measured_REALs[0]
print 'Number of orbits created ',numSamples
print 'chi**2 max = ',chiSquaredMax
print 'number of orbits in the chain = ',len(es)

print '\nKnown/Measured Input Parameters:'
print 'Epoch of observation/image [yrs] = ',epochs[0]
print 'Measured system distance from Earth [PC] = ',Sys_Dist_PC

print '\nRandom Input Parameters means:'
print 'median Longitude of Ascending Node [deg]= ', np.median(longAN_degs)
print 'median e = ', np.median(es)
print 'median T [julian date] = ', np.median(Ts)
print 'median a1 [AU] = ',np.median(a1s)
print 'median a2 [AU] = ',np.median(a2s)
print 'median period [yrs] = ',np.median(periods)
print 'median inclination [deg] = ',np.median(inclination_degs)
print 'median Argument of Pariapsis [deg] = ',np.median(argPeri_degs)
    
print '\nOutput Parameters means:'
print 'median Sep_Dist [AU] = ',np.median(Sep_Dists2)
print 'median True Anomaly [deg] = ', np.median(thetas2)
print 'median Eccentric Anomaly [deg] = ',np.median(Es2)
print 'median Mean Anomaly [deg] = ', np.median(Ms2)
print 'median Mean Motion [rad/yr] = ',np.median(ns2)
print 'median Measured Separation Angle of stars MODEL ["] = ',np.median(Sep_Angle_arcsec_measured_models2)
print 'median Measured Position Angle in image MODEL [deg] = ',np.median(PA_deg_measured_models2)
# plot results
orbElementPlotter(plotFileTitle, longAN_degs, es, Ts, a1s, a2s, periods, inclination_degs,argPeri_degs,\
                  Sep_Dists2, thetas2, Es2, Ms2, ns2, showPlots=showPlots)
# write resultant data to files
dataWriter(plotFileTitle, longAN_degs, es, Ts, a1s, a2s, periods, inclination_degs, argPeri_degs,\
                  Sep_Dists2, thetas2, Es2, Ms2, ns2)
print '\n******************* FINISHED CREATING ORBITS   *************************\n'







