#@Author: Kyle Mede, kylemede@gmail.com, written for the University of Tokyo during UTRIP 2011.
""" This is a Markov Chain Monte Carlo simulator using a Metropolis-Hastings algorithm
    to find the orbital parameters for either a planetary or binary star system 
    with provided separation angle, system distance and position angle.
    This version uses GAUSSIAN distributed random variables as inputs."""

import math
import numpy as np
import random as rand
from orbitToolbox import *
from basicOrbSimulator import basicOrbSim

#*********** In this brick is all you need to play with *************    
numSamples = int(1e5)
verbose = False
numSamplePrints = 50.0
tunningAdjustments = 25 # how many times to adjust a single sigma before moving to the next
tunningSamples = 100 # how many samples to use for calculating acceptance rate
tunningRepeats = 5 # how many times to repeat the tunning sequence
sigmaMAX = 5
triesMAX = 100000
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
PA_deg_measured_REALs = [211.49,210.64] # used to calc kai^2
epochs = [2454781.50,2455194.50]
chiSquaredMax = 50
Sys_Dist_PC = 19.3
binary = False # Calculate for a binary star system (True) or a planet (False)

title = 'mcmcGauss-numSamples_'+str(numSamples)+'_chiSquaredMax_'+str(chiSquaredMax)
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

## Use the basicOrbSim to get initial guess values for the model inputs for first epoch
(longAN_deg_initial, e_initial, T_initial, a1_initial, a2_initial, period_initial, \
 inclination_deg_initial, argPeri_deg_initial, ns_initial, Ms_initial, Es_initial, thetas_initial, Sep_Dists_initial,\
 Sep_Angle_arcsec_measured_models_initial, PA_deg_measured_models_initial) =  \
 basicOrbSim(10000, Sep_Angle_arcsec_measured_REALs[0],\
  PA_deg_measured_REALs[0], chiSquaredMax, epochs[0], Sys_Dist_PC, binary=binary, verbose=verbose)

## Record the initial sigma values to start chain from;
#  these will be tuned throughout the chain.
sigmas = []
sigmas.append([1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0])
sigmaBeingTunned = 0
timesTunned = 0
acceptedCounter = 0
totalCounter = 0
roundsTunned = 0

# store input values to the model
longAN_degs.append(longAN_deg_initial)
es.append(e_initial)
Ts.append(T_initial)
a1s.append(a1_initial)
a2s.append(a2_initial)
periods.append(period_initial)
inclination_degs.append(inclination_deg_initial)
argPeri_degs.append(argPeri_deg_initial)

# Create an initial latest chi**2
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
            
    ## Input variables made with random numbers from a gaussian distribution
    #  Note: These will be stored if they satisfy the chi**2 max
    #        so only the good ones are stored to save processor time
    inRange = False
    tries = 0
    lastFew = []
    while inRange is False:
        tries+=1
        inclination_deg = rand.gauss(inclination_degs[len(es)-1], sigmas[i-1][0])
        if (inclination_deg>=inclination_degMIN)and(inclination_deg<=inclination_degMAX):
            inRange = True
        if tries>triesMAX:
            lastFew.append(inclination_deg)
        if tries>(triesMAX+4):
            print 'inclination_deg while loop hung, lastFew ',repr(lastFew)
            inclination_deg = inclination_degs[len(es)-1]
            break
            
    inRange = False
    tries = 0
    lastFew = []
    while inRange is False:
        tries+=1
        longAN_deg = rand.gauss(longAN_degs[len(es)-1],sigmas[i-1][1])
        if (longAN_deg>=longAN_degMIN)and(longAN_deg<=longAN_degMAX):
            inRange = True
        if tries>triesMAX:
            lastFew.append(longAN_deg)
        if tries>(triesMAX+4):
            print 'longAN_deg while loop hung, lastFew ',repr(lastFew)
            longAN_deg = longAN_degs[len(es)-1]
            break
            
    inRange = False
    tries = 0
    lastFew = []
    while inRange is False:
        tries+=1
        e = rand.gauss(es[len(es)-1],sigmas[i-1][2])
        if (e>=eMIN)and(e<=eMAX):
            inRange = True
        if tries>triesMAX:
            lastFew.append(e)
        if tries>(triesMAX+4):
            print 'e while loop hung, lastFew ',repr(lastFew)
            e = es[len(es)-1]
            break
            
    inRange = False
    tries = 0
    lastFew = []
    while inRange is False:
        tries+=1
        period = rand.gauss(periods[len(es)-1],sigmas[i-1][3])
        if (period>=periodMIN)and(period<=periodMAX):
            inRange = True
        if tries>triesMAX:
            lastFew.append(period)
        if tries>(triesMAX+4):
            print 'period while loop hung, lastFew ',repr(lastFew)
            period = periods[len(es)-1]
            break
            
    inRange = False
    tries = 0
    lastFew = []
    while inRange is False:
        tries+=1
        T = rand.gauss(Ts[len(es)-1],sigmas[i-1][4])
        if (T>=(epochs[0]-period*365.0))and(T<=epochs[0]):
            inRange = True
        if tries>triesMAX:
            lastFew.append(T)
        if tries>(triesMAX+4):
            print 'T while loop hung, lastFew ',repr(lastFew)
            T = Ts[len(es)-1]
            break
            
    if binary:
        # binary star system with primary being the larger star so smaller semi-major
        inRange = False
        tries = 0
        lastFew = []
        while inRange is False:
            tries+=1
            a1 = rand.gauss(a1s[len(es)-1],sigmas[i-1][5])
            if (a1>=a1MIN)and(a1<=a1MAX):
                inRange = True
            if tries>triesMAX:
                lastFew.append(a1)
            if tries>(triesMAX+4):
                print 'a1 while loop hung, lastFew ',repr(lastFew)
                a1 = a1s[len(es)-1]
                break
            
        inRange = False
        tries = 0
        lastFew = []
        while inRange is False:
            a2 = rand.gauss(a2s[len(es)-1],sigmas[i-1][6])
            if (a2>=a1)and(a2<=(a1*a1Multiplier)):
                inRange = True
            if tries>triesMAX:
                lastFew.append(a2)
            if tries>(triesMAX+4):
                print 'a2 while loop hung, lastFew ',repr(lastFew)
                a2 = a2s[len(es)-1]
                break
                
    elif binary is False:
        a1 = 0.0 # always the case for a planet-star system as foci commonly inside star.
        inRange = False
        tries = 0
        lastFew = []
        while inRange is False:
            tries+=1
            a2 = rand.gauss(a2s[len(es)-1],sigmas[i-1][6])
            if (a2>=a2MIN)and(a2<=a2MAX):
                inRange = True
            if tries>triesMAX:
                lastFew.append(a2)
            if tries>(triesMAX+4):
                    print 'a2 while loop hung, lastFew ',repr(lastFew)
                    a2 = a2s[len(es)-1]
                    break
                
    inRange = False
    tries = 0
    lastFew = []
    while inRange is False:
        tries+=1
        argPeri_deg = rand.gauss(argPeri_degs[len(es)-1],sigmas[i-1][7])
        if (argPeri_deg>=argPeri_degMIN)and(argPeri_deg<=argPeri_degMAX):
            inRange = True    
        if tries>triesMAX:
            lastFew.append(argPeri_deg)
        if tries>(triesMAX+4):
            print 'argPeri_deg while loop hung, lastFew ',repr(lastFew)
            argPeri_deg = argPeri_degs[len(es)-1]
            break
    
    (chi_squared_total_cur, ns, Ms, Es, thetas, Sep_Dists, Sep_Angle_arcsec_measured_models, PA_deg_measured_models) =\
     multiEpochOrbCalc(Sep_Angle_arcsec_measured_REALs, PA_deg_measured_REALs, chiSquaredMax, epochs, Sys_Dist_PC,\
                       inclination_deg, longAN_deg, e, T, a1, a2, period, argPeri_deg, verbose=verbose)
    
    # Does current orbit satisfy chi**2?
    if (chi_squared_total_cur != False) and (chi_squared_total_cur<=chiSquaredMax):
        ##### In my original model, I would calc the gauss dist in each direction, BUT
        ##### as the gauss dist is symmetric, the fraction of the two will ALWAYS = 1.
        ##### Thus I am only using a ratio of the probability densities (ie. a likelihood ratio).
        probRatio = (chiSquaredMax-chi_squared_total_cur)/(chiSquaredMax-chi_squared_total_last)
        alpha = 0.9 # rand.uniform(0.0,1.0) #### maybe set this to 1.0 !?!?!?!?!? $$$$$$$$$$$$$
        
        if alpha <= probRatio:
            # keep the NEW value, ie. sample accepted!
            chi_squared_total_last = chi_squared_total_cur
            acceptedCounter = acceptedCounter + 1
                
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

    ##### Sigma tunning block #######
    if (roundsTunned-1)<tunningRepeats:
        # Tunned less times than number of repeats requested, so continue
        
        # copy previous batch of sigmas over to current
        # only change current sigma being adjusted, rest stay same. 
        sigmas.append(sigmas[i-1][:])
        totalCounter = totalCounter + 1
        
        if totalCounter==tunningSamples:
            # it is time to recalculate the acceptance rate and adjust current sigma accordingly 
            acceptRate = float(acceptedCounter)/float(tunningSamples)
            acceptedCounter = 0
            totalCounter = 0
            timesTunned = timesTunned + 1
            
            sigma = sigmas[i-1][sigmaBeingTunned]
            if sigma<=sigmaMAX:
                if acceptRate>0.55:
                    # ie. acceptance rate is too high, so raise current sigma being adjusted
                    sigmas[i][sigmaBeingTunned] = sigma + 0.2
                elif acceptRate<0.45:
                    # ie. acceptance rate is too low, so lower current sigma being adjusted
                    sigmas[i][sigmaBeingTunned] = sigma - 0.2
            else:
                sigmas[i][sigmaBeingTunned] = sigma
        
            if timesTunned>=tunningAdjustments:
                # done tunning this sigma so move to next
                timesTunned = 0
                print 'Last acceptRate, ',acceptRate#####$$$$$$$$$$$$$$$$$$
                print 'Last sigma, ', sigma
                if sigmaBeingTunned==(len(sigmas[:][0])-1):
                    # reached last sigma to tune, so go back to first one
                    sigmaBeingTunned = 0
                    roundsTunned = roundsTunned+1
                    print 'Finished a round of sigma tunning' #$$$$$$$$
                else:
                    sigmaBeingTunned = sigmaBeingTunned + 1
            else:    
                timesTunned = timesTunned + 1
    else:
        sigmas.append(sigmas[i-1][:])
    ##### Sigman tunning block done ##########
        
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
orbElementPlotter(longAN_degs, es, Ts, a1s, a2s, periods, inclination_degs,argPeri_degs,\
                  Sep_Dists2, thetas2, Es2, Ms2, ns2)
print '\n******************* FINISHED CREATING ORBITS   *************************\n'







