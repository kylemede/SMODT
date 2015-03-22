#@Author: Kyle Mede, kylemede@gmail.com, written for the University of Tokyo during UTRIP 2011.
""" This is a Markov Chain Monte Carlo simulator using a Metropolis-Hastings algorithm
    to find the orbital parameters for either a planetary or binary star system 
    with provided separation angle, system distance and position angle.
    This version uses Triangular random variables as inputs.
    
    For a complete description of how this simulator works, please see mcmcOrbSimulatorTriangular2.pdf document.\
    *********** STILL NOT WRITTEN **** SAME AS UNIFORM BUT WITH Triangular RANDOM NUMBERS !!!!"""

import math
import numpy as np
import random as rand
import time
from orbitToolbox2 import *
from basicOrbSimulator2 import basicOrbSim
import gc

def mcmcTriangularOrbSim(numSamples, titleMod, silent=True):

    #*********** In this brick is all you need to play with *************    
    #********************************************************************
    chiSquaredMax = 1000.0 
    verbose = False                                                               
    numSamplePrints = 10.0        
    plotFileTitle = 'TauBoo_'+'mcmcTriangular-'+str(numSamples)+'Samples_'+titleMod
    showPlots = False                                          
    #### Ranges for acceptable random number inputs ######                                        
    longAN_degMIN = 0.001 # [deg]                                                        
    longAN_degMAX = 179.999 # [deg]                                                        
    eMIN = 0.001                                                                        
    eMAX = 0.999                                                                        
    periodMIN = 1.0 # [yrs]                                                                
    periodMAX = 2500.0 # [yrs]                                                            
    inclination_degMIN = 0.001 # [deg]                                                    
    inclination_degMAX = 179.999 # [deg]                                                
    argPeri_degMIN = 0.001 # [deg]                                                        
    argPeri_degMAX = 179.999#89.999 # [deg]                                                    
    #### MEASURED VALUES ########################                                        
    #Sep_Angle_arcsec_measured_REALs = [0.411,0.210,0.326,0.345] # used to calc kai^2 (ACTUALLY NOT IN THIS VERSION, ONLY VERISON 1)(CURRIE)
    #SA_mean_errors = [0.008,0.027,0.013,0.012] (CURRIE)
    #PA_deg_measured_REALs = [31.7,211.49,210.64,209.8] # used to calc kai^2 (CURRIE)
    #PA_mean_errors = [1.3,1.9, 1.2,0.8] (CURRIE)
    #epochs = [2452953.50,2454781.50,2455194.50,2455275.50] (CURRIE)
    #Sys_Dist_PC = 19.3 (CURRIE)
    #Mass1 = 1
    #Mass2 = 1
    Sep_Angle_arcsec_measured_REALs = [2.71, 2.87,2.82,1.93]
    SA_mean_errors = [0.05, 0.03, 0.04, 0.02]                                         
    PA_deg_measured_REALs = [31.3, 30.85, 33.2, 53.1]
    PA_mean_errors = [0.5, 0.03, 1.0, 0.6]
    epochs = [2451945.5, 2451143.5, 2451711.5, 2455590.5]
    Sys_Dist_PC = 15.0
    Mass1 = 1.3
    Mass2 = 0.4
    #********************************************************************
    #********************************************************************
    if not silent:
        print '\nMCMC: $$$$$$$$$$$$$$$$$$$  STARTING THE MARKOV CHAIN  $$$$$$$$$$$$$$$$$'
        print 'Number of sample orbits being created = ',numSamples
    
    # init input lists
    longAN_degs = []
    es = []
    Ts = []
    periods = []
    inclination_degs = []
    argPeri_degs = []
    
    # init output lists
    a1s2 = []
    a2s2 = []
    chiSquareds = []
    Sep_Dists2 = []
    
    ns2 = []
    Ms2 = []
    Es2 = []
    thetas2 = []
    PA_deg_measured_models2 = []
    
    # record the time the chain started
    startTime = time.clock()
    
    ## Use the basicOrbSim to get initial guess values for the model inputs for first epoch
    if not silent:
        print 'MCMC: ** Using basicOrbSim to find the initial parameter set for the chain **'
    
    (longAN_deg_initial, e_initial, T_initial,  period_initial, \
     inclination_deg_initial, argPeri_deg_initial, ns_initial, Ms_initial, Es_initial, thetas_initial, Sep_Dists_initial,\
      PA_deg_measured_models_initial, a1_initial, a2_initial,) =  \
     basicOrbSim(10000, Sep_Angle_arcsec_measured_REALs[0], SA_mean_errors[0], PA_deg_measured_REALs[0],\
      PA_mean_errors[0], chiSquaredMax, epochs[0], Sys_Dist_PC, Mass1=Mass1,Mass2=Mass2, verbose=False)
    
    # Save all the 'initial' values the their respective lists.
    longAN_degs.append(longAN_deg_initial)
    es.append(e_initial) 
    Ts.append(T_initial) 
    periods.append(period_initial)
    inclination_degs.append(inclination_deg_initial)
    argPeri_degs.append(argPeri_deg_initial)  
    
    if not silent: 
        print 'MCMC: ** Initial parameter set found. Starting chain now **'
    
    ## Make an initial chi**2 value to start chain from
    chi_squared_total_last = chiSquaredMax-0.1
    
    # variables for the success rate print block in chain loop
    printTime = numSamples/numSamplePrints
    printCount = 0
    printsDone = 0
    probRatio = 0
    
    ## Start to create samples and check if they are are accepted
    ## ie. Start the chain
    for i in range(1,numSamples+1):
        
        # block to control printing success rate to screen
        printCount = printCount + 1
        if printCount==printTime:
            printsDone = printsDone+1
            printCount = 0
            if not silent:
                print 'latest probRatio = ',probRatio
                print 'latest chi-square = ',chi_squared_total_last
            print str(len(es))+' successful samples out of '+str(i)+'. '+str(printsDone)+'/'+str(int(numSamplePrints))+' completed.'

        ## Input variables made with TRIANGULAR random numbers
        inclination_deg = rand.triangular(inclination_degMIN, inclination_degMAX, inclination_degs[len(es)-1])
        inc_prob1 = TriangProbDensFunc(inclination_degMIN, inclination_degMAX, inclination_degs[len(es)-1],inclination_deg)
        inc_prob2 = TriangProbDensFunc(inclination_degMIN, inclination_degMAX, inclination_deg, inclination_degs[len(es)-1])
        inc_prob_ratio = inc_prob1/inc_prob2
        
        longAN_deg = rand.triangular(longAN_degMIN, longAN_degMAX, longAN_degs[len(es)-1])
        long_prob1 = TriangProbDensFunc(longAN_degMIN, longAN_degMAX, longAN_degs[len(es)-1], longAN_deg)
        long_prob2 = TriangProbDensFunc(longAN_degMIN, longAN_degMAX, longAN_deg, longAN_degs[len(es)-1])
        long_prob_ratio = long_prob1/long_prob2
        
        e = rand.triangular(eMIN, eMAX, es[len(es)-1])
        e_prob1 = TriangProbDensFunc(eMIN, eMAX, es[len(es)-1], e)
        e_prob2 = TriangProbDensFunc(eMIN, eMAX, e, es[len(es)-1])
        e_prob_ratio = e_prob1/e_prob2
        
        period = rand.triangular(periodMIN, periodMAX, periods[len(es)-1]) # [yrs]
        period_prob1 = TriangProbDensFunc(periodMIN, periodMAX, periods[len(es)-1], period)
        period_prob2 = TriangProbDensFunc(periodMIN, periodMAX, period, periods[len(es)-1])
        period_prob_ratio = period_prob1/period_prob2
        
        T = rand.triangular(epochs[-1]-period*365.0, epochs[-1], Ts[len(es)-1]) # thus between a full period ago and now
        T_prob1 = TriangProbDensFunc(epochs[-1]-period*365.0, epochs[-1], Ts[len(es)-1], T)
        T_prob2 = TriangProbDensFunc(epochs[-1]-period*365.0, epochs[-1], T, Ts[len(es)-1])
        T_prob_ratio = T_prob1/T_prob2
        
        argPeri_deg = rand.triangular(argPeri_degMIN, argPeri_degMAX, argPeri_degs[len(es)-1])
        argPeri_prob1 = TriangProbDensFunc(argPeri_degMIN, argPeri_degMAX, argPeri_degs[len(es)-1], argPeri_deg)
        argPeri_prob2 = TriangProbDensFunc(argPeri_degMIN, argPeri_degMAX, argPeri_deg, argPeri_degs[len(es)-1])
        argPeri_prob_ratio = argPeri_prob1/argPeri_prob2
        
        # Calculate the average probability density ratio to use in acceptance formula
        average_prob_ratio = (inc_prob_ratio+long_prob_ratio+e_prob_ratio+period_prob_ratio+T_prob_ratio+argPeri_prob_ratio)/6.0
        
        
        (chi_squared_total_cur, ns, Ms, Es, thetas, Sep_Dists, PA_deg_measured_models, a1s, a2s) =\
        multiEpochOrbCalc(Sep_Angle_arcsec_measured_REALs, PA_deg_measured_REALs, PA_mean_errors, epochs,\
                       chiSquaredMax, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, \
                       Mass1=Mass1, Mass2=Mass2, verbose=verbose)
            
        if chi_squared_total_cur==0:
            #print 'chi_squared_total_cur was zero for sample '+str(i)
            silly = True
        else:
            ##### In my original model, I would calc the gauss dist in each direction, BUT
            ##### as the gauss dist is symmetric, the fraction of the two will ALWAYS = 1.
            ##### Thus I am only using a ratio of the probability densities (ie. a likelihood ratio).
            probRatio = chi_squared_total_last/chi_squared_total_cur
            alpha = rand.uniform(0.0,1.0)#rand.triangular(0.0,1.0,0.7)# maybe set this to 1.0 !?!?!?!?!? $$$$$$$$$$$$$
            
            
            if alpha <= (average_prob_ratio*probRatio):
                # keep the NEW value, ie. sample accepted!
                chi_squared_total_last = chi_squared_total_cur
                
                # store input values to the model
                longAN_degs.append(longAN_deg)
                es.append(e)
                Ts.append(T)
                periods.append(period)
                inclination_degs.append(inclination_deg)
                argPeri_degs.append(argPeri_deg)
                # store output orbital elements 
                a1s2.append(a1s)
                a2s2.append(a2s)
                chiSquareds.append(chi_squared_total_last)
                Sep_Dists2.append(Sep_Dists)
                ns2.append(ns)
                Ms2.append(Ms)
                Es2.append(Es)
                thetas2.append(thetas)
                Sep_Dists2.append(Sep_Dists)
                PA_deg_measured_models2.append(PA_deg_measured_models)
            
    if not silent:
        print 'MCMC: $$$$$$$$$$ Markov Chain complete $$$$$$$$$$'
    
    # record the time the chain finished
    endTime = time.clock()
    totalTime = (endTime-startTime)/60
    if not silent:
        print 'MCMC: The chain took '+str(totalTime)+' minutes to complete.\n'
        
    print 'MCMC for '+titleMod+' took '+str(totalTime)+' minutes to complete.\n' ##### only print in 'silent' mode to track time on night runs
    ## write resultant data to files
    dataWriter(plotFileTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                  ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds)
        
    ## plot results
    orbElementPlotter(plotFileTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
          ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds, showPlots=showPlots, summaryOnly=False, verbose=verbose)
        
    if not silent:
        print '\nMCMC: $$$$$$$$$$$$$ FINISHED MCMCtriangular2 Sim $$$$$$$$$$$$$$$\n'
    
    #$$$$$$$$$ DELETING ALL OBJECTS AND RUNNING GARBAGE COLLECTION TO FREE UP MEMORY $$$$$$$$$$$$
    #$$$$$$$$$ THIS STILL ISN'T HELPING MY RAM ISSUE, BUT IT DOESN'T HURT SO...      $$$$$$$$$$$$
    del plotFileTitle
    del longAN_degs
    del es
    del Ts
    del periods
    del inclination_degs
    del argPeri_degs
    del ns2
    del Ms2
    del Es2
    del thetas2
    del Sep_Dists2
    del a1s2
    del a2s2
    del chiSquareds    
    
    #****** DONE!! ********






