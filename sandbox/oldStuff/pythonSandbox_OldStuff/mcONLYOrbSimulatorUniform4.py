#@Author: Kyle Mede, kylemede@gmail.com, written for the University of Tokyo during UTRIP 2011.
""" This is a Markov Chain Monte Carlo simulator using a Metropolis-Hastings algorithm
    to find the orbital parameters for either a planetary or binary star system 
    with provided separation angle, system distance and position angle.
    This version uses UNIFORM random variables as inputs.
    
    For a complete description of how this simulator works, please see mcmcOrbSimulatorUniform2.pdf document.
    
    !!!! THIS VERSION IS JUST TO FIND THE LISTS OF ACCEPTED PARAMETERS and returns them, NO FILE WRITING OR PLOTTING. !!!!
    """

import math
import numpy as np
import random as rand
import time
from orbitToolbox2 import *
from basicOrbSimulator2 import basicOrbSim

def mcONLYUniformOrbSim(numSamples, silent=True):

    #*********** In this brick is all you need to play with *************    
    #********************************************************************
    chiSquaredMax = 500.0 
    verbose = False                                                               
    numSamplePrints = 10.0        
    
    #### Ranges for acceptable random number inputs ######                                        
    longAN_degMIN = 0.001 # [deg]                                                        
    longAN_degMAX = 179.999 # [deg]                                                        
    eMIN = 0.001                                                                        
    eMAX = 0.999                                                                        
    periodMIN = 1.0 # [yrs]                                                                
    periodMAX = 30.0 # [yrs]                                                            
    inclination_degMIN = 0.001 # [deg]                                                    
    inclination_degMAX = 179.999 # [deg]                                                
    argPeri_degMIN = 0.001 # [deg]                                                        
    argPeri_degMAX = 359.999 #179.999#89.999 # [deg]                                                    
    #### MEASURED VALUES ########################                                        
    ### Currie, betaPic data
    #Sep_Angle_arcsec_measured_REALs = [0.411,0.210,0.326,0.345] # used to calc kai^2 (ACTUALLY NOT IN THIS VERSION, ONLY VERISON 1)(CURRIE)
    #SA_mean_errors = [0.008,0.027,0.013,0.012] (CURRIE)
    #PA_deg_measured_REALs = [31.7,211.49,210.64,209.8] # used to calc kai^2 (CURRIE)
    #PA_mean_errors = [1.3,1.9, 1.2,0.8] (CURRIE)
    #epochs = [2452953.50,2454781.50,2455194.50,2455275.50] (CURRIE)
    #Sys_Dist_PC = 19.3 (CURRIE)
    #Mass1 = 1
    #Mass2 = 1
    #### Liu, 2MAS J1534... data
    Sep_Angle_arcsec_measured_REALs = [0.0628,0.2113,0.199,0.1906,0.1537,0.1144,0.1020]
    SA_mean_errors = [0.0012,0.0015,0.0011,0.0003,0.0004,0.0011,0.0004]                                         
    PA_deg_measured_REALs = [357.1,14.1,14.5,15.43,15.53,21.5,20.4]
    PA_mean_errors = [0.8,0.3,0.6,0.12,0.13,0.9,1.5]
    epochs = [2451774.5,2453491.5,2453754.5,2453860.5,2454212.5,2454480.5,2454557.5]
    Sys_Dist_PC = 13.5
    Mass1 = 1
    Mass2 = 1
    #### tauBoo data
    #Sep_Angle_arcsec_measured_REALs = [2.71, 2.87,2.82,1.93]
    #SA_mean_errors = [0.05, 0.03, 0.04, 0.02]                                         
    #PA_deg_measured_REALs = [31.3, 30.85, 33.2, 53.1]
    #PA_mean_errors = [0.5, 0.03, 1.0, 0.6]
    #epochs = [2451945.5, 2451143.5, 2451711.5, 2455590.5]
    #Sys_Dist_PC = 15.0
    #Mass1 = 1.3
    #Mass2 = 0.4
    #********************************************************************
    #********************************************************************
    if not silent:
        print '\nsimulator: $$$$$$$$$$$$$$$$$$$  STARTING THE MC ONLY sim  $$$$$$$$$$$$$$$$$'
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
        print 'simulator: ** Using basicOrbSim to find the initial parameter set for the chain **'
    
    (longAN_deg_initial, e_initial, T_initial,  period_initial, \
     inclination_deg_initial, argPeri_deg_initial, ns_initial, Ms_initial, Es_initial, thetas_initial, Sep_Dists_initial,\
      PA_deg_measured_models_initial, a1_initial, a2_initial,) =  \
     basicOrbSim(10000, Sep_Angle_arcsec_measured_REALs[0], SA_mean_errors[0], PA_deg_measured_REALs[0],\
      PA_mean_errors[0], chiSquaredMax, epochs[0], Sys_Dist_PC, Mass1=Mass1,Mass2=Mass2, verbose=False)
     
    if not silent: 
        print 'simulator: ** Initial parameter set found. Starting chain now **'
    
    # variables for the success rate print block in chain loop
    printTime = numSamples/numSamplePrints
    printCount = 0
    printsDone = 0
    
    ## Start to create samples and check if they are are accepted
    ## ie. Start the chain
    for i in range(1,numSamples+1):
        
        # block to control printing success rate to screen
        printCount = printCount + 1
        if printCount==printTime:
            printsDone = printsDone+1
            printCount = 0
            if not silent:
                print str(len(es))+' successful samples out of '+str(i)+'. '+str(printsDone)+'/'+str(int(numSamplePrints))+' completed.'
                print 'latest chi-square = ',chi_squared_total_last
                
        ## Input variables made with uniform random numbers
        #  Note: These will be stored if they satisfy the chi**2 max
        #        so only the good ones are stored to save processor time
        inclination_deg = rand.uniform(inclination_degMIN, inclination_degMAX)
        longAN_deg = rand.uniform(longAN_degMIN, longAN_degMAX)
        e = rand.uniform(eMIN, eMAX)
        period = rand.uniform(periodMIN, periodMAX) # [yrs]
        T = rand.uniform(epochs[-1]-period*365.0, epochs[-1]) # thus between a full period ago and now
        argPeri_deg = rand.uniform(argPeri_degMIN, argPeri_degMAX)
        
        (chi_squared_total_cur, ns, Ms, Es, thetas, Sep_Dists, PA_deg_measured_models, a1s, a2s) =\
        multiEpochOrbCalc(Sep_Angle_arcsec_measured_REALs, PA_deg_measured_REALs, PA_mean_errors, epochs,\
                       chiSquaredMax, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, \
                       Mass1=Mass1, Mass2=Mass2, verbose=verbose)
            
        if chi_squared_total_cur==0:
            #print 'chi_squared_total_cur was zero for sample '+str(i)
            silly = True
        else:
           
            if chiSquaredMax >= chi_squared_total_cur:
                # keep the NEW value, ie. sample accepted!
        
                
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
                chiSquareds.append(chi_squared_total_cur)
                Sep_Dists2.append(Sep_Dists)
                ns2.append(ns)
                Ms2.append(Ms)
                Es2.append(Es)
                thetas2.append(thetas)
                Sep_Dists2.append(Sep_Dists)
                PA_deg_measured_models2.append(PA_deg_measured_models)  #not sure why I was saving this as it isn't used after this??
            
    if not silent:
        print 'simulator: $$$$$$$$$$ Monte Carlo simulator complete $$$$$$$$$$'
    
    # record the time the chain finished and print
    endTime = time.clock()
    totalTime = (endTime-startTime) # in seconds
    if totalTime<60:
        print 'simulator: The took '+str(int(totalTime))+' seconds to complete.\n'  ##### only print in 'silent' mode to track time
    elif totalTime>3600:
        totalTime = totalTime/3600.0
        print 'simulator: The took '+str(int(totalTime))+' hours and '+str(int(60*((totalTime)-int(totalTime))))+' minutes to complete.\n'  ##### only print in 'silent' mode to track time
    else:
        totalTime = totalTime/60.0
        print 'simulator: The took '+str(int(totalTime))+' minutes and '+str(int(60*((totalTime)-int(totalTime))))+' seconds to complete.\n'  ##### only print in 'silent' mode to track time

    if not silent:
        print '\nsimulator: $$$$$$$$$$$$$ FINISHED MCONLYuniform4 Sim $$$$$$$$$$$$$$$\n'
  
    return (longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds)
    #****** DONE!! ********






