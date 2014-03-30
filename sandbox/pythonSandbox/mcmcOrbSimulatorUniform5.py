# @Author: Kyle Mede, kylemede@gmail.com, written while at the University of Tokyo 2011.

""" This is a Markov Chain Monte Carlo simulator using a Metropolis-Hastings algorithm
    to find the orbital parameters for either a planetary or binary star system 
    with provided separation angle, system distance and position angle.
    
    NOTE:
    This version runs on UNIFORM random numbers within a variance range to perturb the current 
    parameter values to a new proposal state, rather than simply random numbers for the parameters
     alone.  The variance value will be tuned periodically to attempt an optimal acceptance rate of
     25-35%.  This version also uses the orbitCalculator from orbitToolbox1 rather than version 2.
    The acceptance will in the Metropolis-Hastings algorithm will also be based on a function that is
    half a normalized gaussian centered at 1.0 and ignoring all values less than that; the variance of
    the gaussian will be set so the wings die out around a chiSquared of ~100 and might be changed later
    if this is found to be a poor choice. The initial values of the semi-random input parameters to 
    the model will also be chosen, rather than generated with the basicOrbSim.
    
    !! REQUIRES AT LEAST 4 EPOCHS OF OBSERVATIONS !!
    
    !!!! THIS VERSION IS JUST TO FIND THE LISTS OF ACCEPTED PARAMETERS and returns them, NO FILE WRITING OR PLOTTING. !!!!
"""

import math
import numpy as np
import random as rand
import time
from orbitToolbox import multiEpochOrbCalc, gaussianDist

def mcmcUniformOrbSim(numSamples, silent=True):
    
    #*********** In this brick is all you need to play with *************    
    #********************************************************************
    chiSquaredMax = 2000.0 
    verbose = False                                                               
    numSamplePrints = 10.0        
    targetVariance = 500.0 # 20 has the bell curve dying out around 50
    tuningSamples = 1000
    sigmaMAXpercent = 50.0 
    #### Ranges for acceptable random number inputs ######                                        
    longAN_degMIN = 0.001 # [deg]                                                        
    longAN_degMAX = 179.999 # [deg]                                                        
    eMIN = 0.001                                                                        
    eMAX = 0.999                                                                        
    periodMIN = 1.0 # [yrs]                                                                
    periodMAX = 50.0 # [yrs]                                                            
    inclination_degMIN = 0.001 # [deg]                                                    
    inclination_degMAX = 179.999 # [deg]                                                
    argPeri_degMIN = 0.001 # [deg]                                                        
    argPeri_degMAX = 179.999#89.999 # [deg]          
    a_totalMIN = 0.1
    a_totalMAX = 10.0                                          
    #### MEASURED VALUES ########################                                        
    ### Currie, betaPic data
    #SA_arcsec_measured_REALs = [0.411,0.210,0.326,0.345] # used to calc kai^2 (ACTUALLY NOT IN THIS VERSION, ONLY VERISON 1)(CURRIE)
    #SA_mean_errors = [0.008,0.027,0.013,0.012] (CURRIE)
    #PA_deg_measured_REALs = [31.7,211.49,210.64,209.8] # used to calc kai^2 (CURRIE)
    #PA_mean_errors = [1.3,1.9, 1.2,0.8] (CURRIE)
    #epochs = [2452953.50,2454781.50,2455194.50,2455275.50] (CURRIE)
    #Sys_Dist_PC = 19.3 (CURRIE)
    #Mass1 = 1
    #Mass2 = 1
    #### Liu, 2MAS J1534... data
    SA_arcsec_measured_REALs = [0.0628,0.2113,0.199,0.1912,0.1906,0.1580,0.1537,0.1144,0.1020]
    SA_mean_errors = [0.0012,0.0015,0.0011,0.0011,0.0003,0.0006,0.0004,0.0011,0.0004] # mean simply implies the mean of the + and - uncertainties                                         
    PA_deg_measured_REALs = [357.1,14.1,14.5,15.5,15.43,17.5,15.53,21.5,20.4]
    PA_mean_errors = [0.8,0.3,0.6,0.4,0.12,0.2,0.13,0.9,1.5] # mean simply implies the mean of the + and - uncertainties 
    epochs = [2451774.5,2453491.5,2453754.5,2453826.5,2453860.5,2454185.5,2454212.5,2454480.5,2454557.5]
    Sys_Dist_PC = 13.5
    Mass1 = 1
    Mass2 = 1
    ## Initial values that were found in the past #### maybe create a better method than basicOrbSim to make these.
    longAN_degInitial = 178.0
    eInitial = 0.24
    TInitial = 2456024.5
    periodInitial = 15.2
    inclination_degInitial = 84.3
    argPeri_degInitial = 13.0
    a_totalInitial = 2.3
    #### tauBoo data
    #SA_arcsec_measured_REALs = [2.71, 2.87,2.82,1.93]
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
        print '\nMCMC: $$$$$$$$$$$$$$$$$$$  STARTING THE MARKOV CHAIN  $$$$$$$$$$$$$$$$$'
        print 'Number of sample orbits being created = ',numSamples
    
    # init input lists
    longAN_degs = []
    es = []
    Ts = []
    periods = []
    inclination_degs = []
    argPeri_degs = []
    a_totals = []
    
    # init output lists
    a1s2 = []
    a2s2 = []
    chiSquareds = []
    Sep_Dists2 = []
    ns2 = []
    Ms2 = []
    Es2 = []
    thetas2 = []
    SA_deg_measured_models2 = []
    PA_deg_measured_models2 = []
    
    # Load up input lists with initial values set at the top
    longAN_degs.append(longAN_degInitial)
    es.append(eInitial)
    Ts.append(TInitial)
    periods.append(periodInitial)
    inclination_degs.append(inclination_degInitial)
    argPeri_degs.append(argPeri_degInitial) 
    a_totals.append(a_totalInitial)
    
    # load up output lists with zeros for their first element
    # to match the length of the input lists
    numEpochs = len(epochs)
    a1s2.append(np.zeros((0,numEpochs)))
    a2s2.append(np.zeros((0,numEpochs)))
    chiSquareds.append(0)
    Sep_Dists2.append(np.zeros((0,numEpochs)))
    ns2.append(np.zeros((0,numEpochs)))
    Ms2.append(np.zeros((0,numEpochs)))
    Es2.append(np.zeros((0,numEpochs)))
    thetas2.append(np.zeros((0,numEpochs)))
    SA_deg_measured_models2.append(np.zeros((0,numEpochs)))
    PA_deg_measured_models2.append(np.zeros((0,numEpochs)))
    
    # record the time the chain started
    startTime = time.clock()
    
    # variables for the success rate print block in chain loop
    printTime = numSamples/numSamplePrints
    printCount = 0
    printsDone = 0
    probRatio = 0
    probLast = 0.5 # just a start value
    
    ## Record the initial sigma values to start chain from;
    #  these will be tuned throughout the chain.
    sigmas = []
    sigmas.append(1.0)
    timesTuned = 0
    acceptedCounter = 0
    totalCounter = 0
    roundsTuned = 0
    reduced_chi_squared_total_cur=0 # just an initial value for printing, will be replaced to real value below
    
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
                print 'latest probRatio = ',probRatio
                print 'latest chi-square = ',reduced_chi_squared_total_cur


        ## propose a new set of input parameters by perturbing the current ones 
        
        # inclination
        # sigma is a percentage of total range, so convert its current percentage to a 
        # range value
        sigma = (sigmas[-1]/100.0)*(inclination_degMAX-inclination_degMIN)
        if (inclination_degs[-1]+sigma)>=inclination_degMAX:
            max = inclination_degMAX
            min = inclination_degMAX-(2.0*sigma)
        elif (inclination_degs[-1]-sigma)<=inclination_degMIN:
            max = inclination_degMIN+(2.0*sigma)
            min = inclination_degMIN
        else:
            max = inclination_degs[-1] + sigma
            min = inclination_degs[-1] - sigma
        inclination_deg = rand.uniform(min, max)
        
        # long of acending node
        # sigma is a percentage of total range, so convert its current percentage to a 
        # range value
        sigma = (sigmas[-1]/100.0)*(longAN_degMAX-longAN_degMIN)
        if (longAN_degs[-1]+sigma)>=longAN_degMAX:
            max = longAN_degMAX
            min = longAN_degMAX-(2.0*sigma)
        elif (longAN_degs[-1]-sigma)<=longAN_degMIN:
            max = longAN_degMIN+(2.0*sigma)
            min = longAN_degMIN
        else:
            max = longAN_degs[-1] + sigma
            min = longAN_degs[-1] - sigmas[-1]
        longAN_deg = rand.uniform(min, max)
        
        #eccentricity
        # sigma is a percentage of total range, so convert its current percentage to a 
        # range value
        sigma = (sigmas[-1]/100.0)*(eMAX-eMIN)
        if (es[-1]+sigma)>=eMAX:
            max = eMAX
            min = eMAX-(2.0*sigma)
        elif (es[-1]-sigma)<=eMIN:
            max = eMIN+(2.0*sigma)
            min = eMIN
        else:
            max = es[-1] + sigma
            min = es[-1] - sigma
        e = rand.uniform(min, max)
        if e>1.0:
            print 'e created greater than 1.0, = '+str(e)
            print 'max = '+str(max)
            print 'min = '+str(min)
            print 'sigma = '+str(sigma)
        
        #period
        # sigma is a percentage of total range, so convert its current percentage to a 
        # range value
        sigma = (sigmas[-1]/100.0)*(periodMAX-periodMIN)
        if (periods[-1]+sigma)>=periodMAX:
            max = periodMAX
            min = periodMAX-(2.0*sigma)
        elif (periods[-1]-sigma)<=periodMIN:
            max = periodMIN+(2.0*sigma)
            min = periodMIN
        else:
            max = periods[-1] + sigma
            min = periods[-1] - sigma
        period = rand.uniform(min, max)
        
        # time of last periapsis
        TMAX = epochs[-1]-period*365.0
        TMIN = epochs[-1]
        # sigma is a percentage of total range, so convert its current percentage to a 
        # range value
        sigma = (sigmas[-1]/100.0)*(TMAX-TMIN)
        if (Ts[-1]+sigma)>=TMAX:
            max = TMAX
            min = TMAX-(2.0*sigma)
        elif (Ts[-1]-sigma)<=TMIN:
            max = TMIN+(2.0*sigma)
            min = TMIN
        else:
            max = Ts[-1] + sigma
            min = Ts[-1] - sigma
        T = rand.uniform(min, max)
        
        # argument of perigie
        # sigma is a percentage of total range, so convert its current percentage to a 
        # range value
        sigma = (sigmas[-1]/100.0)*(argPeri_degMAX-argPeri_degMIN)
        if (argPeri_degs[-1]+sigma)>=argPeri_degMAX:
            max = argPeri_degMAX
            min = argPeri_degMAX-(2.0*sigma)
        elif (argPeri_degs[-1]-sigma)>=argPeri_degMIN:
            max = argPeri_degMIN+(2.0*sigma)
            min = argPeri_degMIN
        else:
            max = argPeri_degs[-1] + sigma
            min = argPeri_degs[-1] - sigma
        argPeri_deg = rand.uniform(min, max)
        
        # Total semi-major axis
        # sigma is a percentage of total range, so convert its current percentage to a 
        # range value
        sigma = (sigmas[-1]/100.0)*(a_totalMAX-a_totalMIN)
        if (a_totals[-1]+sigma)>=a_totalMAX:
            max = a_totalMAX
            min = a_totalMAX-(2.0*sigma)
        elif (a_totals[-1]-sigma)>=a_totalMIN:
            max = a_totalMIN+(2.0*sigma)
            min = a_totalMIN
        else:
            max = a_totals[-1] + sigma
            min = a_totals[-1] - sigma
        a_total = rand.uniform(min, max)

        ## send random parameters along with known ones to multi-epoch orbit calculator
        (chi_squared_total_cur, ns, Ms, Es, thetas, Sep_Dists, SA_deg_measured_models, PA_deg_measured_models, a1s, a2s) =\
        multiEpochOrbCalc(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors, chiSquaredMax,\
                       epochs, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total,\
                       Mass1=Mass1, Mass2=Mass2, verbose=verbose)

        # check if a useful set of values were returned
        if chi_squared_total_cur<1.0:

            # or it was less than 1.0 which is not realistic
            skip=True  
        else:
            # values returned were useful so continue
            
            ## Apply the Metropolis-Hastings Algorithm ##
            
            # normalize the chiSquared to give the reduce chiSquared value
            reduced_chi_squared_total_cur = (1/((2*(numEpochs))-6.0))*chi_squared_total_cur
            # get current probability from the target distribution (ie half a gaussian)
            # actually a whole gaussian, but since chiSquared is always positive,
            # it will only sample the positive half.
            probCur = gaussianDist(reduced_chi_squared_total_cur, 1.0, targetVariance)
            # calculate probability ratio (or likelihood ratio)
            probRatio = probCur/probLast
            # find min of probRatio and 1
            if probRatio<1.0:
                rightSide = 1
            else:
                rightSide = probRatio
            #generate a random number between 0 and 1
            alpha = rand.uniform(0.0, 1.0)
            
            if alpha <= rightSide:
                # keep the NEW value, ie. sample accepted!
                print 'probCur '+repr(probCur)
                print 'probLast '+repr(probLast)
                probLast = probCur
                acceptedCounter = acceptedCounter + 1
                
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
                chiSquareds.append(reduced_chi_squared_total_cur)
                ns2.append(ns) #not sure why I was saving this as it isn't used after this??
                Ms2.append(Ms) #not sure why I was saving this as it isn't used after this??
                Es2.append(Es) #not sure why I was saving this as it isn't used after this??
                thetas2.append(thetas) #not sure why I was saving this as it isn't used after this??
                Sep_Dists2.append(Sep_Dists) #not sure why I was saving this as it isn't used after this??
                SA_deg_measured_models2.append(SA_deg_measured_models) #not sure why I was saving this as it isn't used after this??
                PA_deg_measured_models2.append(PA_deg_measured_models)  #not sure why I was saving this as it isn't used after this??
                
            # increment tuning counter
            totalCounter = totalCounter + 1
        
            if totalCounter==tuningSamples:
                # it is time to recalculate the acceptance rate and adjust current sigma accordingly 
                acceptRate = float(acceptedCounter)/float(tuningSamples)
                acceptedCounter = 0
                totalCounter = 0
                timesTuned = timesTuned + 1 # not used, but maybe useful later
                
#                if sigmas[-1]<=sigmaMAXpercent:
#                    if acceptRate<0.25:
#                        # ie. acceptance rate is too high, so raise current sigma being adjusted
#                        sigmas.append(sigmas[-1] + 0.2)
#                    elif acceptRate>0.35:
#                        # ie. acceptance rate is too low, so lower current sigma being adjusted
#                        if (sigmas[-1]-0.2)<0.0:
#                            # all ready below close t 0, so leave it
#                            sigmas.append(sigmas[-1])
#                        else:   
#                            sigmas.append(sigmas[-1] - 0.2)
#                else:
#                    # sigma value is all ready at max, so leave it.
#                    sigmas.append(sigmas[-1])
                sigmas.append(sigmas[-1])
                #print 'Finished a round of sigma tuning' #$$$$$$$$
                
        # DONE sample loop
        
    if not silent:
        print 'MCMC: $$$$$$$$$$ Markov Chain complete $$$$$$$$$$'
    
    # record the time the chain finished and print
    endTime = time.clock()
    totalTime = (endTime-startTime) # in seconds
    if totalTime<60:
        print 'MCMC: The took '+str(int(totalTime))+' seconds to complete.\n'  ##### only print in 'silent' mode to track time
    elif totalTime>3600:
        totalTime = totalTime/3600.0
        print 'MCMC: The took '+str(int(totalTime))+' hours and '+str(int(60*((totalTime)-int(totalTime))))+' minutes to complete.\n'  ##### only print in 'silent' mode to track time
    else:
        totalTime = totalTime/60.0
        print 'MCMC: The took '+str(int(totalTime))+' minutes and '+str(int(60*((totalTime)-int(totalTime))))+' seconds to complete.\n'  ##### only print in 'silent' mode to track time

    if not silent:
        print '\nMCMC: $$$$$$$$$$$$$ FINISHED MCMCuniform4 Sim $$$$$$$$$$$$$$$\n'
  
    return (longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds)
    #****** DONE!! ********
                



