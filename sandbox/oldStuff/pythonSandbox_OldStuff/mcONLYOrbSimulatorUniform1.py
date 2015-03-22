# @Author: Kyle Mede, kylemede@gmail.com, written while at the University of Tokyo 2011-2013.
"""
This is a Non-Markov Chain version of mcmcOrbSimulatorUniform5.py.
It uses uniform random number inputs and simple chiSquare<chiSquareMax acceptance.
"""

import math
import numpy as np
import random as rand
import time
from orbitToolbox import multiEpochOrbCalc

def mcONLYUniformOrbSim(numSamples, silent=True, lowRAM=True):
    
     #*********** In this brick is all you need to play with *************    
    #********************************************************************
    chiSquaredMax = 20
    verbose = True                                                               
    numSamplePrints = 20.0        
    #### Ranges for acceptable random number inputs ######                                        
    longAN_degMIN = 164.0 # [deg]                                                        
    longAN_degMAX = 167.5 # [deg]                                                        
    eMIN = 0.60                                                                       
    eMAX = 0.725                                                                      
    periodMIN = 380.0 # [yrs]                                                                
    periodMAX = 505 # [yrs]                                                            
    inclination_degMIN = 52.5 # [deg]                                                    
    inclination_degMAX = 67.0 # [deg]                                                
    argPeri_degMIN = 351.0 # [deg]                                                        
    argPeri_degMAX = 373.0#89.999 # [deg]          
    a_totalMIN = 86.0
    a_totalMAX = 92.0                                      
    #### MEASURED VALUES ########################                                        
    
    #### TEST system data (HD 130948BC Dupuy2009)
#    SA_arcsec_measured_REALs = [0.0993, 0.0946, 0.0961, 0.0978, 0.0568, 0.0564, 0.0573, 0.1117, 0.1091, 0.0979, 0.0722, 0.0574, 0.0573, 0.0584, 0.0517]
#    SA_mean_errors = [0.0017, 0.0017, 0.0023, 0.0023, 0.0009, 0.0012, 0.0009, 0.0015, 0.0011, 0.0009, 0.0006, 0.0005, 0.0009, 0.0010, 0.0006] # mean simply implies the mean of the + and - uncertainties                                         
#    PA_deg_measured_REALs = [307.9,306.9,308.9,307.0,144.9,146.4,148.6,132.6,132.28,130.73,127.6,124.7,124.6,124.1,123.9]
#    PA_mean_errors = [1.1,1.0,1.6,1.4,0.5,0.6,0.6,0.4,0.13,0.17,0.3,0.4,0.6,0.7,0.5]# mean simply implies the mean of the + and - uncertainties 
#    epochs = [2452519.702,2452519.706,2452519.709,2452519.719,2453425.152,2453425.219,2453425.282,2454127.085,2454185.028,2454306.737,2454481.192,2454555.068,2454555.09,2454555.157,2454584.106]
#    Sys_Dist_PC = 18.17
#    Mass1 = 1
#    Mass2 = 1
    #### LP 349-25AB
#    SA_arcsec_measured_REALs = [0.125, 0.107, 0.1056, 0.1166, 0.1189, 0.1233, 0.1282, 0.1373, 0.1152, 0.1029, 0.10245, 0.0711, 0.0834, 0.1126]
#    SA_mean_errors = [0.01, 0.01, 0.0003, 0.0006, 0.0006, 0.0004, 0.0007, 0.0004, 0.0005, 0.0002, 0.00019, 0.0003, 0.0003, 0.0004]
#    PA_deg_measured_REALs = [12.7, 7.1, 247.3, 240.4, 238.9, 236.0, 232.7, 207.92, 194.64, 189.71, 187.46, 98.3, 81.6, 59.76]
#    PA_mean_errors = [2.0, 0.5, 0.2, 0.2, 0.3, 0.2, 0.3, 0.09, 0.13, 0.11, 0.1, 0.6, 0.4, 0.13]
#    epochs = [2453190.5, 2453275.5, 2453930.5, 2454000.5, 2454018.5, 2454051.5, 2454094.5, 2454482.5, 2454648.5, 2454699.5, 2454719.5, 2455103.5, 2455181.5, 2455339.5]
#    Sys_Dist_PC = 13.2
#    Mass1 = 1
#    Mass2 = 1
    #### tauBoo data
    ######## ALL but P&P  #########################
#    SA_arcsec_measured_REALs = [3.401, 2.865, 2.82, 2.71, 2.181, 1.9421, 1.9324]
#    SA_mean_errors = [0.06, 0.04, 0.05, 0.03, 0.02, 0.01, 0.01]                                         
#    PA_deg_measured_REALs = [20.65, 30.85, 33.2, 31.3,46.26, 52.65, 53.68]
#    PA_mean_errors = [0.5, 0.4, 1.0, 1.0, 0.42, 0.16, 0.07]
#    epochs = [2448337.5, 2451143.5, 2451711.5, 2451945.5, 2454579.5, 2455589.5, 2455704.5]
    ######## ALL ####################
    SA_arcsec_measured_REALs = [10.26, 9.03, 8.91, 8.71, 8.9, 8.89, 8.72, 8.41, 8.45, 8.54, 8.66, 7.28, 7.35, 6.96, 7.06, 6.4, 6.15, 6.09, 5.87, 5.4, 5.32, 5.14, 4.94, 5.17,      3.401, 2.865, 2.82, 2.71, 2.181, 1.9421, 1.9324]
    SA_mean_errors = [1.46, 0.12, 0.04, 0.2, 0.08, 0.1, 0.05, 0.22, 0.13, 0.04, 0.22, 0.89, 0.37, 0.66, 0.37, 0.37, 0.45, 0.42, 0.25, 0.33, 0.22, 0.07, 0.03, 0.38,         0.06, 0.04, 0.05, 0.03, 0.02, 0.01, 0.01]                                         
    PA_deg_measured_REALs = [347.8, 348.7, 350.8, 351.9, 354.3, 352.4, 353.6, 356, 355.1, 354.8, 355.8, 356.8, 359.1, 359.5, 1.7, 3.5, 5.1, 5.1, 6, 6.8, 7.8, 11.4, 8, 13.5,               20.65, 30.85, 33.2, 31.3,46.26, 52.65, 53.68]
    PA_mean_errors = [0.34, 1.88, 0.63, 0.17, 1.79, 0.38, 0.31, 1.43, 0.19, 0.16, 0.45, 0.21, 0.1, 0.09, 1.51, 0.93, 1.34, 1.04, 0.28, 0.75, 0.68, 1.2, 3.69, 1.08,           0.5, 0.4, 1.0, 1.0, 0.42, 0.16, 0.07]
    epochs = [2396557.5, 2403065.5, 2405359.5, 2407064.5, 2408219.5, 24708927.5, 2411846.5, 2413509.5, 2414351.5, 2414465.5, 2415404.5, 2419263.5, 2423443.5, 2424260.5, 2425728.5, 2429753.5, 2431535.5, 2431959.5, 2434136.5, 2436256.5, 2437227.5, 2438866.5, 2440132.5, 2440707.5,      2448337.5, 2451143.5, 2451711.5, 2451945.5, 2454579.5, 2455589.5, 2455704.5]
    Sys_Dist_PC = 15.62
    Mass1 = 1 #1.3
    Mass2 = 1 #0.4
    #********************************************************************
    #********************************************************************
    
    if not silent:
        print '\nMCONLY: $$$$$$$$$$$$$$$$$$$  STARTING THE SIMULATOR  $$$$$$$$$$$$$$$$$'
        numProcesses = 1
        numSamplesString = str(numSamples/int(1e6))+'-million'
        print 'Number of sample orbits being created = '+numSamplesString
    
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
    if lowRAM==False:
        Sep_Dists2 = []
        ns2 = []
        Ms2 = []
        Es2 = []
        thetas2 = []
        SA_deg_measured_models2 = []
        PA_deg_measured_models2 = []
    
    numEpochs = len(epochs)
    
     # variables for the success rate print block in chain loop
    printTime = numSamples/numSamplePrints
    printCount = 0
    printsDone = 0
    acceptedCounter = 0
    chiSquareMin = chiSquaredMax
    bestOrbit = 0
    
     # record the time the chain started
    startTime = time.clock()
    
    ## Start to create samples and check if they are are accepted
    ## ie. Start the chain
    for i in range(1,numSamples+1):
        
        # block to control printing success rate to screen
        printCount = printCount + 1
        if printCount==printTime:
            printsDone = printsDone+1
            printCount = 0
            if not silent:
                timeStr = time.strftime("%H:%M:%S", time.localtime())
                chiSquareMinStr = ". Lowest chiSquare so far = "+str(chiSquareMin)
                print str(len(es))+' successful out of '+str(i)+'. '+str(printsDone)+'/'+str(int(numSamplePrints))+' completed at '+timeStr+chiSquareMinStr
                #print 'latest chi-square = ',reduced_chi_squared_total_cur
    
    
        ## Input variables made with uniform random numbers
        #  Note: These will be stored if they satisfy the chi**2 max
        #        so only the good ones are stored to save processor time
        inclination_deg = rand.uniform(inclination_degMIN,inclination_degMAX)
        longAN_deg = rand.uniform(longAN_degMIN,longAN_degMAX)
        argPeri_deg = rand.uniform(argPeri_degMIN,argPeri_degMAX)
        e = rand.uniform(eMIN,eMAX)
        period = rand.uniform(periodMIN,periodMAX) # [yrs]
        T = rand.uniform(epochs[0]-period*365.0,epochs[0]) # thus between a full period before first observation and the time of first observation
        a_total = rand.uniform(a_totalMIN,a_totalMAX)    
    
        ## send random parameters along with known ones to multi-epoch orbit calculator
        (chi_squared_total_cur, ns, Ms, Es, thetas, Sep_Dists, SA_deg_measured_models, PA_deg_measured_models, a1s, a2s) =\
        multiEpochOrbCalc(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors,\
                       epochs, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total,\
                       Mass1=Mass1, Mass2=Mass2, verbose=False)
    
        # normalize the chiSquared to give the reduce chiSquared value
        reduced_chi_squared_total_cur = (1.0/((2.0*(numEpochs))-7.0))*chi_squared_total_cur
    
        if reduced_chi_squared_total_cur<chiSquaredMax:
            
            acceptedCounter = acceptedCounter + 1
            
            if reduced_chi_squared_total_cur<chiSquareMin:
                chiSquareMin = reduced_chi_squared_total_cur
                bestOrbit = acceptedCounter-1
            
            # store input values to the model
            longAN_degs.append(longAN_deg)
            es.append(e)
            Ts.append(T)
            periods.append(period)
            inclination_degs.append(inclination_deg)
            argPeri_degs.append(argPeri_deg)
            
            # store output orbital elements of model
            a1s2.append(a1s)
            a2s2.append(a2s)
            chiSquareds.append(reduced_chi_squared_total_cur)
            if lowRAM==False:
                ns2.append(ns) #not sure why I was saving this as it isn't used after this??
                Ms2.append(Ms) #not sure why I was saving this as it isn't used after this??
                Es2.append(Es) #not sure why I was saving this as it isn't used after this??
                thetas2.append(thetas) #not sure why I was saving this as it isn't used after this??
                Sep_Dists2.append(Sep_Dists) #not sure why I was saving this as it isn't used after this??
                SA_deg_measured_models2.append(SA_deg_measured_models) #not sure why I was saving this as it isn't used after this??
                PA_deg_measured_models2.append(PA_deg_measured_models)  #not sure why I was saving this as it isn't used after this??
            
    # DONE sample loop
    print 'MCONLY: '+str(len(es))+' were found during the mcONLY simulator'
    
    if not silent:
        print 'MCONLY: $$$$$$$$$$ SIMULATOR complete $$$$$$$$$$\n'
        print 'Best orbit found:'
        print "chiSquaredMin = ",chiSquareMin
        print "LongAN = ",longAN_degs[bestOrbit]
        print "e = ",es[bestOrbit]
        print "To = ",Ts[bestOrbit]
        print "period = ",periods[bestOrbit]
        print "inclination = ",inclination_degs[bestOrbit]
        print "argPeri = ",argPeri_degs[bestOrbit]
        print "a1 = ",np.mean(a1s2[bestOrbit])
        print "a2 = ",np.mean(a2s2[bestOrbit])
    
    print 'MCONLY: Completed at '+time.strftime("%H:%M:%S", time.localtime())
    
    # record the time the chain finished and print
    endTime = time.clock()
    totalTime = (endTime-startTime) # in seconds
    if totalTime<60:
        print 'MCONLY: It took '+str(int(totalTime))+' seconds.\n'  ##### only print in 'silent' mode to track time
    elif totalTime>3600:
        totalTime = totalTime/3600.0
        print 'MCONLY: It took '+str(int(totalTime))+' hours and '+str(int(60*((totalTime)-int(totalTime))))+' minutes.\n'  ##### only print in 'silent' mode to track time
    else:
        totalTime = totalTime/60.0
        print 'MCONLY: It took '+str(int(totalTime))+' minutes and '+str(int(60*((totalTime)-int(totalTime))))+' seconds.\n'  ##### only print in 'silent' mode to track time

    if not silent:
        print '\nMCONLY: $$$$$$$$$$$$$ FINISHED mcONLYUniform1 Sim $$$$$$$$$$$$$$$\n'
  
    return (longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds)
    #****** DONE!! ********        
            
            
            
            
            
            
            
            
            
            
            
            
    
    