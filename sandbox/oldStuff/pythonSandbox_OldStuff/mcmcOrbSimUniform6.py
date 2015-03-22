# @Author: Kyle Mede, kylemede@gmail.com, written while at the University of Tokyo 2011.
"""
This is the latest version of the MCMC simulator modified from mcmcOrbSimulatorUniform5.py
with the changes understood after meeting with Tim Brant in January 2012.

    THIS VERSION HAS THE SIGMA TUNING TO FIND THE BEST VALUE THAT CAUSES THE OPTIMAL 
    ACCEPTANCE RATE.
"""

import math
import numpy as np
import os
import random as rand
import time
from orbitToolbox import multiEpochOrbCalc

#def mcmcUniformOrbSim(numSamples, filename, silent=True):
    
#*********** In this brick is all you need to play with *************    
#********************************************************************
##### Overall simulation parameters ########
numSamples = int(1e6)
filename = '../data/NewSimTest'
silent = False
chiSquaredMax = 1000.0
verbose = True                                                               
numSamplePrints = 100.0     
temperature = 50000.0
##### sigma tuning parameters #########
acceptRate_MAX = 0.35
acceptRate_MIN = 0.2
sigmaPercent_MAX = 5.0 # ie. normal percentage units with max of 100%
sigmaPercent_MIN = 0.25
sigmaTuningIncrement = 0.25
tuningSamples = 1000
#### Ranges for acceptable random number inputs ######                                        
longAN_degMIN = 1.0 # [deg]                                                        
longAN_degMAX = 359.0 # [deg]                                                        
eMIN = 0.01                                                                     
eMAX = 0.05                                                                       
periodMIN = 1.0 # [yrs]                                                                
periodMAX = 20.0 # [yrs]                                                            
inclination_degMIN = 90.01 # [deg]                                                    
inclination_degMAX = 179.0 # [deg]                                   
argPeri_degMIN = 1.0 # [deg]                                                        
argPeri_degMAX = 179.0 #179.999#89.999 # [deg]          
a_totalMIN = 1.0
a_totalMAX = 5.0                                         
#### MEASURED VALUES ########################                                        

#### TEST system data (HD 130948BC Dupuy2009)
SA_arcsec_measured_REALs = [0.0993, 0.0946, 0.0961, 0.0978, 0.0568, 0.0564, 0.0573, 0.1117, 0.1091, 0.0979, 0.0722, 0.0574, 0.0573, 0.0584, 0.0517]
SA_mean_errors = [0.0017, 0.0017, 0.0023, 0.0023, 0.0009, 0.0012, 0.0009, 0.0015, 0.0011, 0.0009, 0.0006, 0.0005, 0.0009, 0.0010, 0.0006] # mean simply implies the mean of the + and - uncertainties                                         
PA_deg_measured_REALs = [307.9,306.9,308.9,307.0,144.9,146.4,148.6,132.6,132.28,130.73,127.6,124.7,124.6,124.1,123.9]
PA_mean_errors = [1.1,1.0,1.6,1.4,0.5,0.6,0.6,0.4,0.13,0.17,0.3,0.4,0.6,0.7,0.5]# mean simply implies the mean of the + and - uncertainties 
epochs = [2452519.702,2452519.706,2452519.709,2452519.719,2453425.152,2453425.219,2453425.282,2454127.085,2454185.028,2454306.737,2454481.192,2454555.068,2454555.09,2454555.157,2454584.106]
Sys_Dist_PC = 18.17
Mass1 = 1
Mass2 = 1
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
    print '\nMCMC: $$$$$$$$$$$$$$$$$$$  STARTING THE SIMULATOR  $$$$$$$$$$$$$$$$$'
    numProcesses = 1
    numSamplesString = str(numSamples/int(1e6))+'-million'
    print 'Number of sample orbits being created = '+numSamplesString

numEpochs = len(epochs)

# init input lists
longAN_deg_latest = 133.37
e_latest = 0.145
T_latest = 2451124.1
period_latest = 9.74
inclination_deg_latest = 95.9 
sine_inc_latest = math.sin(math.radians(inclination_deg_latest))
argPeri_deg_latest = 73.85
a_total_latest = 2.18
reduced_chiSquared_latest = 174.0
chiSquare_latest = (2.0*numEpochs-7.0)*reduced_chiSquared_latest
sigmaPercent_latest = int((sigmaPercent_MAX + sigmaPercent_MIN)/2.0) # ie. the middle of the range
timesBeenHere = 1

# set max and min for sin(inclination_deg)
if inclination_degMIN>90.0:
    sine_inc_MAX = math.sin(math.radians(inclination_degMIN))        
    sine_inc_MIN = math.sin(math.radians(inclination_degMAX))  
else:
    sine_inc_MIN = math.sin(math.radians(inclination_degMIN))        
    sine_inc_MAX = math.sin(math.radians(inclination_degMAX))

# variables for the success rate print block in chain loop
printTime = numSamples/numSamplePrints
printCount = 0
printsDone = 0
totalAcceptedCounter = 0
tuningAcceptedCounter = 0
reduced_chiSquareMin = chiSquaredMax
bestOrbit = 0
acceptRate = 0
tuningTotalCounter = 0

 # record the time the chain started
startTime = time.clock()

## Open and start the file to write the orbits/models too
if filename[-4:]!='.txt':
    filename = filename+".txt"
if not silent:
    print "Writing data to "+filename
f = open(filename, 'w')
fname = os.path.basename(filename)
f.write(fname+'\n')
if ((Mass1==1) and (Mass2==1)):
    # ie. both set to default and only saving input a_total
    f.write('longAN [deg]      e [N/A]       To [julian date]  period [yrs]   inclination [deg]   argPeri [deg]   a_total [AU]    chiSquared    timesBeenHere\n')
else:
    # ie. non-default, so store individual a1 and a2
    f.write('longAN [deg]      e [N/A]     To [julian date]    period [yrs]   inclination [deg]   argPeri [deg]   a1 [AU]    a2[AU]   chiSquared  timesBeenHere\n')

## Start to create samples and check if they are are accepted
## ie. Start the chain
for curSample in range(1,numSamples+1):
    
    # block to control printing success rate to screen
    printCount = printCount + 1
    if printCount==printTime:
        printsDone = printsDone+1
        printCount = 0
        if not silent:
            timeStr = time.strftime("%H:%M:%S", time.localtime())
            reduced_chiSquareMinStr = ". Lowest reduced_chiSquare so far = "+str(reduced_chiSquareMin)
            print "\n"+str(totalAcceptedCounter)+' successful out of '+str(curSample)+'. '+str(printsDone)+'/'+str(int(numSamplePrints))+' completed at '+timeStr+reduced_chiSquareMinStr
            print "Latest acceptance rate = "+str(acceptRate)+" and sigmaPercent = "+str(sigmaPercent_latest)
            print "Latest likelihood_ratio = "+str(likelihood_ratio)+" and alpha = "+str(alpha)
            print "Latest chiSquare_latest - chi_squared_total_cur = "+str(chiSquare_latest - chi_squared_total_cur)
            print 'latest reduced_chi-square = ',reduced_chi_squared_total_cur
            print "inclination_deg_proposed = ",inclination_deg_proposed


    ## propose a new set of input parameters by perturbing the current ones 
    
    # inclination
    # sigma is a percentage of total range, so convert its current percentage to a 
    # range value
    # NOTE: using uniform in sin(i) rather than i as it then avoids even sampling
    #       of near face-on orbits that are rare or unlikely. 
    sigmaRange = (sigmaPercent_latest/100.0)*(sine_inc_MAX-sine_inc_MIN)
    if (sine_inc_latest+sigmaRange)>=sine_inc_MAX:
        max = sine_inc_MAX
        min = sine_inc_latest-sigmaRange
    elif (sine_inc_latest-sigmaRange)<=sine_inc_MIN:
        max = sine_inc_latest+sigmaRange
        min = sine_inc_MIN
    else:
        max = sine_inc_latest + sigmaRange
        min = sine_inc_latest - sigmaRange
    sine_inc_proposed = rand.uniform(min, max)
    try:
        inclination_deg_proposed = math.degrees(math.asin(sine_inc_proposed))
        # handle case where degrees form of asin(sin_inc) is on wrong side of 90 degrees
        if inclination_degMIN>90.0:
            if inclination_deg_proposed<90.0:
                inclination_deg_proposed = 180.0-inclination_deg_proposed
        elif inclination_degMIN<90.0:
            if inclination_deg_proposed>90.0:
                inclination_deg_proposed = 180.0-inclination_deg_proposed
        
    except:
        print "sine_inc_proposed = ",sine_inc_proposed
        print "sine_inc_MIN = ",sine_inc_MIN
        print "min = ",min
        print "sine_inc_MAX = ",sine_inc_MAX
        print "max = ",max
        print "math.asin(sine_inc_proposed) = ",math.asin(sine_inc_proposed)
    if ((inclination_deg_proposed>inclination_degMAX) or (inclination_deg_proposed<inclination_degMIN)):
        print "\nsine_inc_proposed = ",sine_inc_proposed
        print "sine_inc_MIN = ",sine_inc_MIN
        print "min = ",min
        print "sine_inc_MAX = ",sine_inc_MAX
        print "max = ",max
        print "math.asin(sine_inc_proposed) = ",math.asin(sine_inc_proposed)
        print "inclination_deg_proposed = ",inclination_deg_proposed
        print "inclination_degMIN = ",inclination_degMIN
        print "inclination_degMAX = ",inclination_degMAX
        
    # long of acending node
    # sigma is a percentage of total range, so convert its current percentage to a 
    # range value
    sigmaRange = (sigmaPercent_latest/100.0)*(longAN_degMAX-longAN_degMIN)
    if (longAN_deg_latest+sigmaRange)>=longAN_degMAX:
        max = longAN_degMAX
        min = longAN_deg_latest-sigmaRange
    elif (longAN_deg_latest-sigmaRange)<=longAN_degMIN:
        max = longAN_deg_latest+sigmaRange
        min = longAN_degMIN
    else:
        max = longAN_deg_latest + sigmaRange
        min = longAN_deg_latest - sigmaRange
    longAN_deg_proposed = rand.uniform(min, max)
    
    #eccentricity
    # sigma is a percentage of total range, so convert its current percentage to a 
    # range value
    # NOTE: using f(e)=2e as suggested by Dupeunnoy for orbits with periods over 1000days.
    sigmaRange = (sigmaPercent_latest/100.0)*(eMAX-eMIN)
    if (e_latest+sigmaRange)>=eMAX:
        max = eMAX
        min = e_latest-sigmaRange
    elif (e_latest-sigmaRange)<=eMIN:
        max = e_latest+sigmaRange
        min = eMIN
    else:
        max = e_latest + sigmaRange
        min = e_latest - sigmaRange
    double_e_proposed = rand.uniform(2.0*min, 2.0*max)
    e_proposed = 0.5*double_e_proposed
    
    #period
    # sigma is a percentage of total range, so convert its current percentage to a 
    # range value
    sigmaRange = (sigmaPercent_latest/100.0)*(periodMAX-periodMIN)
    if (period_latest+sigmaRange)>=periodMAX:
        max = periodMAX
        min = period_latest-sigmaRange
    elif (period_latest-sigmaRange)<=periodMIN:
        max = period_latest+sigmaRange
        min = periodMIN
    else:
        max = period_latest + sigmaRange
        min = period_latest - sigmaRange
    period_proposed = rand.uniform(min, max)
    
    # time of last periapsis
    # sigma is a percentage of total range, so convert its current percentage to a 
    # range value
    # set the max/min values for Time of last periapsis based on first epoch of observation
    TMAX = epochs[-1]-period_proposed*365.0
    TMIN = epochs[-1]
    sigmaRange = (sigmaPercent_latest/100.0)*(TMAX-TMIN)
    if (T_latest+sigmaRange)>=TMAX:
        max = TMAX
        min = T_latest-sigmaRange
    elif (T_latest-sigmaRange)<=TMIN:
        max = T_latest+sigmaRange
        min = TMIN
    else:
        max = T_latest + sigmaRange
        min = T_latest - sigmaRange
    T_proposed = rand.uniform(min, max)

    # argument of perigie
    # sigma is a percentage of total range, so convert its current percentage to a 
    # range value
    sigmaRange = (sigmaPercent_latest/100.0)*(argPeri_degMAX-argPeri_degMIN)
    if (argPeri_deg_latest+sigmaRange)>=argPeri_degMAX:
        max = argPeri_degMAX
        min = argPeri_deg_latest-sigmaRange
    elif (argPeri_deg_latest-sigmaRange)<=argPeri_degMIN:
        max = argPeri_deg_latest+sigmaRange
        min = argPeri_degMIN
    else:
        max = argPeri_deg_latest + sigmaRange
        min = argPeri_deg_latest - sigmaRange
    argPeri_deg_proposed = rand.uniform(min, max)
    
    # Total semi-major axis
    # sigma is a percentage of total range, so convert its current percentage to a 
    # range value
    sigmaRange = (sigmaPercent_latest/100.0)*(a_totalMAX-a_totalMIN)
    if (a_total_latest+sigmaRange)>=a_totalMAX:
        max = a_totalMAX
        min = a_total_latest-sigmaRange
    elif (a_total_latest-sigmaRange)<=a_totalMIN:
        max = a_total_latest+sigmaRange
        min = a_totalMIN
    else:
        max = a_total_latest + sigmaRange
        min = a_total_latest - sigmaRange
    a_total_proposed = rand.uniform(min, max)
    
    ## send random parameters along with known ones to multi-epoch orbit calculator
    (chi_squared_total_cur, ns, Ms, Es, thetas, Sep_Dists, SA_deg_measured_models, PA_deg_measured_models, a1s, a2s) =\
    multiEpochOrbCalc(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors,\
                   epochs, Sys_Dist_PC, inclination_deg_proposed, longAN_deg_proposed, e_proposed, T_proposed, \
                   period_proposed, argPeri_deg_proposed, a_total_proposed, Mass1=Mass1, Mass2=Mass2, verbose=False)

    # calculate the likelihood ratio for the previous(latest) and current (proposed) models/states
    #likelihood_ratio = math.exp((chiSquare_latest - chi_squared_total_cur)/ 2.0)
    likelihood_ratio = math.exp((chiSquare_latest - chi_squared_total_cur)/ (2.0*temperature))
    
    
    # normalize the chiSquared to give the reduce chiSquared value
    reduced_chi_squared_total_cur = (1.0/((2.0*(numEpochs))-7.0))*chi_squared_total_cur
    
    alpha = rand.uniform(0.0, 1.0)
    
    if likelihood_ratio>1.0:
        RHS = 1.0
    else:
        RHS = likelihood_ratio

    if alpha < RHS:
        
        tuningAcceptedCounter += 1
        totalAcceptedCounter += 1
        
        longAN_deg_latest = longAN_deg_proposed
        e_latest = e_proposed
        T_latest = T_proposed
        period_latest = period_proposed
        inclination_deg_latest = inclination_deg_proposed
        argPeri_deg_latest = argPeri_deg_proposed
        a_total_latest = a_total_proposed
        sine_inc_latest = sine_inc_proposed
        chiSquare_latest = chi_squared_total_cur
        
        if reduced_chi_squared_total_cur<reduced_chiSquareMin:
            reduced_chiSquareMin = reduced_chi_squared_total_cur
            longAN_deg_best = longAN_deg_proposed
            e_best = e_proposed
            T_best = T_proposed
            period_best = period_proposed
            inclination_deg_best = inclination_deg_proposed
            argPeri_deg_best = argPeri_deg_proposed
            a_total_best = a_total_proposed
        
        # store output orbital elements of model
        a1 = np.mean(a1s[:])  # saving mean of each epoch to save room from start compared to old versions
        a2 = np.mean(a2s[:])  # saving mean of each epoch to save room from start compared to old versions
        
        # create string line of data to write to file
        line = str(longAN_deg_proposed)
        line = line +'   '+ str(e_proposed)
        line = line +'   '+ str(T_proposed)
        line = line +'   '+ str(period_proposed)
        line = line +'    '+ str(inclination_deg_proposed)
        line = line +'      '+ str(argPeri_deg_proposed)
        # take care of the case where no a1 is calculated (ie is zero)
        # as it throws off the number of spaces and column headers to make output file look pretty
        if a1 != 0.0:
            line = line +'   '+ str(a1)
            line = line +'   '+ str(a2)
        else:
            line = line +'   '+ str(a2) # a_total=a2 if Mass1 or Mass2 are non-default (ie !=1)
        line = line +'   '+ str(chi_squared_total_cur)
        line = line +'      '+ str(timesBeenHere)+'\n'
        # write final string to file
        f.write(line)
        
        # reset the timesBeenHere to 1
        timesBeenHere = 1
     
    # did not pass the Metropolis-Hastings Algorithm acceptance 
    # so we are staying on this model for another round so 
    # increment the timesBeenHere counter
    else:
        timesBeenHere +=1
    
    # increment tuning counter
    tuningTotalCounter += 1
    
    if tuningTotalCounter==tuningSamples:
        # it is time to recalculate the acceptance rate and adjust current sigma accordingly 
        acceptRate = float(tuningAcceptedCounter)/float(tuningSamples)
        tuningAcceptedCounter = 0
        tuningTotalCounter = 0
        
        if (sigmaPercent_latest+sigmaTuningIncrement)<=sigmaPercent_MAX:
            if acceptRate>=acceptRate_MAX:
                # ie. acceptance rate is too high, so raise current sigma being adjusted
                sigmaPercent_latest = sigmaPercent_latest + sigmaTuningIncrement
            elif acceptRate<acceptRate_MIN:
                # ie. acceptance rate is too low, so lower current sigma being adjusted
                if (sigmaPercent_latest-sigmaTuningIncrement)>sigmaPercent_MIN:
                    sigmaPercent_latest = sigmaPercent_latest - sigmaTuningIncrement
                #else:   
                    # all ready below close t 0, so leave it
        #else:
            # all ready at max so do nothing; ie leave sigmaPercent_lates as is
        
        #if not silent:
        #    print 'Finished a round of sigma tuning' #$$$$$$$$
           
# Done looping through all the samples

# Finished all samples, so close the file data is being written to
f.close()

print 'MCMC: '+str(totalAcceptedCounter)+' steps were taken during the MCMC simulator'

if not silent:
    print 'MCMC: $$$$$$$$$$ SIMULATOR complete $$$$$$$$$$\n'
    print 'Best orbit found:'
    print "chiSquaredMin = ",reduced_chiSquareMin
    print "LongAN = ",longAN_deg_best
    print "e = ",e_best
    print "To = ",T_best
    print "period = ",period_best
    print "inclination = ",inclination_deg_best
    print "argPeri = ",argPeri_deg_best
    print "a_total = ",a_total_best

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
    print '\nMCMC: $$$$$$$$$$$$$ FINISHED MCMC-Uniform6 Sim $$$$$$$$$$$$$$$\n'

#****** DONE!! ********        
            
            

    
    