# @Author: Kyle Mede, kylemede@gmail.com, written while at the University of Tokyo 2011-2013.
"""
This is the latest version of the MCMC simulator modified from mcmcOrbSimulatorUniform5.py
with the changes understood after meeting with Tim Brant in January 2012.

THIS VERSION HAS THE FIXED SIGMA VALUE AND VARIES A SINGLE PARAMETER AT A TIME, INSTEAD
 OF VARYING ALL THE PARAMETERS AT TIME WHICH HAD POOR PERFORMANCE.  The parameter 
 being perturbed will change each trial, rather than every 100 or so like before.
"""

import math
import numpy as np
import os
import random as rand
import time
from orbitToolbox import multiEpochOrbCalc
    
#*********** In this brick is all you need to play with *************    
#********************************************************************
##### Overall simulation parameters ########
numSamples = int(1e7)
filename = '../data/MCMC7_1-Tau-Boo-ALL-1'
silent = False
verbose = True                                                               
numSamplePrints = 100.0     
temperature = 10000.0
startRandom = True
##### model input parameter perturbing variables #########
sigmaPercent = 1.1  # fixed sigma value
rateSamples = 1000 # used to calculate acceptance rate only
#### Ranges for acceptable random number inputs ######                                        
longAN_degMIN = 163.0 # [deg]                                                        
longAN_degMAX = 168.0 # [deg]                                                        
eMIN = 0.56                                                                      
eMAX = 0.74                                                                      
periodMIN = 365.0 # [yrs]                                                                
periodMAX = 530.0 # [yrs]                                                            
inclination_degMIN = 52.5 # [deg]                                                    
inclination_degMAX = 68.5 # [deg]                                                
argPeri_degMIN = 350.0 # [deg]                                                        
argPeri_degMAX = 377.0 # [deg]          
a_totalMIN = 86.0
a_totalMAX = 96.0                                        
#### MEASURED VALUES ## MUST BE IN ORDER FROM OLDEST TO NEWEST!! ########################                                        
#### TEST system data (HD 130948BC Dupuy2009)
#SA_arcsec_measured_REALs = [0.0993, 0.0946, 0.0961, 0.0978, 0.0568, 0.0564, 0.0573, 0.1117, 0.1091, 0.0979, 0.0722, 0.0574, 0.0573, 0.0584, 0.0517]
#SA_mean_errors = [0.0017, 0.0017, 0.0023, 0.0023, 0.0009, 0.0012, 0.0009, 0.0015, 0.0011, 0.0009, 0.0006, 0.0005, 0.0009, 0.0010, 0.0006] # mean simply implies the mean of the + and - uncertainties                                         
#PA_deg_measured_REALs = [307.9,306.9,308.9,307.0,144.9,146.4,148.6,132.6,132.28,130.73,127.6,124.7,124.6,124.1,123.9]
#PA_mean_errors = [1.1,1.0,1.6,1.4,0.5,0.6,0.6,0.4,0.13,0.17,0.3,0.4,0.6,0.7,0.5]# mean simply implies the mean of the + and - uncertainties 
#epochs = [2452519.702,2452519.706,2452519.709,2452519.719,2453425.152,2453425.219,2453425.282,2454127.085,2454185.028,2454306.737,2454481.192,2454555.068,2454555.09,2454555.157,2454584.106]
#Sys_Dist_PC = 18.17
#Mass1 = 1
#Mass2 = 1
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
SA_arcsec_measured_REALs = [10.26, 9.03, 8.91, 8.71, 8.9, 8.89, 8.72, 8.41, 8.45, 8.54, 8.66, 7.28, 7.35, 6.96, 7.06, 6.4, 6.15, 6.09, 5.87, 5.4, 5.32, 5.14, 4.94, 5.17,      3.401, 2.865, 2.82, 2.71, 2.181, 1.9421, 1.9324]
SA_mean_errors = [1.46, 0.12, 0.04, 0.2, 0.08, 0.1, 0.05, 0.22, 0.13, 0.04, 0.22, 0.89, 0.37, 0.66, 0.37, 0.37, 0.45, 0.42, 0.25, 0.33, 0.22, 0.07, 0.03, 0.38,         0.06, 0.04, 0.05, 0.03, 0.02, 0.01, 0.01]                                         
PA_deg_measured_REALs = [347.8, 348.7, 350.8, 351.9, 354.3, 352.4, 353.6, 356, 355.1, 354.8, 355.8, 356.8, 359.1, 359.5, 1.7, 3.5, 5.1, 5.1, 6, 6.8, 7.8, 11.4, 8, 13.5,               20.65, 30.85, 33.2, 31.3,46.26, 52.65, 53.68]
PA_mean_errors = [0.34, 1.88, 0.63, 0.17, 1.79, 0.38, 0.31, 1.43, 0.19, 0.16, 0.45, 0.21, 0.1, 0.09, 1.51, 0.93, 1.34, 1.04, 0.28, 0.75, 0.68, 1.2, 3.69, 1.08,           0.5, 0.4, 1.0, 1.0, 0.42, 0.16, 0.07]
epochs = [2396557.5, 2403065.5, 2405359.5, 2407064.5, 2408219.5, 24708927.5, 2411846.5, 2413509.5, 2414351.5, 2414465.5, 2415404.5, 2419263.5, 2423443.5, 2424260.5, 2425728.5, 2429753.5, 2431535.5, 2431959.5, 2434136.5, 2436256.5, 2437227.5, 2438866.5, 2440132.5, 2440707.5,      2448337.5, 2451143.5, 2451711.5, 2451945.5, 2454579.5, 2455589.5, 2455704.5]
Sys_Dist_PC = 15.62
Mass1 = 1 
Mass2 = 1 
#********************************************************************
#********************************************************************

if not silent:
    print '\nMCMC: $$$$$$$$$$$$$$$$$$$  STARTING THE SIMULATOR  $$$$$$$$$$$$$$$$$'
    numProcesses = 1
    numSamplesString = str(numSamples/int(1e6))+'-million'
    print 'Number of sample orbits being created = '+numSamplesString

numEpochs = len(epochs)

# set the max/min values for Time of last periapsis based on first epoch of observation
TMAX = epochs[-1]
TMIN = TMAX-periodMAX*365.25

pars = []
parsProposed = []
# init input parameter values
if startRandom:
    pars.append(rand.uniform(inclination_degMIN, inclination_degMAX))
    pars.append(rand.uniform(eMIN, eMAX))
    pars.append(rand.uniform(longAN_degMIN, longAN_degMAX))
    pars.append(rand.uniform(periodMIN, periodMAX))
    pars.append(rand.uniform(argPeri_degMIN, argPeri_degMAX))
    pars.append(rand.uniform(a_totalMIN, a_totalMAX))
    pars.append(rand.uniform(TMIN, TMAX))
else:
    pars.append(95.9) 
    pars.append(0.145)
    pars.append(133.37)
    pars.append(9.74)
    pars.append(73.85)
    pars.append(2.18)
    pars.append(2451124.1)

reduced_chiSquared_latest = 200.0
chiSquare_latest = (2.0*numEpochs-7.0)*reduced_chiSquared_latest
timesBeenHere = 1
paramBeingVaried = 0

# set max and min for sin(inclination_deg)
if inclination_degMIN>90.0:
    sine_inc_MAX = math.sin(math.radians(inclination_degMIN))        
    sine_inc_MIN = math.sin(math.radians(inclination_degMAX))  
else:
    sine_inc_MIN = math.sin(math.radians(inclination_degMIN))        
    sine_inc_MAX = math.sin(math.radians(inclination_degMAX))

# Load up parMINs and parMAXs lists
parMINs = [inclination_degMIN,eMIN,longAN_degMIN,periodMIN,argPeri_degMIN,a_totalMIN,TMIN]
parMAXs = [inclination_degMAX,eMAX,longAN_degMAX,periodMAX,argPeri_degMAX,a_totalMAX,TMAX]

# sigma is a percentage of total range, so convert its current percentage to a range value
sine_inc_sigma = (sigmaPercent/100.0)*(sine_inc_MAX-sine_inc_MIN)
e_sigma = (sigmaPercent/100.0)*(eMAX-eMIN)
longAN_sigma = (sigmaPercent/100.0)*(longAN_degMAX-longAN_degMIN)
period_sigma = (sigmaPercent/100.0)*(periodMAX-periodMIN)
argPeri_sigma = (sigmaPercent/100.0)*(argPeri_degMAX-argPeri_degMIN)
a_total_sigma = (sigmaPercent/100.0)*(a_totalMAX-a_totalMIN)
T_sigma = (sigmaPercent/100.0)*(TMAX-TMIN)

# Load up sigmas list
sigmas = [sine_inc_sigma,e_sigma,longAN_sigma,period_sigma,argPeri_sigma,a_total_sigma,T_sigma]

# variables for the success rate print block in chain loop
printTime = numSamples/numSamplePrints
printCount = 0
printsDone = 0
totalAcceptedCounter = 0
rateAcceptedCounter = 0
reduced_chiSquareMin = 1e100
bestOrbit = 0
acceptRates = []
acceptRates.append(0)
rateTotalCounter = 0

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
    # Set proposed pars to current pars and update 
    parsProposed = pars
    
    if (parsProposed[1]<0.0) or (parsProposed[1]>1.0):#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        #print 'e proposed out of range at top of sample loop, value = ',parsProposed[1]#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        l=1##$$$$$$$$$$$$$$$
    # block to control printing success rate to screen
    printCount = printCount + 1
    if printCount==printTime:
        printsDone = printsDone+1
        printCount = 0
        if not silent:
            timeStr = time.strftime("%H:%M:%S", time.localtime())
            reduced_chiSquareMinStr = ". Lowest reduced_chiSquare so far = "+str(reduced_chiSquareMin)
            print "\n"+str(totalAcceptedCounter)+' successful out of '+str(curSample)+'. '+str(printsDone)+\
                                '/'+str(int(numSamplePrints))+' completed at '+timeStr+reduced_chiSquareMinStr
            print "Latest acceptance rate = "+str(acceptRates[-1])+" and sigmaPercent = "+str(sigmaPercent)
            print "Latest likelihood_ratio = "+str(likelihood_ratio)+" and alpha = "+str(alpha)
            print "Latest chiSquare_latest - chi_squared_total_cur = "+str(chiSquare_latest - chi_squared_total_cur)
            print 'latest reduced_chi-square = ',reduced_chi_squared_total_cur
            print "inclination_deg_proposed = ",parsProposed[0]


    ## propose a new set of input parameters by perturbing the current ones 
    if paramBeingVaried==0:
        # inclination
        # NOTE: using uniform in sin(i) rather than i as it then avoids even sampling
        #       of near face-on orbits that are rare or unlikely. 
        try:
            increment = rand.uniform(-1.0,1.0)*sigmas[0] # sigma is sine_inc based so increment is sine_inc based
            parsProposed[0] = math.degrees(math.asin(math.sin(math.radians(pars[0]))+increment)) # par is inc based
            initialProposedInc = parsProposed[0]
            #try:
            # handle case where degrees form of asin(sin_inc) is on wrong side of 90 degrees
            if parMINs[0]>90.0:
                if parsProposed[0]<90.0:
                    parsProposed[0] = 180.0-parsProposed[0]
            elif parMAXs[0]<90.0:
                if parsProposed[0]>90.0:
                    parsProposed[0] = 180.0-parsProposed[0]
            
        except:
            print '\nmath.sin(math.radians(pars[0]))+increment = ',math.sin(math.radians(pars[0]))+increment
            print 'increment = ',increment
            print 'pars[0] = ',pars[0]
            print 'math.asin(math.sin(math.radians(pars[0]))+increment) = ',math.asin(math.sin(math.radians(pars[0]))+increment)
#            print "sine_inc_proposed = ",math.sin(math.radians(pars[0]))+increment
#            print "sine_inc_MIN = ",sine_inc_MIN
#            print "sine_inc_MAX = ",sine_inc_MAX
#            print "inclination proposed = ",parsProposed[0]
#        if ((parsProposed[0]>parMAXs[0]) or (parsProposed[0]<parMINs[0])):
#            print "\nsine_inc_proposed = ",math.sin(math.radians(pars[0]))+increment
#            print "sine_inc_MIN = ",sine_inc_MIN
#            print "sine_inc_MAX = ",sine_inc_MAX
#            print "math.asin(sine_inc_proposed) = ",math.asin(math.sin(math.radians(pars[0]))+increment)
#            print "initially proposed inc degs = ",initialProposedInc
#            print "inclination_deg_proposed = ",parsProposed[0]
#            print "inclination_degMIN = ",inclination_degMIN
#            print "inclination_degMAX = ",inclination_degMAX
    
    elif paramBeingVaried==1:
        # eccentricity
        # NOTE: using f(e)=2e as suggested by Dupeunnoy for orbits with periods over 1000days.
        increment = rand.uniform(-1.0,1.0)*sigmas[1]#*2.0 #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        parsProposed[1] = pars[1]+increment
        
        caughtByPrevCheck = False
        if parsProposed[1]<0.0:
            caughtByPrevCheck = True
            #print 'e proposed below zero = ',parsProposed[1]
            #print 'pars[1] = ',pars[1]
            #print 'eMIN = ',parMINs[1]
            #print 'increment = ',increment
        if parsProposed[1]>1.0:
            caughtByPrevCheck = True
            #print 'e proposed above 1 = ',parsProposed[1]
            #print 'pars[1] = ',pars[1]
            #print 'increment = ',increment
    else:
        # for all other parameters.  ie. longAN, period, argPeri, a and T.
        increment = rand.uniform(-1.0,1.0)*sigmas[paramBeingVaried]
        parsProposed[paramBeingVaried] = pars[paramBeingVaried]+increment

    ### Check if all pars were found to be in their allowed ranges
    ALLparsPassed = True
    for par in range(0,len(pars)):
        if (parsProposed[par]>parMAXs[par]) or (parsProposed[par]<parMINs[par]):
            ALLparsPassed = False
    #if ALLparsPassed==False:
    #    print 'par out of range, so skipping orbCalc'
    # ALL proposed pars passed, so continue to calculating orbit
    if ALLparsPassed:
        if (parsProposed[1]<0.0) or (parsProposed[1]>1.0):
            print 'made it through ALLparsPassed but e still out of range'
            print 'status of caughtByPrevCheck is ',str(caughtByPrevCheck)
        ## send random parameters along with known ones to multi-epoch orbit calculator
        (chi_squared_total_cur, ns, Ms, Es, thetas, Sep_Dists, SA_deg_measured_models, PA_deg_measured_models, a1s, a2s) =\
        multiEpochOrbCalc(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors,\
                       epochs, Sys_Dist_PC, parsProposed[0], parsProposed[2], parsProposed[1], parsProposed[6], \
                       parsProposed[3], parsProposed[4], parsProposed[5], Mass1=Mass1, Mass2=Mass2, verbose=False)
    
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
            
            rateAcceptedCounter += 1 # for calculating the acceptance rate
            totalAcceptedCounter += 1 # for counting total accepted overall
            
            # accepted, so proposed becomes new set of parameters to perturb 
            pars = parsProposed
            
            if (parsProposed[1]<0.0) or (parsProposed[1]>1.0):#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                print 'e proposed out of range inside accepted orbit brick, value = ',parsProposed[1]#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            
            if reduced_chi_squared_total_cur<reduced_chiSquareMin:
                reduced_chiSquareMin = reduced_chi_squared_total_cur
                longAN_deg_best = parsProposed[2]
                e_best = parsProposed[1]
                T_best = parsProposed[6]
                period_best = parsProposed[3]
                inclination_deg_best = parsProposed[0]
                argPeri_deg_best = parsProposed[4]
                a_total_best = parsProposed[5]
            
            # store output orbital elements of model
            a1 = np.mean(a1s[:])  # saving mean of each epoch to save room from start compared to old versions
            a2 = np.mean(a2s[:])  # saving mean of each epoch to save room from start compared to old versions
            
            ##************** Writing new line of data to file ******************************************
            # create string line of data to write to file
            line = str(parsProposed[2])
            line = line +'   '+ str(parsProposed[1])
            line = line +'   '+ str(parsProposed[6])
            line = line +'   '+ str(parsProposed[3])
            line = line +'    '+ str(parsProposed[0])
            line = line +'      '+ str(parsProposed[4])
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
            ##************** new line finished being written to file ***********************************
            
            # reset the timesBeenHere to 1
            timesBeenHere = 1
         
        # did not pass the Metropolis-Hastings Algorithm acceptance 
        # so we are staying on this model/set of pars for another round so 
        # increment the timesBeenHere counter
        else:
            timesBeenHere +=1
        
        # increment the paramBeingVaried value
        if paramBeingVaried<6:
            paramBeingVaried += 1
        else:
            paramBeingVaried = 0
                
        # increment rate counter
        rateTotalCounter += 1
        
        if rateTotalCounter==rateSamples:
            # it is time to recalculate the acceptance rate 
            acceptRates.append(float(rateAcceptedCounter)/float(rateSamples))
            rateAcceptedCounter = 0
            rateTotalCounter = 0
           
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
    print "mean acceptance rate = ",np.mean(acceptRates)

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
    print '\nMCMC: $$$$$$$$$$$$$ FINISHED MCMC-Uniform7 Sim $$$$$$$$$$$$$$$\n'

#****** DONE!! ********        
            
            

    
    