# @Author: Kyle Mede, kylemede@gmail.com, written while at the University of Tokyo 2011-2013.
"""
This is the latest version of the MCMC simulator modified from mcmcOrbSimulatorUniform5.py
with the changes understood after meeting with Tim Brant in January 2012.

THIS VERSION HAS THE FIXED SIGMA VALUE AND VARIES A SINGLE PARAMETER AT A TIME, INSTEAD
 OF VARYING ALL THE PARAMETERS AT TIME WHICH HAD POOR PERFORMANCE.  The parameter 
 being perturbed will change each trial, rather than every 100 or so like before.
 TRY 2 at making this work.
"""

import math
import numpy as np
import os
import random as rand
import shutil
import time
from orbitToolbox import multiEpochOrbCalc, timeString, dictToFile, chiSquaredCalc, rv2bodyCalculator3, rv2bodyCalculator4,\
                               rv1bodyCalculator, rv1bodyCalculator2, semiMajorConverter
from paramSettingsDict import paramSettingsDict
    
def mcOrbSimUniform(filename, numSamples, silent=True, verbose=False,\
                        dataTempDir='/home/Kyle/DataTempForSims/', dataFinalDir='/media/Data1/Todai_Work/workspace/Binary-project/data/'):
    #*********** In this brick is all you need to play with *************    
    #********************************************************************
    ##### Overall simulation parameters ########
    #dataTempDir = '/media/Data1/Todai_Work/workspace/Binary-project/data/'#'/home/Kyle/DataTempForSims/' 
    #dataFinalDir = '/media/Data1/Todai_Work/workspace/Binary-project/data/'
    #filename = 'SimAneal7_2-Tau-Boo-ALL-T10-TEMP'  #This is a feature implemented as I have an SSD which is much 
    #                                            #MUCH faster at reading/writing for use during the sim
    #                                            # and the finished file is then moved to the storage dir on 
    #                                            # my standard HDD at the end of the sim.
    #                                            # making the Temp and Final the same negates this if you don't have an SSD.
                                                           

    ##### model input parameter variables #########
    numSamplePrints = paramSettingsDict['General']['numSamplePrints']
    chi_squared_MAX = paramSettingsDict['General']['chi_squared_MAX']
    gaussMassAndDist = paramSettingsDict['General']['gaussMassAndDist'] # calculate the mass and sys_dist values using gaussian dists?
    #NOTE: leave both of these 'False' to run standard 3D simulations
    DIonly = paramSettingsDict['General']['DIonly'] #setting this makes RV chiSquared=0
    RVonly = paramSettingsDict['General']['RVonly'] #setting this makes DI chiSquared=0
    # check both are not set to True
    if (DIonly==True) and (RVonly==True):
        print 'mcOrbSimUni47: WARNING: can not have both RVonly and DIonly set to True!'
    #### Ranges for acceptable random number inputs ######                                        
    longAN_degMIN = paramSettingsDict['MINS']['longAN']#130.0 # [deg]                                                        
    longAN_degMAX = paramSettingsDict['MAXS']['longAN']#179.9 # [deg]                                                        
    eMIN = paramSettingsDict['MINS']['e']#0.35                                                                      
    eMAX = paramSettingsDict['MAXS']['e']#0.95                                                                      
    periodMIN = paramSettingsDict['MINS']['P']#280.0 # [yrs]                                                                
    periodMAX = paramSettingsDict['MAXS']['P']#900.0 # [yrs]                                                            
    inclination_degMIN = paramSettingsDict['MINS']['i']#25.0 # [deg]                                                    
    inclination_degMAX = paramSettingsDict['MAXS']['i']#75.0 # [deg]                                                
    argPeri_degMIN = paramSettingsDict['MINS']['argPeri']#300.0 # [deg]                                                        
    argPeri_degMAX = paramSettingsDict['MAXS']['argPeri']#420.0 # [deg]          
    a_totalMIN = paramSettingsDict['MINS']['a']#80.0
    a_totalMAX = paramSettingsDict['MAXS']['a']#125.0                        
    #### MEASURED VALUES ########################                                        
    # DI data
    SA_arcsec_measured_REALs = paramSettingsDict['Data']['SAs']#[9.03, 8.91, 8.71, 8.9, 8.89, 8.72, 8.41, 8.45, 8.54, 8.66, 7.28, 7.35, 6.96, 7.06, 6.4, 6.15, 6.09, 5.87, 5.4, 5.32, 5.14, 4.94, 5.17,      3.401, 2.865, 2.82, 2.71, 2.181, 1.9421, 1.9324]
    SA_mean_errors = paramSettingsDict['Data']['SAerrors']#[0.12, 0.04, 0.2, 0.08, 0.1, 0.05, 0.22, 0.13, 0.04, 0.22, 0.89, 0.37, 0.66, 0.37, 0.37, 0.45, 0.42, 0.25, 0.33, 0.22, 0.07, 0.03, 0.38,         0.06, 0.04, 0.05, 0.03, 0.02, 0.01, 0.01]                                         
    PA_deg_measured_REALs = paramSettingsDict['Data']['PAs']#[348.7, 350.8, 351.9, 354.3, 352.4, 353.6, 356, 355.1, 354.8, 355.8, 356.8, 359.1, 359.5, 1.7, 3.5, 5.1, 5.1, 6, 6.8, 7.8, 11.4, 8, 13.5,               20.65, 30.85, 33.2, 31.3,46.26, 52.65, 53.68]
    PA_mean_errors = paramSettingsDict['Data']['PAerrors']#[1.88, 0.63, 0.17, 1.79, 0.38, 0.31, 1.43, 0.19, 0.16, 0.45, 0.21, 0.1, 0.09, 1.51, 0.93, 1.34, 1.04, 0.28, 0.75, 0.68, 1.2, 3.69, 1.08,           0.5, 0.4, 1.0, 1.0, 0.42, 0.16, 0.07]
    DI_epochs = paramSettingsDict['Data']['DI_epochs']#[2403065.5, 2405359.5, 2407064.5, 2408219.5, 24708927.5, 2411846.5, 2413509.5, 2414351.5, 2414465.5, 2415404.5, 2419263.5, 2423443.5, 2424260.5, 2425728.5, 2429753.5, 2431535.5, 2431959.5, 2434136.5, 2436256.5, 2437227.5, 2438866.5, 2440132.5, 2440707.5,      2448337.5, 2451143.5, 2451711.5, 2451945.5, 2454579.5, 2455589.5, 2455704.5]
    # RV data
    RV_epochs = paramSettingsDict['Data']['RV_epochs']
    RVerrors = paramSettingsDict['Data']['RVerrors']
    RVs = paramSettingsDict['Data']['RVs']
    # General System Data
    Sys_Dist_PC = paramSettingsDict['Data']['sysDist']#15.62
    Sys_Dist_Error = paramSettingsDict['Data']['sysDistError']#
    Mass1 = paramSettingsDict['Data']['M1']#1 
    Mass1_Error = paramSettingsDict['Data']['M1Error']#
    Mass2 = paramSettingsDict['Data']['M2']#1
    Mass2_Error = paramSettingsDict['Data']['M2Error']#
    #********************************************************************
    #********************************************************************
    
    if not silent:
        print '\nMCONLY: $$$$$$$$$$$$$$$$$$$  STARTING THE SIMULATOR  $$$$$$$$$$$$$$$$$'
        timeStr = time.strftime("%H:%M:%S", time.localtime())
        print 'Starting time for this simulation was, '+timeStr
        numProcesses = 1
        if numSamples<int(1e6):
            numSamplesString = str(numSamples/int(1e3))+'-thousand'
        else:
            numSamplesString = str(numSamples/int(1e6))+'-million'
        print 'Number of sample orbits being created = '+numSamplesString
    
    numEpochs_DI = len(DI_epochs)
    
    ## set the max/min values for Time of last periapsis based on first epoch of observation
    # first find first epoch of observation between RV and DI data
    if DI_epochs[0]<RV_epochs[0][0]:
        firstEpoch = DI_epochs[0]
    else:
        firstEpoch = RV_epochs[0][0]
    TMAX = RV_epochs[1][0]#firstEpoch
    TMIN = TMAX-periodMAX*365.25
    
    DI_epochs[0] = RV_epochs[0][0]#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    # set up min and max for RV origin params
    RV_origin_vel_0_MIN = -250#np.min(RVs[0])
    RV_origin_vel_0_MAX = 250#np.max(RVs[0])
    RV_origin_vel_1_MIN = 233.5#np.min(RVs[1])
    RV_origin_vel_1_MAX = 237#np.max(RVs[1])
    
    # set max and min for sin(inclination_deg)
    if inclination_degMIN>90.0:
        sine_inc_MAX = math.sin(math.radians(inclination_degMIN))        
        sine_inc_MIN = math.sin(math.radians(inclination_degMAX))  
    else:
        sine_inc_MIN = math.sin(math.radians(inclination_degMIN))        
        sine_inc_MAX = math.sin(math.radians(inclination_degMAX))
    
    # variables for the success rate print block in chain loop
    # other params here are just initial values that will be replaced soon in the loop
    printTime = numSamples/numSamplePrints
    printCount = 0
    printsDone = 0
    totalAcceptedCounter = 0
    reduced_chiSquareMin = 10000000.0
    chi_squared_DI_reducedMIN = 10000000.0
    chi_squared_RV_reducedMIN = 10000000.0
    bestOrbit = 0
    resultsDict = {}
    chi_squared_total_reduced = '?'
    inclination_deg_proposed ='?'
    
     # record the time the chain started
    startTime = time.clock()
    
    # start prints line for later writing to results dict
    printLines = ''
    
    ## Open and start the file to write the orbits/models too
    
    if filename[-4:]!='.txt':
        filename = filename+".txt"
    # Make temp filename in data temp dir for fast writing during sim
    filenameTEMP = dataTempDir+filename
    filenameFINAL = dataFinalDir+filename
    if not silent:
        print "Writing data temporarily to "+filenameTEMP
    f = open(filenameTEMP, 'w')
    fname = os.path.basename(filename)
    f.write(fname+'\n')
    #if ((Mass1==1) and (Mass2==1)):
    # ie. both set to default and only saving input a_total
    f.write('longAN [deg]      e [N/A]       To [julian date]  period [yrs]   inclination [deg]   argPeri [deg]   a_total [AU]  chiSquared   RVorigin_0 [m/s]   RVorigin_1 [m/s]    timesBeenHere\n')
    #else:
    #    # ie. non-default, so store individual a1 and a2
    #    f.write('longAN [deg]      e [N/A]     To [julian date]    period [yrs]   inclination [deg]   argPeri [deg]   a1 [AU]    a2[AU]   chiSquared  timesBeenHere\n')
    
    ## Start to create samples and check if they are are accepted
    ## ie. Start the chain
    for curSample in range(1,numSamples+1):
        
        # block to control printing success rate to screen
        printCount = printCount + 1
        if printCount==printTime:
            printsDone = printsDone+1
            printCount = 0
            if not silent:
                curTime = time.clock()
                timeSoFar = (curTime-startTime) # in seconds
                timeSoFarString = timeString(timeSoFar)
                timeStr = time.strftime("%H:%M:%S", time.localtime())
                reduced_chiSquareMinStr = "Lowest reduced chiSquareds so far:  DI = " +str(chi_squared_DI_reducedMIN)+', RV = '+str(chi_squared_RV_reducedMIN)+' , total = '+str(reduced_chiSquareMin)+', MAX = '+str(chi_squared_MAX)
                line = "\n"+str(printsDone)+'/'+str(int(numSamplePrints))+' done at '+timeStr+', running time: '+timeSoFarString+', '+str(totalAcceptedCounter)+' / '+str(curSample)+' successful. '    
                line = line+'\n'+reduced_chiSquareMinStr
                line = line+'\n'+ 'latest reduced_chi-square = '+str(chi_squared_total_reduced)
                line = line+'\n'+ "inclination_deg_proposed = "+str(inclination_deg_proposed)
                print line
                printLines = printLines+line
                
        ## Draw the values that have error but are not Orbital Elements(ie. the sys_dist and obj masses)
        ## then calculate their chiSquareds
        if gaussMassAndDist:
            # Draw a value for the System Distance from a Gaussian centered 
            # on the HIP value and sigma=error in that value
            Sys_Dist_proposed = rand.gauss(Sys_Dist_PC, Sys_Dist_Error/2.0)
            gauss_chi_squared_sys_dist = chiSquaredCalc(Sys_Dist_PC, Sys_Dist_Error, Sys_Dist_proposed)
            # Draw a value for the primary star's mass from a Gaussian 
            # centered on the previously found value and sigma=error in that value
            Mass1_proposed = rand.gauss(Mass1, Mass1_Error/2.0)
            gauss_chi_squared_sys_dist = chiSquaredCalc(Mass1, Mass1_Error, Mass1_proposed)
            # Draw a value for the companion star's (or planet) mass from a Gaussian 
            # centered on the previously found value and sigma=error in that value
            Mass2_proposed = rand.gauss(Mass2, Mass2_Error/2.0)
            gauss_chi_squared_sys_dist = chiSquaredCalc(Mass1, Mass2_Error, Sys_Dist_proposed)
        
            gauss_chi_squared_total = gauss_chi_squared_sys_dist+gauss_chi_squared_sys_dist+gauss_chi_squared_sys_dist
            numGaussParams = 3
        else:
            Mass1_proposed = Mass1
            Mass2_proposed = Mass2
            Sys_Dist_proposed = Sys_Dist_PC
            gauss_chi_squared_total = 0.0
            numGaussParams = 0
            
        # inclination
        # NOTE: using uniform in sin(i) rather than i as it then avoids even sampling
        #       of near face-on orbits that are rare or unlikely. 
        inclination_deg_proposed = rand.uniform(inclination_degMIN, inclination_degMAX) #87.21#
#        sine_inc_proposed = rand.uniform(sine_inc_MIN, sine_inc_MAX)
#        try:
#            inclination_deg_proposed = math.degrees(math.asin(sine_inc_proposed))
#            # handle case where degrees form of asin(sin_inc) is on wrong side of 90 degrees
#            if inclination_degMIN>90.0:
#                if inclination_deg_proposed<90.0:
#                    inclination_deg_proposed = 180.0-inclination_deg_proposed
#            elif inclination_degMIN<90.0:
#                if inclination_deg_proposed>90.0:
#                    inclination_deg_proposed = 180.0-inclination_deg_proposed
#            
#        except:
#            print "sine_inc_proposed = ",sine_inc_proposed
#            print "sine_inc_MIN = ",sine_inc_MIN
#            print "min = ",min
#            print "sine_inc_MAX = ",sine_inc_MAX
#            print "max = ",max
#            print "math.asin(sine_inc_proposed) = ",math.asin(sine_inc_proposed)
#        if ((inclination_deg_proposed>inclination_degMAX) or (inclination_deg_proposed<inclination_degMIN)):
#            print "\nsine_inc_proposed = ",sine_inc_proposed
#            print "sine_inc_MIN = ",sine_inc_MIN
#            print "min = ",min
#            print "sine_inc_MAX = ",sine_inc_MAX
#            print "max = ",max
#            print "math.asin(sine_inc_proposed) = ",math.asin(sine_inc_proposed)
#            print "inclination_deg_proposed = ",inclination_deg_proposed
#            print "inclination_degMIN = ",inclination_degMIN
#            print "inclination_degMAX = ",inclination_degMAX
        
        # long of accending node
        longAN_deg_proposed = rand.uniform(longAN_degMIN, longAN_degMAX)#1#
        
        #eccentricity
        e_proposed = rand.uniform(eMIN,eMAX) #0.68#
      
        #period
        period_proposed = 21.2165/365.25#rand.uniform(periodMIN, periodMAX)#21.2165/365.25#
        

        # time of last periapsis
        T_proposed = 2454777.94#rand.uniform(TMIN, TMAX)#2454777.94#
    
        # argument of perigie
        argPeri_deg_proposed = rand.uniform(argPeri_degMIN, argPeri_degMAX) #121.6#
        
#        print '\n** Proposed params are:'#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#        print 'inclination_deg_proposed = ',inclination_deg_proposed#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#        print 'longAN_deg_proposed = ',longAN_deg_proposed#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#        print 'e_proposed = ',e_proposed#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#        print 'period_proposed = ',period_proposed#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#        print 'T_proposed = ',T_proposed#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#        print 'argPeri_deg_proposed = ',argPeri_deg_proposed#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
      
#        # Total semi-major axis
#        a_total_proposed = rand.uniform(a_totalMIN, a_totalMAX)
        
        ## send random parameters along with known ones to multi-epoch orbit calculator
        if RVonly:
            (chi_squared_total_cur, ns, Ms, Es, thetas, Sep_Dists, SA_deg_measured_models, PA_deg_measured_models, a1s, a2s) =\
            multiEpochOrbCalc([SA_arcsec_measured_REALs[0]], [SA_mean_errors[0]], [PA_deg_measured_REALs[0]], [PA_mean_errors[0]],\
                       [DI_epochs[0]], Sys_Dist_proposed, inclination_deg_proposed, longAN_deg_proposed, e_proposed, T_proposed, \
                       period_proposed, argPeri_deg_proposed, a_total=False, Mass1=Mass1_proposed, Mass2=Mass2_proposed, verbose=verbose)
        else:
            (chi_squared_total_cur, ns, Ms, Es, thetas, Sep_Dists, SA_deg_measured_models, PA_deg_measured_models, a1s, a2s) =\
            multiEpochOrbCalc(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors,\
                       DI_epochs, Sys_Dist_proposed, inclination_deg_proposed, longAN_deg_proposed, e_proposed, T_proposed, \
                       period_proposed, argPeri_deg_proposed, a_total=False, Mass1=Mass1_proposed, Mass2=Mass2_proposed, verbose=verbose)
            
        # calculate a total from a1 plus a2
        a_total_proposed = np.mean(a1s)+np.mean(a2s)
        
        chi_squared_DI = chi_squared_total_cur
        # convert to a reduced chiSquared
        chi_squared_DI_reduced = (1.0/((2.0*(numEpochs_DI))-6.0))*chi_squared_DI
        
        if RVonly:
            chi_squared_DI_reduced = 0.0 
            chi_squared_DI = 0.0
            numEpochs_DI = 0.0
        # if only saving orbits based on DI data, set all RV chiSquared related values to 0.0
        # Else, go ahead and calculate the proper values for the RV data.
        if DIonly:
            chi_squared_RV_reduced = 0.0
            chi_squared_RV = 0.0
            numEpochs_RV = 0.0
            RV_origin_vel_0_proposed = 0
            RV_origin_vel_1_proposed = 0
            totalNumModelParams = 6.0 
        else:
            totalNumModelParams = 8.0 
            # propose RV origin velocity for each dataset
            RV_origin_vel_0_proposed = rand.uniform(RV_origin_vel_0_MIN, RV_origin_vel_0_MAX) #$$$$$$$$$$$$$
            RV_origin_vel_1_proposed = rand.uniform(RV_origin_vel_1_MIN, RV_origin_vel_1_MAX)
            # update RV values in each dataset for input to rv2bodyCalculator below
            RVs_proposed = RVs
            for epoch in range(0,len(RVs_proposed[0])):
                RVs_proposed[0][epoch] = RVs[0][epoch]-RV_origin_vel_0_proposed
            for epoch in range(0,len(RVs_proposed[1])):
                RVs_proposed[1][epoch] = RVs[1][epoch]-RV_origin_vel_1_proposed
        
            # data for planet's orbit that are needed to calculate proper chiSquare for RV fit
            # all below values for the planet orbit and RV data are from Butler2006.
#            planet_K = 461.1 #[m/s]
#            planet_p = 3.31246   #[days]
#            #planet_p_years = planet_p/365.25
#            planet_e = 0.023    
#            planet_argPeri = 188.0   #[deg]
#            planet_To = 2446957.8   #[JD]
            sigma_jitter = 3.0 #15.0    #[m/s]
            
#            ## solution from Narita's paper
#            argPeri_deg_proposed = 121.6
#            T_proposed = 2454777.94
#            e_proposed = 0.68
#            period_proposed = 0.0574949
#            a = 0.1614
#            (a_total, a1, a2) = semiMajorConverter(Mass1_proposed, Mass2_proposed, a_total=a,a1=0.0,a2=0.0)
#            a1s = [a1,a1]
            
            ##NOTE: version '3' uses the individual masses. version '4' uses the semi-major axis.
            ##      Both should be equivalent, but are not seeming so for some reason.
            chi_squared_RV = 0.0
            numEpochs_RV = 0.0
            for RVdataSet in range(0,len(RVs)):
#                chi_squared_RV_curr = rv2bodyCalculator3(RV_epochs[RVdataSet], RVs[RVdataSet], RVerrors[RVdataSet], sigma_jitter, Mass1_proposed, Mass2_proposed, inclination_deg_proposed, \
#                                                            period_proposed, e_proposed, T_proposed, argPeri_deg_proposed, \
#                                                                planet_K, planet_p, planet_e, planet_argPeri, planet_To)
                
#                chi_squared_RV_curr = rv2bodyCalculator4(RV_epochs[RVdataSet], RVs_proposed[RVdataSet], RVerrors[RVdataSet], sigma_jitter, Mass1_proposed, np.mean(a1s), inclination_deg_proposed, \
#                                                                   period_proposed, e_proposed, T_proposed, argPeri_deg_proposed, \
#                                                                       planet_K, planet_p, planet_e, planet_argPeri, planet_To)
                
                chi_squared_RV_curr = rv1bodyCalculator(RV_epochs[RVdataSet], RVs_proposed[RVdataSet], RVerrors[RVdataSet], sigma_jitter,\
                                                         inclination_deg_proposed, period_proposed, e_proposed, T_proposed, argPeri_deg_proposed,\
                                                          np.mean(a1s), verbose=False)
                
#                chi_squared_RV_curr = rv1bodyCalculator2(RV_epochs[RVdataSet], RVs_proposed[RVdataSet], RVerrors[RVdataSet], sigma_jitter,\
#                                                          period_proposed, e_proposed, T_proposed, argPeri_deg_proposed, Mass1_proposed, Mass2_proposed, verbose=verbose)
                
                chi_squared_RV = chi_squared_RV+chi_squared_RV_curr
                numEpochs_RV_curr = len(RV_epochs[RVdataSet])
                numEpochs_RV = numEpochs_RV+numEpochs_RV_curr
    
            # convert to a reduced chiSquared
            chi_squared_RV_reduced = (1.0/((1.0*numEpochs_RV)-totalNumModelParams))*chi_squared_RV
  
        #print 'chi_squared_RV_reduced = ',chi_squared_RV_reduced#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        # sum non-reduced chiSquareds
        chi_squared_total = chi_squared_DI + chi_squared_RV + gauss_chi_squared_total
        
        # Total Reduced Chi Squared
        totalNumDataPoints = 2.0*numEpochs_DI + 1.0*numEpochs_RV
        TotalNu = (totalNumDataPoints-totalNumModelParams)
        chi_squared_total_reduced = (1.0/TotalNu)*chi_squared_total
        #print 'chi_squared_total_reduced = ',chi_squared_total_reduced#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$        
        # decide if the proposed orbit is to be accepted
        if chi_squared_total_reduced < chi_squared_MAX:
            
            totalAcceptedCounter += 1
            #print '\nchi_squared_total_reduced = '+str(chi_squared_total_reduced)+' was found to be lower than chi_squared_MAX = '+str(chi_squared_MAX)#$$$$$$$$$$$$$$$$
            if chi_squared_total_reduced<reduced_chiSquareMin:
                reduced_chiSquareMin = chi_squared_total_reduced
                chi_squared_DI_reducedMIN = chi_squared_DI_reduced
                chi_squared_RV_reducedMIN = chi_squared_RV_reduced
                resultsDict['reduced_chiSquareMin_DI'] = chi_squared_DI_reducedMIN
                resultsDict['reduced_chiSquareMin_RV'] = chi_squared_RV_reducedMIN
                resultsDict['reduced_chiSquareMin'] = reduced_chiSquareMin
                resultsDict['longAN_deg_best'] = longAN_deg_best = longAN_deg_proposed
                resultsDict['e_best'] = e_best = e_proposed
                resultsDict['T_best'] = T_best = T_proposed
                resultsDict['period_best'] = period_best = period_proposed
                resultsDict['inclination_deg_best'] = inclination_deg_best = inclination_deg_proposed
                resultsDict['argPeri_deg_best'] = argPeri_deg_best = argPeri_deg_proposed
                resultsDict['a_total_best'] = a_total_best = a_total_proposed
                resultsDict['RV_origin_vel_0_best'] = RV_origin_vel_0_best = RV_origin_vel_0_proposed
                resultsDict['RV_origin_vel_1_best'] = RV_origin_vel_1_best = RV_origin_vel_1_proposed
            
            ## store output orbital elements of model
        
            # create string line of data to write to file
            line = str(longAN_deg_proposed)
            line = line +'   '+ str(e_proposed)
            line = line +'   '+ str(T_proposed)
            line = line +'   '+ str(period_proposed)
            line = line +'    '+ str(inclination_deg_proposed)
            line = line +'      '+ str(argPeri_deg_proposed)
            line = line +'   '+ str(a_total_proposed) # a_total=a2 if Mass1 or Mass2 are non-default (ie !=1)
            line = line +'   '+ str(chi_squared_total_reduced)
            line = line +'    '+ str(RV_origin_vel_0_proposed)
            line = line +'    '+ str(RV_origin_vel_1_proposed)
            line = line +'      '+ str(1)+'\n'
            # write final string to file
            f.write(line)
        else:
            if chi_squared_total_reduced<reduced_chiSquareMin:
                reduced_chiSquareMin = chi_squared_total_reduced
                chi_squared_DI_reducedMIN = chi_squared_DI_reduced
                chi_squared_RV_reducedMIN = chi_squared_RV_reduced
            #print '\nchi_squared_total_reduced = '+str(chi_squared_total_reduced)+' was NOT found to be lower than chi_squared_MAX = '+str(chi_squared_MAX)#$$$$$$$$$$$$$$$$
    # Done looping through all the samples
    # Finished writing all samples, so close the file data is being written to
    f.close()
    # Move data file from Temp dir to the Final dir 
    shutil.move(filenameTEMP,filenameFINAL)
    
    if not silent:
        print 'MCONLY: $$$$$$$$$$ SIMULATOR complete $$$$$$$$$$\n'
        print 'Best orbit found:'
        print "LongAN = ",resultsDict['longAN_deg_best']
        print "e = ",resultsDict['e_best']
        print "To = ",resultsDict['T_best']
        print "period = ",resultsDict['period_best']
        print "inclination = ",resultsDict['inclination_deg_best']
        print "argPeri = ",resultsDict['argPeri_deg_best']
        print "a_total = ",resultsDict['a_total_best']
        print "Reduced chiSquaredMin = ",resultsDict['reduced_chiSquareMin']
        print 'RV_origin_vel_0_best = ',resultsDict['RV_origin_vel_0_best']
        print 'RV_origin_vel_1_best = ',resultsDict['RV_origin_vel_1_best']
        print 'reduced_chiSquareMin_DI = ',resultsDict['reduced_chiSquareMin_DI']
        print 'reduced_chiSquareMin_RV = ',resultsDict['reduced_chiSquareMin_RV'] 
        startParamStr = '['+str(resultsDict['longAN_deg_best'])+','+str(resultsDict['e_best'])+','+str(resultsDict['T_best'])+\
                        ','+str(resultsDict['period_best'])+','+str(resultsDict['inclination_deg_best'])+','+str(resultsDict['argPeri_deg_best'])+','+str(resultsDict['a_total_best'])+']'
        print startParamStr
    # record the time the chain finished and print
    endTime = time.clock()
    totalTime = (endTime-startTime) # in seconds
    totalTimeString = timeString(totalTime)
    print 'MCONLY: The took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
    resultsDict['totalTime'] = totalTimeString
    
    # write finished results dict to file
    dictToFile(resultsDict, filenameTEMP[:-4]+'-resultsDict.txt')
    # Add summary print lines to results dict
    f = open(filenameTEMP[:-4]+'-resultsDict.txt','a')
    f.write('\n\n*** Summary prints during simulation ***')
    f.write(printLines)
    f.write('\n\n#### Best Orbit as startParam format ####\n\n'+startParamStr)
    f.close
    
    if not silent:
        print '\nMCONLY: $$$$$$$$$$$$$ FINISHED mcONLY3DorbSimUniform Sim $$$$$$$$$$$$$$$\n'
    
#****** DONE!! ********        
            
            

    
    