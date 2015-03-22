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
from orbitToolbox import multiEpochOrbCalcTH_I, timeString, semiMajorConverter, proposalMaxMinsCalc, ABCFG_MaxMins, ABCFG_MaxMins2, ABCFG_values, dictToFile, chiSquaredCalc, rv2bodyCalculator3, rv2bodyCalculator4
from paramSettingsDict import paramSettingsDict
    
def mcmcOrbSimUniformTH_I(filename, numSamples, temperature=1.0, sigmaPercent=1.0, startParams=False, silent=True, verbose=False,\
                        dataTempDir='/home/Kyle/DataTempForSims/', dataFinalDir='/run/media/Kyle/Data1/Todai_Work/workspace/Binary-project/data/'):
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
                                                               
    if not startParams:
        startRandom = True
    else:
        startRandom = False
        
    ##### model input parameter perturbing variables #########
    numSamplePrints = paramSettingsDict['General']['numSamplePrints']
    rateSamples = paramSettingsDict['General']['rateSamples'] # used to calculate acceptance rate only
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
    SA_arcsec_measured_REALs = paramSettingsDict['Data']['SAs']
    SA_mean_errors = paramSettingsDict['Data']['SAerrors']     
    PA_deg_measured_REALs = paramSettingsDict['Data']['PAs']
    PA_mean_errors = paramSettingsDict['Data']['PAerrors']
    DI_epochs = paramSettingsDict['Data']['DI_epochs']
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
    ######## Ranges of Thiele-Innes parameters ####
    #NOTE: must divide a_total values from above by sys_dist to convert from [AU]->["]
    (a_totalMax, a1MAX, a2MAX) = semiMajorConverter(Mass1, Mass2, a_total=(a_totalMAX/Sys_Dist_PC),a1=0.0,a2=0.0)
    (a_totalMIN, a1MIN, a2MIN) = semiMajorConverter(Mass1, Mass2, a_total=(a_totalMIN/Sys_Dist_PC),a1=0.0,a2=0.0)
    (Amax,Amin, Bmax,Bmin,  Cmax,Cmin, Fmax,Fmin, Gmax,Gmin) = \
        ABCFG_MaxMins2(a_totalMax, a_totalMIN, math.radians(argPeri_degMAX), math.radians(argPeri_degMIN), math.radians(longAN_degMAX), math.radians(longAN_degMIN), \
                                                                math.radians(inclination_degMAX), math.radians(inclination_degMIN))
        
    if not silent:
        print '\nMCMC: $$$$$$$$$$$$$$$$$$$  STARTING THE SIMULATOR  $$$$$$$$$$$$$$$$$'
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
    if DI_epochs[0]<RV_epochs[0]:
    	firstEpoch = DI_epochs[0]
    else:
    	firstEpoch = RV_epochs[0]
    TMAX = firstEpoch
    TMIN = TMAX-periodMAX*365.25
    
    # set up min and max for RV origin params
    RV_origin_vel_0_MIN = np.min(RVs[0])
    RV_origin_vel_0_MAX = np.max(RVs[0])
    RV_origin_vel_1_MIN = np.min(RVs[1])
    RV_origin_vel_1_MAX = np.max(RVs[1])
    
    # init input values
    if startRandom:
        longAN_deg_latest = rand.uniform(longAN_degMIN, longAN_degMAX)
        e_latest = rand.uniform(eMIN, eMAX)
        period_latest = rand.uniform(periodMIN, periodMAX)
        T_latest = rand.uniform(TMIN, TMAX)
        inclination_deg_latest = rand.uniform(inclination_degMIN, inclination_degMAX)
        argPeri_deg_latest = rand.uniform(argPeri_degMIN, argPeri_degMAX)
        a_total_latest = rand.uniform(a_totalMIN, a_totalMAX)
        RV_origin_vel_0_latest = rand.uniform(RV_origin_vel_0_MIN, RV_origin_vel_0_MAX)
        RV_origin_vel_1_latest = rand.uniform(RV_origin_vel_1_MIN, RV_origin_vel_1_MAX)
        (A_latest,B_latest,C_latest,F_latest,G_latest) = \
            ABCFG_values((a_total_latest/Sys_Dist_PC), math.radians(argPeri_deg_latest), math.radians(longAN_deg_latest), math.radians(inclination_deg_latest))
    else:
        longAN_deg_latest = startParams[0]
        e_latest = startParams[1]
        T_latest = startParams[2]
        period_latest = startParams[3]
        inclination_deg_latest = startParams[4]
        argPeri_deg_latest = startParams[5]
        a_total_latest = startParams[6]
        RV_origin_vel_0_latest = 0
        RV_origin_vel_1_latest = 0
        (A_latest,B_latest,C_latest,F_latest,G_latest) = \
            ABCFG_values((a_total_latest/Sys_Dist_PC), math.radians(argPeri_deg_latest), math.radians(longAN_deg_latest), math.radians(inclination_deg_latest))
    
    reduced_chiSquared_latest = 20.0
    chiSquare_latest = (2.0*numEpochs_DI-7.0)*reduced_chiSquared_latest
    timesBeenHere = 1
    paramBeingVaried = 0
    
    # sigma is a percentage of total range, so convert its current percentage to a range value
    A_sigma = (sigmaPercent/100.0)*(Amax-Amin)
    B_sigma = (sigmaPercent/100.0)*(Bmax-Bmin)
    C_sigma = (sigmaPercent/100.0)*(Cmax-Cmin)
    F_sigma = (sigmaPercent/100.0)*(Fmax-Fmin)
    G_sigma = (sigmaPercent/100.0)*(Gmax-Gmin)
    e_sigma = (sigmaPercent/100.0)*(eMAX-eMIN)
    T_sigma = (sigmaPercent/100.0)*(TMAX-TMIN)
    RV_0_sigma = (sigmaPercent/100.0)*(RV_origin_vel_0_MIN-RV_origin_vel_0_MAX)
    RV_1_sigma = (sigmaPercent/100.0)*(RV_origin_vel_1_MIN-RV_origin_vel_1_MAX)
#    #$$$$$$$$$$$$$$$$$$$
#    print 'a_total_latest = ',a_total_latest
#    print '(a_total_latest/Sys_Dist_PC) = ',str((a_total_latest/Sys_Dist_PC))
#    print 'A_sigma = ',A_sigma
#    print 'B_sigma = ',B_sigma
#    print 'C_sigma = ',C_sigma
#    print 'F_sigma = ',F_sigma
#    print 'G_sigma = ',G_sigma
    
    # setting initially 'proposed' states equal to initial 'latest' states
    A_proposed = A_latest
    B_proposed = B_latest
    C_proposed = C_latest
    F_proposed = F_latest
    G_proposed = G_latest
    #$$$$$$$$$$$$$$$$$$$
#    print 'A_latest = ',A_latest
#    print 'B_latest = ',B_latest
#    print 'C_latest = ',C_latest
#    print 'F_latest = ',F_latest
#    print 'G_latest = ',G_latest
    longAN_deg_proposed = longAN_deg_latest
    e_proposed = e_latest
    T_proposed = T_latest
    period_proposed = period_latest
    inclination_deg_proposed = inclination_deg_latest
    argPeri_deg_proposed = argPeri_deg_latest
    a_total_proposed = a_total_latest
    RV_origin_vel_0_proposed = RV_origin_vel_0_latest
    RV_origin_vel_1_proposed = RV_origin_vel_1_latest
    
    # variables for the success rate print block in chain loop
    printTime = numSamples/numSamplePrints
    printCount = 0
    printsDone = 0
    totalAcceptedCounter = 0
    rateAcceptedCounter = 0
    reduced_chiSquareMin = 1e12
    chi_squared_DI_reducedMIN = 1e12
    chi_squared_RV_reducedMIN = 1e12
    bestOrbit = 0
    likelihood_ratio = 0
    alpha = 0
    acceptRates = []
    acceptRates.append(0)
    rateTotalCounter = 0
    resultsDict = {}
    writtenToFile = False
    
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
    f.write('longAN [deg]      e [N/A]       To [julian date]  period [yrs]   inclination [deg]   argPeri [deg]   a_total [AU]  RVorigin_0 [m/s]   RVorigin_1 [m/s]   chiSquared    timesBeenHere\n')
    #else:
    #    # ie. non-default, so store individual a1 and a2
    #    f.write('longAN [deg]      e [N/A]     To [julian date]    period [yrs]   inclination [deg]   argPeri [deg]   a1 [AU]    a2[AU]   chiSquared  timesBeenHere\n')
    
    ## Start to create samples and check if they are are accepted
    ## ie. Start the chain
    for curSample in range(1,numSamples+1):
        # reset boolean for current orbit being written to file as False
        # as it hasn't been written yet.
        writtenToFile = False
        
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
                reduced_chiSquareMinStr = "Lowest reduced chiSquareds so far:  DI = " +str(chi_squared_DI_reducedMIN)+', RV = '+str(chi_squared_RV_reducedMIN)+' , total = '+str(reduced_chiSquareMin)
                line = "\n"+str(printsDone)+'/'+str(int(numSamplePrints))+' done at '+timeStr+', running time: '+timeSoFarString+', '+str(totalAcceptedCounter)+' / '+str(curSample)+' successful. '   
                line = line+'\n'+"Latest acceptance rate = "+str(acceptRates[-1])+" and sigmaPercent = "+str(sigmaPercent)+ " and temperature = "+str(temperature)
                line = line+'\n'+reduced_chiSquareMinStr
                line = line+'\n'+"Latest likelihood_ratio = "+str(likelihood_ratio)+" and alpha = "+str(alpha)
                line = line+'\n'+"Latest chiSquare_latest - chi_squared_total_cur = "+str(chiSquare_latest - chi_squared_total_cur)
                line = line+'\n'+'latest reduced_chi-square = '+str(chi_squared_total_reduced)
                line = line+'\n'+"inclination_deg_proposed = "+str(inclination_deg_proposed) +'\n'
                print line
                printLines = printLines+line
                #print "\n"+str(printsDone)+'/'+str(int(numSamplePrints))+' completed at '+timeStr+'. '+str(totalAcceptedCounter)+' successful out of '+str(curSample)+'. '    
                #print reduced_chiSquareMinStr
                #print "Latest acceptance rate = "+str(acceptRates[-1])+" and sigmaPercent = "+str(sigmaPercent)+ " and temperature = "+str(temperature)
                #print "Latest likelihood_ratio = "+str(likelihood_ratio)+" and alpha = "+str(alpha)
                #print "Latest chiSquare_latest - chi_squared_total_cur = "+str(chiSquare_latest - chi_squared_total_cur)
                #print 'latest reduced_chi-square = ',chi_squared_total_reduced
                #print "inclination_deg_proposed = ",inclination_deg_proposed
                
        ## Draw the values that have error but are not Orbital Elements (ie. the sys_dist and obj masses)
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
            
        ## propose a new set of input parameters by perturbing the current ones 
        if paramBeingVaried==0:
            # A
            (min,max) = proposalMaxMinsCalc(A_latest, A_sigma, Amin, Amax)
            A_proposed = rand.uniform(min,max)
        
        elif paramBeingVaried==1:
            # B
            (min,max) = proposalMaxMinsCalc(B_latest, B_sigma, Bmin, Bmax)
            B_proposed = rand.uniform(min,max)
        
        elif paramBeingVaried==2:
            # C
            (min,max) = proposalMaxMinsCalc(C_latest, C_sigma, Cmin, Cmax)
            C_proposed = rand.uniform(min,max)
        
        elif paramBeingVaried==3:
            # F
            (min,max) = proposalMaxMinsCalc(F_latest, F_sigma, Fmin, Fmax)
            F_proposed = rand.uniform(min,max)
        
        elif paramBeingVaried==4:
            # G
            (min,max) = proposalMaxMinsCalc(G_latest, G_sigma, Gmin, Gmax)
            G_proposed = rand.uniform(min,max)
        
        elif paramBeingVaried==5:
            # eccentricity
            (min,max) = proposalMaxMinsCalc(e_latest, e_sigma, eMIN, eMAX)
            e_proposed = rand.uniform(min,max)
        
        elif paramBeingVaried==6:
            # time of last periapsis
            (min,max) = proposalMaxMinsCalc(T_latest, T_sigma, TMIN, TMAX)
            T_proposed = rand.uniform(min, max)        
        
        ## send random parameters along with known ones to multi-epoch orbit calculator
        if RVonly:
            (a_arcsec_proposed, argPeri_rad_proposed, longAN_rad_proposed, inclination_rad_proposed, period_proposed, chi_squared_total_cur) = \
            multiEpochOrbCalcTH_I(e_proposed, T_proposed, A_proposed, B_proposed, C_proposed, F_proposed, G_proposed, Mass1_proposed, \
             Mass2_proposed, Sys_Dist_PC, [SA_arcsec_measured_REALs[0]], [PA_deg_measured_REALs[0]], [DI_epochs[0]], [SA_mean_errors[0]], [PA_mean_errors[0]])
        else:
            (a_arcsec_proposed, argPeri_rad_proposed, longAN_rad_proposed, inclination_rad_proposed, period_proposed, chi_squared_total_cur) = \
            multiEpochOrbCalcTH_I(e_proposed, T_proposed, A_proposed, B_proposed, C_proposed, F_proposed, G_proposed, Mass1_proposed, \
             Mass2_proposed, Sys_Dist_PC, SA_arcsec_measured_REALs, PA_deg_measured_REALs, DI_epochs, SA_mean_errors, PA_mean_errors)
         
        # calculate a total 
        a_total_proposed = a_arcsec_proposed*Sys_Dist_proposed
    	
    	# convert angles in radians to degrees
        argPeri_deg_proposed = math.degrees(argPeri_rad_proposed)
        inclination_deg_proposed = math.degrees(inclination_rad_proposed)
        longAN_deg_proposed = math.degrees(longAN_rad_proposed)
        
        
#        #$$$$$$$$$$$$$$$$$$$
#        print 'longAN_deg_proposed = ',longAN_deg_proposed
#        print 'inclination_deg_proposed = ',inclination_deg_proposed
#        print 'argPeri_deg_proposed = ',argPeri_deg_proposed
#        print 'a_total_proposed = ',a_total_proposed
#        print 'e_proposed = ',e_proposed
#        print 'T_proposed = ',T_proposed
#        print 'A_proposed = ',A_proposed
#        print 'B_proposed = ',B_proposed
#        print 'C_proposed = ',C_proposed
#        print 'F_proposed = ',F_proposed
#        print 'G_proposed = ',G_proposed
        
        
        # Check if the parameters outside of the Thiele-Innes ones are in ranges
        # this is required as I can't get the max/min of Thiele-Innes params to 
        # make sure the outputs are in range, although this is impossible if you 
        # realize that there is too much contributing to the TH_I params to make them
        # properly set output ranges.
        AllInRange = True
#        if (longAN_rad_proposed>longAN_degMAX) or (longAN_rad_proposed<longAN_degMIN):
#        	AllInRange = False
        if (inclination_deg_proposed>inclination_degMAX) or (inclination_deg_proposed<inclination_degMIN):
         	AllInRange = False
#        if (argPeri_deg_proposed>argPeri_degMAX) or (argPeri_deg_proposed<argPeri_degMIN):
#        	AllInRange = False
#        if (a_total_proposed>a_totalMAX) or (a_total_proposed<a_totalMIN):
#        	AllInRange = False
		if (period_proposed>periodMAX) or (period_proposed<periodMIN):
			AllInRange = False
        
        ## if ALL the parameters are in range continue
		## else do nothing and go back to top of loop and try again.
        if AllInRange:
            chi_squared_DI = chi_squared_total_cur
            # convert to a reduced chiSquared
            chi_squared_DI_reduced = (1.0/((2.0*(numEpochs_DI))-6.0))*chi_squared_DI
            # if only saving orbits based on RV data, set all DI chiSquared related values to 0.0
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
                totalNumModelParams = 6.0 
            else:
                totalNumModelParams = 8.0 
                # propose RV origin velocity for each dataset
                (min,max) = proposalMaxMinsCalc(RV_origin_vel_0_latest, RV_0_sigma, RV_origin_vel_0_MIN, RV_origin_vel_0_MAX)
                RV_origin_vel_0_proposed = rand.uniform(min, max)
                
                (min,max) = proposalMaxMinsCalc(RV_origin_vel_1_latest, RV_1_sigma, RV_origin_vel_1_MIN, RV_origin_vel_1_MAX)
                RV_origin_vel_1_proposed = rand.uniform(min, max)
            
                # update RV values in each dataset for input to rv2bodyCalculator below
                RVs_proposed = RVs
                for epoch in range(0,len(RVs_proposed[0])):
                    RVs_proposed[0][epoch] = RVs[0][epoch]+RV_origin_vel_0_proposed
                for epoch in range(0,len(RVs_proposed[1])):
                    RVs_proposed[1][epoch] = RVs[1][epoch]+RV_origin_vel_1_proposed
                
                # data for planet's orbit that are needed to calculate proper chiSquare for RV fit
                # all below values for the planet orbit and RV data are from Butler2006.
                planet_K = 461.1 #[m/s]
                planet_p = 3.31246   #[days]
                #planet_p_years = planet_p/365.25
                planet_e = 0.023    
                planet_argPeri = 188.0   #[deg]
                planet_To = 2446957.8   #[JD]
                sigma_jitter = 15.0    #[m/s]
            
                ##NOTE: version '3' uses the individual masses. version '4' uses the semi-major axis.
                ##      Both should be equivalent, but are not seeming so for some reason.
                # calculate a1 from a_total and masses
                (a_total, a1, a2) = semiMajorConverter(Mass1, Mass2, a_total=a_total_proposed,a1=0.0,a2=0.0)
#                chi_squared_RV = rv2bodyCalculator3(RV_epochs, RVs, RVerrors, sigma_jitter, Mass1_proposed, Mass2_proposed, inclination_deg_proposed, \
#                                                                period_proposed, e_proposed, T_proposed, argPeri_deg_proposed, \
#                                                                    planet_K, planet_p, planet_e, planet_argPeri, planet_To)
                chi_squared_RV = 0.0
                numEpochs_RV = 0.0
                for RVdataSet in range(1,len(RVs)):
                   chi_squared_RV_curr = rv2bodyCalculator4(RV_epochs[RVdataSet], RVs_proposed[RVdataSet], RVerrors[RVdataSet], sigma_jitter, Mass1_proposed, a1, inclination_deg_proposed, \
                                                               period_proposed, e_proposed, T_proposed, argPeri_deg_proposed, \
                                                                   planet_K, planet_p, planet_e, planet_argPeri, planet_To)
                   chi_squared_RV = chi_squared_RV+chi_squared_RV_curr
                   numEpochs_RV_curr = len(RV_epochs[RVdataSet])
                   numEpochs_RV = numEpochs_RV+numEpochs_RV_curr
		
		
		        # convert to a reduced chiSquared
                chi_squared_RV_reduced = (1.0/((1.0*numEpochs_RV)-6.0))*chi_squared_RV
            
            
            
            # sum non-reduced chiSquareds
            #chi_squared_total = chi_squared_DI + chi_squared_RV + gauss_chi_squared_total
            #weighting scheme chi squared total
            chi_squared_total = (chi_squared_DI*3.0) + (chi_squared_RV/3.0)
            
            # Total Reduced Chi Squared
            totalNumDataPoints = 2.0*numEpochs_DI + 1.0*numEpochs_RV 
            TotalNu = (totalNumDataPoints-totalNumModelParams)
            chi_squared_total_reduced = (1.0/TotalNu)*chi_squared_total      
            
            #print 'chi_squared_total_reduced = ',chi_squared_total_reduced #$$$$$$$$$$$$$$$$$$$
            
            ## calculate the likelihood ratio for the previous(latest) and current (proposed) models/states
            likelihood_ratio=0.0
            try:
                likelihood_ratio = math.exp((chiSquare_latest - chi_squared_total_reduced)/ (2.0*temperature))#in the past this was the non-reduced version, not sure if reduced is allowed by I like it better
            except:
                blah=1 #$$$$$$$$$
                #print 'chiSquare_latest = ',chiSquare_latest
                #print 'chi_squared_total_cur = ',chi_squared_total_cur
                #print 'temperature = ',temperature
            
            alpha = rand.uniform(0.0, 1.0)
            
            if likelihood_ratio>1.0:
                RHS = 1.0
            else:
                RHS = likelihood_ratio
            
            if alpha < RHS:
                # update counters for acceptance total and rate acceptance
                rateAcceptedCounter += 1
                totalAcceptedCounter += 1
                
                # update boolean for latest orbit info being written to file
                writtenToFile = True
                
                #update latest values with ones just proposed
                # inputs
                A_latest = A_proposed
                B_latest = B_proposed
                C_latest = C_proposed #$$$$$$$$$$$$$$$$$$$ aren't using this one anymore
                F_latest = F_proposed
                G_latest = G_proposed
                longAN_deg_latest = longAN_deg_proposed
                e_latest = e_proposed
                T_latest = T_proposed
                #outputs
                period_latest = period_proposed
                inclination_deg_latest = inclination_deg_proposed
                argPeri_deg_latest = argPeri_deg_proposed
                a_total_latest = a_total_proposed
                chiSquare_latest = chi_squared_total_reduced #in the past this was the non-reduced version, not sure if reduced is allowed by I like it better               
                RV_origin_vel_0_latest = RV_origin_vel_0_proposed
                RV_origin_vel_1_latest = RV_origin_vel_1_proposed
                
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
                if chi_squared_total_reduced<chi_squared_MAX:  #################################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ NOTE REALLY LEGAL, BUT SAVES MEMORY
                    # create string line of data to write to file
                    line = str(longAN_deg_proposed)
                    line = line +'   '+ str(e_proposed)
                    line = line +'   '+ str(T_proposed)
                    line = line +'   '+ str(period_proposed)
                    line = line +'    '+ str(inclination_deg_proposed)
                    line = line +'      '+ str(argPeri_deg_proposed)
                    line = line +'   '+ str(a_total_proposed) 
                    line = line +'   '+ str(chi_squared_total_reduced)
                    line = line +'    '+ str(RV_origin_vel_0_proposed)
                    line = line +'    '+ str(RV_origin_vel_1_proposed)
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
    ## store output orbital elements of model that was the last one in the loop.
    ## if they were not written above inside the loop.  This happens when the 
    ## last accepted state had further proposed orbits after without any being accepted.
    
    # create string line of data to write to file
    if writtenToFile==False:
        line = str(longAN_deg_proposed)
        line = line +'   '+ str(e_proposed)
        line = line +'   '+ str(T_proposed)
        line = line +'   '+ str(period_proposed)
        line = line +'    '+ str(inclination_deg_proposed)
        line = line +'      '+ str(argPeri_deg_proposed)
        line = line +'   '+ str(a_total_proposed) # a_total=a2 if Mass1 or Mass2 are non-default (ie !=1)
        line = line +'   '+ str(chi_squared_total_cur)
        line = line +'    '+ str(RV_origin_vel_0_proposed)
        line = line +'    '+ str(RV_origin_vel_1_proposed)
        line = line +'      '+ str(timesBeenHere)+'\n'
        # write final string to file
        f.write(line)
        
    ## Finished writing all samples, so close the file data is being written to
    f.close()
    # Move data file from Temp dir to the Final dir 
    shutil.move(filenameTEMP,filenameFINAL)
    
    print 'MCMC: '+str(totalAcceptedCounter)+' steps were taken during the MCMC simulator'
    resultsDict['meanAcceptanceRate'] = np.mean(acceptRates)
    resultsDict['totalAcceptedCounter'] = totalAcceptedCounter
    
    if not silent:
        print 'MCMC: $$$$$$$$$$ SIMULATOR complete $$$$$$$$$$\n'
        print 'Best orbit found:'
        print "LongAN = ",resultsDict['longAN_deg_best']
        print "e = ",resultsDict['e_best']
        print "To = ",resultsDict['T_best']
        print "period = ",resultsDict['period_best']
        print "inclination = ",resultsDict['inclination_deg_best']
        print "argPeri = ",resultsDict['argPeri_deg_best']
        print "a_total = ",resultsDict['a_total_best']
        print "Reduced chiSquaredMin = ",resultsDict['reduced_chiSquareMin']
        print "\nmean acceptance rate = ",resultsDict['meanAcceptanceRate']
        print 'RV_origin_vel_0_best = ',resultsDict['RV_origin_vel_0_best']
        print 'RV_origin_vel_1_best = ',resultsDict['RV_origin_vel_1_best']
        print 'reduced_chiSquareMin_DI = ',resultsDict['reduced_chiSquareMin_DI']
        print 'reduced_chiSquareMin_RV = ',resultsDict['reduced_chiSquareMin_RV'] 
        startParamStr = '['+str(resultsDict['longAN_deg_best'])+','+str(resultsDict['e_best'])+','+str(resultsDict['T_best'])+\
                        ','+str(resultsDict['period_best'])+','+str(resultsDict['inclination_deg_best'])+','+\
                        str(resultsDict['argPeri_deg_best'])+','+str(resultsDict['a_total_best'])+\
                        ','+str(resultsDict['RV_origin_vel_0_best'])+','+str(resultsDict['RV_origin_vel_1_best'])+']'
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
        print '\nMCMC: $$$$$$$$$$$$$ FINISHED MCMC-UniformTH_I Sim $$$$$$$$$$$$$$$\n'
    
#****** DONE!! ********        
            
            

    
    