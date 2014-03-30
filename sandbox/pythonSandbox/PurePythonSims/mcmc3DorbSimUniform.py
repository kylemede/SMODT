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
from orbitToolbox import multiEpochOrbCalc, multiEpochOrbCalc3, timeString, dictToFile, chiSquaredCalc,rv1bodyCalculator, rv2bodyCalculator3, rv2bodyCalculator4, proposalMaxMinsCalc
from paramSettingsDict import paramSettingsDict
    
def mcmcOrbSimUniform(filename, numSamples, temperature=1.0, sigmaPercent=1.0, startParams=False, silent=True, verbose=False,\
                        dataTempDir='/home/Kyle/DataTempForSims/', dataFinalDir='/media/Data1/Todai_Work/workspace/Binary-project/data/'):
    #*********** In this brick is all you need to play with *************    
    #********************************************************************
    ##### Overall simulation parameters ########
    #numSamples = int(1e5)
    #dataTempDir = '/media/Data1/Todai_Work/workspace/Binary-project/data/'#'/home/Kyle/DataTempForSims/' 
    #dataFinalDir = '/media/Data1/Todai_Work/workspace/Binary-project/data/'
    #filename = 'SimAneal7_2-Tau-Boo-ALL-T10-TEMP'  #This is a feature implemented as I have an SSD which is much 
    #                                            #MUCH faster at reading/writing for use during the sim
    #                                            # and the finished file is then moved to the storage dir on 
    #                                            # my standard HDD at the end of the sim.
    #                                            # making the Temp and Final the same negates this if you don't have an SSD.
    #silent = False
    #verbose = True                                                               
    numSamplePrints = 10.0     
    #temperature = 50.0
    if not startParams:
        startRandom = True
    else:
        startRandom = False
    ##### model input parameter perturbing variables #########
    #sigmaPercent = 1.0  # fixed sigma value
    rateSamples = 1000 # used to calculate acceptance rate only
    chi_squared_MAX = paramSettingsDict['General']['chi_squared_MAX']
    gaussMassAndDist = paramSettingsDict['General']['gaussMassAndDist'] # calculate the mass and sys_dist values using gaussian dists?
    #NOTE: leave both of these 'False' to run standard 3D simulations
    DIonly = paramSettingsDict['General']['DIonly'] #setting this makes RV chiSquared=0
    RVonly = paramSettingsDict['General']['RVonly'] #setting this makes DI chiSquared=0
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
    TMIN = TMAX-periodMAX*365.242
    
    # set up min and max for RV origin params
    RV_origin_vel_0_MIN = np.min(RVs[0])
    RV_origin_vel_0_MAX = np.max(RVs[0])
    RV_origin_vel_1_MIN = -250#np.min(RVs[1])
    RV_origin_vel_1_MAX = -100#np.max(RVs[1])
    RV_origin_vel_0s = []
    RV_origin_vel_1s = []
    
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
    else:
        longAN_deg_latest = startParams[0]
        e_latest = startParams[1]
        T_latest = startParams[2]
        period_latest = startParams[3]
        inclination_deg_latest = startParams[4]
        argPeri_deg_latest = startParams[5]
        a_total_latest = startParams[6]
        try:
            RV_origin_vel_0_latest = startParams[7]
            RV_origin_vel_1_latest = startParams[8]
        except:
            RV_origin_vel_0_latest = 0
            RV_origin_vel_1_latest = 0
        
    
    sine_inc_latest = math.sin(math.radians(inclination_deg_latest))
    reduced_chiSquared_latest = 20.0
    chiSquare_latest = (2.0*numEpochs_DI-7.0)*reduced_chiSquared_latest
    timesBeenHere = 1
    paramBeingVaried = 0
    
    # set max and min for sin(inclination_deg)
    if inclination_degMIN>90.0:
        sine_inc_MAX = math.sin(math.radians(inclination_degMIN))        
        sine_inc_MIN = math.sin(math.radians(inclination_degMAX))  
    else:
        sine_inc_MIN = math.sin(math.radians(inclination_degMIN))        
        sine_inc_MAX = math.sin(math.radians(inclination_degMAX))
    
    # sigma is a percentage of total range, so convert its current percentage to a range value
    sine_inc_sigma = (sigmaPercent/100.0)*(sine_inc_MAX-sine_inc_MIN)
    e_sigma = (sigmaPercent/100.0)*(eMAX-eMIN)
    longAN_sigma = (sigmaPercent/100.0)*(longAN_degMAX-longAN_degMIN)
    period_sigma = (sigmaPercent/100.0)*(periodMAX-periodMIN)
    argPeri_sigma = (sigmaPercent/100.0)*(argPeri_degMAX-argPeri_degMIN)
    a_total_sigma = (sigmaPercent/100.0)*(a_totalMAX-a_totalMIN)
    T_sigma = (sigmaPercent/100.0)*(TMAX-TMIN)
    RV_0_sigma = (sigmaPercent/100.0)*(RV_origin_vel_0_MIN-RV_origin_vel_0_MAX)
    RV_1_sigma = (sigmaPercent/100.0)*(RV_origin_vel_1_MIN-RV_origin_vel_1_MAX)
    
    # setting initially 'proposed' states equal to initial 'latest' states
    longAN_deg_proposed = longAN_deg_latest
    e_proposed = e_latest
    T_proposed = T_latest
    period_proposed = period_latest
    inclination_deg_proposed = inclination_deg_latest
    sine_inc_proposed = sine_inc_latest
    argPeri_deg_proposed = argPeri_deg_latest
    a_total_proposed = a_total_latest
    RV_origin_vel_0_proposed = RV_origin_vel_0_latest
    RV_origin_vel_1_proposed = RV_origin_vel_1_latest
        
    # variables for the success rate print block in chain loop
    # other params here are just initial values that will be replaced soon in the loop
    printTime = numSamples/numSamplePrints
    printCount = 0
    printsDone = 0
    totalAcceptedCounter = 0
    rateAcceptedCounter = 0
    reduced_chiSquareMin = 10000000.0
    chi_squared_DI_reducedMIN = 10000000.0
    chi_squared_RV_reducedMIN = 10000000.0
    bestOrbit = 0
    acceptRates = []
    acceptRates.append(0)
    rateTotalCounter = 0
    resultsDict = {}
    writtenToFile = False
    chi_squared_total_reduced = '?' 
    inclination_deg_proposed ='?'
    likelihood_ratio = '?'
    alpha = '?'
    chi_squared_total_cur = chi_squared_DI_reducedMIN
    
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
                line = line+'\n'+reduced_chiSquareMinStr
                line = line+'\n'+"Latest acceptance rate = "+str(acceptRates[-1])+" and sigmaPercent = "+str(sigmaPercent)+ " and temperature = "+str(temperature)
                line = line+'\n'+"Latest likelihood_ratio = "+str(likelihood_ratio)+" and alpha = "+str(alpha)
                line = line+'\n'+"Latest chiSquare_latest - chi_squared_total_cur = "+str(chiSquare_latest - chi_squared_total_cur)
                line = line+'\n'+'latest reduced_chi-square = '+str(chi_squared_total_reduced)
                line = line+'\n'+"inclination_deg_proposed = "+str(inclination_deg_proposed) +'\n'
                print line
                printLines = printLines+line
                
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
            # inclination
            # NOTE: using uniform in sin(i) rather than i as it then avoids even sampling
            #       of near face-on orbits that are rare or unlikely. 
            if (sine_inc_latest+sine_inc_sigma)>=sine_inc_MAX:
                max = sine_inc_MAX
                min = sine_inc_latest-sine_inc_sigma
            elif (sine_inc_latest-sine_inc_sigma)<=sine_inc_MIN:
                max = sine_inc_latest+sine_inc_sigma
                min = sine_inc_MIN
            else:
                max = sine_inc_latest + sine_inc_sigma
                min = sine_inc_latest - sine_inc_sigma
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
        
        elif paramBeingVaried==1:
            # long of accending node
            if (longAN_deg_latest+longAN_sigma)>=longAN_degMAX:
                max = longAN_degMAX
                min = longAN_deg_latest-longAN_sigma
            elif (longAN_deg_latest-longAN_sigma)<=longAN_degMIN:
                max = longAN_deg_latest+longAN_sigma
                min = longAN_degMIN
            else:
                max = longAN_deg_latest + longAN_sigma
                min = longAN_deg_latest - longAN_sigma
            longAN_deg_proposed = rand.uniform(min, max)
        
        elif paramBeingVaried==2:
            #eccentricity
            if (e_latest+e_sigma)>=eMAX:
                max = eMAX
                min = e_latest-e_sigma
            elif (e_latest-e_sigma)<=eMIN:
                max = e_latest+e_sigma
                min = eMIN
            else:
                max = e_latest + e_sigma
                min = e_latest - e_sigma
            e_proposed = rand.uniform(min,max)
        
        elif paramBeingVaried==3:
            #period
#            if (period_latest+period_sigma)>=periodMAX:
#                max = periodMAX
#                min = period_latest-period_sigma
#            elif (period_latest-period_sigma)<=periodMIN:
#                max = period_latest+period_sigma
#                min = periodMIN
#            else:
#                max = period_latest + period_sigma
#                min = period_latest - period_sigma
#            period_proposed = rand.uniform(min, max)
            period_proposed = 21.2165/365.242
            
        elif paramBeingVaried==4:
            # time of last periapsis
#            if (T_latest+T_sigma)>=TMAX:
#                max = TMAX
#                min = T_latest-T_sigma
#            elif (T_latest-T_sigma)<=TMIN:
#                max = T_latest+T_sigma
#                min = TMIN
#            else:
#                max = T_latest + T_sigma
#                min = T_latest - T_sigma
#            T_proposed = rand.uniform(min, max)
            T_proposed = 2454777.94
    
        elif paramBeingVaried==5:
            # argument of perigie
            if (argPeri_deg_latest+argPeri_sigma)>=argPeri_degMAX:
                max = argPeri_degMAX
                min = argPeri_deg_latest-argPeri_sigma
            elif (argPeri_deg_latest-argPeri_sigma)<=argPeri_degMIN:
                max = argPeri_deg_latest+argPeri_sigma
                min = argPeri_degMIN
            else:
                max = argPeri_deg_latest + argPeri_sigma
                min = argPeri_deg_latest - argPeri_sigma
            argPeri_deg_proposed = rand.uniform(min, max)
#        
#        elif paramBeingVaried==6:
#            # Total semi-major axis
#            if (a_total_latest+a_total_sigma)>=a_totalMAX:
#                max = a_totalMAX
#                min = a_total_latest-a_total_sigma
#            elif (a_total_latest-a_total_sigma)<=a_totalMIN:
#                max = a_total_latest+a_total_sigma
#                min = a_totalMIN
#            else:
#                max = a_total_latest + a_total_sigma
#                min = a_total_latest - a_total_sigma
#            a_total_proposed = rand.uniform(min, max)
        
        
        ## send random parameters along with known ones to multi-epoch orbit calculator
#        if RVonly:
#            (chi_squared_total_cur, ns, Ms, Es, thetas, xs, ys, a1s, a2s) =\
#            multiEpochOrbCalc3([SA_arcsec_measured_REALs[0]], [SA_mean_errors[0]], [PA_deg_measured_REALs[0]], [PA_mean_errors[0]],\
#                       [DI_epochs[0]], Sys_Dist_proposed, inclination_deg_proposed, longAN_deg_proposed, e_proposed, T_proposed, \
#                       period_proposed, argPeri_deg_proposed, a_total=False, Mass1=Mass1_proposed, Mass2=Mass2_proposed, verbose=False)
#        else:
#            (chi_squared_total_cur, ns, Ms, Es, thetas, xs, ys, a1s, a2s) =\
#            multiEpochOrbCalc3(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors,\
#                       DI_epochs, Sys_Dist_proposed, inclination_deg_proposed, longAN_deg_proposed, e_proposed, T_proposed, \
#                       period_proposed, argPeri_deg_proposed, a_total=False, Mass1=Mass1_proposed, Mass2=Mass2_proposed, verbose=False)
        if RVonly:
            (chi_squared_total_cur, ns, Ms, Es, thetas, Sep_Dists, SA_arcsec_measured_models, PA_deg_measured_models, a1s, a2s) =\
            multiEpochOrbCalc([SA_arcsec_measured_REALs[0]], [SA_mean_errors[0]], [PA_deg_measured_REALs[0]], [PA_mean_errors[0]],\
                       [DI_epochs[0]], Sys_Dist_proposed, inclination_deg_proposed, longAN_deg_proposed, e_proposed, T_proposed, \
                       period_proposed, argPeri_deg_proposed, a_total=False, Mass1=Mass1_proposed, Mass2=Mass2_proposed, verbose=False)
        else:
            (chi_squared_total_cur, ns, Ms, Es, thetas, Sep_Dists, SA_arcsec_measured_models, PA_deg_measured_models, a1s, a2s) =\
            multiEpochOrbCalc(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors,\
                       DI_epochs, Sys_Dist_proposed, inclination_deg_proposed, longAN_deg_proposed, e_proposed, T_proposed, \
                       period_proposed, argPeri_deg_proposed, a_total=False, Mass1=Mass1_proposed, Mass2=Mass2_proposed, verbose=False)
            
        # calculate a total from a1 plus a2
        a_total_proposed = np.mean(a1s)+np.mean(a2s)
    
        chi_squared_DI = chi_squared_total_cur
        # convert to a reduced chiSquared
        chi_squared_DI_reduced = (1.0/((2.0*(numEpochs_DI))-6.0))*chi_squared_DI
        
        DIweight = 1.0#(2.0/3.0)
        RVweight = 1.0#(1.0/3.0)
        if RVonly:
            chi_squared_DI_reduced = 0.0 
            chi_squared_DI = 0.0
            numEpochs_DI = 0.0
            DIweight = 1.0
            RVweight = 1.0
        # if only saving orbits based on DI data, set all RV chiSquared related values to 0.0
        # Else, go ahead and calculate the proper values for the RV data.
        if DIonly:
            chi_squared_RV_reduced = 0.0
            chi_squared_RV = 0.0
            numEpochs_RV = 0.0
            RV_origin_vel_0_proposed = 0
            RV_origin_vel_1_proposed = 0
            totalNumModelParams = 6.0 
            DIweight = 1.0
            RVweight = 1.0
        else:
            totalNumModelParams = 8.0 
            # propose RV origin velocity for each dataset
            #(min,max) = proposalMaxMinsCalc(RV_origin_vel_0_latest, RV_0_sigma, RV_origin_vel_0_MIN, RV_origin_vel_0_MAX)
            RV_origin_vel_0_proposed = 0#rand.uniform(min, max)
            
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
            #planet_p_years = planet_p/365.242
            planet_e = 0.023    
            planet_argPeri = 188.0   #[deg]
            planet_To = 2446957.8   #[JD]
            sigma_jitter = 15.0    #[m/s]
    
            ##NOTE: version '3' uses the individual masses. version '4' uses the semi-major axis.
            ##      Both should be equivalent, but are not seeming so for some reason.
            chi_squared_RV = 0.0
            numEpochs_RV = 0.0
            for RVdataSet in range(1,len(RVs)):
#                chi_squared_RV_curr = rv2bodyCalculator3(RV_epochs[RVdataSet], RVs_proposed[RVdataSet], RVerrors[RVdataSet], sigma_jitter, \
#                                                Mass1_proposed, Mass2_proposed, inclination_deg_proposed, \
#                                                        period_proposed, e_proposed, T_proposed, argPeri_deg_proposed, \
#                                                            planet_K, planet_p, planet_e, planet_argPeri, planet_To, verbose=verbose)
#                chi_squared_RV = rv2bodyCalculator4(RV_epochs[RVdataSet], RVs_proposed[RVdataSet], RVerrors[RVdataSet], sigma_jitter, \
#                                                Mass1_proposed, np.mean(a1s), inclination_deg_proposed, \
#                                                        period_proposed, e_proposed, T_proposed, argPeri_deg_proposed, \
#                                                            planet_K, planet_p, planet_e, planet_argPeri, planet_To, verbose=verbose)
                chi_squared_RV_curr = rv1bodyCalculator(RV_epochs[RVdataSet], RVs_proposed[RVdataSet], RVerrors[RVdataSet], sigma_jitter,\
                                                         inclination_deg_proposed, period_proposed, e_proposed, T_proposed, argPeri_deg_proposed,\
                                                          np.mean(a1s), verbose=False)
                
                chi_squared_RV = chi_squared_RV+chi_squared_RV_curr
                numEpochs_RV_curr = len(RV_epochs[RVdataSet])
                numEpochs_RV = numEpochs_RV+numEpochs_RV_curr


            # convert to a reduced chiSquared
            chi_squared_RV_reduced = (1.0/((1.0*numEpochs_RV)-totalNumModelParams))*chi_squared_RV
        
        
        
        # sum non-reduced chiSquareds
        #chi_squared_total = chi_squared_DI + chi_squared_RV + gauss_chi_squared_total
        #weighting scheme chi squared total
        chi_squared_total = (chi_squared_DI*DIweight) + (chi_squared_RV*RVweight)
        
        # Total Reduced Chi Squared
        totalNumDataPoints = 2.0*numEpochs_DI + 1.0*numEpochs_RV 
        TotalNu = (totalNumDataPoints-totalNumModelParams)
        chi_squared_total_reduced = (1.0/TotalNu)*chi_squared_total      
        
        # calculate the likelihood ratio for the previous(latest) and current (proposed) models/states
        likelihood_ratio=0.0
        try:
            likelihood_ratio = math.exp((chiSquare_latest - chi_squared_total_reduced)/ (2.0*temperature))#in the past this was the non-reduced version, not if reduced is allowed by I like it better
        except:
            blah=1 #$$$$$$$$$
        
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
            
            longAN_deg_latest = longAN_deg_proposed
            e_latest = e_proposed
            T_latest = T_proposed
            period_latest = period_proposed
            inclination_deg_latest = inclination_deg_proposed
            argPeri_deg_latest = argPeri_deg_proposed
            a_total_latest = a_total_proposed
            sine_inc_latest = sine_inc_proposed
            chiSquare_latest = chi_squared_total_reduced #in the past this was the non-reduced version, not if reduced is allowed by I like it better
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
                line = line +'   '+ str(a_total_proposed) # a_total=a2 if Mass1 or Mass2 are non-default (ie !=1)
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
                        ','+str(resultsDict['period_best'])+','+str(resultsDict['inclination_deg_best'])+','+str(resultsDict['argPeri_deg_best'])+\
                        ','+str(resultsDict['a_total_best'])+','+str(resultsDict['RV_origin_vel_0_best'])+','+str(resultsDict['RV_origin_vel_1_best'])+']'
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
        print '\nMCMC: $$$$$$$$$$$$$ FINISHED MCMC-Uniform7 Sim $$$$$$$$$$$$$$$\n'
    
#****** DONE!! ********        
            
            

    
    