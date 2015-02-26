#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import math
import gc
import numpy as np
import os
import pylab
from math import pi
import generalToolbox as genTools

def TAcalculator(t,e, T, period, T_center=0, verbose=False, debug=False):
    """
    $$$ Copy of exact function from generalToolbox as something odd was happening $$$$
    
    :param float t: Epoch [JD]
    :param float e: Eccentricity 
    :param float T: Time of Last Periapsis [JD]
    :param float period: Period of the orbit [years]
    :param float T_center: Time of Center Transit [JD], default 0 indicates to 
        it isn't required as it only applies to radial velocity data.
    :param bool verbose: Print 
    """
    
    ## calculate the Mean Motion
    n = (2*pi)/period
    if verbose:
        print '#'*50
        print 'Mean Motion [rad/yr]= '+str(n)
    
    ## calculate Mean Anomaly
    period_days = period*365.242
    timeDiff_days = (t- T)-int((t-T)/period_days)*period_days 
    if timeDiff_days<0.0:
        timeDiff_days = timeDiff_days+period_days
    phase = 0.0
    if T_center!=0.0:
        phaseDiff_days = (T_center-T)-int((T_center -T)/period_days)*period_days 
        if T>T_center:
            phaseDiff_days = phaseDiff_days+period_days
        phase = phaseDiff_days/period_days
    if verbose:
        print "Unitless phase calculated to be "+str(phase)+", using T_center = "+str(T_center)+" and To = "+str(T)
        print 'timeDiff_days = '+str(timeDiff_days)+', phaseDiff_days = '+str(phaseDiff_days)+', period_days = '+str(period_days)
    
    M = n*(((timeDiff_days)/365.242)+phase)#+(phase*2.0*pi)
    ## Push M value into 0-2pi range ######
    numCirclesBiggerD = abs(M/(2.0*pi))
    numCirclesBiggerI = int(numCirclesBiggerD)
    if numCirclesBiggerI<1:
        numCirclesBiggerI = 1
    if (M<0.0):
        M_out = (numCirclesBiggerI+1)*2.0*pi + M
        if verbose:
            print "M updated from, "+str(M)+" to "+str(M_out)+", numCirclesBiggerI = "+str(numCirclesBiggerI)+", numCirclesBiggerD = "+str(numCirclesBiggerD)
        
    elif M>(2.0*pi):
        M_out = M-numCirclesBiggerI*2.0*pi
        if verbose:
            print "M updated from, "+str(M)+" to "+str(M_out)+", numCirclesBiggerI = "+str(numCirclesBiggerI)+", numCirclesBiggerD = "+str(numCirclesBiggerD)
    else:
        M_out = M
    M = M_out
    
    M_deg = math.degrees(M) # convert resulting M to degrees
    if verbose:
        print 'Mean Anomaly [deg]= ',M_deg

    ### Performing Newton's Method to get the Eccentric Anomaly value ###
    if verbose:
        print '-'*50
    # initial guess (E_last), will be updated in loop.  
    # Anything works, just takes longer if further from real value. => pi
    E_last = 2*pi
    # stored initial value to be updated in loop
    # this value is always very close to the true value and will minimize the number of loops
    # the inital guess used here is from Argyle "Observing and Measuring Visual Double Stars"
    try:
        E_latest = M+e*math.sin(M)+((e**2.0)/(2.0*M))*math.sin(2.0*M)
    except:
        # not sure why, but there were cases of division by zero
        # thus, for these cases, I will resort to my old predicted start value that doesn't suffer from this.
        E_latest = M+e*math.sin(M)
    M_last = M
    # show input value to 
    if debug:
        print "Inputs to Newton's Method are : \nM [rad]= "+str(M)+"\nE_last [rad]= "+\
        str(E_last)+"\nE_latest [rad] = "+str(E_latest)+"\ne = "+str(e)
        
        print "\nStarting to run Newton's while loop."
    
    count = 0 # a counter to stop inf loops in Newton's method below
    while (abs(E_last-E_latest) > (1.0e-10))and(count<50):
        if debug:
            print 'current E [rad]= ', E_latest
        E_last = E_latest
        M_last = E_last - e*math.sin(E_last)
        E_latest = E_last - ((M_last-M)/(1.0-e*math.cos(E_last)))
        count = count+1

    E_latest_deg = math.degrees(E_latest) # convert resulting E to degrees
    if verbose:
        print "The resultant E value is [deg] = ", E_latest_deg
    # check if the resultant value solves the original equation.
    Mnewton = math.degrees(E_latest-e*math.sin(E_latest))
    if abs(M_deg-Mnewton)>(1.0e-5):
        if verbose:
            print "PROBLEM: This resultant E does not satisfy the original equation, Newton's method FAILED !!!"
            print 'M from this E Equals = '+str(Mnewton)
            print 'M original = '+str(M_deg)
            print 'E initial = '+str(math.degrees(M+e*math.sin(M) ))
            print 'e = '+str(e)
    else:
        if debug:
            print "This resultant E solves the original equation, Newton's Method worked :-)"
            print '-'*50
    ### Newton's loop finished! ###
    
    ## calculate True Anomaly from Eccentric Anomaly
    top = (math.cos(E_latest)-e)
    btm = (1.0-e*math.cos(E_latest))
    TA_rad  = math.acos(top/btm) 
    
        
#    if E_latest<0.0:
#        #print 'top is negative'
#        TA_rad = -1.0*TA_rad
#    if False:
#        print "\nTA_rad calc failed!"
#        print "E_latest = ",E_latest
#        print "e = ",e
#    
#    if E_latest>(2.0*pi):
#        # convert E to be inside one orbit (ie. under 2*PI)
#        numCirclesD = E_latest/(2.0*pi)
#        numCirclesI = int(numCirclesD)
#        E_latest_oneCircle  = E_latest-numCirclesI*2.0*pi
#        if verbose: 
#            print "E_latest found to be "+str(numCirclesI)+" times over 2pi, so made E_latest_circle = "+str(math.degrees(E_latest_oneCircle))
#    else:
#        E_latest_oneCircle  = E_latest

    if (E_latest>pi) or (E_latest<0):
        # this is to take care of the silly fact that after E=180deg,
        # the equation above produces TAs that go down from 180, rather than up.
        ## This has been found to be due to the equation eliminating negative
        ## values because of the cosine, thus the new update of E_latest<0.0:
        ## should fix this problem, but the solutions will not go 0->360
        ## but rather 0->180, -180->-0.
        TA_rad_orig = TA_rad
        TA_rad = 2.0*pi - TA_rad
        if verbose:
            print "E_latest found to be over PI, so changed TA from "+str(math.degrees(TA_rad_orig))+" to "+str(math.degrees(TA_rad))
        
    #TA_deg = math.degrees(TA_rad)#$$$$$$$$$$$$$$$$
    #print 'epoch = '+str(t)+', T = '+str(T)+', timeDiff_days = '+str(timeDiff_days)+', numPeriodDiff = '+str(int((t-T)/period_days))+', TA_deg = '+str(TA_deg)#$$$$$$$$$$$$$$$$$$$$$$$$
#    ## Calculate TA in another way    
#    x = ((1.0-e**2.0)**(1.0/2.0))*math.cos(E_latest/2.0)
#    y = ((1.0+e**2.0)**(1.0/2.0))*math.sin(E_latest/2.0)
#    TA_rad2 = 2.0*math.atan2(y, x)
#    #print 'TA_2 = '+str(math.degrees(TA_rad2))
#    
#    print 'TA = ',math.degrees(TA_rad)
#    print 'TA2 = ',math.degrees(TA_rad2)
    
#    if True:
#        if t<2452381:
#            print "epoch = "+str(t)+", To = "+str(T)+", Tc = "+str(T_center)+", e = "+str(e)+", period = "+str(period)+ ", E = "+str(math.degrees(E_latest))+", TA = "+str(math.degrees(TA_rad))
    
    return (n, M_deg, E_latest_deg,TA_rad)

def DIdataToDict(filename):
    """
    A function to extract the Direct Imaging/Astrometry data from a .dat or .txt file
    and output it as a Python dictionary.
    
    Data must be in the columns:
    obsDate[JD]  PA[deg]  PA_error[deg]  SA[arcsec]  SA_error[arcsec]  
    
    :param str filename: File name, including the full path, to the Astrometry 
        data file to convert into a dictionary.
    :return: dictionary with the following members: obsDates, PAs, PA_errors,
        SAs, SA_errors and numEpochs.
    :rtype: dict
    """
    verbose = False
    
    # checking if it has '.txt' on the end, if not add it
    if (filename[-4:]!='.txt' and filename[-4:]!='.dat'):
        print 'Changing output filename from '+filename+' to '+filename+'.dat'
        filename = filename+'.dat'
    
    # First load first file into memory
    file = open(filename, 'r')
    lines = file.readlines()
    INtitle = lines[0]
    colHeaders = lines[1]
    numColumns = len(lines[2].split())
    numLines = len(lines)
    file.close()
    
    dataDict = {}
    obsDates = []
    PAs = []
    PA_errors = []
    SAs = []
    SA_errors = []
    numEpochs = 0
    
    if (numColumns is not 5):
        print 'WARNING: There must be 5 columns of data to use DIdataToDict!'
        print str(numColumns)+" were found."
    else:
        for lineNum in range(2,numLines):
            vals =  lines[lineNum].split()
            if len(vals)==5:
                try:
                    obsDates.append(float(vals[0]))
                    PAs.append(float(vals[1]))
                    PA_errors.append(float(vals[2]))
                    SAs.append(float(vals[3]))
                    SA_errors.append(float(vals[4]))
                    numEpochs+=1
                except:
                    if verbose:
                        print "WARNING: found line that didn't have good data in DIdataToDict"
                        print "line was:"+lines[lineNum]
                        print "continuing to next line."
    
    # all lines loaded int lists, so convert into dict
    dataDict['DI_epochs'] = obsDates
    dataDict['PAs'] = PAs
    dataDict['PA_errors'] = PA_errors
    dataDict['SAs'] = SAs
    dataDict['SA_errors'] = SA_errors
    dataDict['numEpochs'] = numEpochs
    
    return dataDict   
            
def proposalMaxMinsCalc(latestVal, sigma, Min, Max):
    """
    Find max and min values for proposing a new random parameter value when doing MCMC.
    
    Note: This function is not used anymore as the MCMC process is now done in C++.
        
    :param float latestVal: Latest successful value in a MCMC or Simulated Annealing
        chain.
    :param float sigma: The step size to generate the next proposed step from.
    :param float Min: Minimum allowed value for the parameter being proposed.
    :param float Max: Maximum allowed value for the parameter being proposed.
    """
    if (latestVal+sigma)>=Max:
        max = Max
        min = latestVal-sigma
    elif (latestVal-sigma)<=Min:
        max = latestVal+sigma
        min = Min
    else:
        max = latestVal + sigma
        min = latestVal - sigma
    
    return (min, max)

def multiEpochOrbCalc(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors,\
                       epochs, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total=False,\
                       Mass1=1, Mass2=1, verbose=False):
    """
    USES orbitCalculatorSAPA to compare the SA and PA values.
    
    Note: SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors inputs MUST all have
          the same length.
          
    :param SA_arcsec_measured_REALs:      measured Separation Angles ["]
    :type SA_arcsec_measured_REALs:       list of floats
    :param SA_mean_errors:                errors in the measured Separation Angles ["]
    :type SA_mean_errors:                 list of floats
    :param PA_deg_measured_REALs:         measured Position Angles [deg]
    :type PA_deg_measured_REALs:          list of floats
    :param PA_mean_errors:                error in the measured Position Angles [deg]
    :type PA_mean_errors:                 list of floats
    :param epochs:                        epochs of observation/image [julian date]
    :type epochs:                         list of floats
    :param Sys_Dist_PC:                   measured system distance from Earth [PC]
    :type Sys_Dist_PC:                    float
    :param inclination_deg:               inclination [deg]
    :type inclination_deg:                float
    :param longAN_deg:                    Longitude of Ascending Node [deg]
    :type longAN_deg:                     float
    :param e:                             eccentricity of orbits [unitless]
    :type e:                              float
    :param T:                             Last Periapsis Epoch/time [julian date] 
    :type T:                              float
    :param period:                        period of orbits [yrs]
    :type period:                         float
    :param argPeri_deg:                   Argument of Periapsis in orbital plane [deg]
    :type argPeri_deg:                    float
    :param a_total:                       Total semi-major axis of system [AU]
    :type a_total:                        float
    :param Mass1:                         Mass of primary in system [solar masses] (Default=1, means M2 is a planet)
    :type Mass1:                          float
    :param Mass2:                         Mass of the secondary in system [solar masses] (Default=1, means M2 is a planet)
    :type Mass2:                          float
    :param verbose:                       Send prints to screen? [True/False] (Default = False)
    :type verbose:                        Python boolean (True/False). Default = False
    """
    testing = False
    # init output lists
    ns = []
    Ms = []
    Es = []
    thetas = []
    Sep_Dists = []
    SA_arcsec_measured_models = []
    PA_deg_measured_models = []
    a1s = []
    a2s = []

    # initial values for boolean equality of while statement 
    i = 0
    chi_squared_total = 0.0
    
    while i<len(epochs):    

        SA_arcsec_measured_REAL = SA_arcsec_measured_REALs[i]
        SA_mean_error = SA_mean_errors[i]
        PA_deg_measured_REAL = PA_deg_measured_REALs[i]
        PA_mean_error = PA_mean_errors[i]
        t = epochs[i]  
        
        # call orbitCalculator to take input variables and calc orbital elements
        (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU, SA_arcsec_measured_model, PA_deg_measured_model, a1, a2) = \
        orbitCalculatorSAPA(t, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total=a_total,\
                        Mass1=Mass1, Mass2=Mass2, verbose=verbose)
        
        # store output orbital elements in lists for later plotting
        ns.append(n)
        Ms.append(M_deg)
        Es.append(E_latest_deg)
        thetas.append(TA_deg)
        Sep_Dists.append(Sep_Dist_AU)
        SA_arcsec_measured_models.append(SA_arcsec_measured_model)
        PA_deg_measured_models.append(PA_deg_measured_model)
        a1s.append(a1)
        a2s.append(a2)
            
        ## calc both PA and SA kai's 
        SA_chi_squared = genTools.chiSquaredCalc(SA_arcsec_measured_REAL, SA_mean_error, SA_arcsec_measured_model)
        # we must account for the 360deg boundry, using +- 25deg from it as the region of conflict
        if ((PA_deg_measured_model-25.0)<0.0) and (PA_deg_measured_REAL>335.0):
            #ie real is close to 360 and model is just over 360
            PA_deg_measured_model = PA_deg_measured_model +360.0
        if ((PA_deg_measured_REAL-25.0)<0.0) and (PA_deg_measured_model>335.0):
            #ie model is close to 360 and real is just over 360
            PA_deg_measured_model = PA_deg_measured_model -360.0
        PA_chi_squared = genTools.chiSquaredCalc(PA_deg_measured_REAL, PA_mean_error, PA_deg_measured_model)
        # Add them to get the updated total
        chi_squared_curr = PA_chi_squared + SA_chi_squared
        chi_squared_total = chi_squared_total + chi_squared_curr
        
        if testing:
            print str(i)+", chi_squared_curr = "+str(chi_squared_curr)#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        # increment to next epoch to check these input for 
        i = i + 1   
        
        ##### while loop ends here!
        
    return (chi_squared_total, ns, Ms, Es, thetas, Sep_Dists, SA_arcsec_measured_models, PA_deg_measured_models, a1s, a2s) 

def multiEpochOrbCalc3(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors,\
                       epochs, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total=False,\
                       Mass1=1, Mass2=1, verbose=False):
    """
    Uses orbitCalculatorXY to compare the x and y values.
    
    Note: SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors inputs MUST all have
          the same length.
          
    :param SA_arcsec_measured_REALs:      measured Separation Angles ["]
    :type SA_arcsec_measured_REALs:       list of floats
    :param SA_mean_errors:                errors in the measured Separation Angles ["]
    :type SA_mean_errors:                 list of floats
    :param PA_deg_measured_REALs:         measured Position Angles [deg]
    :type PA_deg_measured_REALs:          list of floats
    :param PA_mean_errors:                error in the measured Position Angles [deg]
    :type PA_mean_errors:                 list of floats
    :param epochs:                        epochs of observation/image [julian date]
    :type epochs:                         list of floats
    :param Sys_Dist_PC:                   measured system distance from Earth [PC]
    :type Sys_Dist_PC:                    float
    :param inclination_deg:               inclination [deg]
    :type inclination_deg:                float
    :param longAN_deg:                    Longitude of Ascending Node [deg]
    :type longAN_deg:                     float
    :param e:                             eccentricity of orbits [unitless]
    :type e:                              float
    :param T:                             Last Periapsis Epoch/time [julian date] 
    :type T:                              float
    :param period:                        period of orbits [yrs]
    :type period:                         float
    :param argPeri_deg:                   Argument of Periapsis in orbital plane [deg]
    :type argPeri_deg:                    float
    :param a_total:                       Total semi-major axis of system [AU]
    :type a_total:                        float
    :param Mass1:                         Mass of primary in system [solar masses] (Default=1, means M2 is a planet)
    :type Mass1:                          float
    :param Mass2:                         Mass of the secondary in system [solar masses] (Default=1, means M2 is a planet)
    :type Mass2:                          float
    :param verbose:                       Send prints to screen? [True/False] (Default = False)
    :type verbose:                        Python boolean (True/False). Default = False
    """
    testing = False
    # init output lists
    ns = []
    Ms = []
    Es = []
    thetas = []
    xs = []
    ys = []
    a1s = []
    a2s = []

    # initial values for boolean equality of while statement 
    i = 0
    chi_squared_total = 0.0
    SAPA_chisquared_total = 0.0
    while i<len(epochs):    

        SA_arcsec_measured_REAL = SA_arcsec_measured_REALs[i]
        SA_mean_error = SA_mean_errors[i]
        PA_deg_measured_REAL = PA_deg_measured_REALs[i]
        PA_mean_error = PA_mean_errors[i]
        t = epochs[i]  
        
        # call orbitCalculator to take random variables and calc orbital elements
        #(n, M_deg, E_latest_deg, TA_deg, x, y, a1, a2) = \
        (n, M_deg, E_latest_deg, TA_deg, x, y, a1, a2, SA_arcsec_RP_model, PA_deg_RP_model)=\
        orbitCalculatorXY(t, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total=a_total,\
                        Mass1=Mass1, Mass2=Mass2, verbose=verbose)
        # store output orbital elements in lists for later plotting
        ns.append(n)
        Ms.append(M_deg)
        Es.append(E_latest_deg)
        thetas.append(TA_deg)
        xs.append(x)
        ys.append(y)
        a1s.append(a1)
        a2s.append(a2)
            
        x_real = SA_arcsec_measured_REAL*math.cos(math.radians(PA_deg_measured_REAL))  
        y_real = SA_arcsec_measured_REAL*math.sin(math.radians(PA_deg_measured_REAL))  
        
        # calc error in x,y data
        tempAa = SA_mean_error*np.cos(np.radians(PA_deg_measured_REAL))
        tempA = tempAa*tempAa
        tempBa = SA_arcsec_measured_REAL*np.sin(np.radians(PA_deg_measured_REAL))*np.radians(PA_mean_error)
        tempB = tempBa*tempBa
        x_data_error = np.sqrt(tempA+tempB)

        tempCa =SA_mean_error*np.sin(np.radians(PA_deg_measured_REAL))
        tempC = tempCa*tempCa
        tempDa = SA_arcsec_measured_REAL*np.cos(np.radians(PA_deg_measured_REAL))*np.radians(PA_mean_error)
        tempD = tempDa*tempDa
        y_data_error = np.sqrt(tempC+tempD)
        
        ## calc both PA and SA kai's 
        x_chi_squared = genTools.chiSquaredCalc(x_real, x_data_error, x)
        y_chi_squared = genTools.chiSquaredCalc(y_real, y_data_error, y)
        
        # Add them to get the updated total
        chi_squared_curr = y_chi_squared + x_chi_squared
        chi_squared_total = chi_squared_total + chi_squared_curr
        
        if testing:
            # we must account for the 360deg boundry, using +- 25deg from it as the region of conflict
            if ((PA_deg_RP_model-25.0)<0.0) and (PA_deg_measured_REAL>335.0):
                #ie real is close to 360 and model is just over 360
                orig = PA_deg_RP_model
                PA_deg_RP_model = PA_deg_RP_model +360.0
                print 'changed PA from '+str(orig)+" to "+str(PA_deg_RP_model)
            if ((PA_deg_measured_REAL-25.0)<0.0) and (PA_deg_RP_model>335.0):
                #ie model is close to 360 and real is just over 360
                orig = PA_deg_RP_model
                PA_deg_RP_model = PA_deg_RP_model -360.0
                print 'changed PA from '+str(orig)+" to "+str(PA_deg_RP_model)
            print '************************************************************************************************'
            print str(i)
            print 'SA_arcsec_measured_REAL = '+str(SA_arcsec_measured_REAL)+', PA_deg_measured_REAL = '+str(PA_deg_measured_REAL)
            print 'SA_arcsec_RP_model = '+str(SA_arcsec_RP_model)+", PA_deg_RP_model = "+str(PA_deg_RP_model)
            SA_chiSquared = genTools.chiSquaredCalc(SA_arcsec_measured_REAL, SA_mean_error, SA_arcsec_RP_model)
            PA_chiSquared = genTools.chiSquaredCalc(PA_deg_measured_REAL, PA_mean_error, PA_deg_RP_model)
            SAPA_chisquared_total = SAPA_chisquared_total+SA_chiSquared+PA_chiSquared
            print 'SA_mean_error = '+str(SA_mean_error)+', PA_mean_error = '+str(PA_mean_error) 
            print 'chiSquared from SA = '+str(SA_chiSquared)+', from PA = '+str(PA_chiSquared)+', total = '+str(SA_chiSquared+PA_chiSquared)
            print 'SAPA chiSquared total = '+str(SAPA_chisquared_total)
            print '\nx_real  = '+str(x_real)+', y_real  = '+str(y_real)
            print 'x_model  = '+str(x)+', y_model  = '+str(y)
            print 'x_real_error = '+str(x_data_error)+', y_real_error = '+str(y_data_error)
            print 'x_chi_squared = '+str(x_chi_squared)+', y_chi_squared = '+str(y_chi_squared)+', total = '+str(x_chi_squared +y_chi_squared)
            print 'xy chiSquared total = '+str(chi_squared_total)
            #print '################################################################################################'
        
        # increment to next epoch to check these input for 
        i = i + 1   
        
        ##### while loop ends here!
        
    return (chi_squared_total, ns, Ms, Es, thetas, xs, ys, a1s, a2s)


def orbitCalculator(t, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total,\
                    Mass1=1, Mass2=1, verbose=False):
    """
    NOTE: This is the old version of the orbitCalculator not used anymore.  
    See the orbitCalculatorSAPA and orbitCalculatorXY being used now.
    
    This function will calculate the output orbital parameters from the inputs for a SINGLE epoch.  Use the 
    multiEpochOrbCalc to do this for multiple epochs at the same time.  
    The t and Sys_Dist_PC parameters are generally known/measured values,
    while the rest are typically random numbers.
    This function is designed to work as the calculator inside a Monte Carlo Simulator.
    All the equations are the result of Kepler's Laws and naturally occurring symmetries of binary star systems.
    
    Note: RP stands for Reference Plane, ie the plane that values are measured in from Earth.
          OP stands for Orbital Plane, ie in the plane that the stars are orbiting in.
    
    :param t:                             epoch of observation/image [julian date]
    :type t:                              float
    :param Sys_Dist_PC:                   measured system distance from Earth [PC]
    :type Sys_Dist_PC:                    float
    :param inclination_deg:               inclination [deg]
    :type inclination_deg:                float
    :param longAN_deg:                    Longitude of Ascending Node [deg]
    :type longAN_deg:                     float
    :param e:                             eccentricity of orbits [unitless]
    :type e:                              float
    :param T:                             Last Periapsis Epoch/time [julian date] 
    :type T:                              float
    :param period:                        period of orbits [yrs]
    :type period:                         float
    :param argPeri_deg:                   Argument of Periapsis in orbital plane [deg]
    :type argPeri_deg:                    float
    :param a_total:                       Total semi-major axis of system [AU]
    :type a_total:                        float
    :param Mass1:                         Mass of primary in system [solar masses] (Default=1, means M2 is a planet)
    :type Mass1:                          float
    :param Mass2:                         Mass of the secondary in system [solar masses] (Default=1, means M2 is a planet)
    :type Mass2:                          float
    :param verbose:                       Send prints to screen? [True/False] (Default = False)
    :type verbose:                        Python boolean (True/False). Default = False
    
    :returns: (Mean Motion [rad/yr], Mean Anomaly [deg], Eccentric Anomaly [deg],
        True Anomaly [deg], Separation Distance in orbital plane [AU], 
        measured Position Angle in image [deg], measured Position Angle in image [deg],
        semi-major axis of M1 [AU], semi-major axis of M2 [AU])
    :rtype: list of floats
    """
    if verbose:
        print '\n'+'*'*50
        print 'Starting to calculate Orbital Parameters\n'
        print 'Input variable values: '+'\nCurrent epoch time [JD] = '+str(t)+\
        '\nSystem Distance from Earth [PC] = '+str(Sys_Dist_PC)+\
        '\nInclination [deg] = '+ str(inclination_deg)+\
        '\nLongitude of Ascending Node [deg] = '+str(longAN_deg)+\
        '\nEccentricity = '+str(e)+'\nLast Periapsis Epoch [JD] = '+\
        str(T)+'\nPeriod [yrs] = '+str(period)+\
        '\nArgument of Periapsis in orbital plane [deg] = '+\
        str(argPeri_deg)+'\nverbose = ',str(verbose)
    
    ## get the TAcalculator to find the TA in radians
    (n, M_deg, E_latest_deg,TA_rad) = TAcalculator(t,e, T, period, verbose=verbose)
    TA_deg = math.degrees(TA_rad)
    if verbose:
        print 'True Anomaly [deg] = ', TA_deg
        
    ## calculate the Separation Distance in orbital plane
    Sep_Dist_AU_OP = a_total*((1.0-e*e))/(1.0+e*math.cos(TA_rad))
    if verbose:
        print 'Separation Distance in orbital plane [AU] = ',Sep_Dist_AU_OP
    
    ## calculate measured Separation Distance in reference plane
    inclination_rad = math.radians(inclination_deg) # convert inclination to radians 
    # angle between ascending node and M2 in orbital plane
    ang_LN_to_M2_deg = argPeri_deg + TA_deg
    ang_LN_to_M2_rad = math.radians(ang_LN_to_M2_deg)
    if verbose:
        print "inclination_deg = ",inclination_deg
        print "inclination_rad = ",inclination_rad
        print "ang_LN_to_M2_rad = ",ang_LN_to_M2_rad
    # find horizontal and vertical WRT Line of Nodes (X and Y) components of Sep_Dist corrected for inclination
    # Note: only the 'vertical (Y)' component will be effected by inclination
    # Note2: on this axis, Z would be the direction of 'line of sight' of the observer.
    Sep_Dist_AU_RP_y = Sep_Dist_AU_OP*math.sin(ang_LN_to_M2_rad)*math.cos(inclination_rad)
    Sep_Dist_AU_RP_x = Sep_Dist_AU_OP*math.cos(ang_LN_to_M2_rad)
    # use Pythagorous theorem to find the final Sep Dist in the Ref Plane
    Sep_Dist_AU_RP = math.sqrt(math.pow(Sep_Dist_AU_RP_x, 2.0)+math.pow(Sep_Dist_AU_RP_y,2.0))
    if verbose:
        print 'Sep_Dist_AU_RP_y = ',Sep_Dist_AU_RP_y
        print 'Sep_Dist_AU_RP_x = ',Sep_Dist_AU_RP_x
        print 'Separation Distance in the reference plane [AU] = ',Sep_Dist_AU_RP

    ## calculate corrected angle between Line of Nodes and M2 (ie. AN-focus-M2)
    # the first atan will produce an angle between the X-axis, rather than only CC from positive X-axis.
    # must take all 4 quadrants into account.  
    ang_LN_to_M2_corrPOSNEG = math.degrees(math.atan(Sep_Dist_AU_RP_y/Sep_Dist_AU_RP_x))
    
    if (Sep_Dist_AU_RP_x>=0.0)and (Sep_Dist_AU_RP_y>=0.0):
        # quadrant 1
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG 
    elif (Sep_Dist_AU_RP_x<0.0)and (Sep_Dist_AU_RP_y>=0.0):
        # quadrant 2
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +180.0
    elif (Sep_Dist_AU_RP_x<0.0)and (Sep_Dist_AU_RP_y<0.0):
        # quadrant 3
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +180.0
    elif (Sep_Dist_AU_RP_x>=0.0)and (Sep_Dist_AU_RP_y<0.0):
        # quadrant 4
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +360.0
    else:
        print ' The value of '+str(ang_LN_to_M2_corrPOSNEG)+' is apparently out of ranges'
    
    if verbose:
        print 'uncorrected angle between line of nodes and M2 [deg] = ', ang_LN_to_M2_deg
        print 'Corrected angle between Line of nodes and M2 [deg] = ',ang_LN_to_M2_corr
    
    ## calculate measured Position Angle 
    # it is measured from M1, only if the value is over 360deg do we need to fix it
    totalAngle_RP = longAN_deg + ang_LN_to_M2_corr
    if verbose:
        print "Total angle in the Reference Plane [deg] = ",totalAngle_RP
    if totalAngle_RP>360.0:
        PA_deg_RP_model = totalAngle_RP -360.0
    else:
        PA_deg_RP_model = totalAngle_RP

    if verbose:
        print 'Position Angle in reference plane [deg] = ',PA_deg_RP_model
        
    ## calculate measured Separation Angle (model)
    SA_arcsec_RP_model = Sep_Dist_AU_RP/Sys_Dist_PC # PC->AU and rad->arcsec cancel
    if verbose:
        print 'Separation Angle measured (model) [arcsec] = ',SA_arcsec_RP_model

    if (Mass1 is not 1) and (Mass2 is not 1):
        # this means these two parameters are non-default values and 
        # the system must then be a binary star system, thus calculate the
        # a1 and a2 values of the system.
        # NOTE: Mass2 MUST < Mass1, ie. primary is bigger.
        
        # first find the mass ratio of system
        X = Mass2/Mass1
        
        # find a1 from the a1=X*a2 and a_total=a1+a2
        a2 = a_total/(1.0+X)
        
        # now use a1=X*a2 to get a2
        a1 = a2*X
        
        if (a1+a2)-a_total>1.0e-5:
            print "obCalc-ln1000: a1 = ",a1
            print "obCalc-ln1001:a2 = ",a2
            print "obCalc-ln1002:a_total = ",a_total

        if verbose: 
            print 'The system must be a binary star system with a1 = '+str(a1)+' and a2 = '+str(a2)
    else:
        # default Mass1 and Mass2, thus a planet system and a_total = a2 and a1=0.0
        # this is because we assume most planet's mass << star's mass, so a1<<a2
        a1 = 0.0
        a2 = a_total
    
    if verbose:
        print '\nFinished calculating Orbital Parameters'
        print '*'*50+'\n'
   
    return (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SA_arcsec_RP_model, PA_deg_RP_model, a1, a2)

def orbitCalculatorSAPA(t, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total=False,\
                    Mass1=1, Mass2=1, verbose=False):
    """
    This function will calculate the output orbital parameters from the inputs for a SINGLE epoch.  Use the 
    multiEpochOrbCalc to do this for multiple epochs at the same time.  
    The t and Sys_Dist_PC parameters are generally known/measured values,
    while the rest are typically random numbers.
    This function is designed to work as the calculator inside a Monte Carlo Simulator.
    All the equations are the result of Kepler's Laws and naturally occurring symmetries of binary star systems.
    
    Note: RP stands for Reference Plane, ie the plane that values are measured in from Earth.
          OP stands for Orbital Plane, ie in the plane that the stars are orbiting in.
    
    :param t:                             epoch of observation/image [julian date]
    :type t:                              float
    :param Sys_Dist_PC:                   measured system distance from Earth [PC]
    :type Sys_Dist_PC:                    float
    :param inclination_deg:               inclination [deg]
    :type inclination_deg:                float
    :param longAN_deg:                    Longitude of Ascending Node [deg]
    :type longAN_deg:                     float
    :param e:                             eccentricity of orbits [unitless]
    :type e:                              float
    :param T:                             Last Periapsis Epoch/time [julian date] 
    :type T:                              float
    :param period:                        period of orbits [yrs]
    :type period:                         float
    :param argPeri_deg:                   Argument of Periapsis in orbital plane [deg]
    :type argPeri_deg:                    float
    :param a_total:                       Total semi-major axis of system [AU]
    :type a_total:                        float if it is to be used, boolean False if to be calculated from Kepler's 3rd
    :param Mass1:                         Mass of primary in system [kg] (Default=1, means M2 is a planet)
    :type Mass1:                          float
    :param Mass2:                         Mass of the secondary in system [kg] (Default=1, means M2 is a planet)
    :type Mass2:                          float
    :param verbose:                       Send prints to screen? [True/False] (Default = False)
    :type verbose:                        Python boolean (True/False). Default = False
    
    :returns: (Mean Motion [rad/yr], Mean Anomaly [deg], Eccentric Anomaly [deg],
        True Anomaly [deg], Separation Distance in orbital plane [AU], 
        measured Position Angle in image [deg], measured Position Angle in image [deg],
        semi-major axis of M1 [AU], semi-major axis of M2 [AU])
    :rtype: list of floats
    """
    if verbose:
        print '\n'+'*'*50
        print 'Starting to calculate Orbital Parameters\n'
        print 'Input variable values: '+'\nCurrent epoch time [JD] = '+str(t)+\
        '\nSystem Distance from Earth [PC] = '+str(Sys_Dist_PC)+\
        '\nInclination [deg] = '+ str(inclination_deg)+\
        '\nLongitude of Ascending Node [deg] = '+str(longAN_deg)+\
        '\nEccentricity = '+str(e)+'\nLast Periapsis Epoch [JD] = '+\
        str(T)+'\nPeriod [yrs] = '+str(period)+\
        '\nArgument of Periapsis in orbital plane [deg] = '+\
        str(argPeri_deg)+'\nverbose = ',str(verbose)
    
    ## check if a_total is to be calculated using Kepler's 3rd
    if a_total==False:
        
        # conversion factors and constants
        SecPerYear = 31557600.0
        G = 6.67300e-11
        MperAU = 149598000000.0
        KGperMsun = 1.98892e30
        
        # convert units of years to seconds
        period_seconds = period*SecPerYear
        
        # convert masses from Msun to kilograms
        Mass1_kg = Mass1*KGperMsun
        Mass2_kg = Mass2*KGperMsun
        # apply K3
        a2 = (((G*(Mass1_kg+Mass2_kg)*(period_seconds**2.0))/(4.0*(pi**2.0))))**(1.0/3.0) # Masses must be in [kg]!!!!! 
        
        # find a1 using mass ratio
        a1 = a2*(Mass2_kg/Mass1_kg)
        
        # sum a's and convert into [AU]
        a_total = (a1+a2)/MperAU
        
    ## get the TAcalculator to find the TA in radians
    (n, M_deg, E_latest_deg,TA_rad) = TAcalculator(t,e, T, period, verbose=verbose)
    TA_deg = math.degrees(TA_rad)
    if verbose:
        print 'True Anomaly [deg] = ', TA_deg
        
    ## calculate the Separation Distance in orbital plane
    Sep_Dist_AU_OP = (a_total*(1.0-e**2.0))/(1.0+e*math.cos(TA_rad))
    if verbose:
        print 'Separation Distance in orbital plane [AU] = ',Sep_Dist_AU_OP
    
    ## calculate measured Separation Distance in reference plane
    inclination_rad = math.radians(inclination_deg) # convert inclination to radians 
    # angle between ascending node and M2 in orbital plane
    ang_LN_to_M2_deg = argPeri_deg + TA_deg
    if ang_LN_to_M2_deg>360.0:
        ang_LN_to_M2_deg = ang_LN_to_M2_deg-360.0
    ang_LN_to_M2_rad = math.radians(ang_LN_to_M2_deg)
    if verbose:
        print "inclination_deg = ",inclination_deg
        print "inclination_rad = ",inclination_rad
        print "ang_LN_to_M2_rad = ",ang_LN_to_M2_rad
    # find horizontal and vertical WRT Line of Nodes (X and Y) components of Sep_Dist corrected for inclination
    # Note: only the 'vertical (Y)' component will be effected by inclination
    # Note2: on this axis, Z would be the direction of 'line of sight' of the observer.
    Sep_Dist_AU_RP_y = Sep_Dist_AU_OP*math.sin(ang_LN_to_M2_rad)*math.cos(inclination_rad)
    Sep_Dist_AU_RP_x = Sep_Dist_AU_OP*math.cos(ang_LN_to_M2_rad)
    # use Pythagorous theorem to find the final Sep Dist in the Ref Plane
    Sep_Dist_AU_RP = ((Sep_Dist_AU_RP_x**2.0)+(Sep_Dist_AU_RP_y**2.0))**(1.0/2.0)
    if verbose:
        print 'Sep_Dist_AU_RP_y = ',Sep_Dist_AU_RP_y
        print 'Sep_Dist_AU_RP_x = ',Sep_Dist_AU_RP_x
        print 'Separation Distance in the reference plane [AU] = ',Sep_Dist_AU_RP

    ## calculate corrected angle between Line of Nodes and M2 (ie. AN-focus-M2)
    # the first atan will produce an angle between the X-axis, rather than only CC from positive X-axis.
    # must take all 4 quadrants into account.  
    ang_LN_to_M2_corrPOSNEG = math.degrees(math.atan(Sep_Dist_AU_RP_y/Sep_Dist_AU_RP_x))
    
    if (Sep_Dist_AU_RP_x>=0.0)and (Sep_Dist_AU_RP_y>=0.0):
        # quadrant 1
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG 
    elif (Sep_Dist_AU_RP_x<0.0)and (Sep_Dist_AU_RP_y>=0.0):
        # quadrant 2
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +180.0
    elif (Sep_Dist_AU_RP_x<0.0)and (Sep_Dist_AU_RP_y<0.0):
        # quadrant 3
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +180.0
    elif (Sep_Dist_AU_RP_x>=0.0)and (Sep_Dist_AU_RP_y<0.0):
        # quadrant 4
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +360.0
    else:
        print ' The value of '+str(ang_LN_to_M2_corrPOSNEG)+' is apparently out of ranges'
    
    #print '\nuncorrected angle between line of nodes and M2 [deg] = ', ang_LN_to_M2_deg
    #print 'Corrected angle between Line of nodes and M2 [deg] = ',ang_LN_to_M2_corr
    if verbose:
        print 'uncorrected angle between line of nodes and M2 [deg] = ', ang_LN_to_M2_deg
        print 'Corrected angle between Line of nodes and M2 [deg] = ',ang_LN_to_M2_corr
    
    ## calculate measured Position Angle 
    # it is measured from M1, only if the value is over 360deg do we need to fix it
    totalAngle_RP = longAN_deg + ang_LN_to_M2_corr
    if verbose:
        print "Total angle in the Reference Plane [deg] = ",totalAngle_RP
    if totalAngle_RP>360.0:
        PA_deg_RP_model = totalAngle_RP -360.0
    else:
        PA_deg_RP_model = totalAngle_RP

    if verbose:
        print 'Position Angle in reference plane [deg] = ',PA_deg_RP_model
            
    ## calculate measured Separation Angle (model)
    SA_arcsec_RP_model = Sep_Dist_AU_RP/Sys_Dist_PC # PC->AU and rad->arcsec cancel
    if verbose:
        print 'Separation Angle measured (model) [arcsec] = ',SA_arcsec_RP_model
            
    if (Mass1 is not 1) and (Mass2 is not 1) and (a_total!=False):
        # this means these two parameters are non-default values and 
        # the system must then be a binary star system, thus calculate the
        # a1 and a2 values of the system.
        # NOTE: Mass2 MUST < Mass1, ie. primary is bigger.
        
        # first find the mass ratio of system
        X = Mass2/Mass1
        # find a1 from the a1=X*a2 and a_total=a1+a2
        a2 = a_total/(1.0+X)
        # now use a1=X*a2 to get a2
        a1 = a2*X
        
        if (a1+a2)-a_total>1.0e-5:
            print "obCalc-ln1000: a1 = ",a1
            print "obCalc-ln1001:a2 = ",a2
            print "obCalc-ln1002:a_total = ",a_total

        if verbose: 
            print 'The system must be a binary star system with a1 = '+str(a1)+' and a2 = '+str(a2)
    else:
        # default Mass1 and Mass2, thus a planet system and a_total = a2 and a1=0.0
        # this is because we assume most planet's mass << star's mass, so a1<<a2
        a1 = 0.0
        a2 = a_total
    
    if verbose:
        print '\nFinished calculating Orbital Parameters'
        print '*'*50+'\n'
   
    return (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SA_arcsec_RP_model, PA_deg_RP_model, a1, a2)

def orbitCalculatorXY(t, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total=False,\
                    Mass1=1, Mass2=1, verbose=False):
    """
    This version will output the x and y values predicted for a certain epoch for comparison to measured x and y
    values in the images.  If you want to compare to SA and PA, use orbitCalculatorSAPA.
    
    This function will calculate the output orbital parameters from the inputs for a SINGLE epoch.  Use the 
    multiEpochOrbCalc to do this for multiple epochs at the same time.  
    The t and Sys_Dist_PC parameters are generally known/measured values,
    while the rest are typically random numbers.
    This function is designed to work as the calculator inside a Monte Carlo Simulator.
    All the equations are the result of Kepler's Laws and naturally occurring symmetries of binary star systems.
    
    Note: RP stands for Reference Plane, ie the plane that values are measured in from Earth.
          OP stands for Orbital Plane, ie in the plane that the stars are orbiting in.
    
    :param t:                             epoch of observation/image [julian date]
    :type t:                              float
    :param Sys_Dist_PC:                   measured system distance from Earth [PC]
    :type Sys_Dist_PC:                    float
    :param inclination_deg:               inclination [deg]
    :type inclination_deg:                float
    :param longAN_deg:                    Longitude of Ascending Node [deg]
    :type longAN_deg:                     float
    :param e:                             eccentricity of orbits [unitless]
    :type e:                              float
    :param T:                             Last Periapsis Epoch/time [julian date] 
    :type T:                              float
    :param period:                        period of orbits [yrs]
    :type period:                         float
    :param argPeri_deg:                   Argument of Periapsis in orbital plane [deg]
    :type argPeri_deg:                    float
    :param a_total:                       Total semi-major axis of system [AU]
    :type a_total:                        float if it is to be used, boolean False if to be calculated from Kepler's 3rd
    :param Mass1:                         Mass of primary in system [kg] (Default=1, means M2 is a planet)
    :type Mass1:                          float
    :param Mass2:                         Mass of the secondary in system [kg] (Default=1, means M2 is a planet)
    :type Mass2:                          float
    :param verbose:                       Send prints to screen? [True/False] (Default = False)
    :type verbose:                        Python boolean (True/False). Default = False
    
    :returns: (Mean Motion [rad/yr], Mean Anomaly [deg], Eccentric Anomaly [deg],
        True Anomaly [deg], Separation Distance in orbital plane [AU], 
        measured Position Angle in image [deg], measured Position Angle in image [deg],
        semi-major axis of M1 [AU], semi-major axis of M2 [AU])
    :rtype: list of floats
    """
    testing = False
    

########################################################################################################
###### UP TO THIS LINE, THIS IS THE EXACT SAME AS ORBITCALCULATOR2, SO JUST CALL IT TO SAVE CODE #######
########################################################################################################
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SA_arcsec_RP_model, PA_deg_RP_model, a1, a2) = \
    orbitCalculatorSAPA(t, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total=a_total,\
                    Mass1=1, Mass2=1, verbose=False)    
    x2 = SA_arcsec_RP_model*math.cos(math.radians(PA_deg_RP_model))
    y2 = SA_arcsec_RP_model*math.sin(math.radians(PA_deg_RP_model))
    if testing:
        print '\n################################################################################################'
        print 'SA_arcsec_RP_model = ',SA_arcsec_RP_model
        print 'PA_deg_RP_model = ',PA_deg_RP_model
        print '\nx2 = SA_arcsec_RP_model*math.cos(math.radians(PA_deg_RP_model)) = ',x2
        print 'y2 = SA_arcsec_RP_model*math.sin(math.radians(PA_deg_RP_model)) = ',y2
    
    if verbose:
        print '\nFinished calculating Orbital Parameters'
        print '*'*50+'\n'
   
    #return (n, M_deg, E_latest_deg, TA_deg, x2, y2, a1, a2)
    return (n, M_deg, E_latest_deg, TA_deg, x2, y2, a1, a2, SA_arcsec_RP_model, PA_deg_RP_model)

def orbitCalculatorTH_I(A, B, C, F, G, Mass1_Msun, Mass2_Msun, Sys_Dist_PC):
    """
    A function to calculate the Orbital Elements from the input Theile-Innes 
    constants.
    
    :returns: (a1_arcsec, argPeri_rad, longAN_rad, inclination_rad, period)
    """
    rad_per_arcsec = 206264.806  # Extra? 
    
    # calc semi-major axis of primary
    a1_arcsec = math.sqrt(A**2.0+B**2.0+C**2.0)
    
    ## calc argPeri of primary
    #$$$ NOTE: might be a more efficient way to do this, but couldn't see it 
    #$$$       see it immediately, so using this verbose method for now.
    temp1 = math.atan((B-F)/(A+G))
    # convert angle to match quadrant it should be
    quadrant = 0
    if ((B-F)>=0.0)and(((B-F)/(A+G))>=0.0):
        quadrant = 1
    elif ((B-F)>=0.0)and(((B-F)/(A+G))<=0.0):
        quadrant = 2
    elif ((B-F)<=0.0)and(((B-F)/(A+G))>=0.0):
        quadrant = 3
    elif ((B-F)<=0.0)and(((B-F)/(A+G))<=0.0):
        quadrant = 4
#    print 'temp1 before = ',math.degrees(temp1) #$$$$$$$$$$$$$$$$$$
    if (quadrant==2) or (quadrant==3):    
        temp1 = temp1 + pi
    elif quadrant==4:
        temp1 = temp1 + 2.0*pi
    
#    print 'temp1 after = ',math.degrees(temp1) #$$$$$$$$$$$$$$$$
        
        
    temp2 = math.atan((-B-F)/(A-G))
    # convert angle to match quadrant it should be
    quadrant = 0
    if ((-B-F)>=0.0)and(((-B-F)/(A-G))>=0.0):
        quadrant = 1
    elif ((-B-F)>=0.0)and(((-B-F)/(A-G))<=0.0):
        quadrant = 2
    elif ((-B-F)<=0.0)and(((-B-F)/(A-G))>=0.0):
        quadrant = 3
    elif ((-B-F)<=0.0)and(((-B-F)/(A-G))<=0.0):
        quadrant = 4
        
    if (quadrant==2) or (quadrant==3):    
        temp2 = temp2 + pi
    elif quadrant==4:
        temp2 = temp2 + 2.0*pi
    
    # put two corrected angles together to calc final argPeri of primary
    argPeri_primary_rad = temp1+temp2
    
#    print 'argPeri_primary_rad before = ',math.degrees(argPeri_primary_rad) #$$$$$$$$$$$$$$$$$$$$$$$$
    # convert angle to be inside 2*pi
    while argPeri_primary_rad>(2.0*pi):
        argPeri_primary_rad = argPeri_primary_rad-2.0*pi
    if argPeri_primary_rad<0.0:
        argPeri_primary_rad = argPeri_primary_rad +2.0*pi
#    print 'argPeri_primary_rad after = ',math.degrees(argPeri_primary_rad) #$$$$$$$$$$$$$$$$$$$$$$$$
    
    # add pi radians to convert primary's argPeri to companion's
    argPeri_rad = argPeri_primary_rad + pi
    
    
#    print 'argPeri_rad before = ',math.degrees(argPeri_rad) #$$$$$$$$$$$$$$$$$$$$$$$$
    #convert angle to be inside 2.0*pi
    while argPeri_rad>(2.0*pi):
        argPeri_rad = argPeri_rad-2.0*pi
#    print 'argPeri_rad after = ',math.degrees(argPeri_rad) #$$$$$$$$$$$$$$$$$$$$$$$$
    
    # Calc longitude of Ascending Node
    longAN_rad = temp1-argPeri_primary_rad
    
    if longAN_rad>pi:
        longAN_rad = longAN_rad -pi
    
    # Calc inclination
    temp3 = (A-G)*math.cos(argPeri_primary_rad+longAN_rad)
    temp4 = (A+G)*math.cos(argPeri_primary_rad-longAN_rad)
#    print 'A-G = ',(A-G)+
#    print 'temp3 = ',temp3
#    print 'A+G = ',(A+G)
#    print 'temp4 = ',temp4
    inclination_rad = 2.0*math.atan(math.sqrt(abs(temp3/temp4)))
    
#    # calculate inclination another way to ensure it is the same
#    print 'C = ',C
#    print 'a1_arcsec ',a1_arcsec
#    print 'math.sin(argPeri_primary_rad) ',str(math.sin(argPeri_primary_rad))
#    print 'a1_arcsec*math.sin(argPeri_primary_rad) = ',str(a1_arcsec*math.sin(argPeri_primary_rad))
#    print 'abs(C/(a1_arcsec*math.sin(argPeri_primary_rad)))-1.0 ',str(abs(C/(a1_arcsec*math.sin(argPeri_primary_rad)))-1.0)
#    
#    if abs(C/(a1_arcsec*math.sin(argPeri_primary_rad)))>1.0:
#        inside = abs(C/(a1_arcsec*math.sin(argPeri_primary_rad)))-1.0
#    if abs(C/(a1_arcsec*math.sin(argPeri_primary_rad)))>2.0:
#        inside = abs(C/(a1_arcsec*math.sin(argPeri_primary_rad)))-2.0
#    
#    print 'inside = ',inside
#    inclination_rad2 = math.asin(inside)
#    
#    # check both incliantions are the same
#    if (inclination_rad-inclination_rad2)>0.00001:
#        print 'orbCalculatorTH_I2023: WARNING! two inclinations calculated differ by more than 0.00001 !!!'
#        print 'inclination_rad = ',inclination_rad
#        print 'inclination_rad2 = ',inclination_rad2
    
    # Calc the period from Kepler's 3rd law
    # NOTE: multiplying by the system distance in [PC] converts the ["] distances into [AU] to avoid
    #       raising a value less than 1.0 to the third power which has the opposite effect than over 1.0.
    #       the units cancel in the end so it just makes sure the ratio works in the right way.
    
    # conversion factors and constants
    SecPerYear = 31557600.0
    G = 6.67300e-11
    MperAU = 149598000000.0
    KGperMsun = 1.98892e30
    
    temp5 = 4.0*pi**2.0*(Sys_Dist_PC*a1_arcsec*MperAU)**3.0*Mass1_Msun**3.0
    temp6 = G*KGperMsun*(Mass1_Msun+Mass2_Msun)*Mass2_Msun**3.0
    period = (math.sqrt(temp5/temp6))/SecPerYear
    
    return (a1_arcsec, argPeri_rad, longAN_rad, inclination_rad, period)

def orbitCalculatorTH_I2(A, B, C, F, G, Mass1_Msun, Mass2_Msun, Sys_Dist_PC):
    """
    A version of orbitCalculatorTH_I using ideas for calculating argPeri, longAN and a from 
    'The Binary Stars' by robert G. Aitken.
    
    A function to calculate the Orbital Elements from the input Theile-Innes 
    constants.
    
    :returns: (a1_arcsec, argPeri_rad, longAN_rad, inclination_rad, period)
    """
    
    ## calc argPeri of primary
    #$$$ NOTE: might be a more efficient way to do this, but couldn't see it 
    #$$$       see it immediately, so using this verbose method for now.
    argPeri_plus_longAN = math.atan((B-F)/(A+G))
    # convert angle to match quadrant it should be
    quadrant = 0
    if ((B-F)>=0.0)and(((B-F)/(A+G))>=0.0):
        quadrant = 1
    elif ((B-F)>=0.0)and(((B-F)/(A+G))<=0.0):
        quadrant = 2
    elif ((B-F)<=0.0)and(((B-F)/(A+G))>=0.0):
        quadrant = 3
    elif ((B-F)<=0.0)and(((B-F)/(A+G))<=0.0):
        quadrant = 4
        
    if (quadrant==2) or (quadrant==3):    
        argPeri_plus_longAN = argPeri_plus_longAN + pi
    elif quadrant==4:
        argPeri_plus_longAN = argPeri_plus_longAN + 2.0*pi
    
        
        
    argPeri_minus_longAN = math.atan((-B-F)/(A-G))
    # convert angle to match quadrant it should be
    quadrant = 0
    if ((-B-F)>=0.0)and(((-B-F)/(A-G))>=0.0):
        quadrant = 1
    elif ((-B-F)>=0.0)and(((-B-F)/(A-G))<=0.0):
        quadrant = 2
    elif ((-B-F)<=0.0)and(((-B-F)/(A-G))>=0.0):
        quadrant = 3
    elif ((-B-F)<=0.0)and(((-B-F)/(A-G))<=0.0):
        quadrant = 4
        
    if (quadrant==2) or (quadrant==3):    
        argPeri_minus_longAN = argPeri_minus_longAN + pi
    elif quadrant==4:
        argPeri_minus_longAN = argPeri_minus_longAN + 2.0*pi
    
    # new way to calculate longAN and arePeri
    double_argPeri = argPeri_plus_longAN + argPeri_minus_longAN
    argPeri_primary_rad = double_argPeri/2.0
    double_longAN = argPeri_plus_longAN - argPeri_minus_longAN
    longAN_rad = double_longAN/2.0
    
    # make corrections for if longAN is out of accepted 0-180deg range
    if longAN_rad>pi:
        longAN_rad = longAN_rad-pi
        argPeri_primary_rad = argPeri_primary_rad-pi
    if (longAN_rad<0.0) and (longAN_rad>(-pi)):
        longAN_rad = longAN_rad+pi
        argPeri_primary_rad = argPeri_primary_rad+pi
#    if argPeri_primary_rad>(2.0*pi):
#        argPeri_primary_rad = argPeri_primary_rad-2.0*pi
        
    # add pi radians to convert primary's argPeri to companion's
    argPeri_rad = argPeri_primary_rad 
    
    # Calc inclination
    temp3 = (A-G)*math.cos(argPeri_primary_rad+longAN_rad)
    temp4 = (A+G)*math.cos(argPeri_primary_rad-longAN_rad)

    inclination_rad = 2.0*math.atan(math.sqrt(abs(temp3/temp4)))
        
#    tempA = (A+G)/math.cos(argPeri_primary_rad+longAN_rad)
#    tempB = (A-G)/math.cos(argPeri_primary_rad-longAN_rad)
#    inside = ((2.0*tempA)/(tempA+tempB))-1.0
#    inclination_rad2 = math.acos(inside)
#    
#    print 'inclination_rad = '+str(inclination_rad)+' , or in deg = '+str(math.degrees(inclination_rad))
#    print 'inclination_rad2 = '+str(inclination_rad2)+' , or in deg = '+str(math.degrees(inclination_rad2))+', and inside = '+str(inside)
    
    # calc semi-major axis of apparent ellipse
    a_arcsec_try1 = (B-F)/(2.0*math.sin(longAN_rad+argPeri_primary_rad)*(math.cos(inclination_rad/2.0))**2.0)
#    a_try1 = a_arcsec_try1#*Sys_Dist_PC
#    print 'a_try1 = ',a_try1
    a_arcsec = a_arcsec_try1
#    
#    a_arcsec_try2 = ((A*G-B*F)/(math.cos(inclination_rad)))**(1.0/2.0)
#    a_try2 = a_arcsec_try2#*Sys_Dist_PC
#    print 'a_try2 = ',a_try2
#    
#    a_arcsec_try3 = ((A**2.0+B**2.0+F**2.0+G**2.0)/(1.0+(math.cos(inclination_rad))**2.0))**(1.0/2.0)
#    a_try3 = a_arcsec_try3#*Sys_Dist_PC
#    print 'a_try3 = ',a_try3
    
    # Calc the period from Kepler's 3rd law
    # NOTE: multiplying by the system distance in [PC] converts the ["] distances into [AU] to avoid
    #       raising a value less than 1.0 to the third power which has the opposite effect than over 1.0.
    #       the units cancel in the end so it just makes sure the ratio works in the right way.
    
    # conversion factors and constants
    SecPerYear = 31557600.0
    G = 6.67300e-11
    MperAU = 149598000000.0
    KGperMsun = 1.98892e30
    
    temp5 = 4.0*pi**2.0*(Sys_Dist_PC*a_arcsec*MperAU)**3.0
    temp6 = G*KGperMsun*(Mass1_Msun+Mass2_Msun)
    period = (math.sqrt(temp5/temp6))/SecPerYear
    
    return (a_arcsec, argPeri_rad, longAN_rad, inclination_rad, period)


def orbitCalculatorTH_I_PP(A, B, C, F, G, Mass1_Msun, Mass2_Msun, Sys_Dist_PC):
    """
    Way orbital elements are found from the Thiele-Innes constants in PP1995
    
    A function to calculate the Orbital Elements from the input Theile-Innes 
    constants.
    
    :returns: (a_arcsec, argPeri_rad, longAN_rad, inclination_rad, period)
    """
    k = A**2.0 + B**2.0 + F**2.0 + G**2.0
    m = A*G-B*F
    j = math.sqrt(k**2.0-m**2.0)
    z = math.atan((B-F)/(A+G))
    r = math.atan((B+F)/(G-A))
    
    a_arcsec = math.sqrt(j+k)
    inclination_rad = math.atan(math.sqrt(a_arcsec**4.0+m**2.0)/(m))
    if inclination_rad<0.0:
        inclination_rad = inclination_rad+pi
    argPeri_rad = (z+r)/2.0 +pi
    longAN_rad = (z-r)/2.0+pi
    
    # conversion factors and constants
    SecPerYear = 31557600.0
    G = 6.67300e-11
    MperAU = 149598000000.0
    KGperMsun = 1.98892e30
    
    temp5 = 4.0*pi**2.0*(Sys_Dist_PC*a_arcsec*MperAU)**3.0*Mass1_Msun**3.0
    temp6 = G*KGperMsun*(Mass1_Msun+Mass2_Msun)*Mass2_Msun**3.0
    period = (math.sqrt(temp5/temp6))/SecPerYear
    
    return (a_arcsec, argPeri_rad, longAN_rad, inclination_rad, period)

def multiEpochOrbCalcTH_I3(e, T, A, B, C, F, G, Mass1_Msun, Mass2_Msun, Sys_Dist_PC, \
                                                    SAs, PAs, epochs, SA_errors, PA_errors):
    """
    A function to calculate the fit to the input astrometry to a set of Thiele-Innes constants.
    
    This utilizes orbitCalculatorTH_I2 to calculate the orbital elements from the 
    given Thiele-Innes constants.  It also uses the most recent equations to calculate the 
    errors in X_data and y_data that causes matching results to that of multiEpochOrbCalc3,
    which also uses X and Y.
    
    :returns: (a_arcsec, argPeri_rad, longAN_rad, inclination_rad, period, chiSquared_total)
    """
    testing = False
    
    ## Call the orbitCalculatorTH_I to give us the params that stay stable throughout the orbit
    (a_arcsec, argPeri_rad, longAN_rad, inclination_rad, period) = \
            orbitCalculatorTH_I2(A, B, C, F, G, Mass1_Msun, Mass2_Msun, Sys_Dist_PC)
    
    
    ## Loop through all epochs of data to calc the predicted (x,y) and use to calc chiSquared
    # set up initial value for chiSquared to be updated in loop
    chiSquared_total = 0
    chiSquared_total2 = 0
    if (len(SAs)==len(PAs)) and (len(epochs)==len(SAs)):
        for i in range(0,len(epochs)):
            
            # use TAcalculator to get orbital params for current epoch
           (n, M_deg, E_latest_deg,TA_rad) = TAcalculator(epochs[i], e, T, period, verbose=False)
           
           # Calc Normalized rectangular coordinates
           X = math.cos(math.radians(E_latest_deg))-e
           Y = math.sqrt(1.0-e**2.0)*math.sin(math.radians(E_latest_deg))
           
           # Calc x,y values on same coord system as plane of sky (same as data)
           x_model = A*X+F*Y
           y_model = B*X+G*Y
           
           
           # Calc x,y values from SA and PA of data
           x_data = SAs[i]*math.cos(math.radians(PAs[i]))
           y_data = SAs[i]*math.sin(math.radians(PAs[i]))
           
           if testing:
               print '\n#########################################################'
               print 'SAs[i] = ',SAs[i]
               print 'PAs[i] = ',PAs[i]
               print '\nx_model = A*X+F*Y = ',x_model
               print 'y_model = B*X+G*Y = ',y_model
               print 'x_data = SAs[i]*math.cos(math.radians(PAs[i])) = ',x_data
               print 'y_data = SAs[i]*math.sin(math.radians(PAs[i])) = ',y_data
               print '#########################################################'
           # calc error in x,y data
           x_data_error = x_data*((SA_errors[i]/SAs[i])+abs(math.radians(PA_errors[i])*math.tan(math.radians(PAs[i]))))
           y_data_error = y_data*((SA_errors[i]/SAs[i])+abs(math.radians(PA_errors[i])/math.tan(math.radians(PAs[i]))))
           
           # calc chiSquared for both x and y
           chiSquared_x = genTools.chiSquaredCalc(x_data, x_data_error, x_model)
           chiSquared_y = genTools.chiSquaredCalc(y_data, y_data_error, y_model)
           
           # add both of those to running chiSquared total
           chiSquared_total = chiSquared_total+chiSquared_x+chiSquared_y
           
           # calc error in x,y data another way
           x_data_error2 = x_data*((SA_errors[i]/SAs[i])+((math.cos( math.radians(PAs[i]+PA_errors[i]))-math.cos( math.radians(PAs[i]) ))/math.cos( math.radians(PAs[i]) ) ))
           y_data_error2 = y_data*((SA_errors[i]/SAs[i])+((math.sin( math.radians(PAs[i]+PA_errors[i]))-math.sin( math.radians(PAs[i]) ))/math.sin( math.radians(PAs[i]) ) ))
           
#           print 'x_data_error = ',x_data_error
#           print 'y_data_error = ',y_data_error
#           print 'x_data_error2 = ',x_data_error2
#           print 'y_data_error2 = ',y_data_error2
#           
#           # calc chiSquared for both x and y
#           chiSquared_x2 = genTools.chiSquaredCalc(x_data, x_data_error2, x_model)
#           chiSquared_y2 = genTools.chiSquaredCalc(y_data, y_data_error2, y_model)
           
           # add both of those to running chiSquared total
           chiSquared_total2 = chiSquared_total2+chiSquared_x+chiSquared_y
           
    else:
        print 'multiOrbcalcTH_I3: WARNING! Number of epochs, SAs and/or PAs data does not match'
#    
#    print 'chiSquared_total = ',chiSquared_total
#    print 'chiSquared_total2 = ',chiSquared_total2
    
    return (a_arcsec, argPeri_rad, longAN_rad, inclination_rad, period, chiSquared_total)

def ABCFG_MaxMins2(aMax, aMin, argPeri_radMax, argPeri_radMin, longAN_radMax, longAN_radMin, \
                                                                inclination_radMax, inclination_radMin):
    """
    A function to calculate the maximum and minimum allowed values for the 
    Thiele-Innes constants from the max/min values of the Orbital Elements.
    
    Note: old version in the sandbox had some issues that were fixed in this newer version.
    """
    # find values for all 8 permutations of variable combinations
    (A1,B1,C1,F1,G1) = ABCFG_values(aMin,argPeri_radMax, longAN_radMax, inclination_radMax)
    (A2,B2,C2,F2,G2) = ABCFG_values(aMin,argPeri_radMax, longAN_radMax, inclination_radMin)
    (A3,B3,C3,F3,G3) = ABCFG_values(aMin,argPeri_radMin, longAN_radMax, inclination_radMax)
    (A4,B4,C4,F4,G4) = ABCFG_values(aMin,argPeri_radMin, longAN_radMax, inclination_radMin)
    (A5,B5,C5,F5,G5) = ABCFG_values(aMin,argPeri_radMax, longAN_radMin, inclination_radMax)
    (A6,B6,C6,F6,G6) = ABCFG_values(aMin,argPeri_radMax, longAN_radMin, inclination_radMin)
    (A7,B7,C7,F7,G7) = ABCFG_values(aMin,argPeri_radMin, longAN_radMin, inclination_radMin)
    (A8,B8,C8,F8,G8) = ABCFG_values(aMin,argPeri_radMin, longAN_radMin, inclination_radMax)
    
    As = [A1,A2,A3,A4,A5,A6,A7,A8]
    Bs = [B1,B2,B3,B4,B5,B6,B7,B8]
    Cs = [C1,C2,C3,C4,C5,C6,C7,C8]
    Fs = [F1,F2,F3,F4,F5,F6,F7,F8]
    Gs = [G1,G2,G3,G4,G5,G6,G7,G8]
    
    Amax = np.max(As)*(aMax/aMin)
    Amin = np.min(As)
    Bmax = np.max(Bs)*(aMax/aMin)
    Bmin = np.min(Bs)
    Cmax = np.max(Cs)*(aMax/aMin)
    Cmin = np.min(Cs)
    Fmax = np.max(Fs)*(aMax/aMin)
    Fmin = np.min(Fs)
    Gmax = np.max(Gs)*(aMax/aMin)
    Gmin = np.min(Gs)
    # make each max and min the absolute max and min 
    # which depends on if they are positive or negative
    if Amax>0.0:
         Amax = np.max(As)*(aMax/aMin)
    else:
        Amax = np.max(As)
    if Amin<0.0:
        Amin = np.min(As)*(aMax/aMin)
    else:
         Amin = np.min(As)
    
    if Bmax>0.0:
         Bmax = np.max(Bs)*(aMax/aMin)
    else:
        Bmax = np.max(Bs)
    if Bmin<0.0:
        Bmin = np.min(Bs)*(aMax/aMin)
    else:
         Bmin = np.min(Bs)
         
    if Cmax>0.0:
         Cmax = np.max(Cs)*(aMax/aMin)
    else:
        Cmax = np.max(Cs)
    if Cmin<0.0:
        Cmin = np.min(Cs)*(aMax/aMin)
    else:
         Cmin = np.min(Cs)
         
    if Fmax>0.0:
         Fmax = np.max(Fs)*(aMax/aMin)
    else:
        Fmax = np.max(Fs)
    if Fmin<0.0:
        Fmin = np.min(Fs)*(aMax/aMin)
    else:
         Fmin = np.min(Fs)
     
    if Gmax>0.0:
         Gmax = np.max(Gs)*(aMax/aMin)
    else:
        Gmax = np.max(Gs)
    if Gmin<0.0:
        Gmin = np.min(Gs)*(aMax/aMin)
    else:
         Gmin = np.min(Gs)
         
#    Amax = -a1Max
#    Amin = -a1Max
#    Bmax = a1Max
#    Bmin = -a1Max
#    Cmax = a1Max
#    Cmin = -a1Max
#    Fmax = a1Max
#    Fmin = -a1Max
#    Gmax = a1Max
#    Gmin = -a1Max
    return (Amax,Amin, Bmax,Bmin,  Cmax,Cmin, Fmax,Fmin, Gmax,Gmin)
   
   
def ABCFG_values(a, argPeri_rad, longAN_rad, inclination_rad):
    """
    Calculates the Thiele-Innes constants for the given Orbital Elements.
    
    :returns: (A,B,C,F,G) 
    """
    if False:
        print '\nin ABCFG_values'
        print 'a1',a
        print 'argPeri_rad',argPeri_rad
        print 'longAN_rad',longAN_rad
        print 'inclination_rad',inclination_rad
    
    A = a*(math.cos(longAN_rad)*math.cos(argPeri_rad)-\
                   math.sin(longAN_rad)*math.sin(argPeri_rad)*math.cos(inclination_rad))

    B = a*(math.sin(longAN_rad)*math.cos(argPeri_rad)+\
                   math.cos(longAN_rad)*math.sin(argPeri_rad)*math.cos(inclination_rad))
    
    F = a*(-math.cos(longAN_rad)*math.sin(argPeri_rad)-\
                   math.sin(longAN_rad)*math.cos(argPeri_rad)*math.cos(inclination_rad))
    
    G = a*(-math.sin(longAN_rad)*math.sin(argPeri_rad)+\
                   math.cos(longAN_rad)*math.cos(argPeri_rad)*math.cos(inclination_rad))
    
    C = a*math.sin(argPeri_rad)*math.sin(inclination_rad)
    
    return (A,B,C,F,G) 

def PASAcalculator(period, t, T, e, inclination_deg, longAN_deg, argPeri_deg, Sys_Dist_PC, a2, a1=0, verbose=True):
    '''
    This function is to predict/calculate the expected Position Angle and Separation Angle for a given set of 
    of orbital parameters and epoch of observation.
    
    NOTE: 
    This function is capable of working for 1 or 2 body systems.  In order to use it for 2 body systems, just
    set the values of a1 % a2.  For 1 body systems, ie, only one of the bodies is moving and the primary is 
    simply stationary at one of the foci, then leave a1=0 and use a2 for the moving body's semi-major axis value.
    
    Inputs:
    
    :param period:              period of orbits [yrs]
    :type period:               float
    :param t:                   epoch of observation/image [julian date]
    :type t:                    float
    :param T:                   Last Periapsis Epoch/time [julian date] 
    :type T:                    float
    :param e:                   eccentricity of orbit
    :type e:                    float
    :param inclination_deg:     orbit's inclination [degrees]
    :type inclination_deg:      float
    :param longAN_deg:          Longitude of the Ascending Node [degrees]
    :type longAN_deg:           float
    :param argPeri_deg:         Argument of periapsis [degrees]
    :type argPeri_deg:          float
    :param Sys_Dist_PC:         Distance to the system from Earth [parsec]
    :type Sys_Dist_PC:          float
    :param a2:                  Semi-major axis of body 2's orbit [AU]
    :type a2:                   float
    :param a1:                  Semi-major axis of body 1's orbit [AU]
    :type a1:                   float, default=0 indicates body 2 is only moving body in system
    :param verbose:             Show progress prints to screen?
    :type verbose:              python boolean (True/False), default=True
    
    Outputs:
    :param PA_deg_RP_model:     Position Angle [degrees]
    :type PA_deg_RP_model:      float
    :param SA_arcsec_RP_model:  Separation Angle [arcsec]
    :type SA_arcsec_RP_model:   float
    '''

    if verbose:
        print '\n'+'*'*50
        print 'Input variable values: '+\
        '\nPeriod [yrs] = '+str(period)+\
        '\nCurrent epoch time [julian date] = '+str(t)+\
        '\nLast Periapsis Epoch [julian date] = '+str(T)+\
        '\nEccentricity = '+str(e)+\
        '\nInclination [deg] = '+ str(inclination_deg)+\
        '\nLongitude of Ascending Node [deg] = '+str(longAN_deg)+\
        '\nArgument of Periapsis in orbital plane [deg] = '+str(argPeri_deg)+\
        '\nSystem Distance from Earth [PC] = '+str(Sys_Dist_PC)+\
        '\nSemi-major of primary [AU] = '+str(a1)+\
        '\nSemi-major of secondary [AU] = '+str(a2)
        print '\n**Starting to calculate PA and SA from Orbital Parameters**\n'
    
    ## calculate the Mean Motion
    n = (2*pi)/period
    if verbose:
        print 'Mean Motion [rad/yr]= '+str(n)

    ## calculate Mean Anomaly
    M = n*((t-T)/365.0)
    M_deg = math.degrees(M) # convert resulting M to degrees
    if verbose:
        print 'Mean Anomaly [deg]= ',M_deg

    ### Performing Newton's Method to get the Eccentric Anomaly value ###
    if verbose:
        print '-'*50    
    
    # initial guess (E_last), will be updated in loop.  
    # Anything works, just takes longer if further from real value. => pi
    E_last = 2*pi
    # stored initial value to be updated in loop
    # this value is always very close to the true value and will minimize the number of loops
    E_latest = M+e*math.sin(M) 
    
    # show input value to 
    if verbose:
        print "Inputs to Newton's Method are : \nM [rad]= "+str(M)+"\nE_last [rad]= "+\
        str(E_last)+"\nE_latest [rad] = "+str(E_latest)+"\ne = "+str(e)
        
        print "\nStarting to run Newton's while loop."
    
    count = 0 # a counter to stop inf loops in Newtons method below
    while (abs(E_last-E_latest) > (1.0e-10))and(count<100):
        if verbose:
            print 'current E [rad]= ', E_latest
        E_last = E_latest
        #E_latest = E_last - ((M-E_last+e*math.sin(E_last))/(e*math.cos(E_last)-1.0))\
        E_latest = E_last - ((E_last-M-e*math.sin(E_last))/(-e*math.cos(E_last)+1.0))
        count = count+1

    E_latest_deg = math.degrees(E_latest) # convert resulting E to degrees
    if verbose:
        print "The resultant E value is [deg] = ", E_latest_deg
    # check if the resultant value solves the original equation.
    Mnewton = math.degrees(E_latest-e*math.sin(E_latest))
    if abs(M_deg-Mnewton)>(1.0e-5):
        if verbose:
            print "PROBLEM: This resultant E does not satisfy the original equation, Newton's method FAILED !!!"
            print 'M from this E Equals = '+str(Mnewton)
            print 'M original = '+str(M_deg)
            print 'E initial = '+str(math.degrees(M+e*math.sin(M) ))
            print 'e = '+str(e)
    else:
        if verbose:
            print "This resultant E solves the original equation, Newton's Method worked :-)"
            print '-'*50
    ### Newton's loop finished! ###
    
    ## calculate True Anomaly from Eccentric Anomaly
    TA_rad  = math.acos((math.cos(E_latest)-e)/(1-e*math.cos(E_latest))) 
    if E_latest>pi:
        # this is to take care of the silly fact that after E=180deg,
        # the equation above produces TAs that go down from 180, rather than up.
        TA_rad = 2.0*pi - TA_rad
    TA_deg = math.degrees(TA_rad)
    if verbose:
        print 'True Anomaly [deg] = ', TA_deg
        
    ## calculate the Separation Distance in orbital plane
    a_total = a1+a2
    Sep_Dist_AU_OP = a_total*((1-e*e))/(1+e*math.cos(TA_rad))
    if verbose: #$$
        print '(1-e*e) = '+str((1-e*e))#$$$
        print '(1+e*cos(TA)) = '+str((1+e*math.cos(TA_rad)))#$$
    if verbose:
        print 'Separation Distance in orbital plane [AU] = ',Sep_Dist_AU_OP
    
    ## calculate measured Separation Distance in reference plane
    inclination_rad = math.radians(inclination_deg) # convert inclination to radians 
    # angle between ascending node and M2 in orbital plane
    ang_LN_to_M2_deg = argPeri_deg + TA_deg
    # find horizontal and vertical WRT Line of Nodes (X and Y) components of Sep_Dist corrected for inclination
    # Note: only the 'vertical (Y)' component will be effected by inclination
    
    Sep_Dist_AU_RP_y = Sep_Dist_AU_OP*math.sin(math.radians(ang_LN_to_M2_deg))*math.cos(inclination_rad)
    Sep_Dist_AU_RP_x = Sep_Dist_AU_OP*math.cos(math.radians(ang_LN_to_M2_deg))
    
#    if verbose:##$$$
#        print 'sin(TA+argPeri) = '+str(math.sin(math.radians(ang_LN_to_M2_deg)))##$$$
#        print 'cos(TA+argPeri) = '+str(math.cos(math.radians(ang_LN_to_M2_deg)))#$$$
#        print 'cos(i) = '+str(math.cos(inclination_rad))#$$$
#        print 'sep_dist_au_op_y = '+str(Sep_Dist_AU_OP*math.sin(math.radians(ang_LN_to_M2_deg)))
#        print 'sep_dist_au_rp_y = '+str(Sep_Dist_AU_RP_y)#$$$
#        print 'sep_dist_au_rp_x = op_x = '+str(Sep_Dist_AU_RP_x)#$$$
        
    # use Pythagorous theorem to find the final Sep Dist in the Ref Plane
    Sep_Dist_AU_RP = math.sqrt(math.pow(Sep_Dist_AU_RP_x, 2.0)+math.pow(Sep_Dist_AU_RP_y,2.0))
    if verbose:
        print 'Separation Distance in the reference plane [AU] = ',Sep_Dist_AU_RP
    
    # calculate corrected angle between Line of Nodes and M2 (ie. AN-focus-M2)
    # the first atan will produce a negative if x or is negative, so must correct that below
    # must take all 4 quadrants into account.
    ang_LN_to_M2_corrPOSNEG = math.degrees(math.atan(Sep_Dist_AU_RP_y/Sep_Dist_AU_RP_x))
    if verbose:
        print ' angle LN to M2 POSNEG [deg] = '+str(ang_LN_to_M2_corrPOSNEG)
    
    if (Sep_Dist_AU_RP_x>=0.0)and (Sep_Dist_AU_RP_y>=0.0):
        # quadrant 1
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG 
    elif (Sep_Dist_AU_RP_x<0.0)and (Sep_Dist_AU_RP_y>=0.0):
        # quadrant 2
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +180.0
    elif (Sep_Dist_AU_RP_x<0.0)and (Sep_Dist_AU_RP_y<0.0):
        # quadrant 3
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +180.0
    elif (Sep_Dist_AU_RP_x>=0.0)and (Sep_Dist_AU_RP_y<0.0):
        # quadrant 4
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +360.0
    else:
        print ' The value of '+str(ang_LN_to_M2_corrPOSNEG)+' is apparently out of ranges'
        
    if verbose:
        print 'uncorrected angle between line of nodes and M2 [deg] = ', ang_LN_to_M2_deg
        print 'Corrected angle between Line of nodes and M2 [deg] = ',ang_LN_to_M2_corr
    
    ## calculate measured Position Angle 
    # it is measured from M1, only if the value is over 360deg do we need to fix it
    totalAngle_RP = longAN_deg + ang_LN_to_M2_corr
    if totalAngle_RP>360.0:
        PA_deg_RP_model = totalAngle_RP -360.0
    else:
        PA_deg_RP_model = totalAngle_RP
#    if Sep_Dist_AU_RP_x>=0.0:
#        PA_deg_RP_model = longAN_deg + 180.0 + ang_LN_to_M2_corr
#        #print '** corrected M2 found in 2nd or 3rd quadrant' #$$$$$$$$
#    elif Sep_Dist_AU_RP_x<0.0:
#        PA_deg_RP_model = longAN_deg + ang_LN_to_M2_corr
#        #print '** corrected M2 found in 1st or 4th quadrant' #$$$$$$$$$$$
    if verbose:
        print 'Position Angle in reference plane [deg] = ',PA_deg_RP_model
        
    ## calculate measured Separation Angle (model)
    SA_arcsec_RP_model = Sep_Dist_AU_RP/Sys_Dist_PC # PC->AU and rad->arcsec cancel
    if verbose:
        print 'Separation Angle measured (model) [arcsec] = ',SA_arcsec_RP_model
    
    return(PA_deg_RP_model, SA_arcsec_RP_model)

def RADECtoPASA(RA, RA_error, DEC, DEC_error):
    """
    Will calculate the Separation and Position Angles for a given RA and DEC, including their errors.
    
    :returns: (PA,PA_error,SA,SA_error)
    """
    
    PA = math.degrees(math.atan(RA/DEC))
    
    # correction factors for flaws of atan
    #NOTE: this issue could be cured with proper use of np.arctan2(x1,x2), where it would normally be atan(x1/x2).
    if (RA<0.0)and(DEC>0.0):
        PA = 360.0+PA
    elif DEC<0.0:
        PA = 180.0+PA
    
    SA = math.sqrt(RA**2.0 + DEC**2.0)
    
    top = (RA/DEC)*math.sqrt((RA_error/RA)**2.0 + (DEC_error/DEC)**2.0)
    btm = 1.0+(RA/DEC)**2.0
    PA_error = abs(math.degrees(top/btm))
    
    top = SA*(abs(RA*RA_error)+abs(DEC*DEC_error))
    btm = RA**2.0+DEC**2.0
    SA_error = abs(top/btm)
    
    return (PA,PA_error,SA,SA_error)

def PASAtoRADEC(PA,PA_error,SA,SA_error):
    """
    Convert provided Position Angle and Separation Angle, and their errors, into 
    RA and DEC with errors.  These are the same equations for calculating 
    x and y in the Thiele-Innes orbit fitting.  Remember that x and y are 
    reversed in that fitting approach due to how Thiele defined the coord 
    system when deriving the equations used.
    
    NOTE: this can also be used to calculate X and Y used in Thiele-Innes
          With RA=y and DEC=x.  
    
    :returns: (RA, RA_error, DEC, DEC_error)
    """
    DEC = SA*math.cos(math.radians(PA))
    RA = SA*math.sin(math.radians(PA))
    
    tempA = (SA_error*math.cos(math.radians(PA)))**2.0
    tempB = (SA*math.sin(math.radians(PA))*math.radians(PA_error))**2.0
    DEC_error = math.sqrt(tempA+tempB)
    
    tempC = (SA_error*math.sin(math.radians(PA)))**2.0
    tempD = (SA*math.cos(math.radians(PA))*math.radians(PA_error))**2.0
    RA_error = math.sqrt(tempC+tempD)
    
    return (RA, RA_error, DEC, DEC_error)


    
    
