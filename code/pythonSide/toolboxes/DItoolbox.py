#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import math
import gc
import numpy as np
import os
import pylab
from math import pi
import generalToolbox as genTools



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

def multiEpochOrbCalc(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors,\
                       epochs, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total=False,\
                       Mass1=1, Mass2=1, verbose=False, useTHI=True):
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
    :param useTHI:                        Use Thiele-Innes approach? [True/False] (Default = True)
    :type useTHI:                         Python boolean (True/False). Default = False
    
    :returns: (chi squared total, Mean Motions [rad/yr], Mean Anomalys [deg], Eccentric Anomalsy [deg],
        True Anomalys [deg], Separation Distances in orbital plane [AU], 
        measured Position Angles in image [deg], measured Position Angles in image [deg],
        measured xs (North) value in image ["], measured ys (East) in the image ["],
        semi-major axis' of M1 [AU], semi-major axis' of M2 [AU])
    :rtype: list of lists of floats
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
    xs = []
    ys = []
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
        if useTHI:
            (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU, SA_arcsec_measured_model, PA_deg_measured_model, x_model, y_model, a1, a2) = \
            orbitCalculatorTH_I(t, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total, Mass1, Mass2, verbose=verbose)
        else:
            (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU, SA_arcsec_measured_model, PA_deg_measured_model, x_model, y_model, a1, a2) = \
            orbitCalculatorFullEQs(t, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total=a_total,\
                        Mass1=Mass1, Mass2=Mass2, verbose=verbose)
        # store output orbital elements in lists for later plotting
        ns.append(n)
        Ms.append(M_deg)
        Es.append(E_latest_deg)
        thetas.append(TA_deg)
        Sep_Dists.append(Sep_Dist_AU)
        SA_arcsec_measured_models.append(SA_arcsec_measured_model)
        PA_deg_measured_models.append(PA_deg_measured_model)
        xs.append(x_model)
        ys.append(y_model)
        a1s.append(a1)
        a2s.append(a2)
            
        ## calc both PA and SA chi squareds
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
        
    return (chi_squared_total, ns, Ms, Es, thetas, Sep_Dists, SA_arcsec_measured_models, PA_deg_measured_models, xs, ys, a1s, a2s) 


def orbitCalculatorFullEQs(t, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total,\
                    Mass1=1, Mass2=1, verbose=False):
    """
    NOTE: This is the old version that calculates the predicted locations for a companion star/planet
    that uses the full equations approach that Kyle Mede personally developed from first principles.  
    While both have been tested and produce the same outputs, orbitCalculatorTH_I is the standard as it is simpler.
    To use this formulation instead fo the Thiele-Innes approach, in the multiEpochOrbCalc, just set the useTHI boolean 
    to False.
    
    NOTE2: Keep in mind that the 'x' and 'y' returned are the North and East values respectively.  Very important.
    
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
        measured x (North) value in image ["], measured y (East) in the image ["],
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
    (n, M_deg, E_latest_deg,TA_rad) = genTools.TAcalculator(t,e, T, period, verbose=verbose)
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
     
    x_model = SA_arcsec_RP_model*math.cos(math.radians(PA_deg_RP_model))
    y_model = SA_arcsec_RP_model*math.sin(math.radians(PA_deg_RP_model))

    if (Mass1 is not 1) and (Mass2 is not 1):
        # this means these two parameters are non-default values and 
        # the system must then be a binary star system, thus calculate the
        # a1 and a2 values of the system.
        # NOTE: Mass2 MUST < Mass1, ie. primary is bigger.
        (a_total, a1, a2, period) = genTools.semiMajorConverter(Mass1, Mass2, a_total,0,0,period)
    else:
        # default Mass1 and Mass2, thus a planet system and a_total = a2 and a1=0.0
        # this is because we assume most planet's mass << star's mass, so a1<<a2
        a1 = 0.0
        a2 = a_total
    
    if verbose:
        print '\nFinished calculating Orbital Parameters'
        print '*'*50+'\n'
   
    return (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SA_arcsec_RP_model, PA_deg_RP_model, x_model, y_model, a1, a2)    
   
def orbitCalculatorTH_I(t, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total,Mass1=1, Mass2=1, verbose=False):
    """
    NOTE: this is the standard function for calculating the predicted location of a companion star/planet.
    The full equations approach developed from first principles by Kyle Mede is also availble, and has
    been tested to verify both approaches give the same results.  This version is based on the Thiele-Innes 
    equations and is simpler, thus prefered.
    
    NOTE2: Keep in mind that the 'x' and 'y' returned are the North and East values respectively.  Very important.
    
    This function will calculate the output orbital parameters from the inputs for a SINGLE epoch.  Use the 
    multiEpochOrbCalc to do this for multiple epochs at the same time.  
    The t and Sys_Dist_PC parameters are generally known/measured values,
    while the rest are typically random numbers.
    This function is designed to work as the calculator inside a Monte Carlo Simulator.
    It is based on the Thiele-Innes equations.
    
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
        measured x (North) value in image ["], measured y (East) in the image ["],
        semi-major axis of M1 [AU], semi-major axis of M2 [AU])
    :rtype: list of floats
    """
    a_arcsec = a_total/Sys_Dist_PC
    longAN_rad = math.radians(longAN_deg)
    argPeri_rad = math.radians(argPeri_deg)
    inclination_rad = math.radians(inclination_deg)
    A = a_arcsec*(math.cos(longAN_rad)*math.cos(argPeri_rad)-\
                   math.sin(longAN_rad)*math.sin(argPeri_rad)*math.cos(inclination_rad))

    B = a_arcsec*(math.sin(longAN_rad)*math.cos(argPeri_rad)+\
                   math.cos(longAN_rad)*math.sin(argPeri_rad)*math.cos(inclination_rad))
    
    F = a_arcsec*(-math.cos(longAN_rad)*math.sin(argPeri_rad)-\
                   math.sin(longAN_rad)*math.cos(argPeri_rad)*math.cos(inclination_rad))
    
    G = a_arcsec*(-math.sin(longAN_rad)*math.sin(argPeri_rad)+\
                   math.cos(longAN_rad)*math.cos(argPeri_rad)*math.cos(inclination_rad))
    ## get the TAcalculator to find the TA in radians
    (n, M_deg, E_latest_deg,TA_rad) = genTools.TAcalculator(t, e, T, period, T_center=0, verbose=verbose, debug=False)
    # Calc Normalized rectangular coordinates
    X = math.cos(math.radians(E_latest_deg))-e
    Y = math.sqrt(1.0-e**2.0)*math.sin(math.radians(E_latest_deg))
    # Calc x,y values on same coord system as plane of sky (same as data)
    x = A*X+F*Y
    y = B*X+G*Y
    
    #try to calc the individual a1 and a2 from a_total
    if (Mass1 is not 1) and (Mass2 is not 1):
        # this means these two parameters are non-default values and 
        # the system must then be a binary star system, thus calculate the
        # a1 and a2 values of the system.
        # NOTE: Mass2 MUST < Mass1, ie. primary is bigger.
        (a_total, a1, a2, period) = genTools.semiMajorConverter(Mass1, Mass2, a_total,0,0,period)
    else:
        # default Mass1 and Mass2, thus a planet system and a_total = a2 and a1=0.0
        # this is because we assume most planet's mass << star's mass, so a1<<a2
        a1 = 0.0
        a2 = a_total
    ## calculate the Separation Distance in orbital plane
    Sep_Dist_AU_OP = a_total*((1.0-e*e))/(1.0+e*math.cos(TA_rad))
    #ALSO calculate the values in SA and PA incase the user needs them instead
    (PA,PA_error,SA,SA_error) = ENtoPASA(y, 0, x, 0)
    
    return (n, M_deg, E_latest_deg, math.degrees(TA_rad), Sep_Dist_AU_OP, SA, PA, x, y, a1, a2)

def ENtoPASA(E, E_error, N, N_error):
    """
    Will calculate the Separation and Position Angles for a given East and North, including their errors.
    PA and error will be in [deg], with SA and error in ["]
    :returns: (PA,PA_error,SA,SA_error)
    """
    PA = math.degrees(np.arctan2(np.array(E),np.array(N)))
#     # correction factors for flaws of atan
#     #NOTE: this issue could be cured with proper use of np.arctan2(x1,x2), where it would normally be atan(x1/x2).
#     if (RA<0.0)and(DEC>0.0):
#         PA = 360.0+PA
#     elif DEC<0.0:
#         PA = 180.0+PA

    # arctan2 allows negative output angles, so make PA a positive angle if not one   
    if PA<0:
        PA = PA+360.0
    
    SA = math.sqrt(E**2.0 + N**2.0)
    
    PA_error=SA_error=0
    if (E_error==0)or(N_error==0):
        if False:
            print "either the E or N error value was zero, so setting the PA and SA return errors to zero!!"
    else:
        top = (E/N)*math.sqrt((E_error/E)**2.0 + (N_error/N)**2.0)
        btm = 1.0+(E/N)**2.0
        PA_error = abs(math.degrees(top/btm))

        top = SA*(abs(E*E_error)+abs(N*N_error))
        btm = E**2.0+N**2.0
        SA_error = abs(top/btm)
    
    return (PA,PA_error,SA,SA_error)

def PASAtoEN(PA,PA_error,SA,SA_error):
    """
    Convert provided Position Angle and Separation Angle, and their errors, into 
    RA and DEC with errors.  These are the same equations for calculating 
    x and y in the Thiele-Innes orbit fitting.  Remember that x and y are 
    flipped in that fitting approach due to how Thiele defined the coord 
    system when deriving the equations used.
    
    NOTE: this can also be used to calculate x and y used in Thiele-Innes
          With East=RA=y and North=DEC=x.  
    
    :returns: (E, E_error, N, N_error)
    """
    N = SA*math.cos(math.radians(PA))
    E = SA*math.sin(math.radians(PA))
    
    E_error=N_error=0
    if (SA_error==0)or(PA_error==0):
        if False:
            print "either the PA and SA error value was zero, so setting the E or N return errors to zero!!"
    else:
        tempA = (SA_error/SA)**2.0
        tempB = ((math.cos(math.radians(PA+PA_error))-math.cos(math.radians(PA)))/math.cos(math.radians(PA)))**2.0
        N_error = abs(N*math.sqrt(tempA+tempB))
        
        # Another way to calculate the error, but the one above is belived to be more currect 
        tempA2 = (SA_error*math.cos(math.radians(PA)))**2.0
        tempB2 = (SA*math.sin(math.radians(PA))*math.radians(PA_error))**2.0
        N_error2 = math.sqrt(tempA2+tempB2)
        
        tempC = (SA_error/SA)**2.0
        tempD = ((math.sin(math.radians(PA+PA_error))-math.sin(math.radians(PA)))/math.sin(math.radians(PA)))**2.0
        E_error = abs(E*math.sqrt(tempC+tempD))
        
        # Another way to calculate the error, but the one above is belived to be more currect 
        tempC2 = (SA_error*math.sin(math.radians(PA)))**2.0
        tempD2 = (SA*math.cos(math.radians(PA))*math.radians(PA_error))**2.0
        E_error2 = math.sqrt(tempC2+tempD2)
        
        if False:
            print 'N_error2-N_error = '+str(N_error2-N_error)
            print 'E_error2-E_error = '+str(E_error2-E_error)
            print 'E_error2 = '+str(E_error2)+', E_error = '+str(E_error)
            print 'N_error2 = '+str(N_error2)+', N_error = '+str(N_error)
    
    return (E, E_error, N, N_error)

