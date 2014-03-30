import math

def orbitCalculator(t, Sep_Angle_arcsec_measured, Sys_Dist_PC, PA_deg_measured, longAN_deg, e, T, a1, a2, period, verbose=False):
    """
    Inputs:
    @param t:                            = epoch of observation/image [yrs]***********
    @param Sep_Angle_arcsec_measured:    = measured Separation Angle of stars ["]
    @param Sys_Dist_PC:                  = measured system distance from Earth [PC]
    @param PA_deg_measured:              = measured Position Angle in image [deg]
    @param longAN_deg:                   = Longitude of Ascending Node [deg]
    @param e:                            = eccentricity of orbits [unitless]
    @param T:                            = Last Periapsis Epoch/time [yrs] ***********
    @param a1:                           = semi-major axis of M1 [AU]
    @param a2:                           = semi-major axis of M2 [AU]
    @param period:                       = period of orbits [yrs]
    @param verbose:                      = Send prints to screen? [True/False] (Default = False)
    
    Outputs:
    @param n:                            = Mean Motion [rad/yr]
    @param M_deg:                        = Mean Anomaly [deg]
    @param E_latest_deg:                 = Eccentric Anomaly [deg]
    @param TA_deg:                       = True Anomaly [deg]
    @param Sep_Dist_AU:                  = Separation Distance in orbital plane [AU]
    @param inclination_deg:              = inclination [deg]
    @param argPeri_deg:                  = Argument of Periapsis in orbital plane [deg]
    
    """
    if verbose:
        print '\n'+'*'*50
        print 'Starting to calculate Orbital Parameters\n'
        print 'Input variable values: '+'\nCurrent epoch time [yrs] = '+str(t)+\
        'Seperation Angle measured [arcsec] = '+ str(Sep_Angle_arcsec_measured)+\
        '\nSystem Distance from Earth [PC] = '+str(Sys_Dist_PC)+\
        '\nPosition Angle measured [deg] = '+str(PA_deg_measured)+\
        '\nEccentricity = '+str(e)+'\nLast Periapsis Epoch [yrs] = '+\
        str(T)+'\na1 [AU] = '+str(a1)+'\na2 [AU] = '+str(a2)+\
        '\nPeriod [yrs] = '+str(period)+'\nverbose = ',str(verbose)
    #----------------------------------------------------------------------

    # initial guess (E_last), will be updated in loop.  
    # Anything works, just takes longer if further from real value. => pi
    E_last = 3.14
    # stored initial value to be updated in loop
    E_latest = 1.0

    ## calculate the Mean Motion
    n = (2*math.pi)/period
    if verbose:
        print 'Mean Motion [rad/yr]= '+str(n)

    ## calculate Mean Anomaly
    M = n*(t-T)
    M_deg = math.degrees(M) # convert resulting M to degrees
    if verbose:
        print 'Mean Anomaly [deg]= ',M_deg

    ### Performing Newton's Method to get the Eccentric Anomaly value ###
    if verbose:
        print '-'*50

    # show input value to 
    if verbose:
        print 'Inputs to Newtons Method are : \nM [rad]= '+str(M)+"\nE_last [rad]= "+str(E_last)+\
        "\nE_latest [rad] = "+str(E_latest)+"\ne = "+str(e)
        
        print '\nStarting to run Newtons while loop.'

    count = 0 # a counter to stop inf loops in Newtons method below
    while (abs(E_last-E_latest) > (1.0e-10))and(count<100):
        if verbose:
            print 'current E [rad]= ', E_latest
        E_last = E_latest
        E_latest = E_last - ((M-E_last+e*math.sin(E_last))/(e*math.cos(E_last)-1.0))
        count = count+1

    E_latest_deg = math.degrees(E_latest) # convert resulting E to degrees
    if verbose:
        print "The resultant E value is [deg] = ", E_latest_deg
        print '-'*50

    ## calculate True Anomaly from Eccentric Anomaly
    TA_rad  = math.acos((math.cos(E_latest)-e)/(1-e*math.cos(E_latest)))
    TA_deg = math.degrees(TA_rad)
    if verbose:
        print '\nTrue Anomaly [deg] = ', TA_deg
        
    ## calculate the Separation Distance in orbital plane
    Sep_Dist_AU = ((a1+a2)*(1-e*e))/(1+e*math.cos(TA_rad))
    if verbose:
        print '\nSeparation Distance in orbital plane [AU]= ',Sep_Dist_AU

    ## calculate inclination angle
    inclination_rad = (Sep_Angle_arcsec_measured*Sys_Dist_PC)/Sep_Dist_AU # PC->AU and rad->arcsec cancel
    inclination_deg = math.degrees(inclination_rad)
    
    ## calculate True Anomaly measured in reference plane
    TA_deg_measured = TA_deg*math.cos(inclination_rad)
    
    ## calculate Argument of Periapsis in reference plane
    argPeri_deg_measured = PA_deg_measured - longAN_deg - TA_deg_measured
    
    ## convert Argument of Periapsis in ref plane to orbital plane
    argPeri_deg = argPeri_deg_measured/math.cos(inclination_rad)
    if verbose:
        print '\nArgument of Periapsis in Orbital Plane [deg] = ',argPeri_deg
    
    if verbose:
        print '\nFinished calculating Orbital Parameters'
        print '*'*50+'\n'
        
    return (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU, inclination_deg, argPeri_deg)

    
   
