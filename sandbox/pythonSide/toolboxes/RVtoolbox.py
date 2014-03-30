#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import math
import gc
import numpy as np
import os
import pylab
from math import pi
from orbitToolboxDuo import *

def TAcalculator2(t,e, T, period, T_center=0, verbose=False, debug=False):
    """
    
    $$$ Copy of exact function from orbitToolboxDuo as something odd was hapening $$$$
    
    This is the same as TAcalculator but with some updates to the Newton's method loop
    solving the Kepler's equation found in Double Stars by Heintz.
    $$$$$$$$$$$$ Later testing of the two versions will prove which is better. $$$$$$$$$$$$$
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

def RVdataToDict(filename):
    """ 
    Data must be in the columns:
    obsDate[JD]  RV[m/s]  RV_error[m/s]
    Jitter can be included as a 4th column value for the first line of each dataset.
    These jitter values will be encorporated into the rv errors using new_err = sqrt(err**2+jitter**2)
    Note: this same process takes place in the C++ code (But don't worry, it isn't done twice in a row :-P ).
    
    Spaces after rows of data indicates the end of a dataset and the following
    data will be loaded as a separate dataset.
    
    
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
    numLines = len(lines)
    file.close()
    
    lineNum = 0
    gotNumColumns = False
    numColumns = 0
    while (lineNum<numLines) and (gotNumColumns==False):
        if lines[lineNum][0]=='#':
            #comment line so do nothing
            gotNumColumns = False
        elif len(lines[lineNum])<=2:
            #blank line so do nothing
            gotNumColumns = False
        else:
            gotNumColumns = True
            numColumns = len(lines[lineNum].split())
        lineNum = lineNum+1
    
    dataDict = {}
    obsDates2 = []
    RVs2 = []
    RV_errors2 = []
    Jitters2 = []
    reachedEndOfDataset = True
    obsDates = []
    RVs = []
    RV_errors = []
    Jitters = [0]
    numEpochs = 0

    if ((numColumns is not 4) and (numColumns is not 3)):
        print 'WARNING: There must be 3 or 4columns of data to use RVdataToDict!'
        print str(numColumns)+" were found."
        print 'The line the columns were counted with was:'+lines[2]
    else:
        reachedEndOfDataset = True
        obsDates = []
        RVs = []
        RV_errors = []
        Jitters = [0]
        for lineNum in range(2,numLines):
            if verbose:
                print 'cur line:'+lines[lineNum]
            # Check if line is a comment
            if lines[lineNum][0]!='#':
                if reachedEndOfDataset:
                    obsDates = []
                    RVs = []
                    RV_errors = []
                    Jitters = []
                    
                vals =  lines[lineNum].split()
                if len(vals)==3 or len(vals)==4:
                    if verbose:
                        print 'Found a dataline on linenum '+str(lineNum)+", and it was:"+lines[lineNum]
                    reachedEndOfDataset = False
                    try:
                        obsDates.append(float(vals[0]))
                        RVs.append(float(vals[1]))
                        RV_errors.append(float(vals[2]))
                        numEpochs+=1
                        if verbose:
                            print str(len(vals))+" values were loaded the data arrays"
                    except:
                        if verbose:
                            print "WARNING: found line that didn't have good data in RVdataToDict"
                            print "line was:"+lines[lineNum]
                            print "The line had "+str(len(vals))+", values extracted"
                            print "continuing to next line."
                            if len(vals)>=3:
                                print 'The first value was:'+str(float(vals[0]))
                                print 'The second value was:'+str(float(vals[1]))
                                print 'The third value was:'+str(float(vals[2]))
                                
                if len(vals)==4:
                    reachedEndOfDataset = False
                    Jitters.append(float(vals[3]))
                if (((len(vals)==0)or(lineNum==(numLines-1)))and(reachedEndOfDataset==False)):
                    # blank line, or last line in file, so it means we are onto next dataset
                    if verbose:
                        print "loading a dataset of length "+str(len(RVs))+' into the ouput data array'
                    reachedEndOfDataset = True
                    obsDates2.append(obsDates)
                    RVs2.append(RVs)
                    RV_errors2.append(RV_errors)
                    Jitters2.append(Jitters)
       
    # all lines loaded int lists, so convert into dict
    dataDict['RV_epochs'] = obsDates2
    dataDict['RVs'] = RVs2
    dataDict['RV_errors'] = RV_errors2
    dataDict['numEpochs'] = numEpochs
    
    if verbose:
        print 'A total of '+str(len(RVs2))+" epochs of RV data were loaded"
    
    ## Update errors with the jitter values if they exist.
    RV_errors2_updated = []
    for dataset in range(0,len(RVs2)):
        jitter = 0
        RV_errors_updated = []
        for epoch in range(0,len(RV_errors2[dataset])):
            if len(Jitters2[dataset])==0:
                Jitters2[dataset] = [0]
            elif len(Jitters2[dataset])==1:
                jitter = Jitters2[dataset][0]
            elif len(Jitters2[dataset])==len(RV_errors2[dataset]):
                jitter = Jitters2[dataset][epoch]   
            errorOut = math.sqrt(RV_errors2[dataset][epoch]**2.0 + jitter**2.0)
            RV_errors_updated.append(errorOut)
            if verbose:
                print 'Input RV_error = ',RV_errors2[dataset][epoch]
                print 'Input Jitter = ',jitter
                print 'Output RV_error = ',errorOut
        RV_errors2_updated.append(RV_errors_updated)
        
    dataDict['RV_errors'] = RV_errors2_updated
    dataDict['jitters'] = Jitters2

    return dataDict
   
def rvResidualWithoutPlanetResidual():
    """
    A totally temporary function to calculate the residual velocity of the primary star
    WITHOUT the residual from the planet.  Thus only having the velocity that could be 
    caused by the secondary star/companion.
    """

    from paramSettingsDict import paramSettingsDict
    RV_epochs = paramSettingsDict['Data']['RV_epochs']
    #RVerrors = paramSettingsDict['Data']['RVerrors']
    RVs = paramSettingsDict['Data']['RVs']
    M1 = paramSettingsDict['Data']['M1']
    # data for planet's orbit that are needed to calculate proper chiSquare for RV fit
    # all below values for the planet orbit and RV data are from Butler2006.
    K_p = 461.1 #[m/s]
    p_p = 3.31246   #[days]
    #planet_p_years = planet_p/365.25
    e_p = 0.023    
    argPeri_p = 188.0   #[deg]
    T_p = 2446957.8   #[JD]
    sigma_jitter = 15.0    #[m/s]
    
    rvWithoutPlanetResiduals = []
    for RVdataSet in range(0,len(RVs)):
        print '\nworking on RVdataSet ',RVdataSet
        for epoch in range(0,len(RV_epochs[RVdataSet])):
            v_r_p = vrCalculatorPlanet(RV_epochs[RVdataSet][epoch], e_p, T_p, p_p, argPeri_p, M1, M2SineI=False, K=K_p, verbose=False)
            RVprimary = RVs[RVdataSet][epoch]
            rvWithoutPlanetResidual = RVprimary - v_r_p
            rvWithoutPlanetResiduals.append(rvWithoutPlanetResidual)
            print v_r_p
 
def rv1bodyCalculator(RV_epochs, RVs, RVerrors, sigma_jitter, i, p, e, T, argPeri, a, verbose=False):
    """
    This is for calculating the RV and resulting chiSquared due to a planet ONLY, thus no companion star.
    Naturally this system would have a planet with mass<< primary star mass.
    
    
    sigma_jitter must be the unsquared version and in units of [m/s]
    i in degrees
    argPeri in degrees
    a in AU
    T in JD
    t in JD
    e unitless
    period in years
    """
    
    chiSquaredTotal = 0.0
    
    if verbose:
        print '\n## Using the semi-major to calculate K for the primary star due to companion star ##'
    
    for epoch in range(0,len(RVs)):
        if verbose:
            print '\n******************* RV values for epoch '+str(epoch+1)+' *****************'
        # calculate the velocity residual due to the companion 
        (v_r_c, K_c) = vrCalculatorStar2(RV_epochs[epoch],e, T, p, argPeri, a, i=i, K=False, verbose=False)
        
        # calculate chiSquared for this epoch
        chiSquaredCurr = chiSquaredCalc(RVs[epoch], (RVerrors[epoch]+sigma_jitter), (v_r_c))
        
        chiSquaredTotal = chiSquaredTotal + chiSquaredCurr
        
        if verbose:
            print 'RV = (calculated)-measured = ('+str(v_r_c)+') - '+str(RVs[epoch])+' = '+str((v_r_c)-RVs[epoch]) 
            print '(RVerrors[epoch]+sigma_jitter) = '+str((RVerrors[epoch]+sigma_jitter))
            print 'chiSquaredCurr = ',chiSquaredCurr
    
    if verbose:
            print '\nK of companion = ',K_c
            
    return chiSquaredTotal

def rv1bodyCalculator2(RV_epochs, RVs, RVerrors, sigma_jitter, p, e, T, argPeri, M1,M2SineI, verbose=False):
    """
    This is for calculating the RV and resulting chiSquared due to a planet ONLY, thus no companion star.
    Naturally this system would have a planet with mass<< primary star mass.
    
    
    sigma_jitter must be the unsquared version and in units of [m/s]
    i in degrees
    argPeri in degrees
    a in AU
    T in JD
    t in JD
    e unitless
    period in years
    """
    
    chiSquaredTotal = 0.0
    
    if verbose:
        print '\n## Using the semi-major to calculate K for the primary star due to companion star ##'
    
    for epoch in range(0,len(RVs)):
        if verbose:
            print '\n*** RV values for epoch '+str(epoch+1)+' ***'
        # calculate the velocity residual due to the companion 
        (v_r_c, K_c) = vrCalculatorPlanet(RV_epochs[epoch],e, T, p, argPeri,M1,M2SineI=M2SineI, K=False, verbose=False)
        # calculate chiSquared for this epoch
        chiSquaredCurr = chiSquaredCalc(RVs[epoch], (RVerrors[epoch]+sigma_jitter), (v_r_c))
        
        chiSquaredTotal = chiSquaredTotal + chiSquaredCurr
        
        if verbose:
            print 'RV = (calculated)-measured = ('+str(v_r_c)+') - '+str(RVs[epoch])+' = '+str((v_r_c)-RVs[epoch]) 
            print '(RVerrors[epoch]+sigma_jitter) = '+str((RVerrors[epoch]+sigma_jitter))
            print 'chiSquaredCurr = ',chiSquaredCurr
    
    if verbose:
            print '\nK of companion = ',K_c
            
    return chiSquaredTotal 
        
def rv2bodyCalculator3(RV_epochs, RVs, RVerrors, sigma_jitter, M1, M2_c, i_c, p_c, e_c, T_c, argPeri_c, \
                                                                    K_p, p_p, e_p, argPeri_p, T_p, verbose=False):
    """
    NOTE: This version uses the K equation with each objects mass.  
        Use rv2bodyCalculator4 to use the semi-major axis instead.
    
    This will calclulate the radial velocity residue based chiSquared considering RV data for the primary star, orbital elements 
    for the companion star (indicated with '_c' parameters), and some of the parameters for the planet around the primary
    star's orbital elements (indicated with the '_p' parameters).
    
    sigma_jitter must be the unsquared version and in units of [m/s]
    M2SineI_c and M1 in kg
    i in degrees
    argPeri in degrees
    a in AU
    T in JD
    t in JD
    e unitless
    period in years
    """
    
    chiSquaredTotal = 0.0
    
    if verbose:
        print '\n## Using the masses to calculate K for the primary star due to companion star ##'
    for epoch in range(0,len(RVs)):
        if verbose:
            print '\n*** RV values for epoch '+str(epoch+1)+' ***'
        # calculate the velocity residual due to the companion star
        (v_r_c, K_c) = vrCalculatorStar(RV_epochs[epoch],e_c, T_c, p_c, argPeri_c, M1, M2=M2_c, i=i_c, K=False, verbose=False)
        # calculate the velocity residual due to the planet around primary
        (v_r_p, K_p) = vrCalculatorPlanet(RV_epochs[epoch], e_p, T_p, p_p, argPeri_p, M1, M2SineI=False, K=K_p, verbose=False)

        # calculate chiSquared for this epoch
        chiSquaredCurr = chiSquaredCalc(RVs[epoch], (RVerrors[epoch]+sigma_jitter), (v_r_c+v_r_p))
        
        chiSquaredTotal = chiSquaredTotal + chiSquaredCurr
        
        if verbose:
            print 'RV = ('+str(v_r_c)+' + '+str(v_r_p)+') - '+str(RVs[epoch])+' = '+str((v_r_c+v_r_p)-RVs[epoch]) 
            print '(RVerrors[epoch]+sigma_jitter) = '+str((RVerrors[epoch]+sigma_jitter))
            print 'chiSquaredCurr = ',chiSquaredCurr
            
    return chiSquaredTotal

def rv2bodyCalculator4(RV_epochs, RVs, RVerrors, sigma_jitter, M1, a1, i_c, p_c, e_c, T_c, argPeri_c, \
                                                                    K_p, p_p, e_p, argPeri_p, T_p, verbose=False):
    """
    NOTE: This version uses the K equation with the semi-major axis instead of each objects masses.
          Use rv2bodyCalculator3 to use the objects masses instead.
    
    This will calclulate the radial velocity residue based chiSquared considering RV data for the primary star, orbital elements 
    for the companion star (indicated with '_c' parameters), and some of the parameters for the planet around the primary
    star's orbital elements (indicated with the '_p' parameters).
    
    sigma_jitter must be the unsquared version and in units of [m/s]
    M2SineI_c and M1 in kg
    i in degrees
    argPeri in degrees
    a in AU
    T in JD
    t in JD
    e unitless
    period in years
    """
    
    chiSquaredTotal = 0.0
    
    if verbose:
        print '\n## Using the semi-major to calculate K for the primary star due to companion star ##'
    
    for epoch in range(0,len(RVs)):
        if verbose:
            print '\n*** RV values for epoch '+str(epoch+1)+' ***'
        # calculate the velocity residual due to the companion star
        (v_r_c, K_c) = vrCalculatorStar2(RV_epochs[epoch],e_c, T_c, p_c, argPeri_c, a1, i=i_c, K=False, verbose=False)
        # calculate the velocity residual due to the planet around primary
        (v_r_p, K_p) = vrCalculatorPlanet(RV_epochs[epoch], e_p, T_p, p_p, argPeri_p, M1, M2SineI=False, K=K_p, verbose=False)
        
        # calculate chiSquared for this epoch
        chiSquaredCurr = chiSquaredCalc(RVs[epoch], (RVerrors[epoch]+sigma_jitter), (v_r_c+v_r_p))
        
        chiSquaredTotal = chiSquaredTotal + chiSquaredCurr
        
        if verbose:
            print 'RV = ('+str(v_r_c)+' + '+str(v_r_p)+') - '+str(RVs[epoch])+' = '+str((v_r_c+v_r_p)-RVs[epoch]) 
            print '(RVerrors[epoch]+sigma_jitter) = '+str((RVerrors[epoch]+sigma_jitter))
            print 'chiSquaredCurr = ',chiSquaredCurr
            
    return chiSquaredTotal

def vrCalculatorPlanet(t,e,T,period,argPeri,M1,T_center=0,M2SineI=False, K=False, verbose=False):
    """
    Version for when the companion is a planet.  Thus it calculates the residual velocity 
    assuming this and the consequence that (M1+M2~=M1).  Returns residual velocity 
    of the primary star!! due to the planet, not the the vel of the planet.
    
    M1 in Msun
    M2SineI in Mjupiter
    argPeri in degrees
    a in AU
    T in JD
    t in JD
    e unitless
    period in days
    
    K=False implies that the provided values will be used to 
    calculate it from scratch, else the provided value will
    be used to calculate the radial velocity residual.  Thus,
    this param must be either False or a float in units of 
    m/s.
    
    returned vel will be in m/s
    """
        
    if (K==False) and (M2SineI==False):
        print 'vrCalcPlanet1225: The value of M2SineI and K cannot both be False, one MUST be defined'
    if (K==False) and (M2SineI!=False):
        if period==0:
            K=0
        else:
            # convert units of days to seconds
            period_seconds = period*86400.0#31557600.0
            M1_kg = M1*1.98892e30
            M2SineI_kg = M2SineI*1.8986e27
            G = 6.67300e-11
            
            # Calc K in parts to see equation better
            A = ((2.0*pi*G)/period_seconds)**(1.0/3.0)
            B = M2SineI_kg/(M1_kg**(2.0/3.0))
            C = 1.0/math.sqrt(1.0-e**2.0)
            # put it all together
            K = A*B*C
    
    if verbose:
        print 'K_planet = ',K
    
    if K==0:
        v_r = 0
    else:
        period_years = period/365.25
        
        ## get the TAcalculator to find the TA in radians
        (n, M_deg, E_latest_deg,TA_rad) = TAcalculator2(t,e, T, period_years, T_center=T_center, verbose=False, debug=False)
        
        v_r = K*(math.cos(math.radians(argPeri)+TA_rad)+e*math.cos(math.radians(argPeri)))

    return (v_r, K)

def vrCalculatorStar(t,e,T,period,argPeri,M1,M2,T_center=0,i=False, K=False, verbose=False):
    """
    NOTE: This version is the one that uses the objects masses to calculate K.
          To use the one that calculates it using the semi-major axis, use vrCalculatorStar2.
    
    M1 and M2 in Msun
    argPeri in degrees
    i in degrees
    T in JD
    t in JD
    e unitless
    period in years
    
    K=False implies that the provided values will be used to 
    calculate it from scratch, else the provided value will
    be used to calculate the radial velocity residule.  Thus,
    this param must be either False or a float in units of 
    m/s.
    
    returned vel will be in m/s
    """
        
    if (K==False) and (i==False):
        print 'vrCalcStar1268: The value of i and K cannot both be False, one MUST be defined'
    if (K==False) and (i!=False):
        if period==0:
            K=0
        else:
            # convert units of years to seconds
            period_seconds = period*31557600.0
            M1_kg = M1*1.98892e30
            M2_kg = M2*1.98892e30
            G = 6.67300e-11
            
            # Calc K in parts to see equation better
            A = ((2.0*pi*G*(M1_kg+M2_kg))/period_seconds)**(1.0/3.0)
            B = M2_kg/(M1_kg+M2_kg)
            C = math.sin(math.radians(i))/math.sqrt(1.0-e**2.0)
            # put it all together
            K = A*B*C

    if verbose:
        print 'K_Star = ',K
        
    if K==0:
        v_r = 0
    else:
        ## get the TAcalculator to find the TA in radians
        (n, M_deg, E_latest_deg,TA_rad) = TAcalculator2(t,e, T, period, T_center=T_center, verbose=False, debug=False)
        
        v_r = K*(math.cos(math.radians(argPeri)+TA_rad)+e*math.cos(math.radians(argPeri)))

    return (v_r, K)

def vrCalculatorStar2(t,e,T,period,argPeri,a1,T_center=0,i=False, K=False, verbose=False):
    """
    NOTE: this is the version which uses the K equation with the semi-major axis of primary's orbit.
          To use the masses instead, use vrCalculatorStar.
    
    argPeri in degrees
    i in degrees
    T in JD
    t in JD
    e unitless
    period in years
    a1 in AU
    
    K=False implies that the provided values will be used to 
    calculate it from scratch, else the provided value will
    be used to calculate the radial velocity residule.  Thus,
    this param must be either False or a float in units of 
    m/s.
    
    returned vel will be in m/s
    """
        
    if (K==False) and (i==False):
        print 'vrCalcStar1268: The value of i and K cannot both be False, one MUST be defined'
    if (K==False) and (i!=False):
        # convert units of years to seconds
        seconds_per_yr = 31557600.0
        period_seconds = period*seconds_per_yr
        #M1_kg = M1*1.98892e30
        #M2_kg = M2*1.98892e30
        meters_per_AU = 149598000000.0
        a1_meters = a1*meters_per_AU
        
#        # Calc K in parts to see equation better
#        A = ((2.0*pi*6.67e-11*(M1_kg+M2_kg))/period_seconds)**(1.0/3.0)
#        B = M2_kg/M1_kg
#        C = math.sin(math.radians(i))/math.sqrt(1.0-e**2)
#        # put it all together
#        K = A*B*C
        
        # Calc K in parts to see equation better
        A = (2.0*pi)/period_seconds
        B = a1_meters
        C = math.sin(math.radians(i))/math.sqrt(1.0-e**2.0)
        # put it all together
        K = A*B*C
        
    if verbose:
        print '\nperiod = ',period
        print 'a1 = ',a1
        print 'e = ',e
        print 'i = ',i
        print 'K_Star = ',K
        
    ## get the TAcalculator to find the TA in radians
    (n, M_deg, E_latest_deg,TA_rad) = TAcalculator2(t,e, T, period, T_center=T_center, verbose=False, debug=False)
    if verbose:
        #print '#######################################################'
        print 'TA = '+str(math.degrees(TA_rad))  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        if False:
            print 'math.degrees(math.radians(argPeri-pi)+TA_rad) = '+str(math.degrees(math.radians(argPeri-pi)+TA_rad))
            print 'math.cos(math.radians(argPeri-pi)+TA_rad) = '+str(math.cos(math.radians(argPeri-pi)+TA_rad))
            print 'e = '+str(e)
            print 'argPeri = '+str(argPeri)
            print 'e*math.cos(math.radians(argPeri-pi)) = '+str(e*math.cos(math.radians(argPeri-pi)))
            print '#######################################################'
    v_r = K*(math.cos(math.radians(argPeri)+TA_rad)+e*math.cos(math.radians(argPeri)))

    return (v_r, K)


