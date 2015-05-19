#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import math
import gc
import numpy as np
import os
import pylab
from math import pi
from DItoolbox import *
import generalToolbox as genTools



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
    
    if os.path.exists(filename):
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
    else:
        print "filename '"+filename+"' does not exist!!!"
   
def subtractPlanetRV(rvDataFilename,e_p,p_p,K_p,argPeri_p,T_p,Tc_p):
    """
    This function will subtract the RV due to a companion planet leaving only the residual
    that could be caused by a secondary star or another planet.  If the system you are 
    investigating has a long period companion star, then this is useful to subtract 
    the RV due to the circum-primary star with known orbital elements and focus the simulation
    on solving only for the companion star.  Calculating both for each proposed set of 
    elements of the companion star during the simulation is un-necesarily CPU taxing.
    """
    rvDataDict = RVdataToDict(rvDataFilename)
    RV_epochs = rvDataDict['RV_epochs']
    #RVerrors = paramSettingsDict['Data']['RVerrors']
    RVs = rvDataDict['RVs']
    
    rvWithoutPlanetResiduals = []
    for RVdataSet in range(0,len(RVs)):
        print '\nworking on RVdataSet ',RVdataSet
        for epoch in range(0,len(RV_epochs[RVdataSet])):
            v_r_p = vrCalculatorPlanetMassType(RV_epochs[RVdataSet][epoch], e_p, T_p, p_p, argPeri_p, M1=0, T_center=Tc_p, M2SineI=False, K=K_p, verbose=False)
            RVprimary = RVs[RVdataSet][epoch]
            rvWithoutPlanetResidual = RVprimary - v_r_p
            rvWithoutPlanetResiduals.append(rvWithoutPlanetResidual)
            print "Raw RV = "+str(RVprimary)+", planet RV = "+str(v_r_p)
            print "Raw-planet => residual = "+str(rvWithoutPlanetResidual)
    return rvWithoutPlanetResiduals

def vrCalculatorPlanetMassType(t,e,T,period,argPeri,M1,T_center=0,M2SineI=False, K=False, verbose=False):
    """
    NOTE: Specific version ONLY for a companion planet, general version is vrCalculatorMassType.
    
    Version ONLY for when the companion is a planet.  Thus it calculates the residual velocity 
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
        print 'vrCalcPlanet: The value of M2SineI and K cannot both be False, one MUST be defined'
    if (K==False) and (M2SineI!=False):
        if period==0:
            K=0
        else:
            # convert units of days to seconds
            period_seconds = period*86400.0#31557600.0
            M1_kg = M1*1.98855e30
            M2SineI_kg = M2SineI*1.8983e27
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
        
        ## get the genTools.TAcalculator to find the TA in radians
        (n, M_deg, E_latest_deg,TA_rad) = genTools.TAcalculator(t,e, T, period_years, T_center=T_center, verbose=False, debug=False)
        
        v_r = K*(math.cos(math.radians(argPeri)+TA_rad)+e*math.cos(math.radians(argPeri)))

    return (v_r, K)

def vrCalculatorMassType(t,e,T,period,argPeri,M1,M2,T_center=0,i=False, K=False, verbose=False):
    """
    NOTE: This version is the one that uses the objects masses to calculate K. 
          It will work for a companion star or planet.  Keep in mind that the M2sin(i) needs to be handled
          before calling this function, providing the resulting M2 not M2sin(i) value.
          To use the one that calculates it using the semi-major axis, use vrCalculatorSemiMajorType.
    
    M1 and M2 in Msun
    argPeri in degrees
    i in degrees
    T in JD
    t in JD
    e unitless
    period in years
    
    K=False implies that the provided values will be used to 
    calculate it from scratch, else the provided value will
    be used to calculate the radial velocity residual.  Thus,
    this param must be either False or a float in units of 
    m/s.
    
    returned vel will be in m/s
    """
        
    if (K==False) and (i==False):
        print 'vrCalculatorMassType: The value of i and K cannot both be False, one MUST be defined'
    if (K==False) and (i!=False):
        if period==0:
            K=0
            print "vrCalculatorStarMassType: Both K and period = zero, thus cannot calculate RV!!"
        else:
            # convert units of years to seconds
            period_seconds = period*31556908.799999997
            M1_kg = M1*1.98855e30
            M2_kg = M2*1.98855e30
            G = 6.67384e-11
            
            # Calc K in parts to see equation better
            A = ((2.0*pi*G*(M1_kg+M2_kg))/period_seconds)**(1.0/3.0)
            B = M2_kg/(M1_kg)
            C = math.sin(math.radians(i))/math.sqrt(1.0-e**2.0)
            # put it all together
            K = A*B*C

    if verbose:
        print 'K = ',K
        
    if K==0:
        v_r = 0
    else:
        ## get the genTools.TAcalculator to find the TA in radians
        (n, M_deg, E_latest_deg,TA_rad) = genTools.TAcalculator(t,e, T, period, T_center=T_center, verbose=False, debug=False)
        
        v_r = K*(math.cos(math.radians(argPeri)+TA_rad)+e*math.cos(math.radians(argPeri)))

    return (v_r, K)

def vrCalculatorSemiMajorType(t,e,T,period,argPeri,a1,T_center=0,i=False, K=False, verbose=False):
    """
    Most recommended way to calculate the radial velocity as it avoids all the mass issues (reduced 
    formula or M2sin(i)).  Although please provide the proper semi-major axis value for either 
    the primary star if you are comparing to the primary's RV values, or the companion's if 
    using its RVs (ie. provide a1 if primary RV's, or a2 if companion's).
    
    NOTE: this is the version which uses the K equation with the semi-major axis of primary's orbit.
          This version will work for both a companion star or planet.
          To use the masses instead, use vrCalculatorStarMassType or vrCalculatorPlanetMassType.
    
    argPeri in degrees
    i in degrees
    T in JD
    t in JD
    e unitless
    period in years
    a1 in AU
    
    K=False implies that the provided values will be used to 
    calculate it from scratch, else the provided value will
    be used to calculate the radial velocity residual.  Thus,
    this param must be either False or a float in units of 
    m/s.
    
    returned vel will be in m/s
    """
        
    if (K==False) and (i==False):
        print 'vrCalcStar1268: The value of i and K cannot both be False, one MUST be defined'
    if (K==False) and (i!=False):
        # convert units of years to seconds
        seconds_per_yr = 31556908.799999997
        period_seconds = period*seconds_per_yr
        meters_per_AU = 149597870700.0
        a1_meters = a1*meters_per_AU
        
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
        
    ## get the genTools.TAcalculator to find the TA in radians
    (n, M_deg, E_latest_deg,TA_rad) = genTools.TAcalculator(t,e, T, period, T_center=T_center, verbose=False, debug=False)
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


