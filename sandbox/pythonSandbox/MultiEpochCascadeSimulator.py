#@Author: Kyle Mede, kylemede@gmail.com, written for the University of Tokyo during UTRIP 2011.
import math
import numpy as np
import random as rand
import pylab 
from orbitToolbox import *
    
#*********** In this brick is all you need to play with *************    
# number of orbits to calculate for each set of inputs    
numOrbs = int(1e6)
chiSquaredMax = 0.5
Sys_Dist_PC = 19.3 # same system, so same system distance for each epoch
verbose = False

Sep_Angle_arcsec_measured_REALs = [0.210,0.326] # used to calc kai^2 
PA_deg_measured_REALs = [211.49,210.64] # used to calc kai^2
epochs = [2454781.50,2455194.50]
#*********************************************************************

print '\n*************  STARTING TO CREATE ORBITS  *****************'
print 'Number of orbits being created = ',numOrbs

# init input lists
longAN_degs = []
es = []
Ts = []
a1s = []
a2s = []
periods = []
inclination_degs = []
argPeri_degs = []

# init output lists
ns2 = []
Ms2 = []
Es2 = []
thetas2 = []
Sep_Dists2 = []
Sep_Angle_arcsec_measured_models2 = []
PA_deg_measured_models2 = []

for int in range(1,numOrbs+1):
           
    ## Input variables made with uniform random numbers
    # Note: These will be stored if they satisfy the chi**2 max
    #       so only the good ones are stored to save processor time
    longAN_deg = rand.uniform(0.0,180.0)
    e = rand.uniform(0.0,1.0)
    period = rand.uniform(1.0,100.0) # [yrs]
    T = rand.uniform(epochs[0]-period*365.0,epochs[0])# thus between (now-one period) and now
    a1 = 0.0#rand.uniform(0.1,20.0) #$$$$$$$$$$$$$$$$$$$$$$
    a2 = rand.uniform(0.1,50.0) #rand.uniform(a1,10.0*a1) #####$$$$$$$$$$$$$$$$$
    inclination_deg = rand.uniform(0.0,180.0)
    argPeri_deg = rand.uniform(0.0,180.0)
    
    (chi_squared_total, ns, Ms, Es, thetas, Sep_Dists, Sep_Angle_arcsec_measured_models, PA_deg_measured_models) =\
    multiEpochOrbCalc(Sep_Angle_arcsec_measured_REALs, PA_deg_measured_REALs, epochs, Sys_Dist_PC,\
                       inclination_deg, longAN_deg, e, T, a1, a2, period, argPeri_deg, verbose=verbose)
        
        ##### while loop ends here!
    if (chi_squared_total != False) and (chi_squared_total<=chiSquaredMax):
        # store input values to the model
        longAN_degs.append(longAN_deg)
        es.append(e)
        Ts.append(T)
        a1s.append(a1)
        a2s.append(a2)
        periods.append(period)
        inclination_degs.append(inclination_deg)
        argPeri_degs.append(argPeri_deg)
        
        # store output orbital elements in lists for later plotting
        ns2.append(ns)
        Ms2.append(Ms)
        Es2.append(Es)
        thetas2.append(thetas)
        Sep_Dists2.append(Sep_Dists)
        Sep_Angle_arcsec_measured_models2.append(Sep_Angle_arcsec_measured_models)
        PA_deg_measured_models2.append(PA_deg_measured_models)
    
#TEMP print mean orbital element values for kicks
print '\n*************** RESULTS SUMMARY  *********************'
print 'Measured Separation Angle of stars REAL ["] = ',Sep_Angle_arcsec_measured_REAL
print 'Measured Position Angle in image REAL [deg] = ',PA_deg_measured_REAL
print 'Number of orbits being created ',numOrbs
print 'chi**2 max = ',chiSquaredMax
print 'number of orbits with chi**2 under max = ',len(Sep_Dists2)

print '\nKnown/Measured Input Parameters:'
print 'Epoch of observation/image [yrs] = ',t
print 'Measured system distance from Earth [PC] = ',Sys_Dist_PC

print '\nRandom Input Parameters means:'
print 'median Longitude of Ascending Node [deg]= ', np.median(longAN_degs)
print 'median e = ', np.median(es)
print 'median T [julian date] = ', np.median(Ts)
print 'median a1 [AU] = ',np.median(a1s)
print 'median a2 [AU] = ',np.median(a2s)
print 'median period [yrs] = ',np.median(periods)
print 'median inclination [deg] = ',np.median(inclination_degs)
print 'median Argument of Pariapsis [deg] = ',np.median(argPeri_degs)
    
print '\nOutput Parameters means:'
print 'median Sep_Dist [AU] = ',np.median (Sep_Dists2)
print 'median True Anomaly [deg] = ', np.median(thetas2)
print 'median Eccentric Anomaly [deg] = ',np.median(Es2)
print 'median Mean Anomaly [deg] = ', np.median(Ms2)
print 'median Mean Motion [rad/yr] = ',np.median(ns2)
print 'median Measured Separation Angle of stars MODEL ["] = ',np.median(Sep_Angle_arcsec_measured_models2)
print 'median Measured Position Angle in image MODEL [deg] = ',np.median(PA_deg_measured_models2)
(n,bins,patches) = pylab.matplotlib.pyplot.hist(inclination_degs2,bins=100)
pylab.show()
print '\n******************* FINISHED CREATING ORBITS   *************************\n'
    
           