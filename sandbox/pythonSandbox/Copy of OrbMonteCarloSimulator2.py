import math
import numpy as np
import random as rand
from OrbCalculator import orbitCalculator

numOrbs = int(10e4)
print '\n*************  STARTING TO CREATE ORBITS  *****************'
print 'Number of orbits being created = ',numOrbs

#$$$$$$ change the obs epoch variable to use julian date or something else. 
#$$$$$$ ditto for last periapsis epoch param in mc
t = 10.0
Sep_Angle_arcsec_measured = 1.5
Sys_Dist_PC = 5e6
PA_deg_measured = 125.0
verbose = False

# init input lists
longAN_degs = []
es = []
Ts = []
a1s = []
a2s = []
periods = []

# init output lists
ns = []
Ms = []
Es = []
thetas = []
Sep_Dists = []
inclination_degs = []
argPerige_degs = []


for int in range(1,numOrbs+1):
    ## Input variables made with uniform random numbers
    longAN_deg = rand.uniform(0.0,180.0)
    longAN_degs.append(longAN_deg)
    e = rand.uniform(0.0,1.0)
    es.append(e)
    T = rand.uniform(0.1,t) #$$$$$$ change to diff date format
    Ts.append(T)
    a1 = rand.uniform(0.1,100.0)
    a1s.append(a1)
    a2 = rand.uniform(a1,20.0*a1)
    a2s.append(a2)
    period = rand.uniform(5.0,100.0)
    periods.append(period)
    
    # call orbitCalculator to take random variables and calc orbital elements
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU, inclination_deg, argPerige_deg) = \
        orbitCalculator(t, Sep_Angle_arcsec_measured, Sys_Dist_PC, PA_deg_measured, longAN_deg, e, T, a1, a2, period, verbose=verbose)
    
    # store output orbital elements in lists for later plotting
    ns.append(n)
    Ms.append(M_deg)
    Es.append(E_latest_deg)
    thetas.append(TA_deg)
    Sep_Dists.append(Sep_Dist_AU)
    inclination_degs.append(inclination_deg)
    argPerige_degs.append(argPerige_deg)
    
    
#TEMP print mean orbital element values for kicks
print '\nKnown/Measured Input Parameters:'
print 'Number of orbits being created ',numOrbs
print 'Epoch of observation/image [yrs] = ',t
print 'Measured Separation Angle of stars ["] = ',Sep_Angle_arcsec_measured
print 'Measured system distance from Earth [PC] = ',Sys_Dist_PC
print 'Measured Position Angle in image [deg] = ',PA_deg_measured

print '\nRandom Input Parameters means:'
print 'mean Longitude of Ascending Node [deg]= ', np.median(longAN_degs)
print 'mean e = ', np.median(es)
print 'mean T [yrs] = ', np.median(Ts)
print 'mean a1 [AU] = ',np.median(a1s)
print 'mean a2 [AU] = ',np.median(a2s)
print 'mean period [yrs] = ',np.median(periods)
    
print '\nOutput Parameters means:'
print 'mean Sep_Dist [AU] = ',np.median (Sep_Dists)
print 'mean True Anomaly [deg] = ', np.median(thetas)
print 'mean Eccentric Anomaly [deg] = ',np.median(Es)
print 'mean Mean Anomaly [deg] = ', np.median(Ms)
print 'mean Mean Motion [rad/yr] = ',np.median(ns)

print '\n******************* FINISHED CREATING ORBITS   *************************'
    