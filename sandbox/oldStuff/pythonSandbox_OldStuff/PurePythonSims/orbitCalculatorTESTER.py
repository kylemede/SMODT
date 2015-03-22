from orbitToolbox import multiEpochOrbCalc
import math

#### MEASURED VALUES ########################                                        
### Currie, betaPic data
#SA_arcsec_measured_REALs = [0.411,0.210,0.326,0.345] # used to calc kai^2 (ACTUALLY NOT IN THIS VERSION, ONLY VERISON 1)(CURRIE)
#SA_mean_errors = [0.008,0.027,0.013,0.012] (CURRIE)
#PA_deg_measured_REALs = [31.7,211.49,210.64,209.8] # used to calc kai^2 (CURRIE)
#PA_mean_errors = [1.3,1.9, 1.2,0.8] (CURRIE)
#epochs = [2452953.50,2454781.50,2455194.50,2455275.50] (CURRIE)
#Sys_Dist_PC = 19.3 (CURRIE)
#Mass1 = 1
#Mass2 = 1
#### Liu, 2MAS J1534... data
SA_arcsec_measured_REALs = [0.199,0.1912,0.1906,0.1580,0.1537,0.1144,0.1020]
SA_mean_errors = [0.0011,0.0011,0.0003,0.0006,0.0004,0.0011,0.0004] # mean simply implies the mean of the + and - uncertainties                                         
PA_deg_measured_REALs = [14.5,15.5,15.43,17.5,15.53,21.5,20.4]
PA_mean_errors = [0.6,0.4,0.12,0.2,0.13,0.9,1.5] # mean simply implies the mean of the + and - uncertainties 
epochs = [2453754.5,2453826.5,2453860.5,2454185.5,2454212.5,2454480.5,2454557.5]
Sys_Dist_PC = 13.5
Mass1 = 1
Mass2 = 1
## Initial values that were found in the past #### maybe create a better method than basicOrbSim to make these.
longAN_degInitial = 178.0
eInitial = 0.24
periodInitial = 15.2
TInitial = 2456024.5 - (periodInitial*365.0)# to convert it to last periapsis vs next periapsis
inclination_degInitial = 84.3
argPeri_degInitial = 13.0
a_totalInitial = 2.3
#### tauBoo data
#SA_arcsec_measured_REALs = [2.71, 2.87,2.82,1.93]
#SA_mean_errors = [0.05, 0.03, 0.04, 0.02]                                         
#PA_deg_measured_REALs = [31.3, 30.85, 33.2, 53.1]
#PA_mean_errors = [0.5, 0.03, 1.0, 0.6]
#epochs = [2451945.5, 2451143.5, 2451711.5, 2455590.5]
#Sys_Dist_PC = 15.0
#Mass1 = 1.3
#Mass2 = 0.4
#********************************************************************

## send random parameters along with known ones to multi-epoch orbit calculator
(chi_squared_total_cur, ns, Ms, Es, thetas, Sep_Dists, SA_deg_measured_models, PA_deg_measured_models, a1s, a2s) =\
multiEpochOrbCalc(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors, 1000,\
               epochs, Sys_Dist_PC, inclination_degInitial, longAN_degInitial, eInitial, TInitial, periodInitial, argPeri_degInitial, a_totalInitial,\
               Mass1=Mass1, Mass2=Mass2, verbose=True)
    
print '\nINPUTS:'
print 'SA_arcsec_measured_REALs '+repr(SA_arcsec_measured_REALs)
print 'PA_deg_measured_REALs '+repr(PA_deg_measured_REALs)
print 'epochs '+repr(epochs)
print 'Sys_Dist_PC '+str(Sys_Dist_PC)
print 'longAN_degInitial '+str(longAN_degInitial)
print 'eInitial '+str(eInitial)
print 'TInitial '+str(TInitial)
print 'periodInitial '+str(periodInitial)
print 'inclination_degInitial '+str(inclination_degInitial)
print 'argPeri_degInitial '+str(argPeri_degInitial)
print 'a_totalInitial '+str(a_totalInitial)


print '\nOUTPUTS:'
print 'chi_squared_total_cur '+str(chi_squared_total_cur)
print 'ns '+repr(ns)
print 'Ms '+repr(Ms)
print 'Es '+repr(Es)
print 'thetas '+repr(thetas)
print 'Sep_Dists '+repr(Sep_Dists)
print 'SA_deg_measured_models '+repr(SA_deg_measured_models)
print 'PA_deg_measured_models '+repr(PA_deg_measured_models)
print 'a1s '+repr(a1s)
print 'a2s '+repr(a2s)


e = eInitial
#M = (math.pi*1.9)
#E_last = 2*math.pi
#E_latest = M+e*math.sin(M)
#
#print "\nInputs to Newton's Method are : \nM [rad]= "+str(M)+"\nE_last [rad]= "+\
#    str(E_last)+"\nE_latest [rad] = "+str(E_latest)+"\ne = "+str(e)
#
#print "\nStarting to run Newton's while loop."
#count = 0 # a counter to stop inf loops in Newtons method below
#while (abs(E_last-E_latest) > (1.0e-10))and(count<100):
#    print 'current E [rad]= ', E_latest
#    E_last = E_latest
#    #E_latest = E_last - ((M-E_last+e*math.sin(E_last))/(e*math.cos(E_last)-1.0))\
#    E_latest = E_last - ((E_last-M-e*math.sin(E_last))/(-e*math.cos(E_last)+1.0))
#    count = count+1
#    
#E_latest_deg = math.degrees(E_latest) # convert resulting E to degrees
#print "The resultant E value is [deg] = ", E_latest_deg
#print '-'*50  
#print 'e = '+str(e)  
#for i in range(50):
#    E_latest = i*((math.pi*2.0)/50.00)
#    print 'E [deg]= '+str(math.degrees(E_latest))
#    TA_rad  = math.acos((math.cos(E_latest)-e)/(1-e*math.cos(E_latest)))
#    if E_latest>math.pi:
#        TA_rad = 2.0*math.pi - TA_rad
#    TA_deg = math.degrees(TA_rad)
#    print 'True Anomaly [deg] = ', TA_deg

