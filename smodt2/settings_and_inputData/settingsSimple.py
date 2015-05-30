#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp

simpleSettingsDict={
# The number of samples orbital parameters to try/draw [int]
'nSamples' : (100000,"Number of MCMC or MC samples"),
# Number of simulation chains to run in parallel, [1,100] [int].  
# NOTE: greater than numCores-1 causes system to slow down!
'nChains' : (20,"Number of chains"),
# Directory where you want the output data folder to go [string, at least 2 chars long]
'outDir' : '/run/media/kmede/Data1/Todai_Work/Data/data_SMODT',
# General filename for the simulation output folder to distinguish between simulation runs [string, at least 2 chars long]
#*************************************************************************************************************************
'outRoot' : "SMODT2-FakeData",
#*************************************************************************************************************************
# data mode, choices {'RV','DI','3D'} [string]
'dataMode' : ('3D',"Data Mode (RV,DI,3D)"),
# mode to run simulation in, choices {'MC','SA','MCMC'} [string]
'symMode' : ('MCMC',"Simulator mode (MC,SA,MCMC)"),
###################################################
# Ranges for acceptable random number inputs ######
###################################################
# Minimum/Maximum allowed value for the mass of the primary body [double][Msun]
'mass1MIN' : 0.1,
'mass1MAX' : 2.0,
# Minimum/Maximum allowed value for the mass of the secondary body [double][Msun]
'mass2MIN' : 0.01,
'mass2MAX' : 0.5,
# Minimum/Maximum allowed value for the system distance from Earth [double][PC]
'distMIN' : 2.0,
'distMAX' : 8.0,
# Minimum/Maximum allowed value for the Longitude of the Ascending Node [double][deg]
'OmegaMIN' : 20.0,
'OmegaMAX' : 110.0,
# Minimum/Maximum allowed value for the Eccentricity [double]
'eMIN' : 0.2,
'eMAX' : 0.6,
# Minimum/Maximum value for the Time of Last Periapsis (or Time of Center Transit) [JD]
#(-1 indicates to use [earliestsEpoch-period,earliestEpoch])
'TMIN' : 2456700,
'TMAX' : 2457350,
# Minimum/Maximum allowed value for the Period [double][yrs]
'PMIN' : 12.0,
'PMAX' : 17.0,
# Minimum/Maximum allowed value for the Inclination [double][deg]
'incMIN' : 1,
'incMAX' : 60.0,
# Minimum/Maximum allowed value for the Argument of Perigee [double][deg]
'omegaMIN' : 80,
'omegaMAX' : 140,
## Minimum/Maximum allowed value for the total semi-major axis [double][AU]{NOTE: only useful for DIonly simulations as RV requires separate a1,a2,M1,M2!}
#a_totMIN :2,
#a_totMAX :7,
# Minimum/Maximum value for Semi-major amplitude of RV curve [m/s]
'KMIN' : 1000,
'KMAX' : 1300,
# Minimum/Maximum values of Radial Velocity Offsets.  
# Must be one per set of RV data in same order as data comes in RVdata.dat, or the a single value to be used by all [comma separated list of doubles]
'vMINs' :[-550],
'vMAXs' :[550],
}


