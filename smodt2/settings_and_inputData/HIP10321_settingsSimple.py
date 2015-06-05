#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp

simpleSettingsDict={
# The number of samples orbital parameters to try/draw [int]
'nSamples' : (1000000,"Number of MCMC or MC samples"),
# Number of simulation chains to run in parallel, [1,100] [int].  
# NOTE: greater than numCores-1 causes system to slow down!
'nChains' : (6,"Number of chains"),
# Directory where you want the output data folder to go [string, at least 2 chars long]
'outDir' : '/run/media/kmede/Data1/Todai_Work/Data/data_SMODT',
# General filename for the simulation output folder to distinguish between simulation runs [string, at least 2 chars long]
#*************************************************************************************************************************
'outRoot' : "SMODT2-HIP10321-RV-omegaStraight",
#*************************************************************************************************************************
# data mode, choices {'RV','DI','3D'} [string]
'dataMode' : ('RV',"Data Mode (RV,DI,3D)"),
# mode to run simulation in, choices {'MC','SA','MCMC'} [string]
'symMode' : ('MC',"Simulator mode (MC,SA,MCMC)"),
###################################################
# Ranges for acceptable random number inputs ######
###################################################
# Minimum/Maximum allowed value for the mass of the primary body [double][Msun]
'mass1MIN' : 0.5,
'mass1MAX' : 1.6,
# Minimum/Maximum allowed value for the mass of the secondary body [double][Msun]
'mass2MIN' : 0.08,
'mass2MAX' : 0.085,
# Minimum/Maximum allowed value for the system distance from Earth [double][PC]
'distMIN' : 26.0,
'distMAX' : 27.0,
# Minimum/Maximum allowed value for the Longitude of the Ascending Node [double][deg]
'OmegaMIN' : 1.0,
'OmegaMAX' : 180.0,
# Minimum/Maximum allowed value for the Eccentricity [double]
'eMIN' : 0.33,
'eMAX' : 0.36,
# Minimum/Maximum value for the Time of Last Periapsis (or Time of Center Transit) [JD]
#(-1 indicates to use [earliestsEpoch-period,earliestEpoch])
'TMIN' : 2440000,
'TMAX' : 2445900,
# Minimum/Maximum allowed value for the Period [double][yrs]
'PMIN' : 19.0,
'PMAX' : 21.0,
# Minimum/Maximum allowed value for the Inclination [double][deg]
'incMIN' : 1.0,
'incMAX' : 89.0,
# Minimum/Maximum allowed value for the Argument of Perigee [double][deg]
'omegaMIN' : 150,
'omegaMAX' : 174.0,
# Minimum/Maximum value for Semi-major amplitude of RV curve [m/s]
'KMIN' : 840,
'KMAX' : 850,
# Minimum/Maximum values of Radial Velocity Offsets.  
# Must be one per set of RV data in same order as data comes in RVdata.dat, or the a single value to be used by all [comma separated list of doubles]
'vMINs' :[4000,-700,5000],
'vMAXs' :[8000,1000,8000],
}


