#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp

simpleSettingsDict={
# The number of samples orbital parameters to try/draw [int]
'nSamples' : (100000,"Number of MCMC or MC samples"),
# Number of simulation chains to run in parallel, [1,100] [int].  
# NOTE: greater than numCores-1 causes system to slow down!
'nChains' : (7,"Number of chains"),
# Directory where you want the output data folder to go [string, at least 2 chars long]
'outDir' : '/run/media/kmede/Data1/Todai_Work/Data/data_SMODT',
# General filename for the simulation output folder to distinguish between simulation runs [string, at least 2 chars long]
#*************************************************************************************************************************
'outRoot' : "SMODT2-HIP10321-RV-omegaStraight-MCMC",
#*************************************************************************************************************************
# data mode, choices {'RV','DI','3D'} [string]
'dataMode' : ('RV',"Data Mode (RV,DI,3D)"),
# mode to run simulation in, choices {'MC','SA','MCMC'} [string]
'symMode' : ('MCMC',"Simulator mode (MC,SA,MCMC)"),
###################################################
# Ranges for acceptable random number inputs ######
###################################################
# Minimum/Maximum allowed value for the mass of the primary body [double][Msun]
'mass1MIN' : 0.01,
'mass1MAX' : 2.2,
# Minimum/Maximum allowed value for the mass of the secondary body [double][Msun]
'mass2MIN' : 0.07,
'mass2MAX' : 0.09,
# Minimum/Maximum allowed value for the system distance from Earth [double][PC]
'distMIN' : 26.4,
'distMAX' : 27.2,
# Minimum/Maximum allowed value for the Longitude of the Ascending Node [double][deg]
'OmegaMIN' : 1.0,
'OmegaMAX' : 180.0,
# Minimum/Maximum allowed value for the Eccentricity [double]
'eMIN' : 0.34,
'eMAX' : 0.38,
# Minimum/Maximum value for the Time of Last Periapsis (or Time of Center Transit) [JD]
#(-1 indicates to use [earliestsEpoch-period,earliestEpoch])
'TMIN' : 2450000,
'TMAX' : 2455900,
# Minimum/Maximum allowed value for the Period [double][yrs]
'PMIN' : 19.0,
'PMAX' : 21.0,
# Minimum/Maximum allowed value for the Inclination [double][deg]
'incMIN' : 1.0,
'incMAX' : 89.0,
# Minimum/Maximum allowed value for the Argument of Perigee [double][deg]
'omegaMIN' : 166,
'omegaMAX' : 174.0,
# Minimum/Maximum value for Semi-major amplitude of RV curve [m/s]
'KMIN' : 835,
'KMAX' : 870,
# Minimum/Maximum values of Radial Velocity Offsets.  
# Must be one per set of RV data in same order as data comes in RVdata.dat, or the a single value to be used by all [comma separated list of doubles]
'vMINs' :[6100,300,6200],
'vMAXs' :[6200,400,6400],
}


