#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp

simpleSettingsDict={
# The number of samples orbital parameters to try/draw [int]
'nSamples' : (1000000,"Number of MCMC or MC samples"),
# Number of simulation chains to run in parallel, [1,100] [int].  
# NOTE: greater than numCores-1 causes system to slow down!
'nChains' : (3,"Number of chains"),
# Directory where you want the output data folder to go [string, at least 2 chars long]
'outDir' : '/run/media/kmede/Data1/Todai_Work/Data/data_SMODT',
# Directory on an SSD for fast reading/writting of temp data files to save on RAM usage.
##NOTE: just set to same as outDir if you don't have an SSD.
'tmpDir': '/run/media/kmede/HOME/tmp',
# General filename for the simulation output folder to distinguish between simulation runs [string, at least 2 chars long]
#*************************************************************************************************************************
'outRoot' : "SMODT2-JUPITER2-3D-10percent-startAtBest-lowEccTrue",
#*************************************************************************************************************************
# data mode, choices {'RV','DI','3D'} [string]
'dataMode' : ('3D',"Data Mode (RV,DI,3D)"),
# mode to run simulation in, choices {'MC','SA','MCMC'} [string]
'symMode' : ('MCMC',"Simulator mode (MC,SA,MCMC)"),
###################################################
# Ranges for acceptable random number inputs ######
###################################################
# Minimum/Maximum allowed value for the mass of the primary body [double][Msun]
# NOTE: For DI only cases, use mass1 values as total mass and set mass2 values to zero.
'mass1MIN' : 0.1,
'mass1MAX' : 3.0,
# Minimum/Maximum allowed value for the mass of the secondary body [double][Msun]
'mass2MIN' : 0.0001,
'mass2MAX' : 0.0030,
# Minimum/Maximum allowed value for the Parallax [double][mas]
'paraMIN' : 40.00,
'paraMAX' : 60.00,
# Minimum/Maximum allowed value for the Longitude of the Ascending Node [double][deg]
'OmegaMIN' : 1.0,
'OmegaMAX' : 180.0,
# Minimum/Maximum allowed value for the Eccentricity [double]
'eMIN' : 0.0001,
'eMAX' : 0.2,
# Minimum/Maximum value for the Time of Last Periapsis (or Time of Center Transit) [JD]
#(-1 indicates to use [earliestsEpoch-period,earliestEpoch])
'TMIN' : 2450000,
'TMAX' : 2451100,
# Minimum/Maximum allowed value for the Period [double][yrs]
'PMIN' : 0.01,
'PMAX' : 20.0,
# Minimum/Maximum allowed value for the Inclination [double][deg]
'incMIN' : 1,
'incMAX' : 80.0,
# Minimum/Maximum allowed value for the Argument of Perigee [double][deg]
'omegaMIN' : -20,
'omegaMAX' : 70,
# Minimum/Maximum value for Semi-major amplitude of RV curve [m/s]
'KMIN' : 0,
'KMAX' : 0,
# Minimum/Maximum values of Radial Velocity Offsets.  
# Must be one per set of RV data in same order as data comes in RVdata.dat, or the a single value to be used by all [comma separated list of doubles]
'vMINs' :[-50],
'vMAXs' :[50],
}


