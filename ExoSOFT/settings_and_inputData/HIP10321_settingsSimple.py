#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp

simpleSettingsDict={
# The number of samples orbital parameters to try/draw [int]
'nSamples' : (5000000,"Number of MCMC or MC samples"),
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
#'outRoot' : "HIP10321-3D-retro-MCMC-all155msJitter-PDMFm1m2-fixedRanges4-SUPERlong",
'outRoot' : "HIP10321-3D-retro-MCMC-all14point6msJitter-1.05m1-PDMFm2-fixedRanges6-DI5point3masErr",
#'outRoot' : "HIP10321-3D--OPEN-MCMC-all14point6msJitter-PDMFm1m2-DI5point3masErr",
#*************************************************************************************************************************
# data mode, choices {'RV','DI','3D'} [string]
'dataMode' : ('3D',"Data Mode (RV,DI,3D)"),
# mode to run simulation in, choices {'MC','SA','MCMC'} [string]
'symMode' : ('MCMC',"Simulator mode (MC,SA,MCMC)"),
###################################################
# Ranges for acceptable random number inputs ######
###################################################
# Minimum/Maximum allowed value for the mass of the primary body [double][Msun]
'mass1MIN' : 1.05,
'mass1MAX' : 0,
# Minimum/Maximum allowed value for the mass of the secondary body [double][Msun]
'mass2MIN' : 0.1,
'mass2MAX' : 0.65,
# Minimum/Maximum allowed value for the Parallax [double][mas]
'paraMIN' : 30,
'paraMAX' : 50.0,
# Minimum/Maximum allowed value for the Longitude of the Ascending Node [double][deg]
'OmegaMIN' : 220.0,
'OmegaMAX' : 290.0,
# Minimum/Maximum allowed value for the Eccentricity, allowed range [0,0.98]. [double]
'eMIN' : 0.3,
'eMAX' : 0.5,
# Minimum/Maximum value for the Time of Last Periapsis (or Time of Center Transit) [JD]
#(-1 indicates to use [earliestsEpoch-period,earliestEpoch])
'TMIN' : 2444300,#2452200,
'TMAX' : 2445800,#2452500,
# Minimum/Maximum allowed value for the Period [double][yrs]
'PMIN' : 18,
'PMAX' : 27,
# Minimum/Maximum allowed value for the Inclination [double][deg]
'incMIN' : 120.0,
'incMAX' : 190.0,
# Minimum/Maximum allowed value for the Argument of Perigee [double][deg]
'omegaMIN' : 330.0,
'omegaMAX' : 380.0,
# Minimum/Maximum value for Semi-major amplitude of RV curve [m/s]
'KMIN' : 800,
'KMAX' : 900,
# Minimum/Maximum values of Radial Velocity Offsets.  
# Must be one per set of RV data in same order as data comes in RVdata.dat, or the a single value to be used by all [comma separated list of doubles]
'vMINs' :[6100,310,6220],
'vMAXs' :[6390,490,6520],
}
