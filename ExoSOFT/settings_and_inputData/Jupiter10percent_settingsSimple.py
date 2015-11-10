#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp

simpleSettingsDict={
# The number of samples orbital parameters to try/draw [int]
'nSamples' : (90000000,"Number of MCMC or MC samples"),
# Number of simulation chains to run in parallel, [1,100] [int].  
# NOTE: greater than numCores-1 causes system to slow down!
'nChains' : (7,"Number of chains"),
# Directory where you want the output data folder to go [string, at least 2 chars long]
'outDir' : '/run/media/kmede/Data1/Todai_Work/Data/data_SMODT',
# Directory on an SSD for fast reading/writting of temp data files to save on RAM usage.
##NOTE: just set to same as outDir if you don't have an SSD.
'tmpDir': '/run/media/kmede/HOME/tmp',
# General filename for the simulation output folder to distinguish between simulation runs [string, at least 2 chars long]
#*************************************************************************************************************************
'outRoot' : "JUPITER2-3D-MCMC-10percent-lowEcc-PDMFm1m2-newParaPrior-Long",
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
'mass1MIN' : 0.4,
'mass1MAX' : 2.35,
# Minimum/Maximum allowed value for the mass of the secondary body [double][Msun]
'mass2MIN' : 0.00045,
'mass2MAX' : 0.002,
# Minimum/Maximum allowed value for the Parallax [double][mas]
'paraMIN' : 37.50,
'paraMAX' : 63.50,
# Minimum/Maximum allowed value for the Longitude of the Ascending Node [double][deg]
'OmegaMIN' : 85.0,
'OmegaMAX' : 115.0,
# Minimum/Maximum allowed value for the Eccentricity [double]
'eMIN' : 0.000005,
'eMAX' : 0.18,
# Minimum/Maximum value for the Time of Last Periapsis (or Time of Center Transit) [JD]
#(-1 indicates to use [earliestsEpoch-period,earliestEpoch])
'TMIN' : 2449300,
'TMAX' : 2452800,
# Minimum/Maximum allowed value for the Period [double][yrs]
'PMIN' : 10,
'PMAX' : 14.5,
# Minimum/Maximum allowed value for the Inclination [double][deg]
'incMIN' : 20,
'incMAX' : 70.0,
# Minimum/Maximum allowed value for the Argument of Perigee [double][deg]
'omegaMIN' : -100,
'omegaMAX' : 190,
# Minimum/Maximum value for Semi-major amplitude of RV curve [m/s]
'KMIN' : 0,
'KMAX' : 0,
# Minimum/Maximum values of Radial Velocity Offsets.  
# Must be one per set of RV data in same order as data comes in RVdata.dat, or the a single value to be used by all [comma separated list of doubles]
'vMINs' :[-2.2],
'vMAXs' :[1.5],
}


