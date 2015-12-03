#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp

simpleSettingsDict={
# The number of samples orbital parameters to try/draw [int]
'nSamples' : (500000,"Number of MCMC or MC samples"),
# Number of simulation chains to run in parallel, [1,100] [int].  
# NOTE: greater than numCores-1 causes system to slow down!
# For MCMC mode this is the number of SA and ST chains.
'nChains' : (7,"Number MC/SA/ST of chains"),
# Number of MCMC chains to run in parallel. ONLY available in 'MCMC' mode. [1,100] [int].  
# NOTE: greater than numCores-1 causes system to slow down!
'nMCMCcns' : (7,"Number MCMC of chains"),
# Directory where you want the output data folder to go [string, at least 2 chars long]
'outDir' : '/run/media/kmede/Data1/Todai_Work/Data/data_SMODT',
# Directory on an SSD for fast reading/writting of temp data files to save on RAM usage.
##NOTE: just set to same as outDir if you don't have an SSD.
'tmpDir': '/run/media/kmede/HOME/tmp',
# General filename for the simulation output folder to distinguish between simulation runs [string, at least 2 chars long]
#*************************************************************************************************************************
'outRoot' : "SMODT2-HD130948BC-DI-MCMC-startAtBest-fixedDIdata-corrLengthTest1",
#*************************************************************************************************************************
# data mode, choices {'RV','DI','3D'} [string]
'dataMode' : ('DI',"Data Mode (RV,DI,3D)"),
# mode to run simulation in, choices {'MC','SA','SAST','SASTMCMC,'MCMC'} [string]
'symMode' : ('SASTMCMC',"Simulator mode (MC,SA,SAST,SASTMCMC,MCMC)"),
############################################
# Starting parameters and sigmas for MCMC  #
# Can be found with prior run in SAST mode #
############################################
# if unknown, set to False!! else [comma separated list of doubles]
'startParams' : False,
# if unknown, set to False!! else [comma separated list of doubles]
'startSigmas' : False,
###################################################
# Ranges for acceptable random number inputs ######
###################################################
# Minimum/Maximum allowed value for the mass of the primary body [double][Msun]
# NOTE: For DI only cases, use mass1 values as total mass and set mass2 values to zero.
'mass1MIN' : 0.095,
'mass1MAX' : 0.12,
# Minimum/Maximum allowed value for the mass of the secondary body [double][Msun]
'mass2MIN' : 0.0,
'mass2MAX' : 0.0,
# Minimum/Maximum allowed value for the Parallax [double][mas]
'paraMIN' : 52.00,
'paraMAX' : 58.00,
# Minimum/Maximum allowed value for the Longitude of the Ascending Node [double][deg]
'OmegaMIN' : 310.0,
'OmegaMAX' : 316.0,
# Minimum/Maximum allowed value for the Eccentricity, allowed range [0,0.98]. [double]
'eMIN' : 0.012,
'eMAX' : 0.21,
# Minimum/Maximum value for the Time of Last Periapsis (or Time of Center Transit) [JD]
#(-1 indicates to use [earliestsEpoch-period,earliestEpoch])
'TMIN' : 2454570,
'TMAX' : 2454720,
# Minimum/Maximum allowed value for the Period [double][yrs]
'PMIN' : 9.5,
'PMAX' : 11.0,
# Minimum/Maximum allowed value for the Inclination [double][deg]
'incMIN' : 95.0,
'incMAX' : 97.0,
# Minimum/Maximum allowed value for the Argument of Perigee [double][deg]
'omegaMIN' : 235.0,
'omegaMAX' : 255.0,
# Minimum/Maximum value for Semi-major amplitude of RV curve [m/s]
'KMIN' : 0,
'KMAX' : 0,
# Minimum/Maximum values of Radial Velocity Offsets.  
# Must be one per set of RV data in same order as data comes in RVdata.dat, or the a single value to be used by all [comma separated list of doubles]
'vMINs' :[],
'vMAXs' :[],
}


