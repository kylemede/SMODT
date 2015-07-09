#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp

simpleSettingsDict={
# The number of samples orbital parameters to try/draw [int]
'nSamples' : (20000000,"Number of MCMC or MC samples"),
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
'outRoot' : "SMODT2-HR7672-3D-MCMC-fixedM1M2Parallax-wideOpenP",
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
'mass1MIN' : 1.079,
'mass1MAX' : 1.081,
# Minimum/Maximum allowed value for the mass of the secondary body [double][Msun]
'mass2MIN' : 0.0655,
'mass2MAX' : 0.0657,
# Minimum/Maximum allowed value for the Parallax [double][mas]
'paraMIN' : 56.27,
'paraMAX' : 56.29,
# Minimum/Maximum allowed value for the Longitude of the Ascending Node [double][deg]
'OmegaMIN' : 145.0,
'OmegaMAX' : 160.0,
# Minimum/Maximum allowed value for the Eccentricity, allowed range [0,0.98]. [double]
'eMIN' : 0.4,
'eMAX' : 0.7,
# Minimum/Maximum value for the Time of Last Periapsis (or Time of Center Transit) [JD]
#(-1 indicates to use [earliestsEpoch-period,earliestEpoch])
'TMIN' : 2456400,
'TMAX' : 2457700,
# Minimum/Maximum allowed value for the Period [double][yrs]
'PMIN' : 30,
'PMAX' :150,
# Minimum/Maximum allowed value for the Inclination [double][deg]
'incMIN' : 90.0,
'incMAX' : 105.0,
# Minimum/Maximum allowed value for the Argument of Perigee [double][deg]
'omegaMIN' : 240.0,
'omegaMAX' : 290.0,
# Minimum/Maximum value for Semi-major amplitude of RV curve [m/s]
'KMIN' : 0,
'KMAX' : 0,
# Minimum/Maximum values of Radial Velocity Offsets.  
# Must be one per set of RV data in same order as data comes in RVdata.dat, or the a single value to be used by all [comma separated list of doubles]
'vMINs' :[200,200,200,200,200,200,200],
'vMAXs' :[550,550,550,550,550,550,550],
}


