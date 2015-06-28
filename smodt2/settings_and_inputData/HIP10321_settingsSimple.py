#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp

simpleSettingsDict={
# The number of samples orbital parameters to try/draw [int]
'nSamples' : (20000000,"Number of MCMC or MC samples"),
# Number of simulation chains to run in parallel, [1,100] [int].  
# NOTE: greater than numCores-1 causes system to slow down!
'nChains' : (7,"Number of chains"),
# Directory where you want the output data folder to go [string, at least 2 chars long]
'outDir' : '/run/media/kmede/Data1/Todai_Work/Data/data_SMODT',
# General filename for the simulation output folder to distinguish between simulation runs [string, at least 2 chars long]
#*************************************************************************************************************************
'outRoot' : "SMODT2-HIP10321-3D-retro-MCMC-lowerDIerrors-M1equal1",
#*************************************************************************************************************************
# data mode, choices {'RV','DI','3D'} [string]
'dataMode' : ('3D',"Data Mode (RV,DI,3D)"),
# mode to run simulation in, choices {'MC','SA','MCMC'} [string]
'symMode' : ('MCMC',"Simulator mode (MC,SA,MCMC)"),
###################################################
# Ranges for acceptable random number inputs ######
###################################################
# Minimum/Maximum allowed value for the mass of the primary body [double][Msun]
'mass1MIN' : 1.00,
'mass1MAX' : 0,
# Minimum/Maximum allowed value for the mass of the secondary body [double][Msun]
'mass2MIN' : 0.14,
'mass2MAX' : 0.55,
# Minimum/Maximum allowed value for the Parallax [double][mas]
'paraMIN' : 35.5,
'paraMAX' : 40.00,
# Minimum/Maximum allowed value for the Longitude of the Ascending Node [double][deg]
'OmegaMIN' : 180.0,
'OmegaMAX' : 360.0,
# Minimum/Maximum allowed value for the Eccentricity, allowed range [0,0.98]. [double]
'eMIN' : 0.32,
'eMAX' : 0.42,
# Minimum/Maximum value for the Time of Last Periapsis (or Time of Center Transit) [JD]
#(-1 indicates to use [earliestsEpoch-period,earliestEpoch])
'TMIN' : 2452200,
'TMAX' : 2452500,
# Minimum/Maximum allowed value for the Period [double][yrs]
'PMIN' : 19.5,
'PMAX' : 22.5,
# Minimum/Maximum allowed value for the Inclination [double][deg]
'incMIN' : 140.0,
'incMAX' : 169.0,
# Minimum/Maximum allowed value for the Argument of Perigee [double][deg]
'omegaMIN' : 345.0,
'omegaMAX' : 360.0,
# Minimum/Maximum value for Semi-major amplitude of RV curve [m/s]
'KMIN' : 800,
'KMAX' : 900,
# Minimum/Maximum values of Radial Velocity Offsets.  
# Must be one per set of RV data in same order as data comes in RVdata.dat, or the a single value to be used by all [comma separated list of doubles]
'vMINs' :[6100,330,6250],
'vMAXs' :[6250,440,6400],
}