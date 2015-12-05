#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp

simpleSettingsDict={
# The number of samples orbital parameters to try/draw [int]
'nSamples' : (40000000,"Number of MCMC or MC samples"),
# Number of simulation chains to run in parallel, [1,100] [int].  
# NOTE: greater than numCores-1 causes system to slow down!
# For MCMC mode this is the number of SA and ST chains.
'nChains' : (7,"Number MC/SA/ST of chains"),
# Number of MCMC chains to run in parallel. ONLY available in 'MCMC' mode. [1,100] [int].  
# NOTE: greater than numCores-1 causes system to slow down!
'nMCMCcns' : (7,"Number MCMC of chains"),
# Directory where you want the output data folder to go [string, at least 2 chars long]
'outDir' : '/run/media/kmede/Data1/Todai_Work/Data/data_SMODT',
# General filename for the simulation output folder to distinguish between simulation runs [string, at least 2 chars long]
#*************************************************************************************************************************
'outRoot' : "SMODT2-SyntheticJUPITER-3D-20percent-startAtBest-lowEccTrue",
#*************************************************************************************************************************
# set level of log messages to screen [int],recommend 50, ignoring critical msgs can cause problems. 
# choices: ('NONE'=100,'CRITICAL'=50,'ERROR'=40,'WARNING'=30,'INFO'=20,'DEBUG'10,'ALL'=0)
'logLevel' : 10,
# data mode, choices {'RV','DI','3D'} [string]
'dataMode' : ('3D',"Data Mode (RV,DI,3D)"),
# Run in Automatic mode? This will perform checks and select the stages to run automatically. [bool]
'autoMode' : (True, 'Run in Automatic mode?'),
# mode to run simulation in, choices {'MC','SA','ST','SAST','SASTMCMC,'MCMC'} [string]
# NOTE: 'ST' and 'MCMC' modes need a full list of parameters for startParams, else they fail!
#       'MCMC' also needs a full list of sigmas in startSigmas.
'stages' : 'MCMC',
############################################
# Starting parameters and sigmas for MCMC  #
# Can be found with prior run in SAST mode #
############################################
# if unknown, set to False!! else [comma separated list of doubles]
'startParams' : False,
# if unknown, set to False!! else [comma separated list of doubles]
'startSigmas' : False,
# If better startParams are found during ExoSOFT, push them and the newest sigmas into this file? [bool]
# NOTE: this can be helpful, but use caution if you do not wish to overwrite the values in here.
"pushToSettFiles":True,
}


