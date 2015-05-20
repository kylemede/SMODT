import numpy as np
import symSettingsSimple
import constants

advancedDict = {
########################
### General Settings ###
########################
# This will set the maximum reduced ChiSquared value to accept and write to the output file.  ONLY for Monte Carlo and Simulated Annealing, not MCMC!! [double]
'chiMAX' : 3000.0,
# set to 'true' to have NOTHING print to screen while running [bool]
'SILENT' : False,
# set to 'false' to receive extra prints from the main simulations progress for testing [bool]
'quiet' : True,
# set to 'true' to receive prints from the functions called by main for testing [bool]
'verbose' : False,
# Simulated Annealing starting temperature [double]
'strtTemp' : 800.0,
# make plot of posterior distributions? [bool]
'pltDists' :True,
# make plots of RV and DI/AM orbit fits [bool]
'pltOrbit' :True,
# Delete chain files after simulation is complete? [bool]
'delChains' :True,
# Delete combined data files after simulation is complete? [bool]
'delCombined' :False,
###$$$$$$$$$$$$$$$$$$$$$$ Keep in final version $$$$$$$$$$$$$$$$$$$$$$$$$$
# Copy output non-data files to a Dropbox folder? [bool]  
'CopyToDB' :False,
'dbFolder' : '/run/media/kmede/Data1/Todai_Work/Dropbox/SMODT-outputCopies/',
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############################
# Settings for MCMC mode ###
############################
# Calculate the length of the burn in for each chain (must be more than 1 chain)? [bool] 
'CalcBurn' :True,
# remove burn-in of output MCMC chains before combining (must be more than 1 chain) (should already be handled by SimAnneal stage2 though...)?
'delBurn' : True,
# Calculate the Correlation lengths and number of effective points of each chain (must be more than 1 chain)? [bool]
'calcCL' :True,
# number of samples to draw for simulated annealing stage [int] 
'nSAsamp' :1000000,
# number of samples to draw for sigma tuning stage [int] 
'nSTsamp' :500000,
# Make plots of MCMC progress plots? [bool]
'pltMCMCprog' :False,
# Make plots of Simulated Annealing progress plots? [bool]
'pltSAprog' :False,
# Calculate the Gelman-Rubin statistic? [bool]
'CalcGR' :True,
# How many times do you want the Gelman-Rubin statistic calculated [int]
'nGRcalc' :10,
#####################################
# Special Settings for the models ###
#####################################
# fit to the primary's RV orbit [bool]
'fitPrime' : False,
# Are the RVs in the RVdata.dat for the Primary star? [bool]
'primeRVs' : True,
# Step through parameter space in Time of Center Transit (Inferior Conjunction)?  [bool]
'TcStep' : False,
# take the time of center transit (inferior conjunction) into account? [bool]
'TcEqualT' : True,
# force adding a value in degrees to argument of periapsis used in RV orbit fit [double]
'omegaPrv' : 0.0,
##################################################
## Special settings DI model:
# force adding a value in degrees to argument of periapsis used in RV orbit fit [double]
'omegaPdi' : 0.0,
######################
# System Information #
######################
#best estimate of primary's mass, and error [double][Msun]
'mass1Est' : 1.0,
'mass1Err' : 0.1,
#best estimate of secondary's mass, and error [double][Msun]
'mass2Est' : 0.2,
'mass2Err' : 0.1,
#best estimate of system's distance from Earth, and error [double][PC]
'distEst' : 5.0,
'distErr' :0.1,

##################################
# Push prior functions into dict #
##################################
'ePrior' :ePrior,
'pPrior' :pPrior,
'incPrior' :incPrior,
'mass1Prior' :mass1Prior,
'mass2Prior' :mass2Prior,
'distPrior' :distPrior,
}

########################################
#Define the priors as python functions #
########################################
#NOTE: only change the code and not the name of the functions or their inputs.
def ePrior(e,P):
    if (P*constants.daysPerYear<1000.0)and(symSettingsSimple.eMAX!=0):
        return 2.0*e
    else:
        return 1.0
def pPrior(P):
    if symSettingsSimple.PMAX!=symSettingsSimple.PMIN!=0:
        return P
    else:
        return 1.0
def incPrior(inc):
    if symSettingsSimple.incMAX!=symSettingsSimple.incMIN!=0:
        return np.sin(inc*(constants.pi/180.0))
    else:
        return 1.0
def mass1Prior(mass):
    if symSettingsSimple.mass1MIN!=symSettingsSimple.mass1MAX!=0:
        return gaussian(mass, mass1Est, mass1Err)
    else:
        return 1.0
def mass2Prior(mass):
    if symSettingsSimple.mass2MIN!=symSettingsSimple.mass2MAX!=0:
        return gaussian(mass, mass2Est, mass2Err)
    else:
        return 1.0
def distPrior(dist):
    if symSettingsSimple.distMIN!=symSettingsSimple.distMAX!=0:
        return gaussian(dist, distEst, distErr)
    else:
        return 1.0

def gaussian(x,mu,sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
