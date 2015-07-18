#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import numpy as np
from settingsSimple import simpleSettingsDict
import constants

########################################
#Define the priors as python functions #
########################################
#NOTE: key max = 8characters, value+comment max = 68 characters, comment Max=47 it seems in testing.
#NOTE: only change the code and not the name of the functions or their inputs.
def ePrior(e,P):
    if advancedDict['lowEcc'][0]==False:
        if (e!=0)and(P!=0):
            if ((P*constants.daysPerYear)>1000.0)and(simpleSettingsDict['eMAX']!=0):
                return 2.0*e
            else:
                return 1.0
    else:
        return 1.0
def pPrior(P):
    if P!=0:
        if simpleSettingsDict['PMAX']!=simpleSettingsDict['PMIN']!=0:
            return P
        else:
            return 1.0
    else:
        return 1.0
def incPrior(inc):
    if inc!=0:
        if simpleSettingsDict['incMAX']!=simpleSettingsDict['incMIN']!=0:
            return np.sin(inc*(constants.pi/180.0))
        else:
            return 1.0
    else:
        return 1.0
def mass1Prior(mass):
    return 1.0
#     if mass!=0:
#         if simpleSettingsDict['mass1MIN']!=simpleSettingsDict['mass1MAX']!=0:
#             return gaussian(mass, advancedDict['mass1Est'][0], advancedDict['mass1Err'][0])
#         else:
#             return 1.0
#     else:
#         return 1.0
def mass2Prior(mass):
    return 1.0
#     if simpleSettingsDict['mass2MIN']!=simpleSettingsDict['mass2MAX']!=0:
#         return gaussian(mass, advancedDict['mass2Est'][0], advancedDict['mass2Err'][0])
#     else:
#         return 1.0
def paraPrior(parallax):
    if parallax!=0:
        if simpleSettingsDict['paraMIN']!=simpleSettingsDict['paraMAX']!=0:
            return gaussian(parallax, advancedDict['paraEst'][0], advancedDict['paraErr'][0])
        else:
            return 1.0
    else:
        return 1.0


advancedDict = {
########################
### General Settings ###
########################
# This will set the maximum reduced ChiSquared value to accept and write to the output file during MC mode. [double]
'chiMAX' : (300.0,"Max reduced chiSquared during MC and SA"),
# maximum allowed best reduced chiSquared out of SA before entering ST [double]
'chiMaxST':(5,'Max reduced chiSquared to enter ST.'),
# maximum allowed best reduced chiSquared out of ST before entering MCMC [double]
'cMaxMCMC':(3.0,'Max reduced chiSquared to enter MCMC.'),
# set level of log messages to screen [int],recommend 50, ignoring critical msgs can cause problems. 
# choices: ('NONE'=100,'CRITICAL'=50,'ERROR'=40,'WARNING'=30,'INFO'=20,'DEBUG'10,'ALL'=0)
'logLevel' : 30,
#number of times to produce a summary log msg during a stage's progress [int]
'nSumry'  :20,
# make plot of posterior distributions? [bool]
'pltDists' : True,
# make plots of RV and DI/AM orbit fits [bool]
'pltOrbit' :True,
# Delete chain files after simulation is complete? [bool]
'delChains' :True,
# Delete combined data files after simulation is complete? [bool]
'delCombined' :False,
# run 'make' on C++/SWIG code to make sure it is up-to-date [bool]
'remake' :False,
###$$$$$$$$$$$$$$$$$$$$$$ Keep in final version? $$$$$$$$$$$$$$$$$$$$$$$$$$
# Copy output non-data files to a Dropbox folder? [bool]  $$$$$ still not coded up $$$
'CopyToDB' :False,
'dbFolder' : '/run/media/kmede/HOME/Dropbox/SMODT-outputCopies/',
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############################
# Settings for MCMC mode ###
############################
# Calculate the length of the burn in for each chain (must be more than 1 chain)? [bool] 
'CalcBurn' :True,
# remove burn-in of output MCMC chains before combining (must be more than 1 chain) (should already be handled by ST though...)? [bool]
'rmBurn' : (True,"Remove Burn-in?"),
# Calculate the Correlation lengths and number of effective points of each chain (must be more than 1 chain)? [bool]
# NOTE: CAUTION, can take a long time for long runs.  Still needs to be sped up somehow.
'calcCL' :False,
# number of samples to draw for simulated annealing stage [int] 
'nSAsamp' :(1000000,"Num SA samples"),
# Simulated Annealing starting temperature [double]
'strtTemp' : (300.0,"SA start temp."),
# Starting sigma size, % of parameter range, recommend [0.05,0.25].  [double]
# After first trial of SA and ST, take ST output and use here.
'strtSig' : (0.15,"start percent param range for SA"),
# Number of samples till temperature drop. [int]
# Allowed vals [1,nSAsamp), Ideal is ~50.
'tempInt'  : (50,"Num steps till temp drops in SA."),
# number of samples to draw for sigma tuning stage [int].
'nSTsamp' :(100000,"Num ST samples"),
# number of steps per varying parameter until calculating the acceptance rate and tuning sigmas. [int]
# Allowed vals [1,nSTsamp), testing shows a value of ~200 works well.
'sigInt': (200,"Num steps/par till calc acc rate/tune sigs."),
# Maximum step size allowed, as a ratio of each parameters range ie. 1.0=100% [double]
'sigMax' :(1.0,'Max ratio of params range,for step size.'),
# Minimum step size allowed, as a ratio of each parameters range ie. 1.0=100% [double]
'sigMin' :(0.02,'Min ratio of params range,for step size.'),
# interval of accepted values between storing in output array (for SA,ST,MCMC, not MC) [int]
# Make sure to save enough that R~1.0 at max, posteriors look smooth, BUT not too much data is saved that you are just wasting disk space.
'saveInt' : (10,"Int between saving params, for all but MC."),
# Interval of saved values before write/dump the data to disk to avoid consuming too much RAM during long runs. They take 11MB/100000.
'dmpInt'   : 100000,
# Start MCMC at the best params from the ST stage? [bool]
'strBest' : (True,"Start MCMC at best fit params from ST"),
## NOTE: progress plots have no code yet, so MUST be False!!!
# Make plots of MCMC progress plots? [bool]$$$$$ still not coded up $$$
'pltMCMCprog' :False,
# Make plots of Simulated Annealing progress plots? [bool]$$$$$ still not coded up $$$
'pltSAprog' :False,
# Calculate the Gelman-Rubin statistic? [bool]
'CalcGR' :True,
# How many times do you want the Gelman-Rubin statistic calculated [int]  $$$$$ still not coded up $$$
'nGRcalc' :10,
#####################################
# Special Settings for the models ###
#####################################
# Operate in low eccenctricity mode? [bool]
# Then step through in sqrt(e)sin(omega) and sqrt(e)cos(omega) instead of e & omega directly
'lowEcc'   : (False,"low eccentricty stepping?"),
# fit to the primary's RV orbit [bool]
'fitPrime' : (False,"Fit primary's orbit?"),
# Are the RVs in the RVdata.dat for the Primary star? [bool]
'primeRVs' : (True,"RVs measured from Primary?"),
# Draw values for K directly, do NOT calculate it [bool]. Kills varying of Inclination.  Only possible in RV only mode.
'Kdirect'  : (False,'Vary K direct, do not calc it'),
# Step through parameter space in Time of Center Transit (Inferior Conjunction)?  [bool]
'TcStep' : (False,"Step in Tc not T?"),
# take the time of center transit (inferior conjunction) into account? [bool]
'TcEqualT' : (True,"Fix Tc=T?"),
# force adding a value in degrees to argument of periapsis used in RV orbit fit [double]
'omegaPrv' : (0.0,"Custom fixed val added to RV omega in model"),
##################################################
## Special settings DI model:
# force adding a value in degrees to argument of periapsis used in RV orbit fit [double]
'omegaPdi' : (0.0,"Custom fixed val added to DI omega in model"),
######################
# System Information #
######################
#best estimate of primary's mass, and error [double][Msun]
'mass1Est' : (1.09,"Primary's estimated mass"),
'mass1Err' : (0.1,"Primary's estimated mass error"),
#best estimate of secondary's mass, and error [double][Msun]
'mass2Est' : (0.0818,"Secondary's estimated mass"),
'mass2Err' : (0.002,"Secondary's estimated mass error"),
#best estimate of parallax, and error [double][mas]
'paraEst' : (37.25,"Estimated parallax"),
'paraErr' : (0.55,"Estimated parallax error"),
##################################
# Push prior functions into dict #
##################################
'ePrior'    :(True,'Use prior for eccentricity?',ePrior),
'pPrior'    :(True,'Use prior for period?',pPrior),
'incPrior'  :(True,'Use prior for inclination?',incPrior),
'M1Prior':(False,'Use prior for M1?',mass1Prior),
'M2Prior':(False,'Use prior for M2?',mass2Prior),
'parPrior' :(True,'Use prior for parallax?',paraPrior),
}


def gaussian(x,mu,sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

######################
# Merge the two dicts#
######################
settingsDict = {}
for key in simpleSettingsDict:
    settingsDict[key]=simpleSettingsDict[key]
for key in advancedDict:
    settingsDict[key]=advancedDict[key]