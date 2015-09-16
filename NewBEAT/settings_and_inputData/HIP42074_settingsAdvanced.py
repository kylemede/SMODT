#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import numpy as np
from settingsSimple import simpleSettingsDict
import constants

########################################
#Define the priors as python functions #
########################################
#NOTE: only change the code and not the name of the functions or their inputs.
def ePriorRatio(eProposed,eLast):
    if (advancedDict['lowEcc'][0]==False)and(simpleSettingsDict['eMAX']!=0):
        if eProposed!=eLast!=0:
            if (simpleSettingsDict['PMIN']*constants.daysPerYear)>1000.0:
                return eProposed/eLast
            else:
                return 1.0
        else:
            return 1.0
    else:
        return 1.0
def pPriorRatio(Pproposed,Plast):
    if simpleSettingsDict['PMAX']!=0:
        if Pproposed!=0:
            return Plast/Pproposed
        else:
            return 1.0
    else:
        return 1.0
def incPriorRatio(incProposed,incLast):
    if simpleSettingsDict['incMAX']!=0:
        if (incLast%90.0)!=0:
            return np.sin(incProposed*(constants.pi/180.0))/np.sin(incLast*(constants.pi/180.0))
        else:
            return 1.0
    else:
        return 1.0
    
def mass1PriorRatio(M1Proposed,M1Last):
    #we are assuming M1>70Mj
    if (simpleSettingsDict['mass1MAX']!=0)and True:
        if M1Proposed!=M1Last!=0:
            if True:#simpleSettingsDict['mass1MIN']>=0.5:
                #From Table 1 of Chabrier2003
                return (M1Proposed**(-2.3))/(M1Last**(-2.3))
            #elif (simpleSettingsDict['mass1MIN']>=0.07)and(simpleSettingsDict['mass1MAX']<=0.5):
            #    return (M1Proposed**(-1.3))/(M1Last**(-1.3))
            else:
                return 1.0
        else:
            return 1.0
    else:
        return 1.0
#     if mass!=0:
#         if simpleSettingsDict['mass1MIN']!=simpleSettingsDict['mass1MAX']!=0:
#             return gaussian(mass, advancedDict['mass1Est'][0], advancedDict['mass1Err'][0])
#         else:
#             return 1.0
#    else:
#         return 1.0
def mass2PriorRatio(M2Proposed,M2Last,aProposed,aLast):
    if (simpleSettingsDict['mass2MAX']!=0)and True:
        if M2Proposed!=M2Last!=0:
            if True:#simpleSettingsDict['mass2MIN']>=0.5:
                #From Table 1 of Chabrier2003
                return (M2Proposed**(-2.3))/(M2Last**(-2.3))
            #elif (simpleSettingsDict['mass2MIN']>=0.07)and(simpleSettingsDict['mass2MAX']<=0.5):
            #    return (M2Proposed**(-1.3))/(M2Last**(-1.3))
            #elif (simpleSettingsDict['mass2MIN']>=0.005)and(simpleSettingsDict['mass2MAX']<=0.07):
            #    if aProposed!=aLast!=0:
            #        return ((M2Proposed**(-0.65))*(aProposed**(-0.85)))/((M2Last**(-0.65))*(aLast**(-0.85)))
            #    else:
            #        return 1.0
            else:
                return 1.0
        else:
            return 1.0
    else:
        return 1.0
#     if simpleSettingsDict['mass2MIN']!=simpleSettingsDict['mass2MAX']!=0:
#         return gaussian(mass, advancedDict['mass2Est'][0], advancedDict['mass2Err'][0])
#     else:
#         return 1.0
def paraPriorRatio(paraProposed,paraLast):
    if paraProposed!=paraLast!=simpleSettingsDict['paraMAX']!=0:
        if False:
            return (paraLast**2.0)/(paraProposed**2.0)
        elif True:
            ## a Gaussian prior centered on hipparcos and width of hipparcos estimated error
            top = gaussian(paraProposed, advancedDict['paraEst'][0], advancedDict['paraErr'][0])
            btm = gaussian(paraLast, advancedDict['paraEst'][0], advancedDict['paraErr'][0])
            return top/btm
        else:
            return 1.0
    else:
        return 1.0


advancedDict = {
#NOTE: key max = 8characters, value+comment max = 68 characters, comment Max=47 it seems in testing.
########################
### General Settings ###
########################
# This will set the maximum reduced ChiSquared value to accept and write to the output file during MC mode. [double]
'chiMAX' : (300.0,"Max reduced chiSquared during MC and SA"),
# maximum allowed best reduced chiSquared out of SA before entering ST [double]
'chiMaxST':(30,'Max reduced chiSquared to enter ST.'),
# maximum allowed best reduced chiSquared out of ST before entering MCMC [double]
'cMaxMCMC':(20,'Max reduced chiSquared to enter MCMC.'),
# set level of log messages to screen [int],recommend 50, ignoring critical msgs can cause problems. 
# choices: ('NONE'=100,'CRITICAL'=50,'ERROR'=40,'WARNING'=30,'INFO'=20,'DEBUG'10,'ALL'=0)
'logLevel' : 10,
#number of times to produce a summary log msg during a stage's progress [int]
'nSumry'  :20,
# make plot of posterior distributions? [bool]
'pltDists' : True,
# make plots of RV and DI/AM orbit fits [bool]
'pltOrbit' :True,
# Delete chain files after simulation is complete? [bool]
'delChains' :True,
# Delete combined data files after simulation is complete? [bool]
'delCombined' :True,
# run 'make' on C++/SWIG code to make sure it is up-to-date [bool]
'remake' :False,
###$$$$$$$$$$$$$$$$$$$$$$ Keep in final version? $$$$$$$$$$$$$$$$$$$$$$$$$$
# Copy output non-data files to a Dropbox folder? [bool]  $$$$$ still not coded up $$$
'CopyToDB' :True,
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
'nSAsamp' :(1500000,"Num SA samples"),
# Simulated Annealing starting temperature [double]
'strtTemp' : (500.0,"SA start temp."),
# Starting sigma size, % of parameter range, recommend [0.05,0.25].  [double]
# After first trial of SA and ST, take ST output and use here.
'strtSig' : (0.15,"start percent param range for SA"),
# Number of samples till temperature drop. [int]
# Allowed vals [1,nSAsamp), Ideal is ~50.
'tempInt'  : (50,"Num steps till temp drops in SA."),
# number of samples to draw for sigma tuning stage [int].
'nSTsamp' :(500000,"Num ST samples"),
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
'dmpInt'   : 1000000,
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
'mass1Est' : (0.0,"Primary's estimated mass"),
'mass1Err' : (0.0,"Primary's estimated mass error"),
#best estimate of secondary's mass, and error [double][Msun]
'mass2Est' : (0.0,"Secondary's estimated mass"),
'mass2Err' : (0.00,"Secondary's estimated mass error"),
#best estimate of parallax, and error [double][mas]
'paraEst' : (47.21,"Estimated parallax"),
'paraErr' : (0.72,"Estimated parallax error"),
##################################
# Push prior functions into dict #
##################################
'ePrior'    :(True,'Use prior for eccentricity?',ePriorRatio),
'pPrior'    :(True,'Use prior for period?',pPriorRatio),
'incPrior'  :(True,'Use prior for inclination?',incPriorRatio),
'M1Prior':(True,'Use prior for M1?',mass1PriorRatio),
'M2Prior':(False,'Use prior for M2?',mass2PriorRatio),
'parPrior' :(True,'Use prior for parallax?',paraPriorRatio),
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
