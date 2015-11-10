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
    
def mass1PriorRatio(MProposed,MLast):
    if (simpleSettingsDict['mass1MAX']!=0)and True:
        if MProposed!=MLast!=0:
            prop = chabrierPrior(MProposed,IMF=False)
            lst = chabrierPrior(MLast,IMF=False)
            return prop/lst
        else:
            return 1.0
    else:
        return 1.0
    
def mass2PriorRatio(MProposed,MLast):
    if (simpleSettingsDict['mass2MAX']!=0)and True:
        if MProposed!=MLast!=0:
            prop = chabrierPrior(MProposed,IMF=False)
            lst = chabrierPrior(MLast,IMF=False)
            return prop/lst
        else:
            return 1.0
    else:
        return 1.0
    
def paraPriorRatio(paraProposed,paraLast):
    if paraProposed!=paraLast!=simpleSettingsDict['paraMAX']!=0:
        ratioA = (paraLast**4.0)/(paraProposed**4.0)
        ratioB = 1.0
        if advancedDict['paraEst'][0]!=0:
            ## a Gaussian prior centered on hipparcos and width of hipparcos estimated error
            top = gaussian(paraProposed, advancedDict['paraEst'][0], advancedDict['paraErr'][0])
            btm = gaussian(paraLast, advancedDict['paraEst'][0], advancedDict['paraErr'][0])
            ratioB = top/btm
        return ratioA*ratioB
    else:
        return 1.0
    
def chabrierPrior(m,IMF=False):
    if IMF==False:
        if m<1.0:
            d = (0.068618528140713786/m)*np.exp((-(np.log10(m)+1.1023729087095586)**2)/0.9521999999999998)
        elif m<3.47:
            d = 0.019108957203743077*(m**(-5.37))
        elif m<18.20:
            d = 0.0065144172285487769*(m**(-4.53))
        else:
            d = 0.00010857362047581295*(m**(-3.11))
    else:
        if m<1.0:
            d = (0.068618528140713786/m)*np.exp((-(np.log10(m)+1.1023729087095586)**2)/0.9521999999999998)
        else:
            d = 0.019239245548314052*(m**(-2.3))
    return d


advancedDict = {
#NOTE: key max = 8characters, value+comment max = 68 characters, comment Max=47 it seems in testing.
########################
### General Settings ###
########################
# This will set the maximum reduced ChiSquared value to accept and write to the output file during MC mode. [double]
'chiMAX' : (300.0,"Max reduced chiSquared during MC"),
# maximum allowed reduced chiSquared out of SA before entering ST [double]
'chiMaxST':(5,'Max reduced chiSquared to enter ST.'),
# maximum allowed reduced chiSquared out of ST before entering MCMC [double]
'cMaxMCMC':(5,'Max reduced chiSquared to enter MCMC.'),
# set level of log messages to screen [int],recommend 50, ignoring critical msgs can cause problems. 
# choices: ('NONE'=100,'CRITICAL'=50,'ERROR'=40,'WARNING'=30,'INFO'=20,'DEBUG'10,'ALL'=0)
'logLevel' : 10,
#number of times to produce a summary log msg during a stage's progress [int]
'nSumry'  :20,
# make plot of posterior distributions? [bool]
'pltDists' :True,
# make plots of RV and DI/AM orbit fits [bool]
'pltOrbit' :True,
# Delete chain files after simulation is complete? [bool]
'delChains' :True,
# Delete combined data files after simulation is complete? [bool]
'delCombined' :False,
# run 'make' on C++/SWIG code to make sure it is up-to-date [bool]
'remake' :False,
###$$$$$$$$$$$$$$$$$$$$$$ Keep in final version $$$$$$$$$$$$$$$$$$$$$$$$$$
# Copy output non-data files to a Dropbox folder? [bool]  
'CopyToDB' :False,
'dbFolder' : '/run/media/kmede/HOME/Dropbox/SMODT-outputCopies/',
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############################
# Settings for MCMC mode ###
############################
# Calculate the length of the burn in for each chain (must be more than 1 chain)? [bool] 
'CalcBurn' :True,
# remove burn-in of output MCMC chains before combining (must be more than 1 chain) (should already be handled by SimAnneal stage2 though...)?
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
# Operate in low eccenctricity mode? [bool]
# Then step through in sqrt(e)sin(omega) and sqrt(e)cos(omega) instead of e & omega directly
'lowEcc'   : (True,"low eccentricty stepping?"),
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
# force adding a value in degrees to argument of periapsis used in DI orbit fit [double]
'omegaPdi' : (0.0,"Custom fixed val added to DI omega in model"),
# Is the data in the DIdata.dat in SA,PA format? else, it is in E,N (x,y) format [bool]
'sapa'     : (False,"Is astrometry data in SA,PA format?"),
######################
# System Information #
######################
#best estimate of primary's mass, and error [double][Msun]
'mass1Est' : (0.0,"Primary's estimated mass"),
'mass1Err' : (0.0,"Primary's estimated mass error"),
#best estimate of secondary's mass, and error [double][Msun]
'mass2Est' : (0.0,"Secondary's estimated mass"),
'mass2Err' : (0.0,"Secondary's estimated mass error"),
#best estimate of parallax, and error [double][mas]
'paraEst' : (50,"Estimated parallax"),
'paraErr' : (2.5,"Estimated parallax error"),
##################################
# Push prior functions into dict #
##################################
'ePrior'    :(True,'Use prior for eccentricity?',ePriorRatio),
'pPrior'    :(True,'Use prior for period?',pPriorRatio),
'incPrior'  :(True,'Use prior for inclination?',incPriorRatio),
'M1Prior':(False,'Use prior for M1?',mass1PriorRatio),
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
