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
'chiMAX' : (300.0,"Max reduced chiSquared during MC and SA"),
# maximum allowed best reduced chiSquared out of SA before entering ST [double]
'chiMaxST':(15,'Max reduced chiSquared to enter ST.'),
# maximum allowed best reduced chiSquared out of ST before entering MCMC [double]
'cMaxMCMC':(10.0,'Max reduced chiSquared to enter MCMC.'),
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
'dmpInt'   : 1000000,
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
# force adding a value in degrees to argument of periapsis used in DI orbit fit [double]
'omegaPdi' : (0.0,"Custom fixed val added to DI omega in model"),
# Is the data in the DIdata.dat in PA,SA format? else, it is in E,N (x,y) format [bool]
'pasa'     : (False,"Is astrometry data in PA,SA format?"),
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
'ePrior'    :(True,'Use prior for eccentricity?',ePriorRatio),
'pPrior'    :(True,'Use prior for period?',pPriorRatio),
'incPrior'  :(True,'Use prior for inclination?',incPriorRatio),
'M1Prior':(False,'Use prior for M1?',mass1PriorRatio),
'M2Prior':(False,'Use prior for M2?',mass2PriorRatio),
'parPrior' :(True,'Use prior for parallax?',paraPriorRatio),
###################################################
# Ranges for acceptable random number inputs ######
###################################################
# Minimum/Maximum allowed value for the mass of the primary body [double][Msun]
'mass1MIN' : 0.5,
'mass1MAX' : 2.0,
# Minimum/Maximum allowed value for the mass of the secondary body [double][Msun]
'mass2MIN' : 0.1,
'mass2MAX' : 0.6,
# Minimum/Maximum allowed value for the Parallax [double][mas]
'paraMIN' : 30,
'paraMAX' : 50.0,
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
'PMIN' : 17.5,
'PMAX' : 22.5,
# Minimum/Maximum allowed value for the Inclination [double][deg]
'incMIN' : 145.0,
'incMAX' : 180.0,
# Minimum/Maximum allowed value for the Argument of Perigee [double][deg]
'omegaMIN' : 340.0,
'omegaMAX' : 360.0,
# Minimum/Maximum value for Semi-major amplitude of RV curve [m/s]
'KMIN' : 800,
'KMAX' : 900,
# Minimum/Maximum values of Radial Velocity Offsets.  
# Must be one per set of RV data in same order as data comes in RVdata.dat, or the a single value to be used by all [comma separated list of doubles]
'vMINs' :[6100,330,6250],
'vMAXs' :[6250,440,6400],
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
