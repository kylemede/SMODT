import numpy as np
import symSettings
import constants

########################
### General Settings ###
########################
# This will set the maximum reduced ChiSquared value to accept and write to the output file.  ONLY for Monte Carlo and Simulated Annealing, not MCMC!! [double]
chiMAX = 3000.0
# set to 'true' to have NOTHING print to screen while running [bool]
SILENT = False
# set to 'false' to receive extra prints from the main simulations progress for testing [bool]
quiet = True
# set to 'true' to receive prints from the functions called by main for testing [bool]
verbose = False
# Simulated Annealing starting temperature [double]
strtTemp = 800.0
# make plot of posterior distributions? [bool]
pltDists =True
# make plots of RV and DI/AM orbit fits [bool]
pltOrbit =True
# Delete chain files after simulation is complete? [bool]
delChains =True
# Delete combined data files after simulation is complete? [bool]
delCombined =False
###$$$$$$$$$$$$$$$$$$$$$$ Keep in final version $$$$$$$$$$$$$$$$$$$$$$$$$$
# Copy output non-data files to a Dropbox folder? [bool]  
CopyToDB =False
dbFolder = '/run/media/kmede/Data1/Todai_Work/Dropbox/SMODT-outputCopies/'
##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############################
# Settings for MCMC mode ###
############################
# Calculate the length of the burn in for each chain (must be more than 1 chain)? [bool] 
CalcBurn =True
# remove burn-in of output MCMC chains before combining (must be more than 1 chain) (should already be handled by SimAnneal stage2 though...)?
delBurn = True
# Calculate the Correlation lengths and number of effective points of each chain (must be more than 1 chain)? [bool]
calcCL =True
# number of samples to draw for simulated annealing stage [int] 
nSAsamp =1000000
# number of samples to draw for sigma tuning stage [int] 
nSTsamp =1000000
# Make plots of MCMC progress plots? [bool]
pltMCMCprog =False
# Make plots of Simulated Annealing progress plots? [bool]
pltSAprog =False
# Calculate the Gelman-Rubin statistic? [bool]
CalcGR =True
# How many times do you want the Gelman-Rubin statistic calculated [int]
nGRcalc =10
#####################################
# Special Settings for the models ###
#####################################
# fit to the primary's RV orbit [bool]
fitPrime = False
# Are the RVs in the RVdata.dat for the Primary star? [bool]
primeRVs = True
# Step through parameter space in Time of Center Transit (Inferior Conjunction)?  [bool]
TcStep = False
# take the time of center transit (inferior conjunction) into account? [bool]
TcEqualT = True
# force adding a value in degrees to argument of periapsis used in RV orbit fit [double]
omegaPrv = 0.0
##################################################
## Special settings DI model:
# force adding a value in degrees to argument of periapsis used in RV orbit fit [double]
omegaPdi = 0.0



#Define the priors as python functions.
#NOTE: only change the code and not the name of the functions or their inputs.
def e_prior(e,P):
    if (P*constants.daysPerYear<1000.0)and(symSettings.eMAX!=0):
        return 2.0*e
    else:
        return 1.0
def P_prior(P):
    if (symSettings.PMAX!=0)and(symSettings.PMIN!=0):
        return P
    else:
        return 1.0
    
def inc_prior(inc):
    if (symSettings.incMAX!=0)and(symSettings.incMIN!=0):
        return np.sin(inc*(constants.pi/180.0))
    else:
        return 1.0
