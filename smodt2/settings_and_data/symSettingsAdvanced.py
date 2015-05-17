import numpy as np
import symSettings
import constants

exampleSetting1 = 10.0

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
