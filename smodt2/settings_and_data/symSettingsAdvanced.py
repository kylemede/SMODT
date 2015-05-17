import numpy as np
import symSettings
import constants

exampleSetting1 = 10.0

def e_prior(e):
    return 2.0*e
def P_prior(P):
    if (symSettings.PMAX!=0)and(symSettings.PMIN!=0):
        return P
    else:
        return 1.0
    
def inc_prior(inc):
    return np.sin(inc*(constants.pi/180.0))
