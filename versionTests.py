import smodt1.code.pythonSide.toolboxes as tools1
import smodt2.tools as tools2
from smodt2.settings_and_inputData import constants
import os
import math
import numpy as np

def calcTest():
    #np.set_printoptions(precision=15)
    dataFilenameRoot = "/run/media/kmede/Data1/Todai_Work/Dropbox/EclipseWorkspaceDB/SMODT/FakeData_"
    diFilename1 = dataFilenameRoot+'DIdata.dat'
    rvFilename1 = dataFilenameRoot+'RVdata.dat'
    log = tools2.getLogger('main',lvl=100)
    realData2 = tools2.loadRealData(dataFilenameRoot+'2')
    modelData2 = np.zeros((realData2.shape[0],3))
    RVdataDict = tools1.RVtoolbox.RVdataToDict(rvFilename1)
    DIdataDict = tools1.DItoolbox.DIdataToDict(diFilename1)
    e = 0.5
    Sys_Dist_PC = 5.0
    Mass1 = 1.0
    Mass2 = 0.5
    Omega = 70.0
    omega = 50.0
    T = 2457000.0
    T_center = 2457000.0
    P = 5.0
    inc =  40.0
    offset = 0.0
    
    a_total = (((constants.Grav*constants.KGperMsun*(Mass1+Mass2)*((P*constants.secPerYear)**2.0))/(4.0*(constants.pi**2.0))))**(1.0/3.0)
    a_total = a_total/constants.MperAU
    a1 = a_total/(1.0+(Mass2/Mass1))
    
    params = np.array([Mass1,Mass2,Sys_Dist_PC,Omega,e,T,T_center,P,inc,omega,a_total,0,0,offset])
    
    Orbit = tools2.cppTools.Orbit()
    Orbit.loadRealData(realData2)
    Orbit.loadConstants(constants.Grav,constants.pi,constants.KGperMsun, constants.daysPerYear,constants.secPerYear,constants.MperAU)
    print 'modelData2 BEFORE = '+repr(modelData2)
    Orbit.calculate(modelData2,params)
    print 'modelData2 AFTER = '+repr(modelData2)
    
    epochs = realData2[:,0]
    RVs1 = []
    xs1 =[]
    ys1 = []
    tits=oK
    for t in epochs:
        (v_r, K)=tools1.RVtoolbox.vrCalculatorSemiMajorType(t,e,T,P,omega,a1,T_center,inc, K=False, verbose=False)
        (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SA, PA, x, y, a1, a2) = tools1.DItoolbox.orbitCalculatorTH_I(t, Sys_Dist_PC, inc, Omega, e, T, P, omega, a_total,Mass1, Mass2, verbose=False)
        RVs1.append(v_r)
        xs1.append(x)
        ys1.append(y)
    
if __name__ == '__main__':
    calcTest()