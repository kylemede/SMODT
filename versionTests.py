import smodt1.code.pythonSide.toolboxes as tools1
import smodt2.tools as tools2
from smodt2.settings_and_inputData import constants
import os
import math
import numpy as np
import timeit

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
    e = 0.4
    Sys_Dist_PC = 5.0
    Mass1 = 1.0
    Mass2 = 0.2
    Omega = 60.0
    omega = 110.0
    T = 2457000.0
    T_center = 2457000.0
    P = 15.0
    inc =  30.0
    offset = 0.0
    
    #############################################
    if False:
        print "orig # epochs = "+str(realData2.shape[0])
        for i in range(10):
            realData2 = np.concatenate((realData2,realData2))
        print "new # epochs = "+str(realData2.shape[0])
        modelData2 = np.zeros((realData2.shape[0],3))
        #print 'new modelData.shape = '+repr(modelData2.shape)
    ############################################
    a_total = (((constants.Grav*constants.KGperMsun*(Mass1+Mass2)*((P*constants.secPerYear)**2.0))/(4.0*(constants.pi**2.0))))**(1.0/3.0)
    a_total = a_total/constants.MperAU
    a1 = a_total/(1.0+(Mass1/Mass2))
    #print 'a_total before = '+str(a_total)
    #print 'a1 before = '+str(a1)
    
    params = np.array([Mass1,Mass2,Sys_Dist_PC,Omega,e,T,T_center,P,inc,omega,a_total,0,0,offset])
    
    #for SMODT2.0 model
    tic=timeit.default_timer()
    Orbit = tools2.cppTools.Orbit()
    Orbit.loadRealData(realData2)
    Orbit.loadConstants(constants.Grav,constants.pi,constants.KGperMsun, constants.daysPerYear,constants.secPerYear,constants.MperAU)
    #print 'modelData2 BEFORE = \n'+repr(modelData2)
    #blah=blah
    Orbit.calculate(modelData2,params)
    toc=timeit.default_timer()
    swigTime=toc-tic
    print "\nfor swig it took: "+str(swigTime)
    #print 'atot model 2 = '+str(params[10])
    #print "model 2, K = "+str(params[12])
    #print 'modelData2 AFTER = \n'+repr(modelData2)
    
    #for SMODT1.0 model
    modelData1 = np.zeros(modelData2.shape)
    i=0
    tic=timeit.default_timer()
    for t in realData2[:,0]:
        (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SA, PA, x, y, a1, a2) = tools1.DItoolbox.orbitCalculatorTH_I(t, Sys_Dist_PC, inc, Omega, e, T, P, omega, a_total,Mass1, Mass2, verbose=False)
        #print 'a1 from DI = '+str(a1)
        #print 'TA_deg V1.0 = '+str(TA_deg)
        (v_r, K)=tools1.RVtoolbox.vrCalculatorSemiMajorType(t,e,T,P,omega,a1,T_center,inc, verbose=False)
        #print "model 1, K = "+str(K)
        modelData1[i,:] =[y,x,v_r] 
        i+=1
    toc=timeit.default_timer()
    pyTime=toc-tic
    print "for python it took: "+str(pyTime)
    print "Thus, SWIG is "+str(int(pyTime/swigTime))+" times faster than Python version."
    #print 'modelData1 = \n'+repr(modelData1)
    ##compare 
    diffAry = modelData2-modelData1
    #print 'modelData2-modelData1 = \n'+repr(diffAry)
    diffAry2 = np.where(diffAry<0,-1.0*diffAry,diffAry)#convert all negs to pos
    print '\ndiff maxs = '+str(np.max(diffAry2[:,0]))+", "+str(np.max(diffAry2[:,1]))+", "+str(np.max(diffAry2[:,2]))
    print 'diff mins = '+str(np.min(diffAry2[:,0]))+", "+str(np.min(diffAry2[:,1]))+", "+str(np.min(diffAry2[:,2]))
    
    
if __name__ == '__main__':
    calcTest()