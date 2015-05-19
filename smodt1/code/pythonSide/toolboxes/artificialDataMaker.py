import numpy as np
import os
from PyAstronomy import pyasl
import DItoolbox as diTools

def calc_orbit():
    """
    This is a test function to produce a Keplerian orbit in RA,Dec,RV to verify SMODT.
    NOTE: the accuracy of the output values is accurate to 'one part in 10 to the 5', 
          due to the simple center differencing used to calculate the velocities from the positions and times.
    
    From his email:
    "
    The program calc_orbit.py gives you an eight column file:
    1. phase
    2. time (years)
    3. x position of secondary (arcsec)
    4. y position of secondary (arcsec)
    5. radial velocity of secondary (km/s)
    6. x position of primary (arcsec)
    7. y position of primary (arcsec)
    8. radial velocity of primary (km/s)
    
    I just computed a Kepler orbit's positions and velocities and rotated them.
    "
    """
    #Computer Directory
    baseSaveDir='/run/media/kmede/Data1/Todai_Work/Data/data_SMODT/'#$$$$$$$$$$$$$$$$$$$$ MAKE SURE THIS IS SET TO MACH YOUR COMPUTER!!! 
    #baseSaveDir = '/run/media/kmede/SharedData/Data/data_SMODT/'
    NumDataPointsOut = 10 #must be much less than 10000.  values between 10-500 are suitable.
    storePrimaryRVs = True
    percentError = 0.0001 #error is set to a percentage of the median
    realizeErrors = True
    overlapEnds = True # will ensure some points near end overlap the beginning of the orbit.

    #System settings
    massratio = 5.0
    M_primary = 1.0 #Solar masses
    distance = 5.0 #parsecs
    km_to_arcsec = 1/149597870.7/distance # convert km to arcsecond
    
    #constants
    G = 6.67384e-8 #cgs
    
    #Orbital Elements
    TimeLastPeri = 2457000.0 #JD
    e = 0.4
    period = 15. # years
    Omega = 170*np.pi/180 # Longitude of ascending node
    omega = 90*np.pi/180 # Argument of periastron
    i = 30*np.pi/180 # Inclination
 
    mu = G*M_primary*1.9884e33*(1 + 1./massratio) #gravitational parameter
    a = (mu*(period*86400*365.2422)**2/4/np.pi**2)**(1./3) #in cm
    a_km = a/1e5 #to km
    a_AU = a_km/149597870.7 #to AU
    a2 = a_km/(massratio + 1.)
    a1 = a_km - a2 # Semimajor axis of the low-mass component (in km)
    
    # print input orbital elements
    print "\n\nOrbital Elements Used:\ne = "+str(e)
    print "period = "+str(period)+" Years"
    print "LongAN = "+str(Omega*180.0/np.pi)+" deg"
    print "ArgPeri = "+str(omega*180.0/np.pi)+" deg"
    print "a_total = "+str(a_AU)+" AU"
    print "inclination = "+str(i*180.0/np.pi)+" deg"
    print "Time of Last Periapsis = "+str(TimeLastPeri)+" JD"
    print "Mass 1 = "+str(M_primary)+" Msun"
    print "Mass 2 = "+str(M_primary/massratio)+" Msun"
    print "System distance = "+str(distance)+" PC "
    if storePrimaryRVs:
        print "Saving RVs of primary star relative to Center of Mass\n"
    else:
        print "Saving RVs of companion relative to Center of Mass\n"
    #settings prints
    print 'Errors were calculated as '+str(percentError)+"% of the median value in the observables"
    if realizeErrors:
        print 'Data values were realized from the errors'
    else:
        print 'Data values are perfect with NO realization of the errors'
    print str(NumDataPointsOut)+" epochs will be calculated and stored\n"
        
        
    # Positions of both components in km relative to center of mass
    ke = pyasl.KeplerEllipse(a1, period, e=e, Omega=0.)
    NptsBIG = 10000
    t = (np.arange(NptsBIG) - 1)/(NptsBIG - 2.)*period
   
    ## Extend t to include a fifth of extra points at the end that overlap the beginning of the orbit
    if overlapEnds:
        NumOverlapPts = NptsBIG//10
        t2 = np.empty((t.shape[0]+NumOverlapPts))
        t2[0:t.size]=t
        t2[t.size:]=t[2]*np.arange(NumOverlapPts+1)[1:]+t[-1]
        t=t2
    pos_A = ke.xyzPos(t)
    pos_B = -pos_A/massratio
    
    # Velocities in km/s using centered differencing
    vel_A = (pos_A[2:] - pos_A[:-2])/(t[2] - t[0])/(86400*365.2422)
    if False:
        vel_A2 = pos_A.copy()
        vel_A2 = (pos_A[1:] - pos_A[:-1])/(t[1] - t[0])/(86400*365.2422)
        print "pos_A = "+repr(pos_A[0:10])+"\n"+repr(pos_A[-10:])+"\n"
        print "Vel_A before = "+repr(vel_A[0:10])+"\n"+repr(vel_A[-10:])+"\n"
        print "Vel_A after = "+repr(vel_A2[0:10])+"\n"+repr(vel_A2[-10:])+"\n"
        print "pos_A[2]-pos_A[0] = "+str(pos_A[2]-pos_A[0])+", t[2] - t[0] = "+str(t[2] - t[0])
        print "(pos_A[2]-pos_A[0])/(t[2] - t[0])/(86400*365.2422) = "+str((pos_A[2]-pos_A[0])/(t[2] - t[0])/(86400*365.24))
        print "\npos_A[1]-pos_A[0] = "+str(pos_A[1]-pos_A[0])+", t[1] - t[0] = "+str(t[1] - t[0])
        print "(pos_A[1]-pos_A[0])/(t[1] - t[0])/(86400*365.2422) = "+str((pos_A[1]-pos_A[0])/(t[1] - t[0])/(86400*365.24))
        print "(t[2] - t[0])*365.24 = "+str((t[2] - t[0])*365.2422)
        print "(t[1] - t[0])*365.24 = "+str((t[1] - t[0])*365.2422)
    pos_A = pos_A[1:-1]
    vel_B = (pos_B[2:] - pos_B[:-2])/(t[2] - t[0])/(86400*365.2422)
    pos_B = pos_B[1:-1]
    t = t[1:-1]

    # Construct rotation matrix (from wikipedia [http://en.wikipedia.org/wiki/Orbital_elements#Euler_angle_transformations])
    x1 = np.cos(Omega)*np.cos(omega) - np.sin(Omega)*np.cos(i)*np.sin(omega)
    x2 = np.sin(Omega)*np.cos(omega) + np.cos(Omega)*np.cos(i)*np.sin(omega)
    x3 = np.sin(i)*np.sin(omega)

    y1 = -np.cos(Omega)*np.sin(omega) - np.sin(Omega)*np.cos(i)*np.cos(omega)
    y2 = -np.sin(Omega)*np.sin(omega) + np.cos(Omega)*np.cos(i)*np.cos(omega)
    y3 = np.sin(i)*np.cos(omega)

    z1 = np.sin(i)*np.sin(Omega)
    z2 = -np.sin(i)*np.cos(Omega)
    z3 = np.cos(i)

    rotmat = np.asarray([[x1, x2, x3], [y1, y2, y3], [z1, z2, z3]])

    # Rotate positions, velocities
    pos_A = np.dot(pos_A, rotmat)
    vel_A = np.dot(vel_A, rotmat)
    pos_B = np.dot(pos_B, rotmat)
    vel_B = np.dot(vel_B, rotmat)
    
    ## re-sample position, vel and time arrays to requested number of samples
    pos_Anew = []
    pos_Bnew = []
    vel_Anew = []
    vel_Bnew = []
    tnew = []
    for i in range(0,len(t),int(len(t)/NumDataPointsOut)):
        pos_Anew.append(pos_A[i])
        pos_Bnew.append(pos_B[i])
        vel_Anew.append(vel_A[i])
        vel_Bnew.append(vel_B[i])
        tnew.append(t[i])
    
    pos_A = np.array(pos_Anew)
    pos_B = np.array(pos_Bnew)
    vel_A = np.array(vel_Anew)
    vel_B = np.array(vel_Bnew)
    t=np.array(tnew)

    data = np.zeros((pos_A.shape[0], 8))
    data[:, 0] = t*2*np.pi/period #1. phase
    data[:, 1] = t # 2. time (years)
    data[:, 2] = pos_A[:, 0]*km_to_arcsec #3. x position of secondary (arcsec)
    data[:, 3] = pos_A[:, 1]*km_to_arcsec #4. y position of secondary (arcsec)
    data[:, 4] = vel_A[:, 2] #5. radial velocity of secondary (km/s)
    data[:, 5] = pos_B[:, 0]*km_to_arcsec #6. x position of primary (arcsec)
    data[:, 6] = pos_B[:, 1]*km_to_arcsec #7. y position of primary (arcsec)
    data[:, 7] = vel_B[:, 2] #8. radial velocity of primary (km/s)
    
    data2 = np.zeros((pos_A.shape[0],5))
    data2[:,0] = data[:, 1]*365.2422+TimeLastPeri #JD 
    data2[:,1] = pos_A[:, 1]*km_to_arcsec - pos_B[:, 1]*km_to_arcsec #Ythi=Xplot=RA  separation between two bodies based on primary being at 0,0 ["]
    data2[:,2] = pos_A[:, 0]*km_to_arcsec - pos_B[:, 0]*km_to_arcsec #Xthi=Yplot=Dec  separation between two bodies based on primary being at 0,0 ["]
    data2[:,3] = vel_B[:, 2]*1000.0 # RV of primary compared to center of mass origin[ m/s]
    data2[:,4] = vel_A[:, 2]*1000.0 # RV of secondary compared to center of mass origin[ m/s]
    
    #calculate error and use it to realize the errors in the RV data if requested
    #calculate error and use it to realize the errors in the DI data if requested
    errorRA = np.median(np.abs(data2[:,1]))*(percentError/100.0)
    #print 'errorRA = '+str(errorRA)
    errorDec = np.median(np.abs(data2[:,2]))*(percentError/100.0)
    #print 'errorDec = '+str(errorDec)
    errorRA = np.max([errorRA,errorDec])
    errorDec = np.max([errorRA,errorDec])
    #print 'Using larger of two DI errors for both = '+str(errorDec)
    if realizeErrors:
        for i in range(pos_A.shape[0]):
            data2[i,1]+=np.random.normal(0,errorRA)
    if realizeErrors:
        for i in range(pos_A.shape[0]):
            data2[i,2] += np.random.normal(0,errorDec)
    errorRVprimary = np.median(np.abs(data2[:,3]))*(percentError/100.0)
    if realizeErrors:
        for i in range(pos_A.shape[0]):
            data2[i,3] += np.random.normal(0,errorRVprimary)
    errorRVsecondary = np.median(np.abs(data2[:,4]))*(percentError/100.0)
    if realizeErrors:
        for i in range(pos_A.shape[0]):
            data2[i,4] += np.random.normal(0,errorRVsecondary)
    
    #load up data for SMODT1.0 format
    #convert RA and Dec into the SA and PA used by SMODT
    data3 = np.empty((pos_A.shape[0],7))
    data3[:,0] = data2[:, 0]
    for i in range(0,pos_A.shape[0]):
        # convert x,y to SA and PA with fake errors
        #NOTE: there is an error in the resulting errors from ENtoPASA, but the PASA values are good.
        (data3[i,1],data3[i,2],data3[i,3],data3[i,4]) = diTools.ENtoPASA(data2[i,1], errorRA, data2[i,2], errorDec)
#         if False:
#             print repr((data2[i,1], null, data2[i,2], null))
#             (u1,u2,u3,u4)=diTools.PASAtoEN(data3[i,1],data3[i,2],data3[i,3],data3[i,4])
    #Averaging out the errors in PASA due to issue with ENtoPASA
    errorPAMean = np.mean(data3[:,2])
    #print 'errorPAMean = '+str(errorPAMean)+" = "+str((errorPAMean/np.median(np.abs(data3[:,1])))*100.0)+"% the median of the data"
    errorPA2 = np.median(np.abs(data3[:,1]))*(percentError/100.0)
    #print 'errorPA2 = '+str(errorPA2)
    #print 'errorPA2*1.7 = '+str(errorPA2*1.7)
    data3[:,2] = errorPA2*1.7
    errorSAMean = np.mean(data3[:,4])
    #print 'errorSAMean = '+str(errorSAMean)+" = "+str((errorSAMean/np.median(np.abs(data3[:,3])))*100.0)+"% the median of the data"
    errorSA2 = np.median(np.abs(data3[:,3]))*(percentError/100.0)
    #print "errorSA2 = "+str(errorSA2)
    data3[:,4] = errorSA2
    
#     #calculate error and use it to realize the errors in the DI data if requested
#     errorPA = np.median(np.abs(data3[:,1]))*(percentError/100.0)
#     data3[:,2] = errorPA*1.7
#     if realizeErrors:
#         for i in range(pos_A.shape[0]):
#             data3[i,1] += np.random.normal(0,errorPA)
#     errorSA = np.median(np.abs(data3[:,3]))*(percentError/100.0)
#     data3[:,4] = errorSA
#     if realizeErrors:
#         for i in range(pos_A.shape[0]):
#             data3[i,3]+=np.random.normal(0,errorSA)
        
    #output data3 has format ->SMODT format with DI and RV combined into one file
    #1. JD
    #2. Position Angle [deg]
    #3. Position Angle ERROR [deg]
    #4. Separation Angle ["]
    #5. Separation Angle ERROR ["]
    #6. RV of primary (or secondary) rel to CofM [m/s]
    #7. RV ERROR [m/s]
    
    
    #load up data for SMODT2.0 format
    data4 = np.empty((pos_A.shape[0],7))
    data4[:,0] = data2[:, 0]#1. JD
    data4[:,1] = data2[:,1]#2. RA (x) ["]
    data4[:,2] = errorRA #3. RA ERROR ["]
    data4[:,3] = data2[:,2]#4. Dec (y) ["]
    data4[:,4] = errorDec#5. Dec ERROR ["]
    #output data4 has format ->SMODT format with DI and RV combined into one file
    #1. JD
    #2. RA (x) ["]
    #3. RA ERROR ["]
    #4. Dec (y) ["]
    #5. Dec ERROR ["]
    #6. RV of primary (or secondary) rel to CofM [m/s]
    #7. RV ERROR [m/s]
    if storePrimaryRVs:
        data3[:,5] = data2[:,3]
        data3[:,6] = errorRVprimary
        data4[:,5] = data2[:,3]
        data4[:,6] = errorRVprimary
    else:
        data3[:,5] = data2[:,4]
        data3[:,6] = errorRVsecondary
        data4[:,5] = data2[:,4]
        data4[:,6] = errorRVsecondary
    print "resulting data files have "+str(len(data3[:,3]))+" epochs"
    ## split these up into two separate files
    dataDI1 = np.empty((pos_A.shape[0],5))
    dataDI1[:,0] = data3[:,0]
    dataDI1[:,1] = data3[:,1]
    dataDI1[:,2] = data3[:,2]
    dataDI1[:,3] = data3[:,3]
    dataDI1[:,4] = data3[:,4]
    dataDI2 = np.empty((pos_A.shape[0],5))
    dataDI2[:,0] = data4[:,0]
    dataDI2[:,1] = data4[:,1]
    dataDI2[:,2] = data4[:,2]
    dataDI2[:,3] = data4[:,3]
    dataDI2[:,4] = data4[:,4]
    dataRV = np.empty((pos_A.shape[0],3))
    dataRV[:,0] = data3[:,0]
    dataRV[:,1] = data3[:,5]
    dataRV[:,2] = data3[:,6]
    
    ##write files to disk
    if False:
        np.savetxt(os.path.join(baseSaveDir,'mockdata.dat'), data, fmt="%.10g")
        #np.savetxt(os.path.join(baseSaveDir,'mockdata-SMODTformat.dat'), data3, fmt="%.10g")
    if True:
        np.savetxt(os.path.join(baseSaveDir,'mockdata-SMODT1format-DI.dat'), dataDI1, fmt="%.10g")
        np.savetxt(os.path.join(baseSaveDir,'mockdata-SMODTformat-RV.dat'), dataRV, fmt="%.10g")
        np.savetxt(os.path.join(baseSaveDir,'mockdata-SMODT2format-DI.dat'), dataDI2, fmt="%.10g")
    print '\nOutput data files written to:\n'+baseSaveDir

if __name__ == "__main__":
    calc_orbit()
