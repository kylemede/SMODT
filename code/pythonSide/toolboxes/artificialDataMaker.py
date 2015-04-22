import numpy as np
from PyAstronomy import pyasl
import DItoolbox as diTools

def calc_orbit():
    """
    This is a test function to produce a Keplerian orbit in RA,Dec,RV to verify SMODT.
    I believe it was written by Tim Brandt.
    
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
    #model settings
    Npts = 50
    storePrimaryRVs = True

    #System settings
    massratio = 2.0
    M_primary = 1.0 #Solar masses
    distance = 5.0 #parsecs
    km_to_arcsec = 1/1.49598e8/distance # convert km to arcsecond
    
    #constants
    G = 6.673e-8 #cgs
    
    #Orbital Elements
    TimeLastPeri = 2457000.0 #JD
    e = 0.5
    period = 5. # years
    Omega = 70*np.pi/180 # Longitude of ascending node
    omega = 50*np.pi/180 # Argument of periastron
    i = 40*np.pi/180 # Inclination
 
    mu = G*M_primary*1.989e33*(1 + 1./massratio) #gravitational parameter
    a = (mu*(period*86400*365.24)**2/4/np.pi**2)**(1./3) #in cm
    a_km = a/1e5 #to km
    a_AU = a_km/149597871. #to AU
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
        print "saving RVs of primary star relative to Center of Mass\n"
    else:
        print "saving RVs of companion relative to Center of Mass\n"
    # Positions of both components in km relative to center of mass

    ke = pyasl.KeplerEllipse(a1, period, e=e, Omega=0.)

    t = (np.arange(Npts) - 1)/(Npts - 2.)*period
    pos_A = ke.xyzPos(t)
    pos_B = -pos_A/massratio
    
    # Velocities in km/s using centered differencing

    vel_A = pos_A.copy()
    vel_A = (pos_A[2:] - pos_A[:-2])/(t[2] - t[0])/(86400*365.24)
    if False:
        vel_A2 = pos_A.copy()
        vel_A2 = (pos_A[1:] - pos_A[:-1])/(t[1] - t[0])/(86400*365.24)
        print "pos_A = "+repr(pos_A[0:10])+"\n"+repr(pos_A[-10:])+"\n"
        print "Vel_A before = "+repr(vel_A[0:10])+"\n"+repr(vel_A[-10:])+"\n"
        print "Vel_A after = "+repr(vel_A2[0:10])+"\n"+repr(vel_A2[-10:])+"\n"
        print "pos_A[2]-pos_A[0] = "+str(pos_A[2]-pos_A[0])+", t[2] - t[0] = "+str(t[2] - t[0])
        print "(pos_A[2]-pos_A[0])/(t[2] - t[0])/(86400*365.24) = "+str((pos_A[2]-pos_A[0])/(t[2] - t[0])/(86400*365.24))
        print "\npos_A[1]-pos_A[0] = "+str(pos_A[1]-pos_A[0])+", t[1] - t[0] = "+str(t[1] - t[0])
        print "(pos_A[1]-pos_A[0])/(t[1] - t[0])/(86400*365.24) = "+str((pos_A[1]-pos_A[0])/(t[1] - t[0])/(86400*365.24))
        print "(t[2] - t[0])*365.24 = "+str((t[2] - t[0])*365.24)
        print "(t[1] - t[0])*365.24 = "+str((t[1] - t[0])*365.24)
    pos_A = pos_A[1:-1]

    vel_B = pos_B.copy()
    vel_B = (pos_B[2:] - pos_B[:-2])/(t[2] - t[0])/(86400*365.24)
    pos_B = pos_B[1:-1]
    
    if False:
        print repr(t[0:10])+"\n"+repr(t[-10:])
    t = t[1:-1]
    if False:
        print "\n"+repr(t[0:10])+"\n"+repr(t[-10:])
    
    if False:
        print "\n"+"*"*75+"\nbefore rotation:"
        print "t[o] = "+str(t[0])
        if True:
            print "pos_A[0:10] = "+repr(pos_A[0:10])
            print "vel_A[0:10] = "+repr(vel_A[0:10])
        else:
            print "pos_A[0] = "+repr(pos_A[0])
            print "pos_B[0] = "+repr(pos_B[0])
            print "vel_A[0] = "+repr(vel_A[0])
            print "vel_B[0] = "+repr(vel_B[0])
            print "\npos_A[1] = "+repr(pos_A[1])
            print "pos_B[1] = "+repr(pos_B[1])
            print "vel_A[1] = "+repr(vel_A[1])
            print "vel_B[1] = "+repr(vel_B[1])
            print "\npos_A[50] = "+repr(pos_A[50])
            print "pos_B[50] = "+repr(pos_B[50])
            print "vel_A[50] = "+repr(vel_A[50])
            print "vel_B[50] = "+repr(vel_B[50])
            print "\npos_A[-1] = "+repr(pos_A[-1])
            print "pos_B[-1] = "+repr(pos_B[-1])
            print "vel_A[-1] = "+repr(vel_A[-1])
            print "vel_B[-1] = "+repr(vel_B[-1])+"\n"+"*"*75

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
    
    if False:
        print "\n"+"*"*75+"\nAfter rotation:"
        print "t[o] = "+str(t[0])
        if True:
            print "pos_A[0:10] = "+repr(pos_A[0:10])
            print "vel_A[0:10] = "+repr(vel_A[0:10])
        else:
            print "pos_A[0] = "+repr(pos_A[0])
            print "pos_B[0] = "+repr(pos_B[0])
            print "vel_A[0] = "+repr(vel_A[0])
            print "vel_B[0] = "+repr(vel_B[0])
            print "\npos_A[1] = "+repr(pos_A[1])
            print "pos_B[1] = "+repr(pos_B[1])
            print "vel_A[1] = "+repr(vel_A[1])
            print "vel_B[1] = "+repr(vel_B[1])
            print "\npos_A[50] = "+repr(pos_A[50])
            print "pos_B[50] = "+repr(pos_B[50])
            print "vel_A[50] = "+repr(vel_A[50])
            print "vel_B[50] = "+repr(vel_B[50])
            print "\npos_A[-1] = "+repr(pos_A[-1])
            print "pos_B[-1] = "+repr(pos_B[-1])
            print "vel_A[-1] = "+repr(vel_A[-1])
            print "vel_B[-1] = "+repr(vel_B[-1])+"\n"+"*"*75

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
    data2[:,0] = data[:, 1]*365.24+TimeLastPeri #JD 
    data2[:,1] = pos_A[:, 1]*km_to_arcsec - pos_B[:, 1]*km_to_arcsec #Ythi=Xplot=RA  separation between two bodies based on primary being at 0,0 ["]
    data2[:,2] = pos_A[:, 0]*km_to_arcsec - pos_B[:, 0]*km_to_arcsec #Xthi=Yplot=Dec  separation between two bodies based on primary being at 0,0 ["]
    data2[:,3] = vel_B[:, 2]*1000.0 # RV of primary compared to center of mass origin[ m/s]
    data2[:,4] = vel_A[:, 2]*1000.0 # RV of secondary compared to center of mass origin[ m/s]
    
    #calculate error
    errorRA = np.median(np.abs(data2[:,1]))*0.05
    errorDec = np.median(np.abs(data2[:,2]))*0.05
    errorRVprimary = np.median(np.abs(data2[:,3]))*0.05
    errorRVsecondary = np.median(np.abs(data2[:,4]))*0.05
    
    data3 = np.empty((pos_A.shape[0],7))
    data3[:,0] = data2[:, 0]
    for i in range(0,pos_A.shape[0]):
        # convert x,y to SA and PA with fake errors
        (data3[i,1],data3[i,2],data3[i,3],data3[i,4]) = diTools.ENtoPASA(data2[i,1], errorRA, data2[i,2], errorDec)
        if False:
            print repr((data2[i,1], errorRA, data2[i,2], errorDec))
            (u1,u2,u3,u4)=diTools.PASAtoEN(data3[i,1],data3[i,2],data3[i,3],data3[i,4])
       
     
    if storePrimaryRVs:
        data3[:,5] = data2[:,3]
        data3[:,6] = errorRVprimary
    else:
        data3[:,5] = data2[:,4]
        data3[:,6] = errorRVsecondary
        
    #output data3 has format
    #1. JD
    #2. Position Angle [deg]
    #3. Position Angle ERROR [deg]
    #4. Separation Angle ["]
    #5. Separation Angle ERROR ["]
    #6. RV of primary (or secondary) rel to CofM [m/s]
    #7. RV ERROR [m/s]

    np.savetxt('/mnt/Data1/Todai_Work/Data/data_SMODT/mockdata.dat', data, fmt="%.10g")
    np.savetxt('/mnt/Data1/Todai_Work/Data/data_SMODT/mockdata-SMODTformat.dat', data3, fmt="%.10g")


if __name__ == "__main__":
    calc_orbit()
