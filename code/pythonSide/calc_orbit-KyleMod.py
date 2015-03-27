import numpy as np
from PyAstronomy import pyasl
from toolboxes import DItoolbox as diTools

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

    e = 0.5
    period = 5. # years
    massratio = 2.0
    M_primary = 1.0 #Solar masses
    distance = 5.0 #parsecs

    Omega = 70*np.pi/180 # Longitude of ascending node
    omega = 50*np.pi/180 # Argument of periastron
    i = 40*np.pi/180 # Inclination
    Npts = 1000

    G = 6.673e-8 #cgs
    mu = G*M_primary*1.989e33*(1 + 1./massratio) #gravitational parameter
    a = (mu*(period*86400*365.24)**2/4/np.pi**2)**(1./3) #in cm
    a_km = a/1e5 #to km
    a_AU = a_km/149597871. #to AU
    a2 = a_km/(massratio + 1.)
    a1 = a_km - a2 # Semimajor axis of the low-mass component (in km)
    
    # print input orbital elements
    print "e = "+str(e)
    print "period = "+str(period)
    print "LongAN = "+str(Omega*180.0/np.pi)
    print "ArgPeri = "+str(omega*180.0/np.pi)
    print "a_total = "+str(a_AU)
    print "inclination = "+str(i*180.0/np.pi)
    print "Mass 1 = "+str(M_primary)
    print "Mass 2 = "+str(M_primary/massratio)
    print "distance in PC = "+str(distance)
    
    # Positions of both components in km relative to center of mass

    ke = pyasl.KeplerEllipse(a1, period, e=e, Omega=0.)

    t = (np.arange(Npts) - 1)/(Npts - 2.)*period
    pos_A = ke.xyzPos(t)
    pos_B = -pos_A/massratio
    
    # Velocities in km/s using centered differencing

    vel_A = pos_A.copy()
    vel_A = (pos_A[2:] - pos_A[:-2])/(t[2] - t[0])/(86400*365.24)
    pos_A = pos_A[1:-1]

    vel_B = pos_B.copy()
    vel_B = (pos_B[2:] - pos_B[:-2])/(t[2] - t[0])/(86400*365.24)
    pos_B = pos_B[1:-1]

    t = t[1:-1]

    # Construct rotation matrix (from wikipedia)

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

    km_to_arcsec = 1/1.49598e8/distance # convert km to arcsecond

    data = np.zeros((pos_A.shape[0], 8))
    data[:, 0] = t*2*np.pi/period #1. phase
    data[:, 1] = t # 2. time (years)
    data[:, 2] = pos_A[:, 0]*km_to_arcsec #3. x position of secondary (arcsec)
    data[:, 3] = pos_A[:, 1]*km_to_arcsec #4. y position of secondary (arcsec)
    data[:, 4] = vel_A[:, 2] #5. radial velocity of secondary (km/s)
    data[:, 5] = pos_B[:, 0]*km_to_arcsec #6. x position of primary (arcsec)
    data[:, 6] = pos_B[:, 1]*km_to_arcsec #7. y position of primary (arcsec)
    data[:, 7] = vel_B[:, 2] #8. radial velocity of primary (km/s)
    
    data2 = np.zeros((pos_A.shape[0],4))
    data2[:,0] = data[:, 1]*365.25+2457000.0 #JD 
    data2[:,1] = pos_A[:, 1]*km_to_arcsec - pos_B[:, 1]*km_to_arcsec #Ythi=Xplot=RA  separation between two bodies based on primary being at 0,0 ["]
    data2[:,2] = pos_A[:, 0]*km_to_arcsec - pos_B[:, 0]*km_to_arcsec #Xthi=Yplot=Dec  separation between two bodies based on primary being at 0,0 ["]
    data2[:,3] = -vel_B[:, 2]*1000.0 # RV of primary compared to center of mass origin[ m/s]
    
    data3 = np.zeros((pos_A.shape[0],7))
    data3[:,0] = data2[:, 0]
    for i in range(0,pos_A.shape[0]):
        (data3[i,1],data3[i,2],data3[i,3],data3[i,4]) = diTools.ENtoPASA(data2[i,1], data2[i,1]*0.05, data2[i,2], data2[i,2]*0.05)
        
    data3[:,5] = data2[:,3]
    data3[:,6] = data2[:,3]*0.05
    
    #output data3 has format
    #1. JD
    #2. Position Angle [deg]
    #3. Position Angle ERROR [deg]
    #4. Separation Angle ["]
    #5. Separation Angle ERROR ["]
    #6. RV of primary rel to CofM [m/s]
    #7. RV ERROR [m/s]

    np.savetxt('mockdata.dat', data, fmt="%.10g")
    np.savetxt('mockdata-SMODTformat.dat', data3, fmt="%.10g")


if __name__ == "__main__":
    calc_orbit()
