import numpy as np
from PyAstronomy import pyasl

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
    massratio = 2
    M_primary = 1 #Solar masses
    distance = 5 #parsecs

    Omega = 70*np.pi/180 # Longitude of ascending node
    omega = 50*np.pi/180 # Argument of periastron
    i = 40*np.pi/180 # Inclination
    Npts = 10000

    G = 6.673e-8 #cgs
    mu = G*M_primary*1.989e33*(1 + 1./massratio) #gravitational parameter
    a = (mu*(period*86400*365.24)**2/4/np.pi**2)**(1./3) #in cm
    a_km = a/1e5 #to km
    a2 = a_km/(massratio + 1.)
    a1 = a_km - a2 # Semimajor axis of the low-mass component (in km)
    
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

    km_to_arcsec = 1/1.49598e8/distance # convert km to arcseconds

    data = np.zeros((pos_A.shape[0], 8))
    data[:, 0] = t*2*np.pi/period
    data[:, 1] = t
    data[:, 2] = pos_A[:, 0]*km_to_arcsec
    data[:, 3] = pos_A[:, 1]*km_to_arcsec
    data[:, 4] = vel_A[:, 2]
    data[:, 5] = pos_B[:, 0]*km_to_arcsec
    data[:, 6] = pos_B[:, 1]*km_to_arcsec
    data[:, 7] = vel_B[:, 2]

    np.savetxt('mockdata.dat', data, fmt="%.5g")


if __name__ == "__main__":
    calc_orbit()
