import numpy as np

Grav = 6.67384e-11 #from physics.nist.gov
pi = np.pi
KGperMsun = 1.9884e30 #from asa.usno.navy.mil/static/files/2014/Astronomical_Constants_2014.pdf
daysPerYear = 365.25 #for Julian Calendar.  It would be 365.2422 for Gregorian, but in Astronomy we always use JD values for epochs, T and Tc. http://pumas.jpl.nasa.gov/files/04_21_97_1.pdf
secPerYear = 60*60*24*daysPerYear
MperAU = 149597870700.0 #from asa.usno.navy.mil/static/files/2014/Astronomical_Constants_2014.pdf
KGperMjupiter = 1.8983e27 #from http://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html