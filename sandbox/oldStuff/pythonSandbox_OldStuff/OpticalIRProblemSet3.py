import math as m

H = 0.2#m
T = 273.0#K
p = 18.5e-6#m
d = 37.888e-3#m
F = 12.0#unitless
f = 12.0 #m
h = 6.626e-34 #Js
c = 299792458 #m/s

print "PROBLEM SET 3 ANSWERS\n"

## PROBLEM 3-1
print "\nAnswers to problem 3-1:\n"
# First for Ks-band
Lambda = 2.15 ##microns
deltaLambda = 0.3 ##microns

Na = (d*m.pow(p,2.0)*deltaLambda*(1.191e8))/(8.0*H*h*c*m.pow(Lambda,4.0)*(m.exp(14387.7/(Lambda*T))-1.0))
Nb = ((d/H)+(m.sqrt(2.0)/F))*m.pow(10.0,-6.0)
N = Na*Nb
print "Number of photons/second in Ks-band= ",N

# Second for K-band
Lambda = 2.2 ##microns
deltaLambda = 0.35 ##microns

Na = (d*m.pow(p,2.0)*deltaLambda*(1.191e8))/(8.0*H*h*c*m.pow(Lambda,4.0)*(m.exp(14387.7/(Lambda*T))-1.0))
Nb = ((d/H)+(m.sqrt(2.0)/F))*m.pow(10.0,-6.0)
N = Na*Nb
print "Number of photons/second in K-band= ",N

#---------------------------------------------------------------------------------

## PROBLEM 3-2
print "\nAnswers to problem 3-2:\n"
# Dividing the Ks-band up into 3 sections [1.85-2.05], [2.05-2.25] and [2.25-2.45]
# These corresponds to deltaLambda = 0.1 and Lambda = 1.95, 2.25 and 2.35.
 
# Calculating N for each band:

Lambda = 1.95 ##microns
deltaLambda = 0.1 ##microns
Na = (d*m.pow(p,2.0)*deltaLambda*(1.191e8))/(8.0*H*h*c*m.pow(Lambda,4.0)*(m.exp(14387.7/(Lambda*T))-1.0))
Nb = ((d/H)+(m.sqrt(2.0)/F))*m.pow(10.0,-6.0)
N = Na*Nb
print "Number of photons/second for [1.85-2.05] of Ks-band = ",N

Lambda = 2.25 ##microns
deltaLambda = 0.1 ##microns
Na = (d*m.pow(p,2.0)*deltaLambda*(1.191e8))/(8.0*H*h*c*m.pow(Lambda,4.0)*(m.exp(14387.7/(Lambda*T))-1.0))
Nb = ((d/H)+(m.sqrt(2.0)/F))*m.pow(10.0,-6.0)
N = Na*Nb
print "Number of photons/second for [2.05-2.25] of Ks-band = ",N

Lambda = 2.35 ##microns
deltaLambda = 0.1 ##microns
Na = (d*m.pow(p,2.0)*deltaLambda*(1.191e8))/(8.0*H*h*c*m.pow(Lambda,4.0)*(m.exp(14387.7/(Lambda*T))-1.0))
Nb = ((d/H)+(m.sqrt(2.0)/F))*m.pow(10.0,-6.0)
N = Na*Nb
print "Number of photons/second for [2.25-2.45] of Ks-band = ",N

#---------------------------------------------------------------------------------

# Dividing up the K-band into 4 sections [1.85-2.025], [2.025-2.2], [2.20-2.375] and [2.375-2.55] 
# These correspond to deltaLambda = 0.0875 and Lambda = 1.9375, 2.1125, 2.2875 and 2.4625

# Calculating N for each band:

Lambda = 1.9375 ##microns
deltaLambda = 0.0875 ##microns
Na = (d*m.pow(p,2.0)*deltaLambda*(1.191e8))/(8.0*H*h*c*m.pow(Lambda,4.0)*(m.exp(14387.7/(Lambda*T))-1.0))
Nb = ((d/H)+(m.sqrt(2.0)/F))*m.pow(10.0,-6.0)
N = Na*Nb
print "\nNumber of photons/second for [1.85-2.025] of K-band = ",N

Lambda = 2.1125 ##microns
deltaLambda = 0.0875 ##microns
Na = (d*m.pow(p,2.0)*deltaLambda*(1.191e8))/(8.0*H*h*c*m.pow(Lambda,4.0)*(m.exp(14387.7/(Lambda*T))-1.0))
Nb = ((d/H)+(m.sqrt(2.0)/F))*m.pow(10.0,-6.0)
N = Na*Nb
print "Number of photons/second for [2.025-2.2] of K-band = ",N

Lambda = 2.2875 ##microns
deltaLambda = 0.0875 ##microns
Na = (d*m.pow(p,2.0)*deltaLambda*(1.191e8))/(8.0*H*h*c*m.pow(Lambda,4.0)*(m.exp(14387.7/(Lambda*T))-1.0))
Nb = ((d/H)+(m.sqrt(2.0)/F))*m.pow(10.0,-6.0)
N = Na*Nb
print "Number of photons/second for [2.20-2.375] of K-band = ",N

Lambda = 2.4625 ##microns
deltaLambda = 0.0875 ##microns
Na = (d*m.pow(p,2.0)*deltaLambda*(1.191e8))/(8.0*H*h*c*m.pow(Lambda,4.0)*(m.exp(14387.7/(Lambda*T))-1.0))
Nb = ((d/H)+(m.sqrt(2.0)/F))*m.pow(10.0,-6.0) 
N = Na*Nb
print "Number of photons/second for [2.375-2.55] of K-band = ",N

#**********************************************************************************

## PROBLEM 2
# Using the same band sections of problem 3-2
print "\nAnswers to problem 3-3:\n"
N = 10.0#K

# First for the Ks-band:

Lambda = 1.95 ##microns
deltaLambda = 0.1 ##microns
insideA = (d*m.pow(p,2.0)*deltaLambda*(1.191e8))/(8.0*N*H*h*c*m.pow(Lambda,4.0))
insideB = ((d/H)+(m.sqrt(2.0)/F))*m.pow(10.0,-6.0)
#NOTE: m.log1p calculates the natural logarithm of x+1, thus the +1 of the derived equation is dropped
Tout = 14387.7/(Lambda*m.log1p(insideA*insideB))
print "Temperature where only 10 photons/second are received for [1.85-2.05] of Ks-band = ",Tout

Lambda = 2.25 ##microns
deltaLambda = 0.1 ##microns
insideA = (d*m.pow(p,2.0)*deltaLambda*(1.191e8))/(8.0*N*H*h*c*m.pow(Lambda,4.0))
insideB = ((d/H)+(m.sqrt(2.0)/F))*m.pow(10.0,-6.0)
Tout = 14387.7/(Lambda*m.log1p(insideA*insideB))
print "Temperature where only 10 photons/second are received for [2.05-2.25] of Ks-band = ",Tout

Lambda = 2.35 ##microns
deltaLambda = 0.1 ##microns
insideA = (d*m.pow(p,2.0)*deltaLambda*(1.191e8))/(8.0*N*H*h*c*m.pow(Lambda,4.0))
insideB = ((d/H)+(m.sqrt(2.0)/F))*m.pow(10.0,-6.0)
Tout = 14387.7/(Lambda*m.log1p(insideA*insideB))
print "Temperature where only 10 photons/second are received for [2.25-2.45] of Ks-band = ",Tout

#---------------------------------------------------------------------------------
# Second for the K-band:

Lambda = 1.9375 ##microns
deltaLambda = 0.0875 ##microns
insideA = (d*m.pow(p,2.0)*deltaLambda*(1.191e8))/(8.0*N*H*h*c*m.pow(Lambda,4.0))
insideB = ((d/H)+(m.sqrt(2.0)/F))*m.pow(10.0,-6.0)
Tout = 14387.7/(Lambda*m.log1p(insideA*insideB))
print "\nTemperature where only 10 photons/second are received for [1.85-2.025] of K-band = ",Tout

Lambda = 2.1125 ##microns
deltaLambda = 0.0875 ##microns
insideA = (d*m.pow(p,2.0)*deltaLambda*(1.191e8))/(8.0*N*H*h*c*m.pow(Lambda,4.0))
insideB = ((d/H)+(m.sqrt(2.0)/F))*m.pow(10.0,-6.0)
Tout = 14387.7/(Lambda*m.log1p(insideA*insideB))
print "Temperature where only 10 photons/second are received for [2.025-2.2] of K-band = ",Tout

Lambda = 2.2875 ##microns
deltaLambda = 0.0875 ##microns
insideA = (d*m.pow(p,2.0)*deltaLambda*(1.191e8))/(8.0*N*H*h*c*m.pow(Lambda,4.0))
insideB = ((d/H)+(m.sqrt(2.0)/F))*m.pow(10.0,-6.0)
Tout = 14387.7/(Lambda*m.log1p(insideA*insideB))
print "Temperature where only 10 photons/second are received for [2.20-2.375] of K-band = ",Tout

Lambda = 2.4625 ##microns
deltaLambda = 0.0875 ##microns
insideA = (d*m.pow(p,2.0)*deltaLambda*(1.191e8))/(8.0*N*H*h*c*m.pow(Lambda,4.0))
insideB = ((d/H)+(m.sqrt(2.0)/F))*m.pow(10.0,-6.0)
Tout = 14387.7/(Lambda*m.log1p(insideA*insideB))
print "Temperature where only 10 photons/second are received for [2.375-2.55] of K-band = ",Tout






