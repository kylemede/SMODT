import math
import numpy as np
import pylab
from math import pi
plt = pylab.matplotlib.pyplot

fig = plt.figure(1, figsize=(35,22) ,dpi=250)
ax = fig.add_subplot(111)

NUs = range(1,100000)
nu = 100000
while nu<1e8:
    nu = nu+1000
    NUs.append(nu)
print 'nu = ',nu
while nu<1e12:
    nu = nu+1000000
    NUs.append(nu)
print 'nu = ',nu
while nu<1e14:
    nu = nu+100000000
    NUs.append(nu)
print 'nu = ',nu
print 'len(NUs) = ',str(len(NUs))
y_vals = []

A = (3.0**(1.0/2.0))/pi
C = 17.7
max = 47.4042746513
tenPercentRange =  []
onePercentRange = []
print 'starting y_val calculation loop'
for nu in NUs:
    B = ((1.0e4)**(3.0/2.0))/nu
    y_val = A*(math.log(B)+C)
    y_val = y_val*nu**0.1
    y_vals.append(y_val)
    if (y_val/max)>=0.9:
        tenPercentRange.append(nu)
    if (y_val/max)>=0.99:
        onePercentRange.append(nu)
    
print 'done calculation loop, starting to plot'

MAX = np.max(y_vals)
print 'MAX = ',MAX
tenPercentMin = np.min(tenPercentRange)
tenPercentMax = np.max(tenPercentRange)
onePercentMin = np.min(onePercentRange)
onePercentMax = np.max(onePercentRange)
print 'onePercent range [ '+str(onePercentMin)+', '+str(onePercentMax)+' ]'
print 'tenPercent range [ '+str(tenPercentMin)+', '+str(tenPercentMax)+' ]'
ax.plot(NUs,y_vals)
#ax.axes.set_ylim((-0.1e8,0.1e8))
#ax.axes.set_xlim((4.85e13,4.89e13))
ax.set_xscale('log')
#plt.show()
plt.close()


