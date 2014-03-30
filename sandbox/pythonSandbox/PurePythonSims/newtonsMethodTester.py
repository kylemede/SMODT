import math
import numpy as np
import random as rand
import pylab
from pylab import pi
plt = pylab.matplotlib.pyplot



def newton(M_max, e_max):
    
    for i in range(int(1e6)):
        M = rand.uniform(0.0, M_max)
        e = rand.uniform(0.0001, e_max)
        
        M_deg = math.degrees(M)
        ### Performing Newton's Method to get the Eccentric Anomaly value ###
        # initial guess (E_last), will be updated in loop.  
        # Anything works, just takes longer if further from real value. => pi
        E_last = 2*math.pi
        # stored initial value to be updated in loop
        # this value is always very close to the true value and will minimize the number of loops
        E_latest = M+e*math.sin(M) 
        
        # show input value to 
    #    print "Inputs to Newton's Method are : \nM [rad]= "+str(M)+"\nE_last [rad]= "+\
    #    str(E_last)+"\nE_latest [rad] = "+str(E_latest)+"\ne = "+str(e)
    #    
        #print "\nStarting to run Newton's while loop."
        
        count = 0 # a counter to stop inf loops in Newtons method below
        while (abs(E_last-E_latest) > (1.0e-10))and(count<100):
            
            #print 'current E [rad]= ', E_latest
            
            E_last = E_latest
            #E_latest = E_last - ((M-E_last+e*math.sin(E_last))/(e*math.cos(E_last)-1.0))\
            E_latest = E_last - ((E_last-M-e*math.sin(E_last))/(1.0-e*math.cos(E_last)))
            count = count+1
        
        E_latest_deg = math.degrees(E_latest) # convert resulting E to degrees
        
    #    print '-'*25
    #    print "The resultant E value is [deg] = ", E_latest_deg
    #    print 'M [deg] = ',M_deg
    #    print '-'*25
        # check if the resultant value solves the original equation.
        Mnewton = math.degrees(E_latest-e*math.sin(E_latest))
        if abs(M_deg-Mnewton)>(1.0e-5):
            print "PROBLEM: This resultant E does not satisfy the original equation, Newton's method FAILED !!!"
            print 'M from this E Equals = '+str(Mnewton)
            print 'M original = '+str(M_deg)
            print 'E initial = '+str(math.degrees(M+e*math.sin(M) ))
            print 'e = '+str(e)
            
            E = E_latest
            RHS = []
            ang = []
            for i in range(1000):
                ang.append(rand.uniform(M-2*pi,M+2*pi))
                RHS.append(E - e*math.sin(ang[-1]))
                
            fig = plt.figure(1)
            p = fig.add_subplot(111)
            p.axes.set_xlim((M-2*pi,M+2*pi))
            p.axes.set_ylim((np.min(RHS),np.max(RHS)))
            p.scatter(ang,RHS)
            p.plot([M-2*pi,M+2*pi],[M,M])
            p.plot([E,E],[np.min(RHS),np.max(RHS)], color='red')
            plt.show()
    #    else:
    #        print "This resultant E solves the original equation, Newton's Method worked :-)"
            ### Newton's loop finished! ###
    print 'ALL random Ms tried and passed'   
     
#    E = []
#    RHS = []    
#    for i in range(100):
#            E.append(rand.uniform(0,2*pi))
#            RHS.append(E[-1]-e*math.sin(E[-1]))
    
#    E = E_latest
#    RHS = []
#    ang = []
#    for i in range(1000):
#        ang.append(rand.uniform(M-2*pi,M+2*pi))
#        RHS.append(E - e*math.sin(ang[-1]))
#        
#    fig = plt.figure(1)
#    plt.subplot(111)
#    plt.scatter(ang,RHS)
#    plt.plot([M-2*pi,M+2*pi],[M,M])
#    plt.plot([E_latest,E_latest],[np.min(RHS),np.max(RHS)], color='red')
#    plt.show()