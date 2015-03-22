#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import re
import glob
import os
import numpy as np
import pylab
import math
from math import pi
plt = pylab.matplotlib.pyplot

def giveMeK():
    period = 11.8618
    M1 = 1.0
    M2 = 0.000954265748
    e = 0.048775
    i = 90
    # convert units of years to seconds
    period_seconds = period*31557600.0
    M1_kg = M1*1.98892e30
    M2_kg = M2*1.98892e30
    G = 6.67300e-11
    
    # Calc K in parts to see equation better
    A = ((2.0*pi*G*(M1_kg+M2_kg))/period_seconds)**(1.0/3.0)
    B = M2_kg/M1_kg
    C = math.sin(math.radians(i))/math.sqrt(1.0-e**2.0)
    # put it all together
    K = A*B*C
    print 'K1 = ',K
    # Calc K in parts to see equation better
    A = ((2.0*pi*G)/period_seconds)**(1.0/3.0)
    B = M2_kg*((M1_kg+M2_kg)**(-2.0/3.0))
    C = math.sin(math.radians(i))#/math.sqrt(1.0-e**2.0)
    # put it all together
    K = A*B*C
    print 'K2 = ',K
    # Calc K in parts to see equation better
    A = ((2.0*pi*G)/period_seconds)**(1.0/3.0)
    B = M2_kg*((M1_kg+M2_kg)**(-2.0/3.0))
    C = math.sin(math.radians(i))/math.sqrt(1.0-e**2.0)
    # put it all together
    K = A*B*C
    print 'K3 = ',K
    
def epochConverter():
    epochsIn = [
                
   2454459.605807,
    2454459.621927,
    2454459.629827,
    2454459.637718,
    2454459.645611,
    2454459.653504,
    2454459.6614,
    2454459.669293,
    2454459.677194,
    2454459.685096,
    2454459.692989,
    2454459.700887,
    2454459.707944






                ]
    P = 21.2165298
    Tc = 2454777.94761
    for t in epochsIn:
        if (abs((t-Tc)/P))>1.0:
            t_out = t-int((t-Tc)/P)*P 
            #print 't = '+str(t)+', int((t-Tc)/P) = '+str(int((t-Tc)/P))+', t_out = '+str(t_out)+', (t_out-Tc) = '+str((t_out-Tc))
            if abs(t_out-Tc)>(P/2.0):
                t_out=t_out+P
        else:
            t_out = t
        print t_out

def omegaTester():
    
    omegasIN = [1,89,91,179,181,269,271,359,361,-1,-89,-91,-179,-181,-269,-271,359,-361]
    for omega in omegasIN:
        top = np.sin(np.radians(omega))
        btm = np.cos(np.radians(omega))
        omegaOut = int(np.degrees(np.arctan2(top,btm)))
        print str(omega)+'  '+str(top)+'  '+str(btm)+'  '+str(omegaOut)
  
def orbitEllipsePlotterDuoOLD(longAN_deg, e, period, inc, argPeri_deg, a, sysDataDict, DIdataDict,\
                            xLabel='?? ["]', yLabel='?? ["]', \
                          plotFilename='', xLim=False, yLim=False, show=True, telescopeView=False ):
    """
    This function will plot the resulting orbit for the input parameters.
    NOTE: If a plotFilename is provided, then the resulting figure will be saved.
    
    IN THIS VERSION the input orbital elements are to be lists, so that multiple orbits are drawn on one 
    plot.  The first orbit provided, that in the 0th element of the lists, will be used for the 
    main orbit which will get its 1/4 orbit sections marked by colored stars and other stuff.
    #******** ALL THESE PARAMS SHOULD BE LISTS OF ELEMENTS ****************
    @param longAN_degs     = Longitude of the Acending Node in degrees
    @type longAN_deg:      = any type of number other than int
    @param argPeri_deg:    = Argument of Periapsis in degrees
    @type argPeri_deg:     = any type of number other than int
    @param a:              = Semi-Major axis in AU
    @type a:               = any type of number other than int
    @param e:              = Eccentricity
    @type e:               = any type of number other than int
    @param inc:              = Inclination in degrees
    @type inc:               = any type of number other than int
    @param period:         = period of orbits [yrs]
    @type period:          = float
    #**********************************************************************
    @param Sys_Dist_PC:    = Distance to the system from Earth [parsec]
    @type Sys_Dist_PC:     = float
    @param xLim:           = range to limit the x axis values to
    @type xLim:            = tuple of two numbers of any type, (min,max)
                             default False indicates to use min and max of 
                             X values for points on the ellipse
                             +5% for white space padding
    @param yLim:           = range to limit the y axis values to
    @type yLim:            = tuple of two numbers of any type, (min,max)
                             default False indicates to use min and max of 
                             Y values for points on the ellipse 
                             +5% for white space padding
    @param PAs:            = Position Angles in [degrees]
    @type PAs:             = list of numbers, same length as SAs
    @param SAs:            = Separation Angles in ["]
    @type SAs:             = list of numbers, same lenght as PAs
    @param sys_dist:       = Distance to System in [pc]
    @type sys_dist:        = any type of number other than int
    """
    
    mas = False
    if mas:
        asConversion = 1000.0
    else:
        asConversion = 1.0
    
    # DI data
    SAs = DIdataDict['SAs']                                   
    PAs = DIdataDict['PAs']
    SAerrors = DIdataDict['SA_errors']                                   
    PAerrors = DIdataDict['PA_errors']
    # General System Data
    sys_dist = sysDataDict['Sys_Dist_PC']
    
    if plotFilename!='':
        datadir = os.path.dirname(plotFilename)
        logFilename = os.path.join(datadir,'processManagerLogFile.txt')
        if plotFilename[-4:]!='.png':
            plotFilename = plotFilename+'.png'  
    else:
        datadir = os.path.curdir
        logFilename = os.path.join(datadir,'TEMPlogFile.txt')   
    # create a plot title from the plotFilename
    plotFileTitle = os.path.basename(plotFilename).split('.')[0]
    
    # figure out logfilename and open log
    log = open(logFilename,'a')
    log.write('\n'+75*'#'+'\n Inside orbitEllipsePlotterDuo \n'+75*'#'+'\n')
    
    #check the orbit element inputs to ensure they are lists, else make them lists
    if type(longAN_deg)!=list:
       longAN_deg = [longAN_deg]
    if type(argPeri_deg)!=list:
       argPeri_deg = [argPeri_deg]
    if type(a)!=list:
       a = [a]
    if type(e)!=list:
       e = [e]
    if type(inc)!=list:
       inc = [inc]
    if type(period)!=list:
       period = [period]
       
    # First convert the given angles to radians
    i_rad = np.radians(inc[0])
    longAN_rad = np.radians(longAN_deg[0])
    argPeri_rad = np.radians(argPeri_deg[0])
    
    # must correct argPeri for the inclination
    argPeri_corr = np.arctan((np.sin(argPeri_rad)*np.cos(i_rad))/np.cos(argPeri_rad))
    argPeri_corr_deg = np.degrees(argPeri_corr)
#    print 'argPeri_deg = '+str(argPeri_deg)
#    print 'argPeri_corr_deg = '+str(argPeri_corr_deg)
#    print 'ang orig = '+str(longAN_deg+argPeri_deg+90.0)
#    print 'ang corr = '+str(longAN_deg+argPeri_corr_deg+90.0)
#    print 'inclination = '+str(inc)

    # angle for ellipse plotter
    ang = longAN_deg[0]+argPeri_deg[0]+90.0
    ang_corr = longAN_deg[0]+argPeri_corr_deg+90.0

     # semi-major in orbital plane
    a_op = a[0]
    
    # semi-minor in orbital plane
#    b_op = a[0]*(np.sqrt(1-(e[0]*e[0])))
    
    # semi-major in reference plane
    a_rp_y = a_op*np.sin(argPeri_rad)*np.cos(i_rad)
    a_rp_x = a_op*np.cos(argPeri_rad)
    a_rp = np.sqrt(np.power(a_rp_x, 2)+np.power(a_rp_y, 2))
    
#    print 'a_rp_x = '+str(a_rp_x)
#    print 'a_rp_y = '+str(a_rp_y)
#    print 'a_rp = '+str(a_rp)
#    
    # semi-minor in reference plane 
    # Note: semi-minor and semi-major are at 90deg to each other.
#    b_rp_y = b_op*np.sin((pi/2.0)-argPeri_rad)*np.cos(i_rad)
#    b_rp_x = b_op*np.cos((pi/2.0)-argPeri_rad)
#    b_rp = np.sqrt(np.power(b_rp_x, 2)+np.power(b_rp_y, 2))
    
#    print 'b_op = '+str(b_op)
#    print 'b_rp_x = '+str(b_rp_x)
#    print 'b_rp_y = '+str(b_rp_y)
#    print 'b_rp = '+str(b_rp)

    # Calculate distance of center-foci assuming center is at (0,0)
    c_foci = a_rp-((a_rp*(1.0-e[0]*e[0]))/(1.0-e[0]))
    # Calculate loction of foci where star would lie
    yStar = -c_foci*np.sin(np.radians(ang_corr))*(asConversion/sys_dist)
    xStar = -c_foci*np.cos(np.radians(ang_corr))*(asConversion/sys_dist)
    
    if telescopeView:
        yStar = -yStar
        xStar = -xStar
#    print 'c_foci = '+str(c_foci)
#    print 'yStar = '+str(yStar)
#    print 'xStar = '+str(xStar) 
        
    # create a primary star polygon
    starPolygon = star((asConversion/1000.0)*2.0*a[0], 0, 0, color='black', N=5, thin = 0.5)
    
    ## calculate the locations of companion for 'numOrbs' locations throughout the orbit to make an orbit ellipse    
    ellipseXs2 = []
    ellipseYs2 = []
    #orbitTAs = []
    #orbitSAs = []
    #orbitPAs = []
    #sep_dists = []
    for orb in range(0,len(longAN_deg)):
        ellipseXs = []
        ellipseYs = []
        numSteps = 5000.0
        periodIncrement = (period[orb]*365.25)/numSteps
        t = 1.0 
        for step in range(0,int(numSteps)):
            T = 0.0
            t = t + periodIncrement
            (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SA, PA, a1, a2) =\
                orbitCalculator2(t, sys_dist, inc[orb], longAN_deg[orb], e[orb], T, period[orb], argPeri_deg[orb], a[orb],\
                                                            Mass1=1, Mass2=1, verbose=False)
            #(PA, SA) = PASAcalculator(period, t, T, e, inc, longAN_deg, argPeri_deg, sys_dist, a, a1=0, verbose=False)
            #orbitTAs.append(TA_deg)
            #orbitPAs.append(PA)
            #orbitSAs.append(SA)
            ellipseX = SA*math.sin(math.radians(PA))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
            ellipseY = -SA*math.cos(math.radians(PA))*asConversion#*sys_dist   # This will be in [mas] instead of [AU]
            if telescopeView:
                ellipseX = -ellipseX
                ellipseY = -ellipseY
            ellipseXs.append(ellipseX)
            ellipseYs.append(ellipseY)
            #sep_dist = math.sqrt(math.pow(ellipseX,2.0)+math.pow(ellipseY,2.0))
            #sep_dists.append(sep_dist)
        ellipseXs2.append(ellipseXs)
        ellipseYs2.append(ellipseYs)
        
    ## Get the locations of 500 points on an ellipse representing the orbit # this isn't working right... thus orbit method above is used now.
    #X,Y = ellipse(a_rp, b_rp, ang_corr, -xStar, -yStar, Nb=500)
    #X,Y = ellipse(a_rp, b_rp, ang, 0, 0, Nb=500)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, a1, a2) =\
            orbitCalculator2(0.0, sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    Xstart = SAstart*math.sin(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    Ystart = -SAstart*math.cos(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    if telescopeView:
        Xstart = -Xstart
        Ystart = -Ystart
    startStar = star((asConversion/1000.0)*0.6*a[0], Xstart, Ystart, color='green', N=5, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, a1, a2) =\
            orbitCalculator2(((period[0]/4.0)*365.25), sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    XoneQuarter = SAstart*math.sin(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    YoneQuarter = -SAstart*math.cos(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    if telescopeView:
        XoneQuarter = -XoneQuarter
        YoneQuarter = -YoneQuarter
    oneQuarterStar = star((asConversion/1000.0)*0.6*a[0], XoneQuarter, YoneQuarter, color='blue', N=5, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAhalf, PAhalf, a1, a2) =\
            orbitCalculator2(((period[0]/2.0)*365.25), sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    Xhalf = SAhalf*math.sin(math.radians(PAhalf))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    Yhalf = -SAhalf*math.cos(math.radians(PAhalf))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    if telescopeView:
        Xhalf = -Xhalf
        Yhalf = -Yhalf
    halfStar = star((asConversion/1000.0)*0.6*a[0], Xhalf, Yhalf, color='yellow', N=5, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, a1, a2) =\
            orbitCalculator2((3.0*(period[0]/4.0)*365.25), sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    XthreeQuarter = SAstart*math.sin(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    YthreeQuarter = -SAstart*math.cos(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    if telescopeView:
        XthreeQuarter = -XthreeQuarter
        YthreeQuarter = -YthreeQuarter
    threeQuarterStar = star((asConversion/1000.0)*0.6*a[0], XthreeQuarter, YthreeQuarter, color='orange', N=5, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, a1, a2) =\
            orbitCalculator2(((period[0])*365.25), sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    Xend = SAstart*math.sin(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    Yend = -SAstart*math.cos(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    if telescopeView:
        Xend = -Xend
        Yend = -Yend
    endStar = star((asConversion/1000.0)*0.6*a[0], Xend, Yend, color='purple', N=5, thin = 0.5)
    
    # set the xLim and yLim if their values are False
    # and pad max and min values found by 10%
    if not xLim:
        min = np.min(ellipseXs2[:])
        max = np.max(ellipseXs2[:])
        Range = abs(max)+abs(min)
        xLim = (min-abs(Range*0.05),max+abs(Range*0.05))
    else:
        if not (type(xLim)==tuple):
            s= 'PROBLEM: xLim is not of type tuple'
            print s
            log.write(s+'\n')
    if not yLim:
        min = np.min(ellipseYs2[:])
        max = np.max(ellipseYs2[:])
        Range = abs(max)+abs(min)
        yLim = (min-abs(Range*0.05),max+abs(Range*0.05))
    else:
        if not (type(yLim)==tuple):
            s= 'PROBLEM: yLim is not of type tuple'
            print s
            log.write(s+'\n')
            
    ## Create figure, axes and start plotting/drawing everything
    fig = plt.figure(1,figsize=(10,10))
    main = fig.add_subplot(111)
    main.axes.set_xlim(xLim)
    main.axes.set_ylim(yLim)
    main.set_xlabel(xLabel)
    main.set_ylabel(yLabel)
    main.set_title(plotFileTitle)
        
    # Draw orbits
    for orb in range(0,len(longAN_deg)):
        main.plot(ellipseXs2[orb],ellipseYs2[orb]) #$$$$ add title, labels, axes
    # Draw Ellipse ## found this didn't work and thus the orbit method above is used now.
    #main.plot(X,Y, c='black')
    
    #draw semi-major
    main.plot([Xstart,Xhalf],[Ystart,Yhalf],'g-')
    
    #draw semi-minor 
    #main.plot(X_b_rp_corr,Y_b_rp_corr,'b-')
    
    #draw 'x' for center of ellipse
    main.plot(xStar,yStar,'rx')
    
    # Draw lines for horizontal and vertical from origin
    main.plot([xLim[0],xLim[1]],[0,0],c='black')
    main.plot([0,0],[yLim[0],yLim[1]],c='black')
    
    # draw lines along diagonals to mare out 45 degree locations from origin ### Useless now
    #main.plot( [xLim[0],xLim[1]], [yLim[0],yLim[1]])
    #main.plot( [xLim[1],xLim[0]], [yLim[0],yLim[1]])
    
    # Draw stars for the location of each 1/4 of the orbit  # These are optional, I would not include them for the final versions, maybe just the periapsis one.
    main.add_patch(startStar)
    main.add_patch(oneQuarterStar)
    main.add_patch(halfStar)
    main.add_patch(threeQuarterStar)
    main.add_patch(endStar)
    
    # Draw larger star for primary star's location
    main.add_patch(starPolygon)
    
    ## call function to calculate, create and return polygons for the 
    ## companion star locations and boxes for their errors
    (errorBoxes, m2starPolygons) = starAndErrorPolysOLD(SAs,SAerrors,PAs,PAerrors,asConversion, main.transData, telescopeView)
    
    # draw the error boxes for the companion start locations in the data
    for errorBox in errorBoxes:
        main.add_patch(errorBox)
        
    # Draw red star patches for each of the companion's locations from the data
    for star2 in m2starPolygons:
        main.add_patch(star2)
        
    # add a legend
    #main.legend(('longAN_deg = '+str(longAN_deg),'e = '+str(e), 'period = '+str(period), 'inc = '+str(inc), 'argPeri_deg = '+str(argPeri_deg), 'a = '+str(a)), loc=0, markerscale=0.0000000000000001)
    paramsLegndStr = 'longAN_deg = '+str(longAN_deg[0])+'\ne = '+str(e[0])+'\nperiod = '+str(period[0])+'\ninc = '+str(inc[0])+'\nargPeri_deg = '+str(argPeri_deg[0])+'\na = '+str(a[0])
    main.text(xLim[0]+abs(xLim[0]*0.02),abs(yLim[1]*0.2),paramsLegndStr,ha='left')
    
    # save plot to file
    if plotFilename!='':
        plt.savefig(plotFilename, dpi=300, orientation='landscape')
    
    # show plot
    if show:
        plt.show()
    
    # close figure/plot
    plt.close()    
    log.write('\n'+75*'#'+'\n Leaving orbitEllipsePlotterDuo \n'+75*'#'+'\n')
    log.close()
    
    
def summaryPlotter(plotFileTitle, chiSquareds, inclination_degsAlls, longAN_degsAlls, argPeri_degsAlls, esAlls,\
            periodsAlls, TsAlls, a1Means, a2Means, confLevels=True, weight=True, normed=True, showPlots=False,\
                                                                                     save=False, verbose=False):
    """
    OLD VERSION!!
    
    This advanced plotting function will plot all the data in a 3x3 grid on a single figure.  The data will be plotted
    in histograms that can be will be normalized to a max of 1.0 and the data can be weighted if desired.  The 
    confidence levels of the data can be calculated and the resulting histograms will have the bars inside the
    68% bin shaded as dark grey and 95% as light grey and the rest will be white.
    """
            
    # check if the passed in value for plotFilename includes '.png'
    if '.png' not in plotFileTitle:
        plotFilename = plotFileTitle+'.png'
    else:
        plotFilename = plotFileTitle
    
    # create a plot title from the filename
    plotFileTitle = os.path.basename(plotFilename).split('.')[0]
    
    if verbose:
        print 'Starting to make Total Summary Plot'
    if len(inclination_degsAlls)>0:
        fig = plt.figure(1, figsize=(25,22) ,dpi=300)
        #plt.title(plotFileTitle+' TOTAL Summary')
        
        plot = fig.add_subplot(331)
        data = inclination_degsAlls
        xlabel = 'Inclination [deg]'
        plot.axes.set_ylabel('Probability')
        plot = histConverter(chiSquareds, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
        incMedian = np.median(inclination_degsAlls)
        numOrbits = len(inclination_degsAlls)
        ##################$$$$$$$$$$$$$ This extra garbage collection might not be needed but I want it for now as a code EX. ######
        #del inclination_degsAlls
        #gc.collect()
        print "done plotting inclination_degsAlls"
        
        plot = fig.add_subplot(332)
        data = longAN_degsAlls
        xlabel = 'Longitude of Ascending Node [deg]'
        plt.title(plotFileTitle+'   TOTAL Summary')
        plot = histConverter(chiSquareds, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
        longANMedian = np.median(longAN_degsAlls)
        print "done plotting longAN_degsAlls"
        
        plot = fig.add_subplot(333)
        data = argPeri_degsAlls
        xlabel = 'Argument of Perigie [deg]'
        plot = histConverter(chiSquareds, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
        argPeriMedian = np.median(argPeri_degsAlls)
        print "done plotting argPeri_degsAlls"
        
        plot = fig.add_subplot(334)
        data = esAlls
        xlabel = 'e'
        plot.axes.set_ylabel('Probability')
        plot = histConverter(chiSquareds, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
        eMedian = np.median(esAlls)
        print "done plotting esAlls"
        
        plot = fig.add_subplot(335)
        data = periodsAlls
        xlabel = 'Period [Years]'
        plot = histConverter(chiSquareds, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
        periodMedian = np.median(periodsAlls)
        print "done plotting periodsAlls"
        
        plot = fig.add_subplot(336)
        Tsmin = int(np.min(TsAlls))
        Tsmax = int(np.max(TsAlls))
        TsNEW = TsAlls
        for i in range(0,len(TsNEW)):
            TsNEW[i]=(TsNEW[i]-Tsmin)/365
        data = TsNEW
        xlabel = 'Time of last Periapsis [years since '+str(Tsmin)+']'
        plot = histConverter(chiSquareds, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
        TMedian = np.median(TsAlls)
        print "done plotting TsAlls"
        
        # handle both the case where there are usefull a1 values and not
        if (a1Means[0]==0.0) and (a1Means[1]==0.0):
            # no useful a1's, ie single body orbit
            plot = fig.add_subplot(337)
            data = a2Means
            xlabel = 'Semi-major Axis [AU]'
            plot.axes.set_ylabel('Probability')
            plot = histConverter(chiSquareds, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
            a2Median = np.median(a2Means)
            print "done plotting inclinations"
        else:
            # useful a1 values, ie two orbits
            # plot a1's
            plot = fig.add_subplot(337)
            data = a1Means
            xlabel = 'Semi-major Axis body 1 [AU]'
            plot.axes.set_ylabel('Probability')
            plot = histConverter(chiSquareds, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
            # plot a2's
            plot = fig.add_subplot(338)
            data = a2Means
            xlabel = 'Semi-major Axis body 2 [AU]'
            plot = histConverter(chiSquareds, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
            
            # Also create a total semi-major axis plot too (completes 3x3 grid)
            aALLmeans = []
            for i in range(0,len(a1Means)):
                aALLmeans.append(a1Means[i]+a2Means[i])
            plot = fig.add_subplot(339)
            data = aALLmeans
            xlabel = 'Total Semi-major Axis [AU]'
            plot = histConverter(chiSquareds, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
            a2Median = np.median(a2Means)
            a1Median = np.median(a1Means)
            aALLMedian = np.median(aALLmeans)
            print "done plotting aALLmeans"
        # Save file if requested.
        if save:
            plt.savefig(plotFilename, dpi=300, orientation='landscape')
        # print median values of data for reference if verbose
        if verbose:
            print 'Finished making Total Summary Plot'
            print str(numOrbits)+' orbits used'
            print 'median inclination_degsAlls =  ',incMedian
            print 'median longAN_degsAlls = ', longANMedian
            print 'median argPeri_degsAlls = ', argPeriMedian
            print 'median esAlls = ',eMedian
            print 'median periodsAlls =  ',periodMedian
            print 'median TsAlls [jd] = ', TMedian
            if (a1Means[0]==0.0) and (a1Means[1]==0.0):
                print 'median a2Means = ',a2Median
            else:
                print 'median a1Means = ', a1Median
                print 'median a2Means = ',a2Median
                print 'median aALLmeans = ',aALLMedian
    else:
        if verbose:
            print 'No Gooders to plot in Total Summary !!'    
    # Show the plots on the screen if asked to    
    if showPlots:
        plt.show()
        plt.close()
    else:
        plt.close()
        
def extraTestPlotter(datafilename):
    """
    a temp plotter for plotting a gauss dist produced by c++
    """
    
    f = open(datafilename,'r')
    lines = f.readlines()
    print str(len(lines))+' lines were found in file'
    f.close()
    data = []
    for line in lines:
        data.append(float(line))
    print str(len(data))+" lines were moved into the data ary"
    fig = plt.figure(1)
    subPlot = fig.add_subplot(111)
    (a,b,c)=subPlot.hist(data, bins=50, normed=True)

    plotFilename = datafilename[:-4]+"_PLOT.png"
    plt.savefig(plotFilename, dpi=300, orientation='landscape')
    print 'Summary plot saved to: '+plotFilename
    plt.show()
    plt.close()
    
def dataWriter(plotFileTitle, longAN_degs, es, Ts, a1s, a2s, periods, inclination_degs, argPeri_degs,\
                  Sep_Dists2, thetas2, Es2, Ms2, ns2, chiSquareds):
    
    numSamples = np.shape(Es2[:][:])[0]
    numEpochs = np.shape(Es2[:][:])[1]
    
    filename = '../data/'+plotFileTitle+'_INS'
    f = open('%s.txt' % filename, 'w')
    f.write(plotFileTitle+'_INS'+'\n')
    f.write('longAN [deg]   es [N/A]   Ts [julian date]   a1s [AU]   a2s [AU]   periods [yrs]   inclination [deg]   argPeri [deg]\n')
    for sample in range(0,numSamples):
        line = str(longAN_degs[sample])
        line = line +' '+ str(es[sample])
        line = line +' '+ str(Ts[sample])
        line = line +' '+ str(a1s[sample])
        line = line +' '+ str(a2s[sample])
        line = line +' '+ str(periods[sample])
        line = line +' '+ str(inclination_degs[sample])
        line = line +' '+ str(argPeri_degs[sample])+'\n'
        f.write(line)
    f.close()
    
    ################# OUTS ####################
    for epoch in range(0,numEpochs):
        filename = '../data/'+plotFileTitle+'_OUTSepoch'+str(epoch+1)
        f = open('%s.txt' % filename, 'w')
        f.write(plotFileTitle+'_OUTSepoch'+str(epoch+1)+'\n')
        f.write('Sep_Dists [AU]    thetas [deg]     Es [deg]     Ms [deg]     ns [rad/yr]\n ')
        for sample in range(0,numSamples):
            line = str(Sep_Dists2[sample][epoch])
            line = line +' '+ str(thetas2[sample][epoch])
            line = line +' '+ str(Es2[sample][epoch])
            line = line +' '+ str(Ms2[sample][epoch])
            line = line +' '+ str(ns2[sample][epoch])+'\n'
            f.write(line)
        f.close()
        
def fileExtRename(relDir,inExt,outExt):
    """
    inEXT should be like ".txt"
    outEXT should be the same, ".dat"
    """
    verbose = True
    
    # get data list and sort it
    if relDir[-1]!='/':
        relDir = relDir+'/'
    framelist = np.sort(glob.glob(relDir + "*"+inExt))
    nframes = len(framelist)
    
    print "In dir: "+relDir+", has "+str(nframes)+" files with the original ext "+inExt
    
    for frame in framelist:
        frameIN = frame
        frameOUT = re.sub(inExt,outExt,frameIN)
        
        if verbose:
            print 'Input file found: '+os.path.basename(frameIN)
            print 'File will be renamed '+os.path.basename(frameOUT)
        if os.path.exists(frameOUT)==False:
            os.rename(frameIN, frameOUT)
            if os.path.exists(frameOUT)==True:
                if verbose:
                    print 'File renamed successfully'
            else:
                print 'Warning: renaming did not take hold it seems!!'
            
    print 'All files were renamed in the directory.'
    
def confidenceLevelsFinderNEW(filename, verbose=False):
    """
    CAUTION, this was designed to be used with the outputs of 100ModDataset simulation outputs.
    NOT a good function to load data from a long simulation with lots of output sets.
    
    PURPOSE: This is to determine the errors that the mod dataset runs were designed to determine.
    
    NOTE: This version is meant for the new headings from mcmcOrbSimUniform6
        longAN [deg]      e [N/A]       To [julian date]  period [yrs]   inclination [deg]   argPeri [deg]   a_total [AU]    chiSquared  RVoffset0...  timesBeenHere
        
        file format must be:
        
        line1: title
        line2: data headings
        line3: space delimited data
        line4: space delimited data
    
    """
    
    outDict = outputDatafileToDict(filename)
    
    chiSquareds = outDict["chiSquareds"]
    longAN_degsCLevels = ConfLevelFunc(chiSquareds,outDict["longANs"])
    print '\nlongANs have conf levels: \n68.3% = '+repr(longAN_degsCLevels[0])+'\n95.5% = '+repr(longAN_degsCLevels[1])+'\n'
    esCLevels = ConfLevelFunc(chiSquareds,outDict["es"])
    print 'es have conf levels: \n68.3% = '+repr(esCLevels[0])+' \n95.5% = '+repr(esCLevels[1])+'\n'
    TsCLevels = ConfLevelFunc(chiSquareds,outDict["Ts"])
    print 'Ts have conf levels: \n68.3% = '+repr(TsCLevels[0])+' \n95.5% = '+repr(TsCLevels[1])+'\n'
    periodsCLevels = ConfLevelFunc(chiSquareds,outDict["periods"])
    print 'periods have conf levels: \n68.3% = '+repr(periodsCLevels[0])+' \n95.5% = '+repr(periodsCLevels[1])+'\n'
    inclination_degsCLevels= ConfLevelFunc(chiSquareds,outDict["inclinations"])
    print 'inclinations have conf levels: \n68.3% = '+repr(inclination_degsCLevels[0])+' \n95.5% = '+repr(inclination_degsCLevels[1])+'\n'
    argPeri_degsCLevels = ConfLevelFunc(chiSquareds,outDict["argPeris"])
    print 'argPeris have conf levels: \n68.3% = '+repr(argPeri_degsCLevels[0])+' \n95.5% = '+repr(argPeri_degsCLevels[1])+'\n'
    asCLevels = ConfLevelFunc(chiSquareds,outDict["a_totals"])
    print 'a_totals have conf levels: \n68.3% = '+repr(asCLevels[0])+' \n95.5% = '+repr(asCLevels[1])+'\n'
    
    if type(outDict["RVoffsets"][0])!=list:
        RVoffsets = [outDict["RVoffsets"]]
    else:
        RVoffsets = outDict["RVoffsets"]
        
    for dataset in range(0,len(RVoffsets)):
        offsetsCurr = RVoffsets[dataset]
        offsetsCurrCLevels = ConfLevelFunc(chiSquareds,offsetsCurr)
        print 'dataset # '+str(dataset+1)+' RV offsets have conf levels: \n68.3% = '+repr(offsetsCurrCLevels[0])+' \n95.5% = '+repr(offsetsCurrCLevels[1])+'\n'
    
def ConfLevelFunc(chiSquareds,param):
    """
    return [[68.3% minimum, 68.3% maximum],[95.5% minimum, 95.5% maximum]]
    """
    chiSquareMin = np.min(chiSquareds)
    
    numSamples = len(chiSquareds)
    
    paramCLevels = []
    
    ## find all the parameter values that are inside the 68.3% confidence level
    param68s = []
    for orbit in range(0,numSamples):
        if chiSquareds[orbit]<=(chiSquareMin+1.0):
            param68s.append(param[orbit])
    if len(param68s)==0:
        param68s.append(0)
    paramCLevels.append([findArrayMin(param68s),findArrayMax(param68s)])
    
    ## find all the parameter values that are inside the 95.5% confidence level
    param95s = []
    for orbit in range(0,numSamples):
        if chiSquareds[orbit]<=(chiSquareMin+4.0):
            param95s.append(param[orbit])
    if len(param95s)==0:
        param95s.append(0)
    paramCLevels.append([findArrayMin(param95s),findArrayMax(param95s)])
    
    return paramCLevels