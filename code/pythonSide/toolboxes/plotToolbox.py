#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import math
import gc
import numpy as np
import os
import pylab
import timeit
from math import pi
import generalToolbox as genTools 
import DItoolbox as diTools
import RVtoolbox as rvTools
plt = pylab.matplotlib.pyplot
patches = pylab.matplotlib.patches


def ellipse(ra,rb,ang,x0,y0,Nb=50):
    '''
    For calculating Nb x and y values that lie on an ellipse 
    defined using ra, rb, ang, xo and yo.
    
    :param ra:          Semi-Major axis length
    :type ra:           any number
    :param rb:          Semi-Minor axis length
    :type rb:           any number
    :param ang:         Angle [degrees] to plot the semi-major axis 
        wrt horizontal, measured from positive x-axis
    :type ang:          any number
    :param x0:          X coordinate of center of ellipse
    :type xo:           any number
    :param y0:          Y coordinate of center of ellipse
    :type y0:           any number
    :param Nb:          Number of x,y points to return 
    :type Nb:           int
    
    NOTE:
    based on matlab code ellipse.m written by D.G. Long,
    Brigham Young University, based on the
    CIRCLES.m original
    written by Peter Blattner, Institute of Microtechnology,
    University of
    Neuchatel, Switzerland, blattner@imt.unine.ch
    '''
    xpos = x0
    ypos = y0
    radm = ra
    radn = rb
    an = np.radians(ang)
    
    co = np.cos(an)
    si = np.sin(an)
    the=linspace(0,2*pi,Nb)
    X=radm*np.cos(the)*co-si*radn*np.sin(the)+xpos
    Y=radm*np.cos(the)*si+co*radn*np.sin(the)+ypos
    return X,Y    

def fixPlotBordersAndLabels(plot):
    """
    This will just increase the boarder thicknesses and label sizes to look better for printing or 
    putting in papers.
    """
    plot.tick_params(axis='both',which='major',labelsize=20)
    #plot.axhline(linewidth=2.0)
    #plot.axvline(linewidth=2.0)
    plot.xaxis.set_tick_params(width=3)
    plot.yaxis.set_tick_params(width=3)
    plot.spines['right'].set_linewidth(2.0)
    plot.spines['bottom'].set_linewidth(2.0)
    plot.spines['top'].set_linewidth(2.0)
    plot.spines['left'].set_linewidth(2.0)
    return plot
    
def histConverter(chiSquareds, data, plot, xlabel, confLevels=False, weight=True, normed=True, nu=1, bestVal=0,logY=False):
    """
    This function is for creating a modified histogram plot to be properly normalized to a max value
    of 1.0.  There is also the option to color the histogram bars according to the 95 and 68% confidence
    levels of the data.
    """
    verbose = False
    debug = False
    
    if debug:
        print 'in histConverter'
    if (type(data)==bool):
        if verbose:
            print 'input data was a bool'
        data=data0=dataMid=dataBack=0.0
    elif (type(data)==float):
        if verbose:
            print 'input data was a float'
        data0=dataMid=dataBack=data
    else:
        if verbose:
            print ' input data was not a list or bool'
        #Check if data has any varying values
        data0 = data[0]
        dataMid = data[len(data)//2]
        dataBack = data[-1]
        
    if debug:
        print 'data0 = ',data0
        print 'dataMid = ',dataMid
        print 'dataBack = ',dataBack
        print '(data0!=dataBack) = '+repr((data0!=dataBack))
        print '(data0!=dataMid) = '+repr((data0!=dataMid))
        
    # if data varies, then make proper histogram
    if (data0!=dataBack)and(data0!=dataMid):
        if debug:
            print 'inside section that will make the hist plots'
        # if requested, calculate the likelihoods which will be the weights for the data points
        if weight:
            theWeights = likelihoodsCalc(chiSquareds, nu)
        else:
            theWeights = np.ones(len(data))
        # make initial histogram and update it below
        if verbose:
            print 'making initial hist, then convert it'
        (n,bins,rectangles)=plot.hist(data, bins=50, normed=normed, weights=theWeights,log=logY)#, fill=False)
        if type(confLevels)==list:
            if verbose:
                print 'changing rectangle colors to match 68confLevels: '+repr(confLevels[0])+', and 95confLevels: '+repr(confLevels[1])
            # Update rectangle patche's colors according to Confidence Levels of the data
            recs2 = []
            for rec in rectangles:
                    x=rec.get_x()
                    c = 'w'
                    if(x>confLevels[1][0])and(x<confLevels[1][1]):
                        c = '0.8'
                    if (x>confLevels[0][0])and(x<confLevels[0][1]):
                        c = '0.5'
                    recs2.append(patches.Rectangle(xy=rec.get_xy(), width=rec.get_width(),height=rec.get_height(),facecolor=c, edgecolor='k'))#color=c))
            # draw updated patches on plot
            for rec in recs2:
                    plot.add_patch(rec)
        if normed:
            #update the y limit and its ticks and tick labels
            if not logY:
                plot.axes.set_ylim([0.0,plot.axes.get_ylim()[1]*1.2])
                locs = []
                lim = plot.axes.get_ylim()[1]
                for i in range(0,6):
                    locs.append((lim/6.0)*i)
                plot.axes.set_yticks(locs)
                plot.axes.set_yticklabels([' ','0.2','0.4','0.6','0.8','1.0'])
        # draw a line up the median of the data
        Median = np.median(data)
        plot.plot([Median, Median],plot.axes.get_ylim(),'k',linewidth=3)
        #convert data and chiSquareds arrays to np arrays if needed
        if type(chiSquareds)!=np.ndarray:
            chiSquareds = np.array(chiSquareds)
        if type(data)!=np.ndarray:
            data = np.array(data)
        if bestVal==0:
            bestVal = data[np.where(chiSquareds==chiSquareds.min())][0]
            if debug:
                print "using calculated bestVal = ",bestVal
        else:
            if debug:
                print 'using provided bestVal = ',bestVal
        plot.plot([bestVal, bestVal],plot.axes.get_ylim(),'g',linewidth=3)
        
        if verbose:
            print "min found to be "+str(np.min(data,axis=0))+", max was "+str(np.max(data,axis=0))
            print "Best value was found to be at "+str(bestVal)+", and median at "+str(Median)
    else:
        if debug:
            print 'data vals all equal so making empty plot for kicks'
        if (type(data)!=float):
            data = data[0]
        # This means the data had no useful values in it, just constant 
        plot.plot([data,data],[0.0,1.2],'b',linewidth=2)
        if debug:
            print 'empty plot made'
        plot.axes.set_ylim([0.0,1.2])
        if debug:
            print 'axes limits set'
        locs = []
        lim = plot.axes.get_ylim()[1]
        if debug:
            print 'limits are '+repr(lim)
        for i in range(0,6):
            locs.append((lim/6.0)*i)
        if debug:
            print 'locs found '+repr(locs)
        plot.axes.set_yticks(locs)
        if debug:
            print 'locs set as y ticks'
        plot.axes.set_yticklabels(['0.0','0.2','0.4','0.6','0.8','1.0'])
        print 'Values for '+xlabel+' were found to be constant, so not making a histogram, just a simple line plot!!'
    # add x label
    plot.axes.set_xlabel(xlabel,fontsize=30)
    plot = fixPlotBordersAndLabels(plot)
    #plot.tick_params(axis='both',which='major',labelsize=20)
    
    return plot

def histAndScatterPlotter(xData, yData, xLabel='', yLabel='', plotFilename='', xLim=False, yLim=False, show=True, save=False):
    """
    This function will create a plot with a scatter plot of the xData and yData in the middle, a horizontal histogram 
    of the yData on the right, and vertical histogram of the xData at the top.  This is a good way to compare the 
    parameter space where the two data sets have their values focus.
    
    NOTE: xData and yData must have the same length.
    
    :param xData:          Data to plot on the x axis
    :type xData:           list of numbers of any type, must be same length as yData
    :param yData:          Data to plot on the y axis
    :type yData:           list of numbers of any type, must be same length as xData
    :param xLabel:         label for the x axis data
    :type xLabel:          string, default is ''
    :param yLabel:         label for the y axis data
    :type yLabel:          string, default is ''
    :param plotFilename:   name of the file to save the figure to
    :type plotFilename:    string, ensuring to include directory path if not cwd
    :param xLim:           range to limit the x axis values to
    :type xLim:            tuple of two numbers of any type, (min,max)
        default False indicates to use min and max of xData
    :param yLim:           range to limit the y axis values to
    :type yLim:            tuple of two numbers of any type, (min,max)
        default False indicates to use min and max of yData
    :param show:           show the resulting figure on the screen?
    :type show:            boolean, default is False
    :param save:           save the resulting figure?
    :type save:            boolean, default is True
    """
    
    # check the data sets have matching lengths
    if not (len(xData)==len(yData)):
        print 'PROBLEM: The x and y data sets do not have matching lengths.'    
        
    # set the xLim and yLim if their values are False
    if not xLim:
        min = np.min(xData)
        max = np.max(xData)
        range = abs(max)+abs(min)
        xLim = (min-abs(range*0.05),max+abs(range*0.05))
    else:
        if not (type(xLim)==tuple):
            print 'PROBLEM: xLim is not of type tuple'
    if not yLim:
        min = np.min(yData)
        max = np.max(yData)
        range = abs(max) + abs(min)
        yLim = (min-abs(range*0.05),max+abs(range*0.05))
    else:
        if not (type(yLim)==tuple):
            print 'PROBLEM: yLim is not of type tuple'        
    
    # Check the plotFilename has .png, else add it
    if save:
        if plotFilename[-4:]!='.png':
            plotFilename = plotFilename+'.png'       
        
    # create figure to put plots in
    fig = plt.figure(1,figsize=(10,10) ,dpi=300)
    
    # set up and plot central scatter plot
    main = fig.add_subplot(221)
    main.set_position([0.2,0.2,0.4,0.4])
    main.axes.set_xlim(xLim)
    main.axes.set_ylim(yLim)
    main.axes.set_xlabel(xLabel)
    main.axes.set_ylabel(yLabel)
    main.scatter(xData, yData, s=2, c='black', edgecolors='none')
    
    # set up and plot top X axis histogram
    X = fig.add_subplot(223)
    X.set_position([0.2,0.6,0.4,0.2])
    X.xaxis.set_ticks_position('top')
    X.axes.set_xlim(xLim)
    X.yaxis.set_ticklabels([]) #this is just a hack way of killing the tick labels
    X.axes.set_xlabel(xLabel)
    X.xaxis.set_label_position('top')
    #X.set_title(os.path.basename(plotFilename[:-4]))
    (n,bins,patches) = X.hist(xData, bins=50, range=xLim, color='blue')
    
    # set up and plot right side Y axis histogram
    Y = fig.add_subplot(224)
    Y.set_position([0.6,0.2,0.2,0.4])
    Y.yaxis.set_ticks_position('right')
    Y.axes.set_ylim(yLim)
    Y.xaxis.set_ticklabels([]) #this is just a hack way of killing the tick labels
    Y.axes.set_ylabel(yLabel)
    Y.yaxis.set_label_position('right')
    (n,bins,patches) = Y.hist(yData, bins=50, range=yLim, orientation='horizontal', color='blue')
    if save:
        plt.savefig(plotFilename, dpi=300, orientation='landscape')
    if show:
        plt.show()

def histAndContourPlotter(xData, yData, chiSquareds, xLabel='', yLabel='', plotFilename='', xLim=False, yLim=False, show=True, save=False):
    """
    This function will create a plot with a scatter plot of the xData and yData in the middle, a horizontal histogram 
    of the yData on the right, and vertical histogram of the xData at the top.  This is a good way to compare the 
    parameter space where the two data sets have their values focus.
    
    NOTE: xData and yData must have the same length.
    
    :param xData:          Data to plot on the x axis
    :type xData:           list of numbers of any type, must be same length as yData
    :param yData:          Data to plot on the y axis
    :type yData:           list of numbers of any type, must be same length as xData
    :param xLabel:         label for the x axis data
    :type xLabel:          string, default is ''
    :param yLabel:         label for the y axis data
    :type yLabel:          string, default is ''
    :param plotFilename:   name of the file to save the figure to
    :type plotFilename:    string, ensuring to include directory path if not cwd
    :param xLim:           range to limit the x axis values to
    :type xLim:            tuple of two numbers of any type, (min,max)
        default False indicates to use min and max of xData
    :param yLim:           range to limit the y axis values to
    :type yLim:            tuple of two numbers of any type, (min,max)
        default False indicates to use min and max of yData
    :param show:           show the resulting figure on the screen?
    :type show:            boolean, default is False
    :param save:           save the resulting figure?
    :type save:            boolean, default is True
    """
    
    # check the data sets have matching lengths
    if not (len(xData)==len(yData)):
        print 'PROBLEM: The x and y data sets do not have matching lengths.'    
        
    # set the xLim and yLim if their values are False
    if not xLim:
        min = np.min(xData)
        max = np.max(xData)
        range = abs(max)+abs(min)
        xLim = (min-abs(range*0.05),max+abs(range*0.05))
    else:
        if not (type(xLim)==tuple):
            print 'PROBLEM: xLim is not of type tuple'
    if not yLim:
        min = np.min(yData)
        max = np.max(yData)
        range = abs(max) + abs(min)
        yLim = (min-abs(range*0.05),max+abs(range*0.05))
    else:
        if not (type(yLim)==tuple):
            print 'PROBLEM: yLim is not of type tuple'        
    
    # Check the plotFilename has .png, else add it
    if save:
        if plotFilename[-4:]!='.png':
            plotFilename = plotFilename+'.png'       
        
    # create figure to put plots in
    fig = plt.figure(1,figsize=(10,10) ,dpi=300)
    
    # set up and plot central scatter plot
    main = fig.add_subplot(221)
    main.set_position([0.2,0.2,0.4,0.4])
    main.axes.set_xlim(xLim)
    main.axes.set_ylim(yLim)
    main.axes.set_xlabel(xLabel)
    main.axes.set_ylabel(yLabel)
    
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    if xData.size>2e6:
        chunkSize = 100
    elif xData.size>1e5:
        chunkSize = 10
    else:
        chunkSize = 1
    numChunks = xData.size//chunkSize
    xChunks = xData[:chunkSize*numChunks].reshape((-1,chunkSize))
    yChunks = yData[:chunkSize*numChunks].reshape((-1,chunkSize))
    chiSquareds = np.array(chiSquareds)
    zChunks = chiSquareds[:chunkSize*numChunks].reshape((-1,chunkSize))
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    
    #Z = np.meshgrid(chiSquareds,chiSquareds)
    Z = np.meshgrid(zChunks,zChunks)
    CS=main.contour(xData, yData,Z)
    main.clabel(CS,fontsize=9,inline=1)
    #main.scatter(xData, yData, s=2, c='black', edgecolors='none')
    
    # set up and plot top X axis histogram
    X = fig.add_subplot(223)
    X.set_position([0.2,0.6,0.4,0.2])
    X.xaxis.set_ticks_position('top')
    X.axes.set_xlim(xLim)
    X.yaxis.set_ticklabels([]) #this is just a hack way of killing the tick labels
    X.axes.set_xlabel(xLabel)
    X.xaxis.set_label_position('top')
    #X.set_title(os.path.basename(plotFilename[:-4]))
    (n,bins,patches) = X.hist(xData, bins=50, range=xLim, color='blue')
    
    # set up and plot right side Y axis histogram
    Y = fig.add_subplot(224)
    Y.set_position([0.6,0.2,0.2,0.4])
    Y.yaxis.set_ticks_position('right')
    Y.axes.set_ylim(yLim)
    Y.xaxis.set_ticklabels([]) #this is just a hack way of killing the tick labels
    Y.axes.set_ylabel(yLabel)
    Y.yaxis.set_label_position('right')
    (n,bins,patches) = Y.hist(yData, bins=50, range=yLim, orientation='horizontal', color='blue')
    if save:
        plt.savefig(plotFilename, dpi=300, orientation='landscape')
    if show:
        plt.show()

def star(R, x0, y0, color='w', N=5, thin = 0.5):
    """
    Returns an N-pointed star of size R at (x0, y0) (matplotlib patch).
    
    NOTE: code base taken from star_patch.py in 'Beginning-Python-Visualization'
    """
    
    from pylab import Polygon
    
    polystar = np.zeros((2*N, 2))
    for i in range(2*N):
        angle = i*pi/N
        r = R*(1-thin*(i%2))
        polystar[i] = [r*np.cos(angle)+x0, r*np.sin(angle)+y0]
    return Polygon(polystar, fc=color, ec=color)    

def starAndErrorPolys(SAs,SAerrors,PAs,PAerrors,asConversion, transData, telescopeView=False):
    """
    Creates a rectangle patch for a given secondary star data point location.
    It will return the patch.
    
    SA in arcsec
    PA in degrees
    
    """
        
    errorBoxes = []
    m2starPolygons = []
    xmax=-1e6
    xmin = 1e6
    ymax=-1e6
    ymin = 1e6
    for i in range(0,len(SAs)):
        
        h = 2.0*SAerrors[i]*math.cos(math.radians(2.0*PAerrors[i]))*asConversion
        w = 2.0*SAerrors[i]*math.sin(math.radians(2.0*PAerrors[i]))*asConversion
        xCent = SAs[i]*math.sin(math.radians(PAs[i]))*asConversion
        yCent = SAs[i]*math.cos(math.radians(PAs[i]))*asConversion
        if False:
            print 'data point: SA = '+str(SAs[i])+', PA = '+str(PAs[i])#$$$$$$$$$$$$$$$$$$$$$$$$$$
            print 'location = ['+str(xCent)+', '+str(yCent)+']\n'#$$$$$$$$$$$$$$$$$$$$$$$$$$
        if telescopeView:
            xCent = -xCent
            yCent = -yCent
        xCorner = xCent-0.5*w
        yCorner = yCent-0.5*h
        
        # see if new min max found for x,y
        if xmax<(xCent+0.5*w):
            xmax = xCent+0.5*w
        if xmin>xCorner:
            xmin = xCorner
        if ymax<(yCent+0.5*h):
            ymax = yCent+0.5*h
        if ymin>yCorner:
            ymin = yCorner
            
        rect = patches.Rectangle((xCorner,yCorner),width=w,height=h,color='red',alpha=0.2)
        t = pylab.matplotlib.transforms.Affine2D().rotate_deg_around(xCent,yCent,-PAs[i]) +transData
        rect.set_transform(t)
        errorBoxes.append(rect)
        
        # determin x and y locations of the observed PA and SA's for companion star/planet
        # then make a star polygon for each, same as for M1 but much smaller
        m2starPolygons.append(star((asConversion/1000.0)*7.0*SAs[0], xCent, yCent, color='red', N=5, thin = 0.8))
        
    return (errorBoxes, m2starPolygons,[xmin,xmax,ymin,ymax])

def starAndErrorPolysOLD(SAs,SAerrors,PAs,PAerrors,asConversion, transData, telescopeView=False):
    """
    Creates a rectangle patch for a given secondary star data point location.
    It will return the patch.
    
    SA in arcsec
    PA in degrees
    
    """
        
    errorBoxes = []
    m2starPolygons = []
    for i in range(0,len(SAs)):
        
        h = 2.0*SAerrors[i]*math.cos(math.radians(2.0*PAerrors[i]))*asConversion
        w = 2.0*SAerrors[i]*math.sin(math.radians(2.0*PAerrors[i]))*asConversion
        xCent = SAs[i]*math.sin(math.radians(PAs[i]))*asConversion
        yCent = -SAs[i]*math.cos(math.radians(PAs[i]))*asConversion
        if telescopeView:
            xCent = -xCent
            yCent = -yCent
        xCorner = xCent-0.5*w
        yCorner = (yCent+0.5*h)
        
        rect = patches.Rectangle((xCorner,yCorner),width=w,height=h,color='red',alpha=0.2)
        t = pylab.matplotlib.transforms.Affine2D().rotate_deg_around(xCent,yCent,PAs[i]) +transData
        rect.set_transform(t)
        errorBoxes.append(rect)
        
        # determin x and y locations of the observed PA and SA's for companion star/planet
        # then make a star polygon for each, same as for M1 but much smaller
        m2starPolygons.append(star((asConversion/1000.0)*7.0*SAs[0], xCent, yCent, color='red', N=5, thin = 0.8))
        
    return (errorBoxes, m2starPolygons)
    
def makeCleanSummaryPlot(outputDataFilename=''):
    """
    This will first find the lowest chisquared for a given run,
    then use readTrimWrite to trim all data that is more than 40 higher than that value.
    This trimmed data file will then be plotted by summaryPlotter.
    """
    verbose = True
   
    ## if no datafilename provided, make one
    if outputDataFilename=='':
        topDir = "/mnt/Data1/Todai_Work/Data/data_Binary/data_Duo/"
        dataFolder = "MCMC-TAUBOO-STAR-PreviousTrialsAUGUST/MCMC-TauBoo-STAR-DIonly-ButlerAndDonati-900to2000P-160to190Omegas-2to90incs-100to380omegas-argPeriPlus180_OFF-43--14-Million-in_Total"
        outputDataFilename = os.path.join(topDir+dataFolder,'outputData-ALL.dat')
        
    ## double check datafile exists and get to work  
    if os.path.exists(outputDataFilename):  
        datadir = os.path.dirname(outputDataFilename)
        
        ## Get value of non-reduced chiSquared minimum
        bestOrbitFilename = os.path.join(datadir,'bestOrbit.txt')
        if verbose:
            print '\nFinding best Chisquared from file: '+bestOrbitFilename
        bestOrbitFile = open(bestOrbitFilename,'r')
        lines = bestOrbitFile.readlines()
        bestOrbitFile.close()
        for line in lines:
            if 'chiSquaredMin'in line:
                chiSquaredMin=float(line.split('=')[1])
        if verbose:
            print '\nBest chiSquared found to be = '+str(chiSquaredMin)
            
        ## get nu value, then calculate chiSquared cut off
        # get log
        logFilename = os.path.join(datadir,'log-chain_1.txt')
        [nu,nuRV,nuDI,printStr] = genTools.findNuFromLog(logFilename)
        maxNonReducedChiSquared = int(chiSquaredMin+30)###$$$$$40
        bestReducedChiSquared = chiSquaredMin/nu
        maxReducedChiSquared = maxNonReducedChiSquared/nu
        
        if verbose:
            print 'Best reduced chiSquared found to be = '+str(bestReducedChiSquared)
            print "maxReducedChiSquared = "+str(maxReducedChiSquared)
            print "Setting maxNonReducedChiSquared = "+str(maxNonReducedChiSquared)
        
        ## trim data
        newDataFilename = dataReadTrimWriteNEW2(outputDataFilename, maxNonReducedChiSquared, verbose=False)
        if verbose:
            print '\nTrimmed data file at: '+newDataFilename
        
        ## make new plotfile name
        newPlotFilename = os.path.join(datadir,'SummaryPlot-reducedChiSquaredMax_'+str(int(maxReducedChiSquared))+'.png')
        if verbose:
            print '\nNew plot for trimmed data will be written to: '+newPlotFilename
        
        ## pass trimmed data to the plotter
        if verbose:
            print '\nCalling summaryPlotter2 to plot new trimmed data'
        summaryPlotter(newDataFilename, newPlotFilename, weight=False, confLevels=True, normalize=True, nu=1, plot4x1=False, TcStepping=False)
        if verbose:
            print '\nTrimmed data all plotted up baby!! :-D\n'
    else:
        print "Filename provided does not exist!!  Provided path was: "+outputDataFilename
        
    if True:
        return newDataFilename
            

def summaryPlotter(outputDataFilename, plotFilename, weight=False, confLevels=True, normalize=True, nu=1, plot4x1=False, TcStepping=False):
    """
    This advanced plotting function will plot all the data in a 3x3 grid (or 3x1) on a single figure.  The data will be plotted
    in histograms that can be will be normalized to a max of 1.0 and the data can be weighted if desired.  The 
    confidence levels of the data can be calculated and the resulting histograms will have the bars inside the
    68% bin shaded as dark grey and 95% as light grey and the rest will be white.
    
    :param nu:       The value of nu to convert chi squared to reduced chi squared.
    :type nu:        float
    :param plot4x1:  Only make plots for K, argPeri, e, and To (or Tc if Tc stepping was used).
    :type plot4x1:   Python boolean
    """
    verbose = False
    ## find number of RV datasets
    if os.path.exists(outputDataFilename):  
        datadir = os.path.dirname(outputDataFilename)
        
        ## get log
        logFilename = os.path.join(datadir,'processManagerLogFile.txt')
        log = open(logFilename,'a')
        log.write('\n'+75*'#'+'\n Inside summaryPlotter2 \n'+75*'#'+'\n')
         
        s= '\nCreating summary plot for file:\n'+outputDataFilename
        s=s+ '\nInput plotfilename:\n'+plotFilename
        print s
        log.write(s+'\n')
        # record the time the chain started
        startTime = timeit.default_timer()
        
        ## find number of RV datasets
        f = open(outputDataFilename,'r')
        plotFileTitle = f.readline()[:-5]
        headings = f.readline()
        line = f.readline()
        dataLineCols = line.split()
        if (len(line)>11):
            numRVdatasets = len(dataLineCols) - 11
        else:
            line = f.readline()
            dataLineCols = line.split()
            if (len(line)>11):
                numRVdatasets = len(dataLineCols) - 11
        s= "\nNumber of RV datasets found in summaryPlotter2 was "+str(numRVdatasets)+"\n"
        if TcStepping==True:
            s=s+"\nTcStepping passed in was True, so plotting Tc instead of To"
        else:
            s=s+"\nTcStepping passed in was False, so plotting To instead of Tc"
        
        print s
        log.write(s+'\n')
        f.close()
           
        # check if the passed in value for plotFilename includes '.png'
        if '.png' not in plotFilename:
            plotFilename = plotFilename+'.png'
        else:
            plotFilename = plotFilename
        
        ## make an advanced title for plot from folder and filename
        titleTop = os.path.dirname(outputDataFilename).split('/')[-1]
        titleBtm = os.path.basename(plotFilename).split('.')[0]+'  TOTAL Posterior Distributions plot'
        plotFileTitle = titleTop+'\n'+titleBtm
        
        s='\nStarting to make Total Summary Plot'
        print s
        log.write(s+'\n')

        NumSamples = 0
        # Create empty figure to be filled up with plots
        if not plot4x1:
            fig = plt.figure(1, figsize=(45,20) ,dpi=300) 
        else:
            fig = plt.figure(1, figsize=(40,20) ,dpi=300)
            
        ## give the figure its title
        plt.suptitle(plotFileTitle,fontsize=30)
        
        # Create sub plot and fill it up for the Argument of Perigie
        if not plot4x1:
            subPlot = fig.add_subplot(243)
        else:
            subPlot = fig.add_subplot(223)
        startTime = timeit.default_timer()
        s='\nStarting to plot hist for argPeri_degsAlls'
        print s
        log.write(s+'\n')
        paramColNum = 6
        xlabel = 'Argument of Perigie [deg]'
        startTime1 = timeit.default_timer()
        (CLevels,data,chiSquareds,bestDataVal) = genTools.confLevelFinderNEWdataVersion(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True,fast=False)
        endTime1 = timeit.default_timer()
        totalTime = (endTime1-startTime1) # in seconds
        totalTimeString = genTools.timeString(totalTime)
        s='Getting data took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
        startTime1 = timeit.default_timer()
        subPlot = histConverter(chiSquareds, data, subPlot, xlabel, confLevels=CLevels, weight=weight, normed=normalize, nu=nu, bestVal=bestDataVal)
        endTime1 = timeit.default_timer()
        totalTime = (endTime1-startTime1) # in seconds
        totalTimeString = genTools.timeString(totalTime)
        s=s+'Plotting took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
            
        if (type(data)!=float)and(NumSamples==0):
            NumSamples=data.size
        #argPeriMedian = np.median(argPeri_degsAlls)
        s=s+ "done plotting argPeri_degsAlls:\n"
        # record the time the chain finished and print
        endTime = timeit.default_timer()
        totalTime = (endTime-startTime) # in seconds
        totalTimeString = genTools.timeString(totalTime)
        s=s+'That took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
        print s
        log.write(s+'\n')
        
        if not plot4x1:
            startTime = timeit.default_timer()
            s='\nStarting to plot hist for inclination_degsAlls:'
            print s
            log.write(s+'\n')
            subPlot = fig.add_subplot(241)
            paramColNum = 5
            xlabel = 'Inclination [deg]'
            startTime1 = timeit.default_timer()
            (CLevels,data,bestDataVal) = genTools.confLevelFinderNEWdataVersion(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=False, returnBestDataVal=True)
            endTime1 = timeit.default_timer()
            totalTime = (endTime1-startTime1) # in seconds
            totalTimeString = genTools.timeString(totalTime)
            s='Getting data took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
            startTime1 = timeit.default_timer()
            subPlot = histConverter(chiSquareds, data, subPlot, xlabel, confLevels=CLevels, weight=weight, normed=normalize, nu=nu, bestVal=bestDataVal)
            endTime1 = timeit.default_timer()
            totalTime = (endTime1-startTime1) # in seconds
            totalTimeString = genTools.timeString(totalTime)
            s=s+'Plotting took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
            subPlot.axes.set_ylabel('Probability',fontsize=30)
            if (type(data)!=float)and(NumSamples==0):
                NumSamples=data.size
            ##################$$$$$$$$$$$$$ This extra garbage collection might not be needed but I want it for now as a code EX. ######
            #del inclination_degsAlls
            #gc.collect()
            s= s+"done plotting inclination_degsAlls\n"
            # record the time the chain finished and print
            endTime = timeit.default_timer()
            totalTime = (endTime-startTime) # in seconds
            totalTimeString = genTools.timeString(totalTime)
            s=s+'That took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
            print s
            log.write(s+'\n')
        
        if not plot4x1:#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            startTime = timeit.default_timer()
            s='\nStarting to plot hist for longAN_degsAlls:'
            print s
            log.write(s+'\n')
            subPlot = fig.add_subplot(242)
            paramColNum = 0
            xlabel = 'Longitude of Ascending Node [deg]'
            (CLevels,data,bestDataVal) = genTools.confLevelFinderNEWdataVersion(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=False, returnBestDataVal=True)
            bestLongAN = bestDataVal
            s=''
            if (bestLongAN==0):
                s=s+'bestLongAN==0, so it was fixed and we might plot Ks in its place\n'
            else:
                subPlot = histConverter(chiSquareds, data, subPlot, xlabel, confLevels=CLevels, weight=weight, normed=normalize, nu=nu, bestVal=bestDataVal)
                if (type(data)!=float)and(NumSamples==0):
                    NumSamples=data.size
            #longANMedian = np.median(longAN_degsAlls)
            s=s+"done plotting longAN_degsAlls\n"
            # record the time the chain finished and print
            endTime = timeit.default_timer()
            totalTime = (endTime-startTime) # in seconds
            totalTimeString = genTools.timeString(totalTime)
            s=s+'That took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
            print s
            log.write(s+'\n')
        if True:#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            # Create sub plot and fill it up for the e
            if not plot4x1:
                subPlot = fig.add_subplot(244)
            else:
                subPlot = fig.add_subplot(221)
            startTime = timeit.default_timer()
            s='\nStarting to plot hist for esAlls:'
            print s
            log.write(s+'\n')
            paramColNum = 1
            xlabel = 'e'
            startTime1 = timeit.default_timer()
            (CLevels,data,bestDataVal) = genTools.confLevelFinderNEWdataVersion(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=False, returnBestDataVal=True)
            endTime1 = timeit.default_timer()
            totalTime = (endTime1-startTime1) # in seconds
            totalTimeString = genTools.timeString(totalTime)
            s='Getting data took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
            startTime1 = timeit.default_timer()
            subPlot = histConverter(chiSquareds, data, subPlot, xlabel, confLevels=CLevels, weight=weight, normed=normalize, nu=nu, bestVal=bestDataVal)
            endTime1 = timeit.default_timer()
            totalTime = (endTime1-startTime1) # in seconds
            totalTimeString = genTools.timeString(totalTime)
            s=s+'Plotting took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
            if (type(data)!=float)and(NumSamples==0):
                NumSamples=data.size
            #eMedian = np.median(esAlls)
            s=s+ "done plotting esAlls"
            # record the time the chain finished and print
            endTime = timeit.default_timer()
            totalTime = (endTime-startTime) # in seconds
            totalTimeString = genTools.timeString(totalTime)
            s=s+'That took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
            print s
            log.write(s+'\n')
        
        
        # Create sub plot and fill it up for the Time of last Periapsis OR Time of Center Transit/inferior conjunction
        if not plot4x1:
            subPlot = fig.add_subplot(245)
        else:
            subPlot = fig.add_subplot(222)
        if TcStepping==True:
            paramColNum = 3
            xlabel = 'Time of Center Transit [JD]'
        else:
            paramColNum = 2
            xlabel = 'Time of last Periapsis [JD]'
        startTime = timeit.default_timer()
        s='\nStarting to plot hist for TsAlls:'
        print s
        log.write(s+'\n')
        startTime1 = timeit.default_timer()
        (CLevels,data,bestDataVal) = genTools.confLevelFinderNEWdataVersion(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=False, returnBestDataVal=True)
        endTime1 = timeit.default_timer()
        totalTime = (endTime1-startTime1) # in seconds
        totalTimeString = genTools.timeString(totalTime)
        s='Getting data took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
        startTime1 = timeit.default_timer()
        subPlot = histConverter(chiSquareds, data, subPlot, xlabel, confLevels=CLevels, weight=weight, normed=normalize, nu=nu, bestVal=bestDataVal)
        endTime1 = timeit.default_timer()
        totalTime = (endTime1-startTime1) # in seconds
        totalTimeString = genTools.timeString(totalTime)
        s=s+'Plotting took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
        if (CLevels[1][1]-CLevels[1][0])>50:
            # shink x axis labels as there are too may big numbers and they merge otherwise
            subPlot.tick_params(axis='x',which='major',labelsize=14)
        subPlot.axes.set_ylabel('Probability',fontsize=30)
        if (type(data)!=float)and(NumSamples==0):
            NumSamples=data.size
        #TMedian = np.median(TsAlls)
        s= s+"done plotting TsAlls\n"
        # record the time the chain finished and print
        endTime = timeit.default_timer()
        totalTime = (endTime-startTime) # in seconds
        totalTimeString = genTools.timeString(totalTime)
        s=s+'That took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
        print s
        log.write(s+'\n')
            
        if True:#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            # Create sub plot and fill it up for the Ks
            if plot4x1: 
                startTime = timeit.default_timer()
                s='\nStarting to plot hist for Ks:'
                print s
                log.write(s+'\n')
                subPlot = fig.add_subplot(224)
                paramColNum = 9
                xlabel = 'K [m/s]'
                startTime1 = timeit.default_timer()
                print 'Getting data'
                
                (CLevels,data,bestDataVal) = genTools.confLevelFinderNEWdataVersion(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=False, returnBestDataVal=True)
                
                endTime1 = timeit.default_timer()
                totalTime = (endTime1-startTime1) # in seconds
                totalTimeString = genTools.timeString(totalTime)
                s='Getting data took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
                startTime1 = timeit.default_timer()
                
                subPlot = histConverter(chiSquareds, data, subPlot, xlabel, confLevels=CLevels, weight=weight, normed=normalize, nu=nu, bestVal=bestDataVal)
                
                endTime1 = timeit.default_timer()
                totalTime = (endTime1-startTime1) # in seconds
                totalTimeString = genTools.timeString(totalTime)
                s=s+'Plotting took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
                if (type(data)!=float)and(NumSamples==0):
                    NumSamples=data.size
                #periodMedian = np.median(periodsAlls)
                s= s+"done plotting Ks\n"
                # record the time the chain finished and print
                endTime = timeit.default_timer()
                totalTime = (endTime-startTime) # in seconds
                totalTimeString = genTools.timeString(totalTime)
                s=s+'\nThat took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
                print s
                log.write(s+'\n')
        if False:#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            #elif (bestLongAN==0):
            startTime = timeit.default_timer()
            s='\nStarting to plot hist for Ks:'
            print s
            log.write(s+'\n')
            subPlot = fig.add_subplot(242)
            paramColNum = 9
            xlabel = 'K [m/s]'
            (CLevels,data,bestDataVal) =genTools.confLevelFinderNEWdataVersion(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=False, returnBestDataVal=True)
            subPlot = histConverter(chiSquareds, data, subPlot, xlabel, confLevels=CLevels, weight=weight, normed=normalize, nu=nu, bestVal=bestDataVal)
            if (type(data)!=float)and(NumSamples==0):
                NumSamples=data.size
            #periodMedian = np.median(periodsAlls)
            s= "done plotting Ks\n"
            # record the time the chain finished and print
            endTime = timeit.default_timer()
            totalTime = (endTime-startTime) # in seconds
            totalTimeString = genTools.timeString(totalTime)
            s=s+'That took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
            print s
            log.write(s+'\n')
        
        if False:
            # Create sub plot and fill it up for the semi-majors
            if not plot4x1: 
                startTime = timeit.default_timer()
                s='\nStarting to plot hist for Semi-Majors:'
                print s
                log.write(s+'\n')
                subPlot = fig.add_subplot(247)
                paramColNum = 7
                xlabel = 'Semi-Major[AU]'
                (CLevels,data,bestDataVal) =genTools.confLevelFinderNEWdataVersion(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=False, returnBestDataVal=True)
                subPlot = histConverter(chiSquareds, data, subPlot, xlabel, confLevels=CLevels, weight=weight, normed=normalize, nu=nu, bestVal=bestDataVal)
                if (type(data)!=float)and(NumSamples==0):
                    NumSamples=data.size
                semiMajorCLevels = CLevels
                semiMajorBest = bestDataVal
                if NumSamples<2e7:
                    semiMajors = data
                #periodMedian = np.median(periodsAlls)
                s="done plotting Semi-Majors\n"
                # record the time the chain finished and print
                endTime = timeit.default_timer()
                totalTime = (endTime-startTime) # in seconds
                totalTimeString = genTools.timeString(totalTime)
                s=s+'That took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
                print s
                log.write(s+'\n')
        
        # Create sub plot and fill it up for the Period
        if not plot4x1: 
            startTime = timeit.default_timer()
            s='\nStarting to plot hist for periodsAlls:'
            print s
            log.write(s+'\n')
            subPlot = fig.add_subplot(246)
            paramColNum = 4
            xlabel = 'Period [Years]'
            (CLevels,data,bestDataVal) =genTools.confLevelFinderNEWdataVersion(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=False, returnBestDataVal=True)
            subPlot = histConverter(chiSquareds, data, subPlot, xlabel, confLevels=CLevels, weight=weight, normed=normalize, nu=nu, bestVal=bestDataVal)
            periodCLevels = CLevels
            periodBest = bestDataVal
            if NumSamples<2e7:
                periods = data
            #periodMedian = np.median(periodsAlls)
            s= "done plotting periodsAlls"
            # record the time the chain finished and print
            endTime = timeit.default_timer()
            totalTime = (endTime-startTime) # in seconds
            totalTimeString = genTools.timeString(totalTime)
            s=s+'That took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
            print s
            log.write(s+'\n')
            
        ## create plot for predicted total mass if possible
            if False:
                if not plot4x1:
                    startTime = timeit.default_timer()
                    s='\nStarting to plot hist for Mtotals:'
                    print s
                    log.write(s+'\n')
                    # conversion factors and constants
                    SecPerYear = 31557600.0
                    G = 6.67300e-11
                    MperAU = 149598000000.0
                    KGperMsun = 1.98892e30
                    
                    #print 'Starting to calculate the total Mass array'
                    #print 'There were '+str(NumSamples)+" samples in total"
                    consts = ((4.0*np.pi**2.0)/G)*((MperAU**3.0)/(KGperMsun*SecPerYear**2.0))
                    
                    #semiMajorCLevels = np.array(semiMajorCLevels)
                    #periodCLevels = np.array(periodCLevels)
                    # calculate the CLevels for the mass
                    #CLevels=[[]]
                    print 'semiMajorCLevels = '+repr(semiMajorCLevels)
                    print 'periodCLevels = '+repr(periodCLevels)
                    CLevelsA = consts*((semiMajorCLevels[0][0]**3.0)/(periodCLevels[0][0]**2.0))
                    CLevelsB = consts*((semiMajorCLevels[0][1]**3.0)/(periodCLevels[0][1]**2.0))
                    CLevelsC = consts*((semiMajorCLevels[1][0]**3.0)/(periodCLevels[1][0]**2.0))
                    CLevelsD = consts*((semiMajorCLevels[1][1]**3.0)/(periodCLevels[1][1]**2.0))
                    CLevels=[[CLevelsA,CLevelsB],[CLevelsC,CLevelsD]]
                    print 'CLevels = '+repr(CLevels)
                    #print 'CLevels normal = '+repr([[CLevelsA,CLevelsB],[CLevelsC,CLevelsD]])
                    bestVal = consts*((semiMajorBest**3.0)/(periodBest**2.0))
                    
                    if NumSamples<2e7:
                    
                        Mtotals = consts*np.divide(np.power(semiMajors,3.0),np.power(periods,2.0)) # in Msun
                        #print "total mass array has "+str(Mtotals.size)+" elements"
                        #print "calculating CLevels for Mtotals"
                        
                        #print "Mtotals CLevels found to be "+repr(CLevels)
                        subPlot = fig.add_subplot(248)
                        xlabel = 'Total Mass [Msun]'
                        #print "starting to plot Mtotals"
                        chiSquareds = np.array(chiSquareds)
                        #bestVal = Mtotals[np.where(chiSquareds==chiSquareds.min())][0]
                        subPlot = histConverter(chiSquareds, Mtotals, subPlot, xlabel, confLevels=CLevels, weight=weight, normed=normalize, nu=nu, bestVal=bestVal)
                        s= "done plotting Mtotals"
                        # record the time the chain finished and print
                        endTime = timeit.default_timer()
                        totalTime = (endTime-startTime) # in seconds
                        totalTimeString = genTools.timeString(totalTime)
                        s=s+'\nThat took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
                        Median = np.median(Mtotals)
                    s =s+ "\n"+"*"*25 
                    s=s+"\nBest Fit value of Mtotals = "+str(bestVal)
                    if NumSamples<2e7:
                        s=s+"\nMedian value of Mtotals = "+str(Median)
                    s=s+"\n68.3% conf levels = "+repr(CLevels[0])
                    s=s+"\n95.4% conf levels = "+repr(CLevels[1])
                    if NumSamples<2e7:
                        s=s+"\ntotal range = "+repr([np.min(Mtotals),np.max(Mtotals)])
                    s =s+"\n"+"*"*25+"\n"
                    print s
                    log.write(s+'\n')
        
        if True:
           # Create sub plot and fill it up for the Period
           if not plot4x1: 
               startTime = timeit.default_timer()
               if TcStepping==True:
                   paramColNum = 2
                   xlabel = 'Time of Periapsis [JD]'
                   s='\nStarting to plot hist for Ts:'
               else:
                   paramColNum = 3
                   xlabel = 'Time of Center Transit [JD]'
                   s='\nStarting to plot hist for Tcs:'
               print s
               log.write(s+'\n')
               subPlot = fig.add_subplot(247)
               (CLevels,data,bestDataVal) =genTools.confLevelFinderNEWdataVersion(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=False, returnBestDataVal=True)
               subPlot = histConverter(chiSquareds, data, subPlot, xlabel, confLevels=CLevels, weight=weight, normed=normalize, nu=nu, bestVal=bestDataVal)
               s= "done plotting Tcs or Ts"
               # record the time the chain finished and print
               endTime = timeit.default_timer()
               totalTime = (endTime-startTime) # in seconds
               totalTimeString = genTools.timeString(totalTime)
               s=s+'That took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
               print s
               log.write(s+'\n')
        
        
        if True:#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            if not plot4x1:
                legendStr = ''
                ## Get value of non-reduced chiSquared minimum
                bestOrbitFilename = os.path.join(datadir,'bestOrbit.txt')
                if verbose:
                    print '\nFinding best Chisquared from file: '+bestOrbitFilename
                bestOrbitFile = open(bestOrbitFilename,'r')
                lines = bestOrbitFile.readlines()
                bestOrbitFile.close()
                for line in lines:
                    legendStr=legendStr+line
                    if 'chiSquaredMin'in line:
                        chiSquaredMin=float(line.split('=')[1])
                if verbose:
                    print '\nBest chiSquared found to be = '+str(chiSquaredMin)
            
                ## get nu value, then calculate chiSquared cut off
                # get log
                logFilename = os.path.join(datadir,'log-chain_1.txt')
                [nu,nuRV,nuDI,printStr] = genTools.findNuFromLog(logFilename)
                bestReducedChiSquared = chiSquaredMin/nu
                legendStr=legendStr+'\nReducedChiSquared = '+str(bestReducedChiSquared)+'\n'
                subPlot = fig.add_subplot(248)
                subPlot.text(0.05,0.05,legendStr,ha='left',fontsize=20)
                subPlot.axes.set_yticklabels([])
                subPlot.axes.set_xticklabels([])
                
                
        # Save file if requested.
        print '\nStarting to save param hist figure:'
        if plotFilename!='':
            plt.savefig(plotFilename, dpi=300, orientation='landscape')
            s= 'Summary plot saved to: '+plotFilename
            print s
            log.write(s+'\n')
        plt.close()
        
        if (numRVdatasets>0)and(False):
            ## Create a second figure of RV offsets. ####
            try:
                 # Create empty figure to be filled up with plots
                 # Create sub plot and fill it up for the Semi-major
                 if numRVdatasets<9:
                     fig = plt.figure(2, figsize=(42,50) ,dpi=250)
                 else:
                     s= 'summaryPlotter2: WARNING!!! Plotter not set up to handle more than 9 RV datasets and '\
                     +str(numRVdatasets)+' were found in the output datafile.'
                     print s
                     log.write(s+'\n')
                 
                 #add figure title
                 plt.suptitle(plotFileTitle,fontsize=30)
                 
                 for dataset in range(1,numRVdatasets+1):
                     startTime = timeit.default_timer()
                     subPlot = fig.add_subplot(330+dataset)
                     paramColNum = 9+dataset
                     print 'Starting to plot '+'RV offset '+str(dataset)+":"
                     xlabel = 'RV offset '+str(dataset)+' [m/s]'
                     (CLevels,data,bestDataVal) =genTools.confLevelFinderNEWdataVersion(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=False, returnBestDataVal=True)
                     subPlot = histConverter(chiSquareds, data, subPlot, xlabel, confLevels=CLevels, weight=weight, normed=normalize, nu=nu, bestVal=bestDataVal)
                     if ((dataset==1)or(dataset==4))or(dataset==7):
                        subPlot.axes.set_ylabel('Probability',fontsize=50)
                     subPlot.tick_params(axis='both',which='major',labelsize=30)
                     subPlot.axes.set_xlabel(xlabel,fontsize=50)
                     s= "\nDone plotting RV offsets for dataset "+str(dataset)
                     # record the time the chain finished and print
                     endTime = timeit.default_timer()
                     totalTime = (endTime-startTime) # in seconds
                     totalTimeString = genTools.timeString(totalTime)
                     s=s+'\nThat took '+totalTimeString+' to complete.\n'  ##### only print in 'silent' mode to track time
                     print s
                     log.write(s+'\n') 
                     
                 # Save file 
                 print '\nStarting to save RVoffsets figure'
                 plotFilename2 = plotFilename[0:-4]+'-RVoffsets.png'
                 plt.savefig(plotFilename2, dpi=300, orientation='landscape')
                 s= 'RV offsets summary figure saved to '+plotFilename2
                 print s
                 log.write(s+'\n')  
                 plt.close()
            except:
                 plt.close()
                 s= 'No RV offsets to plot' 
                 print s
                 log.write(s+'\n')  
                 
        ## Make a chiSquared distribution         
        if True:
            fig = plt.figure(2, figsize=(35,15) ,dpi=250)
            subPlot = fig.add_subplot(111)
            xlabel = 'chiSquare - chiSquare_MIN'
            CLevels = [[0,0],[0,0]]
            useAry = chiSquareds[int(len(chiSquareds)/2.0):]
            useAry = np.array(useAry)
            useAry = useAry-useAry.min()
            subPlot = histConverter(chiSquareds, useAry, subPlot, xlabel, weight=False, bestVal=0.0001,logY=True)
            
            #add reduced chisquared axis labels
            x1,x2 = subPlot.get_xlim()
            ax2 = subPlot.twiny()
            ax2.set_xlim(((1.0/nu)*x1),((1.0/nu)*x2))
            ax2.figure.canvas.draw()
            ax2.set_xlabel("reduced(chiSquare - chiSquare_MIN)")
            subPlot.axes.set_ylabel('Probability')
            
            #subPlot.set_yscale('log')
            print 'Starting to save chiSquared figure:'
            plotFilename3 = plotFilename[0:-4]+'-ChiSquaredDist.png'
            plt.savefig(plotFilename3, dpi=250, orientation='landscape')
            s= 'chiSquared dist summary figure saved to '+plotFilename3
            print s
            log.write(s+'\n')  
            plt.close()
            
        # record the time the chain finished and print
        endTime = timeit.default_timer()
        totalTime = (endTime-startTime) # in seconds
        totalTimeString = genTools.timeString(totalTime)
        s= '\n\nsummaryPlotter2: Plotting took '+totalTimeString+' to complete.\n'
        print s
        log.write(s+'\n')
        log.write('\n'+75*'#'+'\n Leaving summaryPlotter2 \n'+75*'#'+'\n')
        log.close()
    else:
        s= "summaryPlotter2: ERROR!!!! file doesn't exist"
        print s
        log.write(s+'\n')
        log.write('\n'+75*'#'+'\n Leaving summaryPlotter2 \n'+75*'#'+'\n')
        log.close()
        
def summaryPlotter2MCMC(outputDataFilename, plotFilename, nu=1, plot4x1=False, logFilename='',TcStepping=False):
    """
    This will plot a progress plot for each varying parameters in the file.  
    It is called by mcmcProgressPlotter to plot up the progress for each chain
    separately and utilizes summaryPlotter2MCMCfunc to do the repetative task 
    of producing the plot for each parameter.
    
    These plots will be made for the parameter values vs step, 
    to show how it is converging with time.
    
    :param nu:       The value of nu to convert chi squared to reduced chi squared.
    :type nu:        float
    :param plot4x1:  Only make plots for K, argPeri, e, To (or Tc if Tc stepping was used), and chiSquared.
    :type plot4x1:   Python boolean
    """
    verbose = False
    #find chain number and update logFilename with path and number
    s = outputDataFilename
    chainNumStr = s[s.find('chain_')+6:s.find('chain_')+6+1]
    datadir = os.path.dirname(outputDataFilename)
    logFilename = os.path.join(datadir,'log-chain_'+chainNumStr+'.txt')
    log = open(logFilename,'a')
    log.write('\n'+75*'*'+'\n Inside summaryPlotter2MCMC \n'+75*'*'+'\n')
    
    ## find number of RV datasets
    if os.path.exists(outputDataFilename):   
        s=  '\nCreating summary plot for file:\n'+outputDataFilename
        s=s+'\nInput plotfilename:\n'+plotFilename
        print s
        log.write(s+'\n')
        # record the time the chain started
        startTime = timeit.default_timer()
        
        ## find number of RV datasets
        f = open(outputDataFilename,'r')
        plotFileTitle = f.readline()[:-5]
        headings = f.readline()
        line = f.readline()
        dataLineCols = line.split()
        if (len(line)>11):
            numRVdatasets = len(dataLineCols) - 11
        else:
            line = f.readline()
            dataLineCols = line.split()
            if (len(line)>11):
                numRVdatasets = len(dataLineCols) - 11
        s= "\nNumber of RV datasets found in summaryPlotter2MCMC was "+str(numRVdatasets)+"\n"
        print s
        log.write(s+'\n')
        f.close()
           
        # check if the passed in value for plotFilename includes '.png'
        if '.png' not in plotFilename:
            s='input plotFilename: ',plotFilename
            if logFilename!='':
                log.write(s+'\n')
            if verbose:
                print s
            plotFilename = plotFilename+'.png'
            s='plotFilename changed to: ',plotFilename
            log.write(s+'\n')
            if verbose:
                print s 
        else:
            plotFilename = plotFilename
        
        ## make an advanced title for plot from folder and filename
        titleTop = os.path.dirname(outputDataFilename).split('/')[-1]
        titleBtm = os.path.basename(plotFilename).split('.')[0]+" Parameter Progress Plot"
        plotFileTitle = titleTop+'\n'+titleBtm
        
        s='Starting to make Total Summary Plot'
        log.write(s+'\n')
        if verbose:
            print s 

        # Create empty figure to be filled up with plots
        if not plot4x1:
            fig2 = plt.figure(2, figsize=(30,55),dpi=200)
            s= '\n** a 7x1 figure will be made for all 6 orbital parameters in the output datafile **'
            print s
            log.write(s+'\n')
        else:
            fig2 = plt.figure(2,figsize=(30,35),dpi=200)
            s= '\n** a 4x1 figure will be made for all 6 orbital parameters in the output datafile **'
            print s
            log.write(s+'\n')
        
        # add figure title
        plt.suptitle(plotFileTitle,fontsize=30)
        
        ## NOTE: this one is being done in full code to get all the needed values for the later plots to use
        ##       for the function that reduces code.
        # Create sub plot and fill it up for the Argument of Perigie
        if not plot4x1:
            subPlot2 = fig2.add_subplot(816)
        else:
            subPlot2 = fig2.add_subplot(512)
        paramColNum = 6
        xlabel = 'Argument of Perigie [deg]'
        (CLevels,data,chiSquareds,bestDataVal) =genTools.confLevelFinderNEWdataVersion(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True,fast=False)
        log.write('\nCLevels for '+xlabel+':\n'+repr(CLevels)+'\n')
        xs = np.array(range(0,data.size))
        # Wrap the array into a 2D array of chunks, truncating the last chunk if 
        if data.size>100e6:
            chunkSize = 10000
        elif data.size>10e6:
            chunkSize = 1000
        elif data.size>1e6:
            chunkSize = 100
        elif data.size>1e5:
            chunkSize = 10
        else:
            chunkSize = 1
        numChunks = data.size//chunkSize
        yChunks = data[:chunkSize*numChunks].reshape((-1,chunkSize))
        xChunks = xs[:chunkSize*numChunks].reshape((-1,chunkSize))
        # Calculate the max, min, and means of chunksize-element chunks...
        max_env = yChunks.max(axis=1)
        min_env = yChunks.min(axis=1)
        yCenters = yChunks.mean(axis=1)
        xCenters = xChunks.mean(axis=1)
        s=  '\n'+50*'#'+'\nNumber of datapoints TOTAL: '+repr(data.size)
        s=s+'\nNumber of binned, chunkSize='+repr(chunkSize)+", datapoints = "+repr(xCenters.size)
        s=s+'\n'+50*'#'+'\n'
        log.write(s+'\n')
        if verbose:
            print s
        subPlot2.fill_between(xCenters, min_env, max_env, color='gray', edgecolor='none', alpha=0.5)
        subPlot2.plot(xCenters, yCenters)
        subPlot2.plot([xCenters.min(),xCenters.max()],[bestDataVal,bestDataVal],color='green')
        #subPlot2.plot(range(0,data.size),data)
        subPlot2.axes.set_ylabel(xlabel)
        #argPeriMedian = np.median(argPeri_degsAlls)
        s= "done plotting argPeri_degsAlls"
        print s
        log.write(s+'\n')
        
        if not plot4x1:
            subPlot2 = fig2.add_subplot(815)
            paramColNum = 5
            xlabel = 'Inclination [deg]'
            (log,subPlot2,data,bestDataVal)=summaryPlotter2MCMCfunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize)
            ##################$$$$$$$$$$$$$ This extra garbage collection might not be needed but I want it for now as a code EX. ######
            #del inclination_degsAlls
            #gc.collect()
            s= "done plotting inclination_degsAlls"
            print s
            log.write(s+'\n')
        
        if not plot4x1:
            subPlot2 = fig2.add_subplot(811)
            paramColNum = 0
            xlabel = 'Longitude of Ascending Node [deg]'
            (log,subPlot2,data,bestDataVal)=summaryPlotter2MCMCfunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize)
            #longANMedian = np.median(longAN_degsAlls)
            s= "done plotting longAN_degsAlls"
            print s
            log.write(s+'\n')
        
        # Create sub plot and fill it up for the e
        if not plot4x1:
            subPlot2 = fig2.add_subplot(812)
        else:
            subPlot2 = fig2.add_subplot(511)
        paramColNum = 1
        xlabel = 'e'
        (log,subPlot2,data,beste)=summaryPlotter2MCMCfunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize)
        #eMedian = np.median(esAlls)
        s= "done plotting esAlls"
        print s
        log.write(s+'\n')
           
        # Create sub plot and fill it up for the Semi-Major Amplitude
        if plot4x1: 
            subPlot2 = fig2.add_subplot(513)
            paramColNum = 9
            xlabel = 'K [m/s]'
            (log,subPlot2,data,bestDataVal)=summaryPlotter2MCMCfunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize)
            #periodMedian = np.median(periodsAlls)
            s= "done plotting Ks"
            print s
            log.write(s+'\n')
            
        # Create sub plot and fill it up for the Period
        if not plot4x1: 
            subPlot2 = fig2.add_subplot(814)
            paramColNum = 4
            xlabel = 'Period [Years]'
            (log,subPlot2,data,bestDataVal)=summaryPlotter2MCMCfunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize)
            #periodMedian = np.median(periodsAlls)
            s= "done plotting periodsAlls"
            print s
            log.write(s+'\n')
            
        # Create sub plot and fill it up for the semi-majors
        if not plot4x1: 
            subPlot2 = fig2.add_subplot(817)
            paramColNum = 7
            xlabel = 'Semi-Majors [AU]'
            (log,subPlot2,data,bestDataVal)=summaryPlotter2MCMCfunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize)
            #periodMedian = np.median(periodsAlls)
            s= "done plotting semi-majors"
            print s
            log.write(s+'\n')
            
        # Create sub plot and fill it up for the Time of last Periapsis OR center transit
        if beste>0.000005:
            # don't bother plotting To if the eccent is too low
            # as it will be same as argPeri anyway
            if not plot4x1:
                subPlot2 = fig2.add_subplot(813)
            else:
                subPlot2 = fig2.add_subplot(514)
            if TcStepping:
                paramColNum = 3
                xlabel = 'Time of Center Transit [JD]'
            else:
                paramColNum = 2
                xlabel = 'Time of Last Periapsis [JD]'
            (log,subPlot2,data,bestDataVal)=summaryPlotter2MCMCfunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize)
            #TMedian = np.median(TsAlls)
            s= "done plotting TsAlls"
            print s
            log.write(s+'\n')
                
        ## NOTE: this chiSquareds plot is special, so it can't use the function to do the work
        if not plot4x1:        
            subPlot2 = fig2.add_subplot(818)
        else:
            subPlot2 = fig2.add_subplot(515)
        paramColNum = 8
        xlabel = 'ChiSquareds'
        (CLevels,data) =genTools.confLevelFinderNEWdataVersion(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=False)
        log.write('\nCLevels for '+xlabel+':\n'+repr(CLevels)+'\n')
        yChunks = data[:chunkSize*numChunks].reshape((-1,chunkSize))
        yCenters = yChunks.mean(axis=1)
        subPlot2.fill_between(xCenters, min_env, max_env, color='gray',edgecolor='none', alpha=0.5)
        minY2 = data[int(data.size*0.69):-2].min()
        maxY2 = data[int(data.size*0.69):-2].max()
        subPlot2.axes.set_ylim([minY2,maxY2])
        subPlot2.plot(xCenters, yCenters)
        subPlot2.axes.set_ylabel(xlabel)
        #add reduced chisquared axis labels
        x1,x2 = subPlot2.get_ylim()
        ax2 = subPlot2.twinx()
        ax2.set_ylim(((1.0/nu)*x1),((1.0/nu)*x2))
        ax2.figure.canvas.draw()
        ax2.set_ylabel("Reduced chiSquareds")
        #################################################
        s= "done plotting ChiSquareds < "+str(maxY2)
        print s
        log.write(s+'\n')                   
        
        # Save file if requested.
        if plotFilename!='':
            plt.savefig(plotFilename, orientation='landscape')
            s= '\nSummary plot saved to:\n'+plotFilename
            print s
            log.write(s+'\n')
        plt.close()
        
        if False:
            ## separate figure for the chiSquareds ##
            fig3 = plt.figure(1,figsize=(30,10),dpi=200)
            subPlot2 = fig3.add_subplot(211)
            paramColNum = 8
            xlabel = 'ChiSquareds'
            (CLevels,data) =genTools.confLevelFinderNEWdataVersion(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=False)
            log.write('\nCLevels for '+xlabel+':\n'+repr(CLevels)+'\n')
            yChunks = data[:chunkSize*numChunks].reshape((-1,chunkSize))
            max_env = yChunks.max(axis=1)
            min_env = yChunks.min(axis=1)
            yCenters = yChunks.mean(axis=1)
            subPlot2.fill_between(xCenters, min_env, max_env, color='gray',edgecolor='none', alpha=0.5)
            subPlot2.plot(xCenters, yCenters)
            #subPlot2.plot(range(0,data.size),data)
            subPlot2.axes.set_ylabel(xlabel)
            s= "done plotting ChiSquareds full Ranges"
            print s
            log.write(s+'\n')
            
            subPlot2 = fig3.add_subplot(212)
            subPlot2.fill_between(xCenters, min_env, max_env, color='gray',edgecolor='none', alpha=0.5)
            minY2 = data[int(data.size*0.8):-2].min()
            maxY2 = data[int(data.size*0.8):-2].max()
            subPlot2.axes.set_ylim([minY2,maxY2])
            subPlot2.plot(xCenters, yCenters)
            subPlot2.axes.set_ylabel(xlabel)
            s= "done plotting ChiSquareds < "+str(maxY2)
            print s
            log.write(s+'\n')
            
            # Save file if requested.
            if plotFilename!='':
                s='input plotFilename: '+plotFilename
                log.write(s+'\n')
                if verbose:
                    print s
                plotFilename3 = plotFilename[:-4]+'-chiSquareds.png'
                s= 'plotFilename changed to: '+plotFilename3
                log.write(s+'\n')
                if verbose:
                    print s
                    
                plt.savefig(plotFilename3, orientation='landscape')
                s= '\nSummary plot saved to:\n'+plotFilename
                print s
                log.write(s+'\n')
            plt.close()
        
        
        ## Create a second figure of RV offsets. ####
        if numRVdatasets>0:
            try:
                 # Create empty figure to be filled up with plots
                 # Create sub plot and fill it up for the Semi-major
                 if numRVdatasets==1:
                     fig2 = plt.figure(2,figsize=(30,10),dpi=200)
                 elif numRVdatasets==2:
                     fig2 = plt.figure(2,figsize=(30,15),dpi=200)
                 elif numRVdatasets==3:
                     fig2 = plt.figure(2,figsize=(30,20),dpi=200)
                 else:
                     s= '\nsummaryPlotter2MCMC: WARNING!!! Plotter not set up to handle more than 3 RV datasets and '\
                     +str(numRVdatasets)+' were found in the output datafile.'
                     print s
                     log.write(s+'\n')
                 
                 
                 if numRVdatasets>=1:
                     # Create sub plot and fill it up for the Semi-major
                     if numRVdatasets==1:
                         subPlot2 = fig2.add_subplot(111)
                     elif numRVdatasets==2:
                         subPlot2 = fig2.add_subplot(211)
                     elif numRVdatasets==3:
                         subPlot2 = fig2.add_subplot(311)
                     paramColNum = 10
                     xlabel = 'RV offset 1 [m/s]'
                     (log,subPlot2,data,bestDataVal)=summaryPlotter2MCMCfunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize)
                     s= "Done plotting RV offsets for dataset 1"
                     print s
                     log.write(s+'\n')
                 
                 if numRVdatasets>=2:
                     # Create sub plot and fill it up for the Semi-major
                     if numRVdatasets==2:
                         subPlot2 = fig2.add_subplot(212)
                     elif numRVdatasets==3:
                         subPlot2 = fig2.add_subplot(312)
                     paramColNum = 11
                     xlabel = 'RV offset 2 [m/s]'
                     (log,subPlot2,data,bestDataVal)=summaryPlotter2MCMCfunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize)
                     s= "Done plotting RV offsets for dataset 2"
                     print s
                     log.write(s+'\n')
                     
                 if numRVdatasets==3:
                     # Create sub plot and fill it up for the Semi-major
                     subPlot2 = fig2.add_subplot(313)
                     paramColNum = 12
                     xlabel = 'RV offset 3 [m/s]'
                     (log,subPlot2,data,bestDataVal)=summaryPlotter2MCMCfunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize)
                     s= "Done plotting RV offsets for dataset 3"
                     print s
                     log.write(s+'\n')
                         
                 # Save file 
                 s = '\ninput plotFilename:\n'+plotFilename
                 log.write(s+'\n')
                 if verbose:
                    print s
                 plotFilename2 = plotFilename[:-4]+'-RVoffsets.png'
                 s='\nplotFilename changed to:\n'+plotFilename2
                 log.write(s+'\n')
                 if verbose:
                    print s
                 plt.savefig(plotFilename2, orientation='landscape')
                 s= '\nRV offsets summary figure saved to:\n'+plotFilename2
                 print s
                 log.write(s+'\n')                 
                 plt.close()
            except:
                 plt.close()
                 s= '\nNo RV offsets to plot' 
                 print s
                 log.write(s+'\n')     
        # record the time the chain finished and print
        endTime = timeit.default_timer()
        totalTime = (endTime-startTime) # in seconds
        totalTimeString = genTools.timeString(totalTime)
        s= '\n\nsummaryPlotter2MCMC: Plotting took '+totalTimeString+' to complete.\n'
        print s
        log.write(s+'\n')     
        log.write('\n'+75*'*'+'\n Leaving summaryPlotter2MCMC \n'+75*'*'+'\n')
        log.close()
    else:
        s= "\nsummaryPlotter2MCMC: ERROR!!!! file doesn't exist"
        print s
        log.write(s+'\n')
        log.write('\n'+75*'*'+'\n Leaving summaryPlotter2MCMC \n'+75*'*'+'\n')
        log.close()

def summaryPlotter2MCMCfunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize):
    """
    works for summaryPlotter2MCMC to do the plot for each param and reduce code doubling
    """
    (CLevels,data,bestDataVal) =genTools.confLevelFinderNEWdataVersion(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=False, returnBestDataVal=True,fast=True)
    log.write('\nCLevels for '+xlabel+':\n'+repr(CLevels)+'\n')
    if type(data)!=float:
        yChunks = data[:chunkSize*numChunks].reshape((-1,chunkSize))
        max_env = yChunks.max(axis=1)
        min_env = yChunks.min(axis=1)
        yCenters = yChunks.mean(axis=1)
        subPlot2.fill_between(xCenters, min_env, max_env, color='gray',edgecolor='none', alpha=0.5)
        subPlot2.plot(xCenters, yCenters)
        subPlot2.plot([xCenters.min(),xCenters.max()],[bestDataVal,bestDataVal],color='green')
    else:
        subPlot2.plot([xCenters.min(),xCenters.max()],[bestDataVal,bestDataVal],color='red')
    subPlot2.axes.set_ylabel(xlabel)
    
    return (log,subPlot2,data,bestDataVal)

def orbitEllipsePlotter(longAN_deg, e, period, inc, argPeri_deg, a, sysDataDict, DIdataDict,\
                            xLabel='E ["]', yLabel='N ["]', \
                          plotFilename='', xLim=False, yLim=False, show=True, telescopeView=False,To=0,nuDI=1 ):
    """
    This function will plot the resulting orbit for the input parameters.
    NOTE: If a plotFilename is provided, then the resulting figure will be saved.
    
    IN THIS VERSION the input orbital elements are to be lists, so that multiple orbits are drawn on one 
    plot.  The first orbit provided, that in the 0th element of the lists, will be used for the 
    main orbit which will get its 1/4 orbit sections marked by colored stars and other stuff.

    :param longAN_degs      Longitude of the Acending Node in degrees
    :type longAN_deg:       float, in a list
    :param argPeri_deg:     Argument of Periapsis in degrees
    :type argPeri_deg:      float, in a list
    :param a:               Semi-Major axis in AU
    :type a:                float, in a list
    :param e:               Eccentricity
    :type e:                float, in a list
    :param inc:             Inclination in degrees
    :type inc:              float, in a list
    :param period:          period of orbits [yrs]
    :type period:           float, in a list

    :param Sys_Dist_PC:     Distance to the system from Earth [parsec]
    :type Sys_Dist_PC:      float
    :param xLim:            range to limit the x axis values to
    :type xLim:             tuple of two numbers of any type, (min,max),
        default False indicates to use min and max of X values for points on
        the ellipse +5% for white space padding.
    :param yLim:            range to limit the y axis values to
    :type yLim:             tuple of two numbers of any type, (min,max),
        default False indicates to use min and max of Y values for points on 
        the ellipse +5% for white space padding.
    :param PAs:             Position Angles in [degrees]
    :type PAs:              list of numbers, same length as SAs
    :param SAs:             Separation Angles in ["]
    :type SAs:              list of numbers, same lenght as PAs
    :param sys_dist:        Distance to System in [pc]
    :type sys_dist:         float
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
    epochs = DIdataDict['DI_epochs']
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
    
    ## make an advanced title for plot from folder and filename
    titleTop = os.path.dirname(plotFilename).split('/')[-1]
    titleBtm = os.path.basename(plotFilename).split('.')[0]
    plotFileTitle = titleTop+'\n'+titleBtm

    # figure out logfilename and open log
    log = open(logFilename,'a')
    log.write('\n'+75*'#'+'\n Inside orbitEllipsePlotter \n'+75*'#'+'\n')
    
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
    if type(To)!=list:
        To = [To]
       
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
    yStar = c_foci*np.sin(np.radians(ang_corr))*(asConversion/sys_dist)
    xStar = c_foci*np.cos(np.radians(ang_corr))*(asConversion/sys_dist)
    
    if telescopeView:
        yStar = -yStar
        xStar = -xStar
#    print 'c_foci = '+str(c_foci)
#    print 'yStar = '+str(yStar)
#    print 'xStar = '+str(xStar) 
        
    # create a primary star polygon
    starPolygon = star((asConversion/1000.0)*2.0*a[0], 0, 0, color='black', N=10, thin = 0.4)
    
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
                diTools.orbitCalculatorSAPA(t, sys_dist, inc[orb], longAN_deg[orb], e[orb], T, period[orb], argPeri_deg[orb], a[orb],\
                                                            Mass1=1, Mass2=1, verbose=False)
            #(PA, SA) = PASAcalculator(period, t, T, e, inc, longAN_deg, argPeri_deg, sys_dist, a, a1=0, verbose=False)
            #orbitTAs.append(TA_deg)
            #orbitPAs.append(PA)
            #orbitSAs.append(SA)
            ellipseX = SA*math.sin(math.radians(PA))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
            ellipseY = SA*math.cos(math.radians(PA))*asConversion#*sys_dist   # This will be in [mas] instead of [AU]
            if telescopeView:
                ellipseX = -ellipseX
                ellipseY = -ellipseY
            ellipseXs.append(ellipseX)
            ellipseYs.append(ellipseY)
            #sep_dist = math.sqrt(math.pow(ellipseX,2.0)+math.pow(ellipseY,2.0))
            #sep_dists.append(sep_dist)
        ellipseXs2.append(ellipseXs)
        ellipseYs2.append(ellipseYs)
        
    
    ##************************************************************************************
    ##****************** Calculate the predicted location and fixed JDs ******************
    ##************************************************************************************
    ts = [2452847,
          2453166
          ]
    
    for orb in range(0,len(longAN_deg)):
        s= '@'*50+'\n'
        s=s+ "Calculating the predicted location of the companion for a fixed list of epochs:\n\n"
        s2=''
        s3=''
        print 'To[orb] = '+str(To[orb])#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        for t in ts:
            (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SA, PA, a1, a2) =\
                diTools.orbitCalculatorSAPA(t, sys_dist, inc[orb], longAN_deg[orb], e[orb], To[orb], period[orb], argPeri_deg[orb], a[orb],\
                                                            Mass1=1, Mass2=1, verbose=False)
                
            x = SA*math.sin(math.radians(PA))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
            y = SA*math.cos(math.radians(PA))*asConversion#*sys_dist   # This will be in [mas] instead of [AU]
            s=s+ 'Epoch '+str(t)+': x = '+str(x)+' , y = '+str(y)+'\n'
            s2=s2+str(x)+'    '+str(y)+'\n'
            s3=s3+str(SA)+'     '+str(PA)+'\n'
        s=s+ '@'*50+"\n"
        s=s+ '\n       *** Excel Format ***\n       x                  y      \n'
        s=s+ s2
        
        s = s+'\n       SA           PA        \n'
        s=s+s3
        
        print s
        
        log.write(s+'\n')
        
    ##************************************************************************************
    ##************************************************************************************
    ##************************************************************************************
    
    
    ## Get the locations of 500 points on an ellipse representing the orbit # this isn't working right... thus orbit method above is used now.
    #X,Y = ellipse(a_rp, b_rp, ang_corr, -xStar, -yStar, Nb=500)
    #X,Y = ellipse(a_rp, b_rp, ang, 0, 0, Nb=500)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, a1, a2) =\
            diTools.orbitCalculatorSAPA(0.0, sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    Xstart = SAstart*math.sin(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    Ystart = SAstart*math.cos(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    if telescopeView:
        Xstart = -Xstart
        Ystart = -Ystart
    startStar = star((asConversion/1000.0)*0.6*a[0], Xstart, Ystart, color='green', N=20, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, a1, a2) =\
            diTools.orbitCalculatorSAPA(((period[0]/4.0)*365.25), sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    XoneQuarter = SAstart*math.sin(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    YoneQuarter = SAstart*math.cos(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    if telescopeView:
        XoneQuarter = -XoneQuarter
        YoneQuarter = -YoneQuarter
    oneQuarterStar = star((asConversion/1000.0)*0.6*a[0], XoneQuarter, YoneQuarter, color='blue', N=20, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAhalf, PAhalf, a1, a2) =\
            diTools.orbitCalculatorSAPA(((period[0]/2.0)*365.25), sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    Xhalf = SAhalf*math.sin(math.radians(PAhalf))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    Yhalf = SAhalf*math.cos(math.radians(PAhalf))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    if telescopeView:
        Xhalf = -Xhalf
        Yhalf = -Yhalf
    halfStar = star((asConversion/1000.0)*0.6*a[0], Xhalf, Yhalf, color='blue', N=20, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, a1, a2) =\
            diTools.orbitCalculatorSAPA((3.0*(period[0]/4.0)*365.25), sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    XthreeQuarter = SAstart*math.sin(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    YthreeQuarter = SAstart*math.cos(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    if telescopeView:
        XthreeQuarter = -XthreeQuarter
        YthreeQuarter = -YthreeQuarter
    threeQuarterStar = star((asConversion/1000.0)*0.6*a[0], XthreeQuarter, YthreeQuarter, color='blue', N=20, thin = 0.2)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, a1, a2) =\
            diTools.orbitCalculatorSAPA(((period[0])*365.25), sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    Xend = SAstart*math.sin(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    Yend = SAstart*math.cos(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    if telescopeView:
        Xend = -Xend
        Yend = -Yend
    endStar = star((asConversion/1000.0)*0.6*a[0], Xend, Yend, color='purple', N=20, thin = 0.5)
      
    #######################################################################################              
    ## Get the calculated chiSquared fit to the data for these orbital parameters
    legendStr = ''
    for orb in range(0,len(longAN_deg)):
        legendStr = legendStr+"\nFor orbit # "+str(orb)+':\n'
        legendStr = legendStr+"inc[orb] = "+str(inc[orb])+'\n'
        legendStr = legendStr+"longAN_deg[orb] = "+str(longAN_deg[orb])+'\n'
        legendStr = legendStr+"e[orb] = "+str(e[orb])+'\n'
        legendStr = legendStr+"To[orb] = "+str(To[orb])+'\n'
        legendStr = legendStr+"period[orb] = "+str(period[orb])+'\n'
        legendStr = legendStr+"argPeri_deg[orb] = "+str(argPeri_deg[orb])+'\n'
        legendStr = legendStr+"a[orb] = "+str(a[orb])+'\n'
        if False:
            (chi_squared_total, ns, Ms, Es, thetas, Sep_Dists, SA_arcsec_measured_models, PA_deg_measured_models, a1s, a2s) =\
            diTools.multiEpochOrbCalc(SAs, SAerrors, PAs, PAerrors,epochs, sys_dist, inc[orb], longAN_deg[orb],\
                            e[orb], To[orb], period[orb], argPeri_deg[orb], a_total=a[orb], Mass1=1, Mass2=1, verbose=False)
            print "!!!!!!!!!!!!!!!!!!!!\n\n Try1: for orbit #"+str(orb)+", the chiSquared fit to the DI data was = "+str(chi_squared_total)+', or reduced = '+str(chi_squared_total/nuDI)+"\n!!!!!!!!!!!!!!!!!!!!!!!!!"
        (chi_squared_total, ns, Ms, Es, thetas, xs, ys, a1s, a2s) =\
        diTools.multiEpochOrbCalc3(SAs, SAerrors, PAs, PAerrors,epochs, sys_dist, inc[orb], longAN_deg[orb],\
                        e[orb], To[orb], period[orb], argPeri_deg[orb], a_total=a[orb], Mass1=1, Mass2=1, verbose=False)
        chiSquaredStr = "The chiSquared fit to the DI data was = "+str(chi_squared_total)+', or reduced = '+str(chi_squared_total/nuDI)
        legendStr = legendStr+chiSquaredStr+'\n'
        print "\n\n Try2: for orbit #"+str(orb)+", "+chiSquaredStr+"\n"
    #######################################################################################        
    
    ## Create figure, axes and start plotting/drawing everything
    fig = plt.figure(1,figsize=(10,10))
    main = fig.add_subplot(111)
    main.set_xlabel(xLabel, fontsize=30)
    main.set_ylabel(yLabel, fontsize=30)
    plt.suptitle(plotFileTitle, fontsize=10)
        
    # Draw orbits
    for orb in range(0,len(longAN_deg)):
        main.plot(ellipseXs2[orb],ellipseYs2[orb],linewidth=2) #$$$$ add title, labels, axes
    # Draw Ellipse ## found this didn't work and thus the orbit method above is used now.
    #main.plot(X,Y, c='black')
    
    #draw semi-major
    main.plot([Xstart,Xhalf],[Ystart,Yhalf],'g-',linewidth=1)
    
    #draw semi-minor 
    #main.plot(X_b_rp_corr,Y_b_rp_corr,'b-')
    
    #draw 'x' for center of ellipse
    main.plot(-xStar,yStar,'rx',linewidth=1)
    
    # draw lines along diagonals to mare out 45 degree locations from origin ### Useless now
    #main.plot( [xLim[0],xLim[1]], [yLim[0],yLim[1]])
    #main.plot( [xLim[1],xLim[0]], [yLim[0],yLim[1]])
    
    # Draw stars for the location of each 1/4 of the orbit  # These are optional, I would not include them for the final versions, maybe just the periapsis one.
    main.add_patch(startStar)
    if False:
        main.add_patch(oneQuarterStar)
        main.add_patch(halfStar)
        main.add_patch(threeQuarterStar)
        #main.add_patch(endStar)
    
    # Draw larger star for primary star's location
    main.add_patch(starPolygon)
    
    ## call function to calculate, create and return polygons for the 
    ## companion star locations and boxes for their errors
    (errorBoxes, m2starPolygons,[xmin,xmax,ymin,ymax]) = starAndErrorPolys(SAs,SAerrors,PAs,PAerrors,asConversion, main.transData, telescopeView)
    
    dataMaxMins = [xmin,xmax,ymin,ymax]
    
    #print "from starAndErrorPolys [xmin,xmax,ymin,ymax] = "+repr([xmin,xmax,ymin,ymax])
    ## set the xLim and yLim if their values are False
    ## and pad max and min values found by 10%
    if not xLim:
        min = findArrayMin(ellipseXs2[:])
        max = findArrayMax(ellipseXs2[:])
        Range = abs(max)+abs(min)
        xLim = (min-abs(Range*0.05),max+abs(Range*0.05))
        xLim = (findArrayMin([xLim[0],xmin]), findArrayMax([xLim[1],xmax])) 
        #print 'elipseXs2 min = '+str(min)+", max = "+str(max)+", final xLim = "+repr(xLim)
    else:
        if not (type(xLim)==tuple):
            s= 'PROBLEM: xLim is not of type tuple'
            print s
            log.write(s+'\n')
    if not yLim:
        min = findArrayMin(ellipseYs2[:])
        max = findArrayMax(ellipseYs2[:])
        Range = abs(max)+abs(min)
        yLim = (min-abs(Range*0.05),max+abs(Range*0.05))
        yLim = (findArrayMin([yLim[0],ymin]), findArrayMax([yLim[1],ymax])) 
        #print 'elipseYs2 min = '+str(min)+", max = "+str(max)+", final yLim = "+repr(yLim)
    else:
        if not (type(yLim)==tuple):
            s= 'PROBLEM: yLim is not of type tuple'
            print s
            log.write(s+'\n')
    
    # Draw lines for horizontal and vertical from origin
    main.plot([xLim[0],xLim[1]],[0,0],c='black',linewidth=2)
    main.plot([0,0],[yLim[0],yLim[1]],c='black',linewidth=2)
    
    # draw the error boxes for the companion start locations in the data
    for errorBox in errorBoxes:
        main.add_patch(errorBox)
        
    # Draw red star patches for each of the companion's locations from the data
    for star2 in m2starPolygons:
        main.add_patch(star2)
        
    # add a legend
    #main.legend(('longAN_deg = '+str(longAN_deg),'e = '+str(e), 'period = '+str(period), 'inc = '+str(inc), 'argPeri_deg = '+str(argPeri_deg), 'a = '+str(a)), loc=0, markerscale=0.0000000000000001)
    paramsLegndStr = 'longAN_deg = '+str(longAN_deg[0])+'\ne = '+str(e[0])+'\nperiod = '+str(period[0])+'\ninc = '+str(inc[0])+'\nargPeri_deg = '+str(argPeri_deg[0])+'\na = '+str(a[0])
    if False:
        main.text(xLim[0]+abs(xLim[0]*0.02),abs(yLim[1]*0.2),paramsLegndStr,ha='left')
    
    ## FLIP X-AXIS to match backawards Right Ascension definition
    a = main.axis()
    main.axis([a[1],a[0],a[2],a[3]])
    xLim2 = (xLim[1],xLim[0])
    #print "xLim2 = "+repr(xLim2)
    main.axes.set_xlim(xLim2)
    main.axes.set_ylim(yLim)
    main = fixPlotBordersAndLabels(main)
    main.axhline(linewidth=1.0)
    main.axvline(linewidth=1.0)
    
    
    #print "final xlims in plot after flip"+repr(main.axes.get_xlim())
    #main.set_title("HR-8799e  best fit: Face-on")#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # save plot to file
    if plotFilename!='':
        plt.savefig(plotFilename, dpi=300, orientation='landscape')
    
    ## find cropped limits for just data plus 5% pad
    dataMaxMins = [xmin,xmax,ymin,ymax]
    dataXrange = abs(xmax)+abs(xmin)
    dataYrange = abs(ymax)+abs(ymin)
    xlimData = [xmax+0.025*dataXrange,xmin-0.025*dataXrange]
    ylimData = [ymin-0.025*dataYrange,ymax+0.025*dataYrange]
    main.axes.set_xlim(xlimData)
    main.axes.set_ylim(ylimData)
    
    plotFilenameCropped = plotFilename[:-4]+"-CROPPED.png"
    ## saved cropped version to file if requested
    if plotFilename!='':
        plt.savefig(plotFilenameCropped, dpi=300, orientation='landscape')
       
    ## put limits back to full 
    main.axes.set_xlim(xLim2)
    main.axes.set_ylim(yLim)
    # show plot
    if show:
        plt.show()
    
    # close figure/plot
    plt.close()    
    
    if True:
        ## Create figure for writting the sorta legend to
        fig = plt.figure(1,figsize=(10,10))
        #main = fig.add_subplot(111)
        fig.text(0.05,0.05,legendStr,ha='left')
        if plotFilename!='':
            legendFigFilename = plotFilename[:-4]+"-paramInfo.png"
            plt.savefig(legendFigFilename, dpi=300, orientation='landscape')
        plt.close()
    
    # write last log statement and close log
    log.write('\n'+75*'#'+'\n Leaving orbitEllipsePlotter \n'+75*'#'+'\n')
    log.close()
  

def mcmcProgressPlotter(filenames,rootPlotFilename,nu=1, plot4x1=False,TcStepping=False):
    """
    This will produce a time series type view of the paramters in each chain's output
    data file.  It will utilize summaryPlotter2MCMC to produce these plots and figures.
    
    These plots will be made for the parameter values vs step, 
    to show how it is converging with time.
    
    NOTE: Runs on NEW data files only.
    
    :param nu:       The value of nu to convert chi squared to reduced chi squared.
    :type nu:        float
    :param plot4x1:  Only make plots for K, argPeri, e, To (or Tc if Tc stepping was used), and chiSquared.
    :type plot4x1:   Python boolean
    """
    # check if the passed in value for plotFilename includes '.png'
    if '.png' not in rootPlotFilename:
        rootPlotFilename = rootPlotFilename+'.png'
        
    # check filenames parameter
    if type(filenames) is not list:
        print 'the filenames input parameter must be a list of strings'
    
    numFiles = len(filenames)
    
#    ## get settings used for simulation $$$ only need to see if a 3x2 or 3x1 plot should be used for parameter plots
#    inputSettingsFile = os.path.join("/run/media/Kyle/Data1/Todai_Work/Dropbox/workspace/Binary-project/SimSettings_and_InputData",'SimSettings.txt')
#    # make an output settings file name for the duo version output
#    [fName,ext] = os.path.splitext(inputSettingsFile)
#    outputSettingsFile = fName+'_DuoVersion'+ext
#    paramSettingsDict = cFileToSimSettingsDict(outputSettingsFile)
#    
#    print '&& The value of paramSettingsDict["simulate_StarPlanet"] = '+repr(paramSettingsDict["simulate_StarPlanet"])

    for filename in filenames:
        # check if the passed in value for filename includes '.txt'
        if (filename[-4:]!='.txt' and filename[-4:]!='.dat'):
            filename = filename+'.dat'
            print ".dat was added to the end of the filename as it wasn't found to be in the original version"
        #find chain number and update plotFilename with path and number
        s = filename
        chainNumStr = s[s.find('chain_')+6:s.find('chain_')+6+1]
        plotFilename = rootPlotFilename[:-4]+"-chain_"+chainNumStr+'.png'
        datadir = os.path.dirname(filename)
        logFilename = os.path.join(datadir,'log-chain_'+chainNumStr+'.txt')
        log = open(logFilename,'a')
        if os.path.dirname(plotFilename)=='':
            plotFilename = os.path.join(os.path.dirname(filename),plotFilename)
            
        s= '\n###### Making MCMC progress summary plots for chain '+chainNumStr+' ######\n'
        #[nu,printStr] = genTools.findNuFromLog(logFilename)
        s = s#+'\n'+printStr
        print s
        log.write(s+'\n')
        log.close()
        #call plotter for this file
        summaryPlotter2MCMC(filename, plotFilename, nu=nu, plot4x1=plot4x1,TcStepping=TcStepping)

def diPlotterTester(outputDatafile=''):
    """
    A function to allow custom inputs to orbitEllipsePlotter for testing or 
    producing one off plots for a particular set of orbital Elements.
    """
    
    sysDatafilename = os.path.join("/mnt/Data1/Todai_Work/Dropbox/workspace/Binary-project/SimSettings_and_InputData",'SystemData.txt')
    sysDataDict = sysDataToDict(sysDatafilename)
    
    DIdatafilename = os.path.join("/mnt/Data1/Todai_Work/Dropbox/workspace/Binary-project/SimSettings_and_InputData",'DIdata.dat')
    if os.path.exists(DIdatafilename):
        DIdataDict = diTools.DIdataToDict(DIdatafilename)
    if outputDatafile=='':
        outputDatafile = "/mnt/Data1/Todai_Work/Data/data_Binary/data_Duo/MCMC-parmVaryCorrected-HR8799-UniEccentPrior-circular-veryOpenParms-1--28-Million-in_Total/outputData-ALL.dat"
    
    print '#'*50
    bestOrbit = bestOrbitFinderNEW(outputDatafile, printToScreen=True, saveToFile=False, returnAsList=True)
    print '#'*50
#    #Tau Boo planet test params
#    longAN_deg = 148.620748
#    e = 0.187655
#    period = 98.8500326216
#    inc = 44.44006
#    argPeri_deg = 104.51581
#    a = 24.82015
    
    if True:
        
        if True:
            plotFilename=os.path.join(os.path.dirname(outputDatafile),"UpdatedPlot")
            orbitEllipsePlotter(bestOrbit[0],bestOrbit[1],bestOrbit[3],bestOrbit[4],bestOrbit[5],bestOrbit[6], sysDataDict, DIdataDict,\
                            xLabel='E ["]', yLabel='N ["]', \
                          plotFilename=plotFilename, xLim=False, yLim=False, show=True, telescopeView=False, To=bestOrbit[2] )
        if False:
            plotFilename=os.path.join(os.path.dirname(outputDatafile),"All4OrbitsPlot")
            longAns = [bestOrbit[0],28.4,169.96,174.0]
            es = [bestOrbit[1],0.91,0.41,0.76]
            Ps = [bestOrbit[3],2000,389.25,996.0]
            incs = [bestOrbit[4],51.68,70.6,49]
            argPeris = [bestOrbit[5],99.34,0.67,322]
            As = [bestOrbit[6],245,98.41,120]
            Ts = [bestOrbit[2],1727754,2331859,2100788]
            orbitEllipsePlotter(longAns,es,Ps,incs,argPeris,As, sysDataDict, DIdataDict,\
                            xLabel='E ["]', yLabel='N ["]', \
                          plotFilename=plotFilename, xLim=False, yLim=False, show=True, telescopeView=False, To=Ts )
    if False:
        print 'Starting to make a Dual hist, contour plot'
        
        #xlabel = 'Argument of Perigie [deg]'
        paramColNum = 5
        #(CLevels,data,chiSquareds,bestDataVal) =genTools.confLevelFinderNEWdataVersion(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True)
        #print 'Argument of Perigie data loaded'
        
        yLabel = 'Inclination [deg]'
        paramColNum = 4
        (CLevels,yData,chiSquareds,bestDataVal) =genTools.confLevelFinderNEWdataVersion(outputDatafile,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True)
        print 'Inclination data loaded' 
                    
        #xlabel = 'Longitude of Ascending Node [deg]'
        paramColNum = 0
        #(CLevels,data,chiSquareds,bestDataVal) =genTools.confLevelFinderNEWdataVersion(outputDatafile,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True)
        #print 'Longitude of Ascending Node data loaded'
        
        #xlabel = 'e'
        paramColNum = 1
        #(CLevels,data,chiSquareds,bestDataVal) =genTools.confLevelFinderNEWdataVersion(outputDatafile,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True)
        #print 'e data loaded'
        
        xLabel = 'Period [Years]'
        paramColNum = 3
        (CLevels,xData,chiSquareds,bestDataVal) =genTools.confLevelFinderNEWdataVersion(outputDatafile,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True)
        print 'Period data loaded'
        
        #xlabel = 'Time of last Periapsis [JD]'
        paramColNum = 2
        #(CLevels,data,chiSquareds,bestDataVal) =genTools.confLevelFinderNEWdataVersion(outputDatafile,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True)
        #print 'Time of last Periapsis data loaded'
        
        #xlabel = 'K [m/s]'
        paramColNum = 8
        #(CLevels,data,chiSquareds,bestDataVal) =genTools.confLevelFinderNEWdataVersion(outputDatafile,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True)
        #print 'K data loaded'
        
        histAndContourPlotter(xData, yData, chiSquareds, xLabel='', yLabel='', plotFilename='', xLim=False, yLim=False, show=True, save=False)
    
def PostSimCompleteAnalysisFunc(outputDatafile=''):
    """
    A function that performs many of the post simulation completion plotting 
    or Gelman Rubin statistic calculations as done in MCMC_ProcessManager.  This 
    can be used for custom plotting or testing, or to perform any plotting that 
    was turned off in the initial run of a simulation that you might now want.
    """
    if outputDatafile=='':
        sysDatafilename = os.path.join("/mnt/Data1/Todai_Work/Dropbox/workspace/Binary-project/SimSettings_and_InputData",'SystemData.txt')
        RVdatafilename = os.path.join("/mnt/Data1/Todai_Work/Dropbox/workspace/Binary-project/SimSettings_and_InputData",'RVdata.dat')
        inputSettingsFile = os.path.join("/mnt/Data1/Todai_Work/Dropbox/workspace/Binary-project/SimSettings_and_InputData",'SimSettings.txt')
    else:
        sysDatafilename = os.path.join(os.path.dirname(outputDatafile),'code-used/SystemData.txt')
        RVdatafilename = os.path.join(os.path.dirname(outputDatafile),'code-used/RVdata.dat')
        inputSettingsFile = os.path.join(os.path.dirname(outputDatafile),'code-used/SimSettings.txt')
    sysDataDict = sysDataToDict(sysDatafilename)
    RVdataDict = rvTools.RVdataToDict(RVdatafilename)
    paramSettingsDict = cFileToSimSettingsDict(inputSettingsFile)
    
    numModDataSets = 10
    rvPlotFilename = os.path.join(os.path.dirname(outputDatafile),"RVplot-Manual")#"DotaniAndButlerPre1995PlanetOrbitFit"
    diPlotFilename = os.path.join(os.path.dirname(outputDatafile),"DIplot-Manual")
    modDatasetsFilename =  os.path.join("/mnt/Data1/Todai_Work/Data/data_Binary/data_Duo","mod"+str(numModDataSets)+"Dataset.dat")
    TvsEccPlotFilename = os.path.join("/mnt/Data1/Todai_Work/Data/data_Binary/data_Duo","TvsEccentricityPlot")
    
    ## get nu value, then calculate chiSquared cut off
    # get log
    logFilename = os.path.join(os.path.dirname(outputDatafile),'log-chain_1.txt')
    [nu,nuRV,nuDI,printStr] = genTools.findNuFromLog(logFilename)
    
    if outputDatafile=='':
        outputDatafile = "/mnt/Data1/Todai_Work/Data/data_Binary/data_Duo/DotaniAndButlerPre1995-planetOrbit-results/looped5000MCMC-TauBoo-RVonly-ButlerANDdonati-noPre1995data-TrendRemoved-planet-chiSquaredUpdate-withoutFlaggedPoints-nonLogEccePriors-3--700-Thousand-in_Total/outputData-ALL.dat"
        summaryPlotFile = os.path.join("/run/media/Kyle/Data1/Todai_Work/Data/data_Binary/data_Duo","SummaryPlot-Manual")
    else:
        os.path.abspath(outputDatafile)
        summaryPlotFile = os.path.join(os.path.dirname(outputDatafile),"SummaryPlot-Manual")

    if False:
        numChains = 7
        dataFiles = []
        dataFinalFilename = os.path.join(os.path.dirname(outputDatafile),'outputData-ALL.dat')
        for num in range(1,numChains+1):
            fname = os.path.join(os.path.dirname(outputDatafile),'outputData-chain_'+str(num)+'.dat')
            dataFiles.append(fname)
            print 'Adding filename to list: '+fname
        dataFileCombiner(dataFiles, dataFinalFilename)
        
    if True:      
        print '#'*50
        bestOrbit = bestOrbitFinderNEW(outputDatafile, printToScreen=True, saveToFile=True, returnAsList=True)
        print '#'*50
        longAns = bestOrbit[0]
        e = bestOrbit[1]
        period = bestOrbit[4]
        inc = bestOrbit[5]
        argPeri_deg = bestOrbit[6]
        a = bestOrbit[7]
        T = bestOrbit[2]
        Tc = bestOrbit[3]
        K = bestOrbit[8]
        if paramSettingsDict["DIonly"]==False:
            print 'DIonly = '+repr(paramSettingsDict["DIonly"])
            RVoffsets=bestOrbit[9]
            
        ## Now manually choose what you want to plot
        if False:
            ## make RV fit plot
            rvPlotterDuo(e,T,Tc,period,inc,argPeri_deg,a, \
                      sysDataDict,RVdataDict,paramSettingsDict,K=K,RVoffsets=RVoffsets,\
                      nuRV=nuRV,plotFilename=rvPlotFilename+'-FullOrbit', show=False, plotFullOrbit=True)
            if False:
                rvPlotterDuo(e,T,Tc,period,inc,argPeri_deg,a, \
                          sysDataDict,RVdataDict,paramSettingsDict,K=K,RVoffsets=RVoffsets,\
                          plotFilename=rvPlotFilename, show=False, plotFullOrbit=False)
        if False:
            ## Make bootstrap type data for RV data
            rvModDatasetMaker(e, T, period, inc, argPeri_deg, a, sysDataDict, RVdataDict, paramSettingsDict,\
                     RVoffsets=RVoffsets, modDatasetsFilename=modDatasetsFilename, numModDatasets=numModDataSets)
        if False:
            ## Make a DI fit plot
            DIdatafilename = os.path.join(os.path.dirname(outputDatafile),'code-used/DIdata.dat')
            DIdataDict = diTools.DIdataToDict(DIdatafilename)
            if True:
                longANs = [bestOrbit[0],28.4,169.96,174.0]
                argPeris = [bestOrbit[6],99.34,0.67,322.0]
                incs = [bestOrbit[5],51.68,70.6,49.0]
                periods = [bestOrbit[4],2000.0,389.25,996.0]
                a_totals = [bestOrbit[7],245.0,98.41,120.0]
                Ts = [bestOrbit[2],1727754.0,2331859.0,2100788.0]
                es = [bestOrbit[1],0.91,0.4189,0.76]
            else:
                longANs = [bestOrbit[0]]
                argPeris = [bestOrbit[6]]
                incs = [bestOrbit[5]]
                periods = [bestOrbit[4]]
                a_totals = [bestOrbit[7]]
                Ts = [bestOrbit[2]]
                es = [bestOrbit[1]]
           
            orbitEllipsePlotter(longANs,es,periods,incs,argPeris,a_totals,\
                             sysDataDict,DIdataDict,plotFilename=diPlotFilename,show=False,To=Ts, nuDI=nuDI)  
    if False:
        ## make MCMC progress plots
        numChains = 7
        fileList = []
        for num in range(1,numChains+1):
            fileList.append(os.path.join(os.path.dirname(outputDatafile),'outputData-chain_'+str(num)+'.txt'))
        summaryPlotter2MCMC(outputDatafile, summaryPlotFile+"-MCMCprogress", weight=False, confLevels=True, nu=1, SimPlanetStar=paramSettingsDict["simulate_StarPlanet"])
    if True:
        ## Make the posterior prob histograms.
        summaryPlotter2(outputDatafile, summaryPlotFile, weight=False, confLevels=True, nu=1, plot4x1=False)
    if False:
        makeCleanSummaryPlot(outputDatafile)
    
    if False:
        ## make a scatter hist figure ###$$$ not sure if this works right now.
        outDict = outputDatafileToDict(outputDatafile)
        #outDict["Ts"]
        #outDict["es"]
        #esCLevels = ConfLevelFunc(chiSquareds,outDict["es"])
        chiSquareds = outDict["chiSquareds"]
        histAndScatterPlotter(outDict["Ts"], outDict["es"], xLabel='Ts', yLabel='es', plotFilename=TvsEccPlotFilename, xLim=False, yLim=False, show=True, save=False)
    if False:
        ### Perform second round of Gelman-Rubin
        numChains = 7
        fileList = []
        for num in range(1,numChains+1):
            fname = os.path.join(os.path.dirname(outputDatafile),'gelmanRubin-chain_'+str(num)+'.txt')
            fileList.append(fname)
            print 'Adding filename to list: '+fname
        gelmanRubinStage2(fileList)

def rvPlotter(e, T, Tc, period, inc, argPeri_deg, a, sysDataDict, RVdataDict, paramSettingsDict,\
                 K=0, RVoffsets=[0], nuRV=1, plotFilename='', show=True, plotFullOrbit=True):
    """
    create a plot for the RV data and a fit line from the best orbit data params.
    
    NOTES:
    If a star-planet system is being simulated, then the provided inclination will be used with the
    planet's mass in the system data file/ provided dict, to calculate the residual vel.  If the value
    in the dictionary/file is actually, please set inc=0 to tell this func to ignore it and use the file's value.
    """
    verbose = False
    
    ## make string with input values for possible printing
    s= '\nInputs to rvPlotterDuo were:'+'\n'
    s=s+ 'e = '+repr(e)+'\n'
    s=s+ 'To = '+repr(T)+'\n'
    s=s+ 'Tc = '+repr(Tc)+'\n'
    s=s+ 'period = '+repr(period)+'\n'
    s=s+ 'inc = '+repr(inc)+'\n'
    s=s+ 'argPeri_deg = '+repr(argPeri_deg)+'\n'
    s=s+ 'a = '+repr(a)+'\n'
    s=s+ 'K = '+repr(K)+'\n'
    s=s+ 'RVoffsets = '+repr(RVoffsets)+'\n'
    if verbose:
        print s
        
    #check the orbit element inputs to ensure they are lists, else make them lists
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
    if type(T)!=list:
       T = [T]
    if type(Tc)!=list:
       Tc = [Tc]
    if type(K)!=list:
       K = [K]
    
    if type(RVoffsets)!=list:
       RVoffsets = [RVoffsets]
    
    RV_epochs = RVdataDict['RV_epochs']
    RV_errors = RVdataDict['RV_errors']
    RVs = RVdataDict['RVs']
    
    if type(RV_epochs[0])!=list:
        RV_epochs = [RV_epochs]
    if type(RV_errors[0])!=list:
        RV_errors = [RV_errors]
    if type(RVs[0])!=list:
        RVs = [RVs]
        
    Mass1 = sysDataDict['Mass1'] 
    simulate_StarStar = paramSettingsDict["simulate_StarStar"]
    simulate_StarPlanet = paramSettingsDict["simulate_StarPlanet"]
    if (simulate_StarStar is True) and(simulate_StarPlanet is True):
        print "Error: simulate_StarStar and simulate_StarPlanet can NOT BOTH be True!"
#    if simulate_StarPlanet:            
#        T_center = sysDataDict['planet_Tc']
#        if T_center==0:
#            (To,Tc) = eccArgPeri2ToTcCalc(e[0], period[0], argPeri_deg[0], T[0], Tc=0)
#            T_center = Tc
#    else:
#        T_center = T[0]

    if plotFilename!='':
        datadir = os.path.dirname(plotFilename)
        logFilename = os.path.join(datadir,'processManagerLogFile.txt')
        if plotFilename[-4:]!='.png':
            plotFilename = plotFilename+'.png'    
    else:
        datadir = os.path.curdir
        logFilename = os.path.join(datadir,'TEMPlogFile.txt')   
    
    ## make an advanced title for plot from folder and filename
    titleTop = os.path.dirname(plotFilename).split('/')[-1]
    titleBtm = os.path.basename(plotFilename).split('.')[0]+" Best Fit Plot"
    plotFileTitle = titleTop+'\n'+titleBtm
    
    # figure out logfilename and open log
    log = open(logFilename,'a')
    log.write('\n'+75*'#'+'\n Inside rvPlotterDuo \n'+75*'#'+'\n'+s+'\n') 
    
    # make copy of input RVs (not sure if I still need to do this...)
    RV_epochsIN2 = RV_epochs
    s = '\n There were '+str(len(e))+' orbits provided to plot for '+str(len(RV_epochsIN2))+" RV datasets."
    #print 'there were also '+str(len(RV_errors))+" rv error sets"
    print s
    log.write(s+'\n')
    
    ## now convert the new epochs to phases.
    phases3 = []
    for orb in range(0,len(e)):
        s= '\nConverting new epochs to phases for orbit '+str(orb)
        if verbose:
            print s
        log.write(s+'\n')
        phases2 = epochsToPhases(RV_epochsIN2,Tc[orb],period[orb], verbose=False, halfOrbit=True) 
        #phases2 = epochsToPhases(RV_epochsIN2,T[orb],period[orb], verbose=False, halfOrbit=True)         
        phases3.append(phases2)   
    
    RVsIN = RVs
    RVsOUT = []
    numEpochsTotal = 0
    s=''
    for dataset in range(0,len(RVsIN)):
        RVsUSE=[]
        if len(RVoffsets)>=(dataset+1):
            offset = RVoffsets[dataset]
        else:
            offset = 0
        s=s+"offset being subtracted from RV data for data set #"+str(dataset)+" is "+str(offset)+"\n"
        numEpochsTotal+=len(RVsIN[dataset])
        for epoch in range(0,len(RVsIN[dataset])):
            vel = RVsIN[dataset][epoch]-offset
            RVsUSE.append(vel)
        RVsOUT.append(RVsUSE)
    if False:
        s=  '\nlen of input RVs = '+str(len(RVs))
        s=s+'\nlen of output RVs = '+str(len(RVsOUT))
        print s
    log.write(s+'\n')
        
    # data for planet's orbit that are needed to calculate proper chiSquare for RV fit
    # all below values for the planet orbit and RV data are from Butler2006.
    planet_K = sysDataDict["planet_K"] #[m/s]
    if planet_K==0:
        planet_K=False
    planet_P = sysDataDict["planet_P"]   #[days]
    #planet_p_years = planet_p/365.25
    planet_e = sysDataDict["planet_e"]    
    planet_argPeri = sysDataDict["planet_argPeri"]  #[deg]
    planet_T = sysDataDict["planet_T"]   #[JD]
    if verbose:
        print 'planet_T from dict = '+str(planet_T)
    planet_Tc = sysDataDict['planet_Tc']
    if planet_Tc==0:
        (To,Tcent) = eccArgPeri2ToTcCalc(planet_e, planet_P, planet_argPeri, planet_T, Tc=0)
        planet_Tc = Tcent
    planet_MsinI = sysDataDict["planet_MsinI"]   #[Msun]
    if planet_MsinI==0:
        planet_MsinI=False
    
    star_e = sysDataDict["star_e"]
    star_T = sysDataDict["star_T"]
    star_P = sysDataDict["star_P"]
    star_argPeri = sysDataDict["star_argPeri"]
    star_Tc = sysDataDict["star_Tc"]
    if star_Tc==0:
        (To,Tcent) = eccArgPeri2ToTcCalc(star_e, star_P, star_argPeri, star_T, Tc=0)
        star_Tc = Tcent
    star_inc = sysDataDict["star_inc"]
    star_Mass2 = sysDataDict["star_Mass2"]
    
    planet_Ks = []
    planet_Ps = []
    planet_es = []
    planet_argPeris = []
    planet_Ts = []
    planet_Tcs = []
    planet_MsinIs = []
    planet_incs = []
    
    star_es = []
    star_Ts = []
    star_Tcs = []
    star_Ps = []
    star_argPeris = []
    star_incs = []
    K_use = False
    # e, T, period, inc, argPeri_deg, a,
    if simulate_StarPlanet:
        s= 'simulate_StarPlanet==True, so loading up inputs as planet variables'
        log.write(s+'\n')
        if verbose:
            print s
        planet_Ps = period
        planet_es = e
        planet_argPeris = argPeri_deg
        planet_Ts = T
        planet_Tcs = Tc
        Mass2 = planet_MsinI
        planet_incs = inc
        # load up adjusted mass of planet taking inc into account if needed
        for orb in range(0,len(e)):
#            planet_Tc = sysDataDict['planet_Tc']
#            if planet_Tc==0:
#                (To,Tc) = eccArgPeri2ToTcCalc(planet_e, planet_P, planet_argPeri, planet_T, Tc=0)
#                planet_Tc = Tc
#            planet_Tcs.append(planet_Tc)
            if inc[orb]==0:
                planet_MsinIs.append(planet_MsinI)
            else:
                planet_MsinIs.append(planet_MsinI*math.sin(math.radians(inc[orb])))
            if K[orb]==0:
                planet_Ks.append(planet_K)
                K_use = planet_K
            else:
                K_use = K[0]
                planet_Ks.append(K[orb])
            
    else:
        for orb in range(0,len(e)):
            planet_Ps.append(planet_P)
            planet_es.append(planet_e)
            planet_argPeris.append(planet_argPeri)
            planet_Ts.append(planet_T)
            if verbose:
                print 'appending '+str(planet_T)+" to planet_Ts"
            planet_Tcs.append(planet_Tc)
            planet_MsinIs.append(planet_MsinI)
            planet_Ks.append(planet_K)
            planet_incs.append(0)
            
    if simulate_StarStar:
        s= 'simulate_StarStar==True, so loading up inputs as star variables'
        log.write(s+'\n')
        if verbose:
            print s
        star_es = e
        star_Ts = T
        star_Tcs = Tc
        star_Ps = period
        star_argPeris = argPeri_deg
        star_incs = inc
        Mass2 = star_Mass2
    else:
        for orb in range(0,len(e)):
            star_es.append(star_e)
            star_Ts.append(star_T)
            star_Tcs.append(star_Tc)
            star_Ps.append(star_P)
            star_argPeris.append(star_argPeri)
            star_incs.append(star_inc)
            
    residuals3 = []
    planetVRs3 = []
    starVRs3 = []
    for orb in range(0,len(e)):
        residuals2 = []
        planetVRs2 = []
        starVRs2 = []
        chiSquaredTot2 = 0.0
        numEpochs_RV = 0.0
        chi_squared_RV_reducedCur2 = []
        for dataset in range(0,len(RVsOUT)):
            chiSquaredTot = 0.0
            epochs = RV_epochsIN2[dataset]
            rvs = RVsOUT[dataset]
            errors = RV_errors[dataset]                
            residuals = []
            planetVRs = []
            starVRs = []
            (a_total_s, a1_s, a2_s, p_s) = semiMajorConverter(Mass1, star_Mass2, a_total=a[orb],a1=0.0,a2=0.0, period=period[orb],verbose=False)
            (a_total_p, a1_p, a2_p, p_p) = semiMajorConverter(Mass1, planet_MsinIs[orb], a_total=a[orb],a1=0.0,a2=0.0, period=period[orb],verbose=False)
            
            s= '\nFor dataset '+str(dataset)
            s=s+ "\n there are "+str(len(epochs))+" RV epochs in this set"
            s=s+ "\n there are "+str(len(rvs))+" RV vals in this set"
            s=s+ "\n there are "+str(len(errors))+" RV errors in this set\n"
            s=s+ "\nParams:"
            s=s+"\nPLANET:"
            s=s+ "\nplanet_es[orb] = "+str(planet_es[orb])
            s=s+ "\nplanet_Ts[orb] = "+str(planet_Ts[orb])
            s=s+ "\nplanet_Tcs[orb] = "+str(planet_Tcs[orb])
            s=s+ "\nplanet_Ps[orb] = "+str(planet_Ps[orb])
            s=s+ "\nplanet_argPeris[orb] = "+str(planet_argPeris[orb])
            s=s+ "\na1_p = "+str(a1_p)
            s=s+ "\nplanet_incs[orb] = "+str(planet_incs[orb])
            s=s+ "\nplanet_Ks[orb] = "+str(planet_Ks[orb])
            s=s+"\nSTAR:"
            s=s+ "\nstar_es[orb] = "+str(star_es[orb])
            s=s+ "\nstar_Ts[orb] = "+str(star_Ts[orb])
            s=s+ "\nstar_Tcs[orb] = "+str(star_Tcs[orb])
            s=s+ "\nstar_Ps[orb] = "+str(star_Ps[orb])
            s=s+ "\nstar_argPeris[orb] = "+str(star_argPeris[orb])
            s=s+ "\na1_s = "+str(a1_s)
            s=s+ "\nstar_incs[orb] = "+str(star_incs[orb])
            #s=s+ "\nstar_Ks[orb] = "+str(star_Ks[orb])
            log.write(s+'\n')
            if verbose:
                print s
                
            for epoch in range(0,len(epochs)):
                error = errors[epoch]
                if False:
                    s= '\nWorking on epoch ',epoch
                    print s
                    log.write(s+'\n')
                # calculate the velocity residual due to the companion star
                if (star_Ps[orb]==0):
                    (v_r_c,K_s)=(0,0)
                else:
                    (v_r_c,K_s) = rvTools.vrCalculatorStar2(epochs[epoch],star_es[orb],T[orb],star_Ps[orb],star_argPeris[orb],a1_s,T_center=star_Tcs[orb],i=star_incs[orb], K=False, verbose=False)
                # calculate the velocity residual due to the planet around primary
                if (planet_Ps[orb]==0):
                    (v_r_p,K_p)=(0,0)
                else:
                    s=  "\nParams:"
                    s=s+"\nepochs[epoch] = "+str(epochs[epoch])
                    log.write(s+'\n')
                    if False:
                        print s
                        
                    (v_r_p,K_p) = rvTools.vrCalculatorStar2(epochs[epoch],planet_es[orb],planet_Ts[orb],planet_Ps[orb], planet_argPeris[orb],a1_p,T_center=planet_Tcs[orb],i=planet_incs[orb], K=planet_Ks[orb], verbose=False)
                    #print 'K_p being used is = ',planet_Ks[orb]
                    #print "v_r_p = ",v_r_p
                RV =rvs[epoch]- (v_r_c+v_r_p)
                if ((abs(RV)>80) and False):
                    s= 'Bad point found at '+str(rvs[epoch])+", taken at "+str(RV_epochsIN2[dataset][epoch])+", was off by = "+str(RV)
                    print s
                    log.write(s+'\n')
                residuals.append(RV)
                planetVRs.append(v_r_p) 
                starVRs.append(v_r_c) 
                chiSquaredCurr = chiSquaredCalc(rvs[epoch], error, v_r_c+v_r_p)
                chiSquaredTot = chiSquaredTot+chiSquaredCurr
                
                s= '\nRV = '+str(rvs[epoch])+' - ('+str(v_r_c)+' + '+str(v_r_p)+') = '+str(RV)
                #print 'epoch = ',epochs[epoch]
                s=s+ '\nchiSquaredCurr = '+str(chiSquaredCurr)
                s=s+ ', chiSquaredTot = '+str(chiSquaredTot)
                log.write(s+'\n')
                if False:
                    print s
            nuRV_cur = nuRV*(float(len(rvs))/float(numEpochsTotal))
            s=s+ "\nFor dataset"+str(dataset)+": nuRV = "+str(nuRV)+", len(rvs) = "+str(len(rvs))+", numEpochsTotal = "+str(numEpochsTotal)+", nuRV_cur = "+str(nuRV_cur)
            chi_squared_RV_reducedCur = (1.0/nuRV_cur)*chiSquaredTot
            chi_squared_RV_reducedCur2.append(chi_squared_RV_reducedCur)
            chiSquaredTot2 = chiSquaredTot2 + chiSquaredTot
            s='chiSquaredTot2 = '+str(chiSquaredTot2)+', with nuRV_cur = '+str(nuRV_cur)
            log.write(s+'\n')
            if False:
                print s
            residuals2.append(residuals)
            planetVRs2.append(planetVRs) 
            starVRs2.append(starVRs)
        residuals3.append(residuals2)
        planetVRs3.append(planetVRs2) 
        starVRs3.append(starVRs2)
        
    ## Make an updated Radial Velocities for data accounting for planet or star RVs
    RVsOUTupdated3 = []
    # if simulating star's orbit and planet RVs exist sub them
    if (simulate_StarStar and (planetVRs3[0][0][0]!=0)):
        s= '\n Subtracting planet VRs from data!\n'
        if verbose:
            print s
        log.write(s+'\n')
        for orb in range(0,len(e)):
            RVsOUTupdated2 = []
            for dataset in range(0,len(RVsOUT)):
                RVsOUTupdated = []
                for epoch in range(0,len(RVsOUT[dataset])): 
                    RVsOUTupdated.append(RVsOUT[dataset][epoch]-planetVRs3[orb][dataset][epoch]) 
                RVsOUTupdated2.append(RVsOUTupdated)
            RVsOUTupdated3.append(RVsOUTupdated2)
    # if simulating planet's orbit and star RVs exist sub them
    elif (simulate_StarPlanet  and (starVRs3[0][0][0]!=0)):
        s= '\n Subtracting star VRs from data!\n'
        if verbose:
            print s
        log.write(s+'\n')
        for orb in range(0,len(e)):
            RVsOUTupdated2 = []
            for dataset in range(0,len(RVsOUT)):
                RVsOUTupdated = []
                for epoch in range(0,len(RVsOUT[dataset])): 
                    RVsOUTupdated.append(RVsOUT[dataset][epoch]-starVRs3[orb][dataset][epoch])
                RVsOUTupdated2.append(RVsOUTupdated)
            RVsOUTupdated3.append(RVsOUTupdated2) 
    # if simulating either and data of the other doesn't exist, just copy RVs to 3d ary
    else:
        for orb in range(0,len(e)):
            RVsOUTupdated3.append(RVsOUT)
    
    if planetVRs3[0][0][0]!=0:
        s = '\n'+'*'*50+'\nPlanet Radial Velocities:\n'
        for orb in range(0,len(e)):
            for dataset in range(0,len(RVsOUT)):
                s=s+'\nFor dataset # '+str(dataset)+':'
                for epoch in range(0,len(RVsOUT[dataset])): 
                    if planetVRs3[orb][dataset][epoch]!=0:
                        s = s+str(planetVRs3[orb][dataset][epoch])+'\n'
        s =s+ '\n'+'*'*50
        log.write(s+'\n')
        if verbose:
            print s
        
    if starVRs3[0][0][0]!=0:
        s = '\n'+'*'*50+'\nStar Radial Velocities:\n'
        for orb in range(0,len(e)):
            for dataset in range(0,len(RVsOUT)):
                s=s+'\nFor dataset # '+str(dataset)+':'
                for epoch in range(0,len(RVsOUT[dataset])): 
                    if starVRs3[orb][dataset][epoch]!=0:
                        s = s+str(starVRs3[orb][dataset][epoch])+'\n'
        log.write(s+'\n')
        if verbose:
            print s
        
        
    # convert to a reduced chiSquared
    chi_squared_RV_reduced = (1.0/nuRV)*chiSquaredTot2
    s='\nOriginal ChiSquared = '+str(chiSquaredTot2)
    s=s+ '\nNum RV epochs = '+str(numEpochsTotal)
    s=s+ '\none over nu = '+str(1.0/nuRV)
    s=s+ '\nreduced chiSqured RV = '+str(chi_squared_RV_reduced)
    log.write(s+'\n')
    if verbose:
        print s
        
    #Get orbitRVs for best fit plot
    orbitVRs2 = []
    orbitPhases2 = []
    s= '\nabout to calc fit line'
    if verbose:
        print s
    log.write(s+'\n')
    for orb in range(0,len(e)):
        orbitVRs = []
        numSteps = 5000.0
        periodIncrement = (period[orb]*365.242)/numSteps
        t = Tc[orb]-((period[orb]*365.242)/2.0)
        times = []
        (a_total, a1, a2, p) = semiMajorConverter(Mass1, Mass2, a_total=a[orb],a1=0.0,a2=0.0, period=period[orb])
        for step in range(0,int(numSteps)):
            t = t + periodIncrement
            times.append(t)
            # calculate the velocity residual due to the companion 
            (v_r_c,K) = rvTools.vrCalculatorStar2(t,e[orb],T[orb],period[orb],argPeri_deg[orb],a1,T_center=Tc[orb],i=inc[orb], K=K_use, verbose=False)
            orbitVRs.append(v_r_c)
        #print 'times were '+repr(times)
        s= 'Orbit '+str(orb)+" had a K = "+str(K)
        if verbose:
            print s
        log.write(s+'\n')
        orbitPhases = epochsToPhases(times,Tc[orb],period[orb], verbose=False)
        #print 'orbital phases were:\n'+repr(orbitPhases)
        orbitPhases2.append(orbitPhases)
        orbitVRs2.append(orbitVRs)
        #print 'phases: ',repr(orbitPhases)
        #print 'VRs: ',repr(orbitVRs)
        
    ## Create figure, axes and start plotting/drawing everything
    fig = plt.figure(1,figsize=(10,10))
    plt.suptitle(plotFileTitle, fontsize=10)
    residualsPlot = fig.add_subplot(212)
    residualsPlot.set_title("Residuals Plot")
    residualsPlot.axes.set_xlabel("Orbital Phase",fontsize=30)
    residualsPlot.axes.set_ylabel("Residual [m/s]",fontsize=30)
    #colorsList = ['b','m','k','g','y','o','p']
    colorsList =['Blue','BlueViolet','Chartreuse','Fuchsia','Crimson','Aqua','Gold','DarkCyan','OrangeRed','Plum','DarkGreen','Chocolate','SteelBlue ','Teal','Salmon','Brown']
    
    for orb in range(0,len(argPeri_deg)):
        meanMedianChiSquaredSTR = ''
        for dataset in range(0,len(RVs)):
            #meanMedianChiSquaredSTR = meanMedianChiSquaredSTR+'\nmean of residuals for dataset '+str(dataset)+' is = '+str(np.mean(residuals3[orb][dataset]))
            #meanMedianChiSquaredSTR = meanMedianChiSquaredSTR+'\nmedian of residuals for dataset '+str(dataset)+' is = '+str(np.median(residuals3[orb][dataset]))
            meanMedianChiSquaredSTR = meanMedianChiSquaredSTR+'\nchiSquared_reduced for dataset '+str(dataset)+' is = '+str(chi_squared_RV_reducedCur2[dataset])
            residualsPlot.scatter(phases3[orb][dataset],residuals3[orb][dataset],s=35,edgecolor=colorsList[dataset],facecolor=colorsList[dataset])
            
            s= '\n\n*** RMS of residuals '+str(dataset)+'= '+str(np.sqrt(np.mean(np.array(residuals3[orb][dataset])**2)))+' ***\n\n'
            log.write(s+'\n')
            print s
            #plot error bars
            for epoch in range(0,len(phases3[orb][dataset])):
                xs = [phases3[orb][dataset][epoch],phases3[orb][dataset][epoch]]
                ys = [residuals3[orb][dataset][epoch]-abs(RV_errors[dataset][epoch]),residuals3[orb][dataset][epoch]+abs(RV_errors[dataset][epoch])]
                residualsPlot.plot(xs,ys,c='k')
               
    meanMedianChiSquaredSTR = meanMedianChiSquaredSTR+'\nchiSquared_reduced for ALL data is = '+str(chi_squared_RV_reduced)
    s= meanMedianChiSquaredSTR
    if verbose:
        print s
    log.write(s+'\n')
    
    xmin = findArrayMin(phases3)
    xmax = findArrayMax(phases3)
    ymin = findArrayMin(residuals3)
    ymax = findArrayMax(residuals3)

    xLim =[xmin,xmax]
    xrange = xmax-xmin
    xLim2 = (xmin-abs(xrange*0.05), xmax+abs(xrange*0.05))
    yLim = [ymin,ymax]
    yrange = ymax-ymin
    yLim2 = (ymin-abs(yrange*0.05), ymax+abs(yrange*0.05))
    
    ## make plot of fit to data
    fitPlot = fig.add_subplot(211)
    fitXmin = findArrayMin(orbitPhases2)
    if xmin<fitXmin:
        fitXmin = xmin
    fitXmax = findArrayMax(orbitPhases2)
    if xmax>fitXmax:
        fitXmax = xmax
    fitYmin = findArrayMin(orbitVRs2)
    if ymin<fitYmin:
        fitYmin = ymin
    fitYmax = findArrayMax(orbitVRs2)
    if ymax>fitYmax:
        fitYmax = ymax
        
    yMaxRVs = findArrayMax(RVsOUTupdated3)
    yMinRVs = findArrayMin(RVsOUTupdated3)
    yRVsRange = yMaxRVs-yMinRVs
    yLimRVs = [yMinRVs-abs(yRVsRange*0.05),yMaxRVs+abs(yRVsRange*0.05)]
        
    fitXrange = fitXmax-fitXmin
    fitXLim2 = (fitXmin-abs(fitXrange*0.05), fitXmax+abs(fitXrange*0.05))
    fitYrange = fitYmax-fitYmin
    fitYLim2 = (fitYmin-abs(fitYrange*0.05), fitYmax+abs(fitYrange*0.05))
    
    fitYLim3 = [0,0]
    if fitYLim2[0]>yLimRVs[0]:
        fitYLim3[0] = yLimRVs[0]
    else:
        fitYLim3[0] = fitYLim2[0]
    if fitYLim2[1]<yLimRVs[1]:
        fitYLim3[1] = yLimRVs[1]
    else:
        fitYLim3[1] = fitYLim2[1]
    
    if plotFullOrbit:
        fitXLimsUSE = fitXLim2
        fitYLimsUSE = fitYLim3
    else:
        fitXLimsUSE = xLim2
        fitYLimsUSE = yLimRVs
        
    residualsPlot.plot(fitXLimsUSE,[0,0],c='r',linewidth=2.0)
    residualsPlot.axes.set_xlim(fitXLimsUSE)
    residualsPlot.axes.set_ylim(yLim2)
    fitPlot.axes.set_xlim(fitXLimsUSE)
    fitPlot.axes.set_ylim(fitYLimsUSE)
    fitPlot=fixPlotBordersAndLabels(fitPlot)
    fitPlot.axes.set_ylabel("RV [m/s]",fontsize=30)
    residualsPlot=fixPlotBordersAndLabels(residualsPlot)
    for orb in range(0,len(e)):
        fitPlot.plot(orbitPhases2[orb],orbitVRs2[orb],c='r',linewidth=2.0)
        for dataset in range(0,len(RVs)):
            fitPlot.scatter(phases3[orb][dataset],RVsOUTupdated3[orb][dataset],s=35,edgecolor=colorsList[dataset],facecolor=colorsList[dataset])
            for epoch in range(0,len(phases3[orb][dataset])):
                xs = [phases3[orb][dataset][epoch],phases3[orb][dataset][epoch]]
                ys = [RVsOUTupdated3[orb][dataset][epoch]-abs(RV_errors[dataset][epoch]),RVsOUTupdated3[orb][dataset][epoch]+abs(RV_errors[dataset][epoch])]
                fitPlot.plot(xs,ys,c='k',linewidth=2.0)
            
    paramsLegndStr = 'e = '+str(e[0])+\
                '\ninc = '+str(inc[0])+\
                '\na = '+str(a[0])+'\nargPeri_deg = '+str(argPeri_deg[0])+\
                '\nperiod [days] = '+str(period[0]*365.242)+'\nTo = '+str(T[0])+", Tc = "+str(Tc[orb])
    for dataset in range(0,len(RVs)):
        if dataset==0:
            paramsLegndStr=paramsLegndStr+'\n'
        else:
            paramsLegndStr=paramsLegndStr+', '
        paramsLegndStr = paramsLegndStr+'offset '+str(dataset)+'= '+str(RVoffsets[dataset])
    paramsLegndStr = paramsLegndStr+meanMedianChiSquaredSTR
    if verbose:
        print paramsLegndStr
    log.write(paramsLegndStr+'\n')
#    addLegend=True
#    if plotFullOrbit and addLegend:
#        # in bottom right
#        #residualsPlot.text(fitXLimsUSE[1]-abs(fitXLimsUSE[1]*0.002),abs(fitYLimsUSE[1]*0.5),paramsLegndStr,ha='left')
#        # in bottom left
#        paramsLegndY = fitYLimsUSE[0]+abs(fitYLimsUSE[0]*0.02)
#        paramsLegndX = fitXLimsUSE[0]+abs(fitXLimsUSE[0]*0.03)
#        fitPlot.text(paramsLegndX,paramsLegndY,paramsLegndStr,ha='left')
#        s= '\nlegend bottom left corner at [ '+str(paramsLegndX)+' , '+str(paramsLegndY)+' ]'
#        if verbose:
#            print s
#        log.write(s+'\n')
        
    # save plot to file
    if plotFilename!='':
        plt.savefig(plotFilename, dpi=300, orientation='landscape')
        s= '\nFigure saved to:\n'+plotFilename
        if verbose:
            print s
        log.write(s+'\n')
    else: 
        s= '\nWARNING: NO plotFilename provided, so NOT saving it to disk.'
        if verbose:
            print s
        log.write(s+'\n')
    # show plot
    if show:
        plt.show()
        
    plt.close()
    
    ## Another plot to show the residuals over JD time to look for any remaining trend
    fig2 = plt.figure(2,figsize=(10,10))
    residualsPlotTrend = fig2.add_subplot(212)
    RVsPlotTrend = fig2.add_subplot(211)
    RVsPlotTrend.set_title("RVs TREND")
    residualsPlotTrend.set_title("Residuals Plot TREND")
    residualsPlotTrend.axes.set_xlabel("Epoch [JD]",fontsize=30)
    residualsPlotTrend=fixPlotBordersAndLabels(residualsPlotTrend)
    RVsPlotTrend=fixPlotBordersAndLabels(RVsPlotTrend)
    rng = findArrayMax(RV_epochsIN2)-findArrayMin(RV_epochsIN2)
    rngUse= [findArrayMin(RV_epochsIN2)-rng*0.05,findArrayMax(RV_epochsIN2)+rng*0.05]
    for orb in range(0,len(argPeri_deg)):
        for dataset in range(0,len(RVs)):
                residualsPlotTrend.scatter(RV_epochsIN2[dataset],residuals3[orb][dataset],s=35,edgecolor=colorsList[dataset],facecolor=colorsList[dataset])
                RVsPlotTrend.scatter(RV_epochsIN2[dataset],RVsOUT[dataset],s=35,edgecolor=colorsList[dataset],facecolor=colorsList[dataset])
                #plot error bars
                for epoch in range(0,len(RV_epochsIN2[dataset])):
                    xs = [RV_epochsIN2[dataset][epoch],RV_epochsIN2[dataset][epoch]]
                    ys = [residuals3[orb][dataset][epoch]-abs(RV_errors[dataset][epoch]),residuals3[orb][dataset][epoch]+abs(RV_errors[dataset][epoch])]
                    residualsPlotTrend.plot(xs,ys,c='k')
                    xs = [RV_epochsIN2[dataset][epoch],RV_epochsIN2[dataset][epoch]]
                    ys = [RVsOUT[dataset][epoch]-abs(RV_errors[dataset][epoch]),RVsOUT[dataset][epoch]+abs(RV_errors[dataset][epoch])]
                    RVsPlotTrend.plot(xs,ys,c='k')
    residualsPlotTrend.plot(rngUse,[0,0],c='r',linewidth=2.0)
    residualsPlotTrend.axes.set_ylim(yLim2)
    RVsPlotTrend.plot(rngUse,[0,0],c='r',linewidth=2.0)
    RVsPlotTrend.axes.set_ylim([findArrayMin(RVsOUT)-findArrayMax(RV_errors),findArrayMax(RVsOUT)+findArrayMax(RV_errors)])
    
    residualsPlotTrend.axes.set_xlim(rngUse)
    RVsPlotTrend.axes.set_xlim(rngUse)
    
    # save plot to file
    if plotFilename!='':
        plotFilename2 = plotFilename[:-4]+"_TREND.png"
        plt.savefig(plotFilename2, dpi=300, orientation='landscape')
        s= '\nFigure saved to:\n'+plotFilename2
        if verbose:
            print s
        log.write(s+'\n')
    else: 
        s= '\nWARNING: NO plotFilename provided, so NOT saving it to disk.'
        if verbose:
            print s
        log.write(s+'\n')
    # show plot
    if show:
        plt.show()
        
    plt.close()
    
    
    ### New residuals gaussian plot
    residuals3_trimmed = []
    for orb in range(0,len(argPeri_deg)):
        residuals2_trimmed = []
        for dataset in range(0,len(RVs)):
            residual_trimmed = []
            for epoch in range(0,len(residuals3[orb][dataset])):
                if 2449767.0263<=RV_epochsIN2[dataset][epoch]:
                    residual_trimmed.append(residuals3[orb][dataset][epoch]/RV_errors[dataset][epoch])
        residuals2_trimmed.append(residual_trimmed)
    residuals3_trimmed.append(residuals2_trimmed)
    
    fig3 = plt.figure(3,figsize=(10,10))
    residualsGaussian = fig3.add_subplot(111)
    residualsGaussian.hist(residuals3_trimmed,normed=1)
    residualsGaussian=fixPlotBordersAndLabels(residualsGaussian)
    x = np.linspace(-5,5,100)
    y = np.exp(-x**2/2)/np.sqrt(2*pi)
    residualsGaussian.plot(x,y, linewidth=2)
    # save plot to file
    if plotFilename!='':
        plotFilename3 = plotFilename[:-4]+"_DataDist.png"
        plt.savefig(plotFilename3, dpi=300, orientation='landscape')
        s= '\nFigure saved to:\n'+plotFilename3
        if verbose:
            print s
        log.write(s+'\n')
    # show plot and then close it
    if show:
        plt.show()
    plt.close()
    

    if True:
        ## Create figure for writting the sorta legend to
        fig = plt.figure(1,figsize=(10,10))
        #main = fig.add_subplot(111)
        fig.text(0.05,0.05,paramsLegndStr,ha='left')
        if plotFilename!='':
            legendFigFilename = plotFilename[:-4]+"-paramInfo.png"
            plt.savefig(legendFigFilename, dpi=300, orientation='landscape')
        plt.close()
    log.write('\n'+75*'#'+'\n Leaving rvPlotterDuo \n'+75*'#'+'\n') 
    log.close()

def rvModDatasetMaker(e, T_lastPeri, period, inc, argPeri_deg, a, sysDataDict, RVdataDict, paramSettingsDict,\
                 RVoffsets=[0], modDatasetsFilename='', numModDatasets=100):
    """
    First this will create a plot for the RV data and a fit line from the best orbit data params.  
    Then it will create many modified versions of the RV data set by keeping the epoch and 
    radial velocities the same, but randomizing the location of their associated errors.  
    numModDatasets modified data sets will be made and output to a file.  This 
    type of modified data can be used for 'bootstrap' type simulations although a proper MCMC with the 
    original data is the best approach normally.  Finally, a plot of the modified data sets will be 
    made in comparison to the best fit.
    
    NOTES:
    If a star-planet system is being simulated, then the provided inclination will be used with the
    planet's mass in the system data file/ provided dict, to calculate the residual vel.  If the value
    in the dictionary/file is actually, please set inc=0 to tell this func to ignore it and use the file's value.
    """
    
    T = T_lastPeri
    
    #check the orbit element inputs to ensure they are lists, else make them lists
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
    if type(T)!=list:
       T = [T]
    if type(RVoffsets)!=list:
       RVoffsets = [RVoffsets]
    
    RV_epochs = RVdataDict['RV_epochs']
    RV_errors = RVdataDict['RV_errors']
    RVs = RVdataDict['RVs']
    
    if type(RV_epochs[0])!=list:
        RV_epochs = [RV_epochs]
    if type(RV_errors[0])!=list:
        RV_errors = [RV_errors]
    if type(RVs[0])!=list:
        RVs = [RVs]
    
    Mass1 = sysDataDict['Mass1'] 
    simulate_StarStar = paramSettingsDict["simulate_StarStar"]
    simulate_StarPlanet = paramSettingsDict["simulate_StarPlanet"]
    if (simulate_StarStar is True) and(simulate_StarPlanet is True):
        print "Error: simulate_StarStar and simulate_StarPlanet can NOT BOTH be True!"
    if simulate_StarPlanet:
        T_center = sysDataDict['planet_Tc']
    else:
        T_center = T[0]
        
    if modDatasetsFilename!='':
        if ((modDatasetsFilename[-4:]!='.dat')and(modDatasetsFilename[-4:]!='.txt')):
            modDatasetsFilename = modDatasetsFilename+'.dat'     
    
    plotFilename = modDatasetsFilename[:-4]+"-Plot.png"
    if plotFilename!='':
        if plotFilename[-4:]!='.png':
            plotFilename = plotFilename+'.png'     
    ## make an advanced title for plot from folder and filename
    titleTop = os.path.dirname(outputDataFilename).split('/')[-1]
    titleBtm = os.path.basename(plotFilename).split('.')[0]+" Best Fit Plot"
    plotFileTitle = titleTop+'\n'+titleBtm
    
    # make copy of input RVs (not sure if I still need to do this...)
    RV_epochsIN2 = RV_epochs  
    
        ## now convert the new epochs to phases.
    phases3 = []
    for orb in range(0,len(e)):
        print '\nConverting new epochs to phases for orbit '+str(orb)
        phases2 = epochsToPhases(RV_epochsIN2,T_center,period[orb], verbose=False, halfOrbit=True) 
        #phases2 = epochsToPhases(RV_epochsIN2,T[orb],period[orb], verbose=False, halfOrbit=True)         
        phases3.append(phases2)
    
    RVsIN = RVs
    RVsOUT = []
    for dataset in range(0,len(RVsIN)):
        RVsUSE=[]
        if len(RVoffsets)>=(dataset+1):
            offset = RVoffsets[dataset]
        else:
            offset = 0
        for epoch in range(0,len(RVsIN[dataset])):
            vel = RVsIN[dataset][epoch]-offset
            RVsUSE.append(vel)
        RVsOUT.append(RVsUSE)
    if False:
        print 'len of input RVs = ',len(RVs)
        print 'len of output RVs = ',len(RVsOUT)
        
    # data for planet's orbit that are needed to calculate proper chiSquare for RV fit
    # all below values for the planet orbit and RV data are from Butler2006.
    planet_K = sysDataDict["planet_K"] #[m/s]
    if planet_K==0:
        planet_K=False
    planet_P = sysDataDict["planet_P"]   #[days]
    #planet_p_years = planet_p/365.25
    planet_e = sysDataDict["planet_e"]    
    planet_argPeri = sysDataDict["planet_argPeri"]  #[deg]
    planet_T = sysDataDict["planet_T"]   #[JD]
    planet_Tc = sysDataDict['planet_Tc']
    if planet_Tc==0:
        planet_Tc = T[0]
    planet_MsinI = sysDataDict["planet_MsinI"]   #[Msun]
    if planet_MsinI==0:
        planet_MsinI=False
    
    star_e = sysDataDict["star_e"]
    star_T = sysDataDict["star_T"]
    star_P = sysDataDict["star_P"]
    star_argPeri = sysDataDict["star_argPeri"]
    star_inc = sysDataDict["star_inc"]
    star_Mass2 = sysDataDict["star_Mass2"]
    
    planet_Ks = []
    planet_Ps = []
    planet_es = []
    planet_argPeris = []
    planet_Ts = []
    planet_Tcs = []
    planet_MsinIs = []
    planet_incs = []
    
    star_es = []
    star_Ts = []
    star_Ps = []
    star_argPeris = []
    star_incs = []
    K_use = False
    # e, T, period, inc, argPeri_deg, a,
    if simulate_StarPlanet:
        planet_Ps = period
        planet_es = e
        planet_argPeris = argPeri_deg
        planet_Ts = T
        Mass2 = planet_MsinI
        planet_incs = inc
        # load up adjusted mass of planet taking inc into account if needed
        for orb in range(0,len(e)):
            planet_Tcs.append(planet_Tc)
            if inc[orb]==0:
                planet_MsinIs.append(planet_MsinI)
            else:
                planet_MsinIs.append(planet_MsinI*math.sin(math.radians(inc[orb])))
            planet_Ks.append(planet_K)
            K_use = planet_K
    else:
        for orb in range(0,len(e)):
            planet_Ps.append(planet_P)
            planet_es.append(planet_e)
            planet_argPeris.append(planet_argPeri)
            planet_Ts.append(planet_T)
            planet_Tcs.append(planet_Tc)
            planet_MsinIs.append(planet_MsinI)
            planet_Ks.append(planet_K)
            planet_incs.append(0)
            
    if simulate_StarStar:
        star_es = e
        star_Ts = T
        star_Ps = period
        star_argPeris = argPeri_deg
        star_incs = inc
        Mass2 = star_Mass2
    else:
        for orb in range(0,len(e)):
            star_es.append(star_e)
            star_Ts.append(star_T)
            star_Ps.append(star_P)
            star_argPeris.append(star_argPeri)
            star_incs.append(star_inc)
            
    residuals3 = []
    RVsCalcd3 = []
    for orb in range(0,len(e)):
        residuals2 = []
        RVsCalcd2 = []
        chiSquaredTot2 = 0
        numEpochs_RV = 0
        chi_squared_RV_reducedCur2 = []
        for dataset in range(0,len(RVsOUT)):
            chiSquaredTot = 0
            epochs = RV_epochsIN2[dataset]
            rvs = RVsOUT[dataset]
            errors = RV_errors[dataset]                
            residuals = []
            RVsCalcd = []
            (a_total_s, a1_s, a2_s, p_s) = semiMajorConverter(Mass1, star_Mass2, a_total=a[orb],a1=0.0,a2=0.0, period=period[orb],verbose=False)
            (a_total_p, a1_p, a2_p, p_p) = semiMajorConverter(Mass1, planet_MsinIs[orb], a_total=a[orb],a1=0.0,a2=0.0, period=period[orb],verbose=False)
            
            if False:
                print '\nFor dataset '+str(dataset)
            for epoch in range(0,len(epochs)):
                if False:
                    print '\nWorking on epoch ',epoch
                # calculate the velocity residual due to the companion star
                if (star_Ps[orb]==0):
                    (v_r_c,K_s)=(0,0)
                else:
                    (v_r_c,K_s) = rvTools.vrCalculatorStar2(epochs[epoch],star_es[orb],T[orb],star_Ps[orb],star_argPeris[orb],a1_s,T_center=T_center,i=star_incs[orb], K=False, verbose=False)
                # calculate the velocity residual due to the planet around primary
                if (planet_Ps[orb]==0):
                    (v_r_p,K_p)=(0,0)
                else:
                    (v_r_p,K_p) = rvTools.vrCalculatorStar2(epochs[epoch],planet_es[orb],planet_Ts[orb],planet_Ps[orb], planet_argPeris[orb],a1_p,T_center=planet_Tcs[orb],i=planet_incs[orb], K=planet_Ks[orb], verbose=False)
                    #print 'K_p being used is = ',planet_Ks[orb]
                    #print "v_r_p = ",v_r_p
                RV =rvs[epoch]- (v_r_c+v_r_p)
                if (abs(RV)>100):
                    print 'Bad point found at '+str(rvs[epoch])+", "+str(RV_epochsIN2[dataset][epoch])+", was off by = "+str(RV)
                residuals.append(RV)
                if False:
                    print 'Appending RV '+str(RV)+' to the residuals array'
                RVsCalcd.append(v_r_c+v_r_p)
                chiSquaredCurr = chiSquaredCalc(rvs[epoch], errors[epoch], v_r_c+v_r_p)
                chiSquaredTot = chiSquaredTot+chiSquaredCurr
                if False:
                    print 'RV = '+str(rvs[epoch])+'- ('+str(v_r_c)+' + '+str(v_r_p)+') = '+str(RV)
                    print 'chiSquaredCurr = ',chiSquaredCurr
                    print 'chiSquaredTot = ',chiSquaredTot
            chi_squared_RV_reducedCur = (1.0/((1.0*len(rvs))-6.0))*chiSquaredTot
            chi_squared_RV_reducedCur2.append(chi_squared_RV_reducedCur)
            chiSquaredTot2 = chiSquaredTot2 + chiSquaredTot
            if False: 
                print 'chiSquaredTot2 = ',chiSquaredTot2
            if orb==0:
                numEpochs_RV_curr = len(epochs)
                numEpochs_RV = numEpochs_RV+numEpochs_RV_curr
            residuals2.append(residuals)
            RVsCalcd2.append(RVsCalcd)
        residuals3.append(residuals2)
        RVsCalcd3.append(RVsCalcd2)
        
    # convert to a reduced chiSquared
    chi_squared_RV_reduced = (1.0/((1.0*numEpochs_RV)-6.0))*chiSquaredTot2
    if True:
        print '\nOriginal ChiSquared = ',chiSquaredTot2
        print 'Num RV epochs = ',numEpochs_RV
        print 'one over nu = ',(1.0/((1.0*numEpochs_RV)-6.0))
        print 'reduced chiSqured RV = ',chi_squared_RV_reduced
        
    ###### SECTION TO PRODUCE 100 DATASETS FROM THIS BEST FIT
    ###### 3D NOW MEANS [100SETNUM,DATASET,EPOCH]
    RVsCalcdOUT3 = []
    RVerrorsOUT3 = []
    ResidualsModOUT3 = []
    RV_epochsIN3 = []
    #RV_errors3 = []          
    for modSetNum in range(0,numModDatasets):
        RVsCalcdOUT2 = []
        RVerrorsOUT2 = []
        ResidualsModOUT2 = []
        passed = False
        #print 'Making mod data set # ',(modSetNum+1)
        for dataset in range(0,len(RVsCalcd3[0])):
            RVsCalcdOUT = []
            if len(RVoffsets)>=(dataset+1):
                offset = RVoffsets[dataset]
            else:
                offset = 0
            (residuals_random, errors_random) = listRandomizer(residuals3[0][dataset],RV_errors[dataset])
            epochs = RV_epochsIN2[dataset]
            RVscalcd = RVsCalcd3[0][dataset]
            if ((len(residuals_random)==len(epochs))and(len(epochs)==len(RVscalcd))):
                passed = True
                for epoch in range(0,len(epochs)):
                    RVsCalcdOUT.append(residuals_random[epoch] + RVscalcd[epoch] )#+ offset)
                RVsCalcdOUT2.append(RVsCalcdOUT)
                ResidualsModOUT2.append(residuals_random)
                RVerrorsOUT2.append(errors_random)
        if passed:
            #print 'Dataset was appended during modSet making loop'
            #RV_errors3.append(RV_errors)
            RV_epochsIN3.append(RV_epochsIN2)
            RVsCalcdOUT3.append(RVsCalcdOUT2)
            ResidualsModOUT3.append(ResidualsModOUT2)
            RVerrorsOUT3.append(RVerrorsOUT2)
                
    print 'Size of RVsCalcdOUT3 after loop = ',len(RVsCalcdOUT3)
    ## write mod data to file
    filename = modDatasetsFilename
    f = open('%s' % filename, 'w')
    
    f.write('#'+os.path.basename(filename)+'\n')
    #f.write('#JD          RV [m/s]    delta_RV [m/s]\n')
    ## NOTE: by this point the Jitter has all ready been
    ##       added to the RV errors in quadrature
    for modSetNum in range(0,numModDatasets):
        f.write("#Modified Dataset Number "+str(modSetNum+1)+"\n")
        f.write('#JD          RV [m/s]    delta_RV [m/s]     Jitter\n')
        for dataset in range(0,len(RVsOUT)):
            epochs = RV_epochsIN3[modSetNum][dataset]
            for epoch in range(0,len(epochs)):
                line = str(RV_epochsIN3[modSetNum][dataset][epoch])
                line = line +'  '+ str(RVsCalcdOUT3[modSetNum][dataset][epoch])
                line = line +'  '+ str(RVerrorsOUT3[modSetNum][dataset][epoch])
                line = line +'  '+str(0.0)+'\n'
                f.write(line)
        f.write("\n")
            
    f.close()
    print 'Output 100 mod dataset file written to:'+modDatasetsFilename
    
    ## write original RV data with the RVs having the offsets subtracted
    filename2 = modDatasetsFilename[:-4]+"-origRVdataMinusOffsets.dat"
    f2 = open('%s' % filename2, 'w')
    
    f2.write('#'+os.path.basename(filename2)+'\n')
    #f.write('#JD          RV [m/s]    delta_RV [m/s]\n')
    ## NOTE: by this point the Jitter has all ready been
    ##       added to the RV errors in quadrature
    for dataset in range(0,len(RVsOUT)):
        f2.write("#Offset Subtracted Dataset Number "+str(dataset+1)+"\n")
        f2.write('#JD          RV [m/s]    delta_RV [m/s]     Jitter\n')
        epochs = RV_epochsIN2[dataset]
        for epoch in range(0,len(epochs)):
            line = str(RV_epochsIN2[dataset][epoch])
            line = line +'  '+ str(RVsOUT[dataset][epoch])
            line = line +'  '+ str(RV_errors[dataset][epoch])
            line = line +'  '+str(0.0)+'\n'
            f2.write(line)
            #print line
        f2.write("\n")
    f2.close()
    print 'Original RV data minus offsets written to:'+filename2
    
    ## write RV residuals to file
    filename3 = modDatasetsFilename[:-4]+"-residuals.dat"
    f3 = open('%s' % filename2, 'w')
    
    f3.write('#'+os.path.basename(filename2)+'\n')
    #f.write('#JD          RV [m/s]    delta_RV [m/s]\n')
    ## NOTE: by this point the Jitter has all ready been
    ##       added to the RV errors in quadrature
    for dataset in range(0,len(RVsOUT)):
        f3.write("#Residuals for Dataset Number "+str(dataset+1)+"\n")
        f3.write('#JD          RV [m/s]    delta_RV [m/s]     Jitter\n')
        epochs = RV_epochsIN2[dataset]
        for epoch in range(0,len(epochs)):
            line = str(RV_epochsIN2[dataset][epoch])
            line = line +'  '+ str(residuals3[0][dataset][epoch])
            line = line +'  '+ str(RV_errors[dataset][epoch])
            line = line +'  '+str(0.0)+'\n'
            f3.write(line)
            #print line
        f3.write("\n")
    f3.close()
    print 'RV residuals written to:'+filename2
    
    #######################################################################
    ### Extra plot to see ALL mod datasets on screen with original best fit
    #######################################################################
    #Get orbitRVs for best fit plot
    orbitVRs2 = []
    orbitPhases2 = []
    print '\nabout to calc fit line'#$$$$$$$$$$$$$$$$$$$$$$$$$
    for orb in range(0,len(e)):
        orbitVRs = []
        numSteps = 500.0
        periodIncrement = (period[orb]*365.242)/numSteps
        t = T_center-((period[orb]*365.242)/2.0)
        times = []
        (a_total, a1, a2, p) = semiMajorConverter(Mass1, Mass2, a_total=a[orb],a1=0.0,a2=0.0, period=period[orb])
        for step in range(0,int(numSteps)):
            t = t + periodIncrement
            times.append(t)
            # calculate the velocity residual due to the companion 
            (v_r_c,K) = rvTools.vrCalculatorStar2(t,e[orb],T[orb],period[orb],argPeri_deg[orb],a1,T_center=T_center,i=inc[orb], K=planet_K, verbose=False)
            orbitVRs.append(v_r_c)
        #print 'times were '+repr(times)
        print 'Orbit '+str(orb)+" had a K = "+str(K)
        orbitPhases = epochsToPhases(times,T_center,period[orb], verbose=False)
        #print 'orbital phases were:\n'+repr(orbitPhases)
        orbitPhases2.append(orbitPhases)
        orbitVRs2.append(orbitVRs)
        
        
    print 'len(RVs) = ',len(RVs)
    print 'len(phases3) = ',len(phases3)
    print 'len(phases3[0]) = ',len(phases3[0])
    print 'len(RVsCalcdOUT3) = ',len(RVsCalcdOUT3)
    print 'len(RVsCalcdOUT3[0]) = ',len(RVsCalcdOUT3[0])
    
    
    fig = plt.figure(1,figsize=(10,10))
    fitPlot = fig.add_subplot(211)
    residualsPlot = fig.add_subplot(212)
    residualsPlot.set_title("Residuals Plot")
    for modSetNum in range(0,len(RV_epochsIN3)):
        for dataset in range(0,len(RVs)):
            residualsPlot.scatter(phases3[0][dataset],ResidualsModOUT3[modSetNum][dataset],s=8)
            
    for orb in range(0,len(e)):
        fitPlot.plot(orbitPhases2[orb],orbitVRs2[orb],c='r')
    
    for modSetNum in range(0,len(RV_epochsIN3)):
        for dataset in range(0,len(RVs)):
            fitPlot.scatter(phases3[0][dataset],RVsCalcdOUT3[modSetNum][dataset],s=8)
            
    fitPlot.set_title("Best Fit Plot")
    
    # save plot to file
    if plotFilename!='':
        plt.savefig(plotFilename, dpi=300, orientation='landscape')
        print '\nFigure saved to: ',plotFilename
    else: 
        print '\nWARNING: NO plotFilename provided, so NOT saving it to disk.'
    # show plot
    if True:
        plt.show()
        
    plt.close()
    
def epochsToPhases(epochs2,T_center,P_yrs, verbose=False, halfOrbit=False):
    """
    converter to change the epochs of a 1-D or 2-D array to phases.
    Output dimensions of phases2 will be same as input epochs2
    """
    if verbose:
        print '\n *In epochsToPhases*\n'
        
    if type(epochs2[0])!=list:
        epochs2 = [epochs2]
        notAlist = True
    else:
        notAlist = False
        
    phases2 = []
    P_days = P_yrs*365.242 # converted as t and T are in JD
    for epochs in epochs2:
        phases = []
        for epoch in epochs:
            phaseTimeDiff = epoch - int((epoch-T_center)/P_days)*P_days-T_center
            phase = phaseTimeDiff/P_days
            if halfOrbit:
                if phase>0.5:
                    phase = phase-1.0
                elif phase<-0.5:
                    phase = phase+1.0
            phases.append(phase)
            if verbose:
                print '\nepochOut = ',epoch
                print 'period [days] = ',P_days
                print 'phase = ',phase
        phases2.append(phases)
    
    if notAlist:
        epochs2 = epochs2[0]
        phases2 = phases2[0]
    
    if verbose:
        print '\nepochs2 array = '+repr(epochs2)
        print 'phases2 array = '+repr(phases2)
        
    return phases2
            
        
def rvFitPlotter1Body(e, T_lastPeri, period, inc, argPeri_deg, a=0.0, T_transitCenter=0, RVoffsets=[0],\
              plotFilename='', show=True, plotFullOrbit=True):
    """
    create a plot for the RV data and a fit line from the best orbit data params.
    
    NOTE: e, T, period, inc, argPeri_deg, a can all be lists, 
          BUT T_center will be the same for all the orbits plotted.
    
    """
    T = T_lastPeri
    T_center = T_transitCenter
    if T_center==0:
        T_center = T
    
    #check the orbit element inputs to ensure they are lists, else make them lists    
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
    if type(T)!=list:
       T = [T]
    if type(RVoffsets)!=list:
       RVoffsets = [RVoffsets]
       
    cwd = os.getcwd()
    if '/Toolboxes'==cwd[-10:]:
        os.chdir(cwd[:-9])
        print '\nTemporarily changed cwd to '+os.getcwd()
    from Dicts.HARDCODEDsettingsDict import HARDCODEDsettingsDict
    if 'DualLanguageSimCode/Toolboxes' in cwd:
        os.chdir(cwd)
        print '\nChanged cwd back to '+os.getcwd()
    settings_and_InputDataDir = HARDCODEDsettingsDict['settings_and_InputDataDir']
    # get the input settings file name with dir, and RV and system data dicts
    inputSettingsFile = os.path.join(settings_and_InputDataDir,'SimSettings.txt')
    paramSettingsDict = cFileToSimSettingsDict(inputSettingsFile, inputSettingsFile)
    sysDatafilename = os.path.join(settings_and_InputDataDir,'SystemData.txt')
    sysDataDict = sysDataToDict(sysDatafilename)
    RVdatafilename = os.path.join(settings_and_InputDataDir,'RVdata.dat')
    RVdataDict = rvTools.RVdataToDict(RVdatafilename)
    
    ## Pull out required RV data
    RVs = RVdataDict['RVs']
    RV_epochs = RVdataDict['RV_epochs']
    RV_errors = RVdataDict['RV_errors']
    
    # get the companions mass depending on simulation case
    Mass1 = sysDataDict['Mass1']
    simulate_StarStar = paramSettingsDict["simulate_StarStar"]
    simulate_StarPlanet = paramSettingsDict["simulate_StarPlanet"]
    if (simulate_StarStar is True) and(simulate_StarPlanet is True):
        print "Error: simulate_StarStar and simulate_StarPlanet can NOT BOTH be True!"
        
    star_Mass2 = sysDataDict["star_Mass2"]
    planet_MsinI = sysDataDict["planet_MsinI"]   #[Msun]
    if planet_MsinI==0:
        planet_MsinI=False
    if simulate_StarStar:
        Mass2 = star_Mass2
    if simulate_StarPlanet:
        Mass2 = planet_MsinI#/math.sin(math.radians(inc[0]))
    
    if plotFilename!='':
        if plotFilename[-4:]!='.png':
            plotFilename = plotFilename+'.png'     
    
    ## make an advanced title for plot from folder and filename
    titleTop = os.path.dirname(outputDataFilename).split('/')[-1]
    titleBtm = os.path.basename(plotFilename).split('.')[0]+" Best Fit Plot"
    plotFileTitle = titleTop+'\n'+titleBtm
    
    # make copy of input RVs (not sure if I still need to do this...)
    RV_epochsIN2 = RV_epochs
    
    ## now convert the new epochs to phases.
    phases3 = []
    for orb in range(0,len(e)):
        print '\nConverting new epochs to phases for orbit '+str(orb)
        phases2 = epochsToPhases(RV_epochsIN2,T_center,period[orb], verbose=False, halfOrbit=True)         
        phases3.append(phases2)   
    
    RVsIN = RVs
    RVsOUT = []
    for dataset in range(0,len(RVsIN)):
        RVsUSE=[]
        if len(RVoffsets)>=(dataset+1):
            offset = RVoffsets[dataset]
        else:
            offset = 0
        for epoch in range(0,len(RVsIN[dataset])):
            vel = RVsIN[dataset][epoch]-offset
            RVsUSE.append(vel)
        RVsOUT.append(RVsUSE)
    
    #Get residual values for RV scatter plot
    residuals3 = []
    for orb in range(0,len(argPeri_deg)):
        residuals2 = []
        chiSquaredTot2 = 0
        numEpochs_RV = 0
        chi_squared_RV_reducedCur2 = []
        for dataset in range(0,len(RVsOUT)):
            chiSquaredTot = 0
            epochs = RV_epochsIN2[dataset]
            rvs = RVsOUT[dataset]
            errors = RV_errors[dataset]
            if False:
                print 'len(RVsOUT) = ',len(RVsOUT)
                print 'len(RV_epochsOUT3[orb]) = ',len(RV_epochsIN2)
                print 'len(RVsOUT[dataset]) = ',len(RVsOUT[dataset])
                print 'len(RV_epochsOUT3[orb][dataset]) = ',len(RV_epochsIN2[dataset])
            
            residuals = []
            (a_total, a1, a2, p) = semiMajorConverter(Mass1, Mass2, a_total=a[orb],a1=0.0,a2=0.0, period=period[orb],verbose=True)
            
            #print '\nFor dataset '+str(dataset)
            for epoch in range(0,len(epochs)):
                if False:
                    print '\nWorking on epoch ',(epoch+1)
                # calculate the velocity residual due to the companion star
                (v_r_c,K) = rvTools.vrCalculatorStar2(epochs[epoch],e[orb],T[orb],period[orb],argPeri_deg[orb],a1,T_center=T_center, i=inc[orb], K=False, verbose=False)               
                RV = rvs[epoch]-v_r_c
                residuals.append(RV)
                chiSquaredCurr = chiSquaredCalc(rvs[epoch], errors[epoch], v_r_c)
                chiSquaredTot = chiSquaredTot+chiSquaredCurr
                if False:
                    print 'RV = '+str(rvs[epoch])+' - ('+str(v_r_c)+') = '+str(RV)
                    print 'chiSquaredCurr = ',chiSquaredCurr
                    print 'chiSquaredTot = ',chiSquaredTot
            chi_squared_RV_reducedCur = (1.0/((1.0*len(rvs))-6.0))*chiSquaredTot
            chi_squared_RV_reducedCur2.append(chi_squared_RV_reducedCur)
            chiSquaredTot2 = chiSquaredTot2 + chiSquaredTot
            if orb==0:
                numEpochs_RV_curr = len(epochs)
                numEpochs_RV = numEpochs_RV+numEpochs_RV_curr
            residuals2.append(residuals)
        residuals3.append(residuals2) 
        
    # convert to a reduced chiSquared
    chi_squared_RV_reduced = (1.0/((1.0*numEpochs_RV)-6.0))*chiSquaredTot2
    if True:
        print '\nOriginal ChiSquared = ',chiSquaredTot2
        print 'Num RV epochs = ',numEpochs_RV
        print 'one over nu = ',(1.0/((1.0*numEpochs_RV)-6.0))
        print 'reduced chiSqured RV = ',chi_squared_RV_reduced
    
    #Get orbitRVs for best fit plot
    orbitVRs2 = []
    orbitPhases2 = []
    for orb in range(0,len(e)):
        orbitVRs = []
        numSteps = 5000.0
        periodIncrement = (period[orb]*365.242)/numSteps
        t = T_center-((period[orb]*365.242)/2.0)
        times = []
        (a_total, a1, a2, p) = semiMajorConverter(Mass1, Mass2, a_total=a[orb],a1=0.0,a2=0.0, period=period[orb], verbose=True)
        for step in range(0,int(numSteps)):
            t = t + periodIncrement
            times.append(t)
            # calculate the velocity residual due to the companion 
            (v_r_c,K) = rvTools.vrCalculatorStar2(t,e[orb],T[orb],period[orb],argPeri_deg[orb],a1,T_center=T_center,i=inc[orb], K=False, verbose=False)
            orbitVRs.append(v_r_c)
        #print 'times were '+repr(times)
        print 'Orbit '+str(orb)+" had a K = "+str(K)
        orbitPhases = epochsToPhases(times,T_center,period[orb], verbose=False)
        #print 'orbital phases were:\n'+repr(orbitPhases)
        orbitPhases2.append(orbitPhases)
        orbitVRs2.append(orbitVRs)
    
    
    ## Create figure, axes and start plotting/drawing everything
    fig = plt.figure(1,figsize=(10,10))
    residualsPlot = fig.add_subplot(212)
    residualsPlot.set_title("Residuals Plot")
    plt.suptitle(plotFileTitle)
    colorsList = ['b','m','k','g','y']
    
    for orb in range(0,len(argPeri_deg)):
        meanMedianChiSquaredSTR = ''
        residuals3[orb]
        for dataset in range(0,len(RVs)):
            meanMedianChiSquaredSTR = meanMedianChiSquaredSTR+'\nmean of set '+str(dataset)+' = '+str(np.mean(residuals3[orb][dataset]))
            meanMedianChiSquaredSTR = meanMedianChiSquaredSTR+'\nmedian of set '+str(dataset)+' = '+str(np.median(residuals3[orb][dataset]))
            meanMedianChiSquaredSTR = meanMedianChiSquaredSTR+'\nchiSquared_reduced of set '+str(dataset)+' = '+str(chi_squared_RV_reducedCur2[dataset])
            residualsPlot.scatter(phases3[orb][dataset],residuals3[orb][dataset],s=8,edgecolor=colorsList[dataset],facecolor=colorsList[dataset])
            for epoch in range(0,len(phases3[orb][dataset])):
                xs = [phases3[orb][dataset][epoch],phases3[orb][dataset][epoch]]
                ys = [residuals3[orb][dataset][epoch]-abs(RV_errors[dataset][epoch]),residuals3[orb][dataset][epoch]+abs(RV_errors[dataset][epoch])]
                residualsPlot.plot(xs,ys,c='k')
                
    meanMedianChiSquaredSTR = meanMedianChiSquaredSTR+'\nchiSquared_reduced for ALL data is = '+str(chi_squared_RV_reduced)
    print meanMedianChiSquaredSTR
    xmin = findArrayMin(phases3)
    xmax = findArrayMax(phases3)
    ymin = findArrayMin(residuals3)
    ymax = findArrayMax(residuals3)
    print '\nxmin = '+repr(xmin)
    print 'xmax = '+repr(xmax)
    xLim =[xmin,xmax]
    xrange = xmax-xmin
    xLim2 = (xmin-abs(xrange*0.05), xmax+abs(xrange*0.05))
    yLim = [ymin,ymax]
    yrange = ymax-ymin
    yLim2 = (ymin-abs(yrange*0.05), ymax+abs(yrange*0.05))
    
    ## make plot of fit to data
    fitPlot = fig.add_subplot(211)
    fitXmin = findArrayMin(orbitPhases2[0])
    if xmin<fitXmin:
        fitXmin = xmin
    fitXmax = findArrayMax(orbitPhases2[0])
    if xmax>fitXmax:
        fitXmax = xmax
    fitYmin = findArrayMin(orbitVRs2[0])
    if ymin<fitYmin:
        fitYmin = ymin
    fitYmax = findArrayMax(orbitVRs2[0])
    if ymax>fitYmax:
        fitYmax = ymax
        
    yMaxRVs = findArrayMax(RVsOUT)
    yMinRVs = findArrayMin(RVsOUT)
    yRVsRange = yMaxRVs-yMinRVs
    yLimRVs = [yMinRVs-abs(yRVsRange*0.05),yMaxRVs+abs(yRVsRange*0.05)]
        
    fitXrange = fitXmax-fitXmin
    fitXLim2 = (fitXmin-abs(fitXrange*0.5), fitXmax+abs(fitXrange*0.05))
    fitYrange = fitYmax-fitYmin
    fitYLim2 = (fitYmin-abs(fitYrange*0.05), fitYmax+abs(fitYrange*0.05))
    if plotFullOrbit:
        fitXLimsUSE = fitXLim2
        fitYLimsUSE = fitYLim2
    else:
        fitXLimsUSE = xLim2
        fitYLimsUSE = yLimRVs
        
    residualsPlot.plot(fitXLimsUSE,[0,0],c='r')
    residualsPlot.axes.set_xlim(fitXLimsUSE)
    residualsPlot.axes.set_ylim(yLim2)
    fitPlot.axes.set_xlim(fitXLimsUSE)
    fitPlot.axes.set_ylim(fitYLimsUSE)
    
    for orb in range(0,len(e)):
        ## plot RV fit line
        fitPlot.plot(orbitPhases2[orb],orbitVRs2[orb],c='r')
        ## plot RV data
        for dataset in range(0,len(RVs)):
            fitPlot.scatter(phases3[orb][dataset],RVsOUT[dataset],s=8,edgecolor=colorsList[dataset],facecolor=colorsList[dataset])
            for epoch in range(0,len(phases3[orb][dataset])):
                xs = [phases3[orb][dataset][epoch],phases3[orb][dataset][epoch]]
                ys = [RVsOUT[dataset][epoch]-abs(RV_errors[dataset][epoch]),RVsOUT[dataset][epoch]+abs(RV_errors[dataset][epoch])]
                fitPlot.plot(xs,ys,c='k')
                
    paramsLegndStr = 'e = '+str(e[0])+'\nperiod [days] = '+str(period[0]*365.242)+'\nTo = '+str(T[0])+", Tc = "+str(T_center)+\
                '\ninc = '+str(inc[0])+\
                '\nargPeri_deg = '+str(argPeri_deg[0])+'\na = '+str(a[0])
    for dataset in range(0,len(RVs)):
        if dataset==0:
            paramsLegndStr=paramsLegndStr+'\n'
        else:
            paramsLegndStr=paramsLegndStr+', '
        paramsLegndStr = paramsLegndStr+'offset '+str(dataset)+'= '+str(RVoffsets[dataset])

    
    addLegend=True
    if addLegend:
        # in bottom right
        #residualsPlot.text(fitXLimsUSE[1]-abs(fitXLimsUSE[1]*0.002),abs(fitYLimsUSE[1]*0.5),paramsLegndStr,ha='left')
        # in bottom left
        paramsLegndY = fitYLimsUSE[0]+abs(fitYLimsUSE[0]*0.1)
        paramsLegndX = fitXLimsUSE[0]+abs(fitXLimsUSE[0]*0.03)
        residsLegendY = yLim2[0]+abs(yLim2[0]*0.02)
        residsLegendX = paramsLegndX
        fitPlot.text(paramsLegndX,paramsLegndY,paramsLegndStr,ha='left')
        residualsPlot.text(residsLegendX,residsLegendY,meanMedianChiSquaredSTR,ha='left')
        print 'Params legend bottom left corner at [ '+str(paramsLegndX)+' , '+str(paramsLegndY)+' ]'
        print 'Residuals legend bottom left corner at [ '+str(residsLegendX)+' , '+str(residsLegendY)+' ]'
        
    # save plot to file
    if plotFilename!='':
        plt.savefig(plotFilename, dpi=300, orientation='landscape')
        print '\nFigure saved to: ',plotFilename
    else: 
        print '\nWARNING: NO plotFilename provided, so NOT saving it to disk.'
    # show plot
    if show:
        plt.show()
        
    plt.close()
def singleParamReadAndPlotNEW2(filename, paramColNum, subPlot, numBins, weight=True, confLevels=True, normalize=True, drawLine=False, nu=1):
    """
    This is the single column/param plotter for the new Master Hist plotter
    that will call this in a loop to fill up a summary plot.  This is being
    done as there are RAM issues of trying to load all the data for a particular
    column and plot it all at once.  Thus, this one will load the data 
    only one bin of the hist at a time to be the most efficient with the RAM.
    
    file format must follow the 'NEW' style, ie:
    longAN [deg]      e [N/A]       To [julian date]  period [yrs]   inclination [deg]   argPeri [deg]   a_total [AU]    chiSquared    timesBeenHere
    
    All filename checks will occur in the Master function and are thus not needed here.
    
    All 3D data should have its chiSquared values all ready be reduced.
    """
    verboseInternal = True
    # set y axis title
    subPlot.axes.set_ylabel('Probability')
    axs = subPlot.axis()
    
    # get data and conf levels
    (data,confLevels) =genTools.confLevelFinderNEWdataVersion(filename,paramColNum,True)

    # done looping through all bins
    largestBinVal = np.max(binVals)
    #print 'largestBinVal = ',largestBinVal
    #print 'colors: '+repr(colors)
    #print 'binMins: ',repr(binMins)
    #print 'binMaxs: ',repr(binMaxs)
    
    binVals2 = []
    ## convert the bin information into rectangular patches
    for bin in range(0,len(binVals)):
        x = binMins[bin]
        y = 0.0
        width = binWidth
        if normalize:
            try:
                height = float(binVals[bin])/float(largestBinVal)
            except:
                #print 'Using default height as largestBinVal=0.'
                height = 1.0
        else:
            height = binVals[bin]
        binVals2.append(height)
        rec = patches.Rectangle((x,y), width, height, facecolor=colors[bin], fill=True, edgecolor='k')
        #RecPatches.append(rec)
        subPlot.add_patch(rec)
    
    # add a line to plot the top curves of the hist
    if drawLine:
        subPlot.plot(binMids,binVals2, 'r',linewidth=1)
    # add a vertical line at the median value
    med = np.median(binMids)
    subPlot.plot([med, med],subPlot.axes.get_ylim(),'k')
    
    # return fully updated sub plot
    return subPlot
    
def singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight=True, confLevels=True, normalize=True, drawLine=False, nu=1):
    """
    This is the single column/param plotter for the new Master Hist plotter
    that will call this in a loop to fill up a summary plot.  This is being
    done as there are RAM issues of trying to load all the data for a particular
    column and plot it all at once.  Thus, this one will load the data 
    only one bin of the hist at a time to be the most efficient with the RAM.
    
    It is assumed that the chiSquareds in the data file are reduced, so the nu value will de-reduce 
    them to calculated proper likelihoods for weighting.
    
    file format must follow the 'NEW' style, ie:
    longAN [deg]      e [N/A]       To [julian date]  period [yrs]   inclination [deg]   argPeri [deg]   a_total [AU]    chiSquared    timesBeenHere
    
    All filename checks will occur in the Master function and are thus not needed here.
    
    All 3D data should have its chiSquared values all ready be reduced.
    """
    verboseInternal = False
    # set y axis title
    subPlot.axes.set_ylabel('Probability')
    axs = subPlot.axis()
    
    ## First get ranges of param and ChiSquared values
    if os.path.exists(filename):
        print 'Opening and finding ranges for data in column # '+str(paramColNum)
        
        f = open(filename,'r')
        plotFileTitle = f.readline()[:-5]
        headings = f.readline()
        dataMax = 0.1
        dataMin = 10000000.0
        chiSquaredMin = 1000000.0
        bestDataVal = 0
        line = 'asdf'
        totalAccepted = 0
        while line!='':
            line = f.readline()
            if line!='':
                totalAccepted+=1
                dataLineCols = line.split()
                dataValue = float(dataLineCols[paramColNum])
                chiSquared = float(dataLineCols[7])
                if dataValue>dataMax:
                    dataMax = dataValue
                if dataValue<dataMin:
                    dataMin = dataValue
                if chiSquared<chiSquaredMin:
                    chiSquaredMin = chiSquared
                    bestDataVal = dataValue
        f.close()
        
        ## set up axis of plot to match data
        xAxisWidth = dataMax-dataMin
        if verboseInternal:
            print '\nOriginal values:'
            print 'dataMin = ',dataMin
            print 'dataMax = ',dataMax
            print 'xAxisWidth = ',xAxisWidth
        ##########################################
        # special stuff for To plot
        xSubtractValueString = ''
        xSubtractValue = 0
        if dataMin>2e6:
            # assume the data is for time of periapsis
            Divisor=1e9
            if xAxisWidth>10:
                Divisor = 10
            if xAxisWidth>100:
                Divisor = 100
            if xAxisWidth>1000:
                Divisor = 1000
            if xAxisWidth>5000:
                Divisor = 5000
            if xAxisWidth>10000:
                Divisor = 10000        
            xSubtractValue = int(int(dataMin/Divisor)*Divisor)
            xSubtractValueString = " + "+str(xSubtractValue)
            xLabelStr = subPlot.axes.get_xlabel()
            subPlot.axes.set_xlabel(xLabelStr+xSubtractValueString)
            dataMin = dataMin-xSubtractValue
            dataMax = dataMax-xSubtractValue
            xAxisWidth = dataMax-dataMin
            if verboseInternal:
                print 'xLabelStr = ',xLabelStr
                print '\nUpdated values:'
                print 'Divisor = ',Divisor
                print 'dataMin = ',dataMin
                print 'dataMax = ',dataMax
                print 'xAxisWidth = ',xAxisWidth
                print 'xLabelStr+xSubtractValueString = ',xLabelStr+xSubtractValueString            
        ###########################################
        axisXmin = dataMin-(0.05*xAxisWidth)
        axisXmax = dataMax+(0.05*xAxisWidth)
        axisYmin = 0.0
        axisYmax = 1.1
        subPlot.axis([axisXmin, axisXmax, axisYmin, axisYmax])
        
        #print 'chiSquaredMin = ',chiSquaredMin
        
        ## set up data range using a small buffer
        dataRange = (dataMax-dataMin)*1.05
        dataMin2 = dataMin-(dataRange*0.025)
        dataMax2 = dataMax+(dataRange*0.025)
        binWidth = dataRange/numBins
        if verboseInternal:
            print 'binWidth = ',binWidth
            print 'dataMin2 = ',dataMin2
            print 'dataMax2 = ',dataMax2
            
        nextBinMin = dataMin2
        nextBinMax = dataMin2+binWidth
        binMins = []
        binMaxs = []
        binVals = []
        binMids = []
        binMidsALL = []
        colors = []
        RecPatches = []
        
        for bin in range(1,numBins+1):
            if (totalAccepted>1e6):
                print "Looping through to load up bin # "+str(bin)
            # update bin min and maxs
            binMins.append(nextBinMin)
            binMaxs.append(nextBinMax)
            binMids.append((nextBinMax+nextBinMin)/2.0)
            #print 'currBinMin = ',currBinMin
            #print 'currBinMax = ',currBinMax
            binVal = 0
            f = open(filename,'r')
            plotFileTitle = f.readline()[:-5]
            headings = f.readline()
            #curBinChiSquaredMAX = 0.0
            curBinChiSquaredMIN = 10000.0
            line = 'asdf'
            while line!='':
                line = f.readline()
                if line!='':
                    dataLineCols = line.split()
                    dataValue = float(dataLineCols[paramColNum])-xSubtractValue
                    chiSquared = float(dataLineCols[7])
                    timesBeenHere = float(dataLineCols[-1])
                    
                    #print 'chiSquared ',chiSquared #$$$$$$$$$$$
                    likelihood = 0
                    # see if value is inside bin, not including max
                    if bin==(numBins):
                        # take care of last bin including it's max value
                        if (dataValue>=binMins[-1]) and (dataValue<=binMaxs[-1]):
                            if weight:
                                likelihood = math.exp(-nu*chiSquared/2.0)#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                            else:
                                likelihood = 1
                            #if chiSquared>curBinChiSquaredMAX:
                                #print 'updating chiSquaredMAX to ',chiSquared
                                #curBinChiSquaredMAX = chiSquared 
                            if chiSquared<curBinChiSquaredMIN:
                                #print 'updating chiSquaredMIN to ',chiSquared
                                curBinChiSquaredMIN = chiSquared 
                    else:
                        # still not last bin, so don't include max value
                        if (dataValue>=binMins[-1]) and (dataValue<binMaxs[-1]):
                            if weight:
                                likelihood = math.exp(-nu*chiSquared/2.0)#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                            else:
                                likelihood = 1
                            #if chiSquared>curBinChiSquaredMAX:
                                #print 'updating chiSquaredMax to ',chiSquared
                                #curBinChiSquaredMAX = chiSquared 
                            if chiSquared<curBinChiSquaredMIN:
                                #print 'updating chiSquaredMIN to ',chiSquared
                                curBinChiSquaredMIN = chiSquared 
                    # update bin value with which ever case was satisfied
                    binVal = binVal + (likelihood*timesBeenHere)
                    
                    # add bin mid value to total list to find median later
                    if likelihood>0.0:
                        binMidsALL.append(binMids[-1])  ###$$$$$$$$$$$$$$ could be dangerous for long data sets!!!!!!!!!!!!!!!!!!!!!!
            f.close()
            
            if confLevels:
                #print 'curBinChiSquaredMax', curBinChiSquaredMax #$$$$$$$$$$$$$$$$$$
                #print 'chiSquaredMin',chiSquaredMin  #$$$$$$$$$$$$$$$$$$
                # is there anything in the bin?
                if binVal>0.0:
                    c = 'w'
                    ## set color of bin depending on what conf level the data is in
                    # Check if the data is in the 2sigma conf level
                    if (chiSquaredMin+4.0)>curBinChiSquaredMIN:
                        c = '0.8'
                    # Check if the data is in 1sigma conf level
                    if (chiSquaredMin+1.0)>curBinChiSquaredMIN:
                        c = '0.5'
                        
                # nothing in bin, so make it white
                else:
                    c = 'w'
            # don't do confLevels, so make them all blue
            else:
                c = 'b'     
            
            colors.append(c)
            #print 'colors = '+repr(colors)#$$$$$$$$$$$$$$$$
            binVals.append(binVal)
            
            nextBinMin = nextBinMax
            nextBinMax = nextBinMax+binWidth
            #print 'curr binVal = ',binVal
            
            #done looping through the data
            #print 'final binVal = ',binVal
            #print 'curBinChiSquaredMax = ',curBinChiSquaredMax
        
        # done looping through all bins
        largestBinVal = np.max(binVals)
        #print 'largestBinVal = ',largestBinVal
        #print 'colors: '+repr(colors)
        #print 'binMins: ',repr(binMins)
        #print 'binMaxs: ',repr(binMaxs)
        
        binVals2 = []
        ## convert the bin information into rectangular patches
        for bin in range(0,len(binVals)):
            x = binMins[bin]
            y = 0.0
            width = binWidth
            if normalize:
                try:
                    height = float(binVals[bin])/float(largestBinVal)
                except:
                    #print 'Using default height as largestBinVal=0.'
                    height = 1.0
            else:
                height = binVals[bin]
            binVals2.append(height)
            rec = patches.Rectangle((x,y), width, height, facecolor=colors[bin], fill=True, edgecolor='k')
            #RecPatches.append(rec)
            subPlot.add_patch(rec)
        
        # add a line to plot the top curves of the hist
        if drawLine:
            subPlot.plot(binMids,binVals2, 'r',linewidth=1)
        # add a vertical line at the median value
        med = np.median(binMidsALL)
        subPlot.plot([med, med],subPlot.axes.get_ylim(),'k')
        subPlot.plot([bestDataVal,bestDataVal],subPlot.axes.get_ylim(),'g')
        # return fully updated sub plot
        return subPlot
    
    else:
        print 'Filename: '+filename+'\n Does NOT exist!!!'

def dataReadAndPlotNEW3(filename, plotFilename, weight=True, confLevels=True, normalize=True, drawLine=False, nu=1, plot4x1=False):
    """
    Master plotting function for the singleParamReadAndPlotNEW
    
    file format must follow the 'NEW' style, ie:
    longAN [deg]      e [N/A]       To [julian date]  period [yrs]   inclination [deg]   argPeri [deg]   a_total [AU]    chiSquared    timesBeenHere
    
    All 3D data should have its chiSquared values all ready be reduced.
    """
    verbose = False
    
    numBins = 50
    # record the time the chain started
    startTime = timeit.default_timer()
    
    # check if the passed in value for filename includes '.txt'
    if (filename[-4:]!='.txt' and filename[-4:]!='.dat'):
        filename = filename+'.dat'
        print ".dat was added to the end of the filename as it wasn't found to be in the original version"
    
    datadir = os.path.dirname(filename)
    logFilename = os.path.join(datadir,'processManagerLogFile.txt')
    log = open(logFilename,'a')
    log.write('\n'+75*'#'+'\n Inside dataReadAndPlotNEW3 \n'+75*'#'+'\n')
    
    ## find number of RV datasets
    if os.path.exists(filename):
        f = open(filename,'r')
        plotFileTitle = f.readline()[:-5]
        headings = f.readline()
        line = f.readline()
        dataLineCols = line.split()
        if len(line)>7:
            numRVdatasets = len(dataLineCols) - 9
        else:
            line = f.readline()
            dataLineCols = line.split()
            if len(line)>7:
                numRVdatasets = len(dataLineCols) - 9
        s= "\nNumber of RV datasets found in dataReaderAndPlotNEW3 was "+str(numRVdatasets)+"\n"
        print s
        log.write(s+'\n')
        f.close()
        
    # check if the passed in value for plotFilename includes '.png'
    if '.png' not in plotFilename:
        plotFilename = plotFilename+'.png'
    
    ## make an advanced title for plot from folder and filename
    titleTop = os.path.dirname(outputDataFilename).split('/')[-1]
    titleBtm = os.path.basename(plotFilename).split('.')[0]+" Best Fit Plot"
    plotFileTitle = titleTop+'\n'+titleBtm

    # Call the function to find the best orbit values for reference.
    if verbose:
        print '#'*50
        bestOrbitFinderNEW(filename, printToScreen=False, saveToFile=True)
        print '#'*50
    
    s= '\nStarting to create summary plot\n'
    print s
    log.write(s+'\n')
    
    # Create empty figure to be filled up with plots
    if not plot4x1:
        fig = plt.figure(1, figsize=(25,10) ,dpi=300) 
    else:
        fig = plt.figure(1, figsize=(15,20) ,dpi=300)
        
    # Create sub plot and fill it up for the Inclinations
    if not plot4x1:
        subPlot = fig.add_subplot(231)
        paramColNum = 4
        subPlot.axes.set_xlabel('Inclination [deg]')
        subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, nu)
        s= "Done plotting inclination_degsAlls"
        print s
        log.write(s+'\n')
        
    # Create sub plot and fill it up for the Longitude of Ascending Node
    if not plot4x1:
        subPlot = fig.add_subplot(232)
        paramColNum = 0
        subPlot.axes.set_xlabel('Longitude of Ascending Node [deg]')
        subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, nu)
        s= "Done plotting longAN_degsAlls"
        print s
        log.write(s+'\n')
        
    # Create sub plot and fill it up for the Argument of Perigie
    if not plot4x1:
        subPlot = fig.add_subplot(233)
    else:
        subPlot = fig.add_subplot(312)
    paramColNum = 5
    subPlot.axes.set_xlabel('Argument of Perigie [deg]')
    subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, nu)
    s= "Done plotting argPeri_degsAlls"
    print s
    log.write(s+'\n')
    
    # Create sub plot and fill it up for the e
    if not plot4x1:
        subPlot = fig.add_subplot(234)
    else:
        subPlot = fig.add_subplot(311)
    paramColNum = 1
    subPlot.axes.set_xlabel('e')
    subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, nu)
    s= "Done plotting esAlls"
    print s
    log.write(s+'\n')
    
    # Create sub plot and fill it up for the Period
    if not plot4x1:
        subPlot = fig.add_subplot(235)
        paramColNum = 3
        subPlot.axes.set_xlabel('Period [Years]')
        subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, nu)
        s= "Done plotting periodsAlls"
        print s
        log.write(s+'\n')
        
#    # Create sub plot and fill it up for the Semi-major
#    subPlot = fig.add_subplot(336)
#    paramColNum = 6
#    subPlot.axes.set_xlabel('Semi-major Axis [AU]')
#    subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, nu)
#    print "Done plotting Semi-major Axis"
    
    # Create sub plot and fill it up for the Time of last Periapsis
    if not plot4x1:
        subPlot = fig.add_subplot(236)
    else:
        subPlot = fig.add_subplot(313)
    paramColNum = 2
    subPlot.axes.set_xlabel('Time of last Periapsis [JD]')
    subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, nu)
    s= "Done plotting TsAlls"
    print s
    log.write(s+'\n')
        
    # Save file 
    plt.savefig(plotFilename, dpi=300, orientation='landscape')
    s= 'Summary figure saved to '+plotFilename
    print s
    log.write(s+'\n') 
    plt.close()
    
    ## Create a second figure of RV offsets. ####
    try:
        # Create empty figure to be filled up with plots
        # Create sub plot and fill it up for the Semi-major
        if numRVdatasets==1:
            fig = plt.figure(2, figsize=(35,15) ,dpi=250)
        elif numRVdatasets==2:
            fig = plt.figure(2, figsize=(35,22) ,dpi=250)
        elif numRVdatasets==3:
            fig = plt.figure(2, figsize=(35,30) ,dpi=250)
        
        
        if numRVdatasets>=1:
            # Create sub plot and fill it up for the Semi-major
            if numRVdatasets==1:
                subPlot = fig.add_subplot(111)
            elif numRVdatasets==2:
                subPlot = fig.add_subplot(211)
            elif numRVdatasets==3:
                subPlot = fig.add_subplot(311)
            paramColNum = 8
            subPlot.axes.set_xlabel('RV offset 1 [m/s]')
            subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, nu)
            s= "Done plotting RV offsets for dataset 1"
            print s
            log.write(s+'\n')
        if numRVdatasets>=2:
            # Create sub plot and fill it up for the Semi-major
            if numRVdatasets==2:
                subPlot = fig.add_subplot(212)
            elif numRVdatasets==3:
                subPlot = fig.add_subplot(312)
            paramColNum = 9
            subPlot.axes.set_xlabel('RV offset 2 [m/s]')
            subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, nu)
            s= "Done plotting RV offsets for dataset 2"
            print s
            log.write(s+'\n')
        if numRVdatasets==3:
            # Create sub plot and fill it up for the Semi-major
            subPlot = fig.add_subplot(313)
            paramColNum = 10
            subPlot.axes.set_xlabel('RV offset 3 [m/s]')
            subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, nu)
            s= "Done plotting RV offsets for dataset 3"
            print s
            log.write(s+'\n')  
        # Save file 
        plotFilename2 = plotFilename[0:-4]+'-RVoffsets.png'
        plt.savefig(plotFilename2, dpi=300, orientation='landscape')
        s= 'RV offsets summary figure saved to '+plotFilename2
        print s
        log.write(s+'\n') 
        plt.close()
    except:
        plt.close()
        s= 'No RV offsets to plot'
        print s
        log.write(s+'\n')
        
    # record the time the chain finished and print
    endTime = timeit.default_timer()
    totalTime = (endTime-startTime) # in seconds
    totalTimeString = genTools.timeString(totalTime)
    s= '\n\ndataReadAndPlotNEW3: Plotting took '+totalTimeString+' to complete.\n'
    print s
    log.write(s+'\n')
    log.write('\n'+75*'#'+'\n Leaving dataReadAndPlotNEW3 \n'+75*'#'+'\n')
    log.close()

