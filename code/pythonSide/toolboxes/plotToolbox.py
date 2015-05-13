#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import math
#import gc
import numpy as np
import os
import pylab
import timeit
from math import pi
#import dicts
import pickle
from progressbar import ProgressBar
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
    the= np.linspace(0,2*pi,Nb)
    X=radm*np.cos(the)*co-si*radn*np.sin(the)+xpos
    Y=radm*np.cos(the)*si+co*radn*np.sin(the)+ypos
    return X,Y    

def fixPlotBordersAndLabels(plot):
    """
    This will just increase the boarder thicknesses and label sizes to look better for printing or 
    putting in papers.
    """
    plot.tick_params(axis='both',which='major',width=1,length=3,pad=8,direction='in',labelsize=20)
    #plot.axhline(linewidth=2.0)
    #plot.axvline(linewidth=2.0)
    #plot.xaxis.set_tick_params(width=5,length=5,pad=2,direction='in')
    #plot.yaxis.set_tick_params(width=5,length=5,pad=2,direction='in')
    plot.spines['right'].set_linewidth(1.0)
    plot.spines['bottom'].set_linewidth(1.0)
    plot.spines['top'].set_linewidth(1.0)
    plot.spines['left'].set_linewidth(1.0)
    return plot
    
def histMakeAndDump(chiSquareds,data,outFilename='',nbins=50,weight=False, normed=False, nu=1,logY=False,histType='bar'):
    """
    This will make a matplotlib histogram using the input settings, then writing the resulting  
    centers of the bins and number of data points in said bin values to disk, with '.dat' extension
    if not added already.
    
    This function is designed to work with a follow up like histLoadAndPlot_** to produce publication worthy plots.
    """
    if outFilename[-4]!='.dat':
        outFilename=outFilename+'.dat'
        
    if weight:
        theWeights = genTools.likelihoodsCalc(chiSquareds, nu)
    else:
        theWeights = np.ones(len(data))      
    fig = plt.figure(1)
    subPlot = fig.add_subplot(111)
    (n,bins,rectangles)=subPlot.hist(data, bins=nbins, normed=False, weights=theWeights,linewidth=7,histtype=histType,log=logY, fill=False)
    #find center of bins
    if type(bins)!=np.ndarray:
        bins = np.array(bins)
    binCenters = (bins[1:]+bins[:-1])/2.0
    histData=np.zeros((len(n),2))
    histData[:,0]=binCenters
    histData[:,1]=n
    np.savetxt(outFilename,histData)
    if True:
        print "output dat file:\n"+outFilename
    
def histLoadAndPlot_StackedPosteriors(plot,outFilename='',xLabel='X',lineColor='k',xLims=False,latex=False):
    """
    Loads previously plotted histograms that were written to disk by histPlotAndDump, and plot them up 
    in a way that is ready for publication.  This version is to plot a posterior of the same parameter 
    for multiple simulation runs to see how they differ.
    
    It is foreseen that many versions of this function will exist for different specific publication ready plots.
    """
    if outFilename[-4]!='.dat':
        outFilename=outFilename+'.dat'
    histData = np.loadtxt(outFilename)  
    ys=xs=[]
    maxN = np.max(histData[:,1])
    halfBinWidth = (histData[1][0]-histData[0][0])/2.0
    for i in range(0,histData.shape[0]):
        ys.append(histData[i][1]/maxN)
        ys.append(histData[i][1]/maxN)
        xs.append(histData[i][0]-halfBinWidth)
        xs.append(histData[i][0]+halfBinWidth)
    plot.plot(xs,ys,color=lineColor,linewidth=3)
    plot.axes.set_ylim([0.0,1.02])
    if xLims!=False:
        plot.axes.set_xlim(xLims)
    plt.locator_params(axis='x',nbins=5) # maximum number of x labels
    plt.locator_params(axis='y',nbins=6) # maximum number of y labels
    plot.tick_params(axis='both',which='major',width=3,length=4,pad=6,direction='in',labelsize=25)
    plot.spines['right'].set_linewidth(2.0)
    plot.spines['bottom'].set_linewidth(2.0)
    plot.spines['top'].set_linewidth(2.0)
    plot.spines['left'].set_linewidth(2.0)
    # add axes label
    if latex:
        plot.axes.set_ylabel(r'$\displaystyle dp/dx\times(constant)$',fontsize=25)
    else:
        plot.axes.set_ylabel('dp/dx(*constant)',fontsize=25)
    plot.axes.set_xlabel(xLabel,fontsize=20)
    
    return plot
    
def histLoadAndPlot_ShadedPosteriors(plot,outFilename='',confLevels=False,xLabel='X',xLims=False,latex=False,showYlabel=False):
    """
    Loads previously plotted histograms that were written to disk by histPlotAndDump, and plot them up 
    in a way that is ready for publication.  This is the standard plotter used for plotting simple posteriors
    with shaded regions matching the 68% and 95% confidence.
    
    It is foreseen that many versions of this function will exist for different specific publication ready plots.
    """
    if outFilename[-4]!='.dat':
        outFilename=outFilename+'.dat'
    histData = np.loadtxt(outFilename)
    ys=[]
    xs=[]
    maxN = np.max(histData[:,1])
    minSub = 0
    if np.max(histData[:,0])>100000:
        #must be the To or Tc, so subtract int(min) and add to x-axis label
        #doing this as it doesn't go well allowing matplotlib to do it itself formatting wise
        minSub = int(np.min(histData[:,0]))
        histData[:,0]-=minSub
        xLabel = xLabel+" +"+str(minSub)
    halfBinWidth = (histData[1][0]-histData[0][0])/2.0
    # load up list of x,y values for tops of bins
    for i in range(0,histData.shape[0]):
        ys.append(histData[i][1]/maxN)
        ys.append(histData[i][1]/maxN)
        xs.append(histData[i][0]-halfBinWidth)
        xs.append(histData[i][0]+halfBinWidth)
    # load up list of shaded rectangle objects if confidence levels were provided
    recs = []
    if (type(confLevels)==list)or(type(confLevels)==np.ndarray):
        for i in range(0,histData.shape[0]):
            x=histData[i][0]-halfBinWidth+minSub
            c = 'w'
            if (x>confLevels[1][0])and(x<confLevels[1][1]):
                c = '0.8'
            if (x>confLevels[0][0])and(x<confLevels[0][1]):
                c = '0.5'
            recs.append(patches.Rectangle(xy=(histData[i][0]-halfBinWidth,0), width=halfBinWidth*2.0,height=histData[i][1]/maxN,facecolor=c, edgecolor=c))
        # draw updated patches on plot
        for rec in recs:
                plot.add_patch(rec)
            
    # draw the top line of hist
    plot.plot(xs,ys,color='k',linewidth=1)
    plot.axes.set_ylim([0.0,1.02])
    if xLims!=False:
        plot.axes.set_xlim(xLims)
    plt.locator_params(axis='x',nbins=3) # maximum number of x labels
    plt.locator_params(axis='y',nbins=3) # maximum number of y labels
    plot.tick_params(axis='both',which='major',width=1,length=2,pad=4,direction='in',labelsize=10)
    plot.spines['right'].set_linewidth(0.7)
    plot.spines['bottom'].set_linewidth(0.7)
    plot.spines['top'].set_linewidth(0.7)
    plot.spines['left'].set_linewidth(0.7)
    # add axes label
    if showYlabel:
        if latex:
            plot.axes.set_ylabel(r'$\frac{dp}{dx} \times constant $',fontsize=15)
        else:
            plot.axes.set_ylabel('dp/dx(*constant)',fontsize=15)
    else:
        plot.axes.set_yticklabels(['','',''])
    if latex:
        plot.axes.set_xlabel(r''+xLabel,fontsize=10)
    else:
        plot.axes.set_xlabel(xLabel,fontsize=10)
        
    return plot
    
    
def histConverter(chiSquareds, data, plot, xlabel, confLevels=False, weight=False, normed=False, nu=1, bestVal=0,logY=False,lineColor='k'):
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
            theWeights = genTools.likelihoodsCalc(chiSquareds, nu)
        else:
            theWeights = np.ones(len(data))
        # make initial histogram and update it below
        if verbose:
            print 'making initial hist, then convert it'
        if type(confLevels)==list:
            (n,bins,rectangles)=plot.hist(data, bins=50, normed=normed, weights=theWeights,log=logY, fill=False)
        else:
            if len(data)<20000:
                nbins = 50
            elif 50000<len(data)<200000:
                nbins=70
            else: 
                nbins=100
            #data =data/np.max(data)
            (n,bins,rectangles)=plot.hist(data, bins=nbins, normed=False, weights=theWeights,linewidth=7,color=lineColor, histtype = 'step',log=logY, fill=False)
        # Update rectangle patche's colors according to Confidence Levels of the data
        if verbose:
            print 'changing rectangle colors to match 68confLevels: '+repr(confLevels[0])+', and 95confLevels: '+repr(confLevels[1])
        recs2 = []
        if type(confLevels)==list:
            for rec in rectangles:
                    x=rec.get_x()
                    c = 'w'
                    cEdge = 'k'
                    if(x>confLevels[1][0])and(x<confLevels[1][1]):
                        c = '0.8'
                    if (x>confLevels[0][0])and(x<confLevels[0][1]):
                        c = '0.5'
                    recs2.append(patches.Rectangle(xy=rec.get_xy(), width=rec.get_width(),height=rec.get_height(),facecolor=c, edgecolor=cEdge))#color=c))
        # draw updated patches on plot
        for rec in recs2:
                plot.add_patch(rec)
        if False:        
            if type(confLevels)!=list:
                max = np.max(n)
                n=n/max
                #make a line plot based on centers of bins
                if type(bins)!=np.ndarray:
                    bins = np.array(bins)
                binsUse = (bins[1:]+bins[:-1])/2.0
                plot.plot(binsUse,n,color=lineColor,linewidth=4)
                plot.axes.set_ylim([0.0,1.1])
        
        plt.locator_params(axis='x',nbins=5) # maximum number of x labels
        #plt.locator_params(axis='y',nbins=6) # maximum number of y labels
        if True:
            #update the y limit and its ticks and tick labels
            if not logY:
                plot.axes.set_ylim([0.0,plot.axes.get_ylim()[1]*1.005])
                #print 'max = '+str(np.max(n))
                #print '1.2*max = '+str(np.max(n)*1.2)
                locs = []
                lim = np.max(n)*1.2
                for i in range(0,6):
                    locs.append((lim/6.0)*i)
                #print repr(plot.axes.get_yticks())
                #print 'locs = '+repr(locs)
                plot.axes.set_yticks(locs)
                plot.axes.set_yticklabels([' ','0.2','0.4','0.6','0.8','1.0'])
            if False:
                plot.axes.set_ylim([0.0,1.1])
                locs = [0.25,0.5,0.75,1.0]
                locLabels = ['0.25 ','0.5 ','0.75 ','1.0 ']
                plot.axes.set_yticks(locs)
                plot.axes.set_yticklabels(locLabels)

#         if normed:
#             #update the y limit and its ticks and tick labels
#             if not logY:
#                 plot.axes.set_ylim([0.0,plot.axes.get_ylim()[1]*1.3])
#                 locs = []
#                 lim = plot.axes.get_ylim()[1]
#                 for i in range(0,6):
#                     locs.append((lim/6.0)*i)
#                 #plot.axes.set_yticks(locs)
#                 #plot.axes.set_yticklabels([' ','0.2','0.4','0.6','0.8','1.0'])
        # draw a line up the median of the data
        Median = np.median(data)
        if type(confLevels)==list:
            plot.plot([Median, Median],plot.axes.get_ylim(),'k',linewidth=4)
        #convert data and chiSquareds arrays to np arrays if needed
        if type(chiSquareds)!=np.ndarray:
            chiSquareds = np.array(chiSquareds)
        if type(data)!=np.ndarray:
            data = np.array(data)
        if bestVal==0:
            try:
                bestVal = data[np.where(chiSquareds==chiSquareds.min())][0]
                if debug:
                    print "using calculated bestVal = ",bestVal
            except:
                useless=True
        else:
            if debug:
                print 'using provided bestVal = ',bestVal
        if type(confLevels)==list:
            plot.plot([bestVal, bestVal],plot.axes.get_ylim(),'g',linewidth=4)
        
        if verbose:
            print "min found to be "+str(np.min(data,axis=0))+", max was "+str(np.max(data,axis=0))
            print "Best value was found to be at "+str(bestVal)+", and median at "+str(Median)
    else:
        if debug:
            print 'data vals all equal so making empty plot for kicks'
        if (type(data)!=float):
            data = data[0]
        # This means the data had no useful values in it, just constant 
        plot.plot([data,data],[0.0,1.3],'b',linewidth=2)
        if debug:
            print 'empty plot made'
        plot.axes.set_ylim([0.0,1.3])
        if debug:
            print 'axes limits set'
        if False:
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
            #plot.axes.set_yticklabels(['0.0','0.2','0.4','0.6','0.8','1.0'])
        if verbose:
            print 'Values for '+xlabel+' were found to be constant, so not making a histogram, just a simple line plot!!'
    # add x label
    plot.axes.set_xlabel(xlabel,fontsize=8)
    plot.tick_params(axis='both',which='major',width=1,length=3,pad=5,direction='in',labelsize=5)
    #plot.axhline(linewidth=2.0)
    #plot.axvline(linewidth=2.0)
    #plot.xaxis.set_tick_params(width=5,length=5,pad=2,direction='in')
    #plot.yaxis.set_tick_params(width=5,length=5,pad=2,direction='in')
    plot.spines['right'].set_linewidth(1.0)
    plot.spines['bottom'].set_linewidth(1.0)
    plot.spines['top'].set_linewidth(1.0)
    plot.spines['left'].set_linewidth(1.0)
    plot = fixPlotBordersAndLabels(plot)
    
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
    return Polygon(polystar, fc=color, ec='black',linewidth=1.5)    

def addRVdataToPlot(subPlot,RVs,RVerrors,epochs,alf=1.0,color='blue',plotErrorBars=False):
    """
    Add '+' markers for the data locations with respective y axis errors 
    shown as the height of the markers. 
    """
    for epoch in range(0,len(epochs)):
        try:
            xs = [epochs[epoch],epochs[epoch]]
            ys = [RVs[epoch]-abs(RVerrors[epoch]),RVs[epoch]+abs(RVerrors[epoch])]
        except:
            print '\n\nERROR trying to calc x and y for RV data epoch #'+str(epoch)
            if len(epochs)>1:
                print 'epochs[epoch] = '+str(epochs[epoch])
            else:
                print 'len(epochs)<1!!!'
            if len(RVs)>1:
                print 'RVs[epoch] = '+str(RVs[epoch])
            else:
                print 'len(RVs)<1!!!'
            if len(RVerrors)>1:
                print 'RVerrors[epoch] = '+str(RVerrors[epoch])
            else:
                print 'len(RVerrors)<1!!!'
        if plotErrorBars:
            subPlot.plot(xs,ys,c=color,linewidth=2,alpha=alf)
        subPlot.plot(epochs[epoch],RVs[epoch],c='k',marker='.',markersize=4)
    return subPlot

def addDIdataToPlot(subPlot,SAs,SAerrors,PAs,PAerrors,asConversion,telescopeView=False):
    """
    To plot a '+' for each data point with width and height matching the errors converted 
    to x,y coords.
    This is the new version to use for this over the older startAndErrorPolys function.
    """
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
            
        left = xCent-0.5*w
        right = xCent+0.5*w
        top = yCent+0.5*h
        btm = yCent-0.5*h
        
        subPlot.plot([left,right],[yCent,yCent],linewidth=3,color='k',alpha=1.0)
        subPlot.plot([xCent,xCent],[btm,top],linewidth=3,color='k',alpha=1.0)
        # see if new min max found for x,y
        if xmax<right:
            xmax = right
        if xmin>left:
            xmin = left
        if ymax<top:
            ymax = top
        if ymin>btm:
            ymin = btm
            
    return (subPlot,[xmin,xmax,ymin,ymax])
      
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
        if False:
            ### Older version that makes shaded boxes for each
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
            
        rect = patches.Rectangle((xCorner,yCorner),width=w,height=h,facecolor='black',edgecolor='black',alpha=1.0,linewidth=1.5)
        t = pylab.matplotlib.transforms.Affine2D().rotate_deg_around(xCent,yCent,-PAs[i]) +transData
        rect.set_transform(t)
        errorBoxes.append(rect)
        
        # determin x and y locations of the observed PA and SA's for companion star/planet
        # then make a star polygon for each, same as for M1 but much smaller
        m2starPolygons.append(star((asConversion/1000.0)*12.0*SAs[0], xCent, yCent, color='red', N=5, thin = 0.5))
        
    return (errorBoxes, m2starPolygons,[xmin,xmax,ymin,ymax])
    
def makeCleanSummaryPlot(outputDataFilename=''):
    """
    This will first find the lowest chisquared for a given run,
    then use readTrimWrite to trim all data that is more than 40 higher than that value.
    This trimmed data file will then be plotted by summaryPlotter.
    """
    quiet = True
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
        newDataFilename = genTools.dataReadTrimWrite(outputDataFilename, maxNonReducedChiSquared, verbose=False)
        if verbose:
            print '\nTrimmed data file at: '+newDataFilename
        
        ## make new plotfile name
        newPlotFilename = os.path.join(datadir,'SummaryPlot-reducedChiSquaredMax_'+str(int(maxReducedChiSquared))+'.png')
        if verbose:
            print '\nNew plot for trimmed data will be written to: '+newPlotFilename
        
        ## pass trimmed data to the plotter
        if verbose:
            print '\nCalling summaryPlotter to plot new trimmed data'
        summaryPlotter(newDataFilename, newPlotFilename, weight=False, confLevels=True, normalize=True, nu=1, plot4x1=False, TcStepping=False)
        if quiet==False:
            print '\nNew plot for trimmed data written to: '+newPlotFilename
            print '\nTrimmed data all plotted up baby!! :-D\n'
    else:
        if quiet==False:
            print "Filename provided does not exist!!  Provided path was: "+outputDataFilename
        
    if True:
        return newDataFilename
            
def densityPlot(dataDir,sysDataDict, DIdataDict,RVdataDict,paramSettingsDict,nuDI=1,nuRV=1):
    """
    """
    alf = 0.03
    N = 1000
    DI = False
    print 'Starting to make a density plot for the data fit'
    # record the time the chain started
    startTime = timeit.default_timer()
    
    outputDataFile = os.path.join(dataDir,'outputData-ALL.dat')
    (chiSquareds, incs, es, longANs, periods, argPeris, semiMajors, Ts, Tcs, Ks, rvOffsets) = genTools.getEveryNthOrbElements(outputDataFile,N=N)
    print 'getEveryNthOrbElements found '+str(len(incs))+" sets"
    print 'np.median(longANs) = '+str(np.median(longANs))
    print 'np.median(es) = '+str(np.median(es))
    print 'np.median(periods) = '+str(np.median(periods))
    print 'np.median(incs) = '+str(np.median(incs))
    print 'np.median(argPeris) = '+str(np.median(argPeris))
    print 'np.median(Ts) = '+str(np.median(Ts))
    print 'np.median(Tcs) = '+str(np.median(Tcs))
    print 'np.median(semiMajors) = '+str(np.median(semiMajors))
    print 'np.median(Ks) = '+str(np.median(Ks))
    print 'np.median(rvOffsets) = '+str(np.median(rvOffsets))
    if True:
        #update argPeri value to take offset into account
        argPerisUse = []
        for i in range(0,len(argPeris)):
            argPerisUse.append(argPeris[i]+paramSettingsDict['argPeriOffsetDI'])
        orbitEllipsePlotFilename = os.path.join('/mnt/Data1/Todai_Work/Dropbox/SMODT-outputCopies/','DIdensityPlot')
        orbitEllipsePlotter(longANs, es, periods, incs, argPerisUse, semiMajors, Ts,sysDataDict, DIdataDict,\
                              plotFilename=orbitEllipsePlotFilename, xLim=False, yLim=False, show=False,nuDI=1,alf=alf)
    if True:
        #update argPeri value to take offset into account
        argPerisUse = []
        for i in range(0,len(argPeris)):
            argPerisUse.append(argPeris[i]+paramSettingsDict['argPeriOffsetRV'])
        orbitrvPlotFilename = os.path.join('/mnt/Data1/Todai_Work/Dropbox/SMODT-outputCopies/','RVdensityPlot')
        rvPlotter(es, Ts, Tcs, periods, incs, argPerisUse, semiMajors, sysDataDict, RVdataDict, paramSettingsDict,\
                 K=Ks, RVoffsets=rvOffsets, nuRV=1, plotFilename=orbitrvPlotFilename, show=False, plotFullOrbit=True, primaryRVs=True,alf=alf)
    print 'Finished making the density plot for the data fit'
    endTime = timeit.default_timer()
    totalTime = (endTime-startTime) # in seconds
    totalTimeString = genTools.timeString(totalTime)
    print 'That took '+totalTimeString+' to complete.\n'

def stackedPosteriorsPlotterHackStarter():
    outputDataFilenames = []
    if False:
        outputDataFilenames.append('/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-tight-20PercentRealizedError--14-Million-in_Total/outputData-ALL.dat')
        outputDataFilenames.append('/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-tight-10PercentRealizedError--14-Million-in_Total/outputData-ALL.dat')
        outputDataFilenames.append('/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-tight-5PercentRealizedError--14-Million-in_Total/outputData-ALL.dat')
        outputDataFilenames.append('/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-tight-1PercentRealizedError--14-Million-in_Total/outputData-ALL.dat')
    else:
        outputDataFilenames.append('/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-veryOpen-20PercentRealizedError2--50-Million-in_Total/outputData-ALL.dat')
        outputDataFilenames.append('/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-veryOpen-10PercentRealizedError2--50-Million-in_Total/outputData-ALL.dat')
        outputDataFilenames.append('/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-veryOpen-5PercentRealizedError2--50-Million-in_Total/outputData-ALL.dat')
        outputDataFilenames.append('/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-veryOpen-1PercentRealizedError2--50-Million-in_Total/outputData-ALL.dat')
    
    plotFilename = os.path.join(os.path.abspath('/mnt/Data1/Todai_Work/Dropbox/SMODT-outputCopies/'),'stackedPosteriorsPlot')
    stackedPosteriorsPlotterFunc(outputDataFilenames, plotFilename,ALLparams=False)
    #print 'Final stacked plot file written to:\n'+plotFilename
    if True:
        print 'converted to PDF as well'
        os.system("epstopdf "+plotFilename+'.eps')

def stackedPosteriorsPlotterFunc(outputDataFilenames, plotFilename,ALLparams=True):
    """
    will plot a simple posterior distribution for each parameter in the data files
    stacked ontop of eachother for comparison between different runs.
    """
    weight=False
    confLevels=False
    normalize=False
    nu=1
    quiet = True
    verbose = False
    plotFormat = 'eps'
    latex=False
    if ALLparams==False:
        latex=True
    
    #plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rc('font',family='serif')
    if latex:
        plt.rc('text', usetex=True)
        plt.rcParams['text.latex.unicode']=True 
    
    if type(outputDataFilenames)!=list:
        outputDataFilenames = [outputDataFilenames]
    
    colorsList =['Black','Blue','DarkGreen','Red','Purple','Fuchsia','Crimson','Aqua','Gold','OrangeRed','Plum','Chartreuse','Chocolate','SteelBlue ','Teal','Salmon','Brown']
    
    if os.path.exists(outputDataFilenames[0]):  
        s= '\nCreating a simple plot of some key posteriors for files:\n'
        for filename in outputDataFilenames:
            s+=filename+'\n'
        s=s+ '\nInput plotfilename:\n'+plotFilename
        if verbose:
            print s
        # record the time the chain started
        startTime = timeit.default_timer()
        
        ## find number of RV datasets
        f = open(outputDataFilenames[0],'r')
        plotFileTitle = f.readline()[:-5]
        headings = f.readline()
        line = f.readline()
        dataLineCols = line.split()
        numRVdatasets=0
        if (len(line)>10):
            numRVdatasets = len(dataLineCols) - 10
        else:
            line = f.readline()
            dataLineCols = line.split()
            if (len(line)>10):
                numRVdatasets = len(dataLineCols) - 10
        s= "\nNumber of RV datasets found in summaryPlotter was "+str(numRVdatasets)+"\n"
        if verbose:
            print s
        f.close()
        
        # check if the passed in value for plotFilename includes format extension
        if '.'+plotFormat not in plotFilename:
            plotFilename = plotFilename+"."+plotFormat
            print 'updating plotFilename to : '+plotFilename
        else:
            plotFilename = plotFilename
            print 'plotfilename was found to already have the format extension'
        
        s='\nStarting to make Simple Key Posteriors Plot'
        if verbose:
            print s        
        
    if ALLparams:
        fig = plt.figure(1,figsize=(70,55)) 
        #Load up parameter column number and name string arrays
        paramList = [0,1,2,3,4,5,6,7]
        paramStrs = ['Omega [deg]','e','To [JD]', 'Tc [JD]','P [Yrs]','i [deg]','omega [deg]','a_total [AU]']
        perfectVals = [70.0,0.5,2457000.0,2457000.0,5.0,40.0,50.0,3.34718746581]
        if numRVdatasets>0:
            paramList.append(9)
            paramStrs.append('K [m/s]')
            perfectVals.append(4933.18419284)
            for dataset in range(1,numRVdatasets+1):
                paramList.append(9+dataset)
                paramStrs.append('RV offset '+str(dataset)+' [m/s]')
                perfectVals.append(0.0)
    else:
        fig = plt.figure(1,figsize=(9.625,4.375)) 
        positions = [[0.12,0.15,0.36,0.8],[0.60,0.15,0.36,0.8]]
        paramList = [1,4]
        paramStrs = ['Eccentricity','Period [Yrs]']
        perfectVals = [0.50,5.00]
        xLims2 = [[0.3,0.59],[4.51,5.99]]
    
    ## run through all the data files and parameters requested and make root histogram files
    for i in range(0,len(paramList)):
        for outputDataFilename in outputDataFilenames:
            if os.path.exists(os.path.join(os.path.dirname(outputDataFilename),'hist-'+str(i)+'_n.pkl'))==False:
                print 'Initial Plotting for parameter '+str(i)+"/"+str(len(paramList)-1)+" for file:\n"+outputDataFilename
                histDataBaseName = os.path.join(os.path.dirname(outputDataFilename),'hist-'+str(i))
                #(CLevels,data,bestDataVal) = genTools.confLevelFinder(outputDataFilename,paramList[i], returnData=True, returnChiSquareds=False, returnBestDataVal=True,fast=False)
                (log,data,chiSquareds,[bestDataVal,dataMedian,dataValueStart,dataValueMid,dataValueEnd]) = genTools.dataReader(outputDataFilename, paramList[i], returnData=True, returnChiSquareds=False, returnBestDataVal=False, ignoreConstParam=False)
                histMakeAndDump([],data,outFilename=histDataBaseName,weight=False, normed=False, nu=1,logY=False,histType='step')
                
    ## make combined/stacked plot for each parameter in list
    for i in range(0,len(paramList)):
        startTime = timeit.default_timer()
        s='\nStarting to plot combined hist for '+paramStrs[i]
        if verbose:
            print s
        #log.write(s+'\n')
        if ALLparams:
            subPlot = fig.add_subplot(3,4,i+1)
        else:
            subPlot = fig.add_subplot(1,2,i+1)
        colorInt = 0
        for outputDataFilename in outputDataFilenames:
            if os.path.exists(outputDataFilename):
                histDataBaseName = os.path.join(os.path.dirname(outputDataFilename),'hist-'+str(i))
                print 'Loading and re-plotting parameter '+str(i)+"/"+str(len(paramList)-1)+" for file:\n"+outputDataFilename
                xLim=False
                if ALLparams==False:
                    xLim = xLims2[i]  
                subPlot = histLoadAndPlot_StackedPosteriors(subPlot,outFilename=histDataBaseName,xLabel=paramStrs[i],lineColor=colorsList[colorInt],xLims=xLim,latex=latex)
                #subPlot = histConverter([], data, subPlot, paramStrs[i], confLevels=False, weight=weight, normed=normalize, nu=nu, bestVal=perfectVals[i],lineColor=colorsList[colorInt])
                #print 'color = '+colorsList[colorInt]
                #subPlot.axes.set_ylabel('dp/dx (*constant)',fontsize=55)
                if ALLparams==False:
                    subPlot.set_position(positions[i])
                colorInt+=1
            else:
                s= "stackedPosteriorsPlotterFunc: ERROR!!!! file doesn't exist"
                print s
                break
        # record the time the chain finished and print
        endTime = timeit.default_timer()
        totalTime = (endTime-startTime) # in seconds
        totalTimeString = genTools.timeString(totalTime)
        s=s+'That took '+totalTimeString+' to complete.\n'  ##### only print in !'silent' mode to track time
        if quiet==False:
            print s
        #log.write(s+'\n')
        
        
    ## Save file if requested.
    if verbose:
        print '\nStarting to save param hist figure:'
    if plotFilename!='':
        plt.savefig(plotFilename,format=plotFormat)
        s= 'Summary plot saved to: '+plotFilename
        print s
        if quiet==False:
            print s
        #log.write(s+'\n')
    plt.close()
    
    ## record the time the chain finished and print
    endTime = timeit.default_timer()
    totalTime = (endTime-startTime) # in seconds
    totalTimeString = genTools.timeString(totalTime)
    s= '\n\stackedPosteriorsPlotterFunc: Plotting took '+totalTimeString+' to complete.\n'
    if verbose:
        print s
    
def summaryPlotter(outputDataFilename, plotFilename, shadeConfLevels=True, nu=1):
    """
    This advanced plotting function will plot all the data in a 3x3 grid (or 3x1) on a single figure.  The data will be plotted
    in histograms that can be will be normalized to a max of 1.0 and the data can be weighted if desired.  The 
    confidence levels of the data can be calculated and the resulting histograms will have the bars inside the
    68% bin shaded as dark grey and 95% as light grey and the rest will be white.
    
    :param nu:       The value of nu to convert chi squared to reduced chi squared.
    :type nu:        float32
    """
    quiet = True
    verbose = False
    latex=True
    plotFormat = 'eps'
    forceRecalc = False
    
    #plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rc('font',family='serif')
    if latex:
        plt.rc('text', usetex=True)
        plt.rcParams['text.latex.unicode']=True  
        
    ## find number of RV datasets
    if os.path.exists(outputDataFilename):  
        datadir = os.path.dirname(outputDataFilename)
        
        ## get log
        logFilename = os.path.join(datadir,'processManagerLogFile.txt')
        log = open(logFilename,'a')
        log.write('\n'+75*'#'+'\n Inside summaryPlotter \n'+75*'#'+'\n')
         
        s= '\nCreating summary plot for file:\n'+outputDataFilename
        s=s+ '\nInput plotfilename:\n'+plotFilename
        if verbose:
            print s
        log.write(s+'\n')
        
        ## find number of RV datasets
        f = open(outputDataFilename,'r')
        plotFileTitle = f.readline()[:-5]
        headings = f.readline()
        line = f.readline()
        dataLineCols = line.split()
        numRVdatasets=0
        if (len(line)>10):
            numRVdatasets = len(dataLineCols) - 10
        else:
            line = f.readline()
            dataLineCols = line.split()
            if (len(line)>10):
                numRVdatasets = len(dataLineCols) - 10
        s= "\nNumber of RV datasets found in summaryPlotter was "+str(numRVdatasets)+"\n"
        if verbose:
            print s
        log.write(s+'\n')
        f.close()
           
        # check if the passed in value for plotFilename includes format extension
        if '.'+plotFormat not in plotFilename:
            plotFilename = plotFilename+"."+plotFormat
            print 'updating plotFilename to:\n'+plotFilename
        else:
            plotFilename = plotFilename
            print 'plotfilename was found to already have the format extension'
            
        paramList = [0,1,2,3,4,5,6,7]
        paramFileStrs = ['Omega','e','To', 'Tc','P','i','omega','a_total']
        if latex:
            paramStrs = ['$\Omega$ [deg]','$e$','$T_o$ [JD]', '$T_c$ [JD]','$P$ [Yrs]','$i$ [deg]','$\omega$ [deg]','$a_{total}$ [AU]']
        else:
            paramStrs = ['Omega [deg]','e','To [JD]', 'Tc [JD]','P [Yrs]','i [deg]','omega [deg]','a_total [AU]']
        #perfectVals = [70.0,0.5,2457000.0,2457000.0,5.0,40.0,50.0,3.34718746581]
        if numRVdatasets>0:
            paramList.append(9)
            paramFileStrs.append('K')
            if latex:
                paramStrs.append('$K$ [m/s]')
            else:
                paramStrs.append('K [m/s]')
            #perfectVals.append(4933.18419284)
            for dataset in range(1,numRVdatasets+1):
                paramList.append(9+dataset)
                paramFileStrs.append('offset_'+str(dataset))
                if latex:
                    #paramStrs.append("$\alpha$ "+str(dataset)+" [m/s]")
                    paramStrs.append("$offset_{"+str(dataset)+"}$ [m/s]")
                else:
                    paramStrs.append('offset '+str(dataset)+' [m/s]')
                #perfectVals.append(0.0)
            
        ## run through all the data files and parameters requested and make histogram files
        for i in range(0,len(paramList)):
            if (os.path.exists(os.path.join(os.path.dirname(outputDataFilename),'hist-'+str(i)+'.dat'))==False)or forceRecalc:
                if verbose:
                    print 'Initial Plotting for parameter '+str(i+1)+"/"+str(len(paramList))+": "+paramStrs[i]+", for file:\n"+outputDataFilename
                histDataBaseName = os.path.join(os.path.dirname(outputDataFilename),'hist-'+str(i))
                (CLevels,data,bestDataVal) = genTools.confLevelFinder(outputDataFilename,paramList[i], returnData=True, returnChiSquareds=False, returnBestDataVal=True,fast=False)
                #(log,data,chiSquareds,[bestDataVal,dataMedian,dataValueStart,dataValueMid,dataValueEnd]) = genTools.dataReader(outputDataFilename, paramList[i], returnData=True, returnChiSquareds=False, returnBestDataVal=False, ignoreConstParam=False)
                histMakeAndDump([],data,outFilename=histDataBaseName,weight=False, normed=False, nu=1,logY=False,histType='step')
                if (os.path.exists(os.path.join(os.path.dirname(outputDataFilename),'confLevels-'+str(i)+'.dat'))==False)or forceRecalc:
                    np.savetxt(os.path.join(os.path.dirname(outputDataFilename),'confLevels-'+str(i)+'.dat'),CLevels)
                    if verbose:
                        print 'confidence levels data stored to:\n'+os.path.join(os.path.dirname(outputDataFilename),'confLevels-'+str(i)+'.dat')
        
        # Create empty figure to be filled up with plots
        fig = plt.figure(figsize=(10,10)) 
                    
        ## make combined/stacked plot for each parameter in list
        for i in range(0,len(paramList)):
            startTime = timeit.default_timer()
            s='\nStarting to plot shaded hist for '+paramStrs[i]
            if verbose:
                print s
            subPlot = fig.add_subplot(3,4,i+1)
            
            histDataBaseName = os.path.join(os.path.dirname(outputDataFilename),'hist-'+str(i))
            print 'Loading and re-plotting parameter '+str(i+1)+"/"+str(len(paramList))+": "+paramStrs[i]#+" for file:\n"+outputDataFilename
            xLim=False
            CLevels=False
            if shadeConfLevels:
                CLevels=np.loadtxt(os.path.join(os.path.dirname(outputDataFilename),'confLevels-'+str(i)+'.dat'))
            showYlabel=False
            if (i==0)or(i==4)or(i==8):
                showYlabel = True
            subPlot = histLoadAndPlot_ShadedPosteriors(subPlot,outFilename=histDataBaseName,confLevels=CLevels,xLabel=paramStrs[i],xLims=xLim,latex=latex,showYlabel=showYlabel)
            
            # record the time the chain finished and print
            endTime = timeit.default_timer()
            totalTimeString = genTools.timeString(endTime-startTime)
            s=s+'That took '+totalTimeString+' to complete.\n'  ##### only print in !'silent' mode to track time
            if quiet==False:
                print s            
            
        ## Save file if requested.
        if verbose:
            print '\nStarting to save param hist figure:'
        if plotFilename!='':
            plt.savefig(plotFilename,format=plotFormat)
            s= 'Summary plot saved to: '+plotFilename
            print s
            if quiet==False:
                print s
            #log.write(s+'\n')
        plt.close()
        if True:
            print 'converting to PDF as well'
            os.system("epstopdf "+plotFilename)


        
def progessPlotterSingleFile(outputDataFilename, plotFilename, nu=1, plot4x1=False, logFilename='',TcStepping=False):
    """
    NOTE: For a single chain's output data.
    
    This will plot a progress plot for each varying parameters in the file.  
    It is called by mcmcProgressPlotter to plot up the progress for each chain
    separately and utilizes progessPlotterSingleFileFunc to do the repetative task 
    of producing the plot for each parameter.
    
    These plots will be made for the parameter values vs step, 
    to show how it is converging with time.
    
    :param nu:       The value of nu to convert chi squared to reduced chi squared.
    :type nu:        float
    :param plot4x1:  Only make plots for K, argPeri, e, To (or Tc if Tc stepping was used), and chiSquared.
    :type plot4x1:   Python boolean
    """
    quiet = True
    verbose = False
    #find chain number and update logFilename with path and number
    s = outputDataFilename
    chainNumStr = s[s.find('chain_')+6:s.find('chain_')+6+1]
    datadir = os.path.dirname(outputDataFilename)
    logFilename = os.path.join(datadir,'log-chain_'+chainNumStr+'.txt')
    log = open(logFilename,'a')
    log.write('\n'+75*'*'+'\n Inside progessPlotterSingleFile \n'+75*'*'+'\n')
    
    if os.path.exists(outputDataFilename):   
        s=  '\nCreating summary plot for file:\n'+outputDataFilename
        s=s+'\nInput plotfilename:\n'+plotFilename
        if quiet==False:
            print s
        log.write(s+'\n')
        # record the time the chain started
        startTime = timeit.default_timer()
        
        ##*******************************
        ## find number of RV datasets
        ##*******************************
        f = open(outputDataFilename,'r')
        plotFileTitle = f.readline()[:-5]
        headings = f.readline()
        line = f.readline()
        dataLineCols = line.split()
        rvFirstVal = 0
        if (len(line)>10):
            numRVdatasets = len(dataLineCols) - 10
            rvFirstVal = float(dataLineCols[-2])
        else:
            line = f.readline()
            dataLineCols = line.split()
            if (len(line)>10):
                numRVdatasets = len(dataLineCols) - 10
        s= "\nNumber of RV datasets found in progessPlotterSingleFile was "+str(numRVdatasets)+"\n"
        if verbose:
            print s
        log.write(s+'\n')
        f.close()
           
        ##*************************************************************
        # check if the passed in value for plotFilename includes '.png'
        ##*************************************************************
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
        ##*************************************************************
        ## make an advanced title for plot from folder and filename
        ##*************************************************************
        titleTop = os.path.dirname(outputDataFilename).split('/')[-1]
        titleBtm = os.path.basename(plotFilename).split('.')[0]+" Parameter Progress Plot"
        plotFileTitle = titleTop+'\n'+titleBtm
        
        s='Starting to make Total Summary Plot'
        log.write(s+'\n')
        if verbose:
            print s 

        ##*************************************************************
        # Create empty figure to be filled up with plots
        ##*************************************************************
        if not plot4x1:
            fig2 = plt.figure(2, figsize=(30,55),dpi=200)
            s= '\n** a 7x1 figure will be made for all 6 orbital parameters in the output datafile **'
            if verbose:
                print s
            log.write(s+'\n')
        else:
            fig2 = plt.figure(2,figsize=(30,35),dpi=200)
            s= '\n** a 4x1 figure will be made for all 6 orbital parameters in the output datafile **'
            if verbose:
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
        #(CLevels,data,chiSquareds,bestDataVal) = genTools.confLevelFinder(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True,fast=False)
        (log,data,chiSquareds,[bestDataVal,dataMedian,dataValueStart,dataValueMid,dataValueEnd]) = genTools.dataReader(outputDataFilename, paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True, ignoreConstParam=True)
        #log.write('\nCLevels for '+xlabel+':\n'+repr(CLevels)+'\n')
        xs = np.array(range(0,data.size))
        # Wrap the array into a 2D array of chunks, truncating the last chunk if 
        if data.size>90e6:
            chunkSize = 10000
        elif data.size>9e6:
            chunkSize = 1000
        elif data.size>9e5:
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
        try:
            subPlot2.fill_between(xCenters, min_env, max_env, color='gray', edgecolor='none', alpha=0.5)
            subPlot2.plot(xCenters, yCenters)
            subPlot2.plot([xCenters.min(),xCenters.max()],[bestDataVal,bestDataVal],color='green')
            #subPlot2.plot(range(0,data.size),data)
            subPlot2.axes.set_ylabel(xlabel)
            #argPeriMedian = np.median(argPeri_degsAlls)
            s= "done plotting argPeri_degsAlls"
            if quiet==False:
                print s
            log.write(s+'\n')
        except:
            log.write("An error occured while trying to plot summary of argPeri_degs.\n This is most likely due to it being constant.")
        
        if not plot4x1:
            subPlot2 = fig2.add_subplot(815)
            paramColNum = 5
            xlabel = 'Inclination [deg]'
            (log,subPlot2,data,bestDataVal)=progessPlotterSingleFileFunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize,ignoreConstParam=False)
            ##################$$$$$$$$$$$$$ This extra garbage collection might not be needed but I want it for now as a code EX. ######
            #del inclination_degsAlls
            #gc.collect()
            s= "done plotting inclination_degsAlls"
            if quiet==False:
                print s
            log.write(s+'\n')
        
        
        
        # create sub plot and fill it up for the longitude of ascending node
        #longANMedian=0
        if not plot4x1:
            subPlot2 = fig2.add_subplot(811)
            paramColNum = 0
            xlabel = 'Longitude of Ascending Node [deg]'
            (log,subPlot2,data,bestDataVal)=progessPlotterSingleFileFunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize,ignoreConstParam=False)
            #longANMedian = np.median(data)
            s= "done plotting longAN_degsAlls"
            if quiet==False:
                print s
            log.write(s+'\n')
        
        # Create sub plot and fill it up for the e
        if not plot4x1:
            subPlot2 = fig2.add_subplot(812)
        else:
            subPlot2 = fig2.add_subplot(511)
        paramColNum = 1
        xlabel = 'e'
        (log,subPlot2,data,beste)=progessPlotterSingleFileFunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize,ignoreConstParam=False)
        #eMedian = np.median(esAlls)
        s= "done plotting esAlls"
        if quiet==False:
            print s
        log.write(s+'\n')
            
        # Create sub plot and fill it up for the Period
        if not plot4x1: 
            subPlot2 = fig2.add_subplot(814)
            paramColNum = 4
            xlabel = 'Period [Years]'
            (log,subPlot2,data,bestDataVal)=progessPlotterSingleFileFunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize,ignoreConstParam=False)
            #periodMedian = np.median(periodsAlls)
            s= "done plotting periodsAlls"
            if quiet==False:
                print s
            log.write(s+'\n')
            
        # Create sub plot and fill it up for the semi-majors
        semiMajorBest = 0.0
        if not plot4x1: 
            subPlot2 = fig2.add_subplot(817)
            paramColNum = 7
            xlabel = 'Semi-Majors [AU]'
            (log,subPlot2,data,bestDataVal)=progessPlotterSingleFileFunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize,ignoreConstParam=False)
            semiMajorBest = bestDataVal
            #print "\n\n"+"%"*100+"semiMajorBest = "+str(semiMajorBest)+"%"*100+"\n\n"
            s= "done plotting semi-majors"
            if quiet==False:
                print s
            log.write(s+'\n')
            
        # Create sub plot and fill it up for the Semi-Major Amplitude
        # put it over the longAN plot if it is empty
        if not plot4x1: 
            subPlot2 = fig2.add_subplot(817)
        else:
            subPlot2 = fig2.add_subplot(513)
        if semiMajorBest==0.0:
            paramColNum = 9
            xlabel = 'K [m/s]'
            (log,subPlot2,data,bestDataVal)=progessPlotterSingleFileFunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize,ignoreConstParam=False)
            s= "done plotting Ks"
            if quiet==False:
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
            (log,subPlot2,data,bestDataVal)=progessPlotterSingleFileFunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize,ignoreConstParam=False)
            #TMedian = np.median(TsAlls)
            s= "done plotting TsAlls"
            if quiet==False:
                print s
            log.write(s+'\n')
                
        ## NOTE: this chiSquareds plot is special, so it can't use the function to do the work
        try:
            if not plot4x1:        
                subPlot2 = fig2.add_subplot(818)
            else:
                subPlot2 = fig2.add_subplot(515)
            paramColNum = 8
            xlabel = 'ChiSquareds'
            (log,data,chiSquareds,[bestDataVal,dataMedian,dataValueStart,dataValueMid,dataValueEnd]) = genTools.dataReader(outputDataFilename, paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True,ignoreConstParam=False)
            #(CLevels,data) =genTools.confLevelFinder(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=False)
            #log.write('\nCLevels for '+xlabel+':\n'+repr(CLevels)+'\n')
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
            if quiet==False:
                print s
            log.write(s+'\n')                   
            
            # Save file if requested.
            if plotFilename!='':
                plt.savefig(plotFilename, orientation='landscape')
                s= '\nSummary plot saved to:\n'+plotFilename
                if quiet==False:
                    print s
                log.write(s+'\n')
            plt.close()
        except:
            print "Unable to produce a Chi Squared summary distribution plot for some reason, probably not enough data."
            print "Just going to save current version of figure without the chi squared progress plot."
            # Save file if requested.
            if plotFilename!='':
                plt.savefig(plotFilename, orientation='landscape')
                s= '\nSummary plot saved to:\n'+plotFilename
                if quiet==False:
                    print s
                log.write(s+'\n')
            plt.close()
        
        if False:
            ## separate figure for the chiSquareds ##
            fig3 = plt.figure(1,figsize=(30,10),dpi=200)
            subPlot2 = fig3.add_subplot(211)
            paramColNum = 8
            xlabel = 'ChiSquareds'
            (log,data,chiSquareds,[bestDataVal,dataMedian,dataValueStart,dataValueMid,dataValueEnd]) = genTools.dataReader(outputDataFilename, paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True,ignoreConstParam=True)
            #(CLevels,data) =genTools.confLevelFinder(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=False)
            #log.write('\nCLevels for '+xlabel+':\n'+repr(CLevels)+'\n')
            yChunks = data[:chunkSize*numChunks].reshape((-1,chunkSize))
            max_env = yChunks.max(axis=1)
            min_env = yChunks.min(axis=1)
            yCenters = yChunks.mean(axis=1)
            subPlot2.fill_between(xCenters, min_env, max_env, color='gray',edgecolor='none', alpha=0.5)
            subPlot2.plot(xCenters, yCenters)
            #subPlot2.plot(range(0,data.size),data)
            subPlot2.axes.set_ylabel(xlabel)
            s= "done plotting ChiSquareds full Ranges"
            if quiet==False:
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
            if quiet==False:
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
                if quiet==False:
                    print s
                log.write(s+'\n')
            plt.close()
        
        
        ## Create a second figure of RV offsets. ####
        if ((numRVdatasets>0)and(rvFirstVal!=0)):
            try:
                s = "\nTrying to produce posterior histograms of the RV offsets"
                s+="/n ** First RV value was "+str(rvFirstVal)
                if verbose:
                    print s
                log.write(s+'\n')
                # Create empty figure to be filled up with plots
                # Create sub plot and fill it up for the Semi-major
                if numRVdatasets==1:
                    fig2 = plt.figure(2,figsize=(30,10),dpi=200)
                elif numRVdatasets==2:
                    fig2 = plt.figure(2,figsize=(30,15),dpi=200)
                elif numRVdatasets==3:
                    fig2 = plt.figure(2,figsize=(30,20),dpi=200)
                else:
                    s= '\nprogessPlotterSingleFile: WARNING!!! Plotter not set up to handle more than 3 RV datasets and '\
                    +str(numRVdatasets)+' were found in the output datafile.'
                    if quiet==False:
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
                    (log,subPlot2,data,bestDataVal)=progessPlotterSingleFileFunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize,ignoreConstParam=False)
                    s= "Done plotting RV offsets for dataset 1"
                    if quiet==False:
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
                    (log,subPlot2,data,bestDataVal)=progessPlotterSingleFileFunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize,ignoreConstParam=False)
                    s= "Done plotting RV offsets for dataset 2"
                    if quiet==False:
                        print s
                    log.write(s+'\n')
                    
                if numRVdatasets==3:
                    # Create sub plot and fill it up for the Semi-major
                    subPlot2 = fig2.add_subplot(313)
                    paramColNum = 12
                    xlabel = 'RV offset 3 [m/s]'
                    (log,subPlot2,data,bestDataVal)=progessPlotterSingleFileFunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize,ignoreConstParam=False)
                    s= "Done plotting RV offsets for dataset 3"
                    if quiet==False:
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
                if quiet==False:
                    print s
                log.write(s+'\n')                 
                plt.close()
            except:
                plt.close()
                s= '\n!!! A problem occured while trying to make RV offsets to plot !!!\n' 
                if quiet==False:
                    print s
                log.write(s+'\n')     
        else:
            s= '\nNo RV offsets to plot' 
            if quiet==False:
                print s
            log.write(s+'\n') 
            
        # record the time the chain finished and print
        endTime = timeit.default_timer()
        totalTime = (endTime-startTime) # in seconds
        totalTimeString = genTools.timeString(totalTime)
        s= '\n\nprogessPlotterSingleFile: Plotting took '+totalTimeString+' to complete.\n'
        if quiet==False:
            print s
        log.write(s+'\n')     
        log.write('\n'+75*'*'+'\n Leaving progessPlotterSingleFile \n'+75*'*'+'\n')
        log.close()
    else:
        s= "\nprogessPlotterSingleFile: ERROR!!!! file doesn't exist"
        if quiet==False:
            print s
        log.write(s+'\n')
        log.write('\n'+75*'*'+'\n Leaving progessPlotterSingleFile \n'+75*'*'+'\n')
        log.close()

def progessPlotterSingleFileFunc(log,subPlot2,outputDataFilename,xlabel,paramColNum,xCenters,numChunks,chunkSize, ignoreConstParam=False):
    """
    works for progessPlotterSingleFile to do the plot for each param and reduce code doubling
    """
    verboseInternal = False
    #(CLevels,data,bestDataVal) =genTools.confLevelFinder(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=False, returnBestDataVal=True,fast=True,)
    (log,data,chiSquareds,[bestDataVal,dataMedian,dataValueStart,dataValueMid,dataValueEnd]) = genTools.dataReader(outputDataFilename, paramColNum, returnData=True, returnChiSquareds=False, returnBestDataVal=True, ignoreConstParam=ignoreConstParam)
    #log.write('\nCLevels for '+xlabel+':\n'+repr(CLevels)+'\n')
    if type(data)!=float:
        if len(data)>0:
            yChunks = data[:chunkSize*numChunks].reshape((-1,chunkSize))
            max_env = yChunks.max(axis=1)
            min_env = yChunks.min(axis=1)
            yCenters = yChunks.mean(axis=1)
        else:
            yChunks = np.array([dataValueStart]*len(xCenters))
            yCenters = yChunks
#         if paramColNum==0:
#             log.write("len data = "+str(data.size))
#             log.write("\n\nData Vals:  \n")
#             for line in range(0,data.size):
#                 log.write(str(data[line])+'\n')
#             log.write("\n\nyCenters: \n")#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#             for line in range(0,yCenters.size):
#                 log.write(str(yCenters[line])+'\n')#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#             log.write("\n\n")
        if verboseInternal:
            print "\n xLabel = "+xlabel#$$$$$$$$$$$$$$$$$$$$
            print "len(xCenters) = "+str(len(xCenters))#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            print "len(yChunks) = "+str(len(yChunks))#$$$$$$$$$$$$$$$$$$$$$$
            print "min_env = "+str(min_env)#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            print "max_env = "+str(max_env)#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        subPlot2.fill_between(xCenters, min_env, max_env, color='gray',edgecolor='none', alpha=0.5)
        subPlot2.plot(xCenters, yCenters)
        subPlot2.plot([xCenters.min(),xCenters.max()],[bestDataVal,bestDataVal],color='green')
    else:
        subPlot2.plot([xCenters.min(),xCenters.max()],[bestDataVal,bestDataVal],color='red')
    subPlot2.axes.set_ylabel(xlabel)
    
    return (log,subPlot2,data,bestDataVal)

def orbitEllipsePlotter(longAN_deg, e, period, inc, argPeri_deg, a, To, sysDataDict, DIdataDict,\
                            xLabel='E ["]', yLabel='N ["]', \
                          plotFilename='', xLim=False, yLim=False, show=True, telescopeView=False,nuDI=1,alf=1.0):#,chiSquareds=False):
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
    quiet = True
    verboseInternal = False
    colorsList =['Blue','BlueViolet','Chartreuse','Fuchsia','Crimson','Aqua','Gold','DarkCyan','OrangeRed','Plum','DarkGreen','Chocolate','SteelBlue ','Teal','Salmon','Brown']
    mas = False
    if mas:
        asConversion = 1000.0
    else:
        asConversion = 1.0
        
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    #plt.rc('font',family='serif')
    plt.rc('text', usetex=False) 
    
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
    
    ###################################################################################################################
    ## calculate the locations of companion for 'numOrbs' locations throughout the orbit to make an orbit ellipse    
    ###################################################################################################################
    ellipseXs2 = []
    ellipseYs2 = []
    orbitTAs = []
    orbitSAs = []
    orbitPAs = []
    #sep_dists = []
    if len(longAN_deg)>len(colorsList):
        print "\nAbout to produce ellipse for each set of orbital elements"
        P =ProgressBar('red',width=30,block='=',empty='-',lastblock='>')
    for orb in range(0,len(longAN_deg)):
        #print "Producing ellipse for orb # "+str(orb+1)+"/"+str(len(longAN_deg))
        ellipseXs = []
        ellipseYs = []
        numSteps = 1000.0
        if len(longAN_deg)>100:
            numSteps = 1000.0
        periodIncrement = (period[orb]*365.25)/numSteps
        t = 1.0 
        numFailed = 0
        numPassed = 0
        for step in range(0,int(numSteps)):
            T = 0.0
            t = t + periodIncrement
            (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SA, PA, x, y, a1, a2) =\
                diTools.orbitCalculatorTH_I(t, sys_dist, inc[orb], longAN_deg[orb], e[orb], T, period[orb], argPeri_deg[orb], a[orb],\
                                                            Mass1=1, Mass2=1, verbose=False)
            orbitTAs.append(TA_deg)
            orbitPAs.append(PA)
            orbitSAs.append(SA)
            ellipseX = y*asConversion#SA*math.sin(math.radians(PA))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
            ellipseY = x*asConversion#SA*math.cos(math.radians(PA))*asConversion#*sys_dist   # This will be in [mas] instead of [AU]
            if telescopeView:
                ellipseX = -ellipseX
                ellipseY = -ellipseY
            if True:
                ############################$$$$$$$$$$$$$$
                #$$ Crepp2012 test $$$$$$$$$$$$$$$$$$$$$$$
                ##########################$$$$$$$$$$$$$$$$
                rOvera = Sep_Dist_AU_OP/a2
                r = Sep_Dist_AU_OP/sys_dist
                aUse = ((a1+a2)/sys_dist)*asConversion
                x1=r*(math.cos(math.radians(longAN_deg[orb]))*math.cos(math.radians(argPeri_deg[orb]+TA_deg)) -  math.sin(math.radians(longAN_deg[orb]))*math.sin(math.radians(argPeri_deg[orb]+TA_deg))*math.cos(math.radians(inc[orb])))
                y1=r*(math.sin(math.radians(longAN_deg[orb]))*math.cos(math.radians(argPeri_deg[orb]+TA_deg)) +  math.cos(math.radians(longAN_deg[orb]))*math.sin(math.radians(argPeri_deg[orb]+TA_deg))*math.cos(math.radians(inc[orb])))
                if (abs(x1-x)>1e-6)or (abs(y1-y)>1e-6):
                    numFailed+=1
                    print '\n\n\nx CREPP = '+str(x1)+", y CREPP = "+str(y1)
                    print "x TH-I = "+str(x)+", y TH-I = "+str(y)
                    print 'ellipseX = '+str(ellipseX)+", ellipseY = "+str(ellipseY)+"\n\n\n"
                else:
                    numPassed+=1
            ellipseXs.append(ellipseX)#$$$$$$$$$$$$$$$$$$4
            ellipseYs.append(ellipseY)#$$$$$$$$$$$$$$$$$$$$4
            
            #sep_dist = math.sqrt(math.pow(ellipseX,2.0)+math.pow(ellipseY,2.0))
            #sep_dists.append(sep_dist)
        ellipseXs2.append(ellipseXs)
        ellipseYs2.append(ellipseYs)
        if len(longAN_deg)>len(colorsList):
            P.render((orb+1) * 100 // len(longAN_deg), 'Complete so far.')
        if verboseInternal:
            print "\nThe number of epochs that had mismatched for x,y due to TH-I vs CREPP methods were #"+str(numFailed)+"/"+str(numPassed)+"\n"
    ###################### ORBIT DEBUGGING ####################################
    if False:
        print '\n\n'+'*'*75+'\n      TA             PA                 SA'
        for i in range(0,len(orbitTAs)):
            print str(orbitTAs[i])+"   "+str(orbitPAs[i])+"    "+str(orbitSAs[i])
        print '*'*75+"\n\n"
    ###########################################################################
    
    #########################################################################################
    #$$ Fixed value code to plot errors of fit(s).  ie. the 68% or 95% min and maxes $$$$$$$
    #########################################################################################
    ellipseErrorsXs2 = []
    ellipseErrorsYs2 = []
    plotFitErrors = False
    if plotFitErrors:
        incErrors = [33.0,33.1,inc[0],  19.0,19.1,inc[1],  40.0,40.2,inc[2]]
        longAN_degErrors = [125.1,145.3,longAN_deg[0], 92.0,131.0,longAN_deg[1], 137.5,148.5,longAN_deg[2]]
        eErrors = [0,0,0.0,0.0,0,0,0,0,0]
        periodErrors = [54.5,60.9,period[0],  48.5,53.9,period[1],  61.5,67.0,period[2]]
        argPeri_degErrors = [90,90,90.0,90.0,90,90,90,90,90]
        aErrors = [16.7,17.,a[1],  14.99,15.53,a[1],  18.7,19.2,a[2]]

        for orb in range(0,len(longAN_degErrors)):
            ellipseErrorsXs = []
            ellipseErrorsYs = []
            numSteps = 500.0
            periodIncrement = (periodErrors[orb]*365.25)/numSteps
            t = 1.0 
            for step in range(0,int(numSteps)):
                T = 0.0
                t = t + periodIncrement
                (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SA, PA, x, y, a1, a2) =\
                    diTools.orbitCalculatorTH_I(t, sys_dist, incErrors[orb], longAN_degErrors[orb], eErrors[orb], T, periodErrors[orb], argPeri_degErrors[orb], aErrors[orb],\
                                                                Mass1=1, Mass2=1, verbose=False)
                #orbitTAs.append(TA_deg)
                #orbitPAs.append(PA)
                #orbitSAs.append(SA)
                ellipseErrorsX = y*asConversion#SA*math.sin(math.radians(PA))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
                ellipseErrorsY = x*asConversion#SA*math.cos(math.radians(PA))*asConversion#*sys_dist   # This will be in [mas] instead of [AU]
                if telescopeView:
                    ellipseErrorsX = -ellipseErrorsX
                    ellipseErrorsY = -ellipseErrorsY
                ellipseErrorsXs.append(ellipseErrorsX)
                ellipseErrorsYs.append(ellipseErrorsY)
                #sep_dist = math.sqrt(math.pow(ellipseX,2.0)+math.pow(ellipseY,2.0))
                #sep_dists.append(sep_dist)
            ellipseErrorsXs2.append(ellipseErrorsXs)
            ellipseErrorsYs2.append(ellipseErrorsYs)
        
    
    ##************************************************************************************
    ##****************** Calculate the predicted location and fixed JDs ******************
    ##************************************************************************************
    makePredictions = False
    predictedXs = []
    predictedYs = []
    PredictedLocationString=''
    predictedPatches = []
    if makePredictions:
        if "HR7672" in plotFilename:
            ts = [2452144.000000,
                  2452254.000000,
                  2452473.0,
                  2453989.0,
                  2454367.0,
                  2455697.0             
                  ]
        else:
            ts = [2451965,
            2453425.282,
            2454127.085,
            2454185.028,
            2454306.737,
            2454481.192,
            2454555.068,
            2454584.106,
            2455206,
            2455278
            ]

        for orb in range(0,len(longAN_deg)):
            s= '@'*50+'\n'
            s+="for Orbit # "+str(orb)+"\n"
            s=s+ "Calculating the predicted location of the companion for a fixed list of epochs:\n\n"
            s2=''
            s3=''
            print 'To[orb] = '+str(To[orb])#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            tcount = 0
            for t in ts:
                (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SA, PA, x, y, a1, a2) = \
                    diTools.orbitCalculatorTH_I(t, sys_dist, inc[orb], longAN_deg[orb], e[orb], To[orb], period[orb], argPeri_deg[orb], a[orb],\
                                                                Mass1=1, Mass2=1, verbose=False)
                
                x2 = y*asConversion#SA*math.sin(math.radians(PA))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
                y2 = x*asConversion#SA*math.cos(math.radians(PA))*asConversion#*sys_dist   # This will be in [mas] instead of [AU]
                ############################$$$$$$$$$$$$$$
                #$$ Crepp2012 test $$$$$$$$$$$$$$$$$$$$$$$
                ##########################$$$$$$$$$$$$$$$$
                #rOvera = Sep_Dist_AU_OP/a2
                r = Sep_Dist_AU_OP/sys_dist
                #aUse = ((a1+a2)/sys_dist)*asConversion
                x1=r*(math.cos(math.radians(longAN_deg[orb]))*math.cos(math.radians(argPeri_deg[orb]+TA_deg)) - math.sin(math.radians(longAN_deg[orb]))*math.sin(math.radians(argPeri_deg[orb]+TA_deg))*math.cos(math.radians(inc[orb])))
                y1=r*(math.sin(math.radians(longAN_deg[orb]))*math.cos(math.radians(argPeri_deg[orb]+TA_deg)) + math.cos(math.radians(longAN_deg[orb]))*math.sin(math.radians(argPeri_deg[orb]+TA_deg))*math.cos(math.radians(inc[orb])))
                #(A,B,C,F,G) = diTools.ABCFG_values(aUse, math.radians(argPeri_deg[orb]), math.radians(longAN_deg[orb]), math.radians(inc[orb]))
                # Calc Normalized rectangular coordinates
                #X = math.cos(math.radians(E_latest_deg))-e[orb]
                #X2 = math.cos(math.radians(TA_deg))
                #Y = math.sqrt(1.0-e[orb]**2.0)*math.sin(math.radians(E_latest_deg))
                #Y2 = math.sin(math.radians(TA_deg))
                #print "math.cos(math.radians(E_latest_deg))-e[orb] = "+str(X)+", rOvera*math.cos(math.radians(TA_deg)) = "+str(rOvera*X2)
                #print "math.sqrt(1.0-e[orb]**2.0)*math.sin(math.radians(E_latest_deg)) = "+str(Y)+", rOvera*math.sin(math.radians(TA_deg)) = "+str(rOvera*Y2)
                # Calc x,y values on same coord system as plane of sky (same as data)
                #x_model = A*X+F*Y
                #y_model = B*X+G*Y
                xCent = SAs[tcount]*math.sin(math.radians(PAs[tcount]))*asConversion
                yCent = SAs[tcount]*math.cos(math.radians(PAs[tcount]))*asConversion
                #xCent2 = SAs[tcount]*math.cos(math.radians(PAs[tcount]))*asConversion
                #yCent2 = SAs[tcount]*math.sin(math.radians(PAs[tcount]))*asConversion
                ##
                #x3,y3 = diTools.orbitCalculatorTH_Itoxy(aUse, math.radians(argPeri_deg[orb]), math.radians(longAN_deg[orb]), math.radians(inc[orb]), t, e[orb], To[orb], period[orb])
                ##
                tcount+=1
                if False:
                    print "Epoch = "+str(t)
                    print "E_latest_deg = "+str(E_latest_deg)
                if True:
                    print "Epoch = "+str(t)
                    print 'a1 = '+str(a1)
                    #print "E_latest_deg = "+str(E_latest_deg)
                    print 'plot axis x CREPP  = '+str(y1)+", plot axis y CREPP = "+str(x1)
                    #print 'x_THI = '+str(x_model)+', y_THI = '+str(y_model)
                    #print "x_THI-new = "+str(x3)+", y_THI-new = "+str(y3)
                    #print "x_data (SA*COS(PA)) = "+str(xCent2)+", y_data  (SA*SIN(PA)) = "+str(yCent2)
                    #print 'x MINE (SA*cos(PA))= '+str(SA*math.cos(math.radians(PA))*asConversion)+", y MINE (SA*sin(PA))= "+str(SA*math.sin(math.radians(PA))*asConversion)
                    print "x_data (SA*SIN(PA)) = "+str(xCent)+", y_data  (SA*cos(PA)) = "+str(yCent)
                    print 'x MINE (SA*SIN(PA))= '+str(x2)+", y MINE (SA*cos(PA))= "+str(y2)
                    print "\n"
                
                #########################$$$$$$$$$$$$$$$$$
                predictedXs.append(x2)
                predictedYs.append(y2)
                predictedLocPatch = star((asConversion/1000.0)*0.6*a[0], x2, y2, color=colorsList[orb], N=8, thin = 0.5)
                predictedPatches.append(predictedLocPatch)
                s=s+ 'Epoch '+str(t)+': x = '+str(x)+' , y = '+str(y)+'\n'
                s2=s2+str(x2)+'    '+str(y2)+'\n'
                s3=s3+str(SA)+'     '+str(PA)+'\n'
            s=s+ '@'*50+"\n"
            #s=s+ '\n       *** Excel Format ***\n       x                         y      \n'
            #s=s+ s2
            #s = s+'\n       SA                  PA        \n'
            #s=s+s3
            if verboseInternal:
                print s
            log.write(s+'\n')
            PredictedLocationString+= s+"\n"
        
    ##************************************************************************************
    ##************************************************************************************
    ##************************************************************************************
    
    ###########################################################################################
    # Calculate the location of the secondary at the beginning, 1/4, 1/2, 3/4 and end of orbit.
    ###########################################################################################
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, x, y, a1, a2) =\
            diTools.orbitCalculatorTH_I(0.0, sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    Xstart = SAstart*math.sin(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    Ystart = SAstart*math.cos(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    if telescopeView:
        Xstart = -Xstart
        Ystart = -Ystart
    startStar = star((asConversion/1000.0)*0.6*a[0], Xstart, Ystart, color='green', N=20, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, x, y, a1, a2) =\
            diTools.orbitCalculatorTH_I(((period[0]/4.0)*365.25), sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    XoneQuarter = SAstart*math.sin(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    YoneQuarter = SAstart*math.cos(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    if telescopeView:
        XoneQuarter = -XoneQuarter
        YoneQuarter = -YoneQuarter
    oneQuarterStar = star((asConversion/1000.0)*0.6*a[0], XoneQuarter, YoneQuarter, color='blue', N=20, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAhalf, PAhalf, x, y, a1, a2) =\
            diTools.orbitCalculatorTH_I(((period[0]/2.0)*365.25), sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    Xhalf = SAhalf*math.sin(math.radians(PAhalf))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    Yhalf = SAhalf*math.cos(math.radians(PAhalf))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    if telescopeView:
        Xhalf = -Xhalf
        Yhalf = -Yhalf
    halfStar = star((asConversion/1000.0)*0.6*a[0], Xhalf, Yhalf, color='blue', N=20, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, x, y, a1, a2) =\
            diTools.orbitCalculatorTH_I((3.0*(period[0]/4.0)*365.25), sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    XthreeQuarter = SAstart*math.sin(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    YthreeQuarter = SAstart*math.cos(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    if telescopeView:
        XthreeQuarter = -XthreeQuarter
        YthreeQuarter = -YthreeQuarter
    threeQuarterStar = star((asConversion/1000.0)*0.6*a[0], XthreeQuarter, YthreeQuarter, color='blue', N=20, thin = 0.2)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, x, y, a1, a2) =\
            diTools.orbitCalculatorTH_I(((period[0])*365.25), sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
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
    chiSquaredStr=''
    SA_diffs2 = []
    PA_diffs2 = []
    if len(longAN_deg)>len(colorsList):
        print '\nAbout to calculate residuals'
        P =ProgressBar('red',width=30,block='=',empty='-',lastblock='>')
    for orb in range(0,len(longAN_deg)):
        #print 'Finding residuals for orbit # '+str(orb+1)+"/"+str(len(longAN_deg))
#         if len(longAN_deg)<len(colorsList):
#             legendStr = legendStr+"\nFor orbit # "+str(orb)+' (color = '+colorsList[orb]+') :\n'
#             legendStr = legendStr+"inc[orb] = "+str(inc[orb])+'\n'
#             legendStr = legendStr+"longAN_deg[orb] = "+str(longAN_deg[orb])+'\n'
#             legendStr = legendStr+"e[orb] = "+str(e[orb])+'\n'
#             legendStr = legendStr+"To[orb] = "+str(To[orb])+'\n'
#             legendStr = legendStr+"period[orb] = "+str(period[orb])+'\n'
#             legendStr = legendStr+"argPeri_deg[orb] = "+str(argPeri_deg[orb])+'\n'
#             legendStr = legendStr+"a[orb] = "+str(a[orb])+'\n'
        if False:
            (chi_squared_total, ns, Ms, Es, thetas, Sep_Dists, SA_arcsec_measured_models, PA_deg_measured_models, xs, ys, a1s, a2s) =\
            diTools.multiEpochOrbCalc(SAs, SAerrors, PAs, PAerrors,epochs, sys_dist, inc[orb], longAN_deg[orb],\
                            e[orb], To[orb], period[orb], argPeri_deg[orb], a_total=a[orb], Mass1=1, Mass2=1, verbose=False, useTHI=False)
            if verboseInternal or True:
                print "!!!!!!!!!!!!!!!!!!!!\n fullEQs: for orbit #"+str(orb)+", the chiSquared fit to the DI data was = "+\
                str(chi_squared_total)+', or reduced = '+str(chi_squared_total/nuDI)+"\n!!!!!!!!!!!!!!!!!!!!!!!!!"
        (chi_squared_total, ns, Ms, Es, thetas, Sep_Dists, SA_arcsec_measured_models, PA_deg_measured_models, xs, ys, a1s, a2s) =\
        diTools.multiEpochOrbCalc(SAs, SAerrors, PAs, PAerrors,epochs, sys_dist, inc[orb], longAN_deg[orb],\
                        e[orb], To[orb], period[orb], argPeri_deg[orb], a_total=a[orb], Mass1=1, Mass2=1, verbose=False)
        
        if type(SAs)!=np.ndarray:
            SAs = np.array(SAs)
        if type(PAs)!=np.ndarray:
            PAs = np.array(PAs)
        if type(SA_arcsec_measured_models)!=np.ndarray:
            SA_arcsec_measured_models = np.array(SA_arcsec_measured_models)
        if type(PA_deg_measured_models)!=np.ndarray:
            PA_deg_measured_models = np.array(PA_deg_measured_models)
        SA_diffs2.append(SAs-SA_arcsec_measured_models)
        PA_diffs2.append(PAs-PA_deg_measured_models)       
        if len(longAN_deg)<len(colorsList): 
            chiSquaredStr = "The chiSquared fit to the DI data was = "+str(chi_squared_total)+', or reduced = '+str(chi_squared_total/nuDI)
#             legendStr = legendStr+chiSquaredStr+'\n'
        if verboseInternal:
            print "\n\n     TH-I: for orbit #"+str(orb)+", "+chiSquaredStr+"\n"
        if len(longAN_deg)>len(colorsList):
            P.render((orb+1) * 100 // len(longAN_deg), 'Complete so far.')
    #######################################################################################            
    
    ############################################################################################################################################
    ########################################## Starting to Plot ################################################################################
    ############################################################################################################################################

    ## Create figure, axes and start plotting/drawing everything
    fig = plt.figure(1,figsize=(10,9))
    main = fig.add_subplot(111)
        
    # Draw orbit ERRORS !!
    if False:
        for orb in range(0,len(periodErrors)):
            if ((orb==0)or(orb==1)or(orb==3)or(orb==4)or(orb==6)or(orb==7)):
                main.plot(ellipseErrorsXs2[orb],ellipseErrorsYs2[orb],linewidth=1,color=colorsList[int(orb/3.0)],alpha=1.0) 
                if verboseInternal:
                    print str(orb)+"color = "+colorsList[int(orb/3.0)]
            #else:
            #     main.plot(ellipseErrorsXs2[orb],ellipseErrorsYs2[orb],linewidth=1,color=colorsList[int(orb/3.0)],alpha=1.0) 
        if verboseInternal:
            print "error fills"
        for orb in range(0,len(periodErrors)):
            if ((orb==0)or(orb==3)or(orb==6)):
                main.fill(ellipseErrorsXs2[orb:orb+2],ellipseErrorsYs2[orb:orb+2], color=colorsList[int(orb/3.0)], edgecolor=colorsList[int(orb/3.0)], alpha=0.1)
                if verboseInternal:
                    print str(orb)+"color = "+colorsList[int(orb/3.0)]
                #print repr(ellipseErrorsXs2[orb:orb+2])
    
    # Draw orbits
    if len(longAN_deg)>len(colorsList):
            print '\nAbout to draw ellipses'
            P =ProgressBar('red',width=30,block='=',empty='-',lastblock='>')
    for orb in range(0,len(longAN_deg)):
        if False:
            main.plot(ellipseXs2[orb],ellipseYs2[orb],linewidth=2.5,color=colorsList[orb]) 
        else:
            #print 'Plotting ellipse for orbit # '+str(orb)+"/"+str(len(longAN_deg))
            ## For density plot style
            main.plot(ellipseXs2[orb],ellipseYs2[orb],linewidth=2,color='blue',alpha=alf) 
        if len(longAN_deg)>len(colorsList):
            P.render((orb+1) * 100 // len(longAN_deg), 'Complete so far.')
    
    #draw semi-major
    #main.plot([Xstart,Xhalf],[Ystart,Yhalf],'g-',linewidth=2)
    
    #calculate and draw an X at the center of the semi-major axis
    xCOM = (Xhalf+Xstart)/2.0
    yCOM = (Yhalf+Ystart)/2.0
    if telescopeView:
        yCOM = -yCOM
        xCOM = -xCOM
    #main.plot(xCOM,yCOM,'rx',markersize=12,mew=2.5)
    
    # draw lines along diagonals to mark out 45 degree locations from origin ### Useless now
    #main.plot( [xLim[0],xLim[1]], [yLim[0],yLim[1]])
    #main.plot( [xLim[1],xLim[0]], [yLim[0],yLim[1]])
    
    # Draw stars for the location of each 1/4 of the orbit  # These are optional, I would not include them for the final versions, maybe just the periapsis one.
    if False:
        main.add_patch(startStar)
        main.add_patch(oneQuarterStar)
        main.add_patch(halfStar)
        main.add_patch(threeQuarterStar)
        #main.add_patch(endStar)
    if makePredictions:
        for predStar in predictedPatches:
            main.add_patch(predStar)
    
    # Draw larger star for primary star's location
    # create a primary star polygon at center of plot
    starPolygon = star((asConversion/1000.0)*12.0*a[0], 0, 0, color='yellow', N=6, thin = 0.5)
    main.add_patch(starPolygon)
    
    ## call function to calculate, create and return polygons for the 
    ## companion star locations and boxes for their errors
    if False:
        (errorBoxes, m2starPolygons,[xmin,xmax,ymin,ymax]) = starAndErrorPolys(SAs,SAerrors,PAs,PAerrors,asConversion, main.transData, telescopeView)
    else:
        (main,[xmin,xmax,ymin,ymax]) =  addDIdataToPlot(main,SAs,SAerrors,PAs,PAerrors,asConversion,telescopeView=False)
    dataMaxMins = [xmin,xmax,ymin,ymax]
    if verboseInternal:
        print "original dataMaxMins = "+repr(dataMaxMins)
    
    #print "from starAndErrorPolys [xmin,xmax,ymin,ymax] = "+repr([xmin,xmax,ymin,ymax])
    ## set the xLim and yLim if their values are False
    ## and pad max and min values found by 10%
    dataXrange = abs(dataMaxMins[1])+abs(dataMaxMins[0])
    dataYrange = abs(dataMaxMins[3])+abs(dataMaxMins[2])
    xLimData = []
    if not xLim:
        min = genTools.findArrayMin(ellipseXs2[:])
        predictXmin = genTools.findArrayMin(predictedXs)
        s = "min = "+str(min)
        if makePredictions:
            s+=", predictXmin = "+str(predictXmin)
        if verboseInternal:
            print s
        if min>predictXmin:
            min = predictXmin
        if xmin>predictXmin:
            xmin = predictXmin
        max = genTools.findArrayMax(ellipseXs2[:])
        predictXmax = genTools.findArrayMax(predictedXs)
        s = "max = "+str(max)
        if makePredictions:
            s+=", predictXmax = "+str(predictXmax)
        if verboseInternal:
            print s
        if max<predictXmax:
            max = predictXmax
        if xmax<predictXmax:
            xmax = predictXmax
        Range = abs(max)+abs(min)
        if verboseInternal:
            print "new x range = "+str(Range)
        dataXrange = abs(xmax)+abs(xmin)
        xLimData = (xmin-abs(dataXrange*0.05),xmax+abs(dataXrange*0.05))
        xLim = (min-abs(Range*0.05),max+abs(Range*0.05))
        if verboseInternal:
            print  "xLim = "+repr(xLim)
        xLim = (genTools.findArrayMin([xLim[0],xmin]), genTools.findArrayMax([xLim[1],xmax])) 
        #print 'elipseXs2 min = '+str(min)+", max = "+str(max)+", final xLim = "+repr(xLim)
    else:
        if not (type(xLim)==tuple):
            s= 'PROBLEM: xLim is not of type tuple'
            print s
            log.write(s+'\n')
    yLimData = []
    if not yLim:
        min = genTools.findArrayMin(ellipseYs2[:])
        predictYmin = genTools.findArrayMin(predictedYs)
        s = "min = "+str(min)
        if makePredictions:
            s+=", predictYmin = "+str(predictYmin)
        if verboseInternal:
            print s
        if min>predictYmin:
            min = predictYmin
        if ymin>predictYmin:
            ymin = predictYmin
        max = genTools.findArrayMax(ellipseYs2[:])
        predictYmax = genTools.findArrayMax(predictedYs)
        s = "max = "+str(max)
        if makePredictions:
            s+=", predictYmax = "+str(predictYmax)
        if verboseInternal:
            print s
        if max<predictYmax:
            max = predictYmax
        if ymax<predictYmax:
            ymax = predictYmax
        Range = abs(max)+abs(min)
        dataYrange = abs(ymax)+abs(ymin)
        if verboseInternal:
            print "new y range = "+str(Range)
        yLimData = (ymin-abs(dataYrange*0.05),ymax+abs(dataYrange*0.05))
        yLim = (min-abs(Range*0.05),max+abs(Range*0.05))
        if verboseInternal:
            print  "yLim = "+repr(yLim)
        yLim = (genTools.findArrayMin([yLim[0],ymin]), genTools.findArrayMax([yLim[1],ymax])) 
        #print 'elipseYs2 min = '+str(min)+", max = "+str(max)+", final yLim = "+repr(yLim)
    else:
        if not (type(yLim)==tuple):
            s= 'PROBLEM: yLim is not of type tuple'
            print s
            log.write(s+'\n')
    if True:
        dataMaxMins = [xLim[0],xLim[1],yLim[0],yLim[1]]
        if verboseInternal:
            print "new dataMaxMins = "+repr(dataMaxMins)
    
    # Draw lines for horizontal and vertical from origin
    #main.plot([xLim[0],xLim[1]],[0,0],c='black',linewidth=2)
    #main.plot([0,0],[yLim[0],yLim[1]],c='black',linewidth=2)
    
    if False:
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
    # Fix up borders and tick labels to ideal sizes
    #main.axhline(linewidth=2.0,color='k')
    #main.axvline(linewidth=2.0,color='k')
    main.tick_params(axis='both',which='major',width=1,length=3,pad=10,direction='in',labelsize=10)
    #plot.axhline(linewidth=2.0)
    #plot.axvline(linewidth=2.0)
    #plot.xaxis.set_tick_params(width=5,length=5,pad=2,direction='in')
    #plot.yaxis.set_tick_params(width=5,length=5,pad=2,direction='in')
    main.spines['right'].set_linewidth(1.0)
    main.spines['bottom'].set_linewidth(1.0)
    main.spines['top'].set_linewidth(1.0)
    main.spines['left'].set_linewidth(1.0)
    main.set_position([0.19,0.15,0.77,0.80])
    main.set_xlabel(xLabel, fontsize=30)
    main.set_ylabel(yLabel, fontsize=30)
    #plt.suptitle(plotFileTitle, fontsize=20)

    # save plot to file
    if plotFilename!='':
        plt.savefig(plotFilename, dpi=300, orientation='landscape')
    
    ## find cropped limits for just data plus 5% pad
    try:
        xlimData = [xLimData[1],xLimData[0]]
        ylimData = yLimData
    except:
        xlimData = [dataMaxMins[1]+0.0005*dataXrange,dataMaxMins[0]-0.0005*dataXrange]
        ylimData = [dataMaxMins[2]-0.0005*dataYrange,dataMaxMins[3]+0.0005*dataYrange]
    if verboseInternal:
            print "those used in cropped plot xlimData = "+repr(xlimData)+", ylimData = "+repr(ylimData)
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
    
    #Create a figure for the Residuals in SA and PA
    if True:
        fig2 = plt.figure(1,figsize=(30,20))
        saPlot = fig2.add_subplot(211)
        saPlot.set_xlabel('Epochs [JD]', fontsize=30)
        saPlot.set_ylabel('SA Residuals ["]', fontsize=30)
        saPlot.tick_params(axis='both',which='major',width=1,length=3,pad=10,direction='in',labelsize=25)
        paPlot = fig2.add_subplot(212)
        paPlot.set_xlabel('Epochs [JD]', fontsize=30)
        paPlot.set_ylabel("PA Residuals [deg]", fontsize=30)
        paPlot.tick_params(axis='both',which='major',width=1,length=3,pad=10,direction='in',labelsize=25)
        #plt.suptitle(plotFileTitle, fontsize=10)
        if len(longAN_deg)>len(colorsList):
            print '\nAbout to Draw residuals'
            P =ProgressBar('red',width=30,block='=',empty='-',lastblock='>')
        for orb in range(0,len(longAN_deg)):
            if len(longAN_deg)>len(colorsList):
                #print 'Plotting residuals for orbit # '+str(orb+1)+"/"+str(len(longAN_deg))
                saPlot.scatter(epochs,SA_diffs2[orb],color='blue',alpha=alf)
                paPlot.scatter(epochs,PA_diffs2[orb],color='blue',alpha=alf)
            else:
                saPlot.scatter(epochs,SA_diffs2[orb],color=colorsList[int(orb/3.0)])
                paPlot.scatter(epochs,PA_diffs2[orb],color=colorsList[int(orb/3.0)])
            if len(longAN_deg)>len(colorsList):
                P.render((orb+1) * 100 // len(longAN_deg), 'Complete so far.')
        saPlot.plot(saPlot.axes.get_xlim(),[0,0],c='r',linewidth=2.0)
        paPlot.plot(paPlot.axes.get_xlim(),[0,0],c='r',linewidth=2.0)
        plotFilenameOC = plotFilename[:-4]+"-Residuals.png"
        ## saved cropped version to file if requested
        if plotFilename!='':
            plt.savefig(plotFilenameOC, dpi=300, orientation='landscape')
        plt.close()  
        
    
    if True:
        ## Create figure for writting the sorta legend to
        fig = plt.figure(1,figsize=(10,20))
        #main = fig.add_subplot(111)
        if True:
            legendStr = legendStr+"\n\n"+PredictedLocationString
        fig.text(0.05,0.05,legendStr,ha='left')
        if plotFilename!='':
            legendFigFilename = plotFilename[:-4]+"-paramInfo.png"
            plt.savefig(legendFigFilename, dpi=300, orientation='landscape')
        plt.close()
    
    # write last log statement and close log
    log.write('\n'+75*'#'+'\n Leaving orbitEllipsePlotter \n'+75*'#'+'\n')
    log.close()
    
    return chiSquaredStr
  

def mcmcProgressPlotter(filenames,rootPlotFilename,nu=1, plot4x1=False,TcStepping=False):
    """
    NOTE: For a set of output chain's data files.
    
    This will produce a time series type view of the paramters in each chain's output
    data file.  It will utilize progessPlotterSingleFile to produce these plots and figures.
    
    These plots will be made for the parameter values vs step, 
    to show how it is converging with time.
    
    NOTE: Runs on NEW data files only.
    
    :param nu:       The value of nu to convert chi squared to reduced chi squared.
    :type nu:        float
    :param plot4x1:  Only make plots for K, argPeri, e, To (or Tc if Tc stepping was used), and chiSquared.
    :type plot4x1:   Python boolean
    """
    verbose = False
    # check if the passed in value for plotFilename includes '.png'
    if '.png' not in rootPlotFilename:
        rootPlotFilename = rootPlotFilename+'.png'
        
    # check filenames parameter
    if type(filenames) is not list:
        print 'the filenames input parameter must be a list of strings'
    
    #numFiles = len(filenames)
    
#    ## get settings used for simulation $$$ only need to see if a 3x2 or 3x1 plot should be used for parameter plots
#    inputSettingsFile = os.path.join("/run/media/Kyle/Data1/Todai_Work/Dropbox/workspace/Binary-project/SimSettings_and_InputData",'SimSettings.txt')
#    # make an output settings file name for the duo version output
#    [fName,ext] = os.path.splitext(inputSettingsFile)
#    outputSettingsFile = fName+'_DuoVersion'+ext
#    paramSettingsDict = cFileToSimSettingsDict(outputSettingsFile)
#    
#    print '&& The value of paramSettingsDict["simulate_StarPlanetRV"] = '+repr(paramSettingsDict["simulate_StarPlanetRV"])

    for filename in filenames:
        # check if the passed in value for filename includes '.txt'
        if (filename[-4:]!='.txt' and filename[-4:]!='.dat'):
            filename = filename+'.dat'
            if verbose:
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
        if verbose:
            print s
        log.write(s+'\n')
        log.close()
        #call plotter for this file
        progessPlotterSingleFile(filename, plotFilename, nu=nu, plot4x1=plot4x1,TcStepping=TcStepping)

def variousPlotTester(outputDatafile=''):
    """
    NOTE: not a standard function.  A manual/custom function used at times.  
    
    A function to allow custom inputs to orbitEllipsePlotter for testing or 
    producing one off plots for a particular set of orbital Elements.
    There is also some code here for producing a scatter plot of 2 variables.
    """
    
    sysDatafilename = os.path.join("/mnt/Data1/Todai_Work/Dropbox/workspace/Binary-project/SimSettings_and_InputData",'SystemData.txt')
    sysDataDict = genTools.sysDataToDict(sysDatafilename)
    
    DIdatafilename = os.path.join("/mnt/Data1/Todai_Work/Dropbox/workspace/Binary-project/SimSettings_and_InputData",'DIdata.dat')
    if os.path.exists(DIdatafilename):
        DIdataDict = diTools.DIdataToDict(DIdatafilename)
    if outputDatafile=='':
        outputDatafile = "/mnt/Data1/Todai_Work/Data/data_Binary/data_Duo/MCMC-parmVaryCorrected-HR8799-UniEccentPrior-circular-veryOpenParms-1--28-Million-in_Total/outputData-ALL.dat"
    
    print '#'*50
    bestOrbit = genTools.bestOrbitFinder(outputDatafile, printToScreen=True, saveToFile=False, returnAsList=True)
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
        #(CLevels,data,chiSquareds,bestDataVal) =genTools.confLevelFinder(outputDataFilename,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True)
        #print 'Argument of Perigie data loaded'
        
        yLabel = 'Inclination [deg]'
        paramColNum = 4
        (CLevels,yData,chiSquareds,bestDataVal) =genTools.confLevelFinder(outputDatafile,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True)
        print 'Inclination data loaded' 
                    
        #xlabel = 'Longitude of Ascending Node [deg]'
        paramColNum = 0
        #(CLevels,data,chiSquareds,bestDataVal) =genTools.confLevelFinder(outputDatafile,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True)
        #print 'Longitude of Ascending Node data loaded'
        
        #xlabel = 'e'
        paramColNum = 1
        #(CLevels,data,chiSquareds,bestDataVal) =genTools.confLevelFinder(outputDatafile,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True)
        #print 'e data loaded'
        
        xLabel = 'Period [Years]'
        paramColNum = 3
        (CLevels,xData,chiSquareds,bestDataVal) =genTools.confLevelFinder(outputDatafile,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True)
        print 'Period data loaded'
        
        #xlabel = 'Time of last Periapsis [JD]'
        paramColNum = 2
        #(CLevels,data,chiSquareds,bestDataVal) =genTools.confLevelFinder(outputDatafile,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True)
        #print 'Time of last Periapsis data loaded'
        
        #xlabel = 'K [m/s]'
        paramColNum = 8
        #(CLevels,data,chiSquareds,bestDataVal) =genTools.confLevelFinder(outputDatafile,paramColNum, returnData=True, returnChiSquareds=True, returnBestDataVal=True)
        #print 'K data loaded'
        
        histAndContourPlotter(xData, yData, chiSquareds, xLabel='', yLabel='', plotFilename='', xLim=False, yLim=False, show=True, save=False)
    
def PostSimCompleteAnalysisFunc(outputDatafile=''):
    """
    A function that performs many of the post simulation completion plotting 
    or Gelman Rubin statistic calculations as done in ProcessManager.  This 
    can be used for custom plotting or testing, or to perform any plotting that 
    was turned off in the initial run of a simulation that you might now want.
    """
    ####################################################
    ## Determine the files and riectories for input files
    ####################################################
    baseDir = ""
    prepend = ""
    if outputDatafile=='':
        #baseDir = "/mnt/Data1/Todai_Work/Dropbox/EclipseWorkspace/SMODT/settings_and_InputData"
        outputDatafile = "/mnt/Data1/Todai_Work/Data/data_SMODT/FakeData-mcmc-3D-veryOpen-5Percent--500-Thousand-in_Total/outputData-ALL.dat" 
        baseDir = os.path.dirname(outputDatafile)
        prepend = "FakeData_"
    elif outputDatafile!="":
        baseDir = os.path.dirname(outputDatafile)
        prepend = baseDir.split("/")[-1].split("-")[0]+"_"
        print "\n\nUSING the prepend '"+prepend+"'\n\n"
    summaryPlotFile = os.path.join(baseDir,"SummaryPlot-Manual")
    sysDatafilename = os.path.join(baseDir,"code-used/"+prepend+'SystemData.txt')
    RVdatafilename = os.path.join(baseDir,"code-used/"+prepend+'RVdata.dat')    
    inputSettingsFile = os.path.join(baseDir,"code-used/"+prepend+'SimSettings.txt')
    if os.path.exists(inputSettingsFile)==False:
        inputSettingsFile = os.path.join(baseDir,"code-used/"+prepend+'SimSettings_DuoVersion.txt')
    #print "inputSettingsFile found to be "+inputSettingsFile
    ####################################################
    # Get data dicts
    ####################################################
    sysDataDict = genTools.sysDataToDict(sysDatafilename)
    try:
        RVdataDict = rvTools.RVdataToDict(RVdatafilename)
    except: 
        useless=True
    DIdatafilename = os.path.join(baseDir,"code-used/"+prepend+'DIdata.dat')
    try:
        DIdataDict = diTools.DIdataToDict(DIdatafilename)
    except:
        uselss = True
    paramSettingsDict = genTools.cFileToSimSettingsDict(inputSettingsFile)
    ####################################################
    # Make output plot filenames
    ####################################################
    numModDataSets = 10
    rvPlotFilename = os.path.join(baseDir,"RVplot-Manual")
    diPlotFilename = os.path.join(baseDir,"DIplot-Manual")
    modDatasetsFilename =  os.path.join("/mnt/Data1/Todai_Work/Data/data_Binary/data_Duo","mod"+str(numModDataSets)+"Dataset.dat")
    TvsEccPlotFilename = os.path.join("/mnt/Data1/Todai_Work/Data/data_Binary/data_Duo","TvsEccentricityPlot")
    ####################################################
    ## get nu value, then calculate chiSquared cut off
    # get log
    ####################################################
    logFilename = os.path.join(baseDir,'log-chain_1.txt')        
    [nu,nuRV,nuDI,printStr] = genTools.findNuFromLog(logFilename)

    ####################################################
    # combine output chain data files into one $$$$$$$$$ Not sure why you would need this anymore though...$$$
    ####################################################
    if False:
        numChains = 7
        dataFiles = []
        dataFinalFilename = os.path.join(baseDir,'outputData-ALL.dat')
        for num in range(1,numChains+1):
            fname = os.path.join(baseDir,'outputData-chain_'+str(num)+'.dat')
            dataFiles.append(fname)
            print 'Adding filename to list: '+fname
        genTools.dataFileCombiner(dataFiles, dataFinalFilename)
        
    ####################################################
    # Perform orbit plotting, either DI, RV or both
    ####################################################
    #bestOrbit = [61.0,0.5,2456847.0,0,]
    if False:    
        bestOrbit = []  
        try:
            #print '#'*50
            bestOrbit = genTools.bestOrbitFinder(outputDatafile, printToScreen=False, saveToFile=False, returnAsList=True)
            #print '#'*50
            longAN = bestOrbit[0]
            e = bestOrbit[1]
            period = bestOrbit[4]
            inc = bestOrbit[5]
            argPeri_deg = bestOrbit[6]+180.0
            a = bestOrbit[7]
            T = 2457000#bestOrbit[2]
            Tc = 0#bestOrbit[3]
            K = bestOrbit[8]
            if paramSettingsDict["DIonly"]==False:                
                RVoffsets=bestOrbit[9]
        except:
            print 'No data file was provided and no bestOrbit.txt file existed in the directory provided, so the best orbit values could not be found!!!!'
        #################################################### 
        ## Now manually choose what you want to plot
        ####################################################
        ##############   RV plot(s)  #######################
        if False:
            ## make RV fit plot
            rvPlotter(e,T,Tc,period,inc,argPeri_deg,a, \
                      sysDataDict,RVdataDict,paramSettingsDict,K=K,RVoffsets=RVoffsets,\
                      nuRV=nuRV,plotFilename=rvPlotFilename+'-FullOrbit', show=False, plotFullOrbit=True)
            if False:
                rvPlotter(e,T,Tc,period,inc,argPeri_deg,a, \
                          sysDataDict,RVdataDict,paramSettingsDict,K=K,RVoffsets=RVoffsets,\
                          plotFilename=rvPlotFilename, show=False, plotFullOrbit=False)
        if False:
            ## Make bootstrap type data for RV data
            rvModDatasetMaker(e, T, period, inc, argPeri_deg, a, sysDataDict, RVdataDict, paramSettingsDict,\
                     RVoffsets=RVoffsets, modDatasetsFilename=modDatasetsFilename, numModDatasets=numModDataSets)
        ##############   DI plot(s)  #######################
        if False:
            ## Make a DI fit plot
            if False:
                if False:
                    longANs = [136.48,115.03,143.82]
                    argPeris = [90,90,90]
                    incs = [33.04,19.09,40.10]
                    periods = [57.82,51.27,64.21]
                    a_totals = [17.22,15.29,18.99]
                    Ts = [2454635.9,2453560.2,2454996.53]
                    es = [0.0,0.0,0.0]
                    #if 
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
            if True:
                densityPlot(baseDir, sysDataDict, DIdataDict,RVdataDict,paramSettingsDict, nuDI=nuDI,nuRV=nuRV) 
    ####################################################
    # Extra non-orbit plots and statistic calculations
    ####################################################
    ##############   Parameter Summary plot(s)  #######################
    if False:
        ## make MCMC progress plots
        numChains = 7
        fileList = []
        for num in range(1,numChains+1):
            fileList.append(os.path.join(os.path.dirname(outputDatafile),'outputData-chain_'+str(num)+'.txt'))
        progessPlotterSingleFile(outputDatafile, summaryPlotFile+"-MCMCprogress", weight=False, confLevels=True, nu=1, SimPlanetStar=paramSettingsDict["simulate_StarPlanetRV"])
    if False:
        ## Make the posterior prob histograms.
        summaryPlotter(outputDatafile, summaryPlotFile, weight=False, confLevels=True, nu=1, plot4x1=False)
    if True:
        ## Make the posterior prob histograms.
        summaryPlotFile = os.path.join('/mnt/Data1/Todai_Work/Dropbox/SMODT-outputCopies/',"SummaryPlot-Manual")
        summaryPlotter(outputDatafile, summaryPlotFile, shadeConfLevels=True, nu=1)
    if False:
        ## Make the 'cleaned' posterior prob histograms
        makeCleanSummaryPlot(outputDatafile)
    ##############   Scatter plot(s)  #######################
    if False:
        ## make a scatter hist figure ###$$$ not sure if this works right now.
        outDict = genTools.outputDatafileToDict(outputDatafile)
        #esCLevels = ConfLevelFunc(chiSquareds,outDict["es"])
        chiSquareds = outDict["chiSquareds"]
        histAndScatterPlotter(outDict["Ts"], outDict["es"], xLabel='Ts', yLabel='es', plotFilename=TvsEccPlotFilename, xLim=False, yLim=False, show=True, save=False)
    ############### Calc wrap up statistics ####################
    if False:
        ### Perform second round of Gelman-Rubin
        numChains = 7
        fileList = []
        for num in range(1,numChains+1):
            fname = os.path.join(os.path.dirname(outputDatafile),'gelmanRubin-chain_'+str(num)+'.txt')
            fileList.append(fname)
            print 'Adding filename to list: '+fname
        genTools.gelmanRubinStage2(fileList)

def rvPlotter(e, T, Tc, period, inc, argPeri_deg, a, sysDataDict, RVdataDict, paramSettingsDict,\
                 K=0, RVoffsets=[0], nuRV=1, plotFilename='', show=True, plotFullOrbit=True, primaryRVs=True,alf=1.0):
    """
    create a plot for the RV data and a fit line from the best orbit data params.
    
    NOTES:
    #1: If a star-planet system is being simulated, then the provided inclination will be used with the
    planet's mass in the system data file/provided dict, to calculate the residual vel (due to M2sin(i)).
    If the value in the dictionary/file is actually, please set inc=0 to tell this func to ignore it and use the file's value.
    #2: Handling of the argPeriOffsetRV must be done before calling this function.
    #3: primaryRVs parameter in call will set the semi-major value used to a_1 rather than a_2 if K=0, else K provided will be used. 
    Thus, providing K allows the user to ignore this parameter. 
    """
    verbose = False
    addLegend = False
    plotErrorBars = True
    makeTrendPlot = True
    makeDataResidualPlot = False
    
    colorsList =['Blue','BlueViolet','Chartreuse','Fuchsia','Crimson','Aqua','Gold','DarkCyan','OrangeRed','Plum','DarkGreen','Chocolate','SteelBlue ','Teal','Salmon','Brown']
    
    
    ## make string with input values for possible printing
    s= '\nInputs to rvPlotter were:'+'\n'
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
        
    plt.rc('font',family='serif')
    plt.rc('text', usetex=False)
    
    #check the orbit element inputs to ensure they are lists, else make them lists
    if type(K)!=list:
        K = [K]
    try:
        if type(RVoffsets[0])!=list:
            RVoffsets = [RVoffsets]
    except:
        if type(argPeri_deg)!=list:
            RVoffsets = [RVoffsets]
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
    # correct Tc if somehow still zero
    TcUse = []
    for i in range(0,len(Tc)):
        if Tc[i]==0:
            (To,Tcent) = genTools.eccArgPeri2ToTcCalc(e[i], period[i], argPeri_deg[i], T[i], Tc=0)    
            TcUse.append(Tcent)
        else:
            TcUse.append(Tc[i])
    if False:
        print "Tc = "+repr(Tc)
        print "TcUse = "+repr(TcUse)
    Tc = TcUse
    
    RV_epochs = RVdataDict['RV_epochs']
    RV_errors = RVdataDict['RV_errors']
    RVs = RVdataDict['RVs']
    KMperSec = False
    MtoKM = 1.0
    if np.max(RVs)>2000.0:
        KMperSec = True
        MtoKM = 1000.0
    
    if type(RV_epochs[0])!=list:
        RV_epochs = [RV_epochs]
    if type(RV_errors[0])!=list:
        RV_errors = [RV_errors]
    if type(RVs[0])!=list:
        RVs = [RVs]
        
        
    Mass1 = sysDataDict['Mass1'] 
    simulate_StarStarRV = paramSettingsDict["simulate_StarStarRV"]
    simulate_StarPlanetRV = paramSettingsDict["simulate_StarPlanetRV"]
    if (simulate_StarStarRV is True) and(simulate_StarPlanetRV is True):
        print "Error: simulate_StarStarRV and simulate_StarPlanetRV can NOT BOTH be True!"

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
    if verbose:
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
        phases3.append(phases2)   
    #print "phases3 = "+repr(phases3)
    
    #subtract RVoffsets
    RVsIN = RVs
    if type(RVsIN)!=np.ndarray:
        RVsIN = np.array(RVsIN)
    RVsOUT2 = []
    numEpochsTotal = 0
    for orb in range(0,len(e)):
        RVsOUT = []
        for dataset in range(0,len(RVsIN)):
            offset = RVoffsets[orb][dataset]
            s="offset being subtracted from RV data for orb #"+str(orb)+"/"+str(len(e))+", dataset #"+str(dataset)+" is "+str(offset)+"\n"
            numEpochsTotal+=len(RVsIN[dataset])
            RVsOUT.append(RVsIN[dataset]-offset)
        log.write(s+'\n')
        if verbose:
            print s
        RVsOUT2.append(RVsOUT)
        
    # data for planet's orbit that are needed to calculate proper chiSquare for RV fit
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
        (To,Tcent) = genTools.eccArgPeri2ToTcCalc(planet_e, planet_P, planet_argPeri, planet_T, Tc=0)
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
        (To,Tcent) = genTools.eccArgPeri2ToTcCalc(star_e, star_P, star_argPeri, star_T, Tc=0)
        star_Tc = Tcent
    star_K = sysDataDict["star_K"]
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
    star_Ks = []
    star_argPeris = []
    star_incs = []
    K_use = False
    # e, T, period, inc, argPeri_deg, a,
    if simulate_StarPlanetRV:
        s= 'simulate_StarPlanetRV==True, so loading up inputs as planet variables'
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
            #if verbose:
            #    print 'appending '+str(planet_T)+" to planet_Ts"
            planet_Tcs.append(planet_Tc)
            planet_MsinIs.append(planet_MsinI)
            planet_Ks.append(planet_K)
            planet_incs.append(0)
            
    if simulate_StarStarRV:
        s= 'simulate_StarStarRV==True, so loading up inputs as star variables'
        log.write(s+'\n')
        if verbose:
            print s
        star_es = e
        star_Ts = T
        star_Tcs = Tc
        star_Ps = period
        star_argPeris = argPeri_deg
        star_incs = inc
        star_Ks = K
        K_use = K[0]
        Mass2 = star_Mass2
    else:
        for orb in range(0,len(e)):
            star_es.append(star_e)
            star_Ts.append(star_T)
            star_Tcs.append(star_Tc)
            star_Ps.append(star_P)
            star_argPeris.append(star_argPeri)
            star_Ks.append(star_K)
            star_incs.append(star_inc)
            
    residuals3 = []
    planetVRs3 = []
    starVRs3 = []
    print '\nAbout to calculate RVs and residuals'
    if len(argPeri_deg)>len(colorsList):
        P =ProgressBar('red',width=30,block='=',empty='-',lastblock='>')
    for orb in range(0,len(e)):
        residuals2 = []
        planetVRs2 = []
        starVRs2 = []
        K_s2 = []
        K_p2 = []
        chiSquaredTot2 = 0.0
        numEpochs_RV = 0.0
        chi_squared_RV_reducedCur2 = []
        chi_squared_RV_Cur2 = []
        for dataset in range(0,len(RVsOUT2[orb])):
            chiSquaredTot = 0.0
            epochs = RV_epochsIN2[dataset]
            rvs = RVsOUT2[orb][dataset]
            errors = RV_errors[dataset]                
            residuals = []
            planetVRs = []
            starVRs = []
            (a_total_s, a1_s, a2_s, p_s) = genTools.semiMajorConverter(Mass1, star_Mass2, a_total=a[orb],a1=0.0,a2=0.0, period=period[orb],verbose=False)
            (a_total_p, a1_p, a2_p, p_p) = genTools.semiMajorConverter(Mass1, planet_MsinIs[orb], a_total=a[orb],a1=0.0,a2=0.0, period=period[orb],verbose=False)
            
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
            s=s+ "\nstar_Ks[orb] = "+str(star_Ks[orb])
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
                    a_s = a1_s
                    if primaryRVs==False:
                        a_s = a2_s
                    (v_r_c,K_s) = rvTools.vrCalculatorSemiMajorType(epochs[epoch],star_es[orb],star_Ts[orb],star_Ps[orb],star_argPeris[orb],a_s,T_center=star_Tcs[orb],i=star_incs[orb], K=star_Ks[orb], verbose=False)
                    K_s2.append(K_s)
                # calculate the velocity residual due to the planet around primary
                if (planet_Ps[orb]==0):
                    (v_r_p,K_p)=(0,0)
                else:
                    s=  "\nParams:"
                    s=s+"\nepochs[epoch] = "+str(epochs[epoch])
                    log.write(s+'\n')
                    if False:
                        print s
                    a_p = a1_p
                    if primaryRVs==False:
                        a_p = a2_p
                    (v_r_p,K_p) = rvTools.vrCalculatorSemiMajorType(epochs[epoch],planet_es[orb],planet_Ts[orb],planet_Ps[orb], planet_argPeris[orb],a_p,T_center=planet_Tcs[orb],i=planet_incs[orb], K=planet_Ks[orb], verbose=False)
                    K_p2.append(K_p)
                    #print 'K_p being used is = ',planet_Ks[orb]
                    #print "v_r_p = ",v_r_p
                RV =rvs[epoch]- (v_r_c+v_r_p)
                if ((abs(RV)>80) and False):
                    s= 'Bad point found at '+str(rvs[epoch])+", taken at "+str(RV_epochsIN2[dataset][epoch])+", was off by = "+str(RV)
                    print s
                    log.write(s+'\n')
                residuals.append(RV/MtoKM)
                planetVRs.append(v_r_p/MtoKM) 
                starVRs.append(v_r_c/MtoKM) 
                chiSquaredCurr = genTools.chiSquaredCalc(rvs[epoch], error, v_r_c+v_r_p)
                chiSquaredTot = chiSquaredTot+chiSquaredCurr
                
                s= '\nRV = '+str(rvs[epoch])+' - ('+str(v_r_c)+' + '+str(v_r_p)+') = '+str(RV)
                #print 'epoch = ',epochs[epoch]
                s=s+ '\nchiSquaredCurr = '+str(chiSquaredCurr)
                s=s+ ', chiSquaredTot = '+str(chiSquaredTot)
                if False:
                    log.write(s+'\n')
                if False:
                    print s
            nuRV_cur = nuRV*(float(len(rvs))/float(numEpochsTotal))
            s=s+ "\nFor dataset"+str(dataset)+": nuRV = "+str(nuRV)+", len(rvs) = "+str(len(rvs))+", numEpochsTotal = "+str(numEpochsTotal)+", nuRV_cur = "+str(nuRV_cur)
            chi_squared_RV_reducedCur = (1.0/nuRV_cur)*chiSquaredTot
            chi_squared_RV_Cur2.append(chiSquaredTot)
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
        if len(argPeri_deg)>len(colorsList):
            P.render((orb+1) * 100 // len(argPeri_deg), 'Complete so far.')
    if simulate_StarStarRV:
        KcalcUse = K_s2[0]
    else:
        KcalcUse = K_p2[0]
        
    ## Make an updated Radial Velocities, and errors, for data accounting for planet or star RVs
    RVsOUTupdated3 = []
    RVerrorsUpdated3 = []
    for orb in range(0,len(e)):
        RVsOUTupdated2 = []
        RVerrorsUpdated2 = []
        for dataset in range(0,len(RVsOUT)):
            RVsOUTupdated = []
            RVerrorsUpdated = []
            for epoch in range(0,len(RVsOUT[dataset])): 
                RVsOUTupdated.append((RVsOUT[dataset][epoch]/MtoKM))
                RVerrorsUpdated.append(RV_errors[dataset][epoch]/MtoKM)
            RVsOUTupdated2.append(RVsOUTupdated)
            RVerrorsUpdated2.append(RVerrorsUpdated)
        RVsOUTupdated3.append(RVsOUTupdated2)
        RVerrorsUpdated3.append(RVerrorsUpdated2)
    # if simulating star's orbit and planet RVs exist sub them
    if (simulate_StarStarRV and (planetVRs3[0][0][0]!=0)):
        s= '\n Subtracting planet VRs from data!\n'
        if verbose:
            print s
        log.write(s+'\n')
        for orb in range(0,len(e)):
            for dataset in range(0,len(RVsOUT)):
                for epoch in range(0,len(RVsOUT[dataset])): 
                    #print 'RVsOUT[dataset][epoch]/MtoKM = '+str(RVsOUT[dataset][epoch]/MtoKM)+", planetVRs3[orb][dataset][epoch] = "+str(planetVRs3[orb][dataset][epoch])
                    RVsOUTupdated3[orb][dataset][epoch]+=-planetVRs3[orb][dataset][epoch]
    # if simulating planet's orbit and star RVs exist sub them
    elif (simulate_StarPlanetRV  and (starVRs3[0][0][0]!=0)):
        s= '\n Subtracting star VRs from data!\n'
        if verbose:
            print s
        log.write(s+'\n')
        for orb in range(0,len(e)):
            for dataset in range(0,len(RVsOUT)):
                for epoch in range(0,len(RVsOUT[dataset])): 
                    #print 'RVsOUT[dataset][epoch]/MtoKM = '+str(RVsOUT[dataset][epoch]/MtoKM)+", starVRs3[orb][dataset][epoch] = "+str(starVRs3[orb][dataset][epoch])
                    RVsOUTupdated3[orb][dataset][epoch]+=-starVRs3[orb][dataset][epoch]
        
    #load up strings to print out RVs calculate for planet or companion star
    if planetVRs3[0][0][0]!=0:
        s = '\n'+'*'*50+'\nPlanet Radial Velocities:\n'
        for orb in range(0,len(e)):
            for dataset in range(0,len(RVsOUT)):
                s=s+'\nFor dataset # '+str(dataset)+':\n'
                for epoch in range(0,len(RVsOUT[dataset])): 
                    if planetVRs3[orb][dataset][epoch]!=0:
                        s = s+str(planetVRs3[orb][dataset][epoch])+'\n'
        s =s+ '\n'+'*'*50
        if False:
            log.write(s+'\n')
        if verbose:
            print s
        
    if starVRs3[0][0][0]!=0:
        s = '\n'+'*'*50+'\nStar Radial Velocities:\n'
        for orb in range(0,len(e)):
            for dataset in range(0,len(RVsOUT)):
                s=s+'\nFor dataset # '+str(dataset)+':\n'
                for epoch in range(0,len(RVsOUT[dataset])): 
                    if starVRs3[orb][dataset][epoch]!=0:
                        s = s+str(starVRs3[orb][dataset][epoch])+'\n'
        if False:
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
    trendVRs2 = []
    trendTimes2 = []
    orbitPhases2 = []
    times2 = []
    s= '\nabout to calc fit line'
    if verbose:
        print s
    print s
    log.write(s+'\n')
    if len(argPeri_deg)>len(colorsList):
        P =ProgressBar('red',width=30,block='=',empty='-',lastblock='>')
    for orb in range(0,len(e)):
        orbitVRs = []
        numSteps = 1000.0
        periodIncrement = (period[orb]*365.242)/numSteps
        t = Tc[orb]-((period[orb]*365.242)/2.0)
        times = []
        trendVRs = np.zeros(numSteps)
        trendTimes = np.zeros(numSteps)
        (a_total, a1, a2, p) = genTools.semiMajorConverter(Mass1, Mass2, a_total=a[orb],a1=0.0,a2=0.0, period=period[orb])
        for step in range(0,int(numSteps*1.0)):
            t = t + periodIncrement
            times.append(t)
            trendTimes[step]=t+(period[orb]*365.242/2)
            aUse = a1
            if primaryRVs==False:
                aUse = a2
            # calculate the velocity residual due to the companion 
            (v_r_c,Kcalc) = rvTools.vrCalculatorSemiMajorType(t,e[orb],T[orb],period[orb],argPeri_deg[orb],aUse,T_center=Tc[orb],i=inc[orb], K=K_use, verbose=False)
            orbitVRs.append(v_r_c/MtoKM)
            if step<(numSteps//2):
                trendVRs[step+(numSteps//2)]=v_r_c/MtoKM
            else:
                trendVRs[step-(numSteps//2)]=v_r_c/MtoKM
        if len(argPeri_deg)<len(colorsList):
            s= 'Orbit '+str(orb)+" had a K = "+str(Kcalc)
            if verbose:
                print s
            log.write(s+'\n')
        trendTimes2.append(trendTimes)
        trendVRs2.append(trendVRs)
#         if orb==1:
#             for i in range(0,len(orbitVRs)):
#                 print str(trendTimes[i])+',  '+str(trendVRs[i])
        orbitPhases = epochsToPhases(times,Tc[orb],period[orb], verbose=False)
        #print 'orbital phases were:\n'+repr(orbitPhases)
        orbitPhases2.append(orbitPhases)
        orbitVRs2.append(orbitVRs)
        if len(argPeri_deg)>len(colorsList):
            P.render((orb+1) * 100 // len(argPeri_deg), 'Complete so far.')
        #print 'phases: ',repr(orbitPhases)
        #print 'VRs: ',repr(orbitVRs)
    
    ######################################################################################################################################
    ## Start Plotting!!
    ######################################################################################################################################
        
    ############################################################
    ## Create figure, axes and start plotting/drawing everything
    fig = plt.figure(1,figsize=(10,5))
    #plt.suptitle(plotFileTitle, fontsize=10)
    residualsPlot = fig.add_subplot(212)
    residualsPlot.set_position([0.12,0.15,0.83,0.23])
    #residualsPlot.set_title("Residuals Plot")
    residualsPlot.axes.set_xlabel("Orbital Phase",fontsize=20)
    residualsPlot.axes.set_ylabel("Residual",fontsize=15)
    residualsPlot.tick_params(axis='y',which='major',width=3,length=5,pad=10,direction='in',labelsize=8)
    plt.locator_params(axis='y',nbins=5)
    #colorsList = ['b','m','k','g','y','o','p']
    
    #plot RESIDUAL data and fit, plus build up the chiSquaredStr
    print '\nabout to load up residuals plot'
    # create a progress bar object for updates to user of progress
    if len(argPeri_deg)>len(colorsList):
        P =ProgressBar('red',width=30,block='=',empty='-',lastblock='>')
    for orb in range(0,len(argPeri_deg)):
        chiSquaredStr = ''
        c='blue'
        if len(argPeri_deg)<len(colorsList):
            c = colorsList[orb]
        for dataset in range(0,len(RVs)):
            if len(argPeri_deg)<len(colorsList):
                chiSquaredStr = chiSquaredStr+'\ndataset '+str(dataset)+' has chiSquared = '+str(chi_squared_RV_Cur2[dataset])+", or reduced chiSquared = "+str(chi_squared_RV_reducedCur2[dataset])
                chiSquaredStr = chiSquaredStr+'\nchiSquared_reduced for dataset '+str(dataset)+' is = '+str(chi_squared_RV_reducedCur2[dataset])
            s= '\n\n*** RMS of residuals '+str(dataset)+'= '+str(np.sqrt(np.mean(np.array(residuals3[orb][dataset])**2)))+' ***\n\n'
            log.write(s+'\n')
            if verbose:
                print s
            residualsPlot = addRVdataToPlot(residualsPlot,residuals3[orb][dataset],RVerrorsUpdated3[orb][dataset],phases3[orb][dataset],alf=alf,color=c,plotErrorBars=plotErrorBars)
        if len(argPeri_deg)>len(colorsList):
            P.render((orb+1) * 100 // len(argPeri_deg), 'Complete so far.')
               
    chiSquaredStr = chiSquaredStr+'\nchiSquared_reduced for ALL data is = '+str(chi_squared_RV_reduced)
    s= chiSquaredStr
    if verbose:
        print s
    log.write(s+'\n')
    
    xmin = genTools.findArrayMin(phases3)
    xmax = genTools.findArrayMax(phases3)
    ymin = genTools.findArrayMin(residuals3)
    ymax = genTools.findArrayMax(residuals3)

    xLim =[xmin,xmax]
    xrange = xmax-xmin
    xLim2 = (xmin-abs(xrange*0.05), xmax+abs(xrange*0.05))
    yLim = [ymin,ymax]
    yrange = ymax-ymin
    yLim2 = (ymin-abs(yrange*0.05), ymax+abs(yrange*0.05))
    
    ## make plot of fit to data
    fitPlot = fig.add_subplot(211)
    fitPlot.set_position([0.12,0.38,0.83,0.55])
    fitPlot.xaxis.set_ticklabels([])#this is just a hack way of killing the tick labels
    fitXmin = genTools.findArrayMin(orbitPhases2)
    if xmin<fitXmin:
        fitXmin = xmin
    fitXmax = genTools.findArrayMax(orbitPhases2)
    if xmax>fitXmax:
        fitXmax = xmax
    fitYmin = genTools.findArrayMin(orbitVRs2)
    if ymin<fitYmin:
        fitYmin = ymin
    fitYmax = genTools.findArrayMax(orbitVRs2)
    if ymax>fitYmax:
        fitYmax = ymax
        
    yMaxRVs = genTools.findArrayMax(RVsOUTupdated3)
    yMinRVs = genTools.findArrayMin(RVsOUTupdated3)
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
    yLabel = "RV [m/s]"
    if KMperSec:
        yLabel = "RV [km/s]"
    fitPlot.axes.set_ylabel(yLabel,fontsize=20)
    residualsPlot=fixPlotBordersAndLabels(residualsPlot)
    #print 'ln3712: np.max(RVsOUTupdated3) = '+str(np.max(RVsOUTupdated3))
    print "\nAbout to load up fit plot"
    if len(argPeri_deg)>len(colorsList):
        P =ProgressBar('red',width=30,block='=',empty='-',lastblock='>')
    for orb in range(0,len(e)):
        fitPlot.plot(orbitPhases2[orb],orbitVRs2[orb],c='r',linewidth=2.0,alpha=alf)
        c='blue'
        if len(argPeri_deg)<len(colorsList):
            c = colorsList[orb]
        for dataset in range(0,len(RVs)):
            fitPlot = addRVdataToPlot(fitPlot,RVsOUTupdated3[orb][dataset],RVerrorsUpdated3[orb][dataset],phases3[orb][dataset],alf=alf,color=c,plotErrorBars=plotErrorBars)
        if len(argPeri_deg)>len(colorsList):
            P.render((orb+1) * 100 // len(argPeri_deg), 'Complete so far.')
    paramsLegndStr = 'e = '+str(e[0])
    paramsLegndStr +='\ninc = '+str(inc[0])
    paramsLegndStr +='\na = '+str(a[0])
    paramsLegndStr +='\nargPeri_deg = '+str(argPeri_deg[0])
    if period[0]>4.0:
        paramsLegndStr +='\nperiod [Yrs] = '+str(period[0])
    else:
        paramsLegndStr +='\nperiod [days] = '+str(period[0]*365.242)
    paramsLegndStr +='\nTo = '+str(T[0])
    paramsLegndStr +=", Tc = "+str(Tc[0])
    paramsLegndStr+="\nK [m/s] = "+str(KcalcUse)
    for dataset in range(0,len(RVs)):
        if dataset==0:
            paramsLegndStr+='\n'
        else:
            paramsLegndStr+=', '
        paramsLegndStr += 'offset '+str(dataset)+'= '+str(RVoffsets[dataset])
    paramsLegndStr += chiSquaredStr
    if verbose:
        print paramsLegndStr
    log.write(paramsLegndStr+'\n')
    
    ## code to add a legend to the plot, but might need to be tweaked depending on your plot and param values.
    if plotFullOrbit and addLegend:
        # in bottom right
        #residualsPlot.text(fitXLimsUSE[1]-abs(fitXLimsUSE[1]*0.002),abs(fitYLimsUSE[1]*0.5),paramsLegndStr,ha='left')
        # in bottom left
        paramsLegndY = fitYLimsUSE[0]+abs(fitYLimsUSE[0]*0.02)
        paramsLegndX = fitXLimsUSE[0]+abs(fitXLimsUSE[0]*0.03)
        fitPlot.text(paramsLegndX,paramsLegndY,paramsLegndStr,ha='left')
        s= '\nlegend bottom left corner at [ '+str(paramsLegndX)+' , '+str(paramsLegndY)+' ]'
        if verbose:
            print s
        log.write(s+'\n')
        
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
    if makeTrendPlot:       
        print "\nabout to load up trend plot"
        if len(argPeri_deg)>len(colorsList):
            P =ProgressBar('red',width=30,block='=',empty='-',lastblock='>')
        fig2 = plt.figure(2,figsize=(10,5))
        #plt.suptitle(plotFileTitle, fontsize=10)
        residualsPlotTrend = fig2.add_subplot(212)
        residualsPlotTrend.set_position([0.12,0.15,0.83,0.23])
        RVsPlotTrend = fig2.add_subplot(211)
        RVsPlotTrend.set_position([0.12,0.38,0.83,0.55])
        RVsPlotTrend.xaxis.set_ticklabels([])#this is just a hack way of killing the tick labels
        #RVsPlotTrend.set_title("RVs TREND")
        RVsPlotTrend.axes.set_ylabel(yLabel,fontsize=20)
        #residualsPlotTrend.set_title("Residuals Plot TREND")
        residualsPlotTrend.axes.set_ylabel("Residual",fontsize=20)
        residualsPlotTrend.axes.set_xlabel("Epoch [JD]",fontsize=20)
        residualsPlotTrend=fixPlotBordersAndLabels(residualsPlotTrend)
        RVsPlotTrend=fixPlotBordersAndLabels(RVsPlotTrend)
        if False:
            rng = genTools.findArrayMax(RV_epochsIN2)-genTools.findArrayMin(RV_epochsIN2)
            rngUse= [genTools.findArrayMin(RV_epochsIN2)-rng*0.05,genTools.findArrayMax(RV_epochsIN2)+rng*0.05]
        else:
            rng = genTools.findArrayMax(trendTimes2)-genTools.findArrayMin(trendTimes2)
            rngUse= [genTools.findArrayMin(trendTimes2)-rng*0.05,genTools.findArrayMax(trendTimes2)+rng*0.05]
        for orb in range(0,len(argPeri_deg)):
            RVsPlotTrend.plot(trendTimes2[orb],trendVRs2[orb],c='r',linewidth=2.0,alpha=alf)
            c='blue'
            if len(argPeri_deg)<len(colorsList):
                c = colorsList[orb]
            for dataset in range(0,len(RVs)):
                residualsPlotTrend = addRVdataToPlot(residualsPlotTrend,residuals3[orb][dataset],RVerrorsUpdated3[orb][dataset],RV_epochsIN2[dataset],alf=alf,color=c,plotErrorBars=plotErrorBars)
                RVsPlotTrend = addRVdataToPlot(RVsPlotTrend,RVsOUTupdated3[orb][dataset],RVerrorsUpdated3[orb][dataset],RV_epochsIN2[dataset],alf=alf,color=c,plotErrorBars=plotErrorBars)
            if len(argPeri_deg)>len(colorsList):
                P.render((orb+1) * 100 // len(argPeri_deg), 'Complete so far.')
        residualsPlotTrend.plot(rngUse,[0,0],c='r',linewidth=2.0)
        residualsPlotTrend.axes.set_ylim(yLim2)
        RVsPlotTrend.plot(rngUse,[0,0],c='r',linewidth=2.0)
        RVsPlotTrend.axes.set_ylim(fitYLimsUSE)        
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
    if makeDataResidualPlot:
        residuals3_trimmed = []
        for orb in range(0,len(argPeri_deg)):
            residuals2_trimmed = []
            for dataset in range(0,len(RVs)):
                residual_trimmed = []
                for epoch in range(0,len(residuals3[orb][dataset])):
                    residual_trimmed.append(residuals3[orb][dataset][epoch]/RVerrorsUpdated3[orb][dataset][epoch])
            residuals2_trimmed.append(residual_trimmed)
        residuals3_trimmed.append(residuals2_trimmed)
        #print "len(RV_epochsIN2) = "+str(RV_epochsIN2)+"\n"+"len(residuals3) = "+str(residuals3)+"\n"+"len(residual_trimmed) = "+str(residual_trimmed)+"\n"+"len(residuals2_trimmed) = "+str(residuals2_trimmed)+"\n"+"len(residuals3_trimmed) = "+str(residuals3_trimmed)+"\n"   #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
        fig = plt.figure(1,figsize=(10,5))
        #main = fig.add_subplot(111)
        fig.text(0.05,0.05,paramsLegndStr,ha='left')
        if plotFilename!='':
            legendFigFilename = plotFilename[:-4]+"-paramInfo.png"
            plt.savefig(legendFigFilename, dpi=300, orientation='landscape')
        plt.close()
    log.write('\n'+75*'#'+'\n Leaving rvPlotterDuo \n'+75*'#'+'\n') 
    log.close()
    
    return chiSquaredStr

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
    simulate_StarStarRV = paramSettingsDict["simulate_StarStarRV"]
    simulate_StarPlanetRV = paramSettingsDict["simulate_StarPlanetRV"]
    if (simulate_StarStarRV is True) and(simulate_StarPlanetRV is True):
        print "Error: simulate_StarStarRV and simulate_StarPlanetRV can NOT BOTH be True!"
    if simulate_StarPlanetRV:
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
    titleTop = os.path.dirname(plotFilename).split('/')[-1]
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
    if simulate_StarPlanetRV:
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
            
    if simulate_StarStarRV:
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
            (a_total_s, a1_s, a2_s, p_s) = genTools.semiMajorConverter(Mass1, star_Mass2, a_total=a[orb],a1=0.0,a2=0.0, period=period[orb],verbose=False)
            (a_total_p, a1_p, a2_p, p_p) = genTools.semiMajorConverter(Mass1, planet_MsinIs[orb], a_total=a[orb],a1=0.0,a2=0.0, period=period[orb],verbose=False)
            
            if False:
                print '\nFor dataset '+str(dataset)
            for epoch in range(0,len(epochs)):
                if False:
                    print '\nWorking on epoch ',epoch
                # calculate the velocity residual due to the companion star
                if (star_Ps[orb]==0):
                    (v_r_c,K_s)=(0,0)
                else:
                    (v_r_c,K_s) = rvTools.vrCalculatorSemiMajorType(epochs[epoch],star_es[orb],T[orb],star_Ps[orb],star_argPeris[orb],a1_s,T_center=T_center,i=star_incs[orb], K=False, verbose=False)
                # calculate the velocity residual due to the planet around primary
                if (planet_Ps[orb]==0):
                    (v_r_p,K_p)=(0,0)
                else:
                    (v_r_p,K_p) = rvTools.vrCalculatorSemiMajorType(epochs[epoch],planet_es[orb],planet_Ts[orb],planet_Ps[orb], planet_argPeris[orb],a1_p,T_center=planet_Tcs[orb],i=planet_incs[orb], K=planet_Ks[orb], verbose=False)
                    #print 'K_p being used is = ',planet_Ks[orb]
                    #print "v_r_p = ",v_r_p
                RV =rvs[epoch]- (v_r_c+v_r_p)
                if (abs(RV)>100):
                    print 'Bad point found at '+str(rvs[epoch])+", "+str(RV_epochsIN2[dataset][epoch])+", was off by = "+str(RV)
                residuals.append(RV)
                if False:
                    print 'Appending RV '+str(RV)+' to the residuals array'
                RVsCalcd.append(v_r_c+v_r_p)
                chiSquaredCurr = genTools.chiSquaredCalc(rvs[epoch], errors[epoch], v_r_c+v_r_p)
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
            (residuals_random, errors_random) = genTools.listRandomizer(residuals3[0][dataset],RV_errors[dataset])
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
        (a_total, a1, a2, p) = genTools.semiMajorConverter(Mass1, Mass2, a_total=a[orb],a1=0.0,a2=0.0, period=period[orb])
        for step in range(0,int(numSteps)):
            t = t + periodIncrement
            times.append(t)
            # calculate the velocity residual due to the companion 
            (v_r_c,K) = rvTools.vrCalculatorSemiMajorType(t,e[orb],T[orb],period[orb],argPeri_deg[orb],a1,T_center=T_center,i=inc[orb], K=planet_K, verbose=False)
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
                print '\nepoch = ',epoch
                print 'period [days] = ',P_days
                print 'phase = ',phase
                if True:
                    print str(epoch)+" - "+str(int((epoch-T_center)/P_days)*P_days)+" - "+str(T_center)+" = "+str(phaseTimeDiff)
        phases2.append(phases)
    
    if notAlist:
        epochs2 = epochs2[0]
        phases2 = phases2[0]
    
    if verbose:
        print '\nepochs2 array = '+repr(epochs2)
        print 'phases2 array = '+repr(phases2)
        
    return phases2
             
