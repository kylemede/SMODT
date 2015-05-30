#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import numpy as np
import os
import pylab
plt = pylab.matplotlib.pyplot
patches = pylab.matplotlib.patches
import constants as const
import generalTools as genTools
import smodtLogger

log = smodtLogger.getLogger('main.plotTools',lvl=100,addFH=False)

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
    if False:
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
        subPlot.plot(epochs[epoch],RVs[epoch],c='k',marker='.',markersize=6)
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

def summaryPlotter(outputDataFilename, plotFilename,stage='MCMC', shadeConfLevels=True):
    """
    This advanced plotting function will plot all the data in a 3x4 grid on a single figure.  The data will be plotted
    in histograms that will be normalized to a max of 1.0.  The 
    confidence levels of the data can be calculated and the resulting histograms will have the bars inside the
    68% bin shaded as dark grey and 95% as light grey and the rest will be white.
    
    """
    quiet = True
    verbose = False
    latex=True
    plotFormat = 'eps'
    forceRecalc = False
    
    (head,data) = genTools.loadFits(outputDataFilename)
    
    #plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    plt.rc('font',family='serif')
    if latex:
        plt.rc('text', usetex=True)
        plt.rcParams['text.latex.unicode']=True  
    
    ## find number of RV datasets
    if head!=False:  
        log.debug(' Inside summaryPlotter')
         
        s= '\nCreating summary plot for file:\n'+outputDataFilename
        s=s+ '\nInput plotfilename:\n'+plotFilename
        log.info(s)
        
        # check if the passed in value for plotFilename includes format extension
        if '.'+plotFormat not in plotFilename:
            plotFilename = plotFilename+"."+plotFormat
            log.debug('updating plotFilename to:\n'+plotFilename)
        else:
            plotFilename = plotFilename
        
        paramList = genTools.getParInts(head)    
        #print 'paramList = '+repr(paramList) 
        paramFileStrs = ['M1','M2','dist','Omega','e','To', 'Tc','P','i','omega','a_total','chiSquared']
        paramStrs = ['M1 [Msun]','M2 [Msun]','Distance [PC]','Omega [deg]','e','To [JD]', 'Tc [JD]','P [Yrs]','i [deg]','omega [deg]','a_total [AU]','chiSquared']
        if latex:
            paramStrs = ['$M_1$ [$M_{sun}$]','$M_2$ [$M_{sun}$]','$Distance$ [PC]','$\Omega$ [deg]','$e$','$T_o$ [JD]', '$T_c$ [JD]','$P$ [Yrs]','$i$ [deg]','$\omega$ [deg]','$a_{total}$ [AU]','chiSquared']

        #perfectVals = [70.0,0.5,2457000.0,2457000.0,5.0,40.0,50.0,3.34718746581]
        if head["nRVdsets"]>0:
            paramFileStrs.append('K')
            if latex:
                paramStrs.append('$K$ [m/s]')
            else:
                paramStrs.append('K [m/s]')
            #perfectVals.append(4933.18419284)
            for dataset in range(1,head["nRVdsets"]+1):
                paramFileStrs.append('offset_'+str(dataset))
                if latex:
                    paramStrs.append("$offset_{"+str(dataset)+"}$ [m/s]")
                else:
                    paramStrs.append('offset '+str(dataset)+' [m/s]')
                #perfectVals.append(0.0)
        
        ## run through all the data files and parameters requested and make histogram files
        for i in range(0,len(paramList)):
            if (os.path.exists(os.path.join(os.path.dirname(outputDataFilename),'hist-'+stage+"-"+str(paramList[i])+'.dat'))==False)or forceRecalc:
                if verbose:
                    print 'Initial Plotting for parameter '+str(i+1)+"/"+str(len(paramList))+": "+paramStrs[i]+", for file:\n"+outputDataFilename
                histDataBaseName = os.path.join(os.path.dirname(outputDataFilename),'hist-'+stage+"-"+str(paramList[i]))
                (CLevels,data,bestDataVal) = genTools.confLevelFinder(outputDataFilename,paramList[i], returnData=True, returnChiSquareds=False, returnBestDataVal=True,fast=False)
                histMakeAndDump([],data,outFilename=histDataBaseName,weight=False, normed=False, nu=1,logY=False,histType='step')
                if (os.path.exists(os.path.join(os.path.dirname(outputDataFilename),'confLevels-'+stage+"-"+str(paramList[i])+'.dat'))==False)or forceRecalc:
                    np.savetxt(os.path.join(os.path.dirname(outputDataFilename),'confLevels-'+stage+"-"+str(paramList[i])+'.dat'),CLevels)
                    if verbose:
                        print 'confidence levels data stored to:\n'+os.path.join(os.path.dirname(outputDataFilename),'confLevels-'+stage+"-"+str(paramList[i])+'.dat')
        
        # Create empty figure to be filled up with plots
        fig = plt.figure(figsize=(10,10)) 
                    
        ## make combined/stacked plot for each parameter in list
        for i in range(0,len(paramList)):
            s='\nStarting to plot shaded hist for '+paramStrs[i]
            if verbose:
                print s
            subPlot = fig.add_subplot(3,4,i+1)
            
            histDataBaseName = os.path.join(os.path.dirname(outputDataFilename),'hist-'+stage+"-"+str(paramList[i]))
            if quiet==False:
                print 'Loading and re-plotting parameter '+str(i+1)+"/"+str(len(paramList))+": "+paramStrs[paramList[i]]#+" for file:\n"+outputDataFilename
            xLim=False
            CLevels=False
            if shadeConfLevels:
                CLevels=np.loadtxt(os.path.join(os.path.dirname(outputDataFilename),'confLevels-'+stage+"-"+str(paramList[i])+'.dat'))
            showYlabel=False
            if (i==0)or(i==4)or(i==8):
                showYlabel = True
            subPlot = histLoadAndPlot_ShadedPosteriors(subPlot,outFilename=histDataBaseName,confLevels=CLevels,xLabel=paramStrs[paramList[i]],xLims=xLim,latex=latex,showYlabel=showYlabel)         
            
        ## Save file if requested.
        if verbose:
            if quiet==False:
                print '\nStarting to save param hist figure:'
        if plotFilename!='':
            plt.savefig(plotFilename,format=plotFormat)
            plt.close()
            s= 'Summary plot saved to: '+plotFilename
            if quiet==False:
                print s
            log.info(s)
        if True:
            if quiet==False:
                print 'converting to PDF as well'
            os.system("epstopdf "+plotFilename)
            
def star(R, x0, y0, color='w', N=5, thin = 0.5):
    """
    Returns an N-pointed star of size R at (x0, y0) (matplotlib patch).
    NOTE: code base taken from star_patch.py in 'Beginning-Python-Visualization'
    """
    from pylab import Polygon
    polystar = np.zeros((2*N, 2))
    for i in range(2*N):
        angle = i*const.pi/N
        r = R*(1-thin*(i%2))
        polystar[i] = [r*np.cos(angle)+x0, r*np.sin(angle)+y0]
    return Polygon(polystar, fc=color, ec='black',linewidth=1.5)  

def epochsToPhases(epochs2,Tc,P_yrs, halfOrbit=False):
    """
    converter to change the epochs of a 1-D or 2-D array to phases.
    Output dimensions of phases2 will be same as input epochs2
    halfOrbit option shifts phases into [-0.5,0.5] range, 
    rather than [-1,1].
    """    
    verbose=False 
    notAlist = False   
    if type(epochs2[0])!=list:
        epochs2 = [epochs2]
        notAlist = True        
        
    phases2 = []
    P_days = P_yrs*const.daysPerYear
    for epochs in epochs2:
        phases = []
        for epoch in epochs:
            phaseTimeDiff = epoch - int((epoch-Tc)/P_days)*P_days-Tc #phase shifted into [-P,P]
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
                    print str(epoch)+" - "+str(int((epoch-Tc)/P_days)*P_days)+" - "+str(Tc)+" = "+str(phaseTimeDiff)
        phases2.append(phases)
    
    if notAlist:
        epochs2 = epochs2[0]
        phases2 = phases2[0]
    if verbose:
        print '\nepochs2 array = '+repr(epochs2)
        print 'phases2 array = '+repr(phases2)
        
    return phases2

