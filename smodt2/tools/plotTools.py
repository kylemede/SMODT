#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import numpy as np
import os
import pylab
import copy
plt = pylab.matplotlib.pyplot
patches = pylab.matplotlib.patches
import constants as const
import generalTools as genTools
import cppTools
import smodtLogger

log = smodtLogger.getLogger('main.plotTools',lvl=100,addFH=False)  

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

def addRVdataToPlot(subPlot,realData,alf=1.0,color='blue',plotErrorBars=False):
    """
    Add '+' markers for the data locations with respective y axis errors 
    shown as the height of the markers. 
    """
    for i in range(0,realData.shape[0]):
        xs = [realData[i,0],realData[i,0]]
        ys = [realData[i,5]-realData[i,6],realData[i,5]+realData[i,6]]
        if plotErrorBars:
            subPlot.plot(xs,ys,c=color,linewidth=2,alpha=alf)
        subPlot.plot(realData[i,0],realData[i,5],c='k',marker='.',markersize=6)
    return subPlot

def addDIdataToPlot(subPlot,realData,asConversion):
    """
    To plot a '+' for each data point with width and height matching the errors converted 
    to x,y coords.
    """
    xmin = np.min(realData[:,1]-realData[:,2])*asConversion
    xmax = np.max(realData[:,1]+realData[:,2])*asConversion
    ymin = np.min(realData[:,3]-realData[:,4])*asConversion
    ymax = np.max(realData[:,3]+realData[:,4])*asConversion
    for i in range(0,realData.shape[0]):
        xCent = realData[i,1]*asConversion
        yCent = realData[i,3]*asConversion
        left = xCent-realData[i,2]*asConversion
        right = xCent+realData[i,2]*asConversion
        top = yCent+realData[i,4]*asConversion
        btm = yCent-realData[i,4]*asConversion
        subPlot.plot([left,right],[yCent,yCent],linewidth=3,color='k',alpha=1.0)
        subPlot.plot([xCent,xCent],[btm,top],linewidth=3,color='k',alpha=1.0)
    return (subPlot,[xmin,xmax,ymin,ymax])

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
    forceRecalc = True
    
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
            if (os.path.exists(os.path.join(os.path.dirname(outputDataFilename),'hist-'+stage+"-param"+str(paramList[i])+'.dat'))==False)or forceRecalc:
                if verbose:
                    print 'Initial Plotting for parameter '+str(i+1)+"/"+str(len(paramList))+": "+paramStrs[i]+", for file:\n"+outputDataFilename
                histDataBaseName = os.path.join(os.path.dirname(outputDataFilename),'hist-'+stage+"-param"+str(paramList[i]))
                (CLevels,data,bestDataVal) = genTools.confLevelFinder(outputDataFilename,paramList[i], returnData=True, returnChiSquareds=False, returnBestDataVal=True,fast=False)
                histMakeAndDump([],data,outFilename=histDataBaseName,weight=False, normed=False, nu=1,logY=False,histType='step')
                if (os.path.exists(os.path.join(os.path.dirname(outputDataFilename),'confLevels-'+stage+"-param"+str(paramList[i])+'.dat'))==False)or forceRecalc:
                    np.savetxt(os.path.join(os.path.dirname(outputDataFilename),'confLevels-'+stage+"-param"+str(paramList[i])+'.dat'),CLevels)
                    if verbose:
                        print 'confidence levels data stored to:\n'+os.path.join(os.path.dirname(outputDataFilename),'confLevels-'+stage+"-param"+str(paramList[i])+'.dat')
        
        # Create empty figure to be filled up with plots
        sumFig = plt.figure(figsize=(10,10)) 
                    
        ## make combined/stacked plot for each parameter in list
        for i in range(0,len(paramList)):
            s='\nStarting to plot shaded hist for '+paramStrs[i]
            if verbose:
                print s
            subPlot = sumFig.add_subplot(3,4,i+1)
            
            histDataBaseName = os.path.join(os.path.dirname(outputDataFilename),'hist-'+stage+"-param"+str(paramList[i]))
            if quiet==False:
                print 'Loading and re-plotting parameter '+str(i+1)+"/"+str(len(paramList))+": "+paramStrs[paramList[i]]#+" for file:\n"+outputDataFilename
            xLim=False
            CLevels=False
            if shadeConfLevels:
                CLevels=np.loadtxt(os.path.join(os.path.dirname(outputDataFilename),'confLevels-'+stage+"-param"+str(paramList[i])+'.dat'))
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
            s= 'Summary plot saved to: '+plotFilename
            if quiet==False:
                print s
            log.info(s)
        plt.close()
        if True:
            if quiet==False:
                print 'converting to PDF as well'
            try:
                os.system("epstopdf "+plotFilename)
            except:
                log.warning("Seems epstopdf failed.  Check if it is installed properly.")
            
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

def epochsToPhases(epochs,Tc,P_yrs, halfOrbit=False):
    """
    Convert the epochs (from a realData ary) into phase values 
    (ratio of how far from Tc it is), shifted to lie inside [-1,1].
    if 'halfOrbit'=True, the vals will lie inside [-0.5,0.5].
    """    
    verbose=False         
    phases = []
    P_days = P_yrs*const.daysPerYear
    for epoch in epochs:
        phaseTimeDiff = epoch - int((epoch-Tc)/P_days)*P_days-Tc #phase shifted into [-P,P]
        if verbose:
            print str(epoch)+" - "+str(int((epoch-Tc)/P_days)*P_days)+" - "+str(Tc)+" = "+str(phaseTimeDiff)
        phase = phaseTimeDiff/P_days#phase shifted into [-1,1]
        if halfOrbit:
            if phase>0.5:
                phase = phase-1.0#phase shifted into [-0.5,0.5]
            elif phase<-0.5:
                phase = phase+1.0#phase shifted into [-0.5,0.5]
        phases.append(phase)
        if verbose:
            print '\nepoch = ',epoch
            print 'period [days] = ',P_days
            print 'phase = ',phase        
    return phases

def orbitPlotter(orbParams,settingsDict,plotFnameBase=""):
    """
    Make both the DI and RV plots.
    '-DI.png' and/or '-RV.png' will be added to end of plotFnameBase 
    to make the filenames for each type of plot.
    """
    log.debug("Starting to make orbit plots")
    colorsList =['Blue','BlueViolet','Chartreuse','Fuchsia','Crimson','Aqua','Gold','DarkCyan','OrangeRed','Plum','DarkGreen','Chocolate','SteelBlue ','Teal','Salmon','Brown']
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    #plt.rc('font',family='serif')
    plt.rc('text', usetex=False) 
    
    ##get the real data
    realData = genTools.loadRealData(os.path.join(settingsDict['settingsDir'],settingsDict['prepend']),dataMode=settingsDict['dataMode'])
    
    ## Make Orbit cpp obj
    Orbit = cppTools.Orbit()
    Orbit.loadStaticVars(settingsDict['omegaFdi'][0],settingsDict['omegaFrv'][0])
    Orbit.loadConstants(const.Grav,const.pi,const.KGperMsun, const.daysPerYear,const.secPerYear,const.MperAU)
    
    ################
    # Make DI plot #
    ################
    if settingsDict['dataMode'][0]!='RV':
        ##Make model data for 100~1000 points for plotting fit
        nPts = 500
        fakeRealData = np.ones((nPts,8),dtype=np.dtype('d'),order='C')
        for i in range(0,nPts-1):
            fakeRealData[i,0] = orbParams[5]+(const.daysPerYear*orbParams[7]*(i/float(nPts)))
        fakeRealData[nPts-1,0]  = fakeRealData[0,0]+const.daysPerYear*orbParams[7]
        #print 'fakeRealData[:,0] = '+repr(fakeRealData[:,0])
        Orbit.loadRealData(fakeRealData)
        fitDataDI = np.ones((nPts,3),dtype=np.dtype('d'),order='C')
        orbParamsDI = copy.deepcopy(orbParams)
        orbParamsDI[9]+=settingsDict['omegaFdi'][0]
        #Ensuring that params are in required format for SWIG
        paramsDI = []
        for par in orbParamsDI:
            paramsDI.append(par)
        paramsDI=np.array(paramsDI,dtype=np.dtype('d'),order='C')
        Orbit.calculate(fitDataDI,paramsDI)
        if False:
            ## Get 1/4 locations (useful for drawing semi-major axis, and finding loc of COM)
            fakeRealDataQuarter = np.ones((4,8),dtype=np.dtype('d'))
            for i in range(0,4):
                fakeRealDataQuarter[i,0] = orbParams[5]+(const.daysPerYear*orbParams[7]*(i/4.0))
            Orbit.loadRealData(fakeRealDataQuarter)
            fitDataQuarter = np.zeros((4,3),dtype=np.dtype('d'))
            Orbit.calculate(fitDataQuarter,orbParams)
            #find loc of COM for possible use
            xCOM = (fakeRealDataQuarter[3,0]+fakeRealDataQuarter[0,0])/2.0
            yCOM = (fakeRealDataQuarter[3,1]+fakeRealDataQuarter[0,1])/2.0

        diFig = plt.figure(2,figsize=(10,9))
        main = diFig.add_subplot(111)
        #determine if to plot [mas] or ["]
        asConversion=1.0
        unitStr = '"'
        if abs(np.min([np.min(realData[:,1]),np.min(realData[:,3])]))<1.5:
            asConversion = 1000.0
            unitStr = 'mas'
        ## Draw orbit fit
        main.plot(fitDataDI[:,0]*asConversion,fitDataDI[:,1]*asConversion,linewidth=2.5,color='Blue') 
        ## Draw larger star for primary star's location
        starPolygon = star((asConversion/1000.0)*12.0*orbParams[10],0,0,color='yellow',N=6,thin=0.5)
        main.add_patch(starPolygon)
        ## Add DI data to plot
        (main,[xmin,xmax,ymin,ymax]) =  addDIdataToPlot(main,realData,asConversion)
        ## set limits and other basics of plot looks
        ##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ take max/min of fit data into account!!!!!!!!$$$$$$$$$
        xLim = (xmin-(xmax-xmin)*0.05,xmax+(xmax-xmin)*0.05)
        yLim = (ymin-(ymax-ymin)*0.05,ymax+(ymax-ymin)*0.05)
        ## FLIP X-AXIS to match backawards Right Ascension definition
        a = main.axis()
        main.axis([a[1],a[0],a[2],a[3]])
        main.axes.set_xlim((xLim[1],xLim[0]))
        main.axes.set_ylim(yLim)
        main.tick_params(axis='both',which='major',width=1,length=3,pad=10,direction='in',labelsize=15)
        main.spines['right'].set_linewidth(1.0)
        main.spines['bottom'].set_linewidth(1.0)
        main.spines['top'].set_linewidth(1.0)
        main.spines['left'].set_linewidth(1.0)
        main.set_position([0.19,0.15,0.79,0.83])
        main.set_xlabel('RA ['+unitStr+']', fontsize=30)
        main.set_ylabel('Dec ['+unitStr+']', fontsize=30)
        ## save fig to file
        plotFilename = plotFnameBase+'-DI.png'
        if plotFilename!='':
            plt.savefig(plotFilename, dpi=300, orientation='landscape')
            log.info("DI orbit plot saved to:\n"+plotFilename)
        plt.close()
        ## log params used in DI plot
        log.info('\n'+"*"*50+"\nOrbital Elements used in DI plot:\n"+repr(paramsDI)+'\n'+"*"*50+'\n')

    ################
    # Make RV plot #
    ################
    if settingsDict['dataMode'][0]!='DI':        
        realDataRV = copy.deepcopy(realData)
        ##Make model data for 100~1000 points for plotting fit
        nPts = 500
        fakeRealData = np.ones((nPts,8),dtype=np.dtype('d'),order='C')
        #set all RV offsets to zero
        fakeRealData[:,7]=0.0
        for i in range(0,nPts-1):
            fakeRealData[i,0] = orbParams[5]+(const.daysPerYear*orbParams[7]*(i/float(nPts)))
        fakeRealData[nPts-1,0]  = fakeRealData[0,0]+const.daysPerYear*orbParams[7]
        Orbit.loadRealData(fakeRealData)
        fitDataRV = np.ones((nPts,3),dtype=np.dtype('d'),order='C')
        orbParamsRV = copy.deepcopy(orbParams)
        orbParamsRV[9]+=settingsDict['omegaFrv'][0]
        ##Ensuring that params are in required format for SWIG
        paramsRV = []
        for par in orbParamsRV:
            paramsRV.append(par)
        paramsRV=np.array(paramsRV,dtype=np.dtype('d'),order='C')
        Orbit.calculate(fitDataRV,paramsRV)
        ## get the predicted values for the realData epochs
        Orbit.loadRealData(realDataRV)
        modelDataRV = np.ones((realDataRV.shape[0],3),dtype=np.dtype('d'),order='C')
        Orbit.calculate(modelDataRV,paramsRV)
        
        ##convert epochs to phases for plotting
        phasesReal = epochsToPhases(copy.deepcopy(realDataRV[:,0]),paramsRV[6],paramsRV[7], halfOrbit=False)
        phasesFit = epochsToPhases(copy.deepcopy(fitDataRV[:,0]),paramsRV[6],paramsRV[7], halfOrbit=False)
        
        ##Need to subtract RV offsets from the RVs 
        ##The fakeRealData had all offsets set to zero, so realDataRV needs to be "zeroed" to match
        numOffsets = int(len(paramsRV)-13)
        if numOffsets!=np.max(realDataRV[:7]):
            log.critical("Number of RV offsets in params does not match largest value in realData[:,7]")
        else:
            log.debug("There was a matching number of RV offsets in realData and params, = "+str(numOffsets))
            for i in range(0,realDataRV.shape[0]):
                rvBefore = realDataRV[i,5]
                realDataRV[i,5]-=paramsRV[13+int(realDataRV[i,7])]
                print str(rvBefore)+' - '+str(paramsRV[13+int(realDataRV[i,7])])+' = '+str(realDataRV[i,6])
        
            ## determine if to plot [KM/s] or [M/s]
            kmConversion = 1.0/1000.0
            unitStr = '[KM/s]'
            if np.max(realDataRV[:,5])<1500:
                kmConversion = 1.0
                unitStr = '[M/s]'
        
            #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            #######$$$$$$ CHECK ALL THE CODE BELOW THIS AND UPDATE TO 2.0 FORMAT!!!!!!
            #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            ## start making figure
            figRV = plt.figure(3,figsize=(10,5))
            residualsPlot = fig.add_subplot(212)
            residualsPlot.set_position([0.12,0.15,0.83,0.23])
            residualsPlot.axes.set_xlabel("Orbital Phase",fontsize=20)
            residualsPlot.axes.set_ylabel("Residual",fontsize=15)
            residualsPlot.tick_params(axis='y',which='major',width=3,length=5,pad=10,direction='in',labelsize=8)
            plt.locator_params(axis='y',nbins=5) #fix number of y-axis label points
            
            ## real-model=residual
            residualData = copy.deepcopy(realDataRV)
            residualData[:,5]-= modelDataRV[:,5]
            
            residualsPlot = addRVdataToPlot(residualsPlot,residualData,alf=1.0,color='blue',plotErrorBars=False)
            
            ## make plot of fit to data
            fitPlot = fig.add_subplot(211)
            fitPlot.set_position([0.12,0.38,0.83,0.55])
            fitPlot.xaxis.set_ticklabels([])#this is just a hack way of killing the tick labels
            fitPlot.plot(phasesFit,fitDataRV[:,5],c='r',linewidth=2.0,alpha=alf)
            fitPlot.axes.set_ylabel("RV "+unitStr,fontsize=20)
            fitPlot = addRVdataToPlot(fitPlot,realDataRV,alf=1.0,color='blue',plotErrorBars=False)
            
            
            ## set limits and other basics of plot looks
            ##$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ take max/min of fit data into account!!!!!!!!$$$$$$$$$
            residualsPlot.axes.set_xlim(??)
            residualsPlot.axes.set_ylim(??)
            #residualsPlot=fixPlotBordersAndLabels(residualsPlot)
            #plot zero vel line?
            residualsPlot.plot(??,[0,0],c='r',linewidth=2.0)
            fitPlot.axes.set_xlim(??)
            fitPlot.axes.set_ylim(??)
            #fitPlot=fixPlotBordersAndLabels(fitPlot)
            #yLabel = "RV [m/s]"
            #if KMperSec:
            #    yLabel = "RV [km/s]"
            #fitPlot.axes.set_ylabel(yLabel,fontsize=20)
            #residualsPlot=fixPlotBordersAndLabels(residualsPlot)
            
            
            
            ## save fig to file
            plotFilename = plotFnameBase+'-RV.png'
            if plotFilename!='':
                plt.savefig(plotFilename, dpi=300, orientation='landscape')
                log.info("RV orbit plot saved to:\n"+plotFilename)
            ## log params used in RV plot
            log.info('\n'+"*"*50+"\nOrbital Elements used in RV plot:\n"+repr(orbParamsRV)+'\n'+"*"*50+'\n')
    








