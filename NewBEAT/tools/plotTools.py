#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import numpy as np
import os
import pylab
import copy
import glob
import shutil
gridspec =  pylab.matplotlib.gridspec
plt = pylab.matplotlib.pyplot
patches = pylab.matplotlib.patches
MultLoc = pylab.matplotlib.ticker.MultipleLocator
import constants as const
import generalTools as genTools
import cppTools
import smodtLogger
import warnings
warnings.simplefilter("error")

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
        ## use the likelihoods as the weights
        theWeights = np.exp(-chiSquareds/2.0)
    else:
        theWeights = np.ones(len(data))      
    fig = plt.figure(1)
    subPlot = fig.add_subplot(111)
    ##numpy.histogram(a, bins=10, range=None, normed=False, weights=None, density=None) ##with density=True giving the probability distribution back
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


def histLoadAndPlot_StackedPosteriors(plot,outFilename='',xLabel='X',lineColor='k',xLims=False,latex=False,showYlabel=False,parInt=0):
    """
    Loads previously plotted histograms that were written to disk by histPlotAndDump, and plot them up 
    in a way that is ready for publication.  This version is to plot a posterior of the same parameter 
    for multiple simulation runs to see how they differ.
    
    It is foreseen that many versions of this function will exist for different specific publication ready plots.
    NOTE: this is extremely similar to histLoadAndPlot_ShadedPosteriors.  Re-factor to remove doubled code!!
    """
    if outFilename[-4]!='.dat':
        outFilename=outFilename+'.dat'
    histData = np.loadtxt(outFilename)  
    ys=[]
    xs=[]
    maxN = np.max(histData[:,1])
    minSub = 0
    valRange = np.max(histData[:,0])-np.min(histData[:,0])
    ## check if M2 and if it should be in jupiter masses
    if parInt==1:
        if np.max(histData[:,0])<0.02:
            histData[:,0]=histData[:,0]*(const.KGperMsun/const.KGperMjupiter)
            valRange = np.max(histData[:,0])-np.min(histData[:,0])
            xLabel='M2 [Mjupiter]'
            if latex:
                xLabel='$M2$ [$M_{jupiter}$]'
    if (np.max(histData[:,0])>100000) or (valRange<(np.min(histData[:,0])/100.0)):
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
    plot.plot(xs,ys,color=lineColor,linewidth=2)
    plot.axes.set_ylim([0.0,1.02])
    if xLims!=False:
        plot.axes.set_xlim(xLims)
    plot.locator_params(axis='x',nbins=4) # maximum number of x labels
    plot.locator_params(axis='y',nbins=5) # maximum number of y labels
    plot.tick_params(axis='x',which='major',width=0.5,length=3,pad=3,direction='in',labelsize=20)
    plot.tick_params(axis='y',which='major',width=0.5,length=3,pad=3,direction='in',labelsize=20)
    plot.spines['right'].set_linewidth(0.7)
    plot.spines['bottom'].set_linewidth(0.7)
    plot.spines['top'].set_linewidth(0.7)
    plot.spines['left'].set_linewidth(0.7)
    # add axes label
    if showYlabel:
        if latex:
            plot.axes.set_ylabel(r'$\frac{dp}{dx} \times constant $',fontsize=27)
        else:
            plot.axes.set_ylabel('dp/dx(*constant)',fontsize=25)
    else:
        plot.axes.set_yticklabels(['','',''])
    fsize=23
    if xLabel in ['e','$e$']:
        fsize=fsize+10
    if latex:
        plot.axes.set_xlabel(r''+xLabel,fontsize=fsize)
    else:
        plot.axes.set_xlabel(xLabel,fontsize=fsize)
    
    return plot

def histLoadAndPlot_ShadedPosteriors(plot,outFilename='',confLevels=False,xLabel='X',xLims=False,latex=False,showYlabel=False,parInt=0):
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
    valRange = np.max(histData[:,0])-np.min(histData[:,0])
    ## check if M2 and if it should be in jupiter masses
    if parInt==1:
        if np.max(histData[:,0])<0.02:
            histData[:,0]=histData[:,0]*(const.KGperMsun/const.KGperMjupiter)
            valRange = np.max(histData[:,0])-np.min(histData[:,0])
            xLabel='M2 [Mjupiter]'
            if latex:
                xLabel='$M2$ [$M_{jupiter}$]'
            confLevels=confLevels*(const.KGperMsun/const.KGperMjupiter)
    if (np.max(histData[:,0])>100000) or (valRange<(np.min(histData[:,0])/100.0)):
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
    plot.locator_params(axis='x',nbins=3) # maximum number of x labels
    plot.locator_params(axis='y',nbins=5) # maximum number of y labels
    plot.tick_params(axis='x',which='major',width=0.5,length=3,pad=3,direction='in',labelsize=20)
    plot.tick_params(axis='y',which='major',width=0.5,length=3,pad=3,direction='in',labelsize=20)
    plot.spines['right'].set_linewidth(0.7)
    plot.spines['bottom'].set_linewidth(0.7)
    plot.spines['top'].set_linewidth(0.7)
    plot.spines['left'].set_linewidth(0.7)
    # add axes label
    if showYlabel:
        if latex:
            plot.axes.set_ylabel(r'$\frac{dp}{dx} \times constant $',fontsize=27)
        else:
            plot.axes.set_ylabel('dp/dx(*constant)',fontsize=25)
    else:
        plot.axes.set_yticklabels(['','',''])
    fsize=23
    if xLabel in ['e','$e$']:
        fsize=fsize+10
    if latex:
        plot.axes.set_xlabel(r''+xLabel,fontsize=fsize)
    else:
        plot.axes.set_xlabel(xLabel,fontsize=fsize)
        
    return plot

def addRVdataToPlot(subPlot,epochsORphases,RVs,RVerrs,alf=1.0,color='blue',plotErrorBars=False):
    """
    Add '+' markers for the data locations with respective y axis errors 
    shown as the height of the markers. 
    """
    for i in range(0,RVs.shape[0]):
        if RVerrs[i]<1e3:
            xs = [epochsORphases[i],epochsORphases[i]]
            ys = [RVs[i]-RVerrs[i],RVs[i]+RVerrs[i]]
            #print str(RVerrs[i])+", -> ["+str(epochsORphases[i])+", "+str(RVs[i])+']'
            if plotErrorBars:
                subPlot.plot(xs,ys,c=color,linewidth=1,alpha=alf)
            subPlot.plot(epochsORphases[i],RVs[i],c='k',marker='.',markersize=6)
    return subPlot

def addDIdataToPlot(subPlot,realData,asConversion):
    """
    To plot a '+' for each data point with width and height matching the errors converted 
    to x,y coords.
    """
    ## copy realData and kill off parts where DI errors are 1e6
    diData = copy.deepcopy(realData)
    diData = diData[np.where(diData[:,2]<1e6)[0],:]
    xmin = np.min(diData[:,1]-diData[:,2])*asConversion
    xmax = np.max(diData[:,1]+diData[:,2])*asConversion
    ymin = np.min(diData[:,3]-diData[:,4])*asConversion
    ymax = np.max(diData[:,3]+diData[:,4])*asConversion
    for i in range(0,diData.shape[0]):
        xCent = diData[i,1]*asConversion
        yCent = diData[i,3]*asConversion
        #print 'data [x,y] = ['+str(xCent/asConversion)+', '+str(yCent/asConversion)+']'#$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        left = xCent-diData[i,2]*asConversion
        right = xCent+diData[i,2]*asConversion
        top = yCent+diData[i,4]*asConversion
        btm = yCent-diData[i,4]*asConversion
        subPlot.plot([left,right],[yCent,yCent],linewidth=1.5,color='k',alpha=1.0)
        subPlot.plot([xCent,xCent],[btm,top],linewidth=1.5,color='k',alpha=1.0)
    return (subPlot,[xmin,xmax,ymin,ymax])

def stackedPosteriorsPlotter(outputDataFilenames, plotFilename,paramsToPlot=[],xLims=[],stage='MCMC'):
    """
    This will plot a simple posterior distribution for each parameter in the data files
    stacked ontop of each other for comparison between different runs.
    It can only be ran on folders which have their histograms already calculated and saved in the 
    plotData subfolder!
    
    It is called by the custom function stackedPosteriorsPlotterHackStarter in rePlot.py.
    
    Very similar to summaryPlotter, but the loop is first in parameter int, then
    file to ensure all are stacked on same subplot properly.  
    NOTE: might be able to extract doubled code to clean things up...
    """
    log.setStreamLevel(lvl=10)
    latex=True
    plotFormat = 'eps'   
    plt.rcParams['ps.useafm']= True
    plt.rcParams['pdf.use14corefonts'] = True
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    if latex:
        plt.rc('text', usetex=True)
    
    if type(outputDataFilenames)!=list:
        outputDataFilenames = [outputDataFilenames]
    
    colorsList =['Red','DarkGreen','Blue','Black','Purple','Fuchsia','Crimson','Aqua','Gold','OrangeRed','Plum','Chartreuse','Chocolate','SteelBlue ','Teal','Salmon','Brown']
    
    if os.path.exists(outputDataFilenames[0]):  
        log.debug('\nCreating a simple plot of some key posteriors for files:\n'+repr(outputDataFilenames))
        log.debug("writing resulting figure to:\n"+plotFilename)
        
        ## load first data file to get param lists 
        (head,data) = genTools.loadFits(outputDataFilenames[0])
        (paramList,paramStrs,paramFileStrs) = genTools.getParStrs(head,latex=latex)
        (paramList2,paramStrs2,paramFileStrs2) = genTools.getParStrs(head,latex=False)
        ## modify x labels to account for DI only situations where M1=Mtotal
        if np.var(data[:,1])==0:
            paramStrs2[0] = 'M total [Msun]'
            paramStrs[0] = '$M_{total}$ [$M_{\odot}$]'
            paramFileStrs[0] = 'Mtotal'
        ## check if a subset is to be plotted or the whole set
        ## remake lists of params to match subset.
        if len(paramsToPlot)!=0:
            paramStrs2Use = []
            paramStrsUse = []
            paramFileStrsUse = []
            paramListUse = []
            for par in paramsToPlot:
                paramStrs2Use.append(paramStrs2[par])
                paramStrsUse.append(paramStrs[par])
                paramFileStrsUse.append(paramFileStrs[par])
                paramListUse.append(par)
            paramStrs2 = paramStrs2Use
            paramStrs = paramStrsUse
            paramFileStrs = paramFileStrsUse 
            paramList = paramListUse
        ## determine appropriate figure size for number of params to plot
        figSizes =  [(7,4),(8,4),(9,4),(10,4), (12,8),(12,11),(12,13),(12,15)]
        gridSizes = [(1,1),(1,2),(1,3),(1,4),  (2,4),(3,4),(4,4),(5,4)]
        sz = 0
        if len(paramStrs2)>20:
            log.critical("summaryPlotter is only capable of handling up to 20 parameters, "+str(len(paramStrs2))+" were passed in!")
        elif len(paramStrs2)>16:
            sz = 7
        elif len(paramStrs2)>12:
            sz = 6
        elif len(paramStrs2)>8:
            sz = 5
        elif len(paramStrs2)>4:
            sz = 4
        elif len(paramStrs2)>3:
            sz = 3
        elif len(paramStrs2)>2:
            sz = 2
        elif len(paramStrs2)>1:
            sz = 1
        ## Create empty figure to be filled up with plots
        stackedFig = plt.figure(figsize=figSizes[sz])     
        
        ## Go through params and re-load the hist for each files and plot them 
        for i in range(0,len(paramStrs2)):
            colorInt = 0
            log.debug('Starting to plot stacked hist for '+paramStrs2[i])
            subPlot = plt.subplot(gridSizes[sz][0],gridSizes[sz][1],i+1)
            xLim=False
            if len(paramsToPlot)!=0:
                xLim=xLims[i]
            showYlabel=False
            if i in [0,4,8,12]:
                showYlabel = True
            par=0
            try:
                par = paramList[i]
            except:
                log.warning("Parameter "+str(i)+" not in paramList: \n"+repr(paramList))
            ## go through each file and plot this param's hist on same plot
            for outputDataFilename in outputDataFilenames:
                log.debug('Loading and re-plotting parameter '+str(i+1)+"/"+str(len(paramStrs2))+": "+paramStrs2[i]+" for file:\n"+outputDataFilename)
                ## check if plot data dir exists
                plotDataDir = os.path.join(os.path.dirname(outputDataFilename),"plotData")
                if os.path.exists(plotDataDir)==False:      
                    log.critical("PlotDataDir doesn't exist!! at:\n"+plotDataDir)
                else:
                    plotDataDir+='/'
                    histDataBaseName = os.path.join(os.path.dirname(outputDataFilename),'hist-'+stage+"-"+paramFileStrs[i])
                    if os.path.exists(os.path.join(os.path.dirname(plotDataDir),'hist-'+stage+"-"+paramFileStrs[i]+'.dat')):
                        histDataBaseName = os.path.join(os.path.dirname(plotDataDir),'hist-'+stage+"-"+paramFileStrs[i])
                    if os.path.exists(histDataBaseName+'.dat'):
                        log.debug("plotting file:\n"+histDataBaseName)
                        subPlot = histLoadAndPlot_StackedPosteriors(subPlot,outFilename=histDataBaseName,xLabel=paramStrs[i],lineColor=colorsList[colorInt],xLims=xLim,latex=latex,showYlabel=showYlabel,parInt=par)
                    else:
                        log.debug("Not plotting hist for "+paramStrs2[i]+" as its hist file doesn't exist:\n"+histDataBaseName)
                    colorInt+=1
        plt.tight_layout()        
        ## Save file if requested.
        log.debug('\nStarting to save stacked param hist figure:')
        if plotFilename!='':
            plt.savefig(plotFilename+'.'+plotFormat,format=plotFormat)
            s= 'stacked hist plot saved to: '+plotFilename+'.'+plotFormat
            log.info(s)
        plt.close()
        if True:
            log.debug('converting to PDF as well')
            try:
                os.system("epstopdf "+plotFilename+'.'+plotFormat)
            except:
                log.warning("Seems epstopdf failed.  Check if it is installed properly.")
        
def summaryPlotter(outputDataFilename,plotFilename,paramsToPlot=[],xLims=[],stage='MCMC',shadeConfLevels=True,forceRecalc=True):
    """
    This advanced plotting function will plot all the data in a grid on a single figure.  The data will be plotted
    in histograms that will be normalized to a max of 1.0.  The 
    confidence levels of the data can be calculated and the resulting histograms will have the bars inside the
    68% bin shaded as dark grey and 95% as light grey and the rest will be white.
    
    NOTE: currently a maximum of 20 parameters can be plotted, ie. for a system with up to 7 RV data sets.
    """
    latex=True
    plotFormat = 'eps'   
    plt.rcParams['ps.useafm']= True
    plt.rcParams['pdf.use14corefonts'] = True
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    if latex:
        plt.rc('text', usetex=True)
        #plt.rcParams['text.latex.unicode']=True 
        #plt.rcParams['text.latex.preamble'] = '\usepackage{amssymb}' 
        #plt.rcParams['text.latex.preamble'] = '\usepackage{sfmath}' 
    
    ## check if plot data dir exists, else make it
    plotDataDir = os.path.join(os.path.dirname(outputDataFilename),"plotData")
    #print 'plotDataDir = '+plotDataDir
    if os.path.exists(plotDataDir)==False:      
        os.mkdir(plotDataDir)
    plotDataDir+='/'
    
    (head,data) = genTools.loadFits(outputDataFilename)
    if head!=False:  
        log.debug(' Inside summaryPlotter')
         
        s= '\nCreating summary plot for file:\n'+outputDataFilename
        s=s+ '\nInput plotfilename:\n'+plotFilename
        log.info(s)
        
        ## check if the passed in value for plotFilename includes format extension
        if '.'+plotFormat not in plotFilename:
            plotFilename = plotFilename+"."+plotFormat
            log.debug('updating plotFilename to:\n'+plotFilename)
        else:
            plotFilename = plotFilename
                
        (paramList,paramStrs,paramFileStrs) = genTools.getParStrs(head,latex=latex)
        (paramList2,paramStrs2,paramFileStrs2) = genTools.getParStrs(head,latex=False)
        
        ## modify x labels to account for DI only situations where M1=Mtotal
        if np.var(data[:,1])==0:
            paramStrs2[0] = 'M total [Msun]'
            paramStrs[0] = '$M_{total}$ [$M_{\odot}$]'
            paramFileStrs[0] = 'Mtotal'
        
        ## check if a subset is to be plotted or the whole set
        ## remake lists of params to match subset.
        if len(paramsToPlot)!=0:
            paramStrs2Use = []
            paramStrsUse = []
            paramFileStrsUse = []
            paramListUse = []
            for par in paramsToPlot:
                paramStrs2Use.append(paramStrs2[par])
                paramStrsUse.append(paramStrs[par])
                paramFileStrsUse.append(paramFileStrs[par])
                paramListUse.append(par)
            paramStrs2 = paramStrs2Use
            paramStrs = paramStrsUse
            paramFileStrs = paramFileStrsUse 
            paramList = paramListUse
        ## determine appropriate figure size for number of params to plot
        figSizes =  [(5.5,7),(8,7),(9,7),(10,3.5),(11,8),(11,11),(11,14),(11,16)]
        gridSizes = [(1,1),(1,2),(1,3),(1,4),(2,4),(3,4),(4,4),(5,4)]
        sz = 0
        if len(paramStrs2)>20:
            log.critical("summaryPlotter is only capable of handling up to 20 parameters, "+str(len(paramStrs2))+" were passed in!")
        elif len(paramStrs2)>16:
            sz = 7
        elif len(paramStrs2)>12:
            sz = 6
        elif len(paramStrs2)>8:
            sz = 5
        elif len(paramStrs2)>4:
            sz = 4
        elif len(paramStrs2)>3:
            sz = 3
        elif len(paramStrs2)>2:
            sz = 2
        elif len(paramStrs2)>1:
            sz = 1
#         wRatios = []
#         hRatios = []
#         for i in range(0,len(paramStrs2)):
#             wRatios.append(1)
#             hRatios.append(3)
#         print 'gridsize = '+repr(gridSizes[sz])
#         print 'figSize = '+repr(figSizes[sz])
#         gs = gridspec.GridSpec(gridSizes[sz][0],gridSizes[sz][1])#,width_ratios=wRatios,height_ratios=hRatios)    
#         gs.update(wspace=0.1,hspace=0)
        
        ## run through all the data files and parameters requested and make histogram files
        completeCLstr = '-'*22+'\nConfidence Levels are:\n'+'-'*75+'\n'
        for i in range(0,len(paramStrs2)):
            if (os.path.exists(os.path.join(os.path.dirname(plotDataDir),'hist-'+stage+"-"+paramFileStrs[i]+'.dat'))==False)or forceRecalc:
                log.debug('Checking parameter has useful data '+str(i+1)+"/"+str(len(paramStrs2))+": "+paramStrs2[i]+", for file:\n"+outputDataFilename)
                (CLevels,data,bestDataVal,clStr) = genTools.confLevelFinder(outputDataFilename,i, returnData=True, returnChiSquareds=False, returnBestDataVal=True)
                if bestDataVal!=0:
                    completeCLstr+=paramStrs2[i]+clStr+'\n'+'-'*75+'\n'
                    log.debug('Making hist file for parameter '+str(i+1)+"/"+str(len(paramStrs2))+": "+paramStrs2[i]+", for file:\n"+outputDataFilename)
                    histDataBaseName = os.path.join(os.path.dirname(plotDataDir),'hist-'+stage+"-"+paramFileStrs[i])
                    #print 'histDataBaseName = '+histDataBaseName
                    histMakeAndDump([],data,outFilename=histDataBaseName,weight=False, normed=False, nu=1,logY=False,histType='step')
                    if (os.path.exists(os.path.join(os.path.dirname(plotDataDir),'confLevels-'+stage+"-"+paramFileStrs[i]+'.dat'))==False)or forceRecalc:
                        np.savetxt(os.path.join(os.path.dirname(plotDataDir),'confLevels-'+stage+"-"+paramFileStrs[i]+'.dat'),CLevels)
                        log.debug('confidence levels data stored to:\n'+os.path.join(os.path.dirname(plotDataDir),'confLevels-'+stage+"-"+str(i)+'.dat'))
                else:
                    log.debug("Nope! no useful data for "+str(i+1)+"/"+str(len(paramStrs2))+": "+paramStrs2[i]+", in file:\n"+outputDataFilename)
       
        ## Create empty figure to be filled up with plots
        sumFig = plt.figure(figsize=figSizes[sz])       
        ## make shaded posterior for each param
        for i in range(0,len(paramStrs2)):
            histDataBaseName = os.path.join(os.path.dirname(outputDataFilename),'hist-'+stage+"-"+paramFileStrs[i])
            if os.path.exists(os.path.join(os.path.dirname(plotDataDir),'hist-'+stage+"-"+paramFileStrs[i]+'.dat')):
                histDataBaseName = os.path.join(os.path.dirname(plotDataDir),'hist-'+stage+"-"+paramFileStrs[i])
            if os.path.exists(histDataBaseName+'.dat'):
                log.debug('Starting to plot shaded hist for '+paramStrs2[i])
                #print 'gs[i] = '+repr(gs[0,i])+'\n\n'
                subPlot = plt.subplot(gridSizes[sz][0],gridSizes[sz][1],i+1)#gs[i])
                log.debug('Loading and re-plotting parameter '+str(i+1)+"/"+str(len(paramStrs2))+": "+paramStrs2[i])#+" for file:\n"+outputDataFilename
                CLevels=False
                xLim=False
                if len(paramsToPlot)!=0:
                    xLim=xLims[i]
                if shadeConfLevels:
                    clFile = os.path.join(os.path.dirname(outputDataFilename),'confLevels-'+stage+"-"+paramFileStrs[i]+'.dat')
                    if os.path.exists(os.path.join(os.path.dirname(plotDataDir),'confLevels-'+stage+"-"+paramFileStrs[i]+'.dat')):
                        clFile = os.path.join(os.path.dirname(plotDataDir),'confLevels-'+stage+"-"+paramFileStrs[i]+'.dat')
                    CLevels=np.loadtxt(clFile)
                showYlabel=False
                if i in [0,4,8,12]:
                    showYlabel = True
                par=0
                try:
                    par = paramList[i]
                except:
                    log.debug("Parameter "+str(i)+" not in paramList: \n"+repr(paramList))
                subPlot = histLoadAndPlot_ShadedPosteriors(subPlot,outFilename=histDataBaseName,confLevels=CLevels,xLabel=paramStrs[i],xLims=xLim,latex=latex,showYlabel=showYlabel,parInt=par)         
            else:
                log.debug("Not plotting shaded hist for "+paramStrs2[i]+" as its hist file doesn't exist:\n"+histDataBaseName)
        plt.tight_layout()
        ## Save file if requested.
        log.debug('\nStarting to save param hist figure:')
        if plotFilename!='':
            plt.savefig(plotFilename,format=plotFormat)
            s= 'Summary plot saved to: '+plotFilename
            log.info(s)
        plt.close()
        if True:
            log.debug('converting to PDF as well')
            try:
                os.system("epstopdf "+plotFilename)
            except:
                log.warning("Seems epstopdf failed.  Check if it is installed properly.")
        
        return completeCLstr
            
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
    (ratio of how far from Tc it is), shifted to lie inside [0,1].
    if 'halfOrbit'=True, the vals will lie inside [-0.5,0.5].
    """    
    verbose=False         
    phases = []
    P_days = P_yrs*const.daysPerYear
    for epoch in epochs:
        phaseTimeDiff = epoch - int((epoch-Tc)/P_days)*P_days-Tc #phase shifted into [Tc,Tc+P]
        if verbose:
            print str(epoch)+" - "+str(int((epoch-Tc)/P_days)*P_days)+" - "+str(Tc)+" = "+str(phaseTimeDiff)
        phase = phaseTimeDiff/P_days#phase shifted into [0,1]
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

def orbitPlotter(orbParams,settingsDict,plotFnameBase="",format='png',DIlims=[],RVlims=[]):
    """
    Make both the DI and RV plots.
    '-DI.png' and/or '-RV.png' will be added to end of plotFnameBase 
    to make the filenames for each type of plot.
    
    Optional tweaks:
    DIlims=[[[xMin,xMax],[yMin,yMax]],[[xCropMin,xCropMax],[yCropMin,yCropMax]]]
    RVlims=[[yMin,yMax],[yResidMin,yResidMax],[xMin,xMax]]
    """
    latex=True
    plt.rcParams['ps.useafm']= True
    plt.rcParams['pdf.use14corefonts'] = True
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    if latex:
        plt.rc('text', usetex=True)
        plt.rcParams['text.latex.unicode']=True 
        plt.rcParams['text.latex.preamble'] = '\usepackage{amssymb}' 
        #plt.rcParams['text.latex.preamble'] = '\usepackage{sfmath}' 
    else:
        plt.rc('font',family='serif')
        plt.rc('text', usetex=False)
    log.debug("Starting to make orbit plots")
    colorsList =['Blue','BlueViolet','Chartreuse','Fuchsia','Crimson','Aqua','Gold','DarkCyan','OrangeRed','Plum','DarkGreen','Chocolate','SteelBlue ','Teal','Salmon','Brown']

    ## check if plot data dir exists, else make it
    plotDataDir = os.path.join(os.path.dirname(plotFnameBase),"plotData")
    if os.path.exists(plotDataDir)==False:      
        os.mkdir(plotDataDir)
    plotDataDir+='/'
    ##get the real data
    realData = genTools.loadRealData(os.path.join(settingsDict['settingsDir'],settingsDict['prepend']),dataMode=settingsDict['dataMode'][0])
    ## Make Orbit cpp obj
    Orbit = cppTools.Orbit()
    Orbit.loadStaticVars(settingsDict['omegaFdi'][0],settingsDict['omegaFrv'][0],settingsDict['lowEcc'][0])
    Orbit.loadConstants(const.Grav,const.pi,const.KGperMsun, const.daysPerYear,const.secPerYear,const.MperAU)
    ## ensure orbParams are in required format for Orbit
    params = []
    for par in orbParams:
        params.append(par)
    params=np.array(params,dtype=np.dtype('d'),order='C')
    
    ################
    # Make DI plot #
    ################
    if settingsDict['dataMode'][0]!='RV':
        paramsDI = copy.deepcopy(params)
        realDataDI = copy.deepcopy(realData)
        realDataDI = realDataDI[np.where(realDataDI[:,2]<1e6)[0],:]
        ##Make model data for 100~1000 points for plotting fit
        nPts = 500
        fakeRealData = np.zeros((nPts,8),dtype=np.dtype('d'),order='C')
        fakeRealData[:,1:5]=1.0
        fakeRealData[:,6]=1e6
        for i in range(0,nPts-1):
            fakeRealData[i,0] = paramsDI[5]+(const.daysPerYear*paramsDI[7]*(i/float(nPts)))
        fakeRealData[nPts-1,0]  = fakeRealData[0,0]+const.daysPerYear*paramsDI[7]
        #print 'fakeRealData= '+repr(fakeRealData)
        Orbit.loadRealData(fakeRealData)
        fitDataDI = np.ones((nPts,3),dtype=np.dtype('d'),order='C')
        Orbit.calculate(fitDataDI,paramsDI)
        ## calculate the fit locations for the DI epochs to calculate 0-C
        Orbit.loadRealData(realDataDI)
        predictedDataDI = np.ones((realDataDI.shape[0],3),dtype=np.dtype('d'),order='C')
        Orbit.calculate(predictedDataDI,paramsDI)
        ## Get locations of start/end for semi-major axis or COM, and AN/DN for line-of-nodes
        ## Get 1/4 locations (useful for drawing semi-major axis, and finding loc of COM)
        fakeRealDataQuarter = np.ones((4,8),dtype=np.dtype('d'),order='C')
        for i in range(0,4):
            fakeRealDataQuarter[i,0] = paramsDI[5]+(const.daysPerYear*paramsDI[7]*(i/4.0))
        Orbit.loadRealData(fakeRealDataQuarter)
        fitDataQuarter = np.zeros((4,3),dtype=np.dtype('d'),order='C')
        Orbit.calculate(fitDataQuarter,paramsDI)
        ## make semi-major locs
        semiMajorLocs = np.array([[fitDataQuarter[0,0],fitDataQuarter[0,1]],[fitDataQuarter[2,0],fitDataQuarter[2,1]]])
        ## find loc of COM for possible use
        xCOM = (fakeRealDataQuarter[3,0]+fakeRealDataQuarter[0,0])/2.0
        yCOM = (fakeRealDataQuarter[3,1]+fakeRealDataQuarter[0,1])/2.0
        ## Find Ascending and Descending Node locations
        nodeEpochs = nodeEpochsCalc(paramsDI,settingsDict["omegaFdi"][0]) 
        print 'period/2 = '+repr(const.daysPerYear*paramsDI[7]*(1.0/2.0))
        print 'nodeEpochs = '+repr(nodeEpochs)
        lonData = np.ones((2,8),dtype=np.dtype('d'),order='C')
        lonData[:,0]=nodeEpochs
        Orbit.loadRealData(lonData)
        tmpData = np.ones((2,3),dtype=np.dtype('d'),order='C')
        Orbit.calculate(tmpData,paramsDI)
        lonXYs = tmpData[:,:2]#[[tmpData[0,0],tmpData[0,1]]]
        
        ##load resulting data to file for re-plotting by others
        #real [x,xerr,y,yerr]
        outDIdataReal = realDataDI[:,1:5]
        #fit [x,y]
        outDIdataFit = []
        for i in range(0,len(fitDataDI[:,0])):
            outDIdataFit.append([fitDataDI[i,0],fitDataDI[i,1]])
        outPredictedDIdata = []
        for i in range(0,len(predictedDataDI[:,0])):
            outPredictedDIdata.append([predictedDataDI[i,0],predictedDataDI[i,1]])
        fnameBase = os.path.join(os.path.dirname(plotDataDir),'DIplotData')
        np.savetxt(fnameBase+'-real.dat',outDIdataReal)
        np.savetxt(fnameBase+'-predicted.dat',outPredictedDIdata)
        np.savetxt(fnameBase+'-fit.dat',outDIdataFit)

        diFig = plt.figure(2,figsize=(10,9))
        main = diFig.add_subplot(111)
        #determine if to plot [mas] or ["]
        asConversion=1.0
        unitStr = '"'
        if abs(np.min([np.min(realDataDI[:,1]),np.min(realDataDI[:,3])]))<1.5:
            asConversion = 1000.0
            unitStr = '[mas]'
        ## Draw orbit fit
        main.plot(fitDataDI[:,0]*asConversion,fitDataDI[:,1]*asConversion,linewidth=1,color='Blue') 
        ## Draw line-of-nodes
        main.plot(lonXYs[:,0]*asConversion,lonXYs[:,1]*asConversion,'-.',linewidth=1.0,color='Green')
        #print 'lonXYs*asConversion = '+repr(lonXYs*asConversion)
        ## Draw Semi-Major axis
        main.plot(semiMajorLocs[:,0]*asConversion,semiMajorLocs[:,1]*asConversion,'-',linewidth=1.0,color='Green')
        ## Draw larger star for primary star's location
        starPolygon = star(2*paramsDI[10],0,0,color='yellow',N=6,thin=0.5)
        main.add_patch(starPolygon)
        ## Add DI data to plot
        (main,[xmin,xmax,ymin,ymax]) =  addDIdataToPlot(main,realDataDI,asConversion)
        ## set limits and other basics of plot looks
        xLims = (np.min([xmin,np.min(fitDataDI[:,0]*asConversion)]),np.max([xmax,np.max(fitDataDI[:,0]*asConversion)]))
        yLims = (np.min([ymin,np.min(fitDataDI[:,1]*asConversion)]),np.max([ymax,np.max(fitDataDI[:,1]*asConversion)]))
        xLimsFull = (xLims[0]-(xLims[1]-xLims[0])*0.05,xLims[1]+(xLims[1]-xLims[0])*0.05)
        yLimsFull = (yLims[0]-(yLims[1]-yLims[0])*0.05,yLims[1]+(yLims[1]-yLims[0])*0.05)
        xLims = (xmin,xmax)
        yLims = (ymin,ymax)
        xLimsCrop = (xLims[0]-(xLims[1]-xLims[0])*0.05,xLims[1]+(xLims[1]-xLims[0])*0.05)
        yLimsCrop = (yLims[0]-(yLims[1]-yLims[0])*0.05,yLims[1]+(yLims[1]-yLims[0])*0.05)
        ## FLIP X-AXIS to match backawards Right Ascension definition
        a = main.axis()
        main.axis([a[1],a[0],a[2],a[3]])
        ## update lims with custom values if provided
        if len(DIlims)>0:
            #DIlims=[[[xMin,xMax],[yMin,yMax]],[[xCropMin,xCropMax],[yCropMin,yCropMax]]]
            xLimsFull=(DIlims[0][0][0],DIlims[0][0][1])
            yLimsFull=(DIlims[0][1][0],DIlims[0][1][1])
            xLimsCrop=(DIlims[1][0][0],DIlims[1][0][1])
            yLimsCrop=(DIlims[1][1][0],DIlims[1][1][1])
        main.axes.set_xlim((xLimsFull[1],xLimsFull[0]))
        main.axes.set_ylim(yLimsFull)
        plt.minorticks_on()
        main.tick_params(axis='both',which='major',width=1,length=5,pad=10,direction='in',labelsize=25)
        main.tick_params(axis='both',which='minor',width=1,length=2,pad=10,direction='in')
        main.spines['right'].set_linewidth(1.0)
        main.spines['bottom'].set_linewidth(1.0)
        main.spines['top'].set_linewidth(1.0)
        main.spines['left'].set_linewidth(1.0)
        main.set_position([0.15,0.115,0.80,0.84])
        xLabel = 'Relative RA '+unitStr
        yLabel = 'Relative Dec '+unitStr
        if latex:
            xLabel = '$Relative$ $RA$ '+unitStr#'$\alpha$ '+unitStr
            yLabel = '$Relative$ $Dec$ '+unitStr#'$\delta$ '+unitStr
        main.set_xlabel(xLabel, fontsize=25)
        main.set_ylabel(yLabel, fontsize=25)
        #plt.tight_layout()
        ## save fig to file and maybe convert to pdf if format=='eps'
        orientStr = 'landscape'
        if format=='eps':
            orientStr = 'portrait'
        plotFilenameFull = plotFnameBase+'-DI.'+format
        if plotFilenameFull!='':
            plt.savefig(plotFilenameFull, dpi=300, orientation=orientStr)
            log.info("DI orbit plot (Full) saved to:\n"+plotFilenameFull)
        ##crop to limits of data and save
        main.axes.set_xlim((xLimsCrop[1],xLimsCrop[0]))
        main.axes.set_ylim(yLimsCrop)
        ## plot predicted locations for the data points
        if True:
            for i in range(0,len(predictedDataDI[:,0])):
                main.plot(predictedDataDI[i,0]*asConversion,predictedDataDI[i,1]*asConversion,c='red',marker='.',markersize=5)
                #print 'plotted point ['+str(predictedDataDI[i,0]*asConversion)+', '+str(predictedDataDI[i,1]*asConversion)+']'
        plotFilenameCrop = plotFnameBase+'-DI-cropped.'+format
        if plotFilenameCrop!='':
            plt.savefig(plotFilenameCrop, dpi=300, orientation=orientStr)
            log.info("DI orbit plot (cropped) saved to:\n"+plotFilenameCrop)
        plt.close()
        if (format=='eps')and True:
            log.debug('converting to PDF as well')
            try:
                os.system("epstopdf "+plotFilenameFull)
                os.system("epstopdf "+plotFilenameCrop)
            except:
                log.warning("Seems epstopdf failed.  Check if it is installed properly.")
        ## log params used in DI plot
        log.info('\n'+"*"*50+"\nOrbital Elements used in DI plot:\n"+repr(paramsDI))
        log.info("\n with an omega value = "+str(paramsDI[9]+settingsDict["omegaFdi"][0])+'\n'+"*"*50+'\n')

    ################
    # Make RV plot #
    ################
    if settingsDict['dataMode'][0]!='DI':        
        realDataRV = copy.deepcopy(realData)
        realDataRV = realDataRV[np.where(realDataRV[:,6]<1e6)[0],:]
        ##Ensuring that params are in required format for SWIG
        paramsRV = copy.deepcopy(params)
        ##Make model data for 100~1000 points for plotting fit
        nPts = 500
        fakeRealData = np.zeros((nPts-1,8),dtype=np.dtype('d'),order='C')
        fakeRealData[:,5] = 1.0
        #set all RV offsets to zero
        fakeRealData[:,7] = 0.0
        fakeOrbParams = copy.deepcopy(paramsRV)
        fakeOrbParams[13:]=0.0
        #print 'fakeOrbParams = '+repr(fakeOrbParams)
        fakeRealData[:,0] = paramsRV[6]-(const.daysPerYear*paramsRV[7]/2.0)
        for i in range(0,nPts-1):
            fakeRealData[i,0] += const.daysPerYear*paramsRV[7]*((i+1)/float(nPts))
        Orbit.loadRealData(fakeRealData)
        fitDataRV = np.ones((nPts-1,3),dtype=np.dtype('d'),order='C')
        Orbit.calculate(fitDataRV,fakeOrbParams)
        ## get the predicted values for the realData epochsparamsRV
        Orbit.loadRealData(realDataRV)
        modelDataRV = np.ones((realDataRV.shape[0],3),dtype=np.dtype('d'),order='C')
        #print 'paramsRV = '+repr(paramsRV)
        Orbit.calculate(modelDataRV,paramsRV)
        
        ##Need to subtract RV offsets from the RVs 
        ##The fakeRealData had all offsets set to zero, so realDataRV needs to be "zeroed" to match
        numOffsets = int(len(paramsRV)-13)
        if (numOffsets-1)!=np.max(realDataRV[:,7]):
            log.critical("Number of RV offsets in params does not match largest value in realData[:,7]")
            log.critical("# of offsets in params = "+str(numOffsets)+" != # max in realData = "+str(np.max(realDataRV[:,7])+1))
        else:
            log.debug("There was a matching number of RV offsets in realData and params, = "+str(numOffsets))
            #for i in range(0,realDataRV.shape[0]):
            #    if realDataRV[i,6]<1000000.0:
            #        print 'before offset subtract realData = '+str(realDataRV[i,0])+', '+str(realDataRV[i,5])+", "+str(realDataRV[i,6])+", "+str(realDataRV[i,7])
            zeroedRealDataRV = copy.deepcopy(realDataRV)
            for i in range(0,zeroedRealDataRV.shape[0]):
                rvBefore = zeroedRealDataRV[i,5]
                zeroedRealDataRV[i,5]-=paramsRV[13+int(zeroedRealDataRV[i,7])]
                #print str(rvBefore)+' - '+str(paramsRV[13+int(zeroedRealDataRV[i,7])])+" = "+str(zeroedRealDataRV[i,5])
            
            ##convert epochs to phases for plotting
            phasesReal = epochsToPhases(copy.deepcopy(realDataRV[:,0]),paramsRV[6],paramsRV[7], halfOrbit=True)
            phasesFit = epochsToPhases(copy.deepcopy(fakeRealData[:,0]),paramsRV[6],paramsRV[7], halfOrbit=True)            
            
            ## determine if to plot [km/s] or [m/s]
            kmConversion = 1.0/1000.0
            unitStr = '[km/s]'
            if np.max(np.sqrt(zeroedRealDataRV[:,5]**2.0))<1500:
                kmConversion = 1.0
                unitStr = '[m/s]'
            ## start making figure for residual and fit plots
            figRV = plt.figure(3,figsize=(10,5))
            residualsPlot = figRV.add_subplot(212)
            residualsPlot.set_position([0.13,0.15,0.84,0.23])
            fitPlot = figRV.add_subplot(211)
            fitPlot.set_position([0.13,0.38,0.84,0.57])
            xLabel = "Orbital Phase"
            fitYlabel = 'vr '
            residYlabel = 'O-C '
            if latex:
                residYlabel = '$O-C$ '
                fitYlabel = '$vr$ '
                xLabel = "$Orbital$ $Phase$"
            residualsPlot.axes.set_xlabel(xLabel,fontsize=20)
            residualsPlot.axes.set_ylabel(residYlabel,fontsize=20)
            fitPlot.xaxis.set_ticklabels([])#this is just a hack way of killing the tick labels
            fitPlot.axes.set_ylabel(fitYlabel+unitStr,fontsize=20)
            
            ## real-model=residual, then plot it
            residualData = copy.deepcopy(realDataRV)
            residualData[:,5]-= modelDataRV[:,2]
            
            #for i in range(0,residualData.shape[0]):
            #    print 'residual Data = '+str(residualData[i,0])+', '+str(residualData[i,5])+", "+str(residualData[i,6])+", "+str(residualData[i,7])
            #for i in range(0,modelDataRV.shape[0]):
            #    print 'model Data = '+str(modelDataRV[i,2])
            
            ## add real data to plots
            residualsPlot = addRVdataToPlot(residualsPlot,phasesReal,residualData[:,5]*kmConversion,residualData[:,6]*kmConversion,alf=0.1,color='k',plotErrorBars=True)
            fitPlot = addRVdataToPlot(fitPlot,phasesReal,zeroedRealDataRV[:,5]*kmConversion,zeroedRealDataRV[:,6]*kmConversion,alf=0.2,color='k',plotErrorBars=True)
            ##plot fit epochsORphases,RVs,RVerrs
            fitPlot.plot(phasesFit,fitDataRV[:,2]*kmConversion,c='Blue',linewidth=1.0,alpha=1.0)
            
            ## Find and set limits 
            xLims = (np.min([np.min(phasesFit),np.min(phasesReal)]),np.max([np.max(phasesFit),np.max(phasesReal)]))
            xLims = (xLims[0]-(xLims[1]-xLims[0])*.05,xLims[1]+(xLims[1]-xLims[0])*.05)
            fitYlims = (np.min([np.min(fitDataRV[:,2]*kmConversion),np.min(zeroedRealDataRV[:,5]*kmConversion)]),np.max([np.max(fitDataRV[:,2]*kmConversion),np.max(zeroedRealDataRV[:,5]*kmConversion)]))
            fitYrange = fitYlims[1]-fitYlims[0]
            fitYlims = (fitYlims[0]-fitYrange*.05,fitYlims[1]+fitYrange*.05)
            residYlims = (np.min(residualData[:,5]*kmConversion),np.max(residualData[:,5]*kmConversion))
            residYrange = residYlims[1]-residYlims[0]
            residYlims = (residYlims[0]-residYrange*.05,residYlims[1]+residYrange*.05)
            ## update lims with custom values if provided
            if len(RVlims)>0:
                #RVlims=[[yMin,yMax],[yResidMin,yResidMax],[xMin,xMax]]
                fitYlims = (RVlims[0][0],RVlims[0][1])
                residYlims = (RVlims[1][0],RVlims[1][1])
                xLims = (RVlims[2][0],RVlims[2][1])
            residualsPlot.axes.set_xlim(xLims)
            residualsPlot.axes.set_ylim(residYlims)
            #RVlims=[[yMin,yMax],[yResidMin,yResidMax],[xMin,xMax]]
            fitPlot.axes.set_xlim(xLims)
            fitPlot.axes.set_ylim(fitYlims)
            ##plot zero vel line
            residualsPlot.axhline(linewidth=1.0,c='Blue') #adds a x-axis origin line
            
            ##load resulting data to file for re-plotting by others
            #real [phases,JD,offset subtracted RV, residual, dataset#]
            outRVdataReal = []
            for i in range(0,len(phasesReal)):
                outRVdataReal.append([phasesReal[i],realDataRV[i,0],zeroedRealDataRV[i,5],residualData[i,5],realDataRV[i,7]])
            #print repr(outRVdataReal)
            #fit [phases,JD,RV]
            outRVdataFit = []
            for i in range(0,len(phasesFit)):
                outRVdataFit.append([phasesFit[i],fakeRealData[i,0],fitDataRV[i,2]])
            fnameBase = os.path.join(os.path.dirname(plotDataDir),'RVplotData')
            np.savetxt(fnameBase+'-real.dat',outRVdataReal)
            np.savetxt(fnameBase+'-fit.dat',outRVdataFit)
            for i in range(0,int(np.max(residualData[:,7])+1)):
                ary = []
                chiSqr = 0
                print ''
                for j in range(0,len(phasesReal)):
                    if residualData[j,7]==float(i):
                        if residualData[j,5] not in ary:
                            ary.append(residualData[j,5])
                            chiSqr+=(residualData[j,5]**2.0)/(residualData[j,6]**2.0)
                            if np.abs(residualData[j,5])>residualData[j,6]:
                                print 'dataset '+str(i)+', epoch '+str(realDataRV[j,0])+' had residual  '+str(np.abs(residualData[j,5])-residualData[j,6])+' greater than error'
                d = np.array(ary)
                #data = residualData[np.where(residualData[:,7]==float(i))[0],:]
                d = np.abs(d)
                print 'For RV dataset #'+str(i)
                print 'Max abs residual = '+str(np.max(d))
                print 'Min abs residual = '+str(np.min(d))
                print 'mean abs residual = '+str(np.mean(d))
                print 'chiSqr for this set was = '+str(chiSqr)
                print 'chiSqr/numRVs for this set was = '+str(chiSqr/len(d))
            
            ##clean up boarders, axis ticks and such 
            plt.minorticks_on()
            fitPlot.tick_params(axis='both',which='major',width=1,length=5,pad=8,direction='in',labelsize=20)
            fitPlot.tick_params(axis='both',which='minor',width=1,length=2,pad=8,direction='in')
            fitPlot.spines['right'].set_linewidth(1.0)
            fitPlot.spines['bottom'].set_linewidth(1.0)
            fitPlot.spines['top'].set_linewidth(1.0)
            fitPlot.spines['left'].set_linewidth(1.0)
            fitPlot.locator_params(axis='y',nbins=6) #fix number of y-axis label points
            plt.minorticks_on()
            #plot.axhline(linewidth=2.0) #adds a x-axis origin line
            #plot.axvline(linewidth=2.0) #adds a y-axis origin line
            plt.minorticks_on()
            residualsPlot.tick_params(axis='y',which='major',width=1,length=5,pad=8,direction='in',labelsize=15)
            residualsPlot.tick_params(axis='x',which='major',width=1,length=5,pad=8,direction='in',labelsize=20)
            residualsPlot.tick_params(axis='y',which='minor',width=2,length=4,pad=8,direction='in')
            residualsPlot.tick_params(axis='x',which='minor',width=2,bottom='on',length=4,pad=8,direction='in')
            residualsPlot.locator_params(axis='y',nbins=4) #fix number of y-axis label points
            plt.minorticks_on()
            residualsPlot.spines['right'].set_linewidth(1.0)
            residualsPlot.spines['bottom'].set_linewidth(1.0)
            residualsPlot.spines['top'].set_linewidth(1.0)
            residualsPlot.spines['left'].set_linewidth(1.0)
    
            ## save fig to file and maybe convert to pdf if format=='eps'
            #plt.tight_layout()
            orientStr = 'landscape'
            if format=='eps':
                orientStr = 'portrait'
            plotFilename = plotFnameBase+'-RV.'+format
            if plotFilename!='':
                plt.savefig(plotFilename, dpi=300, orientation=orientStr)
                log.info("RV orbit plot saved to:\n"+plotFilename)
            plt.close()
            if (format=='eps')and True:
                log.debug('converting to PDF as well')
                try:
                    os.system("epstopdf "+plotFilename)
                except:
                    log.warning("Seems epstopdf failed.  Check if it is installed properly.")
                ## log params used in RV plot
            log.info('\n'+"*"*50+"\nOrbital Elements used in RV plot:\n"+repr(paramsRV))
            log.info("\n with an omega value = "+str(paramsRV[9]+settingsDict["omegaFrv"][0])+'\n'+"*"*50+'\n')

def nodeEpochsCalc(paramsDI,omegaDIoffset):
    """
    Calculate the epochs for the Ascending and Descending nodes, might be in different orbital 
    periods and AN/DN might be the wrong order, but should work for plotting... I hope...
    """
    taAtNodes = [-1.0*paramsDI[9]+omegaDIoffset,180.0-1.0*paramsDI[9]+omegaDIoffset]
    nodeEpochs = []
    for ta in taAtNodes:
        if ta<0.0:
            ta =ta+360.0
        elif ta>360:
            ta =ta-360.0
        TA_s_rad = ta*(const.pi/180.0)
        top = np.sqrt(1.0-paramsDI[4])*np.sin(TA_s_rad/2.0)   
        btm = np.sqrt(1.0+paramsDI[4])*np.cos(TA_s_rad/2.0) 
        ATAN_rad = np.arctan2(top, btm)
        #NOTE: both math.atan2 and np.arctan2 tried with same results, both produce negatives rather than continuous 0-360 
        #thus, must correct for negative outputs
        if ATAN_rad<0:
            ATAN_rad = ATAN_rad+(2.0*const.pi)
        M_s_rad = ATAN_rad*2.0-paramsDI[4]*np.sin(ATAN_rad*2.0)
        delta_t = (M_s_rad*paramsDI[7]*const.daysPerYear)/(2.0*const.pi)
        nodeEpochs.append(paramsDI[5]+delta_t)
    return nodeEpochs
    

    
    
    
    
    
    
    
    
    
    
    
    
    
#END OF FILE