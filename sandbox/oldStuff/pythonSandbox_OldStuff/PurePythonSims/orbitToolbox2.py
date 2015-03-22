#@Author: Kyle Mede, kylemede@gmail.com, written for the University of Tokyo during 2011/2012.
import math
import gc
import numpy as np
import os
import pylab
from numpy import linspace
from pylab import pi

plt = pylab.matplotlib.pyplot
mpl = pylab.matplotlib.mpl

"""
This toolbox/library is a collection of functions that were used in multiple 
places throughout the code to conduct various types of binary star system simulations.
There are others in here that are for analysis of the simulator's outputs (maybe 
these will be moved to their one library in the future).
"""

def burnInCalc(paramIN):
    """
    This function will calculate the burn in length and return its value.
    
    @param paramIN:    = parameter array including burn in data
    @type paramIN:     = array (list) of doubles
    """
    
    stdALL = np.std(paramIN)
    halfStdALL = stdALL/2.0
    burnInLength = 0
    
    for i in range(0,len(paramIN)):
        stdCur = np.std(paramIN[0:i])
        if stdCur>halfStdALL:
            burnInLength = i+1
            break
    if burnInLength == len(paramIN):
        print "PROBLEM: Param had a burn in length equal to param length, ie. the chain never burned in"
        
    return burnInLength

def chiSquaredCalc(real, error, model):
    """ Just a simple function to calculate chi**2 for a given observed (real) with error
        and experimental (model) value.
        
    INPUTS:
    @param real:    = value of the observed/real parameter
    @type real:     = float
    @param error:   = error in the observed/real parameter.  
                      If the + and - errors are different, use their mean.
    @type error:    = float
    @param model:   = value of the experimental/model parameter matching the observed/real one.
    @type model:    = float
                   
    NOTES: 
        Note:   Ensure the units of both are the same.
        Note2:  For cases where the error of the observation has different errors in the 
                positive and negative from the observed value, just input the mean of these.
                
        output: Resulting chi**2 value following
        
        chi**2 = (obs - model)**2 / error**2
    """
    # since the numerator is the difference, we need to take sign of each into account
    if ((model>=0.0)and(real>=0.0))or((model<0.0)and(real<0.0)):
        # same sign so just subtract, squaring later will clear any resulting negative
        difference = model - real
    if ((model>=0.0)and(real<=0.0))or((model<0.0)and(real>0.0)):
        # different signs so add the abs of each
        difference = abs(model) + abs(real)
    # put it all together to make final chi**2 for each for this model    
    chi_squared = (math.pow(difference,2.0))/math.pow(abs(error),2.0)
    
    return chi_squared

def colorScatterPlotter(xData, yData, zData, xLabel='', yLabel='', zLabel='', plotfiletitle='', \
                                            xLim=False, yLim=False, zLim=False, save=False, show=True):
    """
    This function is to plot an X,Y scatter plot with the color of the dots representing a third 'Z' axis.
    A colorbar legend for the z-axis data will be shown vertically to the right of the scatter plot.
   
   
   
    NOTE: xData and yData must have the same length.
    
    @param xData:         = Data to plot on the x axis
    @type xData:          = list of numbers of any type, must be same length as yData and zData
    @param yData:         = Data to plot on the y axis
    @type yData:          = list of numbers of any type, must be same length as xData and zData
    @param zData:         = Data to plot on the z axis
    @type zData:          = list of numbers of any type, must be same length as xData and yData
    @param xLabel:        = label for the x axis data
    @type xLabel:         = string, default is ''
    @param yLabel:        = label for the y axis data
    @type yLabel:         = string, default is ''
    @param zLabel:        = label for the z axis data
    @type zLabel:         = string, default is ''
    @param plotFilename:  = name of the file to save the figure to
    @type plotFilename:   = string, ensuring to include directory path if not cwd
    @param xLim:          = range to limit the x axis values to
    @type xLim:           = tuple of two numbers of any type, (min,max)
                            default False indicates to use min and max of xData
    @param yLim:          = range to limit the y axis values to
    @type yLim:           = tuple of two numbers of any type, (min,max)
                            default False indicates to use min and max of yData
    @param zLim:          = range to limit the z axis values to
    @type zLim:           = tuple of two numbers of any type, (min,max)
                            default False indicates to use min and max of zData
    @param show:          = show the resulting figure on the screen?
    @type show:           = boolean, default is False
    @param save:          = save the resulting figure?
    @type save:           = boolean, default is True
   
    """ 
    # check all three data sets have matching lengths
    if not (len(xData)==len(yData)==len(zData)):
        print 'PROBLEM: The x, y and z data sets do not all have matching lengths.'    
    
    # set the xLim, yLim and zLim if their values are False
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
        range = abs(max)+abs(min)
        yLim = (min-abs(range*0.05),max+abs(range*0.05))
    else:
        if not (type(yLim)==tuple):
            print 'PROBLEM: yLim is not of type tuple'        
    if not zLim:
        min = np.min(zData)
        max = np.max(zData)
        range = abs(max)+abs(min)
        zLim = (min-abs(range*0.05),max+abs(range*0.05))
    else:
        if not (type(zLim)==tuple):
            print 'PROBLEM: zLim is not of type tuple'            
        
    # create the figure and add a subplot to it
    fig = plt.figure(1, figsize=(10,10), dpi=300)
    main = fig.add_subplot(111)
    main.set_position([0.2,0.2,0.6,0.7])
    main.axes.set_xlim(xLim)
    main.axes.set_ylim(yLim)
    main.axes.set_xlabel(xLabel)
    main.axes.set_ylabel(yLabel)
    
    # make the scatter plot
    main.scatter(xData,yData,s=2,marker='o',c=zData, edgecolors='none')
    
    # make an axes for the colorbar and add it vertically to the right of the scatter plot
    cbarAxes = fig.add_axes([0.85,0.25,0.05,0.6])
    cmap = mpl.cm.spectral
    norm = mpl.colors.Normalize(vmin=np.min(zData), vmax=np.max(zData))
    cbar = mpl.colorbar.ColorbarBase(cbarAxes, cmap=cmap, norm=norm, orientation='vertical')
    cbar.set_label(zLabel)
    
    if save:
        plt.savefig(plotFilename, dpi=300, orientation='landscape')
    if show:
        plt.show()

def confidenceLevelsFinder(INPUTSfilename, OUTPUTSfilenames, verbose=False):
    """
    
    """
    ####$$$$ change to below params being input????? $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    (plotFileTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                  Sep_Dists2, thetas2, Es2, Ms2, ns2, a1s2, a2s2, chiSquareds) = \
                            dataReader(INPUTSfilename, OUTPUTSfilenames, verbose=False)
    
    a1Means = [] #these will all be the same actually, as a_total is now an input to model
    a2Means = []#these will all be the same actually, as a_total is now an input to model
    
    # determine the number of epochs and samples in the data
    numSamples = np.shape(Es2[:][:])[0] 
    numEpochs2 = np.shape(Es2[:][:])[1]
    
    # run through all samples and load mean value into semi-majors.
    # Note: we don't need to do this for the other multi-epoch lists as they are not important
    for orbit in range(0,numSamples): 
        a1Means.append(np.mean(a1s2[orbit][:]))
        a2Means.append(np.mean(a2s2[orbit][:]))
    
    longAN_degsCLevels = CLFunc(chiSquareds,longAN_degs)
    esCLevels = CLFunc(chiSquareds,es)
    TsCLevels = CLFunc(chiSquareds,Ts)
    periodsCLevels = CLFunc(chiSquareds,periods)
    inclination_degsCLevels= CLFunc(chiSquareds,inclination_degs)
    argPeri_degsCLevels = CLFunc(chiSquareds,argPeri_degs)
    a1sCLevels = CLFunc(chiSquareds,a1Means)
    a2sCLevels = CLFunc(chiSquareds,a2Means)
    
    return (longAN_degsCLevels, esCLevels, TsCLevels, periodsCLevels, inclination_degsCLevels, \
                                    argPeri_degsCLevels, a1sCLevels, a2sCLevels)
    
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
    
    paramCLevels.append([np.min(param68s),np.max(param68s)])
    
    ## find all the parameter values that are inside the 95.5% confidence level
    param95s = []
    for orbit in range(0,numSamples):
        if chiSquareds[orbit]<=(chiSquareMin+4.0):
            param95s.append(param[orbit])
    
    paramCLevels.append([np.min(param95s),np.max(param95s)])
    
    return paramCLevels

def ConfLevelXYs(CLevels,n,bins,patches):
    """
    return  [[yminTop68.3%,ymaxTop68.3%],[yminTop95.5%,ymaxTop95.5%]], bins68_and_95
    """
    
    xmin68 = CLevels[0][0]
    xmax68 = CLevels[0][1]
    xmin95 = CLevels[1][0]
    xmax95 = CLevels[1][1]

    # first lets get the Xs and Ys for the 68.3 and 95.5% vertical lines
    #68% Ys and bins
    minFound = False
    maxFound = False
    for i in range(0,len(bins)):
        if (bins[i]>xmin68) and (minFound==False):
            min68bin = i
            yminTop68 = n[i]
            minFound = True
        if (bins[i]<xmax68) and (maxFound==False):
            max68bin = i
            ymaxTop68 = n[i]
            maxFound = True
    #95% Ys and bins
    minFound = False
    maxFound = False
    for i in range(0,len(bins)):
        if (bins[i]>xmin95) and (minFound==False):
            min95bin = i
            yminTop95 = n[i]
            minFound = True
        if (bins[i]<xmax95) and (maxFound==False):
            max95bin = i
            ymaxTop95 = n[i]
            maxFound = True
    
    Ys = [[yminTop68, ymaxTop68],[yminTop95, ymaxTop95]]
    
    # define the float colors matching the 'mpl color spec'
    lightGrey = '0.4'
    darkGrey = '0.7'
    

    # make list of patches of all patches in 95% region, which includes 68% region
    bins68_and_95 = []
    for bin in range(min95bin,max95bin+1):
        bins68_and_95.append(patches[bin].set_color('g'))
    # update color of patches in 68% (middle) section
    for bin in range(min68bin,max68bin+1):
        bins68_and_95[bin].set_color('y')
    
    

    return (Ys, bins68_and_95)

def daisyChainDataFileCombiner(filenameROOT, numRounds, numEpochs, outFilenameROOT, removeIndividuals=False):
    """
    This is a version of the dataFileCombiner to be used after running a 'daisychain' 
    simulation.  It will put all the INS and OUTS files for each round together into
    final INS and OUTS files as if only a single large round was ran. 

    NOTE: The files may either be in the current working directory, or nested.  If not in the 
    cwd, the path must be included in the filenameROOT and outFilenameROOT parameters.    

    @param filenameROOT:        = root of the filenames, both INS and OUTS.
                                    ex. 'mcmcMultiprocessUniform-100million-Round-1_INS.txt'
                                  has the root 'mcmcMultiprocessUniform-100million-',
                                  thus ignore all from 'Round' onward.
    @type filenameROOT:            = string
    @param numRounds:            = number of rounds in the daisy chain
    @type  numRounds:            = int
    @param numEpochs:            = number of epochs in the OUTS
    @type numEpochs:            = int
    @param outFilenameROOT:     = root of output filenames, following same rules as 
                                  filenameROOT ex.  ie put 'mcmcMultiprocessUniform-6billion-'
    @type outFilenameROOT:        = string
    """

    if '.txt' in filenameROOT:
        print 'PROBLEM: there is .txt in the filenameROOT and there should not be'
    
    if '.txt' in outFilenameROOT:
        print 'PROBLEM: there is .txt in the outFilenameROOT and there should not be'
    
    # create the INS file names in a list and combine
    INS = []
    for round in range(1,numRounds+1):
        filename = filenameROOT+'Round-'+str(round)+'_INS.txt'
        INS.append(filename)
    INSoutname = outFilenameROOT+'_INS.txt'
    dataFileCombiner(INS,INSoutname)

    # creat the OUTS file names in a list and combine for each epoch
    for epoch in range(1,numEpochs+1):
        OUTSepoch = []
        for round in range(1,numRounds+1):
            filename = filenameROOT+'Round-'+str(round)+'_OUTSepoch'+str(epoch)+'.txt'
            OUTSepoch.append(filename)
        OUTSepochName = outFilenameROOT+'_OUTSepoch'+str(epoch)+'.txt'
        dataFileCombiner(OUTSepoch, OUTSepochName)
        
#        if removeIndividuals:
#            for file in OUTSepoch:
#                if os.path.exists(file):
#                    os.remove(file)
                    
    print 'DONE combining all INS and OUTS for that daisy chain'


    
def dataReadAndPlot(filenameROOT, numEpochs, chiSquareMaxGooders=10.0, plotFilename='',lowRAM=True, \
                                                            confLevels=True, weight=False, normed=True):
    """
    This function will load the files in based on the filenameROOT and numEpochs
    into memory, then make a summary plot of the parameters for those orbits with a chiSquared
    less than chiSquareMaxGooders.  This will load all the INS and OUTS for each epoch.
    
    @param filenameROOT:         = root of the filenames, both INS and OUTS.
                                    ex. 'mcmcMultiprocessUniform-100million_OUTSepoch1.txt'
                                   has the root 'mcmcMultiprocessUniform-100million',
                                   thus ignore all from 'Round' onward.
    @type filenameROOT:          = string
    @param numEpochs:            = number of epochs in the OUTS
    @type numEpochs:             = int
    @param chiSquareMaxGooders:  = max value of chiSquared for orbits to be returned
    @type chiSquareMaxGooders:   = any number, default is 10.0
    @param plotFilename:         = file name for the output plot
    @type plotFilename:          = string, including directory if not cwd
    """
    
    # check if the passed in value for filenameROOT includes '.txt'
    if '.txt' in filenameROOT:
        print 'PROBLEM: there is .txt in the filenameROOT and there should not be'
    
    # check if the passed in value for plotFilename includes '.png'
    if '.png' not in plotFilename:
        plotFilename = plotFilename+'.png'
    
    plotFilename2 = plotFilename[:-4]+'-ChiSquareMax_'+str(chiSquareMaxGooders)+'.png'
    
    print '\nLoading data into memory and finding orbits with chiSquare < '+str(chiSquareMaxGooders)
    
    if lowRAM:
        (chiSquaredGooders, inclination_degsAllGooders, longAN_degsAllGooders, argPeri_degsAllGooders, esAllGooders,\
            periodsAllGooders, TsAllGooders, a1Means, a2Means, aTotMeans) = \
            dataReadAndReturn(filenameROOT, numEpochs, chiSquareMaxGooders, lowRAM=lowRAM)
    else:
        (chiSquaredGooders, inclination_degsAllGooders, longAN_degsAllGooders, argPeri_degsAllGooders, esAllGooders,\
            periodsAllGooders, TsAllGooders, a1Means, a2Means, aTotMeans) = \
            dataReadAndReturn(filenameROOT, numEpochs, chiSquareMaxGooders,lowRAM=lowRAM)
    
    # use the likelihoods function to convert the chiSquareds to likelihoods for weighting
    likelihoods = likelihoodsCalc(chiSquaredGooders) 
    if len(chiSquaredGooders)==len(inclination_degsAllGooders):#####$$$$$$$$$
        print 'Lengths are the same!!'##$$$$$$$$$$$$$$$$$
        
    print 'Finished loading data, starting to make summary plot'

    summaryPlotter(plotFilename2, chiSquaredGooders, inclination_degsAllGooders, longAN_degsAllGooders, argPeri_degsAllGooders, \
                   esAllGooders, periodsAllGooders, TsAllGooders, a1Means, a2Means, confLevels=confLevels, weight=weight, normed=normed,\
                    showPlots=False, save=True, verbose=True)
    
    print 'Summary plot complete and written to disk '+plotFilename2
   
def dataReadAndPlot2(filenameROOT, numEpochs, chiSquareMaxGooders=10.0, plotFilename='', weight=False):
    """
    This function will load the files in based on the filenameROOT and numEpochs
    into memory, then make a summary plot of the parameters for those orbits with a chiSquared
    less than chiSquareMaxGooders.  This will load all the INS and OUTS for each epoch.
    
    @param filenameROOT:         = root of the filenames, both INS and OUTS.
                                    ex. 'mcmcMultiprocessUniform-100million_OUTSepoch1.txt'
                                   has the root 'mcmcMultiprocessUniform-100million',
                                   thus ignore all from 'Round' onward.
    @type filenameROOT:          = string
    @param numEpochs:            = number of epochs in the OUTS
    @type numEpochs:             = int
    @param chiSquareMaxGooders:  = max value of chiSquared for orbits to be returned
    @type chiSquareMaxGooders:   = any number, default is 10.0
    @param plotFilename:         = file name for the output plot
    @type plotFilename:          = string, including directory if not cwd
    """
    
    # check if the passed in value for filenameROOT includes '.txt'
    if '.txt' in filenameROOT:
        print 'PROBLEM: there is .txt in the filenameROOT and there should not be'
    
    # check if the passed in value for plotFilename includes '.png'
    if '.png' not in plotFilename:
        plotFilename = plotFilename+'.png'
    
    plotFilename2 = plotFilename[:-4]+'-ChiSquareMax_'+str(int(chiSquareMaxGooders))+'.png'
    
    print '\nLoading data into memory and finding orbits with chiSquare < '+str(chiSquareMaxGooders)
    
    (chiSquaredGooders, inclination_degsAllGooders, longAN_degsAllGooders, argPeri_degsAllGooders, esAllGooders,\
            periodsAllGooders, TsAllGooders, a1Means, a2Means, aTotMeans) = \
            dataReadAndReturn(filenameROOT, numEpochs, chiSquareMaxGooders)
    
    print 'Finished loading data, starting to make summary plot'

    summaryPlotter(plotFilename2, chiSquaredGooders, inclination_degsAllGooders, longAN_degsAllGooders, argPeri_degsAllGooders, \
                   esAllGooders, periodsAllGooders, TsAllGooders, a1Means, a2Means, showPlots=False, \
                   save=True, verbose=True)
    
    print 'Summary plot complete and written to disk '+plotFilename2

def dataReader(INPUTSfilename, OUTPUTSfilenames, verbose=False):
    """ 
    This function is to read data files written to disk by the dataWriter function back
    into memory for plotting or statistical analysis. 
    
    @param INPUTSfilename:     = file name of the file containing the 'inputs' file 
    @type INPUTSfilename:      = single string for the filename or a list of filename strings
    @param OUTPUTSfilenames:   = list containing the file names of the 'outputs' files
    @type OUTPUTSfilenames:    = list of filename strings. Order of list should be epoch1, epoch2,...
        
    NOTE:
        file format must be:
        
        line1: title
        line2: data headings
        line3: space delimited data
        line4: space delimited data
        .
        .
        .
    """
    ### INPUTS #####
    if type(INPUTSfilename)==list:
        inFilename = INPUTSfilename[0]
    else:
        inFilename = INPUTSfilename
    if verbose:
        print 'Starting to load the input data'
    INfile = open(inFilename, 'r')
    lines = INfile.readlines()
    INtitle = lines[0]
    plotFileTitle = INtitle[:-5]
    INcolHeaders = lines[1]
    numLines = len(lines)
    
    # initiate input data lists
    longAN_degs = []
    es = []
    Ts = []
    periods = []
    inclination_degs = []
    argPeri_degs = []
    chiSquareds = []
    
    # loop through all the lines and extract data into data lists
    for row in range(2,numLines):
        dataLine = lines[row]
        dataLineCols = dataLine.split()
        dataLineColsFloats = []
        for data in dataLineCols:
            dataLineColsFloats.append(float(data))
        # put resulting float values of each column into its proper list
        longAN_degs.append(dataLineColsFloats[0])
        es.append(dataLineColsFloats[1])
        Ts.append(dataLineColsFloats[2])
        periods.append(dataLineColsFloats[3])
        inclination_degs.append(dataLineColsFloats[4])
        argPeri_degs.append(dataLineColsFloats[5])
        chiSquareds.append(dataLineColsFloats[6])
    if verbose:
        print 'Done loading Inputs from file'
    INfile.close()
    #################################
    
    # OUTPUTS    
    if type(OUTPUTSfilenames)==list:
        # initiate output data lists
        Sep_Dists2 = []
        thetas2 = []
        Es2 = []
        Ms2 = []
        ns2 = []
        a1s2 = []
        a2s2 = []
        
        for file in range(0,len(OUTPUTSfilenames)):
            epoch = file+1
            curOutFilename = OUTPUTSfilenames[file]
            OUTfile = open(curOutFilename, 'r')
            lines = OUTfile.readlines()
            OUTtitle = lines[0]
            OUTcolHeaders = lines[1]
            numLines = len(lines)
    
            # initiate output data lists
            Sep_Dists = []
            thetas = []
            Es = []
            Ms = []
            ns = []
            a1s = []
            a2s = []
            
            if file == 0:
                if verbose:
                    print 'Loading up output data of first epoch'
                # load first epoch's data into the lists
                # then update/append to these lists with later epochs.
                # loop through all the lines and extract data into data lists
                for row in range(2,numLines):
                    dataLine = lines[row]
                    dataLineCols = dataLine.split()
                    dataLineColsFloats = []
                    for data in dataLineCols:
                        dataLineColsFloats.append(float(data))
        
                    Sep_Dists2.append([dataLineColsFloats[0]])
                    thetas2.append([dataLineColsFloats[1]])
                    Es2.append([dataLineColsFloats[2]])
                    Ms2.append([dataLineColsFloats[3]])
                    ns2.append([dataLineColsFloats[4]])
                    a1s2.append([dataLineColsFloats[5]])
                    a2s2.append([dataLineColsFloats[6]])
            elif file>0:
                if verbose:
                    print 'loading up output data of epoch ',epoch 
                for row in range(2,numLines):
                    dataLine = lines[row]
                    dataLineCols = dataLine.split()
                    dataLineColsFloats = []
                    for data in dataLineCols:
                        dataLineColsFloats.append(float(data))
                    
                    # update/append lists with this epoch's data
                    Sep_Dists2[row-2].append(dataLineColsFloats[0])
                    thetas2[row-2].append(dataLineColsFloats[1])
                    Es2[row-2].append(dataLineColsFloats[2])
                    Ms2[row-2].append(dataLineColsFloats[3])
                    ns2[row-2].append(dataLineColsFloats[4])
                    a1s2[row-2].append(dataLineColsFloats[5])
                    a2s2[row-2].append(dataLineColsFloats[6])
            if verbose:    
                print 'Done loading outputs for epoch '+str(epoch)+' from file'
            OUTfile.close()
        if verbose:    
            print '** All data from files loaded **'
        
    else:
        print 'OUTPUTSfilenames must be a list of output file name strings'
        
    return(plotFileTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                  Sep_Dists2, thetas2, Es2, Ms2, ns2, a1s2, a2s2, chiSquareds)

def dataReadAndReturn(filenameROOT, numEpochs, chiSquareMaxGooders=10.0, lowRAM=True): 
    """
    This function will load the files in based on the filenameROOT and numEpochs
    into memory, then return the parameters for those orbits with a chiSquared
    less than chiSquareMaxGooders.  This will load all the INS and OUTS for each epoch.
    
    @param filenameROOT:         = root of the filenames, both INS and OUTS.
                                    ex. 'mcmcMultiprocessUniform-100million_OUTSepoch1.txt'
                                   has the root 'mcmcMultiprocessUniform-100million',
                                   thus ignore all from 'Round' onward.
    @type filenameROOT:          = string
    @param numEpochs:            = number of epochs in the OUTS
    @type numEpochs:             = int
    @param chiSquareMaxGooders:  = max value of chiSquared for orbits to be returned
    @type chiSquareMaxGooders:   = any number, default is 10.0
    """
    # check if the passed in value for filenameROOT includes '.txt'
    if '.txt' in filenameROOT:
        print 'PROBLEM: there is .txt in the filenameROOT and there should not be'
    
    
    # create the INS file name
    INSfilename = filenameROOT+'_INS.txt'
    
    # creat the OUTS file names in a list
    OUTS = []
    for epoch in range(1,numEpochs+1):
        filename = filenameROOT+'_OUTSepoch'+str(epoch)+'.txt'
        OUTS.append(filename)
    
    if lowRAM:
        (plotFileTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                                 a1s2, a2s2, chiSquareds) = dataReaderLowRAM(INSfilename, OUTS)
    else:
        (plotFileTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
            Sep_Dists2, thetas2, Es2, Ms2, ns2, a1s2, a2s2, chiSquareds) = dataReader(INSfilename, OUTS)
                  
    # init lists to load up and return
    chiSquaredGooders = []
    inclination_degsAllGooders = []
    longAN_degsAllGooders = []
    argPeri_degsAllGooders = []
    esAllGooders = []
    periodsAllGooders = []
    TsAllGooders = []
    a1AllGooders = []
    a2AllGooders = []
    a1Means = []
    a2Means = []
    aTotMeans = []
    
    # determine the number of epochs and samples in the data
    numSamples = np.shape(a2s2[:][:])[0] 
    numEpochs2 = np.shape(a2s2[:][:])[1]
    
    # run through all samples and check the orbits chiSquared
    # if it is lower than chiSquaredMaxGooders, store the orbit
    # else, skip it
    minChiSquared = np.min(chiSquareds)
    for orbit in range(0,numSamples):
        if chiSquareds[orbit]<=chiSquareMaxGooders:
            chiSquaredGooders.append(chiSquareds[orbit])
            inclination_degsAllGooders.append(inclination_degs[orbit])
            esAllGooders.append(es[orbit])
            longAN_degsAllGooders.append(longAN_degs[orbit])
            periodsAllGooders.append(periods[orbit])
            TsAllGooders.append(Ts[orbit])
            argPeri_degsAllGooders.append(argPeri_degs[orbit])
            a1Mean = np.mean(a1s2[orbit][:])
            a1Means.append(a1Mean)
            a2Mean = np.mean(a2s2[orbit][:])
            a2Means.append(a2Mean)
            a_total = np.mean(a1Mean+a2Mean)
            aTotMeans.append(a_total)
            if chiSquareds[orbit]==minChiSquared:
                print '*'*25
                print 'Orbit parameters for lowest chiSquared of '+str(minChiSquared)+ ' are:'
                print 'Inclination [deg]: ',inclination_degs[orbit]
                print 'Eccentricity : ',es[orbit]
                print 'Longitude of Ascending Node [deg] : ',longAN_degs[orbit]
                print 'Period [years]: ',periods[orbit]
                print 'Time of Last Periapsis [jd]: ',Ts[orbit]
                print 'Argument of Perigee [deg]: ',argPeri_degs[orbit]
                print 'Total Semi-major Axis [AU]: ',a_total
                print '*'*25
                
    # check if any orbits were found with chiSquares less than chiSquaredMaxGooders
    if len(aTotMeans)>0:
        print str(len(aTotMeans))+' orbits were found with a chiSquared less than '+str(chiSquareMaxGooders)
    else:
        print 'No orbits were found with a chiSquared less than '+str(chiSquareMaxGooders)
    
            
    return(chiSquaredGooders,inclination_degsAllGooders, longAN_degsAllGooders, argPeri_degsAllGooders, esAllGooders,\
            periodsAllGooders, TsAllGooders, a1Means, a2Means, aTotMeans)        

def dataReadTrimWrite(filenameROOT, numEpochs, chiSquareCutOff, verbose=False):
    """
    This function is ideal for the cases where your data sets are too big because of choosing a chiSquareMax that is
    too big when running the simulator, resulting in large data files that are hard to load and plot.  In these cases
    you can just find out what the standard range of chiSquareds are, and choose a lower cut-off to trim the volume
    of data written to disk.  Just open the 'INS' file and skim over the chiSquared column values to choose a new
    cut-off.
    These circumstances really only happen when you have a long 'daisy chain' run with rough settings, so this is 
    not expected to be a commonly used function.
    
    NOTE: output files will have same name with 'chiSquare-cut-off-####' added before '_INS' or 'OUTS..' to 
          show they are the new trimmed versions.
    
    @param filenameROOT:         = root of the filenames, both INS and OUTS.
                                    ex. 'mcmcMultiprocessUniform-100million_OUTSepoch1.txt'
                                   has the root 'mcmcMultiprocessUniform-100million',
                                   thus ignore all from 'Round' onward.
    @type filenameROOT:          = string
    @param numEpochs:            = number of epochs in the OUTS
    @type numEpochs:             = int
    @param chiSquareCutOff:      = max value of chiSquared for orbits to be written
    @type chiSquareCutOff:       = any number
    """
    
    # check if the passed in value for filenameROOT includes '.txt'
    if '.txt' in filenameROOT:
        print 'PROBLEM: there is .txt in the filenameROOT and there should not be'
    
    # create the INS file name
    INSfilename = filenameROOT+'_INS.txt'
    
    # create the OUTS file names in a list
    OUTS = []
    for epoch in range(1,numEpochs+1):
        filename = filenameROOT+'_OUTSepoch'+str(epoch)+'.txt'
        OUTS.append(filename)


    titleMod = '-chiSquare-cut-off-'+str(chiSquareCutOff)
    ### INPUTS #####
    if verbose:
        print 'Starting to load, trim and write the INS data'
    # open INS file and grab title and headings
    INfile = open(INSfilename, 'r')
    plotFileTitle = INfile.readline()[:-5]
    headings = INfile.readline()
    # create output version of INS file and write new 
    # strip off '_INS...txt', add titleMod, then put it back on
    newINSfilename = INSfilename[:-8]+titleMod+INSfilename[-8:]
    INOUTfile = open(newINSfilename, 'w')
    INOUTfile.write(plotFileTitle+titleMod+'\n')
    INOUTfile.write(headings)
    # loop through data lines and only write those that
    # have chiSquare below cut-off and store all 
    # chiSquareds in a list
    chiSquareds = []
    line='asdf'
    totalKept = 0
    while line!='':
        line = INfile.readline()
        if line!='':
            dataLineCols = line.split()
            chiSquareStr = dataLineCols[6]
            chiSquareds.append(float(chiSquareStr))
            if chiSquareds[-1]<chiSquareCutOff:
                INOUTfile.write(line)
                totalKept = totalKept+1
    if verbose:
        print str(totalKept)+' orbits were found to have a chiSquared < '+str(chiSquareCutOff)
        print 'Done loading, trimming and writing INS'
    INfile.close()
    INOUTfile.close()
    #################################
    
    # OUTPUTS    
    epoch = 1
    if verbose:
        print 'Starting to load and trim the OUTS data'
    for curOutFilename in OUTS:
        # open OUTS file and grab title and headings
        OUTfile = open(curOutFilename, 'r')
        plotFileTitle = OUTfile.readline()[:-5]
        headings = OUTfile.readline()
        # create output version of OUTS file and write new
        # strip off '_OUTS...txt', add titleMod, then put it back on
        if epoch<10:
            newOutFilename = curOutFilename[:-15]+titleMod+curOutFilename[-15:]
        else:
            newOutFilename = curOutFilename[:-16]+titleMod+curOutFilename[-16:]
        OUTOUTfile = open(newOutFilename, 'w')
        OUTOUTfile.write(plotFileTitle+titleMod+'\n')
        OUTOUTfile.write(headings)
        # loop through the data lines and write them to new 
        # file if chiSquared is less than cut-off
        i = 0
        line='asdf'
        while line!='':
            line = OUTfile.readline()
            if line!='':
                if chiSquareds[i]<chiSquareCutOff:
                    OUTOUTfile.write(line)
                i=i+1
            
        if verbose:    
            print 'Done loading, trimming and writing outputs for epoch '+str(epoch)
        OUTfile.close()
        OUTOUTfile.close()
        epoch =epoch+1
    if verbose:    
        print '** All data from files loaded **'

def dataWriter(plotFileTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                  ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds):
    
    """
    This function is to write the resulting inputs and outputs of the mcmcOrbSimulatorUniform2.
    The inputs and outputs for each epoch will all be written to separate files.
    
    The function dataReader can later be used to reload these files back into memory for later
    data analysis or plotting.
    
    CAUTION! the saved data will overwrite a previous version of the file if it has the same name!!
    CAUTION!!! currently the directory the file is written to is hardcoded!!! ie. code/data/~~~ $$$ FIX THIS $$$ maybe accept a dataPath parameter
    
    INPUTS:
    @param plotFileTitle:                 = initial part of the file title.  This will be 
                                            have strings appended to it to create individual
                                            input and output file names.
    @type plotFileTitle:                  = string
    @param longAN_degs:                   = Longitude of Ascending Node [deg]
    @type longAN_degs:                    = list of floats
    @param es:                            = eccentricity of orbits [unitless]
    @type es:                             = list of floats
    @param Ts:                            = Last Periapsis Epoch/time [julian date] 
    @type Ts:                             = list of floats
    @param periods:                       = period of orbits [yrs]
    @type periods:                        = list of floats
    @param inclination_degs:              = inclination [deg]
    @type inclination_degs:               = list of floats
    @param argPeri_degs:                  = Argument of Periapsis in orbital plane [deg]
    @type argPeri_degs:                   = list of floats
    
    @param ns:                            = Mean Motion [rad/yr]
    @type ns:                             = list of floats. Output values for each epoch.
    @param Ms:                            = Mean Anomaly [deg]
    @type Ms:                             = list of floats. Output values for each epoch.
    @param Es:                            = Eccentric Anomaly [deg]
    @type Es:                             = list of floats. Output values for each epoch.
    @param thetas:                        = True Anomaly [deg]
    @type thetas:                         = list of floats. Output values for each epoch.
    @param Sep_Dists:                     = Separation Distance in orbital plane [AU]
    @type Sep_Dists:                      = list of floats. Output values for each epoch.
    @param a1s:                           = semi-major axis of M1 [AU]
    @type a1s:                            = list of floats. Output values for each epoch.
    @param a2s:                           = semi-major axis of M2 [AU]
    @type a2s:                            = list of floats. Output values for each epoch.
    @param chiSquareds:                   = total chi squared values for each orbit.
    @type chiSquareds:                    = list of floats
    
    """
    # find how many samples and epochs are to be written to files
    numSamples = np.shape(a1s2[:][:])[0]
    numEpochs = np.shape(a1s2[:][:])[1]
    
    ############## INS ###################
    filename = '../data/'+plotFileTitle+'_INS'
    print "!!!!!!!data file written to ",filename #$$$$$$$$$$$$$$$
    f = open('%s.txt' % filename, 'w')
    f.write(plotFileTitle+'_INS\n')
    f.write('longAN [deg]    es [N/A]    Ts [julian date]    periods [yrs]   inclination [deg]   argPeri [deg]   chiSquareds\n')
    for sample in range(0,numSamples):
        line = str(longAN_degs[sample])
        line = line +'   '+ str(es[sample])
        line = line +'   '+ str(Ts[sample])
        line = line +'   '+ str(periods[sample])
        line = line +'   '+ str(inclination_degs[sample])
        line = line +'   '+ str(argPeri_degs[sample])
        line = line +'   '+ str(chiSquareds[sample])+'\n'
        f.write(line)
    f.close()
    #######################################
    
    ################# OUTS ####################
    for epoch in range(0,numEpochs):
        filename = '../data/'+plotFileTitle+'_OUTSepoch'+str(epoch+1)
        f = open('%s.txt' % filename, 'w')
        f.write(plotFileTitle+'_OUTSepoch'+str(epoch+1)+'\n')
        f.write('Sep_Dists [AU]    thetas [deg]     Es [deg]     Ms [deg]     ns [rad/yr]  a1s [AU]   a2s [AU]\n')
        for sample in range(0,numSamples):
            line = str(Sep_Dists2[sample][epoch])
            line = line +' '+ str(thetas2[sample][epoch])
            line = line +' '+ str(Es2[sample][epoch])
            line = line +' '+ str(Ms2[sample][epoch])
            line = line +' '+ str(ns2[sample][epoch])
            line = line +' '+ str(a1s2[sample][epoch])
            line = line +' '+ str(a2s2[sample][epoch])+'\n'
            f.write(line)
        f.close()
    

def rv2bodyCalculator(RV_epochs, RVs, RVerrors, sigma_jitter, i_c, p_c, e_c, T_c, argPeri_c, a_c, K_p, p_p, e_p, argPeri_p, T_p):
    """
    This will calclulate the radial velocity residue based chiSquared considering RV data for the primary star, orbital elements 
    for the companion star (indicated with '_c' parameters), and some of the parameters for the planet around the primary
    star's orbital elements (indicated with the '_p' parameters).
    
    sigma_jitter must be the unsquared version and in units of [m/s]
    
    i in degrees
    argPeri in degrees
    a in AU
    T in JD
    t in JD
    e unitless
    period in years
    """
    
    chiSquaredTotal = 0.0
    
    for epoch in range(0,len(RVs)):
        # calculate the velocity residual due to the companion star
        v_r_c = vrCalculator(RV_epochs[epoch],e_c, T_c, p_c, a_c, argPeri_c, i=i_c, K=False)
        #print '\nv_r_c = ',v_r_c###$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        # calculate the velocity residual due to the planet around primary
        v_r_p = vrCalculator(RV_epochs[epoch], e_p, T_p, p_p, 0, argPeri_p, i=False, K=K_p)
        #print 'v_r_p = ',v_r_p###$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
        #print 'RVs[epoch] = ',RVs[epoch]
        #print '(RVerrors[epoch]+sigma_jitter) = '+str((RVerrors[epoch]+sigma_jitter))
        # calculate chiSquared for this epoch
        chiSquaredCurr = chiSquaredCalc(RVs[epoch], (RVerrors[epoch]+sigma_jitter), (v_r_c+v_r_p))
        
        chiSquaredTotal = chiSquaredTotal + chiSquaredCurr
        
    return chiSquaredTotal

def rv2bodyCalculator2(RV_epochs, RVs, RVerrors, sigma_jitter, M1, M2SineI_c, p_c, e_c, T_c, argPeri_c, \
                                                                    K_p, p_p, e_p, argPeri_p, T_p):
    """
    This will calclulate the radial velocity residue based chiSquared considering RV data for the primary star, orbital elements 
    for the companion star (indicated with '_c' parameters), and some of the parameters for the planet around the primary
    star's orbital elements (indicated with the '_p' parameters).
    
    sigma_jitter must be the unsquared version and in units of [m/s]
    M2SineI_c and M1 in kg
    i in degrees
    argPeri in degrees
    a in AU
    T in JD
    t in JD
    e unitless
    period in years
    """
    
    chiSquaredTotal = 0.0
    
    for epoch in range(0,len(RVs)):
        # calculate the velocity residual due to the companion star
        v_r_c = vrCalculator2(RV_epochs[epoch],e_c, T_c, p_c, argPeri_c, M1, M2SineI=M2SineI_c, K=False)
        #print '\nv_r_c = ',v_r_c###$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        # calculate the velocity residual due to the planet around primary
        v_r_p = vrCalculator2(RV_epochs[epoch], e_p, T_p, p_p, argPeri_p, M1, M2SineI=False, K=K_p)
        #print 'v_r_p = ',v_r_p###$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
        #print 'RVs[epoch] = ',RVs[epoch]
        #print '(RVerrors[epoch]+sigma_jitter) = '+str((RVerrors[epoch]+sigma_jitter))
        # calculate chiSquared for this epoch
        chiSquaredCurr = chiSquaredCalc(RVs[epoch], (RVerrors[epoch]+sigma_jitter), (v_r_c+v_r_p))
        
        chiSquaredTotal = chiSquaredTotal + chiSquaredCurr
        
    return chiSquaredTotal

def vrCalculator(t,e, T, period, a, argPeri, i=False, K=False):
    """
    i in degrees
    argPeri in degrees
    a in AU
    T in JD
    t in JD
    e unitless
    period in years
    
    K=False implies that the provided values will be used to 
    calculate it from scratch, else the provided value will
    be used to calculate the radial velocity residule.  Thus,
    this param must be either False or a float in units of 
    m/s.
    
    returned vel will be in m/s
    """
    if (K==False) and (i==False):
        print 'vrCalculator: The value of i and K cannot both be False, one MUST be defined'
    if (K==False) and (i!=False):
        # convert units of years to seconds
        period_seconds = period*31557600.0
        # convert units of AU to meters
        a_meters = a*149598000000.0
        # calculate the K (semi amplitude) for the provided values
        K = (2.0*pi*a_meters*math.sin(math.radians(i)))/(period_seconds*math.sqrt(1.0-e**2)) 
    
    ## get the TAcalculator to find the TA in radians
    (n, M_deg, E_latest_deg,TA_rad) = TAcalculator(t,e, T, period, verbose=False)
    
    v_r = K*(math.cos(math.radians(argPeri)+TA_rad)+e*math.cos(math.radians(argPeri)))

    return v_r

def vrCalculator2(t,e,T,period,argPeri,M1,M2SineI=False, K=False):
    """
    M2SineI and M1 in kg
    argPeri in degrees
    a in AU
    T in JD
    t in JD
    e unitless
    period in years
    
    K=False implies that the provided values will be used to 
    calculate it from scratch, else the provided value will
    be used to calculate the radial velocity residule.  Thus,
    this param must be either False or a float in units of 
    m/s.
    
    returned vel will be in m/s
    """
        
    if (K==False) and (M2SineI==False):
        print 'vrCalculator: The value of i and K cannot both be False, one MUST be defined'
    if (K==False) and (M2SineI!=False):
        # convert units of years to seconds
        period_seconds = period*31557600.0
        mass_total = M1 + M2SineI
        K = ((2.0*pi*6.67e-11)/period_seconds)**(1.0/3.0)*((M2SineI)/(math.sqrt(1.0-e**2)*(mass_total)**(2.0/3.0)))

    ## get the TAcalculator to find the TA in radians
    (n, M_deg, E_latest_deg,TA_rad) = TAcalculator(t,e, T, period, verbose=False)
    
    v_r = K*(math.cos(math.radians(argPeri)+TA_rad)+e*math.cos(math.radians(argPeri)))

    return v_r


def histAndScatterPlotter(xData, yData, xLabel='', yLabel='', plotFilename='', xLim=False, yLim=False, show=True, save=False):
    """
    This function will create a plot with a scatter plot of the xData and yData in the middle, a horizontal histogram 
    of the yData on the right, and vertical histogram of the xData at the top.  This is a good way to compare the 
    parameter space where the two data sets have their values focus.
    
    NOTE: xData and yData must have the same length.
    
    @param xData:         = Data to plot on the x axis
    @type xData:          = list of numbers of any type, must be same length as yData
    @param yData:         = Data to plot on the y axis
    @type yData:          = list of numbers of any type, must be same length as xData
    @param xLabel:        = label for the x axis data
    @type xLabel:         = string, default is ''
    @param yLabel:        = label for the y axis data
    @type yLabel:         = string, default is ''
    @param plotFilename:  = name of the file to save the figure to
    @type plotFilename:   = string, ensuring to include directory path if not cwd
    @param xLim:          = range to limit the x axis values to
    @type xLim:           = tuple of two numbers of any type, (min,max)
                            default False indicates to use min and max of xData
    @param yLim:          = range to limit the y axis values to
    @type yLim:           = tuple of two numbers of any type, (min,max)
                            default False indicates to use min and max of yData
    @param show:          = show the resulting figure on the screen?
    @type show:           = boolean, default is False
    @param save:          = save the resulting figure?
    @type save:           = boolean, default is True
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
    X.set_title(plotFilename[:-4])
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

def multiEpochOrbCalc(Sep_Angle_arcsec_measured_REALs, PA_deg_measured_REALs, PA_mean_errors, epochs,\
                       chiSquaredMax, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, Mass1=1, Mass2=1, verbose=False):
    """
    This function will check if the suggested orbital parameter inputs satisfy all the input Position angles with a total 
    chi squared less than chiSquaredMax.  Thus, it will check if they match all the observations of the system.
    
    
    Inputs:
    @param Sep_Angle_arcsec_measured_REALs:  = measured separation angle ["]
    @type Sep_Angle_arcsec_measured_REALs:   = List of floats
    @param PA_deg_measured_REALs:            = measured position angle [deg]
    @type PA_deg_measured_REALs:             = List of floats
    @param PA_mean_errors:                   = error in measured position angle [deg]. Use average if + and - are different.
    @type PA_mean_errors:                    = List of floats
    @param epochs:                           = epoch of observation/image [julian date]
    @type epochs:                            = List of floats
    @param chiSquaredMax:                    = measured system distance from Earth [PC]
    @type chiSquaredMax:                     = float
    @param Sys_Dist_PC:                      = measured system distance from Earth [PC]
    @type Sys_Dist_PC:                       = float
    @param inclination_deg:                  = inclination [deg]
    @type inclination_deg:                   = float
    @param longAN_deg:                       = Longitude of Ascending Node [deg]
    @type longAN_deg:                        = float
    @param e:                                = eccentricity of orbits [unitless]
    @type e:                                 = float
    @param T:                                = Last Periapsis Epoch/time [julian date] 
    @type T:                                 = float
    @param period:                           = period of orbits [yrs]
    @type period:                            = float
    @param argPeri_deg:                      = Argument of Periapsis in orbital plane [deg]
    @type argPeri_deg:                       = float
    @param Mass1:                            = Mass of primary in system [solar masses] (Default=1, means M2 is a planet)
    @type Mass1:                             = float
    @param Mass2:                            = Mass of the secondary in system [solar masses] (Default=1, means M2 is a planet)
    @type Mass2:                             = float
    @param verbose:                          = Send prints to screen? 
    @type verbose:                           = Python boolean (True/False). Default = False
    
    Outputs:
    @param chi_squared_total:                = Total chi squared the model's position angle for all epochs
    @type chi_squared_total:                 = float
    @param ns:                               = Mean Motion [rad/yr]
    @type ns:                                = list of floats. Output values for each epoch.
    @param Ms:                               = Mean Anomaly [deg]
    @type Ms:                                = list of floats. Output values for each epoch.
    @param Es:                               = Eccentric Anomaly [deg]
    @type Es:                                = list of floats. Output values for each epoch.
    @param thetas:                           = True Anomaly [deg]
    @type thetas:                            = list of floats. Output values for each epoch.
    @param Sep_Dists:                        = Separation Distance in orbital plane [AU]
    @type Sep_Dists:                         = list of floats. Output values for each epoch.
    @param PA_deg_measured_models:           = measured Separation Angle of stars ["]
    @type PA_deg_measured_models:            = list of floats. Output values for each epoch.
    @param a1s:                              = semi-major axis of M1 [AU]
    @type a1s:                               = list of floats. Output values for each epoch.
    @param a2s:                              = semi-major axis of M2 [AU]
    @type a2s:                               = list of floats. Output values for each epoch.
    
    Note: Sep_Angle_arcsec_measured_REALs, PA_deg_measured_REALs, PA_mean_errors inputs MUST all have
          the same length.
    """
    
    # init output lists
    ns = []
    Ms = []
    Es = []
    thetas = []
    Sep_Dists = []
    PA_deg_measured_models = []
    a1s = []
    a2s = []

    # initial values for boolean equality of while statement 
    i = 0
    chi_squared_total = 0.0
    
    while (chiSquaredMax>=chi_squared_total)and(i<len(epochs)):    

        Sep_Angle_arcsec_measured_REAL = Sep_Angle_arcsec_measured_REALs[i]
        PA_deg_measured_REAL = PA_deg_measured_REALs[i]
        PA_mean_error = PA_mean_errors[i]
        t = epochs[i]  
        
        # call orbitCalculator to take random variables and calc orbital elements
        (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU, PA_deg_measured_model, a1, a2) = \
        orbitCalculator(t, Sys_Dist_PC, Sep_Angle_arcsec_measured_REAL, inclination_deg, longAN_deg, e, T, period, argPeri_deg, Mass1=Mass1, Mass2=Mass2, verbose=verbose)
        
        # store output orbital elements in lists for later plotting
        ns.append(n)
        Ms.append(M_deg)
        Es.append(E_latest_deg)
        thetas.append(TA_deg)
        Sep_Dists.append(Sep_Dist_AU)
        PA_deg_measured_models.append(PA_deg_measured_model)
        a1s.append(a1)
        a2s.append(a2)
            
        ## calc both PA and SA kai's and check they are less than chiSquaredMax
        PA_chi_squared = chiSquaredCalc(PA_deg_measured_REAL, PA_mean_error, PA_deg_measured_model)
        # Ad them to get the updated total
        chi_squared_total = chi_squared_total + PA_chi_squared # + SA_chi_squared
        
        # set boolean indicating the while made it to the last epoch and still passed
        lastEpochPass = False
        if (i==(len(epochs)-1))and(chiSquaredMax>=chi_squared_total):
            lastEpochPass = True
            
        # increment to next epoch to check these input for 
        i = i + 1   
        
        # while loop ends here! 
        
    if lastEpochPass:
        return (chi_squared_total, ns, Ms, Es, thetas, Sep_Dists, PA_deg_measured_models, a1s, a2s)
    else:
        return (False, None, None, None, None, None, None, None, None)

def normalizer(inList):
    """
    Just a simple function to normalize the elements of an input list and return the normalized version of the list.
    """

    max = np.max(inList)
    print "max of input lists = ",max
    outList = []
    count = 0
    for i in range(0,len(inList)):
        if count<10:
            print "dividing "+str(inList[i])+"/"+str(max)
            count +=1
        try:
            newVal = inList[i]/max
        except:
            print "Tried to divide "+str(inList[i])+"/"+str(max)+" and failed"
        outList.append(newVal)
    
    return outList

def orbitCalculator(t, Sys_Dist_PC, Sep_Angle_Measured_arcsec, inclination_deg, longAN_deg, e, T, period,\
                      argPeri_deg, Mass1=1, Mass2=1, verbose=False):
    """
    This function will calculate the output orbital parameters from the inputs for a SINGLE epoch.  Use the 
    multiEpochOrbCalc to do this for multiple epochs at the same time.  
    The t, Sys_Dist_PC and Sep_Angle_Measured_arcsec parameters are generally known/measured values,
    while the rest are typically random numbers.
    This function is designed to work as the calculator inside a Monte Carlo Simulator.
    All the equations are the result of Kepler's Laws and naturally occurring symmetries of binary star systems.
    
    See 'OrbitCalculator_2.pdf' for a full description of this function and its uses.
    
    Note: RP stands for Reference Plane, ie the plane that values are measured in from Earth.
          OP stands for Orbital Plane, ie in the plane that the stars are orbiting in.
    
    Inputs:
    @param t:                            = epoch of observation/image [julian date]
    @type t:                             = float
    @param Sys_Dist_PC:                  = measured system distance from Earth [PC]
    @type Sys_Dist_PC:                   = float
    @param Sep_Angle_Measured_arcsec:    = measured separation angle ["]
    @type Sep_Angle_Measured_arcsec:     = float
    @param inclination_deg:              = inclination [deg]
    @type inclination_deg:               = float
    @param longAN_deg:                   = Longitude of Ascending Node [deg]
    @type longAN_deg:                    = float
    @param e:                            = eccentricity of orbits [unitless]
    @type e:                             = float
    @param T:                            = Last Periapsis Epoch/time [julian date] 
    @type T:                             = float
    @param period:                       = period of orbits [yrs]
    @type period:                        = float
    @param argPeri_deg:                  = Argument of Periapsis in orbital plane [deg]
    @type argPeri_deg:                   = float
    @param Mass1:                        = Mass of primary in system [solar masses] (Default=1, means M2 is a planet)
    @type Mass1:                         = float
    @param Mass2:                        = Mass of the secondary in system [solar masses] (Default=1, means M2 is a planet)
    @type Mass2:                         = float
    @param verbose:                      = Send prints to screen? [True/False] (Default = False)
    @type verbose:                       = Python boolean (True/False). Default = False
    
    Outputs:
    @param n:                            = Mean Motion [rad/yr]
    @type n:                             = float
    @param M_deg:                        = Mean Anomaly [deg]
    @type M_deg:                         = float
    @param E_latest_deg:                 = Eccentric Anomaly [deg]
    @type E_latest_deg:                  = float
    @param TA_deg:                       = True Anomaly [deg]
    @type TA_deg:                        = float
    @param Sep_Dist_AU:                  = Separation Distance in orbital plane [AU]
    @type Sep_Dist_AU:                   = float
    @param PA_deg_measured_model:        = measured Position Angle in image [deg]
    @type PA_deg_measured_model:         = float
    @param a1:                           = semi-major axis of M1 [AU]
    @type a1:                            = float
    @param a2:                           = semi-major axis of M2 [AU]
    @type a2:                            = float
    
    """
    if verbose:
        print '\n'+'*'*50
        print 'Starting to calculate Orbital Parameters\n'
        print 'Input variable values: '+'\nCurrent epoch time [julian date] = '+str(t)+\
        '\nSystem Distance from Earth [PC] = '+str(Sys_Dist_PC)+\
        '\nInclination [deg] = '+ str(inclination_deg)+\
        '\nSeparation Angle ["] = '+str(Sep_Angle_Measured_arcsec)+\
        '\nLongitude of Ascending Node [deg] = '+str(longAN_deg)+\
        '\nEccentricity = '+str(e)+'\nLast Periapsis Epoch [julian date] = '+\
        str(T)+'\nPeriod [yrs] = '+str(period)+\
        '\nArgument of Periapsis in orbital plane [deg] = '+\
        str(argPeri_deg)+'\nMass of primary [solar masses] = '+str(Mass1)+\
        '\nMass of secondary [solar masses] = '+str(Mass2)+'\nverbose = ',str(verbose)

    # initial guess (E_last), will be updated in loop.  
    # Anything works, just takes longer if further from real value. => pi
    E_last = 3.14
    # stored initial value to be updated in loop
    E_latest = 1.0
    
    #print '*** period received by orbitCalculator = ',period #######$$$$$$$$$$$$$$$$

    ## calculate the Mean Motion
    n = (2*math.pi)/period
    if verbose:
        print 'Mean Motion [rad/yr]= '+str(n)

    ## calculate Mean Anomaly
    M = n*((t-T)/365.0)
    M_deg = math.degrees(M) # convert resulting M to degrees
    if verbose:
        print 'Mean Anomaly [deg]= ',M_deg

    ### Performing Newton's Method to get the Eccentric Anomaly value ###
    if verbose:
        print '-'*50

    # show input value to 
    if verbose:
        print "Inputs to Newton's Method are : \nM [rad]= "+str(M)+"\nE_last [rad]= "+\
        str(E_last)+"\nE_latest [rad] = "+str(E_latest)+"\ne = "+str(e)
        
        print "\nStarting to run Newton's while loop."

    count = 0 # a counter to stop inf loops in Newtons method below
    while (abs(E_last-E_latest) > (1.0e-10))and(count<100):
        if verbose:
            print 'current E [rad]= ', E_latest
        E_last = E_latest
        E_latest = E_last - ((M-E_last+e*math.sin(E_last))/(e*math.cos(E_last)-1.0))
        count = count+1

    E_latest_deg = math.degrees(E_latest) # convert resulting E to degrees
    if verbose:
        print "The resultant E value is [deg] = ", E_latest_deg
        print '-'*50
    ### Newton's loop finished! ###
    
    ## calculate True Anomaly from Eccentric Anomaly
    TA_rad  = math.acos((math.cos(E_latest)-e)/(1-e*math.cos(E_latest)))
    TA_deg = math.degrees(TA_rad)
    if verbose:
        print 'True Anomaly [deg] = ', TA_deg
        
    # calculate measured Separation Distance in reference plane
    Sep_Dist_measured_AU = Sep_Angle_Measured_arcsec*Sys_Dist_PC
    if verbose:
        print 'Separation Distance in the reference plane [AU] = ',Sep_Dist_measured_AU
        
    ## calculate measured Separation Distance in reference plane
    inclination_rad = math.radians(inclination_deg) # convert inclination to radians 
    
    # angle between ascending node and M1
    ang_AN_to_M1_deg = argPeri_deg + TA_deg
    
    ang_AN_to_M1_rad = math.radians(ang_AN_to_M1_deg) # convert angle to radians
    
    # Calculate the separation distance in the orbital plane
    Sep_Dist_AU_OP = Sep_Dist_measured_AU/math.sqrt(math.pow(math.sin(ang_AN_to_M1_rad)*math.cos(inclination_rad),2.0) + math.pow(math.cos(ang_AN_to_M1_rad),2.0))
    if verbose:
        print 'Separation Distance in orbital plane [AU] = ',Sep_Dist_AU_OP
    
    # find horizontal and vertical WRT Line of Nodes (X and Y) components of Sep_Dist corrected for inclination
    # Note: only the 'vertical (Y)' component will be effected by inclination
    Sep_Dist_AU_RP_y = Sep_Dist_AU_OP*math.sin(ang_AN_to_M1_rad)*math.cos(inclination_rad)
    Sep_Dist_AU_RP_x = Sep_Dist_AU_OP*math.cos(ang_AN_to_M1_rad)
    # use Pythagorous theorem to find the final Sep Dist in the Ref Plane
    # this is done only to varify math is right
    Sep_Dist_measured_AU2 = math.sqrt(math.pow(Sep_Dist_AU_RP_x, 2.0)+math.pow(Sep_Dist_AU_RP_y,2.0))
    #if verbose:
    #    print 'Separation Distance in the reference plane 2 [AU] = ',Sep_Dist_measured_AU2
        
    Sep_Dist_DIFF = math.sqrt(math.pow((Sep_Dist_measured_AU-Sep_Dist_measured_AU2),2.0))
    if Sep_Dist_DIFF > 0.0001:
        print 'The separation distances differ by more than 0.0001 !!!!!!'
        print 'Sep_Dist_measured_AU = ', Sep_Dist_measured_AU
        print 'Sep_Dist_measured_AU2 = ',Sep_Dist_measured_AU2
    
    ## calculate corrected angle between Line of Nodes and M1 (NOT AN to M1)
    ang_LN_to_M1_corr = math.degrees(math.atan(Sep_Dist_AU_RP_y/Sep_Dist_AU_RP_x))
    ## calculate measured Position Angle depending on which quadrant corrected M1 position is in
    if Sep_Dist_AU_RP_x>=0.0:
        PA_deg_measured_model = longAN_deg + 180.0 + ang_LN_to_M1_corr
        #print '** corrected M1 found in 2nd or 3rd quadrant' #$$$$$$$$
    elif Sep_Dist_AU_RP_x<0.0:
        PA_deg_measured_model = longAN_deg + ang_LN_to_M1_corr
        #print '** corrected M1 found in 1st or 4th quadrant' #$$$$$$$$$$$
    if verbose:
        print 'Position Angle in reference plane [deg] = ',PA_deg_measured_model
                                                     
    ## calculate the total semi-major axis (sum of both semi-majors for binary stars) in orbital plane
    a_total = (Sep_Dist_AU_OP*(1+e*math.cos(TA_rad)))/(1-e*e)
    if verbose:
        print 'Total semi-major axis (a1+a2) [AU] = ',a_total
    if (Mass1!=1) and (Mass2!=1):
        # this means these two parameters are non-default values and 
        # the system must then be a binary star system, thus calculate the
        # a1 and a2 values of the system.
        # NOTE: Mass2 MUST < Mass1, ie. primary is bigger.
        
        # first find the mass ratio of system
        X = Mass2/Mass1
        
        # find a1 from the a1=X*a2 and a_total=a1+a2
        a2 = a_total/(1.0+X)
        
        # now use a1=X*a2 to get a2
        a1 = a2*X

        if verbose: 
            print 'The system must be a binary star system with a1 = '+str(a1)+' and a2 = '+str(a2)
    else:
        # default Mass1 and Mass2, thus a planet system and a_total = a2 and a1=0.0
        # this is because we assume most planet's mass << star's mass, so a1<<a2
        a1 = 0.0
        a2 = a_total
        
    if verbose:
        print '\nFinished calculating Orbital Parameters'
        print '*'*50+'\n'
   
    return (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, PA_deg_measured_model, a1, a2)

def orbElementPlotter(plotFileTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
                  ns2, Ms2, Es2, thetas2, Sep_Dists2, a1s2, a2s2, chiSquareds, showPlots=False, summaryOnly=True, verbose=False):
    """
    This function is used to plot the resulting orbital parameters of the orbits that satisfied 
    the chi squared checks in a Monte Carlo simulator of the 'binary' code project.
    
    The plots will be be in multiple figures:
    INS: all the standard orbitalCalculator inputs in separate histogram plots of value vs probability
    OUTPUTS: all the standard orbitalCalculator outputs in separate histogram plots of value vs probability
             with one figure for each epoch.
    SUMMARY: a figure summarizing the orbital parameter histograms for the values quoted by Currie et all 2011.
    
    These figures will be saved to the directory '../figures/'.
    CAUTION! the saved figures will overwrite a previous version of the file if it has the same name!!
    
    INPUTS:
    @param plotFileTitle:                 = initial part of the figure and file title.  This will be 
                                            have strings appended to it to create individual
                                            input and output figure and file names.
    @type plotFileTitle:                  = string
    @param longAN_degs:                   = Longitude of Ascending Node [deg]
    @type longAN_degs:                    = list of floats
    @param es:                            = eccentricity of orbits [unitless]
    @type es:                             = list of floats
    @param Ts:                            = Last Periapsis Epoch/time [julian date] 
    @type Ts:                             = list of floats
    @param periods:                       = period of orbits [yrs]
    @type periods:                        = list of floats
    @param inclination_degs:              = inclination [deg]
    @type inclination_degs:               = list of floats
    @param argPeri_degs:                  = Argument of Periapsis in orbital plane [deg]
    @type argPeri_degs:                   = list of floats
    
    @param ns:                            = Mean Motions [rad/yr]
    @type ns:                             = list of floats. Output values for each epoch.
    @param Ms:                            = Mean Anomalys [deg]
    @type Ms:                             = list of floats. Output values for each epoch.
    @param Es:                            = Eccentric Anomalys [deg]
    @type Es:                             = list of floats. Output values for each epoch.
    @param thetas:                        = True Anomalys [deg]
    @type thetas:                         = list of floats. Output values for each epoch.
    @param Sep_Dists:                     = Separation Distances in orbital plane [AU]
    @type Sep_Dists:                      = list of floats. Output values for each epoch.
    @param a1s:                           = semi-major axis' of M1 [AU]
    @type a1s:                            = list of floats. Output values for each epoch.
    @param a2s:                           = semi-major axis' of M2 [AU]
    @type a2s:                            = list of floats. Output values for each epoch.
    @param chiSquareds:                   = total chi squared values for each orbit.
    @type chiSquareds:                    = list of floats
    
    @param showPlots:                     = show the resulting figures?
    @type showPlots:                      = Python boolean (True/False). Default = False.
    @param summaryOnly:                   = Only save the summary figure?
    @type summaryOnly:                    = Python boolean (True/False). Default = True.
    @param verbose:                       = show prints on screen while running?
    @type verbose:                        =  Python boolean (True/False). Default = False.
    
    """
    
    chiSquareMaxGooders =200
    
    inclination_degsAllGooders = []
    longAN_degsAllGooders = []
    argPeri_degsAllGooders = []
    esAllGooder = []
    periodsAllGooders = []
    TsAllGooders = []
    a1AllGooders = []
    a2AllGooders = []
    
    numSamples = np.shape(Es2[:][:])[0] 
    numEpochs = np.shape(Es2[:][:])[1]
    
    for orbit in range(0,numSamples):
        if chiSquareds[orbit]<=chiSquareMaxGooders:
            inclination_degsAllGooders.append(inclination_degs[orbit])
            esAllGooder.append(es[orbit])
            longAN_degsAllGooders.append(longAN_degs[orbit])
            periodsAllGooders.append(periods[orbit])
            TsAllGooders.append(Ts[orbit])
            argPeri_degsAllGooders.append(argPeri_degs[orbit])
                
    
    ########################## INS #####################
    # plot in's on first figure
    if verbose:
        print 'Starting to plot Inputs'
    plt.figure(1,figsize=(20,7) ,dpi=300)
    plt.subplot(231)
    (n,bins,patches) = plt.hist(longAN_degs,bins=50)
    plt.xlabel('longAN_degs')
    plt.ylabel('Probability')
    plt.subplot(232)
    (n,bins,patches) = plt.hist(es,bins=50)
    plt.xlabel('es')
    plt.title(plotFileTitle+' *Inputs*')
    plt.subplot(233)
    Tsmin = int(np.min(Ts))
    Tsmax = int(np.max(Ts))
    TsNEW = Ts
    for i in range(0,len(Ts)):
        TsNEW[i]=(TsNEW[i]-Tsmin)/365
    (n,bins,patches) = plt.hist(TsNEW,bins=25)
    plt.xlabel('Ts [years since '+str(Tsmin)+']')
    plt.subplot(234)
    (n,bins,patches) = plt.hist(periods,bins=50)
    plt.xlabel('periods')
    plt.ylabel('Probability')
    plt.subplot(235)
    (n,bins,patches) = plt.hist(inclination_degs,bins=50)
    plt.xlabel('inclination_degs')
    plt.subplot(236)
    (n,bins,patches) = plt.hist(argPeri_degs,bins=50)
    plt.xlabel('argPeri_degs')
    if summaryOnly is not True:
        plt.savefig('../figures/'+plotFileTitle+'INS.png', dpi=300, orientation='landscape')
    if verbose:
        print 'Finished plotting Inputs'
    
    ################# OUTS ####################
    if verbose:
        print 'numSamples = ',numSamples
        print 'numEpochs = ',numEpochs
    
    a1AllEpoch = []
    a2AllEpoch = []    
    
    for epoch in range(0,numEpochs):
        if verbose:
            print 'Starting to plot outputs for epoch ',str(epoch+1)
        # reshape the output matrices
        Sep_Dists2SingleEpoch = []
        thetas2SingleEpoch = []
        a1SingleEpoch = []
        a1AllEpoch = []
        a2SingleEpoch = []
        a2AllEpoch = []
        Es2SingleEpoch = []
        Ms2SingleEpoch = []
        ns2SingleEpoch = []
        for sample in range(0,numSamples):
            Sep_Dists2SingleEpoch.append(Sep_Dists2[sample][epoch]) 
            thetas2SingleEpoch.append(thetas2[sample][epoch]) 
            a1SingleEpoch.append(a1s2[sample][epoch]) 
            a1AllEpoch.append(a1s2[sample][epoch])
            a2SingleEpoch.append(a2s2[sample][epoch])
            a2AllEpoch.append(a2s2[sample][epoch]) 
            Es2SingleEpoch.append(Es2[sample][epoch]) 
            Ms2SingleEpoch.append(Ms2[sample][epoch]) 
            ns2SingleEpoch.append(ns2[sample][epoch])
            if chiSquareds[sample]<=chiSquareMaxGooders:
                 a1AllGooders.append(a1s2[sample][epoch])
                 a2AllGooders.append(a2s2[sample][epoch])
            
        # plot OUT's on first figure
        plt.figure(epoch+2, figsize=(20,9) ,dpi=300)
        plt.subplot(241)
        (n,bins,patches) = plt.hist(Sep_Dists2SingleEpoch,bins=50)
        plt.xlabel('Sep_Dists [AU]')
        plt.ylabel('Probability')
        plt.subplot(242)
        (n,bins,patches) = plt.hist(a2SingleEpoch,bins=50)
        plt.xlabel('a2s')
        plt.title(plotFileTitle+' *Outputs of simulator for epoch '+str(epoch+1)+'*')
        # only plot the a1s if they aren't all 0.0
        if np.mean(a1SingleEpoch)>0.0:
            plt.subplot(243)
            (n,bins,patches) = plt.hist(a1SingleEpoch,bins=50)
            plt.xlabel('a1s')
        plt.ylabel('Probability')
        plt.subplot(244)
        (n,bins,patches) = plt.hist(thetas2SingleEpoch,bins=50)
        plt.xlabel('thetas')
        plt.subplot(245)
        (n,bins,patches) = plt.hist(Es2SingleEpoch,bins=50)
        plt.xlabel('Es2')
        plt.subplot(246)
        (n,bins,patches) = plt.hist(Ms2SingleEpoch,bins=50)
        plt.xlabel('Ms2')
        plt.ylabel('Probability')
        plt.subplot(247)
        (n,bins,patches) = plt.hist(ns2SingleEpoch,bins=50)
        plt.xlabel('ns2')
        filename = '../figures/'+plotFileTitle+'OUTSepoch'+str(epoch+1)
        if summaryOnly is not True:
            plt.savefig('%s.png' % filename, dpi=300, orientation='landscape')
            
    # make a summary plot for inputs (i, e, longAN) and output (a2), to compare to Currie et al results
    if verbose:
        print 'Starting to make Summary Plot'
    plt.figure(numEpochs, figsize=(16,9) ,dpi=300)
    plt.subplot(231)
    (n,bins,patches) = plt.hist(longAN_degs,bins=50)
    plt.xlabel('long of AN [deg]')
    plt.ylabel('Probability')
    plt.subplot(232)
    (n,bins,patches) = plt.hist(es,bins=50)
    plt.xlabel('e ')
    plt.title(plotFileTitle+'  Currie et al. Comparison Summary')
    plt.subplot(233)
    (n,bins,patches) = plt.hist(a2AllEpoch,bins=50)#, range=(0.0,50.0)) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    plt.xlabel('a1 [AU]')
    plt.subplot(234)
    (n,bins,patches) = plt.hist(a2AllEpoch,bins=50)#, range=(0.0,50.0)) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    plt.xlabel('a2 [AU]')
    plt.ylabel('Probability')
    plt.subplot(235)
    (n,bins,patches) = plt.hist(inclination_degs,bins=50)
    plt.xlabel('inclination [deg]')
    filename = '../figures/'+plotFileTitle+'_Comparison_Summary'
    if summaryOnly is not True:
        plt.savefig('%s.png' % filename, dpi=300, orientation='landscape')
        

    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ##$$$$$$$$$$$ TOTAL SUMMARY PLOT SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    # make a summary plot for inputs (i, e, longAN) and output (a2), to compare to Currie et al results
    if verbose:
        print 'Starting to make Total Summary Plot'
    if len(inclination_degsAllGooders)>0:
        plt.figure(numEpochs, figsize=(25,22) ,dpi=300)
        plt.subplot(331)
        (n,bins,patches) = plt.hist(inclination_degsAllGooders,bins=50)
        plt.xlabel('inclination [deg]')
        plt.ylabel('Probability')
        plt.subplot(332)
        (n,bins,patches) = plt.hist(longAN_degsAllGooders,bins=50)
        plt.xlabel('long of AN [deg]')
        plt.title(plotFileTitle+' max chi-square='+str(chiSquareMaxGooders)+' TOTAL Summary')
        plt.subplot(333)
        (n,bins,patches) = plt.hist(argPeri_degsAllGooders,bins=50)
        plt.xlabel('argPeri_degs')
        plt.subplot(334)
        (n,bins,patches) = plt.hist(esAllGooder,bins=50)
        plt.xlabel('e ')
        plt.ylabel('Probability')
        plt.subplot(335)
        (n,bins,patches) = plt.hist(periodsAllGooders,bins=50)
        plt.xlabel('periods')
        plt.subplot(336)
        Tsmin = int(np.min(TsAllGooders))
        Tsmax = int(np.max(TsAllGooders))
        TsNEW = TsAllGooders
        for i in range(0,len(TsAllGooders)):
            TsNEW[i]=(TsNEW[i]-Tsmin)/365
        (n,bins,patches) = plt.hist(TsNEW,bins=25)
        plt.xlabel('Ts [years since '+str(Tsmin)+']')
        plt.subplot(337)
        (n,bins,patches) = plt.hist(a1AllGooders,bins=50)#, range=(0.0,50.0)) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        plt.xlabel('a1 [AU]')
        plt.ylabel('Probability')
        plt.subplot(338)
        (n,bins,patches) = plt.hist(a2AllGooders,bins=50)#, range=(0.0,50.0)) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        plt.xlabel('a2 [AU]')
        filename = '../figures/Total_Summaries/'+plotFileTitle+'_chiMax-'+str(chiSquareMaxGooders)+'_Total_Summary'
        plt.savefig('%s.png' % filename, dpi=300, orientation='landscape')
        
        if verbose:
            print 'Finished making Total Summary Plot'
            print str(len(inclination_degsAllGooders))+' orbits found with chi-square under '+str(chiSquareMaxGooders)
            print 'median inclination_degsAllGooders =  ',np.median(inclination_degsAllGooders)
            print 'median longAN_degsAllGooders = ', np.median(longAN_degsAllGooders)
            print 'median argPeri_degsAllGooders = ', np.median(argPeri_degsAllGooders)
            print 'median esAllGooder = ',np.median(esAllGooder)
            print 'median periodsAllGooders =  ',np.median(periodsAllGooders)
            print 'median TsAllGooders [jd] = '+ str((np.median(TsAllGooders)*365.5)+Tsmin)#+' OR [yrs since..] = '+str(np.median(TsAllGooders))
            print 'median a1AllGooders = ', np.median(a1AllGooders)
            print 'median a2AllGooders = ',np.median(a2AllGooders)
    else:
        if verbose:
            print 'No Gooders to plot in Total Summary !!'    
    # Show the plots on the screen if asked to    
    if showPlots:
        plt.show()
        plt.close()
    else:
        plt.close()
        
def orbitEllipsePlotter(longAN_deg, argPeri_deg, a, e, i, xLabel='', yLabel='', plotFilename='', xLim=False, yLim=False, PAs=[], SAs=[], sys_dist=1):
    """
    This function will plot the resulting orbit for the input parameters.
    
    @param longAN_deg:     = Longitude of the Acending Node in degrees
    @type longAN_deg:      = any type of number other than int
    @param argPeri_deg:    = Argument of Periapsis in degrees
    @type argPeri_deg:     = any type of number other than int
    @param a:              = Semi-Major axis in AU
    @type a:               = any type of number other than int
    @param e:              = Eccentricity
    @type e:               = any type of number other than int
    @param i:              = Inclination in degrees
    @type i:               = any type of number other than int
    @param xLim:          = range to limit the x axis values to
    @type xLim:           = tuple of two numbers of any type, (min,max)
                            default False indicates to use min and max of 
                            X values for points on the ellipse
                            +5% for white space padding
    @param yLim:          = range to limit the y axis values to
    @type yLim:           = tuple of two numbers of any type, (min,max)
                            default False indicates to use min and max of 
                            Y values for points on the ellipse 
                            +5% for white space padding
    @param PAs:           = Position Angles in [degrees]
    @type PAs:            = list of numbers, same length as SAs
    @param SAs:           = Separation Angles in ["]
    @type SAs:            = list of numbers, same lenght as PAs
    @param sys_dist:      = Distance to System in [pc]
    @type sys_dist:       = any type of number other than int
    """
    
    # First convert the given angles to radians
    i_rad = np.radians(i)
    longAN_rad = np.radians(longAN_deg)
    argPeri_rad = np.radians(argPeri_deg)
    
    # must correct argPeri for the inclination
    argPeri_corr = np.arctan((np.sin(argPeri_rad)*np.cos(i_rad))/np.cos(argPeri_rad))
    argPeri_corr_deg = np.degrees(argPeri_corr)
    print 'argPeri_deg = '+str(argPeri_deg)
    print 'argPeri_corr_deg = '+str(argPeri_corr_deg)
    print 'ang orig = '+str(longAN_deg+argPeri_deg+90.0)
    print 'ang corr = '+str(longAN_deg+argPeri_corr_deg+90.0)
    
    # angle for ellipse plotter
    ang = longAN_deg+argPeri_deg+90.0
    ang_corr = longAN_deg+argPeri_corr_deg+90.0

     # semi-major in orbital plane
    a_op = a
    
    # semi-minor in orbital plane
    b_op = a*(np.sqrt(1-(e*e)))
    
    # semi-major in reference plane
    a_rp_y = a_op*np.sin(argPeri_rad)*np.cos(i_rad)
    a_rp_x = a_op*np.cos(argPeri_rad)
    a_rp = np.sqrt(np.power(a_rp_x, 2)+np.power(a_rp_y, 2))
    
    print 'a_rp_x = '+str(a_rp_x)
    print 'a_rp_y = '+str(a_rp_y)
    print 'a_rp = '+str(a_rp)
    
    # semi-minor in reference plane 
    # Note: semi-minor and semi-major are at 90deg to each other.
    b_rp_y = b_op*np.sin((pi/2.0)-argPeri_rad)*np.cos(i_rad)
    b_rp_x = b_op*np.cos((pi/2.0)-argPeri_rad)
    b_rp = np.sqrt(np.power(b_rp_x, 2)+np.power(b_rp_y, 2))
    
    print 'b_op = '+str(b_op)
    print 'b_rp_x = '+str(b_rp_x)
    print 'b_rp_y = '+str(b_rp_y)
    print 'b_rp = '+str(b_rp)

    # Calculate distance of center-foci assuming center is at (0,0)
    c_foci = a_rp-((a_rp*(1.0-e*e))/(1.0-e))
    # Calculate loction of foci where star would lie
    xStar = c_foci*np.cos(np.radians(ang_corr))
    yStar = c_foci*np.sin(np.radians(ang_corr))
    
    print 'c_foci = '+str(c_foci)
    print 'yStar = '+str(yStar)
    print 'xStar = '+str(xStar)
    
    ## semi-major line   
    x_a_rp_corr = a_rp*np.cos(np.radians(ang_corr))
    print 'x_a_rp_corr = '+str(x_a_rp_corr)
    y_a_rp_corr = a_rp*np.sin(np.radians(ang_corr))
    print 'y_a_rp_corr = '+str(y_a_rp_corr)
    # take care of the two cases for >< 90deg inclination
    if x_a_rp_corr>0.0:
        if i>90.0:
            #print 'in >90 , neg x case'
            k = 1.0
        else:
            #print 'in <90 , neg x case'
            k = (-1.0)
    if x_a_rp_corr<0.0:
        if i>90.0:
            #print 'in >90 , pos x case'
            k = (-1.0)
        else:
            #print 'in <90 , pos x case'
            k = 1.0
    
    if y_a_rp_corr>0.0:
        if i>90.0:
            #print 'in >90 , neg y case'
            j = 1.0
        else:
            #print 'in <90 , neg y case'
            j = (-1.0)
    if y_a_rp_corr<0.0:
        if i>90.0:
            #print 'in >90 , pos y case'
            j = (-1.0)
        else:
            #print 'in <90 , pos y case'
            j = 1.0
    X_a_rp_corr = [k*abs(x_a_rp_corr)-xStar,(-k)*abs(x_a_rp_corr)-xStar]
    Y_a_rp_corr = [j*abs(y_a_rp_corr)-yStar,(-j)*abs(y_a_rp_corr)-yStar]
    
    print 'x_a_corrs = '+repr(X_a_rp_corr)
    print 'y_a_corrs = '+repr(Y_a_rp_corr)
    
    ## semi-minor line
    x_b_rp_corr = b_rp*np.cos(np.radians(ang_corr+90.0))
    print 'x_b_rp_corr = '+str(x_b_rp_corr)
    y_b_rp_corr = b_rp*np.sin(np.radians(ang_corr+90.0))
    print 'y_b_rp_corr = '+str(y_b_rp_corr)
    # take care of the two cases for >< 90deg inclination
    if x_b_rp_corr>0.0:
        if i>90.0:
            #print 'in >90 , neg x case'
            k = 1.0
        else:
            #print 'in <90 , neg x case'
            k = (-1.0)
    if x_b_rp_corr<0.0:
        if i>90.0:
            #print 'in >90 , pos x case'
            k = (-1.0)
        else:
            #print 'in <90 , pos x case'
            k = 1.0
    
    if y_b_rp_corr>0.0:
        if i>90.0:
            #print 'in >90 , neg y case'
            j = 1.0
        else:
            #print 'in <90 , neg y case'
            j = (-1.0)
    if y_b_rp_corr<0.0:
        if i>90.0:
            #print 'in >90 , pos y case'
            j = (-1.0)
        else:
            #print 'in <90 , pos y case'
            j = 1.0
    X_b_rp_corr = [(-k)*abs(x_b_rp_corr)-xStar,k*abs(x_b_rp_corr)-xStar]
    Y_b_rp_corr = [(-j)*abs(y_b_rp_corr)-yStar,j*abs(y_b_rp_corr)-yStar]
    
    print 'x_b_corrs = '+repr(X_b_rp_corr)
    print 'y_b_corrs = '+repr(Y_b_rp_corr)
        
    # create a star polygon
    starPolygon = star(0.06, 0, 0, color='black', N=5, thin = 0.5)
    
    # determin x and y locations of the observed PA and SA's for companion star/planet
    # then make a star polygon for each, same as for M1 but much smaller
    m2starPolygons = []
    for i in range(0,len(SAs)):
        xStar2 = -SAs[i]*sys_dist*math.sin(math.radians(PAs[i]))
        yStar2 = SAs[i]*sys_dist*math.cos(math.radians(PAs[i]))
        m2starPolygons.append(star(0.03, xStar2, yStar2, color='red', N=5, thin = 0.5))
        
    X,Y = ellipse(a_rp, b_rp, ang_corr, -xStar, -yStar, Nb=500)
    #X,Y = ellipse(a_rp, b_rp, ang, 0, 0, Nb=500)
    
    # set the xLim and yLim if their values are False
    # and pad max and min values found by 10%
    if not xLim:
        min = np.min(X)
        max = np.max(X)
        Range = abs(max)-abs(min)
        xLim = (min-abs(Range*0.05),max+abs(Range*0.05))
    else:
        if not (type(xLim)==tuple):
            print 'PROBLEM: xLim is not of type tuple'
    if not yLim:
        min = np.min(Y)
        max = np.max(Y)
        Range = abs(max)-abs(min)
        yLim = (min-abs(Range*0.05),max+abs(Range*0.05))
    else:
        if not (type(yLim)==tuple):
            print 'PROBLEM: yLim is not of type tuple'
    
    fig = plt.figure(1,figsize=(10,10))
    main = fig.add_subplot(111)
    main.axes.set_xlim(xLim)
    main.axes.set_ylim(yLim)
    main.axes.set_xlabel(xLabel)
    main.axes.set_ylabel(yLabel)
    #draw ellipse
    main.plot(X,Y) #$$$$ add title, labels, axes
    #draw semi-major
    main.plot(X_a_rp_corr,Y_a_rp_corr,'g-')
    #draw semi-minor 
    main.plot(X_b_rp_corr,Y_b_rp_corr,'b-')
    #draw circle for center of ellipse
    main.plot(-xStar,-yStar,'rx')
    main.add_patch(starPolygon)
    for star2 in m2starPolygons:
        main.add_patch(star2)
    plt.show()
    plt.close()
    
    
def orbitEllipsePlotter2(longAN_deg, argPeri_deg, a, e, inc, period, xLabel='delta RA [mas]', yLabel='delta Dec [mas]', plotFilename='', \
                                                        xLim=False, yLim=False, PAs=[], SAs=[], sys_dist=1):
    """
    This function will plot the resulting orbit for the input parameters.
    NOTE: If a plotFilename is provided, then the resulting figure will be saved.
    
    @param longAN_deg:     = Longitude of the Acending Node in degrees
    @type longAN_deg:      = any type of number other than int
    @param argPeri_deg:    = Argument of Periapsis in degrees
    @type argPeri_deg:     = any type of number other than int
    @param a:              = Semi-Major axis in AU
    @type a:               = any type of number other than int
    @param e:              = Eccentricity
    @type e:               = any type of number other than int
    @param i:              = Inclination in degrees
    @type i:               = any type of number other than int
    @param period:         = period of orbits [yrs]
    @type period:          = float
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
    if plotFilename!='':
        if plotFilename[-4:]!='.png':
            plotFilename = plotFilename+'.png'     
    # create a plot title from the filename
    plotFileTitle = os.path.basename(plotFilename).split('.')[0]
    
    # First convert the given angles to radians
    i_rad = np.radians(inc)
    longAN_rad = np.radians(longAN_deg)
    argPeri_rad = np.radians(argPeri_deg)
    
    # must correct argPeri for the inclination
    argPeri_corr = np.arctan((np.sin(argPeri_rad)*np.cos(i_rad))/np.cos(argPeri_rad))
    argPeri_corr_deg = np.degrees(argPeri_corr)
#    print 'argPeri_deg = '+str(argPeri_deg)
#    print 'argPeri_corr_deg = '+str(argPeri_corr_deg)
#    print 'ang orig = '+str(longAN_deg+argPeri_deg+90.0)
#    print 'ang corr = '+str(longAN_deg+argPeri_corr_deg+90.0)
#    print 'inclination = '+str(inc)

    # angle for ellipse plotter
    ang = longAN_deg+argPeri_deg+90.0
    ang_corr = longAN_deg+argPeri_corr_deg+90.0

     # semi-major in orbital plane
    a_op = a
    
    # semi-minor in orbital plane
    b_op = a*(np.sqrt(1-(e*e)))
    
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
    b_rp_y = b_op*np.sin((pi/2.0)-argPeri_rad)*np.cos(i_rad)
    b_rp_x = b_op*np.cos((pi/2.0)-argPeri_rad)
    b_rp = np.sqrt(np.power(b_rp_x, 2)+np.power(b_rp_y, 2))
    
#    print 'b_op = '+str(b_op)
#    print 'b_rp_x = '+str(b_rp_x)
#    print 'b_rp_y = '+str(b_rp_y)
#    print 'b_rp = '+str(b_rp)

    # Calculate distance of center-foci assuming center is at (0,0)
    c_foci = a_rp-((a_rp*(1.0-e*e))/(1.0-e))
    # Calculate loction of foci where star would lie
    xStar = c_foci*np.cos(np.radians(ang_corr))*(1000.0/sys_dist)
    yStar = c_foci*np.sin(np.radians(ang_corr))*(1000.0/sys_dist)
    
#    print 'c_foci = '+str(c_foci)
#    print 'yStar = '+str(yStar)
#    print 'xStar = '+str(xStar)
    
        
    # create a primary star polygon
    starPolygon = star(0.9*a, 0, 0, color='black', N=5, thin = 0.5)
    
    # determin x and y locations of the observed PA and SA's for companion star/planet
    # then make a star polygon for each, same as for M1 but much smaller
    m2starPolygons = []
    for i in range(0,len(SAs)):
        xStar2 = -SAs[i]*math.sin(math.radians(PAs[i]))*1000.0#*sys_dist   # This will be in [mas] instead of [AU]
        yStar2 = SAs[i]*math.cos(math.radians(PAs[i]))*1000.0#*sys_dist   # This will be in [mas] instead of [AU]
        m2starPolygons.append(star(0.4*a, xStar2, yStar2, color='red', N=5, thin = 0.5))
        
    ## calculate the locations of companion for 'numOrbs' locations throughout the orbit to make an orbit ellipse    
    ellipseXs = []
    ellipseYs = []
    orbitTAs = []
    orbitSAs = []
    orbitPAs = []
    sep_dists = []
    numOrbs = 5000.0
    periodIncrement = (period*365.25)/numOrbs
    t = 1.0 #3.0*((period*365.25)/4.0)
    for orb in range(0,int(numOrbs)):
        T = 0.0
        t = t + periodIncrement
        (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SA, PA, a1, a2) =\
            orbitCalculator(t, sys_dist, inc, longAN_deg, e, T, period, argPeri_deg, a,\
                                                        Mass1=1, Mass2=1, verbose=False)
        #(PA, SA) = PASAcalculator(period, t, T, e, inc, longAN_deg, argPeri_deg, sys_dist, a, a1=0, verbose=False)
        orbitTAs.append(TA_deg)
        orbitPAs.append(PA)
        orbitSAs.append(SA)
        ellipseX = -SA*math.sin(math.radians(PA))*1000.0#*sys_dist  # This will be in [mas] instead of [AU]
        ellipseY = SA*math.cos(math.radians(PA))*1000.0#*sys_dist   # This will be in [mas] instead of [AU]
        ellipseXs.append(ellipseX)
        ellipseYs.append(ellipseY)
        sep_dist = math.sqrt(math.pow(ellipseX,2.0)+math.pow(ellipseY,2.0))
        sep_dists.append(sep_dist)
        
    ## Get the locations of 500 points on an ellipse representing the orbit # this isn't working right... thus orbit method above is used now.
    #X,Y = ellipse(a_rp, b_rp, ang_corr, -xStar, -yStar, Nb=500)
    #X,Y = ellipse(a_rp, b_rp, ang, 0, 0, Nb=500)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, a1, a2) =\
            orbitCalculator(0.0, sys_dist, inc, longAN_deg, e, 0.0, period, argPeri_deg, a,\
                                                        Mass1=1, Mass2=1, verbose=False)
    Xstart = -SAstart*math.sin(math.radians(PAstart))*1000.0#*sys_dist  # This will be in [mas] instead of [AU]
    Ystart = SAstart*math.cos(math.radians(PAstart))*1000.0#*sys_dist  # This will be in [mas] instead of [AU]
    startStar = star(0.6*a, Xstart, Ystart, color='green', N=5, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, a1, a2) =\
            orbitCalculator(((period/4.0)*365.25), sys_dist, inc, longAN_deg, e, 0.0, period, argPeri_deg, a,\
                                                        Mass1=1, Mass2=1, verbose=False)
    XoneQuarter = -SAstart*math.sin(math.radians(PAstart))*1000.0#*sys_dist  # This will be in [mas] instead of [AU]
    YoneQuarter = SAstart*math.cos(math.radians(PAstart))*1000.0#*sys_dist  # This will be in [mas] instead of [AU]
    oneQuarterStar = star(0.6*a, XoneQuarter, YoneQuarter, color='blue', N=5, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAhalf, PAhalf, a1, a2) =\
            orbitCalculator(((period/2.0)*365.25), sys_dist, inc, longAN_deg, e, 0.0, period, argPeri_deg, a,\
                                                        Mass1=1, Mass2=1, verbose=False)
    Xhalf = -SAhalf*math.sin(math.radians(PAhalf))*1000.0#*sys_dist  # This will be in [mas] instead of [AU]
    Yhalf = SAhalf*math.cos(math.radians(PAhalf))*1000.0#*sys_dist  # This will be in [mas] instead of [AU]
    halfStar = star(0.6*a, Xhalf, Yhalf, color='yellow', N=5, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, a1, a2) =\
            orbitCalculator((3.0*(period/4.0)*365.25), sys_dist, inc, longAN_deg, e, 0.0, period, argPeri_deg, a,\
                                                        Mass1=1, Mass2=1, verbose=False)
    XthreeQuarter = -SAstart*math.sin(math.radians(PAstart))*1000.0#*sys_dist  # This will be in [mas] instead of [AU]
    YthreeQuarter = SAstart*math.cos(math.radians(PAstart))*1000.0#*sys_dist  # This will be in [mas] instead of [AU]
    threeQuarterStar = star(0.6*a, XthreeQuarter, YthreeQuarter, color='orange', N=5, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, a1, a2) =\
            orbitCalculator(((period)*365.25), sys_dist, inc, longAN_deg, e, 0.0, period, argPeri_deg, a,\
                                                        Mass1=1, Mass2=1, verbose=False)
    Xend = -SAstart*math.sin(math.radians(PAstart))*1000.0#*sys_dist  # This will be in [mas] instead of [AU]
    Yend = SAstart*math.cos(math.radians(PAstart))*1000.0#*sys_dist  # This will be in [mas] instead of [AU]
    endStar = star(0.6*a, Xend, Yend, color='purple', N=5, thin = 0.5)
    
    # set the xLim and yLim if their values are False
    # and pad max and min values found by 10%
    if not xLim:
        min = np.min(ellipseXs)
        max = np.max(ellipseXs)
        Range = abs(max)+abs(min)
        xLim = (min-abs(Range*0.5),max+abs(Range*0.5))
    else:
        if not (type(xLim)==tuple):
            print 'PROBLEM: xLim is not of type tuple'
    if not yLim:
        min = np.min(ellipseYs)
        max = np.max(ellipseYs)
        Range = abs(max)+abs(min)
        yLim = (min-abs(Range*0.5),max+abs(Range*0.5))
    else:
        if not (type(yLim)==tuple):
            print 'PROBLEM: yLim is not of type tuple'
            
    ## Create figure, axes and start plotting/drawing everything
    fig = plt.figure(1,figsize=(10,10))
    main = fig.add_subplot(111)
    main.axes.set_xlim(xLim)
    main.axes.set_ylim(yLim)
    main.set_xlabel(xLabel)
    main.set_ylabel(yLabel)
    main.set_title(plotFileTitle)
    # Draw orbit
    main.plot(ellipseXs,ellipseYs) #$$$$ add title, labels, axes
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
    
    # Draw red star patches for each of the companion's locations from the data
    for star2 in m2starPolygons:
        main.add_patch(star2)
    if plotFilename!='':
        plt.savefig(plotFilename, dpi=300, orientation='landscape')
    plt.show()
    plt.close()    
    
def PASApredictor(period, t, T, e, inclination_deg, longAN_deg, argPeri_deg, Sys_Dist_PC, a1, a2, verbose=True):
    
    if verbose:
        print '\n'+'*'*50
        print 'Input variable values: '+\
        '\nPeriod [yrs] = '+str(period)+\
        '\nCurrent epoch time [julian date] = '+str(t)+\
        '\nLast Periapsis Epoch [julian date] = '+str(T)+\
        '\nEccentricity = '+str(e)+\
        '\nInclination [deg] = '+ str(inclination_deg)+\
        '\nLongitude of Ascending Node [deg] = '+str(longAN_deg)+\
        '\nArgument of Periapsis in orbital plane [deg] = '+str(argPeri_deg)+\
        '\nSystem Distance from Earth [PC] = '+str(Sys_Dist_PC)+\
        '\nSemi-major of primary [AU] = '+str(a1)+\
        '\nSemi-major of secondary [AU] = '+str(a2)
        print '\n**Starting to calculate PA and SA from Orbital Parameters**\n'
    #----------------------------------------------------------------------
    
    # initial guess (E_last), will be updated in loop.  
    # Anything works, just takes longer if further from real value. => pi
    E_last = 3.14
    # stored initial value to be updated in loop
    E_latest = 1.0

    ## calculate the Mean Motion
    n = (2*math.pi)/period
    if verbose:
        print 'Mean Motion [rad/yr]= '+str(n)

    ## calculate Mean Anomaly
    M = n*((t-T)/365.0)
    M_deg = math.degrees(M) # convert resulting M to degrees
    if verbose:
        print 'Mean Anomaly [deg]= ',M_deg

    ### Performing Newton's Method to get the Eccentric Anomaly value ###
    if verbose:
        print '-'*50

    ## show input value to 
    #if verbose:
    #    print "Inputs to Newton's Method are : \nM [rad]= "+str(M)+"\nE_last [rad]= "+\
    #    str(E_last)+"\nE_latest [rad] = "+str(E_latest)+"\ne = "+str(e)
        
        print "Starting to run Newton's while loop to find the Eccentric Anomaly."

    count = 0 # a counter to stop inf loops in Newtons method below
    while (abs(E_last-E_latest) > (1.0e-10))and(count<100):
        #if verbose:
        #    print 'current E [rad]= ', E_latest
        E_last = E_latest
        E_latest = E_last - ((M-E_last+e*math.sin(E_last))/(e*math.cos(E_last)-1.0))
        count = count+1

    E_latest_deg = math.degrees(E_latest) # convert resulting E to degrees
    if verbose:
        print "Finished and the resultant E value is [deg] = ", E_latest_deg
        print '-'*50
    ### Newton's loop finished! ###
    
    ## calculate True Anomaly from Eccentric Anomaly
    TA_rad  = math.acos((math.cos(E_latest)-e)/(1-e*math.cos(E_latest)))
    TA_deg = math.degrees(TA_rad)
    if verbose:
        print 'True Anomaly [deg] = ', TA_deg
    
    ## calculate the Separation Distance in orbital plane
    Sep_Dist_AU_OP = ((a1+a2)*(1-e*e))/(1+e*math.cos(TA_rad))
    if verbose:
        print 'Separation Distance in orbital plane [AU] = ',Sep_Dist_AU_OP
        
    ## calculate measured Separation Distance in reference plane
    inclination_rad = math.radians(inclination_deg) # convert inclination to radians 
    
    # angle between ascending node and M1
    ang_AN_to_M1_deg = argPeri_deg + TA_deg
    
    ang_AN_to_M1_rad = math.radians(ang_AN_to_M1_deg) # convert angle to radians
        
    # Calculate the separation distance in the orbital plane
    Sep_Dist_AU_RP = Sep_Dist_AU_OP*math.sqrt(math.pow(math.sin(ang_AN_to_M1_rad)*math.cos(inclination_rad),2.0) + math.pow(math.cos(ang_AN_to_M1_rad),2.0))
    if verbose:
        print 'Separation Distance in reference plane [AU] = ',Sep_Dist_AU_RP
        
    ## calculate measured Separation Angle (model)
    Sep_Angle_arcsec_measured_model = Sep_Dist_AU_RP/Sys_Dist_PC # PC->AU and rad->arcsec cancel
    print '\nSeparation Angle measured (model) [arcsec] = ',Sep_Angle_arcsec_measured_model    
        
    # find horizontal and vertical WRT Line of Nodes (X and Y) components of Sep_Dist corrected for inclination
    # Note: only the 'vertical (Y)' component will be effected by inclination
    Sep_Dist_AU_RP_y = Sep_Dist_AU_OP*math.sin(ang_AN_to_M1_rad)*math.cos(inclination_rad)
    Sep_Dist_AU_RP_x = Sep_Dist_AU_OP*math.cos(ang_AN_to_M1_rad)
    
    ## calculate corrected angle between Line of Nodes and M1 (NOT AN to M1)
    ang_LN_to_M1_corr = math.degrees(math.atan(Sep_Dist_AU_RP_y/Sep_Dist_AU_RP_x))
    ## calculate measured Position Angle depending on which quadrant corrected M1 position is in
    if Sep_Dist_AU_RP_x>=0.0:
        PA_deg_measured_model = longAN_deg + 180.0 + ang_LN_to_M1_corr
        #print '** corrected M1 found in 2nd or 3rd quadrant' #$$$$$$$$
    elif Sep_Dist_AU_RP_x<0.0:
        PA_deg_measured_model = longAN_deg + ang_LN_to_M1_corr
        #print '** corrected M1 found in 1st or 4th quadrant' #$$$$$$$$$$$
    print 'Position Angle in reference plane [deg] = ',PA_deg_measured_model
    print '\n'+'*'*50    
    return(PA_deg_measured_model, Sep_Angle_arcsec_measured_model)
    
def TAcalculator(t,e, T, period, verbose=False):
    
    ## calculate the Mean Motion
    n = (2*pi)/period
    if verbose:
        print 'Mean Motion [rad/yr]= '+str(n)
    
    ## calculate Mean Anomaly
    M = n*((t-T)/365.25)
    M_deg = math.degrees(M) # convert resulting M to degrees
    if verbose:
        print 'Mean Anomaly [deg]= ',M_deg

    ### Performing Newton's Method to get the Eccentric Anomaly value ###
    if verbose:
        print '-'*50
        
    # initial guess (E_last), will be updated in loop.  
    # Anything works, just takes longer if further from real value. => pi
    E_last = 2*pi
    # stored initial value to be updated in loop
    # this value is always very close to the true value and will minimize the number of loops
    E_latest = M+e*math.sin(M) 
    # show input value to 
    if verbose:
        print "Inputs to Newton's Method are : \nM [rad]= "+str(M)+"\nE_last [rad]= "+\
        str(E_last)+"\nE_latest [rad] = "+str(E_latest)+"\ne = "+str(e)
        
        print "\nStarting to run Newton's while loop."
    
    count = 0 # a counter to stop inf loops in Newton's method below
    while (abs(E_last-E_latest) > (1.0e-10))and(count<100):
        if verbose:
            print 'current E [rad]= ', E_latest
        E_last = E_latest
        #E_latest = E_last - ((M-E_last+e*math.sin(E_last))/(e*math.cos(E_last)-1.0))\
        E_latest = E_last - ((E_last-M-e*math.sin(E_last))/(-e*math.cos(E_last)+1.0))
        count = count+1

    E_latest_deg = math.degrees(E_latest) # convert resulting E to degrees
    if verbose:
        print "The resultant E value is [deg] = ", E_latest_deg
    # check if the resultant value solves the original equation.
    Mnewton = math.degrees(E_latest-e*math.sin(E_latest))
    if abs(M_deg-Mnewton)>(1.0e-5):
        if verbose:
            print "PROBLEM: This resultant E does not satisfy the original equation, Newton's method FAILED !!!"
            print 'M from this E Equals = '+str(Mnewton)
            print 'M original = '+str(M_deg)
            print 'E initial = '+str(math.degrees(M+e*math.sin(M) ))
            print 'e = '+str(e)
    else:
        if verbose:
            print "This resultant E solves the original equation, Newton's Method worked :-)"
            print '-'*50
    ### Newton's loop finished! ###
    
    ## calculate True Anomaly from Eccentric Anomaly
    try:
        TA_rad  = math.acos((math.cos(E_latest)-e)/(1.0-e*math.cos(E_latest))) 
        
    except:
        print "\nTA_rad calc failed!"
        print "E_latest = ",E_latest
        print "e = ",e
    
    ## Calculate TA in another way    
    x = ((1.0-e**2.0)**(1.0/2.0))*math.cos(E_latest/2.0)
    y = ((1.0+e**2.0)**(1.0/2.0))*math.sin(E_latest/2.0)
    TA_rad2 = 2.0*math.atan2(y, x)
    #print 'TA_2 = '+str(math.degrees(TA_rad2))

    if E_latest>pi:
        # this is to take care of the silly fact that after E=180deg,
        # the equation above produces TAs that go down from 180, rather than up.
        TA_rad = 2.0*pi - TA_rad
    
    return (n, M_deg, E_latest_deg,TA_rad2)

        
def TriangProbDensFunc(min, max, mode, x):
    '''
    This function will calculate the normalized probability density (NPD) of a triangular 
    random number produced using the built in 'rand' package of Python.
    
    '''
    worked=False
    if x<min:
        f = 0
        worked=True
    elif (min<=x)and(x<mode):
        f = (x-min)/(mode-min)
        worked=True
    elif x==mode:
        f = 1
        worked=True
    elif (mode<x)and(x<=max):
        f = (max-x)/(max-mode)
        worked=True
    elif max<x:
        f = 0
        worked=True
        
    if not worked:
        print 'min='+str(min)+', max='+str(max)+' mode='+str(mode)+' x='+str(x)
        
    normal = f
        
    return normal
    
def dataWriterLowRAM(plotFileTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs, Sep_Dists2, a1s2, a2s2, chiSquareds):
    
    """
    *****************************************************************************************************
    *** THIS IS THE LOW RAM VERSION of dataWriter for laptops or computers with less than 16GB of RAM ***
    *****************************************************************************************************
    
    This function is to write the resulting inputs and outputs of the mcmcOrbSimulatorUniform2.
    The inputs and outputs for each epoch will all be written to separate files.
    
    The function dataReader can later be used to reload these files back into memory for later
    data analysis or plotting.
    
    CAUTION! the saved data will overwrite a previous version of the file if it has the same name!!
    
    INPUTS:
    @param plotFileTitle:                 = initial part of the file title.  This will be 
                                            have strings appended to it to create individual
                                            input and output file names.
    @type plotFileTitle:                  = string
    @param longAN_degs:                   = Longitude of Ascending Node [deg]
    @type longAN_degs:                    = list of floats
    @param es:                            = eccentricity of orbits [unitless]
    @type es:                             = list of floats
    @param Ts:                            = Last Periapsis Epoch/time [julian date] 
    @type Ts:                             = list of floats
    @param periods:                       = period of orbits [yrs]
    @type periods:                        = list of floats
    @param inclination_degs:              = inclination [deg]
    @type inclination_degs:               = list of floats
    @param argPeri_degs:                  = Argument of Periapsis in orbital plane [deg]
    @type argPeri_degs:                   = list of floats
    
    @param Sep_Dists2:                     = Separation Distance in orbital plane [AU]
    @type Sep_Dists2:                      = list of floats. Output values for each epoch.
    @param a1s:                           = semi-major axis of M1 [AU]
    @type a1s:                            = list of floats. Output values for each epoch.
    @param a2s:                           = semi-major axis of M2 [AU]
    @type a2s:                            = list of floats. Output values for each epoch.
    @param chiSquareds:                   = total chi squared values for each orbit.
    @type chiSquareds:                    = list of floats
    
    """
    # find how many samples and epochs are to be written to files
    numSamples = np.shape(a1s2[:][:])[0]
    numEpochs = np.shape(a1s2[:][:])[1]
    
    ############## INS ###################
    filename = '../data/'+plotFileTitle+'_INS'
    f = open('%s.txt' % filename, 'w')
    f.write(plotFileTitle+'_INS'+'\n')
    f.write('longAN [deg]    es [N/A]    Ts [julian date]    periods [yrs]   inclination [deg]   argPeri [deg]   chiSquareds\n')
    for sample in range(0,numSamples):
        line = str(longAN_degs[sample])
        line = line +'   '+ str(es[sample])
        line = line +'   '+ str(Ts[sample])
        line = line +'   '+ str(periods[sample])
        line = line +'   '+ str(inclination_degs[sample])
        line = line +'   '+ str(argPeri_degs[sample])
        line = line +'   '+ str(chiSquareds[sample])+'\n'
        f.write(line)
    f.close()
    #######################################
    
    ################# OUTS ####################
    for epoch in range(0,numEpochs):
        filename = '../data/'+plotFileTitle+'_LowRAM_OUTSepoch'+str(epoch+1)
        f = open('%s.txt' % filename, 'w')
        f.write(plotFileTitle+'_LowRAM_OUTSepoch'+str(epoch+1)+'\n')
        f.write('a1s [AU]   a2s [AU]\n')
        for sample in range(0,numSamples):
            #line = str(Sep_Dists2[sample][epoch])
            #line = line +' '+ str(thetas2[sample][epoch])
            #line = line +' '+ str(Es2[sample][epoch])
            #line = line +' '+ str(Ms2[sample][epoch])
            #line = line +' '+ str(ns2[sample][epoch])
            line = str(a1s2[sample][epoch])
            line = line +' '+ str(a2s2[sample][epoch])+'\n'
            f.write(line)
        f.close()
      
def dataReaderLowRAM(INPUTSfilename, OUTPUTSfilenames, verbose=False):
    """ 
    This function is to read data files written to disk by the dataWriter function back
    into memory for plotting or statistical analysis. 
    
    @param INPUTSfilename:     = file name of the file containing the 'inputs' file 
    @type INPUTSfilename:      = single string for the filename or a list of filename strings
    @param OUTPUTSfilenames:   = list containing the file names of the 'outputs' files
    @type OUTPUTSfilenames:    = list of filename strings. Order of list should be epoch1, epoch2,...
        
    NOTE:
        file format must be:
        
        line1: title
        line2: data headings
        line3: space delimited data
        line4: space delimited data
        .
        .
        .
    """
    ### INPUTS #####
    if type(INPUTSfilename)==list:
        inFilename = INPUTSfilename[0]
    else:
        inFilename = INPUTSfilename
    if verbose:
        print 'Starting to load the input data'
    INfile = open(inFilename, 'r')
    lines = INfile.readlines()
    INtitle = lines[0]
    plotFileTitle = INtitle[:-5]
    INcolHeaders = lines[1]
    numLines = len(lines)
    
    # initiate input data lists
    longAN_degs = []
    es = []
    Ts = []
    periods = []
    inclination_degs = []
    argPeri_degs = []
    chiSquareds = []
    
    # loop through all the lines and extract data into data lists
    for row in range(2,numLines):
        dataLine = lines[row]
        dataLineCols = dataLine.split()
        dataLineColsFloats = []
        for data in dataLineCols:
            dataLineColsFloats.append(float(data))
        # put resulting float values of each column into its proper list
        longAN_degs.append(dataLineColsFloats[0])
        es.append(dataLineColsFloats[1])
        Ts.append(dataLineColsFloats[2])
        periods.append(dataLineColsFloats[3])
        inclination_degs.append(dataLineColsFloats[4])
        argPeri_degs.append(dataLineColsFloats[5])
        chiSquareds.append(dataLineColsFloats[6])
    if verbose:
        print 'Done loading Inputs from file'
    INfile.close()
    #################################
    
    # OUTPUTS    
    if type(OUTPUTSfilenames)==list:
        # initiate output data lists
        a1s2 = []
        a2s2 = []
        
        for file in range(0,len(OUTPUTSfilenames)):
            epoch = file+1
            curOutFilename = OUTPUTSfilenames[file]
            OUTfile = open(curOutFilename, 'r')
            lines = OUTfile.readlines()
            OUTtitle = lines[0]
            OUTcolHeaders = lines[1]
            numLines = len(lines)
    
            # initiate output data lists
            a1s = []
            a2s = []
            
            if file == 0:
                if verbose:
                    print 'Loading up output data of first epoch'
                # load first epoch's data into the lists
                # then update/append to these lists with later epochs.
                # loop through all the lines and extract data into data lists
                for row in range(2,numLines):
                    dataLine = lines[row]
                    dataLineCols = dataLine.split()
                    dataLineColsFloats = []
                    for data in dataLineCols:
                        dataLineColsFloats.append(float(data))
        
                    a1s2.append([dataLineColsFloats[0]])
                    a2s2.append([dataLineColsFloats[1]])
            elif file>0:
                if verbose:
                    print 'loading up output data of epoch ',epoch 
                for row in range(2,numLines):
                    dataLine = lines[row]
                    dataLineCols = dataLine.split()
                    dataLineColsFloats = []
                    for data in dataLineCols:
                        dataLineColsFloats.append(float(data))
                    
                    # update/append lists with this epoch's data
                    a1s2[row-2].append(dataLineColsFloats[0])
                    a2s2[row-2].append(dataLineColsFloats[1])
            if verbose:    
                print 'Done loading outputs for epoch '+str(epoch)+' from file'
            OUTfile.close()
        if verbose:    
            print '** All data from files loaded **'
        
    else:
        print 'OUTPUTSfilenames must be a list of output file name strings'
        
    return(plotFileTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs, a1s2, a2s2, chiSquareds)

#def orbElementPlotterLowRAM(plotFileTitle, longAN_degs, es, Ts, periods, inclination_degs, argPeri_degs,\
#                            Sep_Dists2, a1s2, a2s2, chiSquareds, showPlots=False, summaryOnly=True, verbose=False):
#    """
#    ************************************************************************************************************
#    *** THIS IS THE LOW RAM VERSION of orbElementPlotter for laptops or computers with less than 16GB of RAM ***
#    ************************************************************************************************************
#    
#    This function is used to plot the resulting orbital parameters of the orbits that satisfied 
#    the chi squared checks in a Monte Carlo simulator of the 'binary' code project.
#    
#    The plots will be be in multiple figures:
#    INS: all the standard orbitalCalculator inputs in separate histogram plots of value vs probability
#    OUTPUTS: all the standard orbitalCalculator outputs in separate histogram plots of value vs probability
#             with one figure for each epoch.
#    SUMMARY: a figure summarizing the orbital parameter histograms for the values quoted by Currie et all 2011.
#    
#    These figures will be saved to the directory '../figures/'.
#    CAUTION! the saved figures will overwrite a previous version of the file if it has the same name!!
#    
#    INPUTS:
#    @param plotFileTitle:                 = initial part of the figure and file title.  This will be 
#                                            have strings appended to it to create individual
#                                            input and output figure and file names.
#    @type plotFileTitle:                  = string
#    @param longAN_degs:                   = Longitude of Ascending Node [deg]
#    @type longAN_degs:                    = list of floats
#    @param es:                            = eccentricity of orbits [unitless]
#    @type es:                             = list of floats
#    @param Ts:                            = Last Periapsis Epoch/time [julian date] 
#    @type Ts:                             = list of floats
#    @param periods:                       = period of orbits [yrs]
#    @type periods:                        = list of floats
#    @param inclination_degs:              = inclination [deg]
#    @type inclination_degs:               = list of floats
#    @param argPeri_degs:                  = Argument of Periapsis in orbital plane [deg]
#    @type argPeri_degs:                   = list of floats
#    
#    @param Sep_Dists:                     = Separation Distances in orbital plane [AU]
#    @type Sep_Dists:                      = list of floats. Output values for each epoch.
#    @param a1s:                           = semi-major axis' of M1 [AU]
#    @type a1s:                            = list of floats. Output values for each epoch.
#    @param a2s:                           = semi-major axis' of M2 [AU]
#    @type a2s:                            = list of floats. Output values for each epoch.
#    @param chiSquareds:                   = total chi squared values for each orbit.
#    @type chiSquareds:                    = list of floats
#    
#    @param showPlots:                     = show the resulting figures?
#    @type showPlots:                      = Python boolean (True/False). Default = False.
#    @param summaryOnly:                   = Only save the summary figure?
#    @type summaryOnly:                    = Python boolean (True/False). Default = True.
#    @param verbose:                       = show prints on screen while running?
#    @type verbose:                        =  Python boolean (True/False). Default = False.
#    
#    """    
#    inclination_degsAllGooders = []
#    longAN_degsAllGooders = []
#    argPeri_degsAllGooders = []
#    esAllGooder = []
#    periodsAllGooders = []
#    TsAllGooders = []
#    a1AllGooders = []
#    a2AllGooders = []
#    
#    numSamples = np.shape(a1s2[:][:])[0] 
#    numEpochs = np.shape(a1s2[:][:])[1]
#    
#    for orbit in range(0,numSamples):
#        if chiSquareds[orbit]<=100.0:
#            inclination_degsAllGooders.append(inclination_degs[orbit])
#            esAllGooder.append(es[orbit])
#            longAN_degsAllGooders.append(longAN_degs[orbit])
#            periodsAllGooders.append(periods[orbit])
#            TsAllGooders.append(Ts[orbit])
#            argPeri_degsAllGooders.append(argPeri_degs[orbit])
#                
#    
#    ########################## INS #####################
#    # plot in's on first figure
#    if verbose:
#        print 'Starting to plot Inputs'
#    plt.figure(1,figsize=(20,7) ,dpi=300)
#    plt.subplot(231)
#    (n,bins,patches) = plt.hist(longAN_degs,bins=50)
#    plt.xlabel('longAN_degs')
#    plt.ylabel('Probability')
#    plt.subplot(232)
#    (n,bins,patches) = plt.hist(es,bins=50)
#    plt.xlabel('es')
#    plt.title(plotFileTitle+' *Inputs*')
#    plt.subplot(233)
#    Tsmin = int(np.min(Ts))
#    Tsmax = int(np.max(Ts))
#    TsNEW = Ts
#    for i in range(0,len(Ts)):
#        TsNEW[i]=(TsNEW[i]-Tsmin)/365
#    (n,bins,patches) = plt.hist(TsNEW,bins=25)
#    plt.xlabel('Ts [years since '+str(Tsmin)+']')
#    plt.subplot(234)
#    (n,bins,patches) = plt.hist(periods,bins=50)
#    plt.xlabel('periods')
#    plt.ylabel('Probability')
#    plt.subplot(235)
#    (n,bins,patches) = plt.hist(inclination_degs,bins=50)
#    plt.xlabel('inclination_degs')
#    plt.subplot(236)
#    (n,bins,patches) = plt.hist(argPeri_degs,bins=50)
#    plt.xlabel('argPeri_degs')
#    if summaryOnly is not True:
#        plt.savefig('../figures/'+plotFileTitle+'INS.png', dpi=300, orientation='landscape')
#    if verbose:
#        print 'Finished plotting Inputs'
#    
#    ################# OUTS ####################
#    if verbose:
#        print 'numSamples = ',numSamples
#        print 'numEpochs = ',numEpochs
#    
#    a1AllEpoch = []
#    a2AllEpoch = []    
#    
#    for epoch in range(0,numEpochs):
#        if verbose:
#            print 'Starting to plot outputs for epoch ',epoch
#        # reshape the output matrices
#        a1SingleEpoch = []
#        a1AllEpoch = []
#        a2SingleEpoch = []
#        a2AllEpoch = []
#    
#        for sample in range(0,numSamples):
#            a1SingleEpoch.append(a1s2[sample][epoch]) 
#            a1AllEpoch.append(a1s2[sample][epoch])
#            a2SingleEpoch.append(a2s2[sample][epoch])
#            a2AllEpoch.append(a2s2[sample][epoch]) 
#            if chiSquareds[sample]<=100.0:
#                 a1AllGooders.append(a1s2[sample][epoch])
#                 a2AllGooders.append(a2s2[sample][epoch])
#            
#        # plot OUT's on first figure
#        plt.figure(epoch+2, figsize=(16,9) ,dpi=300)
#        plt.subplot(122)
#        (n,bins,patches) = plt.hist(a1SingleEpoch,bins=50)
#        plt.xlabel('a1s')
#        plt.subplot(121)
#        (n,bins,patches) = plt.hist(a2SingleEpoch,bins=50)
#        plt.xlabel('a2s')
#        plt.ylabel('Probability')
#        plt.title(plotFileTitle+' *Outputs of simulator for epoch '+str(epoch+1)+'*')
#       
#        filename = '../figures/'+plotFileTitle+'_LowRAM_OUTSepoch'+str(epoch+1)
#        if summaryOnly is not True:
#            plt.savefig('%s.png' % filename, dpi=300, orientation='landscape')
#    
#    # make a summary plot for inputs (i, e, longAN) and output (a2), to compare to Currie et al results
#    if verbose:
#        print 'Starting to make Summary Plot'
#    plt.figure(numEpochs, figsize=(16,14) ,dpi=300)
#    plt.subplot(231)
#    (n,bins,patches) = plt.hist(longAN_degs,bins=50)
#    plt.xlabel('long of AN [deg]')
#    plt.ylabel('Probability')
#    plt.subplot(232)
#    (n,bins,patches) = plt.hist(es,bins=50)
#    plt.xlabel('e ')
#    plt.title(plotFileTitle+'  Currie et al. Comparison Summary')
#    (n,bins,patches) = plt.hist(a2AllEpoch,bins=50)#, range=(0.0,50.0)) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#    plt.subplot(233)
#    plt.xlabel('a1 [AU]')
#    plt.subplot(234)
#    (n,bins,patches) = plt.hist(a2AllEpoch,bins=50)#, range=(0.0,50.0)) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#    plt.xlabel('a2 [AU]')
#    plt.ylabel('Probability')
#    plt.subplot(235)
#    (n,bins,patches) = plt.hist(inclination_degs,bins=50)
#    plt.xlabel('inclination [deg]')
#    filename = '../figures/'+plotFileTitle+'_Comparison_Summary'
#    plt.savefig('%s.png' % filename, dpi=300, orientation='landscape')
#    if verbose:
#        print 'Finished making Summary Plot'
#        print 'median longAN =  ',np.median(longAN_degs)
#        print 'median e = ', np.median(es)
#        print 'median a2s = ', np.median(a2AllEpoch)
#        print 'median i = ',np.median(inclination_degs)
#
#    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#    ##$$$$$$$$$$$ TOTAL SUMMARY PLOT SECTION $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#    #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#    # make a summary plot for inputs (i, e, longAN) and output (a2), to compare to Currie et al results
#    if verbose:
#        print 'Starting to make Total Summary Plot'
#    plt.figure(numEpochs, figsize=(25,20) ,dpi=300)
#    plt.subplot(331)
#    (n,bins,patches) = plt.hist(inclination_degsAllGooders,bins=50)
#    plt.xlabel('inclination [deg]')
#    plt.ylabel('Probability')
#    plt.subplot(332)
#    (n,bins,patches) = plt.hist(longAN_degsAllGooders,bins=50)
#    plt.xlabel('long of AN [deg]')
#    plt.title(plotFileTitle+'  TOTAL CLEAN Summary')
#    plt.subplot(333)
#    (n,bins,patches) = plt.hist(argPeri_degsAllGooders,bins=50)
#    plt.xlabel('argPeri_degs')
#    plt.subplot(334)
#    (n,bins,patches) = plt.hist(esAllGooder,bins=50)
#    plt.xlabel('e ')
#    plt.ylabel('Probability')
#    plt.subplot(335)
#    (n,bins,patches) = plt.hist(periodsAllGooders,bins=50)
#    plt.xlabel('periods')
#    plt.subplot(336)
#    Tsmin = int(np.min(TsAllGooders))
#    Tsmax = int(np.max(TsAllGooders))
#    TsNEW = TsAllGooders
#    for i in range(0,len(TsAllGooders)):
#        TsNEW[i]=(TsNEW[i]-Tsmin)/365
#    (n,bins,patches) = plt.hist(TsNEW,bins=25)
#    plt.xlabel('Ts [years since '+str(Tsmin)+']')
#    plt.subplot(337)
#    (n,bins,patches) = plt.hist(a1AllGooders,bins=50)#, range=(0.0,50.0)) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#    plt.xlabel('a1 [AU]')
#    plt.ylabel('Probability')
#    plt.subplot(338)
#    (n,bins,patches) = plt.hist(a2AllGooders,bins=50)#, range=(0.0,50.0)) #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#    plt.xlabel('a2 [AU]')
#    
#    filename = '../figures/'+plotFileTitle+'_Total_Clean_Summary'
#    plt.savefig('%s.png' % filename, dpi=300, orientation='landscape')
#    if verbose:
#        print 'Finished making Total Summary Plot'
#        print 'median inclination_degsAllGooders =  ',np.median(inclination_degsAllGooders)
#        print 'median longAN_degsAllGooders = ', np.median(longAN_degsAllGooders)
#        print 'median argPeri_degsAllGooders = ', np.median(argPeri_degsAllGooders)
#        print 'median esAllGooder = ',np.median(esAllGooder)
#        print 'median periodsAllGooders =  ',np.median(periodsAllGooders)
#        print 'median TsNEW = ', np.median(TsNEW)
#        print 'median a1AllGooders = ', np.median(a1AllGooders)
#        print 'median a2AllGooders = ',np.median(a2AllGooders)
#        
#    # Show the plots on the screen if asked to    
#    if showPlots:
#        plt.show()        


#def orbitCalculatorLowRAM(t, Sys_Dist_PC, Sep_Angle_Measured_arcsec, inclination_deg, longAN_deg, e, T, period,\
#                      argPeri_deg, Mass1=1, Mass2=1, verbose=False):
#    """
#    This function will calculate the output orbital parameters from the inputs for a SINGLE epoch.  Use the 
#    multiEpochOrbCalc to do this for multiple epochs at the same time.  
#    The t, Sys_Dist_PC and Sep_Angle_Measured_arcsec parameters are generally known/measured values,
#    while the rest are typically random numbers.
#    This function is designed to work as the calculator inside a Monte Carlo Simulator.
#    All the equations are the result of Kepler's Laws and naturally occurring symmetries of binary star systems.
#    
#    See 'OrbitCalculator_2.pdf' for a full description of this function and its uses.
#    
#    Inputs:
#    @param t:                            = epoch of observation/image [julian date]
#    @type t:                             = float
#    @param Sys_Dist_PC:                  = measured system distance from Earth [PC]
#    @type Sys_Dist_PC:                   = float
#    @param Sep_Angle_Measured_arcsec:    = measured separation angle ["]
#    @type Sep_Angle_Measured_arcsec:     = float
#    @param inclination_deg:              = inclination [deg]
#    @type inclination_deg:               = float
#    @param longAN_deg:                   = Longitude of Ascending Node [deg]
#    @type longAN_deg:                    = float
#    @param e:                            = eccentricity of orbits [unitless]
#    @type e:                             = float
#    @param T:                            = Last Periapsis Epoch/time [julian date] 
#    @type T:                             = float
#    @param period:                       = period of orbits [yrs]
#    @type period:                        = float
#    @param argPeri_deg:                  = Argument of Periapsis in orbital plane [deg]
#    @type argPeri_deg:                   = float
#    @param Mass1:                        = Mass of primary in system [solar masses] (Default=1, means M2 is a planet)
#    @type Mass1:                         = float
#    @param Mass2:                        = Mass of the secondary in system [solar masses] (Default=1, means M2 is a planet)
#    @type Mass2:                         = float
#    @param verbose:                      = Send prints to screen? [True/False] (Default = False)
#    @type verbose:                       = Python boolean (True/False). Default = False
#    
#    Outputs:
#    @param Sep_Dist_AU_OP:               = Separation Distance of stars in orbital plane [AU]
#    @type Sep_Dist_AU_OP:                = float
#    @param PA_deg_measured:              = measured Position Angle in image [deg]
#    @type PA_deg_measured:               = float
#    @param a1:                           = semi-major axis of M1 [AU]
#    @type a1:                            = float
#    @param a2:                           = semi-major axis of M2 [AU]
#    @type a2:                            = float
#    
#    """
#    if verbose:
#        print '\n'+'*'*50
#        print 'Starting to calculate Orbital Parameters\n'
#        print 'Input variable values: '+\
#        '\nCurrent epoch time [julian date] = '+str(t)+\
#        '\nSystem Distance from Earth [PC] = '+str(Sys_Dist_PC)+\
#        '\nInclination [deg] = '+ str(inclination_deg)+\
#        '\nSeparation Angle ["] = '+str(Sep_Angle_Measured_arcsec)+\
#        '\nLongitude of Ascending Node [deg] = '+str(longAN_deg)+\
#        '\nEccentricity = '+str(e)+'\nLast Periapsis Epoch [julian date] = '+\
#        str(T)+'\nPeriod [yrs] = '+str(period)+\
#        '\nArgument of Periapsis in orbital plane [deg] = '+\
#        str(argPeri_deg)+'\nMass of primary [solar masses] = '+str(Mass1)+\
#        '\nMass of secondary [solar masses] = '+str(Mass2)
#    #----------------------------------------------------------------------
#
#    # initial guess (E_last), will be updated in loop.  
#    # Anything works, just takes longer if further from real value. => pi
#    E_last = 3.14
#    # stored initial value to be updated in loop
#    E_latest = 1.0
#    
#    #print '*** period recieved by orbitCalculator = ',period #######$$$$$$$$$$$$$$$$
#
#    ## calculate the Mean Motion
#    n = (2*math.pi)/period
#    if verbose:
#        print 'Mean Motion [rad/yr]= '+str(n)
#
#    ## calculate Mean Anomaly
#    M = n*((t-T)/365.0)
#    M_deg = math.degrees(M) # convert resulting M to degrees
#    if verbose:
#        print 'Mean Anomaly [deg]= ',M_deg
#
#    ### Performing Newton's Method to get the Eccentric Anomaly value ###
#    if verbose:
#        print '-'*50
#
#    # show input value to 
#    if verbose:
#        print "Inputs to Newton's Method are : \nM [rad]= "+str(M)+"\nE_last [rad]= "+\
#        str(E_last)+"\nE_latest [rad] = "+str(E_latest)+"\ne = "+str(e)
#        
#        print "\nStarting to run Newton's while loop."
#
#    count = 0 # a counter to stop inf loops in Newtons method below
#    while (abs(E_last-E_latest) > (1.0e-10))and(count<100):
#        if verbose:
#            print 'current E [rad]= ', E_latest
#        E_last = E_latest
#        E_latest = E_last - ((M-E_last+e*math.sin(E_last))/(e*math.cos(E_last)-1.0))
#        count = count+1
#
#    E_latest_deg = math.degrees(E_latest) # convert resulting E to degrees
#    if verbose:
#        print "The resultant E value is [deg] = ", E_latest_deg
#        print '-'*50
#    ### Newton's loop finished! ###
#    
#    ## calculate True Anomaly from Eccentric Anomaly
#    TA_rad  = math.acos((math.cos(E_latest)-e)/(1-e*math.cos(E_latest)))
#    TA_deg = math.degrees(TA_rad)
#    if verbose:
#        print 'True Anomaly [deg] = ', TA_deg
#        
#    # calculate measured Separation Distance in reference plane
#    Sep_Dist_measured_AU = Sep_Angle_Measured_arcsec*Sys_Dist_PC
#    if verbose:
#        print 'Separation Distance in the reference plane [AU] = ',Sep_Dist_measured_AU
#        
#    ## calculate measured Separation Distance in reference plane
#    inclination_rad = math.radians(inclination_deg) # convert inclination to radians 
#    
#    # angle between ascending node and M1
#    ang_AN_to_M1_deg = argPeri_deg + TA_deg
#    
#    ang_AN_to_M1_rad = math.radians(ang_AN_to_M1_deg) # convert angle to radians
#    
#    # Calculate the separation distance in the orbital plane
#    Sep_Dist_AU_OP = Sep_Dist_measured_AU/math.sqrt(math.pow(math.sin(ang_AN_to_M1_rad)*math.cos(inclination_rad),2.0) + math.pow(math.cos(ang_AN_to_M1_rad),2.0))
#    if verbose:
#        print 'Separation Distance in orbital plane [AU] = ',Sep_Dist_AU_OP
#    
#    # find horizontal and vertical WRT Line of Nodes (X and Y) components of Sep_Dist corrected for inclination
#    # Note: only the 'vertical (Y)' component will be effected by inclination
#    Sep_Dist_AU_RP_y = Sep_Dist_AU_OP*math.sin(ang_AN_to_M1_rad)*math.cos(inclination_rad)
#    Sep_Dist_AU_RP_x = Sep_Dist_AU_OP*math.cos(ang_AN_to_M1_rad)
#    # use Pythagorous theorem to find the final Sep Dist in the Ref Plane
#    # this is done only to varify math is right
#    Sep_Dist_measured_AU2 = math.sqrt(math.pow(Sep_Dist_AU_RP_x, 2.0)+math.pow(Sep_Dist_AU_RP_y,2.0))
#    Sep_Dist_DIFF = math.sqrt(math.pow((Sep_Dist_measured_AU-Sep_Dist_measured_AU),2.0))
#    if Sep_Dist_DIFF > 0.0001:
#        print 'The separation distances differ by more than 0.0001 !!!!!!'
#        print 'Sep_Dist_measured_AU = ', Sep_Dist_measured_AU
#        print 'Sep_Dist_measured_AU2 = ',Sep_Dist_measured_AU2
#    
#    ## calculate corrected angle between Line of Nodes and M1 (NOT AN to M1)
#    ang_LN_to_M1_corr = math.degrees(math.atan(Sep_Dist_AU_RP_y/Sep_Dist_AU_RP_x))
#    ## calculate measured Position Angle depending on which quadrant corrected M1 position is in
#    if Sep_Dist_AU_RP_x>=0.0:
#        PA_deg_measured_model = longAN_deg + 180.0 + ang_LN_to_M1_corr
#        #print '** corrected M1 found in 2nd or 3rd quadrant' #$$$$$$$$
#    elif Sep_Dist_AU_RP_x<0.0:
#        PA_deg_measured_model = longAN_deg + ang_LN_to_M1_corr
#        #print '** corrected M1 found in 1st or 4th quadrant' #$$$$$$$$$$$
#    if verbose:
#        print 'Position Angle in reference plane [deg] = ',PA_deg_measured_model
#                                                     
#    ## calculate the total semi-major axis (sum of both semi-majors for binary stars) in orbital plane
#    a_total = (Sep_Dist_AU_OP*(1+e*math.cos(TA_rad)))/(1-e*e)
#    if verbose:
#        print 'Total semi-major axis (a1+a2) [AU] = ',a_total
#    if Mass1!=1 and Mass2!=1:
#        # this means these two parameters are non-default values and 
#        # the system must then be a binary star system, thus calculate the
#        # a1 and a2 values of the system.
#        # NOTE: Mass2 MUST < Mass1, ie. primary is bigger.
#        
#        # first find the mass ratio of system
#        X = Mass2/Mass1
#        
#        # find a1 from the a1=X*a2 and a_total=a1+a2
#        a2 = a_total/(1.0+X)
#        
#        # now use a1=X*a2 to get a2
#        a1 = a2*X
#
#        if verbose: 
#            print 'The system must be a binary star system with a1 = '+str(a1)+' and a2 = '+str(a2)
#    else:
#        # default Mass1 and Mass2, thus a planet system and a_total = a2 and a1=0.0
#        # this is because we assume most planet's mass << star's mass, so a1<<a2
#        a1 = 0.0
#        a2 = a_total
#        
#    if verbose:
#        print '\nFinished calculating Orbital Parameters'
#        print '*'*50+'\n'
#   
#    return (Sep_Dist_AU_OP, PA_deg_measured_model, a1, a2)


#def multiEpochOrbCalcLowRAM(Sep_Angle_arcsec_measured_REALs, PA_deg_measured_REALs, PA_mean_errors, epochs,\
#                        Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, Mass1=1, Mass2=1, verbose=False):
#    """
#    This function will check if the suggested orbital parameter inputs satisfy all the input Position angles with a total 
#    chi squared less than chiSquaredMax.  Thus, it will check if they match all the observations of the system.
#    
#    
#    Inputs:
#    @param Sep_Angle_arcsec_measured_REALs:  = measured separation angle ["]
#    @type Sep_Angle_arcsec_measured_REALs:   = List of floats
#    @param PA_deg_measured_REALs:            = measured position angle [deg]
#    @type PA_deg_measured_REALs:             = List of floats
#    @param PA_mean_errors:                   = error in measured position angle [deg]. Use average if + and - are different.
#    @type PA_mean_errors:                    = List of floats
#    @param epochs:                           = epoch of observation/image [julian date]
#    @type epochs:                            = List of floats
#    @param Sys_Dist_PC:                      = measured system distance from Earth [PC]
#    @type Sys_Dist_PC:                       = float
#    @param inclination_deg:                  = inclination [deg]
#    @type inclination_deg:                   = float
#    @param longAN_deg:                       = Longitude of Ascending Node [deg]
#    @type longAN_deg:                        = float
#    @param e:                                = eccentricity of orbits [unitless]
#    @type e:                                 = float
#    @param T:                                = Last Periapsis Epoch/time [julian date] 
#    @type T:                                 = float
#    @param period:                           = period of orbits [yrs]
#    @type period:                            = float
#    @param argPeri_deg:                      = Argument of Periapsis in orbital plane [deg]
#    @type argPeri_deg:                       = float
#    @param Mass1:                            = Mass of primary in system [solar masses] (Default=1, means M2 is a planet)
#    @type Mass1:                             = float
#    @param Mass2:                            = Mass of the secondary in system [solar masses] (Default=1, means M2 is a planet)
#    @type Mass2:                             = float
#    @param verbose:                          = Send prints to screen? 
#    @type verbose:                           = Python boolean (True/False). Default = False
#    
#    Outputs:
#    @param chi_squared_total:                = Total chi squared the model's position angle for all epochs
#    @type chi_squared_total:                 = float
#    @param Sep_Dists:                        = Separation Distance in orbital plane [AU]
#    @type Sep_Dists:                         = list of floats. Output values for each epoch.
#    @param PA_deg_measured_models:           = measured Separation Angle of stars ["]
#    @type PA_deg_measured_models:            = list of floats. Output values for each epoch.
#    @param a1s:                              = semi-major axis of M1 [AU]
#    @type a1s:                               = list of floats. Output values for each epoch.
#    @param a2s:                              = semi-major axis of M2 [AU]
#    @type a2s:                               = list of floats. Output values for each epoch.
#    
#    Note: Sep_Angle_arcsec_measured_REALs, PA_deg_measured_REALs, PA_mean_errors inputs MUST all have
#          the same length.
#    """
#    
#    # init output lists
#    Sep_Dists = []
#    PA_deg_measured_models = []
#    a1s = []
#    a2s = []
#
#    # initial values for boolean equality of while statement 
#    chi_squared_total = 0.0
#    
#    for i in range(0,len(epochs)):
#
#        Sep_Angle_arcsec_measured_REAL = Sep_Angle_arcsec_measured_REALs[i]
#        PA_deg_measured_REAL = PA_deg_measured_REALs[i]
#        PA_mean_error = PA_mean_errors[i]
#        t = epochs[i]  
#        if verbose:
#            print '$$$$$$ epoch for counter '+str(i)+' out of '+str(len(epochs))+' is '+str(epochs[i])
#        
#        # call orbitCalculator to take random variables and calc orbital elements
#        (Sep_Dist_AU, PA_deg_measured_model, a1, a2) = \
#        orbitCalculatorLowRAM(t, Sys_Dist_PC, Sep_Angle_arcsec_measured_REAL, inclination_deg, longAN_deg, e, T, period, argPeri_deg, Mass1=Mass1, Mass2=Mass2, verbose=verbose)
#        
#        # store output orbital elements in lists for later plotting
#        Sep_Dists.append(Sep_Dist_AU)
#        PA_deg_measured_models.append(PA_deg_measured_model)
#        a1s.append(a1)
#        a2s.append(a2)
#            
#        ## calc both PA and SA kai's and check they are less than chiSquaredMax
#        PA_chi_squared = chiSquaredCalc(PA_deg_measured_REAL, PA_mean_error, PA_deg_measured_model)
#        # Ad them to get the updated total
#        chi_squared_total = chi_squared_total + PA_chi_squared # + SA_chi_squared
#        
#    return (chi_squared_total, Sep_Dists, PA_deg_measured_models, a1s, a2s)