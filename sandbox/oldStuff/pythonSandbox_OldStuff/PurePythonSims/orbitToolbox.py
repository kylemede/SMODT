#@Author: Kyle Mede, kylemede@gmail.com, written for the University of Tokyo during UTRIP 2011.
import math
import gc
import numpy as np
import os
import pylab
from math import pi
plt = pylab.matplotlib.pyplot
patches = pylab.matplotlib.patches

"""
This toolbox is a collection of the calculator type functions that were used in multiple 
places throughout the code to conduct various types of binary star system simulations.
"""
    
def bestOrbitFinderNEW(filename, printToScreen=True, saveToFile=True, returnAsList=False):
    """
    Just a simple function to find the parameters for the lowest chiSquared orbit in a 
    data file of the NEW format.
    """
  
    lowestChiSquared = 100000
    inclinationBest = 0
    eBest = 0
    longANBest = 0
    periodBest = 0
    argPeriBest = 0
    aBest = 0
    TBest = 0
    rv_origin_0_best = 0
    rv_origin_1_best = 0
    
    offsets = False
    
    # strip off the .txt part to make the plot version of the filename
    print "\nWorking on file: "+os.path.basename(filename)
    f = open(filename, 'r')
    plotFileTitle = f.readline()[:-5]
    headings = f.readline()
    line = 'asdf'
    while line!='':
        line = f.readline()
        if line!='':
            dataLineCols = line.split()
            chiSquared = float(dataLineCols[7])
            
            if chiSquared<lowestChiSquared:
                lowestChiSquared = chiSquared
                incBest = float(dataLineCols[4])
                eBest = float(dataLineCols[1])
                longANBest = float(dataLineCols[0])
                periodBest = float(dataLineCols[3])
                argPeriBest = float(dataLineCols[5])
                aBest = float(dataLineCols[6])
                TBest = float(dataLineCols[2])
                if len(dataLineCols)>9:
                    offsets = True
                    rv_origin_0_best = float(dataLineCols[8])
                    rv_origin_1_best = float(dataLineCols[9])   
        
    # print the values for the best orbit
    line= '\nBest orbit found:'
    line=line+ "\nLongAN = "+str(longANBest)
    line=line+ "\ne = "+str(eBest)
    line=line+ "\nTo = "+str(TBest)
    line=line+ "\nperiod = "+str(periodBest)
    line=line+ "\ninclination = "+str(incBest)
    line=line+ "\nargPeri = "+str(argPeriBest)
    line=line+ "\na_total = "+str(aBest)
    if offsets:
        line = line+ "\nRV offset 0 ="+str(rv_origin_0_best)
        line = line+ "\nRV offset 1 ="+str(rv_origin_1_best)
    line=line+ "\nchiSquaredMin = "+str(lowestChiSquared)
    
    if printToScreen:
        print line
    
    if saveToFile:
        outFilename = os.path.dirname(filename)+'/bestOrbit.txt'
        if os.path.exists(outFilename):
            print '\nbestOrbFinderNEW68: Warning: '+outFilename+' all ready exists, so not overwriting.'
        else:
            f = open(outFilename, 'w')
            f.write(line)
            f.close()
            print '\nbestOrbFinderNEW68: Best orbit values saved to: '+outFilename
            
    if returnAsList:
        if offsets:
            list = [longANBest, eBest, TBest, periodBest, incBest, argPeriBest, aBest, rv_origin_0_best, rv_origin_1_best]
        else:
            list = [longANBest, eBest, TBest, periodBest, incBest, argPeriBest, aBest]
        return list 
    
def burnInCalc3(chiSquareds, medianALLchains):
    """
    This function will calculate the burn in length and return its value.
    This can only truly be done for a proper MCMC that was started
    at a random location in the parameter space, as anything else will
    skew the chiSquareds trend.  ie. they will not roughly decrease as the 
    simulation runs and the burn in length will thus be more or less random and 
    useless information.
    
    @param paramIN:    = parameter array including burn in data
    @type paramIN:     = array (list) of doubles
    """
    burnInLength = 0
    
    if len(paramIN)>=100000:
        jump = 100
    else:
        jump = 10
    
    for i in range(0,int(len(paramIN)/jump)):
        # check if current chiSquared is less than the median of all chains yet
        chiSquareCur = chiSquareds[i*jump]
        if stdCur<medianALLchains:
            # less than median, so do all in last jump to find precise location
            for j in range((i-1)*jump, i*jump):
                #print 'inside second loop'#$$$$$$$$$$$$$$$$$$$$
                chiSquareCur2 = chiSquareds[j]
                if chiSquareCur2<medianALLchains:
                    burnInLength = j+1
                    break
            break
    if burnInLength == len(paramIN):
        print "PROBLEM: Param had a burn in length equal to param length, ie. the chain never burned in"
        
    return burnInLength

def CorrLengthCalc(paramIN):
    """
    This function will calculate the correlation length and return its value.
    
    @param paramIN:    = parameter array after burn in data stripped
    @type paramIN:     = array (list) of doubles
    """
    
    try:
        stdALL = np.std(paramIN)
    except:
        useless=0
    halfStdALL = stdALL/2.0
    CorrLength = 0
    
    if len(paramIN)>=100000:
        jump = 100
    else:
        jump = 10
    
    for i in range(1,int(len(paramIN)/jump)):
        #check std at each jump to see if over halfstd yet
        try:
            stdCur = np.std(paramIN[0:i*jump])
        except:
            useless=1
        if stdCur>halfStdALL:
            # over halfStd, so do all in last jump to find precise location
            for j in range((i-1)*jump, i*jump):
                try:
                    stdCur2 = np.std(paramIN[0:j])
                except:
                    useless=2
                if stdCur2>halfStdALL:
                    CorrLength = j+1
                    break
            break
    if CorrLength == len(paramIN):
        print "PROBLEM: Param had a correlation length equal to param length, ie. the chain never burned in"
        
    return CorrLength

def burnInStripper(fullFilename, burnInLength, burnInStrippedFilename):
    """
    A function to read and strip the burn-in from a file and write
    an output file with the burn-in stripped off the data.
    """
    print "\nStarting to strip burn-in off file: "+os.path.basename(fullFilename)
    # open input file
    inputFile = open(fullFilename, 'r')
    # prepare output file
    outputFile = open(burnInStrippedFilename,'w')
    # strip off the .txt part to make the plot version of the filename
    plotFileTitle = inputFile.readline()[:-5]
    outputFile.write(plotFileTitle+' with its max burn-in length of '+str(burnInLength)+' stripped off\n')
    headings = inputFile.readline()
    outputFile.write(headings)
    
    # go through all lines in input file 
    line = 'asdf'
    currLength = 0
    reachedBurnIn = False
    while line!='':
        line = inputFile.readline()
        if reachedBurnIn is False:
            # still haven't reached the burn-in length, so check length and write
            # to output file if this particular line is the first to reach the it
            dataLineCols = line.split()
            timesBeenHere = int(dataLineCols[-1])
            if timesBeenHere>1:
                # This loop will go through the number of times 
                # the simulator stayed at that step/orbit
                for i in range(0,timesBeenHere): 
                    currLength +=1
                    if currLength>burnInLength:
                        # create a new version of the line that has a timesBeenHere=1
                        line = dataLineCols[0]
                        for col in range(1,(len(dataLineCols)-1)):
                            line = line+'    '+dataLineCols[col]
                        line = line+'      '+str(1)+'\n'
                        # write updated line to output file
                        outputFile.write(line)
                        reachedBurnIn = True
            else:
                # Thus this line has a timesBeenHere=1, so line doesn't need to be looped or modified
                if currLength>burnInLength:
                    outputFile.write(line)
                    reachedBurnIn = True 
                else:
                    currLength +=1                
        else:
            # burn-in length reached, so just write line to output file
            outputFile.write(line)
                
    #Finished writing burn-in stripped data so close files
    inputFile.close()
    outputFile.close()
    
    print 'The burn-in was stripped off and final data written to '+burnInStrippedFilename

def chiSquaredCalc(real, error, model):
    """ Just a simple function to calculate chi**2 for a given observed (real) with error
        and experimental (model) value.
        
        Note: Ensure the units of both are the same.
        Note2:  For cases where the error of the observation has different errors in the 
                positive and negative from the observed value, just input the mean of these.
        Note3: There should never be negative values input into this function for any
                of the three inputs, but it has been made capable of handling those cases as well.
                
        output: Resulting chi**2 value following
        
        chi**2 = (obs - model)**2 / error**2
    """
    # since the numerator is the difference, we need to take sign of each into account
    # NOTE: all these issues is now handled by the diffCalc function.
    #difference = diffCalc(model, real)
    difference = model-real
    
    #print "difference = ",difference
    #print "error = ",error
    #print "1/err**2",(1.0/(error**2.0))
    
    # put it all together to make final chi**2 for each for this model    
    chi_squared = (difference**2.0)/(error**2.0)
    
    return chi_squared

def dataFileCombiner(filenames,outFilename):
    '''
    This function is to combine files written to disk by the dataWriter function and 
    combine them into one large data file.  These files must all be of the same epoch, 
    or all be 'input' files.  
    
    Note:
    output file will have the 1st line be the same as the outFilename parameter and
     the second line will be the same headers as the first file in the 
     filenames parameter has.
     All files must have the same number of columns.
     
     
    @param filenames:    = list of filenames to be put together
    @type filenames:     = list of strings of form 'blahblah.txt'
    @param outFilenames: = new filename for combined output file
    @type outFilenames:  = string of form 'blahblahout.txt'
    '''
    
    # check filenames parameter
    if type(filenames) is not list:
        print 'the filenames input parameter must be a list of strings'
    
    # check the outFilename parameter
    if type(outFilename) is not str:
        print 'the outFilename input parameter must be a single string'
    else:
        # checking if it has '.txt' on the end, if not add it
        if outFilename[-4:]!='.txt':
            print 'Changing output filename from '+outFilename+' to '+outFilename+'.txt'
            outFilename = outFilename+'.txt'
    
    numFiles = len(filenames)
   
    # First load first file into memory
    fileOne = open(filenames[0], 'r')
    lines = fileOne.readlines()
    INtitle = lines[0]
    colHeaders = lines[1]
    numColumns = len(lines[2].split(' '))
    numLines = len(lines)
    fileOne.close()
    
    # replace first files title header and rest into a new output line list
    outLines = lines
    outLines[0] = os.path.basename(outFilename)+'\n' # all 
    
    # start output file and load up with outLines so far
    fOUT = open(outFilename, 'w')
    for line in outLines:
        fOUT.write(line)
        
    # Now load the other files into memory and add their lines, 
    # without the first two header lines
    for file in filenames[1:]:
        f = open(file,'r')
        lines = f.readlines()
        if numColumns!=len(lines[2].split(' ')):
            print 'PROBLEM: One of the files has a different number of columns than the first file'
            break
        else:
            for line in lines[2:]:
                fOUT.write(line)
        f.close()
    
    # close outfile as it is now loaded up
    fOUT.close()     
    
    print 'Done! '+str(numFiles)+' files combined and written to "'+outFilename+'"'

def dataReaderNEW(filename, column=0):
    """ filename must be a string of format: 'blahblah.txt'
        column = which column you want to have the data for, zero indexed.
        
        NOTE: This version is meant for the new headings from mcmcOrbSimUniform6
        longAN [deg]      e [N/A]       To [julian date]  period [yrs]   inclination [deg]   argPeri [deg]   a_total [AU]    chiSquared    timesBeenHere
        
        file format must be:
        
        line1: title
        line2: data headings
        line3: space delimited data
        line4: space delimited data
        .
        .
        .
    """
    ## instantiate data list
    data = []
   
    print "\nWorking on file: "+os.path.basename(filename)
    f = open(filename, 'r')
    # strip off the .txt part to make the plot version of the filename
    plotFileTitle = f.readline()[:-5]
    headings = f.readline()
    line = 'asdf'
    while line!='':
        line = f.readline()
        if line!='':
            dataLineCols = line.split()
            dataValue = dataLineCols[column]
            timesBeenHere = int(dataLineCols[-1])
        
            # This loop will append the value to the list based on the number of times 
            # the simulator stayed at that step/orbit
            for i in range(0,timesBeenHere): 
                data.append(float(dataValue))
                
    print "number of values in data column = "+str(len(data))
    return data

def dataReadAndPlotNEW3(filename, plotFilename, weight=True, confLevels=True, normalize=True, drawLine=False, numEpochs=1):
    """
    Master plotting function for the singleParamReadAndPlotNEW
    
    file format must follow the 'NEW' style, ie:
    longAN [deg]      e [N/A]       To [julian date]  period [yrs]   inclination [deg]   argPeri [deg]   a_total [AU]    chiSquared    timesBeenHere
    
    All 3D data should have its chiSquared values all ready be reduced.
    """
    numBins = 50
    
    # check if the passed in value for filename includes '.txt'
    if '.txt' not in filename:
        filename = filename+'.txt'
        print ".txt was added to the end of the filename as it wasn't found to be in the original version"
    
    # check if the passed in value for plotFilename includes '.png'
    if '.png' not in plotFilename:
        plotFilename = plotFilename+'.png'
    
    # create a plot title from the filename
    plotFileTitle = os.path.basename(plotFilename).split('.')[0]

    # Call the function to find the best orbit values for reference.
    print '#'*50
    bestOrbitFinderNEW(filename, printToScreen=True, saveToFile=True)
    print '#'*50
    
    print '\nStarting to create summary plot\n'
    
    # Create empty figure to be filled up with plots
    fig = plt.figure(1, figsize=(25,10) ,dpi=300) #plt.figure(1, figsize=(25,22) ,dpi=300)
        
    # Create sub plot and fill it up for the Inclinations
    subPlot = fig.add_subplot(331)
    paramColNum = 4
    subPlot.axes.set_xlabel('Inclination [deg]')
    subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, numEpochs)
    print "Done plotting inclination_degsAlls"
    
    # Create sub plot and fill it up for the Longitude of Ascending Node
    subPlot = fig.add_subplot(332)
    paramColNum = 0
    subPlot.axes.set_xlabel('Longitude of Ascending Node [deg]')
    subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, numEpochs)
    print "Done plotting longAN_degsAlls"
    
    # Create sub plot and fill it up for the Argument of Perigie
    subPlot = fig.add_subplot(333)
    paramColNum = 5
    subPlot.axes.set_xlabel('Argument of Perigie [deg]')
    subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, numEpochs)
    print "Done plotting argPeri_degsAlls"
    
    # Create sub plot and fill it up for the e
    subPlot = fig.add_subplot(334)
    paramColNum = 1
    subPlot.axes.set_xlabel('e')
    subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, numEpochs)
    print "Done plotting esAlls"
    
    # Create sub plot and fill it up for the Period
    subPlot = fig.add_subplot(335)
    paramColNum = 3
    subPlot.axes.set_xlabel('Period [Years]')
    subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, numEpochs)
    print "Done plotting periodsAlls"
    
    # Create sub plot and fill it up for the Semi-major
    subPlot = fig.add_subplot(336)
    paramColNum = 6
    subPlot.axes.set_xlabel('Semi-major Axis [AU]')
    subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, numEpochs)
    print "Done plotting Semi-major Axis"
    
    ## Create sub plot and fill it up for the Time of last Periapsis
    #subPlot = fig.add_subplot(337)
    #paramColNum = 2
    #subPlot.axes.set_xlabel('Time of last Periapsis')
    #subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, numEpochs)
    #print "Done plotting TsAlls"
    
    # Save file 
    plt.savefig(plotFilename, dpi=300, orientation='landscape')
    print 'Summary figure saved to '+plotFilename
        
    plt.close()
    
    ## Create a second figure of RV offsets. ####
    try:
        # Create empty figure to be filled up with plots
        fig = plt.figure(2, figsize=(35,22) ,dpi=250)
        
        # Create sub plot and fill it up for the Semi-major
        subPlot = fig.add_subplot(211)
        paramColNum = 8
        subPlot.axes.set_xlabel('RV offset 0 [m/s]')
        subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, numEpochs)
        print "Done plotting RV offsets for dataset 0"
        
        # Create sub plot and fill it up for the Semi-major
        subPlot = fig.add_subplot(212)
        paramColNum = 9
        subPlot.axes.set_xlabel('RV offset 1 [m/s]')
        subPlot = singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight, confLevels, normalize, drawLine, numEpochs)
        print "Done plotting RV offsets for dataset 1"
        
        # Save file 
        plotFilename2 = plotFilename[0:-4]+'-RVoffsets.png'
        plt.savefig(plotFilename2, dpi=300, orientation='landscape')
        print 'RV offsets summary figure saved to '+plotFilename2
            
        plt.close()
    except:
        plt.close()
        print 'No RV offsets to plot'
    
def diffCalc(A, B):
    """
    Simply returns the difference between two numbers (floats) as a positive number (float).
    """
    
    # to take sign of each into account
    #if ((A>=0.0)and(B>=0.0))or((A<0.0)and(B<0.0)):
    if (A*B)>=0.0:
        # same sign so just subtract, and take abs to clear any resulting negative
        difference = abs(A - B)
    #if ((A>=0.0)and(B<=0.0))or((A<0.0)and(B>0.0)):
    elif (A*B)<0.0:
        # different signs so add the abs of each
        difference = abs(A) + abs(B)
        
    return difference

def proposalMaxMinsCalc(latestVal, sigma, Min, Max):
    """
    Find max and min values for proposing a new random parameter value when doing MCMC.
    
    OLD simple version:
    
    if (A_latest+A_sigma)>=Amax:
        max = Amax
        min = A_latest-A_sigma
    elif (A_latest-A_sigma)<=Amin:
        max = A_latest+A_sigma
        min = Amin
    else:
        max = A_latest + A_sigma
        min = A_latest - A_sigma
    """
    if (latestVal+sigma)>=Max:
        max = Max
        min = latestVal-sigma
    elif (latestVal-sigma)<=Min:
        max = latestVal+sigma
        min = Min
    else:
        max = latestVal + sigma
        min = latestVal - sigma
    
    return (min, max)

def singleParamReadAndPlotNEW(filename, paramColNum, subPlot, numBins, weight=True, confLevels=True, normalize=True, drawLine=False, numEpochs=1):
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
    # set y axis title
    subPlot.axes.set_ylabel('Probability')
    axs = subPlot.axis()
    
    f = open(filename,'r')
    plotFileTitle = f.readline()[:-5]
    headings = f.readline()
    dataMax = 0.1
    dataMin = 10000.0
    chiSquaredMin = 1000000.0
    line = 'asdf'
    while line!='':
        line = f.readline()
        if line!='':
            dataLineCols = line.split()
            dataValue = float(dataLineCols[paramColNum])
            chiSquared = float(dataLineCols[7])
            if dataValue>dataMax:
                dataMax = dataValue
            if dataValue<dataMin:
                dataMin = dataValue
            if chiSquared<chiSquaredMin:
                chiSquaredMin = chiSquared
    f.close()
    
    ## set up axis of plot to match data
    xAxisWidth = dataMax-dataMin
    axisXmin = dataMin-(0.05*xAxisWidth)
    axisXmax = dataMax+(0.05*xAxisWidth)
    axisYmin = 0.0
    axisYmax = 1.1
    subPlot.axis([axisXmin, axisXmax, axisYmin, axisYmax])
    
    # convert the chiSquaredMin found to the reduced version if not done yet.
    if numEpochs>1:
        # ie, the input data uses non-reduced chiSquared values so fix that
        chiSquaredMin = (1.0/((2.0*numEpochs)-7.0))*chiSquaredMin
    
    #print 'dataMin = ',dataMin
    #print 'dataMax = ',dataMax
    #print 'chiSquaredMin = ',chiSquaredMin
    
    ## set up data range using a small buffer
    dataRange = (dataMax-dataMin)*1.05
    dataMin2 = dataMax-(dataRange*1.025)
    dataMax2 = dataMin+(dataRange*1.025)
    binWidth = dataRange/numBins
    #print 'binWidth = ',binWidth
    
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
        curBinChiSquaredMAX = 0.0
        curBinChiSquaredMIN = 10000.0
        
        line = 'asdf'
        while line!='':
            line = f.readline()
            if line!='':
                dataLineCols = line.split()
                dataValue = float(dataLineCols[paramColNum])
                chiSquared = float(dataLineCols[7])
                
                #print 'chiSquared = ',chiSquared
                if numEpochs>1:
                    # ie, the input data uses non-reduced chiSquared values so fix that
                    chiSquared = (1.0/((2.0*numEpochs)-7.0))*chiSquared
                #print 'chiSquared ',chiSquared #$$$$$$$$$$$
                likelihood = 0
                # see if value is inside bin, not including max
                if bin==(numBins):
                    # take care of last bin including it's max value
                    if (dataValue>=binMins[-1]) and (dataValue<=binMaxs[-1]):
                        if weight:
                            likelihood = math.exp(-chiSquared/2.0)#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                        else:
                            likelihood = 1
                        if chiSquared>curBinChiSquaredMAX:
                            #print 'updating chiSquaredMAX to ',chiSquared
                            curBinChiSquaredMAX = chiSquared 
                        if chiSquared<curBinChiSquaredMIN:
                            #print 'updating chiSquaredMIN to ',chiSquared
                            curBinChiSquaredMIN = chiSquared 
                else:
                    # still not last bin, so don't include max value
                    if (dataValue>=binMins[-1]) and (dataValue<binMaxs[-1]):
                        if weight:
                            likelihood = math.exp(-chiSquared/2.0)#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                        else:
                            likelihood = 1
                        if chiSquared>curBinChiSquaredMAX:
                            #print 'updating chiSquaredMax to ',chiSquared
                            curBinChiSquaredMAX = chiSquared 
                        if chiSquared<curBinChiSquaredMIN:
                            #print 'updating chiSquaredMIN to ',chiSquared
                            curBinChiSquaredMIN = chiSquared 
                # update bin value with which ever case was satisfied
                binVal = binVal + likelihood
                
                # add bin mid value to total list to find median later
                if likelihood>0.0:
                    binMidsALL.append(binMids[-1])  ###$$$$$$$$$$$$$$ could be dangerous for long data sets!!!!!!!!!!!!!!!!!!!!!!
                    
        if confLevels:
            #print 'curBinChiSquaredMax', curBinChiSquaredMax #$$$$$$$$$$$$$$$$$$
            #print 'chiSquaredMin',chiSquaredMin  #$$$$$$$$$$$$$$$$$$
            # is there anything in the bin?
            if binVal>0.0:
                c = 'w'
                ## set color of bin depending on what conf level the data is in
                # Check if the data is in the 2sigma conf level
                if (chiSquaredMin+4.0)>curBinChiSquaredMIN:
                    c = '0.5'
                # Check if the data is in 1sigma conf level
                if (chiSquaredMin+1.0)>curBinChiSquaredMIN:
                    c = '0.3'
                    
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
    
    # return fully updated sub plot
    return subPlot
    
def dataReadTrimWriteNEW(filename, chiSquareCutOff, verbose=False):
    """
    This function is ideal for the cases where your data sets are too big because of choosing a chiSquareMax that is
    too big when running the simulator, resulting in large data files that are hard to load and plot.  In these cases
    you can just find out what the standard range of chiSquareds are, and choose a lower cut-off to trim the volume
    of data written to disk.  Just open the file and skim over the chiSquared column values to choose a new
    cut-off.
    These circumstances really only happen when you have a long 'daisy chain' run with rough settings, so this is 
    not expected to be a commonly used function.
    
    NOTE: output files will have same name with 'chiSquare-cut-off-####' added to 
          show it is the new trimmed version.
    
    @param filename:             = filename
    @type filename:              = string
    @param chiSquareCutOff:      = max value of chiSquared for orbits to be written
    @type chiSquareCutOff:       = any number
    """
    
    # check if the passed in value for filenameROOT includes '.txt'
    if '.txt' not in filename:
        INSfilename = filename+'.txt'
    else:
        INSfilename = filename

    titleMod = '-chiSquare-cut-off-'+str(chiSquareCutOff)
    # create output version of INS file and write new 
    # strip off '.txt', add titleMod, then put it back on
    newFilename = INSfilename[:-4]+titleMod+'.txt'
    
    if verbose:
        print 'Starting to load, trim and write the data'
    # open INS file and grab title and headings
    INfile = open(INSfilename, 'r')
    plotFileTitle = INfile.readline()[:-5]
    headings = INfile.readline()
    
    OUTfile = open(newFilename, 'w')
    OUTfile.write(plotFileTitle+titleMod+'\n')
    OUTfile.write(headings)
    # loop through data lines and only write those that
    # have chiSquare below cut-off and store all 
    # chiSquareds in a list
    chiSquareds = []
    line='asdf'#just a value that isn't '' to start, will be replaced below
    totalKept = 0
    while line!='':
        line = INfile.readline()
        if line!='':
            dataLineCols = line.split()
            if len(dataLineCols)==10:
                chiSquareStr = dataLineCols[8]
            elif len(dataLineCols)==9:
                chiSquareStr = dataLineCols[7]
            chiSquared = float(chiSquareStr)
            if chiSquared<chiSquareCutOff:
                OUTfile.write(line)
                totalKept = totalKept+1
    if verbose:
        print str(totalKept)+' orbits were found to have a chiSquared < '+str(chiSquareCutOff)
        print 'Done loading, trimming and writing. Output file = '+newINSfilename
    INfile.close()
    OUTfile.close()    
    
    print 'Final trimmed data written to: '+newFilename
    
    if verbose:    
        print '** All data from file trimmed and written to output file **'      

def dataReadTrimWriteNEW2(filename, chiSquareCutOff, verbose=False):
    """
    ### Same as dataReadTrimWriteNEW BUT with a inclination trimming feature that isn't really needed anymore ######
    
    
    This function is ideal for the cases where your data sets are too big because of choosing a chiSquareMax that is
    too big when running the simulator, resulting in large data files that are hard to load and plot.  In these cases
    you can just find out what the standard range of chiSquareds are, and choose a lower cut-off to trim the volume
    of data written to disk.  Just open the file and skim over the chiSquared column values to choose a new
    cut-off.
    These circumstances really only happen when you have a long 'daisy chain' run with rough settings, so this is 
    not expected to be a commonly used function.
    
    NOTE: output files will have same name with 'chiSquare-cut-off-####' added to 
          show it is the new trimmed version.
    
    @param filename:             = filename
    @type filename:              = string
    @param chiSquareCutOff:      = max value of chiSquared for orbits to be written
    @type chiSquareCutOff:       = any number
    """
    
    # check if the passed in value for filenameROOT includes '.txt'
    if '.txt' not in filename:
        INSfilename = filename+'.txt'
    else:
        INSfilename = filename

    titleMod = '-chiSquare-cut-off-'+str(chiSquareCutOff)
    # create output version of INS file and write new 
    # strip off '.txt', add titleMod, then put it back on
    newFilename = INSfilename[:-4]+titleMod+'.txt'
    
    if verbose:
        print 'Starting to load, trim and write the data'
    # open INS file and grab title and headings
    INfile = open(INSfilename, 'r')
    plotFileTitle = INfile.readline()[:-5]
    headings = INfile.readline()
    
    OUTfile = open(newFilename, 'w')
    OUTfile.write(plotFileTitle+titleMod+'\n')
    OUTfile.write(headings)
    # loop through data lines and only write those that
    # have chiSquare below cut-off and store all 
    # chiSquareds in a list
    chiSquareds = []
    line='asdf'#just a value that isn't '' to start, will be replaced below
    totalKept = 0
    while line!='':
        line = INfile.readline()
        if line!='':
            dataLineCols = line.split()
            if len(dataLineCols)==10:
                chiSquareStr = dataLineCols[8]
            elif len(dataLineCols)==9:
                chiSquareStr = dataLineCols[7]
            chiSquared = float(chiSquareStr)
            inclination = float(dataLineCols[4])
            if chiSquared<chiSquareCutOff:
                if inclination<90.0:
                    OUTfile.write(line)
                    totalKept = totalKept+1
    if verbose:
        print str(totalKept)+' orbits were found to have a chiSquared < '+str(chiSquareCutOff)
        print 'Done loading, trimming and writing. Output file = '+newINSfilename
    INfile.close()
    OUTfile.close()    
    
    print 'Final trimmed data written to: '+newFilename
    
    if verbose:    
        print '** All data from file trimmed and written to output file **'

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

def dictToFile(d, filename):
    """
    This function will take an simple or complex/nested dictionary
    and write it's values to a text file for future reference.
    
    @param d:     = dictionary to be written to file
    @type d:      = standard Python dictionary, either simple or nested 
                    (max 3 layers ie. dict inside a dict inside a dict)
    """
    
    # check if the passed in value for filenameROOT includes '.txt'
    if '.txt' not in filename:
        filename = filename+'.txt'
        #print 'PROBLEM: there is .txt in the filename and there should not be'
        
    f = open(filename, 'w')
    
    for parName in d:
        if type(d[parName])==dict:
            # ie a singly nested dict
            f.write('\n'+parName+':')
            for subParName in d[parName]:
                if type(d[parName][subParName])==dict:
                    #ie a doubly nested dict
                    f.write('\n    '+subParName+':')
                    for subSubParName in d[parName][subParName]:
                        f.write('\n         '+subSubParName+' : '+str(d[parName][subParName][subSubParName]))
                else:
                    f.write('\n    '+subParName+' : '+str(d[parName][subParName]))
            
        else:
            #ie a simple dict
            f.write('\n'+parName+' : '+str(d[parName]))
    
    f.close()
    
def ellipse(ra,rb,ang,x0,y0,Nb=50):
    '''
    For calculating Nb x and y values that lie on an ellipse 
    defined using ra, rb, ang, xo and yo.
    
    @param ra:         = Semi-Major axis length
    @type ra:          = any number
    @param rb:         = Semi-Minor axis length
    @type rb:          = any number
    @param ang:        = Angle [degrees] to plot the semi-major axis 
                         wrt horizontal, measured from positive x-axis
    @type ang:         = any number
    @param x0:         = X coordinate of center of ellipse
    @type xo:          = any number
    @param y0:         = Y coordinate of center of ellipse
    @type y0:          = any number
    @param Nb:         = Number of x,y points to return 
    @type Nb:          = int
    
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
    
def findTop20orbitsNEW(filename):
    """
    This function will hunt through a data file in the NEW format and display the top 20
    orbits in the format that is used in the 'paramSettingsDict'.
    
    NOTE: format of data file must be:
    longAN [deg]      e [N/A]       To [julian date]  period [yrs]   inclination [deg]   argPeri [deg]   a_total [AU]    chiSquared    timesBeenHere
    """

    # check if the passed in value for filename includes '.txt'
    if '.txt' not in filename:
        filename = filename+'.txt'
        print ".txt was added to the end of the filename as it wasn't found to be in the original version"
    
    ## FIRST: get a list of ALL chiSquareds
    chiSquareds = []
    f = open(filename,'r')
    # strip off the .txt part to make the plot version of the filename
    plotFileTitle = f.readline()[:-5]
    headings = f.readline()
    line = 'asdf'
    while line!='':
        line = f.readline()
        if line!='':
            dataLineCols = line.split()
            curChiSquared = float(dataLineCols[7])
            chiSquareds.append(curChiSquared)
    ## reached end of file, so close it
    f.close()
    
    ## sort chiSquared array and trim it if it is over 20 elements long
    ascendingOrder = np.argsort(chiSquareds)
    if len(ascendingOrder)>20:
        ascendingOrder = ascendingOrder[0:20]
    
#    # get top 20 chiSquared values to match later
#    chiSquaredGooders = []
#    for ind in ascendingOrder:
#        chiSquaredGooders.append(chiSquareds[ind])
    
    #kill off extra arrays to save memory (NOT SURE IF THIS HELPS THOUGH!!)
    del chiSquareds
       
    # de-sort order array for file hunting and then re-sort after
    srt = np.argsort(ascendingOrder)
    unSortedAscendingOrder = []
    for ind in srt:
        unSortedAscendingOrder.append(ascendingOrder[ind])
               
    ## go back through file and pull out orbits that are in top 20
    bestOrbits = []
    saveLineChiSquareds = []
    f = open(filename,'r')
    # strip off the .txt part to make the plot version of the filename
    plotFileTitle = f.readline()[:-5]
    headings = f.readline()
    line = 'asdf'
    lineNumber = 0
    indicesNumber = 0
    while line!='':
        line = f.readline()
        if line!='':
            #Check if this line is on the list to save
            if lineNumber==unSortedAscendingOrder[indicesNumber]:
                dataLineCols = line.split()
                curChiSquared = float(dataLineCols[7])
                saveLine = '['+dataLineCols[0]+', '+dataLineCols[1]+', '+dataLineCols[2]+', '+dataLineCols[3]+\
                ', '+dataLineCols[4]+', '+dataLineCols[5]+', '+dataLineCols[6]+', '+ dataLineCols[8]+', '+dataLineCols[9] +']'+'    chiSquared = '+dataLineCols[7]
                bestOrbits.append(saveLine)
                saveLineChiSquareds.append(float(dataLineCols[7]))
                if indicesNumber<19:
                    indicesNumber +=1
                    
            # increment line number
            lineNumber +=1
        
    ## reached end of file, so close it
    f.close()   
    
    ## sort bestOrbit strings array to proper ascending order
    order = np.argsort(saveLineChiSquareds)
    bestOrbitsSorted = []
    for ind in order:
        bestOrbitsSorted.append(bestOrbits[ind])
        
    ## print final best 20 results, which would be at end of list
    print '\nHere you go! ;-)  :\n'
    for i in range(0,20):
        try:
            print bestOrbitsSorted[i]
        except:
            print 'OOPS, only found the top '+str(i)+', sorry about that. \nYou can send Kyle an angry email if you like :-P'
            break
        
def gaussianDist(x, mean, variance): 
    """
    Just a simple function to return the probability value of a value in a gaussian distribution.
    """
    varSquared = np.power(variance,2)
    
    exponent = -(np.power((x-mean),2)/(2*varSquared))
    #front = 1/np.sqrt(2*pi*varSquared) ## multiply the prob by this to give you the normal distribution
    
    prob = np.exp(exponent)
    
    return prob

def histConverter(chiSquareds, data, plot, xlabel, confLevels=True, weight=True, normed=True):
    """
    This function is for creating a modified histogram plot to be properly normalized to a max value
    of 1.0.  There is also the option to color the histogram bars according to the 95 and 68% confidence
    levels of the data.
    """
    # if requested, calculate the likelihoods which will be the weights for the data points
    if weight:
        theWeights = likelihoodsCalc(chiSquareds)
    else:
        theWeights = np.ones(len(data))
    # make initial histogram and update it below
    (n,bins,rectangles)=plot.hist(data, bins=50, normed=normed, weights=theWeights)#, fill=False)
    if confLevels:
        # get Confidence Levels [[68.3% minimum, 68.3% maximum],[95.5% minimum, 95.5% maximum]]
        CLevels = ConfLevelFunc(chiSquareds,data)
        # Update rectangel patche's colors according to Confidence Levels of the data
        recs2 = []
        for rec in rectangles:
                x=rec.get_x()
                c = 'w'
                if(x>CLevels[1][0])and(x<CLevels[1][1]):
                    c = '0.5'
                if (x>CLevels[0][0])and(x<CLevels[0][1]):
                    c = '0.3'
                recs2.append(mpl.patches.Rectangle(xy=rec.get_xy(), width=rec.get_width(),height=rec.get_height(), color=c))
        # draw updated patches on plot
        for rec in recs2:
                plot.add_patch(rec)
    
    if normed:
        #update the y limit and its ticks and tick labels
        plot.axes.set_ylim([0.0,plot.axes.get_ylim()[1]*1.2])
        locs = []
        lim = plot.axes.get_ylim()[1]
        for i in range(0,6):
            locs.append((lim/6.0)*i)
        plot.axes.set_yticks(locs)
        plot.axes.set_yticklabels(['0.0','0.2','0.4','0.6','0.8','1.0'])
    # draw a line up the median of the data
    Median = np.median(data)
    plot.plot([Median, Median],plot.axes.get_ylim(),'k')
    # add x label
    plot.axes.set_xlabel(xlabel)
    
    return plot

def likelihoodsCalc(chiSquareds):
    """
    This function will convert the chiSqureds into likelihoods following 
    L =e^(-chiSquared/2)
    """
    
    likelihoods = []
    for chiSquared in chiSquareds:
        likelihood = math.exp(-chiSquared/2.0)
        likelihoods.append(likelihood)
    
    return likelihoods
    
def multiEpochOrbCalc(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors,\
                       epochs, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total=False,\
                       Mass1=1, Mass2=1, verbose=False):
    """
    
    Note: SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors inputs MUST all have
          the same length.
    """
    
    # init output lists
    ns = []
    Ms = []
    Es = []
    thetas = []
    Sep_Dists = []
    SA_arcsec_measured_models = []
    PA_deg_measured_models = []
    a1s = []
    a2s = []

    # initial values for boolean equality of while statement 
    i = 0
    chi_squared_total = 0.0
    
    while i<len(epochs):    

        SA_arcsec_measured_REAL = SA_arcsec_measured_REALs[i]
        SA_mean_error = SA_mean_errors[i]
        PA_deg_measured_REAL = PA_deg_measured_REALs[i]
        PA_mean_error = PA_mean_errors[i]
        t = epochs[i]  
        
        # call orbitCalculator to take random variables and calc orbital elements
#        (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU, SA_arcsec_measured_model, PA_deg_measured_model, a1, a2) = \
#        orbitCalculator(t, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total,\
#                        Mass1=Mass1, Mass2=Mass2, verbose=verbose)
        (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU, SA_arcsec_measured_model, PA_deg_measured_model, a1, a2) = \
        orbitCalculator2(t, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total=a_total,\
                        Mass1=Mass1, Mass2=Mass2, verbose=verbose)
        # store output orbital elements in lists for later plotting
        ns.append(n)
        Ms.append(M_deg)
        Es.append(E_latest_deg)
        thetas.append(TA_deg)
        Sep_Dists.append(Sep_Dist_AU)
        SA_arcsec_measured_models.append(SA_arcsec_measured_model)
        PA_deg_measured_models.append(PA_deg_measured_model)
        a1s.append(a1)
        a2s.append(a2)
            
        ## calc both PA and SA kai's 
        SA_chi_squared = chiSquaredCalc(SA_arcsec_measured_REAL, SA_mean_error, SA_arcsec_measured_model)
        # we must account for the 360deg boundry, using +- 25deg from it as the region of conflict
        if ((PA_deg_measured_model-25.0)<0.0) and (PA_deg_measured_REAL>335.0):
            #ie real is close to 360 and model is just over 360
            PA_deg_measured_model = PA_deg_measured_model +360.0
        if ((PA_deg_measured_REAL-25.0)<0.0) and (PA_deg_measured_model>335.0):
            #ie model is close to 360 and real is just over 360
            PA_deg_measured_REAL = PA_deg_measured_REAL +360.0
        PA_chi_squared = chiSquaredCalc(PA_deg_measured_REAL, PA_mean_error, PA_deg_measured_model)
        # Add them to get the updated total
        chi_squared_total = chi_squared_total + PA_chi_squared + SA_chi_squared
            
        # increment to next epoch to check these input for 
        i = i + 1   
        
        ##### while loop ends here!
        
    return (chi_squared_total, ns, Ms, Es, thetas, Sep_Dists, SA_arcsec_measured_models, PA_deg_measured_models, a1s, a2s) 

def multiEpochOrbCalc3(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors,\
                       epochs, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total=False,\
                       Mass1=1, Mass2=1, verbose=False):
    """
    
    Note: SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors inputs MUST all have
          the same length.
    """
    testing = False
    # init output lists
    ns = []
    Ms = []
    Es = []
    thetas = []
    xs = []
    ys = []
    a1s = []
    a2s = []

    # initial values for boolean equality of while statement 
    i = 0
    chi_squared_total = 0.0
    
    while i<len(epochs):    

        SA_arcsec_measured_REAL = SA_arcsec_measured_REALs[i]
        SA_mean_error = SA_mean_errors[i]
        PA_deg_measured_REAL = PA_deg_measured_REALs[i]
        PA_mean_error = PA_mean_errors[i]
        t = epochs[i]  
        
        # call orbitCalculator to take random variables and calc orbital elements
        (n, M_deg, E_latest_deg, TA_deg, x, y, a1, a2) = \
        orbitCalculator3(t, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total=a_total,\
                        Mass1=Mass1, Mass2=Mass2, verbose=verbose)
        # store output orbital elements in lists for later plotting
        ns.append(n)
        Ms.append(M_deg)
        Es.append(E_latest_deg)
        thetas.append(TA_deg)
        xs.append(x)
        ys.append(y)
        a1s.append(a1)
        a2s.append(a2)
            
        x_real = SA_arcsec_measured_REAL*math.cos(math.radians(PA_deg_measured_REAL))  
        y_real = SA_arcsec_measured_REAL*math.sin(math.radians(PA_deg_measured_REAL))  
        if testing:
            print '************************************************************************************************'
            print 'SA_arcsec_measured_REAL = ',SA_arcsec_measured_REAL
            print 'PA_deg_measured_REAL = ',PA_deg_measured_REAL
            print '\nx_real =SA_arcsec_measured_REAL*math.cos(math.radians(PA_deg_measured_REAL)) = ',x_real
            print 'y_real =SA_arcsec_measured_REAL*math.sin(math.radians(PA_deg_measured_REAL)) = ',y_real
            print '################################################################################################'
        
        # calc error in x,y data
        x_real_error = x_real*((SA_mean_error/SA_arcsec_measured_REAL)+abs(math.radians(PA_mean_error)*math.tan(PA_deg_measured_REAL)))
        y_real_error = y_real*((SA_mean_error/SA_arcsec_measured_REAL)+abs(math.radians(PA_mean_error)/math.tan(PA_deg_measured_REAL)))
        
        ## calc both PA and SA kai's 
        x_chi_squared = chiSquaredCalc(x_real, x_real_error, x)
        y_chi_squared = chiSquaredCalc(y_real, y_real_error, y)
        # Add them to get the updated total
        chi_squared_total = chi_squared_total + y_chi_squared + x_chi_squared
            
        # increment to next epoch to check these input for 
        i = i + 1   
        
        ##### while loop ends here!
        
    return (chi_squared_total, ns, Ms, Es, thetas, xs, ys, a1s, a2s)

def TAcalculator2(t,e, T, period, verbose=False, debug=False):
    """
    This is the same as TAcalculator but with some updates to the Newton's method loop
    solving the Kepler's equation found in Double Stars by Heintz.
    $$$$$$$$$$$$ Later testing of the two versions will prove which is better. $$$$$$$$$$$$$
    """
    
    ## calculate the Mean Motion
    n = (2*pi)/period
    if verbose:
        print '#'*50
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
    # the inital guess used here is from Argyle "Observing and Measuring Visual Double Stars"
    try:
        E_latest = M+e*math.sin(M)+((e**2.0)/(2.0*M))*math.sin(2.0*M)
    except:
        # not sure why, but there were cases of division by zero
        # thus, for these cases, I will resort to my old predicted start value that doesn't suffer from this.
        E_latest = M+e*math.sin(M)
    M_last = M
    # show input value to 
    if debug:
        print "Inputs to Newton's Method are : \nM [rad]= "+str(M)+"\nE_last [rad]= "+\
        str(E_last)+"\nE_latest [rad] = "+str(E_latest)+"\ne = "+str(e)
        
        print "\nStarting to run Newton's while loop."
    
    count = 0 # a counter to stop inf loops in Newton's method below
    while (abs(E_last-E_latest) > (1.0e-10))and(count<100):
        if debug:
            print 'current E [rad]= ', E_latest
        E_last = E_latest
        M_last = E_last - e*math.sin(E_last)
        E_latest = E_last + ((M-M_last)/(1.0-e*math.cos(E_last)))
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
        if debug:
            print "This resultant E solves the original equation, Newton's Method worked :-)"
            print '-'*50
    ### Newton's loop finished! ###
    
    ## calculate True Anomaly from Eccentric Anomaly
    try:
        top = (math.cos(E_latest)-e)
        btm = (1.0-e*math.cos(E_latest))
        TA_rad  = math.acos(top/btm) 
        if E_latest<0.0:
            #print 'top is negative'
            TA_rad = -1.0*TA_rad
    except:
        print "\nTA_rad calc failed!"
        print "E_latest = ",E_latest
        print "e = ",e
    
    if E_latest>(2.0*pi):
        # convert E to be inside one orbit (ie. under 2*PI)
        numCirclesD = E_latest/(2.0*pi)
        numCirclesI = int(numCirclesD)
        E_latest_oneCircle  = E_latest-numCirclesI*2.0*pi
        if verbose: 
            print "E_latest found to be "+str(numCirclesI)+" times over 2pi, so made E_latest_circle = "+str(math.degrees(E_latest_oneCircle))
    else:
        E_latest_oneCircle  = E_latest
    if (E_latest_oneCircle>pi):
        # this is to take care of the silly fact that after E=180deg,
        # the equation above produces TAs that go down from 180, rather than up.
        ## This has been found to be due to the equation eliminating negative
        ## values because of the cosine, thus the new update of E_latest<0.0:
        ## should fix this problem, but the solutions will not go 0->360
        ## but rather 0->180, -180->-0.
        TA_rad_orig = TA_rad
        TA_rad = 2.0*pi - TA_rad
        if verbose:
            print "E_latest found to be over PI, so changed TA from "+str(math.degrees(TA_rad_orig))+" to "+str(math.degrees(TA_rad))
        
#    ## Calculate TA in another way    
#    x = ((1.0-e**2.0)**(1.0/2.0))*math.cos(E_latest/2.0)
#    y = ((1.0+e**2.0)**(1.0/2.0))*math.sin(E_latest/2.0)
#    TA_rad2 = 2.0*math.atan2(y, x)
#    #print 'TA_2 = '+str(math.degrees(TA_rad2))
#    
#    print 'TA = ',math.degrees(TA_rad)
#    print 'TA2 = ',math.degrees(TA_rad2)
    
    
    return (n, M_deg, E_latest_deg,TA_rad)

def timeString(duration):
    """
    takes a time duration in seconds and returns it in a nice string with info on number of 
    hours, minutes and/or seconds depending on duration.
    Input: float
    Output: string
    """
    if duration<60:
        totalTimeString = str(int(duration))+' seconds'
    elif duration>3600:
        duration = duration/3600.0
        totalTimeString = str(int(duration))+' hours and '+str(int(60*((duration)-int(duration))))+' minutes'
    else:
        duration = duration/60.0
        totalTimeString = str(int(duration))+' minutes and '+str(int(60*((duration)-int(duration))))+' seconds'
        
    return totalTimeString


def rvResidualWithoutPlanetResidual():
    """
    A totally temporary function to calculate the residual velocity of the primary star
    WITHOUT the residual from the planet.  Thus only having the velocity that could be 
    caused by the secondary star/companion.
    """

    from paramSettingsDict import paramSettingsDict
    RV_epochs = paramSettingsDict['Data']['RV_epochs']
    #RVerrors = paramSettingsDict['Data']['RVerrors']
    RVs = paramSettingsDict['Data']['RVs']
    M1 = paramSettingsDict['Data']['M1']
    # data for planet's orbit that are needed to calculate proper chiSquare for RV fit
    # all below values for the planet orbit and RV data are from Butler2006.
    K_p = 461.1 #[m/s]
    p_p = 3.31246   #[days]
    #planet_p_years = planet_p/365.25
    e_p = 0.023    
    argPeri_p = 188.0   #[deg]
    T_p = 2446957.8   #[JD]
    sigma_jitter = 15.0    #[m/s]
    
    rvWithoutPlanetResiduals = []
    for RVdataSet in range(0,len(RVs)):
        print '\nworking on RVdataSet ',RVdataSet
        for epoch in range(0,len(RV_epochs[RVdataSet])):
            v_r_p = vrCalculatorPlanet(RV_epochs[RVdataSet][epoch], e_p, T_p, p_p, argPeri_p, M1, M2SineI=False, K=K_p, verbose=False)
            RVprimary = RVs[RVdataSet][epoch]
            rvWithoutPlanetResidual = RVprimary - v_r_p
            rvWithoutPlanetResiduals.append(rvWithoutPlanetResidual)
            print v_r_p
 
def rv1bodyCalculator(RV_epochs, RVs, RVerrors, sigma_jitter, i, p, e, T, argPeri, a, verbose=False):
    """
    This is for calculating the RV and resulting chiSquared due to a planet ONLY, thus no companion star.
    Naturally this system would have a planet with mass<< primary star mass.
    
    
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
    
    if verbose:
        print '\n## Using the semi-major to calculate K for the primary star due to companion star ##'
    
    for epoch in range(0,len(RVs)):
        if verbose:
            print '\n******************* RV values for epoch '+str(epoch+1)+' *****************'
        # calculate the velocity residual due to the companion 
        (v_r_c, K_c) = vrCalculatorStar2(RV_epochs[epoch],e, T, p, argPeri, a, i=i, K=False, verbose=False)
        
        # calculate chiSquared for this epoch
        chiSquaredCurr = chiSquaredCalc(RVs[epoch], (RVerrors[epoch]+sigma_jitter), (v_r_c))
        
        chiSquaredTotal = chiSquaredTotal + chiSquaredCurr
        
        if verbose:
            print 'RV = (calculated)-measured = ('+str(v_r_c)+') - '+str(RVs[epoch])+' = '+str((v_r_c)-RVs[epoch]) 
            print '(RVerrors[epoch]+sigma_jitter) = '+str((RVerrors[epoch]+sigma_jitter))
            print 'chiSquaredCurr = ',chiSquaredCurr
    
    if verbose:
            print '\nK of companion = ',K_c
            
    return chiSquaredTotal

def rv1bodyCalculator2(RV_epochs, RVs, RVerrors, sigma_jitter, p, e, T, argPeri, M1,M2SineI, verbose=False):
    """
    This is for calculating the RV and resulting chiSquared due to a planet ONLY, thus no companion star.
    Naturally this system would have a planet with mass<< primary star mass.
    
    
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
    
    if verbose:
        print '\n## Using the semi-major to calculate K for the primary star due to companion star ##'
    
    for epoch in range(0,len(RVs)):
        if verbose:
            print '\n*** RV values for epoch '+str(epoch+1)+' ***'
        # calculate the velocity residual due to the companion 
        (v_r_c, K_c) = vrCalculatorPlanet(RV_epochs[epoch],e, T, p, argPeri,M1,M2SineI=M2SineI, K=False, verbose=False)
        # calculate chiSquared for this epoch
        chiSquaredCurr = chiSquaredCalc(RVs[epoch], (RVerrors[epoch]+sigma_jitter), (v_r_c))
        
        chiSquaredTotal = chiSquaredTotal + chiSquaredCurr
        
        if verbose:
            print 'RV = (calculated)-measured = ('+str(v_r_c)+') - '+str(RVs[epoch])+' = '+str((v_r_c)-RVs[epoch]) 
            print '(RVerrors[epoch]+sigma_jitter) = '+str((RVerrors[epoch]+sigma_jitter))
            print 'chiSquaredCurr = ',chiSquaredCurr
    
    if verbose:
            print '\nK of companion = ',K_c
            
    return chiSquaredTotal 
        
def rv2bodyCalculator3(RV_epochs, RVs, RVerrors, sigma_jitter, M1, M2_c, i_c, p_c, e_c, T_c, argPeri_c, \
                                                                    K_p, p_p, e_p, argPeri_p, T_p, verbose=False):
    """
    NOTE: This version uses the K equation with each objects mass.  
        Use rv2bodyCalculator4 to use the semi-major axis instead.
    
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
    
    if verbose:
        print '\n## Using the masses to calculate K for the primary star due to companion star ##'
    for epoch in range(0,len(RVs)):
        if verbose:
            print '\n*** RV values for epoch '+str(epoch+1)+' ***'
        # calculate the velocity residual due to the companion star
        (v_r_c, K_c) = vrCalculatorStar(RV_epochs[epoch],e_c, T_c, p_c, argPeri_c, M1, M2=M2_c, i=i_c, K=False, verbose=True)
        # calculate the velocity residual due to the planet around primary
        (v_r_p, K_p) = vrCalculatorPlanet(RV_epochs[epoch], e_p, T_p, p_p, argPeri_p, M1, M2SineI=False, K=K_p, verbose=True)

        # calculate chiSquared for this epoch
        chiSquaredCurr = chiSquaredCalc(RVs[epoch], (RVerrors[epoch]+sigma_jitter), (v_r_c+v_r_p))
        
        chiSquaredTotal = chiSquaredTotal + chiSquaredCurr
        
        if verbose:
            print 'RV = ('+str(v_r_c)+' + '+str(v_r_p)+') - '+str(RVs[epoch])+' = '+str((v_r_c+v_r_p)-RVs[epoch]) 
            print '(RVerrors[epoch]+sigma_jitter) = '+str((RVerrors[epoch]+sigma_jitter))
            print 'chiSquaredCurr = ',chiSquaredCurr
            
    return chiSquaredTotal

def rv2bodyCalculator4(RV_epochs, RVs, RVerrors, sigma_jitter, M1, a1, i_c, p_c, e_c, T_c, argPeri_c, \
                                                                    K_p, p_p, e_p, argPeri_p, T_p, verbose=False):
    """
    NOTE: This version uses the K equation with the semi-major axis instead of each objects masses.
          Use rv2bodyCalculator3 to use the objects masses instead.
    
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
    
    if verbose:
        print '\n## Using the semi-major to calculate K for the primary star due to companion star ##'
    
    for epoch in range(0,len(RVs)):
        if verbose:
            print '\n*** RV values for epoch '+str(epoch+1)+' ***'
        # calculate the velocity residual due to the companion star
        (v_r_c, K_c) = vrCalculatorStar2(RV_epochs[epoch],e_c, T_c, p_c, argPeri_c, a1, i=i_c, K=False, verbose=False)
        # calculate the velocity residual due to the planet around primary
        (v_r_p, K_p) = vrCalculatorPlanet(RV_epochs[epoch], e_p, T_p, p_p, argPeri_p, M1, M2SineI=False, K=K_p, verbose=False)
        
        # calculate chiSquared for this epoch
        chiSquaredCurr = chiSquaredCalc(RVs[epoch], (RVerrors[epoch]+sigma_jitter), (v_r_c+v_r_p))
        
        chiSquaredTotal = chiSquaredTotal + chiSquaredCurr
        
        if verbose:
            print 'RV = ('+str(v_r_c)+' + '+str(v_r_p)+') - '+str(RVs[epoch])+' = '+str((v_r_c+v_r_p)-RVs[epoch]) 
            print '(RVerrors[epoch]+sigma_jitter) = '+str((RVerrors[epoch]+sigma_jitter))
            print 'chiSquaredCurr = ',chiSquaredCurr
            
    return chiSquaredTotal

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

def starAndErrorPolys(SAs,SAerrors,PAs,PAerrors,asConversion, transData):
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
        yCent = SAs[i]*math.cos(math.radians(PAs[i]))*asConversion
        xCorner = xCent-0.5*w
        yCorner = -(yCent+0.5*h)
        
        rect = patches.Rectangle((xCorner,yCorner),width=w,height=h,color='red',alpha=0.2)
        t = pylab.matplotlib.transforms.Affine2D().rotate_deg_around(xCent,-yCent,PAs[i]) +transData
        rect.set_transform(t)
        errorBoxes.append(rect)
        
        # determin x and y locations of the observed PA and SA's for companion star/planet
        # then make a star polygon for each, same as for M1 but much smaller
        m2starPolygons.append(star((asConversion/1000.0)*7.0*SAs[0], xCent, -yCent, color='red', N=5, thin = 0.8))
        
    return (errorBoxes, m2starPolygons)
    
def summaryPlotter(plotFileTitle, chiSquareds, inclination_degsAlls, longAN_degsAlls, argPeri_degsAlls, esAlls,\
            periodsAlls, TsAlls, a1Means, a2Means, confLevels=True, weight=True, normed=True, showPlots=False,\
                                                                                     save=False, verbose=False):
    """
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
        

def vrCalculatorPlanet(t,e,T,period,argPeri,M1,M2SineI=False, K=False, verbose=False):
    """
    Version for when the companion is a planet.  Thus it calculates the residule velocity 
    assuming this and the consequence that (M1+M2~=M1).  Returns residule velocity 
    of the primary star!! due to the planet, not the the vel of the planet.
    
    M1 in Msun
    M2SineI in Mjupiter
    argPeri in degrees
    a in AU
    T in JD
    t in JD
    e unitless
    period in days
    
    K=False implies that the provided values will be used to 
    calculate it from scratch, else the provided value will
    be used to calculate the radial velocity residule.  Thus,
    this param must be either False or a float in units of 
    m/s.
    
    returned vel will be in m/s
    """
        
    if (K==False) and (M2SineI==False):
        print 'vrCalcPlanet1225: The value of M2SineI and K cannot both be False, one MUST be defined'
    if (K==False) and (M2SineI!=False):
        # convert units of days to seconds
        period_seconds = period*86400.0#31557600.0
        M1_kg = M1*1.98892e30
        M2SineI_kg = M2SineI*1.8986e27
        G = 6.67300e-11
        
        # Calc K in parts to see equation better
        A = ((2.0*pi*G)/period_seconds)**(1.0/3.0)
        B = M2SineI_kg/(M1_kg**(2.0/3.0))
        C = 1.0/math.sqrt(1.0-e**2.0)
        # put it all together
        K = A*B*C
    
    if verbose:
        print 'K_planet = ',K
        
    period_years = period/365.25
    
    ## get the TAcalculator to find the TA in radians
    (n, M_deg, E_latest_deg,TA_rad) = TAcalculator2(t,e, T, period_years, verbose=False)
    
    if verbose:
        print "Planet:"
        print "TA_deg planet = "+str(math.degrees(TA_rad))
        print "T = "+str(T)+", t = "+str(t)+", period_years = "+str(period_years)
    
    v_r = K*(math.cos(math.radians(argPeri)+TA_rad)+e*math.cos(math.radians(argPeri)))

    return (v_r, K)

def vrCalculatorStar(t,e,T,period,argPeri,M1,M2,i=False, K=False, verbose=False):
    """
    NOTE: This version is the one that uses the objects masses to calculate K.
          To use the one that calculates it using the semi-major axis, use vrCalculatorStar2.
    
    M1 and M2 in Msun
    argPeri in degrees
    i in degrees
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
        print 'vrCalcStar1268: The value of i and K cannot both be False, one MUST be defined'
    if (K==False) and (i!=False):
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

    if verbose:
        print 'K_Star = ',K
        
    ## get the TAcalculator to find the TA in radians
    (n, M_deg, E_latest_deg,TA_rad) = TAcalculator2(t,e, T, period, verbose=False)
    
    if verbose:
        print "Star:"
        print "TA_deg star = "+str(math.degrees(TA_rad))
        print "T = "+str(T)+", t = "+str(t)+", period_years = "+str(period)
    
    v_r = K*(math.cos(math.radians(argPeri)+TA_rad)+e*math.cos(math.radians(argPeri)))

    return (v_r, K)

def vrCalculatorStar2(t,e,T,period,argPeri,a1,i=False, K=False, verbose=False):
    """
    NOTE: this is the version which uses the K equation with the semi-major axis of primary's orbit.
          To use the masses instead, use vrCalculatorStar.
    
    argPeri in degrees
    i in degrees
    T in JD
    t in JD
    e unitless
    period in years
    a1 in AU
    
    K=False implies that the provided values will be used to 
    calculate it from scratch, else the provided value will
    be used to calculate the radial velocity residule.  Thus,
    this param must be either False or a float in units of 
    m/s.
    
    returned vel will be in m/s
    """
        
    if (K==False) and (i==False):
        print 'vrCalcStar1268: The value of i and K cannot both be False, one MUST be defined'
    if (K==False) and (i!=False):
        # convert units of years to seconds
        seconds_per_yr = 31557600.0
        period_seconds = period*seconds_per_yr
        #M1_kg = M1*1.98892e30
        #M2_kg = M2*1.98892e30
        meters_per_AU = 149598000000.0
        a1_meters = a1*meters_per_AU
        
#        # Calc K in parts to see equation better
#        A = ((2.0*pi*6.67e-11*(M1_kg+M2_kg))/period_seconds)**(1.0/3.0)
#        B = M2_kg/M1_kg
#        C = math.sin(math.radians(i))/math.sqrt(1.0-e**2)
#        # put it all together
#        K = A*B*C
        
        # Calc K in parts to see equation better
        A = (2.0*pi)/period_seconds
        B = a1_meters
        C = math.sin(math.radians(i))/math.sqrt(1.0-e**2.0)
        # put it all together
        K = A*B*C
        
    if False:
        print '\nperiod = ',period
        print 'a1 = ',a1
        print 'e = ',e
        print 'K_Star = ',K
        
    ## get the TAcalculator to find the TA in radians
    (n, M_deg, E_latest_deg,TA_rad) = TAcalculator2(t,e, T, period, verbose=False, debug=False)
    if verbose:
        #print '#######################################################'
        print 'TA = '+str(math.degrees(TA_rad))  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        if False:
            print 'math.degrees(math.radians(argPeri-pi)+TA_rad) = '+str(math.degrees(math.radians(argPeri-pi)+TA_rad))
            print 'math.cos(math.radians(argPeri-pi)+TA_rad) = '+str(math.cos(math.radians(argPeri-pi)+TA_rad))
            print 'e = '+str(e)
            print 'argPeri = '+str(argPeri)
            print 'e*math.cos(math.radians(argPeri-pi)) = '+str(e*math.cos(math.radians(argPeri-pi)))
            print '#######################################################'
    v_r = K*(math.cos(math.radians(argPeri)+TA_rad)+e*math.cos(math.radians(argPeri)))

    return (v_r, K)

def orbitCalculator(t, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total,\
                    Mass1=1, Mass2=1, verbose=False):
    """
    This function will calculate the output orbital parameters from the inputs for a SINGLE epoch.  Use the 
    multiEpochOrbCalc to do this for multiple epochs at the same time.  
    The t and Sys_Dist_PC parameters are generally known/measured values,
    while the rest are typically random numbers.
    This function is designed to work as the calculator inside a Monte Carlo Simulator.
    All the equations are the result of Kepler's Laws and naturally occurring symmetries of binary star systems.
    
    Note: RP stands for Reference Plane, ie the plane that values are measured in from Earth.
          OP stands for Orbital Plane, ie in the plane that the stars are orbiting in.
    
    Inputs:
    @param t:                            = epoch of observation/image [julian date]
    @type t:                             = float
    @param Sys_Dist_PC:                  = measured system distance from Earth [PC]
    @type Sys_Dist_PC:                   = float
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
    @param a_total:                        = Total semi-major axis of system [AU]
    @type a_total:                         = float
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
    @param Sep_Dist_AU_OP:               = Separation Distance in orbital plane [AU]
    @type Sep_Dist_AU_OP:                = float
    @param SA_arcsec_RP_model:           = measured Position Angle in image [deg]
    @type SA_arcsec_RP_model:            = float
    @param PA_deg_RP_model:              = measured Position Angle in image [deg]
    @type PA_deg_RP_model:               = float
    @param a1:                           = semi-major axis of M1 [AU]
    @type a1:                            = float
    @param a2:                           = semi-major axis of M2 [AU]
    @type a2:                            = float
    """
    if verbose:
        print '\n'+'*'*50
        print 'Starting to calculate Orbital Parameters\n'
        print 'Input variable values: '+'\nCurrent epoch time [JD] = '+str(t)+\
        '\nSystem Distance from Earth [PC] = '+str(Sys_Dist_PC)+\
        '\nInclination [deg] = '+ str(inclination_deg)+\
        '\nLongitude of Ascending Node [deg] = '+str(longAN_deg)+\
        '\nEccentricity = '+str(e)+'\nLast Periapsis Epoch [JD] = '+\
        str(T)+'\nPeriod [yrs] = '+str(period)+\
        '\nArgument of Periapsis in orbital plane [deg] = '+\
        str(argPeri_deg)+'\nverbose = ',str(verbose)
    #----------------------------------------------------------------------

#    ## calculate the Mean Motion
#    n = (2*pi)/period
#    if verbose:
#        print 'Mean Motion [rad/yr]= '+str(n)
#    
#    ## calculate Mean Anomaly
#    M = n*((t-T)/365.25)
#    M_deg = math.degrees(M) # convert resulting M to degrees
#    if verbose:
#        print 'Mean Anomaly [deg]= ',M_deg
#
#    ### Performing Newton's Method to get the Eccentric Anomaly value ###
#    if verbose:
#        print '-'*50
#        
#    # initial guess (E_last), will be updated in loop.  
#    # Anything works, just takes longer if further from real value. => pi
#    E_last = 2*pi
#    # stored initial value to be updated in loop
#    # this value is always very close to the true value and will minimize the number of loops
#    E_latest = M+e*math.sin(M) 
#    # show input value to 
#    if verbose:
#        print "Inputs to Newton's Method are : \nM [rad]= "+str(M)+"\nE_last [rad]= "+\
#        str(E_last)+"\nE_latest [rad] = "+str(E_latest)+"\ne = "+str(e)
#        
#        print "\nStarting to run Newton's while loop."
#    
#    count = 0 # a counter to stop inf loops in Newton's method below
#    while (abs(E_last-E_latest) > (1.0e-10))and(count<100):
#        if verbose:
#            print 'current E [rad]= ', E_latest
#        E_last = E_latest
#        #E_latest = E_last - ((M-E_last+e*math.sin(E_last))/(e*math.cos(E_last)-1.0))\
#        E_latest = E_last - ((E_last-M-e*math.sin(E_last))/(-e*math.cos(E_last)+1.0))
#        count = count+1
#
#    E_latest_deg = math.degrees(E_latest) # convert resulting E to degrees
#    if verbose:
#        print "The resultant E value is [deg] = ", E_latest_deg
#    # check if the resultant value solves the original equation.
#    Mnewton = math.degrees(E_latest-e*math.sin(E_latest))
#    if abs(M_deg-Mnewton)>(1.0e-5):
#        if verbose:
#            print "PROBLEM: This resultant E does not satisfy the original equation, Newton's method FAILED !!!"
#            print 'M from this E Equals = '+str(Mnewton)
#            print 'M original = '+str(M_deg)
#            print 'E initial = '+str(math.degrees(M+e*math.sin(M) ))
#            print 'e = '+str(e)
#    else:
#        if verbose:
#            print "This resultant E solves the original equation, Newton's Method worked :-)"
#            print '-'*50
#    ### Newton's loop finished! ###
#    
#    ## calculate True Anomaly from Eccentric Anomaly
#    try:
#        TA_rad  = math.acos((math.cos(E_latest)-e)/(1.0-e*math.cos(E_latest))) 
#    except:
#        print "\nTA_rad calc failed!"
#        print "E_latest = ",E_latest
#        print "e = ",e
#    if E_latest>pi:
#        # this is to take care of the silly fact that after E=180deg,
#        # the equation above produces TAs that go down from 180, rather than up.
#        TA_rad = 2.0*pi - TA_rad
    
    ## get the TAcalculator to find the TA in radians
    (n, M_deg, E_latest_deg,TA_rad) = TAcalculator2(t,e, T, period, verbose=verbose)
    TA_deg = math.degrees(TA_rad)
    if verbose:
        print 'True Anomaly [deg] = ', TA_deg
        
    ## calculate the Separation Distance in orbital plane
    Sep_Dist_AU_OP = a_total*((1.0-e*e))/(1.0+e*math.cos(TA_rad))
    if verbose:
        print 'Separation Distance in orbital plane [AU] = ',Sep_Dist_AU_OP
    
    ## calculate measured Separation Distance in reference plane
    inclination_rad = math.radians(inclination_deg) # convert inclination to radians 
    # angle between ascending node and M2 in orbital plane
    ang_LN_to_M2_deg = argPeri_deg + TA_deg
    ang_LN_to_M2_rad = math.radians(ang_LN_to_M2_deg)
    if verbose:
        print "inclination_deg = ",inclination_deg
        print "inclination_rad = ",inclination_rad
        print "ang_LN_to_M2_rad = ",ang_LN_to_M2_rad
    # find horizontal and vertical WRT Line of Nodes (X and Y) components of Sep_Dist corrected for inclination
    # Note: only the 'vertical (Y)' component will be effected by inclination
    # Note2: on this axis, Z would be the direction of 'line of sight' of the observer.
    Sep_Dist_AU_RP_y = Sep_Dist_AU_OP*math.sin(ang_LN_to_M2_rad)*math.cos(inclination_rad)
    Sep_Dist_AU_RP_x = Sep_Dist_AU_OP*math.cos(ang_LN_to_M2_rad)
    # use Pythagorous theorem to find the final Sep Dist in the Ref Plane
    Sep_Dist_AU_RP = math.sqrt(math.pow(Sep_Dist_AU_RP_x, 2.0)+math.pow(Sep_Dist_AU_RP_y,2.0))
    if verbose:
        print 'Sep_Dist_AU_RP_y = ',Sep_Dist_AU_RP_y
        print 'Sep_Dist_AU_RP_x = ',Sep_Dist_AU_RP_x
        print 'Separation Distance in the reference plane [AU] = ',Sep_Dist_AU_RP

    ## calculate corrected angle between Line of Nodes and M2 (ie. AN-focus-M2)
    # the first atan will produce an angle between the X-axis, rather than only CC from positive X-axis.
    # must take all 4 quadrants into account.  
    ang_LN_to_M2_corrPOSNEG = math.degrees(math.atan(Sep_Dist_AU_RP_y/Sep_Dist_AU_RP_x))
    
    if (Sep_Dist_AU_RP_x>=0.0)and (Sep_Dist_AU_RP_y>=0.0):
        # quadrant 1
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG 
    elif (Sep_Dist_AU_RP_x<0.0)and (Sep_Dist_AU_RP_y>=0.0):
        # quadrant 2
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +180.0
    elif (Sep_Dist_AU_RP_x<0.0)and (Sep_Dist_AU_RP_y<0.0):
        # quadrant 3
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +180.0
    elif (Sep_Dist_AU_RP_x>=0.0)and (Sep_Dist_AU_RP_y<0.0):
        # quadrant 4
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +360.0
    else:
        print ' The value of '+str(ang_LN_to_M2_corrPOSNEG)+' is apparently out of ranges'
    
    if verbose:
        print 'uncorrected angle between line of nodes and M2 [deg] = ', ang_LN_to_M2_deg
        print 'Corrected angle between Line of nodes and M2 [deg] = ',ang_LN_to_M2_corr
    
    ## calculate measured Position Angle 
    # it is measured from M1, only if the value is over 360deg do we need to fix it
    totalAngle_RP = longAN_deg + ang_LN_to_M2_corr
    if verbose:
        print "Total angle in the Reference Plane [deg] = ",totalAngle_RP
    if totalAngle_RP>360.0:
        PA_deg_RP_model = totalAngle_RP -360.0
    else:
        PA_deg_RP_model = totalAngle_RP

    if verbose:
        print 'Position Angle in reference plane [deg] = ',PA_deg_RP_model
        
    ## calculate measured Separation Angle (model)
    SA_arcsec_RP_model = Sep_Dist_AU_RP/Sys_Dist_PC # PC->AU and rad->arcsec cancel
    if verbose:
        print 'Separation Angle measured (model) [arcsec] = ',SA_arcsec_RP_model

    if (Mass1 is not 1) and (Mass2 is not 1):
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
        
        if (a1+a2)-a_total>1.0e-5:
            print "obCalc-ln1000: a1 = ",a1
            print "obCalc-ln1001:a2 = ",a2
            print "obCalc-ln1002:a_total = ",a_total

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
   
    return (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SA_arcsec_RP_model, PA_deg_RP_model, a1, a2)

def orbitCalculator2(t, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total=False,\
                    Mass1=1, Mass2=1, verbose=False):
    """
    This function will calculate the output orbital parameters from the inputs for a SINGLE epoch.  Use the 
    multiEpochOrbCalc to do this for multiple epochs at the same time.  
    The t and Sys_Dist_PC parameters are generally known/measured values,
    while the rest are typically random numbers.
    This function is designed to work as the calculator inside a Monte Carlo Simulator.
    All the equations are the result of Kepler's Laws and naturally occurring symmetries of binary star systems.
    
    Note: RP stands for Reference Plane, ie the plane that values are measured in from Earth.
          OP stands for Orbital Plane, ie in the plane that the stars are orbiting in.
    
    Inputs:
    @param t:                            = epoch of observation/image [julian date]
    @type t:                             = float
    @param Sys_Dist_PC:                  = measured system distance from Earth [PC]
    @type Sys_Dist_PC:                   = float
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
    @param a_total:                      = Total semi-major axis of system [AU]
    @type a_total:                       = float if it is to be used, boolean False if to be calculated from Kepler's 3rd
    @param Mass1:                        = Mass of primary in system [kg] (Default=1, means M2 is a planet)
    @type Mass1:                         = float
    @param Mass2:                        = Mass of the secondary in system [kg] (Default=1, means M2 is a planet)
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
    @param Sep_Dist_AU_OP:               = Separation Distance in orbital plane [AU]
    @type Sep_Dist_AU_OP:                = float
    @param SA_arcsec_RP_model:           = measured Position Angle in image [deg]
    @type SA_arcsec_RP_model:            = float
    @param PA_deg_RP_model:              = measured Position Angle in image [deg]
    @type PA_deg_RP_model:               = float
    @param a1:                           = semi-major axis of M1 [AU]
    @type a1:                            = float
    @param a2:                           = semi-major axis of M2 [AU]
    @type a2:                            = float
    """
    if verbose:
        print '\n'+'*'*50
        print 'Starting to calculate Orbital Parameters\n'
        print 'Input variable values: '+'\nCurrent epoch time [JD] = '+str(t)+\
        '\nSystem Distance from Earth [PC] = '+str(Sys_Dist_PC)+\
        '\nInclination [deg] = '+ str(inclination_deg)+\
        '\nLongitude of Ascending Node [deg] = '+str(longAN_deg)+\
        '\nEccentricity = '+str(e)+'\nLast Periapsis Epoch [JD] = '+\
        str(T)+'\nPeriod [yrs] = '+str(period)+\
        '\nArgument of Periapsis in orbital plane [deg] = '+\
        str(argPeri_deg)+'\nverbose = ',str(verbose)
    #----------------------------------------------------------------------
    
    ## check if a_total is to be calculated using Kepler's 3rd
    if a_total==False:
        
        # conversion factors and constants
        SecPerYear = 31557600.0
        G = 6.67300e-11
        MperAU = 149598000000.0
        KGperMsun = 1.98892e30
        
        # convert units of years to seconds
        period_seconds = period*SecPerYear
        
        # convert masses from Msun to kilograms
        Mass1_kg = Mass1*KGperMsun
        Mass2_kg = Mass2*KGperMsun
        # apply K3
        a2 = (((G*(Mass1_kg+Mass2_kg)*(period_seconds**2.0))/(4.0*(pi**2.0))))**(1.0/3.0) # Masses must be in [kg]!!!!! 
        
        # find a1 using mass ratio
        a1 = a2*(Mass2_kg/Mass1_kg)
        
        # sum a's and convert into [AU]
        a_total = (a1+a2)/MperAU
        
    ## get the TAcalculator to find the TA in radians
    (n, M_deg, E_latest_deg,TA_rad) = TAcalculator2(t,e, T, period, verbose=verbose)
    TA_deg = math.degrees(TA_rad)
    if verbose:
        print 'True Anomaly [deg] = ', TA_deg
        
    ## calculate the Separation Distance in orbital plane
    Sep_Dist_AU_OP = (a_total*(1.0-e**2.0))/(1.0+e*math.cos(TA_rad))
    if verbose:
        print 'Separation Distance in orbital plane [AU] = ',Sep_Dist_AU_OP
    
    ## calculate measured Separation Distance in reference plane
    inclination_rad = math.radians(inclination_deg) # convert inclination to radians 
    # angle between ascending node and M2 in orbital plane
    ang_LN_to_M2_deg = argPeri_deg + TA_deg
    if ang_LN_to_M2_deg>360.0:
        ang_LN_to_M2_deg = ang_LN_to_M2_deg-360.0
    ang_LN_to_M2_rad = math.radians(ang_LN_to_M2_deg)
    if verbose:
        print "inclination_deg = ",inclination_deg
        print "inclination_rad = ",inclination_rad
        print "ang_LN_to_M2_rad = ",ang_LN_to_M2_rad
    # find horizontal and vertical WRT Line of Nodes (X and Y) components of Sep_Dist corrected for inclination
    # Note: only the 'vertical (Y)' component will be effected by inclination
    # Note2: on this axis, Z would be the direction of 'line of sight' of the observer.
    Sep_Dist_AU_RP_y = Sep_Dist_AU_OP*math.sin(ang_LN_to_M2_rad)*math.cos(inclination_rad)
    Sep_Dist_AU_RP_x = Sep_Dist_AU_OP*math.cos(ang_LN_to_M2_rad)
    # use Pythagorous theorem to find the final Sep Dist in the Ref Plane
    Sep_Dist_AU_RP = ((Sep_Dist_AU_RP_x**2.0)+(Sep_Dist_AU_RP_y**2.0))**(1.0/2.0)
    if verbose:
        print 'Sep_Dist_AU_RP_y = ',Sep_Dist_AU_RP_y
        print 'Sep_Dist_AU_RP_x = ',Sep_Dist_AU_RP_x
        print 'Separation Distance in the reference plane [AU] = ',Sep_Dist_AU_RP

    ## calculate corrected angle between Line of Nodes and M2 (ie. AN-focus-M2)
    # the first atan will produce an angle between the X-axis, rather than only CC from positive X-axis.
    # must take all 4 quadrants into account.  
    ang_LN_to_M2_corrPOSNEG = math.degrees(math.atan(Sep_Dist_AU_RP_y/Sep_Dist_AU_RP_x))
    
    if (Sep_Dist_AU_RP_x>=0.0)and (Sep_Dist_AU_RP_y>=0.0):
        # quadrant 1
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG 
    elif (Sep_Dist_AU_RP_x<0.0)and (Sep_Dist_AU_RP_y>=0.0):
        # quadrant 2
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +180.0
    elif (Sep_Dist_AU_RP_x<0.0)and (Sep_Dist_AU_RP_y<0.0):
        # quadrant 3
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +180.0
    elif (Sep_Dist_AU_RP_x>=0.0)and (Sep_Dist_AU_RP_y<0.0):
        # quadrant 4
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +360.0
    else:
        print ' The value of '+str(ang_LN_to_M2_corrPOSNEG)+' is apparently out of ranges'
    
    #print '\nuncorrected angle between line of nodes and M2 [deg] = ', ang_LN_to_M2_deg
    #print 'Corrected angle between Line of nodes and M2 [deg] = ',ang_LN_to_M2_corr
    if verbose:
        print 'uncorrected angle between line of nodes and M2 [deg] = ', ang_LN_to_M2_deg
        print 'Corrected angle between Line of nodes and M2 [deg] = ',ang_LN_to_M2_corr
    
    ## calculate measured Position Angle 
    # it is measured from M1, only if the value is over 360deg do we need to fix it
    totalAngle_RP = longAN_deg + ang_LN_to_M2_corr
    if verbose:
        print "Total angle in the Reference Plane [deg] = ",totalAngle_RP
    if totalAngle_RP>360.0:
        PA_deg_RP_model = totalAngle_RP -360.0
    else:
        PA_deg_RP_model = totalAngle_RP
    
    #print 'longAN_deg = ',longAN_deg
    #print 'totalAngle_RP = ',totalAngle_RP
    if verbose:
        print 'Position Angle in reference plane [deg] = ',PA_deg_RP_model
    
#    if ang_LN_to_M2_deg>360.0:
#        ang_LN_to_M2_deg = ang_LN_to_M2_deg-360.0
#    
#    tanCos = math.atan(math.tan(math.radians(ang_LN_to_M2_deg))*math.cos(inclination_rad))
#    PA_deg_RP_model2 = math.degrees(tanCos)+longAN_deg
#        
#    if (ang_LN_to_M2_deg>0.0) and(ang_LN_to_M2_deg<90.0):
#        #quadrant 1
#        PA_deg_RP_model2_corr = PA_deg_RP_model2
#    elif (ang_LN_to_M2_deg>90.0) and(ang_LN_to_M2_deg<180.0):
#        #quadrant 2
#        PA_deg_RP_model2_corr = PA_deg_RP_model2+180.0
#    elif (ang_LN_to_M2_deg>180.0) and(ang_LN_to_M2_deg<270.0):
#        #quadrant 3
#        PA_deg_RP_model2_corr = PA_deg_RP_model2+180.0
#    elif (ang_LN_to_M2_deg>270.0) and(ang_LN_to_M2_deg<360.0):
#        #quadrant 4
#        PA_deg_RP_model2_corr = PA_deg_RP_model2+360.0
#        
#    if abs(PA_deg_RP_model2_corr-PA_deg_RP_model)>0.000001:
#        print '\norbCalc2-2153: Warning: Two differently calculated Position Angles DO NOT match!!'
#        print '               Will continue to use only PA_deg_RP_model.'
#        print 'PA_deg_RP_model = ',PA_deg_RP_model
#        print 'PA_deg_RP_model2_corr = ',PA_deg_RP_model2_corr
#        print 'PA_deg_RP_model2 = ',PA_deg_RP_model2
#        print 'ang_LN_to_M2_deg = ',ang_LN_to_M2_deg
#        print 'math.tan(math.radians(ang_LN_to_M2_deg)) = '+str(math.tan(math.radians(ang_LN_to_M2_deg)))
        
    ## calculate measured Separation Angle (model)
    SA_arcsec_RP_model = Sep_Dist_AU_RP/Sys_Dist_PC # PC->AU and rad->arcsec cancel
    if verbose:
        print 'Separation Angle measured (model) [arcsec] = ',SA_arcsec_RP_model
    
#    # calculate Separation Angle a second way (from Double Stars pg 35)
#    squaredSines = (math.sin(math.radians(ang_LN_to_M2_deg)))**2.0*(math.sin(inclination_rad))**2.0
#    SA_arcsec_RP_model2 = (Sep_Dist_AU_OP/Sys_Dist_PC)*(1.0-squaredSines)**(1.0/2.0)
#    
#    # check two differently calculated sep angles match
#    if abs(SA_arcsec_RP_model2-SA_arcsec_RP_model)>0.000001:
#        print '\norbCalc2-2043: Warning: Two differently calculated Separation Angles DO NOT match!!'
#        print '               Will continue and use only SA_arcsec_RP_model.'
#        print 'SA_arcsec_RP_model = ',SA_arcsec_RP_model
#        print 'SA_arcsec_RP_model2 = ',SA_arcsec_RP_model2
        
    if (Mass1 is not 1) and (Mass2 is not 1) and (a_total!=False):
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
        
        if (a1+a2)-a_total>1.0e-5:
            print "obCalc-ln1000: a1 = ",a1
            print "obCalc-ln1001:a2 = ",a2
            print "obCalc-ln1002:a_total = ",a_total

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
   
    return (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SA_arcsec_RP_model, PA_deg_RP_model, a1, a2)

def orbitCalculator3(t, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total=False,\
                    Mass1=1, Mass2=1, verbose=False):
    """
    This function will calculate the output orbital parameters from the inputs for a SINGLE epoch.  Use the 
    multiEpochOrbCalc to do this for multiple epochs at the same time.  
    The t and Sys_Dist_PC parameters are generally known/measured values,
    while the rest are typically random numbers.
    This function is designed to work as the calculator inside a Monte Carlo Simulator.
    All the equations are the result of Kepler's Laws and naturally occurring symmetries of binary star systems.
    
    Note: RP stands for Reference Plane, ie the plane that values are measured in from Earth.
          OP stands for Orbital Plane, ie in the plane that the stars are orbiting in.
    
    Inputs:
    @param t:                            = epoch of observation/image [julian date]
    @type t:                             = float
    @param Sys_Dist_PC:                  = measured system distance from Earth [PC]
    @type Sys_Dist_PC:                   = float
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
    @param a_total:                      = Total semi-major axis of system [AU]
    @type a_total:                       = float if it is to be used, boolean False if to be calculated from Kepler's 3rd
    @param Mass1:                        = Mass of primary in system [kg] (Default=1, means M2 is a planet)
    @type Mass1:                         = float
    @param Mass2:                        = Mass of the secondary in system [kg] (Default=1, means M2 is a planet)
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
    @param Sep_Dist_AU_OP:               = Separation Distance in orbital plane [AU]
    @type Sep_Dist_AU_OP:                = float
    @param SA_arcsec_RP_model:           = measured Position Angle in image [deg]
    @type SA_arcsec_RP_model:            = float
    @param PA_deg_RP_model:              = measured Position Angle in image [deg]
    @type PA_deg_RP_model:               = float
    @param a1:                           = semi-major axis of M1 [AU]
    @type a1:                            = float
    @param a2:                           = semi-major axis of M2 [AU]
    @type a2:                            = float
    """
    testing = False
    
    if verbose:
        print '\n'+'*'*50
        print 'Starting to calculate Orbital Parameters\n'
        print 'Input variable values: '+'\nCurrent epoch time [JD] = '+str(t)+\
        '\nSystem Distance from Earth [PC] = '+str(Sys_Dist_PC)+\
        '\nInclination [deg] = '+ str(inclination_deg)+\
        '\nLongitude of Ascending Node [deg] = '+str(longAN_deg)+\
        '\nEccentricity = '+str(e)+'\nLast Periapsis Epoch [JD] = '+\
        str(T)+'\nPeriod [yrs] = '+str(period)+\
        '\nArgument of Periapsis in orbital plane [deg] = '+\
        str(argPeri_deg)+'\nverbose = ',str(verbose)
    #----------------------------------------------------------------------
    
    ## check if a_total is to be calculated using Kepler's 3rd
    if a_total==False:
        
        # conversion factors and constants
        SecPerYear = 31557600.0
        G = 6.67300e-11
        MperAU = 149598000000.0
        KGperMsun = 1.98892e30
        
        # convert units of years to seconds
        period_seconds = period*SecPerYear
        
        # convert masses from Msun to kilograms
        Mass1_kg = Mass1*KGperMsun
        Mass2_kg = Mass2*KGperMsun
        # apply K3
        a_total = (((G*(Mass1_kg+Mass2_kg)*(period_seconds**2.0))/(4.0*(pi**2.0))))**(1.0/3.0) # Masses must be in [kg]!!!!! 
        
        a2 = a_total/(1.0+(Mass2_kg/Mass1_kg))
        
        # find a1 using mass ratio
        a1 = a2*(Mass2_kg/Mass1_kg)
        
        # sum a's and convert into [AU]
        #a_total = (a1+a2)/MperAU
        
    ## get the TAcalculator to find the TA in radians
    (n, M_deg, E_latest_deg,TA_rad) = TAcalculator2(t,e, T, period, verbose=verbose)
    TA_deg = math.degrees(TA_rad)
    if verbose:
        print 'True Anomaly [deg] = ', TA_deg
        
    ## calculate the Separation Distance in orbital plane
    Sep_Dist_AU_OP = (a_total*(1.0-e**2.0))/(1.0+e*math.cos(TA_rad))
    if verbose:
        print 'Separation Distance in orbital plane [AU] = ',Sep_Dist_AU_OP
    
    ## calculate measured Separation Distance in reference plane
    inclination_rad = math.radians(inclination_deg) # convert inclination to radians 
    # angle between ascending node and M2 in orbital plane
    ang_LN_to_M2_deg = argPeri_deg + TA_deg
    if ang_LN_to_M2_deg>360.0:
        ang_LN_to_M2_deg = ang_LN_to_M2_deg-360.0
    ang_LN_to_M2_rad = math.radians(ang_LN_to_M2_deg)
    if verbose:
        print "inclination_deg = ",inclination_deg
        print "inclination_rad = ",inclination_rad
        print "ang_LN_to_M2_rad = ",ang_LN_to_M2_rad
    # find horizontal and vertical WRT Line of Nodes (X and Y) components of Sep_Dist corrected for inclination
    # Note: only the 'vertical (Y)' component will be effected by inclination
    # Note2: on this axis, Z would be the direction of 'line of sight' of the observer.
    Sep_Dist_AU_RP_y = Sep_Dist_AU_OP*math.sin(ang_LN_to_M2_rad)*math.cos(inclination_rad)
    Sep_Dist_AU_RP_x = Sep_Dist_AU_OP*math.cos(ang_LN_to_M2_rad)
    
##    ## calculate the x and y values of the Ascending Node location in RP
##    M1_to_AN_dist_RP_AU = (a_total*(1.0-e**2.0))/(1.0+e*math.cos(math.radians(-argPeri_deg)))
##    
##    y_AN = M1_to_AN_dist_RP_AU*math.cos(math.radians(longAN_deg))
##    x_AN = M1_to_AN_dist_RP_AU*math.sin(math.radians(longAN_deg))
##    
##    print 
##    x_total = Sep_Dist_AU_RP_x+x_AN
##    y_total = Sep_Dist_AU_RP_y+y_AN

    # use Pythagorous theorem to find the final Sep Dist in the Ref Plane
    Sep_Dist_AU_RP = ((Sep_Dist_AU_RP_x**2.0)+(Sep_Dist_AU_RP_y**2.0))**(1.0/2.0)
    if verbose:
        print 'Sep_Dist_AU_RP_y = ',Sep_Dist_AU_RP_y
        print 'Sep_Dist_AU_RP_x = ',Sep_Dist_AU_RP_x
        print 'Separation Distance in the reference plane [AU] = ',Sep_Dist_AU_RP

    ## calculate corrected angle between Line of Nodes and M2 (ie. AN-focus-M2)
    # the first atan will produce an angle between the X-axis, rather than only CC from positive X-axis.
    # must take all 4 quadrants into account.  
    ang_LN_to_M2_corrPOSNEG = math.degrees(math.atan(Sep_Dist_AU_RP_y/Sep_Dist_AU_RP_x))
    
    if (Sep_Dist_AU_RP_x>=0.0)and (Sep_Dist_AU_RP_y>=0.0):
        # quadrant 1
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG 
    elif (Sep_Dist_AU_RP_x<0.0)and (Sep_Dist_AU_RP_y>=0.0):
        # quadrant 2
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +180.0
    elif (Sep_Dist_AU_RP_x<0.0)and (Sep_Dist_AU_RP_y<0.0):
        # quadrant 3
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +180.0
    elif (Sep_Dist_AU_RP_x>=0.0)and (Sep_Dist_AU_RP_y<0.0):
        # quadrant 4
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +360.0
    else:
        print ' The value of '+str(ang_LN_to_M2_corrPOSNEG)+' is apparently out of ranges'
    
    #print '\nuncorrected angle between line of nodes and M2 [deg] = ', ang_LN_to_M2_deg
    #print 'Corrected angle between Line of nodes and M2 [deg] = ',ang_LN_to_M2_corr
    if verbose:
        print 'uncorrected angle between line of nodes and M2 [deg] = ', ang_LN_to_M2_deg
        print 'Corrected angle between Line of nodes and M2 [deg] = ',ang_LN_to_M2_corr
    
    ## calculate measured Position Angle 
    # it is measured from M1, only if the value is over 360deg do we need to fix it
    totalAngle_RP = longAN_deg + ang_LN_to_M2_corr
    if verbose:
        print "Total angle in the Reference Plane [deg] = ",totalAngle_RP
    if totalAngle_RP>360.0:
        PA_deg_RP_model = totalAngle_RP -360.0
    else:
        PA_deg_RP_model = totalAngle_RP
    
    #print 'longAN_deg = ',longAN_deg
    #print 'totalAngle_RP = ',totalAngle_RP
    if verbose:
        print 'Position Angle in reference plane [deg] = ',PA_deg_RP_model
    
#    if ang_LN_to_M2_deg>360.0:
#        ang_LN_to_M2_deg = ang_LN_to_M2_deg-360.0
#    
#    tanCos = math.atan(math.tan(math.radians(ang_LN_to_M2_deg))*math.cos(inclination_rad))
#    PA_deg_RP_model2 = math.degrees(tanCos)+longAN_deg
#        
#    if (ang_LN_to_M2_deg>0.0) and(ang_LN_to_M2_deg<90.0):
#        #quadrant 1
#        PA_deg_RP_model2_corr = PA_deg_RP_model2
#    elif (ang_LN_to_M2_deg>90.0) and(ang_LN_to_M2_deg<180.0):
#        #quadrant 2
#        PA_deg_RP_model2_corr = PA_deg_RP_model2+180.0
#    elif (ang_LN_to_M2_deg>180.0) and(ang_LN_to_M2_deg<270.0):
#        #quadrant 3
#        PA_deg_RP_model2_corr = PA_deg_RP_model2+180.0
#    elif (ang_LN_to_M2_deg>270.0) and(ang_LN_to_M2_deg<360.0):
#        #quadrant 4
#        PA_deg_RP_model2_corr = PA_deg_RP_model2+360.0
#        
#    if abs(PA_deg_RP_model2_corr-PA_deg_RP_model)>0.000001:
#        print '\norbCalc2-2153: Warning: Two differently calculated Position Angles DO NOT match!!'
#        print '               Will continue to use only PA_deg_RP_model.'
#        print 'PA_deg_RP_model = ',PA_deg_RP_model
#        print 'PA_deg_RP_model2_corr = ',PA_deg_RP_model2_corr
#        print 'PA_deg_RP_model2 = ',PA_deg_RP_model2
#        print 'ang_LN_to_M2_deg = ',ang_LN_to_M2_deg
##        print 'math.tan(math.radians(ang_LN_to_M2_deg)) = '+str(math.tan(math.radians(ang_LN_to_M2_deg)))
#        
    ## calculate measured Separation Angle (model)
    SA_arcsec_RP_model = Sep_Dist_AU_RP/Sys_Dist_PC # PC->AU and rad->arcsec cancel
    
    if verbose:
        print 'Separation Angle measured (model) [arcsec] = ',SA_arcsec_RP_model
    
#    # calculate Separation Angle a second way (from Double Stars pg 35)
#    squaredSines = (math.sin(math.radians(ang_LN_to_M2_deg)))**2.0*(math.sin(inclination_rad))**2.0
#    SA_arcsec_RP_model2 = (Sep_Dist_AU_OP/Sys_Dist_PC)*(1.0-squaredSines)**(1.0/2.0)
#    
#    # check two differently calculated sep angles match
#    if abs(SA_arcsec_RP_model2-SA_arcsec_RP_model)>0.000001:
#        print '\norbCalc2-2043: Warning: Two differently calculated Separation Angles DO NOT match!!'
#        print '               Will continue and use only SA_arcsec_RP_model.'
#        print 'SA_arcsec_RP_model = ',SA_arcsec_RP_model
#        print 'SA_arcsec_RP_model2 = ',SA_arcsec_RP_model2
        
    x2 = SA_arcsec_RP_model*math.cos(math.radians(PA_deg_RP_model))
    y2 = SA_arcsec_RP_model*math.sin(math.radians(PA_deg_RP_model))
    if testing:
        print '\n################################################################################################'
        print 'SA_arcsec_RP_model = ',SA_arcsec_RP_model
        print 'PA_deg_RP_model = ',PA_deg_RP_model
        print '\nx2 = SA_arcsec_RP_model*math.cos(math.radians(PA_deg_RP_model)) = ',x2
        print 'y2 = SA_arcsec_RP_model*math.sin(math.radians(PA_deg_RP_model)) = ',y2
    if verbose:
        print 'x2 = SA_arcsec_RP_model*math.cos(math.radians(PA_deg_RP_model)) = ',x2
        print 'y2 = SA_arcsec_RP_model*math.sin(math.radians(PA_deg_RP_model)) = ',y2
    
    if (Mass1 is not 1) and (Mass2 is not 1) and (a_total!=False):
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
        
        if (a1+a2)-a_total>1.0e-5:
            print "obCalc-ln1000: a1 = ",a1
            print "obCalc-ln1001:a2 = ",a2
            print "obCalc-ln1002:a_total = ",a_total

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
   
    return (n, M_deg, E_latest_deg, TA_deg, x2, y2, a1, a2)

def orbitCalculatorTH_I(A, B, C, F, G, Mass1_Msun, Mass2_Msun, Sys_Dist_PC):
    """
    """
    rad_per_arcsec = 206264.806  # Extra? 
    
    # calc semi-major axis of primary
    a1_arcsec = math.sqrt(A**2.0+B**2.0+C**2.0)
    
    ## calc argPeri of primary
    #$$$ NOTE: might be a more efficient way to do this, but couldn't see it 
    #$$$       see it immediately, so using this verbose method for now.
    temp1 = math.atan((B-F)/(A+G))
    # convert angle to match quadrant it should be
    quadrant = 0
    if ((B-F)>=0.0)and(((B-F)/(A+G))>=0.0):
        quadrant = 1
    elif ((B-F)>=0.0)and(((B-F)/(A+G))<=0.0):
        quadrant = 2
    elif ((B-F)<=0.0)and(((B-F)/(A+G))>=0.0):
        quadrant = 3
    elif ((B-F)<=0.0)and(((B-F)/(A+G))<=0.0):
        quadrant = 4
#    print 'temp1 before = ',math.degrees(temp1) #$$$$$$$$$$$$$$$$$$
    if (quadrant==2) or (quadrant==3):    
        temp1 = temp1 + pi
    elif quadrant==4:
        temp1 = temp1 + 2.0*pi
    
#    print 'temp1 after = ',math.degrees(temp1) #$$$$$$$$$$$$$$$$
        
        
    temp2 = math.atan((-B-F)/(A-G))
    # convert angle to match quadrant it should be
    quadrant = 0
    if ((-B-F)>=0.0)and(((-B-F)/(A-G))>=0.0):
        quadrant = 1
    elif ((-B-F)>=0.0)and(((-B-F)/(A-G))<=0.0):
        quadrant = 2
    elif ((-B-F)<=0.0)and(((-B-F)/(A-G))>=0.0):
        quadrant = 3
    elif ((-B-F)<=0.0)and(((-B-F)/(A-G))<=0.0):
        quadrant = 4
        
    if (quadrant==2) or (quadrant==3):    
        temp2 = temp2 + pi
    elif quadrant==4:
        temp2 = temp2 + 2.0*pi
    
    # put two corrected angles together to calc final argPeri of primary
    argPeri_primary_rad = temp1+temp2
    
#    print 'argPeri_primary_rad before = ',math.degrees(argPeri_primary_rad) #$$$$$$$$$$$$$$$$$$$$$$$$
    # convert angle to be inside 2*pi
    while argPeri_primary_rad>(2.0*pi):
    	argPeri_primary_rad = argPeri_primary_rad-2.0*pi
    if argPeri_primary_rad<0.0:
    	argPeri_primary_rad = argPeri_primary_rad +2.0*pi
#    print 'argPeri_primary_rad after = ',math.degrees(argPeri_primary_rad) #$$$$$$$$$$$$$$$$$$$$$$$$
    
    # add pi radians to convert primary's argPeri to companion's
    argPeri_rad = argPeri_primary_rad + pi
    
    
#    print 'argPeri_rad before = ',math.degrees(argPeri_rad) #$$$$$$$$$$$$$$$$$$$$$$$$
    #convert angle to be inside 2.0*pi
    while argPeri_rad>(2.0*pi):
    	argPeri_rad = argPeri_rad-2.0*pi
#    print 'argPeri_rad after = ',math.degrees(argPeri_rad) #$$$$$$$$$$$$$$$$$$$$$$$$
    
    # Calc longitude of Ascending Node
    longAN_rad = temp1-argPeri_primary_rad
    
    if longAN_rad>pi:
    	longAN_rad = longAN_rad -pi
    
    # Calc inclination
    temp3 = (A-G)*math.cos(argPeri_primary_rad+longAN_rad)
    temp4 = (A+G)*math.cos(argPeri_primary_rad-longAN_rad)
#    print 'A-G = ',(A-G)+
#    print 'temp3 = ',temp3
#    print 'A+G = ',(A+G)
#    print 'temp4 = ',temp4
    inclination_rad = 2.0*math.atan(math.sqrt(abs(temp3/temp4)))
    
#    # calculate inclination another way to ensure it is the same
#    print 'C = ',C
#    print 'a1_arcsec ',a1_arcsec
#    print 'math.sin(argPeri_primary_rad) ',str(math.sin(argPeri_primary_rad))
#    print 'a1_arcsec*math.sin(argPeri_primary_rad) = ',str(a1_arcsec*math.sin(argPeri_primary_rad))
#    print 'abs(C/(a1_arcsec*math.sin(argPeri_primary_rad)))-1.0 ',str(abs(C/(a1_arcsec*math.sin(argPeri_primary_rad)))-1.0)
#    
#    if abs(C/(a1_arcsec*math.sin(argPeri_primary_rad)))>1.0:
#    	inside = abs(C/(a1_arcsec*math.sin(argPeri_primary_rad)))-1.0
#    if abs(C/(a1_arcsec*math.sin(argPeri_primary_rad)))>2.0:
#    	inside = abs(C/(a1_arcsec*math.sin(argPeri_primary_rad)))-2.0
#    
#    print 'inside = ',inside
#    inclination_rad2 = math.asin(inside)
#    
#    # check both incliantions are the same
#    if (inclination_rad-inclination_rad2)>0.00001:
#    	print 'orbCalculatorTH_I2023: WARNING! two inclinations calculated differ by more than 0.00001 !!!'
#    	print 'inclination_rad = ',inclination_rad
#    	print 'inclination_rad2 = ',inclination_rad2
    
    # Calc the period from Kepler's 3rd law
    # NOTE: multiplying by the system distance in [PC] converts the ["] distances into [AU] to avoid
    #       raising a value less than 1.0 to the third power which has the opposite effect than over 1.0.
    #       the units cancel in the end so it just makes sure the ratio works in the right way.
    
    # conversion factors and constants
    SecPerYear = 31557600.0
    G = 6.67300e-11
    MperAU = 149598000000.0
    KGperMsun = 1.98892e30
    
    temp5 = 4.0*pi**2.0*(Sys_Dist_PC*a1_arcsec*MperAU)**3.0*Mass1_Msun**3.0
    temp6 = G*KGperMsun*(Mass1_Msun+Mass2_Msun)*Mass2_Msun**3.0
    period = (math.sqrt(temp5/temp6))/SecPerYear
    
    return (a1_arcsec, argPeri_rad, longAN_rad, inclination_rad, period)

def orbitCalculatorTH_I2(A, B, C, F, G, Mass1_Msun, Mass2_Msun, Sys_Dist_PC):
    """
    version using ideas for calculating argPeri, longAN and a from 
    'The Binary Stars' by robert G. Aitken.
    """
    
    ## calc argPeri of primary
    #$$$ NOTE: might be a more efficient way to do this, but couldn't see it 
    #$$$       see it immediately, so using this verbose method for now.
    argPeri_plus_longAN = math.atan((B-F)/(A+G))
    # convert angle to match quadrant it should be
    quadrant = 0
    if ((B-F)>=0.0)and(((B-F)/(A+G))>=0.0):
        quadrant = 1
    elif ((B-F)>=0.0)and(((B-F)/(A+G))<=0.0):
        quadrant = 2
    elif ((B-F)<=0.0)and(((B-F)/(A+G))>=0.0):
        quadrant = 3
    elif ((B-F)<=0.0)and(((B-F)/(A+G))<=0.0):
        quadrant = 4
        
    if (quadrant==2) or (quadrant==3):    
        argPeri_plus_longAN = argPeri_plus_longAN + pi
    elif quadrant==4:
        argPeri_plus_longAN = argPeri_plus_longAN + 2.0*pi
    
        
        
    argPeri_minus_longAN = math.atan((-B-F)/(A-G))
    # convert angle to match quadrant it should be
    quadrant = 0
    if ((-B-F)>=0.0)and(((-B-F)/(A-G))>=0.0):
        quadrant = 1
    elif ((-B-F)>=0.0)and(((-B-F)/(A-G))<=0.0):
        quadrant = 2
    elif ((-B-F)<=0.0)and(((-B-F)/(A-G))>=0.0):
        quadrant = 3
    elif ((-B-F)<=0.0)and(((-B-F)/(A-G))<=0.0):
        quadrant = 4
        
    if (quadrant==2) or (quadrant==3):    
        argPeri_minus_longAN = argPeri_minus_longAN + pi
    elif quadrant==4:
        argPeri_minus_longAN = argPeri_minus_longAN + 2.0*pi
    
    # new way to calculate longAN and arePeri
    double_argPeri = argPeri_plus_longAN + argPeri_minus_longAN
    argPeri_primary_rad = double_argPeri/2.0
    double_longAN = argPeri_plus_longAN - argPeri_minus_longAN
    longAN_rad = double_longAN/2.0
    
    # make corrections for if longAN is out of accepted 0-180deg range
    if longAN_rad>pi:
        longAN_rad = longAN_rad-pi
        argPeri_primary_rad = argPeri_primary_rad-pi
    if (longAN_rad<0.0) and (longAN_rad>(-pi)):
        longAN_rad = longAN_rad+pi
        argPeri_primary_rad = argPeri_primary_rad+pi
#    if argPeri_primary_rad>(2.0*pi):
#        argPeri_primary_rad = argPeri_primary_rad-2.0*pi
        
    # add pi radians to convert primary's argPeri to companion's
    argPeri_rad = argPeri_primary_rad 
    
    # Calc inclination
    temp3 = (A-G)*math.cos(argPeri_primary_rad+longAN_rad)
    temp4 = (A+G)*math.cos(argPeri_primary_rad-longAN_rad)

    inclination_rad = 2.0*math.atan(math.sqrt(abs(temp3/temp4)))
    
#    # convert inc to be in 0-90 range          #$$$$$$$$$$$$$$$$$$$$
#    if inclination_rad> (pi/2.0):              #$$$$$$$$$$$$$$$$$$$$$$
#        inclination_rad = pi-inclination_rad   #$$$$$$$$$$$$$$$$$$$$$
    
#    tempA = (A+G)/math.cos(argPeri_primary_rad+longAN_rad)
#    tempB = (A-G)/math.cos(argPeri_primary_rad-longAN_rad)
#    inside = ((2.0*tempA)/(tempA+tempB))-1.0
#    inclination_rad2 = math.acos(inside)
#    
#    print 'inclination_rad = '+str(inclination_rad)+' , or in deg = '+str(math.degrees(inclination_rad))
#    print 'inclination_rad2 = '+str(inclination_rad2)+' , or in deg = '+str(math.degrees(inclination_rad2))+', and inside = '+str(inside)
    
    # calc semi-major axis of apparent ellipse
    a_arcsec_try1 = (B-F)/(2.0*math.sin(longAN_rad+argPeri_primary_rad)*(math.cos(inclination_rad/2.0))**2.0)
#    a_try1 = a_arcsec_try1#*Sys_Dist_PC
#    print 'a_try1 = ',a_try1
    a_arcsec = a_arcsec_try1
#    
#    a_arcsec_try2 = ((A*G-B*F)/(math.cos(inclination_rad)))**(1.0/2.0)
#    a_try2 = a_arcsec_try2#*Sys_Dist_PC
#    print 'a_try2 = ',a_try2
#    
#    a_arcsec_try3 = ((A**2.0+B**2.0+F**2.0+G**2.0)/(1.0+(math.cos(inclination_rad))**2.0))**(1.0/2.0)
#    a_try3 = a_arcsec_try3#*Sys_Dist_PC
#    print 'a_try3 = ',a_try3
    
    # Calc the period from Kepler's 3rd law
    # NOTE: multiplying by the system distance in [PC] converts the ["] distances into [AU] to avoid
    #       raising a value less than 1.0 to the third power which has the opposite effect than over 1.0.
    #       the units cancel in the end so it just makes sure the ratio works in the right way.
    
    # conversion factors and constants
    SecPerYear = 31557600.0
    G = 6.67300e-11
    MperAU = 149598000000.0
    KGperMsun = 1.98892e30
    
    temp5 = 4.0*pi**2.0*(Sys_Dist_PC*a_arcsec*MperAU)**3.0
    temp6 = G*KGperMsun*(Mass1_Msun+Mass2_Msun)
    period = (math.sqrt(temp5/temp6))/SecPerYear
    
    return (a_arcsec, argPeri_rad, longAN_rad, inclination_rad, period)


def orbitCalculatorTH_I_PP(A, B, C, F, G, Mass1_Msun, Mass2_Msun, Sys_Dist_PC):
    """
    Way orbital elements are found from the Thiele-Innes constants in PP1995
    """
    k = A**2.0 + B**2.0 + F**2.0 + G**2.0
    m = A*G-B*F
    j = math.sqrt(k**2.0-m**2.0)
    z = math.atan((B-F)/(A+G))
    r = math.atan((B+F)/(G-A))
    
    a_arcsec = math.sqrt(j+k)
    inclination_rad = math.atan(math.sqrt(a_arcsec**4.0+m**2.0)/(m))
    if inclination_rad<0.0:
        inclination_rad = inclination_rad+pi
    argPeri_rad = (z+r)/2.0 +pi
    longAN_rad = (z-r)/2.0+pi
    
    # conversion factors and constants
    SecPerYear = 31557600.0
    G = 6.67300e-11
    MperAU = 149598000000.0
    KGperMsun = 1.98892e30
    
    temp5 = 4.0*pi**2.0*(Sys_Dist_PC*a_arcsec*MperAU)**3.0*Mass1_Msun**3.0
    temp6 = G*KGperMsun*(Mass1_Msun+Mass2_Msun)*Mass2_Msun**3.0
    period = (math.sqrt(temp5/temp6))/SecPerYear
    
    return (a_arcsec, argPeri_rad, longAN_rad, inclination_rad, period)

def multiEpochOrbCalcTH_I(e, T, A, B, C, F, G, Mass1_Msun, Mass2_Msun, Sys_Dist_PC, \
                                                    SAs, PAs, epochs, SA_errors, PA_errors):
    """
    """
    ## Call the orbitCalculatorTH_I to give us the params that stay stable throughout the orbit
    (a_arcsec, argPeri_rad, longAN_rad, inclination_rad, period) = \
            orbitCalculatorTH_I2(A, B, C, F, G, Mass1_Msun, Mass2_Msun, Sys_Dist_PC)
    
    
    ## Loop through all epochs of data to calc the predicted (x,y) and use to calc chiSquared
    # set up initial value for chiSquared to be updated in loop
    chiSquared_total = 0
    if (len(SAs)==len(PAs)) and (len(epochs)==len(SAs)):
        for i in range(0,len(epochs)):
            
            # use TAcalculator to get orbital params for current epoch
           (n, M_deg, E_latest_deg,TA_rad) = TAcalculator(epochs[i], e, T, period, verbose=False)
           
           # Calc Normalized rectangular coordinates
           X = math.cos(math.radians(E_latest_deg))-e
           Y = math.sqrt(1.0-e**2.0)*math.sin(math.radians(E_latest_deg))
           
           # Calc x,y values on same coord system as plane of sky (same as data)
           x_model = A*X+F*Y
           y_model = B*X+G*Y
           
           # Calc x,y values from SA and PA of data
           x_data = SAs[i]*math.cos(math.radians(PAs[i]))
           y_data = SAs[i]*math.sin(math.radians(PAs[i]))
           
           # calc error in x,y data
           #x_data_error = x_data*((SA_errors[i]/SAs[i])+abs(math.radians(PA_errors[i])*math.tan(math.radians(PAs[i]))))
           x_data_error2 = ((SA_errors[i]*math.cos(math.radians(PA_errors[i])))**2.0+(SAs[i]*math.sin(math.radians(PAs[i])*PA_errors[i]))**2.0)**(1.0/2.0)
           #y_data_error = y_data*((SA_errors[i]/SAs[i])+abs(math.radians(PA_errors[i])/math.tan(math.radians(PAs[i]))))
           y_data_error2 = ((SA_errors[i]*math.sin(math.radians(PA_errors[i])))**2.0+(SAs[i]*math.cos(math.radians(PAs[i])*PA_errors[i]))**2.0)**(1.0/2.0)
           # calc chiSquared for both x and y
           chiSquared_x = chiSquaredCalc(x_data, x_data_error2, x_model)
           chiSquared_y = chiSquaredCalc(y_data, y_data_error2, y_model)
           
           # add both of those to running chiSquared total
           chiSquared_total = chiSquared_total+chiSquared_x+chiSquared_y
           
    else:
        print 'multiOrbcalcTH_I: WARNING! Number of epochs, SAs and/or PAs data does not match'
    
    return (a_arcsec, argPeri_rad, longAN_rad, inclination_rad, period, chiSquared_total)

def multiEpochOrbCalcTH_I3(e, T, A, B, C, F, G, Mass1_Msun, Mass2_Msun, Sys_Dist_PC, \
                                                    SAs, PAs, epochs, SA_errors, PA_errors):
    """
    """
    testing = False
    
    ## Call the orbitCalculatorTH_I to give us the params that stay stable throughout the orbit
    (a_arcsec, argPeri_rad, longAN_rad, inclination_rad, period) = \
            orbitCalculatorTH_I2(A, B, C, F, G, Mass1_Msun, Mass2_Msun, Sys_Dist_PC)
    
    
    ## Loop through all epochs of data to calc the predicted (x,y) and use to calc chiSquared
    # set up initial value for chiSquared to be updated in loop
    chiSquared_total = 0
    chiSquared_total2 = 0
    if (len(SAs)==len(PAs)) and (len(epochs)==len(SAs)):
        for i in range(0,len(epochs)):
            
            # use TAcalculator to get orbital params for current epoch
           (n, M_deg, E_latest_deg,TA_rad) = TAcalculator(epochs[i], e, T, period, verbose=False)
           
           # Calc Normalized rectangular coordinates
           X = math.cos(math.radians(E_latest_deg))-e
           Y = math.sqrt(1.0-e**2.0)*math.sin(math.radians(E_latest_deg))
           
           # Calc x,y values on same coord system as plane of sky (same as data)
           x_model = A*X+F*Y
           y_model = B*X+G*Y
           
           
           # Calc x,y values from SA and PA of data
           x_data = SAs[i]*math.cos(math.radians(PAs[i]))
           y_data = SAs[i]*math.sin(math.radians(PAs[i]))
           
           if testing:
               print '\n#########################################################'
               print 'SAs[i] = ',SAs[i]
               print 'PAs[i] = ',PAs[i]
               print '\nx_model = A*X+F*Y = ',x_model
               print 'y_model = B*X+G*Y = ',y_model
               print 'x_data = SAs[i]*math.cos(math.radians(PAs[i])) = ',x_data
               print 'y_data = SAs[i]*math.sin(math.radians(PAs[i])) = ',y_data
               print '#########################################################'
           # calc error in x,y data
           x_data_error = x_data*((SA_errors[i]/SAs[i])+abs(math.radians(PA_errors[i])*math.tan(math.radians(PAs[i]))))
           y_data_error = y_data*((SA_errors[i]/SAs[i])+abs(math.radians(PA_errors[i])/math.tan(math.radians(PAs[i]))))
           
           # calc chiSquared for both x and y
           chiSquared_x = chiSquaredCalc(x_data, x_data_error, x_model)
           chiSquared_y = chiSquaredCalc(y_data, y_data_error, y_model)
           
           # add both of those to running chiSquared total
           chiSquared_total = chiSquared_total+chiSquared_x+chiSquared_y
           
           # calc error in x,y data another way
           x_data_error2 = x_data*((SA_errors[i]/SAs[i])+((math.cos( math.radians(PAs[i]+PA_errors[i]))-math.cos( math.radians(PAs[i]) ))/math.cos( math.radians(PAs[i]) ) ))
           y_data_error2 = y_data*((SA_errors[i]/SAs[i])+((math.sin( math.radians(PAs[i]+PA_errors[i]))-math.sin( math.radians(PAs[i]) ))/math.sin( math.radians(PAs[i]) ) ))
           
#           print 'x_data_error = ',x_data_error
#           print 'y_data_error = ',y_data_error
#           print 'x_data_error2 = ',x_data_error2
#           print 'y_data_error2 = ',y_data_error2
#           
#           # calc chiSquared for both x and y
#           chiSquared_x2 = chiSquaredCalc(x_data, x_data_error2, x_model)
#           chiSquared_y2 = chiSquaredCalc(y_data, y_data_error2, y_model)
           
           # add both of those to running chiSquared total
           chiSquared_total2 = chiSquared_total2+chiSquared_x2+chiSquared_y2
           
    else:
        print 'multiOrbcalcTH_I: WARNING! Number of epochs, SAs and/or PAs data does not match'
#    
#    print 'chiSquared_total = ',chiSquared_total
#    print 'chiSquared_total2 = ',chiSquared_total2
    
    return (a_arcsec, argPeri_rad, longAN_rad, inclination_rad, period, chiSquared_total)
        
def ABCFG_MaxMins(aMax, aMin, argPeri_radMax, argPeri_radMin, longAN_radMax, longAN_radMin, \
                                                                inclination_radMax, inclination_radMin):
    """
    """
    # calc output maxs from input maxs
    (Amax,Bmax,Cmax,Fmax,Gmax) = ABCFG_values(aMax, argPeri_radMax, longAN_radMax, inclination_radMax)
    
    # calc output mins from input mins
    (Amin,Bmin,Cmin,Fmin,Gmin) = ABCFG_values(aMin, argPeri_radMin, longAN_radMin, inclination_radMin)
    
    # Not 100% sure why this happens, but some mins are greater than maxs, so flip them around.
    if Amax<Amin:
        tempH = Amax
        Amax = Amin
        Amin = tempH
    if Bmax<Bmin:
        tempH = Bmax
        Bmax = Bmin
        Bmin = tempH
    if Cmax<Cmin:
        tempH = Cmax
        Cmax = Cmin
        Cmin = tempH
    if Fmax<Fmin:
        tempH = Fmax
        Fmax = Fmin
        Fmin = tempH    
    if Gmax<Gmin:
        tempH = Gmax
        Gmax = Gmin
        Gmin = tempH    

    return (Amax,Amin, Bmax,Bmin,  Cmax,Cmin, Fmax,Fmin, Gmax,Gmin)

def ABCFG_MaxMins2(aMax, aMin, argPeri_radMax, argPeri_radMin, longAN_radMax, longAN_radMin, \
                                                                inclination_radMax, inclination_radMin):
    """
    """
    # find values for all 8 permutations of variable combinations
    (A1,B1,C1,F1,G1) = ABCFG_values(aMin,argPeri_radMax, longAN_radMax, inclination_radMax)
    (A2,B2,C2,F2,G2) = ABCFG_values(aMin,argPeri_radMax, longAN_radMax, inclination_radMin)
    (A3,B3,C3,F3,G3) = ABCFG_values(aMin,argPeri_radMin, longAN_radMax, inclination_radMax)
    (A4,B4,C4,F4,G4) = ABCFG_values(aMin,argPeri_radMin, longAN_radMax, inclination_radMin)
    (A5,B5,C5,F5,G5) = ABCFG_values(aMin,argPeri_radMax, longAN_radMin, inclination_radMax)
    (A6,B6,C6,F6,G6) = ABCFG_values(aMin,argPeri_radMax, longAN_radMin, inclination_radMin)
    (A7,B7,C7,F7,G7) = ABCFG_values(aMin,argPeri_radMin, longAN_radMin, inclination_radMin)
    (A8,B8,C8,F8,G8) = ABCFG_values(aMin,argPeri_radMin, longAN_radMin, inclination_radMax)
    
    As = [A1,A2,A3,A4,A5,A6,A7,A8]
    Bs = [B1,B2,B3,B4,B5,B6,B7,B8]
    Cs = [C1,C2,C3,C4,C5,C6,C7,C8]
    Fs = [F1,F2,F3,F4,F5,F6,F7,F8]
    Gs = [G1,G2,G3,G4,G5,G6,G7,G8]
    
    Amax = np.max(As)*(aMax/aMin)
    Amin = np.min(As)
    Bmax = np.max(Bs)*(aMax/aMin)
    Bmin = np.min(Bs)
    Cmax = np.max(Cs)*(aMax/aMin)
    Cmin = np.min(Cs)
    Fmax = np.max(Fs)*(aMax/aMin)
    Fmin = np.min(Fs)
    Gmax = np.max(Gs)*(aMax/aMin)
    Gmin = np.min(Gs)
    # make each max and min the absolute max and min 
    # which depends on if they are positive or negative
    if Amax>0.0:
     	Amax = np.max(As)*(aMax/aMin)
    else:
    	Amax = np.max(As)
    if Amin<0.0:
    	Amin = np.min(As)*(aMax/aMin)
    else:
     	Amin = np.min(As)
    
    if Bmax>0.0:
     	Bmax = np.max(Bs)*(aMax/aMin)
    else:
    	Bmax = np.max(Bs)
    if Bmin<0.0:
    	Bmin = np.min(Bs)*(aMax/aMin)
    else:
     	Bmin = np.min(Bs)
     	
    if Cmax>0.0:
     	Cmax = np.max(Cs)*(aMax/aMin)
    else:
    	Cmax = np.max(Cs)
    if Cmin<0.0:
    	Cmin = np.min(Cs)*(aMax/aMin)
    else:
     	Cmin = np.min(Cs)
     	
    if Fmax>0.0:
     	Fmax = np.max(Fs)*(aMax/aMin)
    else:
    	Fmax = np.max(Fs)
    if Fmin<0.0:
    	Fmin = np.min(Fs)*(aMax/aMin)
    else:
     	Fmin = np.min(Fs)
     
    if Gmax>0.0:
     	Gmax = np.max(Gs)*(aMax/aMin)
    else:
    	Gmax = np.max(Gs)
    if Gmin<0.0:
    	Gmin = np.min(Gs)*(aMax/aMin)
    else:
     	Gmin = np.min(Gs)
         
#    Amax = -a1Max
#    Amin = -a1Max
#    Bmax = a1Max
#    Bmin = -a1Max
#    Cmax = a1Max
#    Cmin = -a1Max
#    Fmax = a1Max
#    Fmin = -a1Max
#    Gmax = a1Max
#    Gmin = -a1Max
    return (Amax,Amin, Bmax,Bmin,  Cmax,Cmin, Fmax,Fmin, Gmax,Gmin)
   
   
def ABCFG_values(a, argPeri_rad, longAN_rad, inclination_rad):
    
#    print '\nin ABCFG_values'
#    print 'a1',a1
#    print 'argPeri_rad',argPeri_rad
#    print 'longAN_rad',longAN_rad
#    print 'inclination_rad',inclination_rad
    
    A = a*(math.cos(longAN_rad)*math.cos(argPeri_rad)-\
                   math.sin(longAN_rad)*math.sin(argPeri_rad)*math.cos(inclination_rad))

    B = a*(math.sin(longAN_rad)*math.cos(argPeri_rad)+\
                   math.cos(longAN_rad)*math.sin(argPeri_rad)*math.cos(inclination_rad))
    
    F = a*(-math.cos(longAN_rad)*math.sin(argPeri_rad)-\
                   math.sin(longAN_rad)*math.cos(argPeri_rad)*math.cos(inclination_rad))
    
    G = a*(-math.sin(longAN_rad)*math.sin(argPeri_rad)+\
                   math.cos(longAN_rad)*math.cos(argPeri_rad)*math.cos(inclination_rad))
    
    C = a*math.sin(argPeri_rad)*math.sin(inclination_rad)
    
    return (A,B,C,F,G) 
    
def semiMajorConverter(Mass1, Mass2, a_total=0.0,a1=0.0,a2=0.0, period=0.0, verbose=False):
    """
    Notes:
    Masses must be in same units!  Both MUST be provided for this func to work!!
    
    Units of output semi-majors will be same as input semi-major units.
    only one of the semi-majors needs to be input and the other two will be 
    calculate from that one.
    
    period in years.
    
    only set period to a non-zero value if you do not have any of the semi-majors 
    calculated/provided.  It will use Kepler's third law to find the a_total and 
    then the individual values using the mass ratio.
    """
    # conversion factors and constants
    SecPerYear = 31557600.0
    G = 6.67300e-11
    MperAU = 149598000000.0
    KGperMsun = 1.98892e30
    
    if verbose:
        print '\n**  In SemiMajorConverter **'
        print 'Inputs:\nMass1 = '+str(Mass1)+"\nMass2 = "+str(Mass2)
        print 'a_total = '+str(a_total)+'\na1 = '+str(a1)+'\na2 = '+str(a2)
        print 'period = '+str(period)
    if (a1==a2==a_total==0.0)and(period!=0.0):       
        # convert units of years to seconds
        period_seconds = period*SecPerYear
        
        # convert masses from Msun to kilograms
        Mass1_kg = Mass1*KGperMsun
        Mass2_kg = Mass2*KGperMsun
        # apply K3
        a_total = (((G*(Mass1_kg+Mass2_kg)*(period_seconds**2.0))/(4.0*(pi**2.0))))**(1.0/3.0) # Masses must be in [kg]!!!!! 
        a_total = a_total/MperAU
        #print 'a_total just calculated using K3 to equal ',a_total

    if a1!=0.0:
        a2 = a1*(Mass1/Mass2)
        a_total = a1+a2
    if a2!=0.0:
        a1=a2*(Mass2/Mass1)
        a_total = a1+a2
    if a_total!=0.0:
        a2=a_total/(1.0+(Mass2/Mass1))
        a1=a2*(Mass2/Mass1)
    
    if (period==0.0):
        period_seconds = math.sqrt((4.0*pi**2.0*a_total**3.0)/(G*(Mass1_kg+Mass2_kg)))
        period = period_seconds/SecPerYear
    
    if verbose:
        print '\nOutputs:\nMass1 = '+str(Mass1)+"\nMass2 = "+str(Mass2)
        print 'a_total = '+str(a_total)+'\na1 = '+str(a1)+'\na2 = '+str(a2)
        print 'period = '+str(period)
        print '**  DONE SemiMajorConverter **'
    
    return (a_total, a1, a2, period)
    
def PASAcalculator(period, t, T, e, inclination_deg, longAN_deg, argPeri_deg, Sys_Dist_PC, a2, a1=0, verbose=True):
    '''
    This function is to predict/calculate the expected Position Angle and Separation Angle for a given set of 
    of orbital parameters and epoch of observation.
    
    NOTE: 
    This function is capable of working for 1 or 2 body systems.  In order to use it for 2 body systems, just
    set the values of a1 % a2.  For 1 body systems, ie, only one of the bodies is moving and the primary is 
    simply stationary at one of the foci, then leave a1=0 and use a2 for the moving body's semi-major axis value.
    
    Inputs:
    
    @param period:             = period of orbits [yrs]
    @type period:              = float
    @param t:                  = epoch of observation/image [julian date]
    @type t:                   = float
    @param T:                  = Last Periapsis Epoch/time [julian date] 
    @type T:                   = float
    @param e:                  = eccentricity of orbit
    @type e:                   = float
    @param inclination_deg:    = orbit's inclination [degrees]
    @type inclination_deg:     = float
    @param longAN_deg:         = Longitude of the Ascending Node [degrees]
    @type longAN_deg:          = float
    @param argPeri_deg:        = Argument of periapsis [degrees]
    @type argPeri_deg:         = float
    @param Sys_Dist_PC:        = Distance to the system from Earth [parsec]
    @type Sys_Dist_PC:         = float
    @param a2:                 = Semi-major axis of body 2's orbit [AU]
    @type a2:                  = float
    @param a1:                 = Semi-major axis of body 1's orbit [AU]
    @type a1:                  = float, default=0 indicates body 2 is only moving body in system
    @param verbose:            = Show progress prints to screen?
    @type verbose:             = python boolean (True/False), default=True
    
    Outputs:
    @param PA_deg_RP_model:    = Position Angle [degrees]
    @type PA_deg_RP_model:     = float
    @param SA_arcsec_RP_model: = Separation Angle [arcsec]
    @type SA_arcsec_RP_model:  = float
    '''

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
    
    ## calculate the Mean Motion
    n = (2*pi)/period
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
    
    count = 0 # a counter to stop inf loops in Newtons method below
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
    TA_rad  = math.acos((math.cos(E_latest)-e)/(1-e*math.cos(E_latest))) 
    if E_latest>pi:
        # this is to take care of the silly fact that after E=180deg,
        # the equation above produces TAs that go down from 180, rather than up.
        TA_rad = 2.0*pi - TA_rad
    TA_deg = math.degrees(TA_rad)
    if verbose:
        print 'True Anomaly [deg] = ', TA_deg
        
    ## calculate the Separation Distance in orbital plane
    a_total = a1+a2
    Sep_Dist_AU_OP = a_total*((1-e*e))/(1+e*math.cos(TA_rad))
    if verbose: #$$
        print '(1-e*e) = '+str((1-e*e))#$$$
        print '(1+e*cos(TA)) = '+str((1+e*math.cos(TA_rad)))#$$
    if verbose:
        print 'Separation Distance in orbital plane [AU] = ',Sep_Dist_AU_OP
    
    ## calculate measured Separation Distance in reference plane
    inclination_rad = math.radians(inclination_deg) # convert inclination to radians 
    # angle between ascending node and M2 in orbital plane
    ang_LN_to_M2_deg = argPeri_deg + TA_deg
    # find horizontal and vertical WRT Line of Nodes (X and Y) components of Sep_Dist corrected for inclination
    # Note: only the 'vertical (Y)' component will be effected by inclination
    
    Sep_Dist_AU_RP_y = Sep_Dist_AU_OP*math.sin(math.radians(ang_LN_to_M2_deg))*math.cos(inclination_rad)
    Sep_Dist_AU_RP_x = Sep_Dist_AU_OP*math.cos(math.radians(ang_LN_to_M2_deg))
    
#    if verbose:##$$$
#        print 'sin(TA+argPeri) = '+str(math.sin(math.radians(ang_LN_to_M2_deg)))##$$$
#        print 'cos(TA+argPeri) = '+str(math.cos(math.radians(ang_LN_to_M2_deg)))#$$$
#        print 'cos(i) = '+str(math.cos(inclination_rad))#$$$
#        print 'sep_dist_au_op_y = '+str(Sep_Dist_AU_OP*math.sin(math.radians(ang_LN_to_M2_deg)))
#        print 'sep_dist_au_rp_y = '+str(Sep_Dist_AU_RP_y)#$$$
#        print 'sep_dist_au_rp_x = op_x = '+str(Sep_Dist_AU_RP_x)#$$$
        
    # use Pythagorous theorem to find the final Sep Dist in the Ref Plane
    Sep_Dist_AU_RP = math.sqrt(math.pow(Sep_Dist_AU_RP_x, 2.0)+math.pow(Sep_Dist_AU_RP_y,2.0))
    if verbose:
        print 'Separation Distance in the reference plane [AU] = ',Sep_Dist_AU_RP
    
    # calculate corrected angle between Line of Nodes and M2 (ie. AN-focus-M2)
    # the first atan will produce a negative if x or is negative, so must correct that below
    # must take all 4 quadrants into account.
    ang_LN_to_M2_corrPOSNEG = math.degrees(math.atan(Sep_Dist_AU_RP_y/Sep_Dist_AU_RP_x))
    if verbose:
        print ' angle LN to M2 POSNEG [deg] = '+str(ang_LN_to_M2_corrPOSNEG)
    
    if (Sep_Dist_AU_RP_x>=0.0)and (Sep_Dist_AU_RP_y>=0.0):
        # quadrant 1
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG 
    elif (Sep_Dist_AU_RP_x<0.0)and (Sep_Dist_AU_RP_y>=0.0):
        # quadrant 2
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +180.0
    elif (Sep_Dist_AU_RP_x<0.0)and (Sep_Dist_AU_RP_y<0.0):
        # quadrant 3
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +180.0
    elif (Sep_Dist_AU_RP_x>=0.0)and (Sep_Dist_AU_RP_y<0.0):
        # quadrant 4
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +360.0
    else:
        print ' The value of '+str(ang_LN_to_M2_corrPOSNEG)+' is apparently out of ranges'
        
    if verbose:
        print 'uncorrected angle between line of nodes and M2 [deg] = ', ang_LN_to_M2_deg
        print 'Corrected angle between Line of nodes and M2 [deg] = ',ang_LN_to_M2_corr
    
    ## calculate measured Position Angle 
    # it is measured from M1, only if the value is over 360deg do we need to fix it
    totalAngle_RP = longAN_deg + ang_LN_to_M2_corr
    if totalAngle_RP>360.0:
        PA_deg_RP_model = totalAngle_RP -360.0
    else:
        PA_deg_RP_model = totalAngle_RP
#    if Sep_Dist_AU_RP_x>=0.0:
#        PA_deg_RP_model = longAN_deg + 180.0 + ang_LN_to_M2_corr
#        #print '** corrected M2 found in 2nd or 3rd quadrant' #$$$$$$$$
#    elif Sep_Dist_AU_RP_x<0.0:
#        PA_deg_RP_model = longAN_deg + ang_LN_to_M2_corr
#        #print '** corrected M2 found in 1st or 4th quadrant' #$$$$$$$$$$$
    if verbose:
        print 'Position Angle in reference plane [deg] = ',PA_deg_RP_model
        
    ## calculate measured Separation Angle (model)
    SA_arcsec_RP_model = Sep_Dist_AU_RP/Sys_Dist_PC # PC->AU and rad->arcsec cancel
    if verbose:
        print 'Separation Angle measured (model) [arcsec] = ',SA_arcsec_RP_model
    
    
    return(PA_deg_RP_model, SA_arcsec_RP_model)
    
def orbElementPlotter(plotFileTitle, longAN_degs, es, Ts, a1s, a2s, periods, inclination_degs, argPeri_degs,\
                  Sep_Dists2, thetas2, Es2, Ms2, ns2, showPlots=False):
    
    ########################## INS #####################
    # plot in's on first figure
    plt.figure(1,figsize=(20,7) ,dpi=300)
    plt.subplot(241)
    (n,bins,patches) = plt.hist(longAN_degs,bins=50)
    plt.xlabel('longAN_degs')
    plt.ylabel('Probability')
    plt.subplot(242)
    (n,bins,patches) = plt.hist(es,bins=50)
    plt.xlabel('es')
    plt.title(plotFileTitle+' *Inputs*')
    plt.subplot(243)
    Tsmin = int(np.min(Ts))
    Tsmax = int(np.max(Ts))
    TsNEW = Ts
    for i in range(0,len(Ts)):
        TsNEW[i]=(TsNEW[i]-Tsmin)/365
    (n,bins,patches) = plt.hist(TsNEW,bins=25)
    plt.xlabel('Ts [years since '+str(Tsmin)+']')
    # only plot the a1s if they aren't all 0.0
    if np.mean(a1s)>0.0:
        plt.subplot(244)
        (n,bins,patches) = plt.hist(a1s,bins=50)
        plt.xlabel('a1s')
    plt.subplot(245)
    (n,bins,patches) = plt.hist(a2s,bins=50)
    plt.ylabel('Probability')
    plt.xlabel('a2s')
    plt.subplot(246)
    (n,bins,patches) = plt.hist(periods,bins=50)
    plt.xlabel('periods')
    plt.subplot(247)
    (n,bins,patches) = plt.hist(inclination_degs,bins=50)
    plt.xlabel('inclination_degs')
    plt.subplot(248)
    (n,bins,patches) = plt.hist(argPeri_degs,bins=50)
    plt.xlabel('argPeri_degs')
    
    plt.savefig('../figures/'+plotFileTitle+'INS.png', dpi=300, orientation='landscape')
   
    ################# OUTS ####################
    numSamples = np.shape(Es2[:][:])[0]
    numEpochs = np.shape(Es2[:][:])[1]
    
    for epoch in range(0,numEpochs):
        # reshape the output matrices
        Sep_Dists2SingleEpoch = []
        thetas2SingleEpoch = []
        Es2SingleEpoch = []
        Ms2SingleEpoch = []
        ns2SingleEpoch = []
        for sample in range(0,numSamples):
            Sep_Dists2SingleEpoch.append(Sep_Dists2[sample][epoch])
            thetas2SingleEpoch.append(thetas2[sample][epoch])
            Es2SingleEpoch.append(Es2[sample][epoch])
            Ms2SingleEpoch.append(Ms2[sample][epoch])
            ns2SingleEpoch.append(ns2[sample][epoch])
        # plot OUT's on first figure
        plt.figure(epoch+2, figsize=(20,7) ,dpi=300)
        plt.subplot(231)
        (n,bins,patches) = plt.hist(Sep_Dists2SingleEpoch,bins=50)
        plt.xlabel('Sep_Dists [AU]')
        plt.ylabel('Probability')
        plt.subplot(232)
        (n,bins,patches) = plt.hist(thetas2SingleEpoch,bins=50)
        plt.xlabel('thetas')
        plt.title(plotFileTitle+' *Outputs of simulator for epoch '+str(epoch+1)+'*')
        plt.subplot(233)
        (n,bins,patches) = plt.hist(Es2SingleEpoch,bins=50)
        plt.xlabel('Es2')
        plt.subplot(234)
        (n,bins,patches) = plt.hist(Ms2SingleEpoch,bins=50)
        plt.xlabel('Ms2')
        plt.ylabel('Probability')
        plt.subplot(235)
        (n,bins,patches) = plt.hist(ns2SingleEpoch,bins=50)
        plt.xlabel('ns2')
        filename = '../figures/'+plotFileTitle+'OUTSepoch'+str(epoch+1)
        plt.savefig('%s.png' % filename, dpi=300, orientation='landscape')
    
    # Show the plots on the screen if asked to    
    if showPlots:
        plt.show()

def orbitEllipsePlotter3(longAN_deg, e, period, inc, argPeri_deg, a, xLabel='delta RA [mas]', yLabel='delta Dec [mas]', plotFilename='', \
                                                        xLim=False, yLim=False, show=True):
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
    
    ## Measured data (import from paramSettigsDict)
    from paramSettingsDict import paramSettingsDict
    # DI data
    SAs = paramSettingsDict['Data']['SAs']                                   
    PAs = paramSettingsDict['Data']['PAs']
    SAerrors = paramSettingsDict['Data']['SAerrors']                                   
    PAerrors = paramSettingsDict['Data']['PAerrors']
    # General System Data
    sys_dist = paramSettingsDict['Data']['sysDist']#15.62
    
    
    if plotFilename!='':
        if plotFilename[-4:]!='.png':
            plotFilename = plotFilename+'.png'     
    # create a plot title from the filename
    plotFileTitle = os.path.basename(plotFilename).split('.')[0]
    
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
            orbitCalculator(0.0, sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    Xstart = SAstart*math.sin(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    Ystart = -SAstart*math.cos(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    startStar = star((asConversion/1000.0)*0.6*a[0], Xstart, Ystart, color='green', N=5, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, a1, a2) =\
            orbitCalculator(((period[0]/4.0)*365.25), sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    XoneQuarter = SAstart*math.sin(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    YoneQuarter = -SAstart*math.cos(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    oneQuarterStar = star((asConversion/1000.0)*0.6*a[0], XoneQuarter, YoneQuarter, color='blue', N=5, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAhalf, PAhalf, a1, a2) =\
            orbitCalculator(((period[0]/2.0)*365.25), sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    Xhalf = SAhalf*math.sin(math.radians(PAhalf))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    Yhalf = -SAhalf*math.cos(math.radians(PAhalf))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    halfStar = star((asConversion/1000.0)*0.6*a[0], Xhalf, Yhalf, color='yellow', N=5, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, a1, a2) =\
            orbitCalculator((3.0*(period[0]/4.0)*365.25), sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    XthreeQuarter = SAstart*math.sin(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    YthreeQuarter = -SAstart*math.cos(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    threeQuarterStar = star((asConversion/1000.0)*0.6*a[0], XthreeQuarter, YthreeQuarter, color='orange', N=5, thin = 0.5)
    
    (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SAstart, PAstart, a1, a2) =\
            orbitCalculator(((period[0])*365.25), sys_dist, inc[0], longAN_deg[0], e[0], 0.0, period[0], argPeri_deg[0], a[0],\
                                                        Mass1=1, Mass2=1, verbose=False)
    Xend = SAstart*math.sin(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
    Yend = -SAstart*math.cos(math.radians(PAstart))*asConversion#*sys_dist  # This will be in [mas] instead of [AU]
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
            print 'PROBLEM: xLim is not of type tuple'
    if not yLim:
        min = np.min(ellipseYs2[:])
        max = np.max(ellipseYs2[:])
        Range = abs(max)+abs(min)
        yLim = (min-abs(Range*0.05),max+abs(Range*0.05))
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
    (errorBoxes, m2starPolygons) = starAndErrorPolys(SAs,SAerrors,PAs,PAerrors,asConversion, main.transData)
    
    # draw the error boxes for the companion start locations in the data
    for errorBox in errorBoxes:
        main.add_patch(errorBox)
        
    # Draw red star patches for each of the companion's locations from the data
    for star2 in m2starPolygons:
        main.add_patch(star2)
        
    # add a legend
    #main.legend(('longAN_deg = '+str(longAN_deg),'e = '+str(e), 'period = '+str(period), 'inc = '+str(inc), 'argPeri_deg = '+str(argPeri_deg), 'a = '+str(a)), loc=0, markerscale=0.0000000000000001)
    legndStr = 'longAN_deg = '+str(longAN_deg[0])+'\ne = '+str(e[0])+'\nperiod = '+str(period[0])+'\ninc = '+str(inc[0])+'\nargPeri_deg = '+str(argPeri_deg[0])+'\na = '+str(a[0])
    main.text(xLim[0]+abs(xLim[0]*0.02),abs(yLim[1]*0.2),legndStr,ha='left')
    
    # save plot to file
    if plotFilename!='':
        plt.savefig(plotFilename, dpi=300, orientation='landscape')
    
    # show plot
    if show:
        plt.show()
    
    # close figure/plot
    plt.close()    

def rvPlotter(longAN_deg, e, T, period, inc, argPeri_deg, a, RV_origin_vel_0_proposed=0, RV_origin_vel_1_proposed=0,\
              plotFilename='', show=True):
    """
    create a plot for the RV data and a fit line from the best orbit data params.
    """
       
    ## Measured data (import from paramSettigsDict)
    from paramSettingsDict import paramSettingsDict
    RV_epochs = paramSettingsDict['Data']['RV_epochs']
    RVerrors = paramSettingsDict['Data']['RVerrors']
    RVs = paramSettingsDict['Data']['RVs']
    Mass1 = paramSettingsDict['Data']['M1']#1 
    Mass2 = paramSettingsDict['Data']['M2']#1
    
    if plotFilename!='':
        if plotFilename[-4:]!='.png':
            plotFilename = plotFilename+'.png'     
    # create a plot title from the filename
    plotFileTitle = os.path.basename(plotFilename).split('.')[0]
    
    RVsIN = RVs
    RVsOUT = RVs
    #print 'RV for dataset 0'
    dataset0 = []
    dataset1 = []
    for epoch in range(0,len(RVs[0])):
        vel = RVsIN[0][epoch]+RV_origin_vel_0_proposed
        dataset0.append(vel)
        #print str(RVsIN[0][epoch]) +' + '+str(RV_origin_vel_0_proposed) +" = "+str(vel)
    #print 'RV for dataset 1'
    for epoch in range(0,len(RVs[1])):
        vel = RVsIN[1][epoch]+RV_origin_vel_1_proposed
        dataset1.append(vel)
        #print str(RVsIN[1][epoch]) +' + '+str(RV_origin_vel_1_proposed) +" = "+str(vel)
        
    RVs = [dataset0,dataset1]
    
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
    if type(T)!=list:
       T = [T]
    
    # data for planet's orbit that are needed to calculate proper chiSquare for RV fit
    # all below values for the planet orbit and RV data are from Butler2006.
    planet_K = 461.1 #[m/s]
    planet_p = 3.31246   #[days]
    #planet_p_years = planet_p/365.25
    planet_e = 0.023    
    planet_argPeri = 188.0   #[deg]
    planet_To = 2446957.8   #[JD]

#    RVouts2 = []
#    times2 = []
#    for orb in range(0,len(longAN_deg)):
#        RVouts = []
#        times = []
#        numSteps = 5000.0
#        periodIncrement = (period[orb]*365.25)/numSteps
#        t = 1.0 
#        (a_total, a1, a2, p) = semiMajorConverter(Mass1, Mass2, a_total=a[orb],a1=0.0,a2=0.0)
#        for step in range(0,int(numSteps)):
#            
#            t = t + periodIncrement +T[orb]
#            # calculate the velocity residual due to the companion star
#            v_r_c = vrCalculatorStar2(t,e[orb],T[orb],period[orb],argPeri_deg[orb],a1,i=inc[orb], K=False, verbose=False)
#            # calculate the velocity residual due to the planet around primary
#            v_r_p = vrCalculatorPlanet(t, planet_e, planet_To, planet_p, planet_argPeri, Mass1, M2SineI=False, K=planet_K, verbose=False)
#            RV = v_r_c+v_r_p
#            RVouts.append(RV)
#            times.append(t)
#        RVouts2.append(RVouts)
#        times2.append(times)
    
    residuals3 = []
    for orb in range(0,len(longAN_deg)):
        residuals2 = []
        chiSquaredTot2 = 0
        numEpochs_RV = 0
        chi_squared_RV_reducedCur2 = []
        for dataset in range(0,len(RVs)):
            chiSquaredTot = 0
            epochs = RV_epochs[dataset]
            rvs = RVs[dataset]
            errors = RVerrors[dataset]
            residuals = []
            (a_total, a1, a2, p) = semiMajorConverter(Mass1, Mass2, a_total=a[orb],a1=0.0,a2=0.0)
            print '\For dataset '+str(dataset)
            for epoch in range(0,len(epochs)):
                # calculate the velocity residual due to the companion star
                (v_r_c,K) = vrCalculatorStar2(epochs[epoch],e[orb],T[orb],period[orb],argPeri_deg[orb],a1,i=inc[orb], K=False, verbose=False)
                # calculate the velocity residual due to the planet around primary
                (v_r_p,K) = vrCalculatorPlanet(epochs[epoch], planet_e, planet_To, planet_p, planet_argPeri, Mass1, M2SineI=False, K=planet_K, verbose=False)
                RV = v_r_c+v_r_p-rvs[epoch]#diffCalc(v_r_c+v_r_p,rvs[epoch])
                print 'RV = ('+str(v_r_c)+' + '+str(v_r_p)+') - '+str(rvs[epoch])+' = '+str(RV)
                residuals.append(RV)
                chiSquaredTot = chiSquaredTot+chiSquaredCalc(rvs[epoch], errors[epoch]+15, v_r_c+v_r_p)
            chi_squared_RV_reducedCur = (1.0/((1.0*len(RVs[dataset]))-8.0))*chiSquaredTot
            chi_squared_RV_reducedCur2.append(chi_squared_RV_reducedCur)
            chiSquaredTot2 = chiSquaredTot2 + chiSquaredTot
            numEpochs_RV_curr = len(RV_epochs[dataset])
            numEpochs_RV = numEpochs_RV+numEpochs_RV_curr
            residuals2.append(residuals)
        residuals3.append(residuals2) 
        
        
    # convert to a reduced chiSquared
    chi_squared_RV_reduced = (1.0/((1.0*numEpochs_RV)-8.0))*chiSquaredTot
    ## Create figure, axes and start plotting/drawing everything
    fig = plt.figure(1,figsize=(10,10))
    main = fig.add_subplot(111)
    main.set_title(plotFileTitle)
    
    ## draw RV orbits
#    for orb in range(0,len(longAN_deg)):
#        main.plot(times2[orb],RVouts2[orb])
    
    ## plot RV data
#    for dataset in range(0,len(RVs)):
#        main.scatter(RV_epochs[dataset], RVs[dataset], s=5, c='red', edgecolors='none')
    xmin = 1e10
    xmax = 0
    ymin = 1e10
    ymax = -1e10
    for orb in range(0,len(longAN_deg)):
        meanMedianChiSquaredSTR = ''
        residuals3[orb]
        for dataset in range(0,len(RVs)):
            meanMedianChiSquaredSTR = meanMedianChiSquaredSTR+'\nmean of residuals for dataset '+str(dataset)+' is = '+str(np.mean(residuals3[orb][dataset]))
            meanMedianChiSquaredSTR = meanMedianChiSquaredSTR+'\nmedian of residuals for dataset '+str(dataset)+' is = '+str(np.median(residuals3[orb][dataset]))
            meanMedianChiSquaredSTR = meanMedianChiSquaredSTR+'\nchiSquared_reduced for dataset '+str(dataset)+' is = '+str(chi_squared_RV_reducedCur2[dataset])
            main.scatter(RV_epochs[dataset],residuals3[orb][dataset],s=7)
            main.plot([np.min(RV_epochs[dataset]),np.max(RV_epochs[dataset])],[0,0])
            # set up x,y limits
            if np.min(RV_epochs[dataset])<xmin:
                xmin = np.min(RV_epochs[dataset])
            if np.max(RV_epochs[dataset])>xmax:
                xmax = np.max(RV_epochs[dataset])
            if np.min(residuals3[orb][dataset])<ymin:
                ymin = np.min(residuals3[orb][dataset])
            if np.max(residuals3[orb][dataset])>ymax:
                ymax = np.max(residuals3[orb][dataset])
                
    meanMedianChiSquaredSTR = meanMedianChiSquaredSTR+'\nchiSquared_reduced for ALL data is = '+str(chi_squared_RV_reduced)
    print meanMedianChiSquaredSTR
    xLim =[xmin,xmax]
    yLim = [ymin,ymax]
#    print 'legend bottom left corner at [ '+str(xLim[0]+abs(xLim[0]*0.02))+' , '+str(abs(yLim[1]*0.2))+' ]'
#    print 'xLim = '+repr(xLim)
    # add a legend
    #main.legend(('longAN_deg = '+str(longAN_deg),'e = '+str(e), 'period = '+str(period), 'inc = '+str(inc), 'argPeri_deg = '+str(argPeri_deg), 'a = '+str(a)), loc=0, markerscale=0.0000000000000001)
    legndStr = 'longAN_deg = '+str(longAN_deg[0])+'\ne = '+str(e[0])+'\nperiod = '+str(period[0])+'\ninc = '+str(inc[0])+\
                '\nargPeri_deg = '+str(argPeri_deg[0])+'\na = '+str(a[0])+meanMedianChiSquaredSTR
    main.text(xLim[1]-abs(xLim[1]*0.002),abs(yLim[1]*0.5),legndStr,ha='left')
    
    # save plot to file
    if plotFilename!='':
        plt.savefig(plotFilename, dpi=300, orientation='landscape')
    
    # show plot
    if show:
        plt.show()
        
    plt.close()
        
def rvFitPlotter1Body(longAN_deg, e, T, period, inc, argPeri_deg, a=0.0, RVoffsets=[0],\
              plotFilename='', show=True):
    """
    create a plot for the RV data and a fit line from the best orbit data params.
    """
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
    if type(T)!=list:
       T = [T]
       
    ## Measured data (import from paramSettigsDict)
    from paramSettingsDict import paramSettingsDict
    RV_epochs = paramSettingsDict['Data']['RV_epochs']
    RVerrors = paramSettingsDict['Data']['RVerrors']
    RVs = paramSettingsDict['Data']['RVs']
    Mass1 = paramSettingsDict['Data']['M1']#1 
    Mass2 = paramSettingsDict['Data']['M2']#1
    
    if plotFilename!='':
        if plotFilename[-4:]!='.png':
            plotFilename = plotFilename+'.png'     
    # create a plot title from the filename
    plotFileTitle = os.path.basename(plotFilename).split('.')[0]
    
    RV_epochsIN2 = RV_epochs
    RV_epochsOUT3 = []
    RV_epochsOUT2 = []
    RV_epochsOUT = []
    print 'About to convert epochs to all be in a single orbit'
    for orb in range(0,len(T)):
        print 'Working on input Time of last Periasis ',orb
        Tcurr = T[orb]
        p_curr_days = period[orb]*365.242 # converted as t and T are in JD
        for RV_epochsIN in RV_epochsIN2:
            RV_epochsOUT = []
            for epoch in RV_epochsIN:
                if epoch<Tcurr:
                    if (Tcurr-epoch)>p_curr_days:
                        numPeriodsAwayI = int((Tcurr-epoch)/p_curr_days)
                    else:
                        numPeriodsAwayI = 0
                    epochNEW = epoch+((numPeriodsAwayI+1)*p_curr_days)-Tcurr
                    RV_epochsOUT.append(epochNEW)
                    
                if epoch>Tcurr:
                    if (epoch-Tcurr)>p_curr_days:
                        numPeriodsAwayI = int((epoch-Tcurr)/p_curr_days)
                    else:
                        numPeriodsAwayI = 0
                    epochNEW = epoch-((numPeriodsAwayI-1)*p_curr_days)-Tcurr
                    RV_epochsOUT.append(epochNEW)
                if False:
                    print 'RV epoch updated from '+str(epoch)+" to "+str(epochNEW)+", as numPeriodsAwayI = "+str(numPeriodsAwayI)
                    print "T = "+str(Tcurr)
                    print "Period [days] = "+str(p_curr_days)
            RV_epochsOUT2.append(RV_epochsOUT)
        RV_epochsOUT3.append(RV_epochsOUT2)
    
    if False:
        print '\nOriginal RV epochs array = '+repr(RV_epochsIN2)
        print '\nOUTPUT RV epochs array = '+repr(RV_epochsOUT3)
        
    RVsIN = RVs
    RVsOUT = []
    for dataset in range(0,len(RVsIN)):
        RVsUSE=[]
        for epoch in range(0,len(RVsIN[dataset])):
            vel = RVsIN[dataset][epoch]#-RV_origin_vel_0_proposed
            RVsUSE.append(vel)
        RVsOUT.append(RVsUSE)
    
    #Get residual values for RV scatter plot
    residuals3 = []
    for orb in range(0,len(longAN_deg)):
        residuals2 = []
        chiSquaredTot2 = 0
        numEpochs_RV = 0
        chi_squared_RV_reducedCur2 = []
        for dataset in range(0,len(RVsOUT)):
            chiSquaredTot = 0
            epochs = RV_epochsOUT3[orb][dataset]
            rvs = RVsOUT[dataset]
            if False:
                print 'len(RVsOUT) = ',len(RVsOUT)
                print 'len(RV_epochsOUT3[orb]) = ',len(RV_epochsOUT3[orb])
                print 'len(RVsOUT[dataset]) = ',len(RVsOUT[dataset])
                print 'len(RV_epochsOUT3[orb][dataset]) = ',len(RV_epochsOUT3[orb][dataset])
            errors = RVerrors[dataset]
            
            if len(RVoffsets)>=(dataset+1):
                offset = RVoffsets[dataset]
            else:
                offset = 0
            
            residuals = []
            (a_total, a1, a2, p) = semiMajorConverter(Mass1, Mass2, a_total=a[orb],a1=0.0,a2=0.0, period=period[orb],verbose=True)
            #print '\nFor dataset '+str(dataset)
            for epoch in range(0,len(epochs)):
                print '\nWorking on epoch ',epoch
                # calculate the velocity residual due to the companion star
                (v_r_c,K) = vrCalculatorStar2(epochs[epoch],e[orb],0,period[orb],argPeri_deg[orb],a1,i=inc[orb], K=False, verbose=False)               
                RV = v_r_c-rvs[epoch]+offset
                print 'RV = ('+str(v_r_c)+') - '+str(rvs[epoch])+" + "+str(offset)+' = '+str(RV)
                residuals.append(RV)
                chiSquaredCurr = chiSquaredCalc(rvs[epoch]-offset, errors[epoch], v_r_c)
                chiSquaredTot = chiSquaredTot+chiSquaredCurr
                print 'chiSquaredCurr = ',chiSquaredCurr
                print 'chiSquaredTot = ',chiSquaredTot
            chi_squared_RV_reducedCur = (1.0/((1.0*len(RVs[dataset]))-6.0))*chiSquaredTot
            chi_squared_RV_reducedCur2.append(chi_squared_RV_reducedCur)
            chiSquaredTot2 = chiSquaredTot2 + chiSquaredTot
            print 'chiSquaredTot2 = ',chiSquaredTot2
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
    for orb in range(0,len(e)):
        orbitVRs = []
        numSteps = 5000.0
        periodIncrement = (period[orb]*365.25)/numSteps
        t = 1.0 
        times = []
        (a_total, a1, a2, p) = semiMajorConverter(Mass1, Mass2, a_total=a[orb],a1=0.0,a2=0.0, period=period[orb])
        for step in range(0,int(numSteps)):
            Tcurr = 0
            t = Tcurr + t + periodIncrement
            times.append(t)
            # calculate the velocity residual due to the companion 
            (v_r_c,K) = vrCalculatorStar2(t,e[orb],Tcurr,period[orb],argPeri_deg[orb],a1,i=inc[orb], K=False, verbose=False)
            orbitVRs.append(v_r_c)
        print 'Orbit '+str(orb)+" had a K = "+str(K)
        orbitVRs2.append(orbitVRs)
        
    
    ## Create figure, axes and start plotting/drawing everything
    fig = plt.figure(1,figsize=(10,10))
    residualsPlot = fig.add_subplot(212)
    
    ## draw RV orbits
#    for orb in range(0,len(longAN_deg)):
#        main.plot(times2[orb],RVouts2[orb])
    
    ## plot RV data
#    for dataset in range(0,len(RVs)):
#        main.scatter(RV_epochs[dataset], RVs[dataset], s=5, c='red', edgecolors='none')
    xmin = 1e10
    xmax = 0
    ymin = 1e10
    ymax = -1e10
    for orb in range(0,len(longAN_deg)):
        meanMedianChiSquaredSTR = ''
        residuals3[orb]
        for dataset in range(0,len(RVs)):
            meanMedianChiSquaredSTR = meanMedianChiSquaredSTR+'\nmean of residuals for dataset '+str(dataset)+' is = '+str(np.mean(residuals3[orb][dataset]))
            meanMedianChiSquaredSTR = meanMedianChiSquaredSTR+'\nmedian of residuals for dataset '+str(dataset)+' is = '+str(np.median(residuals3[orb][dataset]))
            meanMedianChiSquaredSTR = meanMedianChiSquaredSTR+'\nchiSquared_reduced for dataset '+str(dataset)+' is = '+str(chi_squared_RV_reducedCur2[dataset])
            residualsPlot.scatter(RV_epochsOUT3[orb][dataset],residuals3[orb][dataset],s=7)
            residualsPlot.plot([np.min(RV_epochsOUT3[orb][dataset]),np.max(RV_epochsOUT3[orb][dataset])],[0,0])
            # set up x,y limits
            if np.min(RV_epochsOUT3[orb][dataset])<xmin:
                xmin = np.min(RV_epochs[dataset])
            if np.max(RV_epochsOUT3[orb][dataset])>xmax:
                xmax = np.max(RV_epochsOUT3[orb][dataset])
            if np.min(residuals3[orb][dataset])<ymin:
                ymin = np.min(residuals3[orb][dataset])
            if np.max(residuals3[orb][dataset])>ymax:
                ymax = np.max(residuals3[orb][dataset])
                
    meanMedianChiSquaredSTR = meanMedianChiSquaredSTR+'\nchiSquared_reduced for ALL data is = '+str(chi_squared_RV_reduced)
    print meanMedianChiSquaredSTR
    xLim =[xmin,xmax]
    xrange = xmax-xmin
    xLim2 = (xmin-abs(xrange*0.05), xmax+abs(xrange*0.05))
    yLim = [ymin,ymax]
    yrange = ymax-ymin
    yLim2 = (ymin-abs(yrange*0.05), ymax+abs(yrange*0.05))
#    print 'legend bottom left corner at [ '+str(xLim[0]+abs(xLim[0]*0.02))+' , '+str(abs(yLim[1]*0.2))+' ]'
#    print 'xLim = '+repr(xLim)
    # add a legend
    #main.legend(('longAN_deg = '+str(longAN_deg),'e = '+str(e), 'period = '+str(period), 'inc = '+str(inc), 'argPeri_deg = '+str(argPeri_deg), 'a = '+str(a)), loc=0, markerscale=0.0000000000000001)
    legndStr = 'longAN_deg = '+str(longAN_deg[0])+'\ne = '+str(e[0])+'\nperiod = '+str(period[0])+'\ninc = '+str(inc[0])+\
                '\nargPeri_deg = '+str(argPeri_deg[0])+'\na = '+str(a[0])+meanMedianChiSquaredSTR
    addLegend=False
    if addLegend:
        residualsPlot.text(xLim[1]-abs(xLim[1]*0.002),abs(yLim[1]*0.5),legndStr,ha='left')
    
    ## make plot of fit to data
    fitPlot = fig.add_subplot(211)
    fitXmin = np.min(times[0])
    if xmin<fitXmin:
        fitXmin = xmin
    fitXmax = np.max(times[0])
    if xmax>fitXmax:
        fitXmax = xmax
    fitYmin = np.min(orbitVRs2[0])
    if ymin<fitYmin:
        fitYmin = ymin
    fitYmax = np.max(orbitVRs2[0])
    if ymax>fitYmax:
        fitYmax = ymax
    fitXrange = fitXmax-fitXmin
    fitXLim2 = (fitXmin-abs(fitXrange*0.05), fitXmax+abs(fitXrange*0.05))
    fitYrange = fitYmax-fitYmin
    fitYLim2 = (fitYmin-abs(fitYrange*0.05), fitYmax+abs(fitYrange*0.05))
    fitPlot.axes.set_xlim(fitXLim2)
    fitPlot.axes.set_ylim(fitYLim2)
    for orb in range(0,len(e)):
        fitPlot.plot(times,orbitVRs2[orb])
        for dataset in range(0,len(RVs)):
            fitPlot.scatter(RV_epochsOUT3[orb][dataset],RVs[dataset],s=7)
            
    fitPlot.set_title(plotFileTitle)
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
        
def orbElementTester(longAN_deg, e, period, inclination_deg, argPeri_deg, a_total, T):
    """
    This version is meant to be used on only the DI data and the original orbitCalculator, 
    not orbitCalculator2.
    
    Use orbElementTester2 to test for 3D data.
    """
    #### HD130948BC
#    SA_arcsec_measured_REALs = [0.0993, 0.0946, 0.0961, 0.0978, 0.0568, 0.0564, 0.0573, 0.1117, 0.1091, 0.0979, 0.0722, 0.0574, 0.0573, 0.0584, 0.0517]
#    SA_mean_errors = [0.0017, 0.0017, 0.0023, 0.0023, 0.0009, 0.0012, 0.0009, 0.0015, 0.0011, 0.0009, 0.0006, 0.0005, 0.0009, 0.0010, 0.0006] # mean simply implies the mean of the + and - uncertainties                                         
#    PA_deg_measured_REALs = [307.9,306.9,308.9,307.0,144.9,146.4,148.6,132.6,132.28,130.73,127.6,124.7,124.6,124.1,123.9]
#    PA_mean_errors = [1.1,1.0,1.6,1.4,0.5,0.6,0.6,0.4,0.13,0.17,0.3,0.4,0.6,0.7,0.5]# mean simply implies the mean of the + and - uncertainties 
#    epochs = [2452519.702,2452519.706,2452519.709,2452519.719,2453425.152,2453425.219,2453425.282,2454127.085,2454185.028,2454306.737,2454481.192,2454555.068,2454555.09,2454555.157,2454584.106]
#    Sys_Dist_PC = 18.17
    Mass1 = 1
    Mass2 = 1
    #### Tau Boo
#    SA_arcsec_measured_REALs = [3.401, 2.865, 2.82, 2.71, 1.93]
#    SA_mean_errors = [0.06, 0.04, 0.05, 0.03, 0.03]                                         
#    PA_deg_measured_REALs = [20.65, 30.85, 33.2, 31.3, 53.1]
#    PA_mean_errors = [0.5, 0.4, 1.0, 1.0, 0.6]
#    epochs = [2448337.5, 2451143.5, 2451711.5, 2451945.5, 2455590.5]
#    Sys_Dist_PC = 15.62
    
#    SA_arcsec_measured_REALs = [10.26, 9.03, 8.91, 8.71, 8.9, 8.89, 8.72, 8.41, 8.45, 8.54, 8.66, 7.28, 7.35, 6.96, 7.06, 6.4, 6.15, 6.09, 5.87, 5.4, 5.32, 5.14, 4.94, 5.17,      3.401, 2.865, 2.82, 2.71, 2.181, 1.9421, 1.9324]
#    SA_mean_errors = [1.46, 0.12, 0.04, 0.2, 0.08, 0.1, 0.05, 0.22, 0.13, 0.04, 0.22, 0.89, 0.37, 0.66, 0.37, 0.37, 0.45, 0.42, 0.25, 0.33, 0.22, 0.07, 0.03, 0.38,         0.06, 0.04, 0.05, 0.03, 0.02, 0.01, 0.01]                                         
#    PA_deg_measured_REALs = [347.8, 348.7, 350.8, 351.9, 354.3, 352.4, 353.6, 356, 355.1, 354.8, 355.8, 356.8, 359.1, 359.5, 1.7, 3.5, 5.1, 5.1, 6, 6.8, 7.8, 11.4, 8, 13.5,               20.65, 30.85, 33.2, 31.3,46.26, 52.65, 53.68]
#    PA_mean_errors = [0.34, 1.88, 0.63, 0.17, 1.79, 0.38, 0.31, 1.43, 0.19, 0.16, 0.45, 0.21, 0.1, 0.09, 1.51, 0.93, 1.34, 1.04, 0.28, 0.75, 0.68, 1.2, 3.69, 1.08,           0.5, 0.4, 1.0, 1.0, 0.42, 0.16, 0.07]
#    epochs = [2396557.5, 2403065.5, 2405359.5, 2407064.5, 2408219.5, 24708927.5, 2411846.5, 2413509.5, 2414351.5, 2414465.5, 2415404.5, 2419263.5, 2423443.5, 2424260.5, 2425728.5, 2429753.5, 2431535.5, 2431959.5, 2434136.5, 2436256.5, 2437227.5, 2438866.5, 2440132.5, 2440707.5,      2448337.5, 2451143.5, 2451711.5, 2451945.5, 2454579.5, 2455589.5, 2455704.5]
#    Sys_Dist_PC = 15.62
#
    SA_arcsec_measured_REALs = [10.26, 9.03, 8.91, 8.71, 8.9, 8.89, 8.72, 8.41, 8.45, 8.54, 8.66, 7.28, 7.35, 6.96, 7.06, 6.4, 6.15, 6.09, 5.87, 5.4, 5.32, 5.14, 4.94, 5.17    ]
    SA_mean_errors = [1.46, 0.12, 0.04, 0.2, 0.08, 0.1, 0.05, 0.22, 0.13, 0.04, 0.22, 0.89, 0.37, 0.66, 0.37, 0.37, 0.45, 0.42, 0.25, 0.33, 0.22, 0.07, 0.03, 0.38      ]                                         
    PA_deg_measured_REALs = [347.8, 348.7, 350.8, 351.9, 354.3, 352.4, 353.6, 356, 355.1, 354.8, 355.8, 356.8, 359.1, 359.5, 1.7, 3.5, 5.1, 5.1, 6, 6.8, 7.8, 11.4, 8, 13.5           ]
    PA_mean_errors = [0.34, 1.88, 0.63, 0.17, 1.79, 0.38, 0.31, 1.43, 0.19, 0.16, 0.45, 0.21, 0.1, 0.09, 1.51, 0.93, 1.34, 1.04, 0.28, 0.75, 0.68, 1.2, 3.69, 1.08         ]
    epochs = [2396557.5, 2403065.5, 2405359.5, 2407064.5, 2408219.5, 24708927.5, 2411846.5, 2413509.5, 2414351.5, 2414465.5, 2415404.5, 2419263.5, 2423443.5, 2424260.5, 2425728.5, 2429753.5, 2431535.5, 2431959.5, 2434136.5, 2436256.5, 2437227.5, 2438866.5, 2440132.5, 2440707.5    ]
    Sys_Dist_PC = 15.62
#    Mass1 = 1 #1.3
#    Mass2 = 1 #0.4
    
    print "\n Orbital Parameters Input:",
    print "Inclination [deg] = ",inclination_deg
    print "Longitude of the Ascending Node [deg] = ",longAN_deg
    print "Eccentricity = ",e
    print "Time of last periapsis [jd] = ",T
    print "Period [years] = ",period
    #print "Tmin = "+str(epochs[0]-period*365.0)
    #print "Tmax = "+str(epochs[0])
    print "Argument of perigee [deg] = ",argPeri_deg
    print "Total Semi-major axis [AU] = ",a_total
    
    (chi_squared_total_cur, ns, Ms, Es, thetas, Sep_Dists, SA_deg_measured_models, PA_deg_measured_models, a1s, a2s) =\
        multiEpochOrbCalc(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors,\
                       epochs, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total,\
                       Mass1=Mass1, Mass2=Mass2, verbose=False )
        
    numEpochs = len(epochs)
    
    print " uncorrected Chi Squared = ",chi_squared_total_cur
    # normalize the chiSquared to give the reduce chiSquared value
    one_over_nu = (1.0/((2.0*(numEpochs))-7.0))
    print "one_over_nu = ",one_over_nu
    reduced_chi_squared_total_cur = one_over_nu*chi_squared_total_cur
    print "\n Resulting total reduced chiSquared is ",reduced_chi_squared_total_cur
    
    ## TEMP for comparing to C++ code
    print "\n measured values:"
    print "SA: ",SA_arcsec_measured_REALs[0]
    print "SA error: ",SA_mean_errors[0]
    print "PA: ",PA_deg_measured_REALs[0]
    print "PA error: ",PA_mean_errors[0]
    print "epoch: ",epochs[0]
    print "sys dist: ",Sys_Dist_PC    
    
    print "\noutputs:\n"
    print "a1s: ",a1s[0]
    print "a2s: ",a2s[0]
    print "reduce chi squared: ",reduced_chi_squared_total_cur
    print "ns: ",ns[0]
    print "Ms: ",Ms[0]
    print "Es: ",Es[0]
    print "thetas: ",thetas[0]
    print "sep dists: ",Sep_Dists[0]
    print "SA model: ",SA_deg_measured_models[0]
    print "PA model: ",PA_deg_measured_models[0]

def orbElementTester2(longAN_deg, e, T, period, inclination_deg, argPeri_deg, a_total, RV_origin_vel_0=0, RV_origin_vel_1=0):
    """
    This version is meant to be used on 3D (DI and RV) data located in the paramSettingsDict
    and the orbitCalculator2.
    
    Use orbElementTester to test for DI only data.
    """
    
    ## Measured data (import from paramSettigsDict)
    from paramSettingsDict import paramSettingsDict
    # DI data
    SA_arcsec_measured_REALs = paramSettingsDict['Data']['SAs']
    SA_mean_errors = paramSettingsDict['Data']['SAerrors']                                         
    PA_deg_measured_REALs = paramSettingsDict['Data']['PAs']
    PA_mean_errors = paramSettingsDict['Data']['PAerrors']
    epochs = paramSettingsDict['Data']['DI_epochs']
    # RV data
    RV_epochs = paramSettingsDict['Data']['RV_epochs']
    RVerrors = paramSettingsDict['Data']['RVerrors']
    RVs = paramSettingsDict['Data']['RVs']
    # General System Data
    Sys_Dist_PC = paramSettingsDict['Data']['sysDist']#15.62
    Sys_Dist_Error = paramSettingsDict['Data']['sysDistError']#
    Mass1 = paramSettingsDict['Data']['M1']#1 
    Mass1_Error = paramSettingsDict['Data']['M1Error']#
    Mass2 = paramSettingsDict['Data']['M2']#1
    Mass2_Error = paramSettingsDict['Data']['M2Error']#
    
    print "\n Orbital Parameters Input:",
    print "Inclination [deg] = ",inclination_deg
    print "Longitude of the Ascending Node [deg] = ",longAN_deg
    print "Eccentricity = ",e
    print "Time of last periapsis [jd] = ",T
    print "Period [years] = ",period
    #print "Tmin = "+str(epochs[0]-period*365.0)
    #print "Tmax = "+str(epochs[0])
    print "Argument of perigee [deg] = ",argPeri_deg
    print "Total Semi-major axis [AU] = ",a_total
    
    numEpochs_DI = len(epochs)
           
    ## send random parameters along with known ones to multi-epoch orbit calculator
    (chi_squared_total_cur, ns, Ms, Es, thetas, Sep_Dists, SA_deg_measured_models, PA_deg_measured_models, a1s, a2s) =\
        multiEpochOrbCalc(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors,\
                       epochs, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total,\
                       Mass1=Mass1, Mass2=Mass2, verbose=False)

    # calculate a total from a1 plus a2
    a_total_proposed = np.mean(a1s)+np.mean(a2s)
    
    chi_squared_DI = chi_squared_total_cur
    # convert to a reduced chiSquared
    chi_squared_DI_reduced = (1.0/((2.0*(numEpochs_DI))-6.0))*chi_squared_DI
    #chi_squared_DI_reduced=0.0 #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    RVs_proposed = RVs
    for epoch in range(0,len(RVs_proposed[0])):
        RVs_proposed[0][epoch] = RVs[0][epoch]+RV_origin_vel_0
    for epoch in range(0,len(RVs_proposed[1])):
        RVs_proposed[1][epoch] = RVs[1][epoch]+RV_origin_vel_1
    
    # data for planet's orbit that are needed to calculate proper chiSquare for RV fit
    # all below values for the planet orbit and RV data are from Butler2006.
    planet_K = 461.1 #[m/s]
    planet_p = 3.31246   #[days]
    #planet_p_years = planet_p/365.25
    planet_e = 0.023    
    planet_argPeri = 188.0   #[deg]
    planet_To = 2446957.8   #[JD]
    sigma_jitter = 15.0    #[m/s]
    
    # calculate a1 from a_total and masses
    (a_total, a1, a2, p) = semiMajorConverter(Mass1, Mass2, a_total=a_total_proposed,a1=0.0,a2=0.0)
    chi_squared_RV = 0.0
    numEpochs_RV = 0.0
    for RVdataSet in range(0,len(RVs)):            
#        chi_squared_RV_curr = rv2bodyCalculator3(RV_epochs[RVdataSet], RVs_proposed[RVdataSet], RVerrors[RVdataSet], sigma_jitter, Mass1, Mass2,  \
#                                                    inclination_deg, period, e, T,  \
#                                                    argPeri_deg, planet_K, planet_p, planet_e, planet_argPeri, planet_To)
        chi_squared_RV_curr = rv2bodyCalculator4(RV_epochs[RVdataSet], RVs_proposed[RVdataSet], RVerrors[RVdataSet], sigma_jitter, Mass1, np.mean(a1s), inclination_deg, \
                                                        period, e, T, argPeri_deg, \
                                                            planet_K, planet_p, planet_e, planet_argPeri, planet_To, verbose=False)
        chi_squared_RV = chi_squared_RV+chi_squared_RV_curr
        numEpochs_RV_curr = len(RV_epochs[RVdataSet])
        numEpochs_RV = numEpochs_RV+numEpochs_RV_curr
        
    # convert to a reduced chiSquared
    chi_squared_RV_reduced = (1.0/((1.0*numEpochs_RV)-8.0))*chi_squared_RV
    
    # sum non-reduced chiSquareds
    chi_squared_total = chi_squared_DI + chi_squared_RV #+ gauss_chi_squared_total
    
    # Total Reduced Chi Squared
    totalNumDataPoints = 2.0*numEpochs_DI + 1.0*numEpochs_RV +3.0 # three are for the gauss params
    totalNumModelParams = 6.0  # 6 proposed standard model params
    TotalNu = (totalNumDataPoints-totalNumModelParams)
    chi_squared_total_reduced = (1.0/TotalNu)*chi_squared_total
    
    print "\nOutputs:\n"
    print 'reduced chisquared DI = ',chi_squared_DI_reduced
    print 'reduced chiSquared RV = ',chi_squared_RV_reduced
    print 'Total reduced chiSquared = ',chi_squared_total_reduced 
#    print "a1s: "+repr(a1s)
#    print "a2s: "+repr(a2s)
#    print "ns: "+repr(ns)
#    print "Ms: "+repr(Ms)
#    print "Es: "+repr(Es)
#    print "thetas: "+repr(thetas)
#    print "sep dists: ",repr(Sep_Dists)
#    print "SA model: "+repr(SA_deg_measured_models)
#    print "PA model: "+repr(PA_deg_measured_models)
    
    
def orbElementTester3(longAN_deg, e, T, period, inclination_deg, argPeri_deg, a_total,RV_origin_vel_0_proposed=0,RV_origin_vel_1_proposed=0):
    """
    This version is meant to be used on 3D (DI and RV) data located in the paramSettingsDict
    and the orbitCalculator2.
    
    Use orbElementTester to test for DI only data.
    """
    
    ## Measured data (import from paramSettigsDict)
    from paramSettingsDict import paramSettingsDict
    # DI data
    SA_arcsec_measured_REALs = paramSettingsDict['Data']['SAs']
    SA_mean_errors = paramSettingsDict['Data']['SAerrors']                                         
    PA_deg_measured_REALs = paramSettingsDict['Data']['PAs']
    PA_mean_errors = paramSettingsDict['Data']['PAerrors']
    epochs = paramSettingsDict['Data']['DI_epochs']
    # RV data
    RV_epochs = paramSettingsDict['Data']['RV_epochs']
    RVerrors = paramSettingsDict['Data']['RVerrors']
    RVs = paramSettingsDict['Data']['RVs']
    # General System Data
    Sys_Dist_PC = paramSettingsDict['Data']['sysDist']#15.62
    Sys_Dist_Error = paramSettingsDict['Data']['sysDistError']#
    Mass1 = paramSettingsDict['Data']['M1']#1 
    Mass1_Error = paramSettingsDict['Data']['M1Error']#
    Mass2 = paramSettingsDict['Data']['M2']#1
    Mass2_Error = paramSettingsDict['Data']['M2Error']#
    
    print "\n Orbital Parameters Input:",
    print "Inclination [deg] = ",inclination_deg
    print "Longitude of the Ascending Node [deg] = ",longAN_deg
    print "Eccentricity = ",e
    print "Time of last periapsis [jd] = ",T
    print "Period [years] = ",period
    #print "Tmin = "+str(epochs[0]-period*365.0)
    #print "Tmax = "+str(epochs[0])
    print "Argument of perigee [deg] = ",argPeri_deg
    print "Total Semi-major axis [AU] = ",a_total
    
    numEpochs_DI = len(epochs)
           
    ## send random parameters along with known ones to multi-epoch orbit calculator
    (chi_squared_total, ns, Ms, Es, thetas, xs, ys, a1s, a2s) =\
        multiEpochOrbCalc3(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors,\
                       epochs, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total,\
                       Mass1=Mass1, Mass2=Mass2, verbose=False)

    # calculate a total from a1 plus a2
    a_total_proposed = np.mean(a1s)+np.mean(a2s)
    
    chi_squared_DI = chi_squared_total
    # convert to a reduced chiSquared
    chi_squared_DI_reduced = (1.0/((2.0*(numEpochs_DI))-6.0))*chi_squared_DI
    #chi_squared_DI_reduced=0.0 #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    # data for planet's orbit that are needed to calculate proper chiSquare for RV fit
    # all below values for the planet orbit and RV data are from Butler2006.
    planet_K = 461.1 #[m/s]
    planet_p = 3.31246   #[days]
    #planet_p_years = planet_p/365.25
    planet_e = 0.023    
    planet_argPeri = 188.0   #[deg]
    planet_To = 2446957.8   #[JD]
    sigma_jitter = 15.0    #[m/s]
    
    # calculate a1 from a_total and masses
    (a_total, a1, a2, p) = semiMajorConverter(Mass1, Mass2, a_total=a_total_proposed,a1=0.0,a2=0.0)
    chi_squared_RV = 0.0
    numEpochs_RV = 0.0
    
    
#    for epoch in range(0,len(RVs[0])):
#        RVs[0][epoch] = RVs[0][epoch]+RV_origin_vel_0_proposed
#    for epoch in range(0,len(RVs[1])):
#        RVs[1][epoch] = RVs[1][epoch]+RV_origin_vel_1_proposed
              
    RVsIN = RVs
    RVsOUT = RVs
    #print 'RV for dataset 0'
    dataset0 = []
    dataset1 = []
    for epoch in range(0,len(RVs[0])):
        vel = RVsIN[0][epoch]+RV_origin_vel_0_proposed
        dataset0.append(vel)
        #print str(RVsIN[0][epoch]) +' + '+str(RV_origin_vel_0_proposed) +" = "+str(vel)
    #print 'RV for dataset 1'
    for epoch in range(0,len(RVs[1])):
        vel = RVsIN[1][epoch]+RV_origin_vel_1_proposed
        dataset1.append(vel)
        #print str(RVsIN[1][epoch]) +' + '+str(RV_origin_vel_1_proposed) +" = "+str(vel)
        
    RVs = [dataset0,dataset1]
      
    for RVdataSet in range(1,len(RVs)):
        chi_squared_RV_curr = rv2bodyCalculator3(RV_epochs[RVdataSet], RVs[RVdataSet], RVerrors[RVdataSet], sigma_jitter, Mass1, Mass2,  \
                                                    inclination_deg, period, e, T,  \
                                                    argPeri_deg, planet_K, planet_p, planet_e, planet_argPeri, planet_To, verbose=True)
#        chi_squared_RV_curr = rv2bodyCalculator4(RV_epochs[RVdataSet], RVs[RVdataSet], RVerrors[RVdataSet], sigma_jitter, Mass1, np.mean(a1s), inclination_deg, \
#                                                        period, e, T, argPeri_deg, \
#                                                            planet_K, planet_p, planet_e, planet_argPeri, planet_To, verbose=True)
        chi_squared_RV = chi_squared_RV+chi_squared_RV_curr
        numEpochs_RV_curr = len(RV_epochs[RVdataSet])
        numEpochs_RV = numEpochs_RV+numEpochs_RV_curr
        
    # convert to a reduced chiSquared
    chi_squared_RV_reduced = (1.0/((1.0*numEpochs_RV)-8.0))*chi_squared_RV
    
    # sum non-reduced chiSquareds
    chi_squared_total = chi_squared_DI + chi_squared_RV #+ gauss_chi_squared_total
    
    # Total Reduced Chi Squared
    totalNumDataPoints = 2.0*numEpochs_DI + 1.0*numEpochs_RV #+3.0 # three are for the gauss params
    totalNumModelParams = 8.0  # 6 proposed standard model params
    TotalNu = (totalNumDataPoints-8.0)
    chi_squared_total_reduced = (1.0/TotalNu)*chi_squared_total
    
    print "\nOutputs:\n"
    print 'reduced chisquared DI = ',chi_squared_DI_reduced
    print 'reduced chiSquared RV = ',chi_squared_RV_reduced
    print 'Total reduced chiSquared = ',chi_squared_total_reduced 
#    print "a1s: "+repr(a1s)
#    print "a2s: "+repr(a2s)
#    print "ns: "+repr(ns)
#    print "Ms: "+repr(Ms)
#    print "Es: "+repr(Es)
#    print "thetas: "+repr(thetas)
#    print "sep dists: ",repr(Sep_Dists)
#    print "SA model: "+repr(SA_deg_measured_models)
#    print "PA model: "+repr(PA_deg_measured_models)    

def orbElementTester4(longAN_deg, e, T, period, inclination_deg, argPeri_deg, a_total,RV_origin_vel_0_proposed=0,RV_origin_vel_1_proposed=0):
    """
    This version is meant to be used on RV data of single planet located in the paramSettingsDict.
    
    Use orbElementTester to test for DI only data.
    """
    
    ## Measured data (import from paramSettigsDict)
    from paramSettingsDict import paramSettingsDict
    # DI data
    SA_arcsec_measured_REALs = paramSettingsDict['Data']['SAs']
    SA_mean_errors = paramSettingsDict['Data']['SAerrors']                                         
    PA_deg_measured_REALs = paramSettingsDict['Data']['PAs']
    PA_mean_errors = paramSettingsDict['Data']['PAerrors']
    epochs = paramSettingsDict['Data']['DI_epochs']
    # RV data
    RV_epochs = paramSettingsDict['Data']['RV_epochs']
    RVerrors = paramSettingsDict['Data']['RVerrors']
    RVs = paramSettingsDict['Data']['RVs']
    # General System Data
    Sys_Dist_PC = paramSettingsDict['Data']['sysDist']#15.62
    Sys_Dist_Error = paramSettingsDict['Data']['sysDistError']#
    Mass1 = paramSettingsDict['Data']['M1']#1 
    Mass1_Error = paramSettingsDict['Data']['M1Error']#
    Mass2 = paramSettingsDict['Data']['M2']#1
    Mass2_Error = paramSettingsDict['Data']['M2Error']#
    
    print "\n Orbital Parameters Input:\n",
    print "Inclination [deg] = ",inclination_deg
    print "Longitude of the Ascending Node [deg] = ",longAN_deg
    print "Eccentricity = ",e
    print "Time of last periapsis [jd] = ",T
    print "Period [years] = ",period
    #print "Tmin = "+str(epochs[0]-period*365.0)
    #print "Tmax = "+str(epochs[0])
    print "Argument of perigee [deg] = ",argPeri_deg
    print "Total Semi-major axis [AU] = ",a_total
    
    
    RV_epochsIN2 = RV_epochs
    RV_epochsOUT2 = []
    RV_epochsOUT = []
    print 'About to convert epochs to all be in a single orbit'
    Tcurr = T
    p_curr_days = period*365.242
    for RV_epochsIN in RV_epochsIN2:
        RV_epochsOUT = []
        for epoch in RV_epochsIN:
            if epoch<Tcurr:
                if (Tcurr-epoch)>p_curr_days:
                    numPeriodsAwayI = int((Tcurr-epoch)/p_curr_days)
                else:
                    numPeriodsAwayI = 0
                epochNEW = epoch+((numPeriodsAwayI+1)*p_curr_days)-Tcurr
                RV_epochsOUT.append(epochNEW)
                
            if epoch>Tcurr:
                if (epoch-Tcurr)>p_curr_days:
                    numPeriodsAwayI = int((epoch-Tcurr)/p_curr_days)
                else:
                    numPeriodsAwayI = 0
                epochNEW = epoch-((numPeriodsAwayI-1)*p_curr_days)-Tcurr
                RV_epochsOUT.append(epochNEW)
            if False:
                print 'RV epoch updated from '+str(epoch)+" to "+str(epochNEW)+", as numPeriodsAwayI = "+str(numPeriodsAwayI)
                print "T = "+str(Tcurr)
                print "Period [days] = "+str(p_curr_days)
        RV_epochsOUT2.append(RV_epochsOUT)
    
    if False:
        print '\nOriginal RV epochs array = '+repr(RV_epochsIN2)
        print '\nOUTPUT RV epochs array = '+repr(RV_epochsOUT3)
    
    
    numEpochs_DI = len(epochs)
           
    ## send random parameters along with known ones to multi-epoch orbit calculator
    (chi_squared_total, ns, Ms, Es, thetas, xs, ys, a1s, a2s) =\
        multiEpochOrbCalc3(SA_arcsec_measured_REALs, SA_mean_errors, PA_deg_measured_REALs, PA_mean_errors,\
                       epochs, Sys_Dist_PC, inclination_deg, longAN_deg, e, T, period, argPeri_deg, a_total,\
                       Mass1=Mass1, Mass2=Mass2, verbose=False)

    # calculate a total from a1 plus a2
    a_total_proposed = np.mean(a1s)+np.mean(a2s)
    
    chi_squared_DI = chi_squared_total
    # convert to a reduced chiSquared
    chi_squared_DI_reduced = (1.0/((2.0*(numEpochs_DI))-6.0))*chi_squared_DI
    #chi_squared_DI_reduced=0.0 #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    print 'Semi-major 1 [AU] = ',np.mean(a1s)
    
    # calculate a1 from a_total and masses
    #(a_total, a1, a2, p) = semiMajorConverter(Mass1, Mass2, a_total=a_total_proposed,a1=0.0,a2=0.0)
    chi_squared_RV = 0.0
    numEpochs_RV = 0.0
    
    
#    for epoch in range(0,len(RVs[0])):
#        RVs[0][epoch] = RVs[0][epoch]+RV_origin_vel_0_proposed
#    for epoch in range(0,len(RVs[1])):
#        RVs[1][epoch] = RVs[1][epoch]+RV_origin_vel_1_proposed
              
    RVsIN = RVs
    RVsOUT = []
    for dataset in range(0,len(RVsIN)):
        RVsUSE=[]
        for epoch in range(0,len(RVsIN[dataset])):
            vel = RVsIN[dataset][epoch]#-RV_origin_vel_0_proposed
            RVsUSE.append(vel)
        RVsOUT.append(RVsUSE)
      
    sigma_jitter = 0.0# 3.0 #15.0    #[m/s]
    
    for RVdataSet in range(1,len(RVs)):        
        chi_squared_RV_curr = rv1bodyCalculator(RV_epochsOUT2[RVdataSet], RVsOUT[RVdataSet], RVerrors[RVdataSet], sigma_jitter,\
                                                         inclination_deg, period, e, T, argPeri_deg, np.mean(a1s), verbose=True)
        
        chi_squared_RV = chi_squared_RV+chi_squared_RV_curr
        numEpochs_RV_curr = len(RV_epochs[RVdataSet])
        numEpochs_RV = numEpochs_RV+numEpochs_RV_curr
        
    # convert to a reduced chiSquared
    chi_squared_RV_reduced = (1.0/((1.0*numEpochs_RV)-8.0))*chi_squared_RV
    
    # sum non-reduced chiSquareds
    chi_squared_total = chi_squared_DI + chi_squared_RV #+ gauss_chi_squared_total
    
    # Total Reduced Chi Squared
    totalNumDataPoints = 2.0*numEpochs_DI + 1.0*numEpochs_RV #+3.0 # three are for the gauss params
    totalNumModelParams = 8.0  # 6 proposed standard model params
    TotalNu = (totalNumDataPoints-8.0)
    chi_squared_total_reduced = (1.0/TotalNu)*chi_squared_total
    
    print "\nOutputs:\n"
    print 'reduced chisquared DI = ',chi_squared_DI_reduced
    print 'reduced chiSquared RV = ',chi_squared_RV_reduced
    print 'Total reduced chiSquared = ',chi_squared_total_reduced 
    
def orbElementTester3TH_I(longAN_deg, e, T, period, inclination_deg, argPeri_deg, a_total):
    """
    This version is meant to be used on 3D (DI and RV) data located in the paramSettingsDict
    and the orbitCalculator2.
    
    Use orbElementTester to test for DI only data.
    """
    
    ## Measured data (import from paramSettigsDict)
    from paramSettingsDict import paramSettingsDict
    # DI data
    SA_arcsec_measured_REALs = paramSettingsDict['Data']['SAs']
    SA_mean_errors = paramSettingsDict['Data']['SAerrors']                                         
    PA_deg_measured_REALs = paramSettingsDict['Data']['PAs']
    PA_mean_errors = paramSettingsDict['Data']['PAerrors']
    epochs = paramSettingsDict['Data']['DI_epochs']
    # RV data
    RV_epochs = paramSettingsDict['Data']['RV_epochs']
    RVerrors = paramSettingsDict['Data']['RVerrors']
    RVs = paramSettingsDict['Data']['RVs']
    # General System Data
    Sys_Dist_PC = paramSettingsDict['Data']['sysDist']#15.62
    Sys_Dist_Error = paramSettingsDict['Data']['sysDistError']#
    Mass1 = paramSettingsDict['Data']['M1']#1 
    Mass1_Error = paramSettingsDict['Data']['M1Error']#
    Mass2 = paramSettingsDict['Data']['M2']#1
    Mass2_Error = paramSettingsDict['Data']['M2Error']#
    
    
    
    numEpochs_DI = len(epochs)
    
    # calculate ABCFG from input params
    (A,B,C,F,G) = ABCFG_values(a_total/Sys_Dist_PC, math.radians(argPeri_deg), math.radians(longAN_deg), math.radians(inclination_deg))
    
    ## send random parameters along with known ones to multi-epoch orbit calculator
    (a_arcsec, argPeri_rad, longAN_rad, inclination_rad, period_out, chiSquared_total) =\
        multiEpochOrbCalcTH_I3(e, T, A, B, C, F, G, Mass1, Mass2, Sys_Dist_PC, \
                                                    SA_arcsec_measured_REALs, PA_deg_measured_REALs, epochs, SA_mean_errors, PA_mean_errors)

    # calculate a total from a1 plus a2
    a_total_proposed = a_arcsec*Sys_Dist_PC
    argPeri_deg_out = math.degrees(argPeri_rad)
    longAN_deg_out = math.degrees(longAN_rad)
    inclination_deg_out = math.degrees(inclination_rad)
    
    print "\n Orbital Parameters Input:",
    print "Inclination [deg] = ",inclination_deg
    print "Longitude of the Ascending Node [deg] = ",longAN_deg
    print "Eccentricity = ",e
    print "Time of last periapsis [jd] = ",T
    print "Period [years] = ",period
    print "Argument of perigee [deg] = ",argPeri_deg
    print "Total Semi-major axis [AU] = ",a_total
    
    print "\n Orbital Parameters Output:",
    print "Inclination [deg] = ",inclination_deg_out
    print "Longitude of the Ascending Node [deg] = ",longAN_deg_out
    print "Period [years] = ",period_out
    print "Argument of perigee [deg] = ",argPeri_deg_out
    print "Total Semi-major axis [AU] = ",a_total_proposed
    
    
    chi_squared_DI = chiSquared_total
    # convert to a reduced chiSquared
    chi_squared_DI_reduced = (1.0/((2.0*(numEpochs_DI))-6.0))*chi_squared_DI
    #chi_squared_DI_reduced=0.0 #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
    # data for planet's orbit that are needed to calculate proper chiSquare for RV fit
    # all below values for the planet orbit and RV data are from Butler2006.
    planet_K = 461.1 #[m/s]
    planet_p = 3.31246   #[days]
    #planet_p_years = planet_p/365.25
    planet_e = 0.023    
    planet_argPeri = 188.0   #[deg]
    planet_To = 2446957.8   #[JD]
    sigma_jitter = 15.0    #[m/s]
    
    # calculate a1 from a_total and masses
    (a_total, a1, a2, p) = semiMajorConverter(Mass1, Mass2, a_total=a_total_proposed,a1=0.0,a2=0.0)
    chi_squared_RV = 0.0
    numEpochs_RV = 0.0
    for RVdataSet in range(0,len(RVs)):
        chi_squared_RV_curr = rv2bodyCalculator3(RV_epochs[RVdataSet], RVs[RVdataSet], RVerrors[RVdataSet], sigma_jitter, Mass1, Mass2,  \
                                                    inclination_deg, period, e, T,  \
                                                    argPeri_deg, planet_K, planet_p, planet_e, planet_argPeri, planet_To)
#        chi_squared_RV = rv2bodyCalculator4(RV_epochs, RVs, RVerrors, sigma_jitter, Mass1_proposed, np.mean(a1s), inclination_deg_proposed, \
#                                                        period_proposed, e_proposed, T_proposed, argPeri_deg_proposed, \
#                                                            planet_K, planet_p, planet_e, planet_argPeri, planet_To, verbose=False)
        chi_squared_RV = chi_squared_RV+chi_squared_RV_curr
        numEpochs_RV_curr = len(RV_epochs[RVdataSet])
        numEpochs_RV = numEpochs_RV+numEpochs_RV_curr
        
    # convert to a reduced chiSquared
    chi_squared_RV_reduced = (1.0/((1.0*numEpochs_RV)-8.0))*chi_squared_RV
    
    # sum non-reduced chiSquareds
    chi_squared_total = chi_squared_DI + chi_squared_RV #+ gauss_chi_squared_total
    
    # Total Reduced Chi Squared
    totalNumDataPoints = 2.0*numEpochs_DI + 1.0*numEpochs_RV +3.0 # three are for the gauss params
    totalNumModelParams = 6.0  # 6 proposed standard model params
    TotalNu = (totalNumDataPoints-totalNumModelParams)
    chi_squared_total_reduced = (1.0/TotalNu)*chi_squared_total
    
    print "\nOutputs:\n"
    print 'reduced chisquared DI = ',chi_squared_DI_reduced
    print 'reduced chiSquared RV = ',chi_squared_RV_reduced
    print 'Total reduced chiSquared = ',chi_squared_total_reduced 
#    print "a1s: "+repr(a1s)
#    print "a2s: "+repr(a2s)
#    print "ns: "+repr(ns)
#    print "Ms: "+repr(Ms)
#    print "Es: "+repr(Es)
#    print "thetas: "+repr(thetas)
#    print "sep dists: ",repr(Sep_Dists)
#    print "SA model: "+repr(SA_deg_measured_models)
#    print "PA model: "+repr(PA_deg_measured_models)  
        
        
        