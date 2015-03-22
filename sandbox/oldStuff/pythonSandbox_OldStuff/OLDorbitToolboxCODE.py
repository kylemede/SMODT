def dataReadAndPlotNEW(filename, plotFilename='', numEpochs=1):
    """
    This function will load the file into memory, then make a summary plot of the 
    parameters. 
    THIS VERSION IS FOR THE 'NEW' DATA FILES WITH ONE FILE ONLY PER SIMULATION.
    longAN [deg]      e [N/A]       To [julian date]  period [yrs]   inclination [deg]   argPeri [deg]   a_total [AU]    chiSquared    timesBeenHere
    
    @param filename:         = name of the file where the data is.
    @type filename:          = string
    @param plotFilename:     = file name for the output plot
    @type plotFilename:      = string, including directory if not cwd
    @param numEpochs:        = number of epochs of observation/data
    @type numEpochs:         = integer
    """
    
    # check if the passed in value for filenameROOT includes '.txt'
    if '.txt' in filename:
        print 'PROBLEM: there is .txt in the filename and there should not be'
    
    # check if the passed in value for plotFilename includes '.png'
    if '.png' not in plotFilename:
        plotFilename = plotFilename+'.png'
    
    print '\nLoading data into memory '
    
    longAN_degsAll = dataReaderNEW(filename,0)
    esAll = dataReaderNEW(filename,1)
    TsAll = dataReaderNEW(filename,2)
    periodsAll = dataReaderNEW(filename,3)
    inclination_degsAll = dataReaderNEW(filename,4)
    argPeri_degsAll = dataReaderNEW(filename,5)
    a1Means = [0.0,0.0] # this tricks the summaryPlotter to only plot a2s, ie a_totals
    a2Means = dataReaderNEW(filename,6)
    chiSquareds = dataReaderNEW(filename,7)
    
    # Convert the chiSqureds to reduced chiSquareds for the plotter IF numEpochs>1
    # At the same time the lowest chiSquared value and location in the array will be found
    if numEpochs>1:
        reducedChiSquareds = []
        reducedChiSquaredMin = 1000.0
        bestOrbit = 0
        #print len(chiSquareds)
        for orbit in range(0,len(chiSquareds)):
            reducedChiSquared = (1.0/((2.0*(numEpochs))-7.0))*chiSquareds[orbit]
            reducedChiSquareds.append(reducedChiSquared)
            if reducedChiSquared<reducedChiSquaredMin:
                reducedChiSquaredMin = reducedChiSquared
                #print 'new chiSquare min = ',reducedChiSquaredMin
                bestOrbit = orbit
    else:
        # Find best orbit params
        chiSquaredMin = 1000.0
        bestOrbit = 0
        for orbit in range(0,len(chiSquareds)):
            if chiSquareds[orbit]<chiSquaredMin:
                chiSquaredMin = chiSquareds[orbit]
                bestOrbit = orbit
   
    print 'Best orbit found:'
    if numEpochs==1:
        print "Non-reduced chiSquaredMin = ",chiSquaredMin
    else:
        print "Reduced chiSquaredMin = ",reducedChiSquaredMin
    print "LongAN = ",longAN_degsAll[bestOrbit]
    print "e = ",esAll[bestOrbit]
    print "period = ",periodsAll[bestOrbit]
    print "inclination = ",inclination_degsAll[bestOrbit]
    print "argPeri = ",argPeri_degsAll[bestOrbit]
    print "a_total = ",a2Means[bestOrbit]
    print "To = ",TsAll[bestOrbit]
    
    print 'Finished loading data, starting to make summary plot'
    if True:
        if numEpochs>1:
            summaryPlotter(plotFilename, reducedChiSquareds, inclination_degsAll, longAN_degsAll, argPeri_degsAll, \
                   esAll, periodsAll, TsAll, a1Means, a2Means, confLevels=True, weight=True, normed=True,\
                    showPlots=False, save=True, verbose=True)
        else:
                summaryPlotter(plotFilename, chiSquareds, inclination_degsAll, longAN_degsAll, argPeri_degsAll, \
                   esAll, periodsAll, TsAll, a1Means, a2Means, confLevels=True, weight=True, normed=True,\
                    showPlots=False, save=True, verbose=True)
    
    print 'Summary plot complete and written to disk '+plotFilename

def dataReadAndPlotNEW2(filename, plotFilename='', numEpochs=1):
    """
    This function will load the file into memory, then make a summary plot of the 
    parameters. 
    THIS VERSION IS FOR THE 'NEW' DATA FILES WITH ONE FILE ONLY PER SIMULATION.
    longAN [deg]      e [N/A]       To [julian date]  period [yrs]   inclination [deg]   argPeri [deg]   a_total [AU]    chiSquared    timesBeenHere
    
    @param filename:         = name of the file where the data is.
    @type filename:          = string
    @param plotFilename:     = file name for the output plot
    @type plotFilename:      = string, including directory if not cwd
    @param numEpochs:        = number of epochs of observation/data
    @type numEpochs:         = integer
    """
    weight = False
    confLevels = False
    normed = False
    
    # check if the passed in value for filename includes '.txt'
    if '.txt' not in filename:
        filename = filename+'.txt'
        print ".txt was added to the end of the filename as it wasn't found to be in the original version"
    
    # check if the passed in value for plotFilename includes '.png'
    if '.png' not in plotFilename:
        plotFilename = plotFilename+'.png'
    
    # create a plot title from the filename
    plotFileTitle = os.path.basename(plotFilename).split('.')[0]
    
    print '\nLoading data into memory '
    
    chiSquareds = dataReaderNEW(filename,7)
    print 'chiSquareds loaded'
    gc.collect()
    print 'garbage collected'
    # Convert the chiSqureds to reduced chiSquareds for the plotter IF numEpochs>1
    # At the same time the lowest chiSquared value and location in the array will be found
    if numEpochs>1:
        print 'starting to calculated reduced chiSquareds'
        reducedChiSquareds = []
        reducedChiSquaredMin = 1000.0
        bestOrbit = 0
        #print len(chiSquareds)
        for orbit in range(0,len(chiSquareds)):
            reducedChiSquared = (1.0/((2.0*(numEpochs))-7.0))*chiSquareds[orbit]
            reducedChiSquareds.append(reducedChiSquared)
            if reducedChiSquared<reducedChiSquaredMin:
                reducedChiSquaredMin = reducedChiSquared
                #print 'new chiSquare min = ',reducedChiSquaredMin
                bestOrbit = orbit
        # Set the reduced chiSquareds to be the chiSquaredsUSE for the histograms
        del chiSquareds
        print 'chiSquareds deleted'
        chiSquaredsUSE = reducedChiSquareds
        del reducedChiSquareds
        print 'Reduced chiSquareds calculated and deleted'
    else:
        # Find best orbit params
        chiSquaredMin = 1000.0
        bestOrbit = 0
        for orbit in range(0,len(chiSquareds)):
            if chiSquareds[orbit]<chiSquaredMin:
                chiSquaredMin = chiSquareds[orbit]
                bestOrbit = orbit
        # Set the non-reduced chiSquareds to be the chiSquaredsUSE for the histograms
        chiSquaredsUSE = chiSquareds
        del chiSquareds
    
    gc.collect()
    print '2nd garbage collected'
    
    print '\nloading inclinations'   
    inclination_degsAlls = dataReaderNEW(filename,4)
    print 'finished loading inclinations'
    
    if len(inclination_degsAlls)>0:
        fig = plt.figure(1, figsize=(25,22) ,dpi=300)
        
        plot = fig.add_subplot(331)
        data = inclination_degsAlls
        xlabel = 'Inclination [deg]'
        plot.axes.set_ylabel('Probability')
        plot = histConverter(chiSquaredsUSE, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
        incMedian = np.median(inclination_degsAlls)
        numOrbits = len(inclination_degsAlls)
        incBest = inclination_degsAlls[bestOrbit]
        del inclination_degsAlls, data
        print "Done plotting inclination_degsAlls"
        
        plot = fig.add_subplot(332)
        longAN_degsAlls = dataReaderNEW(filename,0)
        data = longAN_degsAlls
        xlabel = 'Longitude of Ascending Node [deg]'
        plt.title(plotFileTitle+'   TOTAL Summary')
        plot = histConverter(chiSquaredsUSE, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
        longANMedian = np.median(longAN_degsAlls)
        longANBest = longAN_degsAlls[bestOrbit]
        del longAN_degsAlls, data
        print "Done plotting longAN_degsAlls"
        
        plot = fig.add_subplot(333)
        argPeri_degsAlls = dataReaderNEW(filename,5)
        data = argPeri_degsAlls
        xlabel = 'Argument of Perigie [deg]'
        plot = histConverter(chiSquaredsUSE, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
        argPeriMedian = np.median(argPeri_degsAlls)
        argPeriBest = argPeri_degsAlls[bestOrbit]
        del argPeri_degsAlls, data
        print "Done plotting argPeri_degsAlls"
        
        plot = fig.add_subplot(334)
        esAlls = dataReaderNEW(filename,1)
        data = esAlls
        xlabel = 'e'
        plot.axes.set_ylabel('Probability')
        plot = histConverter(chiSquaredsUSE, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
        eMedian = np.median(esAlls)
        eBest = esAlls[bestOrbit]
        del esAlls, data
        print "Done plotting esAlls"
        
        plot = fig.add_subplot(335)
        periodsAlls = dataReaderNEW(filename,3)
        data = periodsAlls
        xlabel = 'Period [Years]'
        plot = histConverter(chiSquaredsUSE, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
        periodMedian = np.median(periodsAlls)
        periodBest = periodsAlls[bestOrbit]
        del periodsAlls, data
        print "Done plotting periodsAlls"
        
        plot = fig.add_subplot(336)
        TsAlls = dataReaderNEW(filename,2)
        Tsmin = int(np.min(TsAlls))
        Tsmax = int(np.max(TsAlls))
        TsNEW = TsAlls
        for i in range(0,len(TsNEW)):
            TsNEW[i]=(TsNEW[i]-Tsmin)/365
        data = TsNEW
        xlabel = 'Time of last Periapsis [years since '+str(Tsmin)+']'
        plot = histConverter(chiSquaredsUSE, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
        TMedian = np.median(TsAlls)
        TBest = TsAlls[bestOrbit]
        del TsAlls, data
        print "Done plotting TsAlls"
        
        a2Means = dataReaderNEW(filename,6)
        plot = fig.add_subplot(337)
        data = a2Means
        xlabel = 'Semi-major Axis [AU]'
        plot.axes.set_ylabel('Probability')
        plot = histConverter(chiSquaredsUSE, data, plot, xlabel, confLevels=confLevels, weight=weight, normed=normed)
        a2Median = np.median(a2Means)
        a2Best = a2Means[bestOrbit]
        del a2Means, data
        print "Done plotting Semi-major Axis"
   
        # Save file 
        plt.savefig(plotFilename, dpi=300, orientation='landscape')
        print 'Summary figure saved to '+plotFilename
        # print median values of data for reference 
        print 'Finished making Total Summary Plot'
        print '#'*50
        print str(numOrbits)+' orbits used'
        print 'median inclination_degsAlls =  ',incMedian
        print 'median longAN_degsAlls = ', longANMedian
        print 'median argPeri_degsAlls = ', argPeriMedian
        print 'median esAlls = ',eMedian
        print 'median periodsAlls =  ',periodMedian
        print 'median TsAlls [jd] = ', TMedian
        print 'median a = ',a2Median
        print '#'*50
            
        # print the values for the best orbit
        print '*'*50
        print 'Best orbit found:'
        if numEpochs==1:
            print "Non-reduced chiSquaredMin = ",chiSquaredMin
        else:
            print "Reduced chiSquaredMin = ",reducedChiSquaredMin
        print "LongAN = ",longANBest
        print "e = ",eBest
        print "period = ",periodBest
        print "inclination = ",incBest
        print "argPeri = ",argPeriBest
        print "a_total = ",a2Best
        print "To = ",TBest
        print '*'*50
        
    else:
        print 'No orbits to plot in Total Summary !!'    
    
    plt.close()
    
    print 'Summary plot complete and written to disk '+plotFilename
    
    
def correlationLengthCalc(filename):
    """
    Finds the largest Correlation Length of a data set following Tegmark2004.
    ie. the number of steps necessary for the variance to equal 1/2 the total chain variance.
    NOTE: taking square root of both sides reduces this relation to: 1/4 std of total chain
        which avoids wasting processor power in squaring the stds to give the variance.
    """
    colHeaderStrs = ["longAN [deg]","e [N/A]","To [julian date]","period [yrs]","inclination [deg]","argPeri [deg]","a_total [AU]"]
    corrLengths = []
    maxCorrLength = 0
    maxCorrLengthCol = 0
    maxEntireVar = 0
    
    for col in range(0,7):
        data = dataReaderNEW(filename,col)
        # instantiate a list of all stds so far
        stds = []
        
        for i in range(1,int(0.4*len(data))):
            stds.append(np.std(data[:i]))
        
        #print "len(stds) = ",len(stds)
        foundIt = False
        i = 0
        entireVar = np.std(data[:])
        print "entireVar for "+colHeaderStrs[col]+" is "+str(entireVar)#$$$$$$$$$$$$$$$$$$
        while foundIt==False:
            curVar = stds[i]
            #print "curVar = "+str(curVar)#$$$$$$$$$$$$$$$$$$$$$
            if curVar>(0.25*entireVar):
                foundIt = True
            # increment i to next value
            i += 1
            
        corrLength = i
        corrLengths.append(corrLength)
    
        if corrLength>maxCorrLength:
            maxCorrLength = corrLength
            maxCorrLengthCol = col
            maxEntireVar = entireVar
    print "CorrLengths "+repr(corrLengths)
    print "Maximum Correlation Length was found to be "+str(maxCorrLength)+" for the "+colHeaderStrs[maxCorrLengthCol]
    print "This gives Effective Length N = "+str(maxCorrLength)+"/"+str(len(data))+" = "+str(maxCorrLength/len(data))
    print "and parameter accuracy = "+str(maxEntireVar)+'/sqrt('+str(maxCorrLength/len(data))+") = "+str(maxEntireVar/math.sqrt(maxCorrLength/len(data)))
    
 