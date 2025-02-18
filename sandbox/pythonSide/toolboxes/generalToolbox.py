#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import math
import gc
import numpy as np
import os
import random as rand
import timeit
import pylab
import shutil
from math import pi

"""
This toolbox is a collection of the calculator type functions that were used in multiple 
places throughout the code to conduct various types of binary star system simulations.
"""     
     
def eccArgPeri2ToTcCalc(e, period, argPeri_deg, To, Tc=0):
    
    verbose = True
    backHalf = False
    
    if (verbose):
        print "\nInputs to eccArgPeri2ToTcCalc:"
        print "e = ",e
        print "argPeri_deg = ",argPeri_deg
        print "period = ",period
        print "Tc = ",Tc
        print "To = ",To
    
    if (e==0):
        if (Tc!=0):
            To = Tc
        else:
            Tc = To
    else:
        TA_s_deg = 90.0- argPeri_deg
        if backHalf:
            TA_s_deg = TA_s_deg+180.0
            
        if TA_s_deg<0.0:
            TA_s_deg =TA_s_deg+360.0
        TA_s_rad = np.radians(TA_s_deg)
        top = np.sqrt(1.0-e)*np.sin(TA_s_rad/2.0)   
        btm = np.sqrt(1.0+e)*np.cos(TA_s_rad/2.0) 
        ATAN_rad = math.atan2(top, btm)
        E_s_rad = ATAN_rad*2.0
        M_s_rad = E_s_rad-e*np.sin(E_s_rad)
        delta_t = (M_s_rad*period*365.242)/(2.0*pi)

        # in the case that the situation is flipped during the 'backHalf' option
        if backHalf:
            if (argPeri_deg>270):
                delta_t =  delta_t-period*365.242#
            elif (argPeri_deg==270):
                delta_t = 0.0
        # else apply standard To determination for Planet's orbit
        else:
            # check if argPeri is inside the back part of orbit where To>Tc
            if (argPeri_deg>90):
                delta_t =  delta_t-period*365.242#
            elif (argPeri_deg==90):
                delta_t = 0.0
        
        if Tc!=0:
            To = Tc-delta_t
        else:
            Tc = To+delta_t
            
        if (verbose):
            print "Outputs:"
            print 'TA_s_deg = ',TA_s_deg
            print "To = ",To
            print "Tc = ",Tc
            print "delta_t = ",delta_t
        
    return (To,Tc)
        
def bestOrbitFinderNEW(filename, printToScreen=True, saveToFile=True, returnAsList=False):
    """
    Just a simple function to find the parameters for the lowest chiSquared orbit in a 
    data file of the NEW format.
    """
  
    lowestChiSquared = 10000000
    inclinationBest = 0
    eBest = 0
    longANBest = 0
    periodBest = 0
    argPeriBest = 0
    aBest = 0
    TBest = 0
    TcBest = 0
    KBest = 0
    rvOffsetsBest=[]
    
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
            chiSquared = float(dataLineCols[8])
            
            if (chiSquared<lowestChiSquared)and(chiSquared>0.00001):
                lowestChiSquared = chiSquared
                incBest = float(dataLineCols[5])
                eBest = float(dataLineCols[1])
                longANBest = float(dataLineCols[0])
                periodBest = float(dataLineCols[4])
                argPeriBest = float(dataLineCols[6])
                aBest = float(dataLineCols[7])
                TBest = float(dataLineCols[2])
                TcBest = float(dataLineCols[3])
                if len(dataLineCols)>11:
                    KBest = float(dataLineCols[9])
                    rvOffsetsBest=[]
                    for dataset in range(0,int(len(dataLineCols) - 11)):
                        rvOffsetsBest.append(float(dataLineCols[10+dataset]))
            
    # print the values for the best orbit
    line= '\nBest orbit found:'
    line=line+ "\nLongAN = "+str(longANBest)
    line=line+ "\ne = "+str(eBest)
    line=line+ "\nTo = "+str(TBest)
    line=line+ "\nTc = "+str(TcBest)
    line=line+ "\nperiod = "+str(periodBest)
    line=line+ "\ninclination = "+str(incBest)
    line=line+ "\nargPeri = "+str(argPeriBest)
    line=line+ "\na_total = "+str(aBest)
    line=line+ "\nK = "+str(KBest)
    for dataset in range(0,len(rvOffsetsBest)):
        line = line+ "\nRV offset "+str(dataset+1)+" ="+str(rvOffsetsBest[dataset])
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
        if len(rvOffsetsBest)>0:
            list = [longANBest, eBest, TBest, TcBest, periodBest, incBest, argPeriBest, aBest, KBest, rvOffsetsBest]
        else:
            list = [longANBest, eBest, TBest, TcBest, periodBest, incBest, argPeriBest, aBest, KBest]
        return list 
    
def burnInCalc3(chiSquareds, medianALLchains,jumpy=True):
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
    verbose = False
    
    if jumpy:
        if len(chiSquareds)>=100000:
            jump = 100
        else:
            jump = 10
    else:
        jump = 1
    
    if verbose:
        print 'jump size = '+str(jump)
        print 'medianALLchains = '+str(medianALLchains)
        
    for i in range(0,int(len(chiSquareds)/jump)):
        # check if current chiSquared is less than the median of all chains yet
        chiSquareCur = chiSquareds[i*jump]
        if chiSquareCur<medianALLchains:
            if verbose:
                print 'found '+str(chiSquareCur)+' < '+str(medianALLchains)+' at step '+str(i*jump)+', so going back to previous jump step and finding exact location'
            if jump>1:
                # less than median, so do all in last jump to find precise location
                for j in range((i-1)*jump, i*jump):
                    chiSquareCur2 = chiSquareds[j]
                    if chiSquareCur2<medianALLchains:
                        burnInLength = j+1
                        if verbose:
                            print 'found '+str(chiSquareCur2)+' < '+str(medianALLchains)+' at step '+str(j)+', so returning burnInLength as '+str(burnInLength)
                        break
            else:
                burnInLength = i+1
            break
    if burnInLength == len(chiSquareds):
        print "PROBLEM: Param had a burn in length equal to param length, ie. the chain never burned in"
        
    return burnInLength

def CorrLengthCalc(paramIN):
    """
    This version uses np.std
    
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

def CorrLengthCalc2(paramIN):
    """
    This version uses np.var
    
    This function will calculate the correlation length and return its value.
    
    @param paramIN:    = parameter array after burn in data stripped
    @type paramIN:     = array (list) of doubles
    """
    verbose = False
    try:
        varALL = np.var(paramIN)
    except:
        useless=0
    halfVarALL = varALL/2.0
    CorrLength = 0
    
    if len(paramIN)>10e6:
        jump=1000
    elif len(paramIN)>1e5:
        jump = 100
    else:
        jump = 10
    
    if paramIN[0]==paramIN[-1]:
        if verbose:
            print 'First and last parameters were the same, so returning a correlation length of 0.'
        CorrLength=0
    else:
        for i in range(1,int(len(paramIN)/jump)):
            #check std at each jump to see if over halfstd yet
            try:
                varCur = np.var(paramIN[0:i*jump])
            except:
                useless=1
            if varCur>halfVarALL:
                # over halfStd, so do all in last jump to find precise location
                for j in range((i-1)*jump, i*jump):
                    try:
                        varCur2 = np.var(paramIN[0:j])
                    except:
                        useless=2
                    if varCur2>halfVarALL:
                        CorrLength = j+1
                        break
                break
        if CorrLength == len(paramIN):
            print "PROBLEM: Param had a correlation length equal to param length, ie. the chain never burned in"
    if verbose:
        print 'Correlation length found to be = ',CorrLength
    
    return CorrLength

def burnInCalcMultiFile(dataFilenames,simAnneal=True):
    """
    """
    verbose=False
    
    # push dataFilenames in to a list if not one already
    if type(dataFilenames)!=list:
        dataFilenames = [dataFilenames]
    
    chiSquaredsALL = np.array([])
    startMCMCsample=0
    
    if simAnneal:
        for filename in dataFilenames:
            if os.path.exists(filename):
                chiSquaredsChain = dataReaderNEW(filename,7)
                if startMCMCsample==0:
                    startMCMCsample = int(0.75*len(chiSquaredsChain))
                    if verbose:
                        print 'startMCMCsample found to be = '+str(startMCMCsample)
                chiSquaredsChain = chiSquaredsChain[startMCMCsample:]
                chiSquaredsALL=np.concatenate((chiSquaredsALL,chiSquaredsChain),axis=0)
    else:
        ALLfilename = os.path.join(os.path.dirname(dataFilenames[0]),'outputData-ALL.dat')
        chiSquaredsALL = dataReaderNEW(ALLfilename,7)
        
    # calculate median of 'all' array
    if type(chiSquaredsALL)!=np.ndarray:
        chiSquaredsALL = np.array(chiSquaredsALL)
        print 'chiSquaredsALL.shape = '+repr(chiSquaredsALL.shape)
    medainALL = np.median(chiSquaredsALL,axis=0)         
    
    if verbose:
        print "medainALL = "+str(medainALL)
    
    for filename in dataFilenames:
        if os.path.exists(filename):
            #find chain number and update logFilename with path and number
            s = filename
            chainNumStr = s[s.find('chain_')+6:s.find('chain_')+6+1]
            datadir = os.path.dirname(filename)
            logFilename = os.path.join(datadir,'log-chain_'+chainNumStr+'.txt')
            log = open(logFilename,'a')
            log.write('\n'+75*'+'+'\n Inside burnInCalc \n'+75*'+'+'\n')
            if simAnneal:
                log.write('Calculating the burn-in for the last 25% of the Simulated Annealing run:\n')
            else:
                log.write("Calculating the burn-in for MCMC run:\n")
            
            chiSquaredsChain = dataReaderNEW(filename,7)
            if simAnneal:
                chiSquaredsChain = chiSquaredsChain[startMCMCsample:-1]
            #medianChain = np.median(chiSquaredsChain)
            burnInLength = burnInCalc3(chiSquaredsChain, medainALL,jumpy=False)
            
            s = 'median value for all chains = '+str(medainALL)
            s = s+"\nTotal number of points in the chain = "+str(len(chiSquaredsChain))+"\n\n"
            s = s+"Burn-in length = "+str(burnInLength)+"\n\n"
            log.write(s+"\n\n")
            if verbose:
                print 'For chain # '+chainNumStr+s
            
    log.close()

def gelmanRubinStage2(dataFilenames):
    """
    This will finalize the calculation of the Gelman-Rubin statistic
    and write the results to a file in the same folder as the input files.  
    #Stage 1 of this is done in C++ and this one will just grab
    the results for each individual chain and combine to calculate the final inter-chain
    value.
    """
    verbose = True
    # push dataFilenames in to a list if not one already
    if type(dataFilenames)!=list:
        dataFilenames = [dataFilenames]
        
    # create arrays for the data from all the chains to go into
    # Lcs are the same for all chains, so just make a 1D array for it
    Lcs = []
    LcsLoaded = False
    means_multiDim_ALLChains = []
    vars_multiDim_ALLChains = []
    numChains = len(dataFilenames)
    # create empty strings for the datadir and headers for later use in writing output file
    datadir=''
    headers = ''
    
    #run through all files to load up their data into the ALLchain arrays
    for filename in dataFilenames:           
        if os.path.exists(filename):
            #find chain number and update logFilename with path and number
            s = filename
            chainNumStr = s[s.find('chain_')+6:s.find('chain_')+6+1]
            datadir = os.path.dirname(filename)
            logFilename = os.path.join(datadir,'log-chain_'+chainNumStr+'.txt')
            log = open(logFilename,'a')
            log.write('\n'+75*'+'+'\n Inside MCMCeffectivePointsCalc \n'+75*'+'+'\n')
            # figure out GR filename and open it
            GRfilename = os.path.join(datadir,'gelmanRubin-chain_'+chainNumStr+'.txt')
            if os.path.exists(GRfilename):
                if verbose:
                    print 'Getting data from file to calc 2nd stage of GR: '+GRfilename
                GRfile = open(GRfilename,'r')
                GRlines = GRfile.readlines()
                if verbose:
                    print " this file had "+str(len(GRlines))+' lines'
                # create arrays for the values in this chain's GR file
                means_multiDim = []
                vars_multiDim = []
                headers = GRlines[0]
                headersSplit = GRlines[0].split('  ')
                numParams = len(headersSplit)-1
                # load up multiDim arrays with the right number of subarrays
                for paramNum in range(0,numParams):
                    vars_multiDim.append([])
                    means_multiDim.append([])
                # run through all data lines in the GR file and load up the arrays & subarrays
                for line in GRlines[1:]:
                    lineSplit = line.split('  ')
                    if LcsLoaded==False:
                        Lcs.append(float(lineSplit[0]))
                    for paramNum in range(0,numParams):
                        vars_multiDim[paramNum].append(float(lineSplit[paramNum+1].split(',')[0]))
                        means_multiDim[paramNum].append(float(lineSplit[paramNum+1].split(',')[1]))
            else:
                print 'gelmanRubinStage2: Filename does not exist!!!: '+GRfilename
            #load values into the array for all the chain outputs
            LcsLoaded=True
            means_multiDim_ALLChains.append(means_multiDim)
            vars_multiDim_ALLChains.append(vars_multiDim)
    ## DONE loading all Lc,var and mean values for each param and each chain
    
    ## open file to write resulting Gelman-Rubin values into
    GRoutputFilename = os.path.join(datadir,'GRvalues.txt')
    GRoutputFile = open(GRoutputFilename,'w')
    GRoutputFile.write(headers+'\n')
    ToutputFilename = os.path.join(datadir,'Tvalues.txt')
    ToutputFile = open(ToutputFilename,'w')
    ToutputFile.write(headers+'\n')
    ## calculate the Gelman-Rubin stat for each itteration.
    numItterations = len(Lcs)
    print "numItterations = "+str(numItterations) #$$$$$$$$$$$$$$$$
    print "numChains = "+str(numChains) #$$$$$$$$$$$$$$$$
    GRs = []
    if numChains>1: 
        unBiasConversion = float(numChains)/float(numChains-1.0)
    else:
        unBiasConversion=1.0
    for itteration in range(0,numItterations):
        lineStr = ''
        lineStr2 = ''
        Lc = Lcs[itteration]
        lineStr=lineStr+str(Lc)
        lineStr2=lineStr2+str(Lc)
        for paramNum in range(0,numParams):
            means = []
            vars = []
            for chainNum in range(0,numChains):
                means.append(means_multiDim_ALLChains[chainNum][paramNum][itteration])
                vars.append(vars_multiDim_ALLChains[chainNum][paramNum][itteration])
            ## calculate values needed for R then calculate R
            B = Lc*unBiasConversion*np.var(means,axis=0)
            W = np.mean(vars,axis=0)
            weightedVar = ((Lc-1.0)/Lc)*W+(1.0/Lc)*B
            R = np.sqrt(weightedVar/W)
            lineStr=lineStr+'  '+str(R)
            # calculate 'T' from pg 26 of Ford2006
            # it is an "estimate of the effective number of independent draws"
            # and therefore good to compare to the correlation length
            T = Lc*numChains*np.min([weightedVar/B,1.0])
            lineStr2=lineStr2+'  '+str(T)
      
        # line loaded up with all the R values, so write it to output file
        GRoutputFile.write(lineStr+'\n')
        ToutputFile.write(lineStr2+'\n')
        
    print 'Output file with all Gelman-Rubin values for each parameter at each itteration written to: \n'+GRoutputFilename
    GRoutputFile.close()
    print 'Output file with all T values for each parameter at each itteration written to: \n'+ToutputFilename
    ToutputFile.close()
            
def MCMCeffectivePointsCalc(dataFilenames,simAnneal=False):
    """
    This function will calculate the correlation length of each parameter in the input file(s).
    If simAnneal=True, it will assume that only the last 25% of the chains are suitible for 
    calculating the correlation lengths.
    """
    verbose=False
    
    # push dataFilenames in to a list if not one already
    if type(dataFilenames)!=list:
        dataFilenames = [dataFilenames]
    
    # find start of MCMC part of simAnneal run if requested
    startMCMCsample=0
    if simAnneal:
        if os.path.exists(dataFilenames[0]):
            chiSquaredsChain = dataReaderNEW(dataFilenames[0],7)
            startMCMCsample = int(0.75*len(chiSquaredsChain))
                
    for filename in dataFilenames:           
        if os.path.exists(filename):
            #find chain number and update logFilename with path and number
            s = filename
            chainNumStr = s[s.find('chain_')+6:s.find('chain_')+6+1]
            datadir = os.path.dirname(filename)
            logFilename = os.path.join(datadir,'log-chain_'+chainNumStr+'.txt')
            log = open(logFilename,'a')
            log.write('\n'+75*'+'+'\n Inside MCMCeffectivePointsCalc \n'+75*'+'+'\n')
            
            ## find conff levels of data that is always outputed from sims
            s= '\nlongANs have:'
            data = dataReaderNEW(filename,0)
            if simAnneal:
                data = data[startMCMCsample:]
            s=s+ '\n[Min,Max] = '+repr([data.min(),data.max()])+", and median = "+str(np.median(data))
            #N_eff = effectivePointsCalcFunc(data)
            #s=s+ '\nEffective number of points = '+repr(N_eff)
            CorrLength = CorrLengthCalc2(data)
            if CorrLength == data.size:
                s=s+"\nPROBLEM: Param had a correlation length equal to param length, ie. the chain never burned in"
            else:
                s=s+'\nCorrelation length found to be = '+str(CorrLength)
            if CorrLength>0:
                s = s+',  and  '+str(data.size)+"/"+str(CorrLength)+" = "+str(int(data.size/CorrLength))
            log.write(s+'\n')
            if verbose:
                print s
            
            s= '\nes have'
            data = dataReaderNEW(filename,1)
            if simAnneal:
                data = data[startMCMCsample:]
            s=s+ '\n[Min,Max] = '+repr([data.min(),data.max()])+", and median = "+str(np.median(data))
            #N_eff = effectivePointsCalcFunc(data)
            #s=s+ '\nEffective number of points = '+repr(N_eff)
            CorrLength = CorrLengthCalc2(data)
            if CorrLength == data.size:
                s=s+"\nPROBLEM: Param had a correlation length equal to param length, ie. the chain never burned in"
            else:
                s=s+'\nCorrelation length found to be = '+str(CorrLength)
            if CorrLength>0:
                s = s+',  and  '+str(data.size)+"/"+str(CorrLength)+" = "+str(int(data.size/CorrLength))
            log.write(s+'\n')
            if verbose:
                print s
            
            s= '\nTs have:'
            data = dataReaderNEW(filename,2)
            if simAnneal:
                data = data[startMCMCsample:]
            s=s+ '\n[Min,Max] = '+repr([data.min(),data.max()])+", and median = "+str(np.median(data))
            #N_eff = effectivePointsCalcFunc(data)
            #s=s+ '\nEffective number of points = '+repr(N_eff)
            CorrLength = CorrLengthCalc2(data)
            if CorrLength == data.size:
                s=s+"\nPROBLEM: Param had a correlation length equal to param length, ie. the chain never burned in"
            else:
                s=s+'\nCorrelation length found to be = '+str(CorrLength)
            if CorrLength>0:
                s = s+',  and  '+str(data.size)+"/"+str(CorrLength)+" = "+str(int(data.size/CorrLength))
            log.write(s+'\n')
            if verbose:
                print s
            
            s= '\nperiods have:'
            data = dataReaderNEW(filename,3)
            if simAnneal:
                data = data[startMCMCsample:]
            s=s+ '\n[Min,Max] = '+repr([data.min(),data.max()])+", and median = "+str(np.median(data))
            #N_eff = effectivePointsCalcFunc(data)
            #s=s+ '\nEffective number of points = '+repr(N_eff)
            CorrLength = CorrLengthCalc2(data)
            if CorrLength == data.size:
                s=s+"\nPROBLEM: Param had a correlation length equal to param length, ie. the chain never burned in"
            else:
                s=s+'\nCorrelation length found to be = '+str(CorrLength)
            if CorrLength>0:
                s = s+',  and  '+str(data.size)+"/"+str(CorrLength)+" = "+str(int(data.size/CorrLength))
            log.write(s+'\n')
            if verbose:
                print s
            
            s= '\ninclinations have:'
            data= dataReaderNEW(filename,4)
            if simAnneal:
                data = data[startMCMCsample:]
            s=s+ '\n[Min,Max] = '+repr([data.min(),data.max()])+", and median = "+str(np.median(data))
            #N_eff = effectivePointsCalcFunc(data)
            #s=s+ '\nEffective number of points = '+repr(N_eff)
            CorrLength = CorrLengthCalc2(data)
            if CorrLength == data.size:
                s=s+"\nPROBLEM: Param had a correlation length equal to param length, ie. the chain never burned in"
            else:
                s=s+'\nCorrelation length found to be = '+str(CorrLength)
            if CorrLength>0:
                s = s+',  and  '+str(data.size)+"/"+str(CorrLength)+" = "+str(int(data.size/CorrLength))
            log.write(s+'\n')
            if verbose:
                print s
            
            s= '\nargPeris have:'
            data = dataReaderNEW(filename,5)
            if simAnneal:
                data = data[startMCMCsample:]
            s=s+ '\n[Min,Max] = '+repr([data.min(),data.max()])+", and median = "+str(np.median(data))
            #N_eff = effectivePointsCalcFunc(data)
            #s=s+ '\nEffective number of points = '+repr(N_eff)
            CorrLength = CorrLengthCalc2(data)
            if CorrLength == data.size:
                s=s+"\nPROBLEM: Param had a correlation length equal to param length, ie. the chain never burned in"
            else:
                s=s+'\nCorrelation length found to be = '+str(CorrLength)
            if CorrLength>0:
                s = s+',  and  '+str(data.size)+"/"+str(CorrLength)+" = "+str(int(data.size/CorrLength))
            log.write(s+'\n')
            if verbose:
                print s
            
            s= "\na_totals have:"
            data = dataReaderNEW(filename,6)
            if simAnneal:
                data = data[startMCMCsample:]
            s=s+ '\n[Min,Max] = '+repr([data.min(),data.max()])+", and median = "+str(np.median(data))
            #N_eff = effectivePointsCalcFunc(data)
            #s=s+ '\nEffective number of points = '+repr(N_eff)
            CorrLength = CorrLengthCalc2(data)
            if CorrLength == data.size:
                s=s+"\nPROBLEM: Param had a correlation length equal to param length, ie. the chain never burned in"
            else:
                s=s+'\nCorrelation length found to be = '+str(CorrLength)
            if CorrLength>0:
                s = s+',  and  '+str(data.size)+"/"+str(CorrLength)+" = "+str(int(data.size/CorrLength))
            log.write(s+'\n')
            if verbose:
                print s
            
            s= "\nKs have:"
            data = dataReaderNEW(filename,8)
            if simAnneal:
                data = data[startMCMCsample:]
            s=s+ '\n[Min,Max] = '+repr([data.min(),data.max()])+", and median = "+str(np.median(data))
            #N_eff = effectivePointsCalcFunc(data)
            #s=s+ '\nEffective number of points = '+repr(N_eff)
            CorrLength = CorrLengthCalc2(data)
            if CorrLength == data.size:
                s=s+"\nPROBLEM: Param had a correlation length equal to param length, ie. the chain never burned in"
            else:
                s=s+'\nCorrelation length found to be = '+str(CorrLength)
            if CorrLength>0:
                s = s+',  and  '+str(data.size)+"/"+str(CorrLength)+" = "+str(int(data.size/CorrLength))
            log.write(s+'\n')
            if verbose:
                print s
        
            ## figure out if there is any RV offsets in output file and find their confLevels 
            f = open(filename, 'r')
            plotFileTitle = f.readline()
            headings = f.readline()
            dataline = f.readline()
            dataLineCols = dataline.split()
            numDataCols = len(dataLineCols)
            f.close()
            if numDataCols==10:
                s= 'There were 9 columns of data found in the datafile, thus no RVoffsets were recorded'
                log.write(s+'\n')
                if verbose:
                    print s
            elif numDataCols>10:
                s= 'There were '+str(numDataCols)+' columns of data, thus '+str(numDataCols - 10)+ ' columns must be RV offsets' 
                log.write(s+'\n')
                if verbose:
                    print s
                numRVdatasets = numDataCols - 10
                for dataset in range(0,numRVdatasets):
                    s= '\ndataset # '+str(dataset+1)+' RV offsets have:'
                    data = dataReaderNEW(filename,dataset+9)
                    if simAnneal:
                        data = data[startMCMCsample:]
                    s=s+ '\n[Min,Max] = '+repr([data.min(),data.max()])+", and median = "+str(np.median(data))
                    #N_eff = effectivePointsCalcFunc(data)
                    #s=s+ '\nEffective number of points = '+repr(N_eff)
                    CorrLength = CorrLengthCalc2(data)
                    if CorrLength == data.size:
                        s=s+"\nPROBLEM: Param had a correlation length equal to param length, ie. the chain never burned in"
                    else:
                        s=s+'\nCorrelation length found to be = '+str(CorrLength)
                    if CorrLength>0:
                        s = s+',  and  '+str(data.size)+"/"+str(CorrLength)+" = "+str(int(data.size/CorrLength))
                    log.write(s+'\n')
                    if verbose:
                        print s
        else:
            s= "confidenceLevelsFinderNEW: ERROR!!!! file doesn't exist"
            print s
            log.write(s+'\n')
        log.write('\n'+75*'+'+'\n Leaving MCMCeffectivePointsCalc \n'+75*'+'+'\n')
        log.close()
        
def effectivePointsCalcFunc(data):
    """
    for use by MCMCeffectivePointsCalc
    """
    # Important: if your chain is 1 million steps, then ecc should have 
    # 1 million elements.  Don't have a separate array indicating how 
    # many times you visited each point.  N_pts is the total number of 
    # steps in the following code.
    
    # This is the total variance in eccentricity in your chain.  It should
    # roughly be the (squared) 1 sigma error in eccentricity.
    
    if type(data)!=np.ndarray:
        data = np.array(data)
    var_obs = np.var(data)  
    N_pts = data.size
    #print 'Number of data points in total chain is = ',N_pts
    N_eff=0
    for di in range(1, N_pts):
        # Make sure the array dimensions are compatible for resizing.
        newdim = (N_pts // di) * di
    
        # This returns an array of variances: the variance of points 
        # 0, ..., 0 + di - 1, then di, ..., 2 * di - 1, then 2 * di, ... etc.
        var_di = np.var(np.reshape(data[:newdim], (N_pts // di, di)), axis=1)
        
        # If the typical sub-variance is at least half the total variance,
        # we have found the effective length of a sub-chain
    
        if np.mean(var_di) > 0.5 * var_obs:
            N_eff = N_pts // di
            break
    
    #print str(N_eff) + ' effective points in the MCMC chain.'
    
    return N_eff
    
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
        if (outFilename[-4:]!='.txt' and outFilename[-4:]!='.dat'):
            print 'Changing output filename from '+outFilename+' to '+outFilename+'.dat'
            outFilename = outFilename+'.dat'
    
    numFiles = len(filenames)
    if numFiles==1:
        print 'Only one data file, so just copying it to the a file with the output filename'
        shutil.copy(filenames[0],outFilename)
    elif numFiles>1:
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
        longAN [deg]      e [N/A]       To [julian date]  period [yrs]   inclination [deg]   argPeri [deg]   a_total [AU]    chiSquared  RVoffset0...  timesBeenHere
        
        file format must be:
        
        line1: title
        line2: data headings
        line3: space delimited data
        line4: space delimited data
        .
        .
        .
    """
    verbose = False
    ## instantiate data list
    data = []
    
    if verbose:
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
      
    if verbose:          
        print "number of values in data column = "+str(len(data))
    
    # convert to numpy array
    data = np.array(data)
    
    return data


def outputDatafileToDict(filename):
    """
    CAUTION, this was designed to be used with the outputs of 100ModDataset simulation outputs.
    NOT a good function to load data from a long simulation with lots of output sets.
    
    NOTE: This version is meant for the new headings from mcmcOrbSimUniform6
        longAN [deg]      e [N/A]       To [julian date]  period [yrs]   inclination [deg]   argPeri [deg]   a_total [AU]  chiSquared   K [m/s]  RVoffset0...  timesBeenHere
        
        file format must be:
        
        line1: title
        line2: data headings
        line3: space delimited data
        line4: space delimited data
        .
        .
        .
    """
    if os.path.exists(filename):
        outDict = {}
        print "\nWorking on file: "+os.path.basename(filename)
        f = open(filename, 'r')
        # strip off the .txt part to make the plot version of the filename
        plotFileTitle = f.readline()[:-5]
        headings = f.readline()
        dataline = f.readline()
        dataLineCols = dataline.split()
        numDataCols = len(dataLineCols)
        
        if numDataCols==9:
            print 'There were 9 columns of data found in the datafile, thus no RVoffsets were recorded'
        elif numDataCols>9:
            print 'There were '+str(numDataCols)+' columns of data, thus '+str(numDataCols - 9)+ ' columns must be RV offsets' 
            outDict["longANs"] = dataReaderNEW(filename, column=0)
            outDict["es"] = dataReaderNEW(filename, column=1)
            outDict["Ts"] = dataReaderNEW(filename, column=2)
            outDict["periods"] = dataReaderNEW(filename, column=3)
            outDict["inclinations"] = dataReaderNEW(filename, column=4)
            outDict["argPeris"] = dataReaderNEW(filename, column=5)
            outDict["a_totals"] = dataReaderNEW(filename, column=6)
            outDict["chiSquareds"] = dataReaderNEW(filename, column=7)
            outDict["Ks"]= dataReaderNEW(filename, column=8)
            # last column should be number of times been here
            outDict["timesBeenHere"] = dataReaderNEW(filename, column=(numDataCols-1))
            
            RVoffsets = []
            numRVdatasets = numDataCols - 9
            for dataset in range(0,numRVdatasets):
                colnum = int(8+dataset)
                RVoffsetsCurr = dataReaderNEW(filename, column=colnum)
                if numRVdatasets==1:
                    RVoffsets = RVoffsetsCurr
                else:
                    RVoffsets.append(RVoffsetsCurr)
            outDict["RVoffsets"] = RVoffsets
        
        print 'There were '+str(len(RVoffsets))+' rows of data loaded into the dictionary'
        return outDict
    else:
        print "outputDatafileToDict: ERROR!!!! file doesn't exist"
        
    
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
    if (filename[-4:]!='.txt' and filename[-4:]!='.dat'):
        INSfilename = filename+'.dat'
    else:
        INSfilename = filename

    titleMod = '-chiSquare-cut-off-'+str(chiSquareCutOff)
    # create output version of INS file and write new 
    # strip off '.txt', add titleMod, then put it back on
    newFilename = INSfilename[:-4]+titleMod+'.dat'
    
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
            chiSquareStr = dataLineCols[7]
            chiSquared = float(chiSquareStr)
            if chiSquared<chiSquareCutOff:
                OUTfile.write(line)
                totalKept = totalKept+1
    if verbose:
        print str(totalKept)+' orbits were found to have a chiSquared < '+str(chiSquareCutOff)
        print 'Done loading, trimming and writing. Output file = '+newFilename
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
    if '.dat' not in filename:
        INSfilename = filename+'.dat'
    else:
        INSfilename = filename

    titleMod = '-chiSquare-cut-off-'+str(chiSquareCutOff)
    # create output version of INS file and write new 
    # strip off '.txt', add titleMod, then put it back on
    newFilename = INSfilename[:-4]+titleMod+'.dat'
    
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
            #if len(dataLineCols)==10:
            chiSquareStr = dataLineCols[8]
#            elif len(dataLineCols)==9:
#                chiSquareStr = dataLineCols[7]
            chiSquared = float(chiSquareStr)
#            inclination = float(dataLineCols[4])
#            if chiSquared<chiSquareCutOff:
#                if inclination<90.0:
            if chiSquared<chiSquareCutOff:
                OUTfile.write(line)
                totalKept = totalKept+1
    if verbose:
        print str(totalKept)+' orbits were found to have a chiSquared < '+str(chiSquareCutOff)
        print 'Done loading, trimming and writing. Output file = '+os.path.basename(INSfilename)
    INfile.close()
    OUTfile.close()    
    
    print 'Final trimmed data written to: '+newFilename
    
    if verbose:    
        print '** All data from file trimmed and written to output file **'
        
    if True:
        return newFilename
        
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

def confLevelFinderNEWchiSquaredVersion(filename, columNum=False, returnData=False):
    """
    A function to find the 68.3 and 95.4% confidence levels in a given output data file's column.
    
    return [[68.3% minimum, 68.3% maximum],[95.5% minimum, 95.5% maximum]]
    
    columnNum must be an int.
    
    file format must be:
        
        line1: title
        line2: data headings
        line3: space delimited data
        line4: space delimited data
    
    """
    verboseInternal = False
    
    ## First get ranges of param and ChiSquared values
    if os.path.exists(filename):
        if verboseInternal:
            print '\nOpening and finding ranges for data in column # '+str(columNum)
        
        totalAccepted = 0
        f = open(filename,'r')
        plotFileTitle = f.readline()[:-5]
        headings = f.readline()
        dataAry = []
        chiSquareds = []
        line = 'asdf'
        chiSquaredMin=1e6
        chiSquaredMax = 0.0
        while line!='':
            line = f.readline()
            if line!='':
                totalAccepted+=1
                dataLineCols = line.split()
                dataValue = float(dataLineCols[columNum])
                dataAry.append(dataValue)
                chiSquared = float(dataLineCols[7])
                chiSquareds.append(chiSquared)
                #if dataValue>dataMax:
                #    dataMax = dataValue
                #if dataValue<dataMin:
                #    dataMin = dataValue
                if chiSquared<chiSquaredMin:
                    chiSquaredMin = chiSquared
                    bestDataVal = dataValue
                if chiSquared>chiSquaredMax:
                    chiSquaredMax = chiSquared
        f.close()
        
        #Convert data array to a numpy array
        dataAry = np.array(dataAry)
        chiSquareds = np.array(chiSquareds)
        # find range of data's values
        chiSquaredRange = chiSquaredMax-chiSquaredMin
        if verboseInternal:
            print '\nTotal number of orbits = '+str(totalAccepted)
            print 'Best value found was '+str(bestDataVal)+", which had a reduced chiSquared = "+str(chiSquaredMin)
            print '[Min,Max] values found for chiSquareds were '+repr([chiSquaredMin,chiSquaredMax])
            
        # make a step size = 1000th the range to ensure accuracy of 0.1%
        numDataSteps = 2000.0
        deltaBestStep = chiSquaredRange/numDataSteps
        
        deltaBestCurr = deltaBestStep
        conf68ValsRough = []
        conf95ValsRough = []
        conf68Vals = []
        conf95Vals = []
        chiSquared68rough=0
        chiSquared95rough=0
        chiSquared68=0
        chiSquared95=0
        # loop through till both 68 and 95 levels are found, or until last step possible reached
        while (((len(conf68Vals)==0) or (len(conf95Vals)==0))and(deltaBestCurr<deltaBestStep*(numDataSteps/2.0))):
            #update max and min for range
            curChiSquaredMax = chiSquaredMin+deltaBestCurr
            #calc number of points within this range
            curSize = np.where(chiSquareds<curChiSquaredMax)[0].size
            # calc current percentage included
            curPercent = float(curSize)/float(totalAccepted)
            if (verboseInternal and False):
                print 'Current max chiSquared = '+repr(curChiSquaredMax)+", contains "+str(curPercent)+'% of the data'
            
            # if inside percent ranges for the 68.3 or 95.4% load up confVals list
            if (((curPercent<0.690)and(curPercent>0.680))and(len(conf68Vals)==0)):
                dataSubAry = dataAry[np.where(chiSquareds<curChiSquaredMax)]
                curMin = dataSubAry.min()
                curMax = dataSubAry.max()
                conf68ValsRough = [curMin,curMax]
                chiSquared68rough = curChiSquaredMax
                if (((curPercent<0.684)and(curPercent>0.682))and(len(conf68Vals)==0)):
                    conf68Vals = [curMin,curMax]
                    chiSquared68 = curChiSquaredMax
                    if verboseInternal:
                        print "\n68.3% range found to be "+repr(conf68Vals)
                        print "for all data inside chiSquaredMax = "+str(curChiSquaredMax)
            elif (((curPercent<0.960)and(curPercent>0.950))and(len(conf95Vals)==0)):
                dataSubAry = dataAry[np.where(chiSquareds<curChiSquaredMax)]
                curMin = dataSubAry.min()
                curMax = dataSubAry.max()
                conf95ValsRough = [curMin,curMax]
                chiSquared95rough = curChiSquaredMax
                if (((curPercent<0.955)and(curPercent>0.953))and(len(conf95Vals)==0)):
                    conf95Vals = [curMin,curMax]
                    chiSquared95 = curChiSquaredMax
                    if verboseInternal:
                        print "\n95.4% range found to be "+repr(conf95Vals)
                        print "for all data inside chiSquaredMax = "+str(curChiSquaredMax)+'\n'
            
            # increase delta size
            deltaBestCurr = deltaBestCurr+deltaBestStep
        
        if ((len(conf68Vals)==0) or (len(conf95Vals)==0)):
            if (len(conf68Vals)==0):
                print 'confLevelFinderNEWchiSquaredVersion: ERROR!!! No FINE 68.3% confidence levels were found'
                if (len(conf68ValsRough)==0):
                    print 'confLevelFinderNEWchiSquaredVersion: ERROR!!! No ROUGH 68% confidence levels were found, so returning [0,0]'
                    conf68Vals = [0,0]
                else:
                    chiSquared68 = chiSquared68rough
                    conf68Vals = conf68ValsRough
                    print "confLevelFinderNEWchiSquaredVersion: Had to use ROUGH 68% [68,69] as no FINE 68.3% was found. So, using range "+repr(conf68Vals)+", for all data inside chiSquaredMax = "+str(chiSquared68)+'\n'
            
            if (len(conf95Vals)==0):
                print 'confLevelFinderNEWchiSquaredVersion: ERROR!!! No FINE 95.4% confidence levels were found'
                if (len(conf95ValsRough)==0):
                    print 'confLevelFinderNEWchiSquaredVersion: ERROR!!! No ROUGH 95% confidence levels were found, so returning [0,0]'
                    conf95Vals = [0,0]
                else:
                    chiSquared95 = chiSquared95rough
                    conf95Vals = conf95ValsRough
                    print "confLevelFinderNEWchiSquaredVersion: Had to use ROUGH 95% [95,96] as no FINE 95.4% was found. So, using range "+repr(conf95Vals)+", for all data inside chiSquaredMax = "+str(chiSquared95)+'\n'
                
        if returnData:
            return [conf68Vals,conf95Vals],dataAry
        else:
            return [conf68Vals,conf95Vals]
    else:
        print "confLevelFinderNEWchiSquaredVersion: ERROR!!!! file doesn't exist"

def confLevelFinderNEWdataVersion(filename, columNum=False, returnData=False, returnChiSquareds=False, returnBestDataVal=False,fast=True):
    """
    A function to find the 68.3 and 95.4% confidence levels in a given output data file's column.
    
    return [[68.3% minimum, 68.3% maximum],[95.5% minimum, 95.5% maximum]]
    
    columnNum must be an int.
    
    file format must be:
        
        line1: title
        line2: data headings
        line3: space delimited data
        line4: space delimited data
    
    """
    verboseInternal = False
    bestCentered = False

    gotLog=True
    ## First get ranges of param and ChiSquared values
    if os.path.exists(filename):
        ## figure out chain number and open log
        s = filename
        datadir = os.path.dirname(filename)
        chainNumStr = s[s.find('chain_')+6:s.find('chain_')+6+1]
        if chainNumStr.isdigit():  
            logFilename = os.path.join(datadir,'log-chain_'+chainNumStr+'.txt')
        else:
            logFilename = os.path.join(datadir,'processManagerLogFile.txt')
        if os.path.exists(logFilename):
            log = open(logFilename,'a')
            log.write('\n'+75*'-'+'\n Inside confLevelFinderNEWdataVersion \n'+75*'-'+'\n')
        else:
            gotLog=False
        
        s= '\nOpening and finding ranges for data in column # '+str(columNum)
        if gotLog:
            log.write(s+'\n')
        if verboseInternal:
            print s
        
        ## Check if file has useful data for that column#
        # first find out how many lines in total
        fp = open(filename,'r')
        
        TotalSamples=1
        for i,line in enumerate(fp):
            if i>2:
                TotalSamples+= float(line.split()[-1])
        fp.close()
        TotalSamples = int(TotalSamples)
        numDataLines =i-2
        # find values at start, mid and end of file
        fp = open(filename,'r')
        lastColLoc = 0
        dataValueStart = dataValueMid = dataValueEnd =0
        for i,line in enumerate(fp):
            if i==(0+2):
                splitAry = line.split()
                lastColLoc = len(splitAry)-1
                dataValueStart = float(splitAry[columNum])
            elif i==((numDataLines//2)+2):
                splitAry = line.split()
                dataValueMid = float(splitAry[columNum])
            elif i==numDataLines:
                splitAry = line.split()
                dataValueEnd = float(splitAry[columNum])
        fp.close()
        
        doesntVary = True
        if ((dataValueStart!=dataValueMid)and(dataValueStart!=dataValueEnd)):
            doesntVary = False
            
        if (doesntVary==False)or(fast==False):  
            s=''
            if True:
                #Old string parsing directly version UPDATED
                startTime2 = timeit.default_timer()
                totalAccepted = 0
                dataAry = [None]*TotalSamples
                chiSquareds = [None]*TotalSamples
                i = 0
                chiSquaredMin=1e6
                bestOrbit = 0
                lineNum=0
                for line in open(filename,'r'):
                    lineNum+=1
                    if lineNum>2:
                        try:
                            dataLineCols = line.split()
                            ## this should never happen, but it is a check for a double decimal value
                            decimalSplit = dataLineCols[columNum].split('.')
                            if len(decimalSplit)>2:
                                dataValue = float(decimalSplit[0]+'.'+decimalSplit[1])
                            else:
                                dataValue = float(dataLineCols[columNum])
                            chiSquared = float(dataLineCols[8])
                            timesBeenHere = float(dataLineCols[-1])
                            if (chiSquared==0)or(dataValue==0):
                                if verboseInternal:
                                    print line
                            for j in range(0,int(timesBeenHere)):
                                totalAccepted+=1
                                dataAry[i]=dataValue
                                if returnChiSquareds:
                                    chiSquareds[i]=chiSquared
                                i+=1
                            #if dataValue>dataMax:
                            #    dataMax = dataValue
                            #if dataValue<dataMin:
                            #    dataMin = dataValue
                            if chiSquared<chiSquaredMin:
                                chiSquaredMin = chiSquared
                                bestDataVal = dataValue
                                bestOrbit=lineNum      
                        except:
                            print 'Failed for line: '+line 
                endTime2 = timeit.default_timer()
                totalTime = (endTime2-startTime2) # in seconds
                totalTimeString = timeString(totalTime)
                s=s+ '\nUPDATED Direct data loading took '+totalTimeString+' to complete.\n' 
                s=s+'The resulting arrays had '+str(totalAccepted)+' elements, with best value '+str(bestDataVal)+', and minChiSquared '+str(chiSquaredMin)
#            if True:
#                #Old string parsing directly version UPDATED
#                startTime2 = timeit.default_timer()
#                totalAccepted = 0
#                dataAry = [None]*TotalSamples
#                chiSquareds = [None]*TotalSamples
#                i = 0
#                chiSquaredMin=1e6
#                lineNum=0
#                for line in open(filename,'r'):
#                    lineNum+=1
#                    if lineNum>2:
#                        dataLineCols = line.split()
#                        ## this should never happen, but it is a check for a double decimal value
#                        decimalSplit = dataLineCols[columNum].split('.')
#                        if len(decimalSplit)>2:
#                            dataValue = float(decimalSplit[0]+'.'+decimalSplit[1])
#                        else:
#                            dataValue = float(dataLineCols[columNum])
#                        chiSquared = float(dataLineCols[8])
#                        timesBeenHere = int(float(dataLineCols[-1]))
#                        if (chiSquared==0)or(dataValue==0):
#                            if verboseInternal:
#                                print line
#                        totalAccepted+=timesBeenHere
#                        dataAry[i:i+timesBeenHere]=dataValue
#                        if returnChiSquareds:
#                            chiSquareds[i:i+timesBeenHere]=chiSquared
#                        i+=timesBeenHere
#                        if chiSquared<chiSquaredMin:
#                            chiSquaredMin = chiSquared
#                            bestDataVal = dataValue
#                endTime2 = timeit.default_timer()
#                totalTime = (endTime2-startTime2) # in seconds
#                totalTimeString = timeString(totalTime)
#                s=s+ '\nUPDATED2 Direct data loading took '+totalTimeString+' to complete.\n' 
#                s=s+'The resulting arrays had '+str(totalAccepted)+' elements, with best value '+str(bestDataVal)+', and minChiSquared '+str(chiSquaredMin)
#            if False:
#                #Old string parsing directly version UPDATED WITH NUMPY suggestions from stackoverflow
#                startTime2 = timeit.default_timer()
#                totalAccepted = 0
#                dataAry = np.empty(TotalSamples)
#                chiSquareds = np.empty(TotalSamples)
#                i = 0
#                chiSquaredMin=1e6
#                lineNum=0
#                for line in open(filename,'r'):
#                        lineNum+=1
#                        if lineNum>2:
#                            dataLineCols = line.split()
#                            ## this should never happen, but it is a check for a double decimal value
#                            decimalSplit = dataLineCols[columNum].split('.')
#                            if len(decimalSplit)>2:
#                                dataValue = float(decimalSplit[0]+'.'+decimalSplit[1])
#                            else:
#                                dataValue = float(dataLineCols[columNum])
#                            chiSquared = float(dataLineCols[8])
#                            timesBeenHere = float(dataLineCols[-1])
#                            if (chiSquared==0)or(dataValue==0):
#                                if verboseInternal:
#                                    print line
#                            totalAccepted+=timesBeenHere
#                            dataAry[i:i+timesBeenHere]=dataValue
#                            if returnChiSquareds:
#                                chiSquareds[i:i+timesBeenHere]=chiSquared
#                            i+=timesBeenHere
#                            if chiSquared<chiSquaredMin:
#                                chiSquaredMin = chiSquared
#                                bestDataVal = dataValue
#                #f.close()
#                endTime2 = timeit.default_timer()
#                totalTime = (endTime2-startTime2) # in seconds
#                totalTimeString = timeString(totalTime)
#                s=s+ '\nUPDATED WITH NUMPY Direct data loading took '+totalTimeString+' to complete.\n' 
#                s=s+'The resulting arrays had '+str(totalAccepted)+' elements, with best value '+str(bestDataVal)+', and minChiSquared '+str(chiSquaredMin)
            print s
            #Convert data array to a numpy array
            dataAry = np.sort(dataAry)
            # find range of data's values
            dataMax = dataAry[0]
            dataMin = dataAry[-1]
            size = dataAry.size
            if (size%2)==0:
                dataMedian = (dataAry[size / 2 - 1] + dataAry[size / 2]) / 2
            else:
                dataMedian = dataAry[size / 2]
            s=  '\nTotal number of orbits = '+str(totalAccepted)
            s=s+'\nBest value found was '+str(bestDataVal)+", at line Number "+str(bestOrbit)+", and had a chiSquared = "+str(chiSquaredMin)
            s=s+'\nMedian value = '+str(dataMedian)
            s=s+'\n[Min,Max] values found for data were '+repr([dataMin,dataMax])
            if gotLog:
                log.write(s+'\n')
            if verboseInternal:
                print s
        
            if bestCentered:
                mid = np.where(dataAry==bestDataVal)[0][0]
            else:
                mid=size//2
                
            minLoc68=mid-int(float(size)*0.683)//2
            if minLoc68<0:
                minLoc68 = 0
            maxLoc68 = mid+int(float(size)*0.683)//2
            if maxLoc68>(size-1):
                maxLoc68 = size
            minLoc95=mid-int(float(size)*0.958)//2
            if minLoc95<0:
                minLoc95 = 0
            maxLoc95= mid+int(float(size)*0.958)//2
            if maxLoc95>(size-1):
                maxLoc95 = size
            
            conf68Vals = [dataAry[minLoc68],dataAry[maxLoc68]]
            conf95Vals = [dataAry[minLoc95],dataAry[maxLoc95]]
            conf68ValsRough=[]
            conf95ValsRough=[]
            
            if ((len(conf68Vals)==0) or (len(conf95Vals)==0)):
                if (len(conf68Vals)==0):
                    s= 'confLevelFinderNEWdataVersion: ERROR!!! No FINE 68.3% confidence levels were found'
                    if gotLog:
                        log.write(s+'\n')
                    if verboseInternal:
                        print s
                    if (len(conf68ValsRough)==0):
                        s= 'confLevelFinderNEWdataVersion: ERROR!!! No ROUGH 68% confidence levels were found, so returning [0,0]'
                        if gotLog:
                            log.write(s+'\n')
                        if verboseInternal:
                            print s
                        conf68Vals = [0,0]
                    else:
                        conf68Vals = conf68ValsRough
                        s= "confLevelFinderNEWdataVersion: Had to use ROUGH 68% [68,69] as no FINE 68.3% was found. So, using range "+repr(conf68Vals)+'\n'
                        if gotLog:
                            log.write(s+'\n')
                        if verboseInternal:
                            print s
                
                if (len(conf95Vals)==0):
                    s= 'confLevelFinderNEWdataVersion: ERROR!!! No FINE 95.4% confidence levels were found'
                    if gotLog:
                        log.write(s+'\n')
                    if verboseInternal:
                        print s
                    if (len(conf95ValsRough)==0):
                        s= 'confLevelFinderNEWdataVersion: ERROR!!! No ROUGH 95% confidence levels were found, so returning [0,0]'
                        if gotLog:
                            log.write(s+'\n')
                        if verboseInternal:
                            print s
                        conf95Vals = [0,0]
                    else:
                        conf95Vals = conf95ValsRough
                        s= "confLevelFinderNEWdataVersion: Had to use ROUGH 95% [95,96] as no FINE 95.4% was found. So, using range "+repr(conf95Vals)+'\n'
                        if gotLog:
                            log.write(s+'\n')
                        if verboseInternal:
                            print s
            
        else:
            ## There was no useful data, so return values indicating that
            
            dataAry=bestDataVal=dataMedian=dataValueStart
            conf68Vals = [dataValueStart,dataValueStart]
            conf95Vals = [dataValueStart,dataValueStart]
            chiSquareds = 0
            
            
        if gotLog:
            s= "Final 68% range values are: "+repr(conf68Vals)+'\n'
            s=s+"Final 95% range values are: "+repr(conf95Vals)+'\n'
            if bestCentered:
                s=s+ "\nerror is centered on best \n"
                s=s+"68.3% error level = "+str(bestDataVal-conf68Vals[0])+'\n'
                s=s+" =>   "+str(dataMedian)+'  +/-  '+str(bestDataVal-conf68Vals[0])+'\n'
            else:
                s=s+ "\nerror is centered on Median \n"
                s=s+"68.3% error level = "+str(dataMedian-conf68Vals[0])
                s=s+" =>   "+str(dataMedian)+'  +/-  '+str(dataMedian-conf68Vals[0])+'\n'
            
            s=s+'\n'+75*'-'+'\n Leaving confLevelFinderNEWdataVersion \n'+75*'-'+'\n'
            log.write(s)
            log.close()
        if verboseInternal:
            print 'returnData = '+repr(returnData)+', returnChiSquareds = '+repr(returnChiSquareds)+', returnBestDataVal = '+repr(returnBestDataVal)
        if (returnData and returnChiSquareds and (returnBestDataVal==False)):
            if verboseInternal:
                print 'returning first 3'
            returnList =  ([conf68Vals,conf95Vals],dataAry, chiSquareds)
        if (returnData and returnChiSquareds and returnBestDataVal):
            if verboseInternal:
                print 'returning all 4'
            returnList =   ([conf68Vals,conf95Vals],dataAry, chiSquareds, bestDataVal)
        if (returnData and (returnChiSquareds==False)and (returnBestDataVal==False)):
            if verboseInternal:
                print 'returning data only'
            returnList =   ([conf68Vals,conf95Vals],dataAry)
        if (returnData and (returnChiSquareds==False) and returnBestDataVal):
            if verboseInternal:
                print 'returning data and bestval'
            returnList =   ([conf68Vals,conf95Vals],dataAry, bestDataVal)
        if ((returnData==False) and returnChiSquareds):
            if verboseInternal:
                print 'returning just chiSquareds'
            returnList =   ([conf68Vals,conf95Vals], chiSquareds)
        if ((returnData==False) and returnChiSquareds and returnBestDataVal):
            if verboseInternal:
                print 'returning chiSquareds and bestval'
            returnList =   ([conf68Vals,conf95Vals], chiSquareds, bestDataVal)
        if ((returnData==False)and(returnChiSquareds==False) and (returnBestDataVal==False)):
            if verboseInternal:
                print 'returning only CLevels'
            returnList =   [conf68Vals,conf95Vals]
        
        return returnList 
    else:
        s= "confLevelFinderNEWdataVersion: ERROR!!!! file doesn't exist"
        print s
           
        
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
    if os.path.exists(filename):
        ## find conff levels of data that is always outputed from sims
        longAN_degsCLevels = confLevelFinderNEWdataVersion(filename,0)
        print '\nlongANs have conf levels: \n68.3% = '+repr(longAN_degsCLevels[0])+'\n95.4% = '+repr(longAN_degsCLevels[1])+'\n'
        esCLevels = confLevelFinderNEWdataVersion(filename,1)
        print 'es have conf levels: \n68.3% = '+repr(esCLevels[0])+' \n95.4% = '+repr(esCLevels[1])+'\n'
        TsCLevels = confLevelFinderNEWdataVersion(filename,2)
        print 'Ts have conf levels: \n68.3% = '+repr(TsCLevels[0])+' \n95.4% = '+repr(TsCLevels[1])+'\n'
        periodsCLevels = confLevelFinderNEWdataVersion(filename,3)
        print 'periods have conf levels: \n68.3% = '+repr(periodsCLevels[0])+' \n95.4% = '+repr(periodsCLevels[1])+'\n'
        inclination_degsCLevels= confLevelFinderNEWdataVersion(filename,4)
        print 'inclinations have conf levels: \n68.3% = '+repr(inclination_degsCLevels[0])+' \n95.4% = '+repr(inclination_degsCLevels[1])+'\n'
        argPeri_degsCLevels = confLevelFinderNEWdataVersion(filename,5)
        print 'argPeris have conf levels: \n68.3% = '+repr(argPeri_degsCLevels[0])+' \n95.4% = '+repr(argPeri_degsCLevels[1])+'\n'
        asCLevels = confLevelFinderNEWdataVersion(filename,6)
        print 'a_totals have conf levels: \n68.3% = '+repr(asCLevels[0])+' \n95.4% = '+repr(asCLevels[1])+'\n'
        KsCLevels = confLevelFinderNEWdataVersion(filename,8)
        print 'Ks have conf levels: \n68.3% = '+repr(asCLevels[0])+' \n95.4% = '+repr(asCLevels[1])+'\n'
    
        ## figure out if there is any RV offsets in output file and find their confLevels 
        f = open(filename, 'r')
        plotFileTitle = f.readline()
        headings = f.readline()
        dataline = f.readline()
        dataLineCols = dataline.split()
        numDataCols = len(dataLineCols)
        f.close()
        if numDataCols==9:
            print 'There were 9 columns of data found in the datafile, thus no RVoffsets were recorded'
        elif numDataCols>9:
            print 'There were '+str(numDataCols)+' columns of data, thus '+str(numDataCols - 9)+ ' columns must be RV offsets' 
            numRVdatasets = numDataCols - 9
            for dataset in range(0,numRVdatasets):
                offsetsCurrCLevels = confLevelFinderNEWdataVersion(filename,dataset+8)
                print 'dataset # '+str(dataset+1)+' RV offsets have conf levels: \n68.3% = '+repr(offsetsCurrCLevels[0])+' \n95.4% = '+repr(offsetsCurrCLevels[1])+'\n'
    else:
        print "confidenceLevelsFinderNEW: ERROR!!!! file doesn't exist"    
    
def cFileToSimSettingsDict(inputSettingsFile, outputSettingsFile=""):
    """
    This will load in a settings file written in the format for C++ originally 
    into a simulation settings dict for use in the wrapping Python start up
    (and process managers).
    """
    verbose = False # extra prints for testing
    silentInternal = False # Minimum prints about input/output files
    replaceInputFile = False # replace input file with output one?
    
#    ## find the number of available cpus
#    import multiprocessing
#    numCores = multiprocessing.cpu_count()
    
    # start return dict to be loaded up below
    returnDict = {}
    
    if silentInternal==False:
        print "Loading up Python settings dict from file: "+inputSettingsFile
    
    #loop through file and load up dict
    if os.path.exists(inputSettingsFile):
        fIn = open(inputSettingsFile,'r')
        lines = fIn.readlines()
        outLines = lines
        lineNum = -1
        fIn.close()
        for line in lines:
            lineNum+=1
                           
            if verbose:
                line = line.replace('\n','')
                print '\nOriginal line:'+line             
            if len(line)<=1: 
                if verbose:
                    print 'Blank line found'
            elif line[0]=='#':
                if verbose:
                    line = line.replace('\n','')
                    print 'Comment line found: '+line
            
            else:
                # parameter line
                try:
                    [key,val] = line.split('=')
                    # clean up val to make things easier below
                    val = strToStr(val)
                    valUse=val
                    if 'mcONLY'in key:
                        valUse=returnDict['mcONLY']=strToBool(val,True)
                        if verbose:
                            print 'mcONLY found to be = '+str(returnDict['mcONLY'])
                    elif 'numSamples_SimAnneal'in key:
                        try:
                            valUse = int(val)
                        except:
                            print "Value found in settings file for numSamples, '"+val+"', was invalid."
                            print "Using default value of 1e6."
                            valUse = 1e6
                        if (valUse<=0)or(valUse>1e10): 
                            print 'Value found in settings file for numSamples, '+val+\
                                ', was out of range [1,'+str(1e10)+'].'
                            print "Using default value of 1e6."  
                            valUse = 1e6
                        returnDict['numSamples_SimAnneal'] = valUse
                        if verbose:
                            print 'numSamples_SimAnneal found to be = '+str(returnDict['numSamples_SimAnneal'])
                    elif 'numSamples'in key:
                        try:
                            valUse = int(val)
                        except:
                            print "Value found in settings file for numSamples, '"+val+"', was invalid."
                            print "Using default value of 1e6."
                            valUse = 1e6
                        if (valUse<=0)or(valUse>1e10): 
                            print 'Value found in settings file for numSamples, '+val+\
                                ', was out of range [1,'+str(1e10)+'].'
                            print "Using default value of 1e6."  
                            valUse = 1e6
                        returnDict['numSamples'] = valUse
                        if verbose:
                            print 'numSamples found to be = '+str(returnDict['numSamples'])
                    elif 'useMultiProcessing'in key:
                        valUse=returnDict['useMultiProcessing'] = strToBool(val,False)
                        if verbose:
                            print 'useMultiProcessing found to be = '+str(returnDict['useMultiProcessing'])
                    elif 'silent'in key:
                        valUse=returnDict['silent'] = strToBool(val,True)
                        if verbose:
                            print 'silent found to be = '+str(returnDict['silent'])
                    elif 'verbose'in key:
                        returnDict['verbose'] = strToBool(val,False)
                        if verbose:
                            print 'verbose found to be = '+str(returnDict['verbose'])
                    elif 'CopyToDrobox'in key:
                        returnDict['CopyToDrobox'] = strToBool(val,False)
                        if verbose:
                            print 'CopyToDrobox found to be = '+str(returnDict['CopyToDrobox'])      
                    elif 'TcStepping'in key:
                        VALorig = strToBool(val,True)
                        VAL=False
                        if returnDict['DIonly']:
                            VAL = False
                        else:
                            VAL = VALorig
                        returnDict['TcStepping'] = VAL
                        if verbose:
                            print 'TcStepping found to be = '+str(returnDict['TcStepping'])
                            if VALorig!=VAL:
                                print 'TcStepping was initially = '+str(VALorig)+", but it was changed to False because DIonly was True."
                    elif 'outputData_dir'in key:
                        if len(val)<2:
                            default  = '/run/media/Kyle/Data1/Todai_Work/Data/data_Binary/data_Duo/'
                            print 'Value for data_dir, '+val+', was invalid.'
                            print 'Using default value of:'+default
                            valUse = default
                        returnDict['outputData_dir'] = valUse
                        if verbose:
                            print 'outputData_dir found to be = '+returnDict['outputData_dir']
                    elif 'outputData_filenameRoot'in key:
                        if len(val)<2:
                            default  = 'TempRootFilename.dat'
                            print 'Value for filenameRoot, '+val+', was invalid.'
                            print 'Using default value of:'+default
                            valUse = default
                        returnDict['outputData_filenameRoot'] = valUse
                        if verbose:
                            print 'outputData_filenameRoot found to be = '+returnDict['outputData_filenameRoot']
                    elif 'CalcBurnIn'in key:
                        valUse=returnDict['CalcBurnIn'] = strToBool(val,False)
                        if verbose:
                            print 'CalcBurnIn found to be = '+str(returnDict['CalcBurnIn'])
                    elif 'calcCorrLengths'in key:
                        returnDict['calcCorrLengths'] = strToBool(val,False)
                        if verbose:
                            print 'calcCorrLengths found to be = '+str(returnDict['calcCorrLengths'])
                    elif 'makeMCMCprogPlots'in key:
                        returnDict['makeMCMCprogPlots'] = strToBool(val,True)
                        if verbose:
                            print 'makeMCMCprogPlots found to be = '+str(returnDict['makeMCMCprogPlots'])
                    elif 'makeSimAnnealProgPlots'in key:
                        returnDict['makeSimAnnealProgPlots'] = strToBool(val,True)
                        if verbose:
                            print 'makeSimAnnealProgPlots found to be = '+str(returnDict['makeSimAnnealProgPlots'])
                    elif 'delChainsAfter'in key:
                        returnDict['delChainsAfter'] = strToBool(val,False)
                        if verbose:
                            print 'delChainsAfter found to be = '+str(returnDict['delChainsAfter'])
                    elif 'delCombinedDataAfter'in key:
                        returnDict['delCombinedDataAfter'] = strToBool(val,False)
                        if verbose:
                            print 'delCombinedDataAfter found to be = '+str(returnDict['delCombinedDataAfter'])
                    elif 'CalcGelmanRubin'in key:
                        VAL_orig = strToBool(val,True) 
                        if returnDict['useMultiProcessing']==False:
                            VAL = False
                        else:
                            VAL = VAL_orig
                        returnDict['CalcGelmanRubin'] = VAL       
                        if verbose:
                            print 'CalcGelmanRubin found to be = '+str(returnDict['CalcGelmanRubin'])
                            if VAL!=VAL_orig:
                                print 'useMultiProcessing==False, so CalcGelmanRubin changed from '+repr(VAL_orig)+" to "+repr(VAL)
                    elif 'makeOrbitPlots'in key:
                        returnDict['makeOrbitPlots'] = strToBool(val,True)
                        if verbose:
                            print 'makeOrbitPlots found to be = '+str(returnDict['makeOrbitPlots'])
                    elif 'makePosteriorsPlot'in key:
                        returnDict['makePosteriorsPlot'] = strToBool(val,True)
                        if verbose:
                            print 'makePosteriorsPlot found to be = '+str(returnDict['makePosteriorsPlot'])
                    elif 'startTemp'in key:
                        valUse=returnDict['startTemp'] = float(val)
                        if verbose:
                            print 'startTemp found to be = '+str(returnDict['startTemp'])
                    elif 'DIonly'in key:
                        valUse=returnDict['DIonly'] = strToBool(val,False)
                        if verbose:
                            print 'DIonly found to be = '+str(returnDict['DIonly'])
                    elif 'RVonly'in key:
                        valUse=returnDict['RVonly'] = strToBool(val,False)
                        if verbose:
                            print 'RVonly found to be = '+str(returnDict['RVonly'])
                    elif 'simAnneal'in key:
                        valUse=returnDict['simAnneal'] = strToBool(val,False)
                        if verbose:
                            print 'simAnneal found to be = '+str(returnDict['simAnneal'])
                    elif 'loopedMCMC'in key:
                        valUse=returnDict['loopedMCMC'] = strToBool(val,False)
                        if verbose:
                            print 'loopedMCMC found to be = '+str(returnDict['loopedMCMC'])
                    elif 'simulate_StarStar'in key:
                        valUse=returnDict['simulate_StarStar'] = strToBool(val,False)
                        if verbose:
                            print 'simulate_StarStar found to be = '+str(returnDict['simulate_StarStar'])
                    elif 'simulate_StarPlanet'in key:
                        valUse=returnDict['simulate_StarPlanet'] = strToBool(val,False)
                        if verbose:
                            print 'simulate_StarPlanet found to be = '+str(returnDict['simulate_StarPlanet'])
                    elif 'inclination_degMAX'in key:
                        valUse=returnDict['inclination_degMAX'] = float(val)
                        if verbose:
                            print 'inclination_degMAX found to be = '+str(returnDict['inclination_degMAX'])
                    elif 'periodMAX'in key:
                        valUse=returnDict['periodMAX'] = float(val)
                        if verbose:
                            print 'periodMAX found to be = '+str(returnDict['periodMAX'])
                    elif 'longAN_degMAX'in key:
                        valUse=returnDict['longAN_degMAX'] = float(val)
                        if verbose:
                            print 'longAN_degMAX found to be = '+str(returnDict['longAN_degMAX'])
                    # create output version of this line to write to output version of settings file
                    lineOut = key+'='+str(valUse)+'\n'
                    outLines[lineNum]=lineOut
                except:
                    print 'failed to split line:'+line
                    print '[key,val] = '+repr(line.split('='))
                    print 'failed line had length:'+str(len(line))
    else:
        print 'ERROR: Settings file, '+inputSettingsFile+', does NOT exist!'
        
    if outputSettingsFile!="":
        # write output version of lines to an output settings file
        fOut = open(outputSettingsFile,'w')
        fOut.writelines(outLines)
        fOut.close()
        
        if silentInternal==False:
            print 'Output settings file written to: '+outputSettingsFile
        
        # Replace input file with output one with new default values and such
        if (os.path.exists(inputSettingsFile) and os.path.exists(outputSettingsFile) and replaceInputFile):
            # both are on disk, so kill input one and rename output one to that name
            os.remove(inputSettingsFile)
            os.rename(outputSettingsFile,inputSettingsFile)
            if verbose:
                print 'Input settings file: '+inputSettingsFile +' was deleted from disk and '\
                +'output file was renamed to be the same as the input settings file.'
    
    return returnDict

def sysDataToDict(filename):
    """
    Extract needed data of the binary system from the general C++/Python data file 
    in the 'SimSettings_and_InputData' folder for Duo simulations.
    """
    
    verbose = False # extra prints for testing
    silentInternal = False
    
    # start return dict to be loaded up below
    returnDict = {}
    
    if silentInternal==False:
        print "Loading up Python System data dict from file: "+filename
    
    #loop through file and load up dict
    if os.path.exists(filename):
        fIn = open(filename,'r')
        lines = fIn.readlines()
        fIn.close()
        for line in lines:
            if verbose:
                line = line.replace('\n','')
                print '\nOriginal line:'+line
            if len(line)<=1: 
                if verbose:
                    print 'Blank line found'
            elif line[0]=='#':
                if verbose:
                    line = line.replace('\n','')
                    print 'Comment line found: '+line
            else:
                # parameter line
                try:
                    [key,val] = line.split('=')
                    # clean up val to make things easier below
                    val = strToStr(val)
                    valUse=val
                    #load up general system data and primary star mass
                    if ('Sys_Dist_PC'in key) and ('error'not in key):
                        valUse=returnDict['Sys_Dist_PC']=float(val)
                        if verbose:
                            print 'Sys_Dist_PC found to be = '+str(returnDict['Sys_Dist_PC'])
                    elif ('Mass1'in key) and ('error' not in key):
                        valUse=returnDict['Mass1']=float(val)
                        if verbose:
                            print 'Mass1 found to be = '+str(returnDict['Mass1'])
                    # load up companion planet data
                    elif ('planet_K'in key) and ('error' not in key):
                        valUse=returnDict['planet_K']=float(val)
                        if verbose:
                            print 'planet_K found to be = '+str(returnDict['planet_K'])
                    elif ('planet_e'in key) and ('error' not in key):
                        valUse=returnDict['planet_e']=float(val)
                        if verbose:
                            print 'planet_e found to be = '+str(returnDict['planet_e'])
                    elif ('planet_T'in key) and ('Tc' not in key) and ('error' not in key):
                        valUse=returnDict['planet_T']=float(val)
                        if verbose:
                            print 'planet_T found to be = '+str(returnDict['planet_T'])
                    elif ('planet_Tc'in key) and ('error' not in key):
                        valUse=returnDict['planet_Tc']=float(val)
                        if verbose:
                            print 'planet_Tc found to be = '+str(returnDict['planet_Tc'])
                    elif ('planet_P'in key) and ('error' not in key):
                        valUse=returnDict['planet_P']=float(val)
                        if verbose:
                            print 'planet_P found to be = '+str(returnDict['planet_P'])
                    elif ('planet_MsinI'in key) and ('error' not in key):
                        valUse=returnDict['planet_MsinI']=float(val)
                        if verbose:
                            print 'planet_MsinI found to be = '+str(returnDict['planet_MsinI'])
                    elif ('planet_argPeri'in key) and ('error' not in key):
                        valUse=returnDict['planet_argPeri']=float(val)
                        if verbose:
                            print 'planet_argPeri found to be = '+str(returnDict['planet_argPeri'])
                    # load up companion star data
                    elif ('star_e'in key) and ('error' not in key):
                        valUse=returnDict['star_e']=float(val)
                        if verbose:
                            print 'star_e found to be = '+str(returnDict['star_e'])
                    elif ('star_Tc'in key):
                        valUse=returnDict['star_Tc']=float(val)
                        if verbose:
                            print 'star_Tc found to be = '+str(returnDict['star_Tc'])
                    elif (('star_T'in key) and ('error' not in key))and('Tc' not in key):
                        valUse=returnDict['star_T']=float(val)
                        if verbose:
                            print 'star_T found to be = '+str(returnDict['star_T'])
                    elif ('star_P'in key) and ('error' not in key):
                        valUse=returnDict['star_P']=float(val)
                        if verbose:
                            print 'star_P found to be = '+str(returnDict['star_P'])
                    elif ('star_Mass2'in key) and ('error' not in key):
                        valUse=returnDict['star_Mass2']=float(val)
                        if verbose:
                            print 'star_Mass2 found to be = '+str(returnDict['star_Mass2'])
                    elif ('star_argPeri'in key) and ('error' not in key):
                        valUse=returnDict['star_argPeri']=float(val)
                        if verbose:
                            print 'star_argPeri found to be = '+str(returnDict['star_argPeri'])
                    elif ('star_inc'in key) and ('error' not in key):
                        valUse=returnDict['star_inc']=float(val)
                        if verbose:
                            print 'star_inc found to be = '+str(returnDict['star_inc'])
                            
                except:
                    print 'failed to split line:'+line
                    print 'failed line had length:'+str(len(line))
                    
    return returnDict

def strToStr(strIn):
    """
    To clean up a input string.  ie. remove any ", ', or \n values in the string
    """
    strOut = strIn
    strOut = strOut.replace('"','')
    strOut = strOut.replace("'",'')
    strOut = strOut.replace('\n','')
    strOut = strOut.replace(' ','')
    return strOut

def strToBool(strIn, default):
    """
    Convert a string of with 'true','True','False' or 'false' in it into a
    bool accordingly.
    """
    if ('true'in strIn)or('True'in strIn):
        b = True
    elif ('false'in strIn)or('False'in strIn):
        b = False
    else:
        b = default
    return b
    
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

def findArrayMin(ary):
    """
    Find the minimum value inside an array up to 4D.
    """
    min = 1e9
    
    if type(ary)==list:
        for subEl in ary:
            if type(subEl)==list:
                for subsubEl in subEl:
                    if type(subsubEl)==list:
                        for subsubsubEl in subsubEl:
                            if min>subsubsubEl:
                                min = subsubsubEl
                                #print 'new min = ',min
                    else:
                        if min>subsubEl:
                            min = subsubEl
                            #print 'new min = ',min
            else:
                if min>subEl:
                    min = subEl
                    #print 'new min = ',min
    else:
        print 'NOTE: What was provided to findArrayMin was not an array, so just returning value.'

    return min
        
def findArrayMax(ary):
    """
    Find the minimum value inside an array up to 4D.
    """
    
    max = -1e9
    
    if type(ary)==list:
        for subEl in ary:
            if type(subEl)==list:
                for subsubEl in subEl:
                    if type(subsubEl)==list:
                        for subsubsubEl in subsubEl:
                            if max<subsubsubEl:
                                max = subsubsubEl
                                #print 'new max = ',max
                    else:
                        if max<subsubEl:
                            max = subsubEl
                            #print 'new max = ',max
            else:
                if max<subEl:
                    max = subEl
                    #print 'new max = ',max
    else:
        print 'NOTE: What was provided to findArrayMin was not an array, so just returning value.'
    
    return max

def findNuFromLog(logFilename = ''):
    """
    This will pull the value for 'nu' (the number of degrees of freedom) out of the log 
    file written by C++ as that is where nu was calculated and written.  Saves time calculating 
    it again in Python.
    """
    if True:
        print "\n\nFinding nu from log file: "+logFilename+'\n\n'
    log = open(logFilename,'r')
    lines = log.readlines()
    log.close()
    i = 0
    for line in lines:
        i=i+1
        if line.find("one_over_nu_DI")>=0:
            ss = line.split('=')
            onOverNu = float(ss[-1])
            nuDI = 1.0/onOverNu
        elif line.find("one_over_nu_RV")>=0:
            ss = line.split('=')
            onOverNu = float(ss[-1])
            nuRV = 1.0/onOverNu
        elif line.find("one_over_nu_TOTAL")>=0:
            ss = line.split('=')
            onOverNu = float(ss[-1])
            nu = 1.0/onOverNu
            # we can break loop after this is found as it is the 3rd one printed and  
            # would thus be the last found of the 3 so done looping after it is found
            break
        
    printStr =  "Found Nu = "+str(nu)+", NuRV = "+str(nuRV)+", and NuDI = "+str(nuDI)+" on line number "+str(i)
    if False:
        print printStr
    return [nu,nuRV,nuDI,printStr] 

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



def likelihoodsCalc(chiSquareds, nu=1):
    """
    This function will convert the chiSqureds into likelihoods following 
    L =e^(-nu*chiSquared/2)
    """
    
    likelihoods = []
    for chiSquared in chiSquareds:
        likelihood = math.exp(-nu*chiSquared/2.0)
        likelihoods.append(likelihood)
    
    return likelihoods

def listRandomizer(listIN, listIN2=False):
    """
    Returns randomized version of a list.  If two lists pas
    """
    import random as rand
    
    if listIN2==False:
        listIN2 = listIN
    else:
        if len(listIN)!=len(listIN2):
            print 'listRandomizer: ERROR! input lists are different sizes!!'
    
    listOUT = []
    listOUT2 = []
    numValsLeft = len(listIN)
    locationsCurr = range(0,len(listIN))
    locationsOUT = []
    
    while numValsLeft>=1:
        location_proposed = int(rand.uniform(0, len(listIN)))
        if False:
            print 'Current number of values left in input array = ',numValsLeft
        if location_proposed in locationsCurr:
            location = location_proposed
            locationsCurr[location] = -1
            locationsOUT.append(location)
            numValsLeft  = numValsLeft-1    
    
    for loc in locationsOUT:
        listOUT.append(listIN[loc])    
        listOUT2.append(listIN2[loc])
    
    if len(listIN)!=len(listOUT):
        print 'ERROR: length of input list '+str(len(listIN))+' != len of output list '+str(len(listOUT))
        print 'locations remaining array: '+repr(locationsCurr)
    
    if listIN2==False:
        return listOUT
    else:
        return listOUT,listOUT2

def TAcalculator2(t,e, T, period, T_center=0, verbose=False, debug=False):
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
    period_days = period*365.242
    timeDiff_days = (t- T)-int((t-T)/period_days)*period_days 
    if timeDiff_days<0.0:
        timeDiff_days = timeDiff_days+period_days
    phase = 0.0
    if T_center!=0.0:
        phaseDiff_days = (T_center-T)-int((T_center -T)/period_days)*period_days 
        if T>T_center:
            phaseDiff_days = phaseDiff_days+period_days
        phase = phaseDiff_days/period_days
    if verbose:
        print "Unitless phase calculated to be "+str(phase)+", using T_center = "+str(T_center)+" and To = "+str(T)
        print 'timeDiff_days = '+str(timeDiff_days)+', phaseDiff_days = '+str(phaseDiff_days)+', period_days = '+str(period_days)
    
    M = n*(((timeDiff_days)/365.242)+phase)#+(phase*2.0*pi)
    ## Push M value into 0-2pi range ######
    numCirclesBiggerD = abs(M/(2.0*pi))
    numCirclesBiggerI = int(numCirclesBiggerD)
    if numCirclesBiggerI<1:
        numCirclesBiggerI = 1
    if (M<0.0):
        M_out = (numCirclesBiggerI+1)*2.0*pi + M
        if verbose:
            print "M updated from, "+str(M)+" to "+str(M_out)+", numCirclesBiggerI = "+str(numCirclesBiggerI)+", numCirclesBiggerD = "+str(numCirclesBiggerD)
        
    elif M>(2.0*pi):
        M_out = M-numCirclesBiggerI*2.0*pi
        if verbose:
            print "M updated from, "+str(M)+" to "+str(M_out)+", numCirclesBiggerI = "+str(numCirclesBiggerI)+", numCirclesBiggerD = "+str(numCirclesBiggerD)
    else:
        M_out = M
    M = M_out
    
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
    while (abs(E_last-E_latest) > (1.0e-10))and(count<50):
        if debug:
            print 'current E [rad]= ', E_latest
        E_last = E_latest
        M_last = E_last - e*math.sin(E_last)
        E_latest = E_last - ((M_last-M)/(1.0-e*math.cos(E_last)))
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
    top = (math.cos(E_latest)-e)
    btm = (1.0-e*math.cos(E_latest))
    TA_rad  = math.acos(top/btm) 
    
        
#    if E_latest<0.0:
#        #print 'top is negative'
#        TA_rad = -1.0*TA_rad
#    if False:
#        print "\nTA_rad calc failed!"
#        print "E_latest = ",E_latest
#        print "e = ",e
#    
#    if E_latest>(2.0*pi):
#        # convert E to be inside one orbit (ie. under 2*PI)
#        numCirclesD = E_latest/(2.0*pi)
#        numCirclesI = int(numCirclesD)
#        E_latest_oneCircle  = E_latest-numCirclesI*2.0*pi
#        if verbose: 
#            print "E_latest found to be "+str(numCirclesI)+" times over 2pi, so made E_latest_circle = "+str(math.degrees(E_latest_oneCircle))
#    else:
#        E_latest_oneCircle  = E_latest

    if (E_latest>pi) or (E_latest<0):
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
        
    #TA_deg = math.degrees(TA_rad)#$$$$$$$$$$$$$$$$
    #print 'epoch = '+str(t)+', T = '+str(T)+', timeDiff_days = '+str(timeDiff_days)+', numPeriodDiff = '+str(int((t-T)/period_days))+', TA_deg = '+str(TA_deg)#$$$$$$$$$$$$$$$$$$$$$$$$
#    ## Calculate TA in another way    
#    x = ((1.0-e**2.0)**(1.0/2.0))*math.cos(E_latest/2.0)
#    y = ((1.0+e**2.0)**(1.0/2.0))*math.sin(E_latest/2.0)
#    TA_rad2 = 2.0*math.atan2(y, x)
#    #print 'TA_2 = '+str(math.degrees(TA_rad2))
#    
#    print 'TA = ',math.degrees(TA_rad)
#    print 'TA2 = ',math.degrees(TA_rad2)
    
#    if True:
#        if t<2452381:
#            print "epoch = "+str(t)+", To = "+str(T)+", Tc = "+str(T_center)+", e = "+str(e)+", period = "+str(period)+ ", E = "+str(math.degrees(E_latest))+", TA = "+str(math.degrees(TA_rad))
    
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
        
def makeArtificialData(longAN_deg, e, T, Tc, period, inc, argPeri_deg, a_total, sys_dist, Mass1, Mass2, numDataPoints=100):    
    
    epochs = []
    SAs = []
    SA_errors = []
    PAs = []
    PA_errors = []
    VRs = []
    VR_errors = []
    
    if False:
        argPeriDI=argPeri_deg+180.0
    else:
        argPeriDI =argPeri_deg
        
    #rand.uniform(T-(period*365.242*0.5),T+(period*365.242*0.5))
    delta_epoch = (period*365.242)/numDataPoints
    for i in range(0,numDataPoints):
        SA_errors.append(abs(rand.normalvariate(0.015,0.005)))
        PA_errors.append(abs(rand.normalvariate(0.5,0.05)))
        VR_errors.append(abs(rand.normalvariate(15.0,1)))
        #epochs.append(T+(i*delta_epoch))
        epochs.append(rand.uniform(T,T+(period*365.242)))
    
    (a_total, a1IN, a2, p_s) = semiMajorConverter(Mass1, Mass2, a_total,a1=0.0,a2=0.0, period=period,verbose=False)
    for epoch in range(0,numDataPoints):
        (v_r,K) = vrCalculatorStar2(epochs[epoch],e,T,period,argPeri_deg,a1IN,T_center=Tc,i=inc, K=False, verbose=False)
        
        (n, M_deg, E_latest_deg, TA_deg, Sep_Dist_AU_OP, SA, PA, a1, a2) =\
                    orbitCalculator2(epochs[epoch], sys_dist, inc, longAN_deg, e, T, period, argPeriDI, a_total,\
                                                                Mass1=1, Mass2=1, verbose=False)
        #sigma = 1.80
        #SAs.append(SA+SA_errors[epoch]*rand.uniform(-sigma,sigma))
        #PAs.append(PA+PA_errors[epoch]*rand.uniform(-sigma,sigma))
        #VRs.append(v_r+VR_errors[epoch]*rand.uniform(-sigma,sigma))
        SAs.append(rand.normalvariate(SA,0.0341))
        PAs.append(rand.normalvariate(PA,0.181))
        VRs.append(rand.normalvariate(v_r,14.69))

    if True:
        print '\n\nDI data output:\n'
        for epoch in range(0,numDataPoints):
             print str(epochs[epoch])+'    '+str(PAs[epoch])+"    "+str(PA_errors[epoch])+"    "+str(SAs[epoch])+"    "+str(SA_errors[epoch])
        
        print '\n\nRV data output:\n'
        for epoch in range(0,numDataPoints):
             print str(epochs[epoch])+'    '+str(VRs[epoch])+"    "+str(VR_errors[epoch])
         
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
    (a_total, a1, a2) = semiMajorConverter(Mass1, Mass2, a_total=a_total_proposed,a1=0.0,a2=0.0)
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
    (a_total, a1, a2) = semiMajorConverter(Mass1, Mass2, a_total=a_total_proposed,a1=0.0,a2=0.0)
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
    (a_total, a1, a2) = semiMajorConverter(Mass1, Mass2, a_total=a_total_proposed,a1=0.0,a2=0.0)
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
        
        
        