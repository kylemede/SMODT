#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
#import numpy as np
import exoSOFTlogger
import cppTools
import constants as const
import copy
import os
import shutil
import timeit
import glob
import numpy as np
import sys
import pyfits
import warnings
import readWriteTools as rwTools
warnings.simplefilter("error")

log = exoSOFTlogger.getLogger('main.genTools',lvl=100,addFH=False)  
    
def mcmcEffPtsCalc(outputDataFilename):
    """
    Calculate average correlation length and the number of effective steps for each parameter 
    that was varying during the simulation.  The results are put into the log.
    
    Using Correlation Length based on Tegmark2004.  This is the number of steps until the 
    variance (aka correlation) is half that of the total chain.  We perform this in a "box car"
    style that starts calculating it again on the step just after the previously found correlation 
    length step.  This way it produces an average value that is more reliable.
    """
    log.info("Starting to calculate correlation lengths")
    (head,data) = rwTools.loadFits(outputDataFilename)
    numSteps = data.shape[0]
    (paramList,paramStrs,paramFileStrs) = getParStrs(head,latex=False)
    completeStr=""
    try:
        ## Make post tools cpp obj
        PostCTools = cppTools.PostCtools()
        dataC = []
        for i in range(0,numSteps):
            dataC.append(data[i,:])
        dataC=np.array(dataC,dtype=np.dtype('d'),order='C')
        PostCTools.loadParamData(dataC)
        
        completeStr+= '\n'+'-'*45+'\nThe mean correlation lengths of all params are:\n'+'-'*45+'\nparam #, param name, mean correlation length'
        completeStr+= ' -> total # of steps/mean correlation length = number of effective points\n'
        for i in range(0,len(paramList)):
            log.debug( "*"*60+"\n"+'starting to mean calculate corr length for '+paramStrs[i]+' with CPP')
            #Call CPP tool to calculate average correlation length.
            meanCorrLength = PostCTools.corrLenCalc(paramList[i])
            currParamStr = str(paramList[i])+', '+paramStrs[i]+", "+str(meanCorrLength)
            currParamStr+=    ' -> '+str(numSteps)+'/'+str(meanCorrLength)+' = '+str(numSteps/meanCorrLength)+'\n'
            completeStr+=currParamStr
            log.debug(currParamStr)
        log.debug(completeStr)
    except:
        log.error("An error occurred when trying to calculate the correlation lengths, # effective steps")
    return completeStr

def burnInCalc(mcmcFnames,combinedFname):
    """
    NOTE: ExoSOFT was designed to start the full MCMC chain from the last point of the 
        Sigma Tuning stage.  As this stage effectively acts as a form of burn-in period
        the burn-in value found from the pure MCMC tends to be very short.
 
    Calculate the burn in for a set of MCMC chains following the formulation of Tegmark.
     
    Burn-in is defined as the first point in a chain where the likelihood is greater than 
    the median value of all the chains.  Thus, there MUST be more than 1 chain to perform this calculation.
    """
    log.info("Starting to calculate burn-in.")
     
    chiSquaredsALL = np.array([])
    burnInLengths = []
    # calculate median of combined data ary
    (head0,data0) = rwTools.loadFits(combinedFname)
    #nu = float(head0['NU'])
    chiSqsALL = data0[:,11]
    if type(chiSqsALL)!=np.ndarray:
        chiSqsALL = np.array(chiSqsALL)
    likelihoodsALL = np.exp(-chiSqsALL/2.0)
    log.debug("likelihoodsALL min = "+repr(np.min(likelihoodsALL)))
    log.debug("likelihoodsALL max = "+repr(np.max(likelihoodsALL)))
    medainALL = np.median(likelihoodsALL)         
    log.debug("medainALL = "+str(medainALL))
    s =21*'-'+'\nBurn-In lengths were:\n'+21*'-'
    s+='\nmedian value for all chains = '+str(medainALL)
    ## calculate location of medianALL in each chain
    for filename in mcmcFnames:
        if os.path.exists(filename):
            (head,data) = rwTools.loadFits(filename)
            chiSqsChain = data[:,11]
            likelihoodsChain = np.exp(-chiSqsChain/2.0)
            #medianChain = np.median(chiSquaredsChain)
            burnInLength = len(likelihoodsChain)
            i=0
            while i<(len(likelihoodsChain)-1):
                i+=1
                if likelihoodsChain[i]>medainALL:
                    #print 'chiSqsChain[i] = '+str(chiSqsChain[i])
                    burnInLength = i+1
                    break
            burnInLengths.append(burnInLength)
            s2 = "\nfor chain #"+str(head['chainNum'])
            s2 += "\nTotal number of points in the chain = "+str(len(chiSqsChain))+"\n"
            s2 += "Burn-in length = "+str(burnInLength)+"\n"
            s+=s2
            log.debug(s2)
    log.debug(s)
    return (s,burnInLengths)

def burnInStripper(mcmcFnames,burnInLengths):
    """
    Strip the initial burn-in off each chain.
    """
    newFnames=[]
    for i in range(0,len(mcmcFnames)):
        filename = mcmcFnames[i]
        burnIn = burnInLengths[i]
        if os.path.exists(filename):
            (head,data) = rwTools.loadFits(filename)
            ##strip burn-in and write to new fits             
            hdu = pyfits.PrimaryHDU(data[burnIn:,:])
            hdulist = pyfits.HDUList([hdu])
            newHead = hdulist[0].header
            for key in head:
                newHead[key]=(head[key],head.comments[key])
            n = os.path.splitext(os.path.basename(filename))
            newFname = os.path.join(os.path.dirname(mcmcFnames[0]),n[0]+'_BIstripped.fits')
            hdulist.writeto(newFname)
            log.info("output file written to:below\n"+newFname)
            hdulist.close()
            newFnames.append(newFname)
            log.debug("burn-in stripped file written to:\n"+newFname)
    return newFnames

def gelmanRubinCalc(mcmcFileList,nMCMCsamp=1):
    """
    Calculate Gelman-Rubin statistic for each varying param.
    Input MUST be the list of more than one MCMC chain files.
    """
    GRs=[]
    Ts = []
    grStr = '\n'+'-'*30+"\nGelman-Rubin Results:\n"+'-'*30+'\n'
    try:
        Lcfloat = float(nMCMCsamp)
        if os.path.exists(mcmcFileList[0]):
            log.info("Starting to calculate R&T")
            ###########################################################
            ## stage 1 ->  load up values for each param in each chain.
            ## allStg1vals = [chain#, param#, (mean,variance,Lc)]
            ## stage 2 ->  Use them to compare between chains 
            ##             then calc R and T.
            ###########################################################
            (head,data) = rwTools.loadFits(mcmcFileList[0])
            (paramList,paramStrs,paramFileStrs) = getParStrs(head,latex=False)
            
            Nc = len(mcmcFileList)
            ##start stage 1
            allStg1vals=np.zeros((Nc,len(paramList),3))
            for i in range(0,len(mcmcFileList)):
                log.debug("Starting to calc chain #"+str(i)+' GR values')
                (head,data) = rwTools.loadFits(mcmcFileList[i])
                allStg1vals[i,:,2]=data.shape[0]
                for j in range(0,len(paramList)):
                    log.debug("calculating stage 1 of GR for chain #"+str(i)+", param: "+paramStrs[j])
                    allStg1vals[i,j,0]=np.mean(data[:,paramList[j]])
                    allStg1vals[i,j,1]=np.var(data[:,paramList[j]])
            ##start stage 2         
            rHighest = 0
            tLowest = 1e9       
            rHighStr = ''
            tLowStr = ''
            for j in range(0,len(paramList)):
                try:
                    log.debug("Starting stage 2 for param: ("+str(paramList[j])+"/"+str(len(paramList))+"), "+paramStrs[j])
                    Ncfloat = float(Nc)
                    ##calc R
                    W = 0
                    for i in range(0,Nc):
                        #Lcfloat = float(allStg1vals[i,j,2])
                        W+=(Lcfloat/(Lcfloat-1.0))*allStg1vals[i,j,1]
                    W=W/Ncfloat
                    V = np.mean(allStg1vals[:,j,1])+(Ncfloat/(Ncfloat-1.0))*np.var(allStg1vals[:,j,0])
                    R=np.NaN
                    if W!=0:
                        R = np.sqrt(V/W)
                    GRs.append(R)
                    ##calc T
                    #B = (Ncfloat/(Ncfloat-1.0))*np.var(allStg1vals[:,j,0]*allStg1vals[:,j,2])
                    B = (Ncfloat/(Ncfloat-1.0))*np.var(allStg1vals[:,j,0])*Lcfloat
                    #Uses the mean of Lc values
                    T = np.NAN
                    if B!=0:
                        #T = np.mean(allStg1vals[:,j,2])*Ncfloat*np.min([(V/B),1.0])
                        T = Lcfloat*Ncfloat*np.min([(V/B),1.0])
                    Ts.append(T)       
                    grStr+=paramStrs[j]+" had R = "+str(R)+", T = "+str(T)+'\n'
                    if T<tLowest:
                        tLowest=T
                        tLowStr="Lowest T = "+str(T)+' for '+paramStrs[j]+'\n'
                    if R>rHighest:
                        rHighest=R
                        rHighStr="Highest R = "+str(R)+' for '+paramStrs[j]+'\n'
                except:
                    log.critical("an error occured while performing stage 2 of GR calc for param "+paramStrs[j])
            grStr+='\nWorst R&T values were:\n'+rHighStr+tLowStr+'\n'
        else:
            log.critical("Gelman-Rubin stat can NOT be calculated as file does not exist!!:\n"+chainDataFileList[0])
    except:
        log.critical("Gelman-Rubin stat FAILED to be calculated for some reason")
    return (GRs,Ts,grStr)

def timeStrMaker(deltaT):
    """
    Convert a time in seconds into a nicer string.
    """
    timeStr = ''
    if deltaT>60:
        if deltaT>(60*60):
            hr = int(deltaT//(60*60))
            min = int((deltaT-hr*60*60)/60.0)
            timeStr = str(hr)+" hours and "+str(min)+" minutes"
        else:
            min = int(deltaT//(60))
            timeStr = str(min)+" minutes and "+str(int(deltaT-min*60))+" seconds"
    else:
        timeStr = str(int(deltaT))+' seconds'
    return timeStr

def getParStrs(head,latex=True,getALLpars=False):
    """
    Return matching paramList, paramStrs, paramFileStrs for provided header.
    
    latex=True will return latex formated strs for use in Python code.
    getALLpars=True will signal to get all params, not just the ones that were varying.
    """
    paramList = getParInts(head)    
    paramFileStrs = ['m1','m2','parallax','Omega','e','To', 'Tc','P','i','omega','a_total','chiSquared','K']
    paramStrs = ['m1 [Msun]','m2 [Msun]','Parallax [mas]','Omega [deg]','e','To [JD]', 'Tc [JD]','P [Yrs]','i [deg]','omega [deg]','a_total [AU]','chiSquared','K [m/s]']
    if latex:
        paramStrs = [r'$m_1{\rm [M}_{\odot}{\rm ]}$',r'$m_2{\rm [M}_{\odot}{\rm ]}$',r'$\varpi{\rm [mas]}$',r'$\Omega{\rm [deg]}$',r'$e$',r'$T_o{\rm  [JD]}$', r'$T_c{\rm  [JD]}$',r'$P{\rm  [Yrs]}$',r'$i{\rm  [deg]}$',r'$\omega{\rm  [deg]}$',r'$a_{{\rm total}} {\rm [AU]}$',r'$\chi^2$',r'$K{\rm  [m/s]}$']
    if head["nRVdsets"]>0:
        for dataset in range(1,head["nRVdsets"]+1):
            paramFileStrs.append('offset_'+str(dataset))
            if latex:
                paramStrs.append(r"$\gamma_{{\rm "+str(dataset)+"}}{\\rm [m/s]}$")
            else:
                paramStrs.append('offset '+str(dataset)+' [m/s]')        
    ## clean up lists if not returning ALL   
    ## else set paramList to contain its for ALL params
    if (len(paramFileStrs)>len(paramList))and(getALLpars==False):
            paramStrsUse = []
            paramFileStrsUse = []
            for par in paramList:
                paramStrsUse.append(paramStrs[par])
                paramFileStrsUse.append(paramFileStrs[par])
            paramStrs = paramStrsUse
            paramFileStrs = paramFileStrsUse 
    elif getALLpars:
        paramList=np.arange(0,len(paramStrs))
    return (paramList,paramStrs,paramFileStrs)
        
def getSimpleDictVal(dict,key):
    
    """
    Get the value for a key in the settings dictionary.
    This will handle the values that are tuples and not
    returning the value.
    """
    try:
        if type(dict[key])==tuple:
            return dict[key][0]
        else:
            return dict[key]
    except:
        return False
    
def cleanUp(settingsDict,stageList,allFname):
    """
    Clean up final directory after simulation completes
    """
    #os.mkdir(settingsDict['finalFolder'])
    
    delFiles = []
    
    fnames = glob.glob(os.path.join(settingsDict['finalFolder'],"pklTemp-*"))
    for i in range(0,len(fnames)):
        delFiles.append(fnames[i])
    ##get chain data filenames to delete
    if settingsDict["delChains"]:
        for stage in stageList:
            fnames = glob.glob(os.path.join(settingsDict['finalFolder'],"outputData"+stage+"*.fits"))
            for i in range(0,len(fnames)):
                delFiles.append(fnames[i])
    ##get combined data filename to delete
    if settingsDict["delCombined"]:
        delFiles.append(allFname)
    if (settingsDict['rmBurn'][0])and(settingsDict['nChains'][0]>1):
        ##the burn-in was stripped in the final file, so kill the non-stripped version if it exists
        nm = os.path.join(os.path.dirname(allFname),'combinedMCMCdata.fits')
        if os.path.exists(nm):
            delFiles.append(nm)
            
    ##try to delete files
    for fname in delFiles:
        try:
            log.debug('Deleting file: '+os.path.basename(fname))
            os.remove(fname) 
        except:
            log.error('Failed to delete file: '+os.path.basename(fname))


    
def summaryFile(settingsDict,stageList,finalFits,clStr,burnInStr,bestFit,grStr,effPtsStr,allTime,postTime,durationStrings):
    """
    Make a txt file that summarizes the results nicely.
    """
    summaryFname = os.path.join(settingsDict['finalFolder'],'RESULTS.txt')
    if os.path.exists(summaryFname):
        f = open(summaryFname,'a')
    else:
        f = open(summaryFname,'w')
    (head,data) = rwTools.loadFits(finalFits)
    totalSamps = head['NSAMPLES']
    (paramList,paramStrs,paramFileStrs) = getParStrs(head,latex=False,getALLpars=True)
    (paramListCleaned,paramStrsCleaned,paramFileStrsCleaned) = getParStrs(head,latex=False)
    f.write("\n"+"*"*80+"\noutRoot:  "+settingsDict['outRoot']+"\n"+"*"*80+"\n")
    f.write('\n'+'-'*7+'\nBasics:\n'+'-'*7)
    f.write('\nparamList:\n'+repr(paramListCleaned))
    f.write('\nparamStrs:\n'+repr(paramStrsCleaned))
    f.write('\nparamFileStrs:\n'+repr(paramFileStrsCleaned))
    try:
        ## try to make and write the more advanced summary strings to the file
        nusStr = "\nnu values were: [total,DI,RV] = ["+str(head['NU'])+", "+str(head['NUDI'])+", "+str(head['NURV'])+"]\n"
        f.write(nusStr)
        stgNsampStrDict = {"MC":"nSamples","SA":"nSAsamp","ST":"nSTsamp","MCMC":"nSamples"}
        numFilesStr = '\nTotal # of files that finished each stage were:\n'
        chiSquaredsStr = '\nBest Reduced Chi Squareds for each stage were:\n'
        for stage in stageList:
            fnames = np.sort(glob.glob(os.path.join(settingsDict['finalFolder'],"outputData"+stage+"*.fits")))
            if (stage=="MCMC")and settingsDict["rmBurn"][0]:
                fnames = np.sort(glob.glob(os.path.join(settingsDict['finalFolder'],"outputData"+stage+"*BIstripped.fits")))
            numFilesStr+=stage+' = '+str(len(fnames))+", each with "+str(settingsDict[stgNsampStrDict[stage]][0])+" samples\n"
            if len(fnames)>0:
                chiSquaredsStr+=stage+" = ["
                for fname in fnames: 
                    try:
                        bestFit2 = findBestOrbit(fname,bestToFile=False,findAgain=True)
                        chiSquaredsStr+=str(bestFit2[11]/float(head['NU']))+', '
                        #print fname
                        #print str(bestFit2[11])+"/"+str(float(head['NU']))
                    except:
                        log.error("A problem occurred while trying to find best fit of:\n"+fname)
                chiSquaredsStr = chiSquaredsStr[:-2]+']\n'
        numFilesStr+="\n"+"*"*65+"\nThe final combined file was for a total of "+str(totalSamps)+" samples\n"+"*"*65+'\n'
        f.write(numFilesStr)
        f.write(chiSquaredsStr)
        bestStr = '\n'+'-'*21+'\nBest fit values were:\n'+'-'*21+'\n'
        ############################################
        ## calculate chi squareds for the best fit #
        ############################################
        ##get the real data
        realData = rwTools.loadRealData(os.path.join(settingsDict['settingsDir'],settingsDict['prepend']),dataMode=settingsDict['dataMode'][0])
        ## Make Orbit cpp obj
        Orbit = cppTools.Orbit()
        try:
            pasa = settingsDict["pasa"][0]
        except:
            pasa = False
        Orbit.loadStaticVars(settingsDict['omegaFdi'][0],settingsDict['omegaFrv'][0],settingsDict['lowEcc'][0],pasa)
        Orbit.loadConstants(const.Grav,const.pi,const.KGperMsun, const.daysPerYear,const.secPerYear,const.MperAU)
        ## ensure bestFit are in required format for Orbit
        params = []
        for par in bestFit:
            params.append(par)
        params=np.array(params,dtype=np.dtype('d'),order='C')
        Orbit.loadRealData(realData)
        modelData = np.ones((realData.shape[0],3),dtype=np.dtype('d'),order='C')
        Orbit.calculate(modelData,params)
        (raw3D, reducedDI, reducedRV, reduced3D) = chiSquaredCalc3D(realData,modelData,head['NUDI'],head['NURV'],head['NU'],False)
        for i in range(len(bestFit)):
            if i==2:
                bestStr+=paramStrs[2]+" = "+str(bestFit[2])
                if (bestFit[2]!=0):
                    bestStr+=", OR  "+str(1.0/(bestFit[2]/1000.0))+'[PC]\n'
                else:
                    bestStr+='\n'
            elif i==1:
                if bestFit[i]<0.1:
                    mJupMult=(const.KGperMsun/const.KGperMjupiter)
                    bestStr+=paramStrs[i]+" = "+str(bestFit[i])+", OR "+str(bestFit[i]*mJupMult)+' in [Mjupiter]\n'
                else:
                    bestStr+=paramStrs[i]+" = "+str(bestFit[i])+'\n'
            else:
                bestStr+=paramStrs[i]+" = "+str(bestFit[i])+'\n'
        bestStr+='\n'+'*'*100+'\nBEST REDUCED CHISQUAREDS: [total,DI,RV] = ['+str(reduced3D)+", "+str(reducedDI)+", "+str(reducedRV)+"]\n"+'*'*100
        f.write(bestStr)
    except:
        log.critical("A problem occured while trying to produce advanced summary strings.")
    f.write('\n'+clStr)
    f.write('\n'+burnInStr)
    f.write('\n'+grStr)
    f.write(effPtsStr)
    f.write('\n\n'+'-'*24+'\nDurations of each stage:\n'+'-'*24+'\n'+durationStrings)
    f.write('\nPost-Processing took: '+timeStrMaker(postTime)+'\n')
    f.write('Total simulation took: '+timeStrMaker(allTime)+'\n')
    f.write('\n\nEND OF RESULTS :-D')
    f.close()
    log.info("Summary file written to:\n"+summaryFname)  
    
def keplersThird(p=0,atot=0,mtot=0):
    """
    Kepler's Third rule.  
    Find the missing value provided you know 2 of them.
    """
    if (atot==0)and(mtot!=0)and(p!=0):
        atot = (((P**2)*(const.secPerYear**2)*const.Grav*const.KGperMsun*mtot)/(4.0*const.pi**2))**(1.0/3.0)
        atot = atot/const.MperAU
    elif (atot!=0)and(mtot==0)and(p!=0):
        mtot = ((((atot*const.MperAU)**3)*4.0*const.pi**2)/((p**2)*(const.secPerYear**2)*const.Grav*const.KGperMsun))
    elif (atot!=0)and(mtot!=0)and(p==0):
        p = np.sqrt((((atot*const.MperAU)**3)*4.0*(const.pi**2))/(const.Grav*mtot*(const.secPerYear**2)*const.KGperMsun))
    else:
        log.critical('More than 1 parameter was zero, so I can not calc K3')
    
    return (p,atot,mtot)
    
def recheckFit3D(orbParams,settingsDict,finalFits='',nus=[]):
    if finalFits!='':
        (head,data) = rwTools.loadFits(finalFits)
        nu = head['NU']
        nuDI = head['NUDI']
        nuRV = head['NURV']
    elif len(nus)>1:
        nu = nus[0]
        nuDI = nus[1]
        nuRV = nus[2]
    else:
        log.error("nus and finalFits not defined, so cannont calc reduced chiSquareds")
        nu = 1.0
        nuDI = 1.0
        nuRV = 1.0
        
    ##get the real data
    realData = rwTools.loadRealData(os.path.join(settingsDict['settingsDir'],settingsDict['prepend']),dataMode=settingsDict['dataMode'][0])
    ## Make Orbit cpp obj
    Orbit = cppTools.Orbit()
    try:
        pasa = settingsDict["pasa"][0]
    except:
        pasa = False
    Orbit.loadStaticVars(settingsDict['omegaFdi'][0],settingsDict['omegaFrv'][0],settingsDict['lowEcc'][0],pasa)
    Orbit.loadConstants(const.Grav,const.pi,const.KGperMsun, const.daysPerYear,const.secPerYear,const.MperAU)
    ## ensure orbParams are in required format for Orbit
    params = []
    for par in orbParams:
        params.append(par)
    params=np.array(params,dtype=np.dtype('d'),order='C')
    Orbit.loadRealData(realData)
    predictedData = np.ones((realData.shape[0],3),dtype=np.dtype('d'),order='C')
    Orbit.calculate(predictedData,params)
    ## Calculate chi squareds for 3D,DI,RV and update bestPars and bestSumStr if this is better than the best
    (raw3D, reducedDI, reducedRV, reduced3D) = chiSquaredCalc3D(realData,predictedData,nuDI,nuRV,nu)
    print '(raw3D, reducedDI, reducedRV, reduced3D) = '+repr((raw3D, reducedDI, reducedRV, reduced3D))
    
def predictLocation(orbParams,settingsDict,epochs=[]):
    
    ##get the real data
    realData = rwTools.loadRealData(os.path.join(settingsDict['settingsDir'],settingsDict['prepend']),dataMode=settingsDict['dataMode'][0])
    ## Make Orbit cpp obj
    Orbit = cppTools.Orbit()
    try:
        pasa = settingsDict["pasa"][0]
    except:
        pasa = False
    Orbit.loadStaticVars(settingsDict['omegaFdi'][0],settingsDict['omegaFrv'][0],settingsDict['lowEcc'][0],pasa)
    Orbit.loadConstants(const.Grav,const.pi,const.KGperMsun, const.daysPerYear,const.secPerYear,const.MperAU)
    ## ensure orbParams are in required format for Orbit
    params = []
    for par in orbParams:
        params.append(par)
    params=np.array(params,dtype=np.dtype('d'),order='C')
    fakeData = np.ones((len(epochs),7),dtype=np.dtype('d'),order='C')
    fakeData[:,0]=epochs[:]
    Orbit.loadRealData(fakeData)
    predictedData = np.ones((len(epochs),3),dtype=np.dtype('d'),order='C')
    print "fakeData are:\n"+repr(fakeData)
    print "predicted epochs data before are:\n"+repr(predictedData)
    Orbit.calculate(predictedData,params)
    print "predicted epochs data are:\n"+repr(predictedData)
    
def chiSquaredCalc3D(realData,modelData,nuDI,nuRV,nu3D,pasa=False): 
    """
    Based on definition, chiSquared=sum((modelVal_i-dataVal_i)^2/(dataError_i^2)) over all values of 'i'.
    This function will do so for DI, RV and 3D sets of data and provide the reduced chi squared for each.
    The raw 3D value will also be returned.
    NOTES: realData is the standard 7 parameter format from loadRealData function, and modelData is 
           the standard 3 param format.
           
    returned (raw3D, reducedDI, reducedRV, reduced3D)
    """   
    ## convert E,N to SA,PA?
    if pasa:
        (PA,PA_error,SA,SA_error) = ENtoPASA(modelData[:,0], 0, modelData[:,1], 0)
        modelData[:,0] = PA
        modelData[:,1] = SA
    diffs = np.concatenate(((realData[:,1]-modelData[:,0]),(realData[:,3]-modelData[:,1]),(realData[:,5]-modelData[:,2])))
    errors = np.concatenate((realData[:,2],realData[:,4],realData[:,6]))
    raw3D = np.sum((diffs**2)/(errors**2))
    diffsDI = np.concatenate(((realData[:,1]-modelData[:,0]),(realData[:,3]-modelData[:,1])))
    errorsDI = np.concatenate((realData[:,2],realData[:,4]))
    diffsRV = (realData[:,5]-modelData[:,2])
    errorsRV = realData[:,6][np.where(diffsRV!=0)]
    rawDI = np.sum((diffsDI[np.where(diffsDI!=0)]**2)/(errorsDI[np.where(diffsDI!=0)]**2))
    rawRV = np.sum((diffsRV[np.where(diffsRV!=0)]**2)/(errorsRV**2))
    return (raw3D,rawDI/nuDI,rawRV/nuRV,raw3D/nu3D)

def copyToDB(settingsDict):
    """
    Copy vital results files to Dropbox.
    """
    fnamesALL = []
    for extension in ['pdf','txt','log']:
        fnames = glob.glob(os.path.join(settingsDict['finalFolder'],"*."+extension))
        log.debug('found files to copy to DB:\n'+repr(fnames))
        for name in fnames:
            fnamesALL.append(name)
    dbDir = os.path.join(settingsDict['dbFolder'],settingsDict['outRoot'])
    os.mkdir(dbDir)
    log.debug('DB dir is:\n'+repr(dbDir))
    for f in fnamesALL:
        try:
            log.debug('trying to copy file:\n'+repr(os.path.basename(f)))
            shutil.copy(f,os.path.join(dbDir,os.path.basename(f)))
        except:
            log.error('failed to move file:\n'+f+'\nintto DB folder:\n'+dbDir)
    log.info("vital results files copied to DB folder:\n"+dbDir)
    
def getParInts(head):
    """
    convert string version of paramInts into a list again.
    """
    s = head['parInts'] 
    ints = s.split("[")[1].split("]")[0].split(',')
    parInts = []
    for i in ints:
        parInts.append(int(i))  
    return parInts        
        
def confLevelFinder(filename, colNum=False, returnData=False, returnChiSquareds=False, returnBestDataVal=False):
    """
    A function to find the 68.3 and 95.4% confidence levels in a given output data file's column 
    centered on the median value.
    
    return [[68.3% minimum, 68.3% maximum],[95.5% minimum, 95.5% maximum]]
    
    columnNum must be an int.    
    
    NOTE: This function calculates the confidence levels based on the median, while the
          more common standard is to centere it on the mean.  Doing that though requires a 
          time consuming loop over the data to total up that around the mean...
          Can't think of any faster way, so not doing it for now.
    """
    verboseInternal = False
    log.debug('Inside confLevelFinder')
    outStr=''
    if os.path.exists(filename):
        (dataAry,chiSquareds,[bestDataVal,dataMedian,dataValueStart,dataValueMid,dataValueEnd]) = rwTools.dataReader(filename, colNum)
    
        if len(dataAry>0) or (dataValueStart!=dataValueMid!=dataValueEnd):
            #Convert data array to a sorted numpy array
            dataAry = np.sort(dataAry)
            size = dataAry.size
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
                    log.error('confLevelFinder: ERROR!!! No FINE 68.3% confidence levels were found')
                    if (len(conf68ValsRough)==0):
                        log.error('confLevelFinder: ERROR!!! No ROUGH 68% confidence levels were found, so returning [0,0]')
                        conf68Vals = [0,0]
                    else:
                        conf68Vals = conf68ValsRough
                        log.error("confLevelFinder: Had to use ROUGH 68% [68,69] as no FINE 68.3% was found. So, using range "+repr(conf68Vals))                
                if (len(conf95Vals)==0):
                    log.error('confLevelFinder: ERROR!!! No FINE 95.4% confidence levels were found')
                    if (len(conf95ValsRough)==0):
                        log.error('confLevelFinder: ERROR!!! No ROUGH 95% confidence levels were found, so returning [0,0]')
                        conf95Vals = [0,0]
                    else:
                        conf95Vals = conf95ValsRough
                        log.error("confLevelFinder: Had to use ROUGH 95% [95,96] as no FINE 95.4% was found. So, using range "+repr(conf95Vals))
        else:
            ## There was no useful data, so return values indicating that
            dataAry=bestDataVal=dataMedian=dataValueStart
            conf68Vals = [dataValueStart,dataValueStart]
            conf95Vals = [dataValueStart,dataValueStart]
            chiSquareds = 0
            
        mJupMult=(const.KGperMsun/const.KGperMjupiter)
        s = "\nFinal Range values:\nTOTAL "+repr([dataAry[0],dataAry[-1]])
        s+= "\n68%   "+repr(conf68Vals)+'\n95%   '+repr(conf95Vals)+'\n'
        s+= "\nerror is centered on Median \n"
        s+="68.3% error level = "+str(dataMedian-conf68Vals[0])
        s+=" ->   "+str(dataMedian)+'  +/-  '+str(dataMedian-conf68Vals[0])+'\n'
        if (colNum==1) and (dataMedian<0.1):
            s=s+"Or in Mjup: ->   "+str(dataMedian*mJupMult)+'  +/-  '+str(dataMedian*mJupMult-conf68Vals[0]*mJupMult)+'\n'
        outStr+=s
        s=s+'\n'+75*'-'+'\n Leaving confLevelFinder \n'+75*'-'+'\n'
        log.debug(s)
        
        if verboseInternal:
            print 'returnData = '+repr(returnData)+', returnChiSquareds = '+repr(returnChiSquareds)+', returnBestDataVal = '+repr(returnBestDataVal)
        ## return requested form of results
        if (returnData and returnChiSquareds and (returnBestDataVal==False)):
            if verboseInternal:
                print 'returning first 3'
            returnList =  ([conf68Vals,conf95Vals],dataAry, chiSquareds)
        elif (returnData and returnChiSquareds and returnBestDataVal):
            if verboseInternal:
                print 'returning all 4'
            returnList =   ([conf68Vals,conf95Vals],dataAry, chiSquareds, bestDataVal,outStr) ##MODIFIED
        elif (returnData and (returnChiSquareds==False)and (returnBestDataVal==False)):
            if verboseInternal:
                print 'returning data only'
            returnList =   ([conf68Vals,conf95Vals],dataAry)
        elif (returnData and (returnChiSquareds==False) and returnBestDataVal):
            if verboseInternal:
                print 'returning data and bestval'
            returnList =   ([conf68Vals,conf95Vals],dataAry, bestDataVal,outStr) ##MODIFIED
        elif ((returnData==False) and returnChiSquareds):
            if verboseInternal:
                print 'returning just chiSquareds'
            returnList =   ([conf68Vals,conf95Vals], chiSquareds)
        elif ((returnData==False) and returnChiSquareds and returnBestDataVal):
            if verboseInternal:
                print 'returning chiSquareds and bestval'
            returnList =   ([conf68Vals,conf95Vals], chiSquareds, bestDataVal)
        elif ((returnData==False)and(returnChiSquareds==False) and returnBestDataVal):
            if verboseInternal:
                print 'returning CLevels and bestval'
            returnList = ([conf68Vals,conf95Vals], bestDataVal)
        elif ((returnData==False)and(returnChiSquareds==False) and (returnBestDataVal==False)):
            if verboseInternal:
                print 'returning only CLevels'
            returnList =   [conf68Vals,conf95Vals]
        
        return returnList 
        
    else:
        log.critical( "confLevelFinder: ERROR!!!! file doesn't exist")            
                                         
def findBestOrbit(filename,bestToFile=True,findAgain=False):        
    """
    Find the orbital elements for the best fit in a ExoSOFT format fits file.
    """             
    bestFname = os.path.join(os.path.dirname(filename),'bestOrbitParams.txt')  
    gotIt = False
    if os.path.exists(bestFname)and(findAgain==False):
        try:
            orbBest = np.loadtxt(bestFname,delimiter=',')
            log.debug("Using previously found best orbit in file:\n"+bestFname)   
            gotIt=True
        except:
            log.error("Tried to load previously found best orbit from file, but failed, so will find it from data again.")
    if gotIt==False:
        log.debug("trying to find best orbit in file:\n"+filename)   
        (head,data) = rwTools.loadFits(filename)
        chiBest = np.min(data[:,11])
        loc = np.where(data[:,11]==chiBest)
        orbBest = data[loc[0][0],:]
        log.info("Best fit found to be:\n"+repr(orbBest))
        if bestToFile:
            s=''
            for val in orbBest:
                s+=str(val)+","
            f = open(bestFname,'w')
            f.write(s[:-1])
            f.close()
            log.info("Best fit params written to :\n"+bestFname)
    return orbBest
                               
def copytree(src, dst):
    """
    Recursively copy a directory and its contents to another directory.
    
    WARNING: this is not advised for higher level folders as it can also copy subfolders 
    thus leading to a very large copy command if not careful.
    
    Code taken and simplified from:
    http://stackoverflow.com/questions/1868714/how-do-i-copy-an-entire-directory-of-files-into-an-existing-directory-using-pyth
    """
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            try:
                shutil.copytree(s, d)
                log.debug("Copying:\n "+repr(s)+'\nto:\n'+repr(d))
            except:
                log.error('FAILED while copying:\n'+repr(s)+'\nto:\n'+repr(d))
        else:
            try:
                shutil.copy2(s, d)
                log.debug("Copying:\n "+repr(s)+'\nto:\n'+repr(d))
            except:
                log.error('FAILED while copying:\n'+repr(s)+'\nto:\n'+repr(d))                        
    
def ENtoPASA(E, E_error, N, N_error):
    """
    Will calculate the Separation and Position Angles for a given East and North, including their errors.
    PA and error will be in [deg], with SA and error in ["]
    :returns: (PA,PA_error,SA,SA_error)
    """
    verbose = False
    PA = np.degrees(np.arctan2(E,N))
    #NOTE: both math.atan2 and np.arctan2 tried with same results, both produce negatives rather than continuous 0-360 
    #thus, must correct for negative outputs
    if type(PA)!=np.ndarray:
        if PA<0:
            PA = PA+360.0
    else:
        for i in range(0,len(PA)):
            if PA[i]<0:
                PA[i]=PA[i]+360.0
    
    SA = np.sqrt(E**2.0 + N**2.0)
    
    PA_error=0
    SA_error=0
    if (E_error==0)or(N_error==0):
        if False:
            print "either the E or N error value was zero, so setting the PA and SA return errors to zero!!"
    else:
        top = abs(E/N)*np.sqrt((E_error/E)**2.0 + (N_error/N)**2.0)
        btm = 1.0+(E/N)**2.0
        PA_error = abs(np.degrees(top/btm))

        top = SA*(abs(E*E_error)+abs(N*N_error))
        btm = E**2.0+N**2.0
        SA_error = abs(top/btm)
    if verbose:
        print repr((E, E_error, N, N_error))+" -> "+repr((PA,PA_error,SA,SA_error))    
    return (PA,PA_error,SA,SA_error)

def PASAtoEN(PA,PA_error,SA,SA_error):
    """
    Convert provided Position Angle and Separation Angle, and their errors, into 
    RA and DEC with errors.  These are the same equations for calculating 
    x and y in the Thiele-Innes orbit fitting.  Remember that x and y are 
    flipped in that fitting approach due to how Thiele defined the coord 
    system when deriving the equations used.
    
    NOTE: this can also be used to calculate x and y used in Thiele-Innes
          With East=RA=y and North=DEC=x.  
    
    :returns: (E, E_error, N, N_error)
    """
    verbose = False
    N = SA*np.cos(np.radians(PA))
    E = SA*np.sin(np.radians(PA))
    
    E_error=0
    N_error=0
    if (SA_error==0)or(PA_error==0):
        if verbose:
            print "either the PA and SA error value was zero, so setting the E or N return errors to zero!!"
    else:
        tempA = (SA_error/SA)**2.0
        tempB = ((np.cos(np.radians(PA+PA_error))-np.cos(np.radians(PA))) / np.cos(np.radians(PA)))**2.0
        N_error = abs(N*np.sqrt(tempA+tempB))
        
        # Another way to calculate the error, but the one above is belived to be more correct 
        tempA2 = (SA_error*np.cos(np.radians(PA)))**2.0
        tempB2 = (SA*np.sin(np.radians(PA))*np.radians(PA_error))**2.0
        N_error2 = np.sqrt(tempA2+tempB2)
        
        tempC = (SA_error/SA)**2.0
        tempD = ((np.sin(np.radians(PA+PA_error))-np.sin(np.radians(PA))) / np.sin(np.radians(PA)))**2.0
        E_error = abs(E*np.sqrt(tempC+tempD))
        
        # Another way to calculate the error, but the one above is belived to be more correct 
        tempC2 = (SA_error*np.sin(np.radians(PA)))**2.0
        tempD2 = (SA*np.cos(np.radians(PA))*np.radians(PA_error))**2.0
        E_error2 = np.sqrt(tempC2+tempD2)
        
        if verbose:
            print 'N_error2-N_error = '+str(N_error2-N_error)
            print 'E_error2-E_error = '+str(E_error2-E_error)
            print 'E_error2 = '+str(E_error2)+', E_error = '+str(E_error)
            print 'N_error2 = '+str(N_error2)+', N_error = '+str(N_error)+"\n"
    
    return (E, E_error, N, N_error)