import numpy as np
import os
import copy
from scipy.constants.codata import precision
#np.set_printoptions(precision=15)
import tools
import timeit
from settings_and_inputData import constants



class Simulator(object):
    """
    This is the Simulator parent class.  
    It contains the functions to perform basic 'shotgun' Monte Carlo, 
    Simulated Annealing, Sigma Tunning, and pure MCMC simulations.
    """
    def __init__(self,settingsDict):
        self.settingsDict = settingsDict
        self.log = tools.getLogger('main.simulator',lvl=100,addFH=False)
        tools.logSystemInfo(self.log)
        self.Orbit = tools.cppTools.Orbit()
        self.Orbit.loadomegaOffsets(self.dictVal('omegaFdi'),self.dictVal('omegaFrv'))
        self.realData = tools.loadRealData(os.path.join(self.dictVal('settingsDir'),self.dictVal('prepend')))
        self.Orbit.loadRealData(self.realData)
        self.Orbit.loadConstants(constants.Grav,constants.pi,constants.KGperMsun, constants.daysPerYear,constants.secPerYear,constants.MperAU)
        (self.rangeMaxs,self.rangeMins,self.starterSigmas,self.paramInts,self.nu,self.nuDI,self.nuRV) = self.starter()
        self.seed = int(timeit.default_timer())
        self.log.info("random number seed = "+str(self.seed))
        np.random.seed(self.seed)
        self.paramsLast = 0
        self.bestSumStr = ''
        self.latestSumStr = ''
        self.bestRedChiSqr = 1e6
        self.stgNsampDict = {'SA':'nSAsamp','ST':'nSTsamp','MC':'nSamples','MCMC':'nSamples'}
        self.acceptBoolAry = []
        self.parIntVaryAry = []
        self.sigPercentMAX = 0.5 #$$ Put in settings dict???
        self.sigPercentMIN = 0.01 #$$ Put in settings dict???
        #Load real data here? doesn't change so might as well
        #could also make the empty model data array here too
        ##Examples
        #(memTrackProc,memLogFilename) = self.starter()
        #self._memLogFilename = memLogFilename
        #self._memTrackProc = memTrackProc
           
    def starter(self):
        """
        Things needed by all simulator modes that can be internal, but 
        need code to load them up.
        """
        self.log.info("In Simulator.starter")
        ##check there are matching number of RV datasets and provided min/max vals for offsets
        if np.max(self.realData[:,-1])!=(len(self.dictVal('vMINs'))-1):
            self.log.error("THE NUMBER OF vMINs/vMAXs DOES NOT MATCH THE NUMBER OF RV DATASETS!!!\n"+\
                           "please check the vMINs/vMAXs arrays in the simple settings file\n"+\
                           "to make sure they have matching lengths to the number of RV datasets.")
        ##load up range min,max and sigma arrayS
        rangeMaxs = [self.dictVal('mass1MAX'),\
               self.dictVal('mass2MAX'),\
               self.dictVal('distMAX'),\
               self.dictVal('OmegaMAX'),\
               self.dictVal('eMAX'),\
               self.dictVal('TMAX'),\
               self.dictVal('TMAX'),\
               self.dictVal('PMAX'),\
               self.dictVal('incMAX'),\
               self.dictVal('omegaMAX'),\
               0,\
               0,\
               self.dictVal('KMAX')]
        rangeMins = [self.dictVal('mass1MIN'),\
               self.dictVal('mass2MIN'),\
               self.dictVal('distMIN'),\
               self.dictVal('OmegaMIN'),\
               self.dictVal('eMIN'),\
               self.dictVal('TMIN'),\
               self.dictVal('TMIN'),\
               self.dictVal('PMIN'),\
               self.dictVal('incMIN'),\
               self.dictVal('omegaMIN'),\
               0,\
               0,\
               self.dictVal('KMIN')]
        sigmas = [0.1*(self.dictVal('mass1MAX')-self.dictVal('mass1MIN')),\
               0.1*(self.dictVal('mass2MAX')-self.dictVal('mass2MIN')),\
               0.1*(self.dictVal('distMAX')-self.dictVal('distMIN')),\
               0.1*(self.dictVal('OmegaMAX')-self.dictVal('OmegaMIN')),\
               0.1*(self.dictVal('eMAX')-self.dictVal('eMIN')),\
               0.1*(self.dictVal('TMAX')-self.dictVal('TMIN')),\
               0.1*(self.dictVal('TMAX')-self.dictVal('TMIN')),\
               0.1*(self.dictVal('PMAX')-self.dictVal('PMIN')),\
               0.1*(self.dictVal('incMAX')-self.dictVal('incMIN')),\
               0.1*(self.dictVal('omegaMAX')-self.dictVal('omegaMIN')),\
               0,\
               0,\
               0.1*(self.dictVal('KMAX')-self.dictVal('KMIN'))]
        for i in range(0,len(self.dictVal('vMINs'))):
            sigmas.append(0.1*(self.dictVal('vMAXs')[i]-self.dictVal('vMINs')[i]))
            rangeMins.append(self.dictVal('vMINs')[i])
            rangeMaxs.append(self.dictVal('vMAXs')[i])
        
        #figure out which parameters are varying in this run.
        #don't vary atot or chiSquared ever, and take care of TcEqualT cases
        paramInts = []
        for i in range(0,len(rangeMins)):
            if (i!=10)and(i!=11):
                if rangeMaxs[i]!=0:
                    if self.dictVal('TcEqualT'):
                        if self.dictVal('TcStep'):
                            if i!=5:
                                paramInts.append(i)
                        else:
                            if i!=6:
                                paramInts.append(i)
                    else:
                        paramInts.append(i)
        #find total number of RV and DI epochs in real data
        nDIepochs = np.sum(np.where(self.realData[:,1]!=0))
        nRVepochs = np.sum(np.where(self.realData[:,5]!=0))
        nEpochs = len(self.realData[:,0])
        nDIvars = np.sum(np.where(paramInts<10))
        nRVvars = np.sum(np.where(paramInts!=3))
        nVars = len(paramInts)
        nu = nDIepochs*2+nRVepochs-nVars
        nuDI = nDIepochs*2-nDIvars
        nuRV = nRVepochs-nRVvars
        if True:
            print 'rangeMins = '+repr(rangeMins)
            print 'rangeMaxs = '+repr(rangeMaxs)
            print 'sigmas = '+repr(sigmas)
            print 'paramInts = '+repr(paramInts)
        return (rangeMaxs,rangeMins,sigmas,paramInts,nu,nuDI,nuRV)
    
    def dictVal(self,key):
        """
        Get the value for a key in the settings dictionary.
        This will handle the values that are tuples and not
        returning the value.
        """
        if type(self.settingsDict[key])==tuple:
            return self.settingsDict[key][0]
        else:
            return self.settingsDict[key]
    
    def increment(self,pars=[],sigs=[],stage=''):
        """
        Randomly increment one of the parameters
        """
        parsOut = copy.deepcopy(pars)
        varyInt=0
        ## vary all the params if mcONLY
        if stage=='MC':
            for i in range(0,len(pars)):
                if i in self.paramInts:
                    parsOut[i]=np.random.uniform(self.rangeMins[i],self.rangeMaxs[i])
        ## vary a random param if SA, ST or MCMC
        else:
            varyInt = self.paramInts[np.random.randint(0,len(self.paramInts))]
            self.parIntVaryAry.append(varyInt)
            #print '\nbefore params[varyInt] = '+repr(pars[varyInt])
            #print 'np.random.uniform(params[varyInt]-sigmas[varyInt],params[varyInt]+sigmas[varyInt])   ='+repr(np.random.uniform(params[varyInt]-sigmas[varyInt],params[varyInt]+sigmas[varyInt]))
            parsOut[varyInt]=np.random.uniform(pars[varyInt]-sigs[varyInt],pars[varyInt]+sigs[varyInt]) 
            #print 'after params[varyInt] = '+repr(parsOut[varyInt])           
        ## if TcEqualT, push the varied one into the other
        if self.dictVal('TcEqualT'):
            if self.dictVal('TcStep'):
                parsOut[5]=parsOut[6]
            else:
                parsOut[6]=parsOut[5]
        if False:
            #print 'params = '+repr(params)
            print 'varyInt = '+str(varyInt)
            print 'input params[varyInt] = '+repr(pars[varyInt])
            print 'sigmas = '+repr(sigs)
            print 'sigmas[varyInt] = '+repr(sigs[varyInt])
            print 'parsOut[varyInt] = '+repr(parsOut[varyInt])
            print 'parsOut-params = '+repr(parsOut-pars)
            print 'after2 params[varyInt] = '+repr(parsOut[varyInt])
        return parsOut
    
    def rangeCheck(self,pars,numAccepted=0,stage=''):
        """
        Check if values inside allowed range
        """
        inRange=True
        paramsOut = copy.deepcopy(pars)
        for i in range(0,len(pars)):
            if i in self.paramInts:
                if (self.rangeMins[i]>pars[i])or(pars[i]>self.rangeMaxs[i]):
                    inRange=False
                    #print 'out of range param was # '+str(i)+" val = "+repr(pars[i])
        if ((inRange==False)and(numAccepted==0))and(stage=='SA'):
            #jump as starting position was not in dead space.
            paramsOut = self.increment(pars,np.zeros(pars.shape),stage='MC')
            self.log.debug("Nothing accepted yet, so jumping to new starting position.")
            inRange=True
        return (paramsOut,inRange)
    
    def accept(self,sample,pars,modelData,numAccepted=0,temp=1.0,stage=''):
        """
        First this will calculate chi squared for model vs real data.
        
        For mcOnly it performs simple chisquared cut-off acceptance 
        based on 'chiMAX' value in settingsDict.
        Else, it will calculate the priors and accept based on 
        the Metropolis-Hastings algorithm. The temp factor will 
        be set to 1.0 for MCMC and Sigma Tuning, and should be provided 
        for Simulated Annealing.
        """
        paramsOut = copy.deepcopy(pars)
        diffs = np.concatenate(((self.realData[:,1]-modelData[:,0]),(self.realData[:,3]-modelData[:,1]),(self.realData[:,5]-modelData[:,2])))
        errors = np.concatenate((self.realData[:,2],self.realData[:,4],self.realData[:,6]))
        paramsOut[11] = np.sum((diffs**2)/(errors**2))
        diffsDI = np.concatenate(((self.realData[:,1]-modelData[:,0]),(self.realData[:,3]-modelData[:,1])))
        errorsDI = np.concatenate((self.realData[:,2],self.realData[:,4]))
        diffsRV = (self.realData[:,5]-modelData[:,2])
        errorsRV = self.realData[:,6][np.where(diffsRV!=0)]
        chiSquaredDI = np.sum((diffsDI[np.where(diffsDI!=0)]**2)/(errorsDI[np.where(diffsDI!=0)]**2))
        chiSquaredRV = np.sum((diffsRV[np.where(diffsRV!=0)]**2)/(errorsRV**2))
        if (paramsOut[11]/self.nu)<self.bestRedChiSqr:
            self.bestRedChiSqr=(paramsOut[11]/self.nu)
            self.bestSumStr = 'BEST reduced chiSquareds so far: [total,DI,RV] = ['+str(paramsOut[11]/self.nu)+", "+str(chiSquaredDI/self.nuDI)+", "+str(chiSquaredRV/self.nuRV)+"]"
            if self.latestSumStr=='':
                self.latestSumStr='Nothing accepted yet below chi squared max = '+str(self.dictVal('chiMAX'))
        if False:
            print ''
            print 'reals = \n'+repr(self.realData[:,[1,3,5]])
            print 'modelData = \n'+repr(modelData)
            print 'real-model = \n'+repr(self.realData[:,[1,3,5]]-modelData)
            print 'diffs = \n'+repr(diffs)
            print 'errors = \n'+repr(errors)
            print '(diffs**2)/(errors**2) = \n'+repr((diffs**2)/(errors**2))
            print 'chiSquared = '+str(paramsOut[11])
        accept = False
        if stage=='MC':
            if (paramsOut[11]/self.nu)<self.dictVal('chiMAX'):
                accept=True
        else:
            #handle case where doing SA and nothing accepted yet
            if (temp!=0)and(numAccepted==0):
                if (paramsOut[11]/self.nu)<self.dictVal('chiMAX'):
                    #Forces an acceptance of current parameters.  Right idea???$$$$$$$$$$$$$$$
                    accept=True
            ## For SA after first sample, MCMC, and ST
            else:
                likelihoodRatio = np.e**((self.paramsLast[11] - paramsOut[11])/ (2.0*temp))
                ###### put all prior funcs together in dict??
                priorsRatio = (self.dictVal('ePrior')(paramsOut[4],paramsOut[7])/self.dictVal('ePrior')(self.paramsLast[4],self.paramsLast[7]))
                priorsRatio*= (self.dictVal('pPrior')(paramsOut[7])/self.dictVal('pPrior')(self.paramsLast[7]))
                priorsRatio*= (self.dictVal('incPrior')(paramsOut[8])/self.dictVal('incPrior')(self.paramsLast[8]))
                priorsRatio*= (self.dictVal('mass1Prior')(paramsOut[0])/self.dictVal('mass1Prior')(self.paramsLast[0]))
                priorsRatio*= (self.dictVal('mass2Prior')(paramsOut[1])/self.dictVal('mass2Prior')(self.paramsLast[1]))
                priorsRatio*= (self.dictVal('distPrior')(paramsOut[2])/self.dictVal('distPrior')(self.paramsLast[2])) 
                if np.random.uniform(0.0, 1.0)<=(priorsRatio*likelihoodRatio):
                    accept = True
        if accept:
            self.acceptBoolAry.append(1)
            self.paramsLast=paramsOut
            self.latestSumStr = 'Latest accepted reduced chiSquareds: [total,DI,RV] = ['+str(paramsOut[11]/self.nu)+", "+str(chiSquaredDI/self.nuDI)+", "+str(chiSquaredRV/self.nuRV)+"]"
        else:
            self.acceptBoolAry.append(0)
            ## jumping done in another place too, kill this jump?? $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            if ((numAccepted==False)and(stage=='SA'))and(sample%(self.dictVal('nSAsamp')//self.dictVal('nTemps'))==0):
                #jump as starting position was not in dead space. $$$$$ this is second jump, needed???
                paramsOut = self.increment(pars,np.zeros(params.shape),stage='MC')
                self.log.info("Nothing accepted after "+str(sample)+" samples, so jumping to new starting position.")
        ## Write a short summary to log
        #print "sample%(self.dictVal(self.stgNsampDict[stage])//self.dictVal('nSumry') = "+repr(sample%(self.dictVal(self.stgNsampDict[stage])//self.dictVal('nSumry')))
        if sample%(self.dictVal(self.stgNsampDict[stage])//self.dictVal('nSumry'))==0:
            ##Log some summary for this step
            self.log.info("\nAccepted: "+str(numAccepted)+", Finished: "+str(sample)+"/"+str(self.dictVal(self.stgNsampDict[stage]))+"\n"+self.latestSumStr+'\n'+self.bestSumStr+'\n')
        return (paramsOut,accept)
    
    def tempDrop(self,sample,temp,stage=''):
        """
        Determine if it is time to drop the temp, and drop if it is.
        Total temperature range is [strtTemp,0.1), so the minimum 
        temperature is actually <1.0 meaning the last few temperature drops 
        will really push the currently found minimum towards its peak.
        There will be a fixed number of temperature steps = 'nTemps'.
        """
        if stage=='SA':
            if sample%(self.dictVal('nSAsamp')//self.dictVal('nTemps'))==0:
                temp-=(self.dictVal('strtTemp')+0.9)*(1.0/self.dictVal('nTemps'))
                ##Log some summary for this step??
        return temp
    
    def sigTune(self,sample,sigs=[],stage=''):
        """
        Check if it is time to calculate the acceptance rate.
        If stage is 'ST'== Sigma Tuning, then it will also tune the 
        sigmas.  If 'MCMC' then it will just calculate the acceptance 
        rate.  In both cases, if it the sample (i) is at 10% of the total number 
        of samples, a summary message will also be written to the log.
        """
        sigmasOut = copy.deepcopy(sigs)
        if (stage=='ST')or(stage=='MCMC'):
            acceptStr = ''
            shiftStr = ''
            if sample%(self.dictVal(self.stgNsampDict[stage])//10)==0:
                for i in self.paramInts:
                    ##calculate acceptance rate for each param
                    nAcc = np.sum(self.acceptBoolAry[np.where(self.parIntVaryAry==i)])
                    nTot = len(np.where(self.parIntVaryAry==i)[0])
                    acceptStr+= 'parameter # '+str(i)+' acceptance = '+str(nAcc/nTot)+'\n'
                    if stage=='ST':
                        ##check each rate to choose up/down shift and do so and update shiftStr
                        shiftStr+= 'parameter # '+str(i)+" shifting sigma "+str(sigs[i])
                        if ((nAcc/nTot)>0.35)and(sigs[i]<self.sigPercentMAX):
                            sigmasOut[i]+=self.sigPercentMIN
                        elif ((nAcc/nTot)<0.25)and(sigs[i]>self.sigPercentMIN):
                            sigmasOut[i]-=self.sigPercentMIN
                        shiftStr+=str(sigmasOut[i])+"\n"
                self.acceptBoolAry = []
                self.parIntVaryAry = []
                ##Log some summary for this step
                self.log.info(acceptStr+shiftStr)
        return sigmasOut
    
    def simulatorFunc(self,stage='',startParams=[],startSigmas=[]):
        """
        The core function to perform the requested stage of the simulation ('MC','SA','ST','MCMC').
        If stage is SA or ST: final (params,sigmas) are returned, else nothing.
        """
        tic=timeit.default_timer()
        self.log.info("In Simulator.simulatorFunc")
        self.log.info("Trying "+str(self.dictVal(self.stgNsampDict[stage]))+" samples")
        print "Trying "+str(self.dictVal(self.stgNsampDict[stage]))+" samples"
        bar = tools.ProgressBar('green',width=30,block='=',empty='-',lastblock='>')
        modelData = np.zeros((len(self.realData),3))
        acceptedParams = []
        temp = 1.0
        if (stage=='SA')or(stage=='MC'):
            ## get starting params and sigmas as these two stages start at a random point
            sigmas = copy.deepcopy(self.starterSigmas)
            proposedPars = self.increment(np.zeros((len(self.rangeMins))),sigmas,stage='MC')
            
            if stage=='SA':
                temp=self.dictVal('strtTemp')
        else:
            proposedPars = copy.deepcopy(startParams)
            sigmas = copy.deepcopy(startSigmas)
            
        if False:
            print 'VALS AT START OF SIM:'
            print 'params = '+repr(proposedPars)
            print 'rangeMins = '+repr(self.rangeMins)
            print 'rangeMaxs = '+repr(self.rangeMaxs)
            print 'sigmas = '+repr(sigmas)
            print 'paramInts = '+repr(self.paramInts)
            print '-'*25
        for sample in range(0,self.dictVal(self.stgNsampDict[stage])):
            
            (proposedPars,inRange)=self.rangeCheck(proposedPars,len(acceptedParams),stage)
            #print 'sample = '+str(sample)+", temp = "+str(temp)+", number accepted = "+str(len(acceptedParams))+', inRange = '+repr(inRange)
            if inRange:
                self.Orbit.calculate(modelData,proposedPars)
                (params,accept) = self.accept(sample,proposedPars,modelData,len(acceptedParams),temp,stage)
                if accept:
                    acceptedParams.append(params)
                    #print 'acceptedParams[-1] = '+repr(acceptedParams[-1])
                    #print 'params just accepted are: '+repr(params)
            if len(acceptedParams)>0:
                latestPars=acceptedParams[-1]
                #print 'latest accepted = '+repr(latestPars)
                #print 'acceptedParams[-1] = '+repr(acceptedParams[-1])
                #print 'reset params to latest accepted'
            else:
                latestPars = proposedPars
                #print 'reset params to latest proposed as NONE accepted yet'
            proposedPars = self.increment(latestPars,sigmas,stage)
            #print 'latest proposed = '+repr(proposedPars)
            #print 'latestPars-proposed = '+repr(latestPars-proposedPars)
            temp = self.tempDrop(sample,temp,stage)
            sigmas = self.sigTune(sample,sigmas,stage)
            if (True)and(sample%(self.dictVal(self.stgNsampDict[stage])//20)==0):#self.dictVal('SILENT')
                bar.render(sample * 100 // self.dictVal(self.stgNsampDict[stage]), 'Complete so far.')
        if True:#self.dictVal('SILENT')
            bar.render(100, 'Complete so far.')
        toc=timeit.default_timer()
        print "\n"+stage+" it took: "+str(toc-tic)#$$$$$ need time format function still $$$$$$$$$$$$$$$$$$$$$$$$$$$
        #print '\nmodelData = \n'+repr(acceptedParams)
        ## SUMMARY of RUN log msg!!!!!!!!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        print '# Accepted = '+str(len(acceptedParams))
        tools.writeFits('outputData'+stage+'.fits',acceptedParams,self.settingsDict)
        if (stage=='SA')or(stage=='ST'):
            return (params,sigmas)
        
        

        