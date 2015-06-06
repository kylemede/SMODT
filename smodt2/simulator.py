#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import numpy as np
import os
import copy
#from scipy.constants.codata import precision
import tools
import timeit
from tools import constants as const

class Simulator(object):
    """
    This is the Simulator parent class.  
    It contains the functions to perform basic 'shotgun' Monte Carlo, 
    Simulated Annealing, Sigma Tunning, and pure MCMC simulations.
    """
    def __init__(self,settingsDict):
        self.paramsLast = 0
        self.paramsBest = 0
        self.bestSumStr = ''
        self.latestSumStr = ''
        self.bestRedChiSqr = 1e6
        self.stgNsampDict = {'SA':'nSAsamp','ST':'nSTsamp','MC':'nSamples','MCMC':'nSamples'}
        self.acceptCount = 0
        self.acceptStr = ''
        self.shiftStr = ''
        self.acceptBoolAry = []
        self.parIntVaryAry = []
        self.chainNum =0
        self.settingsDict = settingsDict
        self.log = tools.getLogger('main.simulator',lvl=100,addFH=False)
        tools.logSystemInfo(self.log)
        self.realData = tools.loadRealData(os.path.join(self.dictVal('settingsDir'),self.dictVal('prepend')),dataMode=self.dictVal('dataMode'))
        (self.rangeMaxs,self.rangeMins,self.starterSigmas,self.paramInts,self.nu,self.nuDI,self.nuRV) = self.starter() 
        self.Orbit = tools.cppTools.Orbit()
        self.Orbit.loadStaticVars(self.dictVal('omegaFdi'),self.dictVal('omegaFrv'))
        self.Orbit.loadRealData(self.realData)
        self.Orbit.loadConstants(const.Grav,const.pi,const.KGperMsun, const.daysPerYear,const.secPerYear,const.MperAU)
        self.seed = int(timeit.default_timer())#reset in resetTracked()
        np.random.seed(self.seed)
        
        
    def starter(self):
        """
        Things needed by all simulator modes that can be internal, but 
        need code to load them up.
        """
        ##check there are matching number of RV datasets and provided min/max vals for offsets
        if np.max(self.realData[:,-1])!=(len(self.dictVal('vMINs'))-1):
            self.log.error("THE NUMBER OF vMINs DOES NOT MATCH THE NUMBER OF RV DATASETS!!!\n"+\
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
        ##start with uniform sigma values
        sigSize = self.dictVal('strtSig')
        sigmas = [sigSize,sigSize,sigSize,sigSize,sigSize,sigSize,sigSize,sigSize,sigSize,sigSize,0,0,sigSize]
        if len(self.dictVal('vMINs'))!=len(self.dictVal('vMAXs')):
            self.log.critical("THE NUMBER OF vMINs NOT EQUAL TO NUMBER OF vMAXs!!!")
        for i in range(0,len(self.dictVal('vMINs'))):
            sigmas.append(sigSize)
            rangeMins.append(self.dictVal('vMINs')[i])
            rangeMaxs.append(self.dictVal('vMAXs')[i])
        ##figure out which parameters are varying in this run.
        ##don't vary atot or chiSquared ever, 
        ##and take care of TcEqualT and Kdirect cases
        paramInts = []
        for i in range(0,len(rangeMins)):
            if (i!=10)and(i!=11):
                if (i>12):
                    if self.dictVal('dataMode')!='DI':
                        if rangeMaxs[i]!=0:
                            paramInts.append(i) 
                elif (i==8)or(i==12):
                    if (self.dictVal('dataMode')!='RV')or(self.dictVal('Kdirect')==False):
                        if (rangeMaxs[8]!=0)and(i==8):
                            paramInts.append(8)
                    elif self.dictVal('Kdirect'):
                        if (rangeMaxs[12]!=0)and(i==12):
                            paramInts.append(12)                                           
                elif (i==3)or(i==0)or(i==1):
                    if (self.dictVal('dataMode')!='RV'):
                        if(rangeMaxs[i]!=0):
                            paramInts.append(i)
                    elif self.dictVal('Kdirect')==False:
                        if(rangeMaxs[i]!=0):
                            paramInts.append(i)
                elif rangeMaxs[i]!=0:
                    if (i==5)or(i==6):
                        if self.dictVal('TcStep')and(i!=5):
                                paramInts.append(i)
                        elif (self.dictVal('TcStep')==False)and(i!=6):
                                paramInts.append(i)
                    else:
                        paramInts.append(i)
        paramInts = np.array(paramInts)
        ##find total number of RV and DI epochs in real data
        nDIepochs = np.sum(np.where(self.realData[:,2]<1e6,1,0))
        nRVepochs = np.sum(np.where(self.realData[:,6]<1e6,1,0))
        nEpochs = len(self.realData[:,0])
        ##Take mass1, dist, inc and period from those include in nu calcs
        ##as they have clear priors.
        paramIntsClean = copy.deepcopy(paramInts)
        notInNuInts = [0,2,7,8]      
        for val in notInNuInts:
            paramIntsClean=paramIntsClean[np.where(paramIntsClean!=val)]
        nDIvars = np.sum(np.where(paramIntsClean<10,1,0))
        self.log.debug('DIvars = '+repr(paramIntsClean[np.where(paramIntsClean<10)]))
        self.log.debug('RVvars = '+repr(paramIntsClean[np.where(paramIntsClean!=3)]))
        nRVvars = np.sum(np.where(paramIntsClean!=3,1,0))
        if nDIepochs==0:
            nVars = nRVvars
        elif nRVepochs==0:
            nVars = nDIvars
        else:
            nVars = len(paramInts)
        nu = nDIepochs*2+nRVepochs-nVars
        self.log.debug('[nEpochs, nDIepochs, nRVepochs] = ['+str(nEpochs)+', '+str(nDIepochs)+', '+str(nRVepochs)+']')
        self.log.debug('[nVars, nDIvars, nRVvars] = ['+str(nVars)+', '+str(nDIvars)+', '+str(nRVvars)+']')
        nuDI = 1
        nuRV=1
        if nDIepochs>0:
            nuDI = nDIepochs*2-nDIvars
        if nRVepochs>0:
            nuRV = nRVepochs-nRVvars
        #load these into settings dict
        self.settingsDict["nRVdsets"] = (len(self.dictVal('vMINs')),"Number of RV data sets")
        self.settingsDict['nDIepoch'] = (nDIepochs,"Number of DI epochs")
        self.settingsDict['nRVepoch'] = (nRVepochs,"Number of RV epochs")
        self.settingsDict['n3Depoch'] = (nEpochs,"Number of 3D epochs")
        self.settingsDict['nu'] = (nu,"Total nu")
        self.settingsDict['nuDI'] = (nuDI,"nu for DI")
        self.settingsDict['nuRV'] = (nuDI,"nu for RV")
        paramIntsStr = repr(paramInts).replace(' ','')
        self.settingsDict['parInts'] = (paramIntsStr,"Varried params")
        self.settingsDict['chainNum'] = (self.chainNum,"chain number")

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
        Increment all varying parameters if MC.
        Else, just increment one of them at random.
        """
        parsOut = copy.deepcopy(pars)
        varyInt=0
        sig = 0
        ## vary all the params if MC, 
        ##(or special cases at beginning of SA where this func is called pretending to be MC.)
        if  ('MCMC'not in stage) and (stage=='MC'):
            for i in range(0,len(pars)):
                if i in self.paramInts:
                    parsOut[i]=np.random.uniform(self.rangeMins[i],self.rangeMaxs[i])
        ## vary a random param if SA, ST or MCMC
        else:
            varyInt = self.paramInts[np.random.randint(0,len(self.paramInts))]
            self.parIntVaryAry.append(varyInt)
            sig = sigs[varyInt]*(self.rangeMaxs[varyInt]-self.rangeMins[varyInt])
            parsOut[varyInt]=np.random.uniform(pars[varyInt]-sig,pars[varyInt]+sig)        
        ## if TcEqualT, push the varied one into the other
        if self.dictVal('TcEqualT'):
            if self.dictVal('TcStep'):
                parsOut[5]=parsOut[6]
            else:
                parsOut[6]=parsOut[5]
        ## if Kdirect not set, then inclination varys.
        ## then K=0 going into Orbit so Orbit will calc it
        if 8 in self.paramInts:
            parsOut[12] = 0
        return parsOut
    
    def rangeCheck(self,pars,numAccepted=0,stage=''):
        """
        Check if values inside allowed range
        """
        inRange=True
        paramsOut = copy.deepcopy(pars)
        for i in range(0,len(pars)):
            if i in self.paramInts:
                if (i==6)and(self.dictVal('TcStep')):
                    if (self.rangeMins[i]>pars[i])or(pars[i]>self.rangeMaxs[i]):
                        inRange=False
                elif (i==5)and(self.dictVal('TcStep')==False):
                    if (self.rangeMins[i]>pars[i])or(pars[i]>self.rangeMaxs[i]):
                        inRange=False
                elif (i!=5) and (i!=6):
                    if (self.rangeMins[i]>pars[i])or(pars[i]>self.rangeMaxs[i]):
                        inRange=False
        if (numAccepted==0)and(stage=='SA'):
            ##Jump as starting position was in poor part of param space. for SA only.
            paramsOut = self.increment(pars,np.zeros(pars.shape),stage='MC')
            self.log.debug("Chain #"+str(self.chainNum)+" Nothing accepted yet, so jumping to new starting position.")
            inRange=True
        return (paramsOut,inRange)
    
    def accept(self,sample,pars,modelData,nSaved=0,temp=1.0,stage=''):
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
            self.bestSumStr = stage+" chain #"+str(self.chainNum)+' BEST reduced chiSquareds so far: [total,DI,RV] = ['+str(paramsOut[11]/self.nu)+", "+str(chiSquaredDI/self.nuDI)+", "+str(chiSquaredRV/self.nuRV)+"]"
            self.paramsBest = paramsOut
            if self.latestSumStr=='':
                self.latestSumStr=stage+" chain #"+str(self.chainNum)+' Nothing accepted yet below chi squared max = '+str(self.dictVal('chiMAX'))
        ## check if this step is accepted
        accept = False
        if stage=='MC':
            if (paramsOut[11]/self.nu)<self.dictVal('chiMAX'):
                accept=True
        else:
            ## For SA after first sample, MCMC, and ST
            likelihoodRatio = np.e**((self.paramsLast[11] - paramsOut[11])/(2.0*temp))
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
            self.acceptCount+=1
            self.acceptBoolAry.append(1)
            self.paramsLast=paramsOut
            self.latestSumStr = stage+" chain #"+str(self.chainNum)+' Latest accepted reduced chiSquareds: [total,DI,RV] = ['+str(paramsOut[11]/self.nu)+", "+str(chiSquaredDI/self.nuDI)+", "+str(chiSquaredRV/self.nuRV)+"]"
        else:
            self.acceptBoolAry.append(0)
        if sample%(self.dictVal(self.stgNsampDict[stage])//self.dictVal('nSumry'))==0:
            ##log a summary
            sumStr = "below\n"+stage+" chain #"+str(self.chainNum)+", # Accepted: "+str(self.acceptCount)+", # Saved: "+str(nSaved)+", Finished: "+str(sample)+"/"+str(self.dictVal(self.stgNsampDict[stage]))+", Current Temp = "+str(temp)+"\n"
            sumStr+=self.latestSumStr+'\n'+self.bestSumStr+'\n'
            self.log.debug(sumStr)
        if False:
            print 'priorsRatio = '+str(priorsRatio)
            print 'inc priorsRatio = '+str(self.dictVal('incPrior')(paramsOut[8])/self.dictVal('incPrior')(self.paramsLast[8]))
            print 'likelihoodRatio = '+str(likelihoodRatio)
        return (paramsOut,accept)
    
    def tempDrop(self,sample,temp,stage=''):
        """
        Determine if it is time to drop the temp, and drop if it is.
        Total temperature range is [strtTemp,0.1), so the minimum 
        temperature is actually <1.0 meaning the last few temperature drops 
        will really push the currently found minimum towards its peak.
        There will be a fixed number of temperature steps = 'nTmpStps'.
        """
        if stage=='SA':
            if sample%(self.dictVal('nSAsamp')//self.dictVal('nTmpStps'))==0:
                temp-=(self.dictVal('strtTemp')-0.1)*(1.0/self.dictVal('nTmpStps'))
        return temp
    
    def sigTune(self,sample,sigs=[],stage=''):
        """
        Check if it is time to calculate the acceptance rate.
        If stage is ST, then it will also tune the sigmas.
        If MCMC then it will just calculate the acceptance 
        rate.  In both cases, 'nSigStps' throughout the simulation 
        a summary message will also be written to the log.
        """
        sigmasOut = copy.deepcopy(sigs)
        if (stage=='ST')or(stage=='MCMC'):
            if (sample%(self.dictVal(self.stgNsampDict[stage])//self.dictVal('nSigStps'))==0)and(self.acceptCount>1):
                self.acceptStr = '\n'+stage+" chain #"+str(self.chainNum)+'\n'
                self.shiftStr = ''
                self.parIntVaryAry = np.array(self.parIntVaryAry)
                self.acceptBoolAry = np.array(self.acceptBoolAry)
                self.acceptStr+="Number of steps used to calculate acceptance rate = "+repr(len(self.acceptBoolAry))+'\n'
                for i in self.paramInts:
                    ##calculate acceptance rate for each param
                    nAcc = np.sum(np.where(self.parIntVaryAry==i,self.acceptBoolAry,0))
                    nTot = len(np.where(self.parIntVaryAry==i)[0])
                    self.acceptStr+= 'parameter # '+str(i)+' acceptance = '+str(float(nAcc)/float(nTot))+'\n'
                    if stage=='ST':
                        ##check each rate to choose up/down shift and do so and update shiftStr
                        self.shiftStr+= '\n'+stage+" chain #"+str(self.chainNum)+'\nparameter # '+str(i)+" shifting sigma "+str(sigs[i])+" -> "
                        if ((float(nAcc)/float(nTot))>0.35)and(sigs[i]<self.dictVal('sigMax')):
                            sigmasOut[i]+=self.dictVal('sigMin')
                        elif ((float(nAcc)/float(nTot))<0.25)and(sigs[i]>self.dictVal('sigMin')):
                            sigmasOut[i]-=self.dictVal('sigMin')
                        self.shiftStr+=str(sigmasOut[i])+"\n"
                self.acceptBoolAry = []
                self.parIntVaryAry = []
                self.log.debug(self.acceptStr+self.shiftStr)
        return sigmasOut
    
    def endSummary(self,totSaved,temp,sigmas,stage=''):
        """
        Make a final summary of important statistics for the chain.
        """
        sumStr = "END OF "+stage+" CHAIN #"+str(self.chainNum)+" SUMMARY below\nFinalTemp = "
        sumStr+= str(temp)+"\nTotal number of steps accepted = "+str(self.acceptCount)+"\n"
        sumStr+= "Average acceptance rate = "
        sumStr+=str(float(self.acceptCount)/float(self.dictVal(self.stgNsampDict[stage])))+"\n"
        sumStr+= "Total number of steps stored = "+str(totSaved)+"\n"
        sumStr+=self.latestSumStr+'\n'+self.bestSumStr+'\n'
        sumStr+="Last params = "+repr(self.paramsLast)+'\n'
        sumStr+="Best params = "+repr(self.paramsBest)+'\n'
        if (stage=="ST")or(stage=="MCMC"):
            if stage=='ST':
                sumStr+="Final Sigmas = "+repr(sigmas)+'\n'
            sumStr+=self.acceptStr+self.shiftStr
        sumStr+='\n'+'='*70+'\n'
        self.log.debug(sumStr)
    
    def resetTracked(self):
        """
        Reset the internal strings, arys and counters.
        """
        self.bestSumStr = ''
        self.latestSumStr = ''
        self.acceptStr = ''
        self.shiftStr = ''
        self.bestRedChiSqr = 1e6
        self.acceptCount = 0
        self.acceptBoolAry = []
        self.parIntVaryAry = []
        self.seed = int(timeit.default_timer()/(self.chainNum+1))
        self.log.debug("Chain# "+str(self.chainNum)+" has random number seed = "+str(self.seed))
        np.random.seed(self.seed)
    
    def startSummary(self,pars,sigs,stage):
        startStr="Chain #"+str(self.chainNum)+' VALS AT START OF '+stage+' SIM:\n'
        startStr+= 'params = '+repr(pars)+'\n'
        startStr+= 'rangeMins = '+repr(self.rangeMins)+'\n'
        startStr+= 'rangeMaxs = '+repr(self.rangeMaxs)+'\n'
        startStr+= 'sigmas = '+repr(sigs)+'\n'
        startStr+= 'paramInts = '+repr(self.paramInts)+'\n'
        startStr+= '[nu,nuDI,nuRV] = ['+str(self.nu)+', '+str(self.nuDI)+', '+str(self.nuRV)+']\n'
        self.log.debug(startStr)
        
    def simulatorFunc(self,stage='',chainNum=1,startParams=[],startSigmas=[]):
        """
        The core function to perform the requested stage of the simulation ('MC','SA','ST','MCMC').
        If stage is SA or ST: final (params,sigmas) are returned, else nothing.
        """
        tic=timeit.default_timer()
        self.log.debug("Trying "+str(self.dictVal(self.stgNsampDict[stage]))+" samples for chain #"+str(chainNum)+" in "+stage+" mode.")
        self.chainNum = chainNum
        self.resetTracked()
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
        latestPars = copy.deepcopy(proposedPars)
        self.paramsLast = proposedPars
        self.startSummary(proposedPars,sigmas,stage)
        ##loop through each sample 
        ##Follows these steps: in Range?,calc model,accept?,Store?,increment,lower temp?,tune sigmas? DONE ->write output data
        for sample in range(0,self.dictVal(self.stgNsampDict[stage])):
            (proposedPars,inRange)=self.rangeCheck(proposedPars,len(acceptedParams),stage)
            if inRange:
                self.Orbit.calculate(modelData,proposedPars)
                (params,accept) = self.accept(sample,proposedPars,modelData,len(acceptedParams),temp,stage)
                if accept:
                    latestPars = copy.deepcopy(params)
                    if ('MCMC' not in stage)and(stage=='MC'):
                        acceptedParams.append(params)
                    elif (self.acceptCount%self.dictVal('saveInt'))==0:
                        acceptedParams.append(params)                   
            else:
                self.acceptBoolAry.append(0)
            proposedPars = self.increment(latestPars,sigmas,stage)
            temp = self.tempDrop(sample,temp,stage)
            sigmas = self.sigTune(sample,sigmas,stage)
            if (self.dictVal('logLevel')<50)and(sample%(self.dictVal(self.stgNsampDict[stage])//20)==0):
                bar.render(sample * 100 // self.dictVal(self.stgNsampDict[stage]), stage+str(chainNum)+' complete so far.')
        if self.dictVal('logLevel')<50:
            bar.render(100,stage+str(chainNum)+' complete!\n')
        toc=timeit.default_timer()
        self.log.debug(stage+" it took: "+str(int(toc-tic))+' seconds')#$$$$$ need time format function still $$$$$$$$$$$$$$
        self.endSummary(len(acceptedParams),temp,sigmas,stage)
        outFname = tools.writeFits('outputData'+stage+str(chainNum)+'.fits',acceptedParams,self.settingsDict)
        if stage=='ST':
            return (latestPars,sigmas,self.bestRedChiSqr)
        elif stage=='SA':
            #start ST at the best location with tight sigmas, and it will tune to ideal sigmas
            return (self.paramsBest,np.ones(np.array(sigmas).shape)*0.01,self.bestRedChiSqr)
        elif(stage=='MC')or(stage=='MCMC'):
            return outFname
        
#END OF FILE      