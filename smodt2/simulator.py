#@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
import numpy as np
import os
import gc
import copy
import sys
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
        self.nSaved = 0
        self.nSavedPeriodic = 0
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
        (self.rangeMaxsRaw,self.rangeMinsRaw,self.rangeMaxs,self.rangeMins,self.starterSigmas,self.paramInts,self.nu,self.nuDI,self.nuRV) = self.starter() 
        self.Orbit = tools.cppTools.Orbit()
        self.Orbit.loadStaticVars(self.dictVal('omegaFdi'),self.dictVal('omegaFrv'),self.dictVal('lowEcc'))
        self.Orbit.loadRealData(self.realData)
        self.Orbit.loadConstants(const.Grav,const.pi,const.KGperMsun, const.daysPerYear,const.secPerYear,const.MperAU)
        #Just initial seed val, reset in resetTracked() to be unique for each chain.
        self.seed = int(timeit.default_timer())
        np.random.seed(self.seed)
        self.tmpDataFile = os.path.join(self.dictVal('tmpDir'),self.dictVal('outRoot')+"-"+str(self.chainNum)+".npy")
        
    def starter(self):
        """
        Things needed by all simulator modes that can be internal, but 
        need code to load them up.
        """
        ##check there are matching number of RV datasets and provided min/max vals for offsets
        if np.min(self.realData[:,6])<1e6:
            numVmins=len(self.dictVal('vMINs'))
            if numVmins==0:
                numVmins=1
            if np.max(self.realData[:,7])!=(numVmins-1):
                self.log.error("THE NUMBER OF vMINs DOES NOT MATCH THE NUMBER OF RV DATASETS!!!\n"+\
                               "please check the vMINs/vMAXs arrays in the simple settings file\n"+\
                               "to make sure they have matching lengths to the number of RV datasets.")
        
        if (self.dictVal('TMAX')==-1)and(self.dictVal('TMIN')==-1):
            ## set T range to [earliest Epoch-max period,earliest epoch]
            self.settingsDict['TMAX']=np.min(self.realData[:,0])
            self.settingsDict['TMIN']=np.min(self.realData[:,0])-self.dictVal('PMAX')*const.daysPerYear
        ##load up range min,max and sigma arrayS
        rangeMaxs = [self.dictVal('mass1MAX'),\
               self.dictVal('mass2MAX'),\
               self.dictVal('paraMAX'),\
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
               self.dictVal('paraMIN'),\
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
        rangeMaxs = np.array(rangeMaxs)
        rangeMins = np.array(rangeMins)
        ##For lowEcc case, make Raw min/max vals for param drawing during MC mode
        rangeMaxsRaw = copy.deepcopy(rangeMaxs)
        rangeMinsRaw = copy.deepcopy(rangeMins)
        if self.dictVal('lowEcc'):
            ## run through the possible numbers for e and omega to find min/max for RAW versions
            fourMin=1e6
            fourMax=-1e6
            nineMin=1e6
            nineMax=-1e6
            for omeg in range(int(rangeMins[9]*10),int(rangeMaxs[9]*10),1):
                omega = float(omeg)/10.0
                for e in range(int(rangeMins[4]*100),int(rangeMaxs[4]*100),1):
                    ecc = float(e)/100.0
                    four = np.sqrt(ecc)*np.sin((np.pi/180.0)*omega)
                    nine = np.sqrt(ecc)*np.cos((np.pi/180.0)*omega)
                    if four>fourMax:
                        fourMax = four
                        #a= 'fourMax: e = '+str(ecc)+', omega = '+str(omega)
                    if four<fourMin:
                        fourMin = four
                        #b= 'fourMin: e = '+str(ecc)+', omega = '+str(omega)
                    if nine>nineMax:
                        nineMax = nine
                        #c= 'nineMax: e = '+str(ecc)+', omega = '+str(omega)
                    if nine<nineMin:
                        nineMin = nine
                        #d= 'nineMin: e = '+str(ecc)+', omega = '+str(omega)
            rangeMaxsRaw[9] = nineMax
            rangeMaxsRaw[4] = fourMax
            rangeMinsRaw[9] = nineMin
            rangeMinsRaw[4] = fourMin
            #print 'rangeMins[4] = '+str(rangeMins[4])+", rangeMaxs[4] = "+str(rangeMaxs[4])
            #print 'rangeMins[9] = '+str(rangeMins[9])+", rangeMaxs[9] = "+str(rangeMaxs[9])
            #print a
            #print b
            #print c
            #print d
            #print 'rangeMinsRaw[4] = '+str(rangeMinsRaw[4])+", rangeMaxsRaw[4] = "+str(rangeMaxsRaw[4])
            #print 'rangeMinsRaw[9] = '+str(rangeMinsRaw[9])+", rangeMaxsRaw[9] = "+str(rangeMaxsRaw[9])
            #sys.exit(0)
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
                elif (i==2)or(i==3)or(i==0)or(i==1):
                    if (self.dictVal('dataMode')!='RV'):
                        if(rangeMaxs[i]!=0):
                            paramInts.append(i)
                    elif (self.dictVal('Kdirect')==False)and((i!=3)and(i!=2)):
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
        notInNuInts = [2,7,8]      
        for val in notInNuInts:
            paramIntsClean=paramIntsClean[np.where(paramIntsClean!=val)]
        nDIvars = np.sum(np.where(paramIntsClean<10,1,0))
        self.log.debug('DIvars cleaned = '+repr(paramIntsClean[np.where(paramIntsClean<10)]))
        self.log.debug('RVvars cleaned = '+repr(paramIntsClean[np.where(paramIntsClean!=3)]))
        nRVvars = np.sum(np.where(paramIntsClean!=3,1,0))
        if nDIepochs==0:
            nVars = nRVvars
        elif nRVepochs==0:
            nVars = nDIvars
        else:
            nVars = len(paramInts)
        nu = nDIepochs*2+nRVepochs-nVars
        self.log.debug("vars = "+repr(paramInts))
        self.log.debug('[nEpochs, nDIepochs, nRVepochs] = ['+str(nEpochs)+', '+str(nDIepochs)+', '+str(nRVepochs)+']')
        self.log.debug('[nVars, nDIvars, nRVvars] = ['+str(nVars)+', '+str(nDIvars)+', '+str(nRVvars)+']')
        nuDI = 1
        nuRV=1
        if nDIepochs>0:
            nuDI = nDIepochs*2-nDIvars
        if nRVepochs>0:
            nuRV = nRVepochs-nRVvars
        self.log.debug('[nu, nuDI, nuRV] = ['+str(nu)+', '+str(nuDI)+', '+str(nuRV)+']')
        #load these into settings dict
        self.settingsDict["nRVdsets"] = (len(self.dictVal('vMINs')),"Number of RV data sets")
        self.settingsDict['nDIepoch'] = (nDIepochs,"Number of DI epochs")
        self.settingsDict['nRVepoch'] = (nRVepochs,"Number of RV epochs")
        self.settingsDict['n3Depoch'] = (nEpochs,"Number of 3D epochs")
        self.settingsDict['nu'] = (nu,"Total nu")
        self.settingsDict['nuDI'] = (nuDI,"nu for DI")
        self.settingsDict['nuRV'] = (nuRV,"nu for RV")
        paramIntsStr = repr(paramInts).replace(' ','')
        self.settingsDict['parInts'] = (paramIntsStr,"Varried params")
        self.settingsDict['chainNum'] = (self.chainNum,"chain number")
        return (rangeMaxsRaw,rangeMinsRaw,rangeMaxs,rangeMins,sigmas,paramInts,nu,nuDI,nuRV)
    
    def combinedPriors(self,paramsCurr):
        """
        """
        
    #print "likelihoodRatio = "+repr(likelihoodRatio)
    ###### put all prior funcs together in dict??
    #print 'paramsOut[4] = '+str(paramsOut[4])
    #print 'paramsOut[7] = '+str(paramsOut[7])
    #print 'paramsOut[4] = '+str(self.paramsLast[4])
    #print 'paramsOut[7] = '+str(self.paramsLast[7])
    priorsRatio = (self.dictVal('ePrior')(paramsOut[4],paramsOut[7])/self.dictVal('ePrior')(self.paramsLast[4],self.paramsLast[7]))
    #print "a = "+repr((self.dictVal('ePrior')(paramsOut[4],paramsOut[7])/self.dictVal('ePrior')(self.paramsLast[4],self.paramsLast[7])))
    priorsRatio*= (self.dictVal('pPrior')(paramsOut[7])/self.dictVal('pPrior')(self.paramsLast[7]))
    #print "b = "+repr((self.dictVal('pPrior')(paramsOut[7])/self.dictVal('pPrior')(self.paramsLast[7])))
    priorsRatio*= (self.dictVal('incPrior')(paramsOut[8])/self.dictVal('incPrior')(self.paramsLast[8]))
    #print "c = "+repr((self.dictVal('incPrior')(paramsOut[8])/self.dictVal('incPrior')(self.paramsLast[8])))
    priorsRatio*= (self.dictVal('mass1Prior')(paramsOut[0])/self.dictVal('mass1Prior')(self.paramsLast[0]))
    #print "d = "+repr((self.dictVal('mass1Prior')(paramsOut[0])/self.dictVal('mass1Prior')(self.paramsLast[0])))
    priorsRatio*= (self.dictVal('mass2Prior')(paramsOut[1])/self.dictVal('mass2Prior')(self.paramsLast[1]))
    #print "e = "+repr((self.dictVal('mass2Prior')(paramsOut[1])/self.dictVal('mass2Prior')(self.paramsLast[1])))
    priorsRatio*= (self.dictVal('paraPrior')(paramsOut[2])/self.dictVal('paraPrior')(self.paramsLast[2])) 
    #print "f = "+repr((self.dictVal('paraPrior')(paramsOut[2])/self.dictVal('paraPrior')(self.paramsLast[2])) )
    #print 'priorsRatio = '+str(priorsRatio)
    
    
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
                    parsOut[i]=np.random.uniform(self.rangeMinsRaw[i],self.rangeMaxsRaw[i])
        ## vary a random param if SA, ST or MCMC
        else:
            varyInt = self.paramInts[np.random.randint(0,len(self.paramInts))]
            self.parIntVaryAry.append(varyInt)
            sig = sigs[varyInt]*(self.rangeMaxsRaw[varyInt]-self.rangeMinsRaw[varyInt])
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
    
    def rangeCheck(self,pars,stage=''):
        """
        Check if values inside allowed range
        """
        inRange=True
        paramsOut = copy.deepcopy(pars)
        ## convert from Raw form if in lowEcc mode
        self.Orbit.convertParsFromRaw(paramsOut)
        for i in range(0,len(pars)):
            if i in self.paramInts:
                if (i==6)and(self.dictVal('TcStep')):
                    if (self.rangeMins[i]>paramsOut[i])or(paramsOut[i]>self.rangeMaxs[i]):
                        inRange=False
                elif (i==5)and(self.dictVal('TcStep')==False):
                    if (self.rangeMins[i]>paramsOut[i])or(paramsOut[i]>self.rangeMaxs[i]):
                        inRange=False
                elif (i!=5) and (i!=6):
                    if (self.rangeMins[i]>paramsOut[i])or(paramsOut[i]>self.rangeMaxs[i]):
                        inRange=False
        if (self.acceptCount==0)and(stage=='SA'):
            ##Jump as starting position was in poor part of param space. for SA only.
            paramsOut = self.increment(pars,np.zeros(pars.shape),stage='MC')
            ## convert from Raw form if in lowEcc mode
            self.Orbit.convertParsFromRaw(paramsOut)
            #self.log.debug("Chain #"+str(self.chainNum)+" Nothing accepted yet, so jumping to new starting position.")
            inRange=True
        return (paramsOut,inRange)
    
    def accept(self,sample,pars,modelData,temp=1.0,stage=''):
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
        ## Calculate chi squareds for 3D,DI,RV and update bestPars and bestSumStr if this is better than the best
        diffs = np.concatenate(((self.realData[:,1]-modelData[:,0]),(self.realData[:,3]-modelData[:,1]),(self.realData[:,5]-modelData[:,2])))
        errors = np.concatenate((self.realData[:,2],self.realData[:,4],self.realData[:,6]))
        paramsOut[11] = np.sum((diffs**2)/(errors**2))
        diffsDI = np.concatenate(((self.realData[:,1]-modelData[:,0]),(self.realData[:,3]-modelData[:,1])))
        errorsDI = np.concatenate((self.realData[:,2],self.realData[:,4]))
        diffsRV = (self.realData[:,5]-modelData[:,2])
        errorsRV = self.realData[:,6][np.where(diffsRV!=0)]
        chiSquaredDI = np.sum((diffsDI[np.where(diffsDI!=0)]**2)/(errorsDI[np.where(diffsDI!=0)]**2))
        chiSquaredRV = np.sum((diffsRV[np.where(diffsRV!=0)]**2)/(errorsRV**2))
        if self.bestSumStr=='':
            self.bestSumStr = stage+" chain #"+str(self.chainNum)+' Nothing accepted yet below chi squared max = '+str(self.dictVal('chiMAX'))
            self.latestSumStr="Latest reduced chiSquared : [total,DI,RV] = ["+str(paramsOut[11]/self.nu)+", "+str(chiSquaredDI/self.nuDI)+", "+str(chiSquaredRV/self.nuRV)+"]"
        if (paramsOut[11]/self.nu)<self.bestRedChiSqr:
            self.bestRedChiSqr=(paramsOut[11]/self.nu)
            self.bestSumStr = stage+" chain #"+str(self.chainNum)+' BEST reduced chiSquareds so far: [total,DI,RV] = ['+str(paramsOut[11]/self.nu)+", "+str(chiSquaredDI/self.nuDI)+", "+str(chiSquaredRV/self.nuRV)+"]"
            #self.bestSumStr+='\nbest total non-reduced chiSquared/nu = '+str(paramsOut[11])+'/'+str(self.nu)
            bestPars = copy.deepcopy(paramsOut)
            ## convert back to sqrt(e)sin(omega), sqrt(e)cos(omega) if in lowEcc mode
            self.Orbit.convertParsToRaw(bestPars)
            self.paramsBest = bestPars
        ## check if this step is accepted
        accept = False
        if stage=='MC':
            if (paramsOut[11]/self.nu)<self.dictVal('chiMAX'):
                accept=True
        else:
            ## For SA after first sample, MCMC, and ST
            try:
                #print 'self.paramsLast[11] = '+str(self.paramsLast[11])
                #print 'paramsOut[11] = '+str(paramsOut[11])
                likelihoodRatio = np.e**((self.paramsLast[11] - paramsOut[11])/(2.0*temp))
                #print "likelihoodRatio = "+repr(likelihoodRatio)
                ###### put all prior funcs together in dict??
                #print 'paramsOut[4] = '+str(paramsOut[4])
                #print 'paramsOut[7] = '+str(paramsOut[7])
                #print 'paramsOut[4] = '+str(self.paramsLast[4])
                #print 'paramsOut[7] = '+str(self.paramsLast[7])
                priorsRatio = (self.dictVal('ePrior')(paramsOut[4],paramsOut[7])/self.dictVal('ePrior')(self.paramsLast[4],self.paramsLast[7]))
                #print "a = "+repr((self.dictVal('ePrior')(paramsOut[4],paramsOut[7])/self.dictVal('ePrior')(self.paramsLast[4],self.paramsLast[7])))
                priorsRatio*= (self.dictVal('pPrior')(paramsOut[7])/self.dictVal('pPrior')(self.paramsLast[7]))
                #print "b = "+repr((self.dictVal('pPrior')(paramsOut[7])/self.dictVal('pPrior')(self.paramsLast[7])))
                priorsRatio*= (self.dictVal('incPrior')(paramsOut[8])/self.dictVal('incPrior')(self.paramsLast[8]))
                #print "c = "+repr((self.dictVal('incPrior')(paramsOut[8])/self.dictVal('incPrior')(self.paramsLast[8])))
                priorsRatio*= (self.dictVal('mass1Prior')(paramsOut[0])/self.dictVal('mass1Prior')(self.paramsLast[0]))
                #print "d = "+repr((self.dictVal('mass1Prior')(paramsOut[0])/self.dictVal('mass1Prior')(self.paramsLast[0])))
                priorsRatio*= (self.dictVal('mass2Prior')(paramsOut[1])/self.dictVal('mass2Prior')(self.paramsLast[1]))
                #print "e = "+repr((self.dictVal('mass2Prior')(paramsOut[1])/self.dictVal('mass2Prior')(self.paramsLast[1])))
                priorsRatio*= (self.dictVal('paraPrior')(paramsOut[2])/self.dictVal('paraPrior')(self.paramsLast[2])) 
                #print "f = "+repr((self.dictVal('paraPrior')(paramsOut[2])/self.dictVal('paraPrior')(self.paramsLast[2])) )
                #print 'priorsRatio = '+str(priorsRatio)
                if np.random.uniform(0.0, 1.0)<=(priorsRatio*likelihoodRatio):
                    accept = True
            except:
                accept = False
        if accept:
            self.acceptCount+=1
            self.acceptBoolAry.append(1)
            self.paramsLast=paramsOut
            self.latestSumStr = stage+" chain #"+str(self.chainNum)+' Latest accepted reduced chiSquareds: [total,DI,RV] = ['+\
            str(paramsOut[11]/self.nu)+", "+str(chiSquaredDI/self.nuDI)+", "+str(chiSquaredRV/self.nuRV)+"]"
        else:
            self.acceptBoolAry.append(0)
        ##log a status summary?
        if self.dictVal('nSumry')>0:
            if sample%(self.dictVal(self.stgNsampDict[stage])//self.dictVal('nSumry'))==0:
                sumStr = "below\n"+stage+" chain #"+str(self.chainNum)+", # Accepted: "+str(self.acceptCount)+", # Saved: "+\
                str(self.nSaved)+' (curr '+str(self.nSavedPeriodic)+"), Finished: "+str(sample)+"/"+\
                str(self.dictVal(self.stgNsampDict[stage]))+", Current Temp = "+str(temp)+"\n"
                sumStr+=self.latestSumStr+'\n'+self.bestSumStr+'\n'
                self.log.debug(sumStr)
                #print 'priorsRatio = '+repr(priorsRatio)
                #print 'likelihood ratio = '+rerp(likelihoodRatio)
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
                temp-=(self.dictVal('strtTemp')-0.01)*(1.0/self.dictVal('nTmpStps'))
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
                ##log a status summary?
                if self.dictVal('nSumry')>0:
                    if sample%(self.dictVal(self.stgNsampDict[stage])//self.dictVal('nSumry'))==0:
                        self.log.debug(self.acceptStr+self.shiftStr)
        else:
            #No need to track these, so reset them to empty to save any RAM they are using
            self.acceptBoolAry = []
            self.parIntVaryAry = []
        return sigmasOut
    
    def endSummary(self,temp,sigmas,stage=''):
        """
        Make a final summary of important statistics for the chain.
        """
        sumStr = '\n'+"="*70+"\nEND OF "+stage+" CHAIN #"+str(self.chainNum)+" SUMMARY:\nFinalTemp = "
        sumStr+= str(temp)+"\nTotal number of steps accepted = "+str(self.acceptCount)+"\n"
        sumStr+= "Average acceptance rate = "
        sumStr+=str(float(self.acceptCount)/float(self.dictVal(self.stgNsampDict[stage])))+"\n"
        sumStr+= "Total number of steps stored = "+str(self.nSaved)+"\n"
        sumStr+=self.latestSumStr+'\n'+self.bestSumStr+'\n'
        sumStr+="Last params  = "+repr(self.paramsLast)+'\n'
        sumStr+="Best params (Raw) = "+repr(self.paramsBest)+'\n'
        if (stage=="ST")or(stage=="MCMC"):
            if stage=='ST':
                sumStr+="Final Sigmas = "+repr(sigmas)+'\n'
            sumStr+=self.acceptStr+self.shiftStr
        sumStr+='\n'+'='*70+'\n'
        self.log.info(sumStr)
    
    def resetTracked(self,stage):
        """
        Reset the internal strings, arys and counters.
        """
        self.bestSumStr = ''
        self.latestSumStr = ''
        self.acceptStr = ''
        self.shiftStr = ''
        self.bestRedChiSqr = 1e6
        self.nSaved = 0
        self.acceptCount = 0
        self.nSavedPeriodic = 0
        self.acceptBoolAry = []
        self.parIntVaryAry = []
        if (stage=="MC")or(stage=="SA"):
            t = np.random.uniform(1,1e6)
            self.seed = int((timeit.default_timer()/(self.chainNum+1))/t)
            self.log.debug("Chain# "+str(self.chainNum)+" has random number seed = "+str(self.seed))
            self.settingsDict['chainNum'] = (self.chainNum,"chain number")
            np.random.seed(self.seed)
        self.tmpDataFile = os.path.join(self.dictVal('tmpDir'),self.dictVal('outRoot')+"-"+str(self.chainNum)+".npy")
        if os.path.exists(self.tmpDataFile):
            os.remove(self.tmpDataFile)
            self.log.debug("just removed data file from disk:\n"+self.tmpDataFile)
    
    def startSummary(self,pars,sigs,stage):
        startStr=stage+" chain #"+str(self.chainNum)+' VALS AT START OF '+stage+' SIM:\n'
        startStr+= 'params = '+repr(pars)+'\n'
        startStr+= 'rangeMins = '+repr(self.rangeMins)+'\n'
        startStr+= 'rangeMaxs = '+repr(self.rangeMaxs)+'\n'
        if self.dictVal('lowEcc'):
            startStr+= 'rangeMinsRaw = '+repr(self.rangeMinsRaw)+'\n'
            startStr+= 'rangeMaxsRaw = '+repr(self.rangeMaxsRaw)+'\n'
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
        self.resetTracked(stage)
        bar = tools.ProgressBar('green',width=30,block='=',empty='-',lastblock='>')
        modelData = np.zeros((len(self.realData),3))
        acceptedParams = []
        temp = 1.0
        if (stage=='SA')or(stage=='MC'):
            ## get starting params and sigmas as these two stages start at a random point
            sigmas = copy.deepcopy(self.starterSigmas)
            proposedParsRaw = self.increment(self.rangeMinsRaw,sigmas,stage='MC')
            proposedParsRaw[11]=self.dictVal('chiMAX')*10*self.nu
            if stage=='SA':
                temp=self.dictVal('strtTemp')
        else:
            proposedParsRaw = copy.deepcopy(startParams)
            sigmas = copy.deepcopy(startSigmas)
        latestParsRaw = copy.deepcopy(proposedParsRaw)
        paramsLast = copy.deepcopy(proposedParsRaw)
        ## convert from Raw form if in lowEcc mode
        self.Orbit.convertParsFromRaw(paramsLast)
        self.paramsLast = paramsLast
        #self.startSummary(proposedPars,sigmas,stage)
        ##loop through each sample 
        ##Follows these steps: in Range?,calc model,accept?,Store?,increment,lower temp?,tune sigmas? dump data to disk? DONE ->write output data
        sample=0
        while sample<(self.dictVal(self.stgNsampDict[stage])+1):
            sample+=1
            (proposedPars,inRange)=self.rangeCheck(proposedParsRaw,stage)
            if inRange:
                self.Orbit.calculate(modelData,proposedPars)
                (params,accept) = self.accept(sample,proposedPars,modelData,temp,stage)
                if accept:
                    latestParsRaw = copy.deepcopy(params)
                    ## convert back to sqrt(e)sin(omega), sqrt(e)cos(omega) if in lowEcc mode
                    self.Orbit.convertParsToRaw(latestParsRaw)
                    if ('MCMC' not in stage)and(stage=='MC'):
                        acceptedParams.append(params)
                        self.nSaved+=1
                        self.nSavedPeriodic+=1  
                    elif (self.acceptCount%self.dictVal('saveInt'))==0:
                        acceptedParams.append(params)  
                        self.nSaved+=1   
                        self.nSavedPeriodic+=1                
            else:
                self.acceptBoolAry.append(0)
            proposedParsRaw = self.increment(latestParsRaw,sigmas,stage)
            temp = self.tempDrop(sample,temp,stage)
            sigmas = self.sigTune(sample,sigmas,stage)
            if (self.nSavedPeriodic>0)and((self.nSaved%self.dictVal('dmpInt'))==0):
                ## dump acceptedParams array to disk and collect garbage
                self.log.debug('Dumping data to filename:\n'+self.tmpDataFile+\
                               '\nThe acceptedParams Ary had '+str(np.array(acceptedParams).shape[0])+\
                               ' param sets, and size = '+str(np.array(acceptedParams).nbytes/1.0e6)+' MB')
                tools.periodicDataDump(self.tmpDataFile,np.array(acceptedParams))
                acceptedParams = []
                self.nSavedPeriodic = 0
                self.log.debug('about to collect the garbage')
                gc.collect()
            if (self.dictVal('logLevel')<40)and(sample%(self.dictVal(self.stgNsampDict[stage])//self.dictVal('nSumry'))==0):
                bar.render(sample * 100 // self.dictVal(self.stgNsampDict[stage]), stage+str(chainNum)+' complete so far.')
        if self.dictVal('logLevel')<40:
            bar.render(100,stage+str(chainNum)+' complete!\n')
        self.log.debug(stage+" took: "+tools.timeStrMaker(timeit.default_timer()-tic))
        self.endSummary(temp,sigmas,stage)
        tools.periodicDataDump(self.tmpDataFile,np.array(acceptedParams))
        outFname = tools.writeFits('outputData'+stage+str(chainNum)+'.fits',self.tmpDataFile,self.settingsDict)
        if stage=='ST':
            ##start MCMC at best or end of ST?
            if self.dictVal('strBest'):
                return (self.paramsBest,sigmas,self.bestRedChiSqr)
            else:
                return (latestParsRaw,sigmas,self.bestRedChiSqr)
        elif stage=='SA':
            ##start ST at the best location with tight sigmas, and it will tune to ideal sigmas
            return (self.paramsBest,np.ones(np.array(sigmas).shape)*0.01,self.bestRedChiSqr)
        elif(stage=='MC')or(stage=='MCMC'):
            return outFname
#END OF FILE      