import numpy as np
import os
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
        ##check there are 
        if np.max(self.realData[:,-1])!=len(self.dictVal('vMINs')):
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
        #print 'paramInts = '+repr(paramInts)
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
    
    def increment(self,params,sigmas=[],mcOnly=False):
        """
        Randomly increment one of the parameters
        """
        #params order =[Mass1,Mass2,Sys_Dist_PC,Omega,e,T,T_center,P,inc,omega,aTot,chiSquared,K,offset]
        
        #print 'params IN = '+repr(params)
        ## vary all the params if mcONLY
        if mcOnly:
            for i in range(0,len(params)):
                if i in self.paramInts:
                    params[i]=np.random.uniform(self.rangeMins[i],self.rangeMaxs[i])
        ## vary a random param if SA, ST or MCMC
        else:
            varyInt = self.paramInts[np.random.randint(0,len(self.paramInts))]
            params[varyInt]=np.random.uniform(params[varyInt]-sigmas[varyInt],params[varyInt]+sigmas[varyInt])
            #print 'varyInt = '+str(varyInt)
            
        ## if TcEqualT, push the varied one into the other
        if self.dictVal('TcEqualT'):
            if self.dictVal('TcStep'):
                params[5]=params[6]
            else:
                params[6]=params[5]
                
        #print 'params OUT = '+repr(params)
        return params
    
    def accept(self,params,modelData,temp=1.0,mcOnly=False):
        """
        First this will calculate chi squared for model vs real data.
        
        For mcOnly it performs simple chisquared cut-off acceptance 
        based on 'chiMAX' value in settingsDict.
        Else, it will calculate the priors and accept based on 
        the Metropolis-Hastings algorithm. The temp factor will 
        be set to 1.0 for MCMC and Sigma Tuning, and should be provided 
        for Simulated Annealing.
        """
        diffs = np.concatenate(((self.realData[:,1]-modelData[:,0]),(self.realData[:,3]-modelData[:,1]),(self.realData[:,5]-modelData[:,2])))
        errors = np.concatenate((self.realData[:,2],self.realData[:,4],self.realData[:,6]))
        params[11] = np.sum((diffs**2)/(errors**2))
        diffsDI = np.concatenate(((self.realData[:,1]-modelData[:,0]),(self.realData[:,3]-modelData[:,1])))
        errorsDI = np.concatenate((self.realData[:,2],self.realData[:,4]))
        diffsRV = (self.realData[:,5]-modelData[:,2])
        errorsRV = self.realData[:,6][np.where(diffsRV!=0)]
        #print 'diffsDI = \n'+repr(diffsDI)
        #print 'errorsDI = \n'+repr(errorsDI)
        #print 'diffsRV = \n'+repr(diffsRV)
        #print 'errorsRV = \n'+repr(errorsRV)
        chiSquaredDI = np.sum((diffsDI[np.where(diffsDI!=0)]**2)/(errorsDI[np.where(diffsDI!=0)]**2))
        chiSquaredRV = np.sum((diffsRV[np.where(diffsRV!=0)]**2)/(errorsRV**2))
        #print 'chiSquared = '+str(params[11])
        #print 'chiSquaredDI = '+str(chiSquaredDI)
        #print 'chiSquaredRV = '+str(chiSquaredRV)
        #print 'nu,nuDI,nuRV = '+str(self.nu)+", "+str(self.nuDI)+", "+str(self.nuRV)
        #print 'reduced chiSquared: 3D, DI, RV = '+str(params[11]/self.nu)+", "+str(chiSquaredDI/self.nuDI)+", "+str(chiSquaredRV/self.nuRV)
        #print 'chiSquaredDI reduced = '+str(chiSquaredDI/self.nuDI)
        #print 'chiSquaredRV reduced = '+str(chiSquaredRV/self.nuRV)
        #print 'chiSquared reduced = '+str(params[11]/self.nu)
        
        if False:
            print ''
            print 'reals = \n'+repr(self.realData[:,[1,3,5]])
            print 'modelData = \n'+repr(modelData)
            print 'real-model = \n'+repr(self.realData[:,[1,3,5]]-modelData)
            print 'diffs = \n'+repr(diffs)
            print 'errors = \n'+repr(errors)
            print '(diffs**2)/(errors**2) = \n'+repr((diffs**2)/(errors**2))
            print 'chiSquared = '+str(params[11])
        accept = False
        if mcOnly:
            if (params[11]/self.nu)<self.dictVal('chiMAX'):
                accept=True
        else:
            #handle case where doing SA and nothing accepted yet
            if (temp!=0)and(self.paramsLast==0):
                if (params[11]/self.nu)<self.dictVal('chiMAX'):
                    accept=True
            else:
                likelihoodRatio = np.e**((self.paramsLast[11] - params[11])/ (2.0*temp))
                ###### put all prior funcs together in dict??
                priorsRatio = (self.dictVal(ePrior)(params[4])/self.dictVal(ePrior)(self.paramsLast[4]))
                priorsRatio*= (self.dictVal(pPrior)(params[7])/self.dictVal(pPrior)(self.paramsLast[7]))
                priorsRatio*= (self.dictVal(incPrior)(params[8])/self.dictVal(incPrior)(self.paramsLast[8]))
                priorsRatio*= (self.dictVal(mass1Prior)(params[0])/self.dictVal(mass1Prior)(self.paramsLast[0]))
                priorsRatio*= (self.dictVal(mass2Prior)(params[1])/self.dictVal(mass2Prior)(self.paramsLast[1]))
                priorsRatio*= (self.dictVal(distPrior)(params[2])/self.dictVal(distPrior)(self.paramsLast[2])) 
                if np.random.uniform(0.0, 1.0)<=(priorsRatio*likelihoodRatio):
                    accept = True
        if accept:
            self.paramsLast=params
        #print 'accept: reduced chiSquared: 3D, DI, RV = '+repr(accept)+": "+str(params[11]/self.nu)+", "+str(chiSquaredDI/self.nuDI)+", "+str(chiSquaredRV/self.nuRV)
        return (params,accept)
    def tempDrop(self,i,temp,mode=''):
        """
        Determine if it is time to drop the temp, and drop if it is.
        Total temperature range is [strtTemp,0.1), so the minimum 
        temperature is actually <1.0 meaning the last few temperature drops 
        will really push the currently found minimum towards its peak.
        There will be a fixed number of temperature steps = 'nTemps'.
        """
        if mode=='SA':
            if i%(self.dictVal('nSAsamp')//self.dictVal('nTemps')):
                temp-=(self.dictVal('strtTemp')+0.9)*(1.0/self.dictVal('nTemps'))
        return temp
        
    def monteCarlo(self):
        """
        Performs 'shotgun' Monte Carlo.
        Returns ?? not decided yet!!!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$ 
        """
        self.log.info("In Simulator.monteCarlo")
        self.log.info("Trying "+str(self.dictVal('nSamples'))+" samples")
        print "Trying "+str(self.dictVal('nSamples'))+" samples"
        bar = tools.ProgressBar('green',width=30,block='=',empty='-',lastblock='>')
        modelData = np.zeros((self.realData.shape[0],3))
        acceptedParams = []
        e = 0.4
        Sys_Dist_PC = 5.0
        Mass1 = 1.0
        Mass2 = 0.2
        Omega = 60.0
        omega = 110.0
        T = 2457000.0
        T_center = 2457000.0
        P = 15.0
        inc =  30.0
        offset = 0.0
        params = np.array([Mass1,Mass2,Sys_Dist_PC,Omega,e,T,T_center,P,inc,omega,0,0,0,offset])
        tic=timeit.default_timer()
        nsamp = self.dictVal('nSamples')
        for i in range(0,self.dictVal('nSamples')):
            self.Orbit.calculate(modelData,params)
            (params,accept) = self.accept(params,modelData,mcOnly=True)
            if accept:
                acceptedParams.append(params)
            params = self.increment(params,mcOnly=True)
            if i%(self.dictVal('nSamples')//20)==0:
                bar.render(i * 100 // self.dictVal('nSamples'), 'Complete so far.')
        bar.render(100, 'Complete so far.')
        toc=timeit.default_timer()
        swigTime=toc-tic
        print "\nMC it took: "+str(swigTime)
        acceptedParams = np.array(acceptedParams)
        #print '\nmodelData = \n'+repr(acceptedParams)
        print 'number accepted = '+str(acceptedParams.shape[0])
        baseFilename = 'testDataMC.fits'
        tools.writeFits(baseFilename,acceptedParams,self.settingsDict)
        #self.log.info('\nmodelData = \n'+repr(modelData))
        
        
    def simAnneal(self):
        """
        Performs Simulated Annealing.
        Returns ?? not decided yet!!!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$ 
        """
        self.log.info("In Simulator.simAnneal")
        self.log.info("Trying "+str(self.dictVal('nSamples'))+" samples")
        print "Trying "+str(self.dictVal('nSamples'))+" samples"
        bar = tools.ProgressBar('green',width=30,block='=',empty='-',lastblock='>')
        modelData = np.zeros((self.realData.shape[0],3))
        acceptedParams = []
        e = 0.4
        Sys_Dist_PC = 5.0
        Mass1 = 1.0
        Mass2 = 0.2
        Omega = 60.0
        omega = 110.0
        T = 2457000.0
        T_center = 2457000.0
        P = 15.0
        inc =  30.0
        offset = 0.0
        params = np.array([Mass1,Mass2,Sys_Dist_PC,Omega,e,T,T_center,P,inc,omega,0,0,0,offset])
        tic=timeit.default_timer()
        nsamp = self.dictVal('nSamples')
        temp = self.dictVal('strtTemp')
        for i in range(0,self.dictVal('nSamples')):
            self.Orbit.calculate(modelData,params)
            (params,accept) = self.accept(params,modelData,temp=temp)
            if accept:
                acceptedParams.append(params)
            params = self.increment(params,sigmas=self.starterSigmas)
            temp = self.tempDrop(i, temp, mode='SA')
            if i%(self.dictVal('nSamples')//20)==0:
                bar.render(i * 100 // self.dictVal('nSamples'), 'Complete so far.')
        bar.render(100, 'Complete so far.')
        toc=timeit.default_timer()
        swigTime=toc-tic
        print "\nSA it took: "+str(swigTime)
        acceptedParams = np.array(acceptedParams)
        #print '\nmodelData = \n'+repr(acceptedParams)
        print 'number accepted = '+str(acceptedParams.shape[0])
        baseFilename = 'testDataSA.fits'
        tools.writeFits(baseFilename,acceptedParams,self.settingsDict)
        #self.log.info('\nmodelData = \n'+repr(modelData))

    def sigmaTuning(self):
        """
        Performs Sigma Tuning.
        Returns ?? not decided yet!!!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$ 
        """
        self.log.info("In Simulator.sigmaTuning")
        self.log.info("Trying "+str(self.dictVal('nSamples'))+" samples")
        print "Trying "+str(self.dictVal('nSamples'))+" samples"
        bar = tools.ProgressBar('green',width=30,block='=',empty='-',lastblock='>')
        modelData = np.zeros((self.realData.shape[0],3))
        acceptedParams = []
        e = 0.4
        Sys_Dist_PC = 5.0
        Mass1 = 1.0
        Mass2 = 0.2
        Omega = 60.0
        omega = 110.0
        T = 2457000.0
        T_center = 2457000.0
        P = 15.0
        inc =  30.0
        offset = 0.0
        params = np.array([Mass1,Mass2,Sys_Dist_PC,Omega,e,T,T_center,P,inc,omega,0,0,0,offset])
        tic=timeit.default_timer()
        sigmas=self.starterSigmas
        for i in range(0,self.dictVal('nSamples')):
            self.Orbit.calculate(modelData,params)
            (params,accept) = self.accept(params,modelData)
            if accept:
                acceptedParams.append(params)
            params = self.increment(params)
            ##tune sigma HERE!!!!!!
            if i%(self.dictVal('nSamples')//20)==0:
                bar.render(i * 100 // self.dictVal('nSamples'), 'Complete so far.')
        bar.render(100, 'Complete so far.')
        toc=timeit.default_timer()
        swigTime=toc-tic
        print "\nSA it took: "+str(swigTime)
        acceptedParams = np.array(acceptedParams)
        #print '\nmodelData = \n'+repr(acceptedParams)
        print 'number accepted = '+str(acceptedParams.shape[0])
        baseFilename = 'testDataST.fits'
        tools.writeFits(baseFilename,acceptedParams,self.settingsDict)
        #self.log.info('\nmodelData = \n'+repr(modelData))       
 
    def mcmc(self):
        """
        Performs pure Markov Chain Monte Carlo.
        Returns ?? not decided yet!!!!! $$$$$$$$$$$$$$$$$$$$$$$$$$$ 
        """
        self.log.info("In Simulator.mcmc")
        