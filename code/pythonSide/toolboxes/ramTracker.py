import pylab
import os
import numpy as np
from subprocess import Popen
plt = pylab.matplotlib.pyplot

class RAMtracker:
    """
    An object to start, and wrap up RAM usage logging.
    """
   
    def __init__(self,paramSettingsDict,sleep=1):
        self.paramSettingsDict = paramSettingsDict
        self.sleep = sleep
        (memTrackProc,memLogFilename) = self.starter()
        self._memLogFilename = memLogFilename
        self._memTrackProc = memTrackProc
           
    def starter(self):
        """
        The ramTracker module will start running the ramLogger.sh script as a
        background process to log the RAM usage while SMODT is running. 
        The wrapUp function can then be called to terminate the process, then 
        it calls logCleanAndPlot to clean up the log file and plot the results.
        """
        FNULL=open(os.devnull,'w')
        toolDir = os.path.join(self.paramSettingsDict['pythonCodeDir'],'toolboxes')
        shScriptPath = os.path.join(toolDir,'ramLogger.sh')
        memLogFilename = os.path.join(self.paramSettingsDict['outputData_dir'],"RAMusage.log")
        totSamples = self.paramSettingsDict["numSamples"]*self.paramSettingsDict['numProcesses']
        if False:
            print "memLogCleanerAndPlot starter found totalSamples = "+str(totSamples)
        if totSamples<500000001:
            self.sleep = 1
        elif (totSamples>500000000)and(totSamples<50000000000):
            self.sleep = 6 
        else:
            self.sleep = 60
        memTrackProc = Popen([shScriptPath,memLogFilename,str(self.sleep)],stdout=FNULL)
        return (memTrackProc,memLogFilename)

    def chainsDonePrint(self):#,memLogFilename):
        """
        Add a print to log to locate where C++ chains finished.
        """
        f = open(self._memLogFilename,'a')
        f.write('\n'+"*"*50+"\nSMODT chains ENDED NOW!!!!\n"+"*"*50+'\n')
        
    def wrapUp(self):#,memTracProc,memLogFilename):
        """
        Terminates the background process used by ramLogger.sh, , then 
        it calls logCleanAndPlot to clean up the log file and plot the results.
        """
        self._memTrackProc.terminate()
        maxUse = logCleanAndPlot(self._memLogFilename,self.sleep)
        return maxUse
    
def logCleanAndPlot(filename = '',sleep=1,delOrigLog=True):
    """
    A function I whipped together to clean up a RAM usage log produced with a super simple bash script.
    """ 
    if filename=="":
        filename = "/mnt/Data1/Todai_Work/Data/data_SMODT/RAMusage.log"
    
    print "\n\nRunning RAM usage log clean up and plotter\n"
    print 'Input RAMusage log file: '+filename
    f = open(filename,'readonly')
    fnamOut = os.path.abspath(filename)[:-4]+"_clean.log"
    fOut = open(fnamOut,'w')
    lines = f.readlines()
    f.close()
    fOut.write("total[MB]   Used[%] \n")
    used = []
    mem = []
    mcmcEndCounter = 0
    mcmcEnded = False
    totalRAM = 0
    for line in lines:
        if ("Mem" in line)and(totalRAM==0):
            #print line
            l = line.split()[1:3]
            #print repr(l)
            totalRAM = float(l[0])
        if ("ENDED"in line) and (mcmcEnded==False):
            mcmcEnded=True
            print "SMODT chains ended at iteration # "+str(mcmcEndCounter)
        if (line[0]=="-")and(totalRAM>1):
            if (mcmcEnded==False):
                mcmcEndCounter+=1
            usedRAM = float(line.split(":")[-1].split()[0])
            usedPercent = int((usedRAM/totalRAM)*100.)
            if False:
                print 'line = '+line
                print 'usedRAM = '+str(usedRAM)
                print 'usedPercent = '+str(usedPercent)
            lineOut = "  "+str(totalRAM)+"        "+str(usedPercent)+"\n"
            used.append(usedRAM)
            fOut.write(lineOut)
            mem.append(usedPercent)
            #print repr(lineOut)
    fOut.close()
    
    print "cleaned RAMusage log file written to: "+fnamOut
    if delOrigLog:
        os.remove(filename)
        if os.path.exists(filename):
            print "ERROR occurred while trying to delete :"+filename
        else:
            print "Original RAMusage log deleted: "+filename
    if True:
        multmod=sleep/3600.0
        strmod= "[hrs]"
        if ((sleep*len(mem))/3600.0)<1.0:
            multmod=sleep/60.0
            strmod="[minutes]"
        times = np.arange(len(mem))*multmod
        if False:
            print "len(times) = "+repr(times.size)
            print "repr(times) = "+repr(times)
            print "len(mem) = "+repr(len(mem))
            print "repr(mem) = "+repr(mem)
        used= np.array(used)
        mem = np.array(mem)
        fig = plt.figure(1, figsize=(15,10),dpi=200)
        subPlot = fig.add_subplot(211)
        subPlot.plot(times,mem)
        plt.suptitle("RAM usage Information *during* SMODT Run\nNot necessarily solely due to SMODT", fontsize=20)
        memRange = mem.max()-mem.min()
        yLim = [mem.min()-0.1*memRange,mem.max()+0.1*memRange]
        if False:
            print '\n\nadding plot line for when C++ ended using : \nx = '+repr([times[mcmcEndCounter],times[mcmcEndCounter]])+\
                "\ny = "+repr([yLim[0],yLim[1]])+"\nmemRange = "+repr(memRange)+'\nmem = '+repr(mem)+"\n\n"
        subPlot.plot([times[mcmcEndCounter],times[mcmcEndCounter]],[yLim[0],yLim[1]],color='red')
        subPlot.axes.set_ylim(yLim)
        subPlot.axes.set_ylabel("Total System RAM used [%]",fontsize=15)
        timesRange = times[-1]-times[0]
        textX = times[mcmcEndCounter]
        textY = yLim[1]-(yLim[1]-yLim[0])*0.1
        hereStr = "<----------------- HERE "
        #print "\n\n0.2*timesRange = "+str(0.2*timesRange)+"\n\n"
        if textX>0.2*timesRange:
            dashMult = int(textX-(timesRange*0.025+times[0]))
            textX = timesRange*0.025+times[0]
            hereStr = " HERE ------------"+"-"*dashMult+"----->"
        redTextStr = ' Chains finished and\n Post-processing started\n'+hereStr
        print '\nredTextStr = \n'+redTextStr
        print "dashMult = "+str(dashMult)
        print 'textY = '+str(textY)
        subPlot.text(textX,textY,redTextStr,ha='left',fontsize=17,color="red")
        subPlot2 = fig.add_subplot(212)
        used2 = used-used.min()
        subPlot2.plot(times,used2)
        usedRange = used2.max()-used2.min()
        yLim2 = [used2.min()-0.05*usedRange,used2.max()+0.05*usedRange]
        subPlot2.plot([times[mcmcEndCounter],times[mcmcEndCounter]],[yLim2[0],yLim2[1]],color='red')
        subPlot2.axes.set_ylim(yLim2)
        subPlot2.axes.set_ylabel("(During - Before) RAM usage [MB]",fontsize=12)
        subPlot2.axes.set_xlabel("Time from simulation start in "+strmod,fontsize=15)
        maxUse = used.max()-used.min()
        subPlot2.text(textX,yLim2[1]*0.78,' Max RAM used\n     '+str(maxUse)+" MB",ha='left',fontsize=20)
        plotname = os.path.abspath(filename)[:-4]+"_clean.png"
        plt.savefig(plotname, orientation='landscape')
    print "RAMusage plot written to: "+plotname
    
    print "Max RAM used during simulation was "+str(maxUse)+" MB"
    return maxUse
    
if __name__ == '__main__':
    logCleanAndPlot() 
