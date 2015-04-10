import pylab
import os
import numpy as np
from subprocess import Popen
plt = pylab.matplotlib.pyplot

def starter(paramSettingsDict,sleep=1):
    FNULL=open(os.devnull,'w')
    toolDir = os.path.join(paramSettingsDict['pythonCodeDir'],'toolboxes')
    shScriptPath = os.path.join(toolDir,'ramLogger.sh')
    memLogFilename = os.path.join(paramSettingsDict['outputData_dir'],"RAMusage.log")
    totSamples = paramSettingsDict["numSamples"]*paramSettingsDict['numProcesses']
    if False:
        print "memLogCleanerAndPlot starter found totalSamples = "+str(totSamples)
    sleepUse = sleep
    if totSamples<500000001:
        sleepUse = 1
    elif (totSamples>500000000)and(totSamples<50000000000):
        sleepUse = 6 
    else:
        sleepUse = 60
    memTracProc = Popen([shScriptPath,memLogFilename,str(sleepUse)],stdout=FNULL)
    return (memTracProc,memLogFilename,sleepUse)

def wrapUp(proc,memLogFilename,sleep=6):
    proc.terminate()
    maxUse = memUsageLogCleaner(memLogFilename,sleep)
    return maxUse
    
def memUsageLogCleaner(filename = '',sleep=1,delOrigLog=True):
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
    totalRAM = 0
    for line in lines:
        if (line[0]=="M")and(totalRAM==0):
            #print line
            l = line.split()[1:3]
            #print repr(l)
            totalRAM = float(l[0])
        if (line[0]=="-")and(totalRAM>1):
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
        subPlot.axes.set_ylim(yLim)
        subPlot.axes.set_ylabel("Total System RAM used [%]",fontsize=15)
        subPlot2 = fig.add_subplot(212)
        used2 = used-used.min()
        subPlot2.plot(times,used2)
        usedRange = used2.max()-used2.min()
        yLim2 = [used2.min()-0.05*usedRange,used2.max()+0.05*usedRange]
        subPlot2.axes.set_ylim(yLim2)
        subPlot2.axes.set_ylabel("(During - Before) RAM usage [MB]",fontsize=12)
        subPlot2.axes.set_xlabel("Time from simulation start in "+strmod,fontsize=15)
        maxUse = used.max()-used.min()
        subPlot2.text(0.05,usedRange*0.8,'Max RAM used\n      '+str(maxUse),ha='left',fontsize=20)
        plotname = os.path.abspath(filename)[:-4]+"_clean.png"
        plt.savefig(plotname, orientation='landscape')
    print "RAMusage plot written to: "+plotname
    
    print "Max RAM used during simulation was "+str(maxUse)+" MB"
    return maxUse
    
if __name__ == '__main__':
    memUsageLogCleaner() 
