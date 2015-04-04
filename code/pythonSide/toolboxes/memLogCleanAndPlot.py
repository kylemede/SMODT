import pylab
import os
import numpy as np
plt = pylab.matplotlib.pyplot

def memUsageLogCleaner(filename = ''):
    """
    A function I whiped together to clean up a RAM usage log produced with a super simple bash script.
    """ 
    if filename=="":
        filename = "/mnt/Data1/Todai_Work/Data/data_SMODT/RAMusage.log"
    print 'Input RAMusage log file: '+filename
    f = open(filename,'readonly')
    fnamOut = os.path.abspath(filename)[:-4]+"_clean.log"
    fOut = open(fnamOut,'w')
    lines = f.readlines()
    f.close()
    fOut.write("total[MB]   Used[%] \n")
    mem = []
    for line in lines:
        if line[0]=="M":
            #print line
            l = line.split()[1:3]
            #print repr(l)
            usedPercent = int((float(l[1])/float(l[0]))*100.)
            lineOut = "  "+str(l[0])+"        "+str(usedPercent)+"\n"
            fOut.write(lineOut)
            mem.append(usedPercent)
            #print repr(lineOut)
    fOut.close()
    print "cleaned RAMusage log file written to: "+fnamOut
    if True:
        times = np.arange(len(mem))*(1.0/6.0)
        fig = plt.figure(1, figsize=(15,10),dpi=200)
        subPlot = fig.add_subplot(111)
        subPlot.plot(times,mem)
        subPlot.axes.set_ylabel("Percent RAM usage")
        subPlot.axes.set_xlabel("time [hrs]")
        plotname = os.path.abspath(filename)[:-4]+"_clean.png"
        plt.savefig(plotname, orientation='landscape')
    print "RAMusage plot written to: "+plotname
    
if __name__ == '__main__':
    memUsageLogCleaner() 
