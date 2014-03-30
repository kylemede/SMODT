import threading
from threadTestFunc import threadTestFUNC

threadNum = 0
numRuns = 1

class MyThread ( threading.Thread ):

    def run ( self):
        global threadNum
        global numRuns
        print '\nthread ' + str ( threadNum ) + 'started'
        threadTestFUNC(numRuns)
        print 'thread '+str(threadNum)+' finished'
        
    def setVars(self, threadNum1, numRuns1):
        global threadNum
        global numRuns
        threadNum=threadNum1
        numRuns = numRuns1
        
        
for x in range (0,20 ):
    
   Thread = MyThread()
   Thread.setVars(x, int(1e3))
   MyThread().start()