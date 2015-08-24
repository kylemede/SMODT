import os

def nightRun():
    """
    To run a set of SMODT runs in a series to get more done at night.
    """
    cmdList = [
               "python NewBEAT.py HIP10321_",
               "python NewBEAT.py HIP10321-2_",
               #"python NewBEAT.py Jupiter5percent_",
               #"python NewBEAT.py Jupiter1percent_",
               #"python NewBEAT.py Jupiter10percent_",
               #"python customPost.py",
               #"python NewBEAT.py HIP10321_",
               ]
    for runCmd in cmdList:
        print "About to send command to system:\n"+runCmd
        os.system(runCmd)
        print runCmd+" Finished!! :-)"
        
if __name__ == '__main__':
    nightRun()