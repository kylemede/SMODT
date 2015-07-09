import os

def nightRun():
    """
    To run a set of SMODT runs in a series to get more done at night.
    """
    cmdList = [
               "python SMODT.py Jupiter1percent_",
               "python SMODT.py Jupiter10percent_",
               "python SMODT.py Jupiter5percent_",
               "python customPost.py"
               #"python SMODT.py HIP10321_",
               ]
    for runCmd in cmdList:
        print "About to send command to system:\n"+runCmd
        os.system(runCmd)
        print runCmd+" Finished!! :-)"
        
if __name__ == '__main__':
    nightRun()