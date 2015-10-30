import os

def nightRun():
    """
    To run a set of ExoSOFT runs in a series to get more done at night.
    """
    cmdList = [
               #"python ExoSOFT.py HR7672_",
               #"python ExoSOFT.py HR7672-1_",
               #"python ExoSOFT.py HIP10321_",
               #"python customPost.py HIP10321_",
               #"python ExoSOFT.py HIP10321-1_",
               #"python customPost.py HIP10321-1_",
               #"python ExoSOFT.py HIP10321-2_",
               "python ExoSOFT.py Jupiter5percent_",
               "python ExoSOFT.py Jupiter1percent_",
               "python ExoSOFT.py Jupiter10percent_",
               ]
    for runCmd in cmdList:
        print "About to send command to system:\n"+runCmd
        os.system(runCmd)
        print runCmd+" Finished!! :-)"
        
if __name__ == '__main__':
    nightRun()