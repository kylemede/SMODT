import os

def nightRun():
    """
    To run a set of ExoSOFT runs in a series to get more done at night.
    """
    cmdList = [
               #"python ExoSOFT.py HR7672_",
               #"python ExoSOFT.py HR7672-1_",
               #"python ExoSOFT.py HD129333_",
               #"python ExoSOFT.py HD129333-1_",
               "python ExoSOFT.py GJ1108A_",
               "python ExoSOFT.py GJ1108A-1_",
               #"python ExoSOFT.py HIP10321-1_",
               #"python customPost.py HIP10321-1_",
               #"python ExoSOFT.py HIP10321-2_",
               #"python ExoSOFT.py Jupiter10percent_",
               #"python ExoSOFT.py Jupiter1percent_",
               #"python ExoSOFT.py HIP42074_",
               #"python ExoSOFT.py HIP42074-1_",
               #"python ExoSOFT.py HIP42074-2_",
               #"python ExoSOFT.py JupiterOmegaTest_",
               #"python customPost.py HR8799e_"
               ]
    for runCmd in cmdList:
        print "About to send command to system:\n"+runCmd
        os.system(runCmd)
        print runCmd+" Finished!! :-)"
        
if __name__ == '__main__':
    nightRun()