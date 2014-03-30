import os
import shutil
import glob
import numpy as np
from progressbar import ProgressBar

def main():
    
    data_dir = '/mnt/Data1/Todai_Work/Dropbox/RecoveredFiles/'
    print 'Working on directory: '+data_dir
    
    filelist = np.sort(glob.glob(data_dir + "* (deleted*).*"))
    nfiles = len(filelist)
    
    # create an indexed list of data files
    indexlist = [(i, filelist[i]) for i in range(nfiles)]
    
    # create a progress bar object for updates to user of progress
    p = ProgressBar('red',width=30,block='=',empty='-',lastblock='>')
    
    for i, file in indexlist:
        
        newname = file.split(" (deleted")[0]+os.path.splitext(file)[1]
        print '\nold name: '+file+"\n new name: "+newname+'\n'
        os.rename(file,newname)
        
        # Display progress bar
        p.render(i * 100 // nfiles, 'Complete so far.')
    
    
#############################################################
# end
#############################################################

if __name__ == '__main__':
    main()