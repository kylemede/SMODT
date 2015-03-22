import os
import shutil
import glob
import numpy as np
from progressbar import ProgressBar

def main():
    
    data_dir = '/run/media/Kyle/Data1/Todai_Work/Data/data_GDPS/2B/'
    print 'Working on directory: '+data_dir
    
    frameFolderlist = np.sort(glob.glob(data_dir + "*?logkey=pkg%0D"))
    nframes = len(frameFolderlist)
    
    # create an indexed list of data files
    indexlist = [(i, frameFolderlist[i]) for i in range(nframes)]
    
    # create a progress bar object for updates to user of progress
    p = ProgressBar('red',width=30,block='=',empty='-',lastblock='>')
    
    for i, frameFolder in indexlist:
        
        newFolderName = frameFolder.split("?")[0]+'.fits.gz'
        os.rename(frameFolder,newFolderName)
        
        # Display progress bar
        p.render(i * 100 // nframes, 'Complete so far.')
    
    
#############################################################
# end
#############################################################

if __name__ == '__main__':
    main()