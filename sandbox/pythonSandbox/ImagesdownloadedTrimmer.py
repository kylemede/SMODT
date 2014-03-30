import os
import shutil
import glob
import numpy as np
from progressbar import ProgressBar

def main():
    verbose = False 
    #folderList = ['PBblog-2009-10','PBblog-2009-11','PBblog-2009-12','PBblog-2010-01','PBblog-2010-02','PBblog-2010-03','PBblog-2010-04','PBblog-2010-05','PBblog-2010-06']
    years = ['2009','2010','2011','2012','2013']
    months = ["01",'02','03','04','05','06','07','08','09','10','11','12']
    print 'Starting to clean bad images from directories provided\n'
    for year in years:
        for month in months:
            data_dir = '/run/media/Kyle/HOME/shots/PBblog-'+str(year)+'-'+str(month)+'/'
            if os.path.exists(data_dir):
                print 'Cleaning images from: '+data_dir+'\n'
                
                frameFolderlist = np.sort(glob.glob(data_dir + "*.jpg"))
                nframes = len(frameFolderlist)
                
                # create an indexed list of data files
                indexlist = [(i, frameFolderlist[i]) for i in range(nframes)]
                
                # create a progress bar object for updates to user of progress
                p = ProgressBar('red',width=30,block='=',empty='-',lastblock='>')
                
                for i, filename in indexlist:
                    #os.path.getsize(filename)
                    s=filename
                    if (s[s.rfind('x')-1].isdigit())and(s[s.rfind('x')+1].isdigit()):
                        if verbose:
                            print 'Found a bad file: '+filename+'\n'
                        os.remove(filename)
                    # Display progress bar
                    p.render(i * 100 // nframes, 'Complete so far.')
    
    
#############################################################
# end
#############################################################

if __name__ == '__main__':
    main()