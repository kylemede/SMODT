Packages installed on top of standard Fedora 21 libraries.
-gcc compiler for both c and c++ (although we only use c++...)
($sudo yum install gcc) 
($sudo yum install gcc-c++) 
-psutil
available at (https://pypi.python.org/pypi/psutil)
or use command $sudo pip install psutil (in su)
-scipy 
($sudo pip install scipy)
-swig 
($sudo yum install swig)

NOTE: Need sudo for swig compiling --> '$ sudo make' or '$ sudo make clean' for cpp swig code!!!  I guess this will come into the setup.py area.
      NOT TRUE if the user has read/write privilages for the SMODT dirs/files!!!