#! /usr/bin/env python

# System imports
from distutils.core import *
from distutils      import sysconfig
from distutils.extension import Extension

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

_hello = Extension("_hello",
		    ["hello_wrap.cxx",
		     "hello.cc"],
		     include_dirs = [numpy_include],
		     )
_orbCalc = Extension("_orbCalc",
		    ["orbCalc_wrap.cxx",
		     "orbCalc.cc"],
		     include_dirs = [numpy_include],
		     )

setup(name= "KylesInitialTests",
	description = "Functions that links with c++",
	author      = "Kyle Mede",
	py_modules  = ["hello","orbCalc"],
	ext_modules = [_hello,_orbCalc])
