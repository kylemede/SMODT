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

_orbit = Extension("_orbit",
            ["orbit_wrap.cxx",
             "orbit.cc"],
             include_dirs = [numpy_include],
             )
_postctools = Extension("_postctools",
            ["postctools_wrap.cxx",
             "postctools.cc"],
             include_dirs = [numpy_include],
             )

setup(name= "cppStuffSetUp",
	description = "Object for calculating predicted orbit, and some post-processing tools with c++",
	author      = "Kyle Mede",
	py_modules  = ["orbit","postctools"],
	ext_modules = [_orbit,_postctools])
