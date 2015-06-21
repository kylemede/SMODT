/* file: postCtools.i */
%module postCtools
%{
#include "postCtools.h"
#include "string"
%}

%{
#define SWIG_FILE_WITH_INIT
%}
/*include numpy typemaps*/
%include "numpy.i"
/*initialize module*/
%init %{
import_array();
%}
/*typemaps for different non-arrays*/
//%apply double *OUTPUT {double, *result};
/*typemaps for different arrays*/
//pass in 1D array of dynamic size, that can NOT change in CPP func
apply (double* IN_ARRAY1, int DIM1) {(double *x, int x_n)};
%include "std_string.i"
%include "postCtools.h"
