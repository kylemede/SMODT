//@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
#if !defined(POSTCTOOLS_H)
#define POSTCTOOLS_H 1
#include<string>
#include <math.h>
#include <stdio.h>
#include <iostream>

class PostCtools{
public:
	//internal params
	int i,i_last,j,n;
	double s,ep;
	double ave,var,sum;
	double varALL,halfVarALL,numCorrLengths,corrLengthsTotal;
	//input data ary
	double* data;
	int data_nx;
	//funcs
	void loadAry(double *x, int x_nx);
	void sumCalc(int startPoint,int lastPoint);
	void meanCalc(int startPoint, int lastPoint);
	void varianceCalc(int startPoint, int lastPoint);
	double corrLenCalc();
};
#endif
