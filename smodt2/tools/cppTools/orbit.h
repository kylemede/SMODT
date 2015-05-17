#if !defined(ORBIT_H)
#define ORBIT_H 1
#include<string>
#include <math.h>
#include <stdio.h>
#include <iostream>

class Orbit{
public:
    //obj variables
    double testDouble;
    int NewtonCount;
    //declare REAL static data array
    double* dataRealAry;
    int dataRealAry_nx;
    int dataRealAry_ny;
    //declare MODEL in place data array
	double* dataModelAry;
	int dataModelAry_nx;
	int dataModelAry_ny;
	double* params;
	int params_n;
    //Declare all doubles used inside Orbit class.
    double M, E, Eprime,K,theta,thetaPrime,A,B,F,G,X,Y,atot;
    //Declare input varying parameters
    double P,e,omega,Omega,mass1,mass2,inc,To,Tc,dist;
    //Declare static global constants
    double Grav,pi,KGperMsun,daysPerYear,secPerYear,MperAU;
    
    //funcs
    //for passing in a STATIC 2D data array of Real data into the object
    void loadRealData(double *xx, int xx_nx, int xx_ny);
    //For loading in the global constants
    void loadConstants(double Grav_in,double pi_in,double KGperMsun_in, double daysPerYear_in,double secPerYear_in,double MperAU_in);
    //to calculate the model data and load it into an empty 2d array and 1d params array (both inplace arrays)
    void calculate(double *yy, int yy_nx, int yy_ny, double *y, int y_n);
    //void calculate(double *dataModelAry, int dataModelAry_nx, int dataModelAry_ny, double *params, int params_n);
};
double testFunc(double t);


#endif
