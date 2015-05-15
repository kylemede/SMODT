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
    double* dataRealAry;
    int dataRealAry_x;
    int dataRealAry_y;
    
    //funcs
    //for passing in a STATIC 2D data array of Real data into the object
    void loadRealData(double *xx, int nx, int ny);
    //to calculate the model data and load it into an empty 2d array
    void calculate(double  *yy, int nx, int ny);
};
double testFunc(double t);


#endif
