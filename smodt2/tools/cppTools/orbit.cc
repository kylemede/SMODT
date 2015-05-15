#include "orbit.h"

double testFunc(double t){
    std::cout<<"\nInside testFunc"<<std::endl;
    t = t*300;
    return t;
}

void Orbit::loadRealData(double *xx, int nx, int ny){
    std::cout<<"\nInside loadData function"<<std::endl;
    dataRealAry = xx;
    dataRealAry_x = nx;
    dataRealAry_y = ny;
    std::cout<<"data loaded!"<<std::endl;
};

void Orbit::calculate(double  *yy, int nx, int ny){
    /*
    The calculator function to perform the primary 
    orbit calculations for the C++ Orbit object.
    */

	std::cout<<"\nInside Orbit calculator function"<<std::endl;
    //std::cout<<"testDouble = "<<testDouble<<std::endl;
    testDouble = testDouble*10.0;
    //std::cout<<"testDouble = "<<testDouble<<std::endl;
    std::cout<<"For STATIC dataAry:"<<std::endl;
    for (int i=0; i<dataRealAry_x; i++){
		for (int j=0;j<dataRealAry_y;j++){
	        std::cout<<"[i,j = ["<<i<<","<<j<<"] = "<<dataRealAry[j+i*dataRealAry_y]<<std::endl;
		}
	}
	std::cout<<"For INPLACE yy:"<<std::endl;
	int ints;
	double ints2;
	for (int i=0; i<nx; i++){
		for (int j=0;j<ny;j++){
		ints = i+j;
		//std::cout<<"ints = "<<ints<<std::endl;
		ints2 = (double)ints;
		//std::cout<<"ints2 = "<<ints2<<std::endl;
		yy[j+i*ny]=ints2*1.0+1.0;
	        std::cout<<"[i,j] = ["<<i<<","<<j<<"] = "<<yy[j+i*ny]<<std::endl;
		}
	}
	//std::cout<<"testFunc provided = "<<testDouble<<std::endl;
	double r = testFunc(testDouble);
	std::cout<<"testFunc returned = "<<r<<std::endl;
	
}


