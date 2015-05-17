#include "orbit.h"

double testFunc(double t){
    std::cout<<"\nInside testFunc"<<std::endl;
    t = t*300;
    return t;
};

void Orbit::loadRealData(double *xx, int xx_nx, int xx_ny){
	if (false)
		std::cout<<"\nInside loadData function"<<std::endl;
    dataRealAry = xx;
    dataRealAry_nx = xx_nx;
    dataRealAry_ny = xx_ny;
    if (true)
    	std::cout<<"data loaded!"<<std::endl;
};
void Orbit::loadConstants(double Grav_in,double pi_in,double MsunToKG_in, double daysPerYear_in){
	Grav = Grav_in;
	pi = pi_in;
	MsunToKG = MsunToKG_in;
	daysPerYear = daysPerYear_in;
	if (true)
		std::cout<<"constants loaded!"<<std::endl;
};

void Orbit::calculate(double *yy, int yy_nx, int yy_ny, double *x, int x_n){
    /*
    The calculator function to perform the primary 
    orbit calculations for the C++ Orbit object.
    */
	dataModelAry=yy;
	dataModelAry_nx=yy_nx;
	dataModelAry_ny=yy_ny;
	params = x;
	params_n = x_n;

	std::cout<<"\nInside Orbit calculator function"<<std::endl;
    //std::cout<<"testDouble = "<<testDouble<<std::endl;
    testDouble = testDouble*10.0;
    //std::cout<<"testDouble = "<<testDouble<<std::endl;
    std::cout<<"For STATIC dataAry:"<<std::endl;
    for (int i=0; i<dataRealAry_nx; i++){
		for (int j=0;j<dataRealAry_ny;j++){
	        std::cout<<"[i,j = ["<<i<<","<<j<<"] = "<<dataRealAry[j+i*dataRealAry_ny]<<std::endl;
		}
	}

    std::cout<<"For STATIC paramsArys:"<<std::endl;
    for (int i=0; i<params_n;i++)
    	std::cout<<"["<<i<<"] = "<<params[i]<<std::endl;

	std::cout<<"For INPLACE dataAryModel:"<<std::endl;
	int ints;
	double ints2;
	for (int i=0; i<dataModelAry_nx; i++){
		for (int j=0;j<dataModelAry_ny;j++){
		ints = i+j;
		//std::cout<<"ints = "<<ints<<std::endl;
		ints2 = (double)ints;
		//std::cout<<"ints2 = "<<ints2<<std::endl;
		dataModelAry[j+i*dataModelAry_ny]=ints2*1.0+1.0;
	        std::cout<<"[i,j] = ["<<i<<","<<j<<"] = "<<dataModelAry[j+i*dataModelAry_ny]<<std::endl;
		}
	}
	//std::cout<<"testFunc provided = "<<testDouble<<std::endl;
	double r = testFunc(testDouble);
	std::cout<<"testFunc returned = "<<r<<std::endl;
	
	//------------------------------- 	START REAL CALCULATE STEPS -------------
	for (int i=0;i<dataRealAry_nx; i++){
		M = (2.0*pi*(dataRealAry[0+i*dataRealAry_ny]-2.0*))








	}


}


















