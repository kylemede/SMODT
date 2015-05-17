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
void Orbit::loadConstants(double Grav_in,double pi_in,double KGperMsun_in, double daysPerYear_in,double secPerYear_in,double MperAU_in){
	Grav = Grav_in;
	pi = pi_in;
	KGperMsun = KGperMsun_in;
	daysPerYear = daysPerYear_in;
	secPerYear = secPerYear_in;
	MperAU = MperAU_in;
	if (true)
		std::cout<<"constants loaded!"<<std::endl;
};

void Orbit::calculate(double *yy, int yy_nx, int yy_ny, double *y, int y_n){
    /*
    The calculator function to perform the primary
    orbit calculations for the C++ Orbit object.
    */
	dataModelAry=yy;
	dataModelAry_nx=yy_nx;
	dataModelAry_ny=yy_ny;
	params = y;
	params_n = y_n;

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
//	//std::cout<<"testFunc provided = "<<testDouble<<std::endl;
//	double r = testFunc(testDouble);
//	std::cout<<"testFunc returned = "<<r<<std::endl;
	
	//--------------------------------------------------------------------------
	//------------------------------- 	START REAL CALCULATE STEPS -------------
	//--------------------------------------------------------------------------
	//calc those that are static for each epoch
	atot =pow(((params[7]*params[7]*secPerYear*secPerYear*Grav*KGperMsun*(params[0]+params[1]))/(4.0*pi*pi)),(1.0/3.0));
	params[10]=atot/MperAU;
	if (dataRealAry[5]!=0){
		K = ((2.0*pi*((atot)/(1.0+(params[1]/params[0])))*sin(params[8]*(pi/180.0)))/(params[7]*secPerYear*pow((1.0-params[4]*params[4]),(1.0/2.0))));
		params[12]=K;
	}
	//start loop over each epoch of data
	for (int i=0;i<dataRealAry_nx; i++){
		//------------------
		//Calculate TA and E
		//------------------
		M = (2.0*pi*(dataRealAry[0+i*dataRealAry_ny]-2.0*params[5]+params[6])/(P*daysPerYear));
		if (M>2.0*pi)
			M -= 2.0*pi;
		if (M<0)
			M += 2.0*pi;
		Eprime = M+params[4]*sin(M)+((params[4]*params[4])/(2.0*M))*sin(2.0*M);
		NewtonCount = 0;
		while ( (fabs(E-Eprime)>1.0e-10)&&(NewtonCount<50) ){
			E = Eprime;
			Eprime = E-((E-params[4]*sin(E)-M)/(1.0-params[4]*cos(E)));
			NewtonCount +=1;
		}
		//double check it satisfies the original equation
		if (fabs((E-params[4]*sin(E))-M)>1.0e-5)
			std::cout<<"PROBLEM!! resulting E from Newton's loop isn't within error limit!!!"<<std::endl;
		thetaPrime = acos((cos(E)-params[4])/(1.0-params[4]*cos(E)));
		if (E>pi)
			thetaPrime = 2.0*pi-thetaPrime;
		theta = thetaPrime;
		//--------------------------
		//Calculate RV
		//--------------------------
		if (dataRealAry[5]!=0){
			dataModelAry[2+i*dataRealAry_nx]=K*(cos(theta+params[9]*(pi/180.0))+params[4]*cos(params[9]*(pi/180.0)));
		}
		else
			dataModelAry[2+i*dataRealAry_nx]=0.0;
		std::cout<<"RV = "<<dataModelAry[2+i*dataRealAry_nx] <<std::endl;
		//--------------------------
		//Calculate x,y
		//--------------------------
		if ((dataRealAry[1]!=0)&&(dataRealAry[3]!=0)){
			// calculate all the Thiele-Innes constants in ["]
			A = ((atot/MperAU)/params[2])*(cos(params[3]*(pi/180.0))*cos(params[9]*(pi/180.0))-sin(params[3]*(pi/180.0))*sin(params[9]*(pi/180.0))*cos(params[8]*(pi/180.0)));
			B = ((atot/MperAU)/params[2])*(sin(params[3]*(pi/180.0))*cos(params[9]*(pi/180.0))+cos(params[3]*(pi/180.0))*sin(params[9]*(pi/180.0))*cos(params[8]*(pi/180.0)));
			F = ((atot/MperAU)/params[2])*(-cos(params[3]*(pi/180.0))*sin(params[9]*(pi/180.0))-sin(params[3]*(pi/180.0))*cos(params[9]*(pi/180.0))*cos(params[8]*(pi/180.0)));
			G = ((atot/MperAU)/params[2])*(-sin(params[3]*(pi/180.0))*sin(params[9]*(pi/180.0))+cos(params[3]*(pi/180.0))*cos(params[9]*(pi/180.0))*cos(params[8]*(pi/180.0)));
			// The coordinates of the unit orbital ellipse in the true plane (Binnendijk)
			X = cos(E)-params[4];
			Y = sqrt(1.0-params[4]*params[4])*sin(E);

			// Calculate the predicted x&y in ["]
			dataModelAry[0+i*dataRealAry_nx] = A*X +F*Y;
			dataModelAry[1+i*dataRealAry_nx] = B*X +G*Y;
		}
		else
			dataModelAry[0+i*dataRealAry_nx] = dataModelAry[1+i*dataRealAry_nx] = 0.0;
		std::cout<<"x = "<<dataModelAry[0+i*dataRealAry_nx] <<std::endl;
		std::cout<<"y = "<<dataModelAry[1+i*dataRealAry_nx] <<std::endl;
	}


}


















