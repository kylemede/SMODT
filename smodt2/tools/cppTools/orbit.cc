#include "orbit.h"

double testFunc(double t){
    std::cout<<"\nInside testFunc"<<std::endl;
    t = t*300;
    return t;
};

void Orbit::anomalyCalc(double ecc, double T, double Tc,double P, double epoch){
	//------------------
	//Calculate TA and E
	//Remember that in RV, there is a shift due to the Tc!=T, that doesn't exist in DI.
	//------------------
	//std::cout<<"\necc = "<<ecc<<", T = "<<T<<", Tc = "<<Tc<<", P = "<<P<<", epoch = "<<epoch<<std::endl;
	thetaRV=0;
	EDI=0;
	//for RV
	M = (2.0*pi*(epoch-2.0*T+Tc))/(P*daysPerYear);
	M -= (int)(M/(2.0*pi))*(2.0*pi);//shift into [-360,360]
	if ((M!=0)and(M!=(2.0*pi))){
		Eprime = M+ecc*sin(M)+((ecc*ecc)/(2.0*M))*sin(2.0*M);
		newtonCount = 0;
		while ( (fabs(E-Eprime)>1.0e-10)&&(newtonCount<50) ){
			E = Eprime;
			Eprime = E-((E-ecc*sin(E)-M)/(1.0-ecc*cos(E)));
			newtonCount +=1;
		}
		//double check it satisfies the original equation
		if (fabs((E-ecc*sin(E))-M)>1.0e-5){
			std::cout<<"PROBLEM!! resulting E from Newton's loop isn't within error limit!!!"<<std::endl;
			if (true){
				std::cout<<"M = "<<M <<std::endl;
				std::cout<<"e = "<<ecc<<std::endl;
				std::cout<<"T = "<<T<<std::endl;
				std::cout<<"Tc = "<<Tc<<std::endl;
				std::cout<<"P = "<<P<<std::endl;
				std::cout<<"Eprime = "<<Eprime <<"\n" <<std::endl;
			}
		}
		//std::cout<<"E RV = "<<E<<std::endl;
		thetaPrime = acos((cos(E)-ecc)/(1.0-ecc*cos(E)));
		if (E>pi)
			thetaPrime = 2.0*pi-thetaPrime;
		thetaRV = thetaPrime;
	}
	//for DI
	if (T!=Tc){
		M = (2.0*pi*(epoch-T))/(P*daysPerYear);
		M -= (int)(M/(2.0*pi))*(2.0*pi);//shift into [-360,360]
		if ((M!=0)and(M!=(2.0*pi))){
			Eprime = M+ecc*sin(M)+((ecc*ecc)/(2.0*M))*sin(2.0*M);
			newtonCount = 0;
			while ( (fabs(E-Eprime)>1.0e-10)&&(newtonCount<50) ){
				E = Eprime;
				Eprime = E-((E-ecc*sin(E)-M)/(1.0-ecc*cos(E)));
				newtonCount +=1;
			}
			//double check it satisfies the original equation
			if (fabs((E-ecc*sin(E))-M)>1.0e-5){
				std::cout<<"PROBLEM!! resulting E from Newton's loop isn't within error limit!!!"<<std::endl;
				if (true){
					std::cout<<"M = "<<M <<std::endl;
					std::cout<<"e = "<<ecc<<std::endl;
					std::cout<<"T = "<<T<<std::endl;
					std::cout<<"Tc = "<<Tc<<std::endl;
					std::cout<<"P = "<<P<<std::endl;
					std::cout<<"Eprime = "<<Eprime <<"\n" <<std::endl;
				}
			}
		}
	}
	EDI = E;
	//std::cout<<"\nin anomaly calc: EDI = "<<EDI<<", thetaRV = "<<thetaRV<<std::endl;//$$$$$$$$$$$$$$$$$$
};

void Orbit::loadStaticVars(double omegaoffsetDI,double omegaoffsetRV){
	omegaOffsetDI = omegaoffsetDI;
	omegaOffsetRV = omegaoffsetRV;
};

void Orbit::loadRealData(double *xx, int xx_nx, int xx_ny){
	if (false)
		std::cout<<"\nInside loadData function"<<std::endl;
    dataRealAry = xx;
    dataRealAry_nx = xx_nx;
    dataRealAry_ny = xx_ny;
    if (false)
    	std::cout<<"data loaded!"<<std::endl;
};
void Orbit::loadConstants(double Grav_in,double pi_in,double KGperMsun_in, double daysPerYear_in,double secPerYear_in,double MperAU_in){
	Grav = Grav_in;
	pi = pi_in;
	KGperMsun = KGperMsun_in;
	daysPerYear = daysPerYear_in;
	secPerYear = secPerYear_in;
	MperAU = MperAU_in;
	if (false)
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
	bool verbose=false;
	if (verbose)
		std::cout<<"\nInside Orbit calculator function"<<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
	if (verbose){//$$$$$$$$$$$$$$$$$$$$$$$$$
		std::cout<<"real data inside Orbit c++ code:"<<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
		for (int i=0; i<dataRealAry_nx; i++){//$$$$$$$$$$$$$$$$$$$$$$$$$
			std::cout<<"[";//$$$$$$$$$$$$$$$$$$$$$$$$$
			for (int j=0;j<dataRealAry_ny;j++){//$$$$$$$$$$$$$$$$$$$$$$$$$
				std::cout<<dataRealAry[j+i*dataRealAry_ny]<<", ";//$$$$$$$$$$$$$$$$$$$$$$$$$
			}//$$$$$$$$$$$$$$$$$$$$$$$$$
			std::cout<<"]"<<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
		}//$$$$$$$$$$$$$$$$$$$$$$$$$
	}

	//--------------------------------------------------------------------------
	//------------------------------- 	START REAL CALCULATE STEPS -------------
	//--------------------------------------------------------------------------
	//calc those that are static for each epoch
	atot =pow(((params[7]*params[7]*secPerYear*secPerYear*Grav*KGperMsun*(params[0]+params[1]))/(4.0*pi*pi)),(1.0/3.0));
	params[10]=atot/MperAU;
	if ((dataRealAry[5]!=0)&&(params[12]==0)){
		K = ((2.0*pi*((atot)/(1.0+(params[0]/params[1])))*sin(params[8]*(pi/180.0)))/(params[7]*secPerYear*pow((1.0-params[4]*params[4]),(1.0/2.0))));
		params[12]=K;
	}
	omegaDI = params[9]+omegaOffsetDI;
	omegaRV = params[9]+omegaOffsetRV;
	//Calculate Tc <-> T if needed
	if (params[5]!=params[6]){
		//if T=Tc already, do nothing.
		if ((params[4]==0)||(omegaRV==90.0)){
			//Circular, so just set equal.
			if (params[6]==0)
				params[5]=params[6];
			else
				params[6]=params[5];
		}
		else{
			ta = pi/2.0 - omegaRV*(pi/180.0);
			halfE = atan2(sqrt(1.0-params[4])*sin(ta/2.0),sqrt(1.0+params[4])*cos(ta/2.0));
			mTTc = 2.0*halfE-params[4]*sin(2.0*halfE);
			deltaT = (mTTc*params[7]*daysPerYear)/(2.0*pi);
			if (params[6]==0)
				params[6] = params[5]+deltaT;
			else
				params[5] = params[6]-deltaT;
		}
	}
	//start loop over each epoch of data
	for (int i=0;i<dataModelAry_nx; i++){
		if (verbose)//$$$$$$$$$$$$$$$$$$$$$$$$$
			std::cout<<"\nepoch "<<dataRealAry[0+i*dataRealAry_ny]<<":"<<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
		//Calculate true anomaly for RV and eccentric anomaly for DI
		anomalyCalc(params[4],params[5],params[6],params[7],dataRealAry[0+i*dataRealAry_ny]);
		//std::cout<<"in calc: EDI = "<<EDI<<", thetaRV = "<<thetaRV<<std::endl;//$$$$$$$$$$$$$$$$$$
		//--------------------------
		//Calculate RV
		//--------------------------
		if (dataRealAry[5+i*dataRealAry_ny]!=0){
			dataModelAry[2+i*dataModelAry_ny]=K*(cos(thetaRV+omegaRV*(pi/180.0))+params[4]*cos(omegaRV*(pi/180.0)));
			if (false){
				std::cout<<"theta deg V2.0 = "<<(thetaRV*(180.0/pi))<<std::endl;
				std::cout<<"cos(theta+omegaRV*(pi/180.0)) = "<<cos(thetaRV+omegaRV*(pi/180.0))<<std::endl;
				std::cout<<"cos(omegaRV*(pi/180.0)) = "<< cos(omegaRV*(pi/180.0))<<std::endl;
				std::cout<<"params[4]*cos(omegaRV*(pi/180.0)) = "<< params[4]*cos(omegaRV*(pi/180.0))<<std::endl;
				std::cout<<"dataModelAry[2+i*dataModelAry_ny] = "<<dataModelAry[2+i*dataModelAry_ny] <<std::endl;
			}
		}
		else{
			//std::cout<<"RV real = "<< dataRealAry[5+i*dataRealAry_ny]<<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
			dataModelAry[2+i*dataModelAry_ny]=0.0;
		}
		if (verbose){//$$$$$$$$$$$$$$$$$$$$$$$$$
			std::cout<<"RV = "<<dataModelAry[2+i*dataModelAry_ny] <<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
		}//$$$$$$$$$$$$$$$$$$$$$$$$$
		//--------------------------
		//Calculate x,y
		//--------------------------
		if ((dataRealAry[1+i*dataRealAry_ny]!=0)&&(dataRealAry[3+i*dataRealAry_ny]!=0)){
			// calculate all the Thiele-Innes constants in ["]
			A = ((atot/MperAU)/params[2])*(cos(params[3]*(pi/180.0))*cos(omegaDI*(pi/180.0))-sin(params[3]*(pi/180.0))*sin(omegaDI*(pi/180.0))*cos(params[8]*(pi/180.0)));
			B = ((atot/MperAU)/params[2])*(sin(params[3]*(pi/180.0))*cos(omegaDI*(pi/180.0))+cos(params[3]*(pi/180.0))*sin(omegaDI*(pi/180.0))*cos(params[8]*(pi/180.0)));
			F = ((atot/MperAU)/params[2])*(-cos(params[3]*(pi/180.0))*sin(omegaDI*(pi/180.0))-sin(params[3]*(pi/180.0))*cos(omegaDI*(pi/180.0))*cos(params[8]*(pi/180.0)));
			G = ((atot/MperAU)/params[2])*(-sin(params[3]*(pi/180.0))*sin(omegaDI*(pi/180.0))+cos(params[3]*(pi/180.0))*cos(omegaDI*(pi/180.0))*cos(params[8]*(pi/180.0)));
			// The coordinates of the unit orbital ellipse in the true plane (Binnendijk)
			X = cos(EDI)-params[4];
			Y = sqrt(1.0-params[4]*params[4])*sin(EDI);
			// Calculate the predicted x&y in ["]
			//KEY NOTE: x_TH-I = y_plot = North
			//          y_TH-I = x_plot = East
			//THUS, store x=y_TH-I, y=x_TH-I !!!
			dataModelAry[0+i*dataModelAry_ny] = B*X +G*Y;
			dataModelAry[1+i*dataModelAry_ny] = A*X +F*Y;
		}
		else{
			//std::cout<<"x real = "<< dataRealAry[1+i*dataRealAry_nx]<<", y real = "<< dataRealAry[3+i*dataRealAry_nx]<<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
			dataModelAry[0+i*dataModelAry_ny] = dataModelAry[1+i*dataModelAry_ny] = 0.0;
		}
		if (verbose){//$$$$$$$$$$$$$$$$$$$$$$$$$
			std::cout<<"x = "<<dataModelAry[0+i*dataModelAry_ny] <<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
			std::cout<<"y = "<<dataModelAry[1+i*dataModelAry_ny] <<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
		}//$$$$$$$$$$$$$$$$$$$$$$$$$
	}
	if (verbose){//$$$$$$$$$$$$$$$$$$$$$$$$$
		std::cout<<"\nModel data (after calculating all model epochs) inside Orbit c++ code:"<<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
		for (int i=0; i<dataModelAry_nx; i++){//$$$$$$$$$$$$$$$$$$$$$$$$$
			std::cout<<"[";//$$$$$$$$$$$$$$$$$$$$$$$$$
			for (int j=0;j<dataModelAry_ny;j++){//$$$$$$$$$$$$$$$$$$$$$$$$$
				std::cout<<dataModelAry[j+i*dataModelAry_ny]<<", ";//$$$$$$$$$$$$$$$$$$$$$$$$$
			}//$$$$$$$$$$$$$$$$$$$$$$$$$
			std::cout<<"]"<<std::endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
		}//$$$$$$$$$$$$$$$$$$$$$$$$$
	}


}


















