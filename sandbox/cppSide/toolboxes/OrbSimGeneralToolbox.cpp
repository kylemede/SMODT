#include <iostream>
#include <iomanip>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <sstream>
#include <string>
#include <fstream>
#include <functional> //NEW
#include <algorithm> //NEW
#include <numeric> //NEW
#include "orbToolboxes.h"

using namespace std;

outputDataType odtStart(outputDataType ODT, int numTotalSamples)
{
	//repeat for all vectors in ODT
	ODT.a_totals.reserve(numTotalSamples);
	ODT.longAN_degs.reserve(numTotalSamples);
	ODT.es.reserve(numTotalSamples);
	ODT.Ts.reserve(numTotalSamples);
	ODT.Tcs.reserve(numTotalSamples);
	ODT.periods.reserve(numTotalSamples);
	ODT.inclination_degs.reserve(numTotalSamples);
	ODT.argPeri_degs.reserve(numTotalSamples);
	//RVoffsets
	ODT.Ks.reserve(numTotalSamples);
	ODT.chiSquareds.reserve(numTotalSamples);
	ODT.timesBeenHeres.reserve(numTotalSamples);
	return ODT;
}

outputDataType odtFinish(outputDataType ODT)
{
	//repeat for all vectors in ODT
//	ODT.a_totals.shrink_to_fit();
//	ODT.longAN_degs.shrink_to_fit();
//	ODT.es.shrink_to_fit();
//	ODT.Ts.shrink_to_fit();
//	ODT.periods.shrink_to_fit();
//	ODT.inclination_degs.shrink_to_fit();
//	ODT.argPeri_degs.shrink_to_fit();
//	ODT.Ks.shrink_to_fit();
//	ODT.chiSquareds.shrink_to_fit();
//	ODT.timesBeenHeres.shrink_to_fit();

	return ODT;
}

eccArgPeri2ToTcType eccArgPeri2ToTcCalc(eccArgPeri2ToTcType EATT)
{
	bool verbose = false;
	bool backHalf = false;

	if (EATT.e==0)
	{
		if (EATT.Tc!=0)
			EATT.To = EATT.Tc;
		else
			EATT.Tc = EATT.To;
	}
	else
	{
		double To_IN = EATT.To;
		double Tc_IN = EATT.Tc;

		double TA_s_deg = 90.0 - EATT.argPeri_deg;
		if (backHalf)
			TA_s_deg = TA_s_deg+180.0;

		if (TA_s_deg<0.0)
			TA_s_deg = TA_s_deg+360.0;
		double TA_s_rad = TA_s_deg*(PI/180.0);
		double top = sqrt(1.0-EATT.e)*sin(TA_s_rad/2.0);
		double btm = sqrt(1.0+EATT.e)*cos(TA_s_rad/2.0);
		double ATAN_rad ;
	//	ATAN_rad = atanTopBtm(top,btm);
	//	//fix angles calculated to make sense
	//	if (ATAN_rad>(2.0*PI))
	//		ATAN_rad = ATAN_rad-2.0*PI;
	//	if ((ATAN_rad>0)&&(TA_s_rad<0))
	//		ATAN_rad = ATAN_rad-2.0*PI;

//		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//		//$$$$$$$$$$$$$$$$$$$$$$$$$$ TEST $$$$$$$$$$$$$$$$$$$$$$$$$
//		if (true)
//		{
//			vector<double> argPeri_degs;
//			argPeri_degs.push_back(1.0);
//			argPeri_degs.push_back(89.0);
//			argPeri_degs.push_back(91.0);
//			argPeri_degs.push_back(179.0);
//			argPeri_degs.push_back(181.0);
//			argPeri_degs.push_back(359.0);
//			argPeri_degs.push_back(-1.0);
//			double e = 0.1;
//			cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
//			cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
//			for ( int i=0;i<argPeri_degs.size();i++)
//			{
//
//				double TA_s_deg_TEST1 = 90.0 - argPeri_degs[i];
//				if (backHalf)
//					TA_s_deg_TEST1 = TA_s_deg_TEST1+180.0;
//
//				if (TA_s_deg_TEST1<0.0)
//					TA_s_deg_TEST1 = TA_s_deg_TEST1+360.0;
//				double TA_s_rad_TEST1 = TA_s_deg_TEST1*(PI/180.0);
//				double top_TEST = sqrt(1.0-e)*sin(TA_s_rad_TEST1/2.0);
//				double btm_TEST = sqrt(1.0+e)*cos(TA_s_rad_TEST1/2.0);
//				double ATAN_rad_TEST1 ;
//
//				ATAN_rad_TEST1 = atanTopBtm(top_TEST,btm_TEST);
//				//fix angles calculated to make sense
//				if (ATAN_rad_TEST1>(2.0*PI))
//					ATAN_rad_TEST1 = ATAN_rad_TEST1-2.0*PI;
//				if ((ATAN_rad_TEST1>0)&&(TA_s_rad_TEST1<0))
//					ATAN_rad_TEST1 = ATAN_rad_TEST1-2.0*PI;
//				double E_s_rad_TEST1 = ATAN_rad_TEST1*2.0;
//
//				double 	ATAN_rad_TEST2;
//				ATAN_rad_TEST2 = atan2(top_TEST,btm_TEST);
//				double E_s_rad_TEST2 = ATAN_rad_TEST2*2.0;
//
//				cout<<"\n\nargPeri_degs[i] = "<< argPeri_degs[i] <<endl;
//				cout<<"TA_s_deg_TEST1 = "<< TA_s_deg_TEST1 <<endl;
//				cout<<"top_TEST = "<< top_TEST <<endl;
//				cout<<"btm_TEST = "<<btm_TEST  <<endl;
//				cout<<"ATAN_rad_TEST1 = "<<ATAN_rad_TEST1*(180.0/PI)  <<endl;
//				cout<<"E_s_rad_TEST1 = "<<E_s_rad_TEST1*(180.0/PI)  <<endl;
//				cout<<"\nATAN_rad_TEST2 = "<< ATAN_rad_TEST2*(180.0/PI)<<endl;
//				cout<<"E_s_rad_TEST2 = "<< E_s_rad_TEST2*(180.0/PI) <<endl;
//				cout<<"E_s_rad_TEST1-E_s_rad_TEST2 = "<<(E_s_rad_TEST1-E_s_rad_TEST2)*(180.0/PI)  <<endl;
//				//cout<<" = "<<  <<endl;
//			}
//			cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
//			cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
//		}


		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		ATAN_rad = atan2(top,btm);
		double E_s_rad = ATAN_rad*2.0;
		double M_s_rad = E_s_rad-EATT.e*sin(E_s_rad);
		double delta_t = (M_s_rad*EATT.period*365.242)/(2.0*PI);
		double delta_t_orig = delta_t;

		// in the case that the situation is flipped during the 'backHalf' option
		if (backHalf)
		{
			// check if argPeri is inside the back part of orbit where To>Tc
			if (EATT.argPeri_deg>270)
				delta_t = delta_t-EATT.period*365.242;
			else if (EATT.argPeri_deg==270)
				delta_t = 0.0;
		}
		// else apply standard To calc for planet's orbit
		else
		{
			// check if argPeri is inside the back part of orbit where To>Tc
			if (EATT.argPeri_deg>90)
				delta_t = delta_t-EATT.period*365.242;
			else if (EATT.argPeri_deg==90)
				delta_t = 0.0;
		}

		if (EATT.Tc!=0)
			EATT.To = EATT.Tc-delta_t;
		else
			EATT.Tc = EATT.To+delta_t;
		if (verbose)
		{
			cout<<"\ninputs to eccArgPeri2ToTcCalc:"<<endl;
			cout<<"EATT.e = "<< EATT.e<<endl;
			cout<<"EATT.argPeri_deg = "<< EATT.argPeri_deg<<endl;
			cout<<"EATT.period = "<<EATT.period <<endl;
			cout<<"To_IN = "<< To_IN<<endl;
			cout<<"Tc_IN = "<< Tc_IN<<endl;
			cout<<"Calculated:"<<endl;
			cout<<"TA_s_deg = "<<TA_s_deg <<endl;
			cout<<"ATAN_rad*(180.0/PI) = "<< ATAN_rad*(180.0/PI)<<endl;
			cout<<"E_s_rad*(180.0/PI) = "<< E_s_rad*(180.0/PI)<<endl;
			cout<<"M_s_rad*(180.0/PI) = "<<M_s_rad*(180.0/PI) <<endl;
			cout<<"delta_t_orig= "<<delta_t_orig <<endl;
			cout<<"delta_t = "<< delta_t<<endl;
			cout<<"outputs:"<<endl;
			cout<<"EATT.To = "<< EATT.To<<endl;
			cout<<"EATT.Tc = "<< EATT.Tc<<endl;
		}
	}

	return EATT;
}

double atanTopBtm(double top, double btm)
{
	double TAN = top/btm;
	double ATAN_rad = atan(top/btm);

	int quad = 0;
	if ((top>=0.0)and(TAN>=0.0))
		quad = 1;
	if ((top>=0.0)and(TAN<=0.0))
		quad = 2;
	if ((top<=0.0)and(TAN>=0.0))
		quad = 3;
	if ((top<=0.0)and(TAN<=0.0))
		quad = 4;

	if ((quad==2)||(quad==3))
		ATAN_rad = ATAN_rad+PI;
	else if (quad==4)
		ATAN_rad = ATAN_rad+2.0*PI;
	if (ATAN_rad>(2.0*PI))
		ATAN_rad = ATAN_rad-2.0*PI;

	return ATAN_rad;
}

eccArgPeriCalcType eccArgPeriCalc(eccArgPeriCalcType EACT)
{
	bool verbose=false;
	EACT.e;
	EACT.argPeri_deg;
	EACT.sqrtEsinArgPeri;
	EACT.sqrtEcosArgPeri;

	if (verbose)
	{
		cout<<"-----------------------------------------------------"<<endl;
		cout<< "EACT.sqrtEsinArgPeri = "<< EACT.sqrtEsinArgPeri<<endl;
		cout<< "EACT.sqrtEcosArgPeri = "<< EACT.sqrtEcosArgPeri<<endl;
	}

	double ecc = (EACT.sqrtEsinArgPeri*EACT.sqrtEsinArgPeri)+(EACT.sqrtEcosArgPeri*EACT.sqrtEcosArgPeri);
	if (verbose)
		cout<<"\necc = "<<ecc<<endl;

	EACT.e = ecc;

	double ATAN_rad ;
	ATAN_rad = atanTopBtm(EACT.sqrtEsinArgPeri,EACT.sqrtEcosArgPeri);

	if (verbose)
		cout<<"Before conditionals applied, atan_deg = "<<(ATAN_rad * 180 / PI)<<endl;
//
//	int quad = 0;
////	if ((EACT.sqrtEsinArgPeri>=0.0)and(EACT.sqrtEcosArgPeri>=0.0))
////		quad = 1;
////	if ((EACT.sqrtEsinArgPeri>=0.0)and(EACT.sqrtEcosArgPeri<=0.0))
////		quad = 2;
////	if ((EACT.sqrtEsinArgPeri<=0.0)and(EACT.sqrtEcosArgPeri>=0.0))
////		quad = 3;
////	if ((EACT.sqrtEsinArgPeri<=0.0)and(EACT.sqrtEcosArgPeri<=0.0))
////		quad = 4;
//
//	if ((EACT.sqrtEsinArgPeri>=0.0)and(TAN>=0.0))
//		quad = 1;
//	if ((EACT.sqrtEsinArgPeri>=0.0)and(TAN<=0.0))
//		quad = 2;
//	if ((EACT.sqrtEsinArgPeri<=0.0)and(TAN>=0.0))
//		quad = 3;
//	if ((EACT.sqrtEsinArgPeri<=0.0)and(TAN<=0.0))
//		quad = 4;
//
//	if ((quad==2)||(quad==3))
//		ATAN = ATAN+PI;
//	else if (quad==4)
//		ATAN = ATAN+2.0*PI;
//	if (ATAN>(2.0*PI))
//		ATAN = ATAN-2.0*PI;

	if (verbose)
	{
		cout<<"After conditionals applied, atan_deg = "<<(ATAN_rad * 180 / PI)<<endl;
		cout<<"-----------------------------------------------------"<<endl;
	}

	double atan_deg=ATAN_rad * 180 / PI;
	EACT.argPeri_deg = atan_deg;

	return EACT;
}

TAcalcReturnType TAcalculator(TAcalcInputType TACIT)
{
	TAcalcReturnType TACRT;
	bool verboseOrig = TACIT.verbose;
	//TACIT.verbose = true;

	if (TACIT.verbose==true)
		cout<<"******************** Starting TAcalculator *******************"<<endl;
	//Calculate the mean motion and save to output structure
	double n = (2.0*PI)/TACIT.period; //period in years
	double period_days = TACIT.period*365.242;
	double timeDiff_days = (TACIT.t-TACIT.T)-int((TACIT.t-TACIT.T)/period_days)*period_days;
	if (TACIT.verbose==true)
	{
		cout<<"TACIT.t-TACIT.T = "<<(TACIT.t-TACIT.T) <<", int((TACIT.t-TACIT.T)/period_days)*period_days = "<< (int((TACIT.t-TACIT.T)/period_days)*period_days)<<endl;
		cout<<"timeDiff_days = "<<timeDiff_days<<endl;
	}
	if (timeDiff_days<0.0)
		timeDiff_days = timeDiff_days + period_days;
	//Set phase to zero and update if Tc is provided
	double phase = 0.0;
	if (TACIT.Tc!=0)
	{
		double phaseDiff_days = (TACIT.Tc-TACIT.T)-int((TACIT.Tc-TACIT.T)/period_days)*period_days;
		if (TACIT.verbose==true)
			cout<<"(TACIT.Tc-TACIT.T) = "<<(TACIT.Tc-TACIT.T) <<", int((TACIT.t-TACIT.T)/period_days)*period_days = "<< (int((TACIT.Tc-TACIT.T)/period_days)*period_days)<<endl;
		if (TACIT.T>TACIT.Tc)
			phaseDiff_days = phaseDiff_days + period_days;
		phase = phaseDiff_days/period_days;
	}
	if (TACIT.verbose==true)
	{
		cout<<"unitless phase calculated to be "<<phase<<", using Tc = "<<TACIT.Tc<<" and To = "<<TACIT.T<<endl;
		cout<<"epoch = " <<TACIT.t<<", timeDiff_days = "<<timeDiff_days<<", period_days = "<<period_days<<endl;
	}

	//Calculate the Mean Anomaly and add to output structure
	double M = n*((timeDiff_days/365.242)+phase);//+(phase*2.0*PI);//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	double numCirclesBiggerD = fabs((M/(2.0*PI)));
	int numCirclesBiggerI = int(numCirclesBiggerD);
	if (numCirclesBiggerI<1)
		numCirclesBiggerI = 1;
	double M_out;
	if (M<0)
	{
		M_out = (numCirclesBiggerI+1)*2.0*PI+M;
		if (TACIT.verbose==true)
		{
			cout<<"M negative, so updated from, "<<M<<" to "<<M_out<<", numCirclesBiggerI = "<<numCirclesBiggerI<<", numCirclesBiggerD = "<<numCirclesBiggerD<<endl;
		}
	}
	else if (M>(2.0*PI))
	{
		M_out = M-(numCirclesBiggerI)*2.0*PI;
		if (TACIT.verbose==true)
		{
			cout<<"M bigger than one circle, so updated from, "<<M<<" to "<<M_out<<", numCirclesBiggerI = "<<numCirclesBiggerI<<", numCirclesBiggerD = "<<numCirclesBiggerD<<endl;
		}
	}
	else
		M_out = M;
	M = M_out;
	//convert to degrees
	double M_deg = M*(180.0/PI);
	if (TACIT.verbose==true)
			cout<<"n = "<<n << ", M_deg = " << M_deg<<endl;

	//Set temp input values for E_latest and E_last for initial values of Newton's loops below
	double E_last = 2.0*PI;
	double E_best_intial_guess = M+(TACIT.e*sin(M))+((TACIT.e*TACIT.e)/(2.0*M))*sin(2.0*M);
	//cout<<"e = "<<TACIT.e <<", M = "<< M<<", E_best_intial_guess = "<< E_best_intial_guess <<endl;
	double E_latest = E_best_intial_guess;
	int count = 0;

	//Perform Newton's loop to find the solution to Kepler's equation
	double M_last;
	while ( (fabs(E_last-E_latest)>1.0e-10)&&(count<50) )
	{
		E_last = E_latest;
		M_last = E_last-TACIT.e*sin(E_last);
		E_latest = E_last-((M_last-M)/(1.0-TACIT.e*cos(E_last)));
		count +=1;
	}
	if (E_latest<0.0)
	{
		if (true) //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			cout<<"E_latest found to be negative, = "<<E_latest<<endl;
		E_latest = 2.0*PI-E_latest;
		if (true)//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			cout<<"E_latest updated to = "<<E_latest<<endl;
	}
	//Double check solution satisfies original equation
	//Then save the found Eccentric Anomaly if satisfied
	double Mnewton = (180.0/PI)*(E_latest-TACIT.e*sin(E_latest));
	if ( fabs(M_deg-Mnewton)>1.0e-5 )
	{
		cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
		cout << "PROBLEM: This resultant E does not satisfy the original equation, Newton's method FAILED !!!"<<endl;
		cout << "E_latest found "<<E_latest*(180.0/PI) <<" causes Mnewton "<<Mnewton <<" != M_deg " <<M_deg << " using best initial guess for E = "<< E_best_intial_guess*(180.0/PI)<<endl;
		cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
	}
	else
	{
		if ( TACIT.verbose==true )
			cout<<"This resultant E solves the original equation, Newton's Method worked :-)"<<endl;
		//Save solution to output structure
		TACRT.E_deg = E_latest*(180.0/PI);
	}
	if ( TACIT.verbose==true )
		cout<<"Eccentric Anomaly [deg] = "  << TACRT.E_deg <<endl;

	//Calculate the True Anomaly from the Eccentric Anomaly
	double TA_rad = acos((cos(E_latest)-TACIT.e)/(1.0-TACIT.e*cos(E_latest)));
	double TA_rad_out;
	if (( E_latest>PI ) || (E_latest<0))
		TA_rad_out = 2.0*PI - TA_rad;
	else
		TA_rad_out = TA_rad;
	//Save True Anomaly to output structure
	TACRT.TA_deg = TA_rad_out*(180.0/PI);

	if (false) //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	{
		cout<<"\nTACRT.E_deg = "<< TACRT.E_deg<<endl;
		cout<<"TA_rad*(180.0/PI) = "<<TA_rad*(180.0/PI) <<endl;
		cout<<"TA_rad_out*(180.0/PI) = "<<TA_rad_out*(180.0/PI) <<endl;

	}   //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	if ( TACIT.verbose==true )
		cout<< "True Anomaly [deg] = "<<TACRT.TA_deg<<endl;

	if (TACIT.verbose==true)
			cout<<"******************** Done TAcalculator *******************"<<endl;
	TACIT.verbose = verboseOrig;
	return TACRT;
}

double chiSquaredCalc(double real, double inv_var, double model)
{
	double difference = real - model;

	double chi_squared = (difference*difference)*inv_var;

	if ( (difference ==0)||(chi_squared==0) )
	{
		cout <<"!!!!! difference or chi_squared in chiSquaredCalc = 0"<<endl;
		cout<<"real = "<<real<<", model = "<<model<<endl;//$$$$$$
		cout<<"difference = "<<difference<<endl;//$$$$$
		cout<<"inv_var = "<<inv_var<<endl;//$$$$$$$$$$$$$$$$
		cout<<"chi_squared = "<<chi_squared<<endl;//$$$$$$$$$$$$$$$$$$$$$
	}
	return chi_squared;
}

void fileWriter(outputDataType ODT)
{
	// check inputs filename
	if (!(ODT.data_filename.find(".txt"))&&!(ODT.data_filename.find(".dat")))
		ODT.data_filename = ODT.data_filename+".dat";
	// make an output file for the INS and load it up.
	ofstream file;
	file.open(ODT.data_filename.c_str()) ;     //open the output file
	file << ODT.data_filename <<endl;
	file << "longAN [deg]   e [N/A]   To [julian date]  Tc [julian date]     period [yrs]  ";
	file <<"   inclination [deg]  argPeri [deg]  a_total [AU]    chiSquared    K [m/s]";
	for (int dataset=0;dataset<int(ODT.RVoffsets[0].size());++dataset)
		file<<"  RVorigin_"<<dataset<<" [m/s]  ";
	file <<" timesBeenHere"<<endl;

	for (int sample=0; sample<ODT.numSamplesAccepted; sample++)
	{
		file<< fixed<<std::setprecision(6)<<ODT.longAN_degs[sample];//
		file<< "      "<<ODT.es[sample];
		file<< "    " <<ODT.Ts[sample];
		file<< "    " <<ODT.Tcs[sample];
		file<< "   "<< std::setprecision(15)<<ODT.periods[sample];
		file<< "       "<< std::setprecision(5)<<ODT.inclination_degs[sample];
		file<< "          "<< ODT.argPeri_degs[sample];
		file<< "     "<< ODT.a_totals[sample];
		file<< "     "<< std::setprecision(8)<<ODT.chiSquareds[sample];
		file<< "     "<<ODT.Ks[sample];
		for (int set=0;set<int(ODT.RVoffsets[0].size());++set)
			file<<"      "<<ODT.RVoffsets[sample][set]<<"   ";

		file<< "             "<<ODT.timesBeenHeres[sample]<<endl;
	}//finished writing inputs file
	file.close () ;
	cout<<"\n***************************************************************"<<endl;
	cout<<"Output data file written to: "<<ODT.data_filename<<endl;
	cout<<"***************************************************************\n"<<endl;
}

void logFileWriter(string filename, string LOGlines)
{
	std::stringstream ss;

	int fnameLength = filename.length();
	if (fnameLength==0)
		cout<<"ERROR: length of filename for log was zero!"<<endl;
	else
	{
		string directory;
		const size_t last_slash_idx = filename.rfind('/');
		if (std::string::npos != last_slash_idx)
		{
		    directory = filename.substr(0, last_slash_idx);
		}
		//cout<<"directory was found to be: "<<directory<<endl;
		string logname;
		if (filename.find("chain_"))
		{
			//assuming if it has that str it is named by ***chain_#.dat or .txt

			string chainNumber1 = filename.substr(fnameLength-5,fnameLength-4);
			string chainNumber = chainNumber1.substr(0,chainNumber1.find("."));
			//cout<<"Chain number was found to be: "<<chainNumber<<endl;
			ss<<directory<<"/log-chain_"<<chainNumber<<".txt";
			logname = ss.str();
		}
		else
			logname = directory+"/log.txt";

		//cout<<"using logfilename :"<<logname<<endl;

		ofstream file;
		file.open(logname.c_str(),ios::app) ;     //open the output file
		file<<LOGlines;
		file.close () ;
		cout<<"\n###############################################################"<<endl;
		cout<<"Logfile written to: "<< logname<<endl;
		cout<<"###############################################################"<<endl;
	}

}
DItools DItoolsParamLoadUp(DIdataObj DIdo)
{
	bool verboseInternal = false;
	if (verboseInternal)
		cout<<"Inside DItoolsParamLoadUp"<<endl;
	DItools DIt;
	DIt.SAs_arcsec_observed = DIdo.SAs_arcsec_observed;
	DIt.SA_errors = DIdo.SA_errors ;
	DIt.PAs_deg_observed = DIdo.PAs_deg_observed ;
	DIt.PA_errors = DIdo.PA_errors ;
	DIt.epochs_DI = DIdo.epochs_DI ;
	DIt.Sys_Dist_PC = DIdo.Sys_Dist_PC;
	DIt.Mass1 = DIdo.Mass1;
	if (DIdo.star_Mass2!=0)
		DIt.Mass2  = DIdo.star_Mass2;
	else if (DIdo.planet_MsinI !=0)
		DIt.Mass2  = DIdo.planet_MsinI;

	if (verboseInternal)
		cout<<"Done DItoolsParamLoadUp and returning DIt loaded object"<<endl;
	return DIt;
}

VRcalcStarStar VRcalcStarStarLoadUp(RVdataObj RVdo)
{
	VRcalcStarStar VRCss;
	VRCss.T_s = RVdo.star_T;
	VRCss.Tc_s = RVdo.star_Tc;
	VRCss.argPeri_deg_s = RVdo.star_argPeri;
	VRCss.inclination_deg_s = RVdo.star_inc;
	VRCss.period_s = RVdo.star_P;
	VRCss.e_s = RVdo.star_e;
	VRCss.Mass1 = RVdo.Mass1;
	VRCss.Mass2_s = RVdo.star_Mass2;

	//convert semi majors to make a complete set
	semiMajorType SMTin;
	SMTin.Mass1 = RVdo.Mass1;
	SMTin.Mass2 =RVdo.star_Mass2;
	SMTin.a1 = 0;
	SMTin.a2 = 0;
	SMTin.a_total = 0;
	SMTin.period = RVdo.star_P;
	semiMajorType SMTout;

	SMTout = semiMajorConverter(SMTin);
	VRCss.a1 = SMTout.a1;
	VRCss.a_total = SMTout.a_total;

	return VRCss;
}

VRcalcStarPlanet VRcalcStarPlanetLoadUp(RVdataObj RVdo)
{
	//cout<<"\n### In VRcalcStarPlanetLoadUp!! ###"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	VRcalcStarPlanet VRCsp;
	VRCsp.K_p = RVdo.planet_K ;
	VRCsp.K_p_error = RVdo.planet_K_error;
	VRCsp.e_p = RVdo.planet_e ;
	VRCsp.T_p = RVdo.planet_T ;
	VRCsp.Tc_p = RVdo.planet_Tc;
	VRCsp.period_p = RVdo.planet_P ;
	VRCsp.Mass1 = RVdo.Mass1;
	VRCsp.Mass2sinI_p = RVdo.planet_MsinI;
	VRCsp.argPeri_deg_p = RVdo.planet_argPeri ;
	VRCsp.inclination_deg_p = RVdo.planet_inc;
	VRCsp.a1 = 0.0;
	if (VRCsp.K_p==0)
	{
		//convert semi majors to make a complete set
		semiMajorType SMTin;
		SMTin.Mass1 = RVdo.Mass1;
		SMTin.Mass2 = RVdo.planet_MsinI;
		SMTin.a1 = 0;
		SMTin.a2 = 0;
		SMTin.a_total = 0;
		SMTin.period = RVdo.planet_P;
		semiMajorType SMTout;

		SMTout = semiMajorConverter(SMTin);
		VRCsp.a1 = SMTout.a1;
		VRCsp.a_total = SMTout.a_total;
	}
	else
	{
		if (VRCsp.a_total>1e4)
			VRCsp.a_total=0.0;
	}

	return VRCsp;
}

semiMajorType semiMajorConverter(semiMajorType SMT)
{
	bool verboseInternal = false;
	semiMajorType SMTout;
	SMTout = SMT;
	if (verboseInternal)
	{
		cout<<"\n*In semiMajorConverter*"<<endl;
		cout<<"INPUTS were:"<<endl;
		cout<<"SMT.Mass1 = "<<SMT.Mass1 <<endl;
		cout<<"SMT.Mass2 = "<<SMT.Mass2 <<endl;
		cout<<"SMT.a1 = "<< SMT.a1<<endl;
		cout<<"SMT.a2 = "<<SMT.a2 <<endl;
		cout<<"SMT.a_total = "<< SMT.a_total<<endl;
		cout<<"SMT.period = "<< SMT.period<<endl;
	}
	if ((SMT.a1==0.0)&&(SMT.a2==0.0)&&(SMT.a_total==0.0))
	{
		if (verboseInternal)
			cout<<"calculating a_total from Period"<<endl;
		double temp5 = SMT.period*SMT.period*SecPerYear*SecPerYear*GravConst*KGperMsun*(SMT.Mass1+SMT.Mass2);
		double temp6 = 4.0*PI*PI;
		SMTout.a_total = pow((temp5/temp6),(1.0/3.0))/MperAU;
		if (verboseInternal)
			cout<<"a_total  = "<<SMTout.a_total<<endl;
	}

	if (SMTout.a1!=0.0)
	{
		if (verboseInternal)
			cout<<"a1 found to be non-zero"<<endl;
		SMTout.a2 = SMTout.a1*(SMTout.Mass1/SMTout.Mass2);
		SMTout.a_total = SMTout.a1+SMTout.a2;
	}
	else if (SMTout.a2!=0.0)
	{
		if (verboseInternal)
			cout<<"a2 found to be non-zero"<<endl;
		SMTout.a1 = SMTout.a2*(SMTout.Mass2/SMTout.Mass1);
		SMTout.a_total = SMTout.a1+SMTout.a2;
	}
	else if (SMTout.a_total!=0.0)
	{
		if (verboseInternal)
			cout<<"a_total found to be non-zero = "<<SMTout.a_total<<endl;
		SMTout.a2 = SMTout.a_total/(1.0+(SMTout.Mass2/SMTout.Mass1));
		SMTout.a1 = SMTout.a2*(SMTout.Mass2/SMTout.Mass1);
	}

	if (verboseInternal)
	{
		cout<<"\nOUTPUTS were:"<<endl;
		cout<<"SMTout.a1 = "<< SMTout.a1<<endl;
		cout<<"SMTout.a2 = "<<SMTout.a2 <<endl;
		cout<<"SMTout.a_total = "<< SMTout.a_total<<endl;
		cout<<"SMTout.period = "<< SMTout.period<<endl;
	}

	//check if any of the output values are crazy
	if ( SMTout.a_total>1e4)
		SMTout.a_total = 0;
	if ( SMTout.a2>1e4)
		SMTout.a2= 0;
	if (SMTout.a1 >1e4)
		SMTout.a1= 0;


	return SMTout;
}

double earliestEpochFinder(DIdataObj DIdo, RVdataObj RVdo)
{
	bool verbose = true;
	double earliestDIepoch = 1e9;
	double earliestRVepoch = 1e9;
	double earliestEpoch;
	if (DIdo.epochs_DI.size()>0)
	{
		for (int epoch=0;epoch<DIdo.epochs_DI.size();++epoch)
		{
			if (DIdo.epochs_DI[epoch]<earliestDIepoch)
				earliestDIepoch = DIdo.epochs_DI[epoch];
		}
	}
	if (RVdo.epochs_RV.size()>0)
	{
		for (int dataset=0;dataset<RVdo.epochs_RV.size();++dataset)
		{
			for (int epoch=0; epoch<RVdo.epochs_RV[dataset].size();++epoch)
			{
				if (RVdo.epochs_RV[dataset][epoch]<earliestRVepoch)
					earliestRVepoch = RVdo.epochs_RV[dataset][epoch];
			}
		}
	}
	if (verbose)
	{
		cout<<"earliestDIepoch = "<<earliestDIepoch <<endl;
		cout<<"earliestRVepoch = "<< earliestRVepoch<<endl;
	}
	if (earliestRVepoch<earliestDIepoch)
		earliestEpoch = earliestRVepoch;
	else if (earliestDIepoch<=earliestRVepoch)
		earliestEpoch = earliestDIepoch;

	return earliestEpoch;
}

string CorrelationLengthCalc(vector<double> data, string paramName, vector<int> timesBeenHere)
{
	std::stringstream ss;
	string outputStr;
	ss<<"\n****************************************************************"<<endl;
	//cout<<"\n****************************************************************"<<endl;
	if (data[0]==data.back())
		ss<<"\n"<<paramName<<" had a matching start and end value, so not calculating the correlation length."<<endl;
	else
	{
		vector<double> data2;
		data2 = fillOutDataVector(data,timesBeenHere);
		int corrLength = corrLengthJumpyCalc(data2);
		ss<<paramName<<" had Correlation length: "<<corrLength<<endl;
		ss<<data2.size()<<"/"<<corrLength<<" = "<<(data2.size()/double(corrLength))<<endl;
	}
	outputStr = ss.str();
	ss.clear();
	ss.str(std::string());
	return outputStr;
}

double standardDeviation(vector<double> v)
{
	//old version from info on the web
//	double sum = 0;
//	for (int i=0;i<v.size();i++)
//		sum = sum+v[i];
//	double mean = sum / v.size();
//	vector<double> diff(v.size());
//	std::transform(v.begin(), v.end(), diff.begin(),std::bind2nd(std::minus<double>(), mean));
//	double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
//	double stdev = sqrt(sq_sum / v.size());

	//new version using the Numerical Recipes variance calculator coded up in this toolbox
	double stdev;
	double var;
	var = varianceCalc(v, v.size());
	stdev = sqrt(var);
	return stdev;
}

double meanCalc(vector<double> v, int lastPoint)
{
	bool verbose = false;

	double sum;
	sum=0.0;
	//loop through all data points to get total value
	for (int j=0;j<(lastPoint+1);j++)
	{
		sum+=v[j];
//		if ((false)&&(lastPoint>((v.size()/10)*5)))
//			cout<<"sum in loop = "<<ave<<endl;
	}
	//convert total value into average value
	//ave/=(lastPoint+1);
	double ave;
	ave = sum/(double(lastPoint+1));
	if ((verbose)&&(lastPoint>((v.size()/10)*5)))
		cout<<"sum = "<<fixed<<std::setprecision(15)<<sum<<", = "<<double(lastPoint+1) <<", ave output = "<<ave<<endl;
	return ave;
}

void gelmanRubinStage1(outputDataType ODT,int numTimes)
{
	// This function will calculate the variance and mean numTimes times
	// for each parameter and then write output to a file to be read
	// in by Python function 'gelmanRubinStage2' to finish the calculations
	// and write the final file to disk.

	bool verbose = false;
	if (true)
	{
		cout<<"\nStarting stage 1 of Gelman-Rubin calculation"<<endl;
		cout<<" Calculating the Lc, mean and var for each parameter "<<numTimes<<" times."<<endl;
	}
	// find folder to write output file into
	string folder;
	folder = findDirectory(ODT.data_filename);
	// make output filename
	string outputFilename;
	string chainNumStr;
	chainNumStr=findChainNumberStr(ODT.data_filename);
	outputFilename = folder+"/gelmanRubin-chain_"+chainNumStr+".txt";

	// make an output file for the INS and load it up.
	ofstream file;
	file.open(outputFilename.c_str()) ;     //open the output file
	//write headers to file
	file << "Lc  longAN  e  To  Tc  period  inclination  argPeri  a_total  K";
//	for (int dataset=0;dataset<int(ODT.RVoffsets[0].size());++dataset)
//		file<<"  RVorigin_"<<dataset;
	file<<endl;

	// Clear memory for RVoffsets not used in GR calcs to make room
	// for extending the vectors to full size during thei calculations
//	vector<double> emptyVector;
//	emptyVector.push_back(0.0);
	if (true)
	{
		cout<<"clearing memory for RV offset and chiSquareds vectors!!!"<<endl;
		cout << "- size  before   = " << ODT.RVoffsets.size()     << endl;
		cout << "- capacity before = " << ODT.RVoffsets.capacity() << endl;
		cout << "- memory before = " << int((ODT.RVoffsets.capacity()*sizeof(double)+sizeof(ODT.RVoffsets))/1048576) <<" MB"<< endl;
	}
	ODT.RVoffsets.clear();
	vector<vector<double> >().swap(ODT.RVoffsets);
	ODT.chiSquareds.clear();
	vector<double>().swap(ODT.chiSquareds);

	//Check other vectors to see if they can also be deleted, if so kill them off
	if (ODT.longAN_degs[0]==ODT.longAN_degs.back())
	{
		cout<<"Clearing memory for longAN_degs"<<endl;
		ODT.longAN_degs.clear();
		vector<double>().swap(ODT.longAN_degs);
	}
	if (ODT.es[0]==ODT.es.back())
	{
		cout<<"Clearing memory for es"<<endl;
		ODT.es.clear();
		vector<double>().swap(ODT.es);
	}
	if (ODT.Ts[0]==ODT.Ts.back())
	{
		cout<<"Clearing memory for Ts"<<endl;
		ODT.Ts.clear();
		vector<double>().swap(ODT.Ts);
	}
	//In cases of DIonly, Tc is set to that of To, so they are the same vector
	//and no need to calculate the GR values of both.  Kill Tc and just do it for To
	if ((ODT.Tcs[0]==ODT.Ts.back())||((ODT.Tcs[0]==ODT.Ts[0])&&(ODT.Ts.back()==ODT.Ts.back())))
	{
		cout<<"Clearing memory for Tcs"<<endl;
		ODT.Tcs.clear();
		vector<double>().swap(ODT.Tcs);
	}
	if (ODT.periods[0]==ODT.periods.back())
	{
		cout<<"Clearing memory for periods"<<endl;
		ODT.periods.clear();
		vector<double>().swap(ODT.periods);
	}
	if (ODT.inclination_degs[0]==ODT.inclination_degs.back())
	{
		cout<<"Clearing memory for inclination_degs"<<endl;
		ODT.inclination_degs.clear();
		vector<double>().swap(ODT.inclination_degs);
	}
	if (ODT.argPeri_degs[0]==ODT.argPeri_degs.back())
	{
		cout<<"Clearing memory for argPeri_degs"<<endl;
		ODT.argPeri_degs.clear();
		vector<double>().swap(ODT.argPeri_degs);
	}
	if (ODT.a_totals[0]==ODT.a_totals.back())
	{
		cout<<"Clearing memory for a_totals"<<endl;
		ODT.a_totals.clear();
		vector<double>().swap(ODT.a_totals);
	}
	if (ODT.Ks[0]==ODT.Ks.back())
	{
		cout<<"Clearing memory for Ks"<<endl;
		ODT.Ks.clear();
		vector<double>().swap(ODT.Ks);
	}

	if (true)
	{
		cout << "- size  after   = " << ODT.RVoffsets.size()     << endl;
		cout << "- capacity after = " << ODT.RVoffsets.capacity() << endl;
	}

	vector<int> Lcs;

	GRfuncReturnType GRFRT_longAN_degs;
	if (ODT.longAN_degs.size()>1)
	{
		cout<<"Starting to calculate GR values for longAN_degs"<<endl;
		// load up string vectors for each parameters var,mean
		GRFRT_longAN_degs = gelmanRubinStage1func(ODT.longAN_degs,ODT.timesBeenHeres,numTimes);
		//clear memory from input data vector
		ODT.longAN_degs.clear();
		vector<double>().swap(ODT.longAN_degs);
		if (verbose)
			cout<<"finished performing stage 1 of GR calcs on LongAns"<<endl;
		Lcs.swap(GRFRT_longAN_degs.Lcs);
	}

	GRfuncReturnType GRFRT_es;
	if (ODT.es.size()>1)
	{
		cout<<"Starting to calculate GR values for es"<<endl;
		GRFRT_es = gelmanRubinStage1func(ODT.es,ODT.timesBeenHeres,numTimes);
		ODT.es.clear();
		vector<double>().swap(ODT.es);
		if (verbose)
			cout<<"finished performing stage 1 of GR calcs on es"<<endl;
		if (Lcs.size()<1)
			Lcs.swap(GRFRT_es.Lcs);
	}

	GRfuncReturnType GRFRT_Ts;
	if (ODT.Ts.size()>1)
	{
		cout<<"Starting to calculate GR values for Ts"<<endl;
		GRFRT_Ts = gelmanRubinStage1func(ODT.Ts,ODT.timesBeenHeres,numTimes);
		ODT.Ts.clear();
		vector<double>().swap(ODT.Ts);
		if (verbose)
			cout<<"finished performing stage 1 of GR calcs on Ts"<<endl;
		if (Lcs.size()<1)
			Lcs.swap(GRFRT_Ts.Lcs);
	}
	GRfuncReturnType GRFRT_Tcs;
	if (ODT.Tcs.size()>1)
	{
		cout<<"Starting to calculate GR values for Tcs"<<endl;
		GRFRT_Tcs = gelmanRubinStage1func(ODT.Tcs,ODT.timesBeenHeres,numTimes);
		ODT.Tcs.clear();
		vector<double>().swap(ODT.Tcs);
		if (verbose)
			cout<<"finished performing stage 1 of GR calcs on Tcs"<<endl;
		if (Lcs.size()<1)
			Lcs.swap(GRFRT_Tcs.Lcs);
	}
	GRfuncReturnType GRFRT_periods;
	if (ODT.periods.size()>1)
	{
		cout<<"Starting to calculate GR values for periods"<<endl;
		GRFRT_periods = gelmanRubinStage1func(ODT.periods,ODT.timesBeenHeres,numTimes);
		ODT.periods.clear();
		vector<double>().swap(ODT.periods);
		if (verbose)
			cout<<"finished performing stage 1 of GR calcs on v"<<endl;
		if (Lcs.size()<1)
			Lcs.swap(GRFRT_periods.Lcs);
	}

	GRfuncReturnType GRFRT_inclination_degs;
	if (ODT.inclination_degs.size()>1)
	{
		cout<<"Starting to calculate GR values for inclination_degs"<<endl;
		GRFRT_inclination_degs = gelmanRubinStage1func(ODT.inclination_degs,ODT.timesBeenHeres,numTimes);
		ODT.inclination_degs.clear();
		vector<double>().swap(ODT.inclination_degs);
		if (verbose)
			cout<<"finished performing stage 1 of GR calcs on inclination_degs"<<endl;
		if (Lcs.size()<1)
			Lcs.swap(GRFRT_inclination_degs.Lcs);
	}

	GRfuncReturnType GRFRT_argPeri_degs;
	if (ODT.argPeri_degs.size()>1)
	{
		cout<<"Starting to calculate GR values for argPeri_degs"<<endl;
		GRFRT_argPeri_degs = gelmanRubinStage1func(ODT.argPeri_degs,ODT.timesBeenHeres,numTimes);
		ODT.argPeri_degs.clear();
		vector<double>().swap(ODT.argPeri_degs);
		if (verbose)
			cout<<"finished performing stage 1 of GR calcs on argPeri_degs"<<endl;
		if (Lcs.size()<1)
			Lcs.swap(GRFRT_argPeri_degs.Lcs);
	}

	GRfuncReturnType GRFRT_a_totals;
	if (ODT.a_totals.size()>1)
	{
		cout<<"Starting to calculate GR values for a_totals"<<endl;
		GRFRT_a_totals = gelmanRubinStage1func(ODT.a_totals,ODT.timesBeenHeres,numTimes);
		ODT.a_totals.clear();
		vector<double>().swap(ODT.a_totals);
		if (verbose)
			cout<<"finished performing stage 1 of GR calcs on a_totals"<<endl;
		if (Lcs.size()<1)
			Lcs.swap(GRFRT_a_totals.Lcs);
	}

	GRfuncReturnType GRFRT_Ks;
	if (ODT.Ks.size()>1)
	{
		cout<<"Starting to calculate GR values for Ks"<<endl;
		GRFRT_Ks = gelmanRubinStage1func(ODT.Ks,ODT.timesBeenHeres,numTimes);
		ODT.Ks.clear();
		vector<double>().swap(ODT.Ks);
		if (verbose)
			cout<<"finished performing stage 1 of GR calcs on Ks"<<endl;
		if (Lcs.size()<1)
			Lcs.swap(GRFRT_Ks.Lcs);
	}

	int print10=0;
	for (int itt=0;itt<numTimes;++itt)
	{
		++print10;
		if ((print10==int(numTimes/10))&&verbose)
		{
			print10=0;
			cout<<"printing line # "<<itt<<" to GR file"<<endl;
		}
		file<<Lcs[itt]<<fixed<<std::setprecision(15);
		if (GRFRT_longAN_degs.vars.size()>1)
			file<<"  "<<GRFRT_longAN_degs.vars[itt]<<","<<GRFRT_longAN_degs.means[itt];
		else
			file<<"  "<<"0.0,0.0";
		if (GRFRT_es.vars.size()>1)
			file<<"  "<<GRFRT_es.vars[itt]<<","<<GRFRT_es.means[itt];
		else
			file<<"  "<<"0.0,0.0";
		if (GRFRT_Ts.vars.size()>1)
			file<<"  "<<GRFRT_Ts.vars[itt]<<","<<GRFRT_Ts.means[itt];
		else
			file<<"  "<<"0.0,0.0";
		if (GRFRT_Tcs.vars.size()>1)
			file<<"  "<<GRFRT_Tcs.vars[itt]<<","<<GRFRT_Tcs.means[itt];
		else
			file<<"  "<<"0.0,0.0";
		if (GRFRT_periods.vars.size()>1)
			file<<"  "<<GRFRT_periods.vars[itt]<<","<<GRFRT_periods.means[itt];
		else
			file<<"  "<<"0.0,0.0";
		if (GRFRT_inclination_degs.vars.size()>1)
			file<<"  "<<GRFRT_inclination_degs.vars[itt]<<","<<GRFRT_inclination_degs.means[itt];
		else
			file<<"  "<<"0.0,0.0";
		if (GRFRT_argPeri_degs.vars.size()>1)
			file<<"  "<<GRFRT_argPeri_degs.vars[itt]<<","<<GRFRT_argPeri_degs.means[itt];
		else
			file<<"  "<<"0.0,0.0";
		if (GRFRT_a_totals.vars.size()>1)
			file<<"  "<<GRFRT_a_totals.vars[itt]<<","<<GRFRT_a_totals.means[itt];
		else
			file<<"  "<<"0.0,0.0";
		if (GRFRT_Ks.vars.size()>1)
			file<<"  "<<GRFRT_Ks.vars[itt]<<","<<GRFRT_Ks.means[itt];
		else
			file<<"  "<<"0.0,0.0";
		file<<endl;
	}
	file.close();
	//cout<<"Finished calculating the mean and variance values for each param for stage one of the Gelman-Rubin calculation :-)"<<endl;
	cout<<"\n***************************************************************"<<endl;
	cout<<"Writing file to: "<<outputFilename<<endl;
	cout<<"***************************************************************\n"<<endl;
}

vector<double> fillOutDataVector(vector<double> data,vector<int> timesBeenHere)
{
	// put data into proper vector taking timesBeenHere into account
	vector<double> data2;
	//int timesBeenHereTOTAL = 0;
	for (int i=0;i<data.size();++i)
	{
		if (timesBeenHere[i]>0)
		{
			for (int j=0;j<timesBeenHere[i];++j)
			{
				if (data[i]>1e-10)
				{
					data2.push_back(data[i]);
					//timesBeenHereTOTAL+=1;
				}
			}
		}
	}
	return data2;
}

GRfuncReturnType gelmanRubinStage1func(vector<double> data,vector<int> timesBeenHere,int numTimes)
{
	GRfuncReturnType GRFRT;
	bool verbose = false;

	vector<double> data2;
	data2 = fillOutDataVector(data,timesBeenHere);
//	//$$$$$$$$
//	for (int i=0;i<data2.size();++i)
//	{
//		cout<<data2[i]<<endl;
//	}

	int print10=0;
	for (int itt=0;itt<numTimes;++itt)
	{

			int lastPoint;
		if (false)
			cout<<"For itt+1 "<<itt+1<<", numTimes = "<<numTimes<<", and data2.size() = "<<data2.size()<<": (data2.size()/numTimes)*(itt+1) = "<<(data2.size()/numTimes)*(itt+1)<<endl;
		lastPoint = int((data2.size()/numTimes)*(itt+1));
		if (false)
			cout<<"lastPoint = "<<lastPoint <<endl;
		GRFRT.Lcs.push_back(lastPoint);
		double ave,var;
		if (data2.size()>0)
		{
			if (data2.back()==data2[0])
			{
				ave=0;
				var =0;
			}
			else
			{
				ave = meanCalc(data2, lastPoint);
				if (isnan(ave))
					ave=0;
				var = varianceCalc(data2, lastPoint);
				if (isnan(var))
					var=0;
			}
		}
		else
		{
			ave=0;
			var=0;
		}
		GRFRT.means.push_back(ave);
		GRFRT.vars.push_back(var);
		if (false)
			cout<<"mean = "<<ave <<", var = "<< var <<endl;
		++print10;
		if ((print10==int(numTimes/10))&&(verbose))
		{
			print10=0;
			cout<<"For itt+1 "<<itt+1<<", numTimes = "<<numTimes<<", and data2.size() = "<<data2.size()<<": (data2.size()/numTimes)*(itt+1) = "<<(data2.size()/numTimes)*(itt+1)<<endl;
			cout<<"mean = "<<ave <<", var = "<< var <<endl;
		}
	}
	//clear memory for data2
	data2.clear();
	vector<double>().swap(data2);

	return GRFRT;
}

int corrLengthJumpyCalc(vector<double> data)
{
	double varALL = varianceCalc(data,(data.size()-1));
	double halfVarALL = varALL/2.0;
	int corrLength = 0;

	//find jump size for input ary
	int jump=10;
	if (data.size()>10000000)
		jump = 1000;
	else if (data.size()>100000)
		jump = 100;

	double varCurr;
	double varCurr2;
	for (int i=0;i<int(data.size()/jump);i++)
	{
		varCurr = varianceCalc(data,i*jump);
		if (varCurr>halfVarALL)
		{
			for (int j=(i-1)*jump; j<(i*jump) ; ++j)
			{
				varCurr2 = varianceCalc(data,j);
				if (varCurr>halfVarALL)
					corrLength = j+1;
					break;
			}
			break;
		}
	}
	return corrLength;
}

double varianceCalc(vector<double> data, int lastPoint)
{
	// This will calculate the "bias-corrected sample variance"
	// and uses an advanced "corrected two-pass algorithm" to reduce roundoff error

	// a modified version of function on pg 728 of Numerical Recipes 3rd
	double s,ep;
	int j,n=lastPoint+1;//add 1 to lastPoint to make it 1 indexed, instead of zero indexed
	double ave,var;
	ave = meanCalc(data, lastPoint);
	var=ep=0.0;
	//loop through all points to load up sums needed for "corrected two-pass algorithm", eq 14.1.8 pg 724
	for (j=0;j<n;j++)
	{
		s=data[j]-ave;
		ep+=s;
		var+=s*s;
	}
	//calc var and return it
	var=(var-ep*ep/n)/(n-1);
	return var;
}

string numSamplesStringMaker(int numSamples)
{
	std::stringstream ss;

	string numSamplesString = "NON YET!!";

	if (numSamples>=int(1e9))
	{
		int numMillions = numSamples/int(1e6);
		ss<<numMillions<<"-Billion";
		numSamplesString = ss.str();
		ss.clear();
		ss.str(std::string());
	}
	else if (numSamples>=int(1e6))
	{
		int numMillions = numSamples/int(1e6);
		ss<<numMillions<<"-Million";
		numSamplesString = ss.str();
		ss.clear();
		ss.str(std::string());
	}
	else if (numSamples>=int(1000))
	{
		int numThousands = numSamples/(int(1000));
		ss<<numThousands<<"-Thousand";
		numSamplesString = ss.str();
		ss.clear();
		ss.str(std::string());
	}
	else
	{
		ss<<numSamples;
		numSamplesString = ss.str();
		ss.clear();
		ss.str(std::string());
	}

	return numSamplesString;
}

string findDirectory(string str)
{
	size_t found;
	found=str.find_last_of("/\\");
	string folder;
	folder = str.substr(0,found);
	return folder;
}
string fileBasename(string str)
{
	size_t found;
	found=str.find_last_of("/\\");
	string file;
	file = str.substr(found+1);
	return file;
}

string findChainNumberStr(string str)
{
	size_t found;
	found=str.find_last_of(".");
	string chainNumStr;
	chainNumStr = str.substr(found-1,1);
	return chainNumStr;
}

string filenameEndAppend(string inputString, string endAppendString)
{
	// ******** THIS COULD PROBABLY BE RE-WRITTEN TO USE STR.FIND() AND STR.SUBSTR() ***

	bool verboseInternal = false;
	string outputString = "";
	int periodPos = inputString.find('.');

	if (periodPos>0)
	{
		for (int curPos=0; curPos<inputString.size();curPos+=1)
		{
			if (false)
				cout<<" char is:"<<inputString[curPos]<<" , for curPos "<<curPos<<endl;
			if (inputString[curPos]!='.')
				outputString.push_back(inputString[curPos]);
			else
				outputString += endAppendString+".";
		}
	}
	else
		outputString += inputString+endAppendString;
	if (verboseInternal)
		cout<<"\nFilename modified by filenameEndAppend:: \nInput string: "<<inputString<<"\nOutput string:"<<outputString<<endl;

	return outputString;
}

string filenamePrepend(string inputString, string prependString)
{
	// ******** THIS COULD PROBABLY BE RE-WRITTEN TO USE STR.FIND() AND STR.SUBSTR() ***


	bool verboseInternal = false;
	string outputString = "";
	int lastSlashPos = inputString.rfind('/');

	if (verboseInternal)
		cout<<"\nlastSlashPos = "<<lastSlashPos<<endl;

	if (lastSlashPos>0)
	{
		for (int curPos=0; curPos<inputString.size();curPos+=1)
		{
			if (false)
				cout<<" char is:"<<inputString[curPos]<<" , for curPos "<<curPos<<endl;
			if (curPos!=lastSlashPos)
				outputString.push_back(inputString[curPos]);
			else
				outputString += "/"+prependString;
		}
	}
	else
		outputString += prependString+inputString;

	if (verboseInternal)
		cout<<"\nFilename modified by filenamePrepend:: \nInput string: "<<inputString<<"\nOutput string:"<<outputString<<endl;

	return outputString;

}

string boolToStr(bool boolIN)
{
	string boolStr = "";

	if (boolIN==true)
		boolStr = "true";
	else if (boolIN==false)
		boolStr = "false";

	return boolStr;
}

string timeStr(int timeElapsed)
{
	std::stringstream sss;

	string timeString;
	bool verbose=false;

	if (verbose)
		cout<<" in timeStr func, and will find string for timeElapsed="<<timeElapsed<<endl;
	if ( timeElapsed < 1.0 )
	{
		sss<<"0 seconds";
		if (verbose)
			cout<<"timeElapsed determined to be ZERO seconds"<<timeElapsed<<endl;
	}
	else if ( timeElapsed < 60.0 )
	{
		sss<< timeElapsed << " seconds"<<endl;
		if (verbose)
			cout<<"timeElapsed determined to be seconds"<<timeElapsed<<endl;
	}
	else if ( timeElapsed>3600 )
	{
		double totalTime = ((double)timeElapsed)/3600.0;
		int hrs = (int)totalTime;
		int min = int(60*(totalTime-(double)hrs));
		string s1 = "";
		if (hrs!=1)
			s1="s";
		string s2="";
		if (min!=1)
			s2="s";
		sss<< hrs << " hour"<<s1<<" and "<< min << " minute"<<s2 ;
		if (verbose)
			cout<<"timeElapsed determined to be hours"<<endl;
	}
	else
	{
		double totalTime = ((double)timeElapsed)/60.0;
		int min = (int)totalTime;
		int sec = 60*(totalTime-double(min));
		string s2="";
		if (min!=1)
			s2="s";
		sss <<min <<" minute"<<s2<<" and "<<sec<< " seconds";
		if (verbose)
			cout<<"timeElapsed determined to be minutes"<<endl;
	}

	if (verbose)
		cout<<"finished finding timeString "<<endl;

	timeString = sss.str();
	if (verbose)
		cout<<"returning "<<timeString<<endl;

	return timeString;
}
