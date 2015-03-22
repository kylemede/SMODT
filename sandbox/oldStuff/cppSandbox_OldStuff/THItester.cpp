#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <sstream>
#include <string>
#include <fstream>

#include "Toolboxes/orbToolboxes.h"


int main()
{
	DItools DIt;

	double period = DIt.period = 100.0/365.242;
	double T = DIt.T = 0.0;
	double e = DIt.e = 0.0;
	bool verbose = DIt.verbose = false;
	double inclination_deg = DIt.inclination_deg = 0.0;
	double longAN_deg = DIt.longAN_deg = 0.0;
	double argPeri_deg = DIt.argPeri_deg = 0.0;
	double a_total = DIt.a_total = 100.0;
	double Sys_Dist_PC = DIt.Sys_Dist_PC = 10.0;

	vector<double> epochs_DI;

	for (int i=0;i<100;i++)
		epochs_DI.push_back(double(i)*1.0);

	// prepare TA calculator input structure
	// only epoch/t will be updated in the loop
	TAcalcInputType TACIT;
	TACIT.period  = period;
	TACIT.T =  T;
	TACIT.Tc =  0.0;//default value for DI as no need for Tc, only for RV data
	TACIT.e = e;
	TACIT.verbose = verbose;
	double E_deg;

	for ( int i=0; i<((int) epochs_DI.size()); i++ )
	{
		if (verbose)
			cout << "----------------------------------------------------------------------------" <<endl; //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

		// load TACIT structures with necessary values from the MEOCIT & MEOCRT structures
		TACIT.t  = epochs_DI[i];

		// instantiate the TACRT structure and pass the TACIT structure into
		// TAcalculator to calculate the E and TA for this epoch
		TAcalcReturnType TACRT;
		TACRT = TAcalculator(TACIT);

		// push freshly calculated E into OCIT for use in orbitCalculator
		E_deg = TACRT.E_deg;
		//cout<<"E_rad before orbitCalc = "<<E_deg*(PI/180.0)<<endl;
		DIt.E_deg = E_deg;

		// instantiate the OCRT structure and pass the OCIT structure into
		// orbitCalculator to run the model and pass back the complete OCRT.
		orbitCalcReturnType OCRT;
		OCRT = DIt.orbitCalculator();

		OCRT.x_model;
		OCRT.y_model;

		cout<<"For epoch "<<epochs_DI[i]<<", x = "<<OCRT.x_model<<", y = "<< OCRT.y_model<<endl;

//		// grab necessary values for the chi squared calculation from the MEOCIT structure
//		double SA_arcsec_measured_REAL = SAs_arcsec_observed[i];
//		double SA_mean_error = SA_errors[i];
//		double PA_deg_measured_REAL = PAs_deg_observed[i];
//		double PA_mean_error = PA_errors[i];
//
//		// convert SA & PA of data into x&y
//		double x_data = SA_arcsec_measured_REAL*cos(PA_deg_measured_REAL*(PI/180.0));
//		double y_data = SA_arcsec_measured_REAL*sin(PA_deg_measured_REAL*(PI/180.0));
//
//		if (verbose)
//		{
//			cout<<"\nSA_arcsec_measured_REAL = "<<SA_arcsec_measured_REAL <<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//			cout<<"PA_deg_measured_REAL = "<< PA_deg_measured_REAL<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//			cout<<"cos(PA_deg_measured_REAL*(PI/180.0)) = "<< cos(PA_deg_measured_REAL*(PI/180.0))<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//			cout<<"sin(PA_deg_measured_REAL*(PI/180.0)) = "<<sin(PA_deg_measured_REAL*(PI/180.0)) <<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//			cout<<"x_data = "<<x_data <<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//			cout<<"y_data = "<<y_data<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//
//			cout<<"SA_mean_error = "<< SA_mean_error<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//			cout<<"PA_mean_error = "<<PA_mean_error <<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//		}


	}




}
