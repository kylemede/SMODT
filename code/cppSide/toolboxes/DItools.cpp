#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <sstream>
#include <string>
#include <fstream>
#include "orbToolboxes.h"

using namespace std;
/*! \file */
orbitCalcReturnType DItools::orbitCalculator()
{
	/**
	 * This will calculate the X and Y locations for a given epoch
	 */

	double verboseInternal = false; //set to true for debugging
	orbitCalcReturnType OCRT;

	// Pull in values from OCIT
	double a_arcsec = a_use/Sys_Dist_PC;
	double omega_rad = argPeri_deg*(PI/180.0);
	double Omega_rad = longAN_deg*(PI/180.0);
	double i_rad = inclination_deg*(PI/180.0);
	double E_rad = E_deg*(PI/180.0);

	// calculate all the Thiele-Innes constants in ["]
	double A = a_arcsec*(cos(Omega_rad)*cos(omega_rad)-sin(Omega_rad)*sin(omega_rad)*cos(i_rad));
	double B = a_arcsec*(sin(Omega_rad)*cos(omega_rad)+cos(Omega_rad)*sin(omega_rad)*cos(i_rad));
	double F = a_arcsec*(-cos(Omega_rad)*sin(omega_rad)-sin(Omega_rad)*cos(omega_rad)*cos(i_rad));
	double G = a_arcsec*(-sin(Omega_rad)*sin(omega_rad)+cos(Omega_rad)*cos(omega_rad)*cos(i_rad));

	// calculate X&Y
	// The coordinates of the unit orbital ellipse in the true plane (Binnendijk)
	double X = cos(E_rad)-e;
	double Y = sqrt(1.0-e*e)*sin(E_rad);

	// Calculate the predicted x&y in ["]
	double hackNeg = 1.0;//$$$$$$$$$$$$$$$$$$$$$$$$$$$
	OCRT.x_model = hackNeg*(A*X +F*Y);
	OCRT.y_model = B*X +G*Y;

	if (verboseInternal)
	{
		cout<<"\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		cout<<"period = "<<period<<endl;
		cout<<"a_use = "<<a_use <<endl;
		cout<<"a_arcsec = "<<a_arcsec<<endl;
		cout<<"Sys_Dist_PC = "<<Sys_Dist_PC <<endl;
		cout<<"a_arcsec = "<< a_arcsec<<endl;
		cout<<"omega = "<< argPeri_deg<<endl;
		cout<<"Omega = "<< longAN_deg<<endl;
		cout<<"i = "<< inclination_deg<<endl;
		cout<<"E_deg = "<<E_deg <<endl;
		cout<<"E_rad = "<<E_rad <<endl;
		cout<<"A = "<<A <<endl;
		cout<<"B = "<<B <<endl;
		cout<<"F = "<< F<<endl;
		cout<<"G = "<<G <<endl;
		cout<<"X = "<< X<<endl;
		cout<<"Y = "<<Y <<endl;
		cout<<"x_model = "<< OCRT.x_model<<endl;
		cout<<"y_model = "<< OCRT.y_model<<endl;
	}

	return OCRT;
}

multiEpochOrbCalcReturnType DItools::multiEpochOrbCalc()
{
	/**
	 * This will calculate the (X,Y) and/or (Separation Angle,Position Angle)
	 * for each epoch of data provided.  It takes advantage of the Thiele-Innes
	 * equations.
	 */

	multiEpochOrbCalcReturnType  MEOCRT;
	generalTools GT;
	bool verboseInternal;
	verboseInternal = verbose;
	verboseInternal = false;//$$$$$$$$$$$$$$


	double chi_squared_total = 0.0;
	double chi_squared_total2 = 0.0;

	// add pi to the value of the argument of periapsis to convert it to the
	// value for the star instead of the value input which is for the companion.
	//NOTE: This is not fully implemented into the settings files, so left as false for now.
	if (true)
		argPeri_deg = argPeri_deg+180.0;
	//This next step is to put it inside one circle, but not really needed as it
	//should be handled properly either way.  So, just for safe keeping I guess.
	if (argPeri_deg>360.0)
		argPeri_deg = argPeri_deg-360.0;

	// convert the semi-major axis from the period using K3 and the provided masses
	if (a_total==0)
	{
		double temp5 = period*period*SecPerYear*SecPerYear*GravConst*KGperMsun*(Mass1+Mass2);
		double temp6 = 4.0*PI*PI;
		MEOCRT.a_total = a_total = pow((temp5/temp6),(1.0/3.0))/MperAU;
	}
	else
		MEOCRT.a_total = a_total;
	if (MEOCRT.a_total>1e4)
		MEOCRT.a_total=0;

	a_use = a_total;
//	// NOTE: the a1&a2 values are not being saved, kill it?? $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//	// use the mass ratio to convert the a_total into
//	// its a1 and a2 components if the masses are provided
//	if ( (Mass1!=1)&&(Mass2!=1))
//	{
//		// this means these two parameters are non-default values and
//		// the system must then be a binary star system, thus calculate the
//		// a1 and a2 values of the system.
//		// NOTE: Mass2 MUST < Mass1, ie. primary is bigger.
//		semiMajorType SMT_in;
//		semiMajorType SMT_out;
//		SMT_in.a1 = 0;
//		SMT_in.a2 = 0;
//		SMT_in.a_total = MEOCRT.a_total;
//		SMT_in.period = period;
//		SMT_in.Mass1 = Mass1;
//		SMT_in.Mass2 = Mass2;
//		SMT_out = GT.semiMajorConverter(SMT_in);
//		if (Mass2>0.0124054547)
//			a_use = SMT_out.a_total;
//		else
//			a_use = SMT_out.a1;
//		//cout<<"a_use = "<<a_use<<endl;
//		//cout<<"Mass2 = "<<Mass2<<endl;
//	}
//	else
//		cout<<"ERROR! no values for the two object masses were provided!!!"<<endl;

	// prepare TA calculator input structure
	// only epoch/t will be updated in the loop
	TAcalcInputType TACIT;
	TACIT.period  = period;
	TACIT.T =  T;
	TACIT.Tc =  0.0;//default value for DI as no need for Tc, only for RV data
	TACIT.e = e;
	TACIT.verbose = verbose;

	//verbose=true;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	// loop over all epochs of data and find predicted x&y values
	// to compare to those from the data
	for ( int i=0; i<((int) epochs_DI.size()); i++ )
	{
		if (verboseInternal)
			cout << "----------------------------------------------------------------------------" <<endl; //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

		// load TACIT structures with necessary values from the MEOCIT & MEOCRT structures
		TACIT.t  = epochs_DI[i];

		// instantiate the TACRT structure and pass the TACIT structure into
		// TAcalculator to calculate the E and TA for this epoch
		TAcalcReturnType TACRT;
		TACRT = GT.TAcalculator(TACIT);

		// push freshly calculated E into OCIT for use in orbitCalculator
		E_deg = TACRT.E_deg;

		// instantiate the OCRT structure and pass the OCIT structure into
		// orbitCalculator to run the model and pass back the complete OCRT.
		orbitCalcReturnType OCRT;
		OCRT = orbitCalculator();

		//MEOCRT.x_models.push_back(OCRT.x_model);
		//MEOCRT.y_models.push_back(OCRT.y_model);

		// grab necessary values for the chi squared calculation from the MEOCIT structure
		double SA_arcsec_measured_REAL = SAs_arcsec_observed[i];
		double SA_mean_error = SA_errors[i];
		double PA_deg_measured_REAL = PAs_deg_observed[i];
		double PA_mean_error = PA_errors[i];

		// convert SA & PA of data into x&y
		double x_data = SA_arcsec_measured_REAL*cos(PA_deg_measured_REAL*(PI/180.0));
		double y_data = SA_arcsec_measured_REAL*sin(PA_deg_measured_REAL*(PI/180.0));

		if (verboseInternal)
		{
			//cout<<"\nSA_arcsec_measured_REAL = "<<SA_arcsec_measured_REAL <<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			//cout<<"PA_deg_measured_REAL = "<< PA_deg_measured_REAL<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			//cout<<"cos(PA_deg_measured_REAL*(PI/180.0)) = "<< cos(PA_deg_measured_REAL*(PI/180.0))<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			//cout<<"sin(PA_deg_measured_REAL*(PI/180.0)) = "<<sin(PA_deg_measured_REAL*(PI/180.0)) <<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			cout<<"E_deg = "<<E_deg<<endl;
			//cout<<"x_data = "<<x_data <<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			//cout<<"y_data = "<<y_data<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			//cout<<"SA_mean_error = "<< SA_mean_error<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			//cout<<"PA_mean_error = "<<PA_mean_error <<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		}

		// convert SA & PA errors into x&y errors
		double tempAa = SA_mean_error*cos(PA_deg_measured_REAL*(PI/180.0));
		double tempA = tempAa*tempAa;
		double tempBa = SA_arcsec_measured_REAL*sin(PA_deg_measured_REAL*(PI/180.0))*PA_mean_error*(PI/180.0);
		double tempB = tempBa*tempBa;
		double x_data_error = sqrt(tempA+tempB);
		double x_data_inv_var = 1.0/(x_data_error*x_data_error);

		double tempCa =SA_mean_error*sin(PA_deg_measured_REAL*(PI/180.0));
		double tempC = tempCa*tempCa;
		double tempDa = SA_arcsec_measured_REAL*cos(PA_deg_measured_REAL*(PI/180.0))*PA_mean_error*(PI/180.0);
		double tempD = tempDa*tempDa;
		double y_data_error = sqrt(tempC+tempD);
		double y_data_inv_var = 1.0/(y_data_error*y_data_error);

		// calculated associated chiSquared values for each
		double x_chi_squared = GT.chiSquaredCalc(x_data, x_data_inv_var, OCRT.x_model);
        double y_chi_squared = GT.chiSquaredCalc(y_data, y_data_inv_var, OCRT.y_model);

        // Add them to get the updated total
        chi_squared_total = chi_squared_total + x_chi_squared + y_chi_squared;

        if (verboseInternal)
        {
        	cout << "x_data = "<< x_data<<", x_data_error = " <<x_data_error << ", OCRT.x_model = " << OCRT.x_model<<endl;
        	cout << "y_data = "<< y_data<<", y_data_error = " <<y_data_error << ", OCRT.y_model = " << OCRT.y_model<<endl;
        	//double x2 = -OCRT.x_model;
        	//double y2 = -OCRT.y_model;
        	//cout<< "x2 = "<<x2<<", y2 = "<<y2<<endl;
        	//double x_chi_squared2 = GT.chiSquaredCalc(x_data, x_data_inv_var, x2);
        	//double y_chi_squared2 = GT.chiSquaredCalc(y_data, y_data_inv_var, y2);
        	cout << "x chi_squared = "<< x_chi_squared<<  ", y chi_square = "<< y_chi_squared<< endl;
        	//cout << "x chi_squared2 = "<< x_chi_squared2<<  ", y chi_square2 = "<< y_chi_squared2<< endl;
        	//chi_squared_total2 = chi_squared_total2 + x_chi_squared2 + y_chi_squared2;
        	cout << fixed <<" in loop chi_squared_total = "<<chi_squared_total <<endl;
        	//cout << fixed <<" in loop chi_squared_total2 = "<<chi_squared_total2 <<endl;
        }
		// increment to next epoch
	}
	// Loop done, finished for all epochs

	//cout <<"349: output  chi_squared_total = "<<chi_squared_total <<endl;
	MEOCRT.chi_squared_total = chi_squared_total;

	return MEOCRT;
}
