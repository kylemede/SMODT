#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <fstream>
#include "orbToolboxes.h"
//#include "RVtools.h"  //called in orbToolboxes.h, so don't need to call here again

using namespace std;
/*! \file */
double RVtools::VRcalculatorSemiMajorType()
{
	/**
	 Calculate the residual velocity due to a companion star
	based on equation (49) in my thesis.
	 */
	bool verboseInternal = false;

	//check if modeling secondary RVs or the primary star's
	double aUse = a1;
	//cout<<"argPeri_deg = "<<argPeri_deg<<endl;
	//cout<<"primaryStarRVs = "<<primaryStarRVs<<endl;
	if (primaryStarRVs==false)
		aUse = a_total-a1;
	//Got K yet?
	double Kuse;
	if (fabs(K)>1)
	{
		Kuse = K;
	}
	else
	{

		double tempA = (2.0*PI*aUse*MperAU)/(SecPerYear*period);
		double tempC = sin(inclination_deg*(PI/180.0))/sqrt(1-e*e);
		// Calculate the Semi-major Amplitude for star-star system
		Kuse = tempA*tempC;
		K = Kuse;
		if (verboseInternal)
			cout<<"Just calculated Kuse and loaded into K, Kuse = "<<Kuse<<endl;
	}
	double tempD = cos((TA_deg+argPeri_deg)*(PI/180.0))+e*cos(argPeri_deg*(PI/180.0));
	double VR = K*tempD;

	if (verboseInternal)
	{
		cout<<"\n # In VRcalculatorSemiMajorType #"<<endl;
		//cout<<"period = "<< period<<endl;
		cout<<"inclination_deg = "<<inclination_deg<<endl;
		//cout<<"aUse = "<<aUse<<endl;
		//cout<<"e = "<< e<<endl;
		cout<<"Provided K = "<<K<<endl;
		cout<<"Kuse = "<<Kuse<<endl;
		cout<<"TA_deg = "<<TA_deg<<endl;
		cout<<"VR = "<< VR<<endl;
	}
	return VR;
}

vector<double> VRcalcStarStar::multiEpochCalc()
{
	/**
	 * Load up a True Anomaly Input Type and then for each provided epoch of
	 * data calculate the True Anomaly and then the associated radial velocity
	 * of the primary due to the secondary star.
	 */
	bool verboseInternal = false;
	//prep vector for returned VR vals
	vector<double> ResidualVels_s;
	generalTools GT;

	//instantiate and load up constant values for both TA and VR calc inputs
	TAcalcInputType TACIT;

	TACIT.e = e_s;
	TACIT.period = period_s;
	TACIT.verbose = verbose;
	TACIT.T = T_s;
	TACIT.Tc = Tc_s;

	//Load up planet specific params into parent level to calc VR.
	//If there are no differences between companion star and companion planet VR
	// calcs, then these funcs could be merged when I have time.
	K=K_s;
	period = period_s;
	e = e_s;
	inclination_deg = inclination_deg_s;
	argPeri_deg = argPeri_deg_s;

	if (verboseInternal)
		cout<<"#################### start of RV multiEpochCalc ###########################"<<endl;
	for (int epoch=0; epoch<epochs_s.size(); ++epoch)
	{
		if (verboseInternal)
			cout<< "\nWorking on epoch " <<epoch+1 <<endl;

		//prep varied inputs and call TA calc
		TACIT.t = epochs_s[epoch];
		TAcalcReturnType TACRT;
		TACRT = GT.TAcalculator(TACIT);
		//update value of TA for this epoch and call planet VR calc
		TA_deg_s = TACRT.TA_deg;

		double VRc;
		TA_deg = TA_deg_s;
		VRc = VRcalculatorSemiMajorType();
		ResidualVels_s.push_back(VRc);
		if (verboseInternal)
		{
			cout<<"T ="<<TACIT.T <<", t = "<<TACIT.t<<", period = "<<period_s<<endl;
			cout<<"VRc = "<<VRc<<endl;
		}
	}
	// Load value output by VRcalculator into K_s.  Will be the same if K_s was not zero, else will be the updated calculated version
	K_s = K;
	return ResidualVels_s;
}


vector<double> VRcalcStarPlanet::multiEpochCalc()
{
	/**
	 * Load up a True Anomaly Input Type and then for each provided epoch of
	 * data calculate the True Anomaly and then the associated radial velocity
	 * of the primary due to the companion planet.
	 */
	bool verboseInternal = false;
	if (verboseInternal)
		cout<<"in multipEpochCalc for StarPlanet"<<endl;
	//prep vector for returned VR vals
	vector<double> ResidualVels_p;
	generalTools GT;

	//instantiate and load up constant values for both TA and VR calc inputs
	TAcalcInputType TACIT;
	TACIT.e = e_p;
	TACIT.period = period_p;
	TACIT.verbose = false;
	verbose = false;
	TACIT.T = T_p;
	TACIT.Tc = Tc_p;

	//Load up planet specific params into parent level to calc VR.
	//If there are no differences between companion star and companion planet VR
	// calcs, then these funcs could be merged when I have time.
	K=K_p;
	period = period_p;
	e = e_p;
	inclination_deg = inclination_deg_p;
	argPeri_deg = argPeri_deg_p;

	if (verboseInternal)
		cout<<"#################### start of RV multiEpochCalc ###########################"<<endl;
	if (verboseInternal)
	{
		cout<<"input K = "<<K<<endl;
		cout<<"input inc = "<<inclination_deg<<endl;
	}
	for (int epoch=0; epoch<epochs_p.size(); ++epoch)
	{
		if (verboseInternal)
			cout<< "\nWorking on epoch " <<epoch+1 <<endl;

		//prep varied inputs and call TA calc
		TACIT.t = epochs_p[epoch];
		TAcalcReturnType TACRT;
		TACRT = GT.TAcalculator(TACIT);
		//update value of TA for this epoch and call planet VR calc
		TA_deg_p = TACRT.TA_deg;
		E_deg_p = TACRT.E_deg;

		double VRp;
		if (verboseInternal)
			cout<<"sending in inc = "<<inclination_deg<<endl;
		TA_deg = TA_deg_p;
		VRp = VRcalculatorSemiMajorType();
		ResidualVels_p.push_back(VRp);
		if (verboseInternal)
		{
			cout<<"T ="<<TACIT.T <<", t = "<<TACIT.t<<", period = "<<period_p<<endl;
			cout<<"VRp = "<<VRp<<endl;
		}
	}
	// Load value output by VRcalculator into K_p.  Will be the same if K_p was not zero, else will be the updated calculated version
	K_p = K;
	return ResidualVels_p;
}
