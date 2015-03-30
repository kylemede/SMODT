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

	double Kinternal;
	if (fabs(K)>1)
	{
		Kinternal = K;
	}
	else
	{
		double tempA = (2.0*PI*a1*MperAU)/(SecPerYear*period);
		double tempC = sin(inclination_deg*(PI/180.0))/sqrt(1-e*e);
		// Calculate the Semi-major Amplitude for star-star system
		Kinternal = tempA*tempC;
		K = Kinternal;
		if (verboseInternal)
			cout<<"Just calculated Kinternal and loaded into K, Kinternal = "<<Kinternal<<endl;
	}
	double argPeri_deg_internal = argPeri_deg;//extra that isn't needed anymore I think
	double tempD = cos((TA_deg+argPeri_deg_internal)*(PI/180.0))+e*cos(argPeri_deg_internal*(PI/180.0));
	double VR = K*tempD;

	if (verboseInternal)
	{
		cout<<"\n # In VRcalculatorSemiMajorType #"<<endl;
		//cout<<"period = "<< period<<endl;
		cout<<"inclination_deg = "<<inclination_deg<<endl;
		//cout<<"a1 = "<<a1<<endl;
		//cout<<"e = "<< e<<endl;
		cout<<"Provided K = "<<K<<endl;
		cout<<"Kinternal = "<<Kinternal<<endl;
		cout<<"TA_deg = "<<TA_deg<<endl;
		cout<<"VR = "<< VR<<endl;
	}
	return VR;
}
double VRcalcStarStar::VRcalculatorMassType()
{
	/**
	Calculate the residual velocity due to a companion star
	based on equation (47) in my thesis.
	 */
	double Kinternal;
	if (fabs(K_s)>0)
		Kinternal = K_s;
	else
	{
		double tempA = pow((2.0*PI*GravConst*KGperMsun*(Mass1+Mass2_s))/(SecPerYear*period_s),(1.0/3.0));
		double tempB = Mass2_s/Mass1;
		double tempC = sin(inclination_deg_s*(PI/180.0))/sqrt(1-e_s*e_s);
		// Calculate the Semi-major Amplitude for star-star system
		Kinternal = tempA*tempB*tempC;
		K_s = Kinternal;
	}
	double tempD = cos((TA_deg_s+argPeri_deg_s)*(PI/180.0))+e_s*cos(argPeri_deg_s*(PI/180.0));
	double VRc = Kinternal*tempD;

	if (false)
	{
		cout<<"Kinternal = "<<Kinternal<<endl;
		cout<<"TA_deg_s = "<<TA_deg_s<<endl;
	}
	return VRc;
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

	if (false)
		argPeri_deg_s = argPeri_deg_s+180.0;

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
		if (true)
		{
			TA_deg = TA_deg_s;
			VRc = VRcalculatorSemiMajorType();
		}
		else
			VRc = VRcalculatorMassType();
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

double VRcalcStarPlanet::VRcalculatorMassType()
{
	/**

	Calculate the residual velocity due to a companion planet
	based on equation (48) in my thesis.
	*/
	bool verboseInternal = false;

	double Ksp;
	if (K_p>0)
		Ksp = K_p;
	else
	{
		double tempA = pow((2.0*PI*GravConst)/(SecPerYear*period_p),(1.0/3.0));
		double tempB = (Mass2sinI_p*KGperMsun)/pow(Mass1*KGperMsun,(2.0/3.0));
		double tempC = 1.0/sqrt(1-e_p*e_p);
		// Calculate the Semi-major Amplitude for star-planet system
		Ksp = tempA*tempB*tempC;
	}
	if (K_p==0)
		K_p = Ksp;
	double tempD = cos((TA_deg_p+argPeri_deg_p)*(PI/180.0))+e_p*cos(argPeri_deg_p*(PI/180.0));
	double VRp = Ksp*tempD;

	if (verboseInternal)
	{
		cout<<"\n # In VRcalculator #"<<endl;
		//cout<<"\nPI = "<<PI <<endl;
		//cout<<"GravConst = "<< GravConst<<endl;
		//cout<<"SecPerYear = "<< SecPerYear<<endl;
		//cout<<"period_p = "<< period_p<<endl;
		//cout<<"Mass2sinI_p = "<<Mass2sinI_p <<endl;
		//cout<<"KGperMsun = "<< KGperMjupiter<<endl;
		//cout<<"Mass1 = "<< KGperMsun<<endl;
		//cout<<"e_p = "<< e_p<<endl;
		cout<<"Provided K_p = "<<K_p<<endl;
		cout<<"Ksp = "<<Ksp<<endl;
		cout<<"TA_deg_p = "<<TA_deg_p<<endl;
		cout<<"argPeri_deg_p = "<<argPeri_deg_p <<endl;
		cout<<"cos((TA_deg_p)*(PI/180.0)) = "<<cos((TA_deg_p)*(PI/180.0))<<endl;
		cout<<"cos((TA_deg_p+argPeri_deg_p)*(PI/180.0)) = "<<cos((TA_deg_p+argPeri_deg_p)*(PI/180.0))<<endl;
		cout<<"cos(argPeri_deg_p*(PI/180.0)) = "<<cos(argPeri_deg_p*(PI/180.0)) <<endl;
		cout<<"tempD = "<<tempD <<endl;
		cout<<"VRp = "<< VRp<<endl;
	}

	return VRp;
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

	if (false)
		argPeri_deg_p = argPeri_deg_p-180.0;

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
		if (true)
		{
			if (verboseInternal)
				cout<<"sending in inc = "<<inclination_deg<<endl;
			TA_deg = TA_deg_p;
			VRp = VRcalculatorSemiMajorType();
		}
		else
			VRp = VRcalculatorMassType();
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
