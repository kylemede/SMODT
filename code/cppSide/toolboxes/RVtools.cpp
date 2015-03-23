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
double VRcalcStarStar::VRcalculatorMassType()
{
	/**
	Calculate the residual velocity due to a companion star
	based on equation (47) in my thesis.
	 */
	double tempA = pow((2.0*PI*GravConst*KGperMsun*(Mass1+Mass2_s))/(SecPerYear*period_s),(1.0/3.0));
	double tempB = Mass2_s/Mass1;
	double tempC = sin(inclination_deg_s*(PI/180.0))/sqrt(1-e_s*e_s);
	// Calculate the Semi-major Amplitude for star-star system
	double Kss = tempA*tempB*tempC;
	K_s = Kss;
	double tempD = cos((TA_deg_s+argPeri_deg_s)*(PI/180.0))+e_s*cos(argPeri_deg_s*(PI/180.0));
	double VRc = Kss*tempD;

	if (false)
	{
		cout<<"Kss = "<<Kss<<endl;
		cout<<"TA_deg_s = "<<TA_deg_s<<endl;
	}

	return VRc;
}
double VRcalcStarStar::VRcalculatorSemiMajorType()
{
	/**
	 Calculate the residual velocity due to a companion star
	based on equation (49) in my thesis.
	 */
	bool verboseInternal = false;


	double Kss;
	if (fabs(K_s)>1)
	{
		Kss = K_s;
	}
	else
	{
		double tempA = (2.0*PI*a1*MperAU)/(SecPerYear*period_s);
		double tempC = sin(inclination_deg_s*(PI/180.0))/sqrt(1-e_s*e_s);
		// Calculate the Semi-major Amplitude for star-star system
		Kss = tempA*tempC;
		if (verboseInternal)
		{
			cout<<"\nperiod_s = "<< period_s<<endl;
			cout<<"e_s = "<< e_s<<endl;
			cout<<"inclination_deg = "<<inclination_deg_s<<endl;
			cout<<"a1 = "<<a1<<endl;
			cout<<"tempA = "<<tempA<<endl;
			cout<<"tempC = "<<tempC<<endl;
			cout<<"Kss = "<<Kss<<"\n"<<endl;
		}
		if (false)
		{
			double tempAa = pow((2.0*PI*GravConst*KGperMsun*(Mass1+Mass2_s))/(SecPerYear*period_s),(1.0/3.0));
			double tempBa = Mass2_s/Mass1;
			double tempCa = sin(inclination_deg_s*(PI/180.0))/sqrt(1-e_s*e_s);
			// Calculate the Semi-major Amplitude for star-star system
			double Kss2 = tempAa*tempBa*tempCa;
			double diff = Kss=Kss2;
			if (diff>0.1)
			{
				cout<<"Kss-Kss2 = "<<(Kss-Kss2)<<endl;
				cout<<"Kss = "<<Kss<<endl;
				cout<<"Kss2 = "<<Kss2<<endl;
			}
		}
	}
	if (K_s==0)
	{
		K_s = Kss;
		if (verboseInternal)
			cout<<"Just loaded K_s to match calculated Kss = "<<Kss<<endl;
	}
	double argPeri_deg_s_internal = argPeri_deg_s;//extra that isn't needed anymore I think
	double tempD1 = cos((TA_deg_s+argPeri_deg_s_internal)*(PI/180.0));
	double tempD2 = e_s*cos(argPeri_deg_s_internal*(PI/180.0));
	double tempD = tempD1+tempD2;
	double VRc = Kss*tempD;

	//cout<<"TA_deg_p = "<<TA_deg_p<<endl;//$$$$$$$$$$$$$$$$$$
	if (verboseInternal)
	{
		cout<<"\n # In VRcalculatorSemiMajorType #"<<endl;
		cout<<"period_s = "<< period_s<<endl;
		cout<<"inclination_deg = "<<inclination_deg_s<<endl;
		cout<<"a1 = "<<a1<<endl;
		cout<<"e_s = "<< e_s<<endl;
		cout<<"Provided K_s = "<<K_s<<endl;
		cout<<"Kss = "<<Kss<<endl;
		cout<<"TA_deg_s = "<<TA_deg_s<<endl;
//		cout<<"argPeri_deg_p = "<<argPeri_deg_p <<endl;
//		cout<<"cos((TA_deg_p)*(PI/180.0)) = "<<cos((TA_deg_p)*(PI/180.0))<<endl;
//		cout<<"tempD1 = "<<tempD1<<endl;
//		cout<<"tempD2 = "<<tempD2 <<endl;
//		cout<<"tempD = "<<tempD <<endl;
		cout<<"VRc = "<< VRc<<endl;
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
	//prep vector for returned VR vals
	vector<double> ResidualVels_s;
	generalTools GT;

	if (false)
		argPeri_deg_s = argPeri_deg_s-180.0;


//	double temp5 = period_s*period_s*SecPerYear*SecPerYear*GravConst*KGperMsun*(Mass1+Mass2_s);
//	double temp6 = 4.0*PI*PI;
//	a_total = pow((temp5/temp6),(1.0/3.0))/MperAU;
	//cout<<"a_total = "<<a_total<<endl;

	//instantiate and load up constant values for both TA and VR calc inputs
	TAcalcInputType TACIT;

	TACIT.e = e_s;
	TACIT.period = period_s;
	TACIT.verbose = verbose;
	TACIT.T = T_s;
	TACIT.Tc = Tc_s;

	for (int epoch=0; epoch<epochs_s.size(); ++epoch)
	{
		//if ( verbose==true )
		if (false)
			cout<< "\nWorking on epoch " <<epoch+1 <<endl;

		//prep varied inputs and call TA calc
		TACIT.t = epochs_s[epoch];
		TAcalcReturnType TACRT;
		TACRT = GT.TAcalculator(TACIT);
		//update value of TA for this epoch and call planet VR calc
		TA_deg_s = TACRT.TA_deg;

		double VRs;
		if (true)
			VRs = VRcalculatorSemiMajorType();
		else
			VRs = VRcalculatorMassType();
		ResidualVels_s.push_back(VRs);
		if (false)
		{
			cout<<"T ="<<TACIT.T <<", t = "<<TACIT.t<<", period = "<<period_s<<endl;
			cout<<"VRs = "<<VRs<<endl;
		}
	}
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
double VRcalcStarPlanet::VRcalculatorSemiMajorType()
{
	/**
	Calculate the residual velocity due to a companion star
	based on equation (49) in my thesis.
	*/
	bool verboseInternal = false;

//	if (verboseInternal)
//		cout<<"\n # In VRcalculatorSemiMajorType #"<<endl;
	double Ksp;
	if (K_p>0)
	{
		Ksp = K_p;
//		if (verboseInternal)
//			cout<<"Just loaded Ksp to match provided K_p = "<<K_p<<endl;
	}
	else
	{
		double tempA = (2.0*PI*a1*MperAU)/(SecPerYear*period_p);
		double tempC = sin(inclination_deg_p*(PI/180.0))/sqrt(1-e_p*e_p);
		// Calculate the Semi-major Amplitude for star-star system
		Ksp = tempA*tempC;
		if (verboseInternal)
		{
			cout<<"\nperiod_p = "<< period_p<<endl;
			cout<<"e_p = "<< e_p<<endl;
			cout<<"inclination_deg = "<<inclination_deg_p<<endl;
			cout<<"a1 = "<<a1<<endl;
			cout<<"tempA = "<<tempA<<endl;
			cout<<"tempC = "<<tempC<<endl;
			cout<<"Ksp = "<<Ksp<<"\n"<<endl;
		}

	}
	if (K_p==0)
	{
		K_p = Ksp;
//		if (verboseInternal)
//			cout<<"Just loaded K_p to match calculated Ksp = "<<K_p<<endl;
	}
	double argPeri_deg_p_internal = argPeri_deg_p;//extra that isn't needed anymore I think
	double tempD1 = cos((TA_deg_p+argPeri_deg_p_internal)*(PI/180.0));
	double tempD2 = e_p*cos(argPeri_deg_p_internal*(PI/180.0));
	double tempD = tempD1+tempD2;
	double VRp = Ksp*tempD;

	//cout<<"TA_deg_p = "<<TA_deg_p<<endl;//$$$$$$$$$$$$$$$$$$
	if (verboseInternal)
	{
		cout<<"\n # In VRcalculatorSemiMajorType #"<<endl;
		//cout<<"\nPI = "<<PI <<endl;
		//cout<<"GravConst = "<< GravConst<<endl;
		//cout<<"SecPerYear = "<< SecPerYear<<endl;
		cout<<"period_p = "<< period_p<<endl;
		cout<<"inclination_deg = "<<inclination_deg_p<<endl;
		cout<<"a1 = "<<a1<<endl;
		//cout<<"Mass2sinI_p = "<<Mass2sinI_p <<endl;
		//cout<<"KGperMsun = "<< KGperMjupiter<<endl;
		//cout<<"Mass1 = "<< KGperMsun<<endl;
		cout<<"e_p = "<< e_p<<endl;
		cout<<"Provided K_p = "<<K_p<<endl;
		cout<<"Ksp = "<<Ksp<<endl;
		cout<<"TA_deg_p = "<<TA_deg_p<<endl;
//		cout<<"argPeri_deg_p = "<<argPeri_deg_p <<endl;
//		cout<<"cos((TA_deg_p)*(PI/180.0)) = "<<cos((TA_deg_p)*(PI/180.0))<<endl;
//		cout<<"tempD1 = "<<tempD1<<endl;
//		cout<<"tempD2 = "<<tempD2 <<endl;
//		cout<<"tempD = "<<tempD <<endl;
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
	//cout<<"in multipEpochCalc for StarPlanet"<<endl;
	//prep vector for returned VR vals
	vector<double> ResidualVels_p;
	generalTools GT;

	if (false)
		argPeri_deg_p = argPeri_deg_p-180.0;

	//instantiate and load up constant values for both TA and VR calc inputs
	TAcalcInputType TACIT;
	TACIT.e = e_p;
	//cout<<"VRcalc: e = "<<TACIT.e<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	TACIT.period = period_p;
	TACIT.verbose = false;
	verbose = false;
	TACIT.T = T_p;
	TACIT.Tc = Tc_p;

	for (int epoch=0; epoch<epochs_p.size(); ++epoch)
	{
		//cout<< "\n###### Working on epoch " <<epoch<<" ################" <<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		//if ( verbose==true )
		if (false)
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
			VRp = VRcalculatorSemiMajorType();
		else
			VRp = VRcalculatorMassType();
		ResidualVels_p.push_back(VRp);
		//cout<<"#####################################################"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		if (false)
		{
			cout<<"T ="<<TACIT.T <<", t = "<<TACIT.t<<", period = "<<period_p<<endl;
			cout<<"VRp = "<<VRp<<endl;
		}
	}
	return ResidualVels_p;
}
