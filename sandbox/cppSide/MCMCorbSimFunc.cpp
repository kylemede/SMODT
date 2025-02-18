#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <time.h>
#include "Toolboxes/orbToolboxes.h" //Both the DI and RV tools headers are called in here, so it is an all in one toolbox header call
#include "MCMCorbSimFunc.h"
#include "rnd/kylesGaussRand.h"

using namespace std;

void MCMCorbFuncObj::simulator()
{
	// variables for the success rate print block in chain loop
	int printTime = SSO.numSamples/SSO.numSamplePrints;
	int printCount = 0;
	int printsDone = 0;
	std::stringstream ss;
	std::stringstream SSlog;
	// variables for param changing
	int samplesSinceParmaChange = 1;
	// acceptance rate calcs
	double acceptedCounter = 0;
	double acceptCalcTime = double(SSO.numSamples)/(SSO.numSamplePrints*10.0);
	vector<int> paramsVariedRecentlyAry;
	vector<int> acceptedIntsRecentlyAry;
	string acceptString;
	double samplesTillAcceptRateCalc = 0;
	double latestAcceptRate = 0;
	int timesNONEpassed = 0;
	int timesBeenHere = 1;
	int paramBeingVaried = 2;
	bool latestParamsSaved;
	//double K_p_errorPercent = 0;
	double chiSquaredMin_DI=SSO.chiSquaredMax;
	if (SSO.RVonly==true)
		chiSquaredMin_DI=0;
	double chiSquaredMin_RV=SSO.chiSquaredMax;
	if (SSO.DIonly==true)
		chiSquaredMin_RV=0;
	chiSquaredMin = SSO.chiSquaredMax;
	double DI_chiSquared = 0;
	double RV_chiSquared = 0;
	double TOTAL_chiSquared = 0;

	bestOrbit = 0;
	double inc_prior = 12345;
	double P_prior = 12345;
	double e_prior = 12345;
	double priors_ratio=12345;
	double likelihood_ratio =12345;
	double prior_likelihood_ratio = 12345;
	double alpha = 12345;
	double RHS = 0;
	double chiSquare_latest = 12345678;
	string accepted = "?";

	// Printing block for loop/starting params
	string startParms;
	ss<<"\n****************************************"<<endl;
	ss<<"Important loop counters and params:"<<endl;
	ss<<"printTime = "<< printTime<<endl;
	ss<<"acceptCalcTime = "<<acceptCalcTime <<endl;
	string silentStr = boolToStr(SSO.silent);
	ss<<"silent = "<< silentStr <<endl;
	string verboseStr = boolToStr(SSO.verbose);
	ss<<"verbose = "<<verboseStr<<endl;
	ss<<"****************************************"<<endl;
	startParms = ss.str();
	ss.clear();
	ss.str(std::string());
	cout<<startParms;
	SSlog<< startParms;

	//start/choose generator(s)
	CRandomSFMT1 RanGen(randSeed);
	StochasticLib1 RanGen2(randSeed);
	//this is the most advanced uniform random number generator, which combines SFMT and Mother-Of-All.

	// instantiate and load up a second RVdo for proposed planet vals in simStar=true cases
	RVdataObj RVdo2;
	RVdo2 = RVdo;

	double Kp_calculated=0;
	double Ks_calculated=0;

	// set up starting values for input params
	double inclination_deg_latest = start_inc_deg;
	double longAN_deg_latest = start_longAN;
	double argPeri_deg_latest = start_argPeri;
	double e_latest = start_e;
	double period_latest = start_period;
	double K_latest = start_K;
	double a_total_latest = start_a_total;
	double T_latest = start_T;
	double Tc_latest = start_Tc;
	double sqrtESinomega_latest = sqrt(e_latest)*sin((PI/180.0)*argPeri_deg_latest);
	double sqrtECosomega_latest = sqrt(e_latest)*cos((PI/180.0)*argPeri_deg_latest);
//	// calculate starting Tc from provided To
//	eccArgPeri2ToTcType EATT;
//	EATT.period = period_latest;
//	EATT.argPeri_deg = argPeri_deg_latest;
//	EATT.e = e_latest;
//	EATT.To = T_latest;
//	EATT.Tc=0;
//	EATT = eccArgPeri2ToTcCalc(EATT);
//	double Tc_latest;
//	Tc_latest = EATT.Tc;

	string startParms2;
	ss<<"***********************************************"<<endl;
	ss<<"Starting values for model input parameters:"<<endl;
	ss<<"inclination_deg_latest = "<<inclination_deg_latest <<endl;
	ss<<"longAN_deg_latest = "<< longAN_deg_latest<<endl;
	ss<<"argPeri_deg_latest = "<< argPeri_deg_latest<<endl;
	ss<<"e_latest = "<<e_latest <<endl;
	ss<<"period_latest = "<< period_latest<<endl;
	ss<<"T_latest = "<<fixed <<T_latest<<endl;
	ss<<"Tc_latest = "<<fixed <<Tc_latest<<endl;
	ss<<"a_total_latest = "<<a_total_latest<<endl;
	// set up starting offsets if needed
	vector<double> RVoffsets_latest;
	vector<double> RVoffsets_proposed;
	if (SSO.DIonly==false)
	{
		ss<<"K_latest = "<<K_latest<<endl;
		RVoffsets_latest = start_offsets;
		RVoffsets_proposed = RVoffsets_latest;
		for (int dataset=0; dataset<RVdo.epochs_RV.size();++dataset)
			ss<<"RV dataset "<<dataset<<", RVoffset = "<<RVoffsets_latest[dataset]<<endl;
	}
	ss<<"\nNumber of varying parameters for DI = "<< numDIparams<<", RV = " <<numRVparams <<", 3D = " <<numParams<<endl;
	ss<<"one_over_nu_DI = "<<one_over_nu_DI <<" ,one_over_nu_RV = "<<one_over_nu_RV <<", one_over_nu_TOTAL = "<<one_over_nu_TOTAL <<endl;
	ss<<"***********************************************"<<endl;
	startParms2 = ss.str();
	ss.clear();
	ss.str(std::string());
	cout<<startParms2;
	SSlog<< startParms2;

	// setting initially 'proposed' states equal to initial 'latest' states
	double longAN_deg_proposed = longAN_deg_latest;
	double e_proposed = e_latest;
	double T_proposed = T_latest;
	double Tc_proposed = Tc_latest;
	double period_proposed = period_latest;
	double inclination_deg_proposed = inclination_deg_latest;
	double argPeri_deg_proposed = argPeri_deg_latest;
	double a_total_proposed = a_total_latest;
	double K_proposed = K_latest;
	double sqrtESinomega_proposed = sqrtESinomega_latest;
	double sqrtECosomega_proposed = sqrtECosomega_latest;
	double Sys_Dist_PC_proposed = SYSdo.Sys_Dist_PC;
	double Mass1_proposed = SYSdo.Mass1;
	double planet_MsinI_proposed = SYSdo.planet_MsinI;
	double star_Mass2_proposed = SYSdo.star_Mass2;

	double a_total_curr = 0;

	// Reserve data for the output data vectors to safe time resizing them
	ODT= odtStart(ODT, SSO.numSamples);

	bool ALLpassed;

	// ***** Start the samples loop *****
	int sample;
	for ( sample=1; sample<(SSO.numSamples+1); sample++)
	{
		latestParamsSaved=false;

		// Set proposed values to latest
		longAN_deg_proposed = longAN_deg_latest;
		e_proposed = e_latest;
		T_proposed = T_latest;
		Tc_proposed = Tc_latest;
		period_proposed = period_latest;
		inclination_deg_proposed = inclination_deg_latest;
		argPeri_deg_proposed = argPeri_deg_latest;
		a_total_proposed = a_total_latest;
		Sys_Dist_PC_proposed = SYSdo.Sys_Dist_PC;
		Mass1_proposed = SYSdo.Mass1;
		planet_MsinI_proposed = SYSdo.planet_MsinI;
		star_Mass2_proposed = SYSdo.star_Mass2;
		sqrtESinomega_proposed = sqrtESinomega_latest;
		sqrtECosomega_proposed = sqrtECosomega_latest;
		//double argPeri_deg_proposed_use;

		// block to control printing success rate to screen
		printCount = printCount + 1;
		if ( printCount==printTime )
		{
			printsDone = printsDone+1;
			printCount = 0;
			int sampleUSE = sample;

			// create a time object and find current time
			time_t rawtime;
			struct tm * timeinfo;
			time ( &rawtime );
			timeinfo = localtime ( &rawtime );
			// find number of samples accepted so far
			int len=ODT.es.size();
			// print number of samples accepted so far and current time
			string printLine;
			ss <<"\n"<< int(acceptedCounter)<<"/"<<sample<<" Successful. "<<printsDone<<"/"<<SSO.numSamplePrints<<" completed at ";
			ss << asctime (timeinfo);
			ss << "Latest acceptance rate = "<<latestAcceptRate<<endl;
			ss << "Latest param being varied = "<<paramBeingVaried<<", timesBeenHere = "<<timesBeenHere<<endl;
			ss << "Times NONE of params passed = "<<timesNONEpassed<<endl;
			ss << "Largest allowed reduced chiSquareds: DI = "<<SSO.chiSquaredMax*one_over_nu_DI <<", RV = "<<SSO.chiSquaredMax*one_over_nu_RV  <<", Total = "<<SSO.chiSquaredMax*one_over_nu_TOTAL  <<endl;
			ss << "latest reduced chiSquareds: DI = "<< DI_chiSquared*one_over_nu_DI<<", RV = "<<RV_chiSquared*one_over_nu_RV <<", Total = "<< TOTAL_chiSquared*one_over_nu_TOTAL<<endl;
			ss << "LOWEST reduced chiSquareds: DI = "<< chiSquaredMin_DI*one_over_nu_DI <<", RV = "<< chiSquaredMin_RV*one_over_nu_RV <<", Total = "<< chiSquaredMin*one_over_nu_TOTAL <<endl;
			printLine = ss.str();
			ss.clear();
			ss.str(std::string());
			cout<<"\n#######################################################################################"<<endl;
			cout<<printLine;
			if (true)
				cout<<acceptString;

			string printLine2;
			ss<<"\n----------------------------------------------"<<endl;
			ss<<"inclination_deg_proposed = " <<inclination_deg_proposed  <<endl;
			ss<<"longAN_deg_proposed = " << longAN_deg_proposed  <<endl;
			ss<<"argPeri_deg_proposed = "<< argPeri_deg_proposed <<endl;
			ss<<"e_proposed = "<< e_proposed <<endl;
			ss<<"period_proposed = "<< period_proposed  <<endl;
			ss<<"T_proposed = "<<fixed <<T_proposed  <<endl;
			ss<<"Tc_proposed = "<<fixed <<Tc_proposed  <<endl;
			ss<<"K_proposed = "<<K_proposed<<endl;
			ss<<"a_total_proposed = "<<a_total_proposed<<endl;
			ss<<"----------------------------------------------"<<endl;

			if (true)
			{
				for (int dataset=0; dataset<RVdo.epochs_RV.size();++dataset)
					ss<<"RV dataset "<<dataset<<", RVoffset = "<<RVoffsets_latest[dataset]<<endl;
			}

			printLine2 = ss.str();
			ss.clear();
			ss.str(std::string());
			cout<<printLine2;

			if (false)
			{
				cout<<"-----------------------------------------------"<<endl;
				cout<< "chiSquare_latest = "<<chiSquare_latest <<endl;
				cout<<"likelihood_ratio = "<<likelihood_ratio <<endl;
				cout<<"alpha = "<< alpha <<endl;
				cout<<"RHS = "<< RHS<<endl;
				cout<<"accepted = "<<accepted<<endl;
				cout<<"----------------------------------------------"<<endl;
			}
			if (true)
			{
				ss<<"\nKp_calculated = "<<Kp_calculated<<endl;
				ss<<"Ks_calculated = "<<Ks_calculated<<"\n"<<endl;
			}
			cout<<"#######################################################################################"<<endl;
			SSlog<< printLine;
			SSlog<< printLine2;
		}

		// Generate random numbers in the required ranges for the inputs to the orbCalc
		ALLpassed = true;//just the starting value, set to false in proposal block if out of range
		int dataset;
		if (sample>1)
		{
			// ******* Determine which param to vary *************************
			if (paramBeingVaried==0)
				longAN_deg_proposed = RanGen.UniformRandom(longAN_deg_latest-longAN_deg_sigma,longAN_deg_latest+longAN_deg_sigma);
	//		if (paramBeingVaried==1)
	//			e_proposed = RanGen.UniformRandom(e_latest-e_sigma,e_latest+e_sigma);
	//		else if (paramBeingVaried==1)
	//			sqrtESinomega_proposed = RanGen.UniformRandom(sqrtESinomega_latest-sqrtESinomega_sigma,sqrtESinomega_latest+sqrtESinomega_sigma);
			else if (paramBeingVaried==1)
			{
				if (SSO.eMAX<0.3)
					sqrtESinomega_proposed = RanGen.UniformRandom(sqrtESinomega_latest-sqrtESinomega_sigma,sqrtESinomega_latest+sqrtESinomega_sigma);
				else
					e_proposed = RanGen.UniformRandom(e_latest-e_sigma,e_latest+e_sigma);
			}
			else if (paramBeingVaried==2)
			{
				if (SSO.TcStepping)
				{
					Tc_proposed = RanGen.UniformRandom(Tc_latest-T_sigma,Tc_latest+T_sigma);// thus between a full period before first observation and the time of first observation
				}
				else
				{
					T_proposed = RanGen.UniformRandom(T_latest-T_sigma,T_latest+T_sigma);// thus between a full period before first observation and the time of first observation
				}
			}
			else if (paramBeingVaried==3)
				period_proposed = RanGen.UniformRandom(period_latest-period_sigma,period_latest+period_sigma);//  [yrs]
			else if (paramBeingVaried==4)
				inclination_deg_proposed = RanGen.UniformRandom(inclination_deg_latest-inclination_deg_sigma,inclination_deg_latest+inclination_deg_sigma);
	//		else if (paramBeingVaried==5)
	//			argPeri_deg_proposed = RanGen.UniformRandom(argPeri_deg_latest-argPeri_deg_sigma,argPeri_deg_latest+argPeri_deg_sigma);
	//		else if (paramBeingVaried==5)
	//			sqrtECosomega_proposed = RanGen.UniformRandom(sqrtECosomega_latest-sqrtECosomega_sigma,sqrtECosomega_latest+sqrtECosomega_sigma);
			else if (paramBeingVaried==5)
			{
				if (SSO.eMAX<0.3)
					sqrtECosomega_proposed = RanGen.UniformRandom(sqrtECosomega_latest-sqrtECosomega_sigma,sqrtECosomega_latest+sqrtECosomega_sigma);
				else
					argPeri_deg_proposed = RanGen.UniformRandom(argPeri_deg_latest-argPeri_deg_sigma,argPeri_deg_latest+argPeri_deg_sigma);
			}
			else if (paramBeingVaried==6)
				a_total_proposed = RanGen.UniformRandom(a_total_latest-a_total_sigma,a_total_latest+a_total_sigma);
			else if (paramBeingVaried==7)
				K_proposed = RanGen.UniformRandom(K_latest-K_sigma,K_latest+K_sigma);
			else if (paramBeingVaried>7)
			{
				dataset = paramBeingVaried-8;
				if (SSO.RVoffsetMAXs[dataset]!=0)
				{
					RVoffsets_proposed[dataset] = RanGen.UniformRandom(RVoffsets_latest[dataset]-offset_sigmas[dataset],RVoffsets_latest[dataset]+offset_sigmas[dataset]);
					if (false)
						ss<<"dataset = " <<dataset<<", RVoffsets_proposed[dataset] = "<<RVoffsets_proposed[dataset] <<", RVoffsetMIN = "<<SSO.RVoffsetMINs[dataset]<<", MAX = "<< SSO.RVoffsetMAXs[dataset]<<endl;
				}
			}
		}
		//Convert proposed sqrt(e)cos(omega) and ..sin(omega) into useful e and omega
		if ((SSO.eMAX<0.3)&&(SSO.eMAX!=0))
		{
			e_proposed = sqrtESinomega_proposed*sqrtESinomega_proposed+sqrtECosomega_proposed*sqrtECosomega_proposed;
			argPeri_deg_proposed = (180.0/PI)*atan2(sqrtESinomega_proposed,sqrtECosomega_proposed);
			if (SSO.argPeri_degMAX>180)
			{
				if (argPeri_deg_proposed<0)
					argPeri_deg_proposed = argPeri_deg_proposed+360.0;
			}
		}
		if (SSO.DIonly==false)
		{
			if (TMAX!=0)
			{
				if (SSO.TcStepping)
				{
					//cout<<" values in SYSdo: planet_T = "<< SYSdo.planet_T<<", star_T = "<<SYSdo.star_T <<endl;
					if (SSO.simulate_StarPlanet==true)
						T_proposed = SYSdo.planet_T;
					else
						T_proposed = SYSdo.star_T;
				}
				else
				{
					//cout<<" values in SYSdo: planet_Tc = "<<SYSdo.planet_Tc <<", star_Tc = "<< SYSdo.star_Tc<<endl;
					if (SSO.simulate_StarPlanet==true)
						Tc_proposed = SYSdo.planet_Tc;
					else
						Tc_proposed = SYSdo.star_Tc;
				}
			}
			else
			{
				if (SSO.simulate_StarPlanet==true)
				{
					T_proposed = SYSdo.planet_T;
					Tc_proposed = SYSdo.planet_Tc;
				}
				else
				{
					T_proposed = SYSdo.star_T;
					Tc_proposed = SYSdo.star_Tc;
				}
			}
			//cout<<"T in = "<<T_proposed <<", Tc in = "<< Tc_proposed<<endl;

			//update non-updated T if it was 0 in the dictionary
			if ((T_proposed==0)||(Tc_proposed==0))
			{
				// calculate starting To from provided Tc, or the other way around it TcStepping==False
				// calculate the To value from proposed values (omega, e, P & Tc)
				eccArgPeri2ToTcType EATT;
				EATT.period = period_proposed;
				EATT.argPeri_deg = argPeri_deg_proposed;
				EATT.e = e_proposed;
				EATT.To = T_proposed;
				EATT.Tc=Tc_proposed;
				EATT = eccArgPeri2ToTcCalc(EATT);
				//cout<<"T in = "<<T_proposed <<", Tc in = "<< Tc_proposed<<", T out = "<<EATT.To <<", Tc out = "<< EATT.Tc<<endl;
				T_proposed = EATT.To;
				Tc_proposed = EATT.Tc;
			}
		}
		else
			Tc_proposed = T_proposed;

		// calculate the To value from proposed values (omega, e & Tc)
//		eccArgPeri2ToTcType EATT;
//		EATT.period = period_proposed;
//		EATT.argPeri_deg = argPeri_deg_proposed;
//		EATT.Tc = Tc_proposed;
//		EATT.e = e_proposed;
//		EATT = eccArgPeri2ToTcCalc(EATT);
//		T_proposed = EATT.To;

		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//		inclination_deg_proposed =  92.548065;
//		longAN_deg_proposed =  0;
//		argPeri_deg_proposed = 121.613235;
//		T_proposed = 2454757.00787
		//Tc_proposed = 2454756.73134;
//		e_proposed = 0.679959;
////		period_proposed = 164.072068;
		//Tc_proposed = 2455651.491942;
		//K_proposed = 465.0347;
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

		// Generate Gaussian values for the sys dist and masses
		if (false)
		{
			Sys_Dist_PC_proposed = RanGen2.NormalTrunc(SYSdo.Sys_Dist_PC,0.5*SYSdo.Sys_Dist_PC_error,SYSdo.Sys_Dist_PC_error);
			Mass1_proposed = RanGen2.NormalTrunc(SYSdo.Mass1,0.5*SYSdo.Mass1_error,SYSdo.Mass1_error);
			// load up mass2 with correct value depending on star or planet companion
			if (SSO.simulate_StarPlanet==false)
				star_Mass2_proposed = RanGen2.NormalTrunc(SYSdo.star_Mass2,0.5*SYSdo.star_Mass2_error,SYSdo.star_Mass2_error);
			else
				planet_MsinI_proposed = RanGen2.NormalTrunc(SYSdo.planet_MsinI,0.5*SYSdo.planet_MsinI_error,SYSdo.planet_MsinI_error);
		}
		// **** Done producing 'proposed' versions of all params being varied this round ****

		// ****** Check all are good   *********
		if ((SSO.longAN_degMAX!=0)&&(SSO.RVonly==false))
		{
			if ((longAN_deg_proposed>SSO.longAN_degMAX)||(longAN_deg_proposed<SSO.longAN_degMIN))
				ALLpassed=false;
		}
		if (SSO.eMAX!=0)
		{
			if ((e_proposed>SSO.eMAX)||(e_proposed<SSO.eMIN))
				ALLpassed=false;
		}
		if (TMAX!=0)
		{
			if (SSO.TcStepping)
			{
				if ((Tc_proposed>TMAX)||(Tc_proposed<TMIN))
					ALLpassed=false;
			}
			else
			{
				if ((T_proposed>TMAX)||(T_proposed<TMIN))
					ALLpassed=false;
			}
		}
		if ((SSO.periodMIN!=0)&&(SSO.periodMAX!=0))
		{
			if ((period_proposed>SSO.periodMAX)||(period_proposed<SSO.periodMIN))
				ALLpassed=false;
		}
		if (SSO.inclination_degMAX!=0)
		{
			if ((inclination_deg_proposed>=SSO.inclination_degMAX)||(inclination_deg_proposed<=SSO.inclination_degMIN))
				ALLpassed=false;
		}
		if (SSO.argPeri_degMAX!=0)
		{
			if ((argPeri_deg_proposed>SSO.argPeri_degMAX)||(argPeri_deg_proposed<SSO.argPeri_degMIN))
				ALLpassed=false;
		}
		if (SSO.a_totalMAX!=0)
		{
			if ((a_total_proposed>SSO.a_totalMAX)||(a_total_proposed<SSO.a_totalMIN))
				ALLpassed=false;
		}
		if (SSO.DIonly==false)
		{
			if (vary_K)
			{
				if ((K_proposed>SSO.K_MAX)||(K_proposed<SSO.K_MIN))
					ALLpassed=false;
			}
			for (int dataset=0; dataset<RVoffsets_latest.size();++dataset)
			{
				if (SSO.RVoffsetMAXs[dataset]!=0)
				{
					if ((SSO.RVoffsetMINs[dataset]>RVoffsets_proposed[dataset])||(SSO.RVoffsetMAXs[dataset]<RVoffsets_proposed[dataset]))
						ALLpassed = false;
				}
			}
		}

		// **** Done checking 'proposed' versions of all params being varied this round ****

		// if all are good, move on to calculating orbit.
		if(ALLpassed)
		{
			timesNONEpassed = 0;
			DIt.inclination_deg = inclination_deg_proposed;
			DIt.longAN_deg = longAN_deg_proposed;
//			if (true)//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//				DIt.argPeri_deg = argPeri_deg_proposed+180.0;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//			else//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			DIt.argPeri_deg = argPeri_deg_proposed;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			DIt.e = e_proposed;
			DIt.period = period_proposed; //  [yrs]
			DIt.a_total = a_total_proposed;//[AU]
			DIt.T = T_proposed;
			// load up params drawn from fixed gaussians
			DIt.Sys_Dist_PC = Sys_Dist_PC_proposed ;
			DIt.Mass1 = Mass1_proposed;
			if (SSO.simulate_StarStar==true)
				DIt.Mass2 =  star_Mass2_proposed;
			else
				DIt.Mass2 = planet_MsinI_proposed ;

			if ( SSO.silent==false )
				cout<<"ALL random numbers loaded"<<endl;

			// Brick for generating Mass1, Mass2 & Sys_Dist values from $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			// Gaussian distributions.									$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

			if ( SSO.silent==false )
			{
				string printLine2;
				ss<< "\ninclination_deg = " <<DIt.inclination_deg  <<"\n";
				ss<< "longAN_deg = " << DIt.longAN_deg  <<"\n";
				ss<<  "argPeri_deg = "<< DIt.argPeri_deg <<"\n";
				ss<<  "e = "<< DIt.e <<"\n";
				ss<<  "period = "<< DIt.period  <<"\n";
				ss<<  "T = "<< DIt.T  <<"\n\n";
				printLine2 = ss.str();
				ss.clear();
				ss.str(std::string());
				cout<<printLine2;
			}

			multiEpochOrbCalcReturnType MEOCRT;
			if (SSO.RVonly==false)
			{
				//cout<<"In DI block"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				// #### DO STUFF FOR DI ORBIT CALCNS #####
				if ( SSO.silent==false )
					cout<<"Calculating DI orbit for this round "<<endl;


				// Call the orbCalc to have it apply the model to the inputs and produce outputs
				//DIt.verbose=SSO.verbose;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				MEOCRT = DIt.multiEpochOrbCalc();
				a_total_curr = MEOCRT.a_total;
				if (a_total_curr>1e4)
					cout<<"\n\n!!!!! a_total_curr>1e4 in DI section!!!!\n\n"<<endl;

				// Calculate the reduced chiSquared from the returned chiSquared
				DI_chiSquared = MEOCRT.chi_squared_total;


				// update lowest DI reduced chiSquared if current one is lower
//				if ( DI_chiSquared<chiSquaredMin_DI )
//					chiSquaredMin_DI = DI_chiSquared;
			}
			else
			{
				//cout<<"In DI else block"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				//chiSquaredMin_DI = 0;
			}

			if (SSO.DIonly==false)
			{
				//cout<<"In RV block"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				if ( SSO.silent==false )
					cout<<"Calculating RV residuals for this round"<<endl;
				RV_chiSquared = 0.0;
				// ##### DO STUFF FOR RV CALCNS ######
				vector<vector<double> > VRp_vector2;
				vector<vector<double> > VRs_vector2;

				// load up params drawn from fixed gaussians
				RVdo.Sys_Dist_PC = Sys_Dist_PC_proposed ;
				RVdo.Mass1 = Mass1_proposed ;

				// Load up ss or sp parts of RVdo with current trials
				// param values as needed.
				if (SSO.simulate_StarPlanet==true)
				{
					if ( SSO.silent==false )
						cout<<"loading up input params for star-planet rv calcs"<<endl;
					RVdo.planet_e = DIt.e ;
					RVdo.planet_T = DIt.T ;
					RVdo.planet_Tc = Tc_proposed;
					RVdo.planet_P = DIt.period ;
					//RVdo.planet_MsinI  = DIt.Mass2 ;
					if (vary_K)
						RVdo.planet_K = K_proposed;
					RVdo.planet_argPeri = argPeri_deg_proposed ;
					RVdo.planet_inc = DIt.inclination_deg ;
					RVdo.planet_MsinI = planet_MsinI_proposed ;
				}
				if (SSO.simulate_StarStar==true)
				{
					if ( SSO.silent==false )
						cout<<"loading up input params for star-star rv calcs"<<endl;
					RVdo.star_e = DIt.e ;
					RVdo.star_T = DIt.T ;
					RVdo.star_Tc = Tc_proposed;
					RVdo.star_P = DIt.period ;
					//RVdo.star_Mass2 = DIt.Mass2 ;
					RVdo.star_argPeri = argPeri_deg_proposed ;
					RVdo.star_inc = DIt.inclination_deg ;
					RVdo.star_Mass2 =  star_Mass2_proposed;
				}
				// get residual velocities for companion planet if needed
				if (RVdo.planet_P!=0 )
				{
					if ( SSO.silent==false )
						cout<<"Starting to calculate residual vel for star-planet"<<endl;
					// instantiate S-P calc object and load up its params
					//VRcalcStarPlanet VRCsp;
					//generate latest params for planet VR calcs from Gaussians
					if (false)
					{
						if (SSO.simulate_StarStar==true)
						{
							RVdo2.planet_e = RanGen2.NormalTrunc(RVdo.planet_e,RVdo.planet_e_error,3.0*RVdo.planet_e_error);
							RVdo2.planet_T = RanGen2.NormalTrunc(RVdo.planet_T,RVdo.planet_T_error,3.0*RVdo.planet_T_error);
							RVdo2.planet_P = RanGen2.NormalTrunc(RVdo.planet_P,RVdo.planet_P_error,3.0*RVdo.planet_P_error);
							if (RVdo.planet_K>0)
								RVdo2.planet_K=0;
							else
								RVdo2.planet_K = RanGen2.NormalTrunc(RVdo.planet_K,RVdo.planet_K_error,3.0*RVdo.planet_K_error);
							RVdo2.planet_argPeri = RanGen2.NormalTrunc(RVdo.planet_argPeri,RVdo.planet_argPeri_error,3.0*RVdo.planet_argPeri_error);
							//Also do the same for the planet_inc?????
							//load up VRCsp obj with updated RVdo2 vals
							VRCsp = VRcalcStarPlanetLoadUp(RVdo2);
						}
					}
					else
						VRCsp = VRcalcStarPlanetLoadUp(RVdo);
					//K_p_errorPercent = VRCsp.K_p_error/VRCsp.K_p;
					//cout<<"Back from VRcalcStarPlanetLoadUp"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$
					//cout<<"there were "<<RVdo.epochs_RV.size()<<" datasets found in the RVdata file"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$
					VRCsp.verbose = false;
					// run through all RV data sets and calc residuals for it
					for (int dataset=0; dataset<int(RVdo.epochs_RV.size());++dataset)
					{
						if ( SSO.silent==false )
							cout<<"Calculating P-S residuals for dataset "<<(dataset+1)<<"/"<<int(RVdo.epochs_RV.size())<<endl;
						VRCsp.epochs_p = RVdo.epochs_RV[dataset];
						vector<double> VRp_vector;
						//cout<<"about to call multiEpochCalc for dataset "<<dataset<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
						VRp_vector = VRCsp.multiEpochCalc();
						VRp_vector2.push_back(VRp_vector);
					}
					if (SSO.simulate_StarPlanet==true)
					{
						a_total_curr = VRCsp.a_total;
						if (vary_K==false)
							K_proposed = VRCsp.K_p;
					}
					Kp_calculated = VRCsp.K_p;
					if ( SSO.silent==false )
						cout<<"K_p = "<<VRCsp.K_p<<endl;
				}

				// get residual velocities for companion star if needed
				if (RVdo.star_P!=0)
				{
					if ( SSO.silent==false )
						cout<<"Starting to calculate residual vel for star-star"<<endl;
					// instantiate S-S calc object and load up its params
					//VRcalcStarStar VRCss;

					//generate latest params for planet VR calcs from Gaussians
					if (false)
					{
						if (SSO.simulate_StarStar==true)
						{
							RVdo2.star_e = RanGen2.NormalTrunc(RVdo.star_e,RVdo.star_e_error,3.0*RVdo.star_e_error);
							RVdo2.star_T = RanGen2.NormalTrunc(RVdo.star_T,RVdo.star_T_error,3.0*RVdo.star_T_error);
							RVdo2.star_P = RanGen2.NormalTrunc(RVdo.star_P,RVdo.star_P_error,3.0*RVdo.star_P_error);
							RVdo2.star_argPeri = RanGen2.NormalTrunc(RVdo.star_argPeri,RVdo.star_argPeri_error,3.0*RVdo.star_argPeri_error);
							RVdo2.star_inc = RanGen2.NormalTrunc(RVdo.star_inc,RVdo.star_inc_error,3.0*RVdo.star_inc_error);
							//load up VRCsp obj with updated RVdo2 vals
							VRCss = VRcalcStarStarLoadUp(RVdo2);
						}
					}
					else
						VRCss = VRcalcStarStarLoadUp(RVdo);
					VRCss.verbose = false;
					// run through all RV data sets and calc residuals for it
					for (int dataset=0; dataset<int(RVdo.epochs_RV.size());++dataset)
					{
						if ( SSO.silent==false )
							cout<<"Calculating S-S residuals for dataset "<<(dataset+1)<<"/"<<int(RVdo.epochs_RV.size())<<endl;
						VRCss.epochs_s = RVdo.epochs_RV[dataset];
						vector<double> VRs_vector;
						VRs_vector = VRCss.multiEpochCalc();
						VRs_vector2.push_back(VRs_vector);
					}
					if (SSO.simulate_StarStar==true)
					{
						a_total_curr = VRCss.a_total;
						if (vary_K==false)
							K_proposed = VRCss.K_s;
					}
					Ks_calculated = VRCss.K_s;
					if ( SSO.silent==false )
						cout<<"K_s = "<<VRCss.K_s<<endl;
				}

				// go through all residual vels and calculate the RV chiSquareds
				for (int dataset=0; dataset<RVdo.epochs_RV.size();++dataset)
				{
					if ( SSO.silent==false )
						cout<<"\nStarting to calculate chiSquared from residuals for dataset# "<< dataset<<endl;


//						RVoffsets_proposed[0]=161.5790;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//						RVoffsets_proposed[1]=-0.544084;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					for (int epoch=0; epoch<RVdo.epochs_RV[dataset].size(); ++epoch)
					{
						double planetVR = 0;
						double companionStarVR = 0;
						// replace zeroed values above with real values for each as needed
						if (RVdo.planet_P!=0 )
						{
							planetVR = VRp_vector2[dataset][epoch];
							//cout<<"planetVR for epoch "<<epoch <<" is "<<planetVR <<endl;//$$$$$$$$$$$$$$$$
						}
						if (RVdo.star_P!=0)
						{
							companionStarVR = VRs_vector2[dataset][epoch];
							//cout<<"companionStarVR for epoch "<<epoch <<" is "<<companionStarVR <<endl;//$$$$$$$$$$$$$$$$
						}

						double updatedRV_inv_var = RVdo.RV_inv_var[dataset][epoch];
						//double updatedRV_inv_var = 1.0/((1.0/RVdo.RV_inv_var[dataset][epoch])+(K_p_errorPercent*planetVR)*(K_p_errorPercent*planetVR));
						if (false)
						{
							cout<< "RV_inv_var = "<<RVdo.RV_inv_var[dataset][epoch] <<",planetVR  ="<< planetVR <<endl;//<<", K_p_errorPercent = " << K_p_errorPercent <<endl;
							cout<<"updatedRV_inv_var = "<<updatedRV_inv_var <<", RV_inv_var = "<< RVdo.RV_inv_var[dataset][epoch]<<endl;
						}
						double  RV_chiSquared_cur = chiSquaredCalc((RVdo.RVs[dataset][epoch]-RVoffsets_proposed[dataset]),updatedRV_inv_var,(planetVR+companionStarVR));
//							double  RV_chiSquared_cur = chiSquaredCalc((RVdo.RVs[dataset][epoch]-RVoffsets_proposed[dataset]),RVdo.RV_inv_var[dataset][epoch],(planetVR+companionStarVR));
						RV_chiSquared = RV_chiSquared + RV_chiSquared_cur;
						if ( SSO.silent==false )
						{
							cout<<"\noffset = "<< RVoffsets_proposed[dataset]<<endl;
							cout<<"RVdo.RVs[dataset][epoch] = "<<RVdo.RVs[dataset][epoch]<<", ("<<RVdo.RVs[dataset][epoch]<<" - "<<RVoffsets_proposed[dataset]<<")="<<(RVdo.RVs[dataset][epoch]-RVoffsets_proposed[dataset]) <<", planetVR= "<< planetVR<<", companionStarVR= "<< companionStarVR<<endl;
							cout<<"Difference = "<<RVdo.RVs[dataset][epoch]-RVoffsets_proposed[dataset]-planetVR-companionStarVR<<endl;
							cout<<"ChiSquared for this RV is = "<<RV_chiSquared_cur<<endl;
							cout<<"Total NON-reducedChiSquared so far is = "<<RV_chiSquared<<endl;
						}
					}//End epoch loop
				}//End dataset loop




				// update lowest reduced RV chisquared value found if current one is lower
//				if ( RV_chiSquared<chiSquaredMin_RV )
//					chiSquaredMin_RV = RV_chiSquared;
			}//End RV calc block
			else
			{
				//cout<<"In RV else block"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				RVoffsets_latest.push_back(0);
				//chiSquaredMin_RV=0;
			}


			// Do TOTAL chiSquared value calcs
			TOTAL_chiSquared  = DI_chiSquared+RV_chiSquared;

			//
			// Determine if the orbit should be accepted
			//
			//Calculate priors ratio
			//e_prior = 1.0;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			if ((period_latest*365.242)<1000.0)
				e_prior = 1.0;
			else
				e_prior = e_latest/DIt.e;
			if ((SSO.inclination_degMIN!=0)&&(SSO.inclination_degMAX!=0))
				inc_prior = sin(DIt.inclination_deg*(PI/180.0))/sin(inclination_deg_latest*(PI/180.0));
			else
				inc_prior = 1.0;
			//P_prior = 1.0;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			if (false)//((SSO.periodMIN!=0)&&(SSO.periodMAX!=0))
				P_prior = DIt.period/period_latest;
			else
				P_prior = 1.0;
			priors_ratio = P_prior*e_prior*inc_prior;
			//Calculate likelihood ratio
			likelihood_ratio = exp((chiSquare_latest - TOTAL_chiSquared)/2.0);
			RHS = priors_ratio*likelihood_ratio;
			alpha = RanGen.UniformRandom(0.0, 1.0);
			//alpha=0.0;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

			if ( alpha<=RHS )
			{
				if (alpha>prior_likelihood_ratio)
					cout<<"WARNING: alpha>RHS mess up, "<< alpha<<">"<< RHS<<", but was still accepted."<<endl;

				accepted = "true";
				//if (SSO.silent==false)
				if (false)
				{
					cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
					cout<< "chiSquare_latest = "<<chiSquare_latest <<endl;
					cout<< "TOTAL_chiSquared = "<<TOTAL_chiSquared <<endl;
					cout<<"priors_ratio = "<<priors_ratio <<endl;
					cout<<"e_prior = "<<e_prior <<endl;
					cout<<"inc_prior = "<< inc_prior<<endl;
					cout<<"P_prior = "<< P_prior<<endl;
					cout<<"likelihood_ratio = "<<likelihood_ratio <<endl;
					cout<<"alpha = "<< alpha <<endl;
					cout<<"RHS = "<< RHS<<endl;
					cout<<"accepted = "<<accepted<<endl;
					cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
				}

				if (paramsVariedRecentlyAry.size()<acceptCalcTime)
				{
					paramsVariedRecentlyAry.push_back(paramBeingVaried);
					acceptedIntsRecentlyAry.push_back(1);
				}
				acceptedCounter +=1;
				chiSquare_latest = TOTAL_chiSquared;

				//store location of best orbit out of all accepted
				if ( TOTAL_chiSquared<=chiSquaredMin)
				{
					chiSquaredMin = TOTAL_chiSquared;
					bestOrbit = acceptedCounter-1;
					chiSquaredMin_DI = DI_chiSquared;
					chiSquaredMin_RV = RV_chiSquared;
				}

				latestParamsSaved=true;
				// store inputs
				ODT.longAN_degs.push_back(DIt.longAN_deg);
				ODT.es.push_back(DIt.e);
				ODT.Ts.push_back(DIt.T);
				ODT.Tcs.push_back(Tc_proposed);
				ODT.periods.push_back(DIt.period);
				ODT.inclination_degs.push_back(DIt.inclination_deg);
				ODT.argPeri_degs.push_back(argPeri_deg_proposed);
				// store outputs
				ODT.chiSquareds.push_back(TOTAL_chiSquared);
				if (a_total_curr>1e4)
				{
					cout<<"\n\n!!!!! a_total_curr>1e4 in m-h passed section!!!!\n\n"<<endl;
					a_total_curr=0;
				}
				ODT.a_totals.push_back(a_total_curr);
				ODT.Ks.push_back(K_proposed);
				ODT.RVoffsets.push_back(RVoffsets_proposed);
				ODT.timesBeenHeres.push_back(timesBeenHere);
				//reset timesBeenHere counter
				timesBeenHere = 1;
				//Replace 'latest' values
				inclination_deg_latest = DIt.inclination_deg;
				longAN_deg_latest = DIt.longAN_deg;
				argPeri_deg_latest = argPeri_deg_proposed;
				e_latest = DIt.e;
				period_latest = DIt.period; //  [yrs]
				T_latest = DIt.T;
				Tc_latest = Tc_proposed;
				K_latest = K_proposed;
				RVoffsets_latest = RVoffsets_proposed;
				a_total_latest = a_total_proposed;
				sqrtESinomega_latest = sqrtESinomega_proposed;
				sqrtECosomega_latest = sqrtECosomega_proposed;


			}// Done storing accepted orbit parameters
			else
			{
				if (false)//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				{//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					cout<<"I WAS NOT ACCEPTED!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					if (true)													//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					{															//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
						cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
						cout<< "chiSquare_latest = "<<chiSquare_latest <<endl;		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
						cout<< "TOTAL_chiSquared = "<<TOTAL_chiSquared <<endl;		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
						cout<<"priors_ratio = "<<priors_ratio <<endl;				//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
						cout<<"e_prior = "<<e_prior <<endl;							//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
						cout<<"inc_prior = "<< inc_prior<<endl;						//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
						cout<<"P_prior = "<< P_prior<<endl;							//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
						cout<<"likelihood_ratio = "<<likelihood_ratio <<endl;		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
						cout<<"alpha = "<< alpha <<endl;							//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
						cout<<"RHS = "<< RHS<<endl;									//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
						cout<<"accepted = "<<accepted<<endl;						//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
						cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					}																//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					cout<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				}//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				timesBeenHere+=1;
				accepted = "false";
				if (paramsVariedRecentlyAry.size()<acceptCalcTime)
				{
					paramsVariedRecentlyAry.push_back(paramBeingVaried);
					acceptedIntsRecentlyAry.push_back(0);
				}
				if (false)
				{
					cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
					cout<< "chiSquare_latest = "<<chiSquare_latest <<endl;
					cout<< "TOTAL_chiSquared = "<<TOTAL_chiSquared <<endl;
					cout<<"priors_ratio = "<<priors_ratio <<endl;
					cout<<"likelihood_ratio = "<<likelihood_ratio <<endl;
					cout<<"prior_likelihood_ratio = "<<prior_likelihood_ratio <<endl;
					cout<<"alpha = "<< alpha <<endl;
					cout<<"RHS = "<< RHS<<endl;
					cout<<"accepted = "<<accepted<<endl;
					cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
				}
			}
		}//end of ALLpassed block
		else
		{
			timesBeenHere+=1;
			accepted = "false";
			if (paramsVariedRecentlyAry.size()<acceptCalcTime)
			{
			paramsVariedRecentlyAry.push_back(paramBeingVaried);
			acceptedIntsRecentlyAry.push_back(0);
			}
		}

		//**********************************
		//Choose a random param to vary next
		//**********************************
		int randInt = RanGen.IRandomX(0,numParams-1);
		paramBeingVaried = paramsToVaryIntsAry[randInt];
		if (false)//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			cout<<"\n randInt = "<<randInt<<"-> paramBeingVaried = "<<paramBeingVaried<<endl;

		//time to update sigma and acceptance rate?
		if (samplesTillAcceptRateCalc==acceptCalcTime)
		{
			if (true)
				ss<<"\nCalculating sigmas because samplesTillAcceptRateCalc = "<<samplesTillAcceptRateCalc<<endl;
			samplesTillAcceptRateCalc=0;
			for (int i=0;i<paramsToVaryIntsAry.size();i++)
			{
				if (true)
					ss<<"calculating acceptance rate for param at location "<<i <<endl;
				int paramInt = paramsToVaryIntsAry[i];
				double totalVaried = 0;
				double totalAccepted = 0;
				if (true)
					ss<<"size paramsVariedRecentlyAry = "<<paramsVariedRecentlyAry.size() <<", size acceptedIntsRecentlyAry = "<<acceptedIntsRecentlyAry.size() <<endl;

				for (int j=0; j<paramsVariedRecentlyAry.size();j++)
				{
					if (paramsVariedRecentlyAry[j]==paramInt)
					{
						totalVaried = totalVaried+1.0;
						totalAccepted = totalAccepted+double(acceptedIntsRecentlyAry[j]);
					}
				}
				// calculate latest acceptance ratio for param
				latestAcceptRate = totalAccepted/totalVaried;
				ss<<"For paramInt "<<paramInt<<", (totalAccepted) "<< totalAccepted<<"/"<< totalVaried<< " (totalVaried) = latestAcceptRate = "<< latestAcceptRate<<endl;
			}
			// Reset arys for calculating acceptance rates
			paramsVariedRecentlyAry.clear();
			acceptedIntsRecentlyAry.clear();
			// load up acceptString and clear ss
			acceptString = ss.str();
			ss.clear();
			ss.str(std::string());
		}
		else
			samplesTillAcceptRateCalc +=1;

	}//Done sample loops

	//cout<<"\n\n FINAL SAMPLE NUMBER = "<<sample<<"\n\n"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	//Done sampling, so save last position if not done yet
	if (latestParamsSaved==false)
	{
		if (timesBeenHere>1)
		{
			// store inputs
			ODT.longAN_degs.push_back(longAN_deg_latest);
			ODT.es.push_back(e_latest);
			ODT.Ts.push_back(T_latest);
			ODT.Tcs.push_back(Tc_latest);
			ODT.periods.push_back(period_latest);
			ODT.inclination_degs.push_back(inclination_deg_latest);
			ODT.argPeri_degs.push_back(argPeri_deg_latest);
			// store outputs
			ODT.chiSquareds.push_back(chiSquare_latest);
			ODT.a_totals.push_back(ODT.a_totals.back());
			ODT.Ks.push_back(K_latest);
			ODT.RVoffsets.push_back(RVoffsets_latest);
			ODT.timesBeenHeres.push_back(timesBeenHere-1);
		}
	}

	//shrink final ODT vectors to their actual size
	ODT = odtFinish(ODT);

	//move all log prints to log string
	SSlogStr=SSlog.str();
}
