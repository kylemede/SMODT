#include <iostream>
#include <iomanip>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <time.h>
#include "toolboxes/orbToolboxes.h" //Both the DI and RV tools headers are called in here, so it is an all in one toolbox header call
#include "simAnnealOrbSimFunc.h"

using namespace std;
/*! \file */
void simAnealOrbFuncObj::simulator()
{
	/**
	This is the core function that performs the Simulated Annealing process for
	a single chain of a multi or single chain/thread/process Simulated Annealing
	simulation, or MCMC, started by the Python file BinaryOrbSimStarterDuo.py and
	controlled with the MCMC_ProcessManagerDuo.  The entire single chain process
	is controlled by the simAnnealOrbSimulator.cpp for a purely Simulated
	Annealing simulation, or MCMCorbSimulator.cpp for a MCMC one.
	@author Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
	*/

	generalTools GT;

	cout<<"\n$$$$$ inside simAnealfunc $$$$\n"<<endl;//$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$
	// ** overall tuning params **
	double sigmaPercent_min = sigmaPercent;
	double minTemp = 0.1;
	double sigmaPercent_min_simAnneal = 0.1;
	double sigmaPercent_min_mcmc = 0.0005;
	double sigmaPercent_max = 100.0;
	int dropTempCount = 0;
	double percentSimAnneal = 70.0;
	double percentSimAnnealSigmaDrop = 70.0;
	double percentMCMCSigmaDrop = 80.0;

	// variables for the success rate print block in chain loop
	int printTime = numSamples_SA/SSO.numSamplePrints;
	int printCount = 0;
	int printsDone = 0;
	std::stringstream ss;
	std::stringstream SSlog;

	int timesBeenHere = 1;
	timesBeenHereTotal = 0;
	int timesNONEpassed = 0;
	int paramBeingVaried = 2;
	//double K_p_errorPercent = 0;
	double chiSquaredMin_DI=SSO.chiSquaredMax;
	if (SSO.RVonly==true)
		chiSquaredMin_DI=0;
	double chiSquaredMin_RV=SSO.chiSquaredMax;
	if (SSO.DIonly==true)
		chiSquaredMin_RV=0;
	chiSquaredMin = 10000000;//SSO.chiSquaredMax;
	double DI_chiSquared = 0;
	double RV_chiSquared = 0;
	double numDIepochs = 0;
	double numRVepochs = 0;
	double TOTAL_chiSquared = 0;
	double TOTAL_chiSquared_reduced;
	one_over_nu_RV=1;
	one_over_nu_DI=1;
	one_over_nu_TOTAL=1;
	bestOrbit = 0;
	double inc_prior = 12345;
	double P_prior = 12345;
	double e_prior = 12345;
	double priors_ratio=12345;
	double prior_likelihood_ratio = 12345;
	double likelihood_ratio =12345;
	double alpha = 12345;
	double RHS = 0;
	double chiSquare_latest = 12345678;
	string accepted = "?";

	//start/choose generator(s)
	//this is the most advanced uniform random number generator, which combines SFMT and Mother-Of-All
	CRandomSFMT1 RanGen(randSeed);
	StochasticLib1 RanGen2(randSeed);

	// instantiate and load up a second RVdo for proposed planet vals in simStar=true cases
	RVdataObj RVdo2;
	RVdo2 = RVdo;

	//*****************************************************************************
	// set up starting values for input params
	//*****************************************************************************
	sigmaPercent_latest = sigmaPercent;
	string startParamsGenStr;
	// Determine if K will be a varied parameter
	numDIparams=0;
	numRVparams=0;
	double Kp_calculated=0;
	double Ks_calculated=0;
	vary_K = true;
	double K_latest = 0.0;
	if (SSO.DIonly==true)
	{
		vary_K = false;
		ss<<"\nvary_K = false"<<endl;
	}
	else if (SSO.inclination_degMAX!=0)
	{
		vary_K = false;
		ss<<"\nvary_K = false"<<endl;
	}
	else
	{
		numRVparams+=1;
		K_latest = RanGen.UniformRandom(SSO.K_MIN, SSO.K_MAX);
		paramsToVaryIntsAry.push_back(7);
		ss<<"\nvary_K = true"<<endl;
	}
	double inclination_deg_latest = 0;
	if (SSO.inclination_degMAX==0)
	{
		if (SSO.simulate_StarPlanet==true)
			inclination_deg_latest = RVdo.planet_inc;
		else
			inclination_deg_latest = RVdo.star_inc;
	}
	else
	{
		numRVparams+=1;
		numDIparams+=1;
		inclination_deg_latest = RanGen.UniformRandom(SSO.inclination_degMIN, SSO.inclination_degMAX);
		paramsToVaryIntsAry.push_back(4);
	}
	double longAN_deg_latest = 0;
	if (SSO.RVonly==false)
	{
		if (SSO.longAN_degMAX==0)
		{
			if (SSO.simulate_StarPlanet==true)
				longAN_deg_latest = RVdo.planet_long_AN;
			else
				longAN_deg_latest = RVdo.star_long_AN;
		}
		else
		{
			numDIparams+=1;
			longAN_deg_latest = RanGen.UniformRandom(SSO.longAN_degMIN, SSO.longAN_degMAX);
			paramsToVaryIntsAry.push_back(0);
		}
	}
	double argPeri_deg_latest=90;
	if (SSO.argPeri_degMAX==0)
	{
		if (SSO.simulate_StarPlanet==true)
		{
			if (SSO.DIonly==true)
				argPeri_deg_latest = DIdo.planet_argPeri;
			else
				argPeri_deg_latest = RVdo.planet_argPeri;
		}
		else
		{
			if (SSO.DIonly==true)
				argPeri_deg_latest = DIdo.star_argPeri;
			else
				argPeri_deg_latest = RVdo.star_argPeri;
		}
	}
	else
	{
		numRVparams+=1;
		numDIparams+=1;
		argPeri_deg_latest = RanGen.UniformRandom(SSO.argPeri_degMIN, SSO.argPeri_degMAX);
		paramsToVaryIntsAry.push_back(5);
	}
	double a_total_latest = 0;
	if ((SSO.a_totalMAX!=0)&&(SSO.DIonly==true))//{NOTE: only useful for DIonly simulations as RV requires separate a1,a2,M1,M2!}
	{
		numRVparams+=1;
		numDIparams+=1;
		a_total_latest = RanGen.UniformRandom(SSO.a_totalMIN, SSO.a_totalMAX);
		ss<<"a_totalMAX!=0, so value "<<a_total_latest<<" generated from range ["<< SSO.a_totalMIN<<" , "<< SSO.a_totalMAX<<"]"<<endl;
		paramsToVaryIntsAry.push_back(6);
	}
	double period_latest=0;
	if (SSO.periodMAX==0)
	{
		if (SSO.simulate_StarPlanet==true)
			period_latest = RVdo.planet_P;
		if (SSO.simulate_StarStar==true)
			period_latest = RVdo.star_P;
	}
	else
	{
		numRVparams+=1;
		numDIparams+=1;
		period_latest = RanGen.UniformRandom(SSO.periodMIN, SSO.periodMAX); //  [yrs]
		ss<<"\n\n**********************************\nperiodMIN = "<<SSO.periodMIN<<", periodMAX = "<<SSO.periodMAX<<"\n**********************************\n\n"<<endl;
		paramsToVaryIntsAry.push_back(3);
	}
	double e_latest=0;
	double sqrtESinomega_latest;
	double sqrtECosomega_latest;

	if (SSO.eMAX==0)
	{

		if (SSO.simulate_StarPlanet==true)
			e_latest = RVdo.planet_e;
		else
			e_latest = RVdo.star_e;
		cout<<"eMAX==0, so using value in dict: "<< e_latest<<endl;
	}
	else
	{
		numRVparams+=1;
		numDIparams+=1;
		e_latest = RanGen.UniformRandom(SSO.eMIN, SSO.eMAX);
		paramsToVaryIntsAry.push_back(1);
		if (SSO.eMAX<0.3)
		{
			ss<<"\n\n #### eMAX<0.3, So using sqrt(e)sin(omega),sqrt(e)cos(omega) ####\n\n"<<endl;
			if (SSO.argPeri_degMAX!=0)
			{
				//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				//e_latest = 0.001;
				//argPeri_deg_latest = 181;
				//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				sqrtESinomega_latest = sqrt(e_latest)*sin((PI/180.0)*argPeri_deg_latest);
				sqrtECosomega_latest = sqrt(e_latest)*cos((PI/180.0)*argPeri_deg_latest);
			}
		}
		else
			ss<<"\n\n ### Using DIRECT e and omega ### \n\n"<<endl;
	}

	double Tmin;
	TMIN;
	TMAX;
	// set up max and min values for T
	if ((SSO.T_Min==0)&&(SSO.T_Max==0))
	{
		TMIN=0;
		Tmin=0;
		TMAX=0;
	}
	else
	{
		if ((SSO.T_Min==-1)&&(SSO.T_Max==-1))
		{
			Tmin = earliestEpoch-period_latest*365.242;
			TMIN = earliestEpoch-SSO.periodMAX*365.0;
			TMAX = earliestEpoch;
			ss<<"******  both T_Min and T_Max set to -1  ******"<<endl;
		}
		else
		{
			Tmin = SSO.T_Min;
			TMIN = SSO.T_Min;
			TMAX = SSO.T_Max;
		}
	}
	ss<<fixed<<std::setprecision(6)<<"\n\nTMIN = "<<TMIN<<", TMAX = "<<TMAX<<"\n\n"<<endl;
	//cout<<"line # 247"<<endl;//$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$
	// load initial T and Tc values from system data file
	double T_latest;
	double Tc_latest;
	if (SSO.simulate_StarPlanet==true)
	{
		Tc_latest = SYSdo.planet_Tc;
		T_latest = SYSdo.planet_T;
	}
	else
	{
		Tc_latest = SYSdo.star_Tc;
		T_latest = SYSdo.star_T;
	}

	//replace Tc value?
	if (SSO.TcStepping)
	{
		ss<<"Using TcStepping"<<endl;
		if ((SSO.T_Min!=0)&&(SSO.T_Max!=0))
		{

			Tc_latest = RanGen.UniformRandom(Tmin, TMAX);
			numRVparams+=1;
			numDIparams+=1;
			paramsToVaryIntsAry.push_back(2);
			ss<<"Varying Tc"<<endl;

		}
		else
			ss<<"setting Tc to a constant"<<endl;
		if (T_latest==0)
			cout<<"Setting T to 0 and will calculate it with eccArgPeri2ToTcCalc"<<endl;
		else
			ss<<"Setting T to the constant value in system Data file = "<< T_latest<<endl;
	}
	//replace To value?
	else
	{
		ss<<"NOT using TcStepping"<<endl;
		if ((SSO.T_Min!=0)&&(SSO.T_Max!=0))
		{
			T_latest = RanGen.UniformRandom(Tmin, TMAX);
			numRVparams+=1;
			numDIparams+=1;
			paramsToVaryIntsAry.push_back(2);
			ss<<"Varying T"<<endl;
		}
		else
			ss<<"Setting T to a constant"<<endl;
		cout<<"SSO.T_Min = "<<SSO.T_Min<<",SSO.T_Min = "<<SSO.T_Min<<endl;
		if (Tc_latest==0)
			ss<<"Setting Tc to 0 and will calculate it with eccArgPeri2ToTcCalc"<<endl;
		else
			ss<<"Setting Tc to the constant value in system Data file = "<< Tc_latest<<endl;
	}
	//update non-updated T if it was 0 in the dictionary
	if ((T_latest==0)||(Tc_latest==0))
	{
		ss<<"To or Tc value being updated as it was originally zero."<<endl;
		ss<<"Initial values were: To = "<<T_latest<<", Tc = "<<Tc_latest<<endl;
		// calculate starting To from provided Tc
		eccArgPeri2ToTcType EATT;
		EATT.period = period_latest;
		EATT.argPeri_deg = argPeri_deg_latest;
		EATT.e = e_latest;
		EATT.To = T_latest;
		EATT.Tc = Tc_latest;
		EATT = GT.eccArgPeri2ToTcCalc(EATT);
		T_latest = EATT.To;
		Tc_latest = EATT.Tc;
		ss<<"Updated values are: To = "<<T_latest<<", Tc = "<<Tc_latest<<endl;

	}
	if (SSO.DIonly==true)
		Tc_latest = T_latest;

	startParamsGenStr = ss.str();
	ss.clear();
	ss.str(std::string());
	cout<<startParamsGenStr;
	SSlog<< startParamsGenStr;

	// set up starting offsets if needed
	vector<double> RVoffsets_latest;
	vector<double> RVoffsets_proposed;
	vector<double> offset_sigmaPercents_latest;
	double K_proposed;
	if (SSO.DIonly==false)
	{
		K_proposed = K_latest;
		for (int dataset=0; dataset<RVdo.epochs_RV.size();++dataset)
		{
			double offsetStart=0;
			double sig=0;
			if (SSO.RVoffsetMAXs[dataset]!=0)
			{
				numRVparams+=1;
				offsetStart = RanGen.UniformRandom(SSO.RVoffsetMINs[dataset],SSO.RVoffsetMAXs[dataset]);
				offset_sigmas.push_back((sigmaPercent_latest/100.0)*(SSO.RVoffsetMAXs[dataset]-SSO.RVoffsetMINs[dataset]));
				sig = 0.5;
				paramsToVaryIntsAry.push_back(8+dataset);
			}
			RVoffsets_latest.push_back(offsetStart);
			offset_sigmaPercents_latest.push_back(sig);
			ss<<"RV dataset "<<dataset<<", RVoffset = "<<offsetStart <<endl;
		}
		RVoffsets_proposed = RVoffsets_latest;
	}

	numParams = paramsToVaryIntsAry.size();

	double inc_sigmaPercent_latest = 10;//sigmaPercent_min;
	double longAN_sigmaPercent_latest = 1.0;//sigmaPercent_min;
	double argPeri_sigmaPercent_latest = 0.4;//1.0;
	double e_sigmaPercent_latest = 1.6;//2.0;
	double period_sigmaPercent_latest = 0.5;//sigmaPercent_min*0.5;//5.0;
	double T_sigmaPercent_latest = 0.3;//sigmaPercent_min*0.5;//1.0;
	double K_sigmaPercent_latest = 0.5;//sigmaPercent_min;//3.0;
	double a_total_sigmaPercent_latest = 0.8;//sigmaPercent_min;
	double sqrtESinomega_sigmaPercent_latest = 0.02;//0.035;//sigmaPercent_min*0.1;
	double sqrtECosomega_sigmaPercent_latest = 0.055;//0.039;//sigmaPercent_min*0.1;

	double a_total_curr=0;

	// temp and sigma tunning params
	int numSamplesSimAnneal = int(double(numSamples_SA)*(percentSimAnneal/100.0));
	int numSamplesMCMC = numSamples_SA-numSamplesSimAnneal;
	int dropTempTime = numParams*int((1.0/double(numParams))*(tempStepSizePercent/100.0)*double(numSamplesSimAnneal));;
	int switchToMCMCsample = int(double(numSamples_SA)*(percentSimAnneal/100.0));
	int numTempSteps = int(100.0/tempStepSizePercent);
	double temp = startTemp;
	double tempDropDouble = (temp+(1.0-minTemp))*(tempStepSizePercent/100.0);
	double sigmaPercentDropDouble = (sigmaPercent_min-sigmaPercent_min_simAnneal)/(double(numTempSteps)*(percentSimAnnealSigmaDrop/100.0));
	double sigmaPercentDropDouble_mcmc = (sigmaPercent_min_simAnneal-sigmaPercent_min_mcmc)/(double(numTempSteps)*((100-percentSimAnneal)/100.0)*(percentMCMCSigmaDrop/100.0));
	int tempStepNumber = 0;

	double numSaved = 0;
	bool latestParamsSaved;
	double acceptedCounter = 0;
	double samplesTillAcceptRateCalc = 0;
	double acceptCalcTime = int(double(numSamplesMCMC)/1000.0);
	if (numSamplesMCMC>100000)
	{
		if (acceptCalcTime<500)
			acceptCalcTime=500;
	}
	else
	{
		if (acceptCalcTime<100)
			acceptCalcTime=100;
	}
	double latestAcceptRate = 0;
	vector<int> paramsVariedRecentlyAry;
	vector<int> acceptedIntsRecentlyAry;
	string acceptString;

	// Printing block for loop/starting params
	string startParms;
	ss<<"\n**************************************************************"<<endl;
	ss<<"\nConducting simAnneal on "<<int(percentSimAnneal)<<"% of samples, last "<<int(100.0-percentSimAnneal)<<"% will be MCMC!\n"<<endl;
	ss<<"Important loop counters and params:"<<endl;
	ss<<"numSamples_SA = "<<numSamples_SA <<endl;
	ss<<"numSamplesSimAnneal = "<<numSamplesSimAnneal<<endl;
	ss<<"numSamplesMCMC = "<<numSamplesMCMC<<endl;
	ss<<"tempStepSizePercent = "<<tempStepSizePercent<<endl;
	ss<<"sigmaPercent_min = "<<sigmaPercent_min <<endl;
	ss<<"sigmaPercent_min_mcmc = "<< sigmaPercent_min_mcmc<<endl;
	ss<<"sigmaPercent_max = "<< sigmaPercent_max<<endl;
	ss<<"sigmaPercentDropDouble = "<<sigmaPercentDropDouble<<endl;
	ss<<"sigmaPercentDropDouble_mcmc = "<<sigmaPercentDropDouble_mcmc<<endl;
	ss<<"dropTempTime = "<<dropTempTime <<endl;
	ss<<"numTempSteps = "<< numTempSteps<<endl;
	ss<<"startTemp = "<<startTemp <<endl;
	ss<<"minTemp = "<< minTemp<<endl;
	ss<<"tempDropDouble = "<<tempDropDouble <<endl;
	ss<<"printTime = "<< printTime<<endl;
	ss<<"acceptCalcTime = "<<acceptCalcTime <<endl;
	string silentStr = GT.boolToStr(SSO.silent);
	ss<<"silent = "<< silentStr <<endl;
	string verboseStr = GT.boolToStr(SSO.verbose);
	ss<<"verbose = "<<verboseStr<<endl;

	ss<<"\nStarting values for model input parameters:"<<endl;
	ss<<"inclination_deg_latest = "<<inclination_deg_latest <<endl;
	ss<<"longAN_deg_latest = "<< longAN_deg_latest<<endl;
	ss<<"argPeri_deg_latest = "<< argPeri_deg_latest<<endl;
	ss<<"e_latest = "<<e_latest <<endl;
	ss<<"period_latest = "<< period_latest<<endl;
	if (SSO.eMAX<0.3)
	{
		if (false)
		{
			ss<<"sqrtESinomega_latest = "<< sqrtESinomega_latest<<endl;
			ss<<"sqrtECosomega_latest = "<< sqrtECosomega_latest<<endl;
		}
	}
	ss<<"Tc_latest = "<<fixed<<std::setprecision(8)<< Tc_latest<<endl;
	ss<<"T_latest = "<<fixed<<std::setprecision(8)<< T_latest<<endl;
	if (SSO.DIonly==false)
	{
		ss<<"K_latest = "<<K_latest<<endl;
		for (int dataset=0; dataset<RVdo.epochs_RV.size();++dataset)
			if (SSO.RVoffsetMAXs[dataset]!=0)
				ss<<"RV dataset "<<dataset<<", RVoffset = "<<RVoffsets_latest[dataset]<<endl;
	}
	ss<<"\nNumber of varying parameters for DI = "<< numDIparams<<", RV = " <<numRVparams <<", 3D = " <<numParams<<endl;
	ss<<"**************************************************************"<<endl;
	startParms = ss.str();
	ss.clear();
	ss.str(std::string());
	cout<<startParms;
	SSlog<< startParms;

	// setting initially 'proposed' states equal to initial 'latest' states
	double longAN_deg_proposed = longAN_deg_latest;
	double e_proposed = e_latest;
	double T_proposed=T_latest;
	double Tc_proposed = Tc_latest;
	double period_proposed = period_latest;
	double inclination_deg_proposed = inclination_deg_latest;
	double argPeri_deg_proposed = argPeri_deg_latest;
	double a_total_proposed = a_total_latest;
	double Sys_Dist_PC_proposed = SYSdo.Sys_Dist_PC;
	double Mass1_proposed = SYSdo.Mass1;
	double planet_MsinI_proposed = SYSdo.planet_MsinI;
	double star_Mass2_proposed = SYSdo.star_Mass2;
	double sqrtESinomega_proposed = sqrtESinomega_latest;
	double sqrtECosomega_proposed = sqrtECosomega_latest;

	bool ALLpassed;
	//*****************************************************************************
	// ***** Start the samples loop *****
	//*****************************************************************************
	int sample;
	for ( sample=0; sample<numSamples_SA; ++sample)
	{
		//cout<<" ---------------------  Trying sample  "<<sample<<"  ----------------"<<endl;//$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

		latestParamsSaved=false;

		// print block for end of simAnneal stage
		string endSimAnnealStr;
		if (switchToMCMCsample==sample)
		{
			ss << "\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
			ss << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
			ss << int(percentSimAnneal)<<"% Simulated Annealing stage completed at sample "<<sample<<"\n"<<endl;
			ss << int(acceptedCounter)<<"/"<<sample<<" Successful, "<<int(ODT.es.size())<<" saved. ";
			ss << "Latest acceptance rate = "<<std::setprecision(8)<<latestAcceptRate<<endl;
			ss << "sigmaPercent_min = "<<sigmaPercent_min<<" , sigmaPercent_latest = "<<sigmaPercent_latest<<endl;
			sigmaPercent_min=sigmaPercent_min_simAnneal;
			ss << "Latest param being varied = "<<paramBeingVaried<<", timesBeenHere = "<<timesBeenHere<<endl;
			ss << "Largest allowed reduced chiSquareds: DI = "<<SSO.chiSquaredMax*one_over_nu_DI <<", RV = "<<SSO.chiSquaredMax*one_over_nu_RV  <<", Total = "<<SSO.chiSquaredMax*one_over_nu_TOTAL  <<endl;
			ss << "latest reduced chiSquareds: DI = "<< DI_chiSquared*one_over_nu_DI<<", RV = "<<RV_chiSquared*one_over_nu_RV <<", Total = "<< TOTAL_chiSquared*one_over_nu_TOTAL<<endl;
			ss << "LOWEST reduced chiSquareds: DI = "<< chiSquaredMin_DI*one_over_nu_DI <<", RV = "<< chiSquaredMin_RV*one_over_nu_RV <<", Total = "<< chiSquaredMin*one_over_nu_TOTAL <<endl;
			ss << "\nLast Accepted parameters:"<<endl;
			ss << "inclination_deg_latest = "<< inclination_deg_latest<<endl;
			ss << "longAN_deg_latest = "<<longAN_deg_latest <<endl;
			ss << "argPeri_deg_latest = "<< argPeri_deg_latest<<endl;
			ss << "e_latest = "<<e_latest <<endl;
			ss << "period_latest = "<< period_latest<<endl;
			ss << "T_latest = "<<fixed<<T_latest <<endl;
			ss << "Tc_latest = "<<fixed<<Tc_latest <<endl;
			ss << "a_total_latest = "<< a_total_latest<<endl;
			if (SSO.DIonly==false)
			{
				ss << "K_latest = "<<K_latest<<endl;
				for (int dataset=0; dataset<RVoffsets_latest.size();++dataset)
				{
					if (SSO.RVoffsetMAXs[dataset]!=0)
						ss<<"RV dataset "<<dataset<<", RVoffset = "<<RVoffsets_latest[dataset] <<endl;
				}
			}
			if (false)
			{
				ss<<"\ninc_sigmaPercent_latest = "<< inc_sigmaPercent_latest<<endl;
				ss<<"longAN_sigmaPercent_latest = "<< longAN_sigmaPercent_latest<<endl;
				ss<<"period_sigmaPercent_latest = "<< period_sigmaPercent_latest<<endl;
				ss<<"argPeri_sigmaPercent_latest = "<< argPeri_sigmaPercent_latest<<endl;
				ss<<"e_sigmaPercent_latest = "<< e_sigmaPercent_latest<<endl;
				ss<<"T_sigmaPercent_latest = "<< T_sigmaPercent_latest<<endl;
				ss<<"a_total_sigmaPercent_latest = "<< a_total_sigmaPercent_latest<<endl;
				if (SSO.DIonly==false)
				{
					ss<<"K_sigmaPercent_latest = "<< K_sigmaPercent_latest<<endl;
					for (int dataset=0; dataset<RVoffsets_latest.size();++dataset)
					{
						//if (SSO.RVoffsetMAXs[dataset]!=0)
						ss<<"RV dataset "<<dataset<<", offset_sigmaPercents_latest = "<<offset_sigmaPercents_latest[dataset] <<endl;
					}
				}
			}
			if (false)
			{
				ss<<"\n one_over_nu values:"<<endl;
				ss<<"one_over_nu_DI = "<<one_over_nu_DI <<endl;
				ss<<"one_over_nu_RV = "<< one_over_nu_RV<<endl;
				ss<<"one_over_nu_TOTAL = "<< one_over_nu_TOTAL<<endl;
			}
			ss<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
			ss<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"<<endl;
			endSimAnnealStr = ss.str();
			ss.clear();
			ss.str(std::string());
			cout<<endSimAnnealStr;
			SSlog<< endSimAnnealStr;
		}

		// update temp and it's counters
		dropTempCount = dropTempCount+1;
		if (dropTempCount==dropTempTime)
		{
			//cout<<"trying to drop temp"<<endl;
			dropTempCount = 0;
			tempStepNumber = tempStepNumber+1;
			double tempIn = temp;
			temp = temp-tempDropDouble;
			if (switchToMCMCsample<=sample)
			{
				if (temp<1.0)
					temp = 1.0;

				sigmaPercent_min = sigmaPercent_min-sigmaPercentDropDouble_mcmc;
				if (sigmaPercent_min<sigmaPercent_min_mcmc)
					sigmaPercent_min = sigmaPercent_min_mcmc;
			}
			else
			{
				if (temp<minTemp)
					temp = minTemp;
			}
			//cout<<" dropping temp from "<<tempIn<<" to "<<temp<< " at sampleNumber "<<sample<<endl;//$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

			//*****************************************************************************
			//Emergency jump in param values if nothing is being accepted
			//*****************************************************************************
			if ((chiSquaredMin*one_over_nu_TOTAL)>=SSO.chiSquaredMax)
			{
				ss<<"\n****** WARNING: Chain has been stuck for a temperature step worth of trials so trying a new starting point. ******\n"<<endl;
				//cout<<"\n****** WARNING: Chain has been stuck for a temperature step worth of trials so trying a new starting point. ******\n"<<endl;
				if (SSO.inclination_degMAX!=0)
					inclination_deg_latest = RanGen.UniformRandom(SSO.inclination_degMIN, SSO.inclination_degMAX);
				if ((SSO.longAN_degMAX!=0)&&(SSO.RVonly==false))
					longAN_deg_latest = RanGen.UniformRandom(SSO.longAN_degMIN, SSO.longAN_degMAX);
				if (SSO.argPeri_degMAX!=0)
					argPeri_deg_latest = RanGen.UniformRandom(SSO.argPeri_degMIN, SSO.argPeri_degMAX);
				if (SSO.eMAX!=0)
					e_latest = RanGen.UniformRandom(SSO.eMIN, SSO.eMAX);
				if (vary_K)
					K_latest = RanGen.UniformRandom(SSO.K_MIN, SSO.K_MAX);
				if (SSO.periodMAX!=0)
					period_latest = RanGen.UniformRandom(SSO.periodMIN, SSO.periodMAX); //  [yrs]
				if (TMAX!=0)
				{
					if (SSO.TcStepping)
						Tc_latest = RanGen.UniformRandom(TMIN, TMAX);
					else
						T_latest = RanGen.UniformRandom(TMIN, TMAX);
				}
				if ((SSO.a_totalMAX!=0)&&(SSO.DIonly==true))
					a_total_latest = RanGen.UniformRandom(SSO.a_totalMIN,SSO.a_totalMAX); //[AU]
				if (((SSO.eMAX!=0)&&(SSO.argPeri_degMAX!=0))&&(SSO.eMAX<0.3))
				{
					sqrtESinomega_latest = sqrt(e_latest)*sin((PI/180.0)*argPeri_deg_latest);
					sqrtECosomega_latest = sqrt(e_latest)*cos((PI/180.0)*argPeri_deg_latest);
				}
				// set up starting offsets if needed
				if (SSO.DIonly==false)
				{
					double offsetStart;
					for (int dataset=0; dataset<RVdo.epochs_RV.size();++dataset)
					{
						if (SSO.RVoffsetMAXs[dataset]!=0)
							offsetStart = RanGen.UniformRandom(SSO.RVoffsetMINs[dataset],SSO.RVoffsetMAXs[dataset]);
						else
							offsetStart = 0;
						RVoffsets_latest[dataset] = offsetStart;
					}
					RVoffsets_proposed = RVoffsets_latest;
				}
			}
		}
		//*****************************************************************************
		// block to control printing success rate to screen
		//*****************************************************************************
		printCount = printCount + 1;
		if ( printCount==printTime )
		{
			//cout<<"Trying to print print block"<<endl;//$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$
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
			ss<<"\n############################### Simulated Annealing ######################################"<<endl;
			ss << int(acceptedCounter)<<"/"<<sample<<" Successful, "<<len<<" saved. "<<printsDone<<"/"<<SSO.numSamplePrints<<" completed at ";
			ss << asctime (timeinfo);
			ss <<"Finished "<< tempStepNumber<<"/"<< numTempSteps<< " temp steps, Current Temp = "<<temp<<", sigmaPercent_min = "<<sigmaPercent_min<<endl;
			ss << "Latest acceptance rate = "<<latestAcceptRate<<endl<<endl;
			ss << "Latest param being varied = "<<paramBeingVaried<<", timesBeenHere = "<<timesBeenHere<<endl;
			ss << "Times NONE of params passed = "<<timesNONEpassed<<endl;
			ss << "Largest allowed reduced chiSquareds: DI = "<<SSO.chiSquaredMax*one_over_nu_DI <<", RV = "<<SSO.chiSquaredMax*one_over_nu_RV  <<", Total = "<<SSO.chiSquaredMax*one_over_nu_TOTAL  <<endl;
			ss << "latest reduced chiSquareds: DI = "<< DI_chiSquared*one_over_nu_DI<<", RV = "<<RV_chiSquared*one_over_nu_RV <<", Total = "<< TOTAL_chiSquared*one_over_nu_TOTAL<<endl;
			ss << "LOWEST reduced chiSquareds: DI = "<< chiSquaredMin_DI*one_over_nu_DI <<", RV = "<< chiSquaredMin_RV*one_over_nu_RV <<", Total = "<< chiSquaredMin*one_over_nu_TOTAL <<endl;
			//cout<<"SimAnnealFunc, line #"<<679<<endl;//$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			if (false)
			{
				ss<<"-----------------------------------------------"<<endl;
				ss<< "chiSquare_latest = "<<chiSquare_latest <<endl;
				ss<< "TOTAL_chiSquared = "<<TOTAL_chiSquared <<endl;
				ss<< "TOTAL_chiSquared*one_over_nu_TOTAL = "<<TOTAL_chiSquared*one_over_nu_TOTAL <<endl;
				ss<<"likelihood_ratio = "<<likelihood_ratio <<endl;
				ss<<"alpha = "<< alpha <<endl;
				ss<<"RHS = "<< RHS<<endl;
				ss<<"accepted = "<<accepted<<endl;
				ss<<"----------------------------------------------"<<endl;
			}
			if (false)
			{
				ss<<"\ninc_sigmaPercent_latest = "<< inc_sigmaPercent_latest<<endl;
				ss<<"longAN_sigmaPercent_latest = "<< longAN_sigmaPercent_latest<<endl;
				ss<<"period_sigmaPercent_latest = "<< period_sigmaPercent_latest<<endl;
				ss<<"argPeri_sigmaPercent_latest = "<< argPeri_sigmaPercent_latest<<endl;
				ss<<"e_sigmaPercent_latest = "<< e_sigmaPercent_latest<<endl;
				ss<<"T_sigmaPercent_latest = "<< T_sigmaPercent_latest<<endl;
				ss<<"a_total_sigmaPercent_latest = "<<a_total_sigmaPercent_latest<<endl;
				if (SSO.eMAX<0.3)
				{
					if (true)
					{
						ss<<"sqrtESinomega_sigmaPercent_latest = "<<sqrtESinomega_sigmaPercent_latest <<endl;
						ss<<"sqrtECosomega_sigmaPercent_latest = "<<sqrtECosomega_sigmaPercent_latest <<endl;
					}
				}
				if (SSO.DIonly==false)
				{
					ss<<"K_sigmaPercent_latest = "<< K_sigmaPercent_latest<<endl;
					for (int dataset=0; dataset<RVoffsets_latest.size();++dataset)
					{
						//if (SSO.RVoffsetMAXs[dataset]!=0)
						ss<<"RV dataset "<<dataset<<", offset_sigmaPercents_latest = "<<offset_sigmaPercents_latest[dataset] <<endl;
					}
				}
			}
			//cout<<"line # 699"<<endl;//$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$
			if (false)
			{
				ss<<"\nLatest gaussian proposed general parameters were:"<<endl;
				ss<<"Sys_Dist_PC_proposed = "<<Sys_Dist_PC_proposed <<", peak value is = "<< SYSdo.Sys_Dist_PC<<endl;
				ss<<"Mass1_proposed = "<<Mass1_proposed <<", peak value is = "<< SYSdo.Mass1<<endl;
				ss<<"planet_MsinI_proposed = "<< planet_MsinI_proposed<<", peak value is = "<<SYSdo.planet_MsinI <<endl;
				ss<<"star_Mass2_proposed = "<<star_Mass2_proposed <<", peak value is = "<< SYSdo.star_Mass2<<endl;
			}
			if (false)
			{
				ss<<"\nLatest proposed parameters set as latest: "<<endl;
				ss<<"inclination_deg_proposed = "<< inclination_deg_proposed<<endl;
				ss<<"longAN_deg_proposed = "<< longAN_deg_proposed<<endl;
				ss<<"argPeri_deg_proposed = "<< argPeri_deg_proposed<<endl;
				ss<<"e_proposed = "<< e_proposed<<endl;
				ss<<"period_proposed = "<< period_proposed<<endl;
				ss << "a_total_latest = "<< a_total_latest<<endl;
				ss<<"T_proposed = "<<T_proposed <<endl;
				ss<<"Tc_proposed = "<<Tc_proposed <<endl;
				ss<<"K_proposed = "<<K_proposed <<endl;
				if (SSO.eMAX<0.3)
				{
					if (true)
					{
						ss<<"sqrtESinomega_proposed = "<< sqrtESinomega_proposed<<endl;
						ss<<"sqrtECosomega_proposed = "<< sqrtECosomega_proposed<<endl;
					}
				}
				for (int dataset=0; dataset<RVdo.epochs_RV.size();++dataset)
					cout<<"RVoffsets_proposed for dataset "<<dataset <<" = "<<RVoffsets_proposed[dataset] <<endl;
			}
			//cout<<"SimAnnealFunc, line #"<<751<<endl;//$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$
			if (true)
			{
				ss<<"\nLatest parameters set as latest: "<<endl;
				ss<<"inclination_deg_latest = "<< inclination_deg_latest<<endl;
				ss<<"longAN_deg_latest = "<< longAN_deg_latest<<endl;
				ss<<"argPeri_deg_latest = "<< argPeri_deg_latest<<endl;
				ss<<"e_latest = "<< e_latest<<endl;
				ss<<"period_latest = "<< period_latest<<endl;
				ss << "a_total_latest = "<< a_total_latest<<endl;
				ss<<"T_latest =  "<<T_latest <<endl;
				ss<<"Tc_latest = "<<Tc_latest <<endl;
				ss<<"K_latest = "<<K_latest <<endl;
				if ((SSO.eMAX<0.3)&&(SSO.eMAX!=0))
				{
					if (true)
					{
						ss<<"sqrtESinomega_latest = "<< sqrtESinomega_latest<<endl;
						ss<<"sqrtECosomega_latest = "<< sqrtECosomega_latest<<endl;
					}
				}
				for (int dataset=0; dataset<RVdo.epochs_RV.size();++dataset)
				{
					//if (SSO.RVoffsetMAXs[dataset]!=0)
					ss<<"RVoffsets_latest for dataset "<<dataset <<" = "<<RVoffsets_latest[dataset] <<endl;
				}
			}
			//cout<<"SimAnnealFunc, line #"<<777<<endl;//$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$
			if (acceptedCounter>0)
			{
				if (true)
				{
					//cout<<bestOrbit<<endl;
					ss<<"\nBEST parameters set so far: "<<endl;
					ss<<"inclination_deg = "<< ODT.inclination_degs[bestOrbit]<<endl;
					//cout<<"SimAnnealFunc, line #"<<783<<endl;//$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$
					ss<<"longAN_deg = "<< ODT.longAN_degs[bestOrbit]<<endl;
					ss<<"argPeri_deg = "<< ODT.argPeri_degs[bestOrbit]<<endl;
					ss<<"e = "<< ODT.es[bestOrbit]<<endl;
					ss<<"a_total = "<< ODT.a_totals[bestOrbit]<<endl;
					ss<<"period = "<< ODT.periods[bestOrbit]<<endl;
					ss<<"T =  "<<ODT.Ts[bestOrbit] <<endl;
					ss<<"Tc = "<<ODT.Tcs[bestOrbit] <<endl;
					ss<<"K = "<<ODT.Ks[bestOrbit] <<endl;
					for (int dataset=0; dataset<RVdo.epochs_RV.size();++dataset)
					{
						//if (SSO.RVoffsetMAXs[dataset]!=0)
						ss<<"RVoffsets for dataset "<<dataset <<" = "<<ODT.RVoffsets[bestOrbit][dataset] <<endl;
					}
				}
			}
			else
				ss<<"\n****\n NO SAMPLES ACCEPTED YET!!!\n****"<<endl;
			//cout<<"SimAnnealFunc, line #"<<797<<endl;//$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$
			if (true)
			{
				ss<<"\nKp_calculated = "<<Kp_calculated<<endl;
				ss<<"Ks_calculated = "<<Ks_calculated<<"\n"<<endl;
			}
			//cout<<"SimAnnealFunc, line #"<<803<<endl;//$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$
			ss<<"#######################################################################################"<<endl;
			ss<<acceptString;
			printLine = ss.str();
			ss.clear();
			ss.str(std::string());
			cout<<printLine;
			SSlog<< printLine;

		}
		//cout<<"line # 791"<<endl; //$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$
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
		//planet_MsinI_proposed = SYSdo.planet_MsinI;
		//star_Mass2_proposed = SYSdo.star_Mass2;
		sqrtESinomega_proposed = sqrtESinomega_latest;
		sqrtECosomega_proposed = sqrtECosomega_latest;

		// Generate random numbers in the required ranges for the inputs to the orbCalc
		ALLpassed = true;//just the starting value, set to false in proposal block if out of range
		string ALLpassedStr = "true";
		string ParamThatFailed = "?";
		accepted = "?";
		int dataset;
		//*****************************************************************************
		// ******* Determine which param to vary *************************
		//*****************************************************************************
		//cout<<"About to generate proposed params for sample number "<<sample<<endl;//$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		//cout<<"paramBeingVaried = "<<paramBeingVaried<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		if (paramBeingVaried==0)
		{
			// sigma is a percentage of total range, so convert its current percentage to a range value
			longAN_deg_sigma = (longAN_sigmaPercent_latest/100.0)*(SSO.longAN_degMAX-SSO.longAN_degMIN);
			longAN_deg_proposed = RanGen.UniformRandom(longAN_deg_latest-longAN_deg_sigma,longAN_deg_latest+longAN_deg_sigma);
		}
		else if (paramBeingVaried==1)
		{
			if (SSO.eMAX==0)
				cout<<"WARNING!! eMAX==0, but still trying to vary it in proposal brick!!"<<endl;
			if (SSO.eMAX<0.3)
			{
				sqrtESinomega_sigma = sqrtESinomega_sigmaPercent_latest*1.0;
				sqrtESinomega_proposed = RanGen.UniformRandom(sqrtESinomega_latest-sqrtESinomega_sigma,sqrtESinomega_latest+sqrtESinomega_sigma);
			}
			else
			{
				//sigma is a percentage of total range, so convert its current percentage to a range value
				e_sigma = (e_sigmaPercent_latest/100.0)*(SSO.eMAX-SSO.eMIN);
				e_proposed = RanGen.UniformRandom(e_latest-e_sigma,e_latest+e_sigma);
			}
		}
		else if (paramBeingVaried==2)
		{
			// sigma is a percentage of total range, so convert its current percentage to a range value
			T_sigma = (T_sigmaPercent_latest/100.0)*(TMAX-TMIN);
			//if ((SSO.periodMIN!=0)&&(SSO.periodMAX!=0))
			//	Tmin = earliestEpoch-period_proposed*365.242;
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
		{
			// sigma is a percentage of total range, so convert its current percentage to a range value
			period_sigma = (period_sigmaPercent_latest/100.0)*(SSO.periodMAX-SSO.periodMIN);
			period_proposed = RanGen.UniformRandom(period_latest-period_sigma,period_latest+period_sigma);//  [yrs]
		}
		else if (paramBeingVaried==4)
		{
			// sigma is a percentage of total range, so convert its current percentage to a range value
			inclination_deg_sigma = (inc_sigmaPercent_latest/100.0)*(SSO.inclination_degMAX-SSO.inclination_degMIN);
			inclination_deg_proposed = RanGen.UniformRandom(inclination_deg_latest-inclination_deg_sigma,inclination_deg_latest+inclination_deg_sigma);
		}
		else if (paramBeingVaried==5)
		{
			if (SSO.eMAX==0)
				cout<<"WARNING!! eMAX==0, but still trying to vary argPeri in the proposal brick!!"<<endl;
			if (SSO.eMAX<0.3)
			{
				sqrtECosomega_sigma = sqrtECosomega_sigmaPercent_latest*1.0;
				sqrtECosomega_proposed = RanGen.UniformRandom(sqrtECosomega_latest-sqrtECosomega_sigma,sqrtECosomega_latest+sqrtECosomega_sigma);
			}
			else
			{
				// sigma is a percentage of total range, so convert its current percentage to a range value
				argPeri_deg_sigma = (argPeri_sigmaPercent_latest/100.0)*(SSO.argPeri_degMAX-SSO.argPeri_degMIN);
				argPeri_deg_proposed = RanGen.UniformRandom(argPeri_deg_latest-argPeri_deg_sigma,argPeri_deg_latest+argPeri_deg_sigma);
			}
		}
		else if (paramBeingVaried==6)
		{
			// sigma is a percentage of total range, so convert its current percentage to a range value
			a_total_sigma = (a_total_sigmaPercent_latest/100.0)*(SSO.a_totalMAX-SSO.a_totalMIN);
			a_total_proposed = RanGen.UniformRandom(a_total_latest-a_total_sigma,a_total_latest+a_total_sigma);
		}
		else if (paramBeingVaried==7)
		{
			// sigma is a percentage of total range, so convert its current percentage to a range value
			K_sigma = (K_sigmaPercent_latest/100.0)*(SSO.K_MAX-SSO.K_MIN);
			K_proposed = RanGen.UniformRandom(K_latest-K_sigma,K_latest+K_sigma);
		}
		if (paramBeingVaried>7)
		{
			dataset = paramBeingVaried-8;
			if (SSO.RVoffsetMAXs[dataset]!=0)
			{
				// sigma is a percentage of total range, so convert its current percentage to a range value
				offset_sigmas[dataset] = (offset_sigmaPercents_latest[dataset]/100.0)*(SSO.RVoffsetMAXs[dataset]-SSO.RVoffsetMINs[dataset]);
				RVoffsets_proposed[dataset] = RanGen.UniformRandom(RVoffsets_latest[dataset]-offset_sigmas[dataset],RVoffsets_latest[dataset]+offset_sigmas[dataset]);
				if (false)
					ss<<"dataset = " <<dataset<<", RVoffsets_proposed[dataset] = "<<RVoffsets_proposed[dataset] <<", RVoffsetMIN = "<<SSO.RVoffsetMINs[dataset]<<", MAX = "<< SSO.RVoffsetMAXs[dataset]<<endl;
			}
		}
		//cout<<"line # 907"<<endl; //$$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//		if (true)
//		{
//			cout<<"\n\n Ecc ArgPeri Testing \n"<<endl;
//			double e_TEST = 0.5;
//			double Tc_TEST = 2455652.1;
//			vector<double> argPeris_TEST;
//			argPeris_TEST.push_back(1);
//			argPeris_TEST.push_back(89);
//			argPeris_TEST.push_back(91);
//			argPeris_TEST.push_back(179);
//			argPeris_TEST.push_back(181);
//			argPeris_TEST.push_back(269);
//			argPeris_TEST.push_back(271);
//			argPeris_TEST.push_back(359);
//			argPeris_TEST.push_back(361);
//			argPeris_TEST.push_back(449);
//			argPeris_TEST.push_back(451);
//			argPeris_TEST.push_back(539);
//			argPeris_TEST.push_back(541);
//			argPeris_TEST.push_back(-1);
//			argPeris_TEST.push_back(-89);
//			argPeris_TEST.push_back(-91);
//			argPeris_TEST.push_back(-179);
//			argPeris_TEST.push_back(-181);
//			argPeris_TEST.push_back(-269);
//			argPeris_TEST.push_back(-271);
//			argPeris_TEST.push_back(-359);
//			argPeris_TEST.push_back(-361);
//			argPeris_TEST.push_back(-449);
//			argPeris_TEST.push_back(-451);
//			argPeris_TEST.push_back(-539);
//			argPeris_TEST.push_back(-541);
//			double sqrtESinomega_TEST;
//			double sqrtECosomega_TEST;
//			double e_proposed_TEST1;
//			double argPeri_deg_proposed_TEST1;
//			double e_proposed_TEST2;
//			double argPeri_deg_proposed_TEST2;
//			double T_proposed_TEST2;
//			for (int testInt=0;testInt<argPeris_TEST.size();++testInt)
//			{
//				sqrtESinomega_TEST = sqrt(e_TEST)*sin((PI/180.0)*argPeris_TEST[testInt]);
//				sqrtECosomega_TEST = sqrt(e_TEST)*cos((PI/180.0)*argPeris_TEST[testInt]);
//				// Convert proposed values to eccentricity and argPeri
//				eccArgPeriCalcType EACT;
//				EACT.sqrtEsinArgPeri = sqrtESinomega_TEST;
//				EACT.sqrtEcosArgPeri = sqrtECosomega_TEST;
//				EACT = GT.eccArgPeriCalc(EACT);
//				e_proposed_TEST1 = EACT.e;
//				argPeri_deg_proposed_TEST1 = EACT.argPeri_deg;
//
//				e_proposed_TEST2 = sqrtESinomega_TEST*sqrtESinomega_TEST+sqrtECosomega_TEST*sqrtECosomega_TEST;
//				argPeri_deg_proposed_TEST2 = (180.0/PI)*atan2(sqrtESinomega_TEST,sqrtECosomega_TEST);
//				if (SSO.argPeri_degMAX>180)
//				{
//					if (argPeri_deg_proposed_TEST2<0)
//						argPeri_deg_proposed_TEST2 = argPeri_deg_proposed_TEST2+360.0;
//				}
//
//				cout<<"\n***************************  "<<testInt<<"  *****************************"<<endl;
//				eccArgPeri2ToTcType EATT;
//				EATT.period = period_proposed;
//				EATT.argPeri_deg = argPeri_deg_proposed_TEST2;
//				EATT.Tc = Tc_TEST;
//				EATT.e = e_proposed_TEST2;
//				EATT = GT.eccArgPeri2ToTcCalc(EATT);
//				T_proposed_TEST2 = EATT.To;
//				if (true)
//				{
//				cout<<"e_TEST = "<<e_TEST <<endl;
//				cout<<"argPeris_TEST[testInt] = "<<argPeris_TEST[testInt] <<endl;
//				cout<<"sqrtESinomega_TEST = "<<sqrtESinomega_TEST <<endl;
//				cout<<"sqrtECosomega_TEST = "<< sqrtECosomega_TEST<<endl;
//				cout<<"e_proposed_TEST1 = "<<e_proposed_TEST1 <<endl;
//				cout<<"argPeri_deg_proposed_TEST1 = "<<argPeri_deg_proposed_TEST1 <<endl;
//				cout<<"e_proposed_TEST2 = "<< e_proposed_TEST2<<endl;
//				cout<<"argPeri_deg_proposed_TEST2 = "<< argPeri_deg_proposed_TEST2<<endl;
//				cout<<"e_proposed_TEST2-e_proposed_TEST1 = "<< e_proposed_TEST2-e_proposed_TEST1<<endl;
//				cout<<"argPeri_deg_proposed_TEST2-argPeri_deg_proposed_TEST1 = "<<argPeri_deg_proposed_TEST2-argPeri_deg_proposed_TEST1 <<endl;
//				cout<<"Tc_proposed = "<<Tc_proposed <<endl;
//				cout<<"T_proposed_TEST2 = "<<T_proposed_TEST2 <<endl;
//				}
//			}
//			break;
//
//		}
//		break;

		// Convert proposed values to eccentricity and argPeri
		//eccArgPeriCalcType EACT;
		//EACT.sqrtEsinArgPeri = sqrtESinomega_proposed;
		//EACT.sqrtEcosArgPeri = sqrtECosomega_proposed;
		//EACT = GT.eccArgPeriCalc(EACT);
		//e_proposed = EACT.e;
		//argPeri_deg_proposed = EACT.argPeri_deg;

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
			//fix for when arg peri is negative, which is not what we want.
			if (argPeri_deg_proposed<0.0)
				argPeri_deg_proposed = argPeri_deg_proposed+180.0;
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
				EATT.Tc= Tc_proposed;
				EATT = GT.eccArgPeri2ToTcCalc(EATT);
				//cout<<"T in = "<<T_proposed <<", Tc in = "<< Tc_proposed<<", T out = "<<EATT.To <<", Tc out = "<< EATT.Tc<<endl;
				T_proposed = EATT.To;
				Tc_proposed = EATT.Tc;
			}
		}
		else
			Tc_proposed = T_proposed;
		//cout<<"line # 1066"<<endl; //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//		inclination_deg_proposed =  87.89;
//		longAN_deg_proposed =  142.156752;
//		argPeri_deg_proposed = 121.9;
//		e_proposed = 0.6819;
//		period_proposed = 56.1016336126;
		//Tc_proposed = 2454757.00787;//2453738.33274;
		//Tc_proposed = 2454756.73134;
//		Tc_proposed = 0;
//		//T_proposed=2453738.46731;
//		EATT.period = period_proposed;
//		EATT.argPeri_deg = argPeri_deg_proposed;
//		EATT.Tc = Tc_proposed;
//		EATT.e = e_proposed;
//		EATT.To=T_proposed;
//		EATT = GT.eccArgPeri2ToTcCalc(EATT);
//		T_proposed = EATT.To;
//		Tc_proposed = EATT.Tc;
		//T_proposed = 2454778.6475;
//		K_proposed = 279.8;
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

		//*****************************************************************************
		// Generate Gaussian values for the sys dist and masses
		//*****************************************************************************
		if (false)
		{
			Sys_Dist_PC_proposed = RanGen2.NormalTrunc(SYSdo.Sys_Dist_PC,0.5*SYSdo.Sys_Dist_PC_error,SYSdo.Sys_Dist_PC_error);
			Mass1_proposed = RanGen2.NormalTrunc(SYSdo.Mass1,0.5*SYSdo.Mass1_error,3.0*SYSdo.Mass1_error);
			// load up mass2 with correct value depending on star or planet companion
			if (SSO.simulate_StarPlanet==false)
				star_Mass2_proposed = RanGen2.NormalTrunc(SYSdo.star_Mass2,0.5*SYSdo.star_Mass2_error,SYSdo.star_Mass2_error);
			else
				planet_MsinI_proposed = RanGen2.NormalTrunc(SYSdo.planet_MsinI,0.5*SYSdo.planet_MsinI_error,SYSdo.planet_MsinI_error);
		}

		//*****************************************************************************
		//calculate the proposed a_total using K3 in the case of RVonly and 3D simulations
		//*****************************************************************************
		if (SSO.DIonly==false)
		{
			semiMajorType SMT_in;
			semiMajorType SMT_out;
			SMT_in.a1 = 0;
			SMT_in.a2 = 0;
			SMT_in.a_total = 0;
			SMT_in.period = period_proposed;
			SMT_in.Mass1 = Mass1_proposed;
			if (SSO.simulate_StarStar==true)
				SMT_in.Mass2 = star_Mass2_proposed;
			else
				SMT_in.Mass2 = planet_MsinI_proposed/sin(inclination_deg_proposed*(PI/180.0));
			SMT_out = GT.semiMajorConverter(SMT_in);
			a_total_proposed = SMT_out.a_total;
		}
		//cout<<"Done producing proposed params for sample "<<sample<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		// **** Done producing 'proposed' versions of all params being varied this round ****

		//*****************************************************************************
		// ***** Check all proposed values are good  *******
		//*****************************************************************************
		if ((SSO.longAN_degMAX!=0)&&(SSO.RVonly==false))
		{
			if ((longAN_deg_proposed>SSO.longAN_degMAX)||(longAN_deg_proposed<SSO.longAN_degMIN))
			{
				ParamThatFailed = "longAN_deg";
				ALLpassed=false;
				ALLpassedStr = "false";
			}
			//cout<<"longAN_deg_proposed checked and result was = "<< ALLpassedStr <<endl;//$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		}
		if (SSO.eMAX!=0)
		{
			if ((e_proposed>SSO.eMAX)||(e_proposed<SSO.eMIN))
			{
				ParamThatFailed = "e";
				ALLpassed=false;
				ALLpassedStr = "false";
				//cout<<"e_proposed checked and result was = "<< ALLpassedStr <<endl;//$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$
			}
		}
		if (TMAX!=0)
		{
			if (SSO.TcStepping)
			{
				if ((Tc_proposed>TMAX)||(Tc_proposed<TMIN))
				{
					ParamThatFailed = "Tc";
					ALLpassed=false;
					ALLpassedStr = "false";
					//cout<<"Tc_proposed checked and result was = "<< ALLpassedStr <<endl;//$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				}
			}
			else
			{
				if ((T_proposed>TMAX)||(T_proposed<TMIN))
				{
					ParamThatFailed = "T";
					ALLpassed=false;
					ALLpassedStr = "false";
					//cout<<"T_proposed checked and result was = "<< ALLpassedStr <<endl;//$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				}
			}
		}
		if ((SSO.periodMIN!=0)&&(SSO.periodMAX!=0))
		{
			if ((period_proposed>SSO.periodMAX)||(period_proposed<SSO.periodMIN))
			{
				ParamThatFailed = "period";
				ALLpassed=false;
				ALLpassedStr = "false";
			}
			//cout<<"period_proposed checked and result was = "<< ALLpassedStr <<endl;//$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		}
		if (SSO.inclination_degMAX!=0)
		{
			if ((inclination_deg_proposed>=SSO.inclination_degMAX)||(inclination_deg_proposed<=SSO.inclination_degMIN))
			{
				ParamThatFailed = "inclination_deg";
				ALLpassed=false;
				ALLpassedStr = "false";
			}
			//cout<<"inclination_deg_proposed checked and result was = "<< ALLpassedStr <<endl;//$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		}
		if (SSO.argPeri_degMAX!=0)
		{
			if ((argPeri_deg_proposed>SSO.argPeri_degMAX)||(argPeri_deg_proposed<SSO.argPeri_degMIN))
			{
				ParamThatFailed = "argPeri_deg";
				ALLpassed=false;
				ALLpassedStr = "false";
				//cout<<"argPeri_deg_proposed checked and result was = "<< ALLpassedStr <<endl;//$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			}
		}
		if (SSO.a_totalMAX!=0)
		{
			if ((a_total_proposed>SSO.a_totalMAX)||(a_total_proposed<SSO.a_totalMIN))
			{
				ParamThatFailed = "a_total";
				ALLpassed=false;
				ALLpassedStr = "false";
				//cout<<"a_total_proposed checked and result was = "<< ALLpassedStr <<endl;//$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			}
		}
		//cout<<"finished checking argPeri_deg_proposed"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		if (SSO.DIonly==false)
		{
			if (vary_K)
			{
				if ((K_proposed>SSO.K_MAX)||(K_proposed<SSO.K_MIN))
				{
					ParamThatFailed = "K";
					ALLpassed=false;
					ALLpassedStr = "false";
				}
				//cout<<"K_proposed checked and result was = "<< ALLpassedStr <<endl;//$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			}

			for (int dataset=0; dataset<RVoffsets_latest.size();++dataset)
			{
				if (SSO.RVoffsetMAXs[dataset]!=0)
				{
					if ((SSO.RVoffsetMINs[dataset]>RVoffsets_proposed[dataset])||(SSO.RVoffsetMAXs[dataset]<RVoffsets_proposed[dataset]))
					{
						ALLpassed = false;
						ParamThatFailed = "RVoffset";
						ALLpassedStr = "false";
						if (false)
							ss<<"dataset = " <<dataset<<", RVoffsets_proposed[dataset] = "<<RVoffsets_proposed[dataset] <<", RVoffsetMIN = "<<SSO.RVoffsetMINs[dataset]<<", MAX = "<< SSO.RVoffsetMAXs[dataset]<<endl;
					}
				}
				//cout<<"RVoffsets_proposed checked and result was = "<< ALLpassedStr <<endl;//$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			}
		}
		//cout<<"Finsihed checking all params if in range, resulting in ALLpassed = "<<ALLpassedStr<<endl;//$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$
		// **** Done checking 'proposed' versions of all params being varied this round ****


		if (ALLpassed==false)
		{
			timesNONEpassed +=1;
			if ( SSO.silent==false )
				cout<<"At least one parameter was out of range for sample draw # "<<sample<<endl;
			if (false)
				cout<<"The parameter that failed was: "<<ParamThatFailed<<endl;
		}
		//*****************************************************************************
		// if all are good, move on to calculating orbit.
		//*****************************************************************************
		if(ALLpassed)
		{
			//*****************************************************************************
			//Load up a Direct Imaging tool object for use with both DI and RV model inputs
			//*****************************************************************************
			//cout<<"ALLpassed = True"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			if (e_proposed<0)
				cout<<"e<0!!!!!!!!!!!!!!!!!, ALLpassed = "<< ALLpassedStr<<", paramBeingVaried = "<<paramBeingVaried <<endl;

			timesNONEpassed = 0;
			DIt.inclination_deg = inclination_deg_proposed;
			DIt.longAN_deg = longAN_deg_proposed;
//			if (true)//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//				DIt.argPeri_deg = argPeri_deg_proposed+180.0;//$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//			else//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$
			DIt.argPeri_deg = argPeri_deg_proposed;//$$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			DIt.e = e_proposed;
			DIt.period = period_proposed; //  [yrs]
			DIt.T = T_proposed;
			DIt.a_total = a_total_proposed;
			// load up params drawn from fixed gaussians
			DIt.Sys_Dist_PC = Sys_Dist_PC_proposed ;
			DIt.Mass1 = Mass1_proposed ;
			if (SSO.simulate_StarStar==true)
				DIt.Mass2 =  star_Mass2_proposed;
			else
				DIt.Mass2 = planet_MsinI_proposed;

			if ( SSO.silent==false )
				cout<<"ALL DI random numbers loaded"<<endl;
			//cout<<"ALL DI random numbers loaded"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

			if ( (SSO.silent==false)||(timesNONEpassed>100) )
			{
				string printLine2;
				ss<< "\ninclination_deg = " <<DIt.inclination_deg  <<endl;
				ss<< "longAN_deg = " << DIt.longAN_deg  <<endl;
				ss<<  "argPeri_deg = "<< argPeri_deg_proposed <<endl;
				ss<<  "e = "<< DIt.e <<"\n";
				ss<<  "period = "<< DIt.period <<"\n";
				ss<<"a_total = "<<DIt.a_total<<"\n";
				ss<<"Sys_Dist_PC = "<< DIt.Sys_Dist_PC<<"\n";
				ss<<"Mass1 = "<<DIt.Mass1 <<"\n";
				ss<<"Mass2 = "<< DIt.Mass2<<"\n";
				ss<<  "T = "<< DIt.T  <<endl;
				if (SSO.DIonly==false)
				{
					ss<< "K_proposed = "<<K_proposed<<endl;
					int dataset = paramBeingVaried-7;
					ss<< "RVoffsets_proposed["<< dataset<<"] = "<<RVoffsets_proposed[dataset]<<endl;
				}
				printLine2 = ss.str();
				ss.clear();
				ss.str(std::string());
				cout<<printLine2;
			}

			//*****************************************************************************
			// Calculate Direct Imaging fit if requested
			//*****************************************************************************

			if (SSO.RVonly==false)
			{
				multiEpochOrbCalcReturnType MEOCRT;
				//SSO.silent=false;//$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$
				//cout<<"In DI block"<<endl;//$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$
				// #### DO STUFF FOR DI ORBIT CALCNS #####
				if ( SSO.silent==false )
					cout<<"Calculating DI orbit for this round "<<endl;

				// Call the orbCalc to have it apply the model to the inputs and produce outputs
				MEOCRT = DIt.multiEpochOrbCalc();
				//cout<<"line # 1297"<<endl; //$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				a_total_curr = MEOCRT.a_total;
				//cout<<"line # 1300"<<endl; //$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

				// Calculate the chiSquared from the returned chiSquared
				DI_chiSquared = MEOCRT.chi_squared_total;
				//cout<<"line # 1304"<<endl; //$$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				if (one_over_nu_DI==1)
				{
					numDIepochs = DIdo.numEpochs_DI;
					one_over_nu_DI = (1.0/((2.0*numDIepochs)-numDIparams));
				}
				//SSO.silent=false;//$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$
				if ( SSO.silent==false )
				{
					cout<<"line # 1314"<<endl; //$$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					cout<<"DI_chiSquared = "<<DI_chiSquared<<endl;
					cout<<"one_over_nu_DI = "<<one_over_nu_DI<<endl;
					cout<<"DI_chiSquared_reduced = "<<DI_chiSquared*one_over_nu_DI<<endl;
				}
				//SSO.silent=true;//$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$
				// update lowest DI reduced chiSquared if current one is lower
//				if ( DI_chiSquared<chiSquaredMin_DI )
//					chiSquaredMin_DI = DI_chiSquared;
			}
			else
			{
				//cout<<"In DI else block"<<endl;//$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$
				//chiSquaredMin_DI = 0;
				//numDIepochs=0;
			}
			//cout<<"line # 1331"<<endl; //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$
			//*****************************************************************************
			// Calculate Radial Velocity fit if requested
			//*****************************************************************************
			if (SSO.DIonly==false)
			{
				RV_chiSquared = 0;

				//cout<<"In RV block"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$
				if ( SSO.silent==false )
					cout<<"Calculating RV residuals for this round"<<endl;
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
						cout<<"loading up input params for star-planet RV calcs"<<endl;
					RVdo.planet_e  = DIt.e ;
					RVdo.planet_T  = DIt.T ;
					RVdo.planet_Tc = Tc_proposed;
					RVdo.planet_P  = DIt.period ;
					if (vary_K)
						RVdo.planet_K = K_proposed;
					RVdo.planet_argPeri  = argPeri_deg_proposed ;
					RVdo.planet_inc = DIt.inclination_deg ;
					RVdo.planet_MsinI = DIt.Mass2 ;
				}
				if (SSO.simulate_StarStar==true)
				{
					if ( SSO.silent==false )
						cout<<"loading up input params for star-star RV calcs"<<endl;
					RVdo.star_e  = DIt.e ;
					RVdo.star_T  = DIt.T ;
					RVdo.star_Tc = Tc_proposed;
					RVdo.star_P  = DIt.period ;
					RVdo.star_argPeri  = argPeri_deg_proposed ;
					RVdo.star_inc  = DIt.inclination_deg ;
					RVdo.star_Mass2 =  DIt.Mass2;
				}
				// get RVs velocities for companion planet if needed
				if (RVdo.planet_P!=0 )
				{
					if ( SSO.silent==false )
						cout<<"Starting to calculate RVs for star-planet"<<endl;
					// instantiate S-P calc object and load up its params
					VRcalcStarPlanet VRCsp;
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
							VRCsp = GT.VRcalcStarPlanetLoadUp(RVdo2);
						}
					}
					else
						VRCsp = GT.VRcalcStarPlanetLoadUp(RVdo);

					//K_p_errorPercent = VRCsp.K_p_error/VRCsp.K_p;
					//cout<<"Back from VRcalcStarPlanetLoadUp"<<endl;//$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$
					//cout<<"there were "<<RVdo.epochs_RV.size()<<" datasets found in the RVdata file"<<endl;//$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$
					VRCsp.verbose = false;
					// run through all RV data sets and calc RVs for it
					for (int dataset=0; dataset<int(RVdo.epochs_RV.size());++dataset)
					{

						if ( SSO.silent==false )
							cout<<"Calculating P-S RVs for dataset "<<(dataset+1)<<"/"<<int(RVdo.epochs_RV.size())<<endl;
						VRCsp.epochs_p = RVdo.epochs_RV[dataset];
						vector<double> VRp_vector;
						//cout<<"about to call multiEpochCalc for dataset "<<dataset<<endl;//$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$
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

				// get RVs velocities for companion star if needed
				if (RVdo.star_P!=0 )
				{
					if ( SSO.silent==false )
						cout<<"Starting to calculate RVs for star-star"<<endl;
					// instantiate S-S calc object and load up its params
					VRcalcStarStar VRCss;
					//generate latest params for planet VR calcs from Gaussians
					if (false)
					{
						if (SSO.simulate_StarPlanet==true)
						{
							RVdo2.star_e = RanGen2.NormalTrunc(RVdo.star_e,RVdo.star_e_error,3.0*RVdo.star_e_error);
							RVdo2.star_T = RanGen2.NormalTrunc(RVdo.star_T,RVdo.star_T_error,3.0*RVdo.star_T_error);
							RVdo2.star_P = RanGen2.NormalTrunc(RVdo.star_P,RVdo.star_P_error,3.0*RVdo.star_P_error);
							RVdo2.star_argPeri = RanGen2.NormalTrunc(RVdo.star_argPeri,RVdo.star_argPeri_error,3.0*RVdo.star_argPeri_error);
							RVdo2.star_inc = RanGen2.NormalTrunc(RVdo.star_inc,RVdo.star_inc_error,3.0*RVdo.star_inc_error);
							//load up VRCsp obj with updated RVdo2 vals
							VRCss = GT.VRcalcStarStarLoadUp(RVdo2);
						}
					}
					else
						VRCss = GT.VRcalcStarStarLoadUp(RVdo);
					VRCss.verbose = false;
					// run through all RV data sets and calc RVs for it
					for (int dataset=0; dataset<int(RVdo.epochs_RV.size());++dataset)
					{
						if ( SSO.silent==false )
							cout<<"Calculating S-S RVs for dataset "<<(dataset+1)<<"/"<<int(RVdo.epochs_RV.size())<<endl;
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
				//****************************************************************
				// go through all RVs and calculate the residuals and RV chiSquareds
				//****************************************************************
				for (int dataset=0; dataset<RVdo.epochs_RV.size();++dataset)
				{

					if ( SSO.silent==false )
						cout<<"\nStarting to calculate chiSquared from residuals for dataset# "<< dataset<<endl;

					//RVoffsets_proposed[0]=97.7;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					//RVoffsets_proposed[1]=-0.425167;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					//RVoffsets_proposed[2]=138.0;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					for (int epoch=0; epoch<RVdo.epochs_RV[dataset].size(); ++epoch)
					{
						double planetVR = 0;
						double companionStarVR = 0;
						// replace zeroed values above with real values for each as needed
						if (RVdo.planet_P!=0)
						{
							planetVR = VRp_vector2[dataset][epoch];
							//cout<<"planetVR for epoch "<<epoch <<" is "<<planetVR <<endl;//$$$$$$$$$ DEBUGGING $$$$$$$$$$
						}
						if (RVdo.star_P!=0)
						{
							companionStarVR = VRs_vector2[dataset][epoch];
							//cout<<"companionStarVR for epoch "<<epoch <<" is "<<companionStarVR <<endl;//$$$$$$$$$$ DEBUGGING $$$$$$$$$
						}

						//double updatedRV_inv_var = RVdo.RV_inv_var[dataset][epoch];
						//double updatedRV_inv_var = 1.0/((1.0/RVdo.RV_inv_var[dataset][epoch])+(K_p_errorPercent*planetVR)*(K_p_errorPercent*planetVR));
						//if (false)
						//{
						//	cout<< "RV_inv_var = "<<RVdo.RV_inv_var[dataset][epoch] <<",planetVR  ="<< planetVR <<", K_p_errorPercent = " << K_p_errorPercent <<endl;
						//	cout<<"updatedRV_inv_var = "<<updatedRV_inv_var <<", RV_inv_var = "<< RVdo.RV_inv_var[dataset][epoch]<<endl;
						//}
						//double  RV_chiSquared_cur = GT.chiSquaredCalc((RVdo.RVs[dataset][epoch]-RVoffsets_proposed[dataset]),updatedRV_inv_var,(planetVR+companionStarVR));
						double  RV_chiSquared_cur = GT.chiSquaredCalc((RVdo.RVs[dataset][epoch]-RVoffsets_proposed[dataset]),RVdo.RV_inv_var[dataset][epoch],(planetVR+companionStarVR));
						RV_chiSquared = RV_chiSquared + RV_chiSquared_cur;
						if ( SSO.silent==false )
						{
							cout<<"\noffset = "<< RVoffsets_proposed[dataset]<<endl;
							cout<<"RVdo.RVs[dataset][epoch] = "<<RVdo.RVs[dataset][epoch]<<", ("<<RVdo.RVs[dataset][epoch]<<" - "<<RVoffsets_proposed[dataset]<<")="<<(RVdo.RVs[dataset][epoch]-RVoffsets_proposed[dataset]) <<", planetVR= "<< planetVR<<", companionStarVR= "<< companionStarVR<<endl;
							cout<<"Difference = "<<RVdo.RVs[dataset][epoch]-RVoffsets_proposed[dataset]-planetVR<<endl;
							cout<<"ChiSquared for this RV is = "<<RV_chiSquared_cur<<endl;
							cout<<"Total NON-reducedChiSquared so far is = "<<RV_chiSquared<<endl;
						}
					}//End epoch loop
				}//End dataset loop

				// calculate reduced version of ChiSquared for printing
				if (one_over_nu_RV==1)
				{
					numRVepochs = RVdo.numEpochs_RV;
					one_over_nu_RV = (1.0/((1.0*numRVepochs)-numRVparams));
				}

				if ( SSO.silent==false )
				{
					cout<<"\nnumRVepochs = "<< numRVepochs <<endl;
					cout<<"one_over_nu = "<< one_over_nu_RV <<endl;
				}

				// update lowest RV chisquared value found if current one is lower
//				if ( RV_chiSquared<chiSquaredMin_RV )
//					chiSquaredMin_RV = RV_chiSquared;

			}//End RV calc block
			else
			{
				//So don't calc RV stuff ie. DIonly=true
				RVoffsets_latest.push_back(0);
				//chiSquaredMin_RV=0;
			}
			//cout<<"line # 1539"<<endl; //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$

			// ****  Done all DI and RV calculations so put it all together ********
			// ****  and see if this proposed orbit will be accepted        ********

			// Do TOTAL chiSquared value calcs
			TOTAL_chiSquared  = DI_chiSquared+RV_chiSquared;
			if (one_over_nu_TOTAL==1)
				one_over_nu_TOTAL = 1.0/(2.0*numDIepochs+1.0*numRVepochs-numParams);
			TOTAL_chiSquared_reduced = TOTAL_chiSquared*one_over_nu_TOTAL;

			//SSO.silent=false;//$$$$$$$$$$ DEBUGGING $$$$$$$$$$$
			if ( SSO.silent==false )
			{
				cout<<"\nDI_chiSquared = "<<DI_chiSquared <<endl;
				cout<<"RV_chiSquared = "<< RV_chiSquared<<endl;
				cout<<"numDIepochs = "<<numDIepochs <<endl;
				cout<<"numRVepochs = "<< numRVepochs<<endl;
				cout<<"numParams = "<< numParams<<endl;
				cout<<"TOTAL_chiSquared  = "<< TOTAL_chiSquared  <<endl;
				cout<<"one_over_nu_TOTAL = "<< one_over_nu_TOTAL <<endl;
				cout<<"TOTAL_chiSquared_reduced = "<<  TOTAL_chiSquared_reduced<<endl;
			}
			//SSO.silent=true;//$$$$$$$$ DEBUGGING $$$$$$$$$

			//*******************************************
			// Determine if the orbit should be accepted
			//*******************************************
			//Calculate priors ratio
			if (((period_latest*365.242)<1000.0)||(SSO.eMAX==0))
				e_prior = 1.0;
			else
				e_prior = e_latest/DIt.e;
			if ((SSO.inclination_degMIN!=0)&&(SSO.inclination_degMAX!=0))
				inc_prior = sin(DIt.inclination_deg*(PI/180.0))/sin(inclination_deg_latest*(PI/180.0));
			else
				inc_prior = 1.0;
			if (false)//((SSO.periodMIN!=0)&&(SSO.periodMAX!=0))
				P_prior = DIt.period/period_latest;
			else
				P_prior = 1.0;
			priors_ratio = P_prior*e_prior*inc_prior;
			//Calculate likelihood ratio
			likelihood_ratio = exp((chiSquare_latest - TOTAL_chiSquared )/ (2.0*temp));
			RHS = priors_ratio*likelihood_ratio;
			alpha = RanGen.UniformRandom(0.0, 1.0);
			//alpha=0.0;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

			if (false)//$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			{
				cout<<"&&&&&&&&&&&&&&&&&& BEFORE ALFA<RHS  CHECK &&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
				cout<< "chiSquare_latest = "<<chiSquare_latest <<endl;
				cout<< "TOTAL_chiSquared = "<<TOTAL_chiSquared  <<endl;
				cout<<" TOTAL_chiSquared_reduced = "<<TOTAL_chiSquared_reduced<<endl;
				cout<< "temp = "<< temp<<endl;
				cout<<"P_prior = "<<P_prior <<endl;
				cout<<"e_prior = "<<e_prior <<endl;
				cout<<"inc_prior = "<<inc_prior <<endl;
				cout<<"priors_ratio = "<<priors_ratio <<endl;
				cout<<"likelihood_ratio = "<<likelihood_ratio <<endl;
				cout<<"alpha = "<< alpha <<endl;
				cout<<"RHS = "<< RHS<<endl;
				cout<<"accepted = "<<accepted<<endl;
				cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
			}

			if ( alpha<=RHS )
			{
				//**************************************************************************************
				//proposed orbit was accepted, so load it into output array and update 'latest' values
				//**************************************************************************************
				if (alpha>RHS)
					cout<<"WARNING: alpha>RHS messup, "<< alpha<<">"<< RHS<<", but was still accepted."<<endl;

				accepted = "true";
				//if (SSO.silent==false)
				if (false)
				{
					cout<<"&&&&&&&&&&&&&&&&&& VALUES PASSED ALPHA<RHS CHECK &&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
					cout<<"sample = "<<sample<<endl;
					cout<< "chiSquare_latest = "<<chiSquare_latest <<endl;
					cout<< "TOTAL_chiSquared = "<<TOTAL_chiSquared  <<endl;
					cout<<" TOTAL_chiSquared_reduced = "<<TOTAL_chiSquared_reduced<<endl;
					cout<< "temp = "<< temp<<endl;
					cout<<"P_prior = "<<P_prior <<endl;
					cout<<"e_prior = "<<e_prior <<endl;
					cout<<"inc_prior = "<<inc_prior <<endl;
					cout<<"priors_ratio = "<<priors_ratio <<endl;
					cout<<"likelihood_ratio = "<<likelihood_ratio <<endl;
					cout<<"alpha = "<< alpha <<endl;
					cout<<"RHS = "<< RHS<<endl;
					cout<<"accepted = "<<accepted<<endl;
					cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
				}

				acceptedCounter +=1;
				chiSquare_latest = TOTAL_chiSquared ;
				if (sample>=switchToMCMCsample)
				{
					if (false)
						cout<<"latest param accepted and in MCMC section so putting "<<paramBeingVaried<<" into paramarray and 1 into accepted array"<<endl;
					if (paramsVariedRecentlyAry.size()<acceptCalcTime)
					{
						paramsVariedRecentlyAry.push_back(paramBeingVaried);
						acceptedIntsRecentlyAry.push_back(1);
					}
				}

				//store location of best orbit out of all accepted
				if ( TOTAL_chiSquared <chiSquaredMin)
				{
					chiSquaredMin = TOTAL_chiSquared;
					chiSquaredMin_DI = DI_chiSquared;
					chiSquaredMin_RV = RV_chiSquared;
					bestOrbit = numSaved;
				}

				latestParamsSaved=true;
				numSaved +=1;
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
				timesBeenHereTotal+=timesBeenHere;
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

				if (false)
				{
					cout<<"\nLatest parameters set as latest: "<<endl;
					cout<< "TOTAL_chiSquared = "<<TOTAL_chiSquared <<endl;
					cout<<"inclination_deg_latest = "<< inclination_deg_latest<<endl;
					cout<<"longAN_deg_latest = "<< longAN_deg_latest<<endl;
					cout<<"argPeri_deg_latest = "<< argPeri_deg_latest<<endl;
					cout<<"e_latest = "<< e_latest<<endl;
					cout<<"period_latest = "<< period_latest<<endl;
					cout<<"T_latest = "<<T_latest <<endl;
					cout<<"Tc_latest = "<<Tc_latest <<endl;
					cout<<"K_latest = "<<K_latest <<endl;
					for (int dataset=0; dataset<RVdo.epochs_RV.size();++dataset)
						cout<<"RVoffsets_proposed for dataset "<<dataset <<" = "<<RVoffsets_proposed[dataset] <<endl;
				}

			}// Done storing accepted orbit parameters
			else
			{
				//****************************************************************
				// alpha<=RHS not satisfied, increment timesBeenHere and try again
				//****************************************************************
				timesBeenHere+=1;
				accepted = "false";
				if (sample>=switchToMCMCsample)
				{
					if (false)
						cout<<"failed M-H block and in MCMC section so putting "<<paramBeingVaried<<" into paramarray and 0 into accepted array"<<endl;
					if (paramsVariedRecentlyAry.size()<acceptCalcTime)
					{
						paramsVariedRecentlyAry.push_back(paramBeingVaried);
						acceptedIntsRecentlyAry.push_back(0);
					}
				}
				else
				{
					if (false)
						cout<<"failed M-H block for sample number "<< sample<<endl;
				}
			}
		}//end of ALLpassed block
		else
		{
			//****************************************************************************
			// Proposed parameters did not all pass, increment timesBeenHere and try again
			//****************************************************************************
			timesBeenHere+=1;
			accepted = "false";
			if (sample>=switchToMCMCsample)
			{
				if (false)
					cout<<"latest param not accepted and in MCMC section so putting "<<paramBeingVaried<<" into paramarray and 0 into accepted array"<<endl;
				if (paramsVariedRecentlyAry.size()<acceptCalcTime)
				{
					paramsVariedRecentlyAry.push_back(paramBeingVaried);
					acceptedIntsRecentlyAry.push_back(0);
				}
				if (paramsVariedRecentlyAry.size()>acceptCalcTime)
				{
					cout<<"\nSize of vects too big:"<<endl;
					cout<<"paramsVariedRecentlyAry.size() = "<< paramsVariedRecentlyAry.size()<<endl;
					cout<<"acceptedIntsRecentlyAry.size() = "<<acceptedIntsRecentlyAry.size() <<endl;
				}
			}
			if (false)//$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$
				cout<<"No params passed for sample number "<<sample<<endl;
		}

		//**********************************
		//Choose a random param to vary next
		//**********************************
		int randInt = RanGen.IRandomX(0,numParams-1);
		paramBeingVaried = paramsToVaryIntsAry[randInt];
		if (false)
			cout<<"\n randInt = "<<randInt<<"-> paramBeingVaried = "<<paramBeingVaried<<endl;

		//***********************************************************
		// Finished Simulated Annealing, time to switch to MCMC.
		// This means fixing the temp=1 and allowing the sigma values
		// to vary and be tunned.
		//***********************************************************
		if (sample>=switchToMCMCsample)
		{
			//********************************************
			//time to update sigma and acceptance rate?
			//********************************************
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

					if (true)//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
						ss<<"For paramInt "<<paramInt<<", (totalAccepted) "<< totalAccepted<<"/"<< totalVaried<< " (totalVaried) = latestAcceptRate = "<< latestAcceptRate<<endl;

					// ******* Determine which param's sigma to vary *************************
					if (paramInt==0)
						sigmaPercent_latest = longAN_sigmaPercent_latest ;
					if (paramInt==1)
					{
						if (SSO.eMAX<0.3)
							sigmaPercent_latest =sqrtESinomega_sigmaPercent_latest;
						else
							sigmaPercent_latest = e_sigmaPercent_latest;
					}

					if (paramInt==2)
						sigmaPercent_latest = T_sigmaPercent_latest;
					if (paramInt==3)
						sigmaPercent_latest = period_sigmaPercent_latest;
					if (paramInt==4)
						sigmaPercent_latest = inc_sigmaPercent_latest;
					if (paramInt==5)
					{
						if (SSO.eMAX<0.3)
							sigmaPercent_latest =sqrtECosomega_sigmaPercent_latest;
						else
							sigmaPercent_latest = argPeri_sigmaPercent_latest;
					}
					if (paramInt==6)
						sigmaPercent_latest = a_total_sigmaPercent_latest;
					if (paramInt==7)
						sigmaPercent_latest = K_sigmaPercent_latest;
					if (paramInt>7)
						sigmaPercent_latest = offset_sigmaPercents_latest[paramInt-8];
					//***************************************************************************

					//update sigmaPercent_latest
					double upperLimit = 0.35;
					double lowerLimit = 0.25;
					if (latestAcceptRate >=upperLimit)
					{
						if (sigmaPercent_latest<sigmaPercent_max)
							sigmaPercent_latest = sigmaPercent_latest+sigmaPercent_min;
						if (true)//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
							ss<<"latest acceptance was "<<latestAcceptRate <<" which is >"<<upperLimit<<", so sigmaPercent_latest raised to "<< sigmaPercent_latest<<" for parameter # "<<paramInt<<endl;
					}
					else if (latestAcceptRate<lowerLimit)
					{
						if (sigmaPercent_latest>=(sigmaPercent_min*2.0))
							sigmaPercent_latest = sigmaPercent_latest-sigmaPercent_min;
						if (true)//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
							ss<<"latest acceptance was "<<latestAcceptRate <<" which is <"<<lowerLimit<<", so sigmaPercent_latest lowered to "<< sigmaPercent_latest<<" for parameter # "<<paramInt<<endl;
					}

					// **** Load updated sigma back into its varying param sigma name ****
					if (paramInt==0)
						longAN_sigmaPercent_latest = sigmaPercent_latest ;
					else if (paramInt==1)
					{
						if (SSO.eMAX<0.3)
							sqrtESinomega_sigmaPercent_latest = sigmaPercent_latest;
						else
							e_sigmaPercent_latest = sigmaPercent_latest;
					}
					else if (paramInt==2)
						T_sigmaPercent_latest = sigmaPercent_latest;
					else if (paramInt==3)
						period_sigmaPercent_latest = sigmaPercent_latest;
					else if (paramInt==4)
						inc_sigmaPercent_latest = sigmaPercent_latest;
					else if (paramInt==5)
					{
						if (SSO.eMAX<0.3)
							sqrtECosomega_sigmaPercent_latest = sigmaPercent_latest;
						else
							argPeri_sigmaPercent_latest = sigmaPercent_latest;
					}
					else if (paramInt==6)
						a_total_sigmaPercent_latest = sigmaPercent_latest;
					else if (paramInt==7)
						K_sigmaPercent_latest = sigmaPercent_latest;
					else if (paramInt>7)
						offset_sigmaPercents_latest[paramInt-8] = sigmaPercent_latest ;
				}//done varying calculating accept rate and varying sigma for each param

				// Reset arys for calculating acceptance rates
				paramsVariedRecentlyAry.clear();
				acceptedIntsRecentlyAry.clear();
				// load up acceptString and clear ss
				acceptString = ss.str();
				ss.clear();
				ss.str(std::string());
			}
			//not time to update sigmas or accept rate, so just increment counter
			else
				samplesTillAcceptRateCalc +=1;
		}//end of accept and sigma vary block

	}//Done sample loops

	//******************************************************
	//Done sampling, so save last position if not done yet
	//******************************************************
	if (latestParamsSaved==false)
	{
		SSlog<<"\nlatestParamsSaved==false, so storing last values at sample number "<<sample<<", the timesBeenHere = "<<timesBeenHere<<endl;
		SSlog<<"Before storing: ODT.timesBeenHeres.size() = "<<ODT.timesBeenHeres.size()<<endl;
		if (timesBeenHere>1)
		{
			SSlog<<"storing values now"<<endl;
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
			timesBeenHereTotal+=timesBeenHere-1;
			//cout<<"After storing: ODT.timesBeenHeres.size() = "<<ODT.timesBeenHeres.size()<<endl;
		}
		//cout<<"After storing: ODT.timesBeenHeres.size() = "<<ODT.timesBeenHeres.size()<<endl;
	}

	// final print to let us know it was able to get to end of file
	cout<<"\n\n FINAL SAMPLE NUMBER = "<<sample<<endl;
	SSlog<<"\n\n FINAL SAMPLE NUMBER = "<<sample<<endl;
	cout<<"Leaving SimAnnealOrbSimFunc\n\n"<<endl;
	SSlog<<"Leaving SimAnnealOrbSimFunc\n\n"<<endl;

	//move all log prints to log string
	SSlogStr=SSlog.str();
}
