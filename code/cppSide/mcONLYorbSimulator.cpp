#include <iostream>
#include <iomanip>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <time.h>
#include "toolboxes/orbToolboxes.h" //Both the DI and RV tools headers are called in here, so it is an all in one toolbox header call
#include "rnd/kylesGaussRand.h"

using namespace std;

int main(int argc ,char *argv[])
{
	/**
	 This is the core C++ simulator to perform simple Monte Carlo runs.
	 Thanks to its simplicity, unlike MCMC and Simulated Annealing, there is no C++ manager & core function
	 structure needed.
	 @author Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
	 */

	generalTools GT;
	// print to indicate sim has started
	cout<<"\n*** C++ mcONLY simulation has started ***\n";

	std::stringstream ss;
	std::stringstream SSlog;

	// pull in simulator settings filename and output data filename
	// from command line args
	string settingsFilename;
	string outputDataFilename;
	cout<<"There were, "<<argc<<", arguments provided.  The were:"<<endl;
	if (argc>=2)
	{
		settingsFilename = argv[1];//including directory!
		cout<<"Using settings file: "<<settingsFilename<<endl;
	}
	if (argc==3)
	{
		outputDataFilename = argv[2];//including directory!
		cout<<"Output data will be written to: "<<outputDataFilename<<endl;
	}
	if (argc>3)
	{
		cout<<"mcONLYorbSimulator is only set-up to handle settings";
		cout<<"filename and output data filename command line ";
		cout<<"arguments so far."<<endl;
	}
	else if (argc==1)
	{
		settingsFilename = "/run/media/Kyle/Data1/Todai_Work/Dropbox/workspace/Binary-project/SimSettings_and_InputData/SimSettings.txt";
		cout<<"No settings filename provided, so using default: ";
		cout<<settingsFilename<<endl;
		cout<<"No output data filename provided through command line";
		cout<<" so, value within settings file or default being used."<<endl;
	}

	cout<<"Command line arguments all loaded"<<endl; //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	// instantiate the settings object and load it up
	SimSettingsObj SSO;
	SSO.settingsLoadUp(settingsFilename.c_str());
	//cout<< "Settings obj all loaded up"<<endl; //$$$$$$$$$$$$$$

	// instantiate and load up DI and RV data objects as needed
	string outSettingsDir = GT.findDirectory(settingsFilename.c_str())+"/";
	//string SystemDataFilename = SSO.settings_and_InputDataDir + SSO.SystemDataFilename;
	string SystemDataFilename = outSettingsDir + SSO.SystemDataFilename;
	cout<<"SystemDataFilename = "<< SystemDataFilename <<endl;
	DItools DIt;
	DIdataObj DIdo;
	RVdataObj RVdo;
	DataObj SYSdo;

	// log and maybe print all inputs to C++ call
	string inputsStr;
	inputsStr = ss.str();
	ss.clear();
	ss.str(std::string());
	SSlog<<inputsStr;
	if (SSO.verbose==true)
		cout<< inputsStr;

	// Find out what 'mode' simulator will run in then log and print it to screen
	string simModeStr;
	if (SSO.RVonly==false && SSO.DIonly==true)
		ss<<"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nSimulation running in DIonly mode\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	else if (SSO.RVonly==false && SSO.DIonly==false)
		ss<<"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nSimulation running in 3D mode\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	else
		ss<<"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nSimulation running in RVonly mode\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	simModeStr = ss.str();
	ss.clear();
	ss.str(std::string());
	SSlog<<simModeStr;
	cout<<simModeStr<<endl;

	cout<<"\n$$$$ about to try and load up sys data"<<endl;
	SYSdo.systemDataLoadUp(SystemDataFilename.c_str());
	cout<<"DI and RV data objects instantiated"<<endl; //$$$$$$$$$$$$$$$$$$$$$$$
	if ((SSO.DIonly==true) or (SSO.RVonly==false && SSO.DIonly==false))
	{
		//cout<<"Starting to load up DI data obj"<<endl; //$$$$$$$$$$$$$$$$$
		string DIdataFilename = outSettingsDir + SSO.DIdataFilename;
		DIdo.dataLoadUp(DIdataFilename.c_str());
		DIdo.systemDataLoadUp(SystemDataFilename.c_str());
		// instantiate DI tools obj and load it up with same values
		DIt = GT.DItoolsParamLoadUp(DIdo);
		DIt.verbose=SSO.verbose;
		//cout<<"DI data object loaded up"<<endl; //$$$$$$$$$$$$$$$$
	}
	if ((SSO.RVonly==true) or (SSO.RVonly==false && SSO.DIonly==false))
	{
		//cout<<"starting to load up RV data obj"<<endl; //$$$$$$$$$$$$$$
		string RVdataFilename = outSettingsDir + SSO.RVdataFilename;
		RVdo.dataLoadUp(RVdataFilename.c_str());
		RVdo.systemDataLoadUp(SystemDataFilename.c_str());
		//cout<<"RVdo.planet_argPeri = "<<RVdo.planet_argPeri<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		//cout<<"Done loading up RV data obj"<<endl; //$$$$$$$$$$$$$$$$$$
	}
	// find the earliest epoch for the time of last pariapsis min/max settings
	double earliestEpoch = GT.earliestEpochFinder(DIdo, RVdo);
	//cout<<"\n\nearliestEpoch = "<< earliestEpoch<<"\n\n"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	//Start seed for random number generator
	//create very high precision seed for generator (nanosecond precision)
	timespec time1;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
	int time_nsec=time1.tv_nsec;
	//cout<<"Starting time in nanoseconds for use as a random number generator seed = "<<time_nsec<<endl; //$$$$$$$$$$$$$$

	//choose generator
	CRandomSFMT1 RanGen(time_nsec);//this is the most advanced uniform random number generator, which combines SFMT and Mother-Of-All.
	StochasticLib1 RanGen2(time_nsec);


	// ******* STARTING SIMULATION!!!!  *************
	time_t startTime;
	startTime = time(NULL);

	string starterString;
	string numSamplesString =  GT.numSamplesStringMaker(SSO.numSamples);
	ss<<"\nMCONLY: $$$$$$$$$$$$$$$$$$$  Starting Simple Monte Carlo Simulator $$$$$$$$$$$$$$$$$" <<endl;
	ss<<"Number of sample orbits being created = " << numSamplesString <<endl;
	ss<<"\nRandom Generator Seed = "<<time_nsec<<"\n"<<endl;
	starterString = ss.str();
	ss.clear();
	ss.str(std::string());
	cout<< starterString;
	SSlog<<starterString;


	outputDataType ODT;
	ODT.data_filename = outputDataFilename;
	//cout<<"outputDataType instantiated"<<endl;//$$$$$$$$$$$$$$$$$$$$$

	// variables for the success rate print block in chain loop
    int printTime = SSO.numSamples/SSO.numSamplePrints;
    int printCount = 0;
    int printsDone = 0;
    int acceptedCounter = 0;
    //int numParams = 3;//there are 3 params that MUST always vary, so this is the minimum
    double a_total_curr = 0;
    //double K_p_errorPercent = 0;
    double DI_chiSquared_reduced_lowest = 1e6;
    double RV_chiSquared_reduced_lowest = 1e6;
    double TOTAL_chiSquared_reduced_lowest = 1e6;
    double DI_chiSquared_reduced = 0;
    double RV_chiSquared_reduced = 0;
    double TOTAL_chiSquared_reduced = 0;
    double RV_chiSquared_original = 0.0;
    double DI_chiSquared_original = 0.0;
    double numDIepochs = 0;
    double numRVepochs = 0;
    double one_over_nu_RV=1;
	double one_over_nu_DI=1;
	double one_over_nu_TOTAL=1;
    double chiSquareMin = SSO.chiSquaredMax;
    int bestOrbit = 0;

    //****************************************************************************************
    //Set some of the initial values and determine the  number of params to vary for DI and RV
    //****************************************************************************************
    // Determine if K will be a varied parameter
	int numDIparams=0;
	int numRVparams=0;
	int numParams=0;
    bool vary_K = true;
    double K_proposed = 0.0;
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
    	numParams+=1;
    	ss<<"\nvary_K = true"<<endl;
	}


    double inclination_deg_proposed = 0.0;
    if ((SSO.inclination_degMIN!=0)&&(SSO.inclination_degMAX!=0))
    {
    	numRVparams+=1;
    	numDIparams+=1;
    	numParams+=1;
    }
    else
	{
		if (SSO.simulate_StarPlanet==true)
			inclination_deg_proposed = RVdo.planet_inc;
		else
			inclination_deg_proposed = RVdo.star_inc;
	}
    if (SSO.RVonly==false)
    {
    	if ((SSO.longAN_degMIN!=0)&&(SSO.longAN_degMAX!=0))
    	{
    		numDIparams+=1;
    		numParams+=1;
    	}
    }
    if ((SSO.periodMIN!=0)&&(SSO.periodMAX!=0))
    {
    	numRVparams+=1;
    	numDIparams+=1;
    	numParams+=1;
    }
    double longAN_deg_proposed = 0;
	if (SSO.RVonly==false)
	{
		if (SSO.longAN_degMAX==0)
		{
			if (SSO.simulate_StarPlanet==true)
				longAN_deg_proposed = RVdo.planet_long_AN;
			else
				longAN_deg_proposed = RVdo.star_long_AN;
		}
		else
		{
			numDIparams+=1;
			numParams+=1;
			longAN_deg_proposed = RanGen.UniformRandom(SSO.longAN_degMIN, SSO.longAN_degMAX);
		}
	}
    double argPeri_deg_proposed=90;
    if (SSO.argPeri_degMAX!=0)
    {
    	numRVparams+=1;
    	numDIparams+=1;
    	numParams+=1;
    }
    double a_total_proposed = 0;
    if ((SSO.a_totalMAX!=0)&&(SSO.DIonly==true))//{NOTE: only useful for DIonly simulations as RV requires separate a1,a2,M1,M2!}
    {
		numRVparams+=1;
		numDIparams+=1;
		numParams+=1;
    }
    double period_proposed = 0;
    if (SSO.periodMAX!=0)
    {
    	numRVparams+=1;
    	numDIparams+=1;
    	numParams+=1;
    }
    else
	{
		if (SSO.simulate_StarPlanet==true)
			period_proposed = RVdo.planet_P;
		if (SSO.simulate_StarStar==true)
			period_proposed = RVdo.star_P;
	}
    double e_proposed=0;
	double sqrtESinomega_proposed;
	double sqrtECosomega_proposed;
	double sqrtEccMax = sqrt(SSO.eMAX);
	if (SSO.eMAX==0)
	{
		if (SSO.simulate_StarPlanet==true)
			e_proposed = RVdo.planet_e;
		else
			e_proposed = RVdo.star_e;
		ss<<"eMAX==0, so using value in dict: "<< e_proposed<<endl;
	}
	else
	{
		numRVparams+=1;
		numDIparams+=1;
		numParams+=1;
		if (SSO.eMAX<0.3)
			ss<<"\n\n #### eMAX<0.3, So using sqrt(e)sin(omega),sqrt(e)cos(omega) ####\n\n"<<endl;
		else
			ss<<"\n\n ### Using DIRECT e and omega ### \n\n"<<endl;
	}

	double Tmin;
    double TMIN;
    double TMAX;
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
			Tmin = earliestEpoch-SSO.periodMAX*365.242;
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
	// load initial T and Tc values from system data file
	double T_proposed = 0;
	double Tc_proposed = 0;
	if (SSO.simulate_StarPlanet==true)
	{
		Tc_proposed = SYSdo.planet_Tc;
		T_proposed = SYSdo.planet_T;
	}
	else
	{
		Tc_proposed = SYSdo.star_Tc;
		T_proposed = SYSdo.star_T;
	}
	//replace Tc value?
	if (SSO.TcStepping)
	{
		ss<<"Using TcStepping"<<endl;
		if ((SSO.T_Min!=0)&&(SSO.T_Max!=0))
		{
			numRVparams+=1;
			numDIparams+=1;
			numParams+=1;
			ss<<"Varying Tc"<<endl;
		}
		else
			ss<<"setting Tc to a constant value in system Data file = "<<Tc_proposed<<endl;
		if (T_proposed==0)
			cout<<"Setting T to 0 and will calculate it with eccArgPeri2ToTcCalc"<<endl;
		else
			ss<<"Setting T to the constant value in system Data file = "<< T_proposed<<endl;
	}
	//replace To value?
	else
	{
		ss<<"NOT using TcStepping"<<endl;
		if ((SSO.T_Min!=0)&&(SSO.T_Max!=0))
		{
			numRVparams+=1;
			numDIparams+=1;
			numParams+=1;
			ss<<"Varying T"<<endl;
		}
		else
			ss<<"Setting T to a constant"<<endl;
		//cout<<"SSO.T_Min = "<<SSO.T_Min<<",SSO.T_Min = "<<SSO.T_Min<<endl;
		if (Tc_proposed==0)
			ss<<"Setting Tc to 0 and will calculate it with eccArgPeri2ToTcCalc"<<endl;
		else
			ss<<"Setting Tc to the constant value in system Data file = "<< Tc_proposed<<endl;
	}
    if (SSO.DIonly==false)
    {
    	for (int dataset=0; dataset<RVdo.epochs_RV.size();++dataset)
    	{
    		if (SSO.RVoffsetMAXs[dataset]!=0)
    		{
    			numRVparams+=1;
    			numParams+=1;
    		}
    	}
    }
    //NOTE: early testing show Gaussian drawn values for the next 4 caused problems, so not sure if working properly yet
    double Sys_Dist_PC_proposed = SYSdo.Sys_Dist_PC;
    double Mass1_proposed = SYSdo.Mass1;
    double star_Mass2_proposed = SYSdo.star_Mass2;
    double planet_MsinI_proposed = SYSdo.planet_MsinI;

    ss<<"\nNumber of varying parameters for DI = "<< numDIparams<<", RV = " <<numRVparams <<", 3D = " <<numParams<<endl;
    ss<<"**************************************************************"<<endl;

    //load starting param notes into log and screen
    string startParamsNotes = ss.str();
	ss.clear();
	ss.str(std::string());
	cout<<startParamsNotes;
	SSlog<< startParamsNotes;

	//*****************************************************************************
    // ***** Start the samples loop *****
	//*****************************************************************************
    for ( int sample=1; sample<(SSO.numSamples+1); sample++)
    {
    	//cout<<"mcONLYorbSimulator.cpp, line# "<<375<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    	//*****************************************************************************
		// block to control printing success rate to screen
    	//*****************************************************************************
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
			ss <<"\n"<< acceptedCounter<<"/"<<sample<<" Successful. "<<printsDone<<"/"<<SSO.numSamplePrints<<" completed at ";
			ss << asctime (timeinfo);
			ss << "largest reduced chiSquared total allowed = "<<SSO.chiSquaredMax<<endl;
			ss << "latest reduced chiSquareds: DI = "<< DI_chiSquared_reduced<<", RV = "<<RV_chiSquared_reduced <<", Total = "<< TOTAL_chiSquared_reduced<<endl;
			ss << "lowest reduced chiSquare so far = "<< chiSquareMin <<endl;
			printLine = ss.str();
			ss.clear();
			ss.str(std::string());
			cout<<printLine;
			SSlog<< printLine;
        }
        //cout<<"mcONLYorbSimulator.cpp, line# "<<405<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        //*****************************************************************************
        // Generate random numbers in the required ranges for the inputs to the orbCalc
        //*****************************************************************************
        if ((SSO.inclination_degMIN!=0)&&(SSO.inclination_degMAX!=0))
        	inclination_deg_proposed = RanGen.UniformRandom(SSO.inclination_degMIN, SSO.inclination_degMAX);
        	DIt.inclination_deg = inclination_deg_proposed;
        if ((SSO.longAN_degMIN!=0)&&(SSO.longAN_degMAX!=0))
        	longAN_deg_proposed = RanGen.UniformRandom(SSO.longAN_degMIN, SSO.longAN_degMAX);
        DIt.longAN_deg = longAN_deg_proposed;
        if (SSO.DIonly==true)
        {
			if ((SSO.a_totalMAX!=0)&&(SSO.DIonly==true))
				a_total_proposed = RanGen.UniformRandom(SSO.a_totalMIN, SSO.a_totalMAX);
        }
		DIt.a_total = a_total_proposed;
		//cout<<"mcONLYorbSimulator.cpp, line# "<<421<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        //Take care of proposals for all cases of eMAX
        if ((SSO.eMAX<0.3)&&(SSO.eMAX!=0))
        {
        	sqrtESinomega_proposed = RanGen.UniformRandom(-sqrtEccMax,sqrtEccMax);
        	sqrtECosomega_proposed = RanGen.UniformRandom(-sqrtEccMax,sqrtEccMax);
			//convert those to e and argPeri
        	e_proposed = sqrtESinomega_proposed*sqrtESinomega_proposed+sqrtECosomega_proposed*sqrtECosomega_proposed;
			argPeri_deg_proposed = (180.0/PI)*atan2(sqrtESinomega_proposed,sqrtECosomega_proposed);
			if (SSO.argPeri_degMAX>180)
			{
				if (argPeri_deg_proposed<0)
					argPeri_deg_proposed = argPeri_deg_proposed+360.0;
			}
			//fix for when arg peri is -90, which is not what we want.
			if (argPeri_deg_proposed==-90.0)
				argPeri_deg_proposed = 90.0;
        }
        else
        {
        	if (SSO.eMAX!=0)
        		e_proposed = RanGen.UniformRandom(SSO.eMIN, SSO.eMAX);
        	if (SSO.argPeri_degMAX!=0)
        		argPeri_deg_proposed = RanGen.UniformRandom(SSO.argPeri_degMIN, SSO.argPeri_degMAX);
        }
        // Forcing to 90 if not varying it, indicated mostly with MIN && MAX==0
		if ((SSO.argPeri_degMAX==0)&&(SSO.argPeri_degMIN==0))
			argPeri_deg_proposed = 90.0;
        //cout<<"mcONLYorbSimulator.cpp, line# "<<443<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        DIt.e = e_proposed;
        DIt.argPeri_deg = argPeri_deg_proposed;
        if ((SSO.periodMIN!=0)&&(SSO.periodMAX!=0))
        	period_proposed = RanGen.UniformRandom(SSO.periodMIN, SSO.periodMAX); //  [yrs]
        DIt.period = period_proposed;
        Tmin = earliestEpoch-period_proposed*365.242;
       // cout<<"mcONLYorbSimulator.cpp, line# "<<450<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        if (Tmin<TMIN)
        	Tmin = TMIN;
        if (SSO.TcStepping)
			Tc_proposed = RanGen.UniformRandom(Tmin, TMAX);// thus between a full period before first observation and the time of first observation
        else
        	T_proposed = RanGen.UniformRandom(Tmin, TMAX); // thus between a full period before first observation and the time of first observation
        //Calculate the T from Tc or vice versa
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
        DIt.T = T_proposed;
        //cout<<"mcONLYorbSimulator.cpp, line# "<<454<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        //reset RVoffsets_proposed vector and get current param vals for K and offsets.
        vector<double> RVoffsets_proposed;
        if (SSO.DIonly==false)
        {
        	//cout<<"mcONLYorbSimulator.cpp, line# "<<459<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        	if (vary_K)
        		K_proposed = RanGen.UniformRandom(SSO.K_MIN,SSO.K_MAX);
        	//cout<<"mcONLYorbSimulator.cpp, line# "<<462<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        	for (int dataset=0; dataset<RVdo.epochs_RV.size();++dataset)
        	{
        		double offset_proposed=0;
        		if (SSO.RVoffsetMAXs[dataset]!=0)
        		{
					//cout<<"mcONLYorbSimulator.cpp, line# "<<465<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					//cout<< SSO.RVoffsetMINs[dataset]<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					//cout<< SSO.RVoffsetMAXs[dataset]<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					//cout<< <<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
					offset_proposed = RanGen.UniformRandom(SSO.RVoffsetMINs[dataset],SSO.RVoffsetMAXs[dataset]);
					//cout<<"mcONLYorbSimulator.cpp, line# "<<467<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        		}
        		RVoffsets_proposed.push_back(offset_proposed);
        	}
		}
        //cout<<"mcONLYorbSimulator.cpp, line# "<<469<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        //DIt.inclination_deg = 1;
		//DIt.longAN_deg =  1;
		//DIt.argPeri_deg = 162.204772;
		//DIt.e = 0.016822;
		//DIt.T = 2455652.42247;
		//		DIt.period = 0.058088 ;
		//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//        //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//        DIt.period = (21.2165298/365.242);
//        DIt.T = 2454777.94761;
//        //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        if ( SSO.silent==false )
        	cout<<"ALL random numbers loaded"<<endl;

        // Generate Gaussian values for the sys dist and masses
		if (false)//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		{
			Sys_Dist_PC_proposed = RanGen2.NormalTrunc(SYSdo.Sys_Dist_PC,0.5*SYSdo.Sys_Dist_PC_error,SYSdo.Sys_Dist_PC_error);
			Mass1_proposed = RanGen2.NormalTrunc(SYSdo.Mass1,0.5*SYSdo.Mass1_error,3.0*SYSdo.Mass1_error);
			// load up mass2 with correct value depending on star or planet companion
			if (SSO.simulate_StarPlanet==false)
				star_Mass2_proposed = RanGen2.NormalTrunc(SYSdo.star_Mass2,0.5*SYSdo.star_Mass2_error,SYSdo.star_Mass2_error);
			else
				planet_MsinI_proposed = RanGen2.NormalTrunc(SYSdo.planet_MsinI,0.5*SYSdo.planet_MsinI_error,SYSdo.planet_MsinI_error);
		}
		//Load these up into the DIt object
		DIt.Sys_Dist_PC = Sys_Dist_PC_proposed ;
		DIt.Mass1 = Mass1_proposed ;
		if (SSO.simulate_StarStar==true)
			DIt.Mass2 =  star_Mass2_proposed;
		else
			DIt.Mass2 = planet_MsinI_proposed/sin(DIt.inclination_deg*(PI/180.0));

		//calculate the a_total using K3 in the case of RVonly and 3D simulations
		if (SSO.DIonly==false)
		{
			semiMajorType SMT_in;
			semiMajorType SMT_out;
			SMT_in.a1 = 0;
			SMT_in.a2 = 0;
			SMT_in.a_total = 0;
			SMT_in.period = DIt.period;
			SMT_in.Mass1 = DIt.Mass1;
			SMT_in.Mass2 = DIt.Mass2;
			SMT_out = GT.semiMajorConverter(SMT_in);
			DIt.a_total = SMT_out.a_total;
		}

        if ( SSO.silent==false )
        {
        	string printLine2;
        	ss<< "\ninclination_deg = " <<DIt.inclination_deg  <<"\n";
        	ss<< "longAN_deg = " << DIt.longAN_deg  <<"\n";
        	ss<<  "argPeri_deg = "<< DIt.argPeri_deg <<"\n";
        	ss<<  "e = "<< DIt.e <<"\n";
        	ss<<  "period = "<< DIt.period  <<"\n";
        	ss<<  "a_total = "<<DIt.a_total<<"\n";
        	ss<<"Sys_Dist_PC = "<< DIt.Sys_Dist_PC<<"\n";
        	ss<<"Mass1 = "<<DIt.Mass1 <<"\n";
        	ss<<"Mass2 = "<< DIt.Mass2<<"\n";
        	ss<<"K = "<<K_proposed<<"\n";
        	ss<<"Tc = "<<Tc_proposed<<"\n";
        	ss<<  "T = "<< DIt.T  <<"\n\n";
        	printLine2 = ss.str();
			ss.clear();
			ss.str(std::string());
			cout<<printLine2;
        }

        double DI_chiSquared_original = 0;
        double RV_chiSquared_original = 0;
        double numDIepochsInternal = 0;
        double numRVepochsInternal = 0;

        multiEpochOrbCalcReturnType MEOCRT;
        //*****************************************************************************
        // Pass proposed values into models and calc resultant chiSquareds
        //*****************************************************************************
        if (SSO.RVonly==false)
        {
        	//cout<<"In DI block"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        	//*****************************************************************************
        	// #### DO STUFF FOR DI ORBIT CALCNS #####
        	//*****************************************************************************
        	if ( SSO.silent==false )
        		cout<<"Calculating DI orbit for this round "<<endl;

        	// Call the orbCalc to have it apply the model to the inputs and produce outputs
			MEOCRT = DIt.multiEpochOrbCalc();
			//a_total_curr = MEOCRT.a_total;

			// Calculate the reduced chiSquared from the returned chiSquared
			DI_chiSquared_original = MEOCRT.chi_squared_total;
			if (one_over_nu_DI==1)
			{
				numDIepochs = DIdo.numEpochs_DI;
				one_over_nu_DI = (1.0/((2.0*numDIepochs)-numDIparams));
			}
			DI_chiSquared_reduced = one_over_nu_DI*DI_chiSquared_original;
			//SSO.silent=false;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			if ( SSO.silent==false )
			{
				cout<<"DI_chiSquared_original = "<<DI_chiSquared_original<<endl;
				cout<<"one_over_nu_DI = "<<one_over_nu_DI<<endl;
				cout<<"DI_chiSquared_reduced = "<<DI_chiSquared_reduced<<endl;
			}
			//SSO.silent=true;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			// update lowest DI reduced chiSquared if current one is lower
			if ( DI_chiSquared_reduced<DI_chiSquared_reduced_lowest )
				DI_chiSquared_reduced_lowest = DI_chiSquared_reduced;
        }
        else
        {
        	//cout<<"In DI else block"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        	DI_chiSquared_reduced_lowest = 0;
        	numDIepochsInternal=0;
        }

        if (SSO.DIonly==false)
        {
        	//cout<<"In RV block"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        	if ( SSO.silent==false )
        		cout<<"Calculating RV residuals for this round"<<endl;
        	RV_chiSquared_original = 0.0;
        	//*****************************************************************************
        	// ##### DO STUFF FOR RV CALCNS ######
        	//*****************************************************************************
        	vector<vector<double> > VRp_vector2;
        	vector<vector<double> > VRs_vector2;

        	// load up params drawn from fixed gaussians
			RVdo.Sys_Dist_PC = Sys_Dist_PC_proposed ;
			RVdo.Mass1 = Mass1_proposed ;

        	// Load up ss or sp parts of RVdo with current trial's
        	// param values as needed.
        	if (SSO.simulate_StarPlanet==true)
			{
        		if ( SSO.silent==false )
        			cout<<"loading up input params for star-planet rv calcs"<<endl;
        		RVdo.planet_e  = DIt.e ;
        		//cout<<"e = "<<DIt.e<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        		RVdo.planet_T  = DIt.T ;
        		RVdo.planet_P  = DIt.period ;
        		if (vary_K)
        			RVdo.planet_K = K_proposed;
        		RVdo.planet_MsinI  = DIt.Mass2 ;
        		RVdo.planet_argPeri  = DIt.argPeri_deg ;
        		RVdo.planet_inc = DIt.inclination_deg ;
			}
        	if (SSO.simulate_StarStar==true)
        	{
        		if ( SSO.silent==false )
        			cout<<"loading up input params for star-star rv calcs"<<endl;
        		RVdo.star_e  = DIt.e ;
        		RVdo.star_T  = DIt.T ;
        		RVdo.star_P  = DIt.period ;
        		if (vary_K)
        		  	RVdo.star_K = K_proposed;
        		RVdo.star_Mass2  = DIt.Mass2 ;
        		RVdo.star_argPeri  = DIt.argPeri_deg ;
        		RVdo.star_inc  = DIt.inclination_deg ;
        	}
        	// get residual velocities for companion planet if needed
        	if (RVdo.planet_P!=0 and RVdo.planet_e!=0 )
        	{
        		if ( SSO.silent==false )
        			cout<<"Starting to calculate residual vel for star-planet"<<endl;
        		// instantiate S-P calc object and load up its params
        		VRcalcStarPlanet VRCsp;
        		VRCsp = GT.VRcalcStarPlanetLoadUp(RVdo);
        		//K_p_errorPercent = VRCsp.K_p_error/VRCsp.K_p;
        		//cout<<"K_p_errorPercent = "<<K_p_errorPercent <<endl;
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
				//if ((SSO.simulate_StarPlanet==true)&&(SSO.DIonly==false))
				//	a_total_curr = VRCsp.a_total;
				if ( SSO.silent==false )
				        cout<<"K_p = "<<VRCsp.K_p<<endl;
        	}

        	// get residual velocities for companion star if needed
        	if (RVdo.star_P!=0 and RVdo.star_e!=0)
        	{
        		if ( SSO.silent==false )
        			cout<<"Starting to calculate residual vel for star-star"<<endl;
        		// instantiate S-S calc object and load up its params
        		VRcalcStarStar VRCss;
        		VRCss = GT.VRcalcStarStarLoadUp(RVdo);
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
        		//if ((SSO.simulate_StarStar==true)&&(SSO.DIonly==false))
        		//	a_total_curr = VRCss.a_total;
        		if ( SSO.silent==false )
        			cout<<"K_s = "<<VRCss.K_s<<endl;
        		//cout<<"a_total = "<<VRCss.a_total<<endl;
        	}

        	// go through all residual vels and calculate the RV chiSquareds
        	for (int dataset=0; dataset<RVdo.epochs_RV.size();++dataset)
        	{
        		if ( SSO.silent==false )
        			cout<<"\nStarting to calculate chiSquared from residuals for dataset# "<< dataset<<endl;
//				if (dataset==0)
//					RVoffsets_proposed[0] = 10.486022;
//				if (dataset==1)
//					RVoffsets_proposed[1] = 1.0;
				for (int epoch=0; epoch<RVdo.epochs_RV[dataset].size(); ++epoch)
				{
					double planetVR = 0;
					double companionStarVR = 0;
					// replace zeroed values above with real values for each as needed
					if (RVdo.planet_P!=0 and RVdo.planet_e!=0)
					{
						planetVR = VRp_vector2[dataset][epoch];
						//cout<<"planetVR for epoch "<<epoch <<" is "<<planetVR <<endl;//$$$$$$$$$$$$$$$$
					}
					if (RVdo.star_P!=0 and RVdo.star_e!=0)
					{
						companionStarVR = VRs_vector2[dataset][epoch];
						//cout<<"companionStarVR for epoch "<<epoch <<" is "<<companionStarVR <<endl;//$$$$$$$$$$$$$$$$
					}

					double updatedRV_inv_var = RVdo.RV_inv_var[dataset][epoch];
					//double updatedRV_inv_var= 1.0/((1.0/RVdo.RV_inv_var[dataset][epoch])+(K_p_errorPercent*planetVR)*(K_p_errorPercent*planetVR));
					if (false)
					{
						cout<< "RV_inv_var = "<<RVdo.RV_inv_var[dataset][epoch] <<",planetVR  ="<< planetVR <<endl;//<<", K_p_errorPercent = " << K_p_errorPercent <<endl;
						cout<<"updatedRV_inv_var = "<<updatedRV_inv_var <<", RV_inv_var = "<< RVdo.RV_inv_var[dataset][epoch]<<endl;
					}
					double  RV_chiSquared_cur = GT.chiSquaredCalc((RVdo.RVs[dataset][epoch]-RVoffsets_proposed[dataset]),updatedRV_inv_var,(planetVR+companionStarVR));
					RV_chiSquared_original = RV_chiSquared_original + RV_chiSquared_cur;
					if ( SSO.silent==false )
					{
						cout<<"\nWorking on epoch "<<epoch<<endl;
						//cout<<"\noffset = "<< RVoffsets_proposed[dataset]<<endl;
						cout<<"RVdo.RVs[dataset][epoch] = "<<RVdo.RVs[dataset][epoch]<<", RVoffsets_proposed = "<<RVoffsets_proposed[dataset]<<", -> ("<<RVdo.RVs[dataset][epoch]<<" - "<<RVoffsets_proposed[dataset]<<")="<<(RVdo.RVs[dataset][epoch]-RVoffsets_proposed[dataset]) <<", planetVR= "<< planetVR<<", companionStarVR= "<< companionStarVR<<endl;
						cout<<"Difference = "<<RVdo.RVs[dataset][epoch]-RVoffsets_proposed[dataset]-planetVR-companionStarVR<<endl;
						cout<<"ChiSquared for this RV is = "<<RV_chiSquared_cur<<endl;
						cout<<"Total NON-reducedChiSquared so far is = "<<RV_chiSquared_original<<endl;
					}
				}
        	}
        	// calculate reduced version of ChiSquared for printing
        	if (one_over_nu_RV==1)
			{
				numRVepochs = RVdo.numEpochs_RV;
				one_over_nu_RV = (1.0/((1.0*numRVepochs)-numRVparams));
			}
        	RV_chiSquared_reduced = one_over_nu_RV*RV_chiSquared_original;
        	if ( SSO.silent==false )
        	{
				cout<<"\nnumRVepochs = "<< numRVepochs <<endl;
				cout<<"one_over_nu_RV = "<< one_over_nu_RV <<endl;
				cout<<"RV_chiSquared_original = "<< RV_chiSquared_original<<endl;
				cout<<"RV_chiSquared_reduced = "<< RV_chiSquared_reduced <<endl;
        	}

        	// update lowest reduced RV chisquared value found if current one is lower
        	if ( RV_chiSquared_reduced<RV_chiSquared_reduced_lowest )
        		RV_chiSquared_reduced_lowest = RV_chiSquared_reduced;
        }
        else
        {
        	//cout<<"In RV else block"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        	RVoffsets_proposed.push_back(0);
        	RV_chiSquared_reduced_lowest=0;
        }
        if (SSO.DIonly==false && SSO.RVonly==false)
        {
        	//cout<<"In both false block"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        	// THUS, do both DI and RV parts
        }

        // Do TOTAL chiSquared value calcs
        double chiSquared_TOTAL_original = DI_chiSquared_original+RV_chiSquared_original;
        if (one_over_nu_TOTAL==1)
        	one_over_nu_TOTAL = 1.0/(2.0*numDIepochs+1.0*numRVepochs-numParams);
        TOTAL_chiSquared_reduced = one_over_nu_TOTAL*chiSquared_TOTAL_original;

        //SSO.silent=false;//$$$$$$$$$$$$$$$$$$
		if ( SSO.silent==false )
		{
			cout<<"\nDI_chiSquared_original = "<<DI_chiSquared_original <<endl;
			cout<<"RV_chiSquared_original = "<< RV_chiSquared_original<<endl;
			cout<<"RV_chiSquared_original = "<< RV_chiSquared_original<<endl;
			cout<<"numDIepochs = "<<numDIepochs <<endl;
			cout<<"numRVepochs = "<< numRVepochs<<endl;
			cout<<"numParams = "<< numParams<<endl;
			cout<<"chiSquared_TOTAL_original = "<< chiSquared_TOTAL_original <<endl;
			cout<<"one_over_nu_TOTAL = "<< one_over_nu_TOTAL <<endl;
			cout<<"TOTAL_chiSquared_reduced = "<< TOTAL_chiSquared_reduced <<endl;
			//cout<< "output reduced chi squared = "<< TOTAL_chiSquared_reduced <<endl;
		}
		//SSO.silent=true;//$$$$$$$$$$$$$$

		//*****************************************************************************
		// Determine if the orbit should be accepted
		//*****************************************************************************
		if ( TOTAL_chiSquared_reduced<=SSO.chiSquaredMax )
		{
			acceptedCounter +=1;

			//store location of best orbit out of all accepted
			if ( TOTAL_chiSquared_reduced<chiSquareMin)
			{
				chiSquareMin = TOTAL_chiSquared_reduced;
				bestOrbit = acceptedCounter-1;
			}
			// store inputs
			ODT.longAN_degs.push_back(DIt.longAN_deg);
			ODT.es.push_back(DIt.e);
			ODT.Ts.push_back(DIt.T);
			ODT.Tcs.push_back(Tc_proposed);
			ODT.periods.push_back(DIt.period);
			ODT.inclination_degs.push_back(DIt.inclination_deg);
			ODT.argPeri_degs.push_back(DIt.argPeri_deg);
			// store outputs
			ODT.chiSquareds.push_back(chiSquared_TOTAL_original);
			ODT.a_totals.push_back(DIt.a_total);
			ODT.Ks.push_back(K_proposed);
			ODT.RVoffsets.push_back(RVoffsets_proposed);
			ODT.timesBeenHeres.push_back(1);

		}// Done storing accepted orbit parameters
    }//Done sample loops

    //*****************************************************************************
    // Loop done, so perform wrap-up prints and write output data and log to disk
    //*****************************************************************************
    int totalAccepted = ODT.es.size();
    if (acceptedCounter!=totalAccepted)
    	cout<<"Warning: (acceptedCounter) "<<acceptedCounter<< " != "<<totalAccepted<<" (totalAccepted)"<<endl;

    ODT.numSamplesAccepted = totalAccepted;

    string printLine3="";
    // Get all best orbit values
    ss << "\n$$$$$$$$$$$$$$$ SIMULATOR COMPLETE $$$$$$$$$$$$$$$"<<endl;
    ss<< totalAccepted <<" orbits were accepted during simulation"<<endl;
    ss<< "\nBest orbit found:"<<endl;
    ss<< "chiSquaredMin = "<< chiSquareMin <<endl;
    ss<< "LongAN = "<< ODT.longAN_degs[bestOrbit] <<endl;
    ss<< "e = "<< ODT.es[bestOrbit] <<endl;
    ss<< "To = "<< fixed <<ODT.Ts[bestOrbit] <<endl;
    ss<< "period = "<< ODT.periods[bestOrbit] <<endl;
    ss<< "inclination = "<< ODT.inclination_degs[bestOrbit] <<endl;
    ss<< "argPeri = "<< ODT.argPeri_degs[bestOrbit] <<endl;
    if (SSO.DIonly==false)
    {
    	ss<< "K = "<< ODT.Ks[bestOrbit]<<endl;
		for (int set=0;set<ODT.RVoffsets[bestOrbit].size();++set)
			if (ODT.RVoffsets[bestOrbit][set]!=0)
				ss<<"RVoffset for dataset "<<set<<", was = "<< ODT.RVoffsets[bestOrbit][set]<<endl;
    }

    // put inverse variances into the log
    ss<<"\n one_over_nu values for all three chiSquare calcs:"<<endl;
    ss<<"one_over_nu_DI = "<< one_over_nu_DI <<endl;
    ss<<"one_over_nu_RV = "<< one_over_nu_RV <<endl;
    ss<<"one_over_nu_TOTAL = "<< one_over_nu_TOTAL <<endl;

	// calc how long sim took to run
    time_t endTime;
    endTime = time(NULL);
    int timeElapsed = endTime-startTime;

    if ( timeElapsed < 60.0 )
    	ss<< "\nSimulator took "<<timeElapsed << " seconds to complete"<<endl;
    else if ( timeElapsed>3600 )
    {
    	double totalTime = ((double)timeElapsed)/3600.0;
    	int hrs = (int)totalTime;
    	int min = 60*(totalTime-(double)hrs);
    	ss<< "\nSimulator took "<<hrs << " hours and "<< min << " minutes to complete" <<endl;
    }
    else
    {
    	double totalTime = ((double)timeElapsed)/60.0;
    	int min = (int)totalTime;
    	int sec = 60*(totalTime-double(min));
    	ss <<"\nSimulator took "<<min <<" minutes and "<<sec<< " seconds to complete"<<endl;
    }
    printLine3=ss.str();
	ss.clear();
	ss.str(std::string());
	cout<<printLine3;
	// load up log with all prints in SSlog stringstream
	// then write it to file.
	SSlog<< printLine3;
	string LOGlines;
	LOGlines = SSlog.str();
	GT.logFileWriter(ODT.data_filename, LOGlines);

	// get output filename and write output data to it
    GT.fileWriter(ODT);

	return EXIT_SUCCESS;
}
