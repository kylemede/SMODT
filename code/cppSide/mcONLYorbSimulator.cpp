#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <time.h>
//#include <fstream>
#include "Toolboxes/orbToolboxes.h" //Both the DI and RV tools headers are called in here, so it is an all in one toolbox header call
//#include <rnd/stocc.h> //a library with advanced non-uniform random number generators

using namespace std;

int main(int argc ,char *argv[])
{
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

	//cout<<"Command line arguments all loaded"<<endl; //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	// instantiate the settings object and load it up
	SimSettingsObj SSO;
	SSO.settingsLoadUp(settingsFilename.c_str());

	//cout<< "Settings obj all loaded up"<<endl; //$$$$$$$$$$$$$$

	// instantiate and load up DI and RV data objects as needed
	string SystemDataFilename = SSO.settings_and_InputDataDir + SSO.SystemDataFilename;
	cout<<"SystemDataFilename = "<< SystemDataFilename <<endl;
	DItools DIt;
	DIdataObj DIdo;
	RVdataObj RVdo;
	//cout<<"DI and RV data objects instantiated"<<endl; //$$$$$$$$$$$$$$$$$$$$$$$
	if ((SSO.DIonly==true) or (SSO.RVonly==false && SSO.DIonly==false))
	{
		//cout<<"Starting to load up DI data obj"<<endl; //$$$$$$$$$$$$$$$$$
		string DIdataFilename = SSO.settings_and_InputDataDir + SSO.DIdataFilename;
		DIdo.dataLoadUp(DIdataFilename.c_str());
		DIdo.systemDataLoadUp(SystemDataFilename.c_str());
		// instantiate DI tools obj and load it up with same values
		DIt = DItoolsParamLoadUp( DIdo);
		DIt.verbose=SSO.verbose;
		//cout<<"DI data object loaded up"<<endl; //$$$$$$$$$$$$$$$$
	}
	if ((SSO.RVonly==true) or (SSO.RVonly==false && SSO.DIonly==false))
	{
		//cout<<"starting to load up RV data obj"<<endl; //$$$$$$$$$$$$$$
		string RVdataFilename = SSO.settings_and_InputDataDir + SSO.RVdataFilename;
		RVdo.dataLoadUp(RVdataFilename.c_str());
		RVdo.systemDataLoadUp(SystemDataFilename.c_str());
		//cout<<"Done loading up RV data obj"<<endl; //$$$$$$$$$$$$$$$$$$
	}
	// find the earliest epoch for the time of last pariapsis min/max settings
	double earliestEpoch = earliestEpochFinder(DIdo, RVdo);

	//Start seed for random number generator
	//create very high precision seed for generator (nanosecond precision)
	timespec time1;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
	int time_nsec=time1.tv_nsec;
	//cout<<"Starting time in nanoseconds for use as a random number generator seed = "<<time_nsec<<endl; //$$$$$$$$$$$$$$

	//choose generator
	CRandomSFMT1 RanGen(time_nsec);
	//this is the most advanced uniform random number generator, which combines SFMT and Mother-Of-All.

	// ******* STARTING SIMULATION!!!!  *************
	time_t startTime;
	startTime = time(NULL);

	string starterString;
	string numSamplesString =  numSamplesStringMaker(SSO.numSamples);
	ss<<"\nMCONLY: $$$$$$$$$$$$$$$$$$$  STARTING THE SIMULATOR  $$$$$$$$$$$$$$$$$" <<endl;
	ss<<"Number of sample orbits being created = " << numSamplesString <<endl;
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
    int numParams = 3;//there are 3 params that MUST always vary, so this is the minimum
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
    double one_over_nu;
    double chiSquareMin = SSO.chiSquaredMax;
    int bestOrbit = 0;

    bool vary_k = true;
    double K_proposed = 0;
    if (SSO.DIonly==true)
		vary_K = false;
	else if ((SSO.RVonly==true)&&((SSO.inclination_degMIN!=0)&&(SSO.inclination_degMAX!=0)))
		vary_K = false;
	else
		numParams+=1;

    if ((SSO.inclination_degMIN!=0)&&(SSO.inclination_degMAX!=0))
    	numParams+=1;
    if ((SSO.longAN_degMIN!=0)&&(SSO.longAN_degMAX!=0))
    	numParams+=1;
    else if ((SSO.periodMIN!=0)&&(SSO.periodMAX!=0))
    	numParams+=1;

    double TMIN;
    double TMAX;
    if (SSO.TimePeri_Min==0)
    	TMIN = earliestEpoch-SSO.periodMAX*365.0;//2455651.7;//
    else
    	TMIN = SSO.TimePeri_Min;
    if (SSO.TimePeri_Max==0)
    	TMAX = earliestEpoch;//2455653.3;//
    else
    	TMAX = SSO.TimePeri_Max;

    if (SSO.DIonly==false)
    {
    	for (int dataset=0; dataset<RVdo.epochs_RV.size();++dataset)
    		numParams+=1;

    }
    // ***** Start the samples loop *****
    for ( int sample=1; sample<(SSO.numSamples+1); sample++)
    {

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
			ss <<"\n"<< acceptedCounter<<"/"<<sample<<" Successful. "<<printsDone<<"/"<<SSO.numSamplePrints<<" completed at ";
			ss << asctime (timeinfo);
			ss << "latest reduced chiSquareds: DI = "<< DI_chiSquared_reduced<<", RV = "<<RV_chiSquared_reduced <<", Total = "<< TOTAL_chiSquared_reduced<<endl;
			ss << "lowest chiSquare so far = "<< chiSquareMin <<endl;
			printLine = ss.str();
			ss.clear();
			ss.str(std::string());
			cout<<printLine;
			SSlog<< printLine;
        }

        // Generate random numbers in the required ranges for the inputs to the orbCalc
        if ((SSO.inclination_degMIN!=0)&&(SSO.inclination_degMAX!=0))
        	DIt.inclination_deg = RanGen.UniformRandom(SSO.inclination_degMIN, SSO.inclination_degMAX);
        if ((SSO.longAN_degMIN!=0)&&(SSO.longAN_degMAX!=0))
        	DIt.longAN_deg = RanGen.UniformRandom(SSO.longAN_degMIN, SSO.longAN_degMAX);
        DIt.argPeri_deg = RanGen.UniformRandom(SSO.argPeri_degMIN, SSO.argPeri_degMAX);
        DIt.e = RanGen.UniformRandom(SSO.eMIN, SSO.eMAX);
        if ((SSO.simulate_StarPlanet==true)&&(SSO.fixed_planet_period==true))
        	DIt.period = RVdo.planet_P;
        else if ((SSO.periodMIN!=0)&&(SSO.periodMAX!=0))
        	DIt.period = RanGen.UniformRandom(SSO.periodMIN, SSO.periodMAX); //  [yrs]
        DIt.T = RanGen.UniformRandom(TMIN, TMAX); // thus between a full period before first observation and the time of first observation
        //reset RVoffsets_cur vector and get current param vals for K and offsets.
        vector<double> RVoffsets_cur;
        if (SSO.DIonly==false)
        {
        	if (vary_k)
        		K_proposed = RanGen.UniformRandom(SSO.K_MIN,SSO.K_MAX);

        	for (int dataset=0; dataset<RVoffsets_latest.size();++dataset)
        		RVoffsets_cur[dataset] = RanGen.UniformRandom(SSO.RVoffsetMINs[dataset],SSO.RVoffsetMAXs[dataset]);
		}
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

        double DI_chiSquared_original = 0;
        double RV_chiSquared_original = 0;
        double numDIepochsInternal = 0;
        double numRVepochsInternal = 0;

        multiEpochOrbCalcReturnType MEOCRT;
        if (SSO.RVonly==false)
        {
        	//cout<<"In DI block"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        	// #### DO STUFF FOR DI ORBIT CALCNS #####
        	if ( SSO.silent==false )
        		cout<<"Calculating DI orbit for this round "<<endl;

        	// Call the orbCalc to have it apply the model to the inputs and produce outputs
			MEOCRT = DIt.multiEpochOrbCalc();
			a_total_curr = MEOCRT.a_total;

			// Calculate the reduced chiSquared from the returned chiSquared
			DI_chiSquared_original = MEOCRT.chi_squared_total;
			numDIepochsInternal = DIdo.numEpochs_DI;
			one_over_nu = (1.0/((2.0*numDIepochsInternal)-5.0));
			DI_chiSquared_reduced = one_over_nu*DI_chiSquared_original;
			//SSO.silent=false;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			if ( SSO.silent==false )
			{
				cout<<"DI_chiSquared_original = "<<DI_chiSquared_original<<endl;
				cout<<"one_over_nu = "<<one_over_nu<<endl;
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
        	// ##### DO STUFF FOR RV CALCNS ######
        	vector<vector<double> > VRp_vector2;
        	vector<vector<double> > VRs_vector2;

        	// Load up ss or sp parts of RVdo with current trials
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
        		//RVdo.planet_MsinI  = DIt.Mass2 ;
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
        		//RVdo.star_Mass2  = DIt.Mass2 ;
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
        		VRCsp = VRcalcStarPlanetLoadUp(RVdo);
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
				if ((SSO.simulate_StarPlanet==true)&&(SSO.DIonly==false))
					a_total_curr = VRCsp.a_total;
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
        		if ((SSO.simulate_StarStar==true)&&(SSO.DIonly==false))
        			a_total_curr = VRCss.a_total;
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
//					RVoffsets_cur[0] = 10.486022;
//				if (dataset==1)
//					RVoffsets_cur[1] = 1.0;
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
					double  RV_chiSquared_cur = chiSquaredCalc((RVdo.RVs[dataset][epoch]-RVoffsets_cur[dataset]),updatedRV_inv_var,(planetVR+companionStarVR));
					RV_chiSquared_original = RV_chiSquared_original + RV_chiSquared_cur;
					if ( SSO.silent==false )
					{
						cout<<"\nWorking on epoch "<<epoch<<endl;
						//cout<<"\noffset = "<< RVoffsets_cur[dataset]<<endl;
						cout<<"RVdo.RVs[dataset][epoch] = "<<RVdo.RVs[dataset][epoch]<<", ("<<RVdo.RVs[dataset][epoch]<<" - "<<RVoffsets_cur[dataset]<<")="<<(RVdo.RVs[dataset][epoch]-RVoffsets_cur[dataset]) <<", planetVR= "<< planetVR<<", companionStarVR= "<< companionStarVR<<endl;
						cout<<"Difference = "<<RVdo.RVs[dataset][epoch]-RVoffsets_cur[dataset]-planetVR<<endl;
						cout<<"ChiSquared for this RV is = "<<RV_chiSquared_cur<<endl;
						cout<<"Total NON-reducedChiSquared so far is = "<<RV_chiSquared_original<<endl;
					}
				}
        	}
        	// calculate reduced version of ChiSquared for printing
        	numRVepochsInternal = RVdo.numEpochs_RV;
        	one_over_nu = (1.0/((1.0*numRVepochsInternal)-6.0));
        	RV_chiSquared_reduced = one_over_nu*RV_chiSquared_original;
        	if ( SSO.silent==false )
        	{
				cout<<"\nnumRVepochsInternal = "<< numRVepochsInternal <<endl;
				cout<<"one_over_nu = "<< one_over_nu <<endl;
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
        	RVoffsets_cur.push_back(0);
        	RV_chiSquared_reduced_lowest=0;
        }
        if (SSO.DIonly==false && SSO.RVonly==false)
        {
        	//cout<<"In both false block"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        	// THUS, do both DI and RV parts
        }

        // Do TOTAL chiSquared value calcs
        double chiSquared_TOTAL_original = DI_chiSquared_original+RV_chiSquared_original;
        one_over_nu = 1.0/(2.0*numDIepochsInternal+1.0*numRVepochsInternal-numParams);
        TOTAL_chiSquared_reduced = one_over_nu*chiSquared_TOTAL_original;

        //SSO.silent=false;//$$$$$$$$$$$$$$$$$$
		if ( SSO.silent==false )
		{
			cout<<"\nDI_chiSquared_original = "<<DI_chiSquared_original <<endl;
			cout<<"RV_chiSquared_original = "<< RV_chiSquared_original<<endl;
			cout<<"RV_chiSquared_original = "<< RV_chiSquared_original<<endl;
			cout<<"numDIepochsInternal = "<<numDIepochsInternal <<endl;
			cout<<"numRVepochsInternal = "<< numRVepochsInternal<<endl;
			cout<<"numParams = "<< numParams<<endl;
			cout<<"chiSquared_TOTAL_original = "<< chiSquared_TOTAL_original <<endl;
			cout<<"one_over_nu = "<< one_over_nu <<endl;
			cout<<"TOTAL_chiSquared_reduced = "<< TOTAL_chiSquared_reduced <<endl;
			cout<< "output reduced chi squared = "<< TOTAL_chiSquared_reduced <<endl;
		}
		//SSO.silent=true;//$$$$$$$$$$$$$$

		// Determine if the orbit should be accepted
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
			ODT.periods.push_back(DIt.period);
			ODT.inclination_degs.push_back(DIt.inclination_deg);
			ODT.argPeri_degs.push_back(DIt.argPeri_deg);
			// store outputs
			ODT.chiSquareds.push_back(TOTAL_chiSquared_reduced);
			ODT.a_totals.push_back(a_total_curr);
			ODT.Ks.push_back(K_proposed);
			ODT.RVoffsets.push_back(RVoffsets_cur);
			ODT.timesBeenHeres.push_back(1);

		}// Done storing accepted orbit parameters
    }//Done sample loops

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
    	ss<< "K = "<< ODT.Ks[SAOFO.bestOrbit]<<endl;
		for (int set=0;set<ODT.RVoffsets[bestOrbit].size();++set)
			if (ODT.RVoffsets[bestOrbit][set]!=0)
				ss<<"RVoffset for dataset "<<set<<", was = "<< ODT.RVoffsets[bestOrbit][set]<<endl;
    }

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
	logFileWriter(ODT.data_filename, LOGlines);

	// get output filename and write output data to it
    fileWriter(ODT);

	return EXIT_SUCCESS;
}
