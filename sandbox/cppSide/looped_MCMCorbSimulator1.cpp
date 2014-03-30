#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
//#include <string>
//#include <vector>
#include <math.h>
#include <time.h>
#include "Toolboxes/orbToolboxes.h" //Both the DI and RV tools headers are called in here, so it is an all in one toolbox header call
#include "simAnnealOrbSimFunc.h"
#include "MCMCorbSimFunc.h"
//#include <rnd/stocc.h> //a library with advanced non-uniform random number generators

using namespace std;

int main(int argc ,char *argv[])
{
	// print to indicate sim has started
	cout<<"\n*** C++ MCMC simulation has started ***\n";

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
		cout<<"MCMCorbSimulator is only set-up to handle settings";
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

	// ******* STARTING SIMULATION!!!!  *************
	time_t startTime;
	startTime = time(NULL);

	//###############################################################################
	//###############################################################################
	//Instantiate the MCMC object and load it up
	// ODT for the simAnneal in MCMC loop
	outputDataType ODT;
	ODT.data_filename = outputDataFilename;
	// ODT for the MCMC loop
	outputDataType ODT2;
	ODT2.data_filename = outputDataFilename;
	// ODT for 100 mod set outputs
	outputDataType ODTall;
	ODTall.data_filename = outputDataFilename;
	//cout<<"outputDataType instantiated"<<endl;//$$$$$$$$$$$$$$$$$$$$$

	// LOAD up 100 mod RV dataset
	RVdataObj RVdoMOD100;
	string LOGlinesALL;
	if ((SSO.RVonly==true) or (SSO.RVonly==false && SSO.DIonly==false))
	{
		//cout<<"starting to load up RV data obj"<<endl; //$$$$$$$$$$$$$$
		string RV100modDataFilename = SSO.settings_and_InputDataDir + "mod100Dataset.dat";
		cout<<"Using RV100modDataFilename:"<<RV100modDataFilename<<endl;
		RVdoMOD100.dataLoadUp(RV100modDataFilename.c_str());
		RVdoMOD100.systemDataLoadUp(SystemDataFilename.c_str());
		//cout<<"Done loading up RV data obj"<<endl; //$$$$$$$$$$$$$$$$$$
	}

	ss<<"**************************************************************************"<<endl;
	ss<<"****************** TOTAL 100 MOD DATASET SIMULATION STARTED **************"<<endl;
	ss<<"**************************************************************************"<<endl;

	for (int modDatasetNum = 0;modDatasetNum<100;++modDatasetNum)
	{
		time_t startTimeinternal;
		startTimeinternal = time(NULL);

		//*******************************************************************************
		//****************************** SimAnneal block start **************************
		//*******************************************************************************
		simAnealOrbFuncObj SAOFO;
		SAOFO.SSO = SSO;
		SAOFO.SSO.verbose = false;
		SAOFO.SSO.silent = true;
		SAOFO.SSO.numSamplePrints = 1;
		SAOFO.SSO.chiSquaredMax = 700;
		SAOFO.DIt = DIt;
		SAOFO.DIdo = DIdo;
		SAOFO.RVdo = RVdo;
		SAOFO.earliestEpoch = earliestEpoch;
		SAOFO.randSeed = time_nsec;
		SAOFO.ODT = ODT;
		SAOFO.numSamples_SA = 100000;
		SAOFO.tempStepSizePercent = 0.1;
		SAOFO.startTemp = 100;
		SAOFO.sigmaPercent = 2.45;

		SAOFO.RVdo.epochs_RV[0] = RVdoMOD100.epochs_RV[modDatasetNum] ;
		SAOFO.RVdo.RVs[0] = RVdoMOD100.RVs[modDatasetNum];
		SAOFO.RVdo.RV_inv_var[0] = RVdoMOD100.RV_inv_var[modDatasetNum];
		SAOFO.RVdo.numEpochs_RV = RVdoMOD100.epochs_RV[0].size();

		string starterString;
		string numSamplesString =  numSamplesStringMaker(SAOFO.numSamples_SA);
		ss<<"\nlooped_MCMC: $$$$$$$$$$$$$$$$$$$  Simulated Annealing Starting  $$$$$$$$$$$$$$$$$" <<endl;
		ss<<"Number of sample orbits being created = " << numSamplesString <<endl;
		starterString = ss.str();
		ss.clear();
		ss.str(std::string());
		cout<< starterString;
		SSlog<<starterString;

		//*** Run Simulated Annealing simulation first to get starting params
		if (SSO.silent==false)
				cout<<"Calling SAOFO simulator"<<endl;
		SAOFO.simulator();
		//SSO.silent=false;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		if (SSO.silent==false)
			cout<<"Returned from SAOFO simulator :-)"<<endl;
		//move output log string from simulation run to SSlog
		SSlog<<SAOFO.SSlogStr;

		int totalAccepted = SAOFO.ODT.es.size();
		SAOFO.ODT.numSamplesAccepted = totalAccepted;

		string printLine3="";
		// Get all best orbit values
		ss << "\nlooped_MCMC: $$$$$$$$$$$$$$$ Simulated Annealing COMPLETE $$$$$$$$$$$$$$$"<<endl;
		ss<< totalAccepted <<" orbits were accepted during simulation"<<endl;
		ss<< "\nBest orbit found:"<<endl;
		ss<< "chiSquaredMin = "<< SAOFO.chiSquareMin <<endl;
		ss<< "chiSquaredMin from vector = "<<SAOFO.ODT.chiSquareds[SAOFO.bestOrbit]<<endl;
		ss<< "mean Acceptance rate = "<<(double(totalAccepted)/double(SAOFO.numSamples_SA))<<endl;
		ss<< "LongAN = "<< SAOFO.ODT.longAN_degs[SAOFO.bestOrbit] <<endl;
		ss<< "e = "<< SAOFO.ODT.es[SAOFO.bestOrbit] <<endl;
		ss<< "To = "<< SAOFO.ODT.Ts[SAOFO.bestOrbit] <<endl;
		ss<< "period = "<< SAOFO.ODT.periods[SAOFO.bestOrbit] <<endl;
		ss<< "inclination = "<< SAOFO.ODT.inclination_degs[SAOFO.bestOrbit] <<endl;
		ss<< "argPeri = "<< SAOFO.ODT.argPeri_degs[SAOFO.bestOrbit] <<endl;
		for (int set=0;set<SAOFO.ODT.RVoffsets[SAOFO.bestOrbit].size();++set)
			if (SAOFO.ODT.RVoffsets[SAOFO.bestOrbit][set]!=0)
				ss<<"RVoffset for dataset "<<set<<", was = "<< SAOFO.ODT.RVoffsets[SAOFO.bestOrbit][set]<<endl;

		printLine3=ss.str();
		ss.clear();
		ss.str(std::string());
		cout<<printLine3;
		SSlog<< printLine3;

		if (SAOFO.chiSquareMin>SAOFO.SSO.chiSquaredMax)
			cout<<"looped_MCMC: ERROR: no orbit was found during Simulated Annealing below the chiSquaredMax in order to start a proper MCMC run!"<<endl;

		//*******************************************************************************
		//****************************** SimAnneal block end ****************************
		//*******************************************************************************
		else
		{
			//###############################################################################
			//######################### MCMC block start ####################################
			//###############################################################################
			string starterString2;
			string numSamplesString2 =  numSamplesStringMaker(SSO.numSamples);
			ss<<"\nlooped_MCMC: $$$$$$$$$$$$$$$$$$$  MCMC Simulation Starting  $$$$$$$$$$$$$$$$$" <<endl;
			ss<<"Working on MOD RV data set # "<<(modDatasetNum+1)<<endl;
			ss<<"Number of sample orbits being created = " << numSamplesString2 <<endl;
			starterString2 = ss.str();
			ss.clear();
			ss.str(std::string());
			cout<< starterString2;
			SSlog<<starterString2;

			MCMCorbFuncObj MCMCOFO;

	//		MCMCOFO.start_longAN = SAOFO.ODT.longAN_degs.back();
	//		MCMCOFO.start_e = SAOFO.ODT.es.back();
	//		MCMCOFO.start_T = SAOFO.ODT.Ts.back();
	//		MCMCOFO.start_period = SAOFO.ODT.periods.back();
	//		MCMCOFO.start_inc_deg = SAOFO.ODT.inclination_degs.back();
	//		MCMCOFO.start_argPeri = SAOFO.ODT.argPeri_degs.back();
	//		MCMCOFO.start_offsets = SAOFO.ODT.RVoffsets.back();

			MCMCOFO.start_longAN = SAOFO.ODT.longAN_degs[SAOFO.bestOrbit];
			MCMCOFO.start_e = SAOFO.ODT.es[SAOFO.bestOrbit];
			MCMCOFO.start_T = SAOFO.ODT.Ts[SAOFO.bestOrbit];
			MCMCOFO.start_period = SAOFO.ODT.periods[SAOFO.bestOrbit];
			MCMCOFO.start_inc_deg = SAOFO.ODT.inclination_degs[SAOFO.bestOrbit];
			MCMCOFO.start_argPeri = SAOFO.ODT.argPeri_degs[SAOFO.bestOrbit];
			MCMCOFO.start_offsets = SAOFO.ODT.RVoffsets[SAOFO.bestOrbit];

			MCMCOFO.sigmaPercent = SAOFO.sigmaPercent_latest;
			MCMCOFO.inclination_deg_sigma = SAOFO.inclination_deg_sigma;
			MCMCOFO.e_sigma = SAOFO.e_sigma;
			MCMCOFO.longAN_deg_sigma = SAOFO.longAN_deg_sigma;
			MCMCOFO.period_sigma = SAOFO.period_sigma;
			MCMCOFO.argPeri_deg_sigma = SAOFO.argPeri_deg_sigma;
			MCMCOFO.T_sigma = SAOFO.T_sigma;

			MCMCOFO.SSO = SSO;
			MCMCOFO.DIt = DIt;
			MCMCOFO.DIdo = DIdo;
			MCMCOFO.RVdo = RVdo;
			MCMCOFO.earliestEpoch = earliestEpoch;
			MCMCOFO.randSeed = time_nsec;
			MCMCOFO.ODT = ODT2;

			MCMCOFO.RVdo.epochs_RV[0] = RVdoMOD100.epochs_RV[modDatasetNum] ;
			MCMCOFO.RVdo.RVs[0] = RVdoMOD100.RVs[modDatasetNum];
			MCMCOFO.RVdo.RV_inv_var[0] = RVdoMOD100.RV_inv_var[modDatasetNum];
			MCMCOFO.RVdo.numEpochs_RV = RVdoMOD100.epochs_RV[0].size();

			MCMCOFO.paramChangeStepSizePercent = 0.1;
			MCMCOFO.sigmaPercent = 2.45;

			//###############################################################################
			//######################### MCMC block end ######################################
			//###############################################################################
			if (MCMCOFO.chiSquareMin>MCMCOFO.SSO.chiSquaredMax)
				cout<<"looped_MCMC: ERROR: no orbit was found during Simulated Annealing below the chiSquaredMax in order to start a proper MCMC run!"<<endl;
			else
			{
				//*** Run Simulated Annealing simulation first to get starting params
				if (SSO.silent==false)
						cout<<"Calling MCMCOFO simulator"<<endl;
				MCMCOFO.simulator();
				//SSO.silent=false;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
				if (SSO.silent==false)
					cout<<"Returned from MCMCOFO simulator :-)"<<endl;
				//move output log string from simulation run to SSlog
				SSlog<<MCMCOFO.SSlogStr;
				//###############################################################################
				//###############################################################################

				int totalAccepted2 = MCMCOFO.ODT.es.size();
				MCMCOFO.ODT.numSamplesAccepted = totalAccepted;

				string printLine4="";
				// Get all best orbit values
				ss << "\nlooped_MCMC: $$$$$$$$$$$$$$$ MCMC SIMULATOR COMPLETE $$$$$$$$$$$$$$$"<<endl;
				ss<< totalAccepted2 <<" orbits were accepted during simulation"<<endl;
				ss<< "\nBest orbit found:"<<endl;
				ss<< "chiSquaredMin = "<< MCMCOFO.chiSquareMin <<endl;
				ss<< "chiSquaredMin from vector = "<<MCMCOFO.ODT.chiSquareds[MCMCOFO.bestOrbit]<<endl;
				ss<< "mean Acceptance rate = "<<(double(totalAccepted2)/double(SSO.numSamples))<<endl;
				ss<< "LongAN = "<< MCMCOFO.ODT.longAN_degs[MCMCOFO.bestOrbit] <<endl;
				ss<< "e = "<< MCMCOFO.ODT.es[MCMCOFO.bestOrbit] <<endl;
				ss<< "To = "<< MCMCOFO.ODT.Ts[MCMCOFO.bestOrbit] <<endl;
				ss<< "period = "<< MCMCOFO.ODT.periods[MCMCOFO.bestOrbit] <<endl;
				ss<< "inclination = "<< MCMCOFO.ODT.inclination_degs[MCMCOFO.bestOrbit] <<endl;
				ss<< "argPeri = "<< MCMCOFO.ODT.argPeri_degs[MCMCOFO.bestOrbit] <<endl;
				for (int set=0;set<MCMCOFO.ODT.RVoffsets[MCMCOFO.bestOrbit].size();++set)
					if (MCMCOFO.ODT.RVoffsets[MCMCOFO.bestOrbit][set]!=0)
						ss<<"RVoffset for dataset "<<set<<", was = "<< MCMCOFO.ODT.RVoffsets[MCMCOFO.bestOrbit][set]<<endl;

				//Load up ALL 100 output data
				ODTall.numSamplesAccepted = modDatasetNum+1;
				//inputs
				ODTall.longAN_degs.push_back(MCMCOFO.ODT.longAN_degs[MCMCOFO.bestOrbit]);
				ODTall.es.push_back(MCMCOFO.ODT.es[MCMCOFO.bestOrbit]);
				ODTall.Ts.push_back(MCMCOFO.ODT.Ts[MCMCOFO.bestOrbit]);
				ODTall.periods.push_back(MCMCOFO.ODT.periods[MCMCOFO.bestOrbit]);
				ODTall.inclination_degs.push_back(MCMCOFO.ODT.inclination_degs[MCMCOFO.bestOrbit]);
				ODTall.argPeri_degs.push_back(MCMCOFO.ODT.argPeri_degs[MCMCOFO.bestOrbit]);
				ODTall.a_totals.push_back(MCMCOFO.ODT.a_totals[MCMCOFO.bestOrbit]);
				//outputs
				ODTall.chiSquareds.push_back(MCMCOFO.ODT.chiSquareds[MCMCOFO.bestOrbit]) ;
				ODTall.RVoffsets.push_back(MCMCOFO.ODT.RVoffsets[MCMCOFO.bestOrbit]);

				// calc how long sim took to run and print appropriate duration string
				time_t endTime;
				endTime = time(NULL);
				int timeElapsed = endTime-startTimeinternal;

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
				printLine4=ss.str();
				ss.clear();
				ss.str(std::string());
				cout<<printLine4;
				// load up log with all prints in SSlog stringstream
				// then write it to file.
				SSlog<< printLine4;
				string LOGlines;
				LOGlinesALL = SSlog.str()+"\n\n";
			}
		}
	}
	// calc how long sim took to run and print appropriate duration string
	time_t endTime;
	endTime = time(NULL);
	int timeElapsed = endTime-startTime;

	ss<< "\n\n"<<endl;
	ss<<"**************************************************************************"<<endl;
	ss<<"****************** TOTAL 100 MOD DATASET SIMULATION COMPLETE *************"<<endl;
	ss<<"**************************************************************************"<<endl;
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

	// Write output data of Sim Anneal to file
	string ODTall_data_filename = ODTall.data_filename;
	//ODTall_data_filename =  filenamePrepend(ODTall_data_filename, "MCMC_100ModDatasetBestFits_");
	ODTall.data_filename = ODTall_data_filename;
	fileWriter(ODTall);

	logFileWriter(ODTall_data_filename, LOGlinesALL);


	return EXIT_SUCCESS;
}
