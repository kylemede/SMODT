#include <iostream>
#include <iomanip>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <math.h>
#include <time.h>
#include "toolboxes/orbToolboxes.h" //Both the DI and RV tools headers are called in here, so it is an all in one toolbox header call
#include "simAnnealOrbSimFunc.h"
#include "MCMCorbSimFunc.h"
using namespace std;
/**
 * This is a bunch of test documentation.
 * @author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
 * /*! \file */

int main(int argc ,char *argv[])
{
	/**
	This is the "main" that performs the MCMC or Simulated Annealing process for a single chain
	of a multi or single chain/thread/process simulation started by the
	Python file SMODTstarter.py and controlled with the
	SimulationManager.  It will perform the set up tasks, call the simAnnealOrbSimFunc and
	MCMCorbSimFunc to perform the Simulated Annealing and MCMC stages for the
	chain based on the input settings.
	In the case of a request to perform full MCMC, after the MCMC process is complete
	it will perform the requested wrap up and statistical calculations
	requested.  Finally, it will write all the output files to disk for the wrap up tasks and
	plotting done by Python after all the chains are complete.
	*/
	std::stringstream ss;
	std::stringstream SSlog;
	generalTools GT;
	// print to indicate sim has started
	//ss<<"\n*** C++ MCMC simulation has started ***\n";
	ss<<"\n\n $$$$$ Entered MCMC simulator in C++ $$$\n\n"<<endl;
	// pull in simulator settings filename and output data filename
	// from command line args
	string settingsFilename;
	string outputDataFilename;
	ss<<"There were, "<<argc<<", arguments provided.  The were:"<<endl;
	if (argc>=2)
	{
		settingsFilename = argv[1];//including directory!
		ss<<"Using settings file:\n"<<settingsFilename<<endl;
	}
	if (argc==3)
	{
		outputDataFilename = argv[2];//including directory!
		ss<<"Output data will be written to:\n"<<outputDataFilename<<endl;
	}
	if (argc>3)
	{
		ss<<"MCMCorbSimulator is only set-up to handle settings";
		ss<<"filename and output data filename command line ";
		ss<<"arguments so far."<<endl;
	}
	else if (argc==1)
	{
		settingsFilename = "/mnt/Data1/Todai_Work/EclipseWorkspace/SMODT/settings_and_InputData/SimSettings.txt";
		ss<<"No settings filename provided, so using default: ";
		ss<<settingsFilename<<endl;
		ss<<"No output data filename provided through command line";
		ss<<" so, value within settings file or default being used."<<endl;
	}

	//cout<<"Command line arguments all loaded"<<endl; //$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$

	// instantiate the settings object and load it up
	SimSettingsObj SSO;
	SSO.settingsLoadUp(settingsFilename.c_str());

	//cout<< "Settings obj all loaded up"<<endl; //$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$

	// instantiate and load up DI and RV data objects as needed
	string outSettingsDir = GT.findDirectory(settingsFilename.c_str())+"/";
	string SystemDataFilename = outSettingsDir + SSO.SystemDataFilename;
	ss<<"SystemDataFilename = "<< SystemDataFilename <<endl;
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
	if ((SSO.verbose==true)&&(SSO.SILENT==false))
		cout<< inputsStr;

	// Find out what 'mode' simulator will run in then log and print it to screen
	string simModeStr;
	if (SSO.RVonly==false && SSO.DIonly==true)
		ss<<"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nSimulation running in DIonly mode\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	else if (SSO.RVonly==false && SSO.DIonly==false)
		ss<<"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nSimulation running in 3D mode\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	else
		ss<<"\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nSimulation running in RVonly mode\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
	simModeStr = ss.str();
	ss.clear();
	ss.str(std::string());
	SSlog<<simModeStr;
	if (SSO.SILENT==false)
		cout<<simModeStr<<endl;

	//cout<<"\n$$$$ about to try and load up sys data"<<endl;
	SYSdo.systemDataLoadUp(SystemDataFilename.c_str());
	//cout<<"DI and RV data objects instantiated"<<endl; //$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$
	if ((SSO.DIonly==true) or (SSO.RVonly==false && SSO.DIonly==false))
	{
		//cout<<"Starting to load up DI data obj"<<endl; //$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$
		string DIdataFilename = outSettingsDir + SSO.DIdataFilename;
		DIdo.dataLoadUp(DIdataFilename.c_str());
		DIdo.systemDataLoadUp(SystemDataFilename.c_str());
		// instantiate DI tools obj and load it up with same values
		DIt = GT.DItoolsParamLoadUp(DIdo);
		DIt.verbose=SSO.verbose;
		//cout<<"DI data object loaded up"<<endl; //$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$
	}
	if ((SSO.RVonly==true) or (SSO.RVonly==false && SSO.DIonly==false))
	{
		//cout<<"starting to load up RV data obj"<<endl; //$$$$$$$$$$$$$$
		string RVdataFilename = outSettingsDir + SSO.RVdataFilename;
		RVdo.dataLoadUp(RVdataFilename.c_str());
		RVdo.systemDataLoadUp(SystemDataFilename.c_str());
		//cout<<"RVdo.planet_argPeri = "<<RVdo.planet_argPeri<<endl;//$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$
		//cout<<"Done loading up RV data obj"<<endl; //$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$
	}
	// find the earliest epoch for the time of last pariapsis min/max settings
	double earliestEpoch = GT.earliestEpochFinder(DIdo, RVdo);
	//cout<<"\n\nearliestEpoch = "<< earliestEpoch<<"\n\n"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	//Start seed for random number generator
	//create very high precision seed for generator (nanosecond precision)
	timespec time1;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
	int time_nsec=time1.tv_nsec;
	//cout<<"Starting time in nanoseconds for use as a random number generator seed = "<<time_nsec<<endl; //$$$$$$$ DEBUGGING $$$$$$$$$$$$$$


//	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//	//$$$$$$$$$$$$$$$$ SOME TESTER CODE FOR C++ GAUSSIAN RND NUMBER GENERATION $$$$$$
//	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//	if (false)//$$$$$$$$$$$$ DEBUGGING  $$$$$$
//	{
//		StochasticLib1 RanGen2(time_nsec);
//
//		vector<double> randNumbersAry;
//
//		for(int i=0;i<100000;++i)
//		{
//			//double randNumber = RanGen.UniformRandom(1.0, 100.0);
//			double randNumber = RanGen2.NormalTrunc(0.0,10.0,50.0);
//			if (false)//$$$$$$$$$$$$ DEBUGGING  $$$$$$
//				cout<<"Random number is: "<<randNumber<<endl;
//			randNumbersAry.push_back(randNumber);
//		}
//
//		ofstream file;
//		string data_filename = "/mnt/Data1/Todai_Work/Data/data_SMODT/TestCplusplusTrucNormal.dat";
//		//cout<<"\n***************************************************************"<<endl;
//		//cout<<"Writing file to: "<<data_filename<<endl;
//		//cout<<"***************************************************************\n"<<endl;
//		file.open(data_filename.c_str()) ;     //open the output file
//
//		for (int r=0; r<randNumbersAry.size();r++)
//			file<<randNumbersAry[r]<<endl;
//		file.close();
//	}
//	//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



	// ******* STARTING SIMULATION!!!!  *************
	time_t startTime = time(NULL);
	//###############################################################################
	//###############################################################################
	//Instantiate the simAnneal object and load it up
	//
	outputDataType ODT;
	ODT.data_filename = outputDataFilename;
	//cout<<"outputDataType instantiated"<<endl;//$$$$$$$$$$$$ DEBUGGING $$$$$$$$$$$$$$$$$
	//### Find out how many samples to draw and how often to save accepted steps to output data array ###
	int numSamplesUse = SSO.numSamples_SimAnneal;
	if (SSO.simAnneal)
		numSamplesUse = SSO.numSamples;
	int saveEachInt;
	if (numSamplesUse<300001)
		saveEachInt = 1;
	else if ((numSamplesUse>300001)&&(numSamplesUse<30000001))
		saveEachInt=10;
	else if ((numSamplesUse>30000000)&&(numSamplesUse<300000000))
		saveEachInt=100;
	else
		saveEachInt=1000;

	simAnealOrbFuncObj SAOFO;
	SAOFO.SSO = SSO;
	SAOFO.SSO.verbose = SSO.verbose;
	SAOFO.SSO.quiet = SSO.quiet;
	SAOFO.SSO.numSamplePrints = SSO.numSamplePrints;
	SAOFO.SSO.chiSquaredMax = SSO.chiSquaredMax;
	SAOFO.DIt = DIt;
	SAOFO.DIdo = DIdo;
	SAOFO.RVdo = RVdo;
	SAOFO.SYSdo = SYSdo;
	SAOFO.earliestEpoch = earliestEpoch;
	SAOFO.randSeed = time_nsec;
	SAOFO.ODT = ODT;
	SAOFO.numSamples_SA = numSamplesUse;
	SAOFO.tempStepSizePercent = 0.025;//Fixed value found by suitable through testing.
	SAOFO.startTemp = SSO.startTemp;
	SAOFO.sigmaPercent = 1.0;//Just a starting value that is currently overridden inside the simAnneal func anyway.
	SAOFO.saveEachInt = saveEachInt;

	string starterString;
	string numSamplesString =  GT.numSamplesStringMaker(SAOFO.numSamples_SA);
	ss<<"\nMCMC: $$$$$$$$$$$$$$$$$$  Simulated Annealing Starting  $$$$$$$$$$$$$$$" <<endl;
	ss<<"Number of sample orbits being created = " << numSamplesString <<endl;
	starterString = ss.str();
	ss.clear();
	ss.str(std::string());
	if (SSO.SILENT==false)
		cout<< starterString;
	else
		cout<<"MCMC: $$$$$$$$$$$$$$$$$$  Simulated Annealing Starting  $$$$$$$$$$$$$$$" <<endl;
	SSlog<<starterString;

	//Starting log message that ensures vital inputs are logged.
	if (true)
	{
		string startParmsString;
		ss<<"Number of samples = "<<SSO.numSamples<<endl;
		ss<<"randSeed = "<< time_nsec <<endl;
		ss<<"tempStepSizePercent = "<<SAOFO.tempStepSizePercent <<endl;
		ss<<"startTemp = "<< SAOFO.startTemp<< endl;
		ss<<"sigmaPercent = "<< SAOFO.sigmaPercent<< endl;
		ss<<"saveEachInt = "<<SAOFO.saveEachInt<<endl;
		startParmsString = ss.str();
		ss.clear();
		ss.str(std::string());
		if (SSO.SILENT==false)
			cout<< startParmsString;
		SSlog<<startParmsString;
	}
	//*** Run Simulated Annealing simulation first to get starting params
	if ((SSO.quiet==false)&&(SSO.SILENT==false))
			cout<<"Calling SAOFO function"<<endl;
	SAOFO.simulator();
	//SSO.quiet=false;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	if ((SSO.quiet==false)&&(SSO.SILENT==false))
		cout<<"Returned from SAOFO function :-)"<<endl;
	//move output log string from simulation run to SSlog
	SSlog<<SAOFO.SSlogStr;
	//###############################################################################
	//###############################################################################

    int totalAccepted = SAOFO.ODT.es.size();
    SAOFO.ODT.numSamplesAccepted = totalAccepted;

    string printLine3="";
	// Get all best orbit values
	ss << "\nMCMC: $$$$$$$$$$$$$$$ Simulated Annealing COMPLETE $$$$$$$$$$$$$$$"<<endl;
	ss<< totalAccepted <<" orbits were accepted during simulation"<<endl;
	ss<< "\nBest orbit found at step "<<SAOFO.bestOrbit<<" :"<<endl;
	ss<<"chiSquareMin_reduced = "<<SAOFO.chiSquaredMin*SAOFO.one_over_nu_TOTAL<<endl;
	ss<< "chiSquaredMin = "<< SAOFO.chiSquaredMin <<endl;
	ss<< "chiSquaredMin from vector = "<<SAOFO.ODT.chiSquareds[SAOFO.bestOrbit]<<endl;
	ss<< "One before chiSquaredMin from vector = "<<SAOFO.ODT.chiSquareds[SAOFO.bestOrbit-1]<<endl;
	ss<< "One after chiSquaredMin from vector = "<<SAOFO.ODT.chiSquareds[SAOFO.bestOrbit+1]<<endl;
	ss<< "mean Acceptance rate = "<<(double(totalAccepted)/double(SAOFO.numSamples_SA))<<endl;
	ss<< "LongAN = "<< SAOFO.ODT.longAN_degs[SAOFO.bestOrbit] <<endl;
	ss<< "e = "<< SAOFO.ODT.es[SAOFO.bestOrbit] <<endl;
	ss<< "To = "<< fixed <<SAOFO.ODT.Ts[SAOFO.bestOrbit] <<endl;
	ss<< "Tc = "<< SAOFO.ODT.Tcs[SAOFO.bestOrbit]<<endl;
	ss<< "period = "<< SAOFO.ODT.periods[SAOFO.bestOrbit] <<endl;
	ss<< "inclination = "<< SAOFO.ODT.inclination_degs[SAOFO.bestOrbit] <<endl;
	ss<< "argPeri = "<< SAOFO.ODT.argPeri_degs[SAOFO.bestOrbit] <<endl;
	ss<< "a_total = "<< SAOFO.ODT.a_totals[SAOFO.bestOrbit] <<endl;
	if (SAOFO.SSO.DIonly==false)
	{
		ss<< "K = "<< SAOFO.ODT.Ks[SAOFO.bestOrbit]<<endl;
		for (int set=0;set<SAOFO.ODT.RVoffsets[SAOFO.bestOrbit].size();++set)
			ss<<"RVoffset for dataset "<<set<<", was = "<< SAOFO.ODT.RVoffsets[SAOFO.bestOrbit][set]<<endl;
	}

	if (true)
	{
		ss<< "\nLast orbit found:"<<endl;
		ss<< "chiSquaredMin from vector = "<<SAOFO.ODT.chiSquareds.back() <<endl;
		ss<<"chiSquareMin_reduced = "<<SAOFO.ODT.chiSquareds.back()*SAOFO.one_over_nu_TOTAL<<endl;
		ss<< "LongAN = "<< SAOFO.ODT.longAN_degs.back()<<endl;
		ss<< "e = "<< SAOFO.ODT.es.back()<<endl;
		ss<< "To = "<< SAOFO.ODT.Ts.back()<<endl;
		ss<< "Tc = "<< SAOFO.ODT.Tcs.back() <<endl;
		ss<< "period = "<< SAOFO.ODT.periods.back()<<endl;
		ss<< "inclination = "<<SAOFO.ODT.inclination_degs.back() <<endl;
		ss<< "argPeri = "<<SAOFO.ODT.argPeri_degs.back() <<endl;
		ss<< "a_total = "<<SAOFO.ODT.a_totals.back() <<endl;
		if (SAOFO.SSO.DIonly==false)
		{
			ss<< "K = "<< SAOFO.ODT.Ks.back()<<endl;
			for (int set=0;set<SAOFO.ODT.RVoffsets.back().size();++set)
				ss<<"RVoffset for dataset "<<set<<", was = "<< SAOFO.ODT.RVoffsets.back()[set]<<endl;
		}
	}
	if (true)
	{
		ss<< "\nSecond to last orbit found:"<<endl;
		ss<< "chiSquaredMin from vector = "<<SAOFO.ODT.chiSquareds[totalAccepted-2] <<endl;
		ss<<"chiSquareMin_reduced = "<<SAOFO.ODT.chiSquareds[totalAccepted-2]*SAOFO.one_over_nu_TOTAL<<endl;
		ss<< "LongAN = "<< SAOFO.ODT.longAN_degs[totalAccepted-2]<<endl;
		ss<< "e = "<< SAOFO.ODT.es[totalAccepted-2]<<endl;
		ss<< "To = "<< SAOFO.ODT.Ts[totalAccepted-2]<<endl;
		ss<< "Tc = "<< SAOFO.ODT.Tcs[totalAccepted-2] <<endl;
		ss<< "period = "<< SAOFO.ODT.periods[totalAccepted-2]<<endl;
		ss<< "inclination = "<<SAOFO.ODT.inclination_degs[totalAccepted-2] <<endl;
		ss<< "argPeri = "<<SAOFO.ODT.argPeri_degs[totalAccepted-2] <<endl;
		ss<< "a_total = "<<SAOFO.ODT.a_totals[totalAccepted-2] <<endl;
		if (SAOFO.SSO.DIonly==false)
		{
			ss<< "K = "<< SAOFO.ODT.Ks[totalAccepted-2]<<endl;
			for (int set=0;set<SAOFO.ODT.RVoffsets[totalAccepted-2].size();++set)
				ss<<"RVoffset for dataset "<<set<<", was = "<< SAOFO.ODT.RVoffsets[totalAccepted-2][set]<<endl;
		}
	}
	if (true)
	{
		ss<< "\nFinal sigma values for each parameter:"<<endl;
		double sig = SAOFO.inclination_deg_sigma*(100.0/(SAOFO.SSO.inclination_degMAX-SAOFO.SSO.inclination_degMIN));
	    ss<< "inclination_deg_sigma = "<<SAOFO.inclination_deg_sigma<<" = sigPercent = "<<sig<<endl;
		sig = SAOFO.e_sigma*(100.0/(SAOFO.SSO.eMAX-SAOFO.SSO.eMIN));
	    ss<< "e_sigma = "<<SAOFO.e_sigma<<" = sigPercent = "<<sig<<endl;
	    sig = SAOFO.longAN_deg_sigma*(100.0/(SAOFO.SSO.longAN_degMAX-SAOFO.SSO.longAN_degMIN));
	    ss<< "longAN_deg_sigma = "<<SAOFO.longAN_deg_sigma<<" = sigPercent = "<<sig<<endl;
	    sig = SAOFO.period_sigma*(100.0/(SAOFO.SSO.periodMAX-SAOFO.SSO.periodMIN));
	    ss<< "period_sigma = "<<SAOFO.period_sigma<<" = sigPercent = "<<sig<<endl;
	    sig = SAOFO.argPeri_deg_sigma*(100.0/(SAOFO.SSO.argPeri_degMAX-SAOFO.SSO.argPeri_degMIN));
	    ss<< "argPeri_deg_sigma = "<<SAOFO.argPeri_deg_sigma<<" = sigPercent = "<<sig<<endl;
	    sig = SAOFO.T_sigma*(100.0/(SAOFO.SSO.T_Max-SAOFO.SSO.T_Min));
	    ss<< "T_sigma = "<<fixed<<SAOFO.T_sigma<<" = sigPercent = "<<sig<<endl;
	    sig = SAOFO.a_total_sigma*(100.0/(SAOFO.SSO.a_totalMAX-SAOFO.SSO.a_totalMIN));
	    ss<< "a_total_sigma = "<<fixed<<SAOFO.a_total_sigma<<" = sigPercent = "<<sig<<endl;
	    ss<< "sqrtESinomega_sigma = "<<fixed<<SAOFO.sqrtESinomega_sigma<<" = sigPercent = "<<sig<<endl;
	    ss<< "sqrtECosomega_sigma = "<<fixed<<SAOFO.sqrtECosomega_sigma<<" = sigPercent = "<<sig<<endl;
	    if (SAOFO.SSO.DIonly==false)
		{
			sig = SAOFO.K_sigma*(100.0/(SAOFO.SSO.K_MAX-SAOFO.SSO.K_MIN));
			ss<< "K_sigma = "<<SAOFO.K_sigma<<" = sigPercent = "<<sig<<endl;
			for (int set=0;set<SAOFO.ODT.RVoffsets[0].size();++set)
			{
				if (SAOFO.SSO.RVoffsetMAXs[set]!=0)
				{
					sig = SAOFO.offset_sigmas[set]*(100.0/(SAOFO.SSO.RVoffsetMAXs[set]-SAOFO.SSO.RVoffsetMINs[set]));
					ss<<"RVoffset sigma for dataset "<<set<<", was = "<< SAOFO.offset_sigmas[set]<<" = sigPercent = "<<sig<<endl;
				}
				else
					ss<<"No RVoffset sigma for dataset "<<set<<" as RVoffsetMAXs[set]==0"<<endl;
			}
		}
	}

	ss<<"\n one_over_nu values for all three chiSquare calcs:"<<endl;
	ss<<"one_over_nu_DI = "<< SAOFO.one_over_nu_DI <<endl;
	ss<<"one_over_nu_RV = "<< SAOFO.one_over_nu_RV <<endl;
	ss<<"one_over_nu_TOTAL = "<< SAOFO.one_over_nu_TOTAL <<endl;

	printLine3=ss.str();
	ss.clear();
	ss.str(std::string());
	if (SSO.SILENT==false)
		cout<<printLine3;
	SSlog<< printLine3;

	if (SSO.simAnneal==false)
	{
		// Write output data of Sim Anneal to file
		string SAFO_data_filename = SAOFO.ODT.data_filename;
		SAFO_data_filename =  GT.filenamePrepend(SAFO_data_filename, "SimAnneal_");
		SAOFO.ODT.data_filename = SAFO_data_filename;
		GT.fileWriter(SAOFO.ODT);

		//###############################################################################
		//###############################################################################
		//Instantiate the MCMC object and load it up
		//
		outputDataType ODT2;
		ODT2.data_filename = outputDataFilename;
		//cout<<"outputDataType instantiated"<<endl;//$$$$$$$$$$$$$$$$$$$$$

		string starterString2;
		string numSamplesString2 =  GT.numSamplesStringMaker(SSO.numSamples);
		ss<<"\nMCMC: $$$$$$$$$$$$$$$$$$$$$$$$$  MCMC Starting  $$$$$$$$$$$$$$$$$$$$$$$" <<endl;
		ss<<"Number of sample orbits being created = " << numSamplesString2 <<endl;
		starterString2 = ss.str();
		ss.clear();
		ss.str(std::string());
		if (SSO.SILENT==false)
			cout<< starterString2;
		else
			cout<<"MCMC: $$$$$$$$$$$$$$$$$$$$$$$$$  MCMC Starting  $$$$$$$$$$$$$$$$$$$$$$$" <<endl;
		SSlog<<starterString2;

		MCMCorbFuncObj MCMCOFO;
		//### Find out how often to save accepted steps to output data array ###
		if (SSO.numSamples<300001)
			saveEachInt = 1;
		else if ((SSO.numSamples>300001)&&(SSO.numSamples<30000001))
			saveEachInt=10;
		else if ((SSO.numSamples>30000000)&&(SSO.numSamples<30000000))
			saveEachInt=100;
		else
			saveEachInt=1000;

		// Best fit values
	    MCMCOFO.start_longAN = SAOFO.ODT.longAN_degs[SAOFO.bestOrbit];
	    MCMCOFO.start_e = SAOFO.ODT.es[SAOFO.bestOrbit];
	    MCMCOFO.start_T = SAOFO.ODT.Ts[SAOFO.bestOrbit];
	    MCMCOFO.start_Tc = SAOFO.ODT.Tcs[SAOFO.bestOrbit];
	    MCMCOFO.start_period = SAOFO.ODT.periods[SAOFO.bestOrbit];
	    MCMCOFO.start_inc_deg = SAOFO.ODT.inclination_degs[SAOFO.bestOrbit];
	    MCMCOFO.start_argPeri = SAOFO.ODT.argPeri_degs[SAOFO.bestOrbit];
	    MCMCOFO.start_offsets = SAOFO.ODT.RVoffsets[SAOFO.bestOrbit];
	    MCMCOFO.start_K = SAOFO.ODT.Ks[SAOFO.bestOrbit];
	    MCMCOFO.start_a_total = SAOFO.ODT.a_totals[SAOFO.bestOrbit];

		  // Back values
	//	MCMCOFO.start_longAN = SAOFO.ODT.longAN_degs.back();
	//	MCMCOFO.start_e = SAOFO.ODT.es.back();
	//	MCMCOFO.start_T = SAOFO.ODT.Ts.back();
	//	MCMCOFO.start_Tc = SAOFO.ODT.Tcs.back();
	//	MCMCOFO.start_period = SAOFO.ODT.periods.back();
	//	MCMCOFO.start_inc_deg = SAOFO.ODT.inclination_degs.back();
	//	MCMCOFO.start_argPeri = SAOFO.ODT.argPeri_degs.back();
	//	MCMCOFO.start_offsets = SAOFO.ODT.RVoffsets.back();
	//	MCMCOFO.start_K = SAOFO.ODT.Ks.back();
	//	MCMCOFO.start_a_total = SAOFO.ODT.a_totals.back();

		// Second from back params
	//    MCMCOFO.start_longAN = SAOFO.ODT.longAN_degs[totalAccepted-2];
	//   MCMCOFO.start_e = SAOFO.ODT.es[totalAccepted-2];
	//    MCMCOFO.start_T = SAOFO.ODT.Ts[totalAccepted-2];
	//    MCMCOFO.start_Tc = SAOFO.ODT.Tcs[totalAccepted-2];
	//    MCMCOFO.start_period = SAOFO.ODT.periods[totalAccepted-2];
	//    MCMCOFO.start_inc_deg = SAOFO.ODT.inclination_degs[totalAccepted-2];
	//    MCMCOFO.start_argPeri = SAOFO.ODT.argPeri_degs[totalAccepted-2];
	//    MCMCOFO.start_offsets = SAOFO.ODT.RVoffsets[totalAccepted-2];
	//    MCMCOFO.start_K = SAOFO.ODT.Ks[totalAccepted-2];
	//    MCMCOFO.start_a_total = SAOFO.ODT.a_totals[totalAccepted-2];

		MCMCOFO.sigmaPercent = SAOFO.sigmaPercent_latest;
		MCMCOFO.inclination_deg_sigma = SAOFO.inclination_deg_sigma;
		MCMCOFO.e_sigma = SAOFO.e_sigma;
		MCMCOFO.longAN_deg_sigma = SAOFO.longAN_deg_sigma;
		MCMCOFO.period_sigma = SAOFO.period_sigma;
		MCMCOFO.argPeri_deg_sigma = SAOFO.argPeri_deg_sigma;
		MCMCOFO.T_sigma = SAOFO.T_sigma;
		MCMCOFO.TMIN = SAOFO.TMIN;
		MCMCOFO.TMAX = SAOFO.TMAX;
		MCMCOFO.a_total_sigma = SAOFO.a_total_sigma;
		MCMCOFO.K_sigma = SAOFO.K_sigma;
		MCMCOFO.sqrtESinomega_sigma = SAOFO.sqrtESinomega_sigma;
		MCMCOFO.sqrtECosomega_sigma = SAOFO.sqrtECosomega_sigma;

		MCMCOFO.offset_sigmas = SAOFO.offset_sigmas;
		MCMCOFO.vary_K = SAOFO.vary_K;
		MCMCOFO.numParams = SAOFO.numParams;
		MCMCOFO.numDIparams = SAOFO.numDIparams ;
		MCMCOFO.numRVparams = SAOFO.numRVparams;
		MCMCOFO.paramsToVaryIntsAry = SAOFO.paramsToVaryIntsAry;
		MCMCOFO.one_over_nu_RV = SAOFO.one_over_nu_RV;
		MCMCOFO.one_over_nu_DI = SAOFO.one_over_nu_DI;
		MCMCOFO.one_over_nu_TOTAL = SAOFO.one_over_nu_TOTAL;


		MCMCOFO.SSO = SSO;
		//MCMCOFO.SSO.quiet=false;//$$$$$$$$$$$$$$$$$$$$
		//MCMCOFO.SSO.verbose=true;//$$$$$$$$$$$$$$$$$$$$$$$$$
		MCMCOFO.DIt = DIt;
		MCMCOFO.DIdo = DIdo;
		MCMCOFO.RVdo = RVdo;
		MCMCOFO.SYSdo = SYSdo;
		MCMCOFO.earliestEpoch = earliestEpoch;
		MCMCOFO.randSeed = time_nsec;
		MCMCOFO.ODT = ODT2;
		MCMCOFO.saveEachInt = saveEachInt;

		MCMCOFO.paramChangeStepSizePercent = 0.025;
		//MCMCOFO.sigmaPercent = 2.45;

		if (SAOFO.chiSquaredMin*SAOFO.one_over_nu_TOTAL>SSO.chiSquaredMax)
		{
			cout<<"MCMC: ERROR: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cout<<"MCMC: ERROR: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cout<<"MCMC: ERROR: no orbit was found during Simulated Annealing below the chiSquaredMax in order to start a proper MCMC run!"<<endl;
			cout<<"MCMC: ERROR: Lowest reduced chi squared found was "<<SAOFO.chiSquaredMin*SAOFO.one_over_nu_TOTAL<<", and max allowed "<<SSO.chiSquaredMax<<endl;
			cout<<"MCMC: ERROR: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			cout<<"MCMC: ERROR: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
			// load up a string for entire log and write it to file
			string LOGlines;
			LOGlines = SSlog.str();
			GT.logFileWriter(MCMCOFO.ODT.data_filename, LOGlines);
		}
		else
		{
			//*** Run Full MCMC chain with results of the SimAnnealing and sigma tuning
			if ((SSO.quiet==false)&&(SSO.SILENT==false))
					cout<<"Calling MCMCOFO simulator"<<endl;
			MCMCOFO.simulator();
			//SSO.quiet=false;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			if ((SSO.quiet==false)&&(SSO.SILENT==false))
				cout<<"Returned from MCMCOFO simulator :-)"<<endl;
			//move output log string from simulation run to SSlog
			SSlog<<MCMCOFO.SSlogStr;
			//###############################################################################
			//###############################################################################

			int totalAccepted2 = MCMCOFO.ODT.es.size();
			MCMCOFO.ODT.numSamplesAccepted = totalAccepted2;

			string printLine4="";
			// Get all best orbit values
			ss << "\nMCMC: $$$$$$$$$$$$$$$ MCMC SIMULATOR COMPLETE $$$$$$$$$$$$$$$"<<endl;
			ss<< totalAccepted2 <<" orbits were accepted during simulation"<<endl;
			ss<< "\nBest orbit found at step "<<MCMCOFO.bestOrbit<<" :"<<endl;
			ss<<"chiSquareMin_reduced = "<<MCMCOFO.chiSquaredMin*SAOFO.one_over_nu_TOTAL<<endl;
			ss<< "chiSquaredMin = "<< MCMCOFO.chiSquaredMin <<endl;
			ss<< "chiSquaredMin from vector = "<<MCMCOFO.ODT.chiSquareds[MCMCOFO.bestOrbit]<<endl;
			ss<< "mean Acceptance rate = "<<(double(totalAccepted2)/double(SSO.numSamples))<<endl;
			ss<< "LongAN = "<< MCMCOFO.ODT.longAN_degs[MCMCOFO.bestOrbit] <<endl;
			ss<< "e = "<< MCMCOFO.ODT.es[MCMCOFO.bestOrbit] <<endl;
			ss<< "To = "<< fixed <<MCMCOFO.ODT.Ts[MCMCOFO.bestOrbit] <<endl;
			ss<< "Tc = "<< MCMCOFO.ODT.Tcs[MCMCOFO.bestOrbit]  <<endl;
			ss<< "period = "<< MCMCOFO.ODT.periods[MCMCOFO.bestOrbit] <<endl;
			ss<< "inclination = "<< MCMCOFO.ODT.inclination_degs[MCMCOFO.bestOrbit] <<endl;
			ss<< "argPeri = "<< MCMCOFO.ODT.argPeri_degs[MCMCOFO.bestOrbit] <<endl;
			ss<< "a_total = "<< MCMCOFO.ODT.a_totals[MCMCOFO.bestOrbit] <<endl;
			if (MCMCOFO.SSO.DIonly==false)
			{
				ss<< "K = "<< MCMCOFO.ODT.Ks[MCMCOFO.bestOrbit]<<endl;
				for (int set=0;set<MCMCOFO.ODT.RVoffsets[MCMCOFO.bestOrbit].size();++set)
					ss<<"RVoffset for dataset "<<set<<", was = "<< MCMCOFO.ODT.RVoffsets[MCMCOFO.bestOrbit][set]<<endl;
			}

			// calc how long sim took to run and print appropriate duration string
			time_t endTime;
			endTime = time(NULL);
			int timeElapsed = endTime-startTime;

			//cout<<"starting to print time taken string."<<endl;
			string timeString;
			timeString = GT.timeStr( timeElapsed);
			//cout<<"back from timeStr func"<<endl;
			ss<< "Simulated Annealing, Sigma Tuning and MCMC took "<<timeString<<" to complete"<<endl;

			printLine4=ss.str();
			ss.clear();
			ss.str(std::string());
			if (SSO.SILENT==false)
				cout<<printLine4;
			// load up log with all prints in SSlog stringstream
			// then write it to file.
			SSlog<< printLine4;

			// Write output data of MCMC to file
			GT.fileWriter(MCMCOFO.ODT);
			// load up a string for entire log and write it to file
			string LOGlines;
			LOGlines = SSlog.str();
			GT.logFileWriter(MCMCOFO.ODT.data_filename, LOGlines);

			//perform first stage of Gelman-Rubin calculation if requested
			//NOTE: this function clears the memory of the vector, so it
			//      must be performed AFTER ALL other functions that need
			//	    those vectors.
			if ((SSO.useMultiProcessing)&&(SSO.CalcGelmanRubin))
			{
				cout<<"MCMC: $$$$  Starting to calculate Gelman-Rubin statistic, stage 1  $$$$"<<endl;
				GT.gelmanRubinStage1(MCMCOFO.ODT,SSO.numTimesCalcGR);
			}
		}//end attempt at MCMC after finding SimAnneal went well.
	}//end of MCMC stage
	else//Ran just SimAnneal, so just perform basic wrap ups, no statistic calcs.
	{
		// Write output data of Sim Anneal to file
		//string SAFO_data_filename = SAOFO.ODT.data_filename;
		//SAFO_data_filename =  GT.filenamePrepend(SAFO_data_filename, "SimAnneal_");
		//SAOFO.ODT.data_filename = SAFO_data_filename;
		GT.fileWriter(SAOFO.ODT);

		// calc how long sim took to run and print appropriate duration string
		time_t endTime;
		endTime = time(NULL);
		int timeElapsed = endTime-startTime;

		//cout<<"starting to print time taken string."<<endl;
		string timeString;
		timeString = GT.timeStr(timeElapsed);
		//cout<<"back from timeStr func"<<endl;
		ss<< "\nSimulated Annealing and Sigma Tuning took "<<timeString<<" to complete"<<endl;
		string printLine4="";
		printLine4=ss.str();
		ss.clear();
		ss.str(std::string());
		if (SSO.SILENT==false)
			cout<<printLine4;
		// load up log with all prints in SSlog stringstream
		// then write it to file.
		SSlog<< printLine4;
		string LOGlines;
		LOGlines = SSlog.str();
		GT.logFileWriter(SAOFO.ODT.data_filename, LOGlines);
	}//Done wrap up for SimAnneal

	return EXIT_SUCCESS;
}
