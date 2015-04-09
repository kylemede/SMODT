#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
//#include <string>
//#include <vector>
#include <math.h>
#include <time.h>
//#include <fstream>
#include "toolboxes/orbToolboxes.h" //Both the DI and RV tools headers are called in here, so it is an all in one toolbox header call
#include "simAnnealOrbSimFunc.h"

using namespace std;

int main(int argc ,char *argv[])
{
	/**
	This is the "main" that performs the Simulated Annealing process for a
	single chain of a multi or single chain/thread/process Simulated Annealing simulation
	started by the Python file BinaryOrbSimStarterDuo.py and controlled with the
	MCMC_ProcessManagerDuo.  It will perform the set up tasks, call the
	simAnnealOrbSimFunc to perform the Simulated Annealing for the
	chain based on the input settings.  After the chain process is complete
	it will perform the requested wrap up and statistical calculations
	requested then write all the output files to disk for the wrap up tasks and
	plotting done by Python after all the chains are complete.  The
	@author Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
	*/
	std::stringstream ss;
	std::stringstream SSlog;
	generalTools GT;

	// print to indicate sim has started
	ss<<"\n*** C++ Simulated Annealing simulation has started ***\n";

	// pull in simulator settings filename and output data filename
	// from command line args
	string settingsFilename;
	string outputDataFilename;
	ss<<"There were, "<<argc<<", arguments provided.  The were:"<<endl;
	if (argc>=2)
	{
		settingsFilename = argv[1];//including directory!
		ss<<"Using settings file: "<<settingsFilename<<endl;
	}
	if (argc==3)
	{
		outputDataFilename = argv[2];//including directory!
		ss<<"Output data will be written to: "<<outputDataFilename<<endl;
	}
	if (argc>3)
	{
		ss<<"simAnnealOrbSimulator is only set-up to handle settings";
		ss<<"filename and output data filename command line ";
		ss<<"arguments so far."<<endl;
	}
	else if (argc==1)
	{
		settingsFilename = "/run/media/Kyle/Data1/Todai_Work/Dropbox/workspace/Binary-project/SimSettings_and_InputData/SimSettings.txt";
		ss<<"No settings filename provided, so using default: ";
		ss<<settingsFilename<<endl;
		ss<<"No output data filename provided through command line";
		ss<<" so, value within settings file or default being used."<<endl;
	}

	//cout<<"Command line arguments all loaded"<<endl; //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	// instantiate the settings object and load it up
	SimSettingsObj SSO;
	SSO.settingsLoadUp(settingsFilename.c_str());

	//cout<< "Settings obj all loaded up"<<endl; //$$$$$$$$$$$$$$

	// instantiate and load up DI and RV data objects as needed
	string SystemDataFilename = SSO.settings_and_InputDataDir + SSO.SystemDataFilename;
	ss<<"SystemDataFilename = "<< SystemDataFilename <<endl;
	DItools DIt;
	DIdataObj DIdo;
	RVdataObj RVdo;
	DataObj SYSdo;
	//cout<<"SYSdo object instantiated "<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	// log and maybe print all inputs to C++ call
	string inputsStr;
	inputsStr = ss.str();
	ss.clear();
	ss.str(std::string());
	SSlog<<inputsStr;
	if (SSO.verbose)
		cout<< inputsStr;

	SYSdo.systemDataLoadUp(SystemDataFilename.c_str());
	//cout<<"SYSdo object loaded up "<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	//cout<<"DI and RV data objects instantiated"<<endl; //$$$$$$$$$$$$$$$$$$$$$$$
	if ((SSO.DIonly==true) or (SSO.RVonly==false && SSO.DIonly==false))
	{
		//cout<<"Starting to load up DI data obj"<<endl; //$$$$$$$$$$$$$$$$$
		string DIdataFilename = SSO.settings_and_InputDataDir + SSO.DIdataFilename;
		DIdo.dataLoadUp(DIdataFilename.c_str());
		DIdo.systemDataLoadUp(SystemDataFilename.c_str());
		//cout<<"DIdo sys data loaded"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$
		// instantiate DI tools obj and load it up with same values
		DIt = GT.DItoolsParamLoadUp(DIdo);
		//cout<<"DIt loaded up "<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
	double earliestEpoch = GT.earliestEpochFinder(DIdo, RVdo);

	//Start seed for random number generator
	//create very high precision seed for generator (nanosecond precision)
	timespec time1;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
	int time_nsec=time1.tv_nsec;
	//cout<<"Starting time in nanoseconds for use as a random number generator seed = "<<time_nsec<<endl; //$$$$$$$$$$$$$$

	// ******* STARTING SIMULATION!!!!  *************
	time_t startTime;
	startTime = time(NULL);

	string starterString;
	string numSamplesString =  GT.numSamplesStringMaker(SSO.numSamples);
	ss<<"\nSimAnneal: $$$$$$$$$$$$$$$$$$$  STARTING THE SIMULATOR  $$$$$$$$$$$$$$$$$" <<endl;
	ss<<"Number of sample orbits being created = " << numSamplesString <<endl;
	starterString = ss.str();
	ss.clear();
	ss.str(std::string());
	cout<< starterString;
	SSlog<<starterString;

	outputDataType ODT;
	ODT.data_filename = outputDataFilename;
	//cout<<"outputDataType instantiated"<<endl;//$$$$$$$$$$$$$$$$$$$$$
	//### Find out how often to save accepted steps to output data array ###
	int saveEachInt;
	if (SSO.numSamples<100001)
		saveEachInt=10;
	else if ((SSO.numSamples>100000)&&(SSO.numSamples<10000000))
		saveEachInt=100;
	else
		saveEachInt=1000;

	//###############################################################################
	//
	//Instantiate the simAnneal object and load it up
	//
	simAnealOrbFuncObj SAOFO;
	SAOFO.SSO = SSO;
	SAOFO.DIt = DIt;
	SAOFO.DIdo = DIdo;
	SAOFO.RVdo = RVdo;
	SAOFO.SYSdo = SYSdo;
	SAOFO.earliestEpoch = earliestEpoch;
	SAOFO.randSeed = time_nsec;
	SAOFO.ODT = ODT;
	SAOFO.numSamples_SA = SSO.numSamples;
	SAOFO.tempStepSizePercent = 0.025;
	SAOFO.startTemp = SSO.startTemp;
	SAOFO.sigmaPercent = 1.0;
	SAOFO.saveEachInt = saveEachInt;

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
		cout<< startParmsString;
		SSlog<<startParmsString;
	}
	// Run simulator
	if (SSO.silent==false)
			cout<<"Calling SAOFO simulator"<<endl;
	SAOFO.simulator();
	cout<<"Returned from SAOFO simulator :-)"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	if (SSO.silent==false)
		cout<<"Returned from SAOFO simulator :-)"<<endl;
	//move output log string from simulation run to SSlog
	SSlog<<SAOFO.SSlogStr;
	//###############################################################################


    int totalAccepted = SAOFO.ODT.es.size();
//    if (acceptedCounter!=totalAccepted)
//    	cout<<"Warning: (acceptedCounter) "<<acceptedCounter<< " != "<<totalAccepted<<" (totalAccepted)"<<endl;

    SAOFO.ODT.numSamplesAccepted = totalAccepted;

    string printLine3="";
    // Get all best orbit values
    //cout<<"line # 185"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ss << "\n$$$$$$$$$$$$$$$ SIMULATOR COMPLETE $$$$$$$$$$$$$$$"<<endl;
//    ss<< totalAccepted <<" orbits were accepted during simulation"<<endl;
//    ss<< "timesBeenHereTotal = "<<SAOFO.timesBeenHereTotal<<endl;
//    int sum = GT.sumIntCalc(SAOFO.ODT.timesBeenHeres,SAOFO.ODT.timesBeenHeres.size());
//    ss<< "timesBeenHeres total = "<<sum<<endl;
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

	if (false)
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
//		ss<< "timesBeenHere = "<<SAOFO.ODT.timesBeenHeres.back() <<endl;
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
//		ss<< "timesBeenHere = "<<SAOFO.ODT.timesBeenHeres[totalAccepted-2] <<endl;
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

	//cout<<"line # 277"<<endl;//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	// calc how long sim took to run
    time_t endTime;
    endTime = time(NULL);
    int timeElapsed = endTime-startTime;

    if ( timeElapsed < 60.0 )
    	ss<< "\nSimAnneal: Simulator took "<<timeElapsed << " seconds to complete"<<endl;
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
    	ss <<"\nSimAnneal: Simulator took "<<min <<" minutes and "<<sec<< " seconds to complete"<<endl;
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
	GT.logFileWriter(SAOFO.ODT.data_filename, LOGlines);

	// get output filename and write output data to it
    GT.fileWriter(SAOFO.ODT);

    cout<<"\nSimAnneal: $$$$$ Data file written, finished C++ part of simulation.$$$$$\n\n"<<endl;
	return EXIT_SUCCESS;
}
