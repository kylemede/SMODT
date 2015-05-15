#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include "orbToolboxes.h"
//#include "SimSettingsObj.h" //called in orbToolboxes.h, so don't need to call here again

using namespace std;
/*! \file */

void SimSettingsObj::settingsLoadUp(const char* filename)
{
	/**
	 * This function will read in the file containing the settings for a
	 * Simulated Annealing or MCMC run.  This file has already been checked
	 * within the Python code that initializes a simulation, and at this point
	 * all the settings should be well formated and within realistic values
	 * ranges.
	 */
	bool verboseInternal = false;//for testing set = true
	generalTools GT;

	// Setting to fully open values by default
	// defaults replaced later in this func.
	chiSquaredMax=1e6;
	numSamples=1e6;
	useMultiProcessing = true;
	startTemp = 100.0;
	numSamples_SimAnneal = 1e6;
	numSamplePrints=10;
	SILENT = false;
	quiet = true;
	verbose = false;
	settings_and_InputDataDir = "/mnt/Data1/Todai_Work/Dropbox/workspace/SMODT/settings_and_InputData/";
	SystemDataFilename = "SystemData.txt";
	DIdataFilename = "DIdata.dat";
	RVdataFilename = "RVdata.dat";
	data_dir = "/mnt/Data1/Todai_Work/Data/data_SMODT/";
	data_filenameRoot = "TrempRootFilename.dat";
	RVonly = true;
	DIonly = true;
	mcONLY = true;
	simAnneal = false;
	simulate_StarStarRV = false;
	simulate_StarPlanetRV = false;
	simulate_PrimaryOrbitRV = false;
	calcCorrLengths=false;
	CalcGelmanRubin=true;
	numTimesCalcGR=1000;
	TcStepping = true;
	primaryStarRVs=true;
	TcEqualT=true;
	argPeriPlusRV=false;
	argPeriPlusDI=false;

	// Ranges for acceptable random number inputs ######
	longAN_degMIN = 0.0; // [deg]
	longAN_degMAX = 0.0;// [deg]
	eMIN = 0.0;
	eMAX = 0.0;
	a_totalMIN = 0.0; //[AU]
	a_totalMAX = 0.0; //[AU]
	periodMIN = 0.0; // [yrs]
	periodMAX = 0.0; // [yrs]
	inclination_degMIN = 0.0; // [deg]
	inclination_degMAX = 0; // [deg]
	argPeri_degMIN=0.0; // [deg]
	argPeri_degMAX=0; // [deg]
	T_Min=0;// [JD]
	T_Max=0; //[JD]
	K_MIN=0;// [m/s]
	K_MAX=0;// [m/s]


	std::ifstream infile(filename);
	std::stringstream ss;

	if (verboseInternal)
	{
		cout<<"############################# Loading up settings file #################"<<endl;
		cout<<"\nworking on Simulator Settings file: "<<filename<<"\n"<<endl;
	}

	while(!infile.eof())
	{
		string line;
		getline(infile,line);
		if (line.length()>0)
		{
			if (line[0]=='#')
			{
				// line is a comment line
				if (verboseInternal)
					cout<<"comments/headers were found to be: "<<line<<endl;
			}
			//not a blank line
			int loc;
			loc = line.find_first_of('=');

			if (line[0]!='#' && loc>0)
			{
				if (verboseInternal)
					cout<<"line orig: "<<line<<endl;
				// line is a key and value
				string key = line.substr(0,loc);
				ss<<key;
				ss>>key;
				ss.clear();
				ss.str(std::string());
				string val = line.substr(loc+1,line.length());
				//string boolStr;

				//Run through all current possible input keys
				// First: check if it is a general setting
				if (key.compare("chiSquaredMax")==0)
				{
					ss<<val;
					ss>>chiSquaredMax;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"chiSquaredMax: "<<chiSquaredMax<<endl;
				}
				else if (key.compare("startTemp")==0)
				{
					ss<<val;
					ss>>startTemp;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"startTemp: "<<startTemp<<endl;
				}
				else if (key.compare("numSamples_SimAnneal")==0)
				{
					ss<<val;
					ss>>numSamples_SimAnneal;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"numSamples_SimAnneal: "<<numSamples_SimAnneal<<endl;
				}
				else if (key.compare("numSamples")==0)
				{
					ss<<val;
					ss>>numSamples;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"numSamples: "<<numSamples<<endl;
				}
				else if (key.compare("numTimesCalcGR")==0)
				{
					ss<<val;
					ss>>numTimesCalcGR;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"numTimesCalcGR: "<<numTimesCalcGR<<endl;
				}
				else if (key.compare("numSamplePrints")==0)
				{
					ss<<val;
					ss>>numSamplePrints;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"numSamplePrints: "<<numSamplePrints<<endl;
				}
				else if (key.compare("useMultiProcessing")==0)
				{
					useMultiProcessing = boolParser(val);
					if (verboseInternal)
						cout<<"useMultiProcessing: "<<GT.boolToStr(useMultiProcessing)<<endl;
				}
				else if (key.compare("SILENT")==0)
				{
					SILENT = boolParser(val);
					if (verboseInternal)
						cout<<"SILENT: "<<GT.boolToStr(SILENT)<<endl;
				}
				else if (key.compare("quiet")==0)
				{
					quiet = boolParser(val);
					if (verboseInternal)
						cout<<"quiet: "<<GT.boolToStr(quiet)<<endl;
				}
				else if (key.compare("CalcGelmanRubin")==0)
				{
					CalcGelmanRubin = boolParser(val);
					if (verboseInternal)
						cout<<"CalcGelmanRubin: "<<GT.boolToStr(CalcGelmanRubin)<<endl;
				}
				else if (key.compare("TcStepping")==0)
				{
					TcStepping = boolParser(val);
					if (verboseInternal)
						cout<<"TcStepping: "<<GT.boolToStr(TcStepping)<<endl;
				}
				else if (key.compare("calcCorrLengths")==0)
				{
					calcCorrLengths = boolParser(val);
					if (verboseInternal)
						cout<<"calcCorrLengths: "<<GT.boolToStr(calcCorrLengths)<<endl;
				}
				else if (key.compare("verbose")==0)
				{
					verbose = boolParser(val);
					if (verboseInternal)
						cout<<"verbose: "<<GT.boolToStr(verbose)<<endl;
				}
				else if (key.compare("settings_and_InputDataDir")==0)
				{
					ss<<val;
					ss>>settings_and_InputDataDir;
					ss.clear();
					ss.str(std::string());
					// add '/' if needed to end
					if (settings_and_InputDataDir[settings_and_InputDataDir.length()-1]!='/')
						settings_and_InputDataDir = settings_and_InputDataDir+"/";
					if (verboseInternal)
						cout<<"settings_and_InputDataDir: "<<settings_and_InputDataDir<<endl;
				}
				else if (key.compare("SystemDataFilename")==0)
				{
					ss<<val;
					ss>>SystemDataFilename;
					ss.clear();
					if (verboseInternal)
						cout<<"SystemDataFilename: "<<SystemDataFilename<<endl;
				}
				else if (key.compare("DIdataFilename")==0)
				{
					ss<<val;
					ss>>DIdataFilename;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"DIdataFilename: "<<DIdataFilename<<endl;
				}
				else if (key.compare("RVdataFilename")==0)
				{
					ss<<val;
					ss>>RVdataFilename;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"RVdataFilename: "<<RVdataFilename<<endl;
				}
				else if (key.compare("outputData_dir")==0)
				{
					ss<<val;
					ss>>data_dir;
					ss.clear();
					ss.str(std::string());
					if (data_dir[data_dir.length()-1]!='/')
						data_dir = data_dir+"/";
					if (verboseInternal)
						cout<<"data_dir: "<<data_dir<<endl;
				}
				else if (key.compare("outputData_filenameRoot")==0)
				{
					ss<<val;
					ss>>data_filenameRoot;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"data_filenameRoot: "<<data_filenameRoot<<endl;
				}
				else if (key.compare("RVonly")==0)
				{

					RVonly = boolParser(val);
					if (verboseInternal)
						cout<<"RVonly: "<<GT.boolToStr(RVonly)<<endl;
				}
				else if (key.compare("DIonly")==0)
				{

					DIonly = boolParser(val);
					if (verboseInternal)
						cout<<"DIonly: "<<GT.boolToStr(DIonly)<<endl;
				}
				else if (key.compare("mcONLY")==0)
				{
					mcONLY = boolParser(val);
					if (verboseInternal)
						cout<<"mcONLY: "<<GT.boolToStr(mcONLY)<<endl;
				}
				else if (key.compare("simAnneal")==0)
				{
					simAnneal = boolParser(val);
					if (verboseInternal)
						cout<<"simAnneal: "<<GT.boolToStr(simAnneal)<<endl;
				}
				else if (key.compare("simulate_StarStarRV")==0)
				{
					simulate_StarStarRV = boolParser(val);
					if (verboseInternal)
						cout<<"simulate_StarStarRV: "<<GT.boolToStr(simulate_StarStarRV)<<endl;
				}
				else if (key.compare("simulate_StarPlanetRV")==0)
				{
					simulate_StarPlanetRV = boolParser(val);
					if (verboseInternal)
						cout<<"simulate_StarPlanetRV: "<<GT.boolToStr(simulate_StarPlanetRV)<<endl;
				}
				else if (key.compare("simulate_PrimaryOrbitRV")==0)
				{
					simulate_PrimaryOrbitRV = boolParser(val);
					if (verboseInternal)
						cout<<"simulate_PrimaryOrbitRV: "<<GT.boolToStr(simulate_PrimaryOrbitRV)<<endl;
				}
				else if (key.compare("primaryStarRVs")==0)
				{
					primaryStarRVs = boolParser(val);
					if (verboseInternal)
						cout<<"primaryStarRVs: "<<GT.boolToStr(primaryStarRVs)<<endl;
				}
				else if (key.compare("TcEqualT")==0)
				{
					TcEqualT = boolParser(val);
					if (verboseInternal)
						cout<<"TcEqualT: "<<GT.boolToStr(TcEqualT)<<endl;
				}
				else if (key.compare("argPeriPlusRV")==0)
				{
					ss<<val;
					ss>>argPeriPlusRV;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"argPeriPlusRV: "<<argPeriPlusRV<<endl;
				}
				else if (key.compare("argPeriPlusDI")==0)
				{
					ss<<val;
					ss>>argPeriPlusDI;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"argPeriPlusDI: "<<argPeriPlusDI<<endl;
				}

				// Next: check if it is a min/max setting
				else if (key.compare("longAN_degMIN")==0)
				{
					ss<<val;
					ss>>longAN_degMIN;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"longAN_degMIN: "<<longAN_degMIN<<endl;
				}
				else if (key.compare("longAN_degMAX")==0)
				{
					ss<<val;
					ss>>longAN_degMAX;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"longAN_degMAX: "<<longAN_degMAX<<endl;
				}
				else if (key.compare("eMIN")==0)
				{
					ss<<val;
					ss>>eMIN;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"eMIN: "<<eMIN<<endl;
				}
				else if (key.compare("eMAX")==0)
				{
					ss<<val;
					ss>>eMAX;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"eMAX: "<<eMAX<<endl;
				}
				else if (key.compare("a_totalMIN")==0)
				{
					ss<<val;
					ss>>a_totalMIN;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"a_totalMIN: "<<a_totalMIN<<endl;
				}
				else if (key.compare("a_totalMAX")==0)
				{
					ss<<val;
					ss>>a_totalMAX;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"a_totalMAX: "<<a_totalMAX<<endl;
				}
				else if (key.compare("periodMIN")==0)
				{
					ss<<val;
					ss>>periodMIN;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"periodMIN: "<<periodMIN<<endl;
				}
				else if (key.compare("periodMAX")==0)
				{
					ss<<val;
					ss>>periodMAX;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"periodMAX: "<<periodMAX<<endl;
				}
				else if (key.compare("inclination_degMIN")==0)
				{
					ss<<val;
					ss>>inclination_degMIN;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"inclination_degMIN: "<<inclination_degMIN<<endl;
				}
				else if (key.compare("inclination_degMAX")==0)
				{
					ss<<val;
					ss>>inclination_degMAX;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"inclination_degMAX: "<<inclination_degMAX<<endl;
				}
				else if (key.compare("argPeri_degMIN")==0)
				{
					ss<<val;
					ss>>argPeri_degMIN;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"argPeri_degMIN: "<<argPeri_degMIN<<endl;
				}
				else if (key.compare("argPeri_degMAX")==0)
				{
					ss<<val;
					ss>>argPeri_degMAX;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"argPeri_degMAX: "<<argPeri_degMAX<<endl;
				}
				else if (key.compare("T_Min")==0)
				{
					ss<<val;
					ss>>T_Min;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"T_Min: "<<T_Min<<endl;
				}
				else if (key.compare("T_Max")==0)
				{
					ss<<val;
					ss>>T_Max;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"T_Max: "<<T_Max<<endl;
				}
				else if (key.compare("K_MIN")==0)
				{
					ss<<val;
					ss>>K_MIN;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"K_MIN: "<<K_MIN<<endl;
				}
				else if (key.compare("K_MAX")==0)
				{
					ss<<val;
					ss>>K_MAX;
					ss.clear();
					ss.str(std::string());
					if (verboseInternal)
						cout<<"K_MAX: "<<K_MAX<<endl;
				}
				else if (key.compare("RVoffsetMAXs")==0)
				{
					RVoffsetMAXs = rvOffsetsParser(val);
					if (verboseInternal)
						cout<<"RVoffsetMAXs found to be of length: "<<RVoffsetMAXs.size()<<endl;
				}
				else if (key.compare("RVoffsetMINs")==0)
				{
					RVoffsetMINs = rvOffsetsParser(val);
					if (verboseInternal)
						cout<<"RVoffsetMINs found to be of length: "<<RVoffsetMINs.size()<<endl;
				}
			}//finished loading in this param's value
		}//line loaded as comment or param
	}//reached end of file, close it and end func
	if (verboseInternal)
		cout<<"############################# DONE loading up settings file #################"<<endl;
	//figure out the offsets for the argument of periapsis for both the DI and RV models
	findArgPeriOffsets();
	infile.close();
}
void SimSettingsObj::findArgPeriOffsets()
{
	argPeriOffsetDI = 0.0;
	argPeriOffsetRV = 0.0;
	//first using RV special bools
	if ((primaryStarRVs)&&(simulate_PrimaryOrbitRV))
		argPeriOffsetDI=-180.0;
	else if ((primaryStarRVs)&&(simulate_PrimaryOrbitRV==false))
		argPeriOffsetRV=180.0;
	//now update due to fixed argPeriPlus values
	argPeriOffsetRV+=argPeriPlusRV;
	argPeriOffsetDI+=argPeriPlusDI;
}

vector<double> rvOffsetsParser(string strIn)
{
	std::stringstream ss;

	bool verbose=false; //for testing set = true

	vector<double> offsets_vec;
	int offsetNum = 0;
	bool reachedEnd = false;
	double latestOffsetDbl = 0;
	string latestOffsetStr = "";

	// a crazy nested char iterator style parser as I don't know how
	// to do this more eloquently yet...

	for (int pos=0;pos<int(strIn.length());++pos)
	{
		if (verbose)
			cout<<"Char at "<<pos<<" = "<<strIn[pos]<<endl;

		if (!reachedEnd)
		{
			if (strIn[pos]=='[')
			{
				offsetNum=1;
				if (verbose)
					cout<<"Found first square bracket"<<endl;
			}
			else if (strIn[pos]==',')
			{
				++offsetNum;
				if (verbose)
					cout<<"Found comma"<<endl;
				ss<<latestOffsetStr;
				ss>>latestOffsetDbl;
				ss.clear();
				ss.str(std::string());
				offsets_vec.push_back(latestOffsetDbl);
				latestOffsetStr="";
				if (verbose)
					cout<<"latestOffsetDbl: "<<latestOffsetDbl<<endl;

			}
			else if (strIn[pos]==']')
			{
				reachedEnd = true;
				if (verbose)
					cout<<"Found second square bracket"<<endl;
				ss<<latestOffsetStr;
				ss>>latestOffsetDbl;
				ss.clear();
				ss.str(std::string());
				offsets_vec.push_back(latestOffsetDbl);
				latestOffsetStr="";
				if (verbose)
					cout<<"latestOffsetDbl: "<<latestOffsetDbl<<endl;

			}
			else
				latestOffsetStr = latestOffsetStr+strIn[pos];
		}
		else
			cout<<"Warning: There was a character after the ending bracket ]. It was :"<<strIn[pos]<<endl;
	}

	return offsets_vec;
}

bool boolParser(string sIn)
{
	bool verbose=false; //for testing set = true

	bool returnBool=false;

	if (verbose)
		cout<<"Input string to boolParser was:"<<sIn<<endl;

	if ((sIn.find("true")>=0&&sIn.find("true")<100)||(sIn.find("True")>=0&&sIn.find("True")<100))
	{
		returnBool=true;
		if (verbose)
			cout<<"string found to contain true, so set bool accordingly."<<endl;

	}
	else if((sIn.find("false")>=0&&sIn.find("false")<100)||(sIn.find("False")>=0&&sIn.find("False")<100))
	{
		returnBool=false;
		if (verbose)
			cout<<"string found to contain false, so set bool accordingly."<<endl;
	}

	return returnBool;
}
