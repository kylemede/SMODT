#include <iostream>
#include <iomanip>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include "orbToolboxes.h"
//#include "DataObj.h" //called in orbToolboxes.h, so don't need to call here again
/*! \file */
using namespace std;


void DataObj::systemDataLoadUp(const char* filename)
{
	/**
	 * This function will load up the apriori data for the system being investigated.
	 * It will check that the values are within an realistic range to avoid problems
	 * during the simulation.
	 */
	bool verbose = false;//for testing set = true

	// General system data
	Sys_Dist_PC=0;
	Sys_Dist_PC_error=0;
	Mass1=0;
	Mass1_error=0;
	// Orbital params for planet
	// (for if investigating companion star's orbit)
	planet_K = 0;
	planet_K_error = 0;
	planet_e = 0;
	planet_e_error = 0;
	planet_T = 0;
	planet_T_error = 0;
	planet_Tc = 0;
	planet_Tc_error = 0;
	planet_P = 0;
	planet_P_error = 0;
	planet_MsinI = 0;
	planet_MsinI_error = 0;
	planet_argPeri = 90;
	planet_argPeri_error = 0;
	planet_inc = 0;
	planet_inc_error = 0;
	planet_long_AN=0;

	// Orbital params for companion star
	// (for if investigating orbit of star around primary star)
	star_K = 0;
	star_K_error = 0;
	star_e = 0;
	star_e_error = 0;
	star_T = 0;
	star_T_error = 0;
	star_Tc = 0;
	star_Tc_error = 0;
	star_P = 0;
	star_P_error = 0;
	star_Mass2 = 0;
	star_Mass2_error = 0;
	star_argPeri = 90;
	star_argPeri_error = 0;
	star_inc = 0;
	star_inc_error = 0;
	star_long_AN =0;

	std::ifstream infile(filename);
	std::stringstream ss;

	ss<<std::setprecision(15);
	cout<<std::setprecision(15);

	if (verbose)
		cout<<"\nworking on General data file: "<<filename<<"\n"<<endl;

	while(!infile.eof())
	{
		string line;
		getline(infile,line);
		if (line.length()>0)
		{
			if (line[0]=='#')
			{
				// line is a comment line
				if (verbose)
					cout<<"comments/headers were found to be: "<<line<<endl;
			}
			//not a blank line
			int loc;
			loc = line.find_first_of('=');

			if (line[0]!='#' && loc>0)
			{
				if (verbose)
					cout<<"line orig: "<<line<<endl;
				// line is a key and value
				string key = line.substr(0,loc);
				ss<<key;
				ss>>key;
				ss.clear();
				ss.str(std::string());
				string val = line.substr(loc+1,line.length());

				//Run through all current possible input keys

				// General system params and mass of primary star
				if ((key.compare("Sys_Dist_PC")==0)&&(key.compare("Sys_Dist_PC_error")!=0))
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					// ALL the 0<valOut<1000 blocks below are extra
					// any text that is not integers will be
					// turned into zero by default, but these
					// blocks are just in case.
					if (valOut>1000 or valOut<0)
						valOut = Sys_Dist_PC;
					else
						Sys_Dist_PC = valOut;
					if (verbose)
						cout<<"Sys_Dist_PC: "<<Sys_Dist_PC<<endl;
				}
				else if (key.compare("Sys_Dist_PC_error")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					// ALL the 0<valOut<1000 blocks below are extra
					// any text that is not integers will be
					// turned into zero by default, but these
					// blocks are just in case.
					if (valOut>1000 or valOut<0)
						valOut = Sys_Dist_PC_error;
					else
						Sys_Dist_PC_error = valOut;
					if (verbose)
						cout<<"Sys_Dist_PC_error: "<<Sys_Dist_PC_error<<endl;
				}
				else if (key.compare("Mass1")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = Mass1;
					else
						Mass1 = valOut;
					if (verbose)
						cout<<"Mass1: "<<Mass1<<endl;
				}
				else if (key.compare("Mass1_error")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = Mass1_error;
					else
						Mass1_error = valOut;
					if (verbose)
						cout<<"Mass1_error: "<<Mass1_error<<endl;
				}
				// Orbital params for the well studied possible orbiting planet
				else if ((key.compare("planet_K")==0)&&(key.compare("planet_K_error")!=0))
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = planet_K;
					else
						planet_K = valOut;
					if (verbose)
						cout<<"planet_K: "<<planet_K<<endl;
				}
				else if (key.compare("planet_K_error")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = planet_K_error;
					else
						planet_K_error = valOut;
					if (verbose)
						cout<<"planet_K_error: "<<planet_K_error<<endl;
				}
				else if ((key.compare("planet_e")==0)&&(key.compare("planet_e_error")!=0))
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = planet_e;
					else
						planet_e = valOut;
					if (verbose)
						cout<<"planet_e: "<<planet_e<<endl;
				}
				else if ((key.compare("planet_long_AN")==0)&&(key.compare("planet_long_AN_error")!=0))
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = planet_long_AN;
					else
						planet_long_AN = valOut;
					if (verbose)
						cout<<"planet_long_AN: "<<planet_long_AN<<endl;
				}
				else if (key.compare("planet_e_error")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = planet_e_error;
					else
						planet_e_error = valOut;
					if (verbose)
						cout<<"planet_e_error: "<<planet_e_error<<endl;
				}
				else if ((key.compare("planet_T")==0)&&(key.compare("planet_T_error")!=0))
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>3000000 or valOut<0)
						valOut = planet_T;
					else
						planet_T = valOut;
					if (verbose)
						cout<<"planet_T: "<<planet_T<<endl;
				}
				else if (key.compare("planet_T_error")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>3000000 or valOut<0)
						valOut = planet_T_error;
					else
						planet_T_error = valOut;
					if (verbose)
						cout<<"planet_T_error: "<<planet_T_error<<endl;
				}
				else if ((key.compare("planet_Tc")==0)&&(key.compare("planet_Tc_error")!=0))
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>3000000 or valOut<0)
						valOut = planet_Tc;
					else
						planet_Tc = valOut;
					if (verbose)
						cout<<"planet_Tc: "<<planet_Tc<<endl;
				}
				else if (key.compare("planet_Tc_error")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>3000000 or valOut<0)
						valOut = planet_Tc_error;
					else
						planet_Tc_error = valOut;
					if (verbose)
						cout<<"planet_Tc_error: "<<planet_Tc_error<<endl;
				}
				else if ((key.compare("planet_P")==0)&&(key.compare("planet_P_error")!=0))
				{
					double valOut;
					//double tempOut = 0.009069140970019346;
					//ss.precision(15);
					ss <<val;
					ss>> valOut;
					cout<<std::setprecision(15);//<<tempOut<<endl;
					ss.clear();
					ss.str(std::string());
					if (valOut>100000 or valOut<0.00000001)
						valOut = planet_P;
					else
						planet_P = valOut;
					if (verbose)
						cout<<"planet_P: "<<planet_P<<endl;
				}
				else if (key.compare("planet_P_error")==0)
				{
					double valOut;
					ss <<val;
					ss>> valOut;
					cout<<std::setprecision(15);
					ss.clear();
					ss.str(std::string());
					if (valOut>100000 or valOut<0.00000001)
						valOut = planet_P_error;
					else
						planet_P_error = valOut;
					if (verbose)
						cout<<"planet_P_error: "<<planet_P_error<<endl;
				}
				else if (key.compare("planet_MsinI")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = planet_MsinI;
					else
						planet_MsinI = valOut;
					if (verbose)
						cout<<"planet_MsinI: "<<planet_MsinI<<endl;
				}
				else if (key.compare("planet_MsinI_error")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = planet_MsinI_error;
					else
						planet_MsinI_error = valOut;
					if (verbose)
						cout<<"planet_MsinI_error: "<<planet_MsinI_error<<endl;
				}
				else if ((key.compare("planet_argPeri")==0)&&(key.compare("planet_argPeri_error")!=0))
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = planet_argPeri;
					else
						planet_argPeri = valOut;
					if (verbose)
						cout<<"planet_argPeri: "<<planet_argPeri<<endl;
				}
				else if (key.compare("planet_argPeri_error")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = planet_argPeri_error;
					else
						planet_argPeri_error = valOut;
					if (verbose)
						cout<<"planet_argPeri_error: "<<planet_argPeri_error<<endl;
				}
				else if ((key.compare("planet_inc")==0)&&(key.compare("planet_inc_error")!=0))
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = planet_inc;
					else
						planet_inc = valOut;
					if (verbose)
						cout<<"planet_inc: "<<planet_inc<<endl;
				}
				else if (key.compare("planet_inc_error")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = planet_inc_error;
					else
						planet_inc_error = valOut;
					if (verbose)
						cout<<"planet_inc_error: "<<planet_inc_error<<endl;
				}
				// Orbital params for the possible well studied companion star
				else if ((key.compare("star_K")==0)&&(key.compare("star_K_error")!=0))
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>100000 or valOut<0)
						valOut = star_K;
					else
						star_K = valOut;
					if (verbose)
						cout<<"star_K: "<<star_K<<endl;
				}
				else if (key.compare("star_K_error")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>100000 or valOut<0)
						valOut = star_K_error;
					else
						star_K_error = valOut;
					if (verbose)
						cout<<"star_K_error: "<<star_K_error<<endl;
				}
				else if ((key.compare("star_e")==0)&&(key.compare("star_e_error")!=0))
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = star_e;
					else
						star_e = valOut;
					if (verbose)
						cout<<"star_e: "<<star_e<<endl;
				}
				else if (key.compare("star_e_error")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = star_e_error;
					else
						star_e_error = valOut;
					if (verbose)
						cout<<"star_e_error: "<<star_e_error<<endl;
				}
				else if ((key.compare("star_T")==0)&&(key.compare("star_T_error")!=0))
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>3000000 or valOut<0)
						valOut = star_T;
					else
						star_T = valOut;
					if (verbose)
						cout<<"star_T: "<<star_T<<endl;
				}
				else if (key.compare("star_T_error")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>3000000 or valOut<0)
						valOut = star_T_error;
					else
						star_T_error = valOut;
					if (verbose)
						cout<<"star_T_error: "<<star_T_error<<endl;
				}
				else if ((key.compare("star_P")==0)&&(key.compare("star_P_error")!=0))
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>100000 or valOut<0)
						valOut = star_P;
					else
						star_P = valOut;
					if (verbose)
						cout<<"star_P: "<<star_P<<endl;
				}
				else if (key.compare("star_P_error")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>100000 or valOut<0)
						valOut = star_P_error;
					else
						star_P_error = valOut;
					if (verbose)
						cout<<"star_P_error: "<<star_P_error<<endl;
				}
				else if ((key.compare("star_Mass2")==0)&&(key.compare("star_Mass2_error")!=0))
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = star_Mass2;
					else
						star_Mass2 = valOut;
					if (verbose)
						cout<<"star_Mass2: "<<star_Mass2<<endl;
				}
				else if (key.compare("star_Mass2_error")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = star_Mass2_error;
					else
						star_Mass2_error = valOut;
					if (verbose)
						cout<<"star_Mass2_error: "<<star_Mass2_error<<endl;
				}
				else if ((key.compare("star_argPeri")==0)&&(key.compare("star_argPeri_error")!=0))
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = star_argPeri;
					else
						star_argPeri = valOut;
					if (verbose)
						cout<<"star_argPeri: "<<star_argPeri<<endl;
				}
				else if (key.compare("star_argPeri_error")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = star_argPeri_error;
					else
						star_argPeri_error = valOut;
					if (verbose)
						cout<<"star_argPeri_error: "<<star_argPeri_error<<endl;
				}
				else if ((key.compare("star_long_AN")==0)&&(key.compare("star_long_AN_error")!=0))
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = star_long_AN;
					else
						star_long_AN = valOut;
					if (verbose)
						cout<<"star_long_AN: "<<star_long_AN<<endl;
				}
				else if ((key.compare("star_inc")==0)&&(key.compare("star_inc_error")!=0))
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = star_inc;
					else
						star_inc = valOut;
					if (verbose)
						cout<<"star_inc: "<<star_inc<<endl;
				}
				else if (key.compare("star_inc_error")==0)
				{
					double valOut;
					ss<<val;
					ss>>valOut;
					ss.clear();
					ss.str(std::string());
					if (valOut>1000 or valOut<0)
						valOut = star_inc_error;
					else
						star_inc_error = valOut;
					if (verbose)
						cout<<"star_inc_error: "<<star_inc_error<<endl;
				}
			}//finished loading in this param's value
		}//line loaded as comment or param
	}//reached end of file, close it and end func

	if (verbose)
		cout<<"Done loading up system data into obj, so closing file and returning"<<endl;
	infile.close();
}

void DIdataObj::dataLoadUp(const char* filename)
{
	// Data must be in the collumns:
    // obsDate[JD]  PA[deg]  PA_error[deg]  SA[arcsec]  SA_error[arcsec]
	bool verbose = false;//for testing set = true
	numEpochs_DI = 0;

	std::ifstream infile(filename);
	std::stringstream ss;

	if (verbose)
		cout<<"\nworking on DI data file: "<<filename<<"\n"<<endl;

	while(!infile.eof())
	{
		string line;
		getline(infile,line);
		if (line.length()==0)
		{
			// line is blank
			if (verbose)
				cout<<"Found a blank line"<<endl;
		}
		else if (line[0]=='#')
		{
			// line is a comment line
			if (verbose)
				cout<<"headers were found to be: "<<line<<endl;
		}
		else
		{
			double obsDate;
			double PA;
			double PA_error;
			double SA;
			double SA_error;

			// line has data, so convert it
			ss<<line;
			ss>>obsDate>>PA>>PA_error>>SA>>SA_error;
			if (verbose)
				cout<<"\n"<<obsDate<<" "<<PA<<" "<<PA_error<< " "<<SA<<" "<<SA_error<<endl;
			ss.clear(); //clean out stream for next time
			ss.str(std::string());

			if (verbose)
			{
				cout<<"obsDate = "<< obsDate<<endl;
				cout<<"PA = "<<PA <<endl;
				cout<<"PA_error = "<<PA_error <<endl;
				cout<<"SA = "<< SA<<endl;
				cout<<"SA_error = "<<SA_error <<endl;
			}

			// push these data into the data vectors
			// and increment number of epoch count
			epochs_DI.push_back(obsDate);
			PAs_deg_observed.push_back(PA);
			PA_errors.push_back(PA_error);
			SAs_arcsec_observed.push_back(SA);
			SA_errors.push_back(SA_error);
			++numEpochs_DI;
		}
	}//reached end of file, close it and end func
	infile.close();

	if (verbose)
	{
		cout<<"*Finished loading all the data from the file*"<<endl;
		cout<<"\nResulting inputs were:"<<endl;
		// print final resulting vector sizes
		cout<<"\nepochs_DI size: "<<epochs_DI.size()<<endl;
		cout<<"PAs_deg_observed size: "<<PAs_deg_observed.size()<<endl;
		cout<<"PA_errors size: "<<PA_errors.size()<<endl;
		cout<<"SAs_arcsec_observed size: "<<SAs_arcsec_observed.size()<<endl;
		cout<<"SAs_arcsec_observed size: "<<SAs_arcsec_observed.size()<<endl;
		cout<<"The total number of epochs found was "<<numEpochs_DI<<"\n"<<endl;
	}
}//end dataLoadUp



void RVdataObj::dataLoadUp(const char* filename)
{
	// Data must be in the columns:
	// obsDate[JD]  RV[m/s]  RV_error[m/s]
	// jitter can be included as a 4th column value for the first line of each dataset.

	bool verbose = false;//for testing set = true
	bool silent = true;//for level 2 testing set = false
	int datasetNum=0;
	int lineNumber = 0;
	std::ifstream infile(filename);
	std::stringstream ss;
	double jitter;
	if (verbose)
		cout<<"\nworking on RV data file: "<<filename<<"\n"<<endl;

	while(!infile.eof())
	{
		++lineNumber;
		string line;
		getline(infile,line);
		if (line.length()==0)
		{
			jitter=0;
			// line is blank
			if (!silent)
				cout<<"Found a blank line"<<endl;
		}
		else if (line[0]=='#')
		{
			jitter=0;
			// line is a comment line
			//cout<<"headers were found to be: "<<line<<endl;

		}
		else if (line[0]!='#' && line.length()>1)
		{
			// line has data, so convert it
			int numEpochs_cur=0;
			++datasetNum;
			if (!silent)
				cout<<"Now starting to run through data set # "<<datasetNum<<endl;
			//Reached a data line, so read it in and then repeat that till
			//the end of this set of data is indicated by a blank or comment line
			vector<double> epochs_RV_curr;
			vector<double> RVs_cur;
			vector<double> RV_inv_var_cur;
			double obsDate;
			double RV;
			double RV_error;
			double errSquared;

			// NOTE: first data line of all data sets might include a jitter term,
			// NOTE: so it has to added in quadrature with the error.
			// NOTE: If it isn't there, it will be assumed to be zero.

			ss<<line;
			ss>>obsDate>>RV>>RV_error>>jitter;
			if (!silent)
				cout<<obsDate<<" "<<RV<<" "<<RV_error<<" "<<jitter<<endl;
			ss.clear(); //clean out stream for next time
			ss.str(std::string());

			// push these data into the data vectors
			// and increment number of epoch count
			if (obsDate>1)
			{
				epochs_RV_curr.push_back(obsDate);
				RVs_cur.push_back(RV);
				errSquared = RV_error*RV_error+jitter*jitter;
				double inv_var = 1.0/errSquared;
				if (!silent)
				{
					cout<<"jitter found to be = "<<jitter<<endl;
					cout<<"errSquared calculated to be = "<<errSquared<<endl;
					cout<<"inv_var calculated to be = "<<inv_var<<endl;
				}
				RV_inv_var_cur.push_back(inv_var);
				if (errSquared==0 or inv_var==0)
				{
					cout<<"!! errSquared or inv_var is zero!!"<<endl;
					cout<<"Line number "<<lineNumber<<" was: "<<line<<endl;
					cout<<"This was parsed as: obsDate = "<< obsDate<<", RV = "<<RV<<", RV_error = "<<RV_error<<endl;
					cout<<"These were calculated into: errSquared = "<<errSquared<<", and inv_var = "<<inv_var<<endl;
				}
				++numEpochs_cur;
			}

			bool isdataline = true;
			// run through following lines till end of this data set is reached
			while((!infile.eof())&&(isdataline))
			{
				++lineNumber;
				// perform second round of line testing
				//string line;
				getline(infile,line);
				if (line.length()==0)
				{
					jitter=0;
					// line is blank
					if (!silent)
						cout<<"Found a blank line"<<endl;
					isdataline = false;
				}
				else if (line[0]=='#')
				{
					jitter=0;
					// line is a comment line
					if (!silent)
						cout<<"headers were found to be: "<<line<<endl;
					isdataline = false;
				}
				// OK, got a line of data so load it up
				else if (line[0]!='#' && line.length()>1)
				{
					// line has data, so convert it
					ss<<line;
					ss>>obsDate>>RV>>RV_error>>jitter;
					if (!silent)
						cout<<obsDate<<" "<<RV<<" "<<RV_error<<" "<<jitter<<endl;
					ss.clear(); //clean out stream for next time
					ss.str(std::string());

					// push these data into the data vectors
					// and increment number of epoch count
					if (obsDate>1)
					{
						epochs_RV_curr.push_back(obsDate);
						RVs_cur.push_back(RV);
						errSquared = RV_error*RV_error+jitter*jitter;

						double inv_var = 1.0/errSquared;
						if (!silent)
						{
							cout<<"jitter found to be = "<<jitter<<endl;
							cout<<"errSquared calculated to be = "<<errSquared<<endl;
							cout<<"inv_var calculated to be = "<<inv_var<<endl;
						}
						RV_inv_var_cur.push_back(inv_var);
						if (errSquared==0 or inv_var==0)
						{
							cout<<"!! errSquared or inv_var is zero!!"<<endl;
							cout<<"Line number "<<lineNumber<<" was: "<<line<<endl;
							cout<<"This was parsed as: obsDate = "<< obsDate<<", RV = "<<RV<<", RV_error = "<<RV_error<<endl;
							cout<<"These were calculated into: errSquared = "<<errSquared<<", and inv_var = "<<inv_var<<endl;
						}
						++numEpochs_cur;
					}
				}
			}// end of while for current data set

			// done loading up this data set
			// so, push it into the double vector params
			numEpochs_RV = numEpochs_RV+numEpochs_cur;
			epochs_RV.push_back(epochs_RV_curr);
			RVs.push_back(RVs_cur);
			RV_inv_var.push_back(RV_inv_var_cur);
			if (verbose)
				cout<<"There were "<< numEpochs_cur<<" epoch(s) found in data set # "<<datasetNum<<endl;
		}//totally done reading in this data set, move to next(or eof)
	}//reached end of file, close it and end func
	infile.close();

	if (verbose)
	{
		cout<<"*Finished loading all the data from the file*"<<endl;
				cout<<"\nResulting inputs were:"<<endl;
		// print final resulting vector sizes
		cout<<"\nepochs_RV size: "<<epochs_RV.size()<<endl;
		cout<<"RVs size: "<<RVs.size()<<endl;
		cout<<"RV_inv_var size: "<<RV_inv_var.size()<<endl;
		cout<<"The total number of epochs found was "<<numEpochs_RV<<"\n"<<endl;
	}
}//end dataLoadUp
