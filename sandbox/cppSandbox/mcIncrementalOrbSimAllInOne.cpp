#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <sstream>
#include <time.h>
#include <fstream>
#include "rnd/randomc.h" //a library with advanced uniform random number generators
#include "rnd/sfmt.cpp"
//#include <rnd/stocc.h> //a library with advanced non-uniform random number generators

#define PI 3.14159265359

using namespace std;

struct orbitCalcReturnType
{
	double n;
	double M_deg;
	double E_latest_deg;
	double TA_deg;
	double Sep_Dist_AU_OP;
	double SA_arcsec_RP_model;
	double PA_deg_RP_model;
	double a1;
	double a2;
};

struct orbitCalcInputType
{
	double t;
	double Sys_Dist_PC;
	double inclination_deg;
	double longAN_deg;
	double e;
	double T;
	double period;
	double argPeri_deg;
	double a_total;
	double Mass1;
	double Mass2;
	bool verbose;
};

struct multiEpochOrbCalcInputType
{
	vector<double> SA_arcsec_measured_REALs;
	vector<double> SA_mean_errors;
	vector<double> PA_deg_measured_REALs;
	vector<double> PA_mean_errors;
	vector<double> epochs;
	double Sys_Dist_PC;
	double inclination_deg;
	double longAN_deg;
	double e;
	double T;
	double period;
	double argPeri_deg;
	double a_total;
	double Mass1;
	double Mass2;
	bool verbose;
};

struct multiEpochOrbCalcReturnType
{
	double chi_squared_total;
	vector<double> ns;
	vector<double> Ms;
	vector<double> Es;
	vector<double> thetas;
	vector<double> Sep_Dists;
	vector<double> SA_arcsec_measured_models;
	vector<double> PA_deg_measured_models;
	vector<double> a1s;
	vector<double> a2s;

};

struct simDataType
{
	string filename;
	int numEpochs;
	int numSamplesAccepted;
	//inputs
	vector<double> longAN_degs;
	vector<double> es;
	vector<double> Ts;
	vector<double> periods;
	vector<double> inclination_degs;
	vector<double> argPeri_degs;
	vector<double> a_totals;
	//outputs
	vector<double> chiSquareds ;
	vector<vector<double> > a1s2 ;
	vector<vector<double> > a2s2 ;
	vector<vector<double> > Sep_Dists2 ;
	vector<vector<double> > ns2 ;
	vector<vector<double> > Ms2 ;
	vector<vector<double> > Es2 ;
	vector<vector<double> > thetas2 ;
	vector<vector<double> > SA_deg_measured_models2 ;
	vector<vector<double> > PA_deg_measured_models2 ;
};

double diff(double a, double b)
{
	/*
	 This function just quickly calculates the difference between two numbers taking negatives into account.
	 This does not use the abs() function as it is always returning 0.0 for some reason...??? don't know why.
	 */

	double difference;
	// since the numerator is the difference, we need to take sign of each into account
	if ( ((a>=0.0)&&(b>=0.0)) || ((a<0.0)&&(b<0.0)) )
		// same sign so just subtract, squaring later will clear any resulting negative
		if (a>b)
			difference = a - b;
		else
			difference = b - a;
	else if ( (a>=0.0)&&(b<=0.0) )
		// real is negative, so subtract it
		difference = a - b; //NOTE: tried using abs() like in Python, but it causes result to =0.0 in all cases.
	else if ( (a<0.0)&&(b>0.0) )
		// model is negative, so subtract it
		difference = b - a;
	//cout << "In diff: "<< "a ="<<a << ", b ="<<b << ", difference ="<<difference<<endl;
	return difference;
}

double chiSquaredCalc(double real, double error, double model)
{
	double difference = diff(real, model);//using diff function above.

	// put it all together to make final chi**2 for each for this model
	double chi_squared = (pow(difference,2.0))/pow(error,2.0);

	if ( (difference ==0)||(chi_squared==0) )
	{
		cout <<"!!!!! difference or chi_squared in chiSquaredCalc = 0"<<endl;
	}

	/*
	cout <<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
	cout << "*** test: abs(-0.99) = "<<abs(-0.99)<< ", (double)abs((double)-0.99) = "<<(double)abs((double(-0.99)))<<endl;
	cout << scientific <<"real = "<<real<< ", model = "<<model<< ", error "<<error<<endl;
	cout << "difference = "<< difference<<", chi_square = "<< chi_squared <<endl;
	cout <<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
	*/
	return chi_squared;
}

orbitCalcReturnType orbitCalculator(orbitCalcInputType OCIT)
{

	orbitCalcReturnType OCRT;

	//Calculate the mean motion and save to output structure
	OCRT.n = (2.0*PI)/OCIT.period;

	//Calculate the Mean Anomaly and add to output structure
	double M = OCRT.n*((OCIT.t-OCIT.T)/365.25);
	//convert to degrees
	OCRT.M_deg = M*(180.0/PI);

	//Set temp input values for E_latest and E_last for initial values of Newton's loops below
	double E_last = 2.0*PI;
	double E_best_intial_guess = M+(OCIT.e*sin(M))+((OCIT.e*OCIT.e)/(2.0*M))*sin(2.0*M);
	double E_latest = E_best_intial_guess;
	int count = 0;

	//Perform Newton's loop to find the solution to Kepler's equation
	while ( (diff(E_last,E_latest)>1.0e-10)&&(count<100) )
	{
		E_last = E_latest;
		E_latest = E_last-((E_last-M-OCIT.e*sin(E_last))/(1.0-OCIT.e*cos(E_last)));
		count +=1;
	}

	//Double check solution satisfies original equation
	//Then save the found Eccentric Anomaly if satisfied
	double Mnewton = (180.0/PI)*(E_latest-OCIT.e*sin(E_latest));
	if ( diff(OCRT.M_deg,Mnewton)>1.0e-5 )
	{
		cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
		cout << "PROBLEM: This resultant E does not satisfy the original equation, Newton's method FAILED !!!"<<endl;
		cout << "E_latest found "<<E_latest*(180.0/PI) <<" causes Mnewton "<<Mnewton <<" != M_deg " <<OCRT.M_deg << " using best initial guess for E = "<< E_best_intial_guess*(180.0/PI)<<endl;
		cout << fixed << "eccentricity = "<<OCIT.e<<endl; //#$$$$$$$
		cout << "period = "<< OCIT.period<<endl; //##$$$$$$$
		cout << "epoch = "<< OCIT.t <<endl; //$$$$$$$$$
		cout << "time of last periapsis = "<<OCIT.T <<endl;//$$$$$$$$$$$$$
		cout << "n = "<<OCRT.n<< endl;
		if ( count>98 )
			cout << " The Newtons loop went up to 100 count!!"<<endl;
		cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
	}
	else
	{
		if ( OCIT.verbose==true )
			cout<<"This resultant E solves the original equation, Newton's Method worked :-)"<<endl;
		//Save solution to output structure
		OCRT.E_latest_deg = E_latest*(180.0/PI);
	}

	//Calculate the True Anomaly from the Eccentric Anomaly
	double TA_rad = acos((cos(E_latest)-OCIT.e)/(1.0-OCIT.e*cos(E_latest)));
	if ( E_latest>PI )
		TA_rad = 2.0*PI - TA_rad;
	//Save True Anomaly to output structure
	OCRT.TA_deg = TA_rad*(180.0/PI);

	if ( OCIT.verbose==true )
		cout<< "True Anomaly [deg] = "<<OCRT.TA_deg<<endl;

	//Calculate the Separation Distance in the orbital plane and save to output structure
	OCRT.Sep_Dist_AU_OP = OCIT.a_total*((1.0-OCIT.e*OCIT.e)/(1.0+OCIT.e*cos(TA_rad)));

	if ( OCIT.verbose==true )
		cout<<  "Separation Distance in orbital plane [AU] = "<< OCRT.Sep_Dist_AU_OP<<endl;

	////Calculate the Separation Distance in the reference plane
	double inclination_rad = (PI/180.0)*OCIT.inclination_deg;
	double ang_LN_to_M2_deg = (OCIT.argPeri_deg+OCRT.TA_deg);
	double ang_LN_to_M2_rad = (PI/180.0)*ang_LN_to_M2_deg;
	if ( OCIT.verbose==true )
	{
		cout<< "inclination_deg = "<< OCIT.inclination_deg<<endl;
		cout<< "inclination_rad = "<< inclination_rad <<endl;
		cout<< "ang_LN_to_M2_rad = "<<ang_LN_to_M2_rad <<endl;
	}
	//Break it up into it's X and Y components WRT Line of Nodes
	double Sep_Dist_AU_RP_y = OCRT.Sep_Dist_AU_OP*sin(ang_LN_to_M2_rad)*cos(inclination_rad);
	double Sep_Dist_AU_RP_x = OCRT.Sep_Dist_AU_OP*cos(ang_LN_to_M2_rad);
	//Calculate the hypotenuse of this, which is the Separation Distance in the reference plane
	double Sep_Dist_AU_RP = sqrt(pow(Sep_Dist_AU_RP_x,2.0)+pow(Sep_Dist_AU_RP_y,2.0));

	if ( OCIT.verbose==true )
	{
		cout<<"Sep_Dist_AU_RP_y = " <<Sep_Dist_AU_RP_y <<endl;
		cout<< "Sep_Dist_AU_RP_x = "<<Sep_Dist_AU_RP_x <<endl;
		cout << "Separation Distance in the reference plane [AU] = " <<Sep_Dist_AU_RP << endl;
	}

	//Calculate the angle between the Ascending Node and M2 in reference plane
	double ang_LN_to_M2_corrPOSNEG = (180.0/PI)*atan(Sep_Dist_AU_RP_y/Sep_Dist_AU_RP_x);
	double ang_LN_to_M2_corr;
	//the atan will produce an angle between the X-axis in each quadrant, rather than only CC from positive X-axis.
	//must take all 4 quadrants into account.
	if ( (Sep_Dist_AU_RP_x>=0.0)&&(Sep_Dist_AU_RP_y>=0.0) )
	{
		// quadrant 1
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG;
	}
    else if ( (Sep_Dist_AU_RP_x<0.0)&&(Sep_Dist_AU_RP_y>=0.0) )
    {
        // quadrant 2
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +180.0;
    }
     else if ( (Sep_Dist_AU_RP_x<0.0)&&(Sep_Dist_AU_RP_y<0.0) )
     {
        // quadrant 3
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +180.0;
     }
     else if ( (Sep_Dist_AU_RP_x>=0.0)&&(Sep_Dist_AU_RP_y<0.0) )
     {
        // quadrant 4
        ang_LN_to_M2_corr = ang_LN_to_M2_corrPOSNEG +360.0;
     }
     else
    	 cout << " the value of " << ang_LN_to_M2_corrPOSNEG << " is apparently out of range"<<endl;

	if ( OCIT.verbose==true )
	{
		cout<< "uncorrected angle between line of nodes and M2 [deg] = "<<  ang_LN_to_M2_deg <<endl;
		cout<< "Corrected angle between Line of nodes and M2 [deg] = "<< ang_LN_to_M2_corr <<endl;
	}

	////Calculate the Position angle in reference plane
	// it is measured from M1, only if the value is over 360deg do we need to fix it
	double totalAngle_RP = OCIT.longAN_deg + ang_LN_to_M2_corr;
	if ( OCIT.verbose==true )
		cout<< "Total angle in the Reference Plane [deg] = "<< totalAngle_RP<<endl;

	if ( totalAngle_RP>360.0)
		OCRT.PA_deg_RP_model = totalAngle_RP - 360.0;
	else
		OCRT.PA_deg_RP_model = totalAngle_RP;

	if ( OCIT.verbose==true )
	{
		cout<< "Position Angle in reference plane [deg] = "<< OCRT.PA_deg_RP_model<<endl;
	}

	OCRT.SA_arcsec_RP_model = Sep_Dist_AU_RP/OCIT.Sys_Dist_PC;

	if ( OCIT.verbose==true )
	{
		cout<<"Separation Angle measured (model) [arcsec] = " << OCRT.SA_arcsec_RP_model<<endl;
	}


	if ( (OCIT.Mass1!=1)&&(OCIT.Mass2!=1))
	{
		// this means these two parameters are non-default values and
		// the system must then be a binary star system, thus calculate the
		// a1 and a2 values of the system.
		// NOTE: Mass2 MUST < Mass1, ie. primary is bigger.

		//first find the mass ratio of system
		double X = OCIT.Mass2/OCIT.Mass1;

		//find a1 from the a1=X*a2 and a_total=a1+a2
		OCRT.a2 = OCIT.a_total/(1.0+X);

		//now use a1=X*a2 to get a2
		OCRT.a1 = OCRT.a2*X;
	}

	else
	{
		OCRT.a1 = 0.0;
		OCRT.a2 = OCIT.a_total;
	}

	return OCRT;
}

multiEpochOrbCalcReturnType multiEpochOrbCalc(multiEpochOrbCalcInputType MEOCIT)
{
	multiEpochOrbCalcReturnType  MEORCT;

	//int i = 0;
	double chi_squared_total = 0.0;

	for ( int i=0; i<((int) MEOCIT.epochs.size()); i++ )
	{
		orbitCalcInputType OCIT;


		// load OCIT structure with necessary values from the MEOCIT structure
		OCIT.t = MEOCIT.epochs[i];
		OCIT.Mass1 = MEOCIT.Mass1;
		OCIT.Mass2 = MEOCIT.Mass2;
		OCIT.Sys_Dist_PC = MEOCIT.Sys_Dist_PC;
		OCIT.T = MEOCIT.T;
		OCIT.a_total = MEOCIT.a_total;
		OCIT.longAN_deg = MEOCIT.longAN_deg;
		OCIT.argPeri_deg = MEOCIT.argPeri_deg;
		OCIT.e = MEOCIT.e;
		OCIT.inclination_deg = MEOCIT.inclination_deg;
		OCIT.period = MEOCIT.period;
		OCIT.verbose = MEOCIT.verbose;

		if ( MEOCIT.verbose==true )
			cout << "----------------------------------------------------------------------------" <<endl; //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

		// instantiate the OCRT structure and pass the OCIT structure into
		// orbitCalculator to run the model and pass back the complete OCRT.
		orbitCalcReturnType OCRT;
		OCRT = orbitCalculator(OCIT);

		MEORCT.ns.push_back(OCRT.n);
		MEORCT.Ms.push_back(OCRT.M_deg);
		MEORCT.Es.push_back(OCRT.E_latest_deg);
		MEORCT.thetas.push_back(OCRT.TA_deg);
		MEORCT.Sep_Dists.push_back(OCRT.Sep_Dist_AU_OP);
		MEORCT.SA_arcsec_measured_models.push_back(OCRT.SA_arcsec_RP_model);
		MEORCT.PA_deg_measured_models.push_back(OCRT.PA_deg_RP_model);
		MEORCT.a1s.push_back(OCRT.a1);
		MEORCT.a2s.push_back(OCRT.a2);

		// grab necessary values for the chi squared calculation from the MEOCIT structure
		double SA_arcsec_measured_REAL = MEOCIT.SA_arcsec_measured_REALs[i];
		double SA_mean_error = MEOCIT.SA_mean_errors[i];
		double PA_deg_measured_REAL = MEOCIT.PA_deg_measured_REALs[i];
		double PA_mean_error = MEOCIT.PA_mean_errors[i];

		double SA_chi_squared = chiSquaredCalc(SA_arcsec_measured_REAL, SA_mean_error, OCRT.SA_arcsec_RP_model);
		// we must account for the 360deg boundry, using +- 25deg from it as the region of conflict
		double PA_deg_measured_model = OCRT.PA_deg_RP_model;
        if ( ((OCRT.PA_deg_RP_model-25.0)<0.0) && (PA_deg_measured_REAL>335.0) )
            //ie real is close to 360 and model is just over 360
            PA_deg_measured_model = PA_deg_measured_model +360.0;
        else if ( ((OCRT.PA_deg_RP_model-25.0)<0.0) && (PA_deg_measured_model>335.0) )
            //ie model is close to 360 and real is just over 360
            PA_deg_measured_REAL = PA_deg_measured_REAL +360.0;
        double PA_chi_squared = chiSquaredCalc(PA_deg_measured_REAL, PA_mean_error, PA_deg_measured_model);
        // Add them to get the updated total
        chi_squared_total = chi_squared_total + PA_chi_squared + SA_chi_squared;
        if ( MEOCIT.verbose==true )
        {
        	cout << "SA_arcsec_measured_REAL = "<< SA_arcsec_measured_REAL<<", SA_mean_error = " <<SA_mean_error << ", OCRT.SA_arcsec_RP_model = " << OCRT.SA_arcsec_RP_model<<endl;
        	cout << " SA chi_squared = "<< SA_chi_squared<<  ", PA chi_square = "<< PA_chi_squared<< endl;
        	cout << fixed <<" in loop chi_squared_total = "<<chi_squared_total <<endl;
        }
		// increment to next epoch to check these input for
	}

	//cout <<"349: output  chi_squared_total = "<<chi_squared_total <<endl;
	MEORCT.chi_squared_total = chi_squared_total;

	return MEORCT;
}

double gaussianDist(double x, double mean, double variance)
{
	//Just a simple function to return the probability value of a value in a gaussian distribution.

	double varSquared = pow(variance,2.0);

	double exponent = -1*pow((x-mean),2.0)/(2.0*varSquared);
	//double front = 1/sqrt(2.0*PI*varSquared) //multiply the prob by this to give the normal distribution value

	double prob = exp(exponent);

	return prob;
}

int fileWriter(simDataType SDT)
{
	//find the number of epochs and samples saved in the SDT
		//SDT.numEpochs;
		//SDT.numSamplesAccepted;

		//find if the input SDT data is of the lowRAM type or not
		bool lowRAM = false;
		if ( (SDT.Es2[0][0] == 0.0) && ( SDT.Ms2[0][0] == 0.0) )
		{
			cout << "It was found that this is lowRAM data being written to file"<<endl;
			lowRAM = true;
		}
		// make an output file for the INS and load it up.
		ofstream INSfile;
		string INSfilename;
		string rootname(SDT.filename);
		string nameEnding("_INS.txt");
		INSfilename = rootname + nameEnding;
		cout<<"Writing file to: "<<INSfilename.c_str()<<endl;
		INSfile.open(INSfilename.c_str()) ;     //open the output file
		INSfile << INSfilename.c_str() <<endl;
		INSfile << "longAN [deg]  es [N/A]  Ts [julian date]  periods [yrs]  inclination [deg]  argPeri [deg] chiSquareds"<<endl;
		for (int sample=0; sample<SDT.numSamplesAccepted; sample++)
		{
			INSfile<<fixed<< SDT.longAN_degs[sample];
			INSfile<< "        "<<SDT.es[sample];
			INSfile<< "     "<< SDT.Ts[sample];
			INSfile<< "     "<< SDT.periods[sample];
			INSfile<< "         "<< SDT.inclination_degs[sample];
			INSfile<< "         "<< SDT.argPeri_degs[sample];
			INSfile<< "         "<< SDT.chiSquareds[sample]<<endl;
		}//finished writing inputs file
		INSfile.close () ;

		// make an output file for each of the OUTS epochs
		for (int epoch=0; epoch<SDT.numEpochs; epoch++)
		{
			ofstream OUTSfile;
			stringstream OUTSfilename;
			OUTSfilename << SDT.filename<<"_OUTSepoch"<<(epoch+1)<<".txt";
			cout<<"Writing output file to: "<<OUTSfilename.str()<<endl;
			string outputname(OUTSfilename.str());
			OUTSfile.open(outputname.c_str()) ;     //open the output file
			OUTSfile << OUTSfilename.str() << endl;
			OUTSfile<< "Sep_Dists [AU]  thetas [deg]   Es [deg]   Ms [deg]   ns [rad/yr]  a1s [AU]   a2s [AU]"<<endl;
			for (int sample=0; sample<SDT.numSamplesAccepted; sample++)
			{
				//Take care of the lowRAM option for parameter values
				if ( lowRAM==false )
				{
					OUTSfile<< SDT.Sep_Dists2[sample][epoch];
					OUTSfile<< "          "<< SDT.thetas2[sample][epoch] ;
					OUTSfile<< "       "<< SDT.Es2[sample][epoch];
					OUTSfile<< "     "<< SDT.Ms2[sample][epoch];
					OUTSfile<< "     "<< SDT.ns2[sample][epoch];
					OUTSfile<< "     "<< SDT.a1s2[sample][epoch];
					OUTSfile<< "        "<< SDT.a2s2[sample][epoch]<<endl;
				}
				else
				{
					//use zeros to fill up the columns with not saved parameters due to lowRAM feature
					OUTSfile<< 0;//for Sep_dists
					OUTSfile<< "                   "<< 0;//for thetas
					OUTSfile<< "            "<< 0;//for Es
					OUTSfile<< "          "<< 0;//for Ms
					OUTSfile<< "          "<< 0;//for ns
					OUTSfile<< "             "<< SDT.a1s2[sample][epoch];
					OUTSfile<< "         "<< SDT.a2s2[sample][epoch]<<endl;
				}
			}//finished writing current epoch output file
			OUTSfile.close();
		}//finished writing all epoch output files

		return EXIT_SUCCESS;
}

double uniformRandNumber(double min, double max)
{
	//create very high precision seed for generator (nanosecond precision)
	timespec time1;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);//even though this is apparently a 'bug', it works so leaving it.
	int time_nsec=time1.tv_nsec;
	//cout<<time_nsec<<endl;

	//choose generator
	CRandomSFMT1 RanGen(time_nsec);
	//this is the most advanced uniform random number generator, which combines SFMT and Mother-Of-All.

	// now generate random number of double precision between min and max
	double randNumber;
	double rand_zero_to_one = RanGen.Random();
	randNumber = (rand_zero_to_one*(max-min))+min;

	if ( (randNumber>max)||(randNumber<min) )
	{
		cout<< "*** PROBLEM with random number generator!!"<<endl;
		cout << "random number generated "<<randNumber<< ", is outside the range ["<<min<<", "<<max<<"]"<<endl;
		cout<< fixed << "min = "<<min<< ", max = "<<max<<", rand [0,1] = "<<rand_zero_to_one<<", return random number = "<<randNumber<<endl;//$$$$$$$$$$$
	}
	return randNumber;
}

int main()
{
	int numSamples = 1e6;
	bool silent = false;
	bool silent2 = true;
	bool lowRAM = true;//this boolean will cause low priority variables to not be saved, thus saving memory.
	string filename = "/run/media/Kyle/Data1/Todai_Work/Data/data_Binary/data_C++/BLAH";

	//*********** In this brick is all you need to play with *************
    //********************************************************************

    double numSamplePrints = 10.0   ;
    int tuningSamples = (int)(numSamples/5.0);
    double sigma_latest = 0.5;//ie. the whole range
    double chiSquaredMax_latest = 1000.0 ;// $$$$$$$$$ Maybe temp!!! $$$$$$$$$$$$$$$

    // Ranges for acceptable random number inputs ######
    double longAN_degMIN = 0.01; // [deg]
    double longAN_degMAX = 179.0 ;// [deg]
    double eMIN = 0.1 ;
    double eMAX = 0.5 ;
    double periodMIN = 1.0; // [yrs]
    double periodMAX = 20.0; // [yrs]
    double inclination_degMIN = 0.01; // [deg]
    double inclination_degMAX = 180.0; // [deg]
    double argPeri_degMIN = 0.01; // [deg]
    double argPeri_degMAX = 179.0; // [deg]
    double a_totalMIN = 1.0;
    double a_totalMAX = 5.0;
    // MEASURED VALUES ########################

    // TEST system data (HD 130948BC Dupuy2009)
	multiEpochOrbCalcInputType MEOCIT;
	double SA_arcsec_measured_REALs[] = {0.0993, 0.0946, 0.0961, 0.0978, 0.0568, 0.0564, 0.0573, 0.1117, 0.1091, 0.0979, 0.0722, 0.0574, 0.0573, 0.0584, 0.0517};
	double SA_mean_errors[]= {0.0017, 0.0017, 0.0023, 0.0023, 0.0009, 0.0012, 0.0009, 0.0015, 0.0011, 0.0009, 0.0006, 0.0005, 0.0009, 0.0010, 0.0006};
	double PA_deg_measured_REALs[]= {307.9,306.9,308.9,307.0,144.9,146.4,148.6,132.6,132.28,130.73,127.6,124.7,124.6,124.1,123.9};
	double PA_mean_errors[]= {1.1,1.0,1.6,1.4,0.5,0.6,0.6,0.4,0.13,0.17,0.3,0.4,0.6,0.7,0.5};
	double epochs[]= {2452519.702,2452519.706,2452519.709,2452519.719,2453425.152,2453425.219,2453425.282,2454127.085,2454185.028,2454306.737,2454481.192,2454555.068,2454555.09,2454555.157,2454584.106};
	MEOCIT.Sys_Dist_PC = 18.17;
	MEOCIT.Mass1 = 1;
	MEOCIT.Mass2 = 1;
	// ************************************************************************
	// ***********************************************************************



	// ******* STARTING SIMULATION!!!!  *************

	time_t startTime;
	startTime = time(NULL);

	//load MEOCIT struct vectors with above list values
	int epoch = 0;
	int numEpochs=sizeof SA_arcsec_measured_REALs/sizeof(double);

	if ( silent==false )
		MEOCIT.verbose = false;//####$$$$$$$$$$$$$$$$$$$$$ TEMP
	else
		MEOCIT.verbose = false;

	while ( epoch<numEpochs )
	{
		MEOCIT.SA_arcsec_measured_REALs.push_back(SA_arcsec_measured_REALs[epoch]);
		MEOCIT.SA_mean_errors.push_back(SA_mean_errors[epoch]);
		MEOCIT.PA_deg_measured_REALs.push_back(PA_deg_measured_REALs[epoch]);
		MEOCIT.PA_mean_errors.push_back(PA_mean_errors[epoch]);
		MEOCIT.epochs.push_back(epochs[epoch]);
		epoch +=1;
	}
	//cout << "********************************************************"<<endl;
	if ( silent==false )
	{
		cout<< "\nMCONLY: $$$$$$$$$$$$$$$$$$$  STARTING THE SIMULATOR  $$$$$$$$$$$$$$$$$" <<endl;
		int numMillions = numSamples/int(1e6);
		stringstream numSamplesString;
		numSamplesString << numMillions << "-million";
		cout << "Number of sample orbits being created = " << numSamplesString.str() <<endl;
	}

	// instantiate Simulator Data Type
	simDataType SDT;

	//make initial values for the input parameters to model for the sigma based random number generation
	double longAN_degs_latest_min = (longAN_degMAX-longAN_degMIN)/2.0;
	double es_latest_min = (eMAX-eMIN)/2.0;
	double TMAX = epochs[0]-periodMAX*365.0;
	double TMIN = epochs[0];
	double Ts_latest_min = (TMAX-TMIN)/2.0;
	double periods_latest_min = (periodMAX-periodMIN)/2.0;
	double inclination_degs_latest_min = (inclination_degMAX-inclination_degMIN)/2.0;
	double argPeri_degs_latest_min = (argPeri_degMAX-argPeri_degMIN)/2.0;
	double a_totals_latest_min = (a_totalMAX-a_totalMIN)/2.0;

	// load filename defined above into the SDT structure
	SDT.filename = filename;
	SDT.numEpochs = numEpochs;
	// variables for the success rate print block in chain loop
    int printTime = numSamples/numSamplePrints;
    int printCount = 0;
    int printsDone = 0;
    int acceptedCounter = 0;
    double reduced_chiSquared_total_cur;
    double original_chiSquared_total_cur;
    double latest_accepted_chiSquared = 0;
    double one_over_nu;
    double chiSquareMin = chiSquaredMax_latest;
    int bestOrbit = 0;
    int totalCounter = 0;
    double sigma, max, min;


    //Store values of zero in each of the low priority variables to save RAM
    // rather than loading them up.  Later these zeros will trigger the
    // subsequent functions to act accordingly.
    if ( lowRAM==true )
    {
    	vector<double> tempVec;
		tempVec.push_back(0.0);
		SDT.ns2.push_back(tempVec);
		SDT.Ms2.push_back(tempVec);
		SDT.Es2.push_back(tempVec);
		SDT.thetas2.push_back(tempVec);
		SDT.Sep_Dists2.push_back(tempVec);
		SDT.SA_deg_measured_models2.push_back(tempVec);
		SDT.PA_deg_measured_models2.push_back(tempVec);
    }

    if ( silent==false )
    	cout<< "     Starting Sample Loop!   "<<endl;

    // ***** Start the samples loop *****
    for ( int sample=1; sample<(numSamples+1); sample++)
    {

		// block to control printing success rate to screen
        printCount = printCount + 1;
        if ( printCount==printTime )
        {
            printsDone = printsDone+1;
            printCount = 0;
            if ( silent==false )
            {
            	// create a time object and find current time
				time_t rawtime;
				struct tm * timeinfo;
				time ( &rawtime );
				timeinfo = localtime ( &rawtime );
				// find number of samples accepted so far
				int len=SDT.es.size();
				// print number of samples accepted so far and current time
				cout <<fixed<< "latest reduced chi-square = " << reduced_chiSquared_total_cur <<endl;
				cout << "latest accepted ch-square = "<<latest_accepted_chiSquared<<endl;
				cout << "lowest chiSquare so far = "<< chiSquareMin <<endl;
				cout << "latest sigma = "<< sigma_latest<<endl;
                cout << len<<" successful out of "<<sample<<". "<<printsDone<<"/"<<(int)numSamplePrints<<" completed at ";
                cout << asctime (timeinfo) <<endl;
            }
        }

        //// Generate random numbers in the required ranges for the inputs to the orbCalc
        //// by perturbing the current ones

        // Inclination
        // sigma is a percentage of total range, so convert its current percentage to a
        // range value
        sigma = sigma_latest*(inclination_degMAX-inclination_degMIN);
        if ( (inclination_degs_latest_min+sigma)>=inclination_degMAX )
        {
        	max = inclination_degMAX;
        	min = inclination_degMAX-(sigma*2.0);
        }
        else if ( (inclination_degs_latest_min-sigma)<=inclination_degMIN )
        {
        	max = inclination_degMIN+(sigma*2.0);
			min = inclination_degMIN;
        }
        else
        {
        	max = inclination_degs_latest_min + sigma;
			min = inclination_degs_latest_min - sigma;
        }
        MEOCIT.inclination_deg = uniformRandNumber(min, max);


        // Longitude of Ascending Node
		// sigma is a percentage of total range, so convert its current percentage to a
		// range value
		sigma = sigma_latest*(longAN_degMAX-longAN_degMIN);
		if ( (longAN_degs_latest_min+sigma)>=longAN_degMAX )
		{
			max = longAN_degMAX;
			min = longAN_degMAX-(sigma*2.0);
		}
		else if ( (longAN_degs_latest_min-sigma)<=longAN_degMIN )
		{
			max = longAN_degMIN+(sigma*2.0);
			min = longAN_degMIN;
		}
		else
		{
			max = longAN_degs_latest_min + sigma;
			min = longAN_degs_latest_min - sigma;
		}
		MEOCIT.longAN_deg = uniformRandNumber(min, max);


		// Argument of Periapsis
		// sigma is a percentage of total range, so convert its current percentage to a
		// range value
		sigma = sigma_latest*(argPeri_degMAX-argPeri_degMIN);
		if ( (argPeri_degs_latest_min+sigma)>=argPeri_degMAX )
		{
			max = argPeri_degMAX;
			min = argPeri_degMAX-(sigma*2.0);
		}
		else if ( (argPeri_degs_latest_min-sigma)<=argPeri_degMIN )
		{
			max = argPeri_degMIN+(sigma*2.0);
			min = argPeri_degMIN;
		}
		else
		{
			max = argPeri_degs_latest_min + sigma;
			min = argPeri_degs_latest_min - sigma;
		}
		MEOCIT.argPeri_deg = uniformRandNumber(min, max);

		// Eccentricity
		// sigma is a percentage of total range, so convert its current percentage to a
		// range value
		sigma = sigma_latest*(eMAX-eMIN);
		if ( (es_latest_min+sigma)>=eMAX )
		{
			max = eMAX;
			min = eMAX-(sigma*2.0);
		}
		else if ( (es_latest_min-sigma)<=eMIN )
		{
			max = eMIN+(sigma*2.0);
			min = eMIN;
		}
		else
		{
			max = es_latest_min + sigma;
			min = es_latest_min - sigma;
		}
		MEOCIT.e = uniformRandNumber(min, max);


		// Argument of Periapsis
		// sigma is a percentage of total range, so convert its current percentage to a
		// range value
		sigma = sigma_latest*(periodMAX-periodMIN);
		if ( (periods_latest_min+sigma)>=periodMAX )
		{
			max = periodMAX;
			min = periodMAX-(sigma*2.0);
		}
		else if ( (periods_latest_min-sigma)<=periodMIN )
		{
			max = periodMIN+(sigma*2.0);
			min = periodMIN;
		}
		else
		{
			max = periods_latest_min + sigma;
			min = periods_latest_min - sigma;
		}
		MEOCIT.period = uniformRandNumber(min, max);

		// Time of Last Periapsis
		// sigma is a percentage of total range, so convert its current percentage to a
		// range value
		// thus between a full period before first observation and the time of first observation
		double Tmin = epochs[0]-MEOCIT.period*365.0;
		double Tmax = epochs[0];
		sigma = sigma_latest*(Tmax-Tmin);
		if ( (Ts_latest_min+sigma)>=Tmax )
		{
			max = Tmax;
			min = Tmax-(sigma*2.0);
		}
		else if ( (Ts_latest_min-sigma)<=Tmin )
		{
			max = Tmin+(sigma*2.0);
			min = Tmin;
		}
		else
		{
			max = Ts_latest_min + sigma;
			min = Ts_latest_min - sigma;
		}
		MEOCIT.T = uniformRandNumber(min, max);

		// Total Semi-Major Axis
		// sigma is a percentage of total range, so convert its current percentage to a
		// range value
		sigma = sigma_latest*(a_totalMAX-a_totalMIN);
		if ( (a_totals_latest_min+sigma)>=a_totalMAX )
		{
			max = a_totalMAX;
			min = a_totalMAX-sigma;
		}
		else if ( (a_totals_latest_min-sigma)<=a_totalMIN )
		{
			max = a_totalMIN+sigma;
			min = a_totalMIN;
		}
		else
		{
			max = a_totals_latest_min + sigma;
			min = a_totals_latest_min - sigma;
		}
		MEOCIT.a_total = uniformRandNumber(min, max);

        if ( silent2==false )
        {
        	cout<< "inclination_deg = " <<MEOCIT.inclination_deg  <<endl;
        	cout<< "longAN_deg = " << MEOCIT.longAN_deg  <<endl;
        	cout<<  "argPeri_deg = "<< MEOCIT.argPeri_deg <<endl;
        	cout<<  "e = "<< MEOCIT.e <<endl;
        	cout<<  "period = "<<  MEOCIT.period<<endl;
        	cout<<  "T = "<< MEOCIT.T  <<endl;
        	cout<<  "a_total = "<< MEOCIT.a_total <<endl;
        }

        multiEpochOrbCalcReturnType MEOCRT;

        //// Call the orbCalc to have it apply the model to the inputs and produce outputs
		MEOCRT = multiEpochOrbCalc(MEOCIT);

		//// Calculate the Metropolis-Hastings Algorithm

		// Calculate the reduced chiSquared from the returned chiSquared
		original_chiSquared_total_cur = MEOCRT.chi_squared_total;
		one_over_nu = (1.0/((2.0*(double)numEpochs)-7.0));
		reduced_chiSquared_total_cur = one_over_nu*original_chiSquared_total_cur;
		if ( silent2==false )
			cout<< "output reduced chi squared = "<< reduced_chiSquared_total_cur <<endl;

		//// Determine if the orbit should be accepted
		if ( reduced_chiSquared_total_cur <= chiSquaredMax_latest )
		{
			acceptedCounter += 1;
			latest_accepted_chiSquared = reduced_chiSquared_total_cur;

			if ( reduced_chiSquared_total_cur<chiSquareMin)
			{
				chiSquareMin = reduced_chiSquared_total_cur;
				bestOrbit = acceptedCounter-1;
			}

			// store inputs
			SDT.longAN_degs.push_back(MEOCIT.longAN_deg);
			SDT.es.push_back(MEOCIT.e);
			SDT.Ts.push_back(MEOCIT.T);
			SDT.periods.push_back(MEOCIT.period);
			SDT.inclination_degs.push_back(MEOCIT.inclination_deg);
			SDT.argPeri_degs.push_back(MEOCIT.argPeri_deg);

			// store outputs
			SDT.a1s2.push_back(MEOCRT.a1s);
			SDT.a2s2.push_back(MEOCRT.a2s);
			SDT.chiSquareds.push_back(reduced_chiSquared_total_cur);
			// Only save these below if you have plenty of RAM, as they are not really needed.
			// above they are loaded with zeros to trigger later functions to act accordingly
			// if lowRAM=true, so we don't need to do anything here in that case.
			if ( lowRAM==false )
			{
				SDT.ns2.push_back(MEOCRT.ns);
				SDT.Ms2.push_back(MEOCRT.Ms);
				SDT.Es2.push_back(MEOCRT.Es);
				SDT.thetas2.push_back(MEOCRT.thetas);
				SDT.Sep_Dists2.push_back(MEOCRT.Sep_Dists);
				SDT.SA_deg_measured_models2.push_back(MEOCRT.SA_arcsec_measured_models);
				SDT.PA_deg_measured_models2.push_back(MEOCRT.PA_deg_measured_models);
			}
		}// Done storing accepted orbit parameters

		// Increment tuning counter
		totalCounter +=1;

		//// Sigma tuning block
		if ( totalCounter==tuningSamples)
		{
			//cout<< "Entered sigma tunning block!!"<<endl; //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			totalCounter = 0;
			// update the sigma and chiSquredMax by cutting them both in half
			sigma_latest = sigma_latest/2.0;
			chiSquaredMax_latest = chiSquareMin*2.0;
			if ( silent == false )
			{
				cout<< " Sigma updated to be = "<< sigma_latest<<endl;
				cout<<" ChiSquaredMax updated to be = "<<chiSquaredMax_latest<<endl;
			}
			//update min value for each parameter to be centered on for next round
			longAN_degs_latest_min = SDT.longAN_degs[bestOrbit];
			es_latest_min = SDT.es[bestOrbit];
			Ts_latest_min = SDT.Ts[bestOrbit];
			periods_latest_min = SDT.periods[bestOrbit];
			inclination_degs_latest_min = SDT.inclination_degs[bestOrbit];
			argPeri_degs_latest_min = SDT.argPeri_degs[bestOrbit];
		}

    }//Done sample loops

    int totalAccepted = SDT.es.size();
    SDT.numSamplesAccepted = totalAccepted;
    //cout << "********************************************************"<<endl;
    cout<< totalAccepted <<" orbits were accepted during simulation"<<endl;

    if ( silent==false )
    {
    	cout << "$$$$$$$$$$$$$$$ SIMULATOR COMPLETE $$$$$$$$$$$$$$$"<<endl;
    	cout<< "\nBest orbit found:"<<endl;
    	cout<< "chiSquaredMin = "<< chiSquareMin <<endl;
    	cout<< "LongAN = "<< SDT.longAN_degs[bestOrbit] <<endl;
    	cout<< "e = "<< SDT.es[bestOrbit] <<endl;
    	cout<< "To = "<< SDT.Ts[bestOrbit] <<endl;
    	cout<< "period = "<< SDT.periods[bestOrbit] <<endl;
    	cout<< "inclination = "<< SDT.inclination_degs[bestOrbit] <<endl;
    	cout<< "argPeri = "<< SDT.argPeri_degs[bestOrbit] <<endl;
    }

    time_t endTime;
    endTime = time(NULL);
    int timeElapsed = endTime-startTime;

    if ( timeElapsed < 60.0 )
    	cout<< "Simulator took "<<timeElapsed << " seconds to complete"<<endl;
    else if ( timeElapsed>3600 )
    {
    	double totalTime = ((double)timeElapsed)/3600.0;
    	int hrs = (int)totalTime;
    	int min = 60*(totalTime-(double)hrs);
    	cout<< "Simulator took "<<hrs << " hours and "<< min << " minutes to complete" <<endl;
    }
    else
    {
    	double totalTime = ((double)timeElapsed)/60.0;
    	int min = (int)totalTime;
    	int sec = 60*(totalTime-double(min));
    	cout <<"Simulator took "<<min <<" minutes and "<<sec<< " seconds to complete"<<endl;
    }



    fileWriter(SDT);

	return EXIT_SUCCESS;
}
