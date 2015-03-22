#include <iostream>
#include <stdio.h>
#include <stdlib.h>
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
	double chiSquaredMax;
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
	vector<double> SA_arcsec_measured_modelsl;
	vector<double> PA_deg_measured_models;
	vector<double> a1s;
	vector<double> a2s;

};

struct simDataType
{
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
	double E_best_intial_guess = M+(OCIT.e*sin(M));
	double E_latest = E_best_intial_guess;
	int count = 0;

	//Perform Newton's loop to find the solution to Kepler's equation
	while ( (diff(E_last,E_latest)>1.0e-10)&&(count<100) )
	{
		E_last = E_latest;
		E_latest = E_last-((E_last-M-OCIT.e*sin(E_last))/(-OCIT.e*cos(E_last)+1.0));
		count +=1;
	}

	//Double check solution satisfies original equation
	//Then save the found Eccentric Anomaly if satisfied
	double Mnewton = (180.0/PI)*(E_latest-OCIT.e*sin(E_latest));
	if ( abs(OCRT.M_deg-Mnewton)>1.0e-5 )
	{
		cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
		cout << "PROBLEM: This resultant E does not satisfy the original equation, Newton's method FAILED !!!"<<endl;
		cout << "E_latest found "<<E_latest*(180.0/PI) <<" causes Mnewton "<<Mnewton <<" != M_deg " <<OCRT.M_deg << " using best initial guess for E = "<< E_best_intial_guess*(180.0/PI)<<endl;
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
    	 cout << " the value of " << ang_LN_to_M2_corrPOSNEG << " is aparently out of range"<<endl;

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
		MEORCT.SA_arcsec_measured_modelsl.push_back(OCRT.SA_arcsec_RP_model);
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

int fileWriter(simDataType SDT)
{
	ofstream outFile;
	outFile.open ( "outFile.txt" ) ;     //open the output file
	outFile << "The donations, sorted in ascending order are:  \n";
	outFile.close () ;
	return EXIT_SUCCESS;
}

double uniformRandNumber(double min, double max)
{
	//create very high precision seed for generator (nanosecond precision)
	timespec time1;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
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
	int numSamples = 1;
	bool silent = false;

	//*********** In this brick is all you need to play with *************
    //********************************************************************
	//double chiSquaredMax = 10000.0 ;
    //double numSamplePrints = 10.0   ;

    // MEASURED VALUES ########################

	// TEST system data (HD 130948BC Dupuy2009)
	multiEpochOrbCalcInputType MEOCIT;
	double SA_arcsec_measured_REALs[] = {0.0993, 0.0946, 0.0961, 0.0978, 0.0568, 0.0564, 0.0573, 0.1117, 0.1091, 0.0979, 0.0722, 0.0574, 0.0573, 0.0584, 0.0517};
	double SA_mean_errors[]= {0.0011,0.0011,0.0017,0.0017,0.0006,0.0009,0.0006,0.0008,0.0004,0.0003,0.0002,0.0002,0.0006,0.0006,0.0003};
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
		MEOCIT.verbose = true;
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

	simDataType SDT;

	// variables for the success rate print block in chain loop
    //int printTime = numSamples/numSamplePrints;
    //int printCount = 0;
    //int printsDone = 0;
    //int acceptedCounter = 0;
    double reduced_chiSquared_total_cur;
    double original_chiSquared_total_cur;
    double one_over_nu;

    // ***** Start the samples loop *****
    for ( int sample=1; sample<(numSamples+1); sample++)
    {

        // Generate random numbers in the required ranges for the inputs to the orbCalc
        MEOCIT.inclination_deg = 95.90;
        MEOCIT.longAN_deg = 133.37;
        MEOCIT.argPeri_deg = 73.85;
        MEOCIT.e = 0.145;
        MEOCIT.period = 9.74; //  [yrs]
        double Tmin = epochs[0]-MEOCIT.period*365.0;
        double Tmax = epochs[0];
        cout<< fixed <<"\n Tmin = "<<Tmin<< endl;
        cout<<"Tmax = "<<Tmax<< endl;
        MEOCIT.T = 2451124.1; // thus between a full period before first observation and the time of first observation
        MEOCIT.a_total = 2.18;

        cout<<fixed <<"\n INPUTS: "<<endl;
        cout<< "inclination: " << MEOCIT.inclination_deg<<endl;
        cout<< "LongAN: " << MEOCIT.longAN_deg<<endl;
        cout<< "Eccentricity: " << MEOCIT.e <<endl;
        cout<< "Time of last periapsis: "<< MEOCIT.T<<endl;
        cout<< "Period: "<< MEOCIT.period<<endl;
        cout<< "ArgPeri: "<< MEOCIT.argPeri_deg<<endl;
        cout<< "total semi-major"<< MEOCIT.a_total<<endl;

        cout<< "\n measured values:"<<endl;
        cout<< "SA: "<<MEOCIT.SA_arcsec_measured_REALs[0] <<endl;
        cout<< "SA Error: "<< MEOCIT.SA_mean_errors[0]<<endl;
        cout<< "PA: "<< MEOCIT.PA_deg_measured_REALs[0]<<endl;
        cout<< "PA Error: "<< MEOCIT.PA_mean_errors[0]<<endl;
        cout<< "epoch: "<< MEOCIT.epochs[0]<<endl;
        cout<< "sys dist: "<< MEOCIT.Sys_Dist_PC<<endl;


        multiEpochOrbCalcReturnType MEOCRT;

        // Call the orbCalc to have it apply the model to the inputs and produce outputs
		MEOCRT = multiEpochOrbCalc(MEOCIT);

		// Calculate the reduced chiSquared from the returned chiSquared
		original_chiSquared_total_cur = MEOCRT.chi_squared_total;
		one_over_nu = (1.0/((2.0*(double)numEpochs)-7.0));
		if ( silent==false )
		{
			cout<<" original_chiSquared_total_cur = "<<original_chiSquared_total_cur <<endl;
			cout<<"one_over_nu = " << one_over_nu<<endl;
		}
		reduced_chiSquared_total_cur = one_over_nu*original_chiSquared_total_cur;

		cout<< "\n OUPUTS: "<<endl;
		cout<< "a1s: " << MEOCRT.a1s[0] <<endl;
		cout<< "a2s: " << MEOCRT.a2s[0] <<endl;
		cout<< "reduced chi squared: " << reduced_chiSquared_total_cur <<endl;
		cout<< "ns: " << MEOCRT.ns[0] <<endl;
		cout<< "Ms: " << MEOCRT.Ms[0] <<endl;
		cout<< "Es: " << MEOCRT.Es[0] <<endl;
		cout<< "thetas: " << MEOCRT.thetas[0] <<endl;
		cout<< "sep dists: " << MEOCRT.Sep_Dists[0] <<endl;
		cout<< "SA model: " << MEOCRT.SA_arcsec_measured_modelsl[0] <<endl;
		cout<< "PA model: " << MEOCRT.PA_deg_measured_models[0] <<endl;

    }//Done sample loops

    int totalAccepted = SDT.es.size();
    //cout << "********************************************************"<<endl;
    cout<< totalAccepted <<" orbits were accepted during simulation"<<endl;

    if ( silent==false )
    	cout << "$$$$$$$$$$$$$$$ SIMULATOR COMPLETE $$$$$$$$$$$$$$$"<<endl;

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
    	cout<< "Simulator took "<<hrs << " hours and"<< min << " minutes to complete" <<endl;
    }
    else
    {
    	double totalTime = ((double)timeElapsed)/60.0;
    	int min = (int)totalTime;
    	int sec = 60*(totalTime-double(min));
    	cout <<"Simulator took "<<min <<" minutes and "<<sec<< " seconds to complete"<<endl;
    }


    // $$$$$$$$$$$$$$$$$$ fill in with fileWritter call $$$$$$$$$$$$$$$$$$$$$$$$$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	return EXIT_SUCCESS;
}
