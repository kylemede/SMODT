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
using namespace std;

#define PI 3.14159265359


struct orbitCalcInputType
{
	double t;
	double Sys_Dist_PC;
	double inclination_deg;
	double longAN_deg;
	double e;
	double M;
	double T;
	double period;
	double argPeri_deg;
	double a_total;
	double Mass1;
	double Mass2;
	bool verbose;
};

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

orbitCalcReturnType orbitCalculator(orbitCalcInputType OCIT)
{

	orbitCalcReturnType OCRT;
	cout <<"\n Starting to run the Newtons loop \n"<<endl;
	//convert to degrees
	double M = OCIT.M;
	OCRT.M_deg = M*(180.0/PI);

	//Set temp input values for E_latest and E_last for initial values of Newton's loops below
	//double E_last = 2.0*PI;
	double E_best_intial_guess1 = M+(OCIT.e*sin(M))+((OCIT.e*OCIT.e)/(2.0*M))*sin(2.0*M);
	double E_best_intial_guess2 = M+(OCIT.e*sin(M));
	double E_latest = E_best_intial_guess2;
	double E_last = E_latest+10.0;
	int count = 0;

	if (OCIT.verbose==true)
	{
		cout << " using best initial guess for E1 = "<< E_best_intial_guess1*(180.0/PI)<<endl;
		cout << " using best initial guess for E2 = "<< E_best_intial_guess2*(180.0/PI)<<endl;
	}

	cout << "about to enter loop"<<endl;
	//Perform Newton's loop to find the solution to Kepler's equation
	while ( (diff(E_last,E_latest)>1.0e-10)&&(count<100) )
	{
		E_last = E_latest;
		E_latest = E_last-((E_last-M-OCIT.e*sin(E_last))/(1.0-OCIT.e*cos(E_last)));
		count +=1;
		if (OCIT.verbose==true)
		{
			double Mnewton = (180.0/PI)*(E_latest-OCIT.e*sin(E_latest));
			cout << "\n$$$$$$$$$$$$$$$$$$$$$$$$$"<< "  Outputs for round "<< count <<"  $$$$$$$$$$$$$$$$$$$$$$"<<endl;
			cout << "E_latest found "<<E_latest*(180.0/PI) <<" causes Mnewton "<<Mnewton <<" != M_deg " <<OCRT.M_deg<<endl;
			cout << "current diff is = "<<diff(E_last,E_latest) <<endl;
			cout <<endl;
		}
	}
	cout<<"Just after the loop "<<endl;

	//Double check solution satisfies original equation
	//Then save the found Eccentric Anomaly if satisfied

	double Mnewton = (180.0/PI)*(E_latest-OCIT.e*sin(E_latest));
	if ( diff(OCRT.M_deg,Mnewton)>1.0e-5 )
	{
		cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
		cout << "PROBLEM: This resultant E does not satisfy the original equation, Newton's method FAILED !!!"<<endl;
		cout << "E_latest found "<<E_latest*(180.0/PI) <<" causes Mnewton "<<Mnewton <<" != M_deg " <<OCRT.M_deg <<endl;
		cout << fixed << "eccentricity = "<<OCIT.e<<endl;
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
	
	return OCRT;
}

int main()
{
	cout << "Hello World\n";

	orbitCalcInputType OCIT;

	// load OCIT structure
	OCIT.M = 2;
	OCIT.e = 0.77;
	OCIT.verbose = true;

	// instantiate the OCRT structure and pass the OCIT structure into
	// orbitCalculator to run the model and pass back the complete OCRT.
	orbitCalcReturnType OCRT;
	OCRT = orbitCalculator(OCIT);
}
