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

#define PI 3.14159265

using namespace std;


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

double mean(vector <double> ary)
{
	double sum = 0;
	int N = ary.size();

	for ( int i = 0; i<N; i++)
	{
		sum = sum + ary[i];
	}
	double m = sum/(double)N;
	return m;
}


int main()
{
	int numSamples = 100000;

	double min = 2448964.602000;
	double max = 2452519.702000;
	vector <double> numbers;
	for ( int num=0;num<numSamples ;num++ )
	{
		double randNum = uniformRandNumber(min,max);
		numbers.push_back(randNum);
		//cout<< fixed << randNum <<endl;

	}
	cout << fixed<<"Mean of rand Numbers generated = "<< mean(numbers)<<endl;

	return EXIT_SUCCESS;

}
