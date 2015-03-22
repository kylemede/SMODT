#include <iostream>
#include <vector>
#include <math.h>
#include <time.h>
#include <sstream>
#include <string>
#include <fstream>

#include "headerTest.hpp"
//#include "../rnd/randomc.h" //a library with advanced uniform random number generators
//#include "../rnd/sfmt.cpp"
#include "rnd/stocc.h"
//#include "rnd/stoc1.cpp"

using namespace std;


void func(int A,int B,int C)
{

	cout<< "A = "<<A<< endl;
	cout<< "B = "<<B<<endl;
	cout<< "C = "<<C<< endl;
}
//
//void func2(CRandomSFMT1 RanGen)
//{
//	for(int i=0;i<10;++i)
//	{
//		double randNumber = RanGen.UniformRandom(1.0, 100.0);
//		cout<<"Random number is: "<<randNumber<<endl;
//	}
//
//}

string filenameEndAppend(string inputString, string endAppendString)
{
	bool verboseInternal = true;
	string outputString = "";
	int periodPos = inputString.find('.');

	if (periodPos>0)
	{
		for (int curPos=0; curPos<inputString.size();curPos+=1)
		{
			if (false)
				cout<<" char is:"<<inputString[curPos]<<" , for curPos "<<curPos<<endl;
			if (inputString[curPos]!='.')
				outputString.push_back(inputString[curPos]);
			else
				outputString += endAppendString+".";
		}
	}
	else
		outputString += inputString+endAppendString;
	if (verboseInternal)
		cout<<"\nFilename modified by filenameEndAppend:: \nInput string: "<<inputString<<"\nOutput string:"<<outputString<<endl;

	return outputString;
}

string filenamePrepend(string inputString, string prependString)
{
	bool verboseInternal = true;
	string outputString = "";
	int lastSlashPos = inputString.rfind('/');

	if (verboseInternal)
		cout<<"\nlastSlashPos = "<<lastSlashPos<<endl;

	if (lastSlashPos>0)
	{
		for (int curPos=0; curPos<inputString.size();curPos+=1)
		{
			if (false)
				cout<<" char is:"<<inputString[curPos]<<" , for curPos "<<curPos<<endl;
			if (curPos!=lastSlashPos)
				outputString.push_back(inputString[curPos]);
			else
				outputString += "/"+prependString;
		}
	}
	else
		outputString += prependString+inputString;

	if (verboseInternal)
		cout<<"\nFilename modified by filenamePrepend:: \nInput string: "<<inputString<<"\nOutput string:"<<outputString<<endl;

	return outputString;

}
int main(int argc ,char *argv[])
{
	cout << "Hello World\n";

	cout<< "about to enter random number generation block\n"<<endl;
	timespec time1;
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
	int time_nsec=time1.tv_nsec;
	//cout<<time_nsec<<endl;

	//choose generator
	//CRandomSFMT1 RanGen(time_nsec);
	StochasticLib1 RanGen2(time_nsec);

	vector<double> randNumbersAry;

	for(int i=0;i<100;++i)
	{
		//double randNumber = RanGen.UniformRandom(1.0, 100.0);
		double randNumber = RanGen2.NormalTrunc(0.0,10.0,100.0);
		if (true)
			cout<<"Random number is: "<<randNumber<<endl;
		randNumbersAry.push_back(randNumber);
	}

	ofstream file;
	string data_filename = "/run/media/Kyle/Data1/Todai_Work/Data/data_Binary/data_Duo/TestCplusplusTrucNormal.dat";
	cout<<"\n***************************************************************"<<endl;
	cout<<"Writing file to: "<<data_filename<<endl;
	cout<<"***************************************************************\n"<<endl;
	file.open(data_filename.c_str()) ;     //open the output file

	for (int r=0; r<randNumbersAry.size();r++)
	{
		file<<randNumbersAry[r]<<endl;
		file.close();
	}


//
//	func2(RanGen);
//
//	cout<< "Done random number generation block\n"<<endl;
//	// Arguments from the command line are white space delimited!!
//
//	cout <<argc<<" arguments were passed in."<<endl;
//	cout << "They were: "<<endl;
//	for (int nArg=0; nArg < argc; nArg++)
//	{
//		cout << nArg << " " << argv[nArg] << endl;
//		cout << "arg as string: "<<string(argv[nArg])<<endl;
//	}
//
//	vector<vector<double> > lst;
//	for (int i=0; i<5; ++i)
//	{
//		vector<double> lstCur;
//
//		for (int j = 0; j<3;++j)
//			lstCur.push_back(1.234);
//
//		lst.push_back(lstCur);
//	}
//
//	//double temptest[] = {{1,2,3},{4,5,6}};
//	//int sizeoftemptest = int(sizeof(temptest)/sizeof(double));
//	//cout<< "sizeoftemptest = "<<sizeoftemptest<<endl;
////	cout<< "lst.size() = "<<lst.size()<<endl;
////	cout<< "lst[0].size() = "<<lst[0].size()<<endl;
//	std::stringstream ss;
//
////	double something;
////	ss<<'notADouble';
////	ss>>something;
////	cout<<"The value of something is:"<<something<<endl;
//
//	//int Times  = 4;
//	//int t;
//	//t= headerTest(Times);
//	//t = a.tester(Times);
//	//cout<<t<<endl;
//
//	cout<<"\n\n **** String manipulation test ****\n"<<endl;
//
//	string outStr = filenameEndAppend("stringTEST.txt","_addThisText_");
//	string outStr2 = filenameEndAppend("stringTEST2","_addThisText_");
//	string outStr3 = filenamePrepend("this/is/the/dir/stringTEST.txt", "_addThisText_");
//	string outStr4 = filenamePrepend("stringTEST.txt", "_addThisText_");
//
//	//func(1,2,3);
}


