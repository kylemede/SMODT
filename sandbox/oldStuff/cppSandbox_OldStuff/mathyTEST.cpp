#include <iostream>
#include <math.h>
using namespace std;

int main(int argc ,char *argv[])
{

	int numRuns = int(argv[0]);
	//int numRuns = int(1e6);
	int printCounter = int(double(numRuns)/10.0);
	int CurPrintCount = 0;
	int printPercent = 0;
	cout<< "**starting to run mathy loop "<< numRuns<< " many times"<<endl;
	for ( int i=0; i<numRuns ;++i )
	{
		for ( int j=0; j<numRuns;++j)
			sin(0.3);

		if ( i==CurPrintCount )
		{
			cout<< printPercent <<"% complete     "<<endl;
			printPercent = printPercent+10;
			CurPrintCount = CurPrintCount+printCounter;
		}
	}

	cout<<"**Done mathy loop" <<endl;
}


