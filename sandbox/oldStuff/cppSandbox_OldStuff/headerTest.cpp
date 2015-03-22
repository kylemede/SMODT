#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "headerTest.hpp"
using namespace std;


int headerTest(int Times)
{
	int i = 0;
	while ( i<Times )
	{
		cout<< "i = "<<i<<endl;
		i=i+1;
	}
	return i;
}
