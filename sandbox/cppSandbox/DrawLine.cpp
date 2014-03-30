#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

int main (void)   // this is the line drawing function
{
	vector<int> A ;


	for ( int l=0; l<10; l++)
		A.push_back(l);



    for (int i = 0; i < 10; i++)
    {
    	int u = A[i];

    	cout << "value of "<< i << " element of A is "<< u<< "\n";
    }
    return 0;
}
