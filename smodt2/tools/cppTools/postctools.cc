//@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
#include "postctools.h"

void PostCtools::loadParamData(double *zz, int zz_nx, int zz_ny){
	if (false)
			std::cout<<"\nInside loadData function"<<std::endl;
	data = zz;
	data_nx = zz_nx;
	data_ny = zz_ny;
};

void PostCtools::sumCalc(int startPoint,int lastPoint)
{
	sum = 0.0;
	//loop through all data points to get total value
	for (j=startPoint;j<(lastPoint+1);j++)
		sum+=data[colNum+j*data_ny];
};

void PostCtools::meanCalc(int startPoint, int lastPoint)
{
	/**
	 * A very simple function that simply adds up the total of all the elements
	 * in a vector<double> and divides it by the number of elements to give
	 * the mean value.  The precision of the calculation is checked to make sure
	 * the result is sufficiently precise.
	 */
	ave = 0;
	sumCalc(startPoint,lastPoint);
	ave = sum/(double(lastPoint-startPoint+1));
};

double PostCtools::varianceCalc(int startPoint, int lastPoint)
{
	/**
	 * This will calculate the "bias-corrected sample variance"
	 * and uses an advanced "corrected two-pass algorithm" to reduce roundoff error,
	 * a modified version of function on pg 728 of Numerical Recipes 3rd.
	 */
	n=lastPoint+1;//add 1 to lastPoint to make it 1 indexed, instead of zero indexed
	meanCalc(startPoint,lastPoint);
	var=ep=0.0;
	//loop through all points to load up sums needed for "corrected two-pass algorithm", eq 14.1.8 pg 724
	for (j=startPoint;j<n;j++)
	{
		s=data[colNum+j*data_ny]-ave;
		ep+=s;
		var+=s*s;
	}
	//calc var and return it
	var=(var-ep*ep/n)/(n-1);
	return var;
};

double PostCtools::corrLenCalc(int parNum){
	/**
	 * Calculates the average correlation length
	 * of the input parameter's data in a boxcar style.
	 */
	colNum = parNum;
	varALL = varianceCalc(0,(data_nx-1));
	halfVarALL = varALL/2.0;
	numCorrLengths=0.0;
	i_last=0;
	for (i=0;i<data_nx;i++){
		var = varianceCalc(i_last,i);
		if (var>halfVarALL){
			corrLengthsTotal+=(i-i_last+1);
			i_last=i;
			numCorrLengths+=1.0;
			break;
		}
	}
	return corrLengthsTotal/numCorrLengths;
};

