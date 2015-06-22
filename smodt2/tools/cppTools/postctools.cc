//@Author: Kyle Mede, kylemede@astron.s.u-tokyo.ac.jp
#include "postctools.h"

void PostCtools::loadParamData(double *zz, int zz_nx, int zz_ny){
	//std::cout<<"\nInside loadData function"<<std::endl;
	data = zz;
	data_nx = zz_nx;
	data_ny = zz_ny;
};

void PostCtools::sumCalc(int startPoint,int lastPoint)
{
	sum = 0.0;
	//loop through all data points to get total value
	for (j=startPoint;j<(lastPoint+1);j++){
		if ((isnan(sum))||(sum>1e200))
			break;
		else{
			sum+=data[colNum+j*data_ny];
			//std::cout<<"colNum = "<< colNum<<std::endl;
			//std::cout<<"j = "<< j<<std::endl;
			//std::cout<<"data[colNum+j*data_ny] = "<< data[colNum+j*data_ny]<<std::endl;
			//std::cout<<"sum = "<< sum<<std::endl;
		}
	}
	//std::cout<<"sum = "<< sum<<std::endl;
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
	//std::cout<<"ave = "<< ave<<std::endl;
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
	//std::cout<<"ave = "<< ave<<std::endl;
	var=ep=0.0;
	if (n>1){
		//loop through all points to load up sums needed for "corrected two-pass algorithm", eq 14.1.8 pg 724
		for (j=startPoint;j<n;j++)
		{
			s=data[colNum+j*data_ny]-ave;
			//std::cout<<"s = "<< s<<std::endl;
			ep+=s;
			//std::cout<<"ep = "<< ep<<std::endl;
			var+=s*s;
			//std::cout<<"var = "<< var<<std::endl;
		}
		//calc var and return it
		var=(var-ep*ep/n)/(n-1);
		//std::cout<<"ep = "<< ep<<std::endl;
		//std::cout<<"n = "<< n<<std::endl;
		//std::cout<<"var = "<< var<<std::endl;
	}
	return var;
};

double PostCtools::corrLenCalc(int parNum){
	/**
	 * Calculates the average correlation length
	 * of the input parameter's data in a boxcar style.
	 */
	colNum = parNum;
	//std::cout<<"colNum = "<< colNum<<std::endl;
	//std::cout<<"data_nx = "<< data_nx<<std::endl;
	//std::cout<<"data_ny = "<< data_ny<<std::endl;
	varALL = varianceCalc(0,data_nx);
	std::cout<<"varALL = "<< varALL<<std::endl;
	halfVarALL = varALL/2.0;
	std::cout<<"halfVarALL = "<< halfVarALL<<std::endl;
	numCorrLengths=0.0;
	i_last=0;
	for (i=0;i<data_nx;i++){
		var = varianceCalc(i_last,i);
		if (var>halfVarALL){
			corrLengthsTotal+=double(i-i_last+1);
			i_last=i;
			numCorrLengths+=1.0;
		}
	}
	std::cout<<"corrLengthsTotal = "<< corrLengthsTotal<<std::endl;
	std::cout<<"numCorrLengths = "<< numCorrLengths<<std::endl;
	return corrLengthsTotal/numCorrLengths;
};

