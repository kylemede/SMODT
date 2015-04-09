#include <vector>
#include <string>
#include "rnd/kylesGaussRand.h"
using namespace std;


class simAnealOrbFuncObj
{
public:
	// Input structs
	string SSlogStr;
	SimSettingsObj SSO;
	DItools DIt;
	DIdataObj DIdo;
	RVdataObj RVdo;
	DataObj SYSdo;
	multiEpochOrbCalcReturnType MEOCRT;
	VRcalcStarPlanet VRCsp;
	VRcalcStarStar VRCss;
	double chiSquaredMin;
	double one_over_nu_RV;
	double one_over_nu_DI;
	double one_over_nu_TOTAL;
	int bestOrbit;
	int timesBeenHereTotal;
	int saveEachInt;
	double earliestEpoch;
	int randSeed;
	outputDataType ODT;
	double tempStepSizePercent;
	double startTemp;
	double sigmaPercent;
	double sigmaPercent_latest;
	int numSamples_SA;
	double TMIN;
	double TMAX;
	double inclination_deg_sigma;
	double e_sigma ;
	double a_total_sigma;
	double longAN_deg_sigma ;
	double period_sigma ;
	double argPeri_deg_sigma ;
	double T_sigma ;
	double K_sigma;
	double sqrtESinomega_sigma;
	double sqrtECosomega_sigma;
	vector<double> offset_sigmas;
	bool vary_K;
	int numParams;
	int numDIparams;
	int numRVparams;
	vector<int> paramsToVaryIntsAry;

	// funcs available to this object
	void simulator();


};
