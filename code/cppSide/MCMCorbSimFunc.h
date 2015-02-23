#include <vector>
#include <string>
using namespace std;


class MCMCorbFuncObj
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
	int timesBeenHereTotal;
	int bestOrbit;
	double earliestEpoch;
	int randSeed;
	outputDataType ODT;
	double start_longAN;
	double start_e;
	double start_a_total;
	double start_T;
	double start_Tc;
	double start_period;
	double start_inc_deg;
	double start_argPeri;
	vector<double> start_offsets;
	double start_K;
	int numParams;
	int numDIparams;
	int numRVparams;
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
	vector <int> paramsToVaryIntsAry;

	double sigmaPercent;
	double paramChangeStepSizePercent;
	// funcs available to this object
	void simulator();
};
