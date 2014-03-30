#include <vector>
#include <string>
using namespace std;

struct orbitCalcReturnType
{
	double x_model;
	double y_model;
};

struct multiEpochOrbCalcReturnType
{
	double chi_squared_total;
	double a_total; //must be in AU
	//vector<double> x_models;
	//vector<double> y_models;
};

class DItools
{
public:
	vector<double> SAs_arcsec_observed;
	vector<double> SA_errors;
	vector<double> PAs_deg_observed;
	vector<double> PA_errors;
	vector<double> epochs_DI;
	double Sys_Dist_PC;
	double inclination_deg;
	double longAN_deg;
	double e;
	double T;
	double period;
	double argPeri_deg;
	double Mass1; //must be in Msun
	double Mass2; //must be in Msun
	double E_deg;
	double a_total;
	double a_use;
	bool verbose;

	orbitCalcReturnType orbitCalculator();
	multiEpochOrbCalcReturnType multiEpochOrbCalc();
};
