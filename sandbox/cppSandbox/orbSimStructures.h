#include <vector>
#include <string>
using namespace std;

struct outputDataType
{
	string data_filename;
	int numEpochs_DI;
	int numSamplesAccepted;
	//inputs
	vector<double> longAN_degs;
	vector<double> es;
	vector<double> Ts;
	vector<double> Tcs;
	vector<double> periods;
	vector<double> inclination_degs;
	vector<double> argPeri_degs;
	vector<double> a_totals;
	vector<vector<double> > RVoffsets;
	vector<double> Ks;
	//outputs
	vector<double> chiSquareds ;
	vector<int> timesBeenHeres;
};

struct TAcalcInputType
{
	double t;
	double e;
	double T;
	double Tc;
	double period;
	bool verbose;
};

struct TAcalcReturnType
{
	double E_deg;
	double TA_deg;
};

struct semiMajorType
{
	double a1;
	double a2;
	double a_total;
	double Mass1;
	double Mass2;
	double period;
};

struct eccArgPeriCalcType
{
	double e;
	double argPeri_deg;
	double sqrtEsinArgPeri;
	double sqrtEcosArgPeri;
};

struct eccArgPeri2ToTcType
{
	double e;
	double argPeri_deg;
	double period;
	double To;
	double Tc;
};

struct GRfuncReturnType
{
	vector<int> Lcs;
	vector<double> means;
	vector<double> vars;
};


