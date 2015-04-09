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
	//vector<int> timesBeenHeres;
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

class generalTools
{
public:
	TAcalcReturnType TAcalculator(TAcalcInputType TACIT);
	double chiSquaredCalc(double real, double error, double model);
	void fileWriter(outputDataType ODT);
	DItools DItoolsParamLoadUp(DIdataObj DIdo);
	VRcalcStarStar VRcalcStarStarLoadUp(RVdataObj RVdo);
	VRcalcStarPlanet VRcalcStarPlanetLoadUp(RVdataObj RVdo);
	semiMajorType semiMajorConverter(semiMajorType SMT);
	double earliestEpochFinder(DIdataObj DIdo, RVdataObj RVdo);
	string CorrelationLengthCalc(vector<double> data, string paramName);//, vector<int> timesBeenHere);
	string numSamplesStringMaker(int numSamples);
	void logFileWriter(string filename, string LOGlines);
	string filenameEndAppend(string inputString, string endAppendString);
	string filenamePrepend(string inputString, string prependString);
	string boolToStr(bool boolIN);
	double standardDeviation(vector<double> v);
	string timeStr(int timeElapsed);
	int corrLengthJumpyCalc(vector<double> data);
	double varianceCalc(vector<double> data, int lastPoint);
	eccArgPeri2ToTcType eccArgPeri2ToTcCalc(eccArgPeri2ToTcType EATT);
	double atanTopBtm(double top, double btm);
	double meanCalc(vector<double> v, int lastPoint);
	double sumCalc(vector<double> v,int lastPoint);
	int sumIntCalc(vector<int> v,int lastPoint);
	void gelmanRubinStage1(outputDataType ODT,int numTimes);
	GRfuncReturnType gelmanRubinStage1func(vector<double> data,int numTimes);//vector<int> timesBeenHere,int numTimes);
	string fileBasename(string str);
	string findDirectory(string str);
	string findChainNumberStr(string str);
	//vector<double> fillOutDataVector(vector<double> data);//,vector<int> timesBeenHere);
};
