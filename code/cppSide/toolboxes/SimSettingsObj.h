#include <vector>
#include <string>
using namespace std;


class SimSettingsObj
{
public:
	// define params for obj
	// General simulation settings
	double chiSquaredMax;
	int numSamples;
	int numSamples_SimAnneal;
	int startTemp;
	int numSamplePrints;
	bool silent;
	bool verbose;
	string settings_and_InputDataDir;
	string SystemDataFilename;
	string DIdataFilename;
	string RVdataFilename;
	string data_dir;
	string data_filenameRoot;
	bool RVonly;
	bool DIonly;
	bool mcONLY;
	bool simAnneal;
	// advanced RV settings
	bool simulate_StarStar;
	bool simulate_StarPlanet;
	bool simulate_PrimaryOrbitRV;
	bool calcCorrLengths;
	bool CalcGelmanRubin;
	int numTimesCalcGR;
	bool TcStepping;
	bool primaryStarRVs;
	bool TcEqualT;
	bool argPeriPlusPiRV;
	bool argPeriPlusPiDI;
	double argPeriOffsetDI;
	double argPeriOffsetRV;


	// Ranges for acceptable random number inputs ######
	double longAN_degMIN; // [deg]
	double longAN_degMAX;// [deg]
	double eMIN;
	double eMAX;
	double a_totalMIN;
	double a_totalMAX;
	double periodMIN; // [yrs]
	double periodMAX; // [yrs]
	double inclination_degMIN; // [deg]
	double inclination_degMAX; // [deg]
	double argPeri_degMIN; // [deg]
	double argPeri_degMAX; // [deg]
	double T_Min;// [JD]
	double T_Max; //[JD]
	double K_MIN;// [m/s]
	double K_MAX;// [m/s]
	vector<double> RVoffsetMAXs;
	vector<double> RVoffsetMINs;

	// define funcs for obj
	void settingsLoadUp(const char* filename);
	void findArgPeriOffsets();
};

vector<double> rvOffsetsParser(string strIn);
bool boolParser(string sIn);
