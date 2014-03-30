#include <vector>
#include <string>
using namespace std;

class DataObj
{
public:
	// General system data
	double Sys_Dist_PC;
	double Sys_Dist_PC_error;
	double Mass1;
	double Mass1_error;
	// Orbital params for planet
	// (for if investigating companion star's orbit)
	double planet_K;
	double planet_K_error;
	double planet_e;
	double planet_e_error;
	double planet_T;
	double planet_T_error;
	double planet_Tc;
	double planet_Tc_error;
	double planet_P;
	double planet_P_error;
	double planet_MsinI;
	double planet_MsinI_error;
	double planet_argPeri;
	double planet_argPeri_error;
	double planet_inc;
	double planet_inc_error;
	double planet_long_AN;

	// Orbital params for companion star
	// (for if investigating orbit of planet around primary star)
	double star_e;
	double star_e_error;
	double star_T;
	double star_T_error;
	double star_Tc;
	double star_Tc_error;
	double star_P;
	double star_P_error;
	double star_Mass2;
	double star_Mass2_error;
	double star_argPeri;
	double star_argPeri_error;
	double star_inc;
	double star_inc_error;
	double star_long_AN;
	// functions available to all lower level data objects
	void systemDataLoadUp(const char* filename);
};

class DIdataObj : public DataObj
{
public:
	//define the params/variables for the class
	double numEpochs_DI;
	vector<double> epochs_DI;
	vector<double> SAs_arcsec_observed;
	vector<double> SA_errors;
	vector<double> PAs_deg_observed;
	vector<double> PA_errors;
	//define any functions for this class
	void dataLoadUp(const char* filename);
};

class RVdataObj : public DataObj
{
public:
	//define the params/variables for the class
	//keep datasets separate as they have different RVoffsets
	double numEpochs_RV;
	vector<vector<double> > epochs_RV;
	vector<vector<double> > RVs;
	vector<vector<double> > RV_inv_var; //These should include the jitter, else if jitter is defined it will be added in quadrature (ie. [(inv_var)^-1+jitter^2]^-1
	//double jitter; //must be integrated into the RV errors if not done all ready.  Set=0, means it is all ready included, else added in quadrature
	//define any functions for this class
	void dataLoadUp(const char* filename);

};
