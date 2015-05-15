#include <vector>
#include <string>
using namespace std;

//Define parent class, general RV tools can go in here if any are needed
class RVtools
{
public:
	double Mass1; //in Msun
	double a1; //primary star's semi-major axis in AU
	double a_total;//Primary's semi-major PLUS secondary's semi-major [AU]
	bool verbose;
	double K;
	double period;
	double inclination_deg;
	double e;
	double TA_deg;
	double argPeri_deg;
	bool primaryStarRVs;
	//define public functions
	double VRcalculatorSemiMajorType();
};

class VRcalcStarStar : public RVtools
{
	/**
	This is the class for binary star calculations as it uses the more complex
	form of the RV equation that fully takes the secondary's mass into account.
	*/
public:
	//define public variables
	vector<double> epochs_s;
	double K_s;
	double e_s;
	double T_s;
	double Tc_s;
	double period_s; //in [yrs]
	double Mass2_s; //in Msun
	double argPeri_deg_s;
	double inclination_deg_s;
	double TA_deg_s;
	//define public functions
	vector<double> multiEpochCalc();
};

class VRcalcStarPlanet : public RVtools
{
	/**
	 This is the class for calculating the radial velocity of the primary star
	 due to a companion planet of significantly lower mass, and thus uses a
	 simplified version of the RV equation.
	 */
public:
	//define public variables
	vector<double> epochs_p;
	double K_p;
	double K_p_error;
	double e_p;
	double T_p;
	double Tc_p;
	double period_p; //in [yrs]
	double Mass2sinI_p; //in Msun
	double argPeri_deg_p;
	double inclination_deg_p;
	double TA_deg_p;
	double E_deg_p;
	//define public functions
	vector<double> multiEpochCalc();
};
