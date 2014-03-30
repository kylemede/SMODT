#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "RVtools.h"
#include "DItools.h"

using namespace std;

int main()
{
	VRcalcStarStar VRss;
	VRcalcStarPlanet VRsp;
	DItools DIt;

	VRss.epochs_s.push_back(1.0);
	VRsp.epochs_p.push_back(1.0);
	VRsp.e_p = VRss.e_s = 0.4;
	VRsp.T_p = VRss.T_s = 233443344.4;
	VRsp.period_p = VRss.period_s = 20.09; //in [yrs]
	VRsp.Mass1 = VRss.Mass1 = 1.2; //in Msun
	VRsp.Mass2sinI_p = VRss.Mass2_s = 0.3; //in Msun
	VRsp.argPeri_deg_p = VRss.argPeri_deg_s = 88.0;
	VRss.inclination_deg_s = 30.0;
	VRsp.verbose=VRss.verbose = false;

	DIt.SAs_arcsec_observed.push_back(1.1);
	DIt.SA_errors.push_back(1.1);
	DIt.PAs_deg_observed.push_back(1.1);
	DIt.PA_errors.push_back(1.1);
	DIt.epochs.push_back(1341235234.2);
	DIt.Sys_Dist_PC = 12;
	DIt.inclination_deg = 32;
	DIt.longAN_deg = 12;
	DIt.e = 0.3;
	DIt.T = 123214354;
	DIt.period = 10;
	DIt.argPeri_deg = 23.0;
	DIt.Mass1 = 1.2; //must be in Msun
	DIt.Mass2 = 0.4; //must be in Msun
	DIt.verbose = false;

	cout<<"\n*******************************"<<endl;
	cout<<"\nfirst calling starstar version"<<endl;
	vector<double> vrss = VRss.multiEpochCalc();

	cout<<"\nnext calling starplanet version"<<endl;
	vector<double> vrsp = VRsp.multiEpochCalc();

	cout<<"\nnext calling DI multiepoch calc"<<endl;
	multiEpochOrbCalcReturnType MEOCRT = DIt.multiEpochOrbCalc();

	cout<<"\n*******************************"<<endl;


}
