#include "RVtools.h"
#include "DItools.h"
#include "SimSettingsObj.h"
#include "DataObj.h"
#include "generalTools.h"
#include "rnd/sfmt.h"
#include <vector>
/*! \file */
// conversion factors and constants
#define PI 3.14159265359
#define SecPerYear  31557600.0
#define GravConst  6.67300e-11
#define MperAU  149598000000.0
#define KGperMsun  1.98892e30
#define SecPerDay 86400.0
#define KGperMjupiter 1.8986e27
#define MjupiterPerMsun 0.000954265748

////OrbSimGeneralToolbox declarations
//TAcalcReturnType TAcalculator(TAcalcInputType TACIT);
//double chiSquaredCalc(double real, double error, double model);
//void fileWriter(outputDataType ODT);
//DItools DItoolsParamLoadUp(DIdataObj DIdo);
//VRcalcStarStar VRcalcStarStarLoadUp(RVdataObj RVdo);
//VRcalcStarPlanet VRcalcStarPlanetLoadUp(RVdataObj RVdo);
//semiMajorType semiMajorConverter(semiMajorType SMT);
//double earliestEpochFinder(DIdataObj DIdo, RVdataObj RVdo);
//string CorrelationLengthCalc(vector<double> data, string paramName, vector<int> timesBeenHere);
//string numSamplesStringMaker(int numSamples);
//void logFileWriter(string filename, string LOGlines);
//string filenameEndAppend(string inputString, string endAppendString);
//string filenamePrepend(string inputString, string prependString);
//string boolToStr(bool boolIN);
//double standardDeviation(vector<double> v);
//string timeStr(int timeElapsed);
//int corrLengthJumpyCalc(vector<double> data);
//double varianceCalc(vector<double> data, int lastPoint);
//eccArgPeriCalcType eccArgPeriCalc(eccArgPeriCalcType EACT);
//eccArgPeri2ToTcType eccArgPeri2ToTcCalc(eccArgPeri2ToTcType EATT);
//double atanTopBtm(double top, double btm);
//double meanCalc(vector<double> v, int lastPoint);
//void gelmanRubinStage1(outputDataType ODT,int numTimes);
//GRfuncReturnType gelmanRubinStage1func(vector<double> data,vector<int> timesBeenHere,int numTimes);
//string fileBasename(string str);
//string findDirectory(string str);
//string findChainNumberStr(string str);
//vector<double> fillOutDataVector(vector<double> data,vector<int> timesBeenHere);
//outputDataType odtStart(outputDataType ODT, int numTotalSamples);
//outputDataType odtFinish(outputDataType ODT);
////string	confidenceLevelFinder(vector<double> data,vector<int> timesBeenHere);
