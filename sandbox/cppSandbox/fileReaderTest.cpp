#include<stdio.h>
#include<stdlib.h>
#include <sstream>
#include <string.h>
#include <iostream>
#include <fstream>
#include "Toolboxes/DataObj.h"
#include "Toolboxes/SimSettingsObj.h"

using namespace std;

int main()
{
//	const char* filenameRV = "/run/media/Kyle/Data1/Todai_Work/Dropbox/workspace/Binary-project/C++Stuff/SimSettings_and_InputData/RVdata.dat";
//	const char* filenameDI = "/run/media/Kyle/Data1/Todai_Work/Dropbox/workspace/Binary-project/C++Stuff/SimSettings_and_InputData/DIdata.dat";
//	const char* filenameSys = "/run/media/Kyle/Data1/Todai_Work/Dropbox/workspace/Binary-project/C++Stuff/SimSettings_and_InputData/SystemData.txt";
	const char* filenameSimSet = "/run/media/Kyle/Data1/Todai_Work/Dropbox/workspace/Binary-project/C++Stuff/SimSettings_and_InputData/SimSettings.txt";
//	RVdataObj RVdo;
//	DIdataObj DIdo;
	SimSettingsObj SSo;
//	DataObj Do;
//	RVdo.dataLoadUp(filenameRV);
//	DIdo.dataLoadUp(filenameDI);
//	Do.systemDataLoadUp(filenameSys);
	SSo.settingsLoadUp(filenameSimSet);




}
