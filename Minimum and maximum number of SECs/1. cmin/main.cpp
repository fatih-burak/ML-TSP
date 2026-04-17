#include "main.h"
#include <thread>
#include <windows.h>

int main(int argc, char** argv) {
	string path = "path_to_instance_folder";
	string instancePath = path + "INSTANCES\\SELECTED\\ALL\\";		
	string foundInstancesPath = path + "RESULTS\\C_MIN\\ALL_SOL_FOLDER\\"; //path to the saved solutions
	string pathAndFileout = path + "RESULTS\\C_MIN\\results.txt";		//output path

	string searchPath = foundInstancesPath + "*"; // wildcard to match all files
	WIN32_FIND_DATAA findFileData;
	HANDLE hFind = FindFirstFileA(searchPath.c_str(), &findFileData);

	if (hFind == INVALID_HANDLE_VALUE) {
		cerr << "Failed to open directory: " << foundInstancesPath << std::endl;
		return 1;
	}

	do {
		string fileName = findFileData.cFileName;
		if (fileName == "." || fileName == "..") continue;
		if (!(findFileData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
			cout << fileName << endl;
			size_t pos = fileName.find(".txt");
			string baseName = fileName.substr(0, pos);
			cout << "Base name: " << baseName << std::endl;

			string filePath = instancePath + baseName + ".tsp";

			string allSolsPath = path + "RESULTS\\C_MIN\\ALL_SOL_FOLDER\\" + baseName + ".txt";
			cout << allSolsPath << endl;
			Instance inst = readInstance(filePath, allSolsPath);
			
			Solution sol;
			
			inst.timeLimit = 3600;
			
			C_MIN_MAIN(inst, sol);

			sol.timeT = getCPUTime() - inst.startTimeSaved;

			printInfo(pathAndFileout, inst, sol);


		}

	} while (FindNextFileA(hFind, &findFileData) != 0);

	FindClose(hFind);
	return 0;
}
