#include "main.h"
#include <thread>
#include <windows.h>

int main(int argc, char** argv) {
	string path = "path_to_the_instance_folder";
	string instancePath = path + "INSTANCES\\SELECTED\\ALL\\";
	string foundInstancesPath = path + "RESULTS\\C_MIN\\ALL_SOL_FOLDER\\";
	string pathAndFileout = path + "RESULTS\\C_MIN\\results.txt";

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

			string filePath = instancePath + baseName + ".tsp";
			string allSolsPath = path + "RESULTS\\C_MIN\\ALL_SOL_FOLDER\\" + baseName + ".txt";

			Instance inst = readInstance(filePath, allSolsPath);
			
			Solution sol;
			
			inst.timeLimit = 3600;
			
			C_MAX_MAIN(inst, sol);

			sol.timeT = getCPUTime() - inst.startTimeSaved;

			printInfo(pathAndFileout, inst, sol);


		}

	} while (FindNextFileA(hFind, &findFileData) != 0);

	FindClose(hFind);
	return 0;
}
