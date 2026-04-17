#include "main.h"
#include <thread>

int main(int argc, char **argv){          
    // Read input and output paths
	string filePath = argv[1];	
	string pathAndFileout = argv[2];

	Instance inst = readInstance(filePath);
	Solution sol;
				
	inst.timeLimit = 3600;
			
	TSP_DLYD(inst, sol);

	sol.timeT = getCPUTime() - inst.startTimeSaved;

	printInfo(pathAndFileout, inst, sol);

}
