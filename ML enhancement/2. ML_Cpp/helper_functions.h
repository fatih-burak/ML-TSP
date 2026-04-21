#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

using namespace std;
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream> 
#include <time.h>
//#include <ilcp/cp.h>
#include "gurobi_c++.h"
#include <set>
#include <math.h> 
#include <iomanip>  // For setprecision and fixed
#include <tuple>


const double EPSILON = 0.00001; // small constant

struct Instance
{
	string name;
	string fileName;
	string saveFolder;
	bool isSolved;

	int nbCity; 						// number of nodes
	vector<int> V;					// set of vertices
	vector<double> xCoord;	// x-Coords
	vector<double> yCoord;	 // y-Coords
	vector<vector<double>> dist; //distance between nodes
	double timeLimit = 1200;
	double startTimeSaved;
	int howManySolutions = GRB_MAXINT;
	string cutFolder;
	vector<vector<int>> precuts;

	// CYCLE GEN RELATED
	// MST part
	double tMST_cycle_Gen = 0.0; //+
	int nbMST_cycles = 0;           //+
	int nbMSTmatch = 0;
	double percMSTmatch = 0.0;

	// RC part
	double tLPrelaxation = 0.0;//+
	double tBuildingBlocksFromSkeleton = 0.0;//+
	double tCycleXOR = 0.0;//+
	double tGlobalRCMerge = 0.0;//+
	double tAllRC = 0.0;//+

	int nbRC_cycles_NOT_RED = 0;
	int nbRC_cycles_unique = 0;//+

	int nbRCmatch = 0;
	double percRCmatch = 0.0;

	//Both RC and MST
	double tMergeAndDedupe = 0;
	int nbMatch = 0;
	double percMatch = 0.0;

	//feature related
	// Root-LP x-values (only upper triangle used, [i][j] for i<j)

	std::vector<std::vector<double>> lp_x;

	// Root-LP reduced costs (only upper triangle used, [i][j] for i<j)
	std::vector<std::vector<double>> lp_rc;

	double lp_obj;

	// Distance scale for normalization (median distance or MST-based)
	double dist_scale = 1.0;
};


struct InterSol
{
	double objVal = 0.0; //stores the current objective value
	vector<pair<int, int>> sol_edges; //stores the current solution edges
	vector<vector<pair<int, int>>> MULTI_sol_edges; //stores all incumbent solution edges
	vector<double> MULTI_objVal; //stores all incumbent solution's objVals
	int nbCity; //nb of vertices in the TSP instance
	float perc = 1.0; //theshold for SEC size
	vector<vector<int>> shortest_cycle; //stores the shortest cycle for each incumbent solution
	bool isTimeLimitHit = false;
	bool isShortest = false; // false: regular, true: shortest cut added
	string cutSaveLoc;
};

struct Solution
{
	//General
	int opt = 0, iter = 0;
	double timeT = 0.0;
	double optimalObjVal = 0.0, LB = 0.0;
	set<set<int>> cuts;
	vector<vector<int>> solution;
	int eliminatedSECs = 0;
	int needed = 0, RCadded = 0, MSTadded = 0;
	double tNN_pipeline = 0.0;
};

double getCPUTime();
Instance readInstance(const string& filePath);
vector<int> setDifference(const vector<int>& vertices, const set<int>& cycle);
void printInfo(const string& pathAndFileout, const Instance& inst, const Solution& sol);

#endif 
