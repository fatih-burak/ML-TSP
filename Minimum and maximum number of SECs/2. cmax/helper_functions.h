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
#include "gurobi_c++.h"
#include <set>
#include <math.h> 
#include <iomanip>  // For setprecision and fixed
#include <unordered_map>

const double EPSILON = 0.00001; // small constant

struct Instance
{
	string name;
	int nbCity; 						// number of nodes
	vector<int> V;					// set of vertices
	vector<double> xCoord;	// x-Coords
	vector<double> yCoord;	 // y-Coords
	vector<vector<double>> dist; //distance between nodes
	double timeLimit = 1200;
	double startTimeSaved;
	set<set<int>> uniqueSubsets;
	vector<set<set<int>>>C_of_k; //for a given solution k, subtours in the solution
	vector<set<set<int>>>U_of_k; //for a given solution k, subtours in the solution
};

struct Solution
{
	//General
	int opt = 0, numSol = 0;
	double timeT = 0.0;
	double optimalObjVal = 0.0, LB = 0.0;
	vector<vector<int>> solution;
	vector<vector<vector<pair<int, int>>>> cuts;
};

double getCPUTime();
void printInfo(const string& pathAndFileout, const Instance& inst, const Solution& sol);
Instance readInstance(const string& filePath, const string& allSolsPath);
vector<int> setDifference(const vector<int>& vertices, const vector<int>& cycle);
#endif 
