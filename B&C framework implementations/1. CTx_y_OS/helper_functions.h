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
	string cutFolder;
};


struct InterSol
{
	double objVal = 0.0; //stores the current objective value
	vector<pair<int, int>> sol_edges; //stores the current solution edges
	int nbCity; //nb of vertices in the TSP instance
	float perc = 1.0; //theshold for SEC size
	vector<int> shortest_cycle; //stores the shortest cycle
	bool isTimeLimitHit = false;
	bool isShortest = false; // false: regular, true: shortest cut added
};

struct Solution
{
	//General
	int opt = 0, iter = 0;
	double timeT = 0.0;
	double optimalObjVal = 0.0,  LB = 0.0;
	vector<vector<int>> cuts;
	vector<vector<int>> solution;
};

double getCPUTime();
Instance readInstance(const string& filePath);
vector<int> setDifference(const vector<int>& vertices, const vector<int>& cycle);
void printInfo(const string& pathAndFileout, const Instance& inst, const Solution& sol);
#endif 
