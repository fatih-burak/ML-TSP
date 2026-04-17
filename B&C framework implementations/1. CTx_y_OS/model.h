#ifndef MODEL_H
#define MODEL_H

using namespace std;
#include "gurobi_c++.h"
#include "helper_functions.h"

void TSP_DLYD(Instance& inst, Solution& sol);
InterSol ILP(Instance& inst, Solution& sol, int cutType);

vector<vector<int>> findCuts(InterSol& iSol);



#endif 

