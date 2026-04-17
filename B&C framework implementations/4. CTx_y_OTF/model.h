#ifndef MODEL_H
#define MODEL_H

using namespace std;
#include "gurobi_c++.h"
#include "helper_functions.h"

void TSP_DLYD(Instance& inst, Solution& sol);
void ILP(Instance& inst, Solution& sol, int cutType);

vector<vector<int>> findCuts(InterSol& iSol);

class SolutionCollector : public GRBCallback {
public:
    Instance& inst;  // reference to inst
    Solution& sol; // reference to sol
    vector<vector<GRBVar>>& x;         // reference to variable list

    SolutionCollector(vector<vector<GRBVar>>& xvars , Instance& instObj, Solution& solObj) : x(xvars), inst(instObj), sol(solObj) {}

protected:
    void callback();
};


#endif 

