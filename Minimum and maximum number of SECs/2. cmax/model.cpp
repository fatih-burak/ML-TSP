#include "model.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <set>

int isSubsetInASol(set<int>& subset, vector<set<int>>& C) {
	int index = -1;

	for (size_t i = 0; i < C.size(); ++i) {
		if (C[i] == subset) {
			index = static_cast<int>(i);
			break;
		}
	}
	return index;
}

void C_MAX_MAIN(Instance& inst, Solution& sol) {
	// create a model
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 1);
	GRBModel model = GRBModel(env);

	// initizalization of the variables for the model
	vector<set<int>> C_indices;//transforming to vector<set>. why? to retrieve the subtour with its index (if this were set<set<int>>, there is no ordering in the storage) 
	C_indices.insert(C_indices.end(), inst.uniqueSubsets.begin(), inst.uniqueSubsets.end());
	// declaration of the variables for the model
	vector<vector<GRBVar>> x;
	x.resize(C_indices.size(), vector<GRBVar>(inst.C_of_k.size()));
	vector<pair<int, int>> valid_indices;  // to store indices of created variables
	for (int k = 0; k < inst.C_of_k.size(); k++) {
		for (int s = 0; s < C_indices.size(); s++) {
			if (inst.C_of_k[k].find(C_indices[s]) != inst.C_of_k[k].end()) {
				x[s][k] = model.addVar(0, 1, 0, GRB_BINARY);
				valid_indices.emplace_back(s, k);  // store index
			}
		}
	}
	model.update();


	//set objective
	GRBLinExpr obj;
	// Create the objective function
	for (const auto& p : valid_indices) {
		int s = p.first;
		int k = p.second;
		obj += x[s][k];  // this works only if x[s][k] was initialized!
	}
	model.update();
	model.setObjective(obj, GRB_MAXIMIZE);

	//Add the cover constraints
	vector<GRBLinExpr> coverConstraint(inst.C_of_k.size(), 0);
	for (int k = 0; k < inst.C_of_k.size(); k++) {
		for (set<int> subset: inst.U_of_k[k]) {
			int s = isSubsetInASol(subset, C_indices);
			if (s == -1) continue;
			for (int kp=0; kp <= k; kp++) {
				//inner sum: check whether this s is included in k<=kp
				if (inst.C_of_k[kp].find(subset) != inst.C_of_k[kp].end()) {
					coverConstraint[k] += x[s][kp];
				}
			}
			
		}
	}

	for (int k = 0; k < inst.C_of_k.size(); k++) {
		if(coverConstraint[k].size()==0) continue;
		model.addConstr(coverConstraint[k] >= 1);
	}
	//Other constaint.

	for (int k = 0; k < inst.C_of_k.size(); k++) {
		for (set<int> subset : inst.U_of_k[k]) {
			int s = isSubsetInASol(subset, C_indices);
			if (s == -1) continue;
			for (set<int> subset_sp : inst.C_of_k[k]) {
				int sp = isSubsetInASol(subset_sp, C_indices);
				GRBLinExpr sum_inside;
				for (int kp = 0; kp < k; kp++) {
					if (inst.C_of_k[kp].find(subset) != inst.C_of_k[kp].end()) {
						sum_inside += x[s][kp];
					}
				}
				model.addConstr(x[sp][k] <= 1 - sum_inside);
			}
		}
	}

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	//model.getEnv().set(GRB_DoubleParam_TimeLimit, max(inst.timeLimit - (getCPUTime() - inst.startTimeSaved), EPSILON)); //use the remaining time left

	//model.set(GRB_IntParam_PoolSearchMode, 2); // Find all solutions in the pool
	//model.set(GRB_IntParam_PoolSolutions, 100); // Max number of alternative solutions to store
	//model.set(GRB_DoubleParam_PoolGap, 0.0); // Only include optimal solutions (no deviation allowed)

	model.update();
	// optimize	
	model.optimize();

	int numSolutions = model.get(GRB_IntAttr_SolCount);
	double optimalValue = model.get(GRB_DoubleAttr_ObjVal);
	sol.optimalObjVal = optimalValue;
	sol.numSol = numSolutions;

	//if optimality is proven, extract the edges and the objective value
	vector<vector<int>> solution;
	if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
		cout << "Optimal found" << endl;
		sol.opt = 1;
	}
	else if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
		cout << "Time limit!" << endl; return;
	}
	else {
		cout << "Unknown error!" << endl; return;
	}

	cout << "Found " << numSolutions << " solutions at optimal value: " << optimalValue << "\n";

	for (int i = 0; i < numSolutions; ++i) {
		model.set(GRB_IntParam_SolutionNumber, i); // select solution i
		
		cout << "Solution " << i << ": " << endl;

		for (const auto& p : valid_indices) {
			int s = p.first;
			int k = p.second;
			double val = x[s][k].get(GRB_DoubleAttr_Xn); // get value in solution i
			if (val > 0.5) {
				for (int i : C_indices[s]) cout << i+1 << "(" << k +1 << ")" << ", ";
				cout << endl;
				//for (int i : C_indices[s]) cout << "(" << inst.xCoord[i] << ", " << inst.yCoord[i] << ") -- ";
				//cout << endl;
			}
		}
		cout << "-------------------------------------------------------------------------------------------------------\n";
		break;
	}


}

