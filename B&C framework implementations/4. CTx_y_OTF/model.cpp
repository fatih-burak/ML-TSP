#include "model.h"
constexpr int CT = 2;
constexpr double PERC_CUT = 0.5;
constexpr bool IS_SHORTEST = false;

void ILP(Instance& inst, Solution& sol) {
	// create a model
	GRBEnv env = GRBEnv();
	env.set(GRB_IntParam_OutputFlag, 1);
	GRBModel model = GRBModel(env);

	// declaration of the variables for the model
	vector<vector<GRBVar> > x;
	x.resize(inst.nbCity, vector<GRBVar>(inst.nbCity));

	// initizalization of the variables for the model
	for (int j = 0; j < inst.nbCity; j++) {
		for (int i = 0; i < j; i++) {	
			x[i][j] = model.addVar(0, 1, 0, GRB_BINARY);
		}
	}
	model.update();

	//set objective
	GRBLinExpr obj;
	// Create the objective function
	for (int j = 0; j < inst.nbCity; j++) {
		for (int i = 0; i < j; i++) {
			obj += inst.dist[i][j] * x[i][j];
		}
	}
	model.update();
	model.setObjective(obj, GRB_MINIMIZE);

	//Add the degree-2 constraints
	vector<GRBLinExpr> degreeConstraint(inst.nbCity, 0);
	for (int v = 0; v < inst.nbCity; v++) {
		for (int j = 0; j < inst.nbCity; j++) {
			for (int i = 0; i < j; i++) {
				if (v==i || v==j)  degreeConstraint[v] += x[i][j];
			}
		}
	}
	for (int v = 0; v < inst.nbCity; v++) model.addConstr(degreeConstraint[v] == 2); 
	
	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_DoubleParam_TimeLimit, max(inst.timeLimit - (getCPUTime() - inst.startTimeSaved), EPSILON)); //use the remaining time left
	model.getEnv().set(GRB_IntParam_Seed, 0); // Fixed seed
	model.set(GRB_IntParam_LazyConstraints, 1); 		// indicate that we want to add Lazy Constraints
	cout << "Current seed value: " << model.getEnv().get(GRB_IntParam_Seed) << std::endl;
	model.update();
	// optimize	
	SolutionCollector cb(x, inst, sol); // pass variable reference
	model.setCallback(&cb);  // set callback

	//Optimize
	model.optimize();
	

	//if optimality is proven, extract the edges and the objective value
	if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
		cout << "Optimality proven!" << endl;
		sol.opt = 1;
		sol.optimalObjVal = ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
		sol.LB = ceil(model.get(GRB_DoubleAttr_ObjVal) - EPSILON);
	}
	else if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
		cout << "Time limit!" << endl;
		sol.LB = ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON);
	}
	else {
		cout << "Program should not enter here..." << endl;
	}
	return;
}

void SolutionCollector::callback() {
	if (where == GRB_CB_MIPSOL) { // A new feasible solution found
		int old_no = sol.cuts.size();
		sol.iter += 1; //Increase the number of iterations
		cout << "Iteration = " << sol.iter << " and total nb of cuts added = " << sol.cuts.size() << endl;

		// Create an intermediate solution object
		InterSol inter_sol;
		inter_sol.nbCity = inst.nbCity;

		// Extract solution
		vector<pair<int, int>> edges;
		for (int j = 0; j < x.size(); j++) {
			for (int i = 0; i < j; i++) {
				double val = getSolution(x[i][j]);
				if (val > 0.95) {
					inter_sol.sol_edges.push_back({ i,j });
				}
			}
		}
		inter_sol.objVal = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
		cout << "Incumbent with obj val = " << getDoubleInfo(GRB_CB_MIPSOL_OBJ) << " is found!" << endl;

		// SETTINGS
		int cutType = CT; //cutType: 1->SEC Type 3, 2->SEC Type 5, 3->Both
		inter_sol.perc = PERC_CUT; //include all size SECs.
		inter_sol.isShortest = IS_SHORTEST; // true: shortest cut is added
		
		//Check for cycles
		vector<vector<int>> new_cuts = findCuts(inter_sol);

		if (new_cuts.empty()) {
			//do nothing, this solution might be optimal if solver cannot find any other better solution. 
			cout << "A solution without cycles is found!" << endl;
			return;
		}
		// Add the cycles in the new_cuts to sol.cuts list
		sol.cuts.insert(sol.cuts.end(), new_cuts.begin(), new_cuts.end());

		// Add the cuts to the model
		for (const vector<int>& cycle : new_cuts) {
			if (cycle.size() < 3) continue; //this will not happen in symmetric TSP as edges are only defined for i<j

			if (cutType == 1 || (cutType == 3 && cycle.size() <= (2 * inst.nbCity + 1) / 3)) {
				// SEC Type-1
				GRBLinExpr temp_SEC;
				for (int i : cycle) {
					for (int j : cycle) {
						if (i < j) temp_SEC += x[i][j];
					}
				}
				addLazy(temp_SEC <= cycle.size() - 1);
			}
			else if (cutType == 2 || cutType == 3) {
				// SEC Type-2
				GRBLinExpr temp_SEC;
				vector<int> not_in_cycle = setDifference(inst.V, cycle);
				for (int i : cycle) {
					for (int j : not_in_cycle) {
						temp_SEC += x[min(i, j)][max(i, j)];
					}
				}
				addLazy(temp_SEC >= 2);
			}
		}

		cout << "Nb new cuts = " << sol.cuts.size() - old_no << " same as " << new_cuts.size() << endl;
		cout << "_________________________________________________________________________________" << endl;
	}
}

vector<vector<int>> findCuts(InterSol& iSol) {
	/// <summary>
	/// Detects subtours (that are smaller than nbCity*perc ) in a TSP solution and returns cycles that violate the full tour constraint.
	/// </summary>
	/// <param name="iSol">Intermediate solution object: including sol_edges, nbCity, perc, objVal and shortest_cycle.</param>
	/// <returns>A list of cycles (each a vector of vertex indices) that can be used as SECs.</returns>
	vector<vector<int>> new_cuts; 
	
	bool cycleDetected = true; //if no cycle is detected, it will be switched to false.

	while (!iSol.sol_edges.empty()) {
		pair<int, int> edge = iSol.sol_edges[0];
		int first = edge.first;				//take the first vertex
		vector<int> cycle = { first };		//add it as the beginning of the cycle
		int second = edge.second;			//move along to that edge and find the second vertex

		while (cycle.size() < 2 || cycle.front() != cycle.back()) {
			cycle.push_back(second);		//add the vertex in the cycle
			first = second;					//and make it as the current first vertex

			// Remove the used edge
			iSol.sol_edges.erase(remove(iSol.sol_edges.begin(), iSol.sol_edges.end(), edge), iSol.sol_edges.end());

			for (const auto& candid : iSol.sol_edges) {
				if (candid.first == first || candid.second == first) {	//find another edge including the current first vertex
					edge = candid;
					second = (edge.first == first) ? edge.second : edge.first; //find the second vertex
					break;	//as soon as you find, exit
				}
			}
		}
		
		// Check the cycle
		if (cycle.size() - 1 < iSol.nbCity) { // if equal, that means opt
			if (iSol.shortest_cycle.empty() || cycle.size() - 1 < iSol.shortest_cycle.size()) {
				iSol.shortest_cycle = vector<int>(cycle.begin(), cycle.end() - 1);
			}
			if (!iSol.isShortest) {
				if (cycle.size() - 1 <= iSol.nbCity * iSol.perc) {
					// Add the cycle without the last node (repeated start)
					new_cuts.push_back(vector<int>(cycle.begin(), cycle.end() - 1));
				}
				else {
					cout << "Cut is skipped because it is of size = " << (cycle.size() - 1) << endl;
				}
			}
		}
		else cycleDetected = false;
	}
	if (iSol.isShortest && cycleDetected) new_cuts.push_back(iSol.shortest_cycle); // Add the shortest cycle only
	/*
	if (new_cuts.empty() && !iSol.shortest_cycle.empty()) { //if nothing is added, add the smallest (can only happen if perc is less than 0.5)
		new_cuts.push_back(iSol.shortest_cycle);
	}
	*/
	return new_cuts;
}

void TSP_DLYD(Instance& inst, Solution& sol) {
	//INFO
	if (IS_SHORTEST) cout << "CT" << CT << "_" << "SM" << "SEC" << endl;
	else cout << "CT" << CT << "_" << (int)PERC_CUT * 100 << "SEC" << endl;
	// Call the ILP model
	ILP(inst, sol); 
}

