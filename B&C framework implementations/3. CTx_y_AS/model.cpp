#include "model.h"
constexpr int NB_SOL_STORE = GRB_MAXINT;
constexpr int CT = 2;
constexpr double PERC_CUT = 0.5;
constexpr bool IS_SHORTEST = false;


InterSol ILP_ALLSOL(Instance& inst, Solution& sol, int cutType) {
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
				if (v == i || v == j)  degreeConstraint[v] += x[i][j];
			}
		}
	}
	for (int v = 0; v < inst.nbCity; v++) model.addConstr(degreeConstraint[v] == 2);

	for (const set<int>& cycle : sol.cuts) {
		if (cycle.size() < 3) continue; //this will not happen in symmetric TSP as edges are only defined for i<j

		if (cutType == 1 || (cutType == 3 && cycle.size() <= (2 * inst.nbCity + 1) / 3)) {
			// SEC Type-1
			GRBLinExpr temp_SEC;
			for (int i : cycle) {
				for (int j : cycle) {
					if (i < j) temp_SEC += x[i][j];
				}
			}
			model.addConstr(temp_SEC <= cycle.size() - 1);
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
			model.addConstr(temp_SEC >= 2);
		}
	}

	// change some settings
	model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
	model.getEnv().set(GRB_IntParam_Threads, 1);
	model.getEnv().set(GRB_DoubleParam_TimeLimit, max(inst.timeLimit - (getCPUTime() - inst.startTimeSaved), EPSILON)); //use the remaining time left
	model.getEnv().set(GRB_IntParam_Seed, 0); // Fixed seed
	cout << "Current seed value: " << model.getEnv().get(GRB_IntParam_Seed) << std::endl;
	// Set PoolSearchMode = 0 (default) — Gurobi may find additional solutions incidentally
	model.set(GRB_IntParam_PoolSearchMode, 0);
	// Set the PoolSolutions parameter to the maximum allowed value
	model.set(GRB_IntParam_PoolSolutions, inst.howManySolutions);

	model.update();

	//Optimize
	model.optimize();
	// Number of solutions found
	int solCount = model.get(GRB_IntAttr_SolCount);
	cout << "Number of solutions found: " << solCount << endl;

	//create an intermediate solution object
	InterSol inter_sol;
	inter_sol.nbCity = inst.nbCity;

	//if optimality is proven, extract the edges and the objective value
	if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
		// Loop over all available solutions in the pool
		for (int poolSol = 0; poolSol < solCount; ++poolSol) {
			model.set(GRB_IntParam_SolutionNumber, poolSol);
			double objVal = model.get(GRB_DoubleAttr_PoolObjVal);
			inter_sol.MULTI_objVal.push_back(objVal);
			cout << "Solution " << poolSol << ", Obj = " << objVal << endl;

			vector<pair<int, int>> sol_edges;
			for (int j = 0; j < inst.nbCity; j++) {
				for (int i = 0; i < j; i++) {
					if (x[i][j].get(GRB_DoubleAttr_Xn) > 0.95) {
						sol_edges.push_back({ i,j });
					}
				}
			}

			inter_sol.MULTI_sol_edges.push_back(sol_edges);

		}

	}
	else if (model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
		inter_sol.isTimeLimitHit = true;
		cout << "Time limit!" << endl;
		inter_sol.MULTI_objVal.push_back(ceil(model.get(GRB_DoubleAttr_ObjBound) - EPSILON));

	}
	else {
		cout << "Program should not enter here..." << endl;
	}

	sol.iter += 1; //Increase the number of iterations
	return inter_sol;
}

void saveCycle(const vector<int>& cycle, const double& objVal, const string& filename) {
	// Save the full cycle to the file
	std::ofstream outfile(filename, std::ios::out | std::ios::app);
	if (outfile.is_open()) {
		outfile << objVal << " ";
		for (int val : cycle) {
			outfile << val << " ";
		}
		outfile << endl;
		outfile.close();
	}
	else cerr << "Error: Could not open file " << filename << " for writing.\n";
}

vector<vector<int>> findCuts(InterSol& iSol) {
	/// <summary>
	/// Detects subtours (that are smaller than nbCity*perc ) in a TSP solution and returns cycles that violate the full tour constraint.
	/// </summary>
	/// <param name="iSol">Intermediate solution object: including sol_edges, nbCity, perc, objVal and shortest_cycle.</param>
	/// <returns>A list of cycles (each a vector of vertex indices) that can be used as SECs.</returns>
	vector<vector<int>> new_cuts; 
	bool optimalFound = false; //if no cycle is detected in the optimal solution (out of all the solutions found only checking the optimal one)

	for (int sol_no = 0; sol_no < iSol.MULTI_sol_edges.size(); sol_no++) {
		double objValSol = iSol.MULTI_objVal[sol_no];
		vector<pair<int, int>> a_sol_edges = iSol.MULTI_sol_edges[sol_no];
		vector<int> shortest_cycle_in_this_sol;

		bool cycleDetected = true; //if no cycle is detected, it will be switched to false.

		while (!a_sol_edges.empty()) {
			pair<int, int> edge = a_sol_edges[0];
			int first = edge.first;				//take the first vertex
			vector<int> cycle = { first };		//add it as the beginning of the cycle
			int second = edge.second;			//move along to that edge and find the second vertex

			while (cycle.size() < 2 || cycle.front() != cycle.back()) {
				cycle.push_back(second);		//add the vertex in the cycle
				first = second;					//and make it as the current first vertex

				// Remove the used edge
				a_sol_edges.erase(remove(a_sol_edges.begin(), a_sol_edges.end(), edge), a_sol_edges.end());

				for (const auto& candid : a_sol_edges) {
					if (candid.first == first || candid.second == first) {	//find another edge including the current first vertex
						edge = candid;
						second = (edge.first == first) ? edge.second : edge.first; //find the second vertex
						break;	//as soon as you find, exit
					}
				}
			}

			// Check the cycle
			if (cycle.size() - 1 < iSol.nbCity) { // if equal, that means opt
				if (shortest_cycle_in_this_sol.empty() || cycle.size() - 1 < shortest_cycle_in_this_sol.size()) {
					shortest_cycle_in_this_sol = vector<int>(cycle.begin(), cycle.end() - 1);
				}
				if (!iSol.isShortest) {
					if (cycle.size() - 1 <= iSol.nbCity * iSol.perc) {
						// Add the cycle without the last node (repeated start)
						vector<int> foundCycle = vector<int>(cycle.begin(), cycle.end() - 1);
						new_cuts.push_back(foundCycle);
						saveCycle(foundCycle, objValSol,  iSol.cutSaveLoc);
					}
					else {
						cout << "Cut is skipped because it is of size = " << (cycle.size() - 1) << endl;
					}
				}
			}
			else {
				cycleDetected = false;
				if (iSol.MULTI_objVal[sol_no] <= *min_element(iSol.MULTI_objVal.begin(), iSol.MULTI_objVal.end())) {
					// basically checking whether this solution is the best found one (if that is the case, and it has no cycle: then we should be in the optimal case)
					optimalFound = true;
					break;
				}
			}
		}
		if (optimalFound) break;
		if (iSol.isShortest && cycleDetected) iSol.shortest_cycle.push_back(shortest_cycle_in_this_sol); //add the shortest cycle to the list of shortest cycles.
		/*
		if (new_cuts.empty() && !iSol.shortest_cycle.empty()) { //if nothing is added, add the smallest (can only happen if perc is less than 0.5)
			new_cuts.push_back(iSol.shortest_cycle);
		}
		*/
	}
	if (iSol.isShortest) new_cuts = iSol.shortest_cycle;
	if (optimalFound) new_cuts.clear(); // clears the list 

	return new_cuts;
}

void TSP_DLYD(Instance& inst, Solution& sol) {
	//INFO
	if (IS_SHORTEST) cout << "CT" << CT << "_" << "SM" << "SEC" << endl;
	else cout << "CT" << CT << "_" << (int)PERC_CUT * 100 << "SEC" << endl;
	// SET VALUES AND START
	inst.howManySolutions = NB_SOL_STORE;
	while (true) {
		cout << "Iteration = " << sol.iter + 1 << " and total nb of cuts added = " << sol.cuts.size() << endl;
		int old_no = sol.cuts.size();
		// Call the ILP model
		InterSol iSol = ILP_ALLSOL(inst, sol, CT); //last number indicates, cutType: 1->SEC Type 3, 2->SEC Type 5, 3->Both
		if (iSol.MULTI_objVal[0] > EPSILON) sol.LB = iSol.MULTI_objVal[0];
		cout << "Current LB = " << sol.LB << endl;

		//If the time limit is reached, exit.
		if (iSol.isTimeLimitHit) break;

		//Check for cycles
		iSol.perc = PERC_CUT; //include all size SECs.
		iSol.isShortest = IS_SHORTEST; // true: shortest cut is added
		iSol.cutSaveLoc = inst.cutFolder + inst.name + ".txt"; // Directory for the saved cuts
		vector<vector<int>> new_cuts = findCuts(iSol);
		cout << "Extracted nb of cuts = " << new_cuts.size() << endl;
		//for (auto content : iSol.shortest_cycle) cout << content << endl; //IF YOU WANT TO PRINT THE CONTENT OF THE CUT

		if (new_cuts.empty()) { 
			cout << "Optimal solution found!" << endl;
			sol.optimalObjVal = iSol.MULTI_objVal[0];
			sol.opt = 1; break; 
		}
		//Add the cycles in the new_cuts to sol.cuts list
		for (const auto& cut : new_cuts) {
			set<int> cutSet(cut.begin(), cut.end()); // convert vector<int> to set<int>
			sol.cuts.insert(cutSet);                     // only inserts if it's not a duplicate
		}

		cout << "Nb eliminated duplicates = " << new_cuts.size() - (sol.cuts.size() - old_no) << endl;
		sol.eliminatedSECs += new_cuts.size() - (sol.cuts.size() - old_no);

		cout << "Nb new cuts = " << sol.cuts.size() - old_no << endl;
		cout << "___________________________________________________________________" << endl;
	}
}

