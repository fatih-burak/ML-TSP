#include "model.h"
#include <random>
constexpr double PERC_PRECUTS = 1.0; //PERCENTAGE OF ACTUAL TRUE CUTS
constexpr double PRECUT_MULTIPLIER = 0.0;
constexpr int NB_SOL_STORE = GRB_MAXINT;
constexpr int CT = 3;
constexpr double PERC_CUT = 1.0;
constexpr bool IS_SHORTEST = false;

set<set<int>> generate_RNDM_SECs(int number_secs, Instance& inst) {
	cout << "Generate " << number_secs << " random cuts " << endl;

	// Tranform vector<vector<int>> to set<set<int>>
	std::set<std::set<int>> C_star_Grand;
	for (const auto& vec : inst.precuts) {
		std::set<int> inner(vec.begin(), vec.end());
		C_star_Grand.insert(inner);
	}

	set<set<int>> grand_list_SECs;
	int number_vertices = inst.nbCity;
	for (int i = 0; i < number_secs; i++) {

		while (1) {
			// Random number generator
			random_device rd; mt19937 gen(rd());
			// Pick random size between 3 and number_vertices/2
			uniform_int_distribution<> size_dist(3, number_vertices / 2);
			int size_SEC = size_dist(gen);

			// Create list of all vertices [0, ..., number_vertices-1]
			vector<int> vertices(number_vertices);
			for (int v = 0; v < number_vertices; ++v) vertices[v] = v;

			// Shuffle and pick first size_SEC vertices
			shuffle(vertices.begin(), vertices.end(), gen);

			std::set<int> one_SEC(vertices.begin(), vertices.begin() + size_SEC);
			
			// Check if one_SEC already exists in inst.precuts
			bool already_exists = !C_star_Grand.insert(one_SEC).second;

			if (!already_exists && one_SEC.size() <= number_vertices / 2) { //if it is not included in the actual cuts list (C*) AND size is n/2 (which is automatic)
				grand_list_SECs.insert(one_SEC);
				cout << "Cycle No= " << i+1 << ": ";
				for (auto v : one_SEC) cout << v << " ";
				cout << "\n";
				break;
			}
			else if (one_SEC.size() > number_vertices / 2) {
				cout << "size greater than n/2..." << endl;
			}
			else if (already_exists) {
				cout << "exists..." << endl;
			}
			else {
				cout << "unknown problem" << endl;
			}
		}
	
	
	}
	return grand_list_SECs;
}

vector<pair<int, int>> buildMST(const vector<vector<double>>& dist) {
	// -----------------------------
	// Function: buildMST
	// Builds a Minimum Spanning Tree using Prim's algorithm
	// Returns the edges of the MST as (u,v) pairs
	// -----------------------------
	int n = dist.size();
	vector<double> minEdge(n, 1e18);    // Track the minimum edge cost to connect each node to the MST
	vector<int> parent(n, -1);  // Track parent of each node in the MST
	vector<bool> used(n, false);    // Track which nodes are already included in MST

	minEdge[0] = 0.0;     // Start w/ node 0 (arbitrary choice)
	vector<pair<int, int>> mstEdges;   // Store edges of MST

	for (int i = 0; i < n; i++) {     // Repeat n times to add n nodes
		int v = -1;
		// Pick the unused node with smallest edge cost
		for (int j = 0; j < n; j++) {
			if (!used[j] && (v == -1 || minEdge[j] < minEdge[v])) {
				v = j;
			}
		}

		used[v] = true;  // Mark this node as used (added to MST)
		if (parent[v] != -1) {
			mstEdges.push_back({ v, parent[v] }); // If it has a parent, add that edge to MST
			//cout << v << " " << parent[v] << endl; 
		}
		for (int u = 0; u < n; u++) {   // Update minEdge[] values for neighbors of v
			if (dist[v][u] < minEdge[u]) {
				minEdge[u] = dist[v][u];
				parent[u] = v;
			}
		}
	}
	return mstEdges;
}

bool dfsFindPath(int u, int target, vector<int>& path, vector<bool>& visited, const vector<vector<int>>& adj) {
	// -----------------------------
	// Function: dfsFindPath
	// Finds path between two nodes in MST using DFS
	// Stores path in "path" vector
	// -----------------------------
	if (u == target) { // If the tarfet is reached, push and return true
		path.push_back(u);
		return true;
	}
	visited[u] = true;
	for (int v : adj[u]) {  // Try all neighbors
		if (!visited[v]) {
			// If path found through v, push current node and return
			if (dfsFindPath(v, target, path, visited, adj)) {
				path.push_back(u);
				return true;
			}
		}
	}
	return false; // No path found (should not be possible since we work on complete graph. 
}

set<set<int>> CheckAllMSTCycles(Instance& inst) {

	std::ofstream file(inst.name.c_str(), std::ios::out | std::ios::app);
	
	set<set<int>> grand_list_SECs;
	srand(time(NULL)); // Seed for random edge selection

	int n = inst.nbCity;

	// Step 1: Build MST
	vector<pair<int, int>> mstEdges = buildMST(inst.dist);

	// Build adjacency list for MST
	vector<vector<int>> adj(n);
	set<pair<int, int>> mstSet; // For checking whether an edge is in MST
	for (auto e : mstEdges) {
		adj[e.first].push_back(e.second);
		adj[e.second].push_back(e.first);
		mstSet.insert({ min(e.first,e.second), max(e.first,e.second) });
	}

	// -----------------------------
	// Collect non-MST edges
	// -----------------------------
	vector<tuple<double, int, int>> nonMSTEdges;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (inst.dist[i][j] < 1e9) { // valid edge
				pair<int, int> e = { i,j };
				if (mstSet.find(e) == mstSet.end()) {
					nonMSTEdges.push_back(make_tuple(inst.dist[i][j], i, j));
				}
			}
		}
	}

	// Sort by edge weight ascending
	std::sort(nonMSTEdges.begin(), nonMSTEdges.end());


	int k = 0;
	int totalCycles = 0;
	//cout << "using " << nonMSTEdges.size() << " many non-MST edges." << endl;
	while (totalCycles <nonMSTEdges.size()) {
		if (k >= nonMSTEdges.size()) break; // there is no edge left to add to the MST tree
		while (1) { // GENERATE UNTIL YOU FIND ONE THAT WAS NOT INCLUDED IN PRECUTS
			if (k >= nonMSTEdges.size()) break; // There is no edge left to add to the MST tree
			//cout << "Trying nonMST edge " << k + 1 << "/" << nonMSTEdges.size() << "| ";
			set<int> one_SEC;
			double w;
			int a, b;
			tie(w, a, b) = nonMSTEdges[k];
			//cout << "Extra edge added: (" << a << "," << b << ") weight=" << w << "\n";

			// Step 3: find cycle path in MST between a and b
			vector<int> path;
			vector<bool> visited(n, false);
			dfsFindPath(a, b, path, visited, adj);

			// Add b to complete cycle
			path.push_back(b);
			reverse(path.begin(), path.end());

			int cycleLen = 0;
			// Print cycle vertices
			//cout << "Cycle: ";
			for (int v = 0; v < path.size(); v++) {
				if (v < path.size()-1 && path.size()-1<=n/2){
					int i = path[v];
					int j = path[v + 1];
					cycleLen += inst.dist[i][j];
					//cout << i << "-" << j << ": " << inst.dist[i][j] << endl;
				}
				one_SEC.insert(path[v]);
			}
			k++;
			if (one_SEC.size() > n / 2) {
				//cout << "bigger than n/2" << endl;
				break;
			}
			grand_list_SECs.insert(one_SEC);
			totalCycles++;
			//cout << "Cycle No= " << totalCycles << ": "; 
			//cout << "(" << one_SEC.size() << ") ";
			//for (auto v : one_SEC) cout << v << " ";
			//cout << "\n";
			cout << totalCycles << "\t" << one_SEC.size() << "\t" << cycleLen << endl;
			//file << totalCycles << "\t" << one_SEC.size() << "\t" << cycleLen << endl;
			
		}
	}

	file.close();
	return grand_list_SECs;
}

set<set<int>> Find_K_MST_CYCLES(int K, Instance& inst) {

	// Tranform vector<vector<int>> to set<set<int>>
	std::set<std::set<int>> C_star_Grand;
	for (const auto& vec : inst.precuts) {
		std::set<int> inner(vec.begin(), vec.end());
		C_star_Grand.insert(inner);
	}

	set<set<int>> grand_list_SECs;
	srand(time(NULL)); // Seed for random edge selection

	int n = inst.nbCity;

	// Step 1: Build MST
	vector<pair<int, int>> mstEdges = buildMST(inst.dist);

	// Build adjacency list for MST
	vector<vector<int>> adj(n);
	set<pair<int, int>> mstSet; // For checking whether an edge is in MST
	for (auto e : mstEdges) {
		adj[e.first].push_back(e.second);
		adj[e.second].push_back(e.first);
		mstSet.insert({ min(e.first,e.second), max(e.first,e.second) });
	}

	// -----------------------------
	// Collect non-MST edges
	// -----------------------------
	vector<tuple<double, int, int>> nonMSTEdges;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (inst.dist[i][j] < 1e9) { // valid edge
				pair<int, int> e = { i,j };
				if (mstSet.find(e) == mstSet.end()) {
					nonMSTEdges.push_back(make_tuple(inst.dist[i][j], i, j));
				}
			}
		}
	}

	// Sort by edge weight ascending
	std::sort(nonMSTEdges.begin(), nonMSTEdges.end());

	cout << "Generate " << K << " MST cuts ";
 	// repeat process K times, but ordered by edge weight
	//for (int k = 0; k < K && k < (int)nonMSTEdges.size(); k++) {
	int k = 0;
	int totalCycles = 0;
	cout << "using " << nonMSTEdges.size() << " many non-MST edges." << endl;
	while(totalCycles < min(K, (int)nonMSTEdges.size())){
		if (k >= nonMSTEdges.size()) break; // there is no edge left to add to the MST tree
		while (1) { // GENERATE UNTIL YOU FIND ONE THAT WAS NOT INCLUDED IN PRECUTS
			if (k >=nonMSTEdges.size()) break; // There is no edge left to add to the MST tree
			cout << "Trying nonMST edge " << k+1 << "/" << nonMSTEdges.size() << "| ";
			set<int> one_SEC;
			double w;
			int a, b;
			tie(w, a, b) = nonMSTEdges[k];
			//cout << "Extra edge added: (" << a << "," << b << ") weight=" << w << "\n";

			// Step 3: find cycle path in MST between a and b
			vector<int> path;
			vector<bool> visited(n, false);
			dfsFindPath(a, b, path, visited, adj);

			// Add b to complete cycle
			path.push_back(b);
			reverse(path.begin(), path.end());

			// Print cycle vertices
			//cout << "Cycle: ";
			for (int v : path) 	one_SEC.insert(v);
			k++;
			bool already_exists = !C_star_Grand.insert(one_SEC).second;

			if (!already_exists && one_SEC.size() <= n / 2) { //if it is not included in the actual cuts list (C*) AND size is n/2
				grand_list_SECs.insert(one_SEC);
				totalCycles++;
				cout << "Cycle No= " << totalCycles << ": ";
				for (auto v : one_SEC) cout << v << " ";
				cout << "\n";
				break;
			}
			else if (one_SEC.size() > n / 2) {
				cout << "size greater than n/2..." << endl;
			}
			else if (already_exists) {
				cout << "exists..." << endl;
			}
			else cout << "unknown problem" << endl;
		}
	}
	return grand_list_SECs;
}


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

	int all_except_last = (iSol.MULTI_sol_edges.size() >= 2) ? iSol.MULTI_sol_edges.size() - 1 : iSol.MULTI_sol_edges.size();

	for (int sol_no = 0; sol_no < all_except_last; sol_no++) {
		double objValSol = iSol.MULTI_objVal[sol_no];
		cout << "Included SOL- " << sol_no << " with objective = " << objValSol << endl;
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
						saveCycle(foundCycle, objValSol, iSol.cutSaveLoc);
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
	if (IS_SHORTEST) cout << "CT" << CT << "_" << "SM" << "SEC_ASEF" << endl;
	else cout << "CT" << CT << "_" << (int)PERC_CUT * 100 << "SEC_ASEF" << endl;
	cout << PERC_PRECUTS * 100 << "% TRUE/ " << PRECUT_MULTIPLIER * 100 << "% FALSE " << endl;

	//FIND THE SIZE OF PRECUTS (FULL SIZE, READ FROM A TXT FILE WHICH INCLUDES SECS FOUND USING CT1_SM_OS
	int size_saved_cuts = inst.precuts.size();
	int size_selected = round(PERC_PRECUTS * size_saved_cuts);
	cout << "True cuts need to be added = " << size_selected << "| ";
	// Copy precuts so we can shuffle
	vector<vector<int>> shuffled_precuts = inst.precuts;
	// Shuffle using a random generator
	random_device rd; mt19937 g(rd());
	shuffle(shuffled_precuts.begin(), shuffled_precuts.end(), g);
	// Take only the first size_selected cuts
	for (int i = 0; i < size_selected && i < shuffled_precuts.size(); i++) {
		set<int> cutSet(shuffled_precuts[i].begin(), shuffled_precuts[i].end()); // convert vector<int> to set<int>
		sol.cuts.insert(cutSet);                     // only inserts if it's not a duplicate
	}
	int b4_cuts = sol.cuts.size();
	cout << "cuts added = " << sol.cuts.size() << endl;

	//GENERATE RANDOMLY size_selected*PRECUT_MULTIPLIER MANY CUTS
	set<set<int>> randomlyGenSECs = generate_RNDM_SECs(round(size_selected * PRECUT_MULTIPLIER), inst);
	for (set<int> random_SEC : randomlyGenSECs) sol.cuts.insert(random_SEC);   
	int after_cuts = sol.cuts.size();
	cout << "Random cuts added = " << after_cuts - b4_cuts << endl;

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

