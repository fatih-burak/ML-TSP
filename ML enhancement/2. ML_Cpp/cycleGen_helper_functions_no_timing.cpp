#include "cycleGenHelper_functions.h"
#include <climits>
static string shellQuote(const string& s) {
    // Wraps a string in single quotes so PATHs with spaces or quotes can be passed in command strings. We use it to call the pyhton file
    string out = "'";
    for (char c : s) {
        if (c == '\'') out += "'\\''";
        else out += c;
    }
    out += "'";
    return out;
}

static vector<double> readPredictionsFromFile(const string& predPath) {
    vector<double> preds;
    ifstream in(predPath);
    if (!in.is_open()) {
        throw runtime_error("Could not open prediction file: " + predPath);
    }

    string line;
    while (getline(in, line)) {
        if (line.empty()) continue;
        preds.push_back(stod(line));
    }
    return preds;
}
struct FeatureScratch {
    vector<unsigned char> inS;   // n
    vector<int> outside;         // up to n
    vector<double> minEdge;      // up to n/2
    vector<unsigned char> used;  // up to n/2

    //visited marker for greedy walk (global vertex ids)
    vector<int> visited;
    int visited_stamp = 1;

    explicit FeatureScratch(int n = 0) { reset(n); }

    void reset(int n) {
        inS.assign(n, 0);
        outside.clear();
        minEdge.clear();
        used.clear();
        visited.assign(n, 0);
        visited_stamp = 1;
    }

    void nextVisitedStamp() {
        ++visited_stamp;
        if (visited_stamp == INT_MAX) {
            fill(visited.begin(), visited.end(), 0);
            visited_stamp = 1;
        }
    }
    void markVisited(int v) { visited[v] = visited_stamp; }
    bool isVisited(int v) const { return visited[v] == visited_stamp; }
};


// ================== MARKSET (timestamp membership) ==================
struct MarkSet {
    vector<int> mark;  // size n, init 0
    int stamp = 1;

    explicit MarkSet(int n = 0) : mark(n, 0) {}

    void reset(int n) {
        mark.assign(n, 0);
        stamp = 1;
    }

    void next() {
        ++stamp;
        if (stamp == INT_MAX) {
            fill(mark.begin(), mark.end(), 0);
            stamp = 1;
        }
    }

    void add(int v) { mark[v] = stamp; }
    bool has(int v) const { return mark[v] == stamp; }
};

// ================== Top-K neighbor lists ==================
struct TopKLists {
    int n = 0;
    int K = 0;
    vector<vector<int>> distNN; // K nearest by dist (ascending)
    vector<vector<int>> rcNN;   // K smallest rc (ascending)
    vector<vector<int>> lpNN;   // K largest lp (descending)
};

// assumes x/rc stored upper-tri only, symmetric
static double getXfull(const vector<vector<double>>& x, int i, int j) {
    return (i < j) ? x[i][j] : x[j][i];
}
static double getRCfull(const vector<vector<double>>& rc, int i, int j) {
    return (i < j) ? rc[i][j] : rc[j][i];
}

// Maintain Top-K (fast for small K)
template <class Better>
static void pushTopK(vector<int>& idx, vector<double>& val,
    int j, double v, int K, Better better)
{
    const int sz = (int)idx.size();
    if (sz < K) {
        idx.push_back(j);
        val.push_back(v);
        int p = sz;
        while (p > 0 && better(val[p], val[p - 1])) {
            swap(val[p], val[p - 1]);
            swap(idx[p], idx[p - 1]);
            --p;
        }
        return;
    }
    if (!better(v, val[K - 1])) return;

    idx[K - 1] = j;
    val[K - 1] = v;
    int p = K - 1;
    while (p > 0 && better(val[p], val[p - 1])) {
        swap(val[p], val[p - 1]);
        swap(idx[p], idx[p - 1]);
        --p;
    }
}

static TopKLists buildTopKLists(
    const vector<vector<double>>& dist,
    const vector<vector<double>>& x,
    const vector<vector<double>>& rc,
    int K)
{
    TopKLists out;
    out.n = (int)dist.size();
    out.K = K;
    const int n = out.n;

    out.distNN.assign(n, {});
    out.rcNN.assign(n, {});
    out.lpNN.assign(n, {});

    for (int i = 0; i < n; ++i) {
        vector<int> dIdx; dIdx.reserve(K);
        vector<double> dVal; dVal.reserve(K);

        vector<int> rIdx; rIdx.reserve(K);
        vector<double> rVal; rVal.reserve(K);

        vector<int> xIdx; xIdx.reserve(K);
        vector<double> xVal; xVal.reserve(K);

        for (int j = 0; j < n; ++j) {
            if (j == i) continue;

            pushTopK(dIdx, dVal, j, dist[i][j], K, [](double a, double b) { return a < b; }); // smaller better

            pushTopK(rIdx, rVal, j, getRCfull(rc, i, j), K, [](double a, double b) { return a < b; }); // smaller better

            pushTopK(xIdx, xVal, j, getXfull(x, i, j), K, [](double a, double b) { return a > b; }); // larger better
        }

        out.distNN[i] = move(dIdx);
        out.rcNN[i] = move(rIdx);
        out.lpNN[i] = move(xIdx);
    }

    return out;
}

// ================== Top-K query helpers ==============
static double avg2(double a, double b) { return 0.5 * (a + b); }

// OUT: need 2 vertices NOT in S
static void twoMin_OUT_fromTopK_dist(
    int i, const MarkSet& ms, int n,
    const vector<vector<double>>& dist,
    const vector<int>& nnDist,
    double& out1, double& out2)
{
    int found = 0;
    out1 = out2 = 999999999999999999999999999999999; //large number

    // fast path: scan Top-K
    for (int t = 0; t < (int)nnDist.size() && found < 2; ++t) {
        const int j = nnDist[t];
        if (ms.has(j)) continue;
        const double v = dist[i][j];
        if (found == 0) out1 = v;
        else out2 = v;
        ++found;
    }

    // fallback (we rarely do): full scan
    if (found < 2) {
        out1 = out2 = 999999999999999999999999999999999; //large number
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            if (ms.has(j)) continue;
            const double v = dist[i][j];
            if (v < out1) { out2 = out1; out1 = v; }
            else if (v < out2) { out2 = v; }
        }
    }
}

// ----- FIRST IN / FIRST OUT helpers (Top-K with fallback) -----

static double first_OUT_dist_fromTopK(
    int i, const MarkSet& ms, int n,  const vector<vector<double>>& dist,  const vector<int>& nnDist){
    for (int t = 0; t < (int)nnDist.size(); ++t) {
        int j = nnDist[t];
        if (ms.has(j)) continue;
        return dist[i][j];
    }
    // fallback full scan
    double best = 9999999999999999999999999999999999999999; 
    for (int j = 0; j < n; ++j) {
        if (j == i) continue;
        if (ms.has(j)) continue;
        double v = dist[i][j];
        if (v < best) best = v;
    }
    return best;
}

static double first_IN_dist_fromTopK(
    int i, const MarkSet& ms, const vector<int>& S, const vector<vector<double>>& dist, const vector<int>& nnDist){
    for (int t = 0; t < (int)nnDist.size(); ++t) {
        int j = nnDist[t];
        if (j == i) continue;
        if (!ms.has(j)) continue;
        return dist[i][j];
    }
    // fallbackkk scan within S
    double best = 9999999999999999999999999999999999999999; 
    for (int j : S) {
        if (j == i) continue;

        double v = dist[i][j];
        if (v < best) best = v;
    }
    return best;
}

static double first_OUT_rc_fromTopK(int i, const MarkSet& ms, int n, const vector<vector<double>>& rc, const vector<int>& nnRC){
    for (int t = 0; t < (int)nnRC.size(); ++t) {
        int j = nnRC[t];
        if (ms.has(j)) continue;
        return getRCfull(rc, i, j);
    }
    double best = 9999999999999999999999999999999999999999; 
    for (int j = 0; j < n; ++j) {
        if (j == i) continue;
        if (ms.has(j)) continue;
        double v = getRCfull(rc, i, j);

        if (v < best) best = v;
    }
    return best;
}

static double first_IN_rc_fromTopK(int i, const MarkSet& ms, const vector<int>& S, const vector<vector<double>>& rc, const vector<int>& nnRC){
    for (int t = 0; t < (int)nnRC.size(); ++t) {
        int j = nnRC[t];
        if (j == i) continue;
        if (!ms.has(j)) continue;
        return getRCfull(rc, i, j);
    }
    double best = 9999999999999999999999999999999999999999; 
    for (int j : S) {

        if (j == i) continue;
        double v = getRCfull(rc, i, j);
        if (v < best) best = v;
    }
    return best;
}

static double first_OUT_lp_fromTopK( int i, const MarkSet& ms, int n, const vector<vector<double>>& x, const vector<int>& nnLP){
    for (int t = 0; t < (int)nnLP.size(); ++t) {
        int j = nnLP[t];
        if (ms.has(j)) continue;
        return getXfull(x, i, j);
    }
    double best = -9999999999999999999999999999999999999999; 

    for (int j = 0; j < n; ++j) {
        if (j == i) continue;
        if (ms.has(j)) continue;
        double v = getXfull(x, i, j);
        if (v > best) best = v;
    }
    return best;
}

static double first_IN_lp_fromTopK(int i, const MarkSet& ms, const vector<int>& S, const vector<vector<double>>& x, const vector<int>& nnLP){
    for (int t = 0; t < (int)nnLP.size(); ++t) {
        int j = nnLP[t];
        if (j == i) continue;
        if (!ms.has(j)) continue;
        return getXfull(x, i, j);
    }
    double best = -9999999999999999999999999999999999999999; 

    for (int j : S) {
        if (j == i) continue;
        double v = getXfull(x, i, j);
        if (v > best) best = v;
    }
    return best;
}


static void twoMin_OUT_fromTopK_rc(int i, const MarkSet& ms, int n, const vector<vector<double>>& rc, const vector<int>& nnRC, double& out1, double& out2){
    int found = 0;
    out1 = out2 = 9999999999999999999999999999999999999999; 

    for (int t = 0; t < (int)nnRC.size() && found < 2; ++t) {
        const int j = nnRC[t];
        if (ms.has(j)) continue;
        const double v = getRCfull(rc, i, j);
        if (found == 0) out1 = v;
        else out2 = v;
        ++found;
    }

    if (found < 2) {
        out1 = out2 = 9999999999999999999999999999999999999999; 
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            if (ms.has(j)) continue;
            const double v = getRCfull(rc, i, j);
            if (v < out1) { out2 = out1; out1 = v; }
            else if (v < out2) { out2 = v; }
        }
    }
}

static void twoMax_OUT_fromTopK_lp(int i, const MarkSet& ms, int n, const vector<vector<double>>& x, const vector<int>& nnLP, double& out1, double& out2){
    int found = 0;
    out1 = out2 = -9999999999999999999999999999999999999999; 

    for (int t = 0; t < (int)nnLP.size() && found < 2; ++t) {
        const int j = nnLP[t];
        if (ms.has(j)) continue;
        const double v = getXfull(x, i, j);
        if (found == 0) out1 = v;
        else out2 = v;
        ++found;
    }

    if (found < 2) {
        out1 = out2 = -9999999999999999999999999999999999999999; 
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            if (ms.has(j)) continue;
            const double v = getXfull(x, i, j);
            if (v > out1) { out2 = out1; out1 = v; }
            else if (v > out2) { out2 = v; }
        }
    }
}

// IN: need 2 vertices IN S (excluding i)
static void twoMin_IN_fromTopK_dist(int i, const MarkSet& ms, int n, const vector<vector<double>>& dist, const vector<int>& nnDist, double& in1, double& in2){
    int found = 0;
    in1 = in2 = 999999999999999999999999999999999; //large num

    for (int t = 0; t < (int)nnDist.size() && found < 2; ++t) {
        const int j = nnDist[t];
        if (j == i) continue;
        if (!ms.has(j)) continue;
        const double v = dist[i][j];
        if (found == 0) in1 = v;
        else in2 = v;
        ++found;
    }

    if (found < 2) {
        in1 = in2 = 9999999999999999999999999999999999999999; 
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            if (!ms.has(j)) continue;
            const double v = dist[i][j];
            if (v < in1) { in2 = in1; in1 = v; }
            else if (v < in2) { in2 = v; }
        }
    }
}

static void twoMin_IN_fromTopK_rc( int i, const MarkSet& ms, int n, const vector<vector<double>>& rc, const vector<int>& nnRC, double& in1, double& in2){
    int found = 0;
    in1 = in2 = 9999999999999999999999999999999999999999; 

    for (int t = 0; t < (int)nnRC.size() && found < 2; ++t) {
        const int j = nnRC[t];
        if (j == i) continue;
        if (!ms.has(j)) continue;
        const double v = getRCfull(rc, i, j);
        if (found == 0) in1 = v;
        else in2 = v;
        ++found;
    }

    if (found < 2) {
        in1 = in2 = 9999999999999999999999999999999999999999; 
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            if (!ms.has(j)) continue;
            const double v = getRCfull(rc, i, j);
            if (v < in1) { in2 = in1; in1 = v; }
            else if (v < in2) { in2 = v; }
        }
    }
}

static void twoMax_IN_fromTopK_lp(int i, const MarkSet& ms, int n, const vector<vector<double>>& x, const vector<int>& nnLP, double& in1, double& in2){
    int found = 0;
    in1 = in2 = -9999999999999999999999999999999999999999; 

    for (int t = 0; t < (int)nnLP.size() && found < 2; ++t) {
        const int j = nnLP[t];
        if (j == i) continue;
        if (!ms.has(j)) continue;
        const double v = getXfull(x, i, j);
        if (found == 0) in1 = v;
        else in2 = v;
        ++found;
    }

    if (found < 2) {
        in1 = in2 = -9999999999999999999999999999999999999999; 
        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            if (!ms.has(j)) continue;
            const double v = getXfull(x, i, j);
            if (v > in1) { in2 = in1; in1 = v; }
            else if (v > in2) { in2 = v; }
        }
    }
}


// =================== TESTING CONFIGURATION =================================
constexpr int NOT_REDUCED_AND_SIMPLE = 0;

constexpr int XOR2_ACTIVE = 1;
constexpr int XOR3_ACTIVE = 1;
constexpr int XOR4_ACTIVE = 0;
constexpr int XOR5_ACTIVE = 0;

constexpr double RC_LIM = 0.0;
constexpr double MULTIPLIER_MST_EDGE_INCLUSION = 10.0;


// =================== INTERNAL HELPERS (only in this .cpp) ===================

long long encodeEdge(int u, int v) {
    if (u > v) swap(u, v);
    return ((long long)u << 32) | (unsigned int)v;
}


// Hash & equality for bitset<BSZ> (for unordered_set)
struct BitsetHash {
    size_t operator()(const bitset<BSZ>& b) const noexcept {
        constexpr size_t WORDS = (BSZ + 63) / 64; // check how many 64-bit blocks do we need to hold all the bits?
        const uint64_t* ptr = reinterpret_cast<const uint64_t*>(&b); // ok now treat the memory of this bitset as if it were an array of 64-bit integers.
        size_t h = 0; //initiate the running hash value by starting from 0

        for (size_t i = 0; i < WORDS; ++i) { // go through every 64-bit block in the bitset.
            uint64_t word = ptr[i];
            // standard hash-combine (take the current block and mix it into the running hash.)
            h ^= hash<uint64_t>()(word) + 0x9e3779b97f4a7c15ULL + (h << 6)+ (h >> 2); // 0x9e3779b97f4a7c15ULL is used (classic in hashing to spread valuess)
        }
        return h;
    }
};

struct BitsetEq {
    bool operator()(const bitset<BSZ>& a,
        const bitset<BSZ>& b) const {return a == b;}
};

struct VecHash {
    size_t operator()(const vector<int>& v) const noexcept {
        size_t h = 0; // running hash value

        for (int x : v) { // go through each integer in the vector
            h ^= hash<int>()(x) + 0x9e3779b97f4a7c15ULL + (h << 6)+ (h >> 2); // hash the current int, mixing constant, mix in previous hash shifted left and right
        }
        return h; // final hash for the whole vector
    }
};

struct VecEq {
    bool operator()(const vector<int>& a,
        const vector<int>& b) const noexcept {
        return a == b;
    }
};

// =================== LP RELAXATION ===================
pair<vector<vector<int>>, vector<vector<int>>>
SolveLPRelaxation(Instance& inst)
{
    try {
        GRBEnv env = GRBEnv();
        env.set(GRB_IntParam_OutputFlag, 1);
        GRBModel model = GRBModel(env);

        int n = inst.nbCity;
        vector<vector<GRBVar>> x(n, vector<GRBVar>(n));

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                x[i][j] = model.addVar(0, 1, inst.dist[i][j], GRB_CONTINUOUS);
            }
        }
        model.update();

        for (int v = 0; v < n; v++) {
            GRBLinExpr deg(0);
            for (int i = 0; i < v; i++) deg += x[i][v];
            for (int j = v + 1; j < n; j++) deg += x[v][j];
            model.addConstr(deg == 2);
        }

        model.getEnv().set(GRB_IntParam_Threads, 1);
        model.optimize();

        vector<vector<int>> skeleton, remaining_edges;

        // ---- temp storage for pruning ----
        struct SEdge {
            int i, j;
            double rc;
            double dist;
        };
        vector<SEdge> sk; sk.reserve((size_t)n * (n - 1) / 2);

        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            //save the LP relaxation values
            inst.lp_x.assign(n, vector<double>(n, 0.0)); //save x values
            inst.lp_rc.assign(n, vector<double>(n, 0.0)); // save rc valyues
            for (int i = 0; i < n; ++i) {
                for (int j = i + 1; j < n; ++j) {
                    inst.lp_x[i][j] = x[i][j].get(GRB_DoubleAttr_X);
                    inst.lp_rc[i][j] = x[i][j].get(GRB_DoubleAttr_RC);
                }
            }
            const double lpObj = model.get(GRB_DoubleAttr_ObjVal);
            inst.lp_obj = lpObj;
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    double rc = x[i][j].get(GRB_DoubleAttr_RC);

                    if (rc <= RC_LIM * lpObj) {
                        sk.push_back({ i, j, rc, inst.dist[i][j] });
                    }
                    else {
                        remaining_edges.push_back({ i, j });
                    }
                }
            }

            // ---- apply skeleton-size rule ----
            const double totEdges = (double)n * (n - 1) / 2.0;
            const double limitD = 2.2635 * pow(totEdges, 0.4913); // thr trend line we found
            //const double limitD = (2.2635 + 0.1) * pow(totEdges, (0.4913 - 0.02)); // intense trend line
            cout << sk.size() << " " << limitD << " " << totEdges << endl;

            // round the value
            const size_t limit = (size_t)floor(limitD);

            if (sk.size() > limit) {
                // Only remove among rc == 0 edges, removing highest dist first
                constexpr double RC_ZERO_EPS = 1e-12;

                // indices of candidates in sk
                vector<int> cand;
                cand.reserve(sk.size());

                for (int idx = 0; idx < (int)sk.size(); ++idx) {
                    if (abs(sk[idx].rc) <= RC_ZERO_EPS) {
                        cand.push_back(idx);
                    }
                }

                // sort candidates by dist descending
                sort(cand.begin(), cand.end(), [&](int a, int b) {
                    return sk[a].dist > sk[b].dist;
                    });

                // mark removals
                vector<char> removed(sk.size(), 0);

                size_t needRemove = sk.size() - limit;
                for (int t = 0; t < (int)cand.size() && needRemove > 0; ++t) {
                    removed[cand[t]] = 1;
                    --needRemove;
                }

                // If needRemove > 0 here, that means not enough rc==0 edges existed.
                // You can decide what to do: stop, or continue removing rc>0 edges too.
                // For now, we only remove rc==0
                // 
                // rebuild sk keeping non-removed
                vector<SEdge> sk2;
                sk2.reserve(sk.size());
                for (size_t idx = 0; idx < sk.size(); ++idx) {
                    if (!removed[idx]) sk2.push_back(sk[idx]);
                }
                sk.swap(sk2);
            }

            // finally build skeleton from pruned sk
            skeleton.reserve(sk.size());
            for (auto& e : sk) skeleton.push_back({ e.i, e.j });
        }

        return { skeleton, remaining_edges };
    }
    catch (GRBException& e) {
        cerr << "\n=== Gurobi Exception in SolveLPRelaxation ===\n";
        cerr << "Error code: " << e.getErrorCode() << "\n";
        cerr << "Message:    " << e.getMessage() << "\n";
        cerr << "============================================\n";
        throw;
    }
    catch (exception& e) {
        cerr << "\n=== exception in SolveLPRelaxation ===\n";
        cerr << e.what() << "\n";
        cerr << "===========================================\n";
        throw;
    }
    catch (...) {
        cerr << "\n=== Unknown Exception in SolveLPRelaxation ===\n";
        throw;
    }
}


// =================== SPANNING FOREST =============

// Build spanning forest from edge list
vector<pair<int, int>>
buildSpanningForest(const vector<vector<int>>& skeleton,  int n,  vector<int>& component){
    using namespace std;

    // Initialize component array
    component.assign(n, -1);

    // Build adjacency list from skeleton
    vector<vector<int>> adj(n);
    adj.reserve(n);

    for (auto& e : skeleton) {
        int u = e[0], v = e[1];
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    vector<pair<int, int>> forest;  // output tree edges
    forest.reserve(n - 1);

    vector<char> visited(n, 0);
    int compId = 0;

    vector<int> st;
    st.reserve(n);

    // Explore each connected component using DFS
    for (int s = 0; s < n; ++s) {
        if (visited[s]) continue;

        st.clear();
        st.push_back(s);
        visited[s] = 1;
        component[s] = compId;

        while (!st.empty()) {
            int u = st.back();
            st.pop_back();

            for (int v : adj[u]) {
                if (!visited[v]) {
                    visited[v] = 1;
                    component[v] = compId;

                    // Tree edge
                    forest.emplace_back(u, v);

                    st.push_back(v);
                }
            }
        }
        compId++;
    }
    return forest;
}

// =================== FAST SIMPLE-CYCLE CHECK (BLOCK-LOCAL) ===================

bool isSimpleCycle_fast(const bitset<BSZ>& b,  const vector<int>& U,  const vector<int>& V, const vector<vector<int>>& inc){
    int E = (int)U.size();
    int Vn = (int)inc.size();

    vector<int>  degree(Vn, 0);
    vector<char> used(Vn, 0);

    int edgeCount = 0;

    // degree check
    for (int ei = 0; ei < E; ++ei) {
        if (!b.test(ei)) continue;
        int u = U[ei];
        int v = V[ei];
        degree[u]++;
        degree[v]++;
        used[u] = 1;
        used[v] = 1;
        edgeCount++;
    }

    if (edgeCount == 0) return false;

    int usedV = 0;
    int start = -1;
    for (int i = 0; i < Vn; ++i) {
        if (!used[i]) continue;
        if (degree[i] != 2) return false;
        usedV++;
        if (start == -1) start = i;
    }

    if (edgeCount != usedV) return false;

    // connectivity test via DFS using incidence
    vector<char> visited(Vn, 0);
    stack<int> st;
    st.push(start);

    int cnt = 0;
    while (!st.empty()) {
        int u = st.top(); st.pop();
        if (visited[u]) continue;
        visited[u] = 1;
        cnt++;

        for (int ei : inc[u]) {
            if (!b.test(ei)) continue;
            int w = (U[ei] == u ? V[ei] : U[ei]);
            if (!visited[w]) st.push(w);
        }
    }

    return cnt == usedV;
}

// =================== TARJAN BCC (EDGE-BICONNECTED COMPONENTS) ===================

BCCResult tarjanBCC(const vector<vector<int>>& adj, int n)
{
    using namespace std;
    vector<int> disc(n, -1), low(n, -1);
    vector<pair<int, int>> st;
    int timer = 0;

    BCCResult out;

    function<void(int, int)> dfs = [&](int u, int parent)
    {
        disc[u] = low[u] = ++timer;

        for (int v : adj[u])
        {
            if (disc[v] == -1)
            {
                st.emplace_back(u, v);
                dfs(v, u);
                low[u] = min(low[u], low[v]);

                if (low[v] >= disc[u]) {
                    vector<pair<int, int>> edges;
                    unordered_set<int> verts;
                    while (!st.empty()) {
                        auto e = st.back(); st.pop_back();
                        edges.push_back(e);
                        verts.insert(e.first);
                        verts.insert(e.second);
                        if (e.first == u && e.second == v) break;
                    }
                    out.blocks.push_back(move(edges));
                    out.blockVerts.push_back(
                        vector<int>(verts.begin(), verts.end())
                    );
                }
            }
            else if (v != parent && disc[v] < disc[u]) {
                st.emplace_back(u, v);
                low[u] = min(low[u], disc[v]);
            }
        }
    };

    for (int i = 0; i < n; ++i)
        if (disc[i] == -1)
            dfs(i, -1);

    return out;
}

// =================== FUNDAMENTAL CYCLES PER BLOCK (LCA) ===================
//NOTE THAT YOU ELIMINATE ALREADY FUNDAMENTAL CYCLES THAT HAS SIZE BIGGER THAN N/2

// ---------- DSU (Kruskal) ----------
struct DSU {
    vector<int> p, r;
    explicit DSU(int n) : p(n), r(n, 0) { iota(p.begin(), p.end(), 0); }

    int find(int a) { return p[a] == a ? a : p[a] = find(p[a]); }

    bool unite(int a, int b) {
        a = find(a); b = find(b);
        if (a == b) return false;
        if (r[a] < r[b]) swap(a, b);
        p[b] = a;
        if (r[a] == r[b]) r[a]++;
        return true;
    }
};

static long long encEdgeLL(int u, int v) {
    if (u > v) swap(u, v);
    return ((long long)u << 32) | (unsigned int)v;
}

// ---------- MST fundamental cycles per block ----------
vector<vector<int>> computeFundamentalCycles_Block_MST(
    const vector<pair<int, int>>& blockEdges,     // GLOBAL vertex ids
    const vector<int>& blockVerts,     // GLOBAL vertex ids
    const vector<vector<double>>& dist,           // inst.dist
    const int& maxEdgesPerCycle)
{
    using namespace std;

    int V = (int)blockVerts.size();
    if (V <= 1) return {};

    // ---- Global -> Local map
    unordered_map<int, int> g2l;
    g2l.reserve(V);
    for (int i = 0; i < V; i++) g2l[blockVerts[i]] = i;

    // ---- Build weighted edge list for Kruskal (use original dist)
    struct WEdge { double w; int u; int v; }; // u,v are GLOBAL ids
    vector<WEdge> edges;
    edges.reserve(blockEdges.size());
    for (auto& e : blockEdges) {
        int u = e.first, v = e.second;
        edges.push_back({ dist[u][v], u, v });
    }

    sort(edges.begin(), edges.end(),
        [](const WEdge& a, const WEdge& b) { return a.w < b.w; });

    // ---- Kruskal MST on block vertices
    DSU dsu(V);
    vector<pair<int, int>> mstEdges; // GLOBAL ids
    mstEdges.reserve(V - 1);

    for (auto& we : edges) {
        auto itu = g2l.find(we.u);
        auto itv = g2l.find(we.v);
        if (itu == g2l.end() || itv == g2l.end()) continue;

        int lu = itu->second;
        int lv = itv->second;
        if (dsu.unite(lu, lv)) {
            mstEdges.emplace_back(we.u, we.v);
            if ((int)mstEdges.size() == V - 1) break;
        }
    }

    // If this happens, something is off: the "block" edge list isn't connected.
    if ((int)mstEdges.size() != V - 1) {
        return {};
    }

    // ---- Build MST adjacency (LOCAL indices)
    vector<vector<int>> treeAdj(V);
    for (auto& e : mstEdges) {
        int u = g2l[e.first];
        int v = g2l[e.second];
        treeAdj[u].push_back(v);
        treeAdj[v].push_back(u);
    }

    // ---- BFS on MST to build parent + depth
    vector<int> parent(V, -1), depth(V, 0);
    queue<int> q;
    parent[0] = 0;
    q.push(0);

    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : treeAdj[u]) {
            if (parent[v] != -1) continue;
            parent[v] = u;
            depth[v] = depth[u] + 1;
            q.push(v);
        }
    }
    parent[0] = -1;

    // Safety: ensure BFS reached all vertices (should be true for MST)
    for (int i = 0; i < V; ++i) {
        if (i != 0 && parent[i] == -1) return {};
    }

    // ---- Binary lifting table for LCA
    int LOG = (int)ceil(log2((double)V)) + 2;
    vector<vector<int>> up(V, vector<int>(LOG, -1));
    for (int i = 0; i < V; i++) up[i][0] = parent[i];

    for (int k = 1; k < LOG; k++) {
        for (int i = 0; i < V; i++) {
            int mid = up[i][k - 1];
            up[i][k] = (mid == -1 ? -1 : up[mid][k - 1]);
        }
    }

    auto LCA = [&](int a, int b) {
        if (depth[a] < depth[b]) swap(a, b);
        int diff = depth[a] - depth[b];
        for (int k = 0; k < LOG; k++) {
            if (diff & (1 << k)) a = up[a][k];
        }
        if (a == b) return a;

        for (int k = LOG - 1; k >= 0; k--) {
            if (up[a][k] != -1 && up[a][k] != up[b][k]) {
                a = up[a][k];
                b = up[b][k];
            }
        }
        return parent[a];
    };

    // ---- Tree edge set (GLOBAL ids) for skipping MST edges
    unordered_set<long long> treeSet;
    treeSet.reserve(mstEdges.size() * 2);
    for (auto& e : mstEdges) {
        treeSet.insert(encEdgeLL(e.first, e.second));
    }

    // ---- Build fundamental cycles: for every NON-tree edge, take MST path + chord
    vector<vector<int>> cycles;
    cycles.reserve(blockEdges.size()); // rough

    for (auto& e : blockEdges) {
        long long key = encEdgeLL(e.first, e.second);
        if (treeSet.count(key)) continue; // skip MST tree edges

        int u = g2l[e.first]; // local
        int v = g2l[e.second];
        int l = LCA(u, v);

        vector<int> cyc;
        cyc.reserve(depth[u] + depth[v] + 2);

        // u -> l
        int x = u;
        while (x != l) {
            cyc.push_back(blockVerts[x]); // GLOBAL id
            x = parent[x];
        }
        cyc.push_back(blockVerts[l]);

        // v -> l (reverse)
        vector<int> tail;
        x = v;
        while (x != l) {
            tail.push_back(blockVerts[x]);
            x = parent[x];
        }
        reverse(tail.begin(), tail.end());
        cyc.insert(cyc.end(), tail.begin(), tail.end());

        // fundamental cycle length in edges equals number of vertices in this list
        // (since chord closes it). Your downstream encoder closes cyclically.
        if ((int)cyc.size() <= maxEdgesPerCycle) {
            cycles.push_back(move(cyc));
        }
    }

    return cycles;
}



vector<vector<int>> computeFundamentalCycles_Block(const vector<pair<int, int>>& blockEdges,  const vector<int>& blockVerts, const int& maxEdgesPerCycle){
    using namespace std;

    int V = (int)blockVerts.size();
    if (V <= 1) return {};

    // --- Global -> Local map
    unordered_map<int, int> g2l;
    g2l.reserve(V);
    for (int i = 0; i < V; i++) g2l[blockVerts[i]] = i;

    // --- Build adjacency
    vector<vector<int>> adj(V);
    for (auto& e : blockEdges) {
        int u = g2l[e.first];
        int v = g2l[e.second];
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    // --- BFS to build parent + depth (no recursion!)
    vector<int> parent(V, -1), depth(V, 0);
    int LOG = (int)ceil(log2((double)V)) + 2;
    vector<vector<int>> up(V, vector<int>(LOG, -1));

    queue<int> q;
    q.push(0);
    parent[0] = -1;
    depth[0] = 0;

    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : adj[u]) {
            if (v == parent[u]) continue;
            if (parent[v] != -1) continue;     // already visited
            parent[v] = u;
            depth[v] = depth[u] + 1;
            q.push(v);
        }
    }

    // --- Build binary lifting table
    for (int i = 0; i < V; i++) up[i][0] = parent[i];
    for (int k = 1; k < LOG; k++) {
        for (int i = 0; i < V; i++) {
            if (up[i][k - 1] == -1) up[i][k] = -1;
            else up[i][k] = up[up[i][k - 1]][k - 1];
        }
    }

    // --- LCA
    auto LCA = [&](int a, int b) {
        if (depth[a] < depth[b]) swap(a, b);
        int diff = depth[a] - depth[b];
        for (int k = 0; k < LOG; k++) {
            if (diff & (1 << k)) a = up[a][k];
        }
        if (a == b) return a;

        for (int k = LOG - 1; k >= 0; k--) {
            if (up[a][k] != -1 && up[a][k] != up[b][k]) {
                a = up[a][k];
                b = up[b][k];
            }
        }
        return parent[a];
    };

    // --- Tree edges set (global)
    unordered_set<long long> treeSet;
    treeSet.reserve(V);

    auto enc = [&](int u, int v) {
        if (u > v) swap(u, v);
        return ((long long)u << 32) | (unsigned long long)v;
    };

    for (int v = 0; v < V; v++) {
        if (parent[v] != -1) {
            int gu = blockVerts[v];
            int gp = blockVerts[parent[v]];
            treeSet.insert(enc(gu, gp));
        }
    }

    // --- Build fundamental cycles
    vector<vector<int>> cycles;

    for (auto& e : blockEdges) {
        long long key = enc(e.first, e.second);
        if (treeSet.count(key)) continue; // skip tree edges

        int u = g2l[e.first];
        int v = g2l[e.second];
        int l = LCA(u, v);

        vector<int> cyc;

        int x = u;
        while (x != l) {
            cyc.push_back(blockVerts[x]);
            x = parent[x];
        }
        cyc.push_back(blockVerts[l]);

        vector<int> tail;
        x = v;
        while (x != l) {
            tail.push_back(blockVerts[x]);
            x = parent[x];
        }
        reverse(tail.begin(), tail.end());
        cyc.insert(cyc.end(), tail.begin(), tail.end());

        if (NOT_REDUCED_AND_SIMPLE) cycles.push_back(move(cyc)); // Include all len fundamental cycles.
        else {
            if (cyc.size() <= maxEdgesPerCycle) cycles.push_back(move(cyc)); // Include only fundamental cycles that has len<=n/2
        }
    }

    return cycles;
}

// =================== ENCODE CYCLES AS BITSETS PER BLOCK ===================

vector<bitset<BSZ>> encodeCycles_Block( const vector<vector<int>>& cycles, const vector<pair<int, int>>& blockEdges){
    using namespace std;

    unordered_map<long long, int> edgeMap;
    edgeMap.reserve(blockEdges.size() * 2);

    for (int i = 0; i < (int)blockEdges.size(); ++i) {
        int u = blockEdges[i].first;
        int v = blockEdges[i].second;
        if (u > v) swap(u, v);
        edgeMap[((long long)u << 32) | (unsigned int)v] = i;
    }

    vector<bitset<BSZ>> out;

    for (auto& cyc : cycles) {
        bitset<BSZ> b;
        int m = (int)cyc.size();
        if (m < 2) continue;

        for (int i = 0; i < m; ++i) {
            int u = cyc[i];
            int v = cyc[(i + 1) % m];
            if (u > v) swap(u, v);
            long long key = ((long long)u << 32) | (unsigned int)v;
            auto it = edgeMap.find(key);
            if (it != edgeMap.end()) b.set(it->second);
        }

        if (b.any()) out.push_back(b);
    }

    return out;
}

// =================== XOR-2 / XOR-3 / XOR-4 PER BLOCK (old, debug) ===================

vector<bitset<BSZ>> generateXORCycles_Block_old( const vector<bitset<BSZ>>& fundBits, const vector<int>& U_local, const vector<int>& V_local,   const vector<vector<int>>& incEdges,int maxEdgesPerCycle)
{
    using namespace std;

    unordered_set<bitset<BSZ>, BitsetHash, BitsetEq> seen;

    int F = (int)fundBits.size();
    if (F == 0) return {};

    if (XOR2_ACTIVE) {
        for (int i = 0; i < F; ++i) {
            for (int j = i + 1; j < F; ++j) {
                if ((fundBits[i] & fundBits[j]).none()) continue;

                bitset<BSZ> x = fundBits[i] ^ fundBits[j];
                int cnt = (int)x.count();
                if (cnt == 0 || cnt > maxEdgesPerCycle) continue;

                if (isSimpleCycle_fast(x, U_local, V_local, incEdges))
                    seen.insert(x);
            }
        }
    }

    if (XOR3_ACTIVE) {
        for (int i = 0; i < F; ++i) {
            for (int j = i + 1; j < F; ++j) {
                for (int k = j + 1; k < F; ++k) {
                    int Ai = fundBits[i].count();
                    int Aj = fundBits[j].count();
                    int Ak = fundBits[k].count();

                    int o_ij = (fundBits[i] & fundBits[j]).count();
                    int o_ik = (fundBits[i] & fundBits[k]).count();
                    int o_jk = (fundBits[j] & fundBits[k]).count();

                    int LB = Ai + Aj + Ak - 2 * (o_ij + o_ik + o_jk);
                    if (LB > maxEdgesPerCycle) continue;

                    bitset<BSZ> x = fundBits[i] ^ fundBits[j] ^ fundBits[k];
                    int cnt = (int)x.count();
                    if (cnt == 0 || cnt > maxEdgesPerCycle) continue;

                    if (isSimpleCycle_fast(x, U_local, V_local, incEdges))
                        seen.insert(x);
                }
            }
        }
    }

    if (XOR4_ACTIVE) {
        for (int i = 0; i < F; ++i) {
            for (int j = i + 1; j < F; ++j) {
                for (int k = j + 1; k < F; ++k) {
                    for (int l = k + 1; l < F; ++l) {
                        int Ai = fundBits[i].count();
                        int Aj = fundBits[j].count();
                        int Ak = fundBits[k].count();
                        int Al = fundBits[l].count();

                        int o_ij = (fundBits[i] & fundBits[j]).count();
                        int o_ik = (fundBits[i] & fundBits[k]).count();
                        int o_il = (fundBits[i] & fundBits[l]).count();
                        int o_jk = (fundBits[j] & fundBits[k]).count();
                        int o_jl = (fundBits[j] & fundBits[l]).count();
                        int o_kl = (fundBits[k] & fundBits[l]).count();

                        int LB = Ai + Aj + Ak + Al
                            - 2 * (o_ij + o_ik + o_il + o_jk + o_jl + o_kl);
                        if (LB > maxEdgesPerCycle) continue;

                        bitset<BSZ> x = fundBits[i] ^ fundBits[j] ^ fundBits[k] ^ fundBits[l];
                        int cnt = (int)x.count();
                        if (cnt == 0 || cnt > maxEdgesPerCycle) continue;

                        if (isSimpleCycle_fast(x, U_local, V_local, incEdges))
                            seen.insert(x);
                    }
                }
            }
        }
    }

    return vector<bitset<BSZ>>(seen.begin(), seen.end());
}

vector<bitset<BSZ>> generateXORCycles_Block(const vector<bitset<BSZ>>& fundBits, const vector<int>& U_local, const vector<int>& V_local, const vector<vector<int>>& incEdges, int maxEdgesPerCycle)
{
    using namespace std;

    unordered_set<bitset<BSZ>, BitsetHash, BitsetEq> seen;

    const int F = (int)fundBits.size();
    if (F == 0) return {};

    vector<bitset<BSZ>> xor2Simple;
    vector<bitset<BSZ>> xor3Simple;
    vector<bitset<BSZ>> xor4Simple; // NEW (needed for XOR-5)

    if (XOR2_ACTIVE) {

        for (int i = 0; i < F; ++i) {
            for (int j = i + 1; j < F; ++j) {
                if ((fundBits[i] & fundBits[j]).none())
                    continue;

                bitset<BSZ> x = fundBits[i] ^ fundBits[j];

                int cnt = (int)x.count();
                if (cnt == 0 || cnt > maxEdgesPerCycle)
                    continue;

                if (!isSimpleCycle_fast(x, U_local, V_local, incEdges))
                    continue;

                auto [it, inserted] = seen.insert(x);
                if (inserted)
                    xor2Simple.push_back(x);
            }
        }
    }

    if (XOR3_ACTIVE) {

        //xor3Simple.reserve(xor2Simple.size() * 2); // lets not reserve any memory

        for (const auto& x2 : xor2Simple) {
            for (int k = 0; k < F; ++k) {
                if ((x2 & fundBits[k]).none())
                    continue;

                bitset<BSZ> x3 = x2 ^ fundBits[k];

                int cnt = (int)x3.count();
                if (cnt == 0 || cnt > maxEdgesPerCycle)
                    continue;

                if (!isSimpleCycle_fast(x3, U_local, V_local, incEdges))
                    continue;

                auto [it, inserted] = seen.insert(x3);
                if (inserted)
                    xor3Simple.push_back(x3);
            }
        }
    }

    if (XOR4_ACTIVE) {

        xor4Simple.reserve(xor3Simple.size()); // rough guess
        size_t xor4Count = 0;

        for (const auto& x3 : xor3Simple) {
            for (int k = 0; k < F; ++k) {
                if ((x3 & fundBits[k]).none())
                    continue;

                bitset<BSZ> x4 = x3 ^ fundBits[k];

                int cnt = (int)x4.count();
                if (cnt == 0 || cnt > maxEdgesPerCycle)
                    continue;

                if (!isSimpleCycle_fast(x4, U_local, V_local, incEdges))
                    continue;

                auto [it, inserted] = seen.insert(x4);
                if (inserted) {
                    xor4Count++;
                    xor4Simple.push_back(x4); //  keep for XOR-5 layer
                }
            }
        }
    }

    if (XOR5_ACTIVE) {

        size_t xor5Count = 0;

        for (const auto& x4 : xor4Simple) {
            for (int k = 0; k < F; ++k) {
                if ((x4 & fundBits[k]).none())
                    continue;

                bitset<BSZ> x5 = x4 ^ fundBits[k];

                int cnt = (int)x5.count();
                if (cnt == 0 || cnt > maxEdgesPerCycle)
                    continue;

                if (!isSimpleCycle_fast(x5, U_local, V_local, incEdges))
                    continue;

                auto [it, inserted] = seen.insert(x5);
                if (inserted)
                    xor5Count++;
            }
        }
    }

    return vector<bitset<BSZ>>(seen.begin(), seen.end());
}

// =================== NEW PIPELINE HELPERS ===================

vector<vector<int>> buildSkeletonAdjacency(const vector<vector<int>>& skeleton, int n)
{
    vector<vector<int>> adj(n);
    for (const auto& e : skeleton) {
        int u = e[0];
        int v = e[1];
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    return adj;
}

BlockData buildBlocksFromSkeleton(
    const vector<vector<int>>& skeleton,
    const vector<int>& component,
    int n,
    ofstream& log)
{
    using namespace std;

    BlockData out;

    // Build adjacency
    auto adj = buildSkeletonAdjacency(skeleton, n);

    // Tarjan BCC
    BCCResult bres = tarjanBCC(adj, n);

    int B = static_cast<int>(bres.blocks.size());
    out.blockEdgesAll.reserve(B);
    out.blockVertsAll.reserve(B);
    out.blockCompId.reserve(B);

    for (int bi = 0; bi < B; ++bi) {
        auto& edges = bres.blocks[bi];
        auto& verts = bres.blockVerts[bi];

        if (edges.size() < 3) {
            // no cycle possible in a BCC with fewer than 3 edges
            continue;
        }

        int u0 = verts[0];
        int cid = component[u0];

        out.blockEdgesAll.push_back(edges);
        out.blockVertsAll.push_back(verts);
        out.blockCompId.push_back(cid);
    }

    return out;
}

vector<bitset<BSZ>> generateCyclesForSingleBlock(const vector<pair<int, int>>& blockEdges, const vector<int>& blockVerts, int blockCompId, int maxEdgesPerCycle, Instance& inst){
    using namespace std;

    vector<bitset<BSZ>> result;

    int E = static_cast<int>(blockEdges.size());
    cout << "BlockEdges.size()=" << blockEdges.size() << " BSZ=" << BSZ << endl;
    if (E > BSZ) {
        cerr << "Warning ...\n" << flush;
        cerr << "Warning: block (comp " << blockCompId << ") has " << E << " edges, exceeding BSZ=" << BSZ << ". Skipping XOR-generation for this block.\n";
        return result;
    }

    int Vb = static_cast<int>(blockVerts.size());

    // global -> local vertex index
    unordered_map<int, int> g2l;
    g2l.reserve(Vb);
    for (int i = 0; i < Vb; i++) {
        g2l[blockVerts[i]] = i;
    }

    // U_local, V_local
    vector<int> U_local(E), V_local(E);
    for (int ei = 0; ei < E; ei++) {
        U_local[ei] = g2l[blockEdges[ei].first];
        V_local[ei] = g2l[blockEdges[ei].second];
    }

    // incidence list: vertex -> list of edge indices
    vector<vector<int>> incEdges(Vb);
    for (int ei = 0; ei < E; ei++) {
        incEdges[U_local[ei]].push_back(ei);
        incEdges[V_local[ei]].push_back(ei);
    }

    // 1. fundamental cycles (vertex lists)
    //auto fundCycles = computeFundamentalCycles_Block(blockEdges, blockVerts, maxEdgesPerCycle);
    auto fundCycles = computeFundamentalCycles_Block_MST(blockEdges, blockVerts, inst.dist, maxEdgesPerCycle);

    // 2. encode as bitsets
    auto fundBits = encodeCycles_Block(fundCycles, blockEdges);

    // 3. keep fundamental cycles (they should already be simple cycles)
    unordered_set<bitset<BSZ>, BitsetHash, BitsetEq> blockSet;
    //for (auto& fb : fundBits) blockSet.insert(fb);
    int kept = 0, dropped = 0;
    for (const auto& fb : fundBits)
    {
        int cnt = (int)fb.count();
        if (cnt == 0 || cnt > maxEdgesPerCycle) {
            dropped++;
            continue;
        }
        blockSet.insert(fb);
        kept++;
    }
    cout << "Fundamental cycles kept: " << kept << ", dropped: " << dropped << "\n";

    vector<bitset<BSZ>> xorCycles;
    if (NOT_REDUCED_AND_SIMPLE) {
        // 4. XOR2 / XOR3 / XOR4 cycles
        xorCycles = generateXORCycles_Block_old(
            fundBits,
            U_local,
            V_local,
            incEdges,
            maxEdgesPerCycle
        );
    }
    else
    {
        // 4. XOR2 / XOR3 / XOR4 cycles
        xorCycles = generateXORCycles_Block(
            fundBits,
            U_local,
            V_local,
            incEdges,
            maxEdgesPerCycle
        );
    }

    for (auto& x : xorCycles) {
        blockSet.insert(x);
    }

    // 5. store all cycles of this block
    result.assign(blockSet.begin(), blockSet.end());

    cout << "Block (comp " << blockCompId
        << ") has " << fundBits.size() << " fundamental cycles, "
        << result.size()
        << " simple cycles (fund+XOR2+XOR3+XOR4).\n";

    return result;
}

vector<vector<bitset<BSZ>>> generateCyclesForAllBlocks(const BlockData& blocks, int maxEdgesPerCycle, ofstream& log, Instance& inst){
    using namespace std;

    vector<vector<bitset<BSZ>>> allBlockCyclesBits(
        blocks.blockEdgesAll.size()
    );

    // Count how many components exist
    int maxComp = -1;
    for (int cid : blocks.blockCompId)
        if (cid > maxComp) maxComp = cid;

    int totalComponents = maxComp + 1;
    size_t totalBlocks = blocks.blockEdgesAll.size();

    // Count how many blocks per component
    vector<int> blocksPerComponent(totalComponents, 0);
    for (int cid : blocks.blockCompId)
        blocksPerComponent[cid]++;

    for (size_t bi = 0; bi < totalBlocks; ++bi)
    {
        const auto& edges = blocks.blockEdgesAll[bi];
        const auto& verts = blocks.blockVertsAll[bi];
        int compId = blocks.blockCompId[bi];

        int blocksInThisComponent = blocksPerComponent[compId];

        cout << "\n===== Processing Block "
            << bi << " / " << (totalBlocks - 1)
            << "  (Component " << compId
            << " / " << (totalComponents - 1)
            << ", Blocks in this component: " << blocksInThisComponent
            << ") =====\n";

        cout << "Vertices: " << verts.size()
            << " | Edges: " << edges.size() << "\n";

        auto cycles = generateCyclesForSingleBlock(
            edges,
            verts,
            compId,
            maxEdgesPerCycle,
            inst
        );

        cout << "Cycles found: " << cycles.size() << "";

        allBlockCyclesBits[bi] = move(cycles);
    }

    return allBlockCyclesBits;
}

RCCyclesResult buildGlobalRCCycles(const BlockData& blocks, const vector<vector<bitset<BSZ>>>& allBlockCyclesBits, ofstream& log){
    using namespace std;

    RCCyclesResult res;
    //res.cycles_vertices.reserve(500000); //old estimate

    unordered_set<vector<int>, VecHash, VecEq> globalSeen;
    //globalSeen.reserve(500000); //old estimate

    size_t totalCycles = 0;
    for (const auto& v : allBlockCyclesBits) totalCycles += v.size();

    res.cycles_vertices.reserve(totalCycles);
    globalSeen.reserve(totalCycles);


    for (size_t bi = 0; bi < blocks.blockEdgesAll.size(); ++bi)
    {
        const auto& bitsVec = allBlockCyclesBits[bi];
        const auto& edges = blocks.blockEdgesAll[bi];

        if (bitsVec.empty())
            continue;

        int Eb = static_cast<int>(edges.size());

        for (const auto& bs : bitsVec)
        {
            vector<int> verts;
            //verts.reserve(16); //old estimate
            int cntEdges = (int)bs.count(); //Each bitset cycle has cntEdges = bs.count() edges, and we push two endpoints per edge then dedupe, so the raw push count is 2*cntEdges.
            verts.reserve(2 * cntEdges);

            for (int ei = 0; ei < Eb; ++ei) {
                if (!bs.test(ei)) continue;
                verts.push_back(edges[ei].first);
                verts.push_back(edges[ei].second);
            }
            if (verts.empty()) continue;

            sort(verts.begin(), verts.end());
            verts.erase(unique(verts.begin(), verts.end()), verts.end());

            if (globalSeen.insert(verts).second) {
                res.cycles_vertices.push_back(move(verts));
            }
        }
    }

    return res;
}

// =================== MAIN ENTRY: Find_all_RC_CYCLES ===================

vector<vector<int>> Find_all_RC_CYCLES(Instance& inst)
{
    using namespace std;


    string savePath = "time/RC0_";
    //string Path2 = savePath + inst.name + "_time.txt";
    string Path2 = savePath + "time.txt";
    ofstream log(Path2, ios::out | ios::app);

    int n = inst.nbCity;

    // Step 0: LP & skeleton
    auto result = SolveLPRelaxation(inst);
    vector<vector<int>> skeleton = result.first;
    vector<vector<int>> rem_RC_ordered = result.second; // unused here

    // Step 1: Connected components via spanning forest
    vector<int> component;
    vector<pair<int, int>> forest = buildSpanningForest(skeleton, n, component);

    // optional: per-component vertex lists
    //unordered_map<int, vector<int>> verticesByComponent;
    //for (int v = 0; v < (int)component.size(); ++v)
    //    verticesByComponent[component[v]].push_back(v);

    // Step 2 4: BCC -> blocks
    BlockData blocks = buildBlocksFromSkeleton(
        skeleton,
        component,
        n,
        log
    );

    // Step 5: per-block cycles
    int maxEdgesPerCycle = inst.nbCity / 2;
    auto allBlockCyclesBits = generateCyclesForAllBlocks(
        blocks,
        maxEdgesPerCycle,
        log,
        inst
    );

    inst.nbRC_cycles_NOT_RED = allBlockCyclesBits.size();

    // Step 6: global RC cycles as vertex sets
    RCCyclesResult rcResult = buildGlobalRCCycles(
        blocks,
        allBlockCyclesBits,
        log
    );

    cout << "\n===== Global RC cycles (vertex sets) =====\n";
    cout << "Total unique RC cycles: " << rcResult.cycles_vertices.size() << "\n";
    inst.nbRC_cycles_unique = rcResult.cycles_vertices.size();
    /*
    for (size_t i = 0; i < rcResult.cycles_vertices.size(); ++i)
    {
        cout << "Cycle " << i << ": ";
        for (int v : rcResult.cycles_vertices[i])
            cout << v << " ";
        cout << "\n";
    }
    */

    log.close();

    return rcResult.cycles_vertices;
}

// =================== COMPARE PRECUTS VS RC CYCLES ===================

void CompareCstar(vector<vector<int>>& RC_cycles_vertices, Instance& inst, const string& type)
// Compare inst.precuts vs RC_cycles_vertices
{
    using namespace std;

    // ---- Build RC hash set ----
    unordered_set<string> RC_hash;
    RC_hash.reserve(RC_cycles_vertices.size() * 2);

    for (const auto& cyc : RC_cycles_vertices)
    {
        vector<int> v = cyc;
        sort(v.begin(), v.end());  // canonical form

        string key;
        key.reserve(v.size() * 4);
        for (int x : v) {
            key += to_string(x);
            key.push_back(',');
        }

        RC_hash.insert(move(key));
    }

    // ---- Compare with precuts ----
    int precut_total = (int)inst.precuts.size();
    int precut_match = 0;

    for (const auto& cyc : inst.precuts)
    {
        vector<int> v = cyc;
        sort(v.begin(), v.end());

        string key;
        key.reserve(v.size() * 4);
        for (int x : v) {
            key += to_string(x);
            key.push_back(',');
        }

        if (RC_hash.count(key))
            precut_match++;
    }

    // ---- Report ----
    double percentage =
        (precut_total == 0)
        ? 0.0
        : (100.0 * precut_match / precut_total);

    cout << "\n===== Precut Coverage Report =====\n";
    cout << "Precut cycles total    = " << precut_total << "\n";
    cout << "Precut cycles covered  = " << precut_match << "\n";
    cout << "Coverage percentage    = " << percentage << "%\n\n";

    if (type == "MST") {
        inst.nbMSTmatch = precut_match;
        inst.percMSTmatch = percentage;

    }
    else if (type == "RC") {
        inst.nbRCmatch = precut_match;
        inst.percRCmatch = percentage;
    }
    else if (type == "ALL") {
        inst.nbMatch = precut_match;
        inst.percMatch = percentage;
    }
}



// =================== MST & LCA HELPERS ===================

vector<pair<int, int>> buildMST(const vector<vector<double>>& dist) {
    int n = (int)dist.size();
    vector<double> minEdge(n, 1e18);
    vector<int>    parent(n, -1);
    vector<bool>   used(n, false);

    minEdge[0] = 0.0;
    vector<pair<int, int>> mstEdges;

    for (int i = 0; i < n; i++) {
        int v = -1;
        for (int j = 0; j < n; j++) {
            if (!used[j] && (v == -1 || minEdge[j] < minEdge[v])) {
                v = j;
            }
        }

        used[v] = true;
        if (parent[v] != -1) {
            mstEdges.push_back({ v, parent[v] });
        }
        for (int u = 0; u < n; u++) {
            if (dist[v][u] < minEdge[u]) {
                minEdge[u] = dist[v][u];
                parent[u] = v;
            }
        }
    }
    return mstEdges;
}

int getLCA(int a, int b,
    const vector<vector<int>>& up,
    const vector<int>& depth,
    int LOG)
{
    if (depth[a] < depth[b]) swap(a, b);
    int diff = depth[a] - depth[b];

    for (int i = LOG - 1; i >= 0; i--) {
        if (diff & (1 << i)) a = up[a][i];
    }
    if (a == b) return a;

    for (int i = LOG - 1; i >= 0; i--) {
        if (up[a][i] != -1 && up[a][i] != up[b][i]) {
            a = up[a][i];
            b = up[b][i];
        }
    }
    return up[a][0];
}

static double cliqueAverageDistance(const vector<int>& verts,
    const vector<vector<double>>& dist)
{
    const int m = (int)verts.size();
    if (m < 2) return 0.0;

    long double sum = 0.0L;
    long long pairs = 0;

    for (int a = 0; a < m; ++a) {
        int va = verts[a];
        for (int b = a + 1; b < m; ++b) {
            int vb = verts[b];
            sum += (long double)dist[va][vb];
            ++pairs;
        }
    }
    return pairs ? (double)(sum / (long double)pairs) : 0.0;
}

set<set<int>> Check_K_MSTCycles_Ordered_AvgCliqueLen(const int& needed, Instance& inst, Solution& sol)
{
    // Build "seen" from precuts once
    set<set<int>> seen;
    for (const auto& vec : inst.precuts) {
        seen.emplace(vec.begin(), vec.end());
    }
    seen.insert(sol.cuts.begin(), sol.cuts.end()); //Now RC cycles are also there.

    int n = inst.nbCity;

    // 1) Generate ALL MST cycles (your existing function)
    vector<vector<int>> allCycles = CheckAllMSTCycles(inst);
    cout << "# MST cycles generated = " << allCycles.size() << endl;

    // 2) makeUniqueSorted + dedupe + compute clique-average
    struct ScoredCycle {
        vector<int> verts; // canonical: sorted unique vertices
        double score;           // clique average distance
    };

    // Use hash set to dedupe cycles by vertex set
    unordered_set<vector<int>, VecHash, VecEq> uniq;
    uniq.reserve(allCycles.size() * 2);

    vector<ScoredCycle> scored;
    scored.reserve(allCycles.size());

    for (auto& cyc : allCycles)
    {
        // canonical form: sorted unique vertices
        sort(cyc.begin(), cyc.end());
        cyc.erase(unique(cyc.begin(), cyc.end()), cyc.end());

        // optional: keep length constraints here too (safety)
        int m = (int)cyc.size();
        if (m < 3 || m > n / 2) continue;

        if (!uniq.insert(cyc).second) continue;

        double avg = cliqueAverageDistance(cyc, inst.dist);
        scored.push_back({ move(cyc), avg });
    }

    // 3) Sort by smallest average distance first
    sort(scored.begin(), scored.end(),
        [](const ScoredCycle& a, const ScoredCycle& b) {
            return a.score < b.score;
        });

    // 4) Take first "needed" cycles not already in precuts (or previously added)
    set<set<int>> result;
    int countAdded = 0;

    for (const auto& sc : scored)
    {
        if (countAdded >= needed) break;

        set<int> cycleSet(sc.verts.begin(), sc.verts.end());

        if (seen.insert(cycleSet).second) {
            result.insert(move(cycleSet));
            ++countAdded;
        }
    }

    inst.nbMST_cycles = (int)result.size();

    return result;
}


// =================== MST CYCLES (CheckAllMSTCycles) ===================

vector<vector<int>> CheckAllMSTCycles(Instance& inst) {
    int n = inst.nbCity;


    vector<vector<int>> all_MST_cycles;
    all_MST_cycles.reserve(n * (n - 1) / 2 - (n - 1));

    // Step 1: Build MST
    vector<pair<int, int>> mstEdges = buildMST(inst.dist);

    // Build adjacency list for MST
    vector<vector<int>> adj(n);
    set<pair<int, int>> mstSet;
    for (auto e : mstEdges) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
        mstSet.insert({ min(e.first,e.second), max(e.first,e.second) });
    }

    // Collect non-MST edges
    vector<tuple<double, int, int>> nonMSTEdges;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (inst.dist[i][j] < 1e9) {
                pair<int, int> e = { i,j };
                if (mstSet.find(e) == mstSet.end()) {
                    nonMSTEdges.push_back(make_tuple(inst.dist[i][j], i, j));
                }
            }
        }
    }

    cout << "Non MST edges count = " << nonMSTEdges.size() << endl;

    // LCA prep
    int LOG = (int)ceil(log2((double)n)) + 1;
    vector<vector<int>> up(n, vector<int>(LOG, -1));
    vector<int> depth(n, 0);

    function<void(int, int)> dfs = [&](int u, int p) {
        up[u][0] = p;
        for (int i = 1; i < LOG; i++) {
            if (up[u][i - 1] != -1)
                up[u][i] = up[up[u][i - 1]][i - 1];
        }
        for (int v : adj[u]) {
            if (v == p) continue;
            depth[v] = depth[u] + 1;
            dfs(v, u);
        }
    };

    dfs(0, -1);

    //SORT THE EDGES SMALLEST TO LARGEST
    //sort(nonMSTEdges.begin(), nonMSTEdges.end());

    double threshold = (inst.lp_obj / inst.nbCity) * MULTIPLIER_MST_EDGE_INCLUSION;

    for (int k = 0; k < (int)nonMSTEdges.size(); k++) {
        double w;
        int u, v;
        tie(w, u, v) = nonMSTEdges[k];

        // threshold filter
        if (w > threshold) continue; // if it is too large, skip

        int L = getLCA(u, v, up, depth, LOG);

        int len = depth[u] + depth[v] - 2 * depth[L] + 1;
        if (len < 3 || len > n / 2) continue;

        vector<int> cycle;
        int x = u;
        while (x != L) {
            cycle.push_back(x);
            x = up[x][0];
        }
        cycle.push_back(L);
        vector<int> tail;
        x = v;
        while (x != L) {
            tail.push_back(x);
            x = up[x][0];
        }
        reverse(tail.begin(), tail.end());
        cycle.reserve(cycle.size() + tail.size());
        cycle.insert(cycle.end(), tail.begin(), tail.end());

        all_MST_cycles.push_back(move(cycle));
    }

    inst.nbMST_cycles = all_MST_cycles.size();

    return all_MST_cycles;
}

// =================== MERGE & DEDUPE CYCLES ===================

vector<vector<int>> mergeAndDedupe(const vector<vector<int>>& A, const vector<vector<int>>& B){
    using namespace std;

    unordered_set<vector<int>, VecHash, VecEq> seen;
    seen.reserve(A.size() + B.size());

    vector<vector<int>> result;
    result.reserve(A.size() + B.size());

    auto insertCycle = [&](const vector<int>& cyc)
    {
        vector<int> v = cyc;
        sort(v.begin(), v.end());
        v.erase(unique(v.begin(), v.end()), v.end());

        auto [it, inserted] = seen.insert(v);
        if (inserted)
            result.push_back(move(v));
    };

    for (const auto& cyc : A) insertCycle(cyc);
    for (const auto& cyc : B) insertCycle(cyc);

    return result;
}

set<set<int>> Check_K_RCCycles_Ordered_AvgCliqueLen(const int& needed, Instance& inst)
{
    // 0) Build "seen" from precuts once (same as MST version)
    set<set<int>> seen;
    for (const auto& vec : inst.precuts) {
        seen.emplace(vec.begin(), vec.end());
    }

    // 1) Generate ALL RC cycles (vertex sets / vertex lists)
    //    Find_all_RC_CYCLES already returns vertex sets (sorted unique) in most places,
    //   but I still makeUniqueSorted below jus to be safe.
    vector<vector<int>> allRCCycles = Find_all_RC_CYCLES(inst);

    // 2) makeUniqueSorted + remove duplicates + compute clique-average
    struct ScoredCycle {
        vector<int> verts; // canonical: sorted unique vertices
        double score;           // clique average distance
    };

    unordered_set<vector<int>, VecHash, VecEq> uniq;
    uniq.reserve(allRCCycles.size() * 2);

    vector<ScoredCycle> scored;
    scored.reserve(allRCCycles.size());

    const int n = inst.nbCity;

    for (auto& cyc : allRCCycles)
    {
        // canonical: sorted unique
        sort(cyc.begin(), cyc.end());
        cyc.erase(unique(cyc.begin(), cyc.end()), cyc.end());

        // optional safety constraints (match your general philosophy)
        const int m = (int)cyc.size();
        if (m < 3 || m > n / 2) continue;
        if (!uniq.insert(cyc).second) continue;

        double avg = cliqueAverageDistance(cyc, inst.dist);
        scored.push_back({ move(cyc), avg });
    }

    // 3) Sort by smallest clique-average first
    sort(scored.begin(), scored.end(), [](const ScoredCycle& a, const ScoredCycle& b) {
        if (a.score != b.score) return a.score < b.score;
            // tie-breaker: smaller cycle first (optional, but helps determinism)
            return a.verts.size() < b.verts.size();
    });

    // 4) Take first "needed" not already seen
    set<set<int>> result;
    int countAdded = 0;

    for (const auto& sc : scored)
    {
        if (countAdded >= needed) break;

        set<int> cycleSet(sc.verts.begin(), sc.verts.end());

        // avoid duplicates vs precuts + previously added
        if (seen.insert(cycleSet).second) {
            result.insert(move(cycleSet));
            ++countAdded;
        }
    }

    return result;
}

//FEATURES RELATED CONTENT
// ---------- helpers ----------
// Get undirected LP value x(i,j)
static double getXij(const vector<vector<double>>& x, int i, int j) {
    // assumes x stored only for i<j, symmetric
    if (i < j) return x[i][j];
    if (i > j) return x[j][i];
    return 0.0;
}
// Get undirected reduced-cost value rc(i,j)
static double getRCij(const vector<vector<double>>& rc, int i, int j) {
    if (i < j) return rc[i][j];
    if (i > j) return rc[j][i];
    return 0.0;
}

// MST total weight on subset vertices using Prim, O(m^2) where m=|S|
template <class WeightFn>
static double subsetMSTTotal_scratch_weighted(
    const vector<int>& S,
    WeightFn w, // double w(int gi, int gj)
    vector<double>& minEdge,
    vector<unsigned char>& used)
{
    const int m = (int)S.size();
    const double INF = 1e100;

    if ((int)minEdge.size() < m) minEdge.resize(m);
    if ((int)used.size() < m) used.resize(m);

    for (int i = 0; i < m; ++i) { minEdge[i] = INF; used[i] = 0; }

    minEdge[0] = 0.0;
    double total = 0.0;

    for (int it = 0; it < m; ++it) {
        int v = -1;
        for (int i = 0; i < m; ++i)
            if (!used[i] && (v == -1 || minEdge[i] < minEdge[v])) v = i;

        used[v] = 1;
        total += minEdge[v];

        const int gv = S[v];
        for (int u = 0; u < m; ++u) {
            if (used[u]) continue;
            const int gu = S[u];
            const double ww = w(gv, gu);
            if (ww < minEdge[u]) minEdge[u] = ww;
        }
    }
    return total;
}

struct SetFeatures {
    int m = 0;
    int is_found_by_MST = 0;
    int is_found_by_RC = 0;
    double s_frac = 0.0;

    double dist_in_avg_norm = 0.0;
    double dist_2min_out_avg_norm = 0.0;
    double dist_2min_in_avg_norm = 0.0;

    // OLD
    // double mst_avg = 0.0;

    //  3 MST variants
    double mst_dist_avg_norm = 0.0;   // MST on dist / dist_scale
    double mst_1mx_avg = 0.0;         // MST on (1-x) / x_scale (where x_scale = 1since already in a bounded region)
    double mst_rc_avg_norm = 0.0;     // MST on rc / rc_scale

    double rc_in_avg_norm = 0.0;
    double rc_2min_out_avg_norm = 0.0;
    double rc_2min_in_avg_norm = 0.0;

    double lp_in_avg = 0.0;
    double lp_2max_out_avg = 0.0;
    double lp_2max_in_avg = 0.0;

    double greedy_cycle_dist_avg_norm = 0.0;
    double greedy_cycle_1mx_avg = 0.0;
    double greedy_cycle_rc_avg_norm = 0.0;

    //  (1) avg fraction of KNN outside S
    double dist_knn_out_frac_avg = 0.0;
    double lp_knn_out_frac_avg = 0.0;
    double rc_knn_out_frac_avg = 0.0;

    //  (2) avg first internal neighbor value (scaled where needed)
    double dist_first_in_avg_norm = 0.0;
    double lp_first_in_avg = 0.0;
    double rc_first_in_avg_norm = 0.0;

    //  (3) boundary gap (scaled where needed)
    double dist_boundary_gap_avg_norm = 0.0; // (best_out - best_in)/dist_scale
    double lp_boundary_gap_avg = 0.0; // (best_in_x - best_out_x)
    double rc_boundary_gap_avg_norm = 0.0; // (best_out - best_in)/rc_scale

};

// ---------- Greedy cycle helpers ------

// choose next in S (unvisited) minimizing dist
static  int pickNextInS_minDist(
    int cur,
    const vector<int>& S,
    const MarkSet& inS,
    FeatureScratch& scratch,
    const vector<vector<double>>& dist,
    const vector<int>& nnDist,
    int n)
{
    // fast path: scan Top-K neighbors
    for (int t = 0; t < (int)nnDist.size(); ++t) {
        int j = nnDist[t];
        if (j == cur) continue;
        if (!inS.has(j)) continue;
        if (scratch.isVisited(j)) continue;
        return j;
    }

    // fallback: scan all vertices in S
    int best = -1;
    double bestVal = 9999999999999999999999999999999999999999; 
    for (int j : S) {
        if (j == cur) continue;
        if (scratch.isVisited(j)) continue;
        double v = dist[cur][j];
        if (v < bestVal) { bestVal = v; best = j; }
    }
    return best;
}

// choose next in S (unvisited) minimizing rc
static  int pickNextInS_minRC(
    int cur,
    const vector<int>& S,
    const MarkSet& inS,
    FeatureScratch& scratch,
    const vector<vector<double>>& rc,
    const vector<int>& nnRC,
    int n)
{
    for (int t = 0; t < (int)nnRC.size(); ++t) {
        int j = nnRC[t];
        if (j == cur) continue;
        if (!inS.has(j)) continue;
        if (scratch.isVisited(j)) continue;
        return j;
    }

    int best = -1;
    double bestVal = 9999999999999999999999999999999999999999; 
    for (int j : S) {
        if (j == cur) continue;
        if (scratch.isVisited(j)) continue;
        double v = getRCfull(rc, cur, j);
        if (v < bestVal) { bestVal = v; best = j; }
    }
    return best;
}

// choose next in S (unvisited) minimizing (1 - x) == maximizing x
static int pickNextInS_maxX(
    int cur,
    const vector<int>& S,
    const MarkSet& inS,
    FeatureScratch& scratch,
    const vector<vector<double>>& x,
    const vector<int>& nnLP, // already sorted by x desc in your TopKLists
    int n)
{
    for (int t = 0; t < (int)nnLP.size(); ++t) {
        int j = nnLP[t];
        if (j == cur) continue;
        if (!inS.has(j)) continue;
        if (scratch.isVisited(j)) continue;
        return j;
    }

    int best = -1;
    double bestVal = -9999999999999999999999999999999999999999; 
    for (int j : S) {
        if (j == cur) continue;
        if (scratch.isVisited(j)) continue;
        double v = getXfull(x, cur, j);
        if (v > bestVal) { bestVal = v; best = j; }
    }
    return best;
}

static SetFeatures computeFeaturesForSet_fast(
    const vector<int>& S,   // sorted unique
    int n,
    const vector<vector<double>>& dist,
    const vector<vector<double>>& x,
    const vector<vector<double>>& rc,
    double dist_scale,
    double rc_scale,
    const TopKLists& topk,
    FeatureScratch& scratch,
    MarkSet& ms)
{
    SetFeatures f;
    const int m = (int)S.size();
    f.m = m;
    f.s_frac = (double)m / (double)n;


    // -----------------------------
    // membership mark (timestamp)
    // ---------------------------
    {
        ms.next();
        for (int v : S) ms.add(v);
    }

    // --------------------------
    // Pairwise internal averages
    // -----------------------------
    {

        long long pairs = 0;
        long double sumDist = 0.0L, sumRC = 0.0L, sumLP = 0.0L;

        for (int a = 0; a < m; ++a) {
            const int i = S[a];
            for (int b = a + 1; b < m; ++b) {
                const int j = S[b];
                sumDist += (long double)dist[i][j];
                sumRC += (long double)getRCfull(rc, i, j);
                sumLP += (long double)getXfull(x, i, j);
                ++pairs;
            }
        }

        const double invPairs = 1.0 / (double)pairs;
        f.dist_in_avg_norm = (double)(sumDist * invPairs) / dist_scale;
        f.rc_in_avg_norm = (double)(sumRC * invPairs) / rc_scale;
        f.lp_in_avg = (double)(sumLP * invPairs);
    }

    // -----------------------------
    // Per-vertex 2-min / 2-max (Top-K for OUT and IN)
    // -----------------------------
    long double sum_dist2min_out = 0.0L, sum_dist2min_in = 0.0L;
    long double sum_rc2min_out = 0.0L, sum_rc2min_in = 0.0L;
    long double sum_lp2max_out = 0.0L, sum_lp2max_in = 0.0L;
    {
        for (int i : S) {
            // OUT
            double d1o, d2o, r1o, r2o, x1o, x2o;
            twoMin_OUT_fromTopK_dist(i, ms, n, dist, topk.distNN[i], d1o, d2o);
            twoMin_OUT_fromTopK_rc(i, ms, n, rc, topk.rcNN[i], r1o, r2o);
            twoMax_OUT_fromTopK_lp(i, ms, n, x, topk.lpNN[i], x1o, x2o);

            // IN
            double d1i, d2i, r1i, r2i, x1i, x2i;
            twoMin_IN_fromTopK_dist(i, ms, n, dist, topk.distNN[i], d1i, d2i);
            twoMin_IN_fromTopK_rc(i, ms, n, rc, topk.rcNN[i], r1i, r2i);
            twoMax_IN_fromTopK_lp(i, ms, n, x, topk.lpNN[i], x1i, x2i);

            sum_dist2min_out += (long double)avg2(d1o, d2o);
            sum_dist2min_in += (long double)avg2(d1i, d2i);

            sum_rc2min_out += (long double)avg2(r1o, r2o);
            sum_rc2min_in += (long double)avg2(r1i, r2i);

            sum_lp2max_out += (long double)avg2(x1o, x2o);
            sum_lp2max_in += (long double)avg2(x1i, x2i);
        }
    }

    const double invM = 1.0 / (double)m;

    f.dist_2min_out_avg_norm = (double)(sum_dist2min_out * invM) / dist_scale;
    f.dist_2min_in_avg_norm = (double)(sum_dist2min_in * invM) / dist_scale;

    f.rc_2min_out_avg_norm = (double)(sum_rc2min_out * invM) / rc_scale;
    f.rc_2min_in_avg_norm = (double)(sum_rc2min_in * invM) / rc_scale;

    f.lp_2max_out_avg = (double)(sum_lp2max_out * invM);
    f.lp_2max_in_avg = (double)(sum_lp2max_in * invM);

    // (1) Avg fraction of KNN outside S (Dist/LP/RC)
    {
        long double sumDist = 0.0L, sumLP = 0.0L, sumRC = 0.0L;

        const int Kd = topk.K; // each list has size K (or <=K if n small)
        const double invK = (Kd > 0) ? (1.0 / (double)Kd) : 0.0;

        for (int i : S) {
            int outD = 0, outL = 0, outR = 0;

            // Dist
            for (int j : topk.distNN[i]) {
                if (!ms.has(j)) outD++;
            }
            // LP
            for (int j : topk.lpNN[i]) {
                if (!ms.has(j)) outL++;
            }
            // RC
            for (int j : topk.rcNN[i]) {
                if (!ms.has(j)) outR++;
            }

            sumDist += (long double)outD * invK;
            sumLP += (long double)outL * invK;
            sumRC += (long double)outR * invK;
        }

        f.dist_knn_out_frac_avg = (double)(sumDist * invM);
        f.lp_knn_out_frac_avg = (double)(sumLP * invM);
        f.rc_knn_out_frac_avg = (double)(sumRC * invM);
    }

    // (2) Avg  first internal neighbor  value (Dist/LP/RC)
    {

        if (m <= 1) {
            f.dist_first_in_avg_norm = 0.0;
            f.lp_first_in_avg = 0.0;
            f.rc_first_in_avg_norm = 0.0;
        }
        else {
            long double sumBestDist = 0.0L, sumBestLP = 0.0L, sumBestRC = 0.0L;

            for (int i : S) {
                double bd = first_IN_dist_fromTopK(i, ms, S, dist, topk.distNN[i]);
                double bx = first_IN_lp_fromTopK(i, ms, S, x, topk.lpNN[i]);
                double br = first_IN_rc_fromTopK(i, ms, S, rc, topk.rcNN[i]);

                sumBestDist += (long double)bd;
                sumBestLP += (long double)bx;
                sumBestRC += (long double)br;
            }

            f.dist_first_in_avg_norm = (double)(sumBestDist * invM) / dist_scale;
            f.lp_first_in_avg = (double)(sumBestLP * invM);
            f.rc_first_in_avg_norm = (double)(sumBestRC * invM) / rc_scale;
        }
    }
    // (3) Boundary gap (Dist/LP/RC)
    {

        if (m <= 1) {
            f.dist_boundary_gap_avg_norm = 0.0;
            f.lp_boundary_gap_avg = 0.0;
            f.rc_boundary_gap_avg_norm = 0.0;
        }
        else {
            long double sumGapDist = 0.0L, sumGapLP = 0.0L, sumGapRC = 0.0L;

            for (int i : S) {
                double inD = first_IN_dist_fromTopK(i, ms, S, dist, topk.distNN[i]);
                double outD = first_OUT_dist_fromTopK(i, ms, n, dist, topk.distNN[i]);

                double inR = first_IN_rc_fromTopK(i, ms, S, rc, topk.rcNN[i]);
                double outR = first_OUT_rc_fromTopK(i, ms, n, rc, topk.rcNN[i]);

                double inX = first_IN_lp_fromTopK(i, ms, S, x, topk.lpNN[i]);
                double outX = first_OUT_lp_fromTopK(i, ms, n, x, topk.lpNN[i]);

                sumGapDist += (long double)(outD - inD);
                sumGapRC += (long double)(outR - inR);
                sumGapLP += (long double)(inX - outX);
            }

            f.dist_boundary_gap_avg_norm = (double)(sumGapDist * invM) / dist_scale;
            f.rc_boundary_gap_avg_norm = (double)(sumGapRC * invM) / rc_scale;
            f.lp_boundary_gap_avg = (double)(sumGapLP * invM);
        }
    }

    // -----------------------------
    // MST features: 3 weight definitions
    // -----------------------------
    {
        const double denom = (double)(m - 1);

        // 1) MST on original distances
        {
            double mst = subsetMSTTotal_scratch_weighted(S,[&](int a, int b) { return dist[a][b]; }, scratch.minEdge, scratch.used);
            f.mst_dist_avg_norm = (mst / denom) / dist_scale;
        }

        // 2) MST on (1 - x_ij)
        {
            double mst = subsetMSTTotal_scratch_weighted(
                S,
                [&](int a, int b) { return 1.0 - getXfull(x, a, b); },
                scratch.minEdge, scratch.used
            );
            // If you want this unnormalized like greedy_cycle_1mx_avg, set x_scale = 1.
            // Otherwise define x_scale similarly to dist_scale.
            constexpr double x_scale = 1.0;
            f.mst_1mx_avg = (mst / denom) / x_scale;
        }

        // 3) MST on reduced costs
        {
            double mst = subsetMSTTotal_scratch_weighted(S, [&](int a, int b) { return getRCfull(rc, a, b); }, scratch.minEdge, scratch.used);
            f.mst_rc_avg_norm = (mst / denom) / rc_scale;
        }
    }

    // -----------------------------
    //  Greedy cycle features on S
    // -----------------------------
    {
        // For tiny sets, define as 0 (or you can skip writing these rows)

        // 1) Greedy by DIST
        scratch.nextVisitedStamp();
        int start = S[0];   //  start (S is sorted unique)
        int cur = start;
        scratch.markVisited(cur);

        long double sumDist = 0.0L;

        for (int step = 1; step < m; ++step) {
            int nxt = pickNextInS_minDist(
                cur, S, ms, scratch,
                dist, topk.distNN[cur], n
            );
            if (nxt < 0) break; // should not happen
            sumDist += (long double)dist[cur][nxt];
            cur = nxt;
            scratch.markVisited(cur);
        }
        // close the cycle
        sumDist += (long double)dist[cur][start];

        f.greedy_cycle_dist_avg_norm = (double)(sumDist / (long double)m) / dist_scale;


        // 2) Greedy by (1 - x) ---  maximize x
        scratch.nextVisitedStamp();
        start = S[0];
        cur = start;
        scratch.markVisited(cur);

        long double sum1mx = 0.0L;

        for (int step = 1; step < m; ++step) {
            int nxt = pickNextInS_maxX(
                cur, S, ms, scratch,
                x, topk.lpNN[cur], n
            );
            if (nxt < 0) break;
            double xij = getXfull(x, cur, nxt);
            sum1mx += (long double)(1.0 - xij);
            cur = nxt;
            scratch.markVisited(cur);
        }
        sum1mx += (long double)(1.0 - getXfull(x, cur, start));

        f.greedy_cycle_1mx_avg = (double)(sum1mx / (long double)m);


        // 3) Greedy by RC (minimize rc)
        scratch.nextVisitedStamp();
        start = S[0];
        cur = start;
        scratch.markVisited(cur);

        long double sumRC = 0.0L;

        for (int step = 1; step < m; ++step) {
            int nxt = pickNextInS_minRC(
                cur, S, ms, scratch,
                rc, topk.rcNN[cur], n
            );
            if (nxt < 0) break;
            sumRC += (long double)getRCfull(rc, cur, nxt);
            cur = nxt;
            scratch.markVisited(cur);
        }
        sumRC += (long double)getRCfull(rc, cur, start);

        f.greedy_cycle_rc_avg_norm = (double)(sumRC / (long double)m) / rc_scale;
    }

    return f;
}
// turn sorted unique vertex list into a quoted CSV field
static string vertsToCSVField(const vector<int>& v) {
    string s;
    s.reserve(v.size() * 6);
    s.push_back('"');
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) s.push_back(';');
        s += to_string(v[i]);
    }
    s.push_back('"');
    return s;
}

vector<vector<int>> SelectPrecutsWithNN(
    const vector<vector<int>>& sets,
    const vector<vector<int>>& mstCycles,
    const vector<vector<int>>& rcCycles,
    const vector<vector<double>>& dist,
    const vector<vector<double>>& lp_x,
    const vector<vector<double>>& lp_rc,
    const Instance& inst,
    const string& modelPath,
    double threshold)
{
    using namespace std;

    const int n = (int)dist.size();
    vector<vector<int>> selected;

    if (sets.empty()) {
        return selected;
    }

    FeatureScratch scratch(n);
    MarkSet ms(n);

    const int K = min(
        n - 1,
        max(8, (int)ceil(sqrt((double)n)))
    );

    TopKLists topk = buildTopKLists(dist, lp_x, lp_rc, K);

    double dist_scale = inst.lp_obj / n;
    double rc_scale = dist_scale;

    unordered_set<vector<int>, VecHash, VecEq> mstSeen;
    unordered_set<vector<int>, VecHash, VecEq> rcSeen;

    mstSeen.reserve(mstCycles.size());
    rcSeen.reserve(rcCycles.size());

    auto makeUniqueSorted = [&](vector<int> v) {
        sort(v.begin(), v.end());
        v.erase(unique(v.begin(), v.end()), v.end());
        return v;
    };

    for (const auto& cyc : mstCycles) mstSeen.insert(makeUniqueSorted(cyc));
    for (const auto& cyc : rcCycles)  rcSeen.insert(makeUniqueSorted(cyc));

    const string base = inst.saveFolder + inst.name + "_nn_inference";
    const string featurePath = base + "_features.csv";
    const string predPath = base + "_predictions.txt";

    ofstream out(featurePath);
    if (!out.is_open()) {
        throw runtime_error("Could not open feature temp file: " + featurePath);
    }

    // IMPORTANT:
    // Python should only receive model features, not id/m/verts/label.
    out << "is_found_by_MST,is_found_by_RC,"
        << "s_frac,"
        << "dist_in_avg_norm,dist_2min_out_avg_norm,dist_2min_in_avg_norm,"
        << "rc_in_avg_norm,rc_2min_out_avg_norm,rc_2min_in_avg_norm,"
        << "lp_in_avg,lp_2max_out_avg,lp_2max_in_avg,"
        << "dist_knn_out_frac_avg,lp_knn_out_frac_avg,rc_knn_out_frac_avg,"
        << "dist_first_in_avg_norm,lp_first_in_avg,rc_first_in_avg_norm,"
        << "dist_boundary_gap_avg_norm,lp_boundary_gap_avg,rc_boundary_gap_avg_norm,"
        << "mst_dist_avg_norm,mst_1mx_avg,mst_rc_avg_norm,"
        << "greedy_cycle_dist_avg_norm,greedy_cycle_1mx_avg,greedy_cycle_rc_avg_norm\n";

    vector<vector<int>> canonicalSets;
    canonicalSets.reserve(sets.size());


    for (const auto& rawS : sets) {
        vector<int> S = rawS;
        sort(S.begin(), S.end());
        S.erase(unique(S.begin(), S.end()), S.end());

        canonicalSets.push_back(S);

        SetFeatures f = computeFeaturesForSet_fast(
            S, n, dist, lp_x, lp_rc,
            dist_scale, rc_scale,
            topk, scratch, ms
        );

        f.is_found_by_MST = mstSeen.count(S) ? 1 : 0;
        f.is_found_by_RC = rcSeen.count(S) ? 1 : 0;

        out
            << f.is_found_by_MST << ","
            << f.is_found_by_RC << ","
            << f.s_frac << ","
            << f.dist_in_avg_norm << "," << f.dist_2min_out_avg_norm << "," << f.dist_2min_in_avg_norm << ","
            << f.rc_in_avg_norm << "," << f.rc_2min_out_avg_norm << "," << f.rc_2min_in_avg_norm << ","
            << f.lp_in_avg << "," << f.lp_2max_out_avg << "," << f.lp_2max_in_avg << ","
            << f.dist_knn_out_frac_avg << "," << f.lp_knn_out_frac_avg << "," << f.rc_knn_out_frac_avg << ","
            << f.dist_first_in_avg_norm << "," << f.lp_first_in_avg << "," << f.rc_first_in_avg_norm << ","
            << f.dist_boundary_gap_avg_norm << "," << f.lp_boundary_gap_avg << "," << f.rc_boundary_gap_avg_norm << ","
            << f.mst_dist_avg_norm << "," << f.mst_1mx_avg << "," << f.mst_rc_avg_norm << ","
            << f.greedy_cycle_dist_avg_norm << "," << f.greedy_cycle_1mx_avg << "," << f.greedy_cycle_rc_avg_norm
            << "\n";
    }

    out.close();

    const string pyScript = "/home/fakcay/TSP/ML_PART/NN_Keras/TrainOnce_AllData/trained_model/predict_sets.py"; //direct it to the python script

    string cmd = "python3 " + shellQuote(pyScript) + " " + shellQuote(modelPath) + " "+ shellQuote(featurePath) + " " + shellQuote(predPath);

    int rc = system(cmd.c_str());
    if (rc != 0) {
        throw runtime_error("Python prediction command failed with code " + to_string(rc));
    }

    vector<double> preds = readPredictionsFromFile(predPath);
    if (preds.size() != canonicalSets.size()) {
        throw runtime_error(
            "Prediction count mismatch. preds=" + to_string(preds.size()) +
            ", sets=" + to_string(canonicalSets.size())
        );
    }

    selected.reserve(canonicalSets.size());

    for (size_t i = 0; i < canonicalSets.size(); ++i) {
        if (preds[i] >= threshold) {
            if (0) { //printing
                cout << "Selected set (pred=" << preds[i] << "): ";
                for (int v : canonicalSets[i])
                    cout << v << " ";
                cout << "\n";
            }
            selected.push_back(canonicalSets[i]);
        }
    }

    cout << "NN inference done. Threshold = " << threshold  << ", selected = " << selected.size()  << " / " << canonicalSets.size() << "\n";

    return selected;
}

