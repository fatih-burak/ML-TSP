#include "cycleGenHelper_functions.h"

using std::vector;
using std::pair;
using std::cout;
using std::cerr;
using std::endl;
// =================== TESTING CONFIGURATION =================================
constexpr int NOT_REDUCED_AND_SIMPLE = 0;

constexpr int XOR2_ACTIVE = 1;
constexpr int XOR3_ACTIVE = 1;
constexpr int XOR4_ACTIVE = 0;
constexpr int XOR5_ACTIVE = 0;

constexpr double RC_LIM = 0.0;
constexpr double MULTIPLIER_MST_EDGE_INCLUSION = 10.0;

// =================== INTERNAL HELPERS (only in this .cpp) ===================
// hashing... //

inline long long encodeEdge(int u, int v) {
    if (u > v) std::swap(u, v);
    return ((long long)u << 32) | (unsigned int)v;
}

// Hash & equality for bitset<BSZ> (for unordered_set)
struct BitsetHash {
    std::size_t operator()(const std::bitset<BSZ>& b) const noexcept {
        constexpr size_t WORDS = (BSZ + 63) / 64; // check how many 64-bit blocks do we need to hold all the bits?
        const uint64_t* ptr = reinterpret_cast<const uint64_t*>(&b); // ok now treat the memory of this bitset as if it were an array of 64-bit integers.
        std::size_t h = 0; //initiate the running hash value by starting from 0

        for (size_t i = 0; i < WORDS; ++i) { // go through every 64-bit block in the bitset.
            uint64_t word = ptr[i];
            // standard hash-combine (take the current block and mix it into the running hash.)
            h ^= std::hash<uint64_t>()(word) + 0x9e3779b97f4a7c15ULL + (h << 6)+ (h >> 2); // 0x9e3779b97f4a7c15ULL is used (classic in hashing to spread valuess)
        }
        return h;
    }
};

struct BitsetEq {
    bool operator()(const std::bitset<BSZ>& a,
        const std::bitset<BSZ>& b) const {return a == b;}
};

struct VecHash {
    std::size_t operator()(const std::vector<int>& v) const noexcept {
        std::size_t h = 0; // running hash value

        for (int x : v) { // go through each integer in the vector
            h ^= std::hash<int>()(x) + 0x9e3779b97f4a7c15ULL + (h << 6)+ (h >> 2); // hash the current int, mixing constant, mix in previous hash shifted left and right
        }
        return h; // final hash for the whole vector
    }
};

struct VecEq {
    bool operator()(const std::vector<int>& a,
        const std::vector<int>& b) const noexcept {return a == b;}
};

// =================== LP RELAXATION ===================
std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>>
SolveLPRelaxation(Instance& inst)
{
    try {
        GRBEnv env = GRBEnv();
        env.set(GRB_IntParam_OutputFlag, 1);
        GRBModel model = GRBModel(env);

        int n = inst.nbCity;
        std::vector<std::vector<GRBVar>> x(n, std::vector<GRBVar>(n));

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

        std::vector<std::vector<int>> skeleton, remaining_edges;

        // ---- temp storage for pruning ----
        struct SEdge {
            int i, j;
            double rc;
            double dist;
        };
        std::vector<SEdge> sk; sk.reserve((size_t)n * (n - 1) / 2);

        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
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
                std::sort(cand.begin(), cand.end(), [&](int a, int b) {
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
                // 
                // rebuild sk keeping non-removed
                std::vector<SEdge> sk2;
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
        std::cerr << "\n=== Gurobi Exception in SolveLPRelaxation ===\n";
        std::cerr << "Error code: " << e.getErrorCode() << "\n";
        std::cerr << "Message:    " << e.getMessage() << "\n";
        std::cerr << "============================================\n";
        throw;
    }
    catch (std::exception& e) {
        std::cerr << "\n=== std::exception in SolveLPRelaxation ===\n";
        std::cerr << e.what() << "\n";
        std::cerr << "===========================================\n";
        throw;
    }
    catch (...) {
        std::cerr << "\n=== Unknown Exception in SolveLPRelaxation ===\n";
        throw;
    }
}


// =================== SPANNING FOREST ===================

// Build spanning forest from edge list
std::vector<std::pair<int, int>>
buildSpanningForest(const std::vector<std::vector<int>>& skeleton,
    int n,
    std::vector<int>& component)
{
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

bool isSimpleCycle_fast(const std::bitset<BSZ>& b,
    const std::vector<int>& U,
    const std::vector<int>& V,
    const std::vector<std::vector<int>>& inc)
{
    int E = (int)U.size();
    int Vn = (int)inc.size();

    std::vector<int>  degree(Vn, 0);
    std::vector<char> used(Vn, 0);

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
    std::vector<char> visited(Vn, 0);
    std::stack<int> st;
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

// =================== BCC (EDGE-BICONNECTED COMPONENTS) ===================

BCCResult tarjanBCC(const std::vector<std::vector<int>>& adj, int n)
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
                    std::unordered_set<int> verts;
                    while (!st.empty()) {
                        auto e = st.back(); st.pop_back();
                        edges.push_back(e);
                        verts.insert(e.first);
                        verts.insert(e.second);
                        if (e.first == u && e.second == v) break;
                    }
                    out.blocks.push_back(std::move(edges));
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
//NOTE THAT we ELIMINATE ALREADY FUNDAMENTAL CYCLES THAT HAS SIZE BIGGER THAN N/2

// ---------- DSU (Kruskal) ----------
struct DSU {
    std::vector<int> p, r;
    explicit DSU(int n) : p(n), r(n, 0) { std::iota(p.begin(), p.end(), 0); }

    int find(int a) { return p[a] == a ? a : p[a] = find(p[a]); }

    bool unite(int a, int b) {
        a = find(a); b = find(b);
        if (a == b) return false;
        if (r[a] < r[b]) std::swap(a, b);
        p[b] = a;
        if (r[a] == r[b]) r[a]++;
        return true;
    }
};

static inline long long encEdgeLL(int u, int v) {
    if (u > v) std::swap(u, v);
    return ((long long)u << 32) | (unsigned int)v;
}

// ---------- MST fundamental cycles per block ----------
vector<std::vector<int>> computeFundamentalCycles_Block_MST(
    const std::vector<std::pair<int, int>>& blockEdges,     // GLOBAL vertex ids
    const std::vector<int>& blockVerts,     // GLOBAL vertex ids
    const std::vector<std::vector<double>>& dist,           // inst.dist
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
    int LOG = (int)std::ceil(std::log2((double)V)) + 2;
    vector<vector<int>> up(V, vector<int>(LOG, -1));
    for (int i = 0; i < V; i++) up[i][0] = parent[i];

    for (int k = 1; k < LOG; k++) {
        for (int i = 0; i < V; i++) {
            int mid = up[i][k - 1];
            up[i][k] = (mid == -1 ? -1 : up[mid][k - 1]);
        }
    }

    auto LCA = [&](int a, int b) {
        if (depth[a] < depth[b]) std::swap(a, b);
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
        if ((int)cyc.size() <= maxEdgesPerCycle) {
            cycles.push_back(std::move(cyc));
        }
    }

    return cycles;
}



vector<std::vector<int>> computeFundamentalCycles_Block(
    const vector<std::pair<int, int>>& blockEdges,
    const vector<int>& blockVerts,
    const int& maxEdgesPerCycle)
{
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
    int LOG = (int)std::ceil(std::log2((double)V)) + 2;
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
        if (depth[a] < depth[b]) std::swap(a, b);
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
        if (u > v) std::swap(u, v);
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

        if (NOT_REDUCED_AND_SIMPLE) cycles.push_back(std::move(cyc)); // Include all len fundamental cycles.
        else {
            if (cyc.size() <= maxEdgesPerCycle) cycles.push_back(std::move(cyc)); // Include only fundamental cycles that has len<=n/2
        }
    }

    return cycles;
}

// =================== ENCODE CYCLES AS BITSETS PER BLOCK ===================

vector<std::bitset<BSZ>> encodeCycles_Block(
    const vector<std::vector<int>>& cycles,
    const vector<std::pair<int, int>>& blockEdges)
{
    using namespace std;

    unordered_map<long long, int> edgeMap;
    edgeMap.reserve(blockEdges.size() * 2);

    for (int i = 0; i < (int)blockEdges.size(); ++i) {
        int u = blockEdges[i].first;
        int v = blockEdges[i].second;
        if (u > v) std::swap(u, v);
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
            if (u > v) std::swap(u, v);
            long long key = ((long long)u << 32) | (unsigned int)v;
            auto it = edgeMap.find(key);
            if (it != edgeMap.end()) b.set(it->second);
        }

        if (b.any()) out.push_back(b);
    }

    return out;
}

// =================== XOR-2 / XOR-3 / XOR-4 PER BLOCK ===================

std::vector<std::bitset<BSZ>> generateXORCycles_Block_old(
    const vector<std::bitset<BSZ>>& fundBits,const vector<int>& U_local,const vector<int>& V_local,const vector<std::vector<int>>& incEdges,int maxEdgesPerCycle)
{
    using namespace std;

    unordered_set<bitset<BSZ>, BitsetHash, BitsetEq> seen;

    int F = (int)fundBits.size();
    if (F == 0) return {};

    if (XOR2_ACTIVE) {
        // XOR-2
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
        // XOR-3 with timing breakdown
        using Clock = std::chrono::steady_clock;
        using ns = std::chrono::nanoseconds;

        double t_overlap = 0.0;
        double t_xor = 0.0;
        double t_count = 0.0;
        double t_simple = 0.0;

        for (int i = 0; i < F; ++i) {
            for (int j = i + 1; j < F; ++j) {
                for (int k = j + 1; k < F; ++k) {

                    auto t0 = Clock::now();

                    int Ai = fundBits[i].count();
                    int Aj = fundBits[j].count();
                    int Ak = fundBits[k].count();

                    int o_ij = (fundBits[i] & fundBits[j]).count();
                    int o_ik = (fundBits[i] & fundBits[k]).count();
                    int o_jk = (fundBits[j] & fundBits[k]).count();
                    //int o_ijk = (fundBits[i] & fundBits[j] & fundBits[k]).count();
                    //int LB = Ai + Aj + Ak - 2 * (o_ij + o_ik + o_jk) + o_ijk;

                    int LB = Ai + Aj + Ak - 2 * (o_ij + o_ik + o_jk);
                    auto t1 = Clock::now();

                    if (LB > maxEdgesPerCycle)
                        continue;

                    t_overlap += std::chrono::duration_cast<ns>(t1 - t0).count();

                    t0 = Clock::now();
                    bitset<BSZ> x = fundBits[i] ^ fundBits[j] ^ fundBits[k];
                    t1 = Clock::now();
                    t_xor += std::chrono::duration_cast<ns>(t1 - t0).count();

                    t0 = Clock::now();
                    int cnt = x.count();
                    t1 = Clock::now();
                    t_count += std::chrono::duration_cast<ns>(t1 - t0).count();

                    if (cnt == 0 || cnt > maxEdgesPerCycle)
                        continue;

                    t0 = Clock::now();
                    if (isSimpleCycle_fast(x, U_local, V_local, incEdges))
                        seen.insert(x);
                    t1 = Clock::now();
                    t_simple += std::chrono::duration_cast<ns>(t1 - t0).count();
                }
            }
        }

        std::cout << "\n==== Timing breakdown for XOR-3 generation ====\n";
        std::cout << "Time in overlaps (Ai/Aj/Ak + o_ij/o_ik/o_jk): "
            << t_overlap / 1e9 << " sec\n";
        std::cout << "Time in XOR operations: "
            << t_xor / 1e9 << " sec\n";
        std::cout << "Time in count() on XOR: "
            << t_count / 1e9 << " sec\n";
        std::cout << "Time in isSimpleCycle_fast(): "
            << t_simple / 1e9 << " sec\n";
        std::cout << "Total measured time: "
            << (t_overlap + t_xor + t_count + t_simple) / 1e9 << " sec\n";
        std::cout << "=============================================\n";
    }
    if (XOR4_ACTIVE) {
        // ===============================================================
        // XOR-4 timing breakdown
        // ===============================================================
        using Clock = std::chrono::steady_clock;
        using ns = std::chrono::nanoseconds;

        double t4_overlap = 0.0;
        double t4_xor = 0.0;
        double t4_count = 0.0;
        double t4_simple = 0.0;

        int xor4Inserted = 0;

        for (int i = 0; i < F; ++i) {
            for (int j = i + 1; j < F; ++j) {
                for (int k = j + 1; k < F; ++k) {
                    for (int l = k + 1; l < F; ++l) {

                        // -----------------------------
                        // 1) Overlap-based lower bound
                        // -----------------------------
                        auto t0 = Clock::now();

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

                        // Lower bound without 4-way overlap (good enough)
                        int LB = Ai + Aj + Ak + Al
                            - 2 * (o_ij + o_ik + o_il + o_jk + o_jl + o_kl);

                        auto t1 = Clock::now();
                        t4_overlap += std::chrono::duration_cast<ns>(t1 - t0).count();

                        if (LB > maxEdgesPerCycle)
                            continue;

                        // -----------------------------
                        // 2) XOR-4
                        // -----------------------------
                        t0 = Clock::now();
                        bitset<BSZ> x = fundBits[i] ^ fundBits[j]
                            ^ fundBits[k] ^ fundBits[l];
                        t1 = Clock::now();
                        t4_xor += std::chrono::duration_cast<ns>(t1 - t0).count();

                        // -----------------------------
                        // 3) Count bits
                        // -----------------------------
                        t0 = Clock::now();
                        int cnt = x.count();
                        t1 = Clock::now();
                        t4_count += std::chrono::duration_cast<ns>(t1 - t0).count();

                        if (cnt == 0 || cnt > maxEdgesPerCycle)
                            continue;

                        // -----------------------------
                        // 4) Simple cycle test
                        // -----------------------------
                        t0 = Clock::now();
                        if (isSimpleCycle_fast(x, U_local, V_local, incEdges)) {
                            auto [it, ins] = seen.insert(x);
                            if (ins) xor4Inserted++;
                        }
                        t1 = Clock::now();
                        t4_simple += std::chrono::duration_cast<ns>(t1 - t0).count();
                    }
                }
            }
        }

        std::cout << "\n==== Timing breakdown for XOR-4 generation ====\n";
        std::cout << "Time in overlaps (Ai/Aj/Ak/Al + pair overlaps): "
            << t4_overlap / 1e9 << " sec\n";
        std::cout << "Time in XOR operations: "
            << t4_xor / 1e9 << " sec\n";
        std::cout << "Time in count(): "
            << t4_count / 1e9 << " sec\n";
        std::cout << "Time in isSimpleCycle_fast(): "
            << t4_simple / 1e9 << " sec\n";
        std::cout << "Total measured time: "
            << (t4_overlap + t4_xor + t4_count + t4_simple) / 1e9
            << " sec\n";
        std::cout << "Valid XOR-4 cycles inserted: " << xor4Inserted << "\n";
        std::cout << "=============================================\n";
    }
    return vector<bitset<BSZ>>(seen.begin(), seen.end());
}

vector<std::bitset<BSZ>> generateXORCycles_Block(
    const vector<std::bitset<BSZ>>& fundBits, const vector<int>& U_local, const vector<int>& V_local, const vector<std::vector<int>>& incEdges, int maxEdgesPerCycle)
{
    using namespace std;

    unordered_set<bitset<BSZ>, BitsetHash, BitsetEq> seen;

    const int F = (int)fundBits.size();
    if (F == 0) return {};

    vector<bitset<BSZ>> xor2Simple;
    vector<bitset<BSZ>> xor3Simple;
    vector<bitset<BSZ>> xor4Simple; // NEW (needed for XOR-5)

    if (XOR2_ACTIVE) {
        auto t2_start = chrono::steady_clock::now();

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

        auto t2_end = chrono::steady_clock::now();
        double t_xor2 = chrono::duration_cast<chrono::nanoseconds>(t2_end - t2_start).count() / 1e9;
        cout << "  XOR-2 done: " << xor2Simple.size()
            << " cycles  (" << fixed << setprecision(4) << t_xor2 << " s)\n";
    }

    if (XOR3_ACTIVE) {
        auto t3_start = chrono::steady_clock::now();

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

        auto t3_end = chrono::steady_clock::now();
        double t_xor3 = chrono::duration_cast<chrono::nanoseconds>(t3_end - t3_start).count() / 1e9;
        cout << "  XOR-3 done: " << xor3Simple.size()
            << " cycles  (" << fixed << setprecision(4) << t_xor3 << " s)\n";
    }

    if (XOR4_ACTIVE) {
        auto t4_start = chrono::steady_clock::now();

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
                    xor4Simple.push_back(x4); // NEW: keep for XOR-5 layer
                }
            }
        }

        auto t4_end = chrono::steady_clock::now();
        double t_xor4 = chrono::duration_cast<chrono::nanoseconds>(t4_end - t4_start).count() / 1e9;
        cout << "  XOR-4 done: " << xor4Count
            << " cycles  (" << fixed << setprecision(4) << t_xor4 << " s)\n";
    }

    if (XOR5_ACTIVE) {
        auto t5_start = chrono::steady_clock::now();

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

        auto t5_end = chrono::steady_clock::now();
        double t_xor5 = chrono::duration_cast<chrono::nanoseconds>(t5_end - t5_start).count() / 1e9;
        cout << "  XOR-5 done: " << xor5Count
            << " cycles  (" << fixed << setprecision(4) << t_xor5 << " s)\n";
    }

    return vector<bitset<BSZ>>(seen.begin(), seen.end());
}

// =================== NEW HELPERS ===================

vector<std::vector<int>> buildSkeletonAdjacency(
    const std::vector<std::vector<int>>& skeleton, int n)
{
    vector<std::vector<int>> adj(n);
    for (const auto& e : skeleton) {
        int u = e[0];
        int v = e[1];
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    return adj;
}

BlockData buildBlocksFromSkeleton(
    const std::vector<std::vector<int>>& skeleton,
    const std::vector<int>& component,
    int n,
    std::ofstream& log)
{
    using namespace std;
    using namespace std::chrono;

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

vector<std::bitset<BSZ>> generateCyclesForSingleBlock( const vector<std::pair<int, int>>& blockEdges, const vector<int>& blockVerts,  int blockCompId, int maxEdgesPerCycle, Instance& inst)
{
    using namespace std;

    vector<std::bitset<BSZ>> result;

    int E = static_cast<int>(blockEdges.size());
    cout << "BlockEdges.size()=" << blockEdges.size() << " BSZ=" << BSZ << std::endl;
    if (E > BSZ) {
        std::cerr << "Warning ...\n" << std::flush;
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
    std::vector<int> U_local(E), V_local(E);
    for (int ei = 0; ei < E; ei++) {
        U_local[ei] = g2l[blockEdges[ei].first];
        V_local[ei] = g2l[blockEdges[ei].second];
    }

    // incidence list: vertex -> list of edge indices
    std::vector<std::vector<int>> incEdges(Vb);
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
    std::unordered_set<std::bitset<BSZ>, BitsetHash, BitsetEq> blockSet;
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

    std::vector<std::bitset<BSZ>> xorCycles;
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

    std::cout << "Block (comp " << blockCompId
        << ") has " << fundBits.size() << " fundamental cycles, "
        << result.size()
        << " simple cycles (fund+XOR2+XOR3+XOR4).\n";

    return result;
}

vector<std::vector<std::bitset<BSZ>>> generateCyclesForAllBlocks(  const BlockData& blocks,  int maxEdgesPerCycle, ofstream& log, Instance& inst)
{
    using namespace std;
    using namespace std::chrono;

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
    std::vector<int> blocksPerComponent(totalComponents, 0);
    for (int cid : blocks.blockCompId)
        blocksPerComponent[cid]++;

    for (size_t bi = 0; bi < totalBlocks; ++bi)
    {
        const auto& edges = blocks.blockEdgesAll[bi];
        const auto& verts = blocks.blockVertsAll[bi];
        int compId = blocks.blockCompId[bi];

        int blocksInThisComponent = blocksPerComponent[compId];

        std::cout << "\n===== Processing Block "
            << bi << " / " << (totalBlocks - 1)
            << "  (Component " << compId
            << " / " << (totalComponents - 1)
            << ", Blocks in this component: " << blocksInThisComponent
            << ") =====\n";

        std::cout << "Vertices: " << verts.size()
            << " | Edges: " << edges.size() << "\n";

        auto t0 = std::chrono::steady_clock::now();

        auto cycles = generateCyclesForSingleBlock(
            edges,
            verts,
            compId,
            maxEdgesPerCycle,
            inst
        );

        auto t1 = std::chrono::steady_clock::now();
        double t_sec =
            std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() / 1e9;

        std::cout << "Block cycle generation time: "
            << std::fixed << std::setprecision(4)
            << t_sec << " sec\n"
            << "Cycles found: " << cycles.size() << "\n";

        allBlockCyclesBits[bi] = std::move(cycles);
    }

    return allBlockCyclesBits;
}

RCCyclesResult buildGlobalRCCycles( const BlockData& blocks,  const vector<std::vector<std::bitset<BSZ>>>& allBlockCyclesBits, ofstream& log)
{
    using namespace std;

    RCCyclesResult res;
    //res.cycles_vertices.reserve(500000); //old estimate

    std::unordered_set<std::vector<int>, VecHash, VecEq> globalSeen;
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
            std::vector<int> verts;
            //verts.reserve(16); //old estimate
            int cntEdges = (int)bs.count(); //Each bitset cycle has cntEdges = bs.count() edges, and we push two endpoints per edge then dedupe, so the raw push count is 2*cntEdges.
            verts.reserve(2 * cntEdges);

            for (int ei = 0; ei < Eb; ++ei) {
                if (!bs.test(ei)) continue;
                verts.push_back(edges[ei].first);
                verts.push_back(edges[ei].second);
            }
            if (verts.empty()) continue;

            std::sort(verts.begin(), verts.end());
            verts.erase(std::unique(verts.begin(), verts.end()), verts.end());

            if (globalSeen.insert(verts).second) {
                res.cycles_vertices.push_back(std::move(verts));
            }
        }
    }

    return res;
}

// =================== MAIN ENTRY: Find_all_RC_CYCLES ===================

vector<std::vector<int>> Find_all_RC_CYCLES(Instance& inst)
{
    using namespace std;
    using namespace std::chrono;

    auto startALL = steady_clock::now();

    string savePath = "time/RC0_";
    string Path2 = savePath + inst.name + "_time.txt";
    ofstream log(Path2, ios::out | ios::app);

    int n = inst.nbCity;

    // Step 0: LP & skeleton
    auto t0 = steady_clock::now();
    auto result = SolveLPRelaxation(inst);
    vector<vector<int>> skeleton = result.first;
    vector<vector<int>> rem_RC_ordered = result.second; // unused here
    auto t1 = steady_clock::now();
    cout << "LP relaxation took " << duration_cast<nanoseconds>(t1 - t0).count() / 1e9 << " sec\n";
    inst.tLPrelaxation = duration_cast<nanoseconds>(t1 - t0).count() / 1e9;

    // Step 1: Connected components via spanning forest
    vector<int> component;
    vector<pair<int, int>> forest = buildSpanningForest(skeleton, n, component);

    // optional: per-component vertex lists
    //unordered_map<int, vector<int>> verticesByComponent;
    //for (int v = 0; v < (int)component.size(); ++v)
    //    verticesByComponent[component[v]].push_back(v);

    // Step 2 4: BCC -> blocks
    auto t2 = steady_clock::now();
    BlockData blocks = buildBlocksFromSkeleton(
        skeleton,
        component,
        n,
        log
    );
    auto t3 = steady_clock::now();
    cout << "Building Blocks From Skeleton took "
        << duration_cast<nanoseconds>(t3 - t2).count() / 1e9 << " sec\n";
    inst.tBuildingBlocksFromSkeleton = duration_cast<nanoseconds>(t3 - t2).count() / 1e9;

    // Step 5: per-block cycles
    auto t4 = steady_clock::now();
    int maxEdgesPerCycle = inst.nbCity / 2;
    auto allBlockCyclesBits = generateCyclesForAllBlocks(
        blocks,
        maxEdgesPerCycle,
        log,
        inst
    );
    auto t5 = steady_clock::now();
    cout << "Per-block fundamental + XOR generation took " << duration_cast<nanoseconds>(t5 - t4).count() / 1e9 << " sec\n";
    inst.tCycleXOR = duration_cast<nanoseconds>(t5 - t4).count() / 1e9;

    inst.nbRC_cycles_NOT_RED = allBlockCyclesBits.size();

    // Step 6: global RC cycles as vertex sets
    auto t6 = steady_clock::now();
    RCCyclesResult rcResult = buildGlobalRCCycles(
        blocks,
        allBlockCyclesBits,
        log
    );
    auto t7 = steady_clock::now();
    cout << "Build global RC took " << std::chrono::duration_cast<std::chrono::nanoseconds>(t7 - t6).count() / 1e9 << " sec\n";
    inst.tGlobalRCMerge = duration_cast<std::chrono::nanoseconds>(t7 - t6).count() / 1e9;

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
    auto endALL = steady_clock::now();
    cout << "All operations took " << duration_cast<nanoseconds>(endALL - startALL).count() / 1e9 << " sec\n";
    inst.tAllRC = duration_cast<nanoseconds>(endALL - startALL).count() / 1e9;

    //log.close();

    return rcResult.cycles_vertices;
}

// =================== COMPARE PRECUTS VS RC CYCLES ===================

void CompareCstar(std::vector<std::vector<int>>& RC_cycles_vertices, Instance& inst, const string& type)
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

        RC_hash.insert(std::move(key));
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

        if (RC_hash.count(key)) precut_match++;
    }

    // ---- Report ----
    double percentage = (precut_total == 0) ? 0.0: (100.0 * precut_match / precut_total);

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

std::vector<std::pair<int, int>> buildMST(const std::vector<std::vector<double>>& dist) {
    int n = (int)dist.size();
    vector<double> minEdge(n, 1e18);
    vector<int>    parent(n, -1);
    vector<bool>   used(n, false);

    minEdge[0] = 0.0;
    vector<std::pair<int, int>> mstEdges;

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

int getLCA(int a, int b,const vector<std::vector<int>>& up, const vector<int>& depth, int LOG)
{
    if (depth[a] < depth[b]) std::swap(a, b);
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

static double cliqueAverageDistance(const vector<int>& verts, const vector<vector<double>>& dist)
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

set<std::set<int>> Check_K_MSTCycles_Random(const int& needed, Instance& inst, Solution& sol)
{
    // Build "seen" from precuts once
    set<std::set<int>> seen;
    for (const auto& vec : inst.precuts) {
        seen.emplace(vec.begin(), vec.end());
    }
    seen.insert(sol.cuts.begin(), sol.cuts.end()); // Now RC cycles are also there.

    int n = inst.nbCity;
    using namespace std::chrono;
    auto startALL = steady_clock::now();

    // 1) Generate ALL MST cycles
    std::vector<std::vector<int>> allCycles = CheckAllMSTCycles(inst);

    // 2) Canonicalize + dedupe only
    std::unordered_set<std::vector<int>, VecHash, VecEq> uniq;
    uniq.reserve(allCycles.size() * 2);

    std::vector<std::vector<int>> uniqueCycles;
    uniqueCycles.reserve(allCycles.size());

    for (auto& cyc : allCycles)
    {
        // canonical form: sorted unique vertices
        sort(cyc.begin(), cyc.end());
        cyc.erase(std::unique(cyc.begin(), cyc.end()), cyc.end());

        // optional: keep length constraints here too
        int m = (int)cyc.size();
        if (m < 3 || m > n / 2) continue;

        if (!uniq.insert(cyc).second) continue;

        uniqueCycles.push_back(std::move(cyc));
    }

    // 3) Shuffle randomly instead of ordering
    std::random_device rd;
    std::mt19937 rng(rd());
    std::shuffle(uniqueCycles.begin(), uniqueCycles.end(), rng);

    // 4) Take first "needed" cycles not already in precuts / sol.cuts
    std::set<std::set<int>> result;
    int countAdded = 0;

    for (const auto& cyc : uniqueCycles)
    {
        if (countAdded >= needed) break;

        std::set<int> cycleSet(cyc.begin(), cyc.end());

        if (seen.insert(cycleSet).second) {
            result.insert(std::move(cycleSet));
            ++countAdded;
        }
    }

    auto endALL = steady_clock::now();
    inst.tMST_cycle_Gen = duration_cast<nanoseconds>(endALL - startALL).count() / 1e9;
    inst.nbMST_cycles = (int)result.size();

    return result;
}


// =================== MST CYCLES (CheckAllMSTCycles) ===================

std::vector<std::vector<int>> CheckAllMSTCycles(Instance& inst) {
    int n = inst.nbCity;
    using namespace std::chrono;

    auto startALL = steady_clock::now();

    std::vector<std::vector<int>> all_MST_cycles;
    all_MST_cycles.reserve(n * (n - 1) / 2 - (n - 1));

    // Step 1: Build MST
    std::vector<std::pair<int, int>> mstEdges = buildMST(inst.dist);

    // Build adjacency list for MST
    std::vector<std::vector<int>> adj(n);
    std::set<std::pair<int, int>> mstSet;
    for (auto e : mstEdges) {
        adj[e.first].push_back(e.second);
        adj[e.second].push_back(e.first);
        mstSet.insert({ std::min(e.first,e.second), std::max(e.first,e.second) });
    }

    // Collect non-MST edges
    std::vector<std::tuple<double, int, int>> nonMSTEdges;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (inst.dist[i][j] < 1e9) {
                std::pair<int, int> e = { i,j };
                if (mstSet.find(e) == mstSet.end()) {
                    nonMSTEdges.push_back(std::make_tuple(inst.dist[i][j], i, j));
                }
            }
        }
    }

    // LCA prep
    int LOG = (int)std::ceil(std::log2((double)n)) + 1;
    vector<std::vector<int>> up(n, std::vector<int>(LOG, -1));
    vector<int> depth(n, 0);

    std::function<void(int, int)> dfs = [&](int u, int p) {
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
    //std::sort(nonMSTEdges.begin(), nonMSTEdges.end());

    double threshold = (inst.lp_obj / inst.nbCity) * MULTIPLIER_MST_EDGE_INCLUSION;

    for (int k = 0; k < (int)nonMSTEdges.size(); k++) {
        double w;
        int u, v;
        std::tie(w, u, v) = nonMSTEdges[k];
	
	// threshold filter
        if (w > threshold) continue; // if it is too large, skip

        int L = getLCA(u, v, up, depth, LOG);

        int len = depth[u] + depth[v] - 2 * depth[L] + 1;
        if (len < 3 || len > n / 2) continue;

        std::vector<int> cycle;
        int x = u;
        while (x != L) {
            cycle.push_back(x);
            x = up[x][0];
        }
        cycle.push_back(L);
        std::vector<int> tail;
        x = v;
        while (x != L) {
            tail.push_back(x);
            x = up[x][0];
        }
        std::reverse(tail.begin(), tail.end());
        cycle.reserve(cycle.size() + tail.size());
        cycle.insert(cycle.end(), tail.begin(), tail.end());

        all_MST_cycles.push_back(std::move(cycle));
    }

    auto endALL = steady_clock::now();
    //cout << "Whole process took = " << duration / 1e9 << " sec\n";

    inst.tMST_cycle_Gen = duration_cast<nanoseconds>(endALL - startALL).count() / 1e9;
    inst.nbMST_cycles = all_MST_cycles.size();

    return all_MST_cycles;
}

// =================== MERGE & DELETE DUPLICATE CYCLES ===================

vector<std::vector<int>> mergeAndDedupe(const vector<std::vector<int>>& A, const vector<std::vector<int>>& B)
{
    using namespace std;

    unordered_set<vector<int>, VecHash, VecEq> seen;
    seen.reserve(A.size() + B.size());

    vector<vector<int>> result;
    result.reserve(A.size() + B.size());

    auto insertCycle = [&](const vector<int>& cyc)
    {
        vector<int> v = cyc;
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());

        auto [it, inserted] = seen.insert(v);
        if (inserted)
            result.push_back(std::move(v));
    };

    for (const auto& cyc : A) insertCycle(cyc);
    for (const auto& cyc : B) insertCycle(cyc);

    return result;
}

set<set<int>> Check_K_RCCycles_Random(const int& needed, Instance& inst)
{
    // Build "seen" from precuts once
    std::set<std::set<int>> seen;
    for (const auto& vec : inst.precuts) {
        seen.emplace(vec.begin(), vec.end());
    }

    using namespace std::chrono;
    auto startALL = steady_clock::now();

    // 1) Generate ALL RC cycles
    std::vector<std::vector<int>> allRCCycles = Find_all_RC_CYCLES(inst);

    // 2) sort aand dedupe
    std::unordered_set<std::vector<int>, VecHash, VecEq> uniq;
    uniq.reserve(allRCCycles.size() * 2);

    vector<vector<int>> uniqueCycles;
    uniqueCycles.reserve(allRCCycles.size());

    const int n = inst.nbCity;

    for (auto& cyc : allRCCycles){
        std::sort(cyc.begin(), cyc.end());
        cyc.erase(std::unique(cyc.begin(), cyc.end()), cyc.end());

        const int m = (int)cyc.size();
        if (m < 3 || m > n / 2) continue;

        if (!uniq.insert(cyc).second) continue;

        uniqueCycles.push_back(std::move(cyc));
    }

    // 3) Shuffle randomly
    std::random_device rd;
    std::mt19937 rng(rd());
    std::shuffle(uniqueCycles.begin(), uniqueCycles.end(), rng);

    // 4) Take first "needed" cycles not already seen
    std::set<std::set<int>> result;
    int countAdded = 0;

    for (const auto& cyc : uniqueCycles)
    {
        if (countAdded >= needed) break;

        std::set<int> cycleSet(cyc.begin(), cyc.end());

        if (seen.insert(cycleSet).second) {
            result.insert(std::move(cycleSet));
            ++countAdded;
        }
    }

    auto endALL = steady_clock::now();

    // timing field
    // inst.tRC_cycle_Gen = duration_cast<nanoseconds>(endALL - startALL).count() / 1e9;

    return result;
}
