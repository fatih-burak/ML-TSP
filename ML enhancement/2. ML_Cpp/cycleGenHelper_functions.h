#ifndef CYCLEGENHELPER_FUNCTIONS_H
#define CYCLEGENHELPER_FUNCTIONS_H

#include "helper_functions.h"   // for Instance, GRB*, etc.

#include <vector>
#include <bitset>
#include <utility>
#include <fstream>
#include <numeric>   // iota
#include <random>
#include <queue>
#include <functional>
#include <chrono>
#include <cstdint>
#include <unordered_set>
#include <unordered_map>

#include <stack>
#include <cmath>
#include <algorithm>
#include <iostream>

#include <set>
#include <tuple>
#include <string>
#include <sstream>


// =================== CONSTANTS (shared) ===================
inline constexpr int    BSZ = 2000;      // max edges per block
// =================== DATA STRUCTS (shared) ===================

struct BlockData {
    std::vector<std::vector<std::pair<int, int>>> blockEdgesAll;
    std::vector<std::vector<int>>                 blockVertsAll;
    std::vector<int>                              blockCompId;   // component id per block
};

struct RCCyclesResult {
    std::vector<std::vector<int>> cycles_vertices;
};

struct BCCResult {
    std::vector<std::vector<std::pair<int, int>>> blocks;     // edges
    std::vector<std::vector<int>>                 blockVerts; // vertices
};

// =================== PUBLIC API ===================

// LP relaxation -> skeleton + remaining edges
std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>>
SolveLPRelaxation(Instance& inst);

// Spanning forest and components from skeleton
std::vector<std::pair<int, int>>
buildSpanningForest(const std::vector<std::vector<int>>& skeleton,
    int n,
    std::vector<int>& component);

// Main RC cycle generator
std::vector<std::vector<int>> Find_all_RC_CYCLES(Instance& inst);

// All MST-based cycles
std::vector<std::vector<int>> CheckAllMSTCycles(Instance& inst);

// Merge & dedupe two sets of cycles (vertex sets)
std::vector<std::vector<int>> mergeAndDedupe(
    const std::vector<std::vector<int>>& A,
    const std::vector<std::vector<int>>& B);

// =================== LOWER-LEVEL HELPERS (exposed if you want re-use) ===================

bool isSimpleCycle_fast(const std::bitset<BSZ>& b,
    const std::vector<int>& U,
    const std::vector<int>& V,
    const std::vector<std::vector<int>>& inc);

BCCResult tarjanBCC(const std::vector<std::vector<int>>& adj, int n);

std::vector<std::vector<int>> computeFundamentalCycles_Block(
    const std::vector<std::pair<int, int>>& blockEdges,
    const std::vector<int>& blockVerts,
    const int& maxEdgesPerCycle);

std::vector<std::vector<int>> computeFundamentalCycles_Block_MST(
    const std::vector<std::pair<int, int>>& blockEdges,     // GLOBAL vertex ids
    const std::vector<int>& blockVerts,     // GLOBAL vertex ids
    const std::vector<std::vector<double>>& dist,           // inst.dist
    const int& maxEdgesPerCycle);

std::vector<std::bitset<BSZ>> encodeCycles_Block(
    const std::vector<std::vector<int>>& cycles,
    const std::vector<std::pair<int, int>>& blockEdges);

std::vector<std::bitset<BSZ>> generateXORCycles_Block_old(
    const std::vector<std::bitset<BSZ>>& fundBits,
    const std::vector<int>& U_local,
    const std::vector<int>& V_local,
    const std::vector<std::vector<int>>& incEdges,
    int                                  maxEdgesPerCycle);

std::vector<std::bitset<BSZ>> generateXORCycles_Block(
    const std::vector<std::bitset<BSZ>>& fundBits,
    const std::vector<int>& U_local,
    const std::vector<int>& V_local,
    const std::vector<std::vector<int>>& incEdges,
    int                                  maxEdgesPerCycle);

// Skeleton -> adjacency
std::vector<std::vector<int>> buildSkeletonAdjacency(
    const std::vector<std::vector<int>>& skeleton,
    int n);

// Build blocks from skeleton + component info
BlockData buildBlocksFromSkeleton(
    const std::vector<std::vector<int>>& skeleton,
    const std::vector<int>& component,
    int                                  n,
    std::ofstream& log);

// Per-block cycles
std::vector<std::bitset<BSZ>> generateCyclesForSingleBlock(
    const std::vector<std::pair<int, int>>& blockEdges,
    const std::vector<int>& blockVerts,
    int                                     blockCompId,
    int                                     maxEdgesPerCycle,
    Instance& inst);

std::vector<std::vector<std::bitset<BSZ>>> generateCyclesForAllBlocks(
    const BlockData& blocks,
    int              maxEdgesPerCycle,
    std::ofstream& log,
    Instance& inst);

// Global RC (vertex sets)
RCCyclesResult buildGlobalRCCycles(
    const BlockData& blocks,
    const std::vector<std::vector<std::bitset<BSZ>>>& allBlockCyclesBits,
    std::ofstream& log);

// MST helpers
std::vector<std::pair<int, int>>
buildMST(const std::vector<std::vector<double>>& dist);

std::set<set<int>>
Check_K_MSTCycles(const int& needed, Instance& inst);


std::set<std::set<int>>
Check_K_MSTCycles_Ordered_AvgCliqueLen(const int& needed, Instance& inst, Solution& sol);

std::set<std::set<int>>
Check_K_RCCycles_Ordered_AvgCliqueLen(const int& needed, Instance& inst);

int getLCA(int a, int b,
    const std::vector<std::vector<int>>& up,
    const std::vector<int>& depth,
    int                                  LOG);

std::vector<std::vector<int>> SelectPrecutsWithNN(
    const std::vector<std::vector<int>>& sets,
    const std::vector<std::vector<int>>& mstCycles,
    const std::vector<std::vector<int>>& rcCycles,
    const std::vector<std::vector<double>>& dist,
    const std::vector<std::vector<double>>& lp_x,
    const std::vector<std::vector<double>>& lp_rc,
    const Instance& inst,
    const std::string& modelPath,
    double threshold);


#endif // CYCLEGENHELPER_FUNCTIONS_H
