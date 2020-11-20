#include "CliqueSolver.h"

#include <thread>
#include <unordered_set>


using namespace std;


namespace szx {
namespace tsm {

int loadInput(NextEdge nextEdge, const Arr<Weight> &nodeWeights,
    InitOrdering initOrdering = InitOrdering::DegreeDegeneracy, ID lowerBound = 0);
bool tsmMain(Clique &sln, int msTimeout);


bool solveWeightedMaxClique(Clique &sln, NextEdge nextEdge, const Arr<Weight> &nodeWeights, Millisecond timeout, ID lowerBound) {
    bool r = false;
    thread t([&]() {
        if (loadInput(nextEdge, nodeWeights, InitOrdering::DegreeDegeneracy, lowerBound)) {
            r = tsmMain(sln, timeout);
        }
    });
    t.join();
    return r;
}

bool solveWeightedMaxClique(Clique &sln, NextEdge nextEdge, ID nodeNum, Millisecond timeout, ID lowerBound) {
    return solveWeightedMaxClique(sln, nextEdge, Arr<Weight>(nodeNum, 1), timeout, lowerBound);
}

bool solveWeightedMaxClique(Clique &sln, const AdjList &adjList, const Arr<Weight> &nodeWeights, Millisecond timeout, ID lowerBound) {
    ID i = 0;
    ID j = 0;
    auto nextEdge = [&](ID &src, ID &dst) {
        for (;;) {
            if (j >= adjList[i].size()) {
                j = 0;
                ++i;
            }
            if (i >= adjList.size()) {
                i = 0;
                return false;
            }
            if (i >= adjList[i][j].dst) { ++j; continue; }
            src = i;
            dst = adjList[i][j].dst;
            ++j;
            return true;
        }
    };
    return solveWeightedMaxClique(sln, nextEdge, nodeWeights, timeout, lowerBound);
}

bool solveWeightedMaxClique(Clique &sln, const AdjList &adjList, Millisecond timeout, ID lowerBound) {
    return solveWeightedMaxClique(sln, adjList, Arr<Weight>(adjList.size(), 1), timeout, lowerBound);
}

bool solveWeightedMaxClique(Clique &sln, const EdgeList &edgeList, const Arr<Weight> &nodeWeights, Millisecond timeout, ID lowerBound) {
    auto j = edgeList.begin();
    auto nextEdge = [&](ID &src, ID &dst) {
        for (unordered_set<unsigned long long> edgeSet;;) {
            if (j == edgeList.end()) {
                j = edgeList.begin();
                return false;
            }
            
            unsigned long long edge = min(j->src, j->dst);
            edge <<= (sizeof(ID) * 8);
            edge += max(j->src, j->dst);
            if (edgeSet.find(edge) != edgeSet.end()) { continue; }
            edgeSet.insert(edge);

            src = j->src;
            dst = j->dst;
            ++j;
            return true;
        }
    };
    return solveWeightedMaxClique(sln, nextEdge, nodeWeights, timeout, lowerBound);
}

bool solveWeightedMaxClique(Clique &sln, const EdgeList &edgeList, ID nodeNum, Millisecond timeout, ID lowerBound) {
    return solveWeightedMaxClique(sln, edgeList, Arr<Weight>(nodeNum, 1), timeout, lowerBound);
}

bool solveWeightedMaxClique(Clique &sln, const AdjMat &adjMat, const Arr<Weight> &nodeWeights, Millisecond timeout, ID lowerBound) {
    ID i = 0;
    ID j = 0;
    auto nextEdge = [&](ID &src, ID &dst) {
        for (;;) {
            if (j >= i) {
                j = 0;
                ++i;
            }
            if (i >= adjMat.size1()) {
                i = 0;
                return false;
            }
            if (!adjMat.at(i, j)) {
                ++j;
                continue;
            }
            src = i;
            dst = j;
            ++j;
            return true;
        }
    };
    return solveWeightedMaxClique(sln, nextEdge, nodeWeights, timeout, lowerBound);
}

bool solveWeightedMaxClique(Clique &sln, const AdjMat &adjMat, Millisecond timeout, ID lowerBound) {
    return solveWeightedMaxClique(sln, adjMat, Arr<Weight>(adjMat.size1(), 1), timeout, lowerBound);
}

bool solveWeightedIndependentSet(Clique &sln, const AdjMat &adjMat, const Arr<Weight> &nodeWeights, Millisecond timeout, ID lowerBound) {
    ID i = 0;
    ID j = 0;
    auto nextEdge = [&](ID &src, ID &dst) {
        for (;;) {
            if (j >= i) {
                j = 0;
                ++i;
            }
            if (i >= adjMat.size1()) {
                i = 0;
                return false;
            }
            if (adjMat.at(i, j)) {
                ++j;
                continue;
            }
            src = i;
            dst = j;
            ++j;
            return true;
        }
    };
    return solveWeightedMaxClique(sln, nextEdge, nodeWeights, timeout, lowerBound);
}

bool solveWeightedIndependentSet(Clique &sln, const AdjMat &adjMat, Millisecond timeout, ID lowerBound) {
    return solveWeightedIndependentSet(sln, adjMat, Arr<Weight>(adjMat.size1(), 1), timeout, lowerBound);
}

}
}
