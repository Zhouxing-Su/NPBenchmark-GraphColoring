////////////////////////////////
/// usage : 1.	the interface for weighted max clique solvers.
/// 
/// note  : 1.	
////////////////////////////////

#ifndef SMART_SZX_GOAL_MAX_CLIQUE_LIB_CLIQUE_SOLVER_H
#define SMART_SZX_GOAL_MAX_CLIQUE_LIB_CLIQUE_SOLVER_H


#include <functional>

#include "Graph.h"


namespace szx {
namespace tsm {

using ID = Graph::ID; // node ID counts from 0 in the interface, but counts from 1 inside the LKH.
using Weight = int;

using AdjNode = Graph::AdjNode;
using Edge = Graph::Edge;

using AdjMat = Graph::AdjMat<bool>;
using AdjList = Graph::AdjList<AdjNode>;
using EdgeList = Graph::EdgeList<Edge>;

using NodeList = Graph::List<ID>; // `nodeList` is a list of node IDs in increasing order.
using NodeSet = Graph::List<bool>; // `nodeSet[n]` is true if node `n` is included in the set.

using Clique = Graph::Clique<Weight>;

using NextEdge = std::function<bool(ID &src, ID &dst)>; // auto reset after the end of edges.


static constexpr ID CliqueIdBase = 1;

enum InitOrdering {
    DegreeDegeneracy = 1,
    WeightDegeneracy = 3
};


bool solveWeightedMaxClique(Clique &sln, NextEdge nextEdge, const Arr<Weight> &nodeWeights);
bool solveWeightedMaxClique(Clique &sln, NextEdge nextEdge, ID nodeNum);

bool solveWeightedMaxClique(Clique &sln, const AdjList &adjList, const Arr<Weight> &nodeWeights);
bool solveWeightedMaxClique(Clique &sln, const AdjList &adjList);

bool solveWeightedMaxClique(Clique &sln, const EdgeList &edgeList, const Arr<Weight> &nodeWeights);
bool solveWeightedMaxClique(Clique &sln, const EdgeList &edgeList, ID nodeNum);

bool solveWeightedMaxClique(Clique &sln, const AdjMat &adjMat, const Arr<Weight> &nodeWeights);
bool solveWeightedMaxClique(Clique &sln, const AdjMat &adjMat);

}
}


#endif // SMART_SZX_GOAL_MAX_CLIQUE_LIB_CLIQUE_SOLVER_H
