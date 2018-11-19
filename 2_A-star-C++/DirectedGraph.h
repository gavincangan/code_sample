#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <vector>
#include <unordered_map>
#include <set>
#include <unordered_set>

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        return h1 ^ h2;  
    }
};

template<class T>
using uPair = std::pair<T, T>;

template<class T>
using UnorderedMap_uPair_float = std::unordered_map<uPair<T>, float, pair_hash>;

template<class T>
class Graph
{

};

template<class T>
class DirectedGraph : Graph<T>
{
    public:
    T startNode, goalNode;
    size_t numNodes;

    std::set<T> nodes;
    std::unordered_map<T, std::vector<T> > inNeighbors;
    std::unordered_map<T, std::vector<T> > outNeighbors;
    UnorderedMap_uPair_float<T> edgeCosts;

    void addEdge( uPair<T> edge,  float edgeCost );
    void addNode( T node );
    void printGraphAdjacency();
    void printEdgeCosts();
    uint32_t getNumNodes();
    std::set<T> getNodes();
    std::vector<T> getInNeighbors( T node);
    std::vector<T> getOutNeighbors( T node);
    float getEdgeCost(uPair<T> edge);

    DirectedGraph(size_t _numNodes, T _startNode, T _goalNode);
    DirectedGraph(std::string input_file);
    DirectedGraph();
};