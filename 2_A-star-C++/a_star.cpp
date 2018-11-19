#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <limits>
#include <math.h>
#include <iomanip>
#include <list>
#include <queue> 
#include <functional>

#include "DirectedGraph.cpp"

// template<class T>
// using qPair = std::pair<double, T>;

template<class T>
class qPair
{
    public:
    double valueFunc;
    T node;

    qPair(T _node, double _valueFunc) : valueFunc(_valueFunc), node(_node) {}
    qPair() : valueFunc(0), node(0) {}
    bool operator() (qPair nodeA, qPair nodeB)
    {
        bool ret_val = true;
        if( nodeA.valueFunc < nodeB.valueFunc )
        {
            ret_val = false;
        }
        return ret_val;
    }
};

template<class T>
class ShortestPath
{
    public:
    DirectedGraph<T> inputGraph;
    double (*heuristic)(T, T) = NULL;
    size_t num_iter = 0;
    bool compareNodes(std::pair<double, T> nodeA, std::pair<double, T> nodeB)
    {
        bool ret_val = true;
        if( nodeA.second < nodeB.second )
        {
            ret_val = false;
        }
        return ret_val;
    }

    std::priority_queue< qPair<T>, std::vector<qPair<T>>, qPair<T> > open;
    std::unordered_map<T, double> costToCome;
    std::map<T, T> backPointer;

    ShortestPath( DirectedGraph<T> &_inputGraph, double (*_heuristic)(T, T) )
    : inputGraph(_inputGraph), heuristic(_heuristic)
    {

        open.push( qPair<T>(inputGraph.startNode, 0) );

        for( auto node : inputGraph.getNodes() )
        {
            costToCome[node] = std::numeric_limits<double>::infinity();
        }
        costToCome[inputGraph.startNode] = 0;

        qPair<T> currQPair;
        double newCostToCome;
        bool foundPath = false;
        while( !open.empty() )
        {
            num_iter++;
            currQPair = open.top();
            open.pop();

            // std::cout<<"Current node: "<<currQPair.node<<" "<<currQPair.valueFunc<<std::endl;

            for( auto nborNode : inputGraph.getOutNeighbors(currQPair.node) )
            {
                newCostToCome = currQPair.valueFunc + inputGraph.getEdgeCost( std::make_pair( currQPair.node, nborNode ) ) + heuristic(nborNode, inputGraph.goalNode);
                // std::cout<<"\tNbor node: "<< nborNode << " " << costToCome[nborNode] << " " << newCostToCome <<std::endl;
                if( newCostToCome < costToCome[nborNode] )
                {
                    costToCome[nborNode] = newCostToCome;
                    backPointer[nborNode] = currQPair.node;

                    open.push( qPair<T>(nborNode, newCostToCome) );
                }
            }
            if(currQPair.node == inputGraph.goalNode)
            {
                foundPath = true;
                break;
            }
        }

        if(foundPath)
        {
            printPath();
        }
        else
        {
            std::cout<<"No feasible path from start node to goal node!"<<std::endl;
        }
    }

    void printPath()
    {
        auto tnode = inputGraph.goalNode;
        std::list<T> outputPath;
        // std::cout<<"startPath: from "<<inputGraph.goalNode<<" to "<<inputGraph.startNode<<std::endl;
        while( tnode != inputGraph.startNode )
        {
            // std::cout<< tnode << " ";
            outputPath.push_front(tnode);
            tnode = backPointer[tnode];
        }
        outputPath.push_front(tnode);
        // std::cout<<tnode<<std::endl;

        for(auto node : outputPath )
        {
            std::cout<< node <<" ";
        }
        std::cout<<std::endl;
        std::cout<< costToCome[inputGraph.goalNode] <<'\t'<< num_iter <<std::endl;
    }

};

class eucNode
{
    public:
    size_t index;
    double x;
    double y;
    
    eucNode(size_t _index, double _x, double _y) : index(_index), x(_x), y(_y) {}
    eucNode(size_t _index) : index(_index), x(0.0), y(0.0) {}
    eucNode() : index(0), x(0.0), y(0.0) {}
    
    void setCoords(eucNode &node, double _x, double _y)
    {
        node.x = _x;
        node.y = _y;
    }

    bool operator< (const eucNode node)
    {
        return (this->index < node.index);
    }

    bool operator> (const eucNode node)
    {
        return (this->index > node.index);
    }

    bool operator() (const eucNode node)
    {
        return (node.index);
    }
};

bool operator< (const eucNode nodeA, const eucNode nodeB)
{
    return (nodeA.index < nodeB.index);
}

bool operator> (const eucNode nodeA, const eucNode nodeB)
{
    return (nodeA.index > nodeB.index);
}

bool operator== (eucNode nodeA, eucNode nodeB)
{
    return (nodeA.index == nodeB.index);
}

bool operator!= (eucNode nodeA, eucNode nodeB)
{
    return (nodeA.index != nodeB.index);
}

std::ostream& operator<< (std::ostream& os, const eucNode& node) {
    os << node.index;
    return os;
}

namespace std
{
    template <> struct hash<eucNode>
    {
        size_t operator()(const eucNode &node) const
        {
            return node.index;
        }
    };
}

template<typename T>
double zero_dist(T nodeA, T nodeB)
{
    return 0.0;
}

double euclidean_dist(eucNode nodeA, eucNode nodeB)
{
    double ret_dist = sqrt( pow(nodeA.x - nodeB.x, 2) + pow(nodeA.y - nodeB.y, 2) );
    std::cout<<"Euclidean Distance: "<<ret_dist<<std::endl;
    return ret_dist;
}


// DirectedGraph<T>::DirectedGraph(const std::string input_file)
// {
//     std::ifstream input_data;
//     input_data.open( input_file.c_str() );

//     uint32_t _startNode, _goalNode;
//     if( input_data.is_open() )
//     {
//         input_data >> numNodes;
//         input_data >> _startNode;
//         input_data >> _goalNode;

//         startNode = T(_startNode);
//         goalNode = T(_goalNode);

//         int32_t edgeStartNode, edgeEndNode;
//         float edgeCost;

//         this->addNode( T(startNode) );
//         this->addNode( T(goalNode) );

//         // std::cout<<"Graph:: start:"<<startNode<<" goal: "<<goalNode<<std::endl;
//         while(input_data)
//         {
//             input_data >> edgeStartNode;
//             input_data >> edgeEndNode;
//             input_data >> edgeCost;

//             this->addEdge( std::make_pair( T(edgeStartNode), T(edgeEndNode) ), edgeCost );
//             // std::cout<<"Adding edge: "<<edgeStartNode<<"->"<<edgeEndNode<<" : "<<edgeCost<<std::endl;
//         }
//     }
// };


class EuclideanGraph : public DirectedGraph<eucNode>
{
    public:
    EuclideanGraph( const std::string input_file, const std::string coords_file )
    {
        std::ifstream input_data;
        std::ifstream coords_data;

        double _x, _y;
        size_t nodeCount = 0;
        eucNode* _node;

        coords_data.open( coords_file.c_str() );
        if( coords_data.is_open() )
        {
            while(coords_data)
            {
                coords_data >> _x >> _y;
                nodeCount++;
                _node = new eucNode(nodeCount, _x, _y);
                this->nodes.insert(*_node);
            }
        }

        uint32_t _edgeStartIdx, _edgeEndIdx, _edgeCost;
        uint32_t _startNodeIdx, _goalNodeIdx;

        input_data.open( input_file.c_str() );
        if(input_data.is_open())
        {
            input_data >> this->numNodes;
            input_data >> _startNodeIdx;
            input_data >> _goalNodeIdx;
            this->startNode = *std::next(this->nodes.begin(), _startNodeIdx);
            this->goalNode =  *std::next(this->nodes.begin(), _goalNodeIdx);

            while(input_data)
            {
                input_data >> _edgeStartIdx >> _edgeEndIdx >> _edgeCost;

                std::cout<<"Reading edge data: "<<_edgeStartIdx<<" "<<_edgeEndIdx<<"->"<<_edgeCost<<std::endl;
                auto _edgeStartNode = std::next(this->nodes.begin(), _edgeStartIdx);
                auto _edgeEndNode = std::next(this->nodes.begin(), _edgeEndIdx);
                this->edgeCosts[ std::make_pair( *_edgeStartNode, *_edgeEndNode) ] = _edgeCost;
                this->inNeighbors[*_edgeStartNode].push_back( *_edgeEndNode );
                this->outNeighbors[*_edgeEndNode].push_back( *_edgeStartNode );
            }
            input_data.close();
        }
    }

    // ~EuclideanGraph()
    // {
    //     for(auto node : this->nodes)
    //     {
    //         delete &node;
    //     }
    // }
};

class ShortestEucPath : public ShortestPath<eucNode>
{   
    public:
    ShortestEucPath( EuclideanGraph _inputGraph, double (*_heuristic)(eucNode, eucNode) )
    : ShortestPath(_inputGraph, _heuristic) {}
};

int main()
{
    const std::string input_file = "input_3.txt";
    const std::string coords_file = "coords_3.txt";

    // const std::string output_file = "my_output2.txt";

    DirectedGraph<uint32_t> myGraph(input_file);

    std::cout<<"Dijkstra:"<<std::endl;
    ShortestPath<uint32_t> dijkstra( myGraph, zero_dist );

    EuclideanGraph eucGraph(input_file, coords_file);

    // eucGraph.printGraphAdjacency();

    std::cout<<"A-star:"<<std::endl;
    ShortestEucPath astar( eucGraph, euclidean_dist );

    return 0;
}
