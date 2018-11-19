#include "DirectedGraph.h"

template<class T>
DirectedGraph<T>::DirectedGraph(size_t _numNodes, T _startNode, T _goalNode)
    : numNodes(_numNodes)
    , startNode(_startNode)
    , goalNode(_goalNode)
{
    this->addNode(startNode);
    this->addNode(goalNode);
}

template<class T>
DirectedGraph<T>::DirectedGraph(const std::string input_file)
{
    std::ifstream input_data;
    input_data.open( input_file.c_str() );

    uint32_t _startNode, _goalNode;
    if( input_data.is_open() )
    {
        input_data >> numNodes;
        input_data >> _startNode;
        input_data >> _goalNode;

        startNode = T(_startNode);
        goalNode = T(_goalNode);

        int32_t edgeStartNode, edgeEndNode;
        float edgeCost;

        this->addNode( T(startNode) );
        this->addNode( T(goalNode) );

        // std::cout<<"Graph:: start:"<<startNode<<" goal: "<<goalNode<<std::endl;
        while(input_data)
        {
            input_data >> edgeStartNode;
            input_data >> edgeEndNode;
            input_data >> edgeCost;

            this->addEdge( std::make_pair( T(edgeStartNode), T(edgeEndNode) ), edgeCost );
            // std::cout<<"Adding edge: "<<edgeStartNode<<"->"<<edgeEndNode<<" : "<<edgeCost<<std::endl;
        }
        input_data.close();
    }
}

template<class T>
DirectedGraph<T>::DirectedGraph(){}

template<class T>
void DirectedGraph<T>::addEdge( uPair<T> edge,  float edgeCost )
{
    //If the start node is not in the vertices list, add it
    if( this->nodes.find(edge.first) == this->nodes.end() )
    {
        this->addNode( edge.first );
    }
    //If the end node is not in the vertices list, add it      
    if( this->nodes.find(edge.second) == this->nodes.end() )
    {
        this->addNode( edge.second );
    }
    this->edgeCosts[ edge ] = edgeCost;
    this->outNeighbors[ edge.first ].push_back(edge.second);
    this->inNeighbors[edge.second].push_back(edge.first);
}

template<class T>
void DirectedGraph<T>::addNode( T node )
{
    this->nodes.insert(node);
    this->numNodes++;
}

template<class T>
void DirectedGraph<T>::printGraphAdjacency()
{
    // std::cout<<"Nodes >> outNeighbors\n";
    for(auto n:this->nodes)
    {
        std::cout<<n<<" >> ";
        for(auto nbor : this->outNeighbors[n])
        {
            std::cout<<nbor<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

template<class T>
void DirectedGraph<T>::printEdgeCosts()
{
    std::cout<<std::endl;
    for(uint32_t node : this->nodes)
    {
        for(uint32_t nbor: this->outNeighbors[node])
        {
            std::cout<<node<<"->"<<nbor<<" : "<<this->edgeCosts[ std::make_pair(node, nbor) ]<<std::endl;
        }
    }
}

template<class T>
uint32_t DirectedGraph<T>::getNumNodes()
{
    return this->numNodes;
}

template<class T>
std::set<T> DirectedGraph<T>::getNodes()
{
    return this->nodes;
}

template<class T>
std::vector<T> DirectedGraph<T>::getInNeighbors( T node)
{
    auto nborsIterator = this->inNeighbors.find(node);
    if( nborsIterator != this->inNeighbors.end() )
    {
        return nborsIterator->second;
    }
    else
    {
        return std::vector<T>();
    }
}

template<class T>
std::vector<T> DirectedGraph<T>::getOutNeighbors( T node )
{
    auto nborsIterator = this->outNeighbors.find(node);
    if( nborsIterator != this->outNeighbors.end() )
    {
        return nborsIterator->second;
    }
    else
    {
        return std::vector<T>();
    }
}

template<class T>
float DirectedGraph<T>::getEdgeCost(uPair<T> edge)
{
    auto edgeCostIterator = this->edgeCosts.find(edge);
    if( edgeCostIterator != this->edgeCosts.end() )
    {
        return edgeCostIterator->second;
    }
    else
    {
        return std::numeric_limits<float>::infinity();
    }
}