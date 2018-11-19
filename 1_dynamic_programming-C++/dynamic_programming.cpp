#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <limits>
#include <math.h>
#include <iomanip>

#include "DirectedGraph.cpp"

template<class T>
class DynamicProgramming
{
    public:
    DirectedGraph<T> inputGraph;
    std::map<T, float> valueFunc;
    std::map<T, T> backPointer;

    DynamicProgramming( DirectedGraph<T> inputGraph )
    {
        this->inputGraph = inputGraph;
        for( auto node : inputGraph.getNodes() )
        {
            valueFunc[node] = std::numeric_limits<float>::infinity();
        }
        valueFunc[inputGraph.goalNode] = 0;

        bool updatePerformed = false;
        size_t iterCount = 0;
        float newValue;
        while( iterCount++ < 100 )
        {
            updatePerformed = false;
            for( auto nodeValuePair : valueFunc )
            {
                if( nodeValuePair.second != std::numeric_limits<float>::infinity() )
                {
                    for( auto nborNode : inputGraph.getInNeighbors(nodeValuePair.first) )
                    {
                        newValue = nodeValuePair.second + inputGraph.getEdgeCost( std::make_pair( nodeValuePair.first, nborNode) );
                        
                        if( newValue < valueFunc[nborNode] )
                        {
                            valueFunc[nborNode] = newValue;
                            backPointer[nborNode] = nodeValuePair.first;
                            // std::cout<<"bN: "<<nborNode<<"->"<<nodeValuePair.first<<std::endl;
                            updatePerformed = true;
                        }
                    }
                }
            }
            // std::cout<<iterCount<<" : "<<valueFunc[inputGraph.startNode]<<" --> "<<oldValueStartNode<<std::endl;

            if(!updatePerformed)
            {
                printPath();
                printValueFuncTable();
                break;
            }
        }
    }

    void printPath()
    {
        auto tnode = inputGraph.startNode;
        std::cout<<"startPath: from "<<inputGraph.startNode<<" to "<<inputGraph.goalNode<<std::endl;
        while( tnode != inputGraph.goalNode )
        {
            std::cout<< tnode << " ";
            tnode = backPointer[tnode];
        }
        std::cout<<tnode<<std::endl;
    }

    void printValueFuncTable()
    {
        for( auto nodeValuePair : valueFunc )
        {
            // std::cout<<nodeValuePair.first<<" --> "<<nodeValuePair.second<<std::endl;
            std::cout<< std::fixed << std::setw( 11 ) << std::setprecision( 6 ) << nodeValuePair.second<< " ";
        }
        std::cout<<std::endl;
    }
};

int main()
{
    const std::string input_file = "input1.txt";
    const std::string output_file = "my_output1.txt";


    DirectedGraph<uint32_t> myGraph(input_file);
    // myGraph.printGraphAdjacency();

    DynamicProgramming<uint32_t> myAlg(myGraph);

    // myAlg.printValueFuncTable();

    return 0;
}
