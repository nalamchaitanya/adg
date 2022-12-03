#include <iostream>
#include <math.h>
#include <algorithm>
#include <cassert>
#include <string>

using namespace std;
using std::min;
using std::max;

int* parseInput(string inputFile)
{
    // TODO aditya: malloc the graph here after finding the number of edges and length of array needed
    // Determine the number of colors D needed
}

int main()
{
    if(argc == 1)
    {
        cout << "No input" << endl;
    }
    // TODO by aditya
    int* graph = parseInput(argv[1]);
    
}

__device__ int getColor(int* graph, int* rho, int* C, int v, int D)
{
    bool B[D] = 0;
    // if(v==n) TODO make sure n+1 th entry should be the end index of the array to make sure this works.
    // This is very important as we do not want if statement here as this function gets used a lot of times
    for(int i = graph[v];i<graph[v+1];i++)
    {
        if(rho[graph[i]] > rho[v])
        {
            if(C[graph[i]] == -1)
            {
                return -1;
            }
            else
            {
                B[C[graph[i]]] = true;
            }
        }
    }
    for(int i =0;i<D;i++)
    {
        if(!B[i])
        {
            return i;
        }
    }
    // Should not come here at all
    assert(false);
    return -1;
}

__global__ void jpadg(int* graph, int* rho, int* C)
{
    int u = (blockDim.x * blockIdx.x)+threadIdx.x;
    int minC = -1;
    // A very important change to make while loop to if. We just give one chance if doesn't get color
    // give another chance again later, no point in giving a chance again immediately
    if(C[u] == -1)
    {
        minC = getColor(u);
        if(minC != -1)
        {
            C[u] = minC;
        }
    }
    return;
}