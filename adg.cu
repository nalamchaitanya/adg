#include <iostream>
#include <math.h>
#include <algorithm>
#include <cassert>
#include <string>
#include<ifstream>

using namespace std;
using std::min;
using std::max;




void parseInput(string inputFile, int &n, int &m, int *adjList, int* graph, int &D)
{
    // TODO aditya: malloc the graph here after finding the number of edges and length of array needed
    // Determine the number of colors D needed
    ifstream f;
    f.open(inputFile.c_str(), ios::in);
    f >> n >> m;

    graph = new int[n + 2];
    adjList = new int[2*m];
    graph[0] = m;
    vector<int> g[n + 1];
    int u,v;
    for(int i = 0; i < m; i ++)
    {
        f >> u >> v;
        g[u].push_back(v);
        g[v].push_back(u);
    }
    D = 0;
    int ctr = 0;
    for(int i = 1; i <=n; i ++)
    {
        graph[i] = ctr;
        for(auto &x: g[i])
        {
            adjList[ctr++] = x;
        }
        D = max(D, (int)g[i].size());
    }
    graph[n + 1] = ctr;
    for(int i = 0; i <=n ; i++)
    {
        cout << adjList[i] <<" ";
    }
    cout << endl;
    for(int i = 0; i < m; i ++)
    {
        cout << graph[i] <<" ";
    }
    cout << endl;

}
int* getrho(int* graph, int* adjList, int strategy, int n)
{
    // this has to give the total order permutation on the vertices.
    // based on strategy should give random order or adg order or dec-adg order
    // Should we implement this in Device using parallel programming
    // Maybe but later

    int *rho = new int[n + 1];
    rho[0] = -1;
    for(int i  = 1; i <= n; i ++)
    {
        rho[i] = i;
    }
    random_shuffle(rho + 1, rho + n + 1);
    return rho;

int main()
{
    if(argc == 1)
    {
        cout << "No input" << endl;
    }
    // TODO by aditya
    int n, m, D;
    int *adjList; //This is the adjacency list
    int *graph;
    parseInput(argv[1], n , m, adjList, graph, D);

    int* rho = getRho(graph, 1); // 1= random order or largest degree first
    
    dim3 gridDim((n+1023)/1024,1,1);
    dim3 blockDim(1024,1,1);

    int *d_graph, *d_adjList, *d_rho, *d_C;

    if(cudaMalloc(&d_graph,sizeof(int)*(n+2))!=cudaSuccess)
    {
        cout << "Could not allocate d_graph" << endl;
    }

    if(cudaMalloc(&d_adjList, sizeof(int)*(2*m))!=cudaSuccess)
    {
        cout << "Could not allocate d_adjList" << endl;
    }

    if(cudaMalloc(&d_rho,sizeof(int)*n)!=cudaSuccess)
    {
        cout << "Could not allocate d_graph" << endl;
    }

    if(cudaMalloc(&d_C, sizeof(int)*n)!=cudaSuccess)
    {
        cout << "Could not allocate d_C" << endl;
    }

    if(cudaMemcpy(d_graph,graph,sizeof(int)*(n+2),cudaMemcpyHostToDevice) != cudaSuccess){
        cout<<"Could not copy graph into d_graph"<<endl;
    }

    if(cudaMemcpy(d_adjList,adjList,sizeof(int)*(2*m),cudaMemcpyHostToDevice) != cudaSuccess){
        cout<<"Could not copy adjList into d_adjList"<<endl;
    }

    if(cudaMemcpy(d_rho,rho,sizeof(int)*(n),cudaMemcpyHostToDevice) != cudaSuccess){
        cout<<"Could not copy rho into d_rho"<<endl;
    }

    if(cudaMemset(&d_C, 0, sizeof(int)*n) != cudaSuccess)
    {
        cout << "Could not memset C" << endl;
    }

    while(notAllVerticesColored(C,n))
    {
        // We need not run again for all vertices
        // Run only for uncolored vertices VERY IMPORTANT
        jpadg<<<gridDim, blockDim>>>(d_graph, d_adjList, d_rho, d_C);
        if(cudaMemcpy(C,d_C,sizeof(int)*n,cudaMemcpyDeviceToHost) != cudaSuccess)
        {
            cout << "Could not copy d_C into C" << endl;
        }
    }

    assert(true);
    // assert(checkValidColoring(graph, adjList, C));
    return 0;
}

__device__ int getColor(int* graph, int* rho, int* C, int v, int D)
{
    bool B[D] = 0;
    // if(v==n) TODO make sure n+1 th entry should be the end index of the array to make sure this works.
    // This is very important as we do not want if statement here as this function gets used a lot of times
    for(int i = graph[v]; i < graph[v+1]; i++)
    {
        if(rho[adjList[i]] > rho[v])
        {
            if(C[adjList[i]] == -1)
            {
                return -1;
            }
            else
            {
                B[C[adjList[i]]] = true;
            }
        }
    }
    for(int i =1 ;i <= D; i++)
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
