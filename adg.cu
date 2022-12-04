#include <iostream>
#include <math.h>
#include <algorithm>
#include <cassert>
#include <string>
#include <fstream>
#include <vector>

using namespace std;
using std::min;
using std::max;



bool notAllVerticesColored(int* C, int n, int &count)
{
    bool result = false;
    count = n;
    for(int i = 1; i <= n; i++)
    {
        if(C[i] == 0)
        {
            // return true;
            result = true;
            count--;
        }
    }
    return result;
}

bool checkValidColoring(int* graph, int* adjList, int* C, int n)
{
    cout << "coloring" << endl;
    for(int i=1;i<=n;i++)
    {
        cout << i << " " << C[i] << endl;
        for(int j=graph[i];j<graph[i+1];j++)
        {
            if(C[adjList[j]] == C[i])
            {
                return false;
            }
        }
    }
    return true;
}

void parseInput(char* inputFile, int &n, int &m, int* &graph, int* &adjList, int &D)
{
    // TODO aditya: malloc the graph here after finding the number of edges and length of array needed
    // Determine the number of colors D needed
    fstream f;
    f.open(inputFile, ios::in);
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
    D++;
    graph[n + 1] = ctr;
    // for(int i = 0; i <2*m ; i++)
    // {
    //     cout << adjList[i] <<" ";
    // }
    // cout << endl;
    // for(int i = 0; i < n+2; i ++)
    // {
    //     cout << graph[i] <<" ";
    // }
    // cout << endl;

}
int* getRho(int* graph, int* adjList, int strategy, int n)
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
}

__device__ int getColor(int* graph, int* adjList, int* rho, int* C, int v, int D)
{
    bool* B = new bool[D+1]();
    memset(B,0,sizeof(bool)*(D+1));
    // if(v==n) TODO make sure n+1 th entry should be the end index of the array to make sure this works.
    // This is very important as we do not want if statement here as this function gets used a lot of times
    for(int i = graph[v]; i < graph[v+1]; i++)
    {
        if(rho[adjList[i]] > rho[v])
        {
            if(C[adjList[i]] == 0)
            {
                free(B);
                return 0;
            }
            else
            {
                B[C[adjList[i]]] = true;
            }
        }
    }
    // printf("get color %d\n",v);
    for(int i =1 ;i <= D; i++)
    {
        if(!B[i])
        {
            free(B);
            return i;
        }
    }
    // Should not come here at all
    assert(false);
    free(B);
    return 0;
}

__global__ void jpadg(int* graph, int* adjList, int* rho, int* C, int D, int n)
{
    int u = (blockDim.x * blockIdx.x)+threadIdx.x+1;
    if(u> n || u<1)
    {
        return;
    }
    int minC = 0;
    // A very important change to make while loop to if. We just give one chance if doesn't get color
    // give another chance again later, no point in giving a chance again immediately
    if(C[u] == 0)
    {
        minC = getColor(graph, adjList, rho, C, u, D);
        if(minC != 0)
        {
            C[u] = minC;
        }
    }
    return;
}


int main(int argc, char** argv)
{
    if(argc == 1)
    {
        cout << "No input" << endl;
        return 0;
    }
    // TODO by aditya
    int n, m, D;
    int *adjList = NULL; //This is the adjacency list
    int *graph = NULL;
    cout << "Parse inp" << endl;
    parseInput(argv[1], n , m, graph, adjList, D);
    cout << "Parse input D" << D << " n " << n << " m " << m << endl;
    int* rho = getRho(graph, adjList, 1, n); // 1= random order or largest degree first
    cout << "Get Rho" << endl;
    dim3 gridDim((n+1023)/1024,1,1);
    dim3 blockDim(1024,1,1);

    int* C = new int[n+1]();
    memset(C, 0, sizeof(int)*(n+1));

    int *d_graph, *d_adjList, *d_rho, *d_C;

    if(cudaMalloc(&d_graph,sizeof(int)*(n+2))!=cudaSuccess)
    {
        cout << "Could not allocate d_graph" << endl;
    }

    if(cudaMalloc(&d_adjList, sizeof(int)*(2*m))!=cudaSuccess)
    {
        cout << "Could not allocate d_adjList" << endl;
    }

    if(cudaMalloc(&d_rho,sizeof(int)*(n+1))!=cudaSuccess)
    {
        cout << "Could not allocate d_graph" << endl;
    }

    if(cudaMalloc(&d_C, sizeof(int)*(n+1))!=cudaSuccess)
    {
        cout << "Could not allocate d_C" << endl;
    }

    if(cudaMemcpy(d_graph,graph,sizeof(int)*(n+2),cudaMemcpyHostToDevice) != cudaSuccess){
        cout<<"Could not copy graph into d_graph"<<endl;
    }

    if(cudaMemcpy(d_adjList,adjList,sizeof(int)*(2*m),cudaMemcpyHostToDevice) != cudaSuccess){
        cout<<"Could not copy adjList into d_adjList"<<endl;
    }

    if(cudaMemcpy(d_rho,rho,sizeof(int)*(n+1),cudaMemcpyHostToDevice) != cudaSuccess){
        cout<<"Could not copy rho into d_rho"<<endl;
    }

    if(cudaMemset(d_C, 0, sizeof(int)*(n+1)) != cudaSuccess)
    {
        cout << "Could not memset C" << endl;
    }

    int iter = 0;
    int count = 0;

    while(notAllVerticesColored(C,n,count))
    {
        // We need not run again for all vertices
        // Run only for uncolored vertices VERY IMPORTANT
        cout << "Running iteration " << iter++ << " colored : " << count << "/" << n << endl;
        jpadg<<<gridDim, blockDim>>>(d_graph, d_adjList, d_rho, d_C, D, n);
        auto code = cudaMemcpy(C,d_C,sizeof(int)*(n+1),cudaMemcpyDeviceToHost);
        if (code != cudaSuccess)
        {
            cout << "GPUassert:" << cudaGetErrorName(code) << " " <<  cudaGetErrorString(code) << " " << endl;
        }
        // for(int i = 1;i<=n;i++)
        // {
        //     cout << "color of " << i << " " << C[i] << endl;
        // }
    }

    cudaFree(d_graph);
    cudaFree(d_adjList);
    cudaFree(d_rho);
    cudaFree(d_C);

    free(rho);
    assert(checkValidColoring(graph, adjList, C, n));
    free(graph);
    free(adjList);
    free(C);
    return 0;
}