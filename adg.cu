#include <iostream>
#include <math.h>
#include<curand.h>
#include<curand_kernel.h>
#include <algorithm>
#include <cassert>
#include <string>
#include <fstream>
#include <vector>
#include<cstdio>
#include<map>

using namespace std;
using std::min;
using std::max;
const long scale_1 = 1e16,scale_2 = 1e15;

bool notAllVerticesOrdered(long* ordering, int n, int &count)
{
    bool result = false;
    count = n;
    for(int i = 1; i <= n; i++)
    {
        if(ordering[i] == 0)
        {
            // return true;
            // TODO return here, are we using count somewhere?
            result = true;
            count--;
        }
    }
    return result;
}

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
    f.close();
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
    return;
}

__global__ void setup_kernel(curandState *state){

  int idx = threadIdx.x+blockDim.x*blockIdx.x + 1;
  curand_init(clock64(), idx, 0, &state[idx-1]);
}

__global__ void avgDegree1(int n, int sz, int* degree, long* ordering, int* aux_degree, int* aux_active)
{
    //printf("Entered avg degree\n");
    int vert = threadIdx.x + 1;
    if(vert > n or vert < 1)
        return;
    int shift = blockDim.x;
    int su = 0;
    int su2 = 0;
    for(int i = vert; i <= sz;i += shift)
    {
        su += degree[i];
        su2 += (ordering[i] == 0)?1:0;
    }
    aux_degree[vert] = su;
    aux_active[vert] = su2;
    //printf("The sums are %d and %d\n", su, su2);
    // if(su2 == 0)
    //     *avg = 0;
    // else
    return;
}

__global__ void avgDegree2(int n, int sz, int* degree, int* ordering, double* avg)
{
    //printf("Entered avg degree\n");
    int vert = threadIdx.x + 1;
    if(vert > n or vert < 1)
        return;
    int shift = blockDim.x;
    int su = 0;
    long su2 = 0;
    for(int i = vert; i <= sz;i += shift)
    {
        su += degree[i];
        su2 += ordering[i];
    }
    //printf("The sums are %d and %d\n", su, su2);
    // if(su2 == 0)
    //     *avg = 0;
    // else
    *avg = (double)su/su2;
    return;
}

__global__ void halfSum(int limit, long* arr1, long* arr2)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x + 1;
    if(2*i <= limit)
    {
        arr2[i] = arr1[2*i-1]+ arr1[2*i];
    }
    else if(2*i == limit+1)
    {
        arr2[i] = arr1[2*i-1];
    }
    return;
}

__global__ void halfSum2(int limit, int* arr1, long* arr2)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x + 1;
    if(2*i <= limit)
    {
        arr2[i] = arr1[2*i-1]+ arr1[2*i];
    }
    else if(2*i == limit+1)
    {
        arr2[i] = arr1[2*i-1];
    }
    return;
}

long getDegSum(int n, int* degree)
{
    long* arr1;
    if(cudaMalloc(&arr1,sizeof(long)*((n+1)/2 + 1))!=cudaSuccess) //(n+1)/2 for odd case handling, +1 for 1-index
    {
        cout << "Could not allocate temp_d_degree" << endl;
    }

    int threadCount = (n+1)/2;
    dim3 gridDim((threadCount+1023)/1024, 1,1);
    dim3 blockDim(1024,1,1);
    halfSum2<<<gridDim, blockDim>>>(n,degree,arr1);
    n=(n+1)/2;
    long* arr2;
    if(cudaMalloc(&arr2,sizeof(long)*((n+1)/2 + 1))!=cudaSuccess) //(n+1)/2 for odd case handling, +1 for 1-index
    {
        cout << "Could not allocate temp_d_degree" << endl;
    }
    long* temp;
    while(n>1)
    {
        cudaDeviceSynchronize();
        threadCount = (n+1)/2;
        gridDim.x = (threadCount+1023)/1024;
        halfSum<<<gridDim, blockDim>>>(n,arr1,arr2);
        n=(n+1)/2;
        temp = arr1;
        arr1 = arr2;
        arr2 = temp;
    }
    long sum;
    auto code = cudaMemcpy(&sum,&arr1[1],sizeof(long),cudaMemcpyDeviceToHost);
    if (code != cudaSuccess)
    {
        cout << "GPU: arr1 to sum " << cudaGetErrorName(code) << " " <<  cudaGetErrorString(code) << " " << endl;
    }
    cout << "Sum : " << sum << endl;
    cudaFree(arr1);
    cudaFree(arr2);
    return sum;
}

__global__ void getADG(int n, double eps, double* avg, long* ordering, int* degree, curandState *state, int num_partition, int* temp_degree, int* graph, int* adjList)
{
    int u = blockDim.x * blockIdx.x + threadIdx.x + 1; // vertex id
    double avg_val = *avg;
    if(u > n or u < 1 or ordering[u] != 0)
    {
        return;
    }
    if(((double)degree[u]) <= avg_val*(1 + eps)) //check if need to be in current set
    {
        //Need to include this vertex in the set
        double randf = curand_uniform(&state[u-1]);
        double temp = scale_1 * num_partition + randf *(scale_2 + 0.99999);
        ordering[u] = (long) trunc(temp);
        //ordering[u] = num_partition;
        for(int i = graph[u];i<graph[u+1];i++)
        {
            atomicAdd(temp_degree+adjList[i],-1);
        }
    }
    return;
}

__global__  void updateDegree(long* ordering, int* degree, int* temp_degree, int n)
{
    int u = blockDim.x * blockIdx.x + threadIdx.x + 1; // vertex id
    if(u > n || u < 1)
    {
        return;
    }
    //already ordered
    degree[u] = (ordering[u] == 0 ? temp_degree[u] : 0);
    //Update degrees
    // for(int i = graph[u]; i < graph[u+1]; i ++)
    // {
    //     if(ordering[adjList[i]] == 0)
    //     {
    //         degree[u] += 1;
    //     }
    // }
    return;
}

__device__ int getColor(int* graph, int* adjList, long* rho, int* C, int v, int D)
{
    int deg = graph[v+1]-graph[v];
    bool* B = new bool[deg+2]();
    memset(B,0,sizeof(bool)*(deg+2));
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
                if(C[adjList[i]]<=deg+1)
                {
                    B[C[adjList[i]]] = true;
                }
            }
        }
    }
    // printf("get color %d ended\n",v);
    for(int i =1 ;i <= deg+1; i++)
    {
        if(!B[i])
        {
            // printf("get color %d for %d\n",i,v);
            free(B);
            return i;
        }
    }
    // Should not come here at all
    free(B);
    assert(false);
    return 0;
}

__global__ void jpadg(int* graph, int* adjList, long* rho, int* C, int D, int n)
{
    int u = (blockDim.x * blockIdx.x)+threadIdx.x+1;
    if(u> n || u<1 || C[u]!=0)
    {
        return;
    }
    C[u] = getColor(graph, adjList, rho, C, u, D);
    return;
}

long* getRho(int* graph, int* adjList, int strategy, int n)
{
    // this has to give the total order permutation on the vertices.
    // based on strategy should give random order or adg order or dec-adg order
    // Should we implement this in Device using parallel programming
    // Maybe but later

    long *rho = new long[n + 1];
    rho[0] = -1;
    for(int i  = 1; i <= n; i ++)
    {
        rho[i] = i;
    }   
    random_shuffle(rho + 1, rho + n + 1);
    return rho;
}

__global__ void getDegree(int* graph, int* degree, int n)
{
    int u = blockIdx.x * blockDim.x + threadIdx.x + 1;
    if(u>n || u<1)
    {
        return;
    }
    degree[u] = graph[u+1]-graph[u];
    return;
}

long* getRhoAdg(int* d_graph, int* d_adjList, int strategy, int n, double eps, curandState *d_state)
{
    int *d_degree, *d_auxdegree, *d_auxactive;
    long* d_ordering;
    dim3 gridDim((n+1023)/1024,1,1);
    dim3 blockDim(1024,1,1);

{
    if(cudaMalloc(&d_degree,sizeof(int)*(n+1))!=cudaSuccess)
    {
        cout << "Could not allocate d_degree" << endl;
    }
    if(cudaMalloc(&d_auxdegree,sizeof(int)*(n+1))!=cudaSuccess)
    {
        cout << "Could not allocate d_auxdegree" << endl;
    }
    if(cudaMalloc(&d_auxactive,sizeof(int)*(n+1))!=cudaSuccess)
    {
        cout << "Could not allocate d_auxactive" << endl;
    }

    if(cudaMalloc(&d_ordering,sizeof(long)*(n+1))!=cudaSuccess)
    {
        cout << "Could not allocate d_ordering" << endl;
    }

    if(cudaMemset(d_ordering, 0, sizeof(long)*(n+1)) != cudaSuccess)
    {
        cout << "Could not memset d_ordering" << endl;
    }

    if(cudaMemset(d_auxdegree, 0, sizeof(int)*(n+1)) != cudaSuccess)
    {
        cout << "Could not memset d_auxdegree" << endl;
    }

    if(cudaMemset(d_auxactive, 0, sizeof(int)*(n+1)) != cudaSuccess)
    {
        cout << "Could not memset d_auxactive" << endl;
    }
}
    int count = 0;
    int num_partition = 1;
    long *ordering = new long[n+1]();
    memset(ordering, 0 , sizeof(long)*(n + 1));

    double *d_avg;
    // int *d_active; //number of vertices not yet in the ordering
    auto code = cudaMalloc(&d_avg, sizeof(double));
    if (code != cudaSuccess)
    {
        cout <<"Could not malloc d_avg"<<endl;

        cout << "GPU:" << cudaGetErrorName(code) << " " <<  cudaGetErrorString(code) << " " << endl;
    }
    code = cudaMemset(d_avg, 0, sizeof(double));
    // auto code = cudaMalloc(&d_active, sizeof(int));
    if (code != cudaSuccess)
    {
        cout << "GPU:" << cudaGetErrorName(code) << " " <<  cudaGetErrorString(code) << " " << endl;
    }

    int* temp_d_degree;
    if(cudaMalloc(&temp_d_degree,sizeof(int)*(n+1))!=cudaSuccess)
    {
        cout << "Could not allocate temp_d_degree" << endl;
    }

    getDegree<<<gridDim, blockDim>>>(d_graph, temp_d_degree, n);

    while(notAllVerticesOrdered(ordering, n, count))
    {
        cout <<"Finished ordering" << count << endl;
        double *avg = new double;

        updateDegree<<<gridDim,blockDim>>>(d_ordering, d_degree, temp_d_degree, n);
        cudaDeviceSynchronize();

        avgDegree1<<<1, blockDim.x>>>(n, n,  d_degree, d_ordering, d_auxdegree, d_auxactive);
        cudaDeviceSynchronize();

        *avg = double(getDegSum(n,d_degree))/(n-count);

        // avgDegree2<<<1, 1>>>(n, min(n, blockDim.x), d_auxdegree, d_auxactive, d_avg);
        code = cudaMemcpy(d_avg, avg, sizeof(double),cudaMemcpyHostToDevice);
        if (code != cudaSuccess)
        {
            cout << "GPU: d_avg to avg " << cudaGetErrorName(code) << " " <<  cudaGetErrorString(code) << " " << endl;
        }
        cout <<"avg degree is "<<*avg << endl;
        getADG<<<gridDim, blockDim>>>(n, eps, d_avg, d_ordering, d_degree, d_state, num_partition++, temp_d_degree, d_graph, d_adjList);
        code = cudaMemcpy(ordering,d_ordering,sizeof(long)*(n+1),cudaMemcpyDeviceToHost); // copy from device to host
        if (code != cudaSuccess)
        {
            cout << "GPU d_ordering into ordering " << cudaGetErrorName(code) << " " <<  cudaGetErrorString(code) << " " << endl;
        }
    }
    cout <<"Finished ordering " << count << endl;
    return ordering;
}


int main(int argc, char** argv)
{
    if(argc == 1)
    {
        cout << "No input" << endl;
        return 0;
    }
    int n, m, D;
    int *adjList = NULL; //This is the adjacency list
    int *graph = NULL;
    cout << "Parse inp" << endl;
    parseInput(argv[1], n , m, graph, adjList, D);
    cout << "Parse input D" << D << " n " << n << " m " << m << endl;
    //long* rho = getRho(graph, adjList, 1, n); // 1= random order or largest degree first
    
    cout << "Get Rho" << endl;
    dim3 gridDim((n+1023)/1024,1,1);
    dim3 blockDim(1024,1,1);
    curandState *d_state;
    cudaMalloc(&d_state, blockDim.x * gridDim.x * sizeof(curandState));
    // dim3 gridDim(n,1,1);
    // dim3 blockDim(1,1,1);

    setup_kernel<<<gridDim,blockDim>>>(d_state);
    cout <<"finished random number generation" << endl;

    int* C = new int[n+1](); //change back to int if needed
    memset(C, 0, sizeof(int)*(n+1));

    int *d_graph, *d_adjList,*d_C;
    long* d_rho;

    if(cudaMalloc(&d_graph,sizeof(int)*(n+2))!=cudaSuccess)
    {
        cout << "Could not allocate d_graph" << endl;
    }

    if(cudaMalloc(&d_adjList, sizeof(int)*(2*m))!=cudaSuccess)
    {
        cout << "Could not allocate d_adjList" << endl;
    }

    if(cudaMalloc(&d_rho,sizeof(long)*(n+1))!=cudaSuccess)
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


    if(cudaMemset(d_C, 0, sizeof(int)*(n+1)) != cudaSuccess)
    {
        cout << "Could not memset C" << endl;
    }

    int iter = 0;
    int count = 0;
    const double eps = 0.5;

    cout <<"calling getrhoadg" << endl;
    long *rho = getRhoAdg(d_graph, d_adjList, 0, n, eps, d_state);

    //print ADG ordering
    map<long,long> mymap;
    long flag = 0;
    for(int i = 1; i <= n; i ++)
    {
       // cout <<i<< " : "<< rho[i] << endl;
        while(mymap.find(rho[i]) != mymap.end())
        {
            rho[i]++;
            flag++;
        }
        mymap[rho[i]]++;
    }

    if(cudaMemcpy(d_rho,rho,sizeof(long)*(n+1),cudaMemcpyHostToDevice) != cudaSuccess){
        cout<<"Could not copy rho into d_rho"<<endl;
    }


    cout <<"Number of collisions" << flag << endl;

    while(notAllVerticesColored(C,n,count))
    {
        // We need not run again for all vertices
        // Run only for uncolored vertices VERY IMPORTANT
        cout << "Running iteration " << iter++ << " colored : " << count << "/" << n << endl;
        jpadg<<<gridDim, blockDim>>>(d_graph, d_adjList, d_rho, d_C, D, n);
        auto code = cudaMemcpy(C,d_C,sizeof(int)*(n+1),cudaMemcpyDeviceToHost);
        if (code != cudaSuccess)
        {
            cout << "GPU:" << cudaGetErrorName(code) << " " <<  cudaGetErrorString(code) << " " << endl;
        }
        // for(int i = 1;i<=n;i++)
        // {
        //     cout << "color of " << i << " " << C[i] << endl;
        // }
    }
    cout << "Running iteration " << iter++ << " colored : " << count << "/" << n << endl;

    cudaFree(d_graph);
    cudaFree(d_adjList);
    cudaFree(d_rho);
    cudaFree(d_C);

    free(rho);
    if(!checkValidColoring(graph, adjList, C, n))
    {
        cout << "Assert failed" << endl;
    }
    free(graph);
    free(adjList);
    free(C);
    return 0;
}
