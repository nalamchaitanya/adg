#include <iostream>
#include <math.h>
#include <algorithm>
#include <cassert>

using namespace std;
using std::min;
using std::max;

__device__ int IJtoIdx(int i, int j)
{
    int t_blockIdx_x = i/blockDim.x;
    int t_threadIdx_x = i % blockDim.x;
    int t_blockIdx_y = j/blockDim.y;
    int t_threadIdx_y = j % blockDim.y;
    int topCorner = (t_blockIdx_x * gridDim.y + t_blockIdx_y)*blockDim.x*blockDim.y;
    int idx = topCorner + t_threadIdx_x*blockDim.y + t_threadIdx_y;
    return idx;
}

int Host_IJtoIdx(int i, int j, int b_x, int b_y, int g_x, int g_y)
{
    int t_blockIdx_x = i/b_x;
    int t_threadIdx_x = i % b_x;
    int t_blockIdx_y = j/b_y;
    int t_threadIdx_y = j % b_y;
    int topCorner = (t_blockIdx_x * g_y + t_blockIdx_y)*b_x*b_y;
    int idx = topCorner + t_threadIdx_x*b_y + t_threadIdx_y;
    return idx;
}

__device__ double secondmin(double a, double b, double c, double d)
{
    return max(min(min(a,b), max(c,d)), min(max(a,b), min(c,d)));
}

__global__ void Stencil(double *a, double *b, int n)
{
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;
    int idx = IJtoIdx(i,j);
    if(i == 0 || i >= n-1 || j == 0 || j >= n-1)
    {
        b[idx] = a[idx];
    }
    else
    {
        // Can be optimized later by just writing all this in one function to find second min
        double top = a[IJtoIdx(i-1,j-1)];
        double bottom = a[IJtoIdx(i+1,j+1)];
        double left = a[IJtoIdx(i+1,j-1)];
        double right = a[IJtoIdx(i-1,j+1)];
        b[idx] = a[idx] + secondmin(top,bottom,left,right);
    }
    return;
}

__global__ void sumElements1(double *arr, int sz)
{
    int st = threadIdx.x;
    if(st >= sz)
        return;
    int shift = blockDim.x;
    double su = 0;
    for(int i = st; i < sz;i += shift)
    {
        su += arr[i];
    }
    arr[st] = su;
}

__global__ void sumElements2(double *arr, int num)
{
    if(num == 0)
        return;
    double su = 0;
    for(int i = 0; i < num; i++)
    {
        su += arr[i];
    }
    arr[0] = su;
}

void init(double *A, int n, int g_x, int g_y, int b_x, int b_y)
{
    double temp = 0;
    for(int i = 0;i< n;i++)
    {
        for(int j = 0;j<n;j++)
        {
            temp = 1+cos(2*i)+sin(j);
            A[Host_IJtoIdx(i,j,b_x,b_y,g_x,g_y)] = temp*temp;
        }
    }
    return;
}

int main(int argc, char **argv)
{
    //allocate memory for device variables
    double *d_A,*d_B;

    int A_x = atoi(argv[1]);
    int A_y = atoi(argv[1]);
    
    int b_x = 32;
    int b_y = 32;

    int g_x = (A_x + b_x - 1)/b_x;
    int g_y = (A_y + b_y - 1)/b_y;

    dim3 gridDim(g_x,g_y,1);
    dim3 blockDim(b_x,b_y,1);

    int A_size = g_x*g_y*b_x*b_y;

    // Allocating memory on host
    double *A = new double[A_size];
    memset(A, 0, sizeof(double)*A_size);

    //we can check if the cuda functions fail by seeing if they return a cudaSuccess code
    //you get status codes like cudaSuccess for free when you are compiling with nvcc
    if(cudaMalloc(&d_A,sizeof(double)*A_size) != cudaSuccess){
        cout<<"Could not allocate d_A"<<endl;
    }
    cudaMemset(&d_A ,0,A_size*sizeof(double));
    if(cudaMalloc(&d_B,sizeof(double)*A_size) != cudaSuccess){
        cout<<"Could not allocate d_A"<<endl;
    }
    cudaMemset(&d_B ,0,A_size*sizeof(double));

    // Init<<<gridDim,blockDim>>>(d_A, b_x*b_y, A_x);
    init(A,A_x,g_x,g_y, b_x,b_y);
    if(cudaMemcpy(d_A,A,sizeof(double)*A_size,cudaMemcpyHostToDevice) != cudaSuccess){
        cout<<"Could not copy A into d_A"<<endl;
    }

    //Start cuda timer after init
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start);



    int t = 10;
    for(int iter = 0;iter < t;iter++)
    {
        cudaDeviceSynchronize();
        if(iter%2 == 0)
        {
            Stencil<<<gridDim, blockDim>>>(d_A, d_B, A_x);
        }
        else
        {
            Stencil<<<gridDim, blockDim>>>(d_B, d_A, A_x);
        }
    }

    double sum, val;
    int val_idx;
    if(A_x >= 48)
    {
        val_idx = Host_IJtoIdx(37,47,b_x,b_y,g_x,g_y);
        if(cudaMemcpy(&val,d_A+val_idx,sizeof(double),cudaMemcpyDeviceToHost) != cudaSuccess)
        {
            cout<<"Could not copy d_A into A"<<endl;
        }
    }
    else
    {
        val_idx = -1;
    }

    cudaDeviceSynchronize();
    sumElements1<<<1,1024>>>(d_A, A_size);

    cudaDeviceSynchronize();
    sumElements1<<<1,1>>>(d_A, min(1024, A_size));

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    if(cudaMemcpy(&sum,d_A,sizeof(double),cudaMemcpyDeviceToHost) != cudaSuccess){
        cout<<"Could not copy d_A into A"<<endl;
    }

    cudaFree(d_A);
    cudaFree(d_B);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    if(val_idx == -1)
    {
        cout << "n: " << A_x << " sum: " << sum << " val: N/A time(ms): " << milliseconds << endl;
    }
    else
    {
        cout << "n: " << A_x << " sum: " << sum << " val: " << val << " time(ms): " << milliseconds << endl;
    }

    return 0;
}