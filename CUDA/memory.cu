#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <math.h>

#define N 16
#define threadsPerBlock 4
#define blocksPerGrid 2

//fill the vector v with n random floating point numbers.

void vfill(float* v, int n){
  int i;
  for (i=0; i<n;i++){
   v[i] = (float) rand() / RAND_MAX;
  }
};

//print the vector v

void vprint(float*v , int n){
  int i;
  printf("v = \n");
  for (i=0 ; i < n; i++){
    printf("%7.3f\n",v[i]);
}
  printf("\n");
};


__global__ void dot(float* a, float* b, float* partial_c){

__shared__ float cache[threadsPerBlock];

int tid = threadIdx.x + blockIdx.x * blockDim.x;
int cacheIndex = threadIdx.x;
float temp = 0;

while (tid < N){
  temp += a[tid] + b[tid];
  tid += blockDim.x * gridDim.x;
}

cache[cacheIndex] = temp;

//synchronize
__syncthreads();

int i = blockDim.x/2;

while (i != 0) {

  if (cacheIndex < i) {
      cache[cacheIndex] += cache[cacheIndex + 1];}
  __syncthreads();
  i = i/2;
}

//record the result in partial_c

if (cacheIndex == 0){
  partial_c[blockIdx.x] = cache[0];
}
};

int main(void){
  float c;
  float *a, *b, *partial_c;
  float *a_d, *b_d, *partial_c_d;
  
  a = (float*) malloc(N*sizeof(*a));
  b = (float*) malloc(N*sizeof(*b));
  partial_c = (float*) malloc(blocksPerGrid*sizeof(*partial_c));


  cudaMalloc((float**) &a_d, (N*sizeof(a_d)));
  cudaMalloc((float**) &b_d, (N*sizeof(b_d)));
  cudaMalloc((float**) &partial_c_d, (blocksPerGrid*sizeof(partial_c_d)));

  //fill v_h with N RANDOM floating point numbers
  vfill(a,N);
  vfill(b,N);
  //print v_h to the console
  vprint(a,N);
  vprint(b,N);

  cudaMemcpy(a_d,a,threadsPerBlock*sizeof(float),cudaMemcpyHostToDevice);
  cudaMemcpy(b_d,b,threadsPerBlock*sizeof(float),cudaMemcpyHostToDevice);

  dot<<<blocksPerGrid,threadsPerBlock>>>(a_d,b_d,partial_c_d);

  cudaMemcpy(partial_c,partial_c_d,blocksPerGrid*sizeof(float),cudaMemcpyDeviceToHost);
  c = 0;

  for (int i = 0; i<blocksPerGrid; i++){
      c += partial_c[i];
  }

printf("this is c: %.3f\n",c);
return 0;
}
