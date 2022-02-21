//pairwise-sum the elements of vector v and store the results in v[0]

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>

/* 
 *this program computes the sum of the elements of 
 *vector v using the pairwise (cascading) sum algorithm.
 */

#define N 8 //length of vector v... must be a power of 2

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

__global__ void psum(float* v){
  int t = threadIdx.x; // thread index
  int n = blockDim.x; // should be half the length of v.
                      // number of threads per block in x direction.
  while (n != 0){
    if (t < n){
       v[t] = v[t] + v[t+n];}
    __syncthreads(); 
    n = n/2;
  }
};

int main(void){
  float *v_h, *v_d; //host and device copies of our vector, respectively
  //dynamically allocate memory on the host for v_h

  v_h = (float*) malloc(N * sizeof(*v_h));
  
  // dynamically allocate memory on the device for v_d
  cudaMalloc((float**)&v_d, N * sizeof(*v_d));

  //fill v_h with N RANDOM floating point numbers
  vfill(v_h,N);
  //print v_h to the console
  vprint(v_h,N);
  // write the contents of v_h to v_d
  cudaMemcpy(v_d,v_h,N*sizeof(float),cudaMemcpyHostToDevice);
 //compute the pairwise sum of the elements of v_d and store the result in v_d[0]
 psum <<<1,N/2>>>(v_d);
 //write the pairwise sum, v_d[0], to v_h[0]
 cudaMemcpy(v_h,v_d,sizeof(float),cudaMemcpyDeviceToHost);

 printf("pairwise sum= %7.3f\n",v_h[0]);
 //free dynamically-allocated host memory
 free(v_h);
 //free dynamically-allocated device memory
 cudaFree(v_d);
}
