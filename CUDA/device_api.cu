#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand_kernel.h>

#define N 5

__global__ void setup_kernel(curandState *state){
       int id = threadIdx.x + blockIdx.x * N;
       
       /* Each thread gets same seed, a different sequence number, no
       offset */

       curand_init(1234, id, 0, &state[id]);
}

__global__ void generate_kernel(curandState *state, int *result){
       int id = threadIdx.x + blockIdx.x * N;
       int count = 0;
       unsigned int x;

       //Copy state to local memory for efficiency


       curandState localState = state[id];

       //Generate pseudo-random unsigned ints

       for (int n = 0; n < 100000; n++){
         x = curand(&localState);
         /* Check if odd */
         if (x & 1){
            count ++;
         }
       }

       //copy state back to global 
      state[id] = localState;

       //store results
 
       result[id] += count; 
}


int main(int argc, char *argv[]){
    int i, total;
    int *devResults, *hostResults;
    curandState *devStates;

    // Allocate space for results on host

    hostResults = (int *) calloc(N*N,sizeof(int));
   
   //Allocate space for results on device
 
   cudaMalloc((void**) &devResults,N*N*sizeof(int));
  
   //set results to 0 in the device

   cudaMemset(devResults,0,N*N*sizeof(int));
  
   //allocate space for prng states on device */

   cudaMalloc((void **) &devStates, N * N * sizeof(curandState));

   //set-up prng states

   setup_kernel<<<N,N>>>(devStates);

   /* Generate and use pseudorandom numbers*/

   for(i=0; i < 10; i++){
     generate_kernel<<<N,N>>>(devStates, devResults);
   }

   // copy device memory to host

   cudaMemcpy(hostResults, devResults, N*N*sizeof(int),cudaMemcpyDeviceToHost);

   //show results
   total = 0;
   for(i=0; i < N*N; i++){
      total += hostResults[i];
   }
   printf("fraction odd was %10.13f\n", (float) total / (5.0 * 5.0 * 100000.0 * 10.0));

   //clean-up 

   cudaFree(devStates);
   cudaFree(devResults);
   free(hostResults);

   return 0;
}
