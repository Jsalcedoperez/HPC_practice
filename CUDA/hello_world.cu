#include <stdio.h>


/*__global__ says the function is a kernel
will be executed on the GPU by one or more
simultaneous threads when called.*/

__global__ void myKernel(){
};

int main(){
 /*notes <<<nblocks,nthreadsperblock>>>*/
 myKernel<<<1, 1>>>();
 printf("Hello, World!\n");
 return 0;
}
