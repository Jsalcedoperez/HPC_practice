#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <iostream>

//counting iterators

//constant iterators

int main(void){


    //constant iterators

    thrust::constant_iterator<int> first(10);
    thrust::constant_iterator<int> last = first + 3;
     
    thrust::counting_iterator<int> fcount(10);
    thrust::counting_iterator<int> lastf = fcount+3;

    first[0]; //returns 10
    first[1]; //returns 10
    first[100]; //returns 10

    fcount[0]; //returns 10
    fcount[1]; //return 11
    fcount[100]; //returns 110
    
    std::cout<<"lastf[1] = "<<lastf[2]<<"and lastf[24] = "<<lastf[24]<<std::endl;    // sum of [first, last)

    //sum of [first,last)
    thrust::reduce(first,last); //returns 30 (i.e., 3*10)

    return 0;

}
