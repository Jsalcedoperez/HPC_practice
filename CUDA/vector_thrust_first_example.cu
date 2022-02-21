#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <iostream>

int main(void){
    // initialize all ten integers
    thrust::device_vector<int> D(10,1);

    // set the first seven elements of a vector to 9
 
    thrust::fill(D.begin(),D.begin()+7,9);

    //init host vector with the first five elements of D

    //thrust::host_vector<int>H(D.begin(),D.begin()+5);

    thrust::host_vector<int>H(5,2);

    //copy all of H back to the beginning of D

    //create a sequence

    //thrust::sequence(H.begin(), H.end());
    
    //copy all of H back to the beginning of D

    thrust::copy(H.begin(), H.end(), D.begin());

    //print D

    for(int i=0;i < D.size(); i++){
      std::cout <<"D["<<i<<"] = " << D[i] << std::endl;
     }
    return 0;

}
