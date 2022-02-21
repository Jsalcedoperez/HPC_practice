#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/sequence.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/replace.h>
#include <thrust/functional.h>
#include <thrust/count.h>
int main(void) {

  //allocate three device_vectors with 10 elements

  thrust::device_vector <int> X(10);
  thrust::device_vector <int> Y(10);
  thrust::device_vector <int> Z(10);

  //init X to 1,2 and 3

  thrust::sequence(X.begin(), X.end());
  
  //compute Y = -X
  

  //thrust::transform(X.begin(),X.end(),Y.begin(),thrust::negate<int>());

  //fill Z with twos

  thrust::fill(Z.begin(), Z.begin(), 2);

  //compute Y = X mod 2

  thrust::transform(X.begin(),X.end(),Z.begin(),Y.begin(),thrust::modulus<int>());

  //replace all the ones in Y with tens

  thrust::replace(Y.begin(), Y.end(), 1, 10);

  //print Y
  
  //thrust::copy(Y.begin(), Y.end(), std::ostream_iterator<int>(std::cout,"\n"));

  //doing reductions


  thrust::copy(X.begin(), X.end(), std::ostream_iterator<int>(std::cout,"\n"));

  int sum = thrust::reduce(X.begin(),X.end(), (int) 0, thrust::minus<int>());

  std::cout<<"this is sum: "<<sum<<std::endl;

  //third argument is the starting value of the reduction.

  //fourth argument is the binary operation that defines the kind of reduction.

  //counting

  // put three 1s in a device_vector

  thrust::device_vector<int> vec(5,0);

  vec[1] = 1;
  vec[3] = 1;
  vec[4] = 1;

  //count the 1s

  int result = thrust::count(vec.begin(),vec.end(),1);

  std::cout<<"this is how many 3s we have in vec: "<<result<<std:endl; 

  return 0;

}
  
 
  

