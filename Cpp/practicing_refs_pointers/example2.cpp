#include<iostream>

using namespace std;


const int& func(int& x)

{

  return x;

}

int main()
{
  //int z = 10;
  //int& y = func(z);
  //cout << y << endl;
  //in the original prob we had:
  const auto& y = func(10);
  cout << y << endl;
  // which does not work because we cannot bind a non-const l-value int reference to 
  // an r-value (i.e 10)
  return 0;
}
