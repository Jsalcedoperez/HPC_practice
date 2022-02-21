#include<iostream>
using namespace std;

int main()
{
  //int z = 6;
  //int *ptr = &z;
  int* ptr = NULL;
  //references cannot be NULL
  int &ref = *ptr;
  cout << ref << endl;

}
