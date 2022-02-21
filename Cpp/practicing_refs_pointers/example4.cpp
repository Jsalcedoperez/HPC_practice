#include<iostream>
using namespace std;

int main()
{

  int x = 10;
  int *ptr = &x;
  //should print 10
  cout << "*ptr = "<< *ptr << endl;
  int* &ptr1 = ptr;
  //is the above equal to **ptr1 = ptr?
  //ptr1 appears to also be a  &x;
  cout << "*ptr1 = "<< *ptr1 << endl;
  return 0;
}
