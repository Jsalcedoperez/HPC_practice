#include<iostream>
#include<string>
using namespace std;

void swap(int* a_ptr, int* b_ptr)
{

  int temp = *a_ptr;
  *a_ptr = *b_ptr;
  *b_ptr = temp;
  cout << "a is " << *a_ptr << endl;
  cout << "b is " << *b_ptr << endl;

}

int main()
{

  int a = 2;
  int b = 3;
  int* a_ptr = &a;
  int* b_ptr = &b;
  swap(a_ptr, b_ptr);
  cout << "a is " << *a_ptr << endl;
  cout << "b is " << *b_ptr << endl;
  return 0;

}
