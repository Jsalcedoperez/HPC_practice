#include<iostream>
#include<stdio.h>
using namespace std;

void swap(char* &str1, char* &str2)
{

  std::cout << " 1 from swap, str1 is "<< std::endl;
  char* temp = str1;
  str1 = str2;
  str2 = temp;
  std::cout << "from swap, str1 is "<< std::endl;
  cout << "from swap, str2 is " << endl;
  cout << "from swap, str1 is " << str1 << endl;
  cout << "from swap, str2 is " << str2 << endl;

}

int main()
{

  const char* str1 = "GEEKS";
  const char* str2 = "FOR GEEKS";
  cout << "str1 is " << str1 << endl;
  cout << "str2 is " << str2 << endl;
  swap(str1, str2);
  cout << "str1 is " << str1 << endl;
  cout << "str2 is " << str2 << endl;
  return 0;

}
