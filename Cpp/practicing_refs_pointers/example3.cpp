#include<iostream>
#include<string>
using namespace std;

void swap(string& str1, string& str2)
{

  string temp = str1;
  str1 = str2;
  str2 = temp;
  cout << "from swap, str1 is " << str1 << endl;
  cout << "from swap, str2 is " << str2 << endl;

}

int main()
{

  string str1 = "GEEKS";
  string str2 = "FOR GEEKS";
  swap(str1, str2);
  cout << "str1 is " << str1 << endl;
  cout << "str2 is " << str2 << endl;
  return 0;

}
