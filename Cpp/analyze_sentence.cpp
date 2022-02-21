#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <ios>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>

using std::cin; using std::setprecision;
using std::cout; using std::string;
using std::endl; using std::streamsize;

int main()
{

std::fstream newfile;
//code snippet taken from: https://www.tutorialspoint.com/parsing-a-comma-delimited-std-string-in-cplusplus
newfile.open("data_sentence.txt",std::ios::in); //open a file to perform read operation using file object
if (newfile.is_open()){ //checking whether the file is open
      
      std::string student_data;

      std::map<std::string,std::vector<int>> ret;
      
      while(getline(newfile, student_data)){ //read data from file object and put it into string.
      // Returns first token  
      std::string line;
      int line_number = 0;
      
      std::cout << student_data << "\n";
    
      std::stringstream s_stream(student_data);
      
      std::vector<std::string> sentence;
      ++line_number;
      
      while (s_stream.good())                      
      { 

        std::string substr;
        getline(s_stream,substr,' ');
        
        sentence.push_back(substr);
      }

      for (std::vector<std::string>::const_iterator it = sentence.begin(); it != sentence.end(); ++it)

      {

         ret[*it].push_back(line_number);

      }
      
      }

      newfile.close(); //close the file object.
}

return 0;

}
