#include "student_info.h"

#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <ios>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>


void read(std::vector<Student_info>& students)
{
std::fstream newfile;
//code snippet taken from: https://www.tutorialspoint.com/parsing-a-comma-delimited-std-string-in-cplusplus
newfile.open("data.txt",std::ios::in); //open a file to perform read operation using file object
if (newfile.is_open()){ //checking whether the file is open
      std::string student_data;
      while(getline(newfile, student_data)){ //read data from file object and put it into string.
      // Returns first token  
      std::cout << student_data << "\n";
      Student_info student;
      std::stringstream s_stream(student_data);
      int count = 0;
      while (s_stream.good())                      
      { 
        ++count;

        std::string substr;
        getline(s_stream,substr,',');
        
        if (count == 1)
        {
           student.name = substr;
        }
        
       else if (count == 2)
        {
           student.midterm = std::stod(substr);
        }
        
        else if (count == 3)
        {
           student.fin = std::stod(substr);
        }
        
       else
        {
           student.homework.push_back(std::stod(substr));
        }
      }        
      students.push_back(student);
      }
      newfile.close(); //close the file object.
}

}

double median(std::vector<double>& vec)
{

typedef std::vector<double>::size_type vec_sz;
vec_sz size = vec.size();

if (size == 0) 
{

   throw std::domain_error("median of an empty vector");

}

std::sort(vec.begin(),vec.end());

vec_sz mid = size / 2;

return size % 2 == 0 ? (vec[mid] + vec[mid-1]) / 2 : vec[mid];

}

double grade(double midterm, double fin, double median)
{

return 0.2 * midterm + 0.4 * fin + 0.4 * median;

}

bool compare(const Student_info& s1, const Student_info& s2)
{

return s1.name < s2.name;

}
