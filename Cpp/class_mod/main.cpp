#include "aux_files/Student_info.h"

#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <ios>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

using std::cin; using std::setprecision;
using std::cout; using std::string;
using std::endl; using std::streamsize;

int main()

{
std::vector<Student_info> students;

std::fstream newfile;
//code snippet taken from: https://www.tutorialspoint.com/parsing-a-comma-delimited-std-string-in-cplusplus
newfile.open("data.txt",std::ios::in); //open a file to perform read operation using file object
if (newfile.is_open())
{ //checking whether the file is open
      std::string student_data;
      while(getline(newfile, student_data))
      { //read data from file object and put it into string.
        // Returns first token  
        std::cout << student_data << "\n";
        Student_info record;
        std::stringstream s_stream(student_data);
        record.read(s_stream);
        students.push_back(record);
      }
newfile.close(); //close the file object.

}

//alphabetize the students records

std::sort(students.begin(), students.end(), Student_info::compare);

// write the names and the grades


for (std::vector<Student_info>::size_type i = 0; i != students.size(); ++i)
{

// write the name, padded on the right to maxlen + 1 characters
  //std::string(maxlen + 1 - students[i].name().size(), " "); 
  std::cout <<" Processing "<< students[i].rank() << " student "<< students[i].name() <<" 's grades.." << std::endl;

try 
{  
    double med = students[i].median();
}


catch (std::domain_error)
{
     
     std::cout << std::endl << " You must enter yout grades. "
                               " Please try again." << std::endl;
     return 1;
     
}
try
{

    double final_grade = students[i].grade();
   
    streamsize prec = cout.precision(3);

    std::cout << "your course grade is: " << final_grade << setprecision(prec) << std::endl;
   


}

catch (std::runtime_error)
{
     
     std::cout << std::endl << " The student is not initialized. "
                               << std::endl;
     return 1;
     
}

}
Student_info Baldocopy = students[0];

std::cout << "this is baldocopy's name, rank and grade: " << Baldocopy.name() << ", "<< Baldocopy.rank() << ", " << Baldocopy.grade() << std::endl; 

Baldocopy = students[1];

std::cout << "this is carmelos's assignment to Baldo's copy: name, rank and grade: " << Baldocopy.name() << ", "<< Baldocopy.rank() << ", " << Baldocopy.grade() << std::endl; 

return 0;

}






