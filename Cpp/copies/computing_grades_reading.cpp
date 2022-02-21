#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <ios>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "student_info.h"


using std::cin; using std::setprecision;
using std::cout; using std::string;
using std::endl; using std::streamsize;

int main()
{

//Student_info record;
typedef std::vector<Student_info>::size_type student_type;
std::vector<Student_info> students;
std::string::size_type maxlen = 0;

//read and store all the records, and find the length of the longest.

read(students);

//alphabetize the records
std::sort(students.begin(), students.end(), compare);

for (std::vector<Student_info>::size_type i = 0; i != students.size(); ++i)
{

// write the name, padded on the right to maxlen + 1 characters
// std::string(maxlen + 1 - students[i].name.size(), " "); 


try 
{  
    double med = median(students[i].homework);

    double final_grade = grade(students[i].midterm, students[i].fin, med);
   
    streamsize prec = cout.precision(3);

    std::cout << "your course grade is: " << final_grade << setprecision(prec) << std::endl;
   
}

catch (std::domain_error)
{
     
     std::cout << std::endl << " You must enter yout grades. "
                               " Please try again." << std::endl;
     return 1;
     
}

}

return 0;

}

