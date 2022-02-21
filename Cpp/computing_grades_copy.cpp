#include <algorithm>
#include <iomanip>
#include <ios>
#include <iostream>
#include <string>
#include <vector>

struct Student_info
{
  std::string name;
  double midterm, fin;
  std::vector<double> homework;
};

std::istream& read(std::istream& in, Student_info& student);
std::istream& read_hw(std::istream& in, std::vector<double>& homework);
double median(std::vector<double>& vec);
double grade(double midterm, double fin, double median);
bool compare(const Student_Info& s1, const Student_Info& s2);

using std::cin; using std::setprecision;
using std::cout; using std::string;
using std::endl; using std::streamsize;

int main()
{

Student_info record;
std::vector<Student_info>::size_type student_type;
std::vector<Student_info> students;
std::size_type maxlen = 0;

// ask for and read the student's name
std::cout << "Enter the student's name "<< std::endl;
std::cout << "Enter the student's grade in the midterm "<< std::endl;
std::cout << "Enter the student's grade in the final "<< std::endl;


//read and store all the records, and find the length of the longest.

while(read(cin,record))
{

maxlen = max(maxlen, record.name.size());
students.push_back(record);

}

//alphabetize the records
sort(students.begin(), students.end(), compare);

}

sort(students.begin(),students.end(),compare);

for (student_type i = 0; i != students.size(); ++i)
{

// write the name, padded on the right to maxlen + 1 characters
std::string(maxlen + 1 - students[i].name.size(), " "); 


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

std::istream& read(std::istream& in, Student_info& student)
{

// ask for the name
std::string s_midterm;
std::string s_fin;

in >> student.name;
in >> s_midterm;
in >> s_fin;

student.midterm = std::stod(s_midterm);
student.fin = std::stod(s_fin);

std::cout << "Hi " << student.name << ", thanks for entering the first grades" << std::endl;
std::cout << "Now, please enter your homework grades, one by one. Hit Enter after each one and type 'end' after you have all of them." << std::endl;

//invariant: we have read count grades so far, and sum is the sum
//of the first count grades

read_hw(in,student.homework);

return in;

}

std::istream& read_hw(std::istream& in, std::vector<double>& homework)
{
if (in) {
   
   //get rid of previous content
   homework.clear();
   
   //read homework grades
   double x;
    
   while (in >> x) 
   {
     homework.push_back(x);
   }
   
   //clear the stream so that input will work for the next student
   in.clear();

}   

   return in;

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

bool compare(const Student_Info& s1, const Student_Info& s2)
{

return s1.name < s2.name;

}
