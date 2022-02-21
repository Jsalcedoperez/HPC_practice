#include <string>
#include <vector>

//preposcessor variable

#ifndef GUARD_aux_files_h
#define GUARD_aux_files_h

struct Student_info
{
  std::string name;
  double midterm, fin;
  std::vector<double> homework;
};

void read(std::vector<Student_info>& students);
double median(std::vector<double>& vec);
double grade(double midterm, double fin, double median);
bool compare(const Student_info& s1, const Student_info& s2);

#endif
