#include "Core.h"

#include <string>
#include <vector>
#include <stdexcept>
//preposcessor variable

#ifndef aux_files_Student_info_h
#define aux_files_Student_info_h

class Student_info

{

public:

// constructors and copy control

Student_info(): cp(0) {}
Student_info(std::stringstream& ss): cp(0) {read(ss);}
Student_info(const Student_info&);
Student_info& operator=(const Student_info&);
~Student_info() {delete cp;}

// operations 

void read(std::stringstream& ss);

std::string rank() const 
  {
  
  if (cp) {return cp->get_rank();}
  else { throw std::runtime_error("unitialized student");}

  }

std::string name() const 
  {
  
  if (cp) {return cp->get_name();}
  else { throw std::runtime_error("unitialized student");}

  }
double median() const 
  {
  
  if (cp) {return cp->get_median();}
  else { throw std::runtime_error("unitialized student");}

  }

double grade() const 
{
  if (cp) {return cp->get_grade();}
  else{ throw std::runtime_error("unitialized student");}
}

static bool compare(const Student_info& s1, const Student_info& s2)
{
  return s1.name() < s2.name(); 
}

private:
Core* cp;

};
#endif
