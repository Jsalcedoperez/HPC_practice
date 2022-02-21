#include <string>
#include <vector>
//#include "Student_info.h"

//preposcessor variable

#ifndef aux_files_CORE_h
#define aux_files_CORE_h

class Core

{

public:
friend class Student_info;

Core(): midterm(0), final(0) {}
Core(std::stringstream& ss) { read(ss);}


std::string get_rank() const;

double get_median() const;

std::string get_name() const;

virtual void read(std::stringstream& ss);

virtual double get_grade() const;

protected:
virtual Core* clone() const {return new Core(*this);}
std::string name;
std::string rank;
double midterm, final;
std::vector<double> homework;

};

#endif
