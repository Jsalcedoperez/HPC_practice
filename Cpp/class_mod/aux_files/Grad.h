#include <string>
#include <vector>
#include "Core.h"
//preposcessor variable

#ifndef aux_files_GRAD_h
#define aux_files_GRAD_h

class Grad : public Core

{

public:

Grad(): thesis(0) {};
Grad(std::stringstream& ss) { read(ss);}

virtual void read(std::stringstream& ss);

virtual double get_grade() const;

protected:
Grad* clone() const {return new Grad(*this);}

private:
double thesis;


};

#endif
