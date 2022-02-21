#include "../aux_files/Core.h"
#include "../aux_files/Grad.h"

#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <ios>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>


void Grad::read(std::stringstream& ss)
{
std::fstream newfile;
      int count = 0;
      while (ss.good())                      
      { 
        ++count;

        std::string substr;
        getline(ss,substr,',');
        
        this->rank = "Graduate";
        
        
        if (count == 1)
        {
           this->name = substr;
        }
        
        else if (count == 2)
        {
           this->midterm = std::stod(substr);
        }
        
        else if (count == 3)
        {
           this->final = std::stod(substr);
        }
        
        else if ((count > 3) && (count < 7))
        {
           this->homework.push_back(std::stod(substr));
        }
        
        else if (count == 7)
        {
           this->thesis = std::stod(substr);
        }

        else

        {
            throw std::runtime_error("Not the right number of inputs. Input has to include a rank, name, and the following grades: midterm, final, three hws and thesis.");
        }

      }        
}


double Grad::get_grade() const
{

return (0.2 * this->midterm + 0.2 * this->final + 0.2 * this->thesis + 0.4 * this->get_median());

}
