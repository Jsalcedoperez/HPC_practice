#include "../aux_files/Core.h"

#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <ios>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

void Core::read(std::stringstream& ss)
{
std::fstream newfile;
      int count = 0;
      while (ss.good())                      
      { 
        ++count;

        std::string substr;
        getline(ss,substr,',');
        
        this->rank = "UGraduate";
        
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

        else

        {
            throw std::runtime_error("Not the right number of inputs. Input has to include a rank, name, and the following grades: midterm, final and three hws.");
        }

      }        
}

std::string Core::get_rank() const

{

  return this->rank;

}

std::string Core::get_name() const

{

  return this->name;

}

double Core::get_median() const
{

typedef std::vector<double>::size_type vec_sz;
vec_sz size = this->homework.size();

if (size == 0) 
{

   throw std::domain_error("median of an empty vector");

}

std::vector<double> vec(size);

std::copy(this->homework.begin(),this->homework.end(),vec.begin());

std::sort(vec.begin(),vec.end());

vec_sz mid = size / 2;

return size % 2 == 0 ? (vec[mid] + vec[mid-1]) / 2 : vec[mid];

}

double Core::get_grade() const
{

return (0.2 * this->midterm + 0.4 * this->final + 0.4 * this->get_median());

}


