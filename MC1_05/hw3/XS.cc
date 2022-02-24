#include "XS.hh"

#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <ios>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>

XS::XS(std::string& ss){ load_xs(ss);}


void XS::load_data(std::string& filename)

{
    std::fstream file;
    int column, row;
    std::string substr;
    file.open(filename,std::ios::in); //open a file to perform read operation using file object

    if (file.is_open())
    { //checking whether the file is open
          std::string points;
          std::cout << " READING XS FILE DATA:: " << this->filename << std::endl;
          while(getline(file, points))
          { //read data from file object and put it into string.
            // Returns first token
            std::stringstream s_stream(points);
            row++;
            column = -1;
            while (s_stream.good())
            {
                substr = "";
                getline(s_stream,substr,',');
                column++;
                if (column % 2 == 0)
                {
                  //std::cout << "( " << column << ", " << std::stod(substr) << ")" << std::endl; 
                  this->energy.push_back(std::stod(substr));
                
                }
                else

                {
                  this->microXS.push_back(std::stod(substr));
                } 
            }
          }
    }
 
}


double XS::get_microXS() const

{

  return this->microXS;

}

double XS::get_macroXS() const

{

  return this->macroXS;

}

Vec_Dbl XS::get_energy() const

{

  return this->energy;

}

Vec_Dbl XS::compute_macroXS() 

{
  Vec_Dbl::const_iterator microXS_iter = this->d_microXS.begin();
  for (microXS_iter; microXS_iter != this->d_microXS.end() ; ++microXS_iter)
  {
    auto i_macroXS = this->rho * 6.022E23 * *(microXS_iter) * 1.E-24 / this->A;
    //std::cout << "... " << i_macroXS << std::endl;
    this->macroXS.push_back(i_macroXS);
  }
  return this->macroXS; 
}

