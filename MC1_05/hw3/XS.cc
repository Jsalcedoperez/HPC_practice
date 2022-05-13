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

XS::XS(std::string& ss){ load_data(ss);}

void XS::load_data(std::string& filename)

{
    std::fstream file;
    int column = -1;
    std::string substr;
    file.open(filename,std::ios::in); //open a file to perform read operation using file object

    if (file.is_open())
    { //checking whether the file is open
          std::string points;
          std::cout << " READING XS FILE DATA:: "  << std::endl;
          while(getline(file, points))
          { //read data from file object and put it into string.
            // Returns first token
            std::stringstream s_stream(points);
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


Vec_Dbl XS::get_microXS()

{

  return this->microXS;

}


Vec_Dbl XS::get_energy()

{

  return this->energy;

}

