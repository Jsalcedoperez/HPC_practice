#include "Nuclide.hh"

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

Nuclide::Nuclide(std::stringstream& ss){ read(ss);}

void Nuclide::read(std::stringstream& ss)
{
      int count = -1;
      std::string substr;
      while (ss.good())                      
      { 
        ++count;
        substr = "";
        getline(ss,substr,',');
      
        std::cout << "count " << count << " and substr " << substr << std::endl;
        
        if (count == 0)
        {
           this->name = substr;
        }
        
        else if (count == 1)
        {
           this->A = std::stod(substr);
        }
        
        else if (count == 2)
        {
           this->rho = std::stod(substr);
        }
        
        else if (count == 3)
        {
           this->microXS = std::stod(substr);
        }
        
        else if (count == 4)
        {
           std::string Edependency_str = substr;
           if (Edependency_str == "True")
            {
              this->Edependency = true;

            }
            else{ this->Edependency = false;}
        }
        
        else if (count == 5)
        {
           if ( this->Edependency) 
           {
            this->datadir = substr;
           }
           else {break;}
        }  

        else

        {
            std::cout << " count is " << count  << std::endl;
            throw std::runtime_error("Not the right number of inputs. For each nuclide, the input has to include the following information:  name, atomic number, density, micro cross-section.");
        }

      }        
}

void Nuclide::init()

{

  std::unordered_map<std::string,XS> myxs;

  //load nuclide's MTs)
  
  MTs = load_MTs("data.txt");
  
  //fill myxs
  Vec_int::const_interator MTs_iter = MTs.begin(); 
  

  // read 


  for (MTs_iter; MTs_iter != MTs.end() ; ++MTs_iter)

  {
    
    myxs[*MTs_iter] = load_xs();

  }

}

void Nuclide::load_XS()

{
    
    XS::load_XS(this->datadir);
    std::fstream file;
    int column, row;
    std::string substr;
    if (this->Edependency)
    {
      file.open(filename,std::ios::in); //open a file to perform read operation using file object
      if (file.is_open())
      { //checking whether the file is open
            std::string points;
            //this->d_microXS = 0;
            //this->d_microXS = new double*[MAX_L];
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
                    Xs_energy.push_back(std::stod(substr));
                  
                  }
                  else

                  {
                    d_microXS.push_back(std::stod(substr));
                  } 
              }
            }
      }
    } 
}

std::string Nuclide::get_name() const

{

  return this->name;

}

int Nuclide::get_A() const

{

  return this->A;

}

  
double Nuclide::get_rho() const

{

  return this->rho;

}

double Nuclide::get_microXS() const

{

  return this->microXS;

}

double Nuclide::get_macroXS() const

{

  return this->macroXS;

}

Vec_Dbl Nuclide::get_Xs_energy() const

{

  return this->Xs_energy;

}

bool Nuclide::get_Edependency() const


{

  return this->Edependency;


}


double Nuclide::compute_macroXS() 

{

  return this->macroXS = this->rho * 6.022E23 * this->microXS * 1.E-24 / this->A;

}

Vec_Dbl Nuclide::d_compute_macroXS() 

{
  Vec_Dbl::const_iterator microXS_iter = this->d_microXS.begin();
  for (microXS_iter; microXS_iter != this->d_microXS.end() ; ++microXS_iter)
  {
    auto i_macroXS = this->rho * 6.022E23 * *(microXS_iter) * 1.E-24 / this->A;
    //std::cout << "... " << i_macroXS << std::endl;
    this->d_macroXS.push_back(i_macroXS);
  }
  return this->d_macroXS; 
}






