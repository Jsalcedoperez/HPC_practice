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
           this->xs_org = std::stod(substr);
        }

        else if (count == 4)
        {
           this->datadir = substr;
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

  //https://stackoverflow.com/questions/10093891/opening-a-file-in-c-outside-of-the-working-directory

  std::string f_MT = this->datadir + "/" + this->get_name() + "_MTs.txt";

  //std::ifstream input(f);
  //load nuclide's MTs

  load_MTs(f_MT);

  //fill myxs
  Vec_int::const_iterator MTs_iter;

  // read cross-section

  for (MTs_iter = this->MTs.begin(); MTs_iter != this->MTs.end(); ++MTs_iter)

  {
    std::string xs_file = this->datadir + "/" + this->get_name() + "_"+ std::to_string(*MTs_iter) + ".txt";
    this->load_XS(xs_file);
  }

}

void Nuclide::load_MTs(std::string filename)

{

    std::fstream file;
    file.open(filename, std::ios::in); //open a file to perform read operation using file object

    if (file.is_open())
      { //checking whether the file is open
            std::string MT;
            std::cout << " READING MT FILE DATA:: " << std::endl;
            int counter = -1;
            while(getline(file, MT))
            { //read data from file object and put it into string.
              // Returns first token
              //std::stringstream s_stream(MTs);

              this->MTs.push_back(std::stoi(MT));
              counter += 1;
              this->rx_map.push_back(counter);
            }
      }

}


void Nuclide::load_XS(std::string xs_file)

{
    // NEED to make sure i am writing this correctly.

    XS_uptrs xs = std::make_shared<XS>(xs_file);
    this->XS_.push_back(xs);
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

Vec_Dbl Nuclide::get_microXS(int MT)

{
  int indx_ = this->MTs[MT];
  return this->XS_[indx_]->get_microXS();

}


Vec_Dbl Nuclide::get_energy(int MT)

{
  int indx_ = this->MTs[MT];
  return this->XS_[indx_]->get_energy();

}


Vec_Dbl Nuclide::compute_macroXS(int MT)

{


    int indx_ = this->MTs[MT];
    auto microxs = this->XS_[indx_]->get_microXS();
    //auto microxs =  this->microXS[MT].get_microXS();

    Vec_Dbl::const_iterator microXS_iter;

    Vec_Dbl::const_iterator microXS_iter_end = microxs.end();

    for (microXS_iter = microxs.begin() ; microXS_iter != microXS_iter_end; ++microXS_iter)
    {

         auto macro_xs = this->rho * 6.022E23 * *(microXS_iter) * 1.E-24 / this->A;
         this->macroXS.push_back(macro_xs);
        //std::cout << "... " << i_macroXS << std::endl;

    }

  return this->macroXS;
}






