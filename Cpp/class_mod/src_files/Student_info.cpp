#include "../aux_files/Student_info.h"
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


Student_info::Student_info(const Student_info& s) : cp(0)
{

  if (s.cp) cp = s.cp->clone();

}


Student_info& Student_info::operator=(const Student_info& s)
{

  if (&s != this)

  {

    delete cp;
    if (s.cp)
    {
      cp = s.cp->clone();
    }
    else
    {
      cp = 0;
    }


  }


}
void Student_info::read(std::stringstream& ss)

{
      delete this->cp;
      std::string G("Graduate");
      std::string U("UGraduate");
      while (ss.good())                      
      { 

        std::string substr;
        getline(ss,substr,',');
        
        if (substr == U)
        {

          this->cp = new Core(ss);

        }

        else if (substr == G)

        {
            
          this->cp = new Grad(ss);

        }

        else

        {
          std::runtime_error("wrong student rank");
        }

      }

}
