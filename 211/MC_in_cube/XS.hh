#ifndef XS_class_hh
#define XS_class_hh

#include <fstream>
#include <sstream>
#include <vector>

using std::string;
using std::endl; using std::streamsize;
using Vec_Dbl = std::vector<double>;

class XS
{

public:


  XS(std::string& filename);

  ~XS() {};

  void load_data(std::string& filename);

  Vec_Dbl get_microXS();

  Vec_Dbl get_energy();


private:
  Vec_Dbl energy;

  Vec_Dbl microXS;

};

#endif
