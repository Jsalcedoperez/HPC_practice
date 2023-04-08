#ifndef Nuclide_hh
#define Nuclide_hh

#include "XS.hh"

#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <memory>

using std::string;
using std::endl; using std::streamsize;
using Vec_Dbl = std::vector<double>;
using Vec_str = std::vector<std::string>;
using Vec_int = std::vector<int>;
using XS_type = std::unordered_map<std::string,XS>;
using XS_uptrs = std::shared_ptr<XS>;
using XS_vec = std::vector<XS_uptrs>;
class Nuclide
{

public:

  Nuclide(std::stringstream& ss);

  ~Nuclide() {};

  void read(std::stringstream& ss);

  int get_A() const;

  std::string get_name() const;

  double get_rho() const;

  Vec_Dbl get_microXS(int);

  Vec_Dbl compute_macroXS(int);

  Vec_Dbl get_energy(int);

  void init();

  void load_MTs(std::string);

  void load_XS(std::string);

  XS_vec XS_;

  Vec_int rx_list;


private:

  std::string name;
  Vec_int MTs = Vec_int(102,0);
  std::string datadir;
  std::string xs_name;
  int A;
  double rho; // g/cm3
  double xs_org;
  //XS microXS;
  //double macroXS;
  int total_size;
  Vec_Dbl macroXS;
};

#endif
