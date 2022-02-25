#include "XS.hh"

#include <fstream>
#include <sstream>
#include <vector>

using std::string;
using std::endl; using std::streamsize;
using Vec_Dbl = std::vector<double>;

class Nuclide
{

public:

  Nuclide(std::stringstream& ss);

  ~Nuclide() {};

  void read(std::stringstream& ss);

  int get_A() const;

  std::string get_name() const;

  double get_rho() const;

  double get_microXS() const;
  
  double compute_macroXS();
  
  Vec_Dbl d_compute_macroXS();

  double get_macroXS() const;

  Vec_Dbl get_Xs_energy() const;
  
  bool get_Edependency() const;

  Vec_Dbl Xs_energy;

  Vec_Dbl d_microXS;
  
  Vec_Dbl d_macroXS;
 
  void load_XS();
 
private:
  
  bool Edependency;
  std::string name;
  std::string filename;
  int A;
  XS xs_object;
  double rho; // g/cm3
  //XS microXS;
  //double macroXS;
  int total_size;
};
