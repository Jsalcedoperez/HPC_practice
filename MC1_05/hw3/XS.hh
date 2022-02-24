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

  ~Nuclide() {};

  void load_xs(std::string& filename);

  double get_microXS() const;
  
  double compute_macroXS();
  
  Vec_Dbl d_compute_macroXS();

  double get_macroXS() const;

  Vec_Dbl get_energy() const;

 
private:
  Vec_Dbl energy;

  Vec_Dbl microXS;
  
  Vec_Dbl macroXS;
  Vec_Dbl microXS;  
};
