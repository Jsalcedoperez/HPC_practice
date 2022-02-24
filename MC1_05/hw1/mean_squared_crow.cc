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

namespace definitions

{

  int A = 1; 
  int nNeutrons = 1E6;
  float rho = 0.7;
  double rs = 0.0, sum_rs = 0.0, mean_rs = 0.0;
  int H_microxs = 20;
  double H_numdens = rho * 6.022E23 / A;
  double H_macroxs = H_numdens * H_microxs * 1.E-24;

  double E;
  double x, x_new;
  double y, y_new;
  double z, z_new;
  double sinz, cosx, cosy, cosz;  
  double eta, phi, dis_x, dis_y, dis_z;

}

float rn(unsigned long *);

int main()
{

	// Intialize Seed
	unsigned long seed = 1337 * time(NULL);
  
  using namespace definitions;

  int nNeutron_counter = nNeutrons;

  while (nNeutron_counter >= 1)

  {

    E = 1.E+6;
    x = 0.0; y = 0.0; z = 0.0;
    cosx = 0.0; cosy = 0.0; cosz = 1.0;

    while (E > 1) 
    
    {
      // compare the length of the flight's projections
      // on three axes
      double dist = -log(rn(&seed)) / H_macroxs;
      dis_x = dist * cosx;
      dis_y = dist * cosy;
      dis_z = dist * cosz;

      // compute the new position
      x_new = x + dis_x;
      y_new = y + dis_y;
      z_new = z + dis_z;

      // decide the neutron's new flying direction 
      // (isotropic scatter in LAB system)

      eta = 2 * rn(&seed) - 1;
      phi = 2 * M_PI * rn(&seed);
      cosz = eta; 
      sinz = sqrt(1-pow(eta,2));
      cosx = sinz * cos(phi);
      cosy = sinz * sin(phi);

      // compute the new energy of neutron after scatter

      E = E * (pow(A,2) + 1 + 2*A*eta) / pow((A+1),2);
      
      // update [x,y,z] using the new position

      x = x_new;
      y = y_new;
      z = z_new;

    }

    // reduce counter.
    nNeutron_counter -= 1;
    
    // compute r-squared  (square of crow flight length)
    // for this neutron
    rs = pow(x,2) + pow(y,2) + pow(z,2);
    sum_rs += rs;
      
  }

  mean_rs = sum_rs / nNeutrons;
  double q = rn(&seed);
  std::cout << "this is the mean square crow "<< mean_rs << std::endl; 
  return 0;

}

// Park & Miller LCG from
// Numerical Recipes Vol. 2
float rn(unsigned long * seed)
{
    float ret;
    unsigned long n1;
    unsigned long a = 16807;
    unsigned long m = 2147483647;
    n1 = ( a * (*seed) ) % m;
    //std::cout << "this is n1 "<< n1 << std::endl; 
    *seed = n1;
    ret = (float) n1 / m;
    //std::cout << "this is ret "<< ret << std::endl; 
    return ret;
}
