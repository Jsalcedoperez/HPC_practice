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

#include "Nuclide.hh"

namespace definitions

{
  int A,nScat;
  int nNeutrons = 1E6;
  double rs, sum_rs, mean_rs;

  double E, macroXS;
  double x, x_new;
  double y, y_new;
  double z, z_new;
  double sinz, cosx, cosy, cosz;  
  double eta, phi, dis_x, dis_y, dis_z;
  double muBar, mu, mu_denom, cosxNew, cosyNew, coszNew;
  double temp1, temp2;
}

double interpolate(double Ex, std::vector<double>& Xs_energy ,std::vector<double>& Xs_macros);

template<class It, class T>
typename std::iterator_traits<It>::difference_type lower_bound_index(
  It first, It last, const T& value);


double interpolate_adv(double Ex, std::vector<double>& Xs_energy ,std::vector<double>& Xs_macros);

float rn(unsigned long *);

int main()
{

    std::vector<Nuclide> Nuclides;

    std::fstream newfile;
    //code snippet taken from: https://www.tutorialspoint.com/parsing-a-comma-delimited-std-string-in-cplusplus
    newfile.open("data.txt",std::ios::in); //open a file to perform read operation using file object
    if (newfile.is_open())
    { //checking whether the file is open
          std::string nuclide_data;
          while(getline(newfile, nuclide_data))
          { //read data from file object and put it into string.
            // Returns first token  
            std::cout << nuclide_data << "\n";
            std::stringstream s_stream(nuclide_data);
            Nuclide nuclide(s_stream);
            Nuclides.push_back(nuclide);
          }
    }
    newfile.close(); //close the file object.

	// Intialize Seed
	unsigned long seed = 1337 * time(NULL);
  
  using namespace definitions;

  int nNeutron_counter;

  for (auto& nuc : Nuclides)

  {
  rs = 0.0;
  sum_rs = 0.0;
  mean_rs = 0.0;
  nScat = 0;
  macroXS = nuc.compute_macroXS();
  nNeutron_counter = nNeutrons;
  A = nuc.get_A();
  std::vector<double> d_macroXS;
  std::vector<double> d_Xs_energy;
  if (nuc.get_Edependency())

    {
        nuc.load_XS();
        d_macroXS = nuc.d_compute_macroXS();
        d_Xs_energy = nuc.get_Xs_energy();
    }
  while (nNeutron_counter >= 1)

  {

    E = 1.E+6;
    x = 0.0; y = 0.0; z = 0.0;
    cosx = 0.0; cosy = 1.0; cosz = 0.0;

    while (E > 1) 
    
    {
      // compare the length of the flight's projections
      // on three axes
      if (nuc.get_Edependency())

      {
        macroXS = interpolate(E,d_Xs_energy,d_macroXS);
      }

      double dist = -log(rn(&seed)) / macroXS;
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
      mu_denom = sqrt(pow(A,2)+2*A*eta+1);
      mu = (1.0 + A*eta) / mu_denom;
      
      temp1 = sqrt(1-mu*mu);
      temp2 = sqrt(1 - cosz*cosz);
      cosxNew = mu*cosx + temp1*(cosx*cosz*cos(phi)-cosy*sin(phi)) / temp2;
      cosyNew = mu*cosy + temp1*(cosy*cosz*cos(phi)+cosx*sin(phi)) / temp2;
      coszNew = mu*cosz -temp1*temp2*cos(phi);

      muBar = muBar + cosxNew*cosx + cosyNew*cosy + coszNew*cosz;

      nScat += 1; 

      cosx = cosxNew;
      cosy = cosyNew;
      cosz = coszNew;

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
  muBar = muBar / nScat;
  std::cout << "For nuclide " << nuc.get_name() << ", this is the mean squared crow "<< mean_rs << std::endl; 
  
  }
return 0;

}

double interpolate(double Ex, std::vector<double>& Xs_energy ,std::vector<double>& Xs_macros)
{

    int i = 0;
    while ( Ex > Xs_energy[i])

    {

      i++;
      

    }
    
    double xs = (Xs_macros[i] - Xs_macros[i-1]) * (Ex -Xs_energy[i-1]) / ( Xs_energy[i] - Xs_energy[i-1])  + Xs_macros[i-1];
    return xs;
    
}

double interpolate_adv(double Ex, std::vector<double>& Xs_energy, std::vector<double>& Xs_macros)
{

  double xs = 0.0;
  
  if (Ex >= Xs_macros.front())

  {

    auto i = lower_bound_index(Xs_energy.begin(), Xs_energy.end(), Ex);
    if (i > 0)
    {
    xs = (Xs_macros[i] - Xs_macros[i-1]) * (Ex -Xs_energy[i-1]) / ( Xs_energy[i] - Xs_energy[i-1])  + Xs_macros[i-1];
    }

  }

}

//! Perform binary search

template<class It, class T>
typename std::iterator_traits<It>::difference_type lower_bound_index(
  It first, It last, const T& value)
{
  if (*first == value)
    return 0;
  return std::lower_bound(first, last, value) - first - 1;
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
