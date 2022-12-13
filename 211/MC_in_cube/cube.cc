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
#include <unordered_map>
#include <utility>
#define NBINS 1000
#include "Nuclide.hh"
#include <numeric>
#include <functional>

using Vec_Dbl = std::vector<double>;
using int_pair = std::pair<int,int>;

namespace definitions

{
  using reverse_it = std::vector<double>::reverse_iterator;
  int A,nScat;
  int nNeutrons = 1E2;
  int e_bin;
  double rs, sum_rs, mean_rs;

  double E, microXS, sigt;
  double x, x_new;
  double y, y_new;
  double z, z_new;
  double sinz, cosx, cosy, cosz;
  double eta, phi, dis_x, dis_y, dis_z;
  double muBar, mu, mu_denom, cosxNew, cosyNew, coszNew;
  double temp1, temp2;
  double rxn_th;
  double nuclide_type, rxn_type;
  std::vector<double> flux(NBINS,0);

  //211 specific
  int cube_lengh = 50;
  int max_scatt = 0;
  int MT, plane_index;
  int counter_scattering=0;
  int count_neutron_scatt = 0;
  int count_scatt_ver = 0;
  int total_scattering = 0; // before absorption
  int total_scattering_pre = 0; // before absorption and leakage
  int count_leakage=0;
  int count_abs=0;
  int nuclide,nuclide_id;
  bool anisotropic_scattering = false;
  bool one_over_v = true;
  std::vector<int> anisotropic_flag = {1,1};

  double dist_collision_total=0.0;

  std::vector<double> energy;
  double dist;
  double dist_leakage=0.0;
  double dist_leakage_total=0.0;
  double dist_absorption_total=0.0;
  double dist_absorption=0.0;
  double dist_collission=0.0;
  double dist_collission_total=0.0;
  double mu_iso = 0.0;
  double mean_scatt,mean_dist_scatt,mean_dist_leak, mean_dist_abs, mu_mean, E_mean_abs, E_mean_leak;
  double f_abs,f_leak;
  std::vector<Nuclide> Nuclides;


  double E_abs = 0;
  double E_leak = 0;
  double Eth = 1E-5;

  unsigned long rxn_seed;
  unsigned long nuclide_seed;

  std::vector<double> origin;
  std::vector<double> direction;
  std::vector<double> coords;
  std::vector<double> micro_scatXS;
  int flag=0; //0- scattering, 1- absorption, 2- leakage
              //5- low energy

  std::vector<double> sigma_xs;
  std::vector<double> sigma;
  double Nd;

struct Plane

{

    std::vector<double> point;
    std::vector<double> normal;

};


struct Input

{
   double Eo;
   double low_u;
   double kill;
   double source_E;
   double delta_u;

};

  Input input;

  Plane XPlane;
  Plane nXPlane;
  Plane YPlane;
  Plane nYPlane;
  Plane ZPlane;
  Plane nZPlane;


}

/*

function raytrace(ray::Ray, plane::Plane)
    dist::Float64 = dot( plane.point - ray.origin, plane.normal) / dot( ray.direction, plane.normal)
    # Check if parallel
    if dist < 0 || dist == Inf
        return false, NaN
    end
    return true, dist
end

*/

//I need to write it as 'definitions::Plane'
bool distance(std::vector<double>&, std::vector<double>&, definitions::Plane&);

double interpolate(double Ex, std::vector<double>& Xs_energy, std::vector<double>& Xs_macros);

template<class It, class T>
typename std::iterator_traits<It>::difference_type lower_bound_index(
  It first, It last, const T& value);

double interpolate_adv(double Ex, Vec_Dbl& energy, Vec_Dbl& macro_XS);

float rn(unsigned long *);

void load_data(std::string);

void init();

int find_u_bin(double, definitions::Input);

int find_E_bin(double, definitions::Input);

double find_u_bin_u(int, definitions::Input);

double find_u_bin_E(int, definitions::Input);

double find_E_bin_E(int, definitions::Input);

int main(int argc, char** argv)
{


   //std::vector<std::string> MT = {"total","elastic","fission","capture"};
   std::vector<int> MT = {1,2,18,102};

   using namespace definitions;

   //reads a "character'
   auto system_id = argv[1];

   std::string nuc_file = "";

   std::vector<double> Sig_Scat(2,0);


   init();


   if (system_id == 0)

   {
        nuc_file = "system0.txt";

        //load cross-sections and nuclides.
        load_data(nuc_file);

        sigma_xs = {12.7081,0.1048,0.1400,5.00579e-6};
        sigma = {30.1468,3.9759,0.332127,0.00019};
        Sig_Scat[0] = 1.5051;
        Nd = 2 * 0.7 * 0.6022 / 18;
        Sig_Scat[1] = Nd * sigma[1]/2;
        auto& nuc = Nuclides[0];

        {
            nuc.init();
            auto iMT = MT[1];
            micro_scatXS = nuc.get_microXS(iMT);
            energy = nuc.get_energy(iMT);

        }



   }


  // Intialize Seed
  unsigned long seed = 1337 * time(NULL);

  int nNeutron_counter;


  for (int n = 1; n < nNeutrons; ++n)
  {

    flag = 0;
    count_neutron_scatt = 0;
    E = input.source_E;
    x = 0.0; y = 0.0; z = 0.0;
    direction = {0.,1.,0.};
    cosx = 0.0; cosy = 1.0; cosz = 0.0;
    double min_dist;
    while (flag==0)

    {



      count_scatt_ver += 1;
      min_dist = 1e6;

      auto microXS = interpolate_adv(E,energy,micro_scatXS);
      sigma[0] = microXS;

      if (system_id==0)
      {
        Sig_Scat[0] = sigma[0] * Nd;
        if (one_over_v == 1)
        {
            sigt=Sig_Scat[0]+Sig_Scat[1]+sqrt(0.0253/E)*0.01556;
            dist = -log(rn(&seed))/sigt;
        }
        else
        {
                    sigt=Sig_Scat[0]+Sig_Scat[1]+0.01556;
                    dist = -log(rn(&seed))/sigt;
        }
        }
      else
        {
            if (one_over_v == 1)
            {

                sigt=Sig_Scat[0]+sqrt(0.0253/E)*sigma_xs[1];
                dist = -log(rn(&seed))/sigt;
            }
            else
            {
                    sigt=Sig_Scat[0]+sigma_xs[1];
                    dist = -log(rn(&seed))/sigt;
            }
        }


      double dist = -log(rn(&seed)) / sigt;
      dis_x = dist * cosx;
      dis_y = dist * cosy;
      dis_z = dist * cosz;

      // compute the new position
      x_new = x + dis_x;
      y_new = y + dis_y;
      z_new = z + dis_z;

      std::vector<double> direction = {dis_x, dis_y, dis_z};

      // put direction in a single vector
      origin = {x, y, z};

      //coords = {x_new,y_new,z_new};
      coords = {dis_x,dis_y,dis_z};

      std::vector<Plane> planes = {XPlane, YPlane, ZPlane, nXPlane, nYPlane, nZPlane};

      for (int i= 0; i < 3 ; ++i)
      {

        if (abs(coords[i]) >= 50)
        {

            //positive side
            //check if particle left through face
            if (coords[i] > 50)
            { plane_index = i;}
            else {plane_index = i+2;}

             bool leak = distance(origin, coords,planes[plane_index]);

             if (leak)
             {

                //set the flag to leak.
                 flag = 2;
                 //check all coords that are > 50 and place them
                 //right at the boundary
                for (int j = 0;  j < 3; j++)

                {
                    if (coords[i] > 50)
                    {
                        //place particle at box boundary
                        coords[i] = 50;
                    }

                    else if (coords[i] < -50)
                    {

                        coords[i] = -50;

                    }

                }
               break;
             }


        }
      }

      //std::cout << "neutron " << n << " with flag " << flag << " and nscatt " << count_neutron_scatt << std::endl;
      std::cout << "neutron " << n << " with coords (" << x << "," << y << "," <<z<< ")" << std::endl;

      if (flag != 2)

      {
            // Intialize Seed
            rxn_seed = 7777 * time(NULL);
            rxn_type = rn(&rxn_seed);
            if (system_id == 0)
            {
                nuclide_seed = 666 * time(NULL);
                nuclide_type = rn(&nuclide_seed);
                if (one_over_v == 1)
                {


                    if (E < 1)
                    {
                        // use both oxygen and hydrogen in rxn computation
                    nuclide = 2 ;
                    }

                    else{

                    //still can't remember why I am using this ratio
                    nuclide = 10/11;

                    }
                }

                if (nuclide < nuclide_type)

                {
                    nuclide_id = 1;
                    A=16;

                }

                else
                {

                    nuclide_id = 0;
                    A=1;


                }

        }



        else

        {
                rxn_th = sigma_xs[0]/sigt;
        }

        if (one_over_v)

        {
            if (nuclide == 2)
            {
                rxn_th = Sig_Scat[0] / (Sig_Scat[0] + Sig_Scat[1] + sqrt(0.0253/E) * 0.01556);
            }

            else
            {
                rxn_th = sigma[nuclide_id] / (sigma[nuclide_id] + sqrt(0.0253/E) * sigma[nuclide_id + 2]);
            }
        }

        else
        {

            rxn_th = sigma[nuclide_id] / (sigma[nuclide_id] + sigma[nuclide_id + 2]);

        }


        eta = 2 * rn(&seed) - 1;
        phi = 2 * M_PI * rn(&seed);
        if ( anisotropic_scattering)
            {

                mu_denom = sqrt(pow(A,2)+2*A*eta+1);
                mu = (1.0 + A*eta) / mu_denom;

                temp1 = sqrt(1-mu*mu);
                temp2 = sqrt(1 - cosz*cosz);
                cosxNew = mu*cosx + temp1*(cosx*cosz*cos(phi)-cosy*sin(phi)) / temp2;
                cosyNew = mu*cosy + temp1*(cosy*cosz*cos(phi)+cosx*sin(phi)) / temp2;
                coszNew = mu*cosz -temp1*temp2*cos(phi);


                nScat += 1;

            }


            else

            {

                cosxNew = sqrt(1-pow(eta,2))*cos(phi);
                cosyNew = sqrt(1-pow(eta,2))*sin(phi);
                coszNew = eta;


            }
            e_bin = find_E_bin(E,input);

            flux[e_bin] += 1 / sigt;


            /*


           reverse_it group_reverse = group.rbegin();

            for ( group_reverse ; group_reverse < group.rend() ; group_reverse++)
            {
                if (E <= Etop(*group_reverse))
                {
                    group = *group_reverse;
                    break
                }
            }

            flux(group) += 1 / sigt;
            */

            //absorption rxn
            if (rxn_th < rxn_type)
            {
                flag = 1;
                // why is it starting from the origin and no the pre-collision
                // distance
                dist_absorption = (pow(coords[0],2) + pow(coords[1],2) + pow(coords[2],2));
                count_abs+=1;
                E_abs+=1;


            }

            else

            {
                counter_scattering+=1;
                dist_collission=sqrt(pow((coords[0]-x),2)+pow((coords[1]-y),2)+pow((coords[2]-z),2));
                dist_collission_total+=dist_collission;

                mu_iso += cosxNew*cosx + cosyNew*cosy + coszNew*cosz;

                cosx = cosxNew;
                cosy = cosyNew;
                cosz = coszNew;
            }
      }

      else

      {

            count_leakage+=1;
            dist_leakage = (pow(coords[0],2) + pow(coords[1],2) + pow(coords[2],2));
            E_leak += E;

      }

       E = E * (pow(A,2) + 1 + 2*A*eta) / pow((A+1),2);
       x = coords[0];
       y = coords[1];
       z = coords[2];

       if (E <= Eth)

       {
            flag=5;

       }

        if (flag==0)

        {
            if (count_scatt_ver > 1)
            {

                count_neutron_scatt += 1;


            }

            if (counter_scattering > max_scatt)
            {
                max_scatt=counter_scattering;

            }
            total_scattering+=counter_scattering;
            dist_absorption_total += dist_absorption;
        }
        total_scattering_pre+=counter_scattering;
        if (flag == 2)
        {
            dist_leakage_total+=dist_leakage;

        }
    }
      std::cout << "neutron " << n << " with flag " << flag << " and nscatt " << count_neutron_scatt << std::endl;
  }
mu_mean = mu_iso/total_scattering_pre;
E_mean_abs = E_abs / count_abs;
E_mean_leak = E_leak / count_leakage;
f_abs = count_abs / nNeutrons;
f_leak = 1-f_abs;
mean_dist_leak = dist_leakage_total/count_leakage;
mean_dist_abs = dist_absorption_total/count_abs;
mean_scatt = total_scattering/count_neutron_scatt;
mean_dist_scatt = dist_collission_total/total_scattering_pre;


std::cout << " mean absorption energy " << E_mean_abs << std::endl;
std::cout << " mean leak energy " << E_mean_leak << std::endl;
std::cout << " mean leak dist " << mean_dist_leak << std::endl;
std::cout << " mean abs dist " << mean_dist_abs << std::endl;
std::cout << " mean scat dist " << mean_dist_scatt << std::endl;


return 0;

}

void init()

{

  using namespace definitions;


  XPlane.point = {50,0,0};
  nXPlane.point = {-50,0,0};
  YPlane.point = {0,50,0};
  nYPlane.point = {0,-50,0};
  ZPlane.point = {0,0,50};
  nZPlane.point = {0,0,-50};


  XPlane.normal = {1,0,0};
  nXPlane.normal = {-1,0,0};
  YPlane.normal = {0,1,0};
  nYPlane.normal = {0,-1,0};
  ZPlane.normal = {0,0,1};
  nZPlane.normal = {0,0,-1};

  input.kill = 1.E-7;
  input.Eo= 20.E+6;
  input.source_E = 1.E+6;
  input.low_u = log(input.Eo/input.source_E);
  input.delta_u = log(input.Eo/input.kill) - input.low_u;

}


int find_u_bin(double E, definitions::Input input)
{

    double u = log(input.Eo/E);
    double val = (u - input.low_u) / (input.delta_u);
    int bin = val * NBINS;
    return bin;
}


int find_E_bin(double E, definitions::Input input)

{

    return (E / input.source_E)  * NBINS;

}

double find_u_bin_u(int i, definitions::Input input)

{

    double del = (input.delta_u) / NBINS;
    double bin = input.low_u + i*del;
    return bin;

}

double find_u_bin_E(int i,definitions::Input input)

{
    double u = find_u_bin_u(i,input);
    double E = input.Eo / exp(u);
    return E;
}

double find_E_bin_E(int i, definitions::Input input)
{
    return (input.source_E/NBINS) * i;
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


double interpolate_adv(double Ex, Vec_Dbl& energy, Vec_Dbl& XS_macro)
{

  double xs = 0.0;
  if (Ex >= energy.front())

  {

    auto i = lower_bound_index(energy.begin(), energy.end(), Ex);
    if (i > 0)
    {
    xs = (XS_macro[i] - XS_macro[i-1]) * (Ex -energy[i-1]) / ( energy[i] - energy[i-1])  + XS_macro[i-1];
    }

  }

return xs;
}

//! Perform binary search

template<class It, class T>
typename std::iterator_traits<It>::difference_type lower_bound_index(
  It first, It last, const T& value)
{
  if (*first == value)
    return 0;
  return std::lower_bound(first, last, value) - first ;
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


bool distance(std::vector<double>& origin, std::vector<double>& direction, definitions::Plane& plane)

{

    std::vector<double> point_to_origin_dist;

    std::transform(plane.point.begin(), plane.point.end(), origin.begin(), point_to_origin_dist.begin(), std::minus<int>());

    double dist = std::inner_product(point_to_origin_dist.begin(), point_to_origin_dist.end(),plane.normal.begin(),0)/std::inner_product(plane.point.begin(),plane.point.end(),plane.normal.begin(),0);

    if ((dist == 0) or (dist > 1000000000))
    {
            return false;

    }

    else

    {
            return true;

    }

}



void load_data(std::string data_file)

{

    using namespace definitions;


   std::fstream newfile;
    //will I need to add "using namespace definition" or something along these
    //lines?
   //code snippet taken from: https://www.tutorialspoint.com/parsing-a-comma-delimited-std-string-in-cplusplus
   newfile.open(data_file,std::ios::in); //open a file to perform read operation using file object
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

}
