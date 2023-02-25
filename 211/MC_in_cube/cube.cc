#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <cstdlib>
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
  int nNeutrons = 120000;
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
  double nuclide,nuclide_type, rxn_type;
  std::vector<double> flux(NBINS,0);

  //211 specific
  int max_scatt = 0;
  int MT, plane_index;
  int counter_scattering;
  int count_neutron_scatt = 0;
  int count_scatt_ver = 0;
  double total_scattering = 0.0; // before absorption
  double total_scattering_pre = 0; // before absorption and leakage
  double count_leakage=0.0;
  double count_abs=0.0;
  int nuclide_id;
  bool anisotropic_scattering = true;
  bool one_over_v = true;
  bool interpolate_xs = false;
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


  double E_abs = 0.0;
  double E_leak = 0.0;
  double Eth = 1e-5;


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

double rn(unsigned long *);

void load_data(std::string);

void init();

int find_u_bin(double, definitions::Input);

int find_E_bin(double, definitions::Input);

double find_u_bin_u(int, definitions::Input);

double find_u_bin_E(int, definitions::Input);

double find_E_bin_E(int, definitions::Input);

double tot_scat, tot_abs;

int main(int argc, char** argv)
{


   //std::vector<std::string> MT = {"total","elastic","fission","capture"};
   std::vector<int> MT = {1,2,18,102};

   using namespace definitions;
   using namespace std;
   //reads a "character'
   //int system_id = static_cast<int>(argv[1]);

   auto system_id = argv[1];
   std::string nuc_file = "";

   std::vector<double> Sig_Scat(2,0);
   std::vector<double> Sig_Abs(2,0);


   init();


   if (*system_id=='0')

   {
        nuc_file = "system0.txt";

        //load cross-sections and nuclides.
        load_data(nuc_file);

        //sigma_xs = {12.7081,0.1048,0.1400,5.00579e-6};
        //sigma_xs = {1.4120,0.09311,0.01556,4.450e-6};
        sigma = {30.1468,3.9759,0.332127,0.00019};
        //sigma = {{30.0,3.9,0.36,1.79e-4};
        Nd = 2 * 0.7 * 0.6022 / 18;
        std::vector<double> sigma_xs(4,0);
        sigma_xs[0] = Nd * sigma[0];
        sigma_xs[1] = Nd / 2 * sigma[1];
        sigma_xs[2] = Nd * sigma[2];
        sigma_xs[3] = Nd / 2 * sigma[3];

        Sig_Scat[0] = sigma_xs[0];
        Sig_Scat[1] = sigma_xs[1];
        Sig_Abs[0] = sigma_xs[2];
        Sig_Abs[1]= sigma_xs[3];

        tot_scat = std::accumulate(Sig_Scat.begin(),Sig_Scat.end(),0.0);
        tot_abs = std::accumulate(Sig_Abs.begin(),Sig_Abs.end(),0.0);

        std::cout << "Sig_Scat[0] "<< Sig_Scat[0] << std::endl;
        std::cout << "Sig_Scat[1] "<< Sig_Scat[1] << std::endl;
        std::cout << "Sig_Abs[0] "<< Sig_Abs[0] << std::endl;
        std::cout << "Sig_Abs[1] "<< Sig_Abs[1] << std::endl;
        std::cout << "tot_scat "<< tot_scat << std::endl;
        std::cout << "tot_abs "<< tot_abs << std::endl;


        auto& nuc = Nuclides[0];

        {
            nuc.init();
            auto iMT = MT[1];
            micro_scatXS = nuc.get_microXS(iMT);
            energy = nuc.get_energy(iMT);

        }



   }

  // Intialize Seed
  unsigned long seed = 1337;

  unsigned long rxn_seed = 7777 * time(NULL);

  unsigned long nuclide_seed = 2306 * time(NULL);



  int nNeutron_counter;

  std::cout << " seed " << seed << std::endl;

  for (int n = 1; n < nNeutrons; ++n)
  {

    flag = 0;
    count_scatt_ver = 0;
    counter_scattering = 0;
    E = input.source_E;
    x = 0.0; y = 0.0; z = 0.0;
    cosx = 0.0; cosy = 1.0; cosz = 0.0;
    while (flag==0)

    {

      count_scatt_ver += 1;
      //std::cout << "Energy " << E << std::endl;
      //std::cout << "reading block" << std::endl;


      if (interpolate_xs)


      {

        auto microXS = interpolate_adv(E,energy,micro_scatXS);
        sigma[0] = microXS;


      }





      Sig_Scat[0] = sigma[0] * Nd;
      std::cout << "E = " << E << std::endl;
      std::cout << "sigma[0] = " << sigma[0] << std::endl;
      std::cout << "Sig_Scat[0] = " << Sig_Scat[0] << std::endl;

      if (*system_id=='0')
      {

        //std::cout << "aqui" << std::endl;
        if (one_over_v)

        {
            auto no = 0.0253 / E;
            auto sqno = std::sqrt(no);
            std::cout << "no " << no << std::endl;
            std::cout << "sqno " << sqno << std::endl;
            sigt=tot_scat+sqno*tot_abs;
        }
        else
        {
                    sigt=tot_scat+tot_abs;
        }
      }
      else
      {
            if (one_over_v)
            {

                sigt=Sig_Scat[0]+std::sqrt(0.0253/E)*Sig_Abs[0];
            }
            else
            {
                    sigt=Sig_Scat[0]+Sig_Abs[0];
            }
      }


      double dist = -log(rn(&seed)) / sigt;
      std::cout << "dist " << dist << std::endl;
      std::cout << "sigt " << sigt << std::endl;
      dis_x = dist * cosx;
      dis_y = dist * cosy;
      dis_z = dist * cosz;

      // compute the new position
      x_new = x + dis_x;
      y_new = y + dis_y;
      z_new = z + dis_z;


      direction = {dis_x, dis_y, dis_z};
      // put direction in a single vector
      origin = {x, y, z};

      //coords = {x_new,y_new,z_new};
      coords = {x_new, y_new, z_new};

      //std::cout << "neutron " << n << " with old coords (" << x << "," << y << "," <<z<< ")" << std::endl;

      //std::cout << "neutron " << n << " with coords (" << x_new << "," << y_new << "," <<z_new<< ")" << std::endl;

      std::vector<Plane> planes = {XPlane, YPlane, ZPlane, nXPlane, nYPlane, nZPlane};

      for (int i= 0; i < 3 ; ++i)
      {

        if (abs(coords[i]) >= 50)
        {

            //std::cout << "leak block" << std::endl;
            //positive side
            //check if particle left through face
            if (coords[i] > 50)
            { plane_index = i;}
            else {plane_index = i+2;}

             bool leak = distance(origin, coords,planes[plane_index]);

             if (leak)
             {
                 std::cout << "if leak block" << std::endl;
                 //std::cout << "(" <<coords[0] << ","<<coords[1]<<","<<coords[2]<<")"<< std::endl;

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


      std::cout << "collission block" << std::endl;
      std::cout << "system id "<< *system_id << std::endl;
      if (flag != 2)

      {
            // Intialize Seed
            rxn_type = rn(&seed);
            if (*system_id=='0')
            {
                std::cout << "en flag!=0" << std::endl;
                nuclide_type = rn(&seed);
                if (one_over_v)
                {

                    //nuclide = 10.0/11.0;

                    /*
                    auto den_nuc = Sig_Scat[0] + Sig_Scat[1] + (sigma_xs[3]+sigma_xs[2])*std::sqrt(0.0253/E);

                    nuclide = (Sig_Scat[0] + sigma_xs[2]*std::sqrt(0.0253/E)) / den_nuc ;

                    */

                     auto den_nuc = sigma[0] + sigma[1] + (sigma[3]+sigma[2])*std::sqrt(0.0253/E);

                    nuclide = (sigma[0] + sigma[2]*std::sqrt(0.0253/E)) / den_nuc ;



                    if (E < 1.0)
                    {
                        // use both oxygen and hydrogen in rxn computation
                    nuclide = 2.0 ;
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

        std::cout << "nuclide_id " << nuclide_id << std::endl;
        std::cout << "nuclide " << nuclide << std::endl;
        std::cout << " nuclide_type " << nuclide_type << std::endl;

        if (one_over_v)

        {
            if (nuclide == 2.0)
            {

                auto den = 4*Sig_Scat[0] + Sig_Scat[1] + std::sqrt(0.0253/E) * tot_abs;

                //rxn_th = tot_scat/ den;

                rxn_th = ( 4*Sig_Scat[0] + Sig_Scat[1]) / den;
                std::cout << "tot_scat " << tot_scat << std::endl;
                std::cout << "den " << den << std::endl;
            }

            else
            {
                //std::cout << "nuclide id " << nuclide_id << std::endl;
                auto den = sigma[nuclide_id] + std::sqrt(0.0253/E) * sigma[nuclide_id + 2];
                rxn_th = sigma[nuclide_id] / den;

                std::cout << "A that reacted " << A << std::endl;

                std::cout << "sigma[nuclide_id] " << sigma[nuclide_id] << std::endl;
                std::cout << "den " << den << std::endl;
            }
        }

        else
        {
            //std::cout << "nuclide id " << nuclide_id << std::endl;

            rxn_th = sigma[nuclide_id] / (sigma[nuclide_id] + sigma[nuclide_id + 2]);

        }



        //std::cout << "abs block" << std::endl;
            //absorption rxn
        if (rxn_th < rxn_type)
        {
                std::cout << "desde rxn_th " << std::endl;
                std::cout << "rxn_th " << rxn_th << std::endl;
                std::cout << "rxn_type " << rxn_type << std::endl;
                flag = 1;
                // why is it starting from the origin and no the pre-collision
                // distance?
                dist_absorption = std::pow(coords[0],2) + std::pow(coords[1],2) + std::pow(coords[2],2);
                count_abs+=1;
                E_abs+=E;
                std::cout << "E_abs " << E << std::endl;

        }

        else

        {


        eta = 2 * rn(&seed) - 1;
        phi = 2 * M_PI * rn(&seed);
        if ( anisotropic_scattering)
        {

                mu_denom = std::sqrt(std::pow(A,2)+2*A*eta+1);
                mu = (1.0 + A*eta) / mu_denom;

                temp1 = sqrt(1-mu*mu);
                temp2 = sqrt(1 - cosz*cosz);
                std::cout << "mu " << mu << std::endl;
                std::cout << "x " << cosx << std::endl;
                std::cout << "y " << cosy << std::endl;
                std::cout << "z " << cosz << std::endl;
                std::cout << "phi " << phi << std::endl;
                std::cout << "temp1 " << temp1 << std::endl;
                std::cout << "temp2 " << temp2 << std::endl;

                cosxNew = mu*cosx + temp1*(cosx*cosz*std::cos(phi)-cosy*std::sin(phi)) / temp2;
                cosyNew = mu*cosy + temp1*(cosy*cosz*std::cos(phi)+cosx*std::sin(phi)) / temp2;
                coszNew = mu*cosz - temp1*temp2*std::cos(phi);

                std::cout << "cosxNew " << cosxNew << std::endl;


        }


        else

        {
                std::cout << "eta " << eta << std::endl;
                std::cout << "x " << cosx << std::endl;
                std::cout << "y " << cosy << std::endl;
                std::cout << "z " << cosz << std::endl;
                std::cout << "phi " << phi << std::endl;

                cosxNew = std::sqrt(1-std::pow(eta,2))*std::cos(phi);
                cosyNew = std::sqrt(1-std::pow(eta,2))*std::sin(phi);
                coszNew = eta;
                std::cout << "cos(phi) " << std::cos(phi) << std::endl;

                std::cout << "cosxNew " << cosxNew << std::endl;


        }


                //std::cout << "scattering block" << std::endl;
                counter_scattering+=1;
                dist_collission=std::sqrt(std::pow((coords[0]-x),2)+std::pow((coords[1]-y),2)+std::pow((coords[2]-z),2));
                dist_collission_total+=dist_collission;

                mu_iso += cosxNew*cosx + cosyNew*cosy + coszNew*cosz;

                cosx = cosxNew;
                cosy = cosyNew;
                cosz = coszNew;
                /*
               if (nuclide == 2.0)

               {
                    A=18;
               }
                */

               std::cout << "A " << A << std::endl;
               std::cout << "E1 " << E << std::endl;
               std::cout << "eta " << eta << std::endl;
               E = E * (std::pow(A,2) + 1 + 2*A*eta) / std::pow((A+1),2);
               std::cout << "E2 " << E << std::endl;
               x = coords[0];
               y = coords[1];
               z = coords[2];
               //std::cout << "reset block" << std::endl;
        }



       if (E <= input.kill)

       {
            flag=5;
            total_scattering+=counter_scattering;
            //total_scattering_pre+=counter_scattering;

       }

        if (flag==1)

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
            //total_scattering_pre+=counter_scattering;
        }
        std::cout << "counter scattering " << counter_scattering << std::endl;
        std::cout << "total scattering " << total_scattering << std::endl;


      }

      else

      {

            count_leakage+=1;
            dist_leakage_total += std::sqrt(std::pow(coords[0],2) + std::pow(coords[1],2) + std::pow(coords[2],2));
            E_leak += E;

      }
    }
        std::cout << "last block" << std::endl;
      //std::cout << "neutron " << n << " with flag " << flag << " and nscatt " << count_neutron_scatt << std::endl;
  }

std::cout << "very last block" << std::endl;
std::cout << "total scattering pre: "<< total_scattering << std::endl;
mu_mean = mu_iso/total_scattering;
double  mean_col_before_abs = total_scattering /count_abs;
E_mean_abs = E_abs / count_abs;
E_mean_leak = E_leak / count_leakage;
auto f_abs_auto = count_abs / nNeutrons;
auto f_leak_auto = 1-f_abs;
mean_dist_leak = dist_leakage_total/count_leakage;
mean_dist_abs = dist_absorption_total/count_abs;
mean_scatt = total_scattering/count_neutron_scatt;
mean_dist_scatt = dist_collission_total/total_scattering;
std::cout << "count_abs " << count_abs << std::endl;
std::cout << "count_leakage " << count_leakage << std::endl;
std::cout << "f_abs " << f_abs_auto << std::endl;
std::cout << "mean cosine lab " << mu_mean << std::endl;
std::cout << " fraction of neutrons that leaked " << f_leak << std::endl;
std::cout << " mean absorption energy " << E_mean_abs << std::endl;
std::cout << " mean leak energy " << E_mean_leak << std::endl;
std::cout << " mean leak dist " << mean_dist_leak << std::endl;
std::cout << " mean abs dist " << mean_dist_abs << std::endl;
std::cout << " mean scat dist " << mean_dist_scatt << std::endl;
std::cout << " mean # of collisions before absorption " << mean_col_before_abs <<std::endl;
std::cout << " max # collisions before absorption " << max_scatt << std::endl;

return 0;

}

void init()

{

  using namespace definitions;


  XPlane.point = {50.0,0.0,0.0};
  nXPlane.point = {-50.0,0.0,0.0};
  YPlane.point = {0.0,50.0,0.0};
  nYPlane.point = {0.0,-50.0,0.0};
  ZPlane.point = {0.0,0.0,50.0};
  nZPlane.point = {0.0,0.0,-50.0};


  XPlane.normal = {1,0,0};
  nXPlane.normal = {-1,0,0};
  YPlane.normal = {0,1,0};
  nYPlane.normal = {0,-1,0};
  ZPlane.normal = {0,0,1};
  nZPlane.normal = {0,0,-1};

  input.kill = 1.0e-10;
  input.Eo= 20.0e+6;
  input.source_E = 2.0e+6;
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
double rn(unsigned long * seed)
{
    double ret;
    unsigned long n1;
    unsigned long a = 16807;
    unsigned long m = 2147483647;
    n1 = ( a * (*seed) ) % m;
    //std::cout << "this is n1 "<< n1 << std::endl;
    *seed = n1;
    ret = (double) n1 / m;
    //std::cout << "this is ret "<< ret << std::endl;
    return ret;
}


bool distance(std::vector<double>& origin, std::vector<double>& direction, definitions::Plane& plane)

{

    std::vector<double> point_to_origin_dist(3,0);

    std::transform(plane.point.begin(), plane.point.end(), origin.begin(), point_to_origin_dist.begin(), std::minus<int>());

    std::cout << "after transform..." << std::endl;
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
