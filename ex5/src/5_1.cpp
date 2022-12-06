/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Exercise 5_1
//Sampling probability densities of two hydrogen wave functions with Metropolis algorithm.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <functional>
#include "random.h"

//probability density of the wave function psi 1,0,0
double psi100(std::vector<double> coord, int dim);

//probability density of the wave function psi 2,1,0
double psi210(std::vector<double> coord, int dim);

int main(int argc, char* argv[])
{ 
  constexpr int N_blocks = 100;
  constexpr int N_throws = 10000;
  constexpr int space_dim = 3;
  constexpr double steps_100 = 1.22;
  constexpr double steps_210 = 2.97;
  constexpr double sigma_100 = 0.76;
  constexpr double sigma_210 = 1.87;

  Random rnd(SEED "/Primes", SEED "/seed.in");

  std::ofstream out("data/N_data_5_1.txt", std::ofstream::out);
  if (!out.is_open())
  {
    std::cerr << "Unable to open N_data_5_1.txt" << std::endl;
    return 1;
  }
  out << N_blocks << std::endl;
  out << N_throws << std::endl;
  out.close();

  std::ofstream accepted("data/acc_steps5_1.txt", std::ofstream::out);
  if (!accepted.is_open())
  {
    std::cerr << "Unable to open acc_steps5_1.txt" << std::endl;
    return 1;
  }

  //psi100 with uniform distribution Metropolis algorithm
  std::ofstream out_average("data/5_1_uni_psi100_ave.txt", std::ofstream::out);
  if (!out_average.is_open())
  {
    std::cerr << "Unable to open 5_1_uni_psi100_ave.txt" << std::endl;
    return 1;
  }

  out.open("data/5_1_uni_psi100.txt");
  if (!out.is_open())
  {
    std::cerr << "Unable to open 5_1_uni_psi100.txt" << std::endl;
    return 1;
  }

  std::vector<double> uniform_coord100 = {2.24,2.24,2.24}; 
  int acc_uni100 = 0; 
  
  for (int i = 0; i < N_blocks; ++i)
  {
    double rad_100 = 0.0;
    for (int j = 0; j < N_throws; ++j)
    {
      acc_uni100 += rnd.metropolis_uniform(uniform_coord100, space_dim, steps_100, psi100); 
      rad_100 += std::sqrt(uniform_coord100.at(0)*uniform_coord100.at(0) + 
                           uniform_coord100.at(1)*uniform_coord100.at(1) + 
                           uniform_coord100.at(2)*uniform_coord100.at(2));
      out << uniform_coord100.at(0) << '\t' << uniform_coord100.at(1) << '\t' << uniform_coord100.at(2) << std::endl;
    }
    double ave_rad_100 = rad_100 / N_throws;
    out_average << ave_rad_100 << std::endl;
  }
  out.close();
  out_average.close();
  accepted << (double) acc_uni100 / (N_blocks * N_throws) << std::endl;

  //psi100 with gaussian distribution Metropolis algorithm
  out_average.open("data/5_1_gauss_psi100_ave.txt", std::ofstream::out);
  if (!out_average.is_open())
  {
    std::cerr << "Unable to open 5_1_gauss_psi100_ave.txt" << std::endl;
    return 1;
  }

  out.open("data/5_1_gauss_psi100.txt");
  if (!out.is_open())
  {
    std::cerr << "Unable to open 5_1_gauss_psi100.txt" << std::endl;
    return 1;
  }

  std::vector<double> gauss_coord100 = {2.24,2.24,2.24};
  int acc_gauss100 = 0; 
  
  for (int i = 0; i < N_blocks; ++i)
  {
    double rad_100 = 0.0;
    for (int j = 0; j < N_throws; ++j)
    {
      acc_gauss100 += rnd.metropolis_gaussian(gauss_coord100, space_dim, sigma_100, psi100); 
      rad_100 += std::sqrt(gauss_coord100.at(0)*gauss_coord100.at(0) + 
                           gauss_coord100.at(1)*gauss_coord100.at(1) + 
                           gauss_coord100.at(2)*gauss_coord100.at(2));
      out << gauss_coord100.at(0) << '\t' << gauss_coord100.at(1) << '\t' << gauss_coord100.at(2) << std::endl;
    }
    double ave_rad_100 = rad_100 / N_throws;
    out_average << ave_rad_100 << std::endl;
  }
  out.close();
  out_average.close();
  accepted << (double) acc_gauss100 / (N_blocks * N_throws) << std::endl;

  //psi210 with uniform distribution Metropolis algorithm 
  out_average.open("data/5_1_uni_psi210_ave.txt", std::ofstream::out);
  if (!out_average.is_open())
  {
    std::cerr << "Unable to open 5_1_uni_psi210_ave.txt" << std::endl;
    return 1;
  }

  out.open("data/5_1_uni_psi210.txt");
  if (!out.is_open())
  {
    std::cerr << "Unable to open 5_1_uni_psi210.txt" << std::endl;
    return 1;
  }
  
  std::vector<double> uniform_coord210 = {2.24,2.24,2.24};
  int acc_uni210 = 0; 
  
  for (int i = 0; i < N_blocks; ++i)
  {
    double rad_210 = 0.0;
    for (int j = 0; j < N_throws; ++j)
    {
      acc_uni210 += rnd.metropolis_uniform(uniform_coord210, space_dim, steps_210, psi210); 
      rad_210 += std::sqrt(uniform_coord210.at(0)*uniform_coord210.at(0) + 
                           uniform_coord210.at(1)*uniform_coord210.at(1) + 
                           uniform_coord210.at(2)*uniform_coord210.at(2));
      out << uniform_coord210.at(0) << '\t' << uniform_coord210.at(1) << '\t' << uniform_coord210.at(2) << std::endl;
    }
    double ave_rad_210 = rad_210 / N_throws;
    out_average << ave_rad_210 << std::endl;
  }
  out.close();
  out_average.close();
  accepted << (double) acc_uni210 / (N_blocks * N_throws) << std::endl; 

  //psi210 with gaussian distribution Metropolis algorithm 
  out_average.open("data/5_1_gauss_psi210_ave.txt", std::ofstream::out);
  if (!out_average.is_open())
  {
    std::cerr << "Unable to open 5_1_gauss_psi210_ave.txt" << std::endl;
    return 1;
  }

  out.open("data/5_1_gauss_psi210.txt");
  if (!out.is_open())
  {
    std::cerr << "Unable to open 5_1_gauss_psi210.txt" << std::endl;
    return 1;
  }
  
  std::vector<double> gauss_coord210 = {2.24,2.24,2.24};
  int acc_gauss210 = 0; 
  
  for (int i = 0; i < N_blocks; ++i)
  {
    double rad_210 = 0.0;
    for (int j = 0; j < N_throws; ++j)
    {
      acc_gauss210 += rnd.metropolis_gaussian(gauss_coord210, space_dim, sigma_210, psi210); 
      rad_210 += std::sqrt(gauss_coord210.at(0)*gauss_coord210.at(0) + 
                           gauss_coord210.at(1)*gauss_coord210.at(1) + 
                           gauss_coord210.at(2)*gauss_coord210.at(2));
      out << gauss_coord210.at(0) << '\t' << gauss_coord210.at(1) << '\t' << gauss_coord210.at(2) << std::endl;
    }
    double ave_rad_210 = rad_210 / N_throws;
    out_average << ave_rad_210 << std::endl;
  }
  out.close();
  out_average.close();
  accepted << (double) acc_gauss210 / (N_blocks * N_throws) << std::endl; 

  accepted.close();
  rnd.save_seed(SEED "/seed.out");
  return 0;
}

double psi100(std::vector<double> coord, int dim)
{
  double r = 0.0;
  for (const auto& i : coord) r += i*i;
  r = std::sqrt(r);

  return std::exp(-r*2.0);
}

//probability density of the wave function psi 2,1,0
double psi210(std::vector<double> coord, int dim)
{ 
  double r = 0.0;
  for (const auto& i : coord) r += i*i;
  r = std::sqrt(r);
  
  return coord.back() * coord.back() * std::exp(-r); //coord.back() takes the last coordinate in coord vector, z coordinate
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
