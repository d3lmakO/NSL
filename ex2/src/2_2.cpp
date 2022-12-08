/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Exercise 2_2
//3D random walk on a cubic lattice and in the continuum.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <array>
#include "random.h"

int main(int argc, char* argv[])
{ 
  constexpr int N_blocks = 100;
  constexpr int N_throws = 10000; //number of walks per block
  constexpr int N_steps = 100;

  Random rnd(SEED "/Primes", SEED "/seed.in");

  std::ofstream out("data/N_data_2_2.txt", std::ofstream::out);
  if (!out.is_open())
  {
    std::cerr << "Can't open N_data_2_2.txt" << std::endl;
    return 1;
  }
  out << N_blocks << std::endl;
  out << N_throws << std::endl;
  out << N_steps << std::endl;
  out.close();

  out.open("data/2_2_discreteRW.txt");
  if (!out.is_open())
  {
    std::cerr << "Can't open 2_2_discreteRW.txt" << std::endl;
    return 1;
  }

  for (int i = 0; i < N_blocks; ++i)
  {
    std::array<double, N_steps> rw2_dist = {0.0};

    for (int j = 0; j < N_throws; ++j)
    {
      std::array<int, 3> axes = {0};
      
      for (int k = 0; k < N_steps; ++k)
      {
        int ax = static_cast<int>(rnd.rannyu(0.0, 3.0));
        double forward_backward = rnd.rannyu(-1.0, 1.0);
        axes.at(ax) += (forward_backward >= 0.0 ? 1 : -1);
        rw2_dist.at(k) += axes.at(0)*axes.at(0) + axes.at(1)*axes.at(1) + axes.at(2)*axes.at(2);
      }
    }
    
    for (int n = 0; n < N_steps; ++n)
    {
      out << std::sqrt(rw2_dist.at(n) / N_throws) << '\t';
    }
    out << std::endl;
  } 
  out.close();

  out.open("data/2_2_continuumRW.txt");
  if (!out.is_open())
  {
    std::cerr << "Can't open 2_2_continuumRWRW.txt" << std::endl;
    return 1;
  }
  
  for (int i = 0; i < N_blocks; ++i)
  {
    std::array<double, N_steps> rw2_dist = {0.0};

    for (int j = 0; j < N_throws; ++j)
    {
      std::array<double, 3> axes = {0}; //x,y,z in spherical coordinates with a = 1
      
      for (int k = 0; k < N_steps; ++k)
      {
        double phi = rnd.rannyu(0.0, 2*M_PI); 
        double theta = std::acos(1 - 2 * rnd.rannyu());
        axes.at(0) += std::sin(theta) * std::cos(phi);
        axes.at(1) += std::sin(theta) * std::sin(phi);
        axes.at(2) += std::cos(theta); 
        rw2_dist.at(k) += axes.at(0)*axes.at(0) + axes.at(1)*axes.at(1) + axes.at(2)*axes.at(2);
      }      
    }

    for (int n = 0; n < N_steps; ++n)
    {
      out << std::sqrt(rw2_dist.at(n) / N_throws) << '\t';
    }
    out << std::endl;
  }
  out.close();

  return 0;
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
