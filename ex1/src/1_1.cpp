/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// Exercise 1_1
// Test the Pseudo-Random number generator estimating mean value, variance and chi squared.
//

#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include "random.h"

int main(int argc, char* argv[])
{
  constexpr int N_blocks = 100;
  constexpr int N_throws = 10000;

  //Construct RNG
  Random rnd(SEED "/Primes", SEED "/seed.in");

  std::ofstream out("data/N_data_1_1.txt", std::ofstream::out);
  if (!out.is_open())
  {
    std::cerr << "Can't open N_data_1_1.txt" << std::endl;
    return 1;
  }
  out << N_blocks << std::endl;
  out << N_throws << std::endl;
  out.close();

  out.open("data/1_1_data.txt");
  if (!out.is_open())
  {
    std::cerr << "Can't open 1_1_data.txt" << std::endl;
    return 1;
  }
 
  for (int i = 0; i < N_blocks; ++i)
  {
    double sum = 0.;
    double sum_var = 0.; 
    double chi_sq = 0.;
    std::array<int, N_blocks> m = {};

    for (int j = 0; j < N_throws; ++j)
    {
      double x = rnd.rannyu();
      
      //chi^2
      m.at((int) (x*N_blocks))++;
      
      //calculate total sum and total square sum for mean and variance
      sum += x;
      sum_var += (x - 0.5) * (x - 0.5);
    }

    //chi^2
    for (int k = 0; k < N_blocks; ++k)
    {
      chi_sq += (m.at(k) - (N_throws / N_blocks)) * (m.at(k) - (N_throws / N_blocks));
    }
    chi_sq *= ((double) N_blocks / N_throws);

    double ave_sum = sum / N_throws;
    double ave_var = sum_var / N_throws;


    out << ave_sum << '\t' << ave_var << '\t' << chi_sq << std::endl;
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
