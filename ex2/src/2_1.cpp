/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// Exercise 2_1
// Integral estimated via Monte Carlo sampling from two different distributions:
// (to visualize integral see jupyter notebook in ex2 directory)
// - Uniform distribution
// - using importance sampling (non uniform distribution)
//


#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "random.h"

int main(int argc, char* argv[])
{ 
  constexpr int N_blocks = 100;
  constexpr int N_throws = 10000;

  Random rnd(SEED "/Primes", SEED "/seed.in");

  std::ofstream out("data/N_data_2_1.txt", std::ofstream::out);
  if (!out.is_open())
  {
    std::cerr << "Can't open N_data_2_1.txt" << std::endl;
    return 1;
  }
  out << N_blocks << std::endl;
  out << N_throws << std::endl;
  out.close();

  out.open("data/2_1_uniform.txt");
  if (!out.is_open())
  {
    std::cerr << "Can't open 2_1_uniform.txt" << std::endl;
    return 1;
  }

  //uniform distribution 
  for (int i = 0; i < N_blocks; ++i)
  {
    double sum = 0.0;
    double sum_var = 0.0;

    for (int j = 0; j < N_throws; ++j)
    {
      double x = rnd.rannyu(0., 1.); //see slide lection 2: p(x) = 1
      double f_x = M_PI/2 * std::cos(M_PI/2 * x); //see slide lection 2: g(x)

      sum += f_x;
      sum_var += (f_x - 1.0) * (f_x * 1.0);
    }

    double ave_sum = sum / N_throws;
    double ave_var = sum_var / N_throws;

    out << ave_sum << '\t' << ave_var << std::endl;
  }
  out.close();

  out.open("data/2_1_importance.txt");
  if (!out.is_open())
  {
    std::cerr << "Can't open 2_1_importance.txt" << std::endl;
    return 1;
  }

  //importance sampling
  for (int i = 0; i < N_blocks; ++i)
  {
    double sum = 0.0;
    double sum_var = 0.0;
    
    for (int j = 0; j < N_throws; ++j)
    {
      double x = rnd.rannyu();
      x = 1 + std::sqrt(1 - x);  //inverse of the CDF of 2 * (1 - x) with p(y) = 1
      double g_x = M_PI/2 * std::cos(M_PI/2 * x) / (2 - 2*x);

      sum += g_x;
      sum_var += (g_x - 1.0) * (g_x - 1.0);
    }

    double ave_sum = sum / N_throws;
    double ave_var = sum_var / N_throws;
 
    out << ave_sum << '\t' << ave_var << std::endl;
  }
  out.close();

  rnd.save_seed(SEED "/seed.out");
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
