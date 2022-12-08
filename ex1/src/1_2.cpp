/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// Exercise 1_2
// Test Central Limit Theorem for uniform and exponential distributions.
// It's different for a Cauchy - Lorentz distribution: infinite variance,
// CLT fails.
//

#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include "random.h"

int main(int argc, char* argv[])
{ 
  constexpr int N_throws = 10000;
  constexpr std::array<int, 4> N_sums = {1,2,10,100};

  Random rnd(SEED "/Primes", SEED "/seed.in");
 
  std::ofstream out("data/N_data_1_2.txt", std::ofstream::out);
  if (!out.is_open())
  {
    std::cerr << "Can't open N_data_1_2.txt" << std::endl;
    return 1;
  }
  out << N_throws << std::endl;
  for (const auto& el : N_sums)
  {
    out << el << std::endl;
  }
  out.close();
 
  //Uniform distribution  
  //////////////////////////////////////////////////////////////
  out.open("data/1_2_uniform.txt");
  if (!out.is_open())
  {
    std::cerr << "Can't open 1_2_uniform.txt" << std::endl;
    return 1;
  }

  for (int i = 0; i < N_throws; ++i)
  {
    for (const auto& sums : N_sums)
    {
      int sum = 0;
      for (int k = 0; k < sums; ++k)
      {
        sum += (int) rnd.rannyu(1., 7.);
      }
      out << (double) sum / sums << '\t';
    }    
    out << std::endl;
  }  
  out.close();

  //Exponential distribution
  //////////////////////////////////////////////////////////////
  out.open("data/1_2_exp.txt"); 
  if (!out.is_open())
  {
    std::cerr << "Can't open 1_2_exp.txt" << std::endl;
    return 1;
  }

  for (int i = 0; i < N_throws; ++i)
  {
    for (const auto& sums : N_sums)
    {
      double sum = 0.;
      for (int k = 0; k < sums; ++k)
      {
        sum += (int) rnd.exp(1.);
      }
      out << sum / sums << '\t';
    }    
    out << std::endl;
  }  
  out.close();

  //Cauchy-Lorentz distribution
  //////////////////////////////////////////////////////////////
  out.open("data/1_2_cauchy.txt"); 
  if (!out.is_open())
  {
    std::cerr << "Can't open 1_2_cauchy.txt" << std::endl;
    return 1;
  }

  for (int i = 0; i < N_throws; ++i)
  {
    for (const auto& sums : N_sums)
    {
      double sum = 0.;
      for (int k = 0; k < sums; ++k)
      {
        sum += (int) rnd.cauchy(0., 1.);
      }
      out << sum / sums << '\t';
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
