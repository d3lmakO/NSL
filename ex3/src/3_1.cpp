/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Exercise 3_1
//Compute Europen Call-option price and European Put-option price via geometric brownian motion (GBM)
//- sampling directly the final asset price
//- sampling the discretized GBM path of the asset price 
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
  constexpr int N_steps = 100;
  //parameters
  constexpr double S0 = 100.0;
  constexpr double K = 100.0;
  constexpr double T = 1.0;
  constexpr double r = 0.1;
  constexpr double sigma = 0.25;

  Random rnd(SEED "/Primes", SEED "/seed.in");

  std::ofstream out("data/N_data_3_1.txt", std::ofstream::out);
  if (!out.is_open())
  {
    std::cerr << "Can't open N_data_3_1.txt" << std::endl;
    return 1;
  }
  out << N_blocks << std::endl;
  out << N_throws << std::endl;
  out << N_steps << std::endl;
  out.close();

  //European Call-option direct GBM
  out.open("data/3_1_call_directly.txt");
  if (!out.is_open())
  {
    std::cerr << "Can't open 3_1_call_directly.txt" << std::endl;
    return 1;
  }
  
  for (int i = 0; i < N_blocks; ++i)
  {
    double sum_gbm = 0.0;

    for (int j = 0; j < N_throws; ++j)
    {
      double wiener = rnd.gauss(0.0, T);
      double S_T = S0 * std::exp((r - sigma*sigma/2.0) * T + sigma * wiener);
      sum_gbm += S_T > K ? std::exp(-r*T) * (S_T - K) : 0;
    }

    double ave_gbm = sum_gbm / N_throws;
    out << ave_gbm << std::endl;
  }
  out.close();

  //European Call-option discretized GBM
  out.open("data/3_1_call_discretized.txt");
  if (!out.is_open())
  {
    std::cerr << "Can't open 3_1_call_discretized.txt" << std::endl;
    return 1;
  }
  
  for (int i = 0; i < N_blocks; ++i)
  {
    double sum_gbm = 0.0;

    for (int j = 0; j < N_throws; ++j)
    {
      double dt = T/N_steps;
      double sqr_dt = std::sqrt(dt);
      double S_T = S0; 
      for (int k = 0; k < N_steps; ++k)
      {
        double Z = rnd.gauss(0., 1.);
        S_T *= std::exp((r - sigma*sigma/2.0) * dt + sigma * Z * sqr_dt); 
      }
      sum_gbm += S_T > K ? std::exp(-r*T) * (S_T - K) : 0;
    }
    
    double ave_gbm = sum_gbm / N_throws;
    out << ave_gbm << std::endl;
  }
  out.close();

  //European Put-option direct GBM
  out.open("data/3_1_put_directly.txt");
  if (!out.is_open())
  {
    std::cerr << "Can't open 3_1_put_directly.txt" << std::endl;
    return 1;
  }
  
  for (int i = 0; i < N_blocks; ++i)
  {
    double sum_gbm = 0.0;

    for (int j = 0; j < N_throws; ++j)
    {
      double wiener = rnd.gauss(0.0, T);
      double S_T = S0 * std::exp((r - sigma*sigma/2.0) * T + sigma * wiener);
      sum_gbm += S_T < K ? std::exp(-r*T) * (K - S_T) : 0;
    }

    double ave_gbm = sum_gbm / N_throws;
    out << ave_gbm << std::endl;
  }
  out.close();

  //European Put-option discretized GBM
  out.open("data/3_1_put_discretized.txt");
  if (!out.is_open())
  {
    std::cerr << "Can't open 3_1_put_discretized.txt" << std::endl;
    return 1;
  }
  
  for (int i = 0; i < N_blocks; ++i)
  {
    double sum_gbm = 0.0;

    for (int j = 0; j < N_throws; ++j)
    {
      double dt = T/N_steps;
      double sqr_dt = std::sqrt(dt);
      double S_T = S0; 
      for (int k = 0; k < N_steps; ++k)
      {
        double Z = rnd.gauss(0., 1.);
        S_T *= std::exp((r - sigma*sigma/2.0) * dt + sigma * Z * sqr_dt); 
      } 
      sum_gbm += S_T < K ? std::exp(-r*T) * (K - S_T) : 0;
    }
    
    double ave_gbm = sum_gbm / N_throws;
    out << ave_gbm << std::endl;
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
