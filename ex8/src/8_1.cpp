/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Exercise 8_1
//Variational montecarlo to compute the expectation value of the hamiltonian
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "random.h"

//trial function
double trial_wfunc(double x, double mu, double sigma);
//derivative 
double trial_wfuncD2(double x, double mu, double sigma);
double potential_energy(double x);
double energy(double x, double mu, double sigma);

int main(int argc, char* argv[])
{ 
  constexpr int N_blocks = 100;
  constexpr int N_throws = 100000;

  Random rnd(SEED "/Primes", SEED "/seed.in");

  std::ofstream out("data/N_data_8_1.txt", std::ofstream::out);
  if (!out.is_open())
  {
    std::cerr << "Unable to open N_data_8_1.txt" << std::endl;
    return 1;
  }
  out << N_blocks << std::endl;
  out << N_throws << std::endl;

  std::ofstream data_out("data/8_1_H_expect.txt");
  if (!data_out.is_open())
  {
    std::cerr << "Unable to open 8_1_H_expect.txt" << std::endl;
    return 1;
  }
  
  std::vector<double> x = {1.2};
  constexpr double step = 2.35;
  double mu = 0.5;
  double sigma = 1.0;

  //square modulus of the trial function
  auto trial_wfunc2 = [mu, sigma](std::vector<double> coord, int dim)
  {
    double trial = trial_wfunc(coord[0], mu, sigma);         
    return trial * trial;
  };

  int accept = 0;
  for (int i = 0; i < N_blocks; ++i)
  {

    double hamiltonian = 0.0;
    for (int j = 0; j < N_throws; ++j)
    {
      accept += rnd.metropolis_uniform(x, 1, step, trial_wfunc2);
      hamiltonian += energy(x[0], mu, sigma); 
    }

    double average_ham = hamiltonian / N_throws;
    data_out << average_ham << std::endl;
  } 
  data_out.close();
  out << static_cast<double>(accept) / (N_blocks * N_throws) << std::endl;

  out.close();
  return 0;
}

double trial_wfunc(double x, double mu, double sigma)
{
  double mins = (x - mu) / sigma;
  double plus = (x + mu) / sigma;
  return std::exp(-0.5 * mins * mins) + std::exp(-0.5 * plus * plus);
}

double trial_wfuncD2(double x, double mu, double sigma)
{ 
  double mins = (x - mu) / sigma;
  double plus = (x + mu) / sigma;
  double add_1 = x*x + mu*mu - 2.0*x*mu - sigma*sigma;
  double add_2 = x*x + mu*mu + 2.0*x*mu - sigma*sigma;
  double sigma_4 = sigma*sigma*sigma*sigma;
  return (add_1 * std::exp(-0.5 * mins * mins) + add_2 * std::exp(-0.5 * plus * plus)) / sigma_4;
}

double potential_energy(double x)
{
  return x*x*x*x - 2.5*x*x;
}

double energy(double x, double mu, double sigma)
{
  return -0.5 * trial_wfuncD2(x, mu, sigma) / trial_wfunc(x, mu, sigma) + potential_energy(x); 
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
