/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Exercise 8_2
//Variational montecarlo simulated annealing to compute 
//mu and sigma that minimize hamiltonian expectation value.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <functional>
#include "random.h"

//trial function
double trial_wfunc(double x, double mu, double sigma);
//derivative 
double trial_wfuncD2(double x, double mu, double sigma);
double potential_energy(double x);
double energy(double x, double mu, double sigma);
double boltzmann(double beta, double h_new, double h_old);

int main(int argc, char* argv[])
{ 
  constexpr int N_blocks = 100;
  constexpr int N_throws = 10000;
  constexpr int N_throws_final = 1000000;
  constexpr double t_min = 0.1;
  constexpr double t_max = 1.0;
  constexpr double t_step = (t_max - t_min) / 5000.0;

  if (argc != 2)
  {
    std::cerr << "Use: " << argv[0] << " parameters | minimize" << std::endl; 
    return 1;
  }

  Random rnd(SEED "/Primes", SEED "/seed.in");

  std::ofstream out("data/N_data_8_2.txt", std::ofstream::out);
  if (!out.is_open())
  {
    std::cerr << "Unable to open N_data_8_2.txt" << std::endl;
    return 1;
  }
  out << N_blocks << std::endl;
  out << N_throws << std::endl;
  out << N_throws_final << std::endl;
  out << t_min << std::endl;
  out << t_max << std::endl;
  out << t_step << std::endl;
  out.close();
  //pos and step for acceptance ~50%
  std::vector<double> x = {1.2};
  constexpr double step = 2.35;
  //mu and sigma start values

  double mu = 0.5;
  double sigma = 1.0;

  if (std::string{argv[1]} == "parameters")
  {
    std::ofstream parameters("data/8_2_parameters.txt");
    if (!parameters.is_open())
    {
      std::cerr << "Unable to open 8_2_parameters.txt" << std::endl;
      return 1;
    }

    double mu_old = 0.0;
    double sigma_old = 0.0; 
    double h_old = 0.0;
    bool old = false; 
    int steps = 0;

    for (double t = t_max; t >= t_min; t -= t_step) 
    {  
      double beta = 1.0 / t;
      double h_new = 0.0;
      //square modulus of the trial function
      auto trial_wfunc2 = [mu, sigma](std::vector<double> coord, int dim)
      {
        double trial = trial_wfunc(coord[0], mu, sigma);         
        return trial * trial;
      };

      for (int i = 0; i < N_blocks; ++i)
      { 
        for (int j = 0; j < N_throws; ++j)
        {
          rnd.metropolis_uniform(x, 1, step, trial_wfunc2);
          h_new += energy(x[0], mu, sigma); 
        }

      }

      if (old == false)
      {
        h_old = h_new;
        old = true;
        x = {1.2};
      }
      else
      {
        bool accept = rnd.simulated_annealing(beta, h_new, h_old, boltzmann);
        if (accept == true) 
        {
          h_old = h_new;
          mu_old = mu;
          sigma_old = sigma;
        }
        x = {1.2};
      }

      steps++;

      parameters << h_old / (N_blocks * N_throws) << '\t' << mu_old << '\t' << mu << '\t' 
                 << sigma_old << '\t' << sigma << '\t' << steps << std::endl;

      mu = rnd.rannyu();
      sigma = rnd.rannyu();
    }


    parameters.close();
  }
  else if (std::string{argv[1]} == "minimize")
  { 
    std::ofstream minimize("data/8_2_minimize.txt");
    if (!minimize.is_open())
    {
      std::cerr << "Unable to open 8_2_minimize.txt" << std::endl;
      return 1;
    }

    std::ofstream psi2("data/8_2_psi2.txt");
    if (!psi2.is_open())
    {
      std::cerr << "Unable to open 8_2_psi2.txt" << std::endl;
      return 1;
    }

    mu = 0.796330;
    sigma = 0.615019; 
    
    auto trial_wfunc2 = [mu, sigma](std::vector<double> coord, int dim)
    {
      double trial = trial_wfunc(coord[0], mu, sigma);         
      return trial * trial;
    };

    for (int k = 0; k < N_blocks; ++k)
    { 
      double h = 0.0;
      for (int l = 0; l < N_throws_final; ++l)
      {
        rnd.metropolis_uniform(x, 1, step, trial_wfunc2);
        h += energy(x[0], mu, sigma); 
        psi2 << x[0] << std::endl;
      }
      double ave_h = h / N_throws_final;
      minimize << ave_h << std::endl;
    }
    
    psi2.close(); 
    minimize.close();
  }
  else
  {
    std::cerr << "Use: " << argv[0] << " parameters | minimize" << std::endl; 
    return 1;
  }

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

double boltzmann(double beta, double h_new, double h_old)
{
  return std::exp(-beta * (h_new - h_old));
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
