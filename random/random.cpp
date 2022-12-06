/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include "random.h"

Random::Random(std::string primes_in, std::string seed_in)
{
  int seed[4];
  int p1, p2;
  std::ifstream primes(primes_in);

  if (primes.is_open())
  {
    primes >> p1 >> p2 ;
  } 
  else 
  {
   std::cerr << "PROBLEM: Unable to open " << primes_in << std::endl;
   exit(1);
  }
  primes.close();

  std::ifstream input(seed_in);
  std::string property;
  if (input.is_open())
  {
    while ( !input.eof() )
    {
      input >> property;
      if( property == "RANDOMSEED" )
      {
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        set_random(seed, p1, p2);
      }
    }
    input.close();
  } 
  else 
  {
    std::cerr << "PROBLEM: Unable to open " << seed_in << std::endl;
    exit(1);
  }

}

void Random::save_seed(std::string seed_out)
{
  std::ofstream write_seed;
  write_seed.open(seed_out);
  if (write_seed.is_open())
  {
    write_seed << l1 << " " << l2 << " " << l3 << " " << l4 << std::endl;;
  } 
  else 
  {
    std::cerr << "PROBLEM: Unable to open " << seed_out << std::endl;
  }
  write_seed.close();
  return;
}

double Random::gauss(double mean, double sigma) 
{
  double s=rannyu();
  double t=rannyu();
  double x=std::sqrt(-2.*std::log(1.-s))*std::cos(2.*M_PI*t);
  return mean + x * sigma;
}

double Random::rannyu(double min, double max)
{
  return min+(max-min)*rannyu();
}

double Random::rannyu(void)
{
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random::set_random(int * s, int p1, int p2)
{
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

double Random::exp(double lambda)
{
  return (-std::log(1 - rannyu())) / lambda;
}

double Random::cauchy(double mean, double gamma)
{
  return gamma * std::tan(M_PI * (rannyu() - 0.5)) + mean;
}

bool Random::metropolis_uniform(std::vector<double>& coord, int dim, double step, 
  std::function<double(std::vector<double>, int)> prob)
{
  std::vector<double> n_coord(dim);
  for (int i = 0; i < dim; ++i)
  {
    n_coord.at(i) = rannyu(coord.at(i) - step, coord.at(i) + step);
  }
  
  bool accepted = false;
  double prob_ratio = prob(n_coord, dim) / prob(coord, dim);
  if (prob_ratio > 1.0 || rannyu() <= prob_ratio)
  {  
    for (int i = 0; i < dim; ++i)
    {
      coord.at(i) = n_coord.at(i);
    }
    accepted = true;
  }
  
  return accepted;
}

bool Random::metropolis_gaussian(std::vector<double>& coord, int dim, double sigma, 
  std::function<double(std::vector<double>, int)> prob)
{
  std::vector<double> n_coord(dim);
  for (int i = 0; i < dim; ++i)
  {
    n_coord.at(i) = gauss(coord.at(i), sigma);
  }
  
  bool accepted = false;  
  double prob_ratio = prob(n_coord, dim) / prob(coord, dim);
  if (prob_ratio > 1.0 || rannyu() <= prob_ratio)
  {
    for (int i = 0; i < dim; ++i)
    {
      coord.at(i) = n_coord.at(i);
    }
    accepted = true;
  }
  
  return accepted;
}

bool Random::simulated_annealing(double beta, double h_new, double h_old, 
  std::function<double(double, double, double)> prob)
{
  bool accepted = false;
  double weight_function = prob(beta, h_new, h_old);
  if (weight_function > 1.0 || rannyu() <= weight_function)
  {
    accepted = true;
  }
  return accepted;
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
