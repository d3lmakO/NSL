/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef RANDOM_H
#define RANDOM_H

#include <string>
#include <vector>
#include <functional>

class Random 
{
private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

public:
  // constructor
  Random(std::string primes_in, std::string seed_in);

  // methods
  void set_random(int * , int, int);
  void save_seed(std::string seed_out);
  double rannyu(void);
  double rannyu(double min, double max);
  double gauss(double mean, double sigma);
  double exp(double lambda);
  double cauchy(double mean, double gamma);
  bool metropolis_uniform(std::vector<double>& coord, int dim, double step, 
    std::function<double(std::vector<double>, int)> prob); 
  bool metropolis_gaussian(std::vector<double>& coord, int dim, double sigma, 
    std::function<double(std::vector<double>, int)> prob);
  bool simulated_annealing(double beta, double h_new, double h_old, 
    std::function<double(double, double, double)> prob);

};

#endif // RANDOM_H

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
