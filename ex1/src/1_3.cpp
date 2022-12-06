/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

// Exercise 1_3
// Estimate pi employing Buffon's method.
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
  constexpr double needle_l = 0.5;
  constexpr double lines_dist = 1.0;

  Random rnd(SEED "/Primes", SEED "/seed.in");

  std::ofstream out("data/N_data_1_3.txt", std::ofstream::out);
  if (!out.is_open())
  {
    std::cerr << "Can't open N_data_1_3.txt" << std::endl;
    return 1;
  }
  out << N_blocks << std::endl;
  out << N_throws << std::endl;
  out.close();

  out.open("data/1_3_data.txt");
  if (!out.is_open())
  {
    std::cerr << "Can't open 1_3_data.txt" << std::endl;
    return 1;
  }

  for (int i = 0; i < N_blocks; ++i)
  {
    int n_hit = 0;
    for (int j = 0; j < N_throws; ++j)
    {
      double x1 = rnd.rannyu(0.0, lines_dist); 
      
      double x_p;
      double y_p;
      do
      {
        x_p = rnd.rannyu();
        y_p = rnd.rannyu();
      } 
      while (x_p*x_p + y_p*y_p > lines_dist);

      double theta = std::acos(x_p / std::sqrt(x_p*x_p + y_p*y_p));
      double x2 = x1 + needle_l * std::sin(theta);

      if (x2 > lines_dist) n_hit++;
    }
    
    //estimate PI
    double pi = 2 * needle_l * N_throws / (n_hit * lines_dist);

    out << pi << std::endl;  
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

