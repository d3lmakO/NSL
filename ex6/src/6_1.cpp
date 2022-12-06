/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Exercise 6_1
//Calculate different quantities via Ising model with metropolis or gibbs sampling.
//

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "ising.h"

int main(int argc, char* argv[])
{ 
  constexpr double temp_min = 0.5;
  constexpr double temp_max = 2.0;
  if (argc != 2)
  {
    std::cout << "Use: " << argv[0] << " 0|1" << std::endl;
    std::cout << "0: NO restart old ising configuration." << std::endl;
    std::cout << "1: Restart old ising configuration." << std::endl;
  } 
  
  bool restart = std::stoi(argv[1]);
  Ising ising(restart, temp_max, temp_min);
  ising.input(IN_DIR "/input.dat", "data/config.final");
  ising.simulation("data/");

  ising.confFinal("data/config.final");
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
