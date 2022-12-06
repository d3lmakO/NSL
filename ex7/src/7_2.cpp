/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Exercise 7_2
//Print out instantaneous values of potential energy per particle
//

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>
#include "mdmc.h"

int main(int argc, char* argv[])
{ 
  constexpr int mc_steps = 500000; 

  if (argc != 3)
  {
    std::cout << "Use: " << argv[0] << " 0|1 phase" << std::endl;
    std::cout << "0: equilibration" << std::endl;
    std::cout << "1: measure" << std::endl;
    std::cout << "phase: gas|solid|liquid" << std::endl;
    std::cout << "WARNING! You need to call at least one equilibration!" << std::endl;
  }

  if (std::stoi(argv[1]) == 0)
  {
    Mdmc mdmc(std::stoi(argv[1]));
    mdmc.input(MDMC_DIR "/input.in",MDMC_DIR "/config.in", "data/");
    mdmc.equilibration("equilibration/" + std::string{argv[2]} + "_phase.txt", 
                       "data/config.out", 
                       "data/velocity.out");
  }
  else
  {
    Mdmc mdmc(std::stoi(argv[1]));
    mdmc.input(MDMC_DIR "/input.in",MDMC_DIR "/config.in", "data/");
    for (int i = 0; i < mc_steps; ++i)
    {
      mdmc.move();
      mdmc.measure();
      mdmc.out_istant_values("data/"+std::string{argv[2]}+"_istant_pot_energy.txt");
      std::cout << mdmc.accepted / mdmc.attempted << std::endl;
    }
  }

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
