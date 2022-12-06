/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Exercise 4_2
//We perform molecular dynamics simulations via Lennard Jones model of Argon for 
//three phases: solid , liquid , gas and we calculate energy per particle, 
//pot energy per particle, kin energy per particle, temperature and pressure.
//

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>
#include "mdmc.h"

int main(int argc, char* argv[])
{ 
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
    //int nconf = 1;
    for (int iblk = 1; iblk <= mdmc.nblk; ++iblk)
    {
      mdmc.reset(iblk);
      for (int istep = 1; istep <= mdmc.nstep; ++istep)
      {
        mdmc.move(); 
        mdmc.measure();
        mdmc.accumulate();
        //if(istep%100 == 0)
        //{
        //  mdmc.confXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        //  nconf += 1;
        //}
      }
      mdmc.averages(iblk, std::string{argv[2]});
    }
    mdmc.confFinal("data/config.out", "data/velocity.out");
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
