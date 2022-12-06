/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#ifndef ISING_H
#define ISING_H

#include "random.h"
#include <string>

class Ising
{
public:
  Random rnd;
  
  //parameters, observables
  constexpr static int m_props=1000;
  int n_props,iu,ic,im,ix,ig;
  double nbins;
  double walker[m_props];
  
  // averages
  double blk_av[m_props];
  double blk_norm;
  
  double attempted,accepted;
  
  double glob_av[m_props];
  double glob_av2[m_props];
  double stima_u,stima_c,stima_m,stima_x,stima_g;
  double err_u,err_c,err_m,err_x,err_g;
  
  //configuration
  constexpr static int m_spin=50;
  double s[m_spin];
  
  //restart
  bool rest;
  
  // thermodynamical state
  int nspin;
  double beta,temp,J,h,ntherm;
  double temp_max;
  double temp_min;
  double temp_step;
  
  // simulation
  int nstep, nblk, metro;
   
  //functions
  Ising(bool rest, double temp_max, double temp_min);
  void input(std::string, std::string);
  void restart(std::string);
  void move(int);
  void thermalization();
  void simulation(std::string);
  void reset(int);
  void accumulate();
  void averages(int);
  void average_final(std::string);
  void confFinal(std::string);
  void measure();
  double boltzmann(int);
  int pbc(int);
  double error(double,double,int);

};

#endif //ISING_H

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
