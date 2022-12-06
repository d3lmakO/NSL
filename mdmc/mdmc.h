/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef MDMC_H
#define MDMC_H

#include <string>
#include "random.h"

class Mdmc
{
public:
  Random rnd;
  std::string in_file;
  //parameters, observables
  constexpr static int m_props=1000;
  int n_props, nbins, iv, ik, it, ie, iw, ig;
  double vtail, ptail, bin_size, sd;
  double walker[m_props];
  
  // averages
  double blk_av[m_props], blk_norm, accepted, attempted;
  double glob_av[m_props], glob_av2[m_props];
  double stima_pot, stima_pres, stima_kin, stima_etot, stima_temp;
  double err_pot, err_press, err_kin, err_etot, err_temp;
  
  //configuration
  constexpr static int m_part=108;
  double x[m_part],    y[m_part],    z[m_part];
  double xold[m_part], yold[m_part], zold[m_part];
  double vx[m_part],  vy[m_part],   vz[m_part];
  
  // thermodynamical state
  int npart;
  double beta,temp,energy,vol,rho,box,rcut;
  
  // simulation
  int iNVET, nstep, nblk, restart, nequil;
  double delta;
  
  //pigreco
  constexpr static double pi=3.1415927;
  
  //functions
  Mdmc(int restart);
  void input(std::string input_file, std::string conf_file, std::string rest_dir);
  void reset(int);
  void equilibration(std::string equil_file, std::string conf_final, std::string vel_final);
  void out_istant_values(std::string path);
  void accumulate(void);
  void averages(int iblk, std::string phase);
  void move(void);
  void confFinal(std::string conf_final, std::string vel_final);
  void confXYZ(int);
  void measure(void);
  double boltzmann(double, double, double, int);
  double pbc(double);
  double error(double,double,int);
  double force(int, int);
};

#endif //MDMC_H

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
