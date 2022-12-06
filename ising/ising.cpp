/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/


//  read_input >> nspin;
//  read_input >> J;
//  read_input >> h;
//  read_input >> metro;
//  read_input >> nblk;
//  read_input >> nstep;
//  read_input >> ntherm;

#include <iostream>
#include <fstream>
#include <string>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "ising.h"

Ising::Ising(bool rest, double temp_max, double temp_min) : 
       rnd(SEED "/Primes", SEED "/seed.in"), rest{rest}, temp_max{temp_max}, temp_min{temp_min}
{ 
  temp_step = (temp_max - temp_min)/50.0; 
}

void Ising::input(std::string input_file, std::string conf_f)
{
  std::ifstream read_input;

  std::cout << "Classic 1D Ising model             " << std::endl;
  std::cout << "Monte Carlo simulation             " << std::endl << std::endl;
  std::cout << "Nearest neighbour interaction      " << std::endl << std::endl;
  std::cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << std::endl << std::endl;
  std::cout << "The program uses k_B=1 and mu_B=1 units " << std::endl;
  
  //Read input informations
  read_input.open(input_file);

  //read_input >> temp;
  std::cout << "Temperature range: 0.5 - 2.0" << std::endl;

  read_input >> nspin;
  std::cout << "Number of spins = " << nspin << std::endl;

  read_input >> J;
  std::cout << "Exchange interaction = " << J << std::endl;

  read_input >> h;
  std::cout << "External field = " << h << std::endl << std::endl;
    
  read_input >> metro; // if=1 Metropolis else Gibbs

  read_input >> nblk;

  read_input >> nstep;

  read_input >> ntherm;

  if(metro==1) std::cout << "The program perform Metropolis moves" << std::endl;
  else std::cout << "The program perform Gibbs moves" << std::endl;
  std::cout << "Number of blocks = " << nblk << std::endl;
  std::cout << "Number of steps in one block = " << nstep << std::endl << std::endl;
  std::cout << "Used steps for thermalization = " << ntherm << std::endl << std::endl;
  read_input.close();

  //Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

  if (rest == 1)
  {
    restart(conf_f);
  }
  else
  {
    //initial configuration
    for (int i=0; i<nspin; ++i)
    {
      if(rnd.rannyu() >= 0.5) s[i] = 1.0;
      else s[i] = -1.0;
    } 
  }
}

void Ising::restart(std::string conf_f)
{
  std::ifstream restart_config;
  restart_config.open(conf_f);
  int data;
  int id = 0; 
  while (restart_config >> data)
  {
    s[id] = data;
    id++;
  }
  restart_config.close();
}

void Ising::move(int metro)
{
  int o;

  for(int i=0; i<nspin; ++i)
  {
    //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.rannyu()*nspin); 
    if(metro==1) //Metropolis
    {
      attempted++;
      if (rnd.rannyu() <= boltzmann(o))
      {
        s[o] *= -1;
        accepted++;
      }
      else
      {
        accepted += 0;
      }
    }
    else //Gibbs sampling
    {
      attempted = 1.;
      accepted = 1.;
      s[o] = 1;
      if (rnd.rannyu() > (1.0 / (1.0 + boltzmann(o))))
      {
        s[o] = -1;
      }
    }
  }
}

void Ising::thermalization()
{
  for (int i = 1; i <= ntherm; ++i)
  {
    move(metro);
  }
}

void Ising::simulation(std::string averages_final)
{
  for(double t = temp_min; t <= temp_max; t = t + temp_step)
  {
    beta = 1.0/t;
    temp = t; 
    thermalization();
    std::cout << "T = " << t << std::endl;
    
    for (int iblk = 1; iblk <= nblk; ++iblk)
    {
      reset(iblk); 
      for (int istep = 1; istep <= nstep; ++istep)
      {
        move(metro);
        measure();
        accumulate();  
      }
      averages(iblk);
    }
    average_final(averages_final);
  }
}

double Ising::boltzmann(int i)
{
  double ene = 2 * s[i] * (J * (s[pbc(i - 1)] + s[pbc(i + 1)]) + h);
  return std::exp(-beta*ene);
}

void Ising::measure()
{
  double u = 0.0;
  double m = 0.0;

  //cycle over spins
  for (int i=0; i<nspin; ++i)
  {
    u += -J * s[i] * s[pbc(i+1)] - 0.5 * h * (s[i] + s[pbc(i+1)]);
    m += s[i];
  }
  walker[iu] = u;
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;
}

void Ising::reset(int iblk) //Reset block averages
{ 
  if(iblk == 1)
  {
    for(int i=0; i<n_props; ++i)
    {
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }

  for(int i=0; i<n_props; ++i)
  {
    blk_av[i] = 0;
  }
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}

void Ising::accumulate() //Update block averages
{  
  for(int i=0; i<n_props; ++i)
  {
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
}

void Ising::averages(int iblk) //results for current block
{
  stima_u = blk_av[iu]/blk_norm/(double)nspin; //energy
  glob_av[iu]  += stima_u;
  glob_av2[iu] += stima_u*stima_u;
  err_u = error(glob_av[iu], glob_av2[iu], iblk);

  stima_c = beta*beta*(blk_av[ic]/blk_norm/(double)nspin - (double)nspin*stima_u*stima_u); //heat capacity
  glob_av[ic]  += stima_c;
  glob_av2[ic] += stima_c*stima_c;
  err_c = error(glob_av[ic], glob_av2[ic], iblk);

  stima_m = blk_av[im]/blk_norm/(double)nspin; //magnetization
  glob_av[im]  += stima_m;
  glob_av2[im] += stima_m*stima_m;
  err_m = error(glob_av[im], glob_av2[im], iblk);

  stima_x = beta*blk_av[ix]/blk_norm/(double)nspin; //magnetic susceptibility
  glob_av[ix]  += stima_x;
  glob_av2[ix] += stima_x*stima_x;
  err_x = error(glob_av[ix], glob_av2[ix], iblk);
}

void Ising::average_final(std::string averages_final) //print final averages
{
  std::ofstream ene, heat, magn, chi;
  constexpr int wd = 12;

  std::cout << "Acceptance rate " << accepted/attempted << std::endl << std::endl;

  if (metro == 1) //metropolis
  {
    if (h == 0)
    {
      ene.open(averages_final + "/energy_T_metropolis.out",std::ios::app);
      ene << std::setw(wd) << '\t' << temp << std::setw(wd) << '\t' << glob_av[iu]/(double)nblk << std::setw(wd)
          << '\t' << err_u << std::setw(wd) << std::endl;
      ene.close();

      heat.open(averages_final + "/heat_T_metropolis.out", std::ios::app);
      heat << std::setw(wd) << '\t' << temp << std::setw(wd) << '\t' << glob_av[ic]/(double)nblk << std::setw(wd)
          << '\t' << err_c << std::setw(wd) << std::endl;
      heat.close();

      chi.open(averages_final + "/chi_T_metropolis.out", std::ios::app);
      chi << std::setw(wd) << '\t' << temp << std::setw(wd) << '\t' << glob_av[ix]/(double)nblk << std::setw(wd)
          << '\t' << err_x << std::setw(wd) << std::endl;
      chi.close();
    }
    else if (h != 0)
    {  
      magn.open(averages_final + "/magn_T_metropolis.out", std::ios::app);
      magn << std::setw(wd) << '\t' << temp << std::setw(wd) << '\t' << glob_av[im]/(double)nblk << std::setw(wd)
          << '\t' << err_m << std::setw(wd) << std::endl;
      magn.close();
    }
  }
  else //gibbs
  {
    if (h == 0)
    { 
      ene.open(averages_final + "/energy_T_gibbs.out", std::ios::app);
      ene << std::setw(wd) << '\t' << temp << std::setw(wd) << '\t' << glob_av[iu]/(double)nblk << std::setw(wd)
          << '\t' << err_u << std::setw(wd) << std::endl;
      ene.close();

      heat.open(averages_final + "/heat_T_gibbs.out", std::ios::app);
      heat << std::setw(wd) << '\t' << temp << std::setw(wd) << '\t' << glob_av[ic]/(double)nblk << std::setw(wd)
          << '\t' << err_c << std::setw(wd) << std::endl;
      heat.close();

      chi.open(averages_final + "/chi_T_gibbs.out", std::ios::app);
      chi << std::setw(wd) << '\t' << temp << std::setw(wd) << '\t' << glob_av[ix]/(double)nblk << std::setw(wd)
          << '\t' << err_x << std::setw(wd) << std::endl;
      chi.close();
    }
    else if (h != 0)
    { 
      magn.open(averages_final + "/magn_T_gibbs.out", std::ios::app);
      magn << std::setw(wd) << '\t' << temp << std::setw(wd) << '\t' << glob_av[im]/(double)nblk << std::setw(wd)
          << '\t' << err_m << std::setw(wd) << std::endl;
      magn.close();
    }
  }  
}

void Ising::confFinal(std::string conf_f)
{
  std::ofstream write_conf;

  std::cout << "Print final configuration to file config.final " << std::endl << std::endl;
  write_conf.open(conf_f);
  for (int i=0; i<nspin; ++i)
  {
    write_conf << s[i] << std::endl;
  }
  write_conf.close();
  rnd.save_seed(SEED "/seed.out");
}

int Ising::pbc(int i)  //Algorithm for periodic boundary conditions
{
  if(i >= nspin) i -= nspin;
  else if(i < 0) i += nspin;
  return i;
}

double Ising::error(double sum, double sum2, int iblk)
{
  if(iblk==1) return 0.0;
  else return std::sqrt((sum2/(double)iblk - std::pow(sum/(double)iblk,2))/(double)(iblk-1));
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
