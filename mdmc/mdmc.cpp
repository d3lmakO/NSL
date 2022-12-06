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
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include "mdmc.h"

//ReadInput >> iNVET; 0=MD(NVE) 1=MC(NVT)
//ReadInput >> temp;
//ReadInput >> npart;
//ReadInput >> rho;
//ReadInput >> rcut;
//ReadInput >> delta;
//ReadInput >> nblk;
//ReadInput >> nstep;
//ReadInput >> nbins;

Mdmc::Mdmc(int restart) : rnd(SEED "/Primes", SEED "/seed.in"), restart{restart}
{
  if (restart) 
  {
    int p1, p2; 
    int seed[4]; 
    std::ifstream primes(SEED "/Primes");
    if (!primes.is_open())
    {
      std::cerr << "Unable to open Primes" << std::endl;
      exit(1);
    }
    primes >> p1 >> p2;
    primes.close();
    
    std::ifstream seed_out(SEED "/seed.out"); 
    if (!seed_out.is_open())
    {
      std::cerr << "Unable to open seed.out" << std::endl;
      exit(1);
    }
    seed_out >> seed[0] >> seed[1] >> seed[2] >> seed[3]; 
    seed_out.close();

    rnd.set_random(seed, p1, p2);
  }
}

void Mdmc::input(std::string input_file, std::string conf_file, std::string rest_dir)
{
  std::ifstream readInput, readConf, readVelocity;

  std::cout << "Classic Lennard-Jones fluid        " << std::endl;
  std::cout << "MD(NVE)-->0 / MC(NVT)-->1 simulation       " << std::endl << std::endl;
  std::cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << std::endl << std::endl;
  std::cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << std::endl << std::endl;
  std::cout << "The program uses Lennard-Jones units " << std::endl;

  //Read input informations
  readInput.open(input_file);

  readInput >> iNVET;

  readInput >> temp;
  beta = 1.0/temp;
  std::cout << "Temperature = " << temp << std::endl;

  readInput >> npart;
  std::cout << "Number of particles = " << npart << std::endl;

  readInput >> rho;
  std::cout << "Density of particles = " << rho << std::endl;
  vol = (double)npart/rho;
  box = std::pow(vol,1.0/3.0);
  std::cout << "Volume of the simulation box = " << vol << std::endl;
  std::cout << "Edge of the simulation box = " << box << std::endl;

  readInput >> rcut;
  std::cout << "Cutoff of the interatomic potential = " << rcut << std::endl << std::endl;
    
  readInput >> delta;

  readInput >> nblk;

  readInput >> nstep;
  readInput >> nequil;
  readInput >> nbins;

  bin_size = (box/2.0)/(double)nbins;

  //tail corrections for pot energy and pressure 
  vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));

  std::cout << "The program perform Metropolis moves with uniform translations" << std::endl;
  std::cout << "Moves parameter = " << delta << std::endl;
  std::cout << "Number of blocks = " << nblk << std::endl;
  std::cout << "Number of steps in one block = " << nstep << std::endl;
  std::cout << "Number of steps for Temperature equilibration = " << nequil << std::endl;
  std::cout << "Tail correction for potential energy = " << vtail << std::endl;
  std::cout << "Tail correction for pressure = " << ptail << std::endl << std::endl;

  readInput.close();

  //Prepare arrays for measurements
  iv = 0; //Potential energy
  it = 1; //Temperature
  ik = 2; //Kinetic energy
  ie = 3; //Total energy
  iw = 4; //Pressure
  n_props = 5; //Number of observables

  ig = 5; //Radial distribution function

  n_props = n_props + nbins;

  //Read initial configuration
  std::cout << "Read initial configuration" << std::endl << std::endl;
  if(restart)
  {
    readConf.open(rest_dir + "config.out");
    readVelocity.open(rest_dir + "velocity.out");
    if (!readConf.is_open())
    {
      std::cerr << "Unable to open config.out" << std::endl;
      exit(1);
    }

    if (!readVelocity.is_open())
    { 
      std::cerr << "Unable to open velocity.out" << std::endl;
      exit(1);
    }

    for (int i=0; i<npart; ++i) readVelocity >> vx[i] >> vy[i] >> vz[i];
  }
  else 
  {
    readConf.open(conf_file);
    std::cout << "Prepare velocities with center of mass velocity equal to zero " << std::endl;
    double sumv[3] = {0.0, 0.0, 0.0};
    for (int i=0; i<npart; ++i)
    {
      vx[i] = rnd.gauss(0.,std::sqrt(temp));
      vy[i] = rnd.gauss(0.,std::sqrt(temp));
      vz[i] = rnd.gauss(0.,std::sqrt(temp));
      sumv[0] += vx[i];
      sumv[1] += vy[i];
      sumv[2] += vz[i];
    }
    for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i)
    {
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];
      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;
    fs = std::sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
    std::cout << "velocity scale factor: " << fs << std::endl << std::endl;
    for (int i=0; i<npart; ++i)
    {
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;
    }
  }

  for (int i=0; i<npart; ++i)
  {
    readConf >> x[i] >> y[i] >> z[i];
    x[i] = pbc( x[i] * box );
    y[i] = pbc( y[i] * box );
    z[i] = pbc( z[i] * box );
  }
  readConf.close();

  for (int i=0; i<npart; ++i)
  {
    if(iNVET)
    {
      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];
    }
    else
    {
      xold[i] = pbc(x[i] - vx[i] * delta);
      yold[i] = pbc(y[i] - vy[i] * delta);
      zold[i] = pbc(z[i] - vz[i] * delta);
    }
  }
  
  //Evaluate properties of the initial configuration
  measure();

  //Print initial values for measured properties
  std::cout << "Initial potenzial energy with correction = " << walker[iv]/(double)npart + vtail << std::endl;
  std::cout << "Pressure with correction = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << std::endl;
  //std::cout << "Initial potential energy = " << walker[iv]/(double)npart << std::endl;
  //std::cout << "Initial temperature      = " << walker[it] << std::endl;
  //std::cout << "Initial kinetic energy   = " << walker[ik]/(double)npart << std::endl;
  //std::cout << "Initial total energy     = " << walker[ie]/(double)npart << std::endl;
  //std::cout << "Initial pressure         = " << walker[iw]/(double)npart << std::endl;
  return;
}

void Mdmc::move()
{
  int o;
  double p, energy_old, energy_new;
  double xnew, ynew, znew;

  if(iNVET) // Monte Carlo (NVT) move
  {
    for(int i=0; i<npart; ++i)
    {
      //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
      o = (int)(rnd.rannyu()*npart);

      //Old
      energy_old = boltzmann(x[o],y[o],z[o],o);

      //New
      x[o] = pbc( x[o] + delta*(rnd.rannyu() - 0.5) );
      y[o] = pbc( y[o] + delta*(rnd.rannyu() - 0.5) );
      z[o] = pbc( z[o] + delta*(rnd.rannyu() - 0.5) );

      energy_new = boltzmann(x[o],y[o],z[o],o);

      //Metropolis test
      p = std::exp(beta*(energy_old-energy_new));
      if(p >= rnd.rannyu())  
      {
      //Update
        xold[o] = x[o];
        yold[o] = y[o];
        zold[o] = z[o];
        accepted = accepted + 1.0;
      } else {
        x[o] = xold[o];
        y[o] = yold[o];
        z[o] = zold[o];
      }
      attempted = attempted + 1.0;
    }
  } else // Molecular Dynamics (NVE) move
  {
    double fx[m_part], fy[m_part], fz[m_part];

    for(int i=0; i<npart; ++i)
    { //Force acting on particle i
      fx[i] = force(i,0);
      fy[i] = force(i,1);
      fz[i] = force(i,2);
    }

    for(int i=0; i<npart; ++i)
    { //Verlet integration scheme
      xnew = pbc( 2.0 * x[i] - xold[i] + fx[i] * std::pow(delta,2) );
      ynew = pbc( 2.0 * y[i] - yold[i] + fy[i] * std::pow(delta,2) );
      znew = pbc( 2.0 * z[i] - zold[i] + fz[i] * std::pow(delta,2) );

      vx[i] = pbc(xnew - xold[i])/(2.0 * delta);
      vy[i] = pbc(ynew - yold[i])/(2.0 * delta);
      vz[i] = pbc(znew - zold[i])/(2.0 * delta);

      xold[i] = x[i];
      yold[i] = y[i];
      zold[i] = z[i];

      x[i] = xnew;
      y[i] = ynew;
      z[i] = znew;

      accepted += 1.0;
      attempted += 1.0;
    }
  }
  return;
}

double Mdmc::boltzmann(double xx, double yy, double zz, int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
      // distance ip-i in pbc
      dx = pbc(xx - x[i]);
      dy = pbc(yy - y[i]);
      dz = pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = std::sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/std::pow(dr,12) - 1.0/std::pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}

double Mdmc::force(int ip, int idir)
{ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
      dvec[0] = pbc( x[ip] - x[i] );  // distance iw-i in pbc
      dvec[1] = pbc( y[ip] - y[i] );
      dvec[2] = pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = std::sqrt(dr);

      if(dr < rcut)
      {
        f += dvec[idir] * (48.0/std::pow(dr,14) - 24.0/std::pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
 
  return f;
}

void Mdmc::measure() //Properties measurement
{
  double v = 0.0, p = 0.0;
  //double kin = 0.0;
  double vij, pij;
  double dx, dy, dz, dr;
  double r;

  for (int k = ig; k < ig + nbins; ++k)
  {
    walker[k] = 0.0;
  }

  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {
      //distance i-j in pbc
      dx = pbc(x[i] - x[j]);
      dy = pbc(y[i] - y[j]);
      dz = pbc(z[i] - z[j]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = std::sqrt(dr);

      for (int k = ig; k <= ig + nbins; ++k) 
      {
        r = bin_size * (k - ig);
        if (r < dr && dr < r + bin_size)
        {
          walker[k] += 2;
        }
      }

      if(dr < rcut)
      {
        vij = 1.0/std::pow(dr,12) - 1.0/std::pow(dr,6);
        v += vij;
        //Pressure
        pij = 1.0/std::pow(dr,12) - 0.5/std::pow(dr,6);
        p += pij;
      }
    }          
  }

  //for (int i=0; i<npart; ++i) kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

  walker[iv] = 4.0 * v; // Potential energy
  //walker[ik] = kin; // Kinetic energy
  //walker[it] = (2.0 / 3.0) * kin/(double)npart; // Temperature
  //walker[ie] = 4.0 * v + kin;  // Total energy;
  walker[iw] = 48.0 * p / 3.0; //Pressure

  return;
}

void Mdmc::equilibration(std::string equil_file, std::string conf_final, std::string vel_final)
{ 
  const int wd=12;
  std::ofstream out_equil(equil_file, std::ios::app);
  if (!out_equil.is_open())
  {
    std::cerr << "Unable to open equilibration file." << std::endl;
    exit(1);
  }
  for (int i = 0; i < nequil; ++i)
  {
    move();
    measure();
    //out_equil << std::setw(wd) << walker[ik]/(double)npart << '\t' << std::setw(wd) << walker[iv]/(double)npart << '\t'
    //          << std::setw(wd) << walker[ie]/(double)npart << '\t' << std::setw(wd) << walker[it] << std::endl; 
    out_equil << std::setw(wd) << walker[iv]/(double)npart + vtail << std::endl;
              //<< std::setw(wd) << rho * temp + (walker[iw] + (double)npart * ptail) / vol << std::endl;
    std::cout << accepted / attempted << std::endl;
  }

  out_equil.close();

  confFinal(conf_final, vel_final);
}

void Mdmc::out_istant_values(std::string path)
{ 
  const int wd=12;
  std::ofstream pth(path, std::ios::app);
  if (!pth.is_open())
  {
    std::cerr << "Unable to open istantaneous values file." << std::endl;
    exit(1);
  }

  pth << std::setw(wd) << walker[iv]/(double)npart + vtail << std::endl;
  pth.close();
}

void Mdmc::reset(int iblk) //Reset block averages
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

void Mdmc::accumulate(void) //Update block averages
{
  for(int i=0; i<n_props; ++i)
  {
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
}

void Mdmc::averages(int iblk, std::string phase) //Print results for current block
{  
  double r;
  double stima_gdir[ig + nbins];
  double err_gdir[ig + nbins];
  std::ofstream Epot, Ekin, Etot, Temp, Press, Gdir, Gfin;
  const int wd=12;
    
  std::cout << "Block number " << iblk << std::endl;
  std::cout << "Acceptance rate " << accepted/attempted << std::endl << std::endl;
   
  //Epot.open("data/"+phase+"_output_epot.dat",std::ios::app);
  //Ekin.open("data/"+phase+"_output_ekin.dat",std::ios::app);
  //Temp.open("data/"+phase+"_output_temp.dat",std::ios::app);
  //Etot.open("data/"+phase+"_output_etot.dat",std::ios::app);
  //Press.open("data/"+phase+"_output_press.dat",std::ios::app);
  //Gdir.open("data/"+phase+"_output_gdir.dat", std::ios::app);
  Gfin.open("data/"+phase+"_output_gdir_final.dat", std::ios::app);
   
  //stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
  //glob_av[iv] += stima_pot;
  //glob_av2[iv] += stima_pot*stima_pot;
  //err_pot = error(glob_av[iv],glob_av2[iv],iblk);
   
  //stima_kin = blk_av[ik]/blk_norm/(double)npart; //Kinetic energy
  //glob_av[ik] += stima_kin;
  //glob_av2[ik] += stima_kin*stima_kin;
  //err_kin = error(glob_av[ik],glob_av2[ik],iblk);

  //stima_etot = blk_av[ie]/blk_norm/(double)npart; //Total energy
  //glob_av[ie] += stima_etot;
  //glob_av2[ie] += stima_etot*stima_etot;
  //err_etot = error(glob_av[ie],glob_av2[ie],iblk);

  //stima_temp = blk_av[it]/blk_norm; //Temperature
  //glob_av[it] += stima_temp;
  //glob_av2[it] += stima_temp*stima_temp;
  //err_temp = error(glob_av[it],glob_av2[it],iblk);

  //stima_pres = rho * stima_temp + (blk_av[iw]/blk_norm) / vol; //Pressure
  //glob_av[iw] += stima_pres;
  //glob_av2[iw] += stima_pres*stima_pres;
  //err_press = error(glob_av[iw],glob_av2[iw],iblk);

  //stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
  //glob_av[iw] += stima_pres;
  //glob_av2[iw] += stima_pres*stima_pres;
  //err_press = error(glob_av[iw],glob_av2[iw],iblk);

  for (int i = ig; i < ig + nbins; ++i) 
  {
    r = bin_size * (i - ig);
    stima_gdir[i] = blk_av[i]/(blk_norm*(rho * npart * (4.0 * pi / 3.0)*
                              ((r + bin_size)*(r + bin_size)*(r + bin_size) - r*r*r)));
    glob_av[i] += stima_gdir[i];
    glob_av2[i] += stima_gdir[i] * stima_gdir[i];
    err_gdir[i] = error(glob_av[i], glob_av2[i], iblk);
  }

  //Potential energy per particle
  //Epot << std::setw(wd) << iblk << '\t' << std::setw(wd) << stima_pot << '\t' << std::setw(wd) 
  //     << glob_av[iv]/(double)iblk << '\t' << std::setw(wd) << err_pot << std::endl;
  //Kinetic energy
  //Ekin << std::setw(wd) << iblk << '\t' << std::setw(wd) << stima_kin << '\t' << std::setw(wd) 
  //     << glob_av[ik]/(double)iblk << '\t' << std::setw(wd) << err_kin << std::endl;
  //Total energy
  //Etot << std::setw(wd) << iblk << '\t' << std::setw(wd) << stima_etot << '\t' << std::setw(wd) 
  //     << glob_av[ie]/(double)iblk << '\t' << std::setw(wd) << err_etot << std::endl;
  //Temperature
  //Temp << std::setw(wd) << iblk << '\t' << std::setw(wd) << stima_temp << '\t' << std::setw(wd) 
  //     << glob_av[it]/(double)iblk << '\t' << std::setw(wd) << err_temp << std::endl;
  //Pressure
  //Press << std::setw(wd) << iblk << '\t' << std::setw(wd) << stima_pres << '\t' << std::setw(wd)
  //      << glob_av[iw]/(double)iblk << '\t' << std::setw(wd) << err_press << std::endl;

  //Radial distribution function
  //for (int j = ig; j < ig + nbins; ++j)
  //{
  //  Gdir << std::setw(wd) << iblk << '\t' << std::setw(wd) << j << '\t' << std::setw(wd) 
  //       << stima_gdir[j] << std::endl;
  //}

  //Final radial distribution function with statistical uncertainties
  if (iblk == nblk)
  {
    for (int k = ig; k < ig + nbins; ++k)
    { 
      r = bin_size * (k - ig);  
      Gfin << std::setw(wd) << r << '\t' << std::setw(wd) << glob_av[k]/(double)iblk << '\t' 
                            << std::setw(wd) << err_gdir[k] << std::endl;
    }
  }

  std::cout << "----------------------------" << std::endl << std::endl;

  //Epot.close();
  //Ekin.close();
  //Etot.close();
  //Temp.close();
  //Press.close();
  //Gdir.close();
  Gfin.close();
}

void Mdmc::confFinal(std::string conf_final, std::string vel_final)
{
  std::ofstream writeConf, writeVelocity;

  std::cout << "Print final configuration to file config.out" << std::endl << std::endl;
  writeConf.open(conf_final);
  writeVelocity.open(vel_final);
  for (int i=0; i<npart; ++i)
  {
    writeConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << std::endl;
    writeVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << std::endl;
  }
  writeConf.close();
  writeVelocity.close();

  rnd.save_seed(SEED "/seed.out");
}

void Mdmc::confXYZ(int nconf)
{ //Write configuration in .xyz format
  std::ofstream writeXYZ;

  writeXYZ.open("frames/config_" + std::to_string(nconf) + ".xyz");
  writeXYZ << npart << std::endl;
  writeXYZ << "This is only a comment!" << std::endl;
  for (int i=0; i<npart; ++i)
  {
    writeXYZ << "LJ  " << pbc(x[i]) << "   " <<  pbc(y[i]) << "   " << pbc(z[i]) << std::endl;
  }
  writeXYZ.close();
}

double Mdmc::pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
  return r - box * std::rint(r/box);
}

double Mdmc::error(double sum, double sum2, int iblk)
{
  return std::sqrt(std::fabs(sum2/(double)iblk - std::pow(sum/(double)iblk,2))/(double)iblk);
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
