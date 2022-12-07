/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Exercise 10_1
//TSP with GA with MPI 
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <vector>
#include <array>
#include <string>
#include "random.h"
#include "tsp.h"
#include "mpi.h"

bool check_permutation(std::vector<int> idx, int dim);
bool check_first(std::vector<int> path);
double ave_best_half(std::vector<std::vector<int>> paths, TravelCost &t_cost, int dim);

int main(int argc, char* argv[])
{ 
  
  MPI_Init(&argc, &argv);

  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status status;
  MPI_Request req;

  double tstart = MPI_Wtime();

  constexpr int N_cities = 34;
  constexpr int N_paths = 6000;
  constexpr int N_migrations = 20;
  constexpr int N_generations = 10;
  constexpr double p_selection = 0.1;
  constexpr double p_cross = 0.8;
  constexpr double p_mutation = 0.05;

  if (argc != 2)
  { 
    MPI_Finalize();
    std::cerr << "Use: " << argv[0] << " circle | square" << std::endl; 
    return 1;
  }
   
   
  Random rnd( "primes/primes"+std::to_string(rank), SEED "/seed.in");

  std::array<double, N_cities> coord_x;
  std::array<double, N_cities> coord_y;

  MPI_Barrier(MPI_COMM_WORLD);
  
  if (rank == 0)
  { 
    std::ofstream out("data/N_data_10_1.txt", std::ofstream::out);
    if (!out.is_open())
    { 
      MPI_Finalize();
      std::cerr << "Unable to open N_data_10_1.txt" << std::endl;
      return 1;
    }
    out << N_cities << std::endl;
    out << N_paths << std::endl;
    out << N_generations << std::endl; 
    out << p_selection << std::endl;
    out << p_cross << std::endl;
    out << p_mutation << std::endl;

    out.close();
 
    if (std::string{argv[1]} == "circle")
    {
      std::ofstream out_circle("data/10_1_" + std::string{argv[1]} + "_cities.txt" , std::ofstream::out);
      if (!out_circle.is_open())
      { 
        MPI_Finalize();
        std::cerr << "data/10_1_" + std::string{argv[1]} + "_cities.txt" << std::endl;
        return 1;
      }

      for (int i = 0; i < N_cities; ++i)
      {
        double theta = rnd.rannyu(0., 2 * M_PI);
        coord_x[i] = std::cos(theta);
        coord_y[i] = std::sin(theta);
        out_circle << i << '\t' << std::cos(theta) << '\t' << std::sin(theta) << std::endl;
      } 
      out_circle.close();
    }
    else if (std::string{argv[1]} == "square")
    {
      std::ofstream out_square("data/10_1_" + std::string{argv[1]} + "_cities.txt", std::ofstream::out);
      if (!out_square.is_open())
      { 
        MPI_Finalize();
        std::cerr << "data/10_1_" + std::string{argv[1]} + "_cities.txt" << std::endl;
        return 1;
      }

      for (int j = 0; j < N_cities; ++j)
      {
        double x = rnd.rannyu();
        double y = rnd.rannyu();
        coord_x[j] = x;
        coord_y[j] = y;
        out_square << j << '\t' << x << '\t' << y << std::endl;
      } 
      out_square.close();
    }
    else
    { 
      MPI_Finalize();
      std::cerr << "Use: " << argv[0] << " circle | square" << std::endl; 
      return 1;
    }
  }

  MPI_Bcast(coord_x.data(), N_cities, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(coord_y.data(), N_cities, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);

  std::vector<City> cities;
  for (int k = 0; k < N_cities; ++k)
  {
    cities.push_back(City(coord_x[k], coord_y[k]));
  }
 
  //Find the travel cost for every path 
  TravelCost t_cost(cities);
  //Class mutations: selection, crossover and different types of mutations   
  Genetic gen(rnd, N_paths, N_cities, p_selection, p_cross, p_mutation);  
  //different paths
  std::vector<std::vector<int>> paths;
  paths.reserve(N_paths);

  Ranking rkg(t_cost);
  //fill different paths with random cities except the first
  for (int i = 0; i < N_paths; ++i)
  { 
    std::vector<int> city_indexes;
    city_indexes.reserve(N_cities);
    city_indexes.push_back(static_cast<int>(0));
    for (int j = 1; j < N_cities; ++j)
    {
      for(;;)
      {
        int r_city = static_cast<int>(rnd.rannyu(1, N_cities)); 
        if (std::find(city_indexes.begin(), city_indexes.end(), r_city) == city_indexes.end())
        {
          city_indexes.push_back(r_city);
          break;
        }
      }
    }  
    paths.push_back(city_indexes);
  }
  
  //paths' ranking based on travel cost 
  rkg.ranking(paths); 

  std::ofstream out_paths("data/10_1_" + std::string{argv[1]} + std::to_string(rank) + "_paths_out.txt" , std::ofstream::out);
  if (!out_paths.is_open())
  {  
    MPI_Finalize();
    std::cerr << "data/10_1_" + std::string{argv[1]} + std::to_string(rank) + "_paths_out.txt" << std::endl;
    return 1;
  }

  out_paths << t_cost.travel_cost(paths[paths.size() - 1])  << '\t' 
            << ave_best_half(paths, t_cost, N_paths / 2) << std::endl;

  //running GA for TSP: every N_migrations generations exchange their best individuals randomly
  for (int j = 0; j < N_migrations; ++j) 
  {
    for (int i = 0; i < N_generations; ++i) 
    { 
      std::vector<std::vector<int>> new_paths;
      new_paths.reserve(N_paths);
     
      for (size_t k = 0; k < N_paths; k += 2)
      { 
        std::vector<int> p1_path = paths[gen.selection()];
        std::vector<int> p2_path = paths[gen.selection()];

        std::vector<int> s1_path;
        std::vector<int> s2_path;
        gen.crossover(p1_path, p2_path, s1_path, s2_path);

        gen.mutation(s1_path);
        gen.mutation(s2_path);
        new_paths.push_back(s1_path);
        new_paths.push_back(s2_path);
      } 
      paths = new_paths;
      rkg.ranking(paths);
      out_paths << t_cost.travel_cost(paths[paths.size() - 1])  << '\t' 
                << ave_best_half(paths, t_cost, N_paths / 2) << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    std::vector<int> send_cont(size);

    if (rank == 0)
    {
      send_cont[0] = 0;
      for (int i = 1; i < size; ++i)
      {
        for (;;)
        {
          int cont = static_cast<int>(rnd.rannyu(1.0, static_cast<double>(size)));
          if (std::find(send_cont.begin(), send_cont.end(), cont) == send_cont.end()) 
          {
            send_cont[i] = cont;
            break;
          }
        }
      }
    } 
    MPI_Bcast(send_cont.data(), send_cont.size(), MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> send_path(paths[paths.size() - 1]);
    std::vector<int> recv_path(send_path.size());

    if (size % 2 == 0) //even number of cores
    { 
      for (int i = 0; i < size; i += 2)
      {
        if (rank == send_cont[i])
        {
          MPI_Sendrecv(send_path.data(), send_path.size(), MPI_INT, send_cont[i+1], 0, 
                       recv_path.data(), send_path.size(), MPI_INT, send_cont[i+1], MPI_ANY_TAG, MPI_COMM_WORLD, &status); 

          std::cout << "rank send: " << send_cont[i+1] << " rank receive: " << send_cont[i] << std::endl;
        }
        else if (rank == send_cont[i+1])
        {
          MPI_Sendrecv(send_path.data(), send_path.size(), MPI_INT, send_cont[i], 0, 
                       recv_path.data(), send_path.size(), MPI_INT, send_cont[i], MPI_ANY_TAG, MPI_COMM_WORLD, &status); 

          std::cout << "rank send: " << send_cont[i] << " rank receive: " << send_cont[i+1] << std::endl;
        }
      }
    }
    else //odd number of cores
    {
      int add_cont;
      if (rank == 0)
      {
        for (;;)  
        {
          add_cont = static_cast<int>(rnd.rannyu(1.0, static_cast<double>(size)));
          if (add_cont != send_cont[size-1]) break;
        } 
      }
      MPI_Bcast(&add_cont, 1, MPI_INT, 0, MPI_COMM_WORLD);

      for (int i = 0; i < size; i += 2)
      {
        if (send_cont[i] != send_cont[size-1])
        {
          if (rank == send_cont[i])
          {
            MPI_Sendrecv(send_path.data(), send_path.size(), MPI_INT, send_cont[i+1], 0, 
                         recv_path.data(), send_path.size(), MPI_INT, send_cont[i+1], MPI_ANY_TAG, MPI_COMM_WORLD, &status); 

            std::cout << "rank send: " << send_cont[i+1] << " rank receive: " << send_cont[i] << std::endl;
          }
          else if (rank == send_cont[i+1])
          {
            MPI_Sendrecv(send_path.data(), send_path.size(), MPI_INT, send_cont[i], 0, 
                         recv_path.data(), send_path.size(), MPI_INT, send_cont[i], MPI_ANY_TAG, MPI_COMM_WORLD, &status); 

            std::cout << "rank send: " << send_cont[i] << " rank receive: " << send_cont[i+1] << std::endl;
          }
        }
        else
        {    
          if (rank == send_cont[i])
          {
            MPI_Isend(send_path.data(), send_path.size(), MPI_INT, add_cont, 1, MPI_COMM_WORLD, &req); 
            MPI_Recv(recv_path.data(), send_path.size(), MPI_INT, add_cont, 2, MPI_COMM_WORLD, NULL);
            std::cout << "rank send: " << rank << " rank receive: " << add_cont << std::endl;
          }
          else if (rank == add_cont)
          {
            MPI_Send(send_path.data(), send_path.size(), MPI_INT, send_cont[size-1], 2, MPI_COMM_WORLD); 
            MPI_Recv(recv_path.data(), send_path.size(), MPI_INT, send_cont[size-1], 1, MPI_COMM_WORLD, NULL);
            std::cout << "rank send: " << add_cont << " rank receive: " << send_cont[size-1] << std::endl;
          }
        }
      }
    }
   
    MPI_Barrier(MPI_COMM_WORLD);      
    paths[paths.size() - 1] = recv_path; 
    rkg.ranking(paths);
  }
  out_paths.close();

  std::ofstream best_path("data/10_1_" + std::string{argv[1]} + std::to_string(rank) + "_best_path.txt" , std::ofstream::out);
  if (!best_path.is_open())
  {
    MPI_Finalize();
    std::cerr << "data/10_1_" + std::string{argv[1]} + std::to_string(rank) + "_best_path.txt" << std::endl;
    return 1;
  }
  
  for (const auto& el : paths[paths.size()-1])
  {
    best_path << el << std::endl;
  } 

  best_path.close();

  double tend = MPI_Wtime();
  double dt = tend - tstart;

  std::cout << "Rank number: " << rank << "; Time = " << dt << std::endl; 
  MPI_Finalize();
  return 0;
}

bool check_permutation(std::vector<int> idx, int dim)
{
  std::vector<int> test_idx(dim);
  std::iota(test_idx.begin(), test_idx.end(), 0);
  return std::is_permutation(test_idx.begin(), test_idx.end(), idx.begin());
}

bool check_first(std::vector<int> path)
{
  bool flag = false;
  if (path[0] == 0)
  {
    flag = true;
  }
  return flag;
}

double ave_best_half(std::vector<std::vector<int>> paths, TravelCost &t_cost, int dim)
{
  double acc = 0.0;
  return std::accumulate(paths.begin() + dim, paths.end(), acc, 
    [&t_cost](const double i, const std::vector<int> &path){return i + t_cost.travel_cost(path) ;}) / dim;  
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
