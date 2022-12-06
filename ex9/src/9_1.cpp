/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Exercise 9_1
//traveling salesman problem (TSP) with a genetic algorithm (GA)
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <numeric>
#include <tuple>
#include <vector>
#include <string>
#include "random.h"
#include "tsp.h"

bool check_permutation(std::vector<int> idx, int dim);
bool check_first(std::vector<int> path);
double ave_best_half(std::vector<std::vector<int>> paths, TravelCost &t_cost, int dim);

int main(int argc, char* argv[])
{ 

  constexpr int N_cities = 34;
  constexpr int N_paths = 6000;
  constexpr int N_generations = 200;
  constexpr double p_selection = 0.1;
  constexpr double p_cross = 0.8;
  constexpr double p_mutation = 0.05;

  if (argc != 2)
  {
    std::cerr << "Use: " << argv[0] << " circle | square" << std::endl; 
    return 1;
  }
   
  Random rnd(SEED "/Primes", SEED "/seed.in");

  std::ofstream out("data/N_data_9_1.txt", std::ofstream::out);
  if (!out.is_open())
  {
    std::cerr << "Unable to open N_data_9_1.txt" << std::endl;
    return 1;
  }
  out << N_cities << std::endl;
  out << N_paths << std::endl;
  out << N_generations << std::endl; 
  out << p_selection << std::endl;
  out << p_cross << std::endl;
  out << p_mutation << std::endl;

  out.close();

  std::vector<City> cities;
  cities.reserve(N_cities);
  
  if (std::string{argv[1]} == "circle")
  {
    std::ofstream out_circle("data/9_1_" + std::string{argv[1]} + "_cities.txt" , std::ofstream::out);
    if (!out_circle.is_open())
    {
      std::cerr << "data/9_1_" + std::string{argv[1]} + "_cities.txt" << std::endl;
      return 1;
    }

    for (int i = 0; i < N_cities; ++i)
    {
      double theta = rnd.rannyu(0., 2 * M_PI);
      cities.push_back(City{std::cos(theta), std::sin(theta)});
      out_circle << i << '\t' << std::cos(theta) << '\t' << std::sin(theta) << std::endl;
    } 
    out_circle.close();
  }
  else if (std::string{argv[1]} == "square")
  {
    std::ofstream out_square("data/9_1_" + std::string{argv[1]} + "_cities.txt", std::ofstream::out);
    if (!out_square.is_open())
    {
      std::cerr << "data/9_1_" + std::string{argv[1]} + "_cities.txt" << std::endl;
      return 1;
    }

    for (int j = 0; j < N_cities; ++j)
    {
      double x = rnd.rannyu();
      double y = rnd.rannyu();
      cities.push_back(City{x, y});
      out_square << j << '\t' << x << '\t' << y << std::endl;
    } 
    out_square.close();
  }
  else
  {
    std::cerr << "Use: " << argv[0] << " circle | square" << std::endl; 
    return 1;
  }

  //Find the travel cost for every path 
  TravelCost t_cost(cities);
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

  Genetic gen(rnd, N_paths, N_cities, p_selection, p_cross, p_mutation);  

  std::ofstream out_paths("data/9_1_" + std::string{argv[1]} + "_paths_out.txt" , std::ofstream::out);
  if (!out_paths.is_open())
  {
    std::cerr << "data/9_1_" + std::string{argv[1]} + "_paths_out.txt" << std::endl;
    return 1;
  }

  out_paths << t_cost.travel_cost(paths[paths.size() - 1])  << '\t' 
            << ave_best_half(paths, t_cost, N_paths / 2) << std::endl;
  
  //running GA for TSP 
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

  out_paths.close();

  std::ofstream best_path("data/9_1_" + std::string{argv[1]} + "_best_path.txt" , std::ofstream::out);
  if (!best_path.is_open())
  {
    std::cerr << "data/9_1_" + std::string{argv[1]} + "_best_path.txt" << std::endl;
    return 1;
  }
  
  for (const auto& el : paths[paths.size()-1])
  {
    best_path << el << std::endl;
  } 


  best_path.close();
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
