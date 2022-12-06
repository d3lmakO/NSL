#ifndef TSP_H
#define TSP_H

#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <algorithm>
#include "random.h"

class City
{
public:  
  City() { x = 0.0; y = 0.0; }
  City(double x, double y) : x{x}, y{y} {}
  double distance(City c)
  {
    return ((c.x - x)*(c.x - x) + (c.y - y)*(c.y - y));
  }
private:
  double x;
  double y;
};

class TravelCost
{
public:
  TravelCost(std::vector<City> cities) : cities{cities} {}
  double travel_cost(std::vector<int> city_indexes) 
  {
    double dist = 0.0;
    for (auto& idx : city_indexes)
    {
      if (idx == city_indexes.back())
      {
        dist += cities[idx].distance(cities[city_indexes[0]]);
      }
      else
      {
        dist += cities[*(&idx + 1)].distance(cities[idx]);
      }
    } 
    return dist;
  }
private:
  std::vector<City> cities;
};

class Ranking
{
public:
  Ranking(TravelCost &t_cost) :
          t_cost{t_cost} {}
  void ranking(std::vector<std::vector<int>> &paths)
  { 
    //sorting paths based on their travel cost
    std::sort(paths.begin(), paths.end(), [this]
             (const std::vector<int> &path1, const std::vector<int> &path2)
             {return t_cost.travel_cost(path1) > t_cost.travel_cost(path2);});
  }
private:
  TravelCost &t_cost;
};

class Genetic
{
public:
  Genetic(Random &rnd, int N_paths, int N_cities, 
          double p_selection, double p_cross, double p_mutation) : 
          rnd{rnd}, N_paths{N_paths}, N_cities{N_cities}, 
          p_selection{p_selection}, p_cross{p_cross}, p_mutation{p_mutation} {}
  int selection()
  {
    //selects path : selection probability set to 0.1 so we can select paths with lower travel_cost
    //WARNING! : +1 in formula lead to segmentation fault in vector -->  int(N * r^p) + 1 
    auto sel = [this]()
               {return static_cast<int>(N_paths * std::pow(rnd.rannyu(), p_selection));};
    int sel_path = sel();
    return sel_path;
  } 

  //crossover probability set to 0.8
  void crossover(std::vector<int> p1_path, std::vector<int> p2_path,
                 std::vector<int> &s1_path, std::vector<int> &s2_path)
  { 
    s1_path = p1_path;
    s2_path = p2_path;
    if (rnd.rannyu() < p_cross)
    {
      size_t point_cut = static_cast<size_t>(rnd.rannyu(1, N_cities - 1));
      size_t idx1 = point_cut;
      size_t idx2 = point_cut;
      for (size_t i = 0; i < p2_path.size(); ++i)
      {
        int city1 = p2_path[i]; 
        if (std::find(p1_path.begin(), p1_path.begin() + point_cut, city1) == 
            p1_path.begin() + point_cut)
        {
          s1_path[idx1++] = city1;
        }
      }
      for (size_t j = 0; j < p1_path.size(); ++j)
      {
        int city2 = p1_path[j];
        if (std::find(p2_path.begin(), p2_path.begin() + point_cut, city2) == 
            p2_path.begin() + point_cut)
        {
          s2_path[idx2++] = city2;
        }
      }
    }
  }
 
  //mutation probability set to 0.05
  void mutation(std::vector<int> &path)
  {
    //pair permutation of cities
    if (rnd.rannyu() < p_mutation)
    {
      std::swap(path[static_cast<int>(rnd.rannyu(1, N_cities))], 
                path[static_cast<int>(rnd.rannyu(1, N_cities))]);
    } 
    
    //shift of n positions for m contiguous cities
    if (rnd.rannyu() < p_mutation)
    {
      int n1 = static_cast<int>(rnd.rannyu(1, N_cities));
      int n2 = static_cast<int>(rnd.rannyu(1, N_cities));
      int n3 = static_cast<int>(rnd.rannyu(1, N_cities));
      std::vector<int> srt = {n1, n2, n3};
      std::sort(srt.begin(), srt.end());
      std::rotate(path.begin() + srt[0], 
                  path.begin() + srt[1], 
                  path.begin() + srt[2]);
    }
    

    //permutation among m contiguous cities
    if (rnd.rannyu() < p_mutation)
    {
      int p_begin = static_cast<int>(rnd.rannyu(1, N_cities / 2)); 
      int p_end = static_cast<int>(rnd.rannyu(1, N_cities / 2));
      if (p_begin <= p_end)
      {
        std::rotate(path.begin() + (p_end - p_begin), 
                    path.begin() + p_end, 
                    path.begin() + (p_end + p_begin));
      }
      else
      {
        std::rotate(path.begin() + (p_begin - p_end), 
                    path.begin() + p_begin, 
                    path.begin() + (p_begin + p_end));
      }
    }

    //inversion of m contiguous cities
    if (rnd.rannyu() < p_mutation)
    {
      int i = static_cast<int>(rnd.rannyu(1, N_cities - 1));
      int j = static_cast<int>(rnd.rannyu(i + 1, N_cities));
      std::reverse(path.begin() + i, path.begin() + j);
    }
  }
private:
  Random &rnd;
  int N_paths;
  int N_cities;
  double p_selection;
  double p_cross;
  double p_mutation;
};

#endif //TSP_H
