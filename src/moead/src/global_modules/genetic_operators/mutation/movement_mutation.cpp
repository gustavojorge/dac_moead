#include "../../../../headers/global_modules/generate_initial_population/population.h"
#include "../../../../headers/global_modules/genetic_operators/movement_mutation.h"
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

using namespace std;

static double nearest_neighbor_distance(vector<Turbine> &zone_turbines, int selected){
  double nearest = numeric_limits<double>::max();
  for(int i = 0; i < (int) zone_turbines.size(); i++){
    if(i != selected){
      nearest = min(nearest, hypot(zone_turbines[i].x - zone_turbines[selected].x, zone_turbines[i].y - zone_turbines[selected].y));
    }
  }
  return nearest;
}

void movement_mutation(Solution &solution){

  int zone = rand() % solution.turbines.size();
  vector<Turbine> &zone_turbines = solution.turbines[zone];
  int zone_size = solution.grid[zone].size();

  if(zone_turbines.empty() || (int) zone_turbines.size() >= zone_size){
    return;
  }

  int selected = rand() % zone_turbines.size();
  Turbine &turbine = zone_turbines[selected];
  double radius = nearest_neighbor_distance(zone_turbines, selected);

  vector<int> candidates;
  for(int p = 0; p < zone_size; p++){
    if(solution.grid[zone][p] == 0 && hypot(foundations[zone][p].x - turbine.x, foundations[zone][p].y - turbine.y) <= radius){
      candidates.push_back(p);
    }
  }

  if(candidates.empty()){
    return;
  }

  int destination = candidates[rand() % candidates.size()];

  solution.grid[zone][turbine.index] = 0;
  solution.grid[zone][destination] = turbine.id;
  turbine.index = destination;
  turbine.x = foundations[zone][destination].x;
  turbine.y = foundations[zone][destination].y;

  calculate_cost(solution);
  calculate_power(solution);

  Solution *sol = new Solution(solution);
  pareto->adicionarSol(sol);

  countRevalue++;

  if(countRevalue % 100000 == 0){
    string path = instance + "_" + algorithm + "_" + to_string(countRevalue) + ".txt";

    pareto->printAllSolutions(root_folder + path);

    if(countRevalue >= stop_criteria){
      pareto->printAllSolutionsLayout(root_folder + instance + "_" + algorithm + "_layout.txt");
    }
  }

  delete sol;
}
