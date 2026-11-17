#include "../../../../headers/global_modules/generate_initial_population/population.h"
#include "../../../../headers/global_modules/genetic_operators/worst_turbine_mutation.h"
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

static int worst_cost_turbine(vector<Turbine> &zone_turbines, int zone){
  int worst = 0;
  for(int i = 1; i < (int) zone_turbines.size(); i++){
    if(foundations[zone][zone_turbines[i].index].cost > foundations[zone][zone_turbines[worst].index].cost){
      worst = i;
    }
  }
  return worst;
}

static int worst_power_turbine(vector<Turbine> &zone_turbines){
  int worst = 0;
  for(int i = 1; i < (int) zone_turbines.size(); i++){
    if(zone_turbines[i].generated_power < zone_turbines[worst].generated_power){
      worst = i;
    }
  }
  return worst;
}

static bool has_wake_loss(Turbine &turbine){
  double free_wind = wind;
  return turbine.generated_power < power_produced(free_wind, turbine);
}

void worst_turbine_mutation(Solution &solution, double cost_probability){

  int zone = rand() % solution.turbines.size();
  vector<Turbine> &zone_turbines = solution.turbines[zone];
  int zone_size = solution.grid[zone].size();

  if(zone_turbines.empty() || (int) zone_turbines.size() >= zone_size){
    return;
  }

  int worst = 0;
  bool cost_focus = (static_cast<double>(rand()) / RAND_MAX) < cost_probability;

  if(!cost_focus){
    worst = worst_power_turbine(zone_turbines);
    if(!has_wake_loss(zone_turbines[worst])){
      cost_focus = true;
    }
  }

  if(cost_focus){
    worst = worst_cost_turbine(zone_turbines, zone);
  }

  int destination;
  do{
    destination = rand() % zone_size;
  } while(solution.grid[zone][destination] != 0);

  Turbine &turbine = zone_turbines[worst];
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
