#include "../../../../headers/global_modules/generate_initial_population/population.h"
#include "../../../../headers/global_modules/genetic_operators/block_swap_crossover.h"
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace std;

const int BLOCK_SWAP_MAX_ATTEMPTS = 10;
const double BLOCK_SWAP_MIN_FRACTION = 0.1;
const double BLOCK_SWAP_MAX_FRACTION = 0.5;

static double random_fraction(){
  return BLOCK_SWAP_MIN_FRACTION + (static_cast<double>(rand()) / RAND_MAX) * (BLOCK_SWAP_MAX_FRACTION - BLOCK_SWAP_MIN_FRACTION);
}

static void add_turbine(Solution &offspring_solution, Turbine turbine, int zone, int &id_new){
  turbine.id = id_new;
  offspring_solution.grid[zone][turbine.index] = id_new;
  offspring_solution.turbines[zone].push_back(turbine);
  id_new++;
}

static void block_swap_zone(Solution &base_solution, Solution &donor_solution, Solution &offspring_solution, int zone, int &id_new){
  vector<Foundation> &zone_foundations = foundations[zone];
  offspring_solution.grid[zone].assign(zone_foundations.size(), 0);

  if(!donor_solution.turbines[zone].empty()){
    double min_x = zone_foundations[0].x, max_x = zone_foundations[0].x;
    double min_y = zone_foundations[0].y, max_y = zone_foundations[0].y;

    for(const Foundation &f : zone_foundations){
      min_x = min(min_x, f.x);
      max_x = max(max_x, f.x);
      min_y = min(min_y, f.y);
      max_y = max(max_y, f.y);
    }

    double width = max_x - min_x;
    double height = max_y - min_y;

    for(int attempt = 0; attempt < BLOCK_SWAP_MAX_ATTEMPTS; attempt++){
      Turbine &center = donor_solution.turbines[zone][rand() % donor_solution.turbines[zone].size()];
      double half_width = width * random_fraction() / 2.0;
      double half_height = height * random_fraction() / 2.0;

      auto inside_block = [&](const Turbine &t){
        return fabs(t.x - center.x) <= half_width && fabs(t.y - center.y) <= half_height;
      };

      vector<Turbine> outside_base;
      int inside_base = 0;
      for(const Turbine &t : base_solution.turbines[zone]){
        if(inside_block(t)){
          inside_base++;
        } else {
          outside_base.push_back(t);
        }
      }

      vector<Turbine> inside_donor;
      for(const Turbine &t : donor_solution.turbines[zone]){
        if(inside_block(t)){
          inside_donor.push_back(t);
        }
      }

      if((int) inside_donor.size() < inside_base){
        continue;
      }

      int excess = inside_donor.size() - inside_base;
      for(int i = 0; i < excess; i++){
        outside_base.erase(outside_base.begin() + rand() % outside_base.size());
      }

      for(const Turbine &t : outside_base){
        add_turbine(offspring_solution, t, zone, id_new);
      }
      for(const Turbine &t : inside_donor){
        add_turbine(offspring_solution, t, zone, id_new);
      }
      return;
    }
  }

  for(const Turbine &t : base_solution.turbines[zone]){
    add_turbine(offspring_solution, t, zone, id_new);
  }
}

static bool same_layout(Solution &a, Solution &b){
  for(int z = 0; z < num_zones; z++){
    for(int p = 0; p < (int) a.grid[z].size(); p++){
      if((a.grid[z][p] != 0) != (b.grid[z][p] != 0)){
        return false;
      }
    }
  }
  return true;
}

Solution block_swap_crossover(Solution &base_solution, Solution &donor_solution){

  Solution offspring_solution;
  offspring_solution.grid.resize(num_zones);
  offspring_solution.turbines.resize(num_zones);

  int id_new = fixd.size();

  for(int z = 0; z < num_zones; z++){
    block_swap_zone(base_solution, donor_solution, offspring_solution, z, id_new);
  }

  if(same_layout(offspring_solution, base_solution)){
    return base_solution;
  }

  calculate_cost(offspring_solution);
  calculate_power(offspring_solution);

  Solution * sol = new Solution(offspring_solution);
  pareto->adicionarSol(sol);

  delete sol;

  countRevalue++;

  if(countRevalue % 100000 == 0){
    string path = instance + "_" + algorithm + "_" + to_string(countRevalue) + ".txt";
    pareto->printAllSolutions(root_folder + path);

    if(countRevalue >= stop_criteria){
      pareto->printAllSolutionsLayout(root_folder + instance + "_" + algorithm + "_layout.txt");
    }
  }

  return offspring_solution;
}
