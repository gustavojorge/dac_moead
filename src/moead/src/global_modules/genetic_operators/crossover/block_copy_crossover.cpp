#include "../../../../headers/global_modules/generate_initial_population/population.h"
#include "../../../../headers/global_modules/genetic_operators/block_copy_crossover.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace std;

const double BLOCK_COPY_MIN_FRACTION = 0.1;
const double BLOCK_COPY_MAX_FRACTION = 0.5;

static vector<map<pair<long long, long long>, int>> lattice;
static vector<double> lattice_min_x, lattice_min_y, lattice_step_x, lattice_step_y;

static double random_fraction(){
  return BLOCK_COPY_MIN_FRACTION + (static_cast<double>(rand()) / RAND_MAX) * (BLOCK_COPY_MAX_FRACTION - BLOCK_COPY_MIN_FRACTION);
}

static double smallest_step(vector<double> values){
  sort(values.begin(), values.end());
  double step = 0;
  for(int i = 1; i < (int) values.size(); i++){
    double gap = values[i] - values[i - 1];
    if(gap > 1e-6 && (step == 0 || gap < step)){
      step = gap;
    }
  }
  return step == 0 ? 1.0 : step;
}

static pair<long long, long long> lattice_key(int zone, double x, double y){
  return {llround((x - lattice_min_x[zone]) / lattice_step_x[zone]), llround((y - lattice_min_y[zone]) / lattice_step_y[zone])};
}

static void build_lattice(){
  if((int) lattice.size() == num_zones){
    return;
  }

  lattice.assign(num_zones, {});
  lattice_min_x.assign(num_zones, 0);
  lattice_min_y.assign(num_zones, 0);
  lattice_step_x.assign(num_zones, 1);
  lattice_step_y.assign(num_zones, 1);

  for(int z = 0; z < num_zones; z++){
    if(foundations[z].empty()){
      continue;
    }

    vector<double> xs, ys;
    for(const Foundation &f : foundations[z]){
      xs.push_back(f.x);
      ys.push_back(f.y);
    }

    lattice_min_x[z] = *min_element(xs.begin(), xs.end());
    lattice_min_y[z] = *min_element(ys.begin(), ys.end());
    lattice_step_x[z] = smallest_step(xs);
    lattice_step_y[z] = smallest_step(ys);

    for(int p = 0; p < (int) foundations[z].size(); p++){
      lattice[z][lattice_key(z, foundations[z][p].x, foundations[z][p].y)] = p;
    }
  }
}

static void move_turbine(Turbine &turbine, int zone, int foundation){
  turbine.index = foundation;
  turbine.x = foundations[zone][foundation].x;
  turbine.y = foundations[zone][foundation].y;
}

static void add_turbine(Solution &offspring_solution, Turbine turbine, int zone, int &id_new){
  turbine.id = id_new;
  offspring_solution.grid[zone][turbine.index] = id_new;
  offspring_solution.turbines[zone].push_back(turbine);
  id_new++;
}

static void block_copy_zone(Solution &base_solution, Solution &donor_solution, Solution &offspring_solution, int zone, int &id_new){
  vector<Foundation> &zone_foundations = foundations[zone];
  offspring_solution.grid[zone].assign(zone_foundations.size(), 0);

  if(donor_solution.turbines[zone].empty() || base_solution.turbines[zone].empty() || zone_foundations.size() < 2){
    for(const Turbine &t : base_solution.turbines[zone]){
      add_turbine(offspring_solution, t, zone, id_new);
    }
    return;
  }

  double min_x = zone_foundations[0].x, max_x = zone_foundations[0].x;
  double min_y = zone_foundations[0].y, max_y = zone_foundations[0].y;

  for(const Foundation &f : zone_foundations){
    min_x = min(min_x, f.x);
    max_x = max(max_x, f.x);
    min_y = min(min_y, f.y);
    max_y = max(max_y, f.y);
  }

  Turbine &center = donor_solution.turbines[zone][rand() % donor_solution.turbines[zone].size()];
  double half_width = (max_x - min_x) * random_fraction() / 2.0;
  double half_height = (max_y - min_y) * random_fraction() / 2.0;

  int target;
  do{
    target = rand() % zone_foundations.size();
  } while(target == center.index);

  double shift_x = zone_foundations[target].x - center.x;
  double shift_y = zone_foundations[target].y - center.y;

  vector<bool> occupied(zone_foundations.size(), false);

  vector<Turbine> kept_base;
  for(const Turbine &t : base_solution.turbines[zone]){
    if(fabs(t.x - zone_foundations[target].x) > half_width || fabs(t.y - zone_foundations[target].y) > half_height){
      kept_base.push_back(t);
      occupied[t.index] = true;
    }
  }

  vector<Turbine> pasted;
  for(const Turbine &t : donor_solution.turbines[zone]){
    if(fabs(t.x - center.x) > half_width || fabs(t.y - center.y) > half_height){
      continue;
    }

    auto found = lattice[zone].find(lattice_key(zone, t.x + shift_x, t.y + shift_y));
    if(found == lattice[zone].end() || occupied[found->second]){
      continue;
    }

    Turbine copy = t;
    move_turbine(copy, zone, found->second);
    occupied[found->second] = true;
    pasted.push_back(copy);
  }

  int required = base_solution.turbines[zone].size();

  while((int) (kept_base.size() + pasted.size()) > required){
    int removed = rand() % kept_base.size();
    occupied[kept_base[removed].index] = false;
    kept_base.erase(kept_base.begin() + removed);
  }

  if((int) (kept_base.size() + pasted.size()) < required){
    vector<int> free_foundations;
    for(int p = 0; p < (int) zone_foundations.size(); p++){
      if(!occupied[p]){
        free_foundations.push_back(p);
      }
    }

    while((int) (kept_base.size() + pasted.size()) < required){
      int chosen = rand() % free_foundations.size();
      Turbine filler = base_solution.turbines[zone][0];
      move_turbine(filler, zone, free_foundations[chosen]);
      free_foundations.erase(free_foundations.begin() + chosen);
      kept_base.push_back(filler);
    }
  }

  for(const Turbine &t : kept_base){
    add_turbine(offspring_solution, t, zone, id_new);
  }
  for(const Turbine &t : pasted){
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

Solution block_copy_crossover(Solution &base_solution, Solution &donor_solution){

  build_lattice();

  Solution offspring_solution;
  offspring_solution.grid.resize(num_zones);
  offspring_solution.turbines.resize(num_zones);

  int id_new = fixd.size();

  for(int z = 0; z < num_zones; z++){
    block_copy_zone(base_solution, donor_solution, offspring_solution, z, id_new);
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
