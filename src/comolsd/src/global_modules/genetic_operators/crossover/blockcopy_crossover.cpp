#include "../../../../headers/global_modules/generate_initial_population/population.h"
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>

using namespace std;

static const double BLOCKCOPY_WINDOW_FRACTION = 0.35;
static const double BLOCKCOPY_SNAP_RADIUS = 160.0;
static const int BLOCKCOPY_RANDOM_FREE_ATTEMPTS = 1000;

struct BlockCopyZoneIndex {
  double min_x, max_x, min_y, max_y;
  unordered_map<long long, vector<int>> cells;
};

static vector<BlockCopyZoneIndex> blockcopy_zone_indexes;

static long long blockCopyCellCoord(double value, double origin){
  return (long long) floor((value - origin) / BLOCKCOPY_SNAP_RADIUS);
}

static long long blockCopyCellKey(long long cell_x, long long cell_y){
  return (cell_x << 32) ^ (cell_y & 0xffffffffLL);
}

static void blockCopyBuildZoneIndexes(){
  blockcopy_zone_indexes.assign(num_zones, BlockCopyZoneIndex());

  for(int z = 0; z < num_zones; z++){
    BlockCopyZoneIndex &index = blockcopy_zone_indexes[z];
    vector<Foundation> &zone_foundations = foundations[z];

    if(zone_foundations.empty()){
      continue;
    }

    index.min_x = index.max_x = zone_foundations[0].x;
    index.min_y = index.max_y = zone_foundations[0].y;

    for(Foundation &f : zone_foundations){
      index.min_x = min(index.min_x, f.x);
      index.max_x = max(index.max_x, f.x);
      index.min_y = min(index.min_y, f.y);
      index.max_y = max(index.max_y, f.y);
    }

    for(int p = 0; p < (int) zone_foundations.size(); p++){
      long long cell_x = blockCopyCellCoord(zone_foundations[p].x, index.min_x);
      long long cell_y = blockCopyCellCoord(zone_foundations[p].y, index.min_y);
      index.cells[blockCopyCellKey(cell_x, cell_y)].push_back(p);
    }
  }
}

static bool blockCopyInsideWindow(double x, double y, double center_x, double center_y, double half_width, double half_height){
  return fabs(x - center_x) <= half_width && fabs(y - center_y) <= half_height;
}

static int blockCopyPopRandom(vector<int> &positions){
  int r = rand() % positions.size();
  int p = positions[r];
  positions[r] = positions.back();
  positions.pop_back();
  return p;
}

static int blockCopyNearestFreeFoundation(int zone, double target_x, double target_y, vector<int> &grid){
  BlockCopyZoneIndex &index = blockcopy_zone_indexes[zone];
  long long cell_x = blockCopyCellCoord(target_x, index.min_x);
  long long cell_y = blockCopyCellCoord(target_y, index.min_y);

  int best = -1;
  double best_dist_sq = BLOCKCOPY_SNAP_RADIUS * BLOCKCOPY_SNAP_RADIUS;

  for(long long dx = -1; dx <= 1; dx++){
    for(long long dy = -1; dy <= 1; dy++){
      auto cell = index.cells.find(blockCopyCellKey(cell_x + dx, cell_y + dy));

      if(cell == index.cells.end()){
        continue;
      }

      for(int p : cell->second){
        if(grid[p] != 0){
          continue;
        }

        double diff_x = foundations[zone][p].x - target_x;
        double diff_y = foundations[zone][p].y - target_y;
        double dist_sq = diff_x * diff_x + diff_y * diff_y;

        if(dist_sq < best_dist_sq || (dist_sq == best_dist_sq && (best == -1 || p < best))){
          best_dist_sq = dist_sq;
          best = p;
        }
      }
    }
  }

  return best;
}

static int blockCopyRandomFreeFoundation(int zone, vector<int> &grid, double center_x, double center_y, double half_width, double half_height){
  vector<Foundation> &zone_foundations = foundations[zone];
  int zone_size = zone_foundations.size();

  for(int attempt = 0; attempt < BLOCKCOPY_RANDOM_FREE_ATTEMPTS; attempt++){
    int p = rand() % zone_size;

    if(grid[p] == 0 && !blockCopyInsideWindow(zone_foundations[p].x, zone_foundations[p].y, center_x, center_y, half_width, half_height)){
      return p;
    }
  }

  vector<int> free_outside, free_inside;

  for(int p = 0; p < zone_size; p++){
    if(grid[p] == 0){
      bool inside = blockCopyInsideWindow(zone_foundations[p].x, zone_foundations[p].y, center_x, center_y, half_width, half_height);
      (inside ? free_inside : free_outside).push_back(p);
    }
  }

  return blockCopyPopRandom(free_outside.empty() ? free_inside : free_outside);
}

static void blockCopyZone(Solution &parent_solutionB, Solution &offspring_solution, int zone){
  vector<Foundation> &zone_foundations = foundations[zone];
  vector<Turbine> &donor_turbines = parent_solutionB.turbines[zone];
  vector<int> &grid = offspring_solution.grid[zone];
  int zone_size = zone_foundations.size();
  int target_turbines = turbines_per_zone[zone];

  if(zone_size == 0 || target_turbines == 0 || donor_turbines.empty()){
    return;
  }

  BlockCopyZoneIndex &index = blockcopy_zone_indexes[zone];

  double half_width = BLOCKCOPY_WINDOW_FRACTION * (index.max_x - index.min_x) / 2.0;
  double half_height = BLOCKCOPY_WINDOW_FRACTION * (index.max_y - index.min_y) / 2.0;

  Turbine source_anchor = donor_turbines[rand() % donor_turbines.size()];
  Foundation target_anchor = zone_foundations[rand() % zone_size];

  double shift_x = target_anchor.x - source_anchor.x;
  double shift_y = target_anchor.y - source_anchor.y;

  vector<int> kept, copied;

  for(Turbine &turbine : offspring_solution.turbines[zone]){
    if(blockCopyInsideWindow(turbine.x, turbine.y, target_anchor.x, target_anchor.y, half_width, half_height)){
      grid[turbine.index] = 0;
    } else {
      kept.push_back(turbine.index);
    }
  }

  for(Turbine &turbine : donor_turbines){
    if(!blockCopyInsideWindow(turbine.x, turbine.y, source_anchor.x, source_anchor.y, half_width, half_height)){
      continue;
    }

    int best = blockCopyNearestFreeFoundation(zone, turbine.x + shift_x, turbine.y + shift_y, grid);

    if(best != -1){
      grid[best] = -1;
      copied.push_back(best);
    }
  }

  int active = kept.size() + copied.size();

  while(active > target_turbines){
    int p = blockCopyPopRandom(kept.empty() ? copied : kept);
    grid[p] = 0;
    active--;
  }

  while(active < target_turbines){
    int p = blockCopyRandomFreeFoundation(zone, grid, target_anchor.x, target_anchor.y, half_width, half_height);
    grid[p] = -1;
    kept.push_back(p);
    active++;
  }

  Turbine turbine = donor_turbines[0];
  turbine.zone = zone;

  offspring_solution.turbines[zone].clear();

  for(vector<int> *positions : {&kept, &copied}){
    for(int p : *positions){
      turbine.index = p;
      turbine.x = zone_foundations[p].x;
      turbine.y = zone_foundations[p].y;
      offspring_solution.turbines[zone].push_back(turbine);
    }
  }
}

static void blockCopyRenumberTurbines(Solution &solution){
  int next_id = fixd.size();

  for(int z = 0; z < num_zones; z++){
    for(Turbine &turbine : solution.turbines[z]){
      turbine.id = next_id;
      solution.grid[z][turbine.index] = next_id;
      next_id++;
    }
  }
}

static void blockCopyRegisterEvaluation(Solution &offspring_solution){
  calculate_cost(offspring_solution);
  calculate_power(offspring_solution);

  Solution * sol = new Solution(offspring_solution);
  pareto->adicionarSol(sol);

  delete sol;

  countRevalue++;

  if(countRevalue % 100000 == 0){
    string path = instance + "_" + algorithm + "_" + to_string(countRevalue) + ".txt";
    pareto->printAllSolutions(root_folder + path);
  }
}

Solution block_copy_crossover(Solution &parent_solutionA, Solution &parent_solutionB){
  if((int) blockcopy_zone_indexes.size() != num_zones){
    blockCopyBuildZoneIndexes();
  }

  Solution offspring_solution = parent_solutionA;

  for(int z = 0; z < num_zones; z++){
    blockCopyZone(parent_solutionB, offspring_solution, z);
  }

  blockCopyRenumberTurbines(offspring_solution);

  blockCopyRegisterEvaluation(offspring_solution);

  return offspring_solution;
}
