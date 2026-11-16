#ifndef BLOCK_SWAP_CROSSOVER_H
#define BLOCK_SWAP_CROSSOVER_H

#include "../generate_initial_population/population.h"
#include "../../adaptative_grid/BoundedParetoSet.h"

Solution block_swap_crossover(Solution &base_solution, Solution &donor_solution);

#endif
