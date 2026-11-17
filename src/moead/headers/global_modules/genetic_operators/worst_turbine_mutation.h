#ifndef WORST_TURBINE_MUTATION_H
#define WORST_TURBINE_MUTATION_H

#include "../generate_initial_population/population.h"
#include "../../adaptative_grid/BoundedParetoSet.h"

void worst_turbine_mutation(Solution &solution, double cost_probability);

#endif
