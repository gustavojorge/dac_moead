#ifndef BOUNDEDPARETOSET_H
#define BOUNDEDPARETOSET_H

#include "../global_modules/generate_initial_population/generate_rSolution.h"
#include "ParetoSet.h"
#include <list>
#include <cstdio>
#include <string>
#include <cmath>
#include <algorithm>

// Classe BoundedParetoSet
class BoundedParetoSet : public ParetoSet {
    private:
        // Tamanho máximo do conjunto Pareto
        static const int MAXARCSIZE = 300;

    public:
        // Construtor e destrutor
        BoundedParetoSet();
        ~BoundedParetoSet();

        // Adiciona uma solução ao conjunto Pareto delimitado
        bool adicionarSol(Solution *s);
};

#endif
