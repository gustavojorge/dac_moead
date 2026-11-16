#ifndef PARETOSET_H
#define PARETOSET_H

#include <list>
#include <map>
#include <cmath>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include "param.h"
#include "Grid.h"
#include "../global_modules/generate_initial_population/generate_rSolution.h"

typedef struct {
    double min, max;
} range;

double getObj(Solution *s, int obj);
bool dominatesP(Solution &s1, Solution &s2);
bool equals(Solution &s1, Solution &s2);
void printSolucao(FILE* f, Solution *s);
void printLayout(std::ofstream &file, Solution *s);

class ParetoSet {
    protected:
        std::list<Solution *> sol;
        range rangeNovo[2], rangeAtual[2];
        int frequencia[NUMEROVERTICES][NUMEROVERTICES];
        Grid g;

        int calcularGridPos(Solution &s);
        void updateGrid();
        void reiniciarRanges();

    public:
        ParetoSet();
        virtual ~ParetoSet();

        int getPositionCount(Solution &s);
        int getPositionCount(int p);
        std::list<Solution *> getElementos();
        
        virtual bool adicionarSol(Solution *s);
        void printAllSolutions(std::string path);
        void printAllSolutionsLayout(std::string path);
        int getSize();
        Solution *getSolucao(int p);
        void clear();
        bool confereGrid();
};

#endif
