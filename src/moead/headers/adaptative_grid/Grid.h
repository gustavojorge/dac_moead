#ifndef GRID_H
#define GRID_H

#include "./param.h"
#include <cstring>

class Grid {
private:
    const static int tam = 1 << (NUMOBJETIVOS * PROFUNDIDADEGRID);
    int grid[1 << (NUMOBJETIVOS * PROFUNDIDADEGRID)];

public:
    Grid();  

    int getPositionCount(int p);

    void addGrid(int p);

    void removeGrid(int p);

    void clearGrid();

    int getSize();
};

#endif
