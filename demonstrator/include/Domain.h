//
// Created by Johannes Martin on 21.09.22.
//

#ifndef DEMONSTRATOR_DOMAIN_H
#define DEMONSTRATOR_DOMAIN_H

#include <cmath>
#include <vector>

#include "parameter.h"
#include "Logger.h"

class Domain {

public:
    struct Cell {
        Cell(double *bounds) : minX { bounds[0] }, maxX { bounds[DIM] },
                               minY { bounds[1] }, maxY { bounds[DIM+1] }

#if DIM == 3
                              ,minZ { bounds[2] }, maxZ { bounds[DIM+2] }
#endif
                               {}
        Cell() : minX { 0. }, maxX { 0. },
                 minY { 0. }, maxY { 0. }
#if DIM == 3
                ,minZ { 0.] }, maxZ { 0. }
#endif
                 {}

        double minX;
        double maxX;
        double minY;
        double maxY;
#if DIM == 3
        double minZ;
        double maxZ;
#endif
        std::vector<int> prtcls {};
    };

    Domain(Cell bounds);
    ~Domain();

    int numGridCells { 0 };
    std::vector<Cell> grid;

/*#if PERIODIC_BOUNDARIES
    int numGhostCells { 0 };
    std::vector<Cell> ghostGrid;
    void createGhostGrid();
#endif*/

    void createGrid(const double &kernelSize);
    void getNeighborCells(const int &iCell, int *neighborCells);

    Cell bounds; // global cell
    // number of grid cells in each dimension
    int cellsX { 0 }, cellsY { 0 };
#if DIM == 3
    int cellsZ { 0 };
#endif
    // number of grid cells in each dimension
    double cellSizeX { 0 }, cellSizeY { 0 };
#if DIM == 3
    int cellsSizeZ { 0 };
#endif
private:
    int (*dimIndex)[DIM] { nullptr };
/*#if PERIODIC_BOUNDARIES
    int (*ghostDimIndex)[DIM] { nullptr };
#endif*/
};


#endif //DEMONSTRATOR_DOMAIN_H
