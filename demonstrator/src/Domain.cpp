//
// Created by Johannes Martin on 21.09.22.
//

#include "../include/Domain.h"

Domain::Domain(Cell bounds) : bounds { bounds }{}

void Domain::createGrid(const double &kernelSize){
    cellsX = floor((bounds.maxX - bounds.minX)/kernelSize);
    cellsY = floor((bounds.maxY - bounds.minY)/kernelSize);
#if DIM == 3
    cellsZ = floor((bounds.maxZ - bounds.minZ)/kernelSize);
#endif

#if DIM == 2
    numGridCells = cellsX*cellsY;
#else
    numGridCells = cellsX*cellsY*cellsZ;
#endif

    cellSizeX = (bounds.maxX - bounds.minX)/(double)cellsX;
    cellSizeY = (bounds.maxY - bounds.minY)/(double)cellsY;
    Logger(DEBUG) << "      > cellSizeX = " << cellSizeX << ", cellsX = " << cellsX;
    Logger(DEBUG) << "      > cellSizeY = " << cellSizeY<< ", cellsY = " << cellsY;
#if DIM == 3
    cellSizeZ = (bounds.maxZ - bounds.minZ)/(double)cellsZ;
    Logger(DEBUG) << "      > cellSizeZ = " << cellSizeZ << ", cellsZ = " << cellsZ;
#endif

    grid = std::vector<Cell>(numGridCells);
    dimIndex = new int[numGridCells][DIM];

    for(int iX=0; iX<cellsX; ++iX){
        for(int iY=0; iY<cellsY; ++iY){
#if DIM == 2
            double cellBounds[] = { iX*cellSizeX+bounds.minX, iY*cellSizeY+bounds.minY,
                                    (iX+1)*cellSizeX+bounds.minX, (iY+1)*cellSizeY+bounds.minY };
            grid[iX+iY*cellsX] = Cell(cellBounds);
            dimIndex[iX+iY*cellsX][0] = iX;
            dimIndex[iX+iY*cellsX][1] = iY;
#else
            for(int iZ=0; iZ<cellsZ; ++iZ){
                double cellBounds[] = { iX*cellSizeX+bounds.minX, iY*cellSizeY+bounds.minY, iZ*cellSizeZ+bounds.minZ,
                                        (iX+1)*cellSizeX+bounds.minX, (iY+1)*cellSizeY+bounds.minY, (iZ+1)*cellSizeZ+bounds.minZ };
                grid[iX+iY*cellsX+iZ*cellsX*cellsY] = Cell(cellBounds);
                dimIndex[iX+iY*cellsX+iZ*cellsX*cellsY][0] = iX;
                dimIndex[iX+iY*cellsX+iZ*cellsX*cellsY][1] = iY;
                dimIndex[iX+iY*cellsX+iZ*cellsX*cellsY][2] = iZ;
            }
#endif
        }
    }
}

/*#if PERIODIC_BOUNDARIES
void Domain::createGhostGrid(){
    numGhostCells = DIM*(cellsX+cellsY+DIM);
    ghostGrid = std::vector<Cell>(numGhostCells);
    ghostDimIndex = new int[numGhostCells][DIM];

    int iGhost = 0;
#if DIM ==2
    for(int iY=-1; iY<=cellsY; ++iY){
        for(int iX=-1; iX<=cellsX; ++iX){
            if(iY==-1 || iY==cellsY || iX==-1 || iX==cellsX){
                double cellBounds[] = { iX*cellSizeX+bounds.minX, iY*cellSizeY+bounds.minY,
                                        (iX+1)*cellSizeX+bounds.minX, (iY+1)*cellSizeY+bounds.minY };
                ghostGrid[iGhost] = Cell(cellBounds);
                ghostDimIndex[iGhost][0] = iX;
                ghostDimIndex[iGhost][1] = iY;
                ++iGhost;
            }
        }
    }
#else
    Logger(ERROR) << "Ghost cells not implemented for 3D simulations. - Aborting.";
    exit(2);
#endif
}
#endif*/

void Domain::getNeighborCells(const int &iCell, int *neighborCell){
    // recover iX, iY and iZ
    int iX = dimIndex[iCell][0];
    int iY = dimIndex[iCell][1];
#if DIM == 3
    int iZ = dimIndex[iCell][2];
#endif
    int iNeighbor = 0;
    for(int k=iX-1; k<=iX+1; ++k){
        for(int l=iY-1; l<=iY+1; ++l){
#if DIM == 2
            if (k < 0 || k >= cellsX){
                neighborCell[iNeighbor] = -1; // -1: no cell or ghost
            } else if (l < 0 || l >= cellsY) {
                neighborCell[iNeighbor] = -1;
            } else {
                neighborCell[iNeighbor] = k+l*cellsX;
            }
            ++iNeighbor;
#else
            for(int m=iZ-1; m<=iZ+1; ++m){
                if (k < 0 || k >= cellsX){
                    neighborCell[iNeighbor] = -1; // -1: no cell or ghost
                } else if (l < 0 || l >= cellsY) {
                    neighborCell[iNeighbor] = -1;
                } else if (m < 0 || m >= cellsZ) {
                    neighborCell[iNeighbor] = -1;
                } else {
                    neighborCell[iNeighbor] = k+l*cellsX+m*cellsX*cellsY;
                }
                ++iNeighbor;
            }
#endif
        }
    }
}

void Domain::printout(){
    Logger(INFO) << "Domain > X [" << bounds.minX << ", " << bounds.maxX << "]";
    Logger(INFO) << "         Y [" << bounds.minY << ", " << bounds.maxY << "]";
#if DIM == 3
    Logger(INFO) << "         Z [" << bounds.minZ << ", " << bounds.maxZ << "]";
#endif
}

Domain::~Domain(){
    delete[] dimIndex;
}