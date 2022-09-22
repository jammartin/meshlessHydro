//
// Created by Johannes Martin on 21.09.22.
//

#include "../include/Domain.h"

Domain::Domain(Cell bounds) : bounds { bounds }{}

void Domain::createGrid(double kernelSize){
    cellsX = ceil((bounds.maxX - bounds.minX)/kernelSize);
    cellsY = ceil((bounds.maxY - bounds.minY)/kernelSize);
#if DIM == 3
    cellsZ = ceil((bounds.maxZ - bounds.minZ)/kernelSize);
#endif

#if DIM == 2
    numGridCells = cellsX*cellsY;
#else
    numGridCells = cellsX*cellsY*cellsZ;
#endif

    cellSizeX = (bounds.maxX - bounds.minX)/(double)cellsX;
    cellSizeY = (bounds.maxY - bounds.minY)/(double)cellsY;
    Logger(DEBUG) << "      > cellSizeX = " << cellSizeX;
    Logger(DEBUG) << "      > cellSizeY = " << cellSizeY;
#if DIM == 3
    cellSizeZ = (domain.bounds.maxZ - domain.bounds.minZ)/(double)cellsZ;
    Logger(DEBUG) << "      > cellSizeZ = " << cellSizeZ;
#endif

    grid = std::vector<Cell>(numGridCells);
    dimIndex = new int[numGridCells][DIM];

    for(int iX=0; iX<cellsX; ++iX){
        for(int iY=0; iY<cellsY; ++iY){
#if DIM == 2
            double cellBounds[] = { iX*cellSizeX + (iX+1)*cellSizeX,
                                    iY*cellSizeY + (iY+1)*cellSizeY };
            grid[iX+iY*cellsX] = Cell(cellBounds);
            dimIndex[iX+iY*cellsX][0] = iX;
            dimIndex[iX+iY*cellsX][1] = iY;
#else
            for(int iZ=0; iZ<cellsZ; ++iZ){
                double cellBounds[] = { iX*cellSizeX + (iX+1)*cellSizeX,
                                        iY*cellSizeY + (iY+1)*cellSizeY,
                                        iZ*cellSizeZ + (iZ+1)*cellSizeZ };
                grid[iX+iY*cellsX+iZ*cellsX*cellsY] = Cell(cellBounds);
                dimIndex[iX+iY*cellsX+iZ*cellsX*cellsY][0] = iX;
                dimIndex[iX+iY*cellsX+iZ*cellsX*cellsY][1] = iY;
                dimIndex[iX+iY*cellsX+iZ*cellsX*cellsY][2] = iZ;
            }
#endif
        }
    }
}

void Domain::getNeighborCells(int iCell, int *neighborCell){
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
#else
            for(int m=-1; m<=1; ++m){
                if (k < 0 || k >= cellsX){
                    neighborCell[iNeighbor] = -1; // -1: no cell or ghost
                } else if (l < 0 || l >= cellsY) {
                    neighborCell[iNeighbor] = -1;
                } else if (m < 0 || m >= cellsZ)
                    neighborCell[iNeighbor] = -1;
                } else {
                    neighborCell[iNeighbor] = k+l*cellsX+m*cellsX*cellsY;
                }
            }
#endif
        ++iNeighbor;
        }
    }
}

Domain::~Domain(){
    delete[] dimIndex;
}