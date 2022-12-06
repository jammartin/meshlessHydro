//
// Created by Johannes Martin on 28.09.22.
//

#ifndef MESHLESSHYDRO_RIEMANN_H
#define MESHLESSHYDRO_RIEMANN_H

#include <cmath>
#include <RiemannSolver.hpp>

#include "Helper.h"
#include "parameter.h"
#include "Logger.h"

class Riemann {

public:
    /// WR, WL and Aij must be pre-allocated
    Riemann(double *WR, double *WL, double *Aij, int i);

    /**
     * Solving a one-dimensional Riemann problem for an ideal gas
     *
     * @param[out] Fij flux to be computed
     */
    void exact(double *Fij, const double &gamma);

private:
    int i;
    double *WR, *WL, *Aij; // *Fij;
    double rhoSol, vSol[DIM], PSol;
    double AijNorm, hatAij[DIM];
    double unitX[DIM] = { 1, 0
#if DIM==3
            ,0
#endif
    };

    void rotateAndProjectFluxes2D(double *Fij, const double &gamma);

};


#endif //MESHLESSHYDRO_RIEMANN_H
