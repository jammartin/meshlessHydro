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
    Riemann(double *WR, double *WL, double *Aij);

    /**
     * Solving a one-dimensional Riemann problem for an ideal gas
     *
     * @param[out] Fij flux to be computed
     */
    void exact(double *Fij, const double &gamma);

private:
    double *WR, *WL, *Aij, *Fij;
    double AijNorm, hatAij[DIM], Lambda[DIM*DIM];

};


#endif //MESHLESSHYDRO_RIEMANN_H
