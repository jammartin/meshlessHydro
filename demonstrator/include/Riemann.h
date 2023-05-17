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
#include "HLLC.h"

class Riemann {

public:
    /// WR, WL and Aij must be pre-allocated
    Riemann(double *WR, double *WL, double *vFrame, double *Aij,
#if USE_HLLC
                                                double *nUnit,
#endif
                                                        int i);

    /**
     * Solving a one-dimensional Riemann problem for an ideal gas
     *
     * @param[out] Fij flux to be computed
     */
    void exact(double *Fij, const double &gamma);

    // HLL approximate Riemann solver, as in Toro, chapter 10.3
    // For != ideal gas
    // void HLL();

    // Method solveHLLC: eqivalent to riemann_solve_for_flux in SWIFT, BUT:
    // Takes velocities in the frame of reference of the interface (already boosted)
    // Returns the de-boosted fluxes that are not yet rotated to the interface frame
    // Also, not yet projected onto the faces

    // @param WL the left state vector
    // @param WR the right state vector
    // @param nUnit normal vector of the interface
    // @param vij velocity vector of the interface
    // @param Fij flux vector to store the solution in

    // void HLLC(const double *WL, const double *WR, const double *nUnit,
    //                                    const double *vij, const double &gamma, double *Fij);
    void HLLCFlux(double *Fij, const double &gamma);

private:
#if USE_HLLC
    double *nUnit;
#endif

    int i;
    double *WR, *WL, *vFrame, *Aij;
    double rhoSol, vSol[DIM], PSol;
    double hatAij[DIM];
    double unitX[DIM] = { 1, 0
#if DIM==3
            ,0
#endif
    };
#if DIM==2
    void rotateAndProjectFluxes2D(double *Fij, const double &gamma);
#else
    void rotateAndProjectFluxes3D(double *Fij, const double &gamma);
#endif
};


#endif //MESHLESSHYDRO_RIEMANN_H
