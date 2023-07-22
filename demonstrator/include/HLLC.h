//
// Created by Jakob Sappler on 09.05.23.
//

#ifndef MESHLESS_HYDRO_HLLC
#define MESHLESS_HYDRO_HLLC

#include "parameter.h"
#include "Helper.h"
//#include "riemannhelper.h"

class HLLC {

public:

    //HLLC(double *WR, double *WL, double *vFrame, double *Aij, double *n_unit, int i);

    // Method solve: eqivalent to riemann_solve_for_flux in SWIFT, BUT:
    // Takes data in the frame of reference of the interface (already boosted)
    // Returns the boosted fluxes that are not yet rotated to the interface frame
    // Also, not yet projected onto the faces

    // @param rhoL, uL, PL ... rho, u and P on the left hand side.
    // @param n_unit normal vector of the interface
    // @param Fij flux vector to store the solution in

    // static void solveHLLC(double rhoL, double uL, double PL, double rhoR, double uR, double PR,
    //                     const double &gamma, double *Fij);


    // Second try
    // Method solve: eqivalent to riemann_solve_for_flux in SWIFT, BUT:
    // static void solveHLLC0(double *WL, double *WR, double *n, double *totflux,
    //         const double *vij, const double &hydro_gamma);

    static void solveHLLC1(double *WL, double *WR, double *n, double *totflux,
            const double *vij, const double &hydro_gamma);

    // HLL Solver
    static void HLL(double *WL, double *WR, double *totflux, const double &hydro_gamma);
    
// private:
//     int j;
};
#endif
