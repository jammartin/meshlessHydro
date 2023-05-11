//
// Created by Jakob Sappler on 09.05.23.
//

#ifndef MESHLESS_HYDRO_HLLC
#define MESHLESS_HYDRO_HLLC

#include "parameter.h"
#include "Helper.h"


class HLLC {

public:
    // Method solve: eqivalent to riemann_solve_for_flux in SWIFT, BUT:
    // Takes velocities in the frame of reference of the interface (already boosted)
    // Returns the de-boosted fluxes that are not yet rotated to the interface frame
    // Also, not yet projected onto the faces

    // @param WL the left state vector
    // @param WR the right state vector
    // @param n_unit normal vector of the interface
    // @param vij velocity vector of the interface
    // @param Fij flux vector to store the solution in

    void solve(const double *WL, const double *WR, const double *n_unit,
                                        const double *vij, double *Fij);
}
#endif
