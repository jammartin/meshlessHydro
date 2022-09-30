//
// Created by Johannes Martin on 28.09.22.
//

#include "../include/Riemann.h"

Riemann::Riemann(double *WR, double *WL, double *Aij) : WR { WR }, WL { WL }, Aij { Aij }{

    // compute norm of effective face
    AijNorm = sqrt(Helper::dotProduct(Aij, Aij));

    hatAij[0] = 1./AijNorm*Aij[0];
    hatAij[1] = 1./AijNorm*Aij[1];
#if DIM==3
    hatAij[2] = 1./AijNorm*Aij[2];
#endif

    // compute norm of velocity
    double unitX[DIM] = { 1, 0
#if DIM==3
                          ,0
#endif
    };

#if DIM==3
    Logger(ERROR) << "Rotation of states is not yet implemented in 3D. - Aborting.";
    exit(4);
#endif
    // rotate states to be aligned with hatAij
    Helper::rotationMatrix2D(unitX, hatAij, Lambda);

    WR[2] = Lambda[0]*WR[2]+Lambda[1]*WR[3];
    WR[3] = Lambda[2]*WR[2]+Lambda[3]*WR[3];
    WL[2] = Lambda[0]*WL[2]+Lambda[1]*WL[3];
    WL[3] = Lambda[2]*WL[2]+Lambda[3]*WL[3];

    // TODO: add z-components to make it work for 3D
}

void Riemann::exact(double *Fij, const double &gamma){
    RiemannSolver solver { gamma };

    Logger(DEBUG) << "rhoL = " << WL[0] << ", rhoR = " << WR[0]
                  << ", uL = " << WL[2] << ", uR = " << WR[2]
                  << ", PL = " << WL[1] << ", PR = " << WR[1];

    solver.solve(WL[0], WL[2], WL[1], WR[0], WR[2], WR[1],
                 Fij[0], Fij[2], Fij[1]);
}

