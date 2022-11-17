//
// Created by Johannes Martin on 28.09.22.
//

#include "../include/Riemann.h"

Riemann::Riemann(double *WR, double *WL, double Aij[DIM]) : WR { WR }, WL { WL }, Aij { Aij }{

    // compute norm of effective face
    AijNorm = sqrt(Helper::dotProduct(Aij, Aij));

    hatAij[0] = 1./AijNorm*Aij[0];
    hatAij[1] = 1./AijNorm*Aij[1];
#if DIM==3
    hatAij[2] = 1./AijNorm*Aij[2];
#endif

    //Logger(DEBUG) << "Aij = [" << Aij[0] << ", " << Aij[1] << "], AijNorm = " << AijNorm;

#if DIM==3
    Logger(ERROR) << "Rotation of states is not yet implemented in 3D. - Aborting.";
    exit(4);
#endif
    // rotate states to be aligned with hatAij
    double Lambda[DIM*DIM];
    Helper::rotationMatrix2D(unitX, hatAij, Lambda);

    //Logger(DEBUG) << "hatAij = [" << hatAij[0] << ", " << hatAij[1] << "], |hatAij| = "
    //          << sqrt(Helper::dotProduct(hatAij, hatAij));
    //Logger(DEBUG) << "Lambda = [" << Lambda[0] << ", " << Lambda[1];
    //Logger(DEBUG) << "          " << Lambda[2] << ", " << Lambda[3] << "]";

    WR[2] = Lambda[0]*WR[2]+Lambda[1]*WR[3];
    WR[3] = Lambda[2]*WR[2]+Lambda[3]*WR[3];
    WL[2] = Lambda[0]*WL[2]+Lambda[1]*WL[3];
    WL[3] = Lambda[2]*WL[2]+Lambda[3]*WL[3];

    // TODO: add z-components to make it work for 3D
}

void Riemann::exact(double *Fij, const double &gamma){
    RiemannSolver solver { gamma };

    //if(WR[1] < 0. || WL[1] < 0.){
    //    Logger(ERROR) << "Negative pressure encountered. Very bad :( !!";
    //    Logger(DEBUG) << "rhoL = " << WL[0] << ", rhoR = " << WR[0]
    //                  << ", uL = " << WL[2] << ", uR = " << WR[2]
    //                  << ", PL = " << WL[1] << ", PR = " << WR[1];
    //}

    int flagLR = solver.solve(WL[0], WL[2], WL[1], WR[0], WR[2], WR[1],
                              rhoSol, vSol[0], PSol);

    // Fij[3] = 0.; // explicitly setting vy flux to zero
    //if (Fij[2] < 0.){
    //    Logger(DEBUG) << "Fij = [" << Fij[0] << " (mass), " << Fij[2] << " (vx), "
    //              << Fij[3] << " (vy), " << Fij[1] << " (energy)]";
    //}

    /// CONVERT SOLUTION TO FLUXES
#if DIM==3
    // TODO: Implement 3D case
    Logger(ERROR) << "Flux computation from Riemann problem's solution not implemented in 3D. - Aborting.";
    exit(5);
#else
    if (flagLR == 1){
        // right state sampled
        vSol[1] = WR[3];
    } else { // flagLR == -1
        // left state sampled
        vSol[1] = WL[3];
    }
    rotateAndProjectFluxes2D(Fij, gamma);
#endif

    //Logger(DEBUG) << "Fij = [" << Fij[0] << " (mass), " << Fij[2] << " (vx), "
    //          << Fij[3] << " (vy), " << Fij[1] << " (energy)]";
}

void Riemann::rotateAndProjectFluxes2D(double *Fij, const double &gamma){

    // rotate back velocities to simulation coordinate frame
    double LambdaInv[DIM*DIM];
    Helper::rotationMatrix2D(hatAij, unitX, LambdaInv);

    //Logger(DEBUG) << "LambdaInv = [" << LambdaInv[0] << ", " << LambdaInv[1]
    //          << ", " << LambdaInv[2] << ", " << LambdaInv[3] << "]";

    vSol[0] = LambdaInv[0]*vSol[0]+LambdaInv[1]*vSol[1];
    vSol[1] = LambdaInv[2]*vSol[0]+LambdaInv[3]*vSol[1];

    /// Compute fluxes projected in Aij direction
    Fij[0] = hatAij[0]*rhoSol*vSol[0] + hatAij[1]*rhoSol*vSol[1]; // mass flux
    Fij[2] = hatAij[0]*(rhoSol*vSol[0]*vSol[0]+PSol) + hatAij[1]*rhoSol*vSol[0]*vSol[1]; // vx flux
    Fij[3] = hatAij[0]*rhoSol*vSol[0]*vSol[1] + hatAij[1]*(rhoSol*vSol[1]*vSol[1]+PSol); // vy flux
    Fij[1] = hatAij[0]*vSol[0]*(PSol/(gamma-1.)+rhoSol*.5*Helper::dotProduct(vSol, vSol) + PSol)
            + hatAij[1]*vSol[1]*(PSol/(gamma-1.)+rhoSol*.5*Helper::dotProduct(vSol, vSol) + PSol); // energy flux

}

