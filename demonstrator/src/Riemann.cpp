//
// Created by Johannes Martin on 28.09.22.
//

#include "../include/Riemann.h"

Riemann::Riemann(double *WR, double *WL, double *vFrame, double *Aij, int i) :
                    WR { WR }, WL { WL }, vFrame { vFrame },  Aij { Aij }, i { i }{

    // compute norm of effective face
    double AijNorm = sqrt(Helper::dotProduct(Aij, Aij));

    hatAij[0] = 1./AijNorm*Aij[0];
    hatAij[1] = 1./AijNorm*Aij[1];
#if DIM == 3
    hatAij[2] = 1./AijNorm*Aij[2];
#endif

    //if (i == 6){
    //    Logger(DEBUG) << "Aij = [" << Aij[0] << ", " << Aij[1] << "], AijNorm = " << AijNorm;
    //}

//#if DIM == 3
//    Logger(ERROR) << "Rotation of states is not yet implemented in 3D. - Aborting.";
//    exit(4);
//#endif
    // rotate states to be aligned with hatAij
    double Lambda[DIM*DIM];
    //Helper::rotationMatrix2D(unitX, hatAij, Lambda);

#if DIM == 2
    Helper::rotationMatrix2D(hatAij, unitX, Lambda);

    double vBufR[DIM] = { WR[2], WR[3] };
    double vBufL[DIM] = { WL[2], WL[3] };

    WR[2] = Lambda[0]*vBufR[0]+Lambda[1]*vBufR[1];
    WR[3] = Lambda[2]*vBufR[0]+Lambda[3]*vBufR[1];
    WL[2] = Lambda[0]*vBufL[0]+Lambda[1]*vBufL[1];
    WL[3] = Lambda[2]*vBufL[0]+Lambda[3]*vBufL[1];

    //Logger(DEBUG) << "hatAij = [" << hatAij[0] << ", " << hatAij[1] << "], |hatAij| = "
    //          << sqrt(Helper::dotProduct(hatAij, hatAij));
    //Logger(DEBUG) << "Lambda = [" << Lambda[0] << ", " << Lambda[1];
    //Logger(DEBUG) << "          " << Lambda[2] << ", " << Lambda[3] << "]";

    //Logger(DEBUG) << "vR_unrot = [" << WR[2] << ", " << WR[3] << "], vL_unrot = ["
    //          <<  WL[2] << ", " << WL[3] << "]";

    //if(i == 204){
    //    Logger(DEBUG) <<  "rhoR = " << WR[0] << ", rhoL = " << WL[0]
    //                  << ", vR = [" << WR[2] << ", " << WR[3]
    //                  << "], vL = [" <<  WL[2] << ", " << WL[3]
    //                  << "], PR = " << WR[1] << ", PL = " << WL[1];
    //}

    //if(i == 46){
    //    Logger(DEBUG) <<  "rhoR = " << WR[0] << ", rhoL = " << WL[0]
    //                  << ", vR = [" << WR[2] << ", " << WR[3]
    //                  << "], vL = [" <<  WL[2] << ", " << WL[3]
    //                  << "], PR = " << WR[1] << ", PL = " << WL[1];
    //}
#else
    Helper::rotationMatrix3D(hatAij, unitX, Lambda);

    //Logger(DEBUG) << "Lambda = [" << Lambda[0] << ", " << Lambda[1] << ", " << Lambda[2];
    //Logger(DEBUG) << "          " << Lambda[3] << ", " << Lambda[4] << ", " << Lambda[5];
    //Logger(DEBUG) << "          " << Lambda[6] << ", " << Lambda[7] << ", " << Lambda[8] << "]";

    double vBufR[DIM] = { WR[2], WR[3], WR[4] };
    double vBufL[DIM] = { WL[2], WL[3], WL[4] };

    WR[2] = Lambda[0]*vBufR[0]+Lambda[1]*vBufR[1]+Lambda[2]*vBufR[2];
    WR[3] = Lambda[3]*vBufR[0]+Lambda[4]*vBufR[1]+Lambda[5]*vBufR[2];
    WR[4] = Lambda[6]*vBufR[0]+Lambda[7]*vBufR[1]+Lambda[8]*vBufR[2];

    WL[2] = Lambda[0]*vBufL[0]+Lambda[1]*vBufL[1]+Lambda[2]*vBufL[2];
    WL[3] = Lambda[3]*vBufL[0]+Lambda[4]*vBufL[1]+Lambda[5]*vBufL[2];
    WL[4] = Lambda[6]*vBufL[0]+Lambda[7]*vBufL[1]+Lambda[8]*vBufL[2];
#endif
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
    if (flagLR == 1){
        // right state sampled
        vSol[1] = WR[3];
        vSol[2] = WR[4];
    } else if (flagLR == -1){ // flagLR == -1
        // left state sampled
        vSol[1] = WL[3];
        vSol[2] = WL[4];
    } else { // flagLR == 0
        Logger(WARN) << "  > Vacuum state sampled. This is not expected.";
    }

    rotateAndProjectFluxes3D(Fij, gamma);

#else
    if (flagLR == 1){
        // right state sampled
        vSol[1] = WR[3];
    } else if (flagLR == -1){ // flagLR == -1
        // left state sampled
        vSol[1] = WL[3];
    } else { // flagLR == 0
        Logger(WARN) << "  > Vacuum state sampled. This is not expected.";
    }

    //if (i == 46){
    //    Logger(DEBUG) << "riemann solution = [" << rhoSol << ", " << vSol[0] << ", " << vSol[1] << ", " << PSol << "]";
    //}

    rotateAndProjectFluxes2D(Fij, gamma);

    //if (i == 46){
    //    Logger(DEBUG) << "Fij = [" << Fij[0] << " (mass), " << Fij[2] << " (vx), "
    //                  << Fij[3] << " (vy), " << Fij[1] << " (energy)]";
    //}

#endif
}

void Riemann::rotateAndProjectFluxes2D(double *Fij, const double &gamma){

    // rotate back velocities to simulation coordinate frame
    double LambdaInv[DIM*DIM];
    //Helper::rotationMatrix2D(hatAij, unitX, LambdaInv);
    Helper::rotationMatrix2D(unitX, hatAij, LambdaInv);

    //Logger(DEBUG) << "vSol_unrot = [" << vSol[0] << ", " << vSol[1] << "]";
    //Logger(DEBUG) << "LambdaInv = [" << LambdaInv[0] << ", " << LambdaInv[1]
    //          << ", " << LambdaInv[2] << ", " << LambdaInv[3] << "]";

    double vSolBuf[DIM] = { vSol[0], vSol[1] };
    vSol[0] = LambdaInv[0]*vSolBuf[0]+LambdaInv[1]*vSolBuf[1];
    vSol[1] = LambdaInv[2]*vSolBuf[0]+LambdaInv[3]*vSolBuf[1];

    //if (i == 204){
    //    Logger(DEBUG) << "rhoSol = " << rhoSol << ", vSol_rot = [" << vSol[0] << ", " << vSol[1] << "], PSol = " << PSol;
    //}

    /// Compute fluxes projected onto Aij
    Fij[0] = Aij[0]*rhoSol*vSol[0] + Aij[1]*rhoSol*vSol[1]; // mass flux

    // velocity solution in the lab frame
    double vLab[DIM];
    vLab[0] = vSol[0] + vFrame[0];
    vLab[1] = vSol[1] + vFrame[1];

    Fij[2] = Aij[0]*(rhoSol*vLab[0]*vSol[0]+PSol) + Aij[1]*rhoSol*vLab[0]*vSol[1]; // vx flux
    Fij[3] = Aij[0]*rhoSol*vLab[1]*vSol[0] + Aij[1]*(rhoSol*vLab[1]*vSol[1]+PSol); // vy flux

    Fij[1] = Aij[0]*(vSol[0]*(PSol/(gamma-1.)+rhoSol*.5*Helper::dotProduct(vLab, vLab)) + PSol*vLab[0])
            + Aij[1]*(vSol[1]*(PSol/(gamma-1.)+rhoSol*.5*Helper::dotProduct(vLab, vLab)) + PSol*vLab[1]); // energy flux

    //double vFluxBuf[DIM] = { Fij[2], Fij[3] };
    //Fij[2] = LambdaInv[0]*vFluxBuf[0]+LambdaInv[1]*vFluxBuf[1];
    //Fij[3] = LambdaInv[2]*vFluxBuf[0]+LambdaInv[3]*vFluxBuf[1];

}

void Riemann::rotateAndProjectFluxes3D(double *Fij, const double &gamma){

    // rotate back velocities to simulation coordinate frame
    double LambdaInv[DIM*DIM];
    Helper::rotationMatrix3D(unitX, hatAij, LambdaInv);

    double vSolBuf[DIM] = { vSol[0], vSol[1], vSol[2] };
    vSol[0] = LambdaInv[0]*vSolBuf[0]+LambdaInv[1]*vSolBuf[1]+LambdaInv[2]*vSolBuf[2];
    vSol[1] = LambdaInv[3]*vSolBuf[0]+LambdaInv[4]*vSolBuf[1]+LambdaInv[5]*vSolBuf[2];
    vSol[2] = LambdaInv[6]*vSolBuf[0]+LambdaInv[7]*vSolBuf[1]+LambdaInv[8]*vSolBuf[2];

    //if (i == 204){
    //    Logger(DEBUG) << "rhoSol = " << rhoSol << ", vSol_rot = [" << vSol[0] << ", " << vSol[1] << "], PSol = " << PSol;
    //}

    /// Compute fluxes projected onto Aij
    Fij[0] = Aij[0]*rhoSol*vSol[0] + Aij[1]*rhoSol*vSol[1] + Aij[2]*rhoSol*vSol[2]; // mass flux

    // velocity solution in the lab frame
    double vLab[DIM];
    vLab[0] = vSol[0] + vFrame[0];
    vLab[1] = vSol[1] + vFrame[1];
    vLab[2] = vSol[2] + vFrame[2];

    Fij[2] = Aij[0]*(rhoSol*vLab[0]*vSol[0]+PSol) + Aij[1]*rhoSol*vLab[0]*vSol[1] + Aij[2]*rhoSol*vLab[0]*vSol[2]; // vx flux
    Fij[3] = Aij[0]*rhoSol*vLab[1]*vSol[0] + Aij[1]*(rhoSol*vLab[1]*vSol[1]+PSol) + Aij[2]*rhoSol*vLab[1]*vSol[2]; // vy flux
    Fij[4] = Aij[0]*rhoSol*vLab[2]*vSol[0] + Aij[1]*rhoSol*vLab[2]*vSol[1] + Aij[2]*(rhoSol*vLab[2]*vSol[2]+PSol); // vz flux

    Fij[1] = Aij[0]*(vSol[0]*(PSol/(gamma-1.)+rhoSol*.5*Helper::dotProduct(vLab, vLab)) + PSol*vLab[0])
             + Aij[1]*(vSol[1]*(PSol/(gamma-1.)+rhoSol*.5*Helper::dotProduct(vLab, vLab)) + PSol*vLab[1])
             + Aij[2]*(vSol[2]*(PSol/(gamma-1.)+rhoSol*.5*Helper::dotProduct(vLab, vLab)) + PSol*vLab[2]); // energy flux

    //double vFluxBuf[DIM] = { Fij[2], Fij[3] };
    //Fij[2] = LambdaInv[0]*vFluxBuf[0]+LambdaInv[1]*vFluxBuf[1];
    //Fij[3] = LambdaInv[2]*vFluxBuf[0]+LambdaInv[3]*vFluxBuf[1];

}

