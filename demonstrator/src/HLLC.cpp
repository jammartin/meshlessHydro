//
// Created by Jakob Sappler on 09.05.23
//

#include "../include/HLLC.h"


//HLLC::HLLC()


// Method solve: eqivalent to riemann_solve_for_flux in SWIFT, BUT:
// Takes velocities in the frame of reference of the interface (already boosted)
// Returns the de-boosted fluxes that are not yet rotated to the interface frame
// Also, not yet projected onto the faces

// @param WL the left state vector
// @param WR the right state vector
// @param n_unit normal vector of the interface
// @param vij velocity vector of the interface
// @param Fij flux vector to store the solution in

// Difference from SWIFT: WL/WR already in interface frame: Boosted and rotated
void HLLC::solveHLLC(double *WL, double *WR, double *n,
    double *totflux, const double *vij, const double &hydro_gamma){


    const double hydro_gamma_plus_one = hydro_gamma + 1.0d;
    const double hydro_one_over_gamma = 1.0d / hydro_gamma;
    const double hydro_one_over_gamma_minus_one = hydro_one_over_gamma - 1.0d;

    // Step 0: velocity in interface frame
#if DIM == 3
    const double uL = WL[1] * n[0] + WL[2] * n[1] + WL[3] * n[2];
    const double uR = WR[1] * n[0] + WR[2] * n[1] + WR[3] * n[2];
#else
    const double uL = WL[1] * n[0] + WL[2] * n[1];
    const double uR = WR[1] * n[0] + WR[2] * n[1];
#endif

    const double rhoLinv = (WL[0] > 0.0d) ? 1.0d / WL[0] : 0.0d;
    const double rhoRinv = (WR[0] > 0.0d) ? 1.0d / WR[0] : 0.0d;

#if DIM == 3
    const double aL = sqrtf(hydro_gamma * WL[4] * rhoLinv);
    const double aR = sqrtf(hydro_gamma * WR[4] * rhoRinv);
#else
    const double aL = sqrtf(hydro_gamma * WL[3] * rhoLinv);
    const double aR = sqrtf(hydro_gamma * WR[3] * rhoRinv);
#endif

    /* STEP 1: pressure estimate */
    const double rhobar = WL[0] + WR[0];
    const double abar = aL + aR;
    const double pPVRS =
#if DIM == 3
        0.5d * ((WL[4] + WR[4]) - 0.25d * (uR - uL) * rhobar * abar);
#else
        0.5d * ((WL[3] + WR[3]) - 0.25d * (uR - uL) * rhobar * abar);
#endif
    const double pstar = std::max(0.0d, pPVRS);

    /* STEP 2: wave speed estimates
       all these speeds are along the interface normal, since uL and uR are */
#if DIM == 3
    double qL = 1.0d;
    if (pstar > WL[4] && WL[4] > 0.0d) {
      qL = sqrtf(1.0d + 0.5d * hydro_gamma_plus_one * hydro_one_over_gamma *
                            (pstar / WL[4] - 1.0d));
    }
    double qR = 1.0d;
    if (pstar > WR[4] && WR[4] > 0.0d) {
      qR = sqrtf(1.0d + 0.5d * hydro_gamma_plus_one * hydro_one_over_gamma *
                            (pstar / WR[4] - 1.0d));
    }
    const double SLmuL = -aL * qL;
    const double SRmuR = aR * qR;
    const double Sstar =
        (WR[4] - WL[4] + WL[0] * uL * SLmuL - WR[0] * uR * SRmuR) /
        (WL[0] * SLmuL - WR[0] * SRmuR);
#else
    double qL = 1.0d;
    if (pstar > WL[3] && WL[3] > 0.0d) {
      qL = sqrtf(1.0d + 0.5d * hydro_gamma_plus_one * hydro_one_over_gamma *
                            (pstar / WL[3] - 1.0d));
    }
    double qR = 1.0d;
    if (pstar > WR[3] && WR[3] > 0.0d) {
      qR = sqrtf(1.0d + 0.5d * hydro_gamma_plus_one * hydro_one_over_gamma *
                            (pstar / WR[3] - 1.0d));
    }
    const double SLmuL = -aL * qL;
    const double SRmuR = aR * qR;
    const double Sstar =
        (WR[3] - WL[3] + WL[0] * uL * SLmuL - WR[0] * uR * SRmuR) /
        (WL[0] * SLmuL - WR[0] * SRmuR);
#endif

    /* STEP 3: HLLC flux in a frame moving with the interface velocity */
#if DIM == 3
    if (Sstar >= 0.0d) {
        const double rhoLuL = WL[0] * uL;
        const double v2 = WL[1] * WL[1] + WL[2] * WL[2] + WL[3] * WL[3];
        //const double v2 = uL * uL;
        const double eL =
            WL[4] * rhoLinv * hydro_one_over_gamma_minus_one + 0.5d * v2;
        const double SL = SLmuL + uL;

        /* flux FL */
        totflux[0] = rhoLuL;
        /* these are the actual correct fluxes in the boosted lab frame
           (not rotated to interface frame) */
        totflux[1] = rhoLuL * WL[1] + WL[4] * n[0];
        totflux[2] = rhoLuL * WL[2] + WL[4] * n[1];
        totflux[3] = rhoLuL * WL[3] + WL[4] * n[2];
        totflux[4] = rhoLuL * eL + WL[4] * uL;

        if (SL < 0.0d) {

          const double starfac = SLmuL / (SL - Sstar);
          const double rhoLSL = WL[0] * SL;
          const double SstarmuL = Sstar - uL;
          const double rhoLSLstarfac = rhoLSL * (starfac - 1.0d);
          const double rhoLSLSstarmuL = rhoLSL * SstarmuL * starfac;

          totflux[0] += rhoLSLstarfac;
          totflux[1] += rhoLSLstarfac * WL[1] + rhoLSLSstarmuL * n[0];
          totflux[2] += rhoLSLstarfac * WL[2] + rhoLSLSstarmuL * n[1];
          totflux[3] += rhoLSLstarfac * WL[3] + rhoLSLSstarmuL * n[2];
          totflux[4] += rhoLSLstarfac * eL +
                        rhoLSLSstarmuL * (Sstar + WL[4] / (WL[0] * SLmuL));
        }
      } else {
        const double rhoRuR = WR[0] * uR;
        const double v2 = WR[1] * WR[1] + WR[2] * WR[2] + WR[3] * WR[3];
        const double eR =
            WR[4] * rhoRinv * hydro_one_over_gamma_minus_one + 0.5d * v2;
        const double SR = SRmuR + uR;

        /* flux FR */
        totflux[0] = rhoRuR;
        totflux[1] = rhoRuR * WR[1] + WR[4] * n[0];
        totflux[2] = rhoRuR * WR[2] + WR[4] * n[1];
        totflux[3] = rhoRuR * WR[3] + WR[4] * n[2];
        totflux[4] = rhoRuR * eR + WR[4] * uR;

        if (SR > 0.0d) {

          const double starfac = SRmuR / (SR - Sstar);
          const double rhoRSR = WR[0] * SR;
          const double SstarmuR = Sstar - uR;
          const double rhoRSRstarfac = rhoRSR * (starfac - 1.f);
          const double rhoRSRSstarmuR = rhoRSR * SstarmuR * starfac;

          totflux[0] += rhoRSRstarfac;
          totflux[1] += rhoRSRstarfac * WR[1] + rhoRSRSstarmuR * n[0];
          totflux[2] += rhoRSRstarfac * WR[2] + rhoRSRSstarmuR * n[1];
          totflux[3] += rhoRSRstarfac * WR[3] + rhoRSRSstarmuR * n[2];
          totflux[4] += rhoRSRstarfac * eR +
                        rhoRSRSstarmuR * (Sstar + WR[4] / (WR[0] * SRmuR));
        }
    }
#else
    if (Sstar >= 0.0d) {
        const double rhoLuL = WL[0] * uL;
        const double v2 = WL[1] * WL[1] + WL[2] * WL[2];
        //const double v2 = uL * uL;
        const double eL =
            WL[3] * rhoLinv * hydro_one_over_gamma_minus_one + 0.5d * v2;
        const double SL = SLmuL + uL;

        /* flux FL */
        totflux[0] = rhoLuL;
        /* these are the actual correct fluxes in the boosted lab frame
           (not rotated to interface frame) */
        totflux[1] = rhoLuL * WL[1] + WL[4] * n[0];
        totflux[2] = rhoLuL * WL[2] + WL[4] * n[1];
        totflux[3] = rhoLuL * eL + WL[4] * uL;

        if (SL < 0.0d) {

          const double starfac = SLmuL / (SL - Sstar);
          const double rhoLSL = WL[0] * SL;
          const double SstarmuL = Sstar - uL;
          const double rhoLSLstarfac = rhoLSL * (starfac - 1.0d);
          const double rhoLSLSstarmuL = rhoLSL * SstarmuL * starfac;

          totflux[0] += rhoLSLstarfac;
          totflux[1] += rhoLSLstarfac * WL[1] + rhoLSLSstarmuL * n[0];
          totflux[2] += rhoLSLstarfac * WL[2] + rhoLSLSstarmuL * n[1];
          totflux[3] += rhoLSLstarfac * eL +
                        rhoLSLSstarmuL * (Sstar + WL[3] / (WL[0] * SLmuL));
        }
      } else {
        const double rhoRuR = WR[0] * uR;
        const double v2 = WR[1] * WR[1] + WR[2] * WR[2];
        const double eR =
            WR[3] * rhoRinv * hydro_one_over_gamma_minus_one + 0.5d * v2;
        const double SR = SRmuR + uR;

        /* flux FR */
        totflux[0] = rhoRuR;
        totflux[1] = rhoRuR * WR[1] + WR[4] * n[0];
        totflux[2] = rhoRuR * WR[2] + WR[4] * n[1];
        totflux[3] = rhoRuR * eR + WR[3] * uR;

        if (SR > 0.0d) {

          const double starfac = SRmuR / (SR - Sstar);
          const double rhoRSR = WR[0] * SR;
          const double SstarmuR = Sstar - uR;
          const double rhoRSRstarfac = rhoRSR * (starfac - 1.f);
          const double rhoRSRSstarmuR = rhoRSR * SstarmuR * starfac;

          totflux[0] += rhoRSRstarfac;
          totflux[1] += rhoRSRstarfac * WR[1] + rhoRSRSstarmuR * n[0];
          totflux[2] += rhoRSRstarfac * WR[2] + rhoRSRSstarmuR * n[1];
          totflux[3] += rhoRSRstarfac * eR +
                        rhoRSRSstarmuR * (Sstar + WR[3] / (WR[0] * SRmuR));
        }
    }

#endif // IF DIM == 3

    //TODO: Rotate and project here!

    // Rotating:
//
//     double LambdaInv[DIM*DIM];
// #if DIM == 3
//     Helper::rotationMatrix3D(unitX, hatAij, LambdaInv);
//
//     double vSolBuf[DIM] = { vSol[0], vSol[1] , vSol[2]};
//
//     vSol[0] = LambdaInv[0]*vSolBuf[0]+LambdaInv[1]*vSolBuf[1];
//     vSol[1] = LambdaInv[2]*vSolBuf[0]+LambdaInv[3]*vSolBuf[1];

//     /* deboost to lab frame we add the flux contribution due to the
//     movement of the interface the density flux is unchanged
//     we add the extra velocity flux due to the absolute motion of the fluid
//     similarly, we need to add the energy fluxes due to the absolute motion */
// #if DIM == 3
//     const double v2 = vij[0] * vij[0] + vij[1] * vij[1] + vij[2] * vij[2];
//
//     /* order is important: we first use the momentum fluxes to update the energy
//        flux and then de-boost the momentum fluxes! */
//     totflux[4] += vij[0] * totflux[1] + vij[1] * totflux[2] +
//                   vij[2] * totflux[3] + 0.5f * v2 * totflux[0];
//     totflux[1] += vij[0] * totflux[0];
//     totflux[2] += vij[1] * totflux[0];
//     totflux[3] += vij[2] * totflux[0];
// #else
//     const double v2 = vij[0] * vij[0] + vij[1] * vij[1];
//     /* order is important: we first use the momentum fluxes to update the energy
//        flux and then de-boost the momentum fluxes! */
//     totflux[3] += vij[0] * totflux[1] + vij[1] * totflux[2] +
//                   + 0.5f * v2 * totflux[0];
//     totflux[1] += vij[0] * totflux[0];
//     totflux[2] += vij[1] * totflux[0];
// #endif // If DIM == 3
//



}
