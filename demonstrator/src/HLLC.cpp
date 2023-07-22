//
// Created by Jakob Sappler on 09.05.23
//

#include "../include/HLLC.h"


// void HLLC::solveHLLC0(double *WR, double *WL, double *n,
//     double *totflux, const double *vij, const double &hydro_gamma){


//     const double hydro_gamma_plus_one = hydro_gamma + 1.0d;
//     const double hydro_one_over_gamma = 1.0d / hydro_gamma;
//     const double hydro_one_over_gamma_minus_one = hydro_one_over_gamma - 1.0d;

//     // Step 0: velocity in interface frame


// #if DIM == 3
//     const double uL = WL[2] * n[0] + WL[3] * n[1] + WL[4] * n[2];
//     const double uR = WR[3] * n[0] + WR[3] * n[1] + WR[4] * n[2];
// #else
//     const double uL = WL[2] * n[0] + WL[3] * n[1];
//     const double uR = WR[2] * n[0] + WR[3] * n[1];
// #endif


//     // const double uL = WL[2];
//     // const double uR = WR[2];

//     // Logger(DEBUG) << "WL[2]: " << WL[2] << " WL[3]: " << WL[3];

//     const double rhoLinv = (WL[0] > 0.0d) ? 1.0d / WL[0] : 0.0d;
//     const double rhoRinv = (WR[0] > 0.0d) ? 1.0d / WR[0] : 0.0d;

//     const double aL = sqrtf(hydro_gamma * WL[1] * rhoLinv);
//     const double aR = sqrtf(hydro_gamma * WR[1] * rhoRinv);


//     // Logger(DEBUG) << "aL^2: " << hydro_gamma * WL[1] * rhoLinv;
//     // Logger(DEBUG) << "aL: " << sqrtf(hydro_gamma * WL[1] * rhoLinv);


//     /* STEP 1: pressure estimate */
//     const double rhobar = WL[0] + WR[0];
//     const double abar = aL + aR;
//     const double pPVRS = 0.5d * ((WL[1] + WR[1]) - 0.25d * (uR - uL) * rhobar * abar);

//     const double pstar = std::max(0.0d, pPVRS);

//     /* STEP 2: wave speed estimates
//        all these speeds are along the interface normal, since uL and uR are */
// #if DIM == 3
//     double qL = 1.0d;
//     if (pstar > WL[4] && WL[4] > 0.0d) {
//       qL = sqrtf(1.0d + 0.5d * hydro_gamma_plus_one * hydro_one_over_gamma *
//                             (pstar / WL[4] - 1.0d));
//     }
//     double qR = 1.0d;
//     if (pstar > WR[4] && WR[4] > 0.0d) {
//       qR = sqrtf(1.0d + 0.5d * hydro_gamma_plus_one * hydro_one_over_gamma *
//                             (pstar / WR[4] - 1.0d));
//     }
//     const double SLmuL = -aL * qL;
//     const double SRmuR = aR * qR;
//     const double Sstar =
//         (WR[4] - WL[4] + WL[0] * uL * SLmuL - WR[0] * uR * SRmuR) /
//         (WL[0] * SLmuL - WR[0] * SRmuR);
// #else
//     double qL = 1.0d;
//     if (pstar > WL[1] && WL[1] > 0.0d) {
//       qL = sqrtf(1.0d + 0.5d * hydro_gamma_plus_one * hydro_one_over_gamma *
//                             (pstar / WL[1] - 1.0d));
//     }
//     double qR = 1.0d;
//     if (pstar > WR[1] && WR[1] > 0.0d) {
//       qR = sqrtf(1.0d + 0.5d * hydro_gamma_plus_one * hydro_one_over_gamma *
//                             (pstar / WR[1] - 1.0d));
//     }
//     const double SLmuL = -aL * qL;
//     const double SRmuR = aR * qR;
//     const double Sstar =
//         (WR[1] - WL[1] + WL[0] * uL * SLmuL - WR[0] * uR * SRmuR) /
//         (WL[0] * SLmuL - WR[0] * SRmuR);
// #endif
//     // Logger(DEBUG) <<  "SStar: " << Sstar;
//     // Logger(DEBUG) << "SLmuL: " << SLmuL;
//     // Logger(DEBUG) << "aR: " << aR << "SRmuR: " << SRmuR;

//     /* STEP 3: HLLC flux in a frame moving with the interface velocity */
// #if DIM == 3
//     if (Sstar >= 0.0d) {
//         const double rhoLuL = WL[0] * uL;
//         const double v2 = WL[1] * WL[1] + WL[2] * WL[2] + WL[3] * WL[3];
//         //const double v2 = uL * uL;
//         const double eL =
//             WL[4] * rhoLinv * hydro_one_over_gamma_minus_one + 0.5d * v2;
//         const double SL = SLmuL + uL;

//         /* flux FL */
//         totflux[0] = rhoLuL;
//         /* these are the actual correct fluxes in the boosted lab frame
//            (not rotated to interface frame) */
//         totflux[1] = rhoLuL * WL[1] + WL[4] * n[0];
//         totflux[2] = rhoLuL * WL[2] + WL[4] * n[1];
//         totflux[3] = rhoLuL * WL[3] + WL[4] * n[2];
//         totflux[4] = rhoLuL * eL + WL[4] * uL;

//         if (SL < 0.0d) {

//             const double starfac = SLmuL / (SL - Sstar);
//             const double rhoLSL = WL[0] * SL;
//             const double SstarmuL = Sstar - uL;
//             const double rhoLSLstarfac = rhoLSL * (starfac - 1.0d);
//             const double rhoLSLSstarmuL = rhoLSL * SstarmuL * starfac;

//             totflux[0] += rhoLSLstarfac;
//             totflux[1] += rhoLSLstarfac * WL[1] + rhoLSLSstarmuL * n[0];
//             totflux[2] += rhoLSLstarfac * WL[2] + rhoLSLSstarmuL * n[1];
//             totflux[3] += rhoLSLstarfac * WL[3] + rhoLSLSstarmuL * n[2];
//             totflux[4] += rhoLSLstarfac * eL +
//                         rhoLSLSstarmuL * (Sstar + WL[4] / (WL[0] * SLmuL));
//         }
//       } else {
//         const double rhoRuR = WR[0] * uR;
//         const double v2 = WR[1] * WR[1] + WR[2] * WR[2] + WR[3] * WR[3];
//         const double eR =
//             WR[4] * rhoRinv * hydro_one_over_gamma_minus_one + 0.5d * v2;
//         const double SR = SRmuR + uR;

//         /* flux FR */
//         totflux[0] = rhoRuR;
//         totflux[1] = rhoRuR * WR[1] + WR[4] * n[0];
//         totflux[2] = rhoRuR * WR[2] + WR[4] * n[1];
//         totflux[3] = rhoRuR * WR[3] + WR[4] * n[2];
//         totflux[4] = rhoRuR * eR + WR[4] * uR;

//         if (SR > 0.0d) {

//             const double starfac = SRmuR / (SR - Sstar);
//             const double rhoRSR = WR[0] * SR;
//             const double SstarmuR = Sstar - uR;
//             const double rhoRSRstarfac = rhoRSR * (starfac - 1.f);
//             const double rhoRSRSstarmuR = rhoRSR * SstarmuR * starfac;

//             totflux[0] += rhoRSRstarfac;
//             totflux[1] += rhoRSRstarfac * WR[1] + rhoRSRSstarmuR * n[0];
//             totflux[2] += rhoRSRstarfac * WR[2] + rhoRSRSstarmuR * n[1];
//             totflux[3] += rhoRSRstarfac * WR[3] + rhoRSRSstarmuR * n[2];
//             totflux[4] += rhoRSRstarfac * eR +
//                         rhoRSRSstarmuR * (Sstar + WR[4] / (WR[0] * SRmuR));
//         }
//     }
// #else
//     if (Sstar >= 0.0d) {
//         // Logger(DEBUG) << "SStar >= 0";
//         // Logger(DEBUG) << " a ";
//         const double rhoLuL = WL[0] * uL;
//         const double v2 = WL[2] * WL[2] + WL[3] * WL[3];
//         //const double v2 = uL * uL;
//         const double eL =
//             WL[1] * rhoLinv * hydro_one_over_gamma_minus_one + 0.5d * v2;
//         const double SL = SLmuL + uL;

//         /* flux FL */
//         // WATCH THIS
//         totflux[0] = rhoLuL;
//         /* these are the actual correct fluxes in the boosted lab frame
//            (not rotated to interface frame) */
//         totflux[1] = rhoLuL * eL + WL[1] * uL;
//         totflux[2] = rhoLuL * WL[2] + WL[1] * n[0];
//         totflux[3] = rhoLuL * WL[3] + WL[1] * n[1];

//         // From Toro:
//         // totflux[0] = rhoLuL;
//         // totflux[1] = rhoLuL * eL + WL[1] * uL;
//         // totflux[2] = rhoLuL * WL[2];
//         // totflux[3] = rhoLuL * WL[3];

//         if (SL < 0.0d) {
//             // Logger(DEBUG) << "SL < 0";
//             // Logger(DEBUG) << " b ";
//             const double starfac = SLmuL / (SL - Sstar);
//             const double rhoLSL = WL[0] * SL;
//             const double SstarmuL = Sstar - uL;
//             const double rhoLSLstarfac = rhoLSL * (starfac - 1.0d);
//             const double rhoLSLSstarmuL = rhoLSL * SstarmuL * starfac;

//             totflux[0] += rhoLSLstarfac;
//             totflux[1] += rhoLSLstarfac * eL +
//                 rhoLSLSstarmuL * (Sstar + WL[1] / (WL[0] * SLmuL));
//             totflux[2] += rhoLSLstarfac * WL[2] + rhoLSLSstarmuL * n[0];
//             totflux[3] += rhoLSLstarfac * WL[3] + rhoLSLSstarmuL * n[1];

//             // From Toro:
//             // totflux[0] += rhoLSLstarfac;
//             // totflux[1] += SstarmuL / (SL - Sstar) * SL * eL +
//             //     rhoLSL * starfac * SstarmuL * (Sstar  + WL[1] / (WL[0] * SLmuL));
//             // totflux[2] += WL[0] * starfac * Sstar;
//             // totflux[3] += WL[0] * starfac * WL[3];
//         }
//     } else {
//         // Logger(DEBUG) << "SStar < 0";
//         // Logger(DEBUG) << " c ";
//         const double rhoRuR = WR[0] * uR;
//         const double v2 = WR[2] * WR[2] + WR[3] * WR[3];
//         const double eR =
//             WR[1] * rhoRinv * hydro_one_over_gamma_minus_one + 0.5d * v2;
//         const double SR = SRmuR + uR;

//         /* flux FR */
//         totflux[0] = rhoRuR;
//         totflux[1] = rhoRuR * eR + WR[1] * uR;
//         totflux[2] = rhoRuR * WR[2] + WR[1] * n[0];
//         totflux[3] = rhoRuR * WR[3] + WR[1] * n[1];

//         // totflux[0] = rhoRuR;
//         // totflux[1] = rhoRuR * eR + WR[1] * uR;
//         // totflux[2] = rhoRuR * WR[2];
//         // totflux[3] = rhoRuR * WR[3];

//         if (SR > 0.0d) {
//             // Logger(INFO) << "SR > 0";
//             // Logger(DEBUG) << " d ";
//             const double starfac = SRmuR / (SR - Sstar);
//             const double rhoRSR = WR[0] * SR;
//             const double SstarmuR = Sstar - uR;
//             const double rhoRSRstarfac = rhoRSR * (starfac - 1.d);
//             const double rhoRSRSstarmuR = rhoRSR * SstarmuR * starfac;
//             //Logger(DEBUG) << Sstar;
//             totflux[0] += rhoRSRstarfac;
//             totflux[1] += rhoRSRstarfac * eR +
//                 rhoRSRSstarmuR * (Sstar + WR[1] / (WR[0] * SRmuR));
//             totflux[2] += rhoRSRstarfac * WR[2] + rhoRSRSstarmuR * n[0];
//             totflux[3] += rhoRSRstarfac * WR[3] + rhoRSRSstarmuR * n[1];

//             // // From Toro:
//             // totflux[0] += rhoRSRstarfac;
//             // totflux[1] += SstarmuR / (SR - Sstar) * SR * eR +
//             //     rhoRSR * starfac * SstarmuR * (Sstar  + WR[1] / (WR[0] * SRmuR));
//             // totflux[2] += WR[0] * starfac * Sstar;
//             // totflux[3] += WR[0] * starfac * WR[3];
//         }

//     }



// #endif // IF DIM == 3

//     //TODO: Rotate and project here!

//     // Rotating:
// //
// //     double LambdaInv[DIM*DIM];
// // #if DIM == 3
// //     Helper::rotationMatrix3D(unitX, hatAij, LambdaInv);
// //
// //     double vSolBuf[DIM] = { vSol[0], vSol[1] , vSol[2]};
// //
// //     vSol[0] = LambdaInv[0]*vSolBuf[0]+LambdaInv[1]*vSolBuf[1];
// //     vSol[1] = LambdaInv[2]*vSolBuf[0]+LambdaInv[3]*vSolBuf[1];

//     /* deboost to lab frame we add the flux contribution due to the
//     movement of the interface the density flux is unchanged
//     we add the extra velocity flux due to the absolute motion of the fluid
//     similarly, we need to add the energy fluxes due to the absolute motion */
// #if DIM == 3
//     const double v2 = vij[0] * vij[0] + vij[1] * vij[1] + vij[2] * vij[2];

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
//     totflux[1] += vij[0] * totflux[2] + vij[1] * totflux[3] +
//                   + 0.5f * v2 * totflux[0];
//     totflux[2] += vij[0] * totflux[0];
//     totflux[3] += vij[1] * totflux[0];
// #endif // If DIM == 3



// }

// Riemann HLLC directly and exactly from swift.
// Sole exeption: Structure of the state vector differs,
// with pressure and momentum components exchange. Flux changes accordingly

void HLLC::solveHLLC1(double *WR, double *WL, double *n,
    double *totflux, const double *vij, const double &hydro_gamma){


    const double hydro_gamma_plus_one = hydro_gamma + 1.0d;
    const double hydro_one_over_gamma = 1.0d / hydro_gamma;
    const double hydro_one_over_gamma_minus_one = hydro_one_over_gamma - 1.0d;

    // Step 0: velocity in interface frame


#if DIM == 3
    const double uL = WL[2] * n[0] + WL[3] * n[1] + WL[4] * n[2];
    const double uR = WR[3] * n[0] + WR[3] * n[1] + WR[4] * n[2];
#else
    const double uL = WL[2] * n[0] + WL[3] * n[1];
    const double uR = WR[2] * n[0] + WR[3] * n[1];
#endif

    const double rhoLinv = (WL[0] > 0.0d) ? 1.0d / WL[0] : 0.0d;
    const double rhoRinv = (WR[0] > 0.0d) ? 1.0d / WR[0] : 0.0d;

    const double aL = sqrtf(hydro_gamma * WL[1] * rhoLinv);
    const double aR = sqrtf(hydro_gamma * WR[1] * rhoRinv);

    /* STEP 1: pressure estimate */
    const double rhobar = WL[0] + WR[0];
    const double abar = aL + aR;
    const double pPVRS = 0.5d * ((WL[1] + WR[1]) - 0.25d * (uR - uL) * rhobar * abar);

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
    if (pstar > WL[1] && WL[1] > 0.0d) {
      qL = sqrtf(1.0d + 0.5d * hydro_gamma_plus_one * hydro_one_over_gamma *
                            (pstar / WL[1] - 1.0d));
    }
    double qR = 1.0d;
    if (pstar > WR[1] && WR[1] > 0.0d) {
      qR = sqrtf(1.0d + 0.5d * hydro_gamma_plus_one * hydro_one_over_gamma *
                            (pstar / WR[1] - 1.0d));
    }
    const double SLmuL = -aL * qL;
    const double SRmuR = aR * qR;
    const double Sstar =
        (WR[1] - WL[1] + WL[0] * uL * SLmuL - WR[0] * uR * SRmuR) /
        (WL[0] * SLmuL - WR[0] * SRmuR);
#endif
    // Logger(DEBUG) <<  "SStar: " << Sstar;
    // Logger(DEBUG) << "SLmuL: " << SLmuL;
    // Logger(DEBUG) << "aR: " << aR << "SRmuR: " << SRmuR;

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
        // Logger(DEBUG) << "SStar >= 0";
        // Logger(DEBUG) << " a ";
        const double rhoLuL = WL[0] * uL;
        const double v2 = WL[2] * WL[2] + WL[3] * WL[3];
        //const double v2 = uL * uL;
        const double eL =
            WL[1] * rhoLinv * hydro_one_over_gamma_minus_one + 0.5d * v2;
        const double SL = SLmuL + uL;

        /* flux FL */
        // WATCH THIS
        totflux[0] = rhoLuL;
        /* these are the actual correct fluxes in the boosted lab frame
           (not rotated to interface frame) */
        totflux[1] = rhoLuL * eL + WL[1] * uL;
        totflux[2] = rhoLuL * WL[2] + WL[1] * n[0];
        totflux[3] = rhoLuL * WL[3] + WL[1] * n[1];


        if (SL < 0.0d) {
            // Logger(DEBUG) << "SL < 0";
            // Logger(DEBUG) << " b ";
            const double starfac = SLmuL / (SL - Sstar);
            const double rhoLSL = WL[0] * SL;
            const double SstarmuL = Sstar - uL;
            const double rhoLSLstarfac = rhoLSL * (starfac - 1.0d);
            const double rhoLSLSstarmuL = rhoLSL * SstarmuL * starfac;

            totflux[0] += rhoLSLstarfac;
            totflux[1] += rhoLSLstarfac * eL +
                rhoLSLSstarmuL * (Sstar + WL[1] / (WL[0] * SLmuL));
            totflux[2] += rhoLSLstarfac * WL[2] + rhoLSLSstarmuL * n[0];
            totflux[3] += rhoLSLstarfac * WL[3] + rhoLSLSstarmuL * n[1];

        }
    } else {
        // Logger(DEBUG) << "SStar < 0";
        // Logger(DEBUG) << " c ";
        const double rhoRuR = WR[0] * uR;
        const double v2 = WR[2] * WR[2] + WR[3] * WR[3];
        const double eR =
            WR[1] * rhoRinv * hydro_one_over_gamma_minus_one + 0.5d * v2;
        const double SR = SRmuR + uR;

        /* flux FR */
        totflux[0] = rhoRuR;
        totflux[1] = rhoRuR * eR + WR[1] * uR;
        totflux[2] = rhoRuR * WR[2] + WR[1] * n[0];
        totflux[3] = rhoRuR * WR[3] + WR[1] * n[1];


        if (SR > 0.0d) {
            // Logger(INFO) << "SR > 0";
            // Logger(DEBUG) << " d ";
            const double starfac = SRmuR / (SR - Sstar);
            const double rhoRSR = WR[0] * SR;
            const double SstarmuR = Sstar - uR;
            const double rhoRSRstarfac = rhoRSR * (starfac - 1.d);
            const double rhoRSRSstarmuR = rhoRSR * SstarmuR * starfac;
            //Logger(DEBUG) << Sstar;
            totflux[0] += rhoRSRstarfac;
            totflux[1] += rhoRSRstarfac * eR +
                rhoRSRSstarmuR * (Sstar + WR[1] / (WR[0] * SRmuR));
            totflux[2] += rhoRSRstarfac * WR[2] + rhoRSRSstarmuR * n[0];
            totflux[3] += rhoRSRstarfac * WR[3] + rhoRSRSstarmuR * n[1];

        }

    }



#endif // DIM == 3

#if DIM == 3
    const double v2 = vij[0] * vij[0] + vij[1] * vij[1] + vij[2] * vij[2];

    /* order is important: we first use the momentum fluxes to update the energy
       flux and then de-boost the momentum fluxes! */
    totflux[4] += vij[0] * totflux[1] + vij[1] * totflux[2] +
                  vij[2] * totflux[3] + 0.5f * v2 * totflux[0];
    totflux[1] += vij[0] * totflux[0];
    totflux[2] += vij[1] * totflux[0];
    totflux[3] += vij[2] * totflux[0];

#else

    const double v2 = vij[0] * vij[0] + vij[1] * vij[1];
    /* order is important: we first use the momentum fluxes to update the energy
       flux and then de-boost the momentum fluxes! */
    totflux[1] += vij[0] * totflux[2] + vij[1] * totflux[3] +
                  + 0.5f * v2 * totflux[0];
    totflux[2] += vij[0] * totflux[0];
    totflux[3] += vij[1] * totflux[0];

#endif // If DIM == 3



}




// HLL Solver
void HLLC::HLL(double *WL, double *WR, double *totflux, const double &hydro_gamma){

    // Step 1: Sound speed:

    const double aL = sqrtf(hydro_gamma * WL[1] / WL[0]);
    const double aR = sqrtf(hydro_gamma * WR[1] / WR[0]);

    // Step 2: Wave speed estimate:

    double v2L = pow(WL[2], 2) + pow(WL[3], 2);
    double v2R = pow(WR[2], 2) + pow(WR[3], 2);

    double eL = WL[1] / WL[0] * (1 / hydro_gamma - 1.d)
        + v2L;
    double eR = WR[1] / WR[0] * (1 / hydro_gamma - 1.d)
        + v2R;

#if USE_ROE

    // Enthalpy:
    double HL = (eL - WL[1]) / WL[0];
    double HR = (eR - WR[1]) / WR[0];

    // approximate enthalpy
    double HBar = (sqrtf(WL[0]) * HL + sqrtf(WR[0]) * WR[0])
        / (sqrtf(WL[0]) + sqrtf(WR[0]));

    // Approximate particle speed:
    double uBar = (sqrtf(WL[0]) * WL[2] + sqrtf(WR[0]) * WR[2])
        / (sqrtf(WL[0]) + sqrtf(WR[0]));

    // Approximate sound speed:
    double aBar = sqrtf((hydro_gamma - 1.d) * (HBar - 5.d * pow(uBar, 2)));

    // Wave speed estimate:
    double SL = uBar - aBar;
    double SR = uBar + aBar;

#else // USE_ROE

    double SL = std::min(WL[2] - aL, WR[2] - aR);
    double SR = std::max(WL[2] - aL, WR[2] - aR);

#endif // USE_ROE

    // Step 3: Flux estimates:


    //TODO: Add 3d


    if (SL <= 0){
        if (SR >= 0){
            // F_hll

            totflux[0] = (SR * WL[0] * WL[2] - SL * WR[0] * WR[2]
                + SL * SR * (WR[0] - WL[0]))
                 / (SR - SL);

            totflux[1] = (SR * WL[2] * (eL - WL[1])
                - SL * WR[2] * (eR - WR[1])
                + SL * SR * (eR - eL)) / (SR - SL);

            totflux[2] = (SR * (WL[0] * pow(WL[2], 2) + WL[1])
                - SL * (WR[0] * pow(WR[2], 2) + WR[1])
                + SL * SR * (WR[0] * WR[2] - WL[0] * WL[2]))
                    / (SR - SL);

            totflux[2] = (SR * WL[0] * WL[2] * WL[3]
                - SL * WR[0] * WR[2] * WR[3]
                + SL * SR * (WR[0] * WR[3] - WL[0] * WL[3]))
                    / (SR - SL);

            // TODO: Add 3d
        }
        else{
            // F_R

            totflux[0] = WR[0] * WR[2];
            totflux[1] = WR[2] * (eR + WR[1]);
            totflux[2] = WR[0] * pow(WR[2], 2) + WR[1];
            totflux[3] = WR[0] * WR[2] * WR[3];
        }
    }
    else{
        // F_L

        totflux[0] = WL[0] * WL[2];
        totflux[1] = WL[2] * (eL + WL[1]);
        totflux[2] = WL[0] * pow(WL[2], 2) + WL[1];
        totflux[3] = WL[0] * WL[2] * WL[3];
    }
}
