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

void HLLC::solveHLLC(const double *WL, const double *WR, const double *n_unit, const double *vij, const double &gamma, double *Fij){

    // // Do I need this?
    // /* Handle pure vacuum */
    // if (!WL[0] && !WR[0]) {
    //   Fij[0] = 0.0f;
    //   Fij[1] = 0.0f;
    //   Fij[2] = 0.0f;
    //   Fij[3] = 0.0f;
    //   Fij[4] = 0.0f;
    //   return;
    // }

    // /* STEP 0: obtain velocity in interface frame */
    // const double uL = WL[1] * n_unit[0] + WL[2] * n_unit[1] + WL[3] * n_unit[2];
    // const double uR = WR[1] * n_unit[0] + WR[2] * n_unit[1] + WR[3] * n_unit[2];
    // const double rhoLinv = (WL[0] > 0.0f) ? 1.0f / WL[0] : 0.0f;
    // const double rhoRinv = (WR[0] > 0.0f) ? 1.0f / WR[0] : 0.0f;
    // const double aL = sqrtf(hydro_gamma * WL[4] * rhoLinv);
    // const double aR = sqrtf(hydro_gamma * WR[4] * rhoRinv);
    const double hydro_gamma_plus_one = gamma + 1.0d;
    const double hydro_one_over_gamma = 1.0d / gamma;
    const double hydro_one_over_gamma_minus_one = hydro_one_over_gamma - 1.0d;

    const double uL = WL[0];
    const double uR = WR[0];
    const double rhoLinv = (WL[0] > 0.0d) ? 1.0d / WL[0] : 0.0d;
    const double rhoRinv = (WR[0] > 0.0d) ? 1.0d / WR[0] : 0.0d;

    // For wave speed estimates:
#if DIM == 2
    const double aL = sqrtf(gamma * WR[3] * rhoLinv);
    const double aR = sqrtf(gamma * WR[3] * rhoLinv);
#else
    const double aL = sqrtf(gamma * WR[4] * rhoLinv);
    const double aR = sqrtf(gamma * WR[4] * rhoLinv);
#endif

    // // Do I need this? If so, params do not match. Is it smart to call it here?
    // /* Handle vacuum: vacuum does not require iteration and is always exact */
    // if (Helper::riemann_is_vacuum(WL, WR, uL, uR, aL, aR)) {
    //   Helper::riemann_solve_vacuum_flux(WL, WR, uL, uR, aL, aR, n, vij, Fij);
    //   return;
    // }


    // /* Handle vacuum: vacuum does not require iteration and is always exact */
    // if (riemann_is_vacuum(WL, WR, uL, uR, aL, aR)) {
    //     riemann_solve_vacuum_flux(WL, WR, uL, uR, aL, aR, n_unit, vij, Fij);
    //     return;
    // }

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
    double qL = 1.0d;

#if DIM == 3
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
    if (Sstar >= 0.0d) {
        const double rhoLuL = WL[0] * uL;
        const double v2 = WL[1] * WL[1] + WL[2] * WL[2] + WL[3] * WL[3];
        const double eL =
        WL[4] * rhoLinv * hydro_one_over_gamma_minus_one + 0.5d * v2;
        const double SL = SLmuL + uL;

        /* flux FL */
        Fij[0] = rhoLuL;
        /* these are the actual correct fluxes in the boosted lab frame
        (not rotated to interface frame) */
#if DIM == 3
        Fij[1] = rhoLuL * WL[1] + WL[4] * n_unit[0];
        Fij[2] = rhoLuL * WL[2] + WL[4] * n_unit[1];
        Fij[3] = rhoLuL * WL[3] + WL[4] * n_unit[2];
        Fij[4] = rhoLuL * eL + WL[4] * uL;
#else
        Fij[1] = rhoLuL * WL[1] + WL[3] * n_unit[0];
        Fij[2] = rhoLuL * WL[2] + WL[3] * n_unit[1];
        Fij[3] = rhoLuL * eL + WL[3] * uL;
#endif

        if (SL < 0.0d) {

            const double starfac = SLmuL / (SL - Sstar);
            const double rhoLSL = WL[0] * SL;
            const double SstarmuL = Sstar - uL;
            const double rhoLSLstarfac = rhoLSL * (starfac - 1.0d);
            const double rhoLSLSstarmuL = rhoLSL * SstarmuL * starfac;

            Fij[0] += rhoLSLstarfac;
            Fij[1] += rhoLSLstarfac * WL[1] + rhoLSLSstarmuL * n_unit[0];
            Fij[2] += rhoLSLstarfac * WL[2] + rhoLSLSstarmuL * n_unit[1];
            Fij[3] += rhoLSLstarfac * WL[3] + rhoLSLSstarmuL * n_unit[2];
            Fij[4] += rhoLSLstarfac * eL +
            rhoLSLSstarmuL * (Sstar + WL[4] / (WL[0] * SLmuL));
        }
    } else {
        const double rhoRuR = WR[0] * uR;
        const double v2 = WR[1] * WR[1] + WR[2] * WR[2] + WR[3] * WR[3];
        const double eR =
#if DIM == 3
        WR[4] * rhoRinv * hydro_one_over_gamma_minus_one + 0.5d * v2;
#else
        WR[3] * rhoRinv * hydro_one_over_gamma_minus_one + 0.5d * v2;
#endif
        const double SR = SRmuR + uR;


        /* flux FR */
        Fij[0] = rhoRuR;
#if DIM == 3
        Fij[1] = rhoRuR * WR[1] + WR[4] * n_unit[0];
        Fij[2] = rhoRuR * WR[2] + WR[4] * n_unit[1];
        Fij[3] = rhoRuR * WR[3] + WR[4] * n_unit[2];
        Fij[4] = rhoRuR * eR + WR[4] * uR;
#else
        Fij[1] = rhoRuR * WR[1] + WR[3] * n_unit[0];
        Fij[2] = rhoRuR * WR[2] + WR[3] * n_unit[1];
        Fij[3] = rhoRuR * eR + WR[3] * uR;
#endif

        if (SR > 0.0d) {

            const double starfac = SRmuR / (SR - Sstar);
            const double rhoRSR = WR[0] * SR;
            const double SstarmuR = Sstar - uR;
            const double rhoRSRstarfac = rhoRSR * (starfac - 1.d);
            const double rhoRSRSstarmuR = rhoRSR * SstarmuR * starfac;

            Fij[0] += rhoRSRstarfac;
#if DIM == 3
            Fij[1] += rhoRSRstarfac * WR[1] + rhoRSRSstarmuR * n_unit[0];
            Fij[2] += rhoRSRstarfac * WR[2] + rhoRSRSstarmuR * n_unit[1];
            Fij[3] += rhoRSRstarfac * WR[3] + rhoRSRSstarmuR * n_unit[2];
            Fij[4] += rhoRSRstarfac * eR +
            rhoRSRSstarmuR * (Sstar + WR[4] / (WR[0] * SRmuR));
#else
            Fij[1] += rhoRSRstarfac * WR[1] + rhoRSRSstarmuR * n_unit[0];
            Fij[2] += rhoRSRstarfac * WR[2] + rhoRSRSstarmuR * n_unit[1];
            Fij[3] += rhoRSRstarfac * eR +
            rhoRSRSstarmuR * (Sstar + WR[3] / (WR[0] * SRmuR));
#endif
        }
    }

    /* deboost to lab frame we add the flux contribution due to the
    movement of the interface the density flux is unchanged
    we add the extra velocity flux due to the absolute motion of the fluid
    similarly, we need to add the energy fluxes due to the absolute motion */
#if DIM == 3
    const double v2 = vij[0] * vij[0] + vij[1] * vij[1] + vij[2] * vij[2];
#else
    const double v2 = vij[0] * vij[0] + vij[1] * vij[1];
#endif

    /* order is important: we first use the momentum fluxes to update the energy
    flux and then de-boost the momentum fluxes! */
#if DIM == 3
    Fij[4] += vij[0] * Fij[1] + vij[1] * Fij[2] +
    vij[2] * Fij[3] + 0.5d * v2 * Fij[0];
    Fij[1] += vij[0] * Fij[0];
    Fij[2] += vij[1] * Fij[0];
    Fij[3] += vij[2] * Fij[0];
#else
    Fij[3] += vij[0] * Fij[1] + vij[1] * Fij[2] + 0.5d * v2 * Fij[0];
    Fij[1] += vij[0] * Fij[0];
    Fij[2] += vij[1] * Fij[0];
#endif

// #ifdef SWIFT_DEBUG_CHECKS
//     riemann_check_output(WL, WR, n, vij, Fij);
// #endif
}
