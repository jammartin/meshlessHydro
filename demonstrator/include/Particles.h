//
// Created by Johannes Martin on 05.09.22.
//

#ifndef MESHLESSHYDRO_PARTICLES_H
#define MESHLESSHYDRO_PARTICLES_H

#include <cmath>
#include <limits>
#include <highfive/H5File.hpp>

#include "parameter.h"
#include "Logger.h"
#include "Domain.h"
#include "Helper.h"
#include "Riemann.h"

namespace Kernel {
    double cubicSpline(const double &r, const double &h);
    // For SPH only and 2d:
    // dW(r, h)/dr:
    double dWdr(const double &r, const double &h);
    // dW(r, h)/dh:
    double dWdh(const double &r, const double &h);
}

class Particles {

public:
    Particles(int numParticles, bool ghosts=false);
    ~Particles();

    int N;
    int *matId;
    int *cell; // cell in which particle at index resides
    double *m, *u, *x, *y, *vx, *vy, *rho, *P;

    double (*rhoGrad)[DIM], (*vxGrad)[DIM], (*vyGrad)[DIM], (*vzGrad)[DIM], (*PGrad)[DIM];

    void assignParticlesAndCells(Domain &domain);
    void gridNNS(Domain &domain, const double &kernelSize);

#if RUNSPH
    // For SPH
    double *dEdt, *dn, *drho;
    // double *unimportant;
    // For Artificial Viscocity: sound speed, additional acceleration & energy terms
    double *cs, *axArtVisc, *ayArtVisc, *dudtArtVisc;
#if DIM == 3
    double *azArtVisc;
#endif

    double *ax, *ay;
    double *z, *vz;
    double *az;
    // For SPH:
    // For comparable ICs: This sets the internal energies so that P = 2.5 everywhere
    void setInternalEnergy(const double Pressure, const double gamma);

    // Computes the density via kernel smoothing w/o ghost particles
    void compDensitySPH(const double &kernelSize);

    // Calculate acceleration for each particle, w/o Ghosts
    // c.f. eq 8 in Monaghan: SPH and its diverse Applications, Annu.Rev. Fluid mechanics, 2012
    void compAccSPH(const double &kernelSize);

    // Simple euler integration
    void eulerSPH(const double &dt, const Domain &domain);

    // SPH energy calculation:
    void compuis(const double &dt, const double &kernelSize);

    // Compute omegas, as in eq. 7 in GIZMO paper
    void compOmegas(const double &kernelSize);

    void calcdndrho(const double &kernelSize);

    void calcdE(const Particles &ghostParticles, const double &kernelSize);

#if PERIODIC_BOUNDARIES
    // Computes the density via kernel smoothing w/ ghost particles
    void compDensitySPH(const Particles &ghostParticles, const double &kernelSize);

    // Calculate acceleration for each particle, w/ Ghosts
    // c.f. eq 8 in Monaghan: SPH and its diverse Applications, Annu.Rev. Fluid mechanics, 2012
    void compAccSPH(const Particles &ghostParticles, const double &kernelSize);

    // Compute omegas, as in eq. 7 in GIZMO paper
    void compOmegas(const Particles &ghostParticles, const double &kernelSize);

    void calcdndrho(const Particles &ghostParticles, const double &kernelSize);

    void calcdE(const double &kernelSize);

#endif

#if ARTVISC
    // Artificial Viscocity:
    // Implemented as in: "A smooth particle hydrodynamics code to model collisions between solid, self-gravitating objects", Schäfer et al. 2016

    // Compute speed of sound:
    // No need to update ghost particles after this!
    void compCs(const double gamma);

    // Compute Mu_ij:
    double compMuij(int i, int j, const double &kernelSize);

    // Compute PI_ij:
    double compPIij(int i, int j, const double &kernelSize);

    // Compute additional acceleration term for artificial viscocity
    void compAccArtVisc(const double &kernelSize);

    // Compute additional energy terms for artificial viscocity
    void compUiArtVisc(const double &kernelSize);

#if PERIODIC_BOUNDARIES
    // For artificial viscocity & periodic boundaries:

    // Compute Mu_ij for a ghost Particle as an interaction partner:
    double compMuij(const Particles &ghostParticles, int i, int j, const double &kernelSize);

    // Compute PI_ij for a ghost Particle as an interaction partner:
    double compPIij(const Particles &ghostParticles, int i, int j, const double &kernelSize);

    // Compute additional acceleration term for artificial viscocity w/ ghostPartices
    void compAccArtVisc(const Particles &ghostParticles, const double &kernelSize);

    // Compute additional energy terms for artificial viscocity w/ ghostPartices
    void compUiArtVisc(const Particles &ghostParticles, const double &kernelSize);
#endif
#endif
#endif // RUNSPH

    void compDensity(const double &kernelSize); // also computes omega

    void compPsijTilde(Helper &helper, const double &kernelSize);
    void gradient(double *f, double (*grad)[DIM]);
    void slopeLimiter(const double &kernelSize,
                      Particles *ghostParticles=nullptr);
    void compPressure(const double &gamma);
    void compEffectiveFace();

    double compGlobalTimestep(const double &gamma, const double &kernelSize);
    void compRiemannStatesLR(const double &dt, const double &kernelSize, const double &gamma);

    void solveRiemannProblems(const double &gamma, const Particles &ghostParticles);

    void collectFluxes(Helper &helper, const Particles &ghostParticles);

    void updateStateAndPosition(const double &dt, const Domain &domain);

// //For DEBUGGING:
//    void printDensity(const double &gamma);

    double pairwiseLimiter(double phi_0, double phi_i, double phi_j, double xijxi_abs, double xjxi_abs);


#if PERIODIC_BOUNDARIES
    void createGhostParticles(Domain &domain,
                              Particles &ghostParticles, const double &kernelSize);
    void ghostNNS(Domain &domain, const Particles &ghostParticles, const double &kernelSize);

    void compDensity(const Particles &ghostParticles, const double &kernelSize);

    void compPsijTilde(Helper &helper, const Particles &ghostParticles, const double &kernelSize);
    void gradient(double *f, double (*grad)[DIM], double *fGhost, const Particles &ghostParticles); //TODO: remove ghostParticles argument
    void compEffectiveFace(const Particles &ghostParticles);
    void compRiemannStatesLR(const double &dt, const double &kernelSize, const double &gamma,
                             const Particles &ghostParticles);

    /// functions to copy computed quantities to ghosts needed for further processing
    void updateGhostState(Particles &ghostParticles);
    //void updateGhostPsijTilde(Particles &ghostParticles);
    void updateGhostGradients(Particles &ghostParticles);

    void dumpNNL(std::string filename, const Particles &ghostParticles);

    void printNoi();
#else
    void getDomainLimits(double *domainLimits);
#endif

    /// function to move particles for testing purposes
    //void move(const double &dt, Domain &domain); // TODO: remove

    /// sanity check functions
    double sumVolume();
    double sumMass();
    double sumEnergy();
    double sumMomentumX();
    double sumMomentumY();
#if DIM==3
    double sumMomentumZ();
#endif
    void checkFluxSymmetry(Particles *ghostParticles=nullptr);

    void dump2file(std::string filename, double simTime);

private:
    int *nnl; // nearest neighbor list
    int *noi;  // number of interactions
    double *omega; // store omega to avoid recomputing
    double (*psijTilde_xi)[DIM];
    double (*Aij)[DIM];
    double (*WijL)[DIM+2], (*WijR)[DIM+2]; // DIM velocity components, density and pressure
    double (*Fij)[DIM+2];
    double (*vFrame)[DIM];

    double (*kernel)(const double&, const double&){ &Kernel::cubicSpline };

    /// variables for integration
    double *mF, *eF, (*vF)[DIM];

    void compOmega(int i, const double &kernelSize);
    void slopeLimiter(double *f, double (*grad)[DIM], const double &kernelSize,
                      Particles *ghostParticles, double *fGhost);


#if PERIODIC_BOUNDARIES
    void compOmega(int i, const Particles &ghostParticles, const double &kernelSize);
    int *nnlGhosts;
    int *noiGhosts;
    double (*psijTilde_xiGhosts)[DIM];
    double (*AijGhosts)[DIM];
    double (*WijLGhosts)[DIM+2], (*WijRGhosts)[DIM+2]; // DIM velocity components, density and pressure
    double (*FijGhosts)[DIM+2];
    double (*vFrameGhosts)[DIM];
    int *ghostMap;
    int *parent; // if class holds ghost particles the parent node is stored
#endif

    /// Helper variables for gradient estimation
    double B[DIM*DIM];
    double xi[DIM];
    double xj[DIM];
    double xjGhost[DIM];

    bool ghosts;
};


#endif //MESHLESSHYDRO_PARTICLES_H
