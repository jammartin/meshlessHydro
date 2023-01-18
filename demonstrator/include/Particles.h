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
#if DIM == 3
    double *z, *vz;
#endif
    void assignParticlesAndCells(Domain &domain);
    void gridNNS(Domain &domain, const double &kernelSize);
    void compDensity(const double &kernelSize); // also computes omega
    void compPsijTilde(Helper &helper, const double &kernelSize);
    void gradient(double *f, double (*grad)[DIM]);
    void slopeLimiter(const double &kernelSize,
                      Particles *ghostParticles=nullptr);
    void compPressure(const double &gamma);
    void compEffectiveFace();

    void compRiemannStatesLR(const double &dt, const double &kernelSize, const double &gamma);

    void solveRiemannProblems(const double &gamma, const Particles &ghostParticles);

    void collectFluxes(Helper &helper, const Particles &ghostParticles);

    void updateStateAndPosition(const double &dt, const Domain &domain);


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
