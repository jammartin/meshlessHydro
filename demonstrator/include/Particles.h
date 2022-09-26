//
// Created by Johannes Martin on 05.09.22.
//

#ifndef MESHLESSHYDRO_PARTICLES_H
#define MESHLESSHYDRO_PARTICLES_H

#include <cmath>
#include <highfive/H5File.hpp>

#include "parameter.h"
#include "Domain.h"
#include "Helper.h"

namespace Kernel {
    double cubicSpline(const double &r, const double &h);
}

class Particles {

public:
    Particles(int numParticles);
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
    void compPressure(const double &gamma);

#if PERIODIC_BOUNDARIES
    void createGhostParticles(Domain &domain,
                              Particles &ghostParticles, const double &kernelSize);
    void ghostNNS(Domain &domain, const Particles &ghostParticles, const double &kernelSize);
    void compDensity(const Particles &ghostParticles, const double &kernelSize);
    void compPsijTilde(Helper &helper, const Particles &ghostParticles, const double &kernelSize);
    void gradient(double *f, double (*grad)[DIM], const Particles &ghostParticles);
#endif

    /// function to move particles for testing purposes
    void move(const double &dt, Domain &domain);

    void dump2file(std::string filename);

private:
    int *nnl; // nearest neighbor list
    int *noi;  // number of interactions
    double *omega; // store omega to avoid recomputing
    double (*psijTilde_xi)[DIM];
    void compOmega(int i, const double &kernelSize);
    double (*kernel)(const double&, const double&){ &Kernel::cubicSpline };


#if PERIODIC_BOUNDARIES
    void compOmega(int i, const Particles &ghostParticles, const double &kernelSize);
    int *nnlGhosts;
    int *noiGhosts;
#endif

    /// Helper variables for gradient estimation
    double B[DIM*DIM];
    double xi[DIM];
    double xj[DIM];
};


#endif //MESHLESSHYDRO_PARTICLES_H
