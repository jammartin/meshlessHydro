//
// Created by Johannes Martin on 05.09.22.
//

#ifndef MESHLESSHYDRO_PARTICLES_H
#define MESHLESSHYDRO_PARTICLES_H

#include <cmath>
#include <highfive/H5File.hpp>

#include "parameter.h"
#include "Domain.h"

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
    double *m, *u, *x, *y, *vx, *vy, *rho, *p;
#if DIM == 3
    double *z, *vz;
#endif
    void assignParticlesAndCells(Domain &domain);
    void gridNNS(Domain &domain, const double &kernelSize);
    void compDensity(const double &kernelSize);

#if PERIODIC_BOUNDARIES
    void createGhostParticles(Domain &domain,
                              Particles &ghostParticles, const double &kernelSize);
    void ghostNNS(Domain &domain, const Particles &ghostParticles, const double &kernelSize);
    void compDensity(const Particles &ghostParticles, const double &kernelSize);
#endif

    void dump2file(std::string filename);

private:
    int *nnl; // nearest neighbor list
    int *noi;  // number of interactions
    double compOmega(int i, const double &kernelSize);
    double (*kernel)(const double&, const double&){ &Kernel::cubicSpline };

#if PERIODIC_BOUNDARIES
    double compOmega(int i, const Particles &ghostParticles, const double &kernelSize);
    int *nnlGhosts;
    int *noiGhosts;
#endif


};


#endif //MESHLESSHYDRO_PARTICLES_H
