//
// Created by Johannes Martin on 21.09.22.
//

#include "../include/MeshlessScheme.h"

//TODO: add for 3 dimensions
double Kernel::cubicSpline(const double &r, const double &h) {
    const double sigma = 10./(7.*M_PI*h*h);
    const double q = r/h;
    if (0. <= q && q <= 1.){
        return sigma*(1.-3./2.*q*q*(1.-q/2.));
    } else if (1. < q && q < 2.){
        return sigma/4.*pow(2.-q, 3.);
    } else {
        return 0.;
    }
}

MeshlessScheme::MeshlessScheme(Configuration config, Particles *particles,
                               Domain::Cell bounds) : config { config }, particles { particles },
                                                      domain(bounds){

    Logger(INFO) << "    > Creating grid ... ";

    domain.createGrid(config.kernelSize);

    Logger(INFO) << "    > ... got " << domain.numGridCells << " cells";

    Logger(INFO) << "    > Assigning particles ...";
    particles->assignParticlesAndCells(domain);
    Logger(INFO) << "    > ... done.";
    Logger(INFO) << "    > Preparing simulation ...";
#if PERIODIC_BOUNDARIES
    Logger(INFO) << "    > Creating ghost particles ...";
    //Logger(DEBUG) << "      > Creating ghost grid";
    //domain.createGhostGrid();
    Logger(DEBUG) << "      > Creating ghost particles ... ";
    Particles ghostParticles { particles->N/(DIM*2) };
    particles->createGhostParticles(domain, ghostParticles, config.kernelSize);
    Logger(DEBUG) << "      > ... found " << ghostParticles.N << " ghosts";
    Logger(INFO) << "    > ... done.";

#endif
    Logger(DEBUG) << "      > Nearest neighbor search";
    particles->gridNNS(domain, config.kernelSize);
#if PERIODIC_BOUNDARIES
    particles->ghostNNS(domain, ghostParticles, config.kernelSize);
#endif
    Logger(DEBUG) << "      > Computing density";
#ifdef PERIODIC_BOUNDARIES
    particles->compDensity(ghostParticles, config.kernelSize);
#else
    particles->compDensity(config.kernelSize);
#endif
    Logger(DEBUG) << "      > Writing ICs to output";
    particles->dump2file(config.outDir + std::string("/ic.h5"));
}

void MeshlessScheme::run(){

}

MeshlessScheme::~MeshlessScheme(){}