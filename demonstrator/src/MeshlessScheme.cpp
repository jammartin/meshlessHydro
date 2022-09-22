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
    Logger(DEBUG) << "      > Nearest neighbor search";
    particles->gridNNS(domain, config.kernelSize);
    Logger(DEBUG) << "      > Computing density";
    particles->compDensity(config.kernelSize);
    Logger(DEBUG) << "      > Writing ICs to output";
    particles->dump2file(config.outDir + std::string("/ic.h5"));
}

void MeshlessScheme::run(){

}

MeshlessScheme::~MeshlessScheme(){}