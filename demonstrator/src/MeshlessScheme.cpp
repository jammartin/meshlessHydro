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
                                                      domain(bounds), ghostParticles(particles->N/DIM){

    Logger(INFO) << "    > Creating grid ... ";
    domain.createGrid(config.kernelSize);
    Logger(INFO) << "    > ... got " << domain.numGridCells << " cells";
}

void MeshlessScheme::run(){
    double t = 0;
    int step = 0;
    do {
        std::stringstream stepss;
        stepss << std::setw(6) << std::setfill('0') << step;

        Logger(INFO) << "    > Assigning particles ...";
        particles->assignParticlesAndCells(domain);
        Logger(INFO) << "    > ... done.";
#if PERIODIC_BOUNDARIES
        Logger(INFO) << "    > Creating ghost particles ...";
        //Logger(DEBUG) << "      > Creating ghost grid";
        //domain.createGhostGrid();
        Logger(DEBUG) << "      > Creating ghost particles ... ";
        particles->createGhostParticles(domain, ghostParticles, config.kernelSize);
        Logger(DEBUG) << "      > ... found " << ghostParticles.N << " ghosts";
        Logger(INFO) << "    > ... done.";

#endif
        Logger(INFO) << "    > Nearest neighbor search";
        particles->gridNNS(domain, config.kernelSize);
#if PERIODIC_BOUNDARIES
        Logger(DEBUG) << "      > Ghosts NNS";
        particles->ghostNNS(domain, ghostParticles, config.kernelSize);
#endif
        Logger(INFO) << "    > Computing density";
        particles->compDensity(config.kernelSize);
#if PERIODIC_BOUNDARIES
        particles->compDensity(ghostParticles, config.kernelSize);
#endif
        Logger(INFO) << "    > Computing pressure";
        particles->compPressure(config.gamma);

        Logger(INFO) << "    > Computing gradients";
#if PERIODIC_BOUNDARIES
        particles->updateGhostState(ghostParticles);
        particles->compPsijTilde(helper, ghostParticles, config.kernelSize);
        particles->updateGhostPsijTilde(ghostParticles);
        particles->gradient(particles->rho, particles->rhoGrad, ghostParticles.rho, ghostParticles);
        particles->gradient(particles->vx, particles->vxGrad, ghostParticles.vx, ghostParticles);
        particles->gradient(particles->vy, particles->vyGrad, ghostParticles.vy, ghostParticles);
#if DIM == 3
        particles->gradient(particles->vz, particles->vzGrad, ghostParticles.vz, ghostParticles);
#endif
        particles->gradient(particles->P, particles->PGrad, ghostParticles.P, ghostParticles);
        Logger(DEBUG) << "      > Update ghost gradients";
        particles->updateGhostGradients(ghostParticles);
        // TODO: check how to properly limit gradiens
        Logger(DEBUG) << "      > Limiting slopes";
        particles->slopeLimiter(config.kernelSize, &ghostParticles);
        Logger(DEBUG) << "      > Update limited ghost gradients";
        particles->updateGhostGradients(ghostParticles);
#else
        particles->compPsijTilde(helper, config.kernelSize);
        particles->gradient(particles->rho, particles->rhoGrad);
        particles->gradient(particles->vx, particles->vxGrad);
        particles->gradient(particles->vy, particles->vyGrad);
#if DIM == 3
        particles->gradient(particles->vz, particles->vzGrad);
#endif
        particles->gradient(particles->P, particles->PGrad);
        // TODO: check how to properly limit gradiens
        Logger(DEBUG) << "      > Limiting slopes";
        particles->slopeLimiter(config.kernelSize);
#endif
        Logger(INFO) << "    > Preparing Riemann solver";
        Logger(DEBUG) << "      > Computing effective faces";
        particles->compEffectiveFace();
#if PERIODIC_BOUNDARIES
        particles->compEffectiveFace(ghostParticles);
#endif
        Logger(DEBUG) << "      > Computing fluxes";
        particles->compRiemannFluxes(config.timeStep, config.kernelSize, config.gamma);

#if PERIODIC_BOUNDARIES
        Logger(DEBUG) << "      > Computing ghost fluxes";
        particles->compRiemannFluxes(config.timeStep, config.kernelSize, config.gamma,
                                     ghostParticles);
#endif
        Logger(INFO) << "    > Solving Riemann problems";
        //TODO: continue here with implementation after gradient check
        //particles->solveRiemannProblems(config.gamma);

#if PERIODIC_BOUNDARIES
        //particles->solveRiemannProblems(ghostParticles);
#endif

        Logger(INFO) << "    > Dump particles to file";
        particles->dump2file(config.outDir + "/" + stepss.str() + std::string(".h5"));

        Logger(ERROR) << "Aborting for debugging.";
        exit(6);

        Logger(INFO) << "    > Moving particles";
        particles->move(config.timeStep, domain);

        t += config.timeStep;
        ++step;
    } while(t<=config.timeEnd);
}

MeshlessScheme::~MeshlessScheme(){}