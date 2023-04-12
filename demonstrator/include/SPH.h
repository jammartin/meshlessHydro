//
// Created by Johannes Martin on 21.09.22.
//

#ifndef DEMONSTRATOR_SPH_H
#define DEMONSTRATOR_SPH_H

#include <iomanip>

#include "parameter.h"
#include "InitialDistribution.h"
#include "Logger.h"
#include "Domain.h"
#include "Helper.h"

class SPH {

public:
    struct Configuration {
        std::string initFile;
        std::string outDir;
        double timeStep;
        double timeEnd;
        int h5DumpInterval;
        double periodicBoxLimits[2 * DIM];
        double kernelSize;
        double gamma; // adiabatic index
    };

    SPH(Configuration config, Particles *particles, Domain::Cell domain);
    ~SPH();

    void run();

private:
    Configuration config;
    Particles *particles;
#if PERIODIC_BOUNDARIES
    Particles ghostParticles;
#endif
    Domain domain;
    Helper helper {};
};


#endif //DEMONSTRATOR_SPH_H
