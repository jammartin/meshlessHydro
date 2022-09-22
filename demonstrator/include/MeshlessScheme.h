//
// Created by Johannes Martin on 21.09.22.
//

#ifndef DEMONSTRATOR_MESHLESSSCHEME_H
#define DEMONSTRATOR_MESHLESSSCHEME_H

#include "parameter.h"
#include "InitialDistribution.h"
#include "Logger.h"
#include "Domain.h"

class MeshlessScheme {

public:
    struct Configuration {
        std::string initFile;
        std::string outDir;
        double timeStep;
        double timeEnd;
        int h5DumpInterval;
        double periodicBoxLimits[2 * DIM];
        double kernelSize;
    };

    MeshlessScheme(Configuration config, Particles *particles, Domain::Cell domain);
    ~MeshlessScheme();

    void run();

private:
    Configuration config;
    Particles *particles;
    Domain domain;
};


#endif //DEMONSTRATOR_MESHLESSSCHEME_H
