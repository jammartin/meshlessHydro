//
// Created by Johannes Martin on 23.09.21.
//

#ifndef MESHLESSHYDRO_INITIALDISTRIBUTION_H
#define MESHLESSHYDRO_INITIALDISTRIBUTION_H

#include <highfive/H5File.hpp>

#include "Particles.h"

class InitialDistribution {
public:
    InitialDistribution(const std::string &file);

    int getNumberOfParticles() const { return numberOfParticles; };
    void getAllParticles(Particles &particles);

private:
    // containers to be filled from hdf5 file
    std::vector<double> m {}, u {};
    std::vector<std::vector<double>> x {}, v {};
    std::vector<int> matId {};
    int numberOfParticles { 0 };
};


#endif //MESHLESSHYDRO_INITIALDISTRIBUTION_H
