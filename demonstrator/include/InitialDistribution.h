//
// Created by Johannes Martin on 23.09.21.
//

#ifndef MESHLESSHYDRO_INITIALDISTRIBUTION_H
#define MESHLESSHYDRO_INITIALDISTRIBUTION_H

#include <highfive/H5File.hpp>

#include "Particles.h"

class InitialDistribution {
public:
    InitialDistribution(const std::string &file, bool _traceMaterial=false);

    int getNumberOfParticles() const { return numberOfParticles; };
    void getAllParticles(Particles &particles);
    //void getParticles(Particles &particles, int offset, int amount);

private:
    // containers to be filled from hdf5 file
    std::vector<double> m {};
    std::vector<std::vector<double>> x {}, v {};
    std::vector<int> materialId {};
    int numberOfParticles { 0 };
    bool traceMaterial;
};


#endif //MESHLESSHYDRO_INITIALDISTRIBUTION_H
