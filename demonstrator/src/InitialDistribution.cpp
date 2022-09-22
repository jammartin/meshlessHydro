//
// Created by Johannes Martin on 23.09.21.
//

#include "../include/InitialDistribution.h"

InitialDistribution::InitialDistribution(const std::string &file){
    HighFive::File h5file(file, HighFive::File::ReadOnly);

    // read datasets from file
    HighFive::DataSet mass = h5file.getDataSet("/m");
    HighFive::DataSet pos = h5file.getDataSet("/x");
    HighFive::DataSet vel = h5file.getDataSet("/v");
    HighFive::DataSet materialId = h5file.getDataSet("/materialId");
    HighFive::DataSet energy = h5file.getDataSet("/u");

    // read data into containers
    mass.read(m);
    pos.read(x);
    vel.read(v);
    energy.read(u);
    materialId.read(matId);

    // sanity check
    if (x.size() == v.size() && x.size() == m.size() && x.size() == u.size() && x.size() == matId.size()){
        numberOfParticles = x.size();
    } else {
        throw std::length_error("Length mismatch between mass, position and/or velocity vectors.");
    }
}

void InitialDistribution::getAllParticles(Particles &particles){
    std::vector<std::vector<double>>::iterator xit = x.begin();
    std::vector<std::vector<double>>::iterator vit = v.begin();
    std::vector<double>::iterator mit = m.begin();
    std::vector<double>::iterator uit = u.begin();
    std::vector<int>::iterator matIdIt = matId.begin();

    int pCounter = 0;

    while (xit != x.end()){
        particles.m[pCounter] = *mit;
        particles.u[pCounter] = *uit;
        particles.matId[pCounter] = *matIdIt;
        particles.x[pCounter] = (*xit)[0];
        particles.vx[pCounter] = (*vit)[0];
        particles.y[pCounter] = (*xit)[1];
        particles.vy[pCounter] = (*vit)[1];
#if DIM == 3
        particles.z[pCounter] = (*xit)[2];
        particles.vz[pCounter] = (*vit)[2];
#endif
        ++matIdIt;
        ++uit;
        ++xit;
        ++vit;
        ++mit;
        ++pCounter;
    }
}