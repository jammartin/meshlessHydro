//
// Created by Johannes Martin on 23.09.21.
//

#include "../include/InitialDistribution.h"

InitialDistribution::InitialDistribution(const std::string &file, bool _traceMaterial)
    : traceMaterial { _traceMaterial }{
    HighFive::File h5file(file, HighFive::File::ReadOnly);

    // read datasets from file
    HighFive::DataSet mass = h5file.getDataSet("/m");
    HighFive::DataSet pos = h5file.getDataSet("/x");
    HighFive::DataSet vel = h5file.getDataSet("/v");
    HighFive::DataSet matId;

    // read data into containers
    mass.read(m);
    pos.read(x);
    vel.read(v);

    if(traceMaterial){
        matId = h5file.getDataSet("/materialId");
        matId.read(materialId);
    }

    // sanity check
    if (x.size() == v.size() && x.size() == m.size() ){
        numberOfParticles = x.size();
    } else {
        throw std::length_error("Length mismatch between mass, position and/or velocity vectors.");
    }
}

void InitialDistribution::getAllParticles(Particle *&particles){
    std::vector<std::vector<double>>::iterator xit = x.begin();
    std::vector<std::vector<double>>::iterator vit = v.begin();
    std::vector<double>::iterator mit = m.begin();
    std::vector<int>::iterator matIdIt = materialId.begin();

    int pCounter = 0;

    while (xit != x.end()){
        particles[pCounter].m = *mit;
        for (int d=0; d<global::dim; ++d){
            particles[pCounter].x[d] = (*xit)[d];
            particles[pCounter].v[d] = (*vit)[d];
        }
        if (traceMaterial){
            particles[pCounter].materialId = *matIdIt;
            ++matIdIt;
        }
        ++xit;
        ++vit;
        ++mit;
        ++pCounter;
    }
}

void InitialDistribution::getParticles(Particle *&particles, int offset, int amount){
    std::vector<std::vector<double>>::iterator xit = x.begin() + offset;
    std::vector<std::vector<double>>::iterator vit = v.begin() + offset;
    std::vector<double>::iterator mit = m.begin() + offset;
    std::vector<int>::iterator matIdIt = materialId.begin() + offset;

    int pCounter = 0;

    while (pCounter < amount){
        particles[pCounter].m = *mit;
        for (int d=0; d<global::dim; ++d){
            particles[pCounter].x[d] = (*xit)[d];
            particles[pCounter].v[d] = (*vit)[d];
        }
        if (traceMaterial){
            particles[pCounter].materialId = *matIdIt;
            ++matIdIt;
        }
        ++xit;
        ++vit;
        ++mit;
        ++pCounter;
    }
}