//
// Created by Johannes Martin on 05.09.22.
//

#include "../include/Particles.h"

Particles::Particles(int numParticles) : N { numParticles }{
    // allocate memory
    matId = new int[numParticles];
    cell = new int[numParticles];
    m = new double[numParticles];
    u = new double[numParticles];
    rho = new double[numParticles];
    p = new double[numParticles];
    x = new double[numParticles];
    y = new double[numParticles];
    vx = new double[numParticles];
    vy = new double[numParticles];
#if DIM == 3
    z = new double[numParticles];
    vz = new double[numParticles];
#endif
    nnl = new int[numParticles*MAX_NUM_INTERACTIONS];
    noi = new int[numParticles];
}

Particles::~Particles(){
    delete[] matId;
    delete[] cell;
    delete[] m;
    delete[] u;
    delete[] p;
    delete[] rho;
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;
#if DIM == 3
    delete[] z;
    delete[] vz;
#endif
    delete[] nnl;
    delete[] noi;
}

void Particles::assignParticlesAndCells(Domain &domain){
    for(int i=0; i<N; ++i){
        int iGrid = floor((x[i]-domain.bounds.minX)/domain.cellSizeX)
                    + floor((y[i]-domain.bounds.minY)/domain.cellSizeY) * domain.cellsX
#if DIM == 3
        + floor((z[i]-domain.bounds.minZ)/domain.cellSizeZ) * domain.cellsX * domain.cellsY
#endif
        ;
        //Logger(DEBUG) << "      > Assigning particle " << i << " to cell " << iGrid;
        domain.grid[iGrid].prtcls.push_back(i);
        cell[i] = iGrid; // assign cells to particles
    }
}

void Particles::gridNNS(Domain &domain, double kernelSize){
    // loop over particles
    for(int i=0; i<N; ++i){
        int numSearchCells = pow(3, DIM);
        // search for nearest neighbors in the particle cell and neighbor cells
        int cells[numSearchCells];
        // neighboring cells (including particles cell)
        domain.getNeighborCells(cell[i], cells);
        // do nearest neighbor search
        int noiBuf = 0;
        double hSqr = kernelSize * kernelSize;
        for (int iNeighbor=0; iNeighbor<numSearchCells; ++iNeighbor){
            // loop over particle indices in all
            if (cells[iNeighbor] < 0){
                // handle ghost cells
            } else {
                for(auto const &iPrtcl : domain.grid[cells[iNeighbor]].prtcls){
                    double dSqr = pow(x[iPrtcl] - x[i], 2)
                                  + pow(y[iPrtcl] - y[i], 2);
#if DIM == 3
                    dSqr += pow(z[iPrtcl] - z[i], 2);
#endif
                    if (dSqr < hSqr){
                        if(noiBuf >= MAX_NUM_INTERACTIONS){
                            Logger(ERROR) << "MAX_NUM_INTERACTIONS exeeded for particle "
                                               << i << " - Aborting.";
                            exit(1);
                        }
                        nnl[noiBuf+i*MAX_NUM_INTERACTIONS] = iPrtcl;
                        ++noiBuf;
                    }
                }
            }
        }
        noi[i] = noiBuf;
    }
}

void Particles::compDensity(double kernelSize){
    for(int i=0; i<N; ++i){
        rho[i] = m[i]/compVolume(i, kernelSize);
    }
}

double Particles::compVolume(int i, const double &kernelSize){
    double omega = 0;
    for (int j=0; j<noi[i]; ++j){
        double dSqr = pow(x[i] - x[nnl[j+i*MAX_NUM_INTERACTIONS]], 2)
                    + pow(y[i] - y[nnl[j+i*MAX_NUM_INTERACTIONS]], 2);
#if DIM == 3
        dSqr += pow(z[i] - z[nnl[j+i*MAX_NUM_INTERACTIONS]], 2);
#endif
        double r = sqrt(dSqr);
        omega += kernel(r, kernelSize);
    }
    return 1./omega;
}

void Particles::dump2file(std::string filename){
    // open output file
    HighFive::File h5File { filename, HighFive::File::ReadWrite |
                                      HighFive::File::Create |
                                      HighFive::File::Truncate };

    // dimensions for datasets containing vectors
    std::vector<size_t> dataSpaceDims(2);
    dataSpaceDims[0] = std::size_t(N); // number of particles
    dataSpaceDims[1] = DIM;

    // create datasets
    HighFive::DataSet rhoDataSet = h5File.createDataSet<double>("/rho", HighFive::DataSpace(N));
    HighFive::DataSet posDataSet = h5File.createDataSet<double>("/x", HighFive::DataSpace(dataSpaceDims));

    // containers for particle data
    std::vector<double> rhoVec(rho, rho+N);
    std::vector<std::vector<double>> posVec(N);

    // fill containers with data
    for(int i=0; i<N; ++i){
        std::vector<double> posBuf(DIM);
        posBuf[0] = x[i];
        posBuf[1] = y[i];
#if DIM == 3
        posBuf[2] = z[i];
#endif
        posVec[i] = posBuf;
    }

    // write data
    rhoDataSet.write(rhoVec);
    posDataSet.write(posVec);
}

