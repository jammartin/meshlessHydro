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
    P = new double[numParticles];
    x = new double[numParticles];
    y = new double[numParticles];
    vx = new double[numParticles];
    vy = new double[numParticles];
    rhoGrad = new double[numParticles][DIM];
    vxGrad = new double[numParticles][DIM];
    vyGrad = new double[numParticles][DIM];
    vzGrad = new double[numParticles][DIM];
    PGrad = new double[numParticles][DIM];
#if DIM == 3
    z = new double[numParticles];
    vz = new double[numParticles];
#endif
    nnl = new int[numParticles*MAX_NUM_INTERACTIONS];
    noi = new int[numParticles];
    omega = new double[numParticles];
    psijTilde_xi = new double[numParticles][DIM];
#if PERIODIC_BOUNDARIES
    // estimated memory allocation
    nnlGhosts = new int[numParticles*MAX_NUM_GHOST_INTERACTIONS];
    noiGhosts = new int[numParticles];
#endif
}

Particles::~Particles(){
    delete[] matId;
    delete[] cell;
    delete[] m;
    delete[] u;
    delete[] P;
    delete[] rho;
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;
    delete[] rhoGrad;
    delete[] vxGrad;
    delete[] vyGrad;
    delete[] vzGrad;
    delete[] PGrad;
#if DIM == 3
    delete[] z;
    delete[] vz;
#endif
    delete[] nnl;
    delete[] noi;
    delete[] omega;
    delete[] psijTilde_xi;
#if PERIODIC_BOUNDARIES
    delete[] nnlGhosts;
    delete[] noiGhosts;
#endif
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

void Particles::gridNNS(Domain &domain, const double &kernelSize){
    // loop over particles
    for(int i=0; i<N; ++i){
        int numSearchCells = pow(3, DIM);
        // search for nearest neighbors in the particle cell and neighbor cells
        int cells[numSearchCells];
        // neighboring cells (including particles cell)
        domain.getNeighborCells(cell[i], cells);
        // do nearest neighbor search
        int noiBuf = 0;
        const double hSqr = kernelSize * kernelSize;
        for (int iNeighbor=0; iNeighbor<numSearchCells; ++iNeighbor){
            // loop over particle indices in all
            if (cells[iNeighbor] < 0){
                // handle ghost cells in external function
            } else {
                for(auto const &iPrtcl : domain.grid[cells[iNeighbor]].prtcls){
                    double dSqr = pow(x[iPrtcl] - x[i], 2)
                                  + pow(y[iPrtcl] - y[i], 2);
#if DIM == 3
                    dSqr += pow(z[iPrtcl] - z[i], 2);
#endif
                    if (dSqr < hSqr){
                        if(noiBuf >= MAX_NUM_INTERACTIONS){
                            Logger(ERROR) << "MAX_NUM_INTERACTIONS exceeded for particle "
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


void Particles::compDensity(const double &kernelSize){
    for(int i=0; i<N; ++i){
        compOmega(i, kernelSize);
        rho[i] = m[i]*omega[i];
    }
}

void Particles::compOmega(int i, const double &kernelSize){
    double omg = 0;
    for (int j=0; j<noi[i]; ++j){
        double dSqr = pow(x[i] - x[nnl[j+i*MAX_NUM_INTERACTIONS]], 2)
                    + pow(y[i] - y[nnl[j+i*MAX_NUM_INTERACTIONS]], 2);
#if DIM == 3
        dSqr += pow(z[i] - z[nnl[j+i*MAX_NUM_INTERACTIONS]], 2);
#endif
        double r = sqrt(dSqr);
        omg += kernel(r, kernelSize);
    }
    omega[i] = omg;
}

void Particles::compPsijTilde(Helper &helper, const double &kernelSize){

    for (int i=0; i<N; ++i){

        // reset buffer
        for (int k=0; k<DIM*DIM; ++k){
            B[k] = 0.;
        }

        for (int alpha = 0; alpha < DIM; ++alpha) {
            psijTilde_xi[i][alpha] = 0.;
        }

        xi[0] = x[i];
        xi[1] = y[i];
#if DIM==3
        xi[2] = z[i];
#endif

        for (int j=0; j<noi[i]; ++j){

            double dSqr = pow(x[i] - x[nnl[j+i*MAX_NUM_INTERACTIONS]], 2)
                          + pow(y[i] - y[nnl[j+i*MAX_NUM_INTERACTIONS]], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[nnl[j+i*MAX_NUM_INTERACTIONS]], 2);
#endif
            double r = sqrt(dSqr);
            double psij_xi = kernel(r, kernelSize)/omega[i];

            xj[0] = x[nnl[j+i*MAX_NUM_INTERACTIONS]];
            xj[1] = y[nnl[j+i*MAX_NUM_INTERACTIONS]];
#if DIM==3
            xj[2] = z[nnl[j+i*MAX_NUM_INTERACTIONS]];
#endif

            for (int alpha=0; alpha<DIM; ++alpha){
                for(int beta=0; beta<DIM; ++beta){
                    B[alpha+DIM*beta] += (xj[alpha] - xi[alpha])*(xj[beta] - xi[beta]) * psij_xi;
                }
            }
        }

        helper.inverseMatrix(B, DIM);

        for (int j=0; j<noi[i]; ++j) {

            double dSqr = pow(x[i] - x[nnl[j + i * MAX_NUM_INTERACTIONS]], 2)
                          + pow(y[i] - y[nnl[j + i * MAX_NUM_INTERACTIONS]], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[nnl[j+i*MAX_NUM_INTERACTIONS]], 2);
#endif
            double r = sqrt(dSqr);
            double psij_xi = kernel(r, kernelSize) / omega[i];

            for (int alpha = 0; alpha < DIM; ++alpha) {
                psijTilde_xi[i][alpha] = 0.;
                for (int beta = 0; beta < DIM; ++beta) {
                    psijTilde_xi[i][alpha] += B[alpha + beta * DIM] * (xj[beta] - xi[beta]) * psij_xi;
                }
            }
        }
    }
}

void Particles::gradient(double *f, double (*grad)[DIM]){
    for (int i=0; i<N; ++i) {
        for (int alpha = 0; alpha < DIM; ++alpha) {
            grad[i][alpha] = 0;
        }

        for (int j = 0; j < noi[i]; ++j) {
            for (int alpha = 0; alpha < DIM; ++alpha) {
                grad[i][alpha] += (f[i] - f[nnl[j + i * MAX_NUM_INTERACTIONS]]) * psijTilde_xi[i][alpha];
            }
        }
    }
}

void Particles::compPressure(const double &gamma){
    for (int i=0; i<N; ++i){
        P[i] = (gamma-1.)*rho[i]*u[i];
    }
}

#if PERIODIC_BOUNDARIES
void Particles::createGhostParticles(Domain &domain, Particles &ghostParticles,
                                     const double &kernelSize){
    int iGhost = 0;
    for(int i=0; i<N; ++i) {
        bool foundGhost = false;

        if (x[i] < domain.bounds.minX + kernelSize) {
            ghostParticles.x[iGhost] = domain.bounds.maxX + (x[i] - domain.bounds.minX);
            foundGhost = true;
        } else if (domain.bounds.maxX - kernelSize <= x[i]) {
            ghostParticles.x[iGhost] = domain.bounds.minX - (domain.bounds.maxX - x[i]);
            foundGhost = true;
        } else {
            ghostParticles.x[iGhost] = x[i];
        }
        if (y[i] < domain.bounds.minY + kernelSize) {
            ghostParticles.y[iGhost] = domain.bounds.maxY + (y[i] - domain.bounds.minY);
            foundGhost = true;
        } else if (domain.bounds.maxY - kernelSize <= y[i]) {
            ghostParticles.y[iGhost] = domain.bounds.minY - (domain.bounds.maxY - y[i]);
            foundGhost = true;
        } else {
            ghostParticles.y[iGhost] = y[i];
        }
        if (foundGhost) ++iGhost;
#if DIM == 3
        Logger(ERROR) << "Ghost cells not implemented for 3D simulations. - Aborting.";
        exit(2);
#endif
    }
    ghostParticles.N = iGhost;
}

void Particles::ghostNNS(Domain &domain, const Particles &ghostParticles, const double &kernelSize){
    const double hSqr = kernelSize * kernelSize;
    for(int i=0; i<N; ++i){
        int noiBuf = 0;
        for(int iGhost=0; iGhost<ghostParticles.N; ++iGhost){
            double dSqr = pow(ghostParticles.x[iGhost] - x[i], 2)
                          + pow(ghostParticles.y[iGhost] - y[i], 2);
#if DIM == 3
            Logger(ERROR) << "Ghost cells not implemented for 3D simulations. - Aborting.";
        exit(2);
#endif
            if (dSqr < hSqr){
                if(noiBuf >= MAX_NUM_GHOST_INTERACTIONS){
                    Logger(ERROR) << "MAX_NUM_GHOST_INTERACTIONS exceeded for particle "
                                  << i << " - Aborting.";
                    exit(3);
                }
                nnlGhosts[noiBuf+i*MAX_NUM_GHOST_INTERACTIONS] = iGhost;
                ++noiBuf;
            }
        }
        noiGhosts[i] = noiBuf;
    }
}

void Particles::compDensity(const Particles &ghostParticles, const double &kernelSize){
    for(int i=0; i<N; ++i){
        compOmega(i, ghostParticles, kernelSize);
        rho[i] = m[i]*omega[i];
    }
}

void Particles::compOmega(int i, const Particles &ghostParticles, const double &kernelSize){
    double omg = omega[i];
    for (int j=0; j<noiGhosts[i]; ++j){
        double dSqr = pow(x[i] - ghostParticles.x[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2)
                      + pow(y[i] - ghostParticles.y[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2);
#if DIM == 3
        dSqr += pow(z[i] - ghostParticles.z[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2);
#endif
        double r = sqrt(dSqr);
        omg += kernel(r, kernelSize);
    }
    omega[i] = omg;
}

void Particles::compPsijTilde(Helper &helper, const Particles &ghostParticles, const double &kernelSize){
    for (int i=0; i<N; ++i){

        // reset buffer
        for (int k=0; k<DIM*DIM; ++k){
            B[k] = 0.;
        }

        for (int alpha = 0; alpha < DIM; ++alpha) {
            psijTilde_xi[i][alpha] = 0.;
        }

        xi[0] = x[i];
        xi[1] = y[i];
#if DIM==3
        xi[2] = z[i];
#endif

        for (int j=0; j<noi[i]; ++j){

            double dSqr = pow(x[i] - x[nnl[j+i*MAX_NUM_INTERACTIONS]], 2)
                          + pow(y[i] - y[nnl[j+i*MAX_NUM_INTERACTIONS]], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[nnl[j+i*MAX_NUM_INTERACTIONS]], 2);
#endif
            double r = sqrt(dSqr);
            double psij_xi = kernel(r, kernelSize)/omega[i];

            xj[0] = x[nnl[j+i*MAX_NUM_INTERACTIONS]];
            xj[1] = y[nnl[j+i*MAX_NUM_INTERACTIONS]];
#if DIM==3
            xj[2] = z[nnl[j+i*MAX_NUM_INTERACTIONS]];
#endif

            for (int alpha=0; alpha<DIM; ++alpha){
                for(int beta=0; beta<DIM; ++beta){
                    B[alpha+DIM*beta] += (xj[alpha] - xi[alpha])*(xj[beta] - xi[beta]) * psij_xi;
                }
            }
        }

        for (int j=0; j<noiGhosts[i]; ++j){

            double dSqr = pow(x[i] - x[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2)
                          + pow(y[i] - y[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2);
#endif
            double r = sqrt(dSqr);
            double psij_xi = kernel(r, kernelSize)/omega[i];

            xj[0] = x[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]];
            xj[1] = y[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]];
#if DIM==3
            xj[2] = z[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]];
#endif

            for (int alpha=0; alpha<DIM; ++alpha){
                for(int beta=0; beta<DIM; ++beta){
                    B[alpha+DIM*beta] += (xj[alpha] - xi[alpha])*(xj[beta] - xi[beta]) * psij_xi;
                }
            }
        }

        helper.inverseMatrix(B, DIM);

        for (int j=0; j<noi[i]; ++j) {

            double dSqr = pow(x[i] - x[nnl[j + i * MAX_NUM_INTERACTIONS]], 2)
                          + pow(y[i] - y[nnl[j + i * MAX_NUM_INTERACTIONS]], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[nnl[j+i*MAX_NUM_INTERACTIONS]], 2);
#endif
            double r = sqrt(dSqr);
            double psij_xi = kernel(r, kernelSize) / omega[i];

            for (int alpha = 0; alpha < DIM; ++alpha) {
                psijTilde_xi[i][alpha] = 0.;
                for (int beta = 0; beta < DIM; ++beta) {
                    psijTilde_xi[i][alpha] += B[alpha + beta * DIM] * (xj[beta] - xi[beta]) * psij_xi;
                }
            }
        }

        for (int j=0; j<noiGhosts[i]; ++j) {

            double dSqr = pow(x[i] - x[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]], 2)
                          + pow(y[i] - y[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2);
#endif
            double r = sqrt(dSqr);
            double psij_xi = kernel(r, kernelSize) / omega[i];

            for (int alpha = 0; alpha < DIM; ++alpha) {
                for (int beta = 0; beta < DIM; ++beta) {
                    psijTilde_xi[i][alpha] += B[alpha + beta * DIM] * (xj[beta] - xi[beta]) * psij_xi;
                }
            }
        }
    }
}

void Particles::gradient(double *f, double (*grad)[DIM], const Particles &ghostParticles){
    for (int i=0; i<N; ++i) {
        for (int alpha = 0; alpha < DIM; ++alpha) {
            grad[i][alpha] = 0;
        }

        for (int j = 0; j < noi[i]; ++j) {
            for (int alpha = 0; alpha < DIM; ++alpha) {
                grad[i][alpha] += (f[i] - f[nnl[j + i * MAX_NUM_INTERACTIONS]]) * psijTilde_xi[i][alpha];
            }
        }

        for (int j = 0; j < noiGhosts[i]; ++j) {
            for (int alpha = 0; alpha < DIM; ++alpha) {
                grad[i][alpha] +=
                        (f[i] - f[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]]) * psijTilde_xi[i][alpha];
            }
        }
    }
}

#endif

void Particles::move(const double &dt, Domain &domain){
    for(int i=0; i<N; ++i) {
        x[i] = x[i] + vx[i] * dt;
        y[i] = y[i] + vy[i] * dt;
#if DIM == 3
        z[i] = z[i] +vz[i] * dt;
#endif
#if PERIODIC_BOUNDARIES
        if (x[i] <= domain.bounds.minX) {
            x[i] = domain.bounds.maxX - (domain.bounds.minX - x[i]);
        } else if (domain.bounds.maxX < x[i]) {
            x[i] = domain.bounds.minX + (x[i] - domain.bounds.maxX);
        }
        if (y[i] <= domain.bounds.minY) {
            y[i] = domain.bounds.maxY - (domain.bounds.minY - y[i]);
        } else if (domain.bounds.maxY < y[i]) {
            y[i] = domain.bounds.minY + (y[i] - domain.bounds.maxY);
        }
#if DIM ==3
        if (z[i] <= domain.bounds.minZ) {
            z[i] = domain.bounds.maxZ - (domain.bounds.minZ - z[i]);
        } else if (domain.bounds.maxZ < z[i]) {
            z[i] = domain.bounds.minZ + (z[i] - domain.bounds.maxZ);
        }
#endif
#endif
    }
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

