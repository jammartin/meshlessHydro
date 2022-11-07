//
// Created by Johannes Martin on 05.09.22.
//

#include "../include/Particles.h"

//TODO: add for 3 dimensions
double Kernel::cubicSpline(const double &r, const double &h) {

    // TODO: remove this
    double h2 = h/2.;

    const double sigma = 10./(7.*M_PI*h2*h2);
    const double q = r/h2;
    if (0. <= q && q <= 1.){
        return sigma*(1.-3./2.*q*q*(1.-q/2.));
    } else if (1. < q && q < 2.){
        return sigma/4.*pow(2.-q, 3.);
    } else {
        return 0.;
    }
}

Particles::Particles(int numParticles, bool ghosts) : N { numParticles }, ghosts { ghosts }{
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


    //TODO: check if this is needed as array
    //B = new double[numParticles][DIM*DIM];
#if DIM == 3
    z = new double[numParticles];
    vz = new double[numParticles];
#endif
    omega = new double[numParticles];
    if (!ghosts){
        nnl = new int[numParticles*MAX_NUM_INTERACTIONS];
        noi = new int[numParticles];
        psijTilde_xi = new double[numParticles*MAX_NUM_INTERACTIONS][DIM];
        Aij = new double[numParticles*MAX_NUM_INTERACTIONS][DIM];
        WijL = new double[numParticles*MAX_NUM_INTERACTIONS][DIM+2];
        WijR = new double[numParticles*MAX_NUM_INTERACTIONS][DIM+2];
        Fij = new double[numParticles*MAX_NUM_INTERACTIONS][DIM+2]; // TODO: this buffer should not be needed
        vFrame = new double[numParticles*MAX_NUM_INTERACTIONS][DIM];

        mF = new double[numParticles];
        vF = new double[numParticles][DIM];
        eF = new double[numParticles];

#if PERIODIC_BOUNDARIES
        // estimated memory allocation
        nnlGhosts = new int[numParticles*MAX_NUM_GHOST_INTERACTIONS];
        noiGhosts = new int[numParticles];
        ghostMap = new int[numParticles*(DIM+1)]; // TODO: this is only applicable for DIM==2
        AijGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM];
        WijLGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM+2];
        WijRGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM+2];
        FijGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM+2]; // TODO: this buffer should not be needed
        vFrameGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM];
#endif
    }
#if PERIODIC_BOUNDARIES
    psijTilde_xiGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM];
#endif
}

Particles::~Particles() {
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
    delete[] omega;

    if (!ghosts) {
        delete[] psijTilde_xi;
        delete[] WijL;
        delete[] WijR;
        delete[] Fij;
        delete[] vFrame;
        delete[] mF;
        delete[] vF;
        delete[] eF;
#if DIM == 3
        delete[] z;
        delete[] vz;
#endif
        delete[] nnl;
        delete[] noi;
        delete[] Aij;
#if PERIODIC_BOUNDARIES
        delete[] nnlGhosts;
        delete[] noiGhosts;
        delete[] AijGhosts;
        delete[] WijLGhosts;
        delete[] WijRGhosts;
        delete[] ghostMap;
        delete[] FijGhosts;
        delete[] vFrameGhosts;
#endif
    }
#if PERIDOIC_BOUNDARIES
    delete[] psijTilde_xiGhosts;
#endif
}

void Particles::assignParticlesAndCells(Domain &domain){

    // reset particles assigned to grid cells
    for(int iGrid=0; iGrid<domain.numGridCells; ++iGrid){
        domain.grid[iGrid].prtcls = std::vector<int>();
    }

    for(int i=0; i<N; ++i){

        int floorX = floor((x[i]-domain.bounds.minX)/domain.cellSizeX);
        int floorY = floor((y[i]-domain.bounds.minY)/domain.cellSizeY);

        if (x[i] == domain.bounds.maxX){
            floorX -= 1;
        }
        if (y[i] == domain.bounds.maxY){
            floorY -= 1;
        }

        int iGrid = floorX + floorY * domain.cellsX
#if DIM == 3
        + floor((z[i]-domain.bounds.minZ)/domain.cellSizeZ) * domain.cellsX * domain.cellsY
#endif
        ;

        //Logger(DEBUG) << "floor x = " << floorX
        //          << ", floor y = " << floorY;

        //Logger(DEBUG) << "      > Assigning particle@" << i << " = [" << x[i] << ", " << y[i]
        //          << "] to cell " << iGrid;
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
                    if (iPrtcl != i){
                        double dSqr = pow(x[iPrtcl] - x[i], 2)
                                      + pow(y[iPrtcl] - y[i], 2);
#if DIM == 3
                        dSqr += pow(z[iPrtcl] - z[i], 2);
#endif
                        if (dSqr < hSqr) {
                            if (noiBuf >= MAX_NUM_INTERACTIONS) {
                                Logger(ERROR) << "MAX_NUM_INTERACTIONS exceeded for particle "
                                              << i << " - Aborting.";
                                exit(1);
                            }
                            nnl[noiBuf + i * MAX_NUM_INTERACTIONS] = iPrtcl;
                            ++noiBuf;
                        }
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
        if(rho[i] <= 0.){
            Logger(WARN) << "Zero or negative density @" << i;
        }
    }
}

void Particles::compOmega(int i, const double &kernelSize){
    double omg = 0.;
    for (int j=0; j<noi[i]; ++j){
        double dSqr = pow(x[i] - x[nnl[j+i*MAX_NUM_INTERACTIONS]], 2)
                    + pow(y[i] - y[nnl[j+i*MAX_NUM_INTERACTIONS]], 2);
#if DIM == 3
        dSqr += pow(z[i] - z[nnl[j+i*MAX_NUM_INTERACTIONS]], 2);
#endif
        double r = sqrt(dSqr);
        omg += kernel(r, kernelSize);

        //Logger(DEBUG) << "x[" << i << "] = [" << x[i] << ", " <<  y[i] << "], x["
        //          << nnl[j+i*MAX_NUM_INTERACTIONS] << "] = [" << x[nnl[j+i*MAX_NUM_INTERACTIONS]] << ", "
        //          << y[nnl[j+i*MAX_NUM_INTERACTIONS]] << "]";

    }
    omega[i] = omg + kernel(0., kernelSize); // add self interaction to normalization factor
}

void Particles::compPsijTilde(Helper &helper, const double &kernelSize){
    for (int i=0; i<N; ++i){

        // reset buffer
        for (int k=0; k<DIM*DIM; ++k){
            B[k] = 0.;
        }

        //for (int alpha = 0; alpha < DIM; ++alpha) {
        //    psijTilde_xi[i][alpha] = 0.;
        //}

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
                    B[DIM*alpha+beta] += (xj[alpha] - xi[alpha])*(xj[beta] - xi[beta]) * psij_xi;
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

            xj[0] = x[nnl[j+i*MAX_NUM_INTERACTIONS]];
            xj[1] = y[nnl[j+i*MAX_NUM_INTERACTIONS]];

            for (int alpha = 0; alpha < DIM; ++alpha) {
                psijTilde_xi[j + i * MAX_NUM_INTERACTIONS][alpha] = 0.;
                for (int beta = 0; beta < DIM; ++beta) {
                    psijTilde_xi[j + i * MAX_NUM_INTERACTIONS][alpha] += B[DIM * alpha + beta] * (xj[beta] - xi[beta]) * psij_xi;
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
                grad[i][alpha] += (f[nnl[j + i * MAX_NUM_INTERACTIONS]] - f[i])
                                  * psijTilde_xi[j + i * MAX_NUM_INTERACTIONS][alpha];
            }
        }
    }
}

void Particles::compPressure(const double &gamma){
    for (int i=0; i<N; ++i){
        P[i] = (gamma-1.)*rho[i]*u[i];
        if(P[i] <= 0.){
            Logger(WARN) << "Zero or negative pressure @" << i;
        }
    }
}

void Particles::compEffectiveFace(){
    for (int i=0; i<N; ++i){
        for (int j=0; j<noi[i]; ++j){
            int ji = nnl[i*MAX_NUM_INTERACTIONS+j]; // index i of particle j
            // search neighbor i in nnl[] of j
            int ij = 0;
            for(ij=0; ij<noi[ji]; ++ij){
                if (nnl[ij+ji*MAX_NUM_INTERACTIONS] == i) break;
            }
            for (int alpha=0; alpha<DIM; ++alpha){
                Aij[i*MAX_NUM_INTERACTIONS+j][alpha] = 1./omega[i]*psijTilde_xi[i*MAX_NUM_INTERACTIONS+j][alpha]
                        - 1./omega[ji]*psijTilde_xi[ij+ji*MAX_NUM_INTERACTIONS][alpha];
            }
        }
    }
}

void Particles::slopeLimiter(const double &kernelSize,
                             Particles *ghostParticles){
    double *rhoGhost = nullptr, *vxGhost = nullptr, *vyGhost = nullptr,
#if DIM==3
        *vzGhost = nullptr,
#endif
    *PGhost = nullptr;
#if PERIODIC_BOUNDARIES
    rhoGhost = ghostParticles->rho;
    vxGhost = ghostParticles->vx;
    vyGhost = ghostParticles->vy;
#if DIM==3
    vzGhost = ghostParticles->vz;
#endif
    PGhost = ghostParticles->P;
#endif
    slopeLimiter(rho, rhoGrad, kernelSize, ghostParticles, rhoGhost);
    slopeLimiter(vx, vxGrad, kernelSize, ghostParticles, vxGhost);
    slopeLimiter(vy, vyGrad, kernelSize, ghostParticles, vyGhost);
#if DIM==3
    slopeLimiter(vz, vzGrad, kernelSize, ghostParticles, vzGhost);
#endif
    slopeLimiter(P, PGrad, kernelSize, ghostParticles, PGhost);
}

void Particles::slopeLimiter(double *f, double (*grad)[DIM], const double &kernelSize,
                             Particles *ghostParticles, double *fGhost){
    double xij[DIM], xijxi[DIM];

    for (int i=0; i<N; ++i){
        double psiMaxNgb { std::numeric_limits<double>::min() };
        double psiMinNgb { std::numeric_limits<double>::max() };
        double psiMaxMid { std::numeric_limits<double>::min() };
        double psiMinMid { std::numeric_limits<double>::max() };

        for (int jn=0; jn<noi[i]; ++jn) {
            int j = nnl[i * MAX_NUM_INTERACTIONS + jn];

            xij[0] = x[i] + kernelSize / 2. * (x[j] - x[i]);
            xij[1] = y[i] + kernelSize / 2. * (y[j] - y[i]);
            xijxi[0] = xij[0] - x[i];
            xijxi[1] = xij[1] - y[i];
#if DIM == 3
            xij[2] = z[i] + kernelSize/2. * (z[j] - z[i]);
            xijxi[2] = xij[2] - z[i];
#endif
            if (psiMaxNgb < f[j]) psiMaxNgb = f[j];
            if (psiMinNgb > f[j]) psiMinNgb = f[j];
            // reconstruct states at effective face
            double fij = f[i] + Helper::dotProduct(grad[i], xijxi);
            if (psiMaxMid < fij) psiMaxMid = fij;
            if (psiMinMid > fij) psiMinMid = fij;
        }
#if PERIODIC_BOUNDARIES
        for (int jn=0; jn<noiGhosts[i]; ++jn) {
            int j = nnlGhosts[i * MAX_NUM_GHOST_INTERACTIONS + jn];

            xij[0] = x[i] + kernelSize / 2. * (ghostParticles->x[j] - x[i]);
            xij[1] = y[i] + kernelSize / 2. * (ghostParticles->y[j] - y[i]);
            xijxi[0] = xij[0] - x[i];
            xijxi[1] = xij[1] - y[i];
#if DIM == 3
            xij[2] = z[i] + kernelSize/2. * (ghostParticles->z[j] - z[i]);
            xijxi[2] = xij[2] - z[i];
#endif
            if (psiMaxNgb < fGhost[j]) psiMaxNgb = fGhost[j];
            if (psiMinNgb > fGhost[j]) psiMinNgb = fGhost[j];
            // reconstruct states at effective face
            double fij = f[i] + Helper::dotProduct(grad[i], xijxi);
            if (psiMaxMid < fij) psiMaxMid = fij;
            if (psiMinMid > fij) psiMinMid = fij;
        }
#endif
        // update gradients if necessary
        double alphaMax = abs((psiMaxNgb - f[i]) / (psiMaxMid - f[i]));
        double alphaMin = abs((f[i] - psiMinNgb) / (f[i] - psiMinMid));

        /*if(i==4){
            Logger(DEBUG) << "alphaMin = " << alphaMin << ", alphaMax = " << alphaMax
                        << ", psiMinNgb = " << psiMinNgb << ", psiMaxNgb = " << psiMaxNgb
                        << ", psiMinNgb = " << psiMinMid << ", psiMaxNgb = " << psiMaxMid
                        << ", f[i] = " << f[i];
        }*/

        if (alphaMin <= alphaMax && BETA*alphaMin < 1.){
            Logger(DEBUG) << "        > Limiting gradient with alphaMin@" << i << ", alpha = "
                          << alphaMin << ", grad[0] = " << grad[i][0] << ", grad[1] = " << grad[i][1];
            grad[i][0] *= alphaMin;
            grad[i][1] *= alphaMin;
#if DIM==3
            grad[i][2] *= alphaMin;
#endif
        } else if (alphaMax <= alphaMin && BETA*alphaMax < 1.){
            Logger(DEBUG) << "        > Limiting gradient with alphaMax@" << i << ", alpha = "
                        << alphaMax << ", grad[0] = " << grad[i][0] << ", grad[1] = " << grad[i][1];
            grad[i][0] *= alphaMax;
            grad[i][1] *= alphaMax;
#if DIM==3
            grad[i][2] *= alphaMax;
#endif
        }
    }
}

void Particles::compRiemannFluxes(const double &dt, const double &kernelSize, const double &gamma){
    for (int i=0; i<N; ++i){
        double xij[DIM];
        //double vFrame[DIM];
        // helper vectors
        double xijxi[DIM], xjxi[DIM], xijxj[DIM];

        for (int jn=0; jn<noi[i]; ++jn){
            int j = nnl[i*MAX_NUM_INTERACTIONS+jn];

            xjxi[0] = x[j] - x[i];
            xjxi[1] = y[j] - y[i];

            xij[0] = x[i] + kernelSize/2. * (xjxi[0]);
            xij[1] = y[i] + kernelSize/2. * (xjxi[1]);

            xijxi[0] = xij[0] - x[i];
            xijxi[1] = xij[1] - y[i];

            xijxj[0] = xij[0] - x[j];
            xijxj[1] = xij[1] - y[j];
#if DIM==3
            xjxi[2] = z[j] - z[i];
            xij[2] = z[i] + kernelSize/2. * (xjxi[2]);
            xijxi[2] = xij[2] - z[i];
            xijxj[2] = xij[2] - z[j];
#endif
            double dotProd = Helper::dotProduct(xijxi, xjxi);
            double dSqr = Helper::dotProduct(xjxi, xjxi);

            int iW = i*MAX_NUM_INTERACTIONS+jn;

            vFrame[iW][0] = vx[i] + (vx[j]-vx[i]) * dotProd/dSqr;
            vFrame[iW][1] = vy[i] + (vy[j]-vy[i]) * dotProd/dSqr;
#if DIM==3
            vFrame[iW][2] = vz[i] + (vz[j]-vz[i]) * dotProd/dSqr;
#endif
            // boost frame to effective face
            WijR[iW][0] = rho[i];
            WijL[iW][0] = rho[j];
            WijR[iW][1] = P[i];
            WijL[iW][1] = P[j];
            WijR[iW][2] = vx[i] - vFrame[iW][0];
            WijL[iW][2] = vx[j] - vFrame[iW][0];
            WijR[iW][3] = vy[i] - vFrame[iW][1];
            WijL[iW][3] = vy[j] - vFrame[iW][1];
#if DIM == 3
            WijR[iW][4] = vz[i] - vFrame[iW][2];
            WijL[iW][4] = vz[j] - vFrame[iW][2];
#endif

            /*if(i==4){
                Logger(DEBUG) << "vFrame[0] = " << vFrame[0]
                            << ", vFrame[1] = " << vFrame[1]
                            << ", rhoGrad[i][0] = " << rhoGrad[i][0]
                            << ", rhoGrad[i][1] = " << rhoGrad[i][1]
                            << ", rho[i] = " << rho[i]
                            << ", xijxi = " << xijxi[0] << ", " << xijxi[1]
                            << ", xij = " << xij[0] << ", " << xij[1]
                            << ", xj = " << x[j] << ", " << x[j] << " @" << j;
                //exit(5);
            }*/
            // reconstruction at effective face
            WijR[iW][0] += Helper::dotProduct(rhoGrad[i], xijxi);
            WijL[iW][0] += Helper::dotProduct(rhoGrad[j], xijxj);
            WijR[iW][1] += Helper::dotProduct(PGrad[i], xijxi);
            WijL[iW][1] += Helper::dotProduct(PGrad[j], xijxj);
            WijR[iW][2] += Helper::dotProduct(vxGrad[i], xijxi);
            WijL[iW][2] += Helper::dotProduct(vxGrad[j], xijxj);
            WijR[iW][3] += Helper::dotProduct(vyGrad[i], xijxi);
            WijL[iW][3] += Helper::dotProduct(vyGrad[j], xijxj);
#if DIM==3
            WijR[iW][4] += Helper::dotProduct(vzGrad[i], xijxi);
            WijL[iW][4] += Helper::dotProduct(vzGrad[j], xijxj);
#endif
            // predict half a timestep
            double viDiv = vxGrad[i][0] + vyGrad[i][1];
            double vjDiv = vxGrad[j][0] + vyGrad[j][1];
#if DIM==3
            viDiv += vzGrad[i][2];
            vjDiv += vzGrad[j][2];
#endif
            WijR[iW][0] -= dt/2. * rho[i] * viDiv;
            WijL[iW][0] -= dt/2. * rho[j] * vjDiv;
            WijR[iW][1] -= dt/2. * gamma*P[i] * viDiv;
            WijL[iW][1] -= dt/2. * gamma*P[j] * vjDiv;
            WijR[iW][2] -= dt/2. * PGrad[i][0]/rho[i];
            WijL[iW][2] -= dt/2. * PGrad[j][0]/rho[j];
            WijR[iW][3] -= dt/2. * PGrad[i][1]/rho[i];
            WijL[iW][3] -= dt/2. * PGrad[j][1]/rho[j];
#if DIM==3
            WijR[iW][4] -= dt/2. * PGrad[i][2]/rho[i];
            WijL[iW][4] -= dt/2. * PGrad[j][2]/rho[j];
#endif
        }
    }
}


void Particles::solveRiemannProblems(const double &gamma){
    for (int i=0; i<N; ++i){

        if (!(i % (N/VERBOSITY_PARTICLES))){
            Logger(DEBUG) << "        > i = " << i;
        }

        for (int j=0; j<noi[i]; ++j){
            int ii = i*MAX_NUM_INTERACTIONS+j; // interaction index
            Riemann solver { WijR[ii], WijL[ii], Aij[ii] };
            solver.exact(Fij[ii], gamma);
        }

#if PERIODIC_BOUNDARIES
        for (int j=0; j<noiGhosts[i]; ++j){
            int ii = i*MAX_NUM_GHOST_INTERACTIONS+j; // interaction index
            Riemann solver { WijRGhosts[ii], WijLGhosts[ii], AijGhosts[ii] };
            solver.exact(FijGhosts[ii], gamma);
        }
#endif
    }
}

void Particles::collectFluxes(Helper &helper){
    for (int i=0; i<N; ++i){
        mF[i] = 0.;

        vF[i][0] = 0.;
        vF[i][1] = 0.;
#if DIM == 3
        vF[i][2] = 0.;
#endif
        eF[i] = 0.;

        for(int j=0; j<noi[i]; ++j){
            int ii = j+i*MAX_NUM_INTERACTIONS;
            double AijNorm = sqrt(Helper::dotProduct(Aij[ii], Aij[ii]));

            /// MASS FLUXES
            mF[i] += AijNorm*Fij[ii][0];

            /// VELOCITY FLUXES
            // TODO: encapsulate block below into function
            double hatAij[DIM];

            hatAij[0] = 1./AijNorm*Aij[ii][0];
            hatAij[1] = 1./AijNorm*Aij[ii][1];
#if DIM==3
            hatAij[2] = 1./AijNorm*Aij[ii][2];
#endif

            // compute norm of velocity
            double unitX[DIM] = { 1, 0
#if DIM==3
                    ,0
#endif
            };

#if DIM==3
            Logger(ERROR) << "Rotation of states is not yet implemented in 3D. - Aborting.";
            exit(4);
#endif

            // rotate states back to reference frame

            double LambdaInv[DIM*DIM];
            Helper::rotationMatrix2D(unitX, hatAij, LambdaInv);
            helper.inverseMatrix(LambdaInv, DIM);
            //TODO: end block

            double vFrotate[DIM];
            vFrotate[0] = LambdaInv[0]*Fij[ii][2]+LambdaInv[1]*Fij[ii][3];
            vFrotate[1] = LambdaInv[2]*Fij[ii][2]+LambdaInv[3]*Fij[ii][3];

            // allocate buffer for energy update
            double Fv[DIM];
            Fv[0] = Helper::dotProduct(hatAij, vFrotate) + Fij[ii][0]*vFrame[ii][0];
            Fv[1] = Helper::dotProduct(hatAij, vFrotate) + Fij[ii][0]*vFrame[ii][1];

            vF[i][0] += AijNorm*Fv[0];
            vF[i][1] += AijNorm*Fv[1];

            // TODO: implement z-component to work for 3D

            /// ENERGY FLUXES
            eF[i] += AijNorm*(Fij[ii][1] + .5*Helper::dotProduct(vFrame[ii], vFrame[ii])*Fij[ii][0]
                                + Helper::dotProduct(vFrame[ii], Fv));

        }

#if PERIODIC_BOUNDARIES
        for(int j=0; j<noiGhosts[i]; ++j){
            int ii = j+i*MAX_NUM_GHOST_INTERACTIONS;
            double AijNorm = sqrt(Helper::dotProduct(AijGhosts[ii], AijGhosts[ii]));

            /// MASS FLUXES
            mF[i] += AijNorm*FijGhosts[ii][0];
            // TODO: encapsulate block below into function
            double hatAij[DIM];

            hatAij[0] = 1./AijNorm*AijGhosts[ii][0];
            hatAij[1] = 1./AijNorm*AijGhosts[ii][1];
#if DIM==3
            hatAij[2] = 1./AijNorm*AijGhosts[ii][2];
#endif

            // compute norm of velocity
            double unitX[DIM] = { 1, 0
#if DIM==3
                    ,0
#endif
            };

#if DIM==3
            Logger(ERROR) << "Rotation of states is not yet implemented in 3D. - Aborting.";
    exit(4);
#endif

            // rotate states back to reference frame

            double LambdaInv[DIM*DIM];
            Helper::rotationMatrix2D(unitX, hatAij, LambdaInv);
            helper.inverseMatrix(LambdaInv, DIM);
            //TODO: end block

            double vFrotate[DIM];
            vFrotate[0] = LambdaInv[0]*FijGhosts[ii][2]+LambdaInv[1]*FijGhosts[ii][3];
            vFrotate[1] = LambdaInv[2]*FijGhosts[ii][2]+LambdaInv[3]*FijGhosts[ii][3];

            // allocate buffer for energy update
            double Fv[DIM];
            Fv[0] = Helper::dotProduct(hatAij, vFrotate) + FijGhosts[ii][0]*vFrameGhosts[ii][0];
            Fv[1] = Helper::dotProduct(hatAij, vFrotate) + FijGhosts[ii][0]*vFrameGhosts[ii][1];

            vF[i][0] += AijNorm*Fv[0];
            vF[i][1] += AijNorm*Fv[1];

            // TODO: implement z-component to work for 3D

            /// ENERGY FLUXES
            eF[i] += AijNorm*(FijGhosts[ii][1] + .5*Helper::dotProduct(vFrameGhosts[ii], vFrameGhosts[ii])*FijGhosts[ii][0]
                              + Helper::dotProduct(vFrame[ii], Fv));
        }
#endif
    }
}

void Particles::updateStateAndPosition(const double &dt){

    for(int i=0; i<N; ++i){

        // store velocity for second order position update
        double vxi = vx[i], vyi = vy[i];
#if DIM==3
        double vzi = vz[i];
#endif

        // create state vector without mass
        double Q[DIM+1];
#if DIM==3
        Q[0] = m[i]*u[i] + .5*m[i]*(vxi*vxi+vyi*vyi+vzi*vzi);
#else
        Q[0] = m[i]*u[i] + .5*m[i]*(vxi*vxi+vyi*vyi);
#endif
        Q[1] = m[i]*vxi;
        Q[2] = m[i]*vyi;
#if DIM==3
        Q[3] = m[i]*vzi;
#endif

        // UPDATE MASS
        m[i] -= dt*mF[i];

        // UPDATE VELOCITY
        Q[1] -= dt*vF[i][0];
        Q[2] -= dt*vF[i][1];
#if DIM == 3
        Q[3] -= dt*vF[i][2];
#endif
        vx[i] = Q[1]/m[i];
        vy[i] = Q[2]/m[i];
#if DIM==3
        vz[i] = Q[3]/m[i];
#endif
        // UPDATE ENERGY
        Q[0] -= dt*eF[i];
#if DIM==3
        u[i] = (Q[0]-.5*m[i]*(vxi*vxi+vyi*vyi+vzi*vzi))/m[i];
#else
        u[i] = (Q[0]-.5*m[i]*(vxi*vxi+vyi*vyi))/m[i];
#endif
    }
}

#if PERIODIC_BOUNDARIES
void Particles::createGhostParticles(Domain &domain, Particles &ghostParticles,
                                     const double &kernelSize){
    int iGhost = 0;
    for(int i=0; i<N; ++i) {
        bool foundGhostX = false, foundGhostY = false;
        ghostMap[i*(DIM+1)] = -1;
        ghostMap[i*(DIM+1)+1] = -1;
        ghostMap[i*(DIM+1)+2] = -1;

        // x-direction
        if (x[i] <= domain.bounds.minX + kernelSize){ // && x[i] > domain.bounds.minX) {
            ghostParticles.x[iGhost] = domain.bounds.maxX + (x[i] - domain.bounds.minX);
            foundGhostX = true;
        } else if (domain.bounds.maxX - kernelSize < x[i]){ // && x[i] < domain.bounds.maxX) {
            ghostParticles.x[iGhost] = domain.bounds.minX - (domain.bounds.maxX - x[i]);
            foundGhostX = true;
        } else {
            ghostParticles.x[iGhost] = x[i];
        }

        // y-direction
        if (y[i] <= domain.bounds.minY + kernelSize){ // && y[i] > domain.bounds.minY) {
            ghostParticles.y[iGhost] = domain.bounds.maxY + (y[i] - domain.bounds.minY);
            foundGhostY = true;
        } else if (domain.bounds.maxY - kernelSize < y[i]){ // && y[i] < domain.bounds.maxY) {
            ghostParticles.y[iGhost] = domain.bounds.minY - (domain.bounds.maxY - y[i]);
            foundGhostY = true;
        } else {
            ghostParticles.y[iGhost] = y[i];
        }

        // 'corner' particle first if both are true
        if (foundGhostX || foundGhostY) {
            ghostMap[i*(DIM+1)] = iGhost;
            //Logger(DEBUG) << "particle@" << i << " = [" << x[i] << ", " << y[i] << "] makes "
            //          << "ghost@" << iGhost << " = [" << ghostParticles.x[iGhost] << ", "
            //          << ghostParticles.y[iGhost] << "]";
            ++iGhost;
        }

        // create DIM extra normal particles
        if (foundGhostX && foundGhostY){
            ghostParticles.x[iGhost] = x[i];
            if (y[i] <= domain.bounds.minY + kernelSize){ // && y[i] > domain.bounds.minY) {
                ghostParticles.y[iGhost] = domain.bounds.maxY + (y[i] - domain.bounds.minY);
            } else if (domain.bounds.maxY - kernelSize < y[i]){ // && y[i] < domain.bounds.maxY) {
                ghostParticles.y[iGhost] = domain.bounds.minY - (domain.bounds.maxY - y[i]);
            }
            ghostMap[i*(DIM+1)+1] = iGhost;
            //Logger(DEBUG) << "particle@" << i << " = [" << x[i] << ", " << y[i] << "] makes "
            //              <<"ghost@" << iGhost << " = [" << ghostParticles.x[iGhost] << ", "
            //              << ghostParticles.y[iGhost] << "]";
            ++iGhost;
            if (x[i] <= domain.bounds.minX + kernelSize){ // && x[i] > domain.bounds.minX) {
                ghostParticles.x[iGhost] = domain.bounds.maxX + (x[i] - domain.bounds.minX);
            } else if (domain.bounds.maxX - kernelSize < x[i]){ // && x[i] < domain.bounds.maxX) {
                ghostParticles.x[iGhost] = domain.bounds.minX - (domain.bounds.maxX - x[i]);
            }
            ghostParticles.y[iGhost] = y[i];
            ghostMap[i*(DIM+1)+2] = iGhost;
            //Logger(DEBUG) << "particle@" << i << " = [" << x[i] << ", " << y[i] << "] makes "
            //              << "ghost@" << iGhost << " = [" << ghostParticles.x[iGhost] << ", "
            //              << ghostParticles.y[iGhost] << "]";
            ++iGhost;
        }

        //foundGhostX = false;
        //foundGhostY = false;

#if DIM == 3
        Logger(ERROR) << "Ghost cells not implemented for 3D simulations. - Aborting.";
        exit(2);
#endif
    }
    ghostParticles.N = iGhost;
}

void Particles::updateGhostState(Particles &ghostParticles){
    for (int i=0; i<N*(DIM+1); ++i){
        if (ghostMap[i] >= 0){
            ghostParticles.rho[ghostMap[i]] = rho[i/(DIM+1)];
            ghostParticles.P[ghostMap[i]] = P[i/(DIM+1)];
            ghostParticles.omega[ghostMap[i]] = omega[i/(DIM+1)];
            ghostParticles.vx[ghostMap[i]] = vx[i/(DIM+1)];
            ghostParticles.vx[ghostMap[i]] = vy[i/(DIM+1)];
#if DIM==3
            ghostParticles.vz[ghostMap[i]] = vz[i/(DIM+1)];
#endif
        }
    }
}

void Particles::updateGhostGradients(Particles &ghostParticles){
    for (int i=0; i<N*(DIM+1); ++i){
        if (ghostMap[i] >= 0){
            for (int alpha=0; alpha<DIM; ++alpha){
                ghostParticles.rhoGrad[ghostMap[i]][alpha] = rhoGrad[i/(DIM+1)][alpha];
                ghostParticles.vxGrad[ghostMap[i]][alpha] = vxGrad[i/(DIM+1)][alpha];
                ghostParticles.vyGrad[ghostMap[i]][alpha] = vyGrad[i/(DIM+1)][alpha];
#if DIM == 3
                ghostParticles.vzGrad[ghostMap[i]][alpha] = vzGrad[i/(DIM+1)][alpha];
#endif
                ghostParticles.PGrad[ghostMap[i]][alpha] = PGrad[i/(DIM+1)][alpha];
            }
        }
    }
}

void Particles::updateGhostPsijTilde(Particles &ghostParticles){
    for (int i=0; i<N*(DIM+1); ++i){
        if (ghostMap[i] >= 0){
            // below should not be needed
            //for(int j=0; j<noi[i]; ++j){
            //    for (int alpha=0; alpha<DIM; ++alpha){
            //        ghostParticles.psijTilde_xi[j+ghostMap[i]*MAX_NUM_INTERACTIONS][alpha] = psijTilde_xi[j+i/(DIM+1)*MAX_NUM_INTERACTIONS][alpha];
            //    }
            //}
            for (int j=0; j<noiGhosts[i/(DIM+1)]; ++j){
                for (int alpha=0; alpha<DIM; ++alpha){
                    ghostParticles.psijTilde_xiGhosts[j+ghostMap[i]*MAX_NUM_GHOST_INTERACTIONS][alpha]
                                    = psijTilde_xiGhosts[j+i/(DIM+1)*MAX_NUM_GHOST_INTERACTIONS][alpha];
                }
            }
        }
    }
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
        if(rho[i] <= 0.){
            Logger(WARN) << "Zero or negative density @" << i;
        }
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

        //Logger(DEBUG) << "x[" << i << "] = [" << x[i] << ", " <<  y[i] << "], x["
        //              << nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS] << "] = [" << ghostParticles.x[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]] << ", "
        //              << ghostParticles.y[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]] << "]";
    }
    omega[i] = omg;
    //Logger(DEBUG) << "V[" << i << "] = " << 1./omega[i] << ", noiTot = " << noi[i] + noiGhosts[i]
    //          << " noiGhosts = " << noiGhosts[i];
}

void Particles::compPsijTilde(Helper &helper, const Particles &ghostParticles, const double &kernelSize){
    for (int i=0; i<N; ++i){

        // reset buffer
        for (int k=0; k<DIM*DIM; ++k){
            B[k] = 0.;
        }

        //for (int alpha = 0; alpha < DIM; ++alpha) {
        //    psijTilde_xi[i][alpha] = 0.;
        //}

        xi[0] = x[i];
        xi[1] = y[i];
#if DIM==3
        xi[2] = z[i];
#endif

        //Logger(DEBUG) << "In compPsijTilde(): noi[" << i << "] = " << noi[i]
        //        << ", noiGhosts[" << i << "] = " << noiGhosts[i];

        for (int j=0; j<noi[i]; ++j){

            double dSqr = pow(x[i] - x[nnl[j+i*MAX_NUM_INTERACTIONS]], 2)
                          + pow(y[i] - y[nnl[j+i*MAX_NUM_INTERACTIONS]], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[nnl[j+i*MAX_NUM_INTERACTIONS]], 2);
#endif
            double r = sqrt(dSqr);
            //double psij_xi = kernel(r, kernelSize)/omega[nnl[j + i * MAX_NUM_INTERACTIONS]];
            double psij_xi = kernel(r, kernelSize)/omega[i];

            xj[0] = x[nnl[j+i*MAX_NUM_INTERACTIONS]];
            xj[1] = y[nnl[j+i*MAX_NUM_INTERACTIONS]];
#if DIM==3
            xj[2] = z[nnl[j+i*MAX_NUM_INTERACTIONS]];
#endif

            //if (i == 7){
            //    Logger(DEBUG) << "psij_xi = " << psij_xi << ", xj = [" << xj[0] << ", " << xj[1] << "]"
            //                  << "; xi = [" << xi[0] << ", " << xi[1] << "]";
            //}

            for (int alpha=0; alpha<DIM; ++alpha){
                for(int beta=0; beta<DIM; ++beta){
                    B[DIM*alpha+beta] += (xj[alpha] - xi[alpha])*(xj[beta] - xi[beta]) * psij_xi;
                }
            }
        }

        for (int j=0; j<noiGhosts[i]; ++j){

            double dSqr = pow(x[i] - ghostParticles.x[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2)
                          + pow(y[i] - ghostParticles.y[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2);
#if DIM == 3
            dSqr += pow(z[i] - ghostParticles.z[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2);
#endif
            double r = sqrt(dSqr);
            //double psij_xi = kernel(r, kernelSize)/ghostParticles.omega[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]];
            double psij_xi = kernel(r, kernelSize)/omega[i];

            xjGhost[0] = ghostParticles.x[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]];
            xjGhost[1] = ghostParticles.y[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]];
#if DIM==3
            xj[2] = ghostParticles.z[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]];
#endif
            //if (i == 7){
            //    Logger(DEBUG) << "psij_xi = " << psij_xi << ", xjGhost = [" << xjGhost[0] << ", " << xjGhost[1] << "]"
            //                << "; xi = [" << xi[0] << ", " << xi[1] << "]";
            //}

            for (int alpha=0; alpha<DIM; ++alpha){
                for(int beta=0; beta<DIM; ++beta){
                    B[DIM*alpha+beta] += (xjGhost[alpha] - xi[alpha])*(xjGhost[beta] - xi[beta]) * psij_xi;
                }
            }
        }

        //Logger(DEBUG) << "E = " << "[" << B[0] << ", " << B[1] << ", " << B[2] << ", " << B[3] << "]";

        // needed for sanity check of matrix E
        //double normE = 0;
        //for (int alpha=0; alpha<DIM; ++alpha){
        //    for (int beta=0; beta<DIM; ++beta){
        //        normE += B[alpha*DIM+beta]*B[alpha*DIM+beta];
        //    }
        //}

        helper.inverseMatrix(B, DIM);

        //double normB = 0;
        //for (int alpha=0; alpha<DIM; ++alpha){
        //    for (int beta=0; beta<DIM; ++beta){
        //        normB += B[alpha*DIM+beta]*B[alpha*DIM+beta];
        //    }
        //}

        // Check whether Matrix E is ill-conditioned
        //double Ncond = 1./DIM * sqrt(normE*normB);
        //Logger(DEBUG) << "Ncond@" << i << " = " << Ncond;

        //if (i == 7) {
        //    Logger(DEBUG) << "B = " << "[" << B[0] << ", " << B[1] << ", " << B[2] << ", " << B[3] << "]";
        //}
        //Logger(DEBUG) << "noi[" << i << "] = " << noi[i] << ", noiGhosts[" << i << "] = " << noiGhosts[i];
        //exit(7);

        for (int j=0; j<noi[i]; ++j) {

            double dSqr = pow(x[i] - x[nnl[j + i * MAX_NUM_INTERACTIONS]], 2)
                          + pow(y[i] - y[nnl[j + i * MAX_NUM_INTERACTIONS]], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[nnl[j+i*MAX_NUM_INTERACTIONS]], 2);
#endif
            double r = sqrt(dSqr);
            //double psij_xi = kernel(r, kernelSize) / omega[nnl[j + i * MAX_NUM_INTERACTIONS]];
            double psij_xi = kernel(r, kernelSize) / omega[i];

            xj[0] = x[nnl[j+i*MAX_NUM_INTERACTIONS]];
            xj[1] = y[nnl[j+i*MAX_NUM_INTERACTIONS]];

            for (int alpha = 0; alpha < DIM; ++alpha) {
                //psijTilde_xi[nnl[j + i * MAX_NUM_INTERACTIONS]+i*MAX_NUM_INTERACTIONS][alpha] = 0.;
                psijTilde_xi[j + i * MAX_NUM_INTERACTIONS][alpha] = 0.;
                for (int beta = 0; beta < DIM; ++beta) {
                    //psijTilde_xi[nnl[j + i * MAX_NUM_INTERACTIONS]+i*MAX_NUM_INTERACTIONS][alpha] += B[alpha * DIM + beta] * (xj[beta] - xi[beta]) * psij_xi;
                    psijTilde_xi[j + i * MAX_NUM_INTERACTIONS][alpha] += B[alpha * DIM + beta] * (xj[beta] - xi[beta]) * psij_xi;
                }

                //if(i == 86){
                //    Logger(DEBUG) << "psijTilde_xi[" << alpha << "]@" << nnl[j + i * MAX_NUM_INTERACTIONS] << " = " << psijTilde_xi[nnl[j + i * MAX_NUM_INTERACTIONS]][alpha];
                //}
            }
        }

        for (int j=0; j<noiGhosts[i]; ++j) {

            double dSqr = pow(x[i] - ghostParticles.x[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]], 2)
                          + pow(y[i] - ghostParticles.y[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]], 2);
#if DIM == 3
            dSqr += pow(z[i] - ghostParticles.z[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]], 2);
#endif
            double r = sqrt(dSqr);
            //double psij_xi = kernel(r, kernelSize) / ghostParticles.omega[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]];
            double psij_xi = kernel(r, kernelSize) / omega[i];

            xjGhost[0] = ghostParticles.x[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]];
            xjGhost[1] = ghostParticles.y[nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS]];

            for (int alpha = 0; alpha < DIM; ++alpha) {
                //psijTilde_xiGhosts[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]+i*MAX_NUM_GHOST_INTERACTIONS][alpha] = 0.;
                psijTilde_xiGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS][alpha] = 0.;
                for (int beta = 0; beta < DIM; ++beta) {
                    //psijTilde_xiGhosts[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]+i*MAX_NUM_GHOST_INTERACTIONS][alpha] += B[alpha * DIM + beta] * (xjGhost[beta] - xi[beta]) * psij_xi;
                    psijTilde_xiGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS][alpha] += B[alpha * DIM + beta] * (xjGhost[beta] - xi[beta]) * psij_xi;
                }
                //if(i == 86){
                //    Logger(DEBUG) << "psijTildeGhost_xi[" << alpha << "]@" << nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS] << " = " << ghostParticles.psijTilde_xi[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]][alpha];
                //}
            }
        }
    }
}

void Particles::gradient(double *f, double (*grad)[DIM], double *fGhost, const Particles &ghostParticles){
    for (int i=0; i<N; ++i) {
        for (int alpha = 0; alpha < DIM; ++alpha) {
            grad[i][alpha] = 0;
        }

        //Logger(DEBUG) << "      > noi[" << i << "] = " << noi[i]
        //              << ", noiGhosts[" << i << "] = " << noiGhosts[i];

        //if(i == 86) Logger(DEBUG) << "fi[" << i << "] = " << f[i];

        for (int j = 0; j < noi[i]; ++j) {
            for (int alpha = 0; alpha < DIM; ++alpha) {
                grad[i][alpha] += (f[nnl[j + i * MAX_NUM_INTERACTIONS]] - f[i])
                                  //* psijTilde_xi[nnl[j + i * MAX_NUM_INTERACTIONS]+i*MAX_NUM_INTERACTIONS][alpha];
                                  * psijTilde_xi[j + i * MAX_NUM_INTERACTIONS][alpha];

                //if(i == 86){
                //    Logger(DEBUG) << "psijTilde_xi[" << alpha << "]@" << nnl[j + i * MAX_NUM_INTERACTIONS]
                //              << " = " << psijTilde_xi[nnl[j + i * MAX_NUM_INTERACTIONS]][alpha]
                //              << ", f[" << nnl[j + i * MAX_NUM_INTERACTIONS] << "] = " << f[nnl[j + i * MAX_NUM_INTERACTIONS]];
                //}
            }
        }

        for (int j = 0; j < noiGhosts[i]; ++j) {
            for (int alpha = 0; alpha < DIM; ++alpha) {
                grad[i][alpha] += (fGhost[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]] - f[i])
                                  //* psijTilde_xiGhosts[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]+i*MAX_NUM_GHOST_INTERACTIONS][alpha];
                                  * psijTilde_xiGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS][alpha];
                // if(i == 86){
                //    Logger(DEBUG) << "psijTildeGhost_xi[" << alpha << "]@" << nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]
                //              << " = " << ghostParticles.psijTilde_xi[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]][alpha]
                //              << ", fGhost[" << nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS] << "] = " << fGhost[nnlGhosts[j + i * MAX_NUM_GHOST_INTERACTIONS]];
                //}
                            }
        }

        //for (int alpha=0; alpha < DIM; ++alpha ){
        //    if(i == 86) Logger(DEBUG) << "grad[" << i << "][" << alpha << "] = " << grad[i][alpha];
        //}
    }
}

void Particles::compEffectiveFace(const Particles &ghostParticles){
    for (int i=0; i<N; ++i){
        for (int j=0; j<noiGhosts[i]; ++j){
            int ji = nnlGhosts[i*MAX_NUM_GHOST_INTERACTIONS+j]; // index i of particle j
            // search neighbor i in nnlGhosts[] of j
            int ij = 0;
            for(ij=0; ij<noiGhosts[ji]; ++ij){
                for(int iGhost=i*(DIM+1); iGhost<i*(DIM+1)+3; ++iGhost){
                    if (nnlGhosts[ij+ji*MAX_NUM_GHOST_INTERACTIONS] == ghostMap[iGhost]) break;
                }
            }
            for (int alpha=0; alpha<DIM; ++alpha){
                AijGhosts[i*MAX_NUM_GHOST_INTERACTIONS+j][alpha] = 1./omega[i]*psijTilde_xiGhosts[i*MAX_NUM_GHOST_INTERACTIONS+j][alpha]
                                                       - 1./ghostParticles.omega[ji]*ghostParticles.psijTilde_xiGhosts[ij+ji*MAX_NUM_GHOST_INTERACTIONS][alpha];
            }
        }
    }
}

void Particles::compRiemannFluxes(const double &dt, const double &kernelSize, const double &gamma,
                                  const Particles &ghostParticles){
    for (int i=0; i<N; ++i){
        double xij[DIM];
        //double vFrame[DIM];
        // helper vectors
        double xijxi[DIM], xjxi[DIM], xijxj[DIM];

        for (int jn=0; jn<noiGhosts[i]; ++jn){

            int j = nnlGhosts[i*MAX_NUM_GHOST_INTERACTIONS+jn];

            xjxi[0] = ghostParticles.x[j] - x[i];
            xjxi[1] = ghostParticles.y[j] - y[i];

            xij[0] = x[i] + kernelSize/2. * (xjxi[0]);
            xij[1] = y[i] + kernelSize/2. * (xjxi[1]);

            xijxi[0] = xij[0] - x[i];
            xijxi[1] = xij[1] - y[i];

            xijxj[0] = xij[0] - ghostParticles.x[j];
            xijxj[1] = xij[1] - ghostParticles.y[j];
#if DIM==3
            xjxi[2] = ghostParticles.z[j] - z[i];
            xij[2] = z[i] + kernelSize/2. * (xjxi[2]);
            xijxi[2] = xij[2] - z[i];
            xijxj[2] = xij[2] - ghostParticles.z[j];
#endif
            double dotProd = Helper::dotProduct(xijxi, xjxi);
            double dSqr = Helper::dotProduct(xjxi, xjxi);

            int iW = i*MAX_NUM_GHOST_INTERACTIONS+jn;

            vFrameGhosts[iW][0] = vx[i] + (ghostParticles.vx[j]-vx[i]) * dotProd/dSqr;
            vFrameGhosts[iW][1] = vy[i] + (ghostParticles.vy[j]-vy[i]) * dotProd/dSqr;
#if DIM==3
            vFrameGhosts[iW][2] = vz[i] + (ghostParticles.vz[j]-vz[i]) * dotProd/dSqr;
#endif

            /*if(i==4){
                Logger(DEBUG) << "vFrame[0] = " << vFrame[0]
                              << ", vFrame[1] = " << vFrame[1]
                              << ", rhoGrad[i][0] = " << rhoGrad[i][0]
                              << ", rhoGrad[i][1] = " << rhoGrad[i][1]
                              << ", rho[i] = " << rho[i]
                              << ", xijxi = " << xijxi[0] << ", " << xijxi[1]
                              << ", xij = " << xij[0] << ", " << xij[1]
                              << ", xj = " << ghostParticles.x[j] << ", " << ghostParticles.y[j] << " @" << j;
                exit(5);
            }*/
            // boost frame to effective face
            WijRGhosts[iW][0] = rho[i];
            WijLGhosts[iW][0] = ghostParticles.rho[j];
            WijRGhosts[iW][1] = P[i];
            WijLGhosts[iW][1] = ghostParticles.P[j];
            WijRGhosts[iW][2] = vx[i] - vFrameGhosts[iW][0];
            WijLGhosts[iW][2] = ghostParticles.vx[j] - vFrameGhosts[iW][0];
            WijRGhosts[iW][3] = vy[i] - vFrameGhosts[iW][1];
            WijLGhosts[iW][3] = ghostParticles.vy[j] - vFrameGhosts[iW][1];
#if DIM == 3
            WijRGhosts[iW][4] = vz[i] - vFrameGhosts[iW][2];
            WijLGhosts[iW][4] = ghostParticles.vz[j] - vFrameGhosts[iW][2];
#endif

            // reconstruction at effective face
            WijRGhosts[iW][0] += Helper::dotProduct(rhoGrad[i], xijxi);
            WijLGhosts[iW][0] += Helper::dotProduct(ghostParticles.rhoGrad[j], xijxj);
            WijRGhosts[iW][1] += Helper::dotProduct(PGrad[i], xijxi);
            WijLGhosts[iW][1] += Helper::dotProduct(ghostParticles.PGrad[j], xijxj);
            WijRGhosts[iW][2] += Helper::dotProduct(vxGrad[i], xijxi);
            WijLGhosts[iW][2] += Helper::dotProduct(ghostParticles.vxGrad[j], xijxj);
            WijRGhosts[iW][3] += Helper::dotProduct(vyGrad[i], xijxi);
            WijLGhosts[iW][3] += Helper::dotProduct(ghostParticles.vyGrad[j], xijxj);
#if DIM==3
            WijRGhosts[iW][4] += Helper::dotProduct(vzGrad[i], xijxi);
            WijLGhosts[iW][4] += Helper::dotProduct(ghostParticles.vzGrad[j], xijxj);
#endif
            // predict half a timestep
            double viDiv = vxGrad[i][0] + vyGrad[i][1];
            double vjDiv = ghostParticles.vxGrad[j][0] + ghostParticles.vyGrad[j][1];
#if DIM==3
            viDiv += vzGrad[i][2];
            vjDiv += ghostParticles.vzGrad[j][2];
#endif
            WijRGhosts[iW][0] -= dt/2. * rho[i] * viDiv;
            WijLGhosts[iW][0] -= dt/2. * ghostParticles.rho[j] * vjDiv;
            WijRGhosts[iW][1] -= dt/2. * gamma*P[i] * viDiv;
            WijLGhosts[iW][1] -= dt/2. * gamma*ghostParticles.P[j] * vjDiv;
            WijRGhosts[iW][2] -= dt/2. * PGrad[i][0]/rho[i];
            WijLGhosts[iW][2] -= dt/2. * ghostParticles.PGrad[j][0]/rho[j];
            WijRGhosts[iW][3] -= dt/2. * PGrad[i][1]/rho[i];
            WijLGhosts[iW][3] -= dt/2. * ghostParticles.PGrad[j][1]/rho[j];
#if DIM==3
            WijRGhosts[iW][4] -= dt/2. * PGrad[i][2]/rho[i];
            WijLGhosts[iW][4] -= dt/2. * ghostParticles.PGrad[j][2]/rho[j];
#endif
        }
    }
}

/// debugging function dumping nearest neighbors to file
void Particles::dumpNNL(std::string filename, const Particles &ghostParticles){
    // open output file
    HighFive::File h5File { filename, HighFive::File::ReadWrite |
                                      HighFive::File::Create |
                                      HighFive::File::Truncate };

    for (int i=0; i<N; ++i){

        int noiTot = noi[i] + noiGhosts[i];
        std::vector<std::vector<double>> nnlPrtcls {};

        std::vector<size_t> dataSpaceDims(2);
        dataSpaceDims[0] = std::size_t(noiTot);
        dataSpaceDims[1] = DIM;

        for (int j=0; j<noi[i]; ++j){
            nnlPrtcls.push_back(std::vector<double>(DIM));
            nnlPrtcls[j][0] = x[nnl[i*MAX_NUM_INTERACTIONS+j]];
            nnlPrtcls[j][1] = y[nnl[i*MAX_NUM_INTERACTIONS+j]];
#if DIM == 3
            nnlPrtcls[j][2] = z[nnl[i*MAX_NUM_INTERACTIONS+j]];
#endif
        }

        for (int j=0; j<noiGhosts[i]; ++j){
            nnlPrtcls.push_back(std::vector<double>(DIM));
            nnlPrtcls[j+noi[i]][0] = ghostParticles.x[nnlGhosts[i*MAX_NUM_GHOST_INTERACTIONS+j]];
            nnlPrtcls[j+noi[i]][1] = ghostParticles.y[nnlGhosts[i*MAX_NUM_GHOST_INTERACTIONS+j]];
#if DIM == 3
            nnlPrtcls[j+noi[i]][2] = ghostParticles.z[nnlGhosts[i*MAX_NUM_GHOST_INTERACTIONS+j]];
#endif
        }

        HighFive::DataSet nnlDataSet = h5File.createDataSet<double>("/nnlPrtcls" + std::to_string(i),
                                                                    HighFive::DataSpace(dataSpaceDims));
        nnlDataSet.write(nnlPrtcls);
    }
}


#endif

void Particles::move(const double &dt, Domain &domain){

    for(int i=0; i<N; ++i) {

        //Logger(DEBUG) << "v@" << i <<  " = [" << vx[i] << ", " << vy[i] << "]"
        //            << ", x[n] = [" << x[i] << ", " << y[i] << "]";

        x[i] = x[i] + vx[i] * dt;
        y[i] = y[i] + vy[i] * dt;
#if DIM == 3
        z[i] = z[i] +vz[i] * dt;
#endif
#if PERIODIC_BOUNDARIES
        if (x[i] < domain.bounds.minX) {
            x[i] = domain.bounds.maxX - (domain.bounds.minX - x[i]);
        } else if (domain.bounds.maxX <= x[i]) {
            x[i] = domain.bounds.minX + (x[i] - domain.bounds.maxX);
        }
        if (y[i] < domain.bounds.minY) {
            y[i] = domain.bounds.maxY - (domain.bounds.minY - y[i]);
        } else if (domain.bounds.maxY <= y[i]) {
            y[i] = domain.bounds.minY + (y[i] - domain.bounds.maxY);
        }
#if DIM ==3
        if (z[i] < domain.bounds.minZ) {
            z[i] = domain.bounds.maxZ - (domain.bounds.minZ - z[i]);
        } else if (domain.bounds.maxZ <= z[i]) {
            z[i] = domain.bounds.minZ + (z[i] - domain.bounds.maxZ);
        }
#endif
#endif
        //Logger(DEBUG) << "                               x[n+1] = ["
        //          << x[i] << ", " << y[i] << "]";
    }
}

/// Sanity check functions
double Particles::sumVolume(){
    double V = 0.;
    for (int i=0; i<N; ++i){
        V += 1./omega[i];
    }
    return V;
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
    HighFive::DataSet rhoGradDataSet = h5File.createDataSet<double>("/rhoGrad", HighFive::DataSpace(dataSpaceDims));

    // containers for particle data
    std::vector<double> rhoVec(rho, rho+N);
    std::vector<std::vector<double>> posVec(N);
    std::vector<std::vector<double>> rhoGradVec(N);

    // fill containers with data
    std::vector<double> posBuf(DIM);
    std::vector<double> rhoGradBuf(DIM);
    for(int i=0; i<N; ++i){
        //Logger(DEBUG) << "      > Dumping particle @"  << i;
        // position
        posBuf[0] = x[i];
        posBuf[1] = y[i];
#if DIM == 3
        posBuf[2] = z[i];
#endif
        posVec[i] = posBuf;

        // density gradient
        rhoGradBuf[0] = rhoGrad[i][0];
        rhoGradBuf[1] = rhoGrad[i][1];
#if DIM == 3
        rhoGradBuf[2] = rhoGrad[i][2];
#endif
        rhoGradVec[i] = rhoGradBuf;
    }
    // write data
    rhoDataSet.write(rhoVec);
    posDataSet.write(posVec);
    rhoGradDataSet.write(rhoGradVec);
}

