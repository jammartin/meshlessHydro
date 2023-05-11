//
// Created by Johannes Martin on 05.09.22.
//

#include "../include/Particles.h"

double Kernel::cubicSpline(const double &r, const double &h) {

    // TODO: remove this
    double h2 = h/2.;
#if DIM == 2
    const double sigma = 10./(7.*M_PI*h2*h2);
#else // DIM == 3
    const double sigma = 1./(M_PI*h2*h2*h2);
#endif
    const double q = r/h2;
    if (0. <= q && q <= 1.){
        return sigma*(1.-3./2.*q*q*(1.-q/2.));
    } else if (1. < q && q < 2.){
        return sigma/4.*pow(2.-q, 3.);
    } else {
        return 0.;
    }
}

// // dW(r, h)/dr. Scalar. For needs to be multiplied with (x_i - x_j)/r
// // For 2d SPH only
// double Kernel::dWdr(const double &r, const double &h){
//     const double sigma = 10./(7.*M_PI*h*h*h);
//     const double q = r/h;
//     if (0. <= q && q <= 1./2.){
//         return sigma * (- 3 * q + 9./4. * q * q);
//     } else if (1./2. < q && q < 1.){
//         return sigma * -1 * 0.75 * pow((2 - q), 2);
//     } else {
//         return 0.;
//     }
// };

// dW(r, h)/dr. Scalar. For needs to be multiplied with (x_i - x_j)/r
// For 2d SPH only
double Kernel::dWdr(const double &r, const double &h){
    const double sigma = 40./(7.*M_PI);
    const double q = r/h;
    if (0. <= q && q < 1./2.){
        return 6 * sigma / pow(h, DIM + 1) * (3*pow(q, 2) - 2 * q);
    } else if (1./2. <= q && q <= 1.){
        return 6 * sigma / pow(h, DIM + 1) * -1 * pow((1 - q), 2);
    } else {
        return 0.;
    }
};

// dW/dh
double Kernel::dWdh(const double &r, const double &h){
    const double sigma = 40./(7.*M_PI*h*h*h);
    const double q = r/h;
    if (0. <= q && q <= 1./2.){
        return 2 * sigma / pow(h, 6) * (12 * h * r*r - pow(h, 3) - 15 * pow(r, 3));
    } else if (1./2. < q && q < 1.){
        return 2 * sigma / pow(h, 6) * pow((h - r), 2) * (5 * r - 2 * h);
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
    PGrad = new double[numParticles][DIM];


#if RUNSPH
    ax = new double[numParticles];
    ay = new double[numParticles];

    // ArtVisc thingy
    cs = new double[numParticles];
    dudtArtVisc = new double[numParticles];

    axArtVisc = new double[numParticles];
    ayArtVisc = new double[numParticles];
#if DIM == 3
    az = new double[numParticles];
    azArtVisc = new double[numParticles];
#endif


	//unimportant = new double[numParticles];
	dEdt = new double[numParticles];
	dn = new double[numParticles];
	drho = new double[numParticles];
#endif // RUNSPH

    //TODO: check if this is needed as array
    //B = new double[numParticles][DIM*DIM];
#if DIM == 3
    z = new double[numParticles];
    vz = new double[numParticles];
    vzGrad = new double[numParticles][DIM];
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
        psijTilde_xiGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM];
        AijGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM];
        WijLGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM+2];
        WijRGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM+2];
        FijGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM+2]; // TODO: this buffer should not be needed
        vFrameGhosts = new double[numParticles*MAX_NUM_GHOST_INTERACTIONS][DIM];
    } else {
        parent = new int[numParticles]; // store index of original node
#endif
    }
#if DEBUG_LVL > 1
    if (ghosts){
        // This is necessary for dumping ghosts to file
        noi = new int[numParticles];
    }
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
    delete[] PGrad;
    delete[] omega;

#if RUNSPH
	//delete[] unimportant;
	delete[] dEdt;
	delete[] dn;
	delete[] drho;

    delete[] ax;
    delete[] ay;

    delete[] cs;
    delete[] dudtArtVisc;

    delete[] axArtVisc;
    delete[] ayArtVisc;
#if DIM == 3
    delete[] az;
    delete[] azArtVisc;
#endif
#endif // RUNSPH

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
        delete[] vzGrad;
#endif
        delete[] nnl;
        delete[] noi;
        delete[] Aij;
#if PERIODIC_BOUNDARIES
        delete[] nnlGhosts;
        delete[] noiGhosts;
        delete[] psijTilde_xiGhosts;
        delete[] AijGhosts;
        delete[] WijLGhosts;
        delete[] WijRGhosts;
        delete[] ghostMap;
        delete[] FijGhosts;
        delete[] vFrameGhosts;
    } else {
        delete[] parent;
#endif
    }
#if DEBUG_LVL > 1
    if (ghosts){
        // This is necessary for dumping ghosts to file
        delete[] noi;
    }
#endif
}

#if !PERIODIC_BOUNDARIES
void Particles::getDomainLimits(double *domainLimits){

    double minX { std::numeric_limits<double>::max() };
    double maxX { std::numeric_limits<double>::min() };
    double minY { std::numeric_limits<double>::max() };
    double maxY { std::numeric_limits<double>::min() };
#if DIM == 3
    double minZ { std::numeric_limits<double>::max() };
    double maxZ { std::numeric_limits<double>::min() };
#endif

    for(int i=0; i<N; ++i){
        if (x[i] < minX){
            minX = x[i];
        } else if (x[i] > maxX){
            maxX = x[i];
        }
        if (y[i] < minY){
            minY = y[i];
        } else if (y[i] > maxY){
            maxY = y[i];
        }
#if DIM == 3
        if (z[i] < minZ){
            minZ = z[i];
        } else if (z[i] > maxZ){
            maxZ = z[i];
        }
#endif
    }

    domainLimits[0] = minX;
    domainLimits[DIM] = maxX;
    domainLimits[1] = minY;
    domainLimits[DIM+1] = maxY;
#if DIM == 3
    domainLimits[2] = minZ;
    domainLimits[DIM+2] = maxZ;
#endif
}
#endif // !PERIODIC_BOUNDARIES

void Particles::assignParticlesAndCells(Domain &domain){

    // reset particles assigned to grid cells
    for(int iGrid=0; iGrid<domain.numGridCells; ++iGrid){
        domain.grid[iGrid].prtcls = std::vector<int>();
    }

    for(int i=0; i<N; ++i){

        int floorX = floor((x[i]-domain.bounds.minX)/domain.cellSizeX);
        int floorY = floor((y[i]-domain.bounds.minY)/domain.cellSizeY);

        // This rarely happens when x or y is really close to the domain bounds
        if (floorX == domain.cellsX){
            floorX -= 1;
        }
        if (floorY == domain.cellsY){
            floorY -= 1;
        }

#if DIM==3
        int floorZ = floor((z[i]-domain.bounds.minZ)/domain.cellSizeZ);
        if (floorZ == domain.cellsZ){
            floorZ -= 1;
        }
#endif


        int iGrid = floorX + floorY * domain.cellsX
#if DIM == 3
        + floorZ * domain.cellsX * domain.cellsY
#endif
        ;

        // if (floorX >= domain.cellsX || floorY >= domain.cellsY){
        //    Logger(ERROR) << "  > Particle i = " << i << " cannot be properly assigned to search grid.";
        //    Logger(ERROR) << "    > x = " << x[i] << " -> floorX = " << floorX
        //                  << ", y = " << y[i] << " -> floorY = " << floorY;
        // }
        //
        // Logger(DEBUG) << "floor x = " << floorX
        //          << ", floor y = " << floorY;

        //Logger(DEBUG) << "      > Assigning particle@" << i << " = [" << x[i] << ", " << y[i]
//#if DIM == 3
//                      << ", " << z[i]
//#endif
//                      << "] to cell " << iGrid;

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
        if (noiBuf == 0){
            Logger(WARN) << "No neighbors for particle " << i << ". Caution.";
        }
        noi[i] = noiBuf;
    }
}

#if RUNSPH
// For comparable ICs: This sets the internal energies so that P = 2.5 everywhere
void Particles::setInternalEnergy(double Pressure, const double gamma){
    for (int i = 0; i < N; i++){
        u[i] = Pressure / ((gamma - 1) * rho[i]);
    }
}

// Computes the density via kernel smoothing w/o ghost particles
void Particles::compDensitySPH(const double &kernelSize){
    double dnst;
    double dSqr;
    double r;
    int iP;
    for (int i = 0; i < N; i++){
        dnst = 0;
        dSqr = 0;
        r = 0.;
        for(int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            dSqr = pow(x[i] - x[iP], 2)
                        + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            r = sqrt(dSqr);
            dnst += m[i] * kernel(r, kernelSize);
            // For normalization, as in compOmega: Add self interaction
            rho[i] = dnst + m[j] * kernel(0., kernelSize);
        }
    }
}


void Particles::compAccSPH(const double &kernelSize){
    int iP;
    double dSqr;
    double r;
    double PRhoHost;
    double PRhoTarget;
    double PIij;
    for (int i = 0; i < N; i++){
        ax[i] = 0;
        ay[i] = 0;
        #if DIM == 3
        az[i] = 0;
        #endif
        PRhoHost = P[i] / pow(rho[i], 2);
        for (int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];

            dSqr = pow(x[i] - x[iP], 2)
            + pow(y[i] - y[iP], 2);
            #if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
            #endif
            r = sqrt(dSqr);
            PRhoTarget =  P[j] / pow(rho[j], 2);
            #if ARTVISC
            PIij = compPIij(i, iP, kernelSize);
            ax[i] += m[iP] * (PRhoHost + PRhoTarget + PIij)  * (x[iP] - x[i])/r * Kernel::dWdr(r, kernelSize);
            ay[i] += m[iP] * (PRhoHost + PRhoTarget + PIij)  * (y[iP] - y[i])/r * Kernel::dWdr(r, kernelSize);

            //ax[i] -= m[j] * (PRhoHost + PRhoTarget) * (x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
            //ay[i] -= m[j] * (PRhoHost + PRhoTarget) * (y[i] - y[j])/r * Kernel::dWdr(r, kernelSize);
            // Logger(DEBUG) << i << " " << j << " " << m[j]*(PRhoHost+PRhoTarget)*(x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
            #if DIM == 3
            az[i] -= m[j]*(PRhoHost+PRhoTarget+PIij)*(z[i] - z[j])/r * Kernel::dWdr(r, kernelSize);
            #endif
            #else
            ax[i] += m[iP] * (PRhoHost+PRhoTarget)  * (x[iP] - x[i])/r * Kernel::dWdr(r, kernelSize);
            ay[i] += m[iP] * (PRhoHost+PRhoTarget)  * (y[iP] - y[i])/r * Kernel::dWdr(r, kernelSize);

            //ax[i] -= m[j] * (PRhoHost + PRhoTarget) * (x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
            //ay[i] -= m[j] * (PRhoHost + PRhoTarget) * (y[i] - y[j])/r * Kernel::dWdr(r, kernelSize);
            // Logger(DEBUG) << i << " " << j << " " << m[j]*(PRhoHost+PRhoTarget)*(x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
            #if DIM == 3
            az[i] -= m[j]*(PRhoHost+PRhoTarget)*(z[i] - z[j])/r * Kernel::dWdr(r, kernelSize);
            #endif
            ax[i] -= m[i]*(PRhoHost+PRhoTarget)*(x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
            ay[i] -= m[i]*(PRhoHost+PRhoTarget)*(y[i] - y[j])/r * Kernel::dWdr(r, kernelSize);
            #if DIM == 3
            az[i] -= m[i]*(PRhoHost+PRhoTarget)*(z[i] - z[j])/r * Kernel::dWdr(r, kernelSize);
            #endif
            #endif
        }
    }
}


// Euler intgration for SPH
void Particles::eulerSPH(const double &dt, const Domain &domain){
    for (int i = 0; i < N; i++){
        vx[i] += ax[i]*dt;
        vy[i] += ay[i]*dt;
#if DIM == 3
        vz[i] += dt * az[i];
#endif

        x[i] += vx[i]*dt;
        y[i] += vy[i]*dt;
#if DIM == 3
        z[i] += vz[i]*dt;
#endif

#if PERIODIC_BOUNDARIES
        if (x[i] < domain.bounds.minX) {
            x[i] = domain.bounds.maxX - (domain.bounds.minX - x[i]);
        }
        else if (domain.bounds.maxX <= x[i]) {
            //Logger(DEBUG) << "X is " << x[i];
            x[i] = domain.bounds.minX + (x[i] - domain.bounds.maxX);
        }
        if (y[i] < domain.bounds.minY) {
            y[i] = domain.bounds.maxY - (domain.bounds.minY - y[i]);
        }
        else if (domain.bounds.maxY <= y[i]) {
            y[i] = domain.bounds.minY + (y[i] - domain.bounds.maxY);
        }
#if DIM ==3
        if (z[i] < domain.bounds.minZ) {
            z[i] = domain.bounds.maxZ - (domain.bounds.minZ - z[i]);
        }
        else if (domain.bounds.maxZ <= z[i]) {
            z[i] = domain.bounds.minZ + (z[i] - domain.bounds.maxZ);
        }
#endif
#endif

        // Check if particles are out of bounds
        // if (x[i] > domain.bounds.maxX){
            //     Logger(DEBUG) << "X is too big, out of bounds";
            // }
    }
}


// SPH energy computation functions:
void Particles::compuis(const double &dt, const double &kernelSize){
    //double u_i_old;
    for (int i = 0; i < N; i++){
        //u_i_old = u[i];
        u[i] += dEdt[i] * dt;
#if ARTVISC
        u[i] += dudtArtVisc[i] * dt;
#endif
        //std::cout << "u_" << i << ": " << u[i] << ", was " << u_i_old << std::endl;;
        if (u[i] < 0.){
            Logger(WARN) << "+++DANGER+++ negative specific internal energy encountered.";
        }
    }
}


// Computing all omegas:
void Particles::compOmegas(const double &kernelSize){
    for (int i = 0; i < N; i++){
        compOmega(i, kernelSize);
    }
}

// Implementing equations F5 and F6 in Hopkins` GIZMO Paper
void Particles::calcdndrho(const double &kernelSize){
    double r, dSqr, tmp;
    int iP;
    for (int i = 0; i < N; i++){
        dn[i] = 0;
        drho[i] = 0;
        for (int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            dSqr = pow(x[i] - x[iP], 2)
                        + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            r = sqrt(dSqr);
            dn[i] -= 1 / kernelSize * (DIM * Kernel::cubicSpline(r, kernelSize) + r*Kernel::dWdr(r, kernelSize));
            drho[i] -= m[iP]*1 / kernelSize * (DIM * Kernel::cubicSpline(r, kernelSize) + r*Kernel::dWdr(r, kernelSize));
        }
    }
}

// Implementing the sum in equation F3 in Hopkins` GIZMO Paper
void Particles::calcdE(const double &kernelSize){
    double tmp, dSqr, r, fij;
    int iP;
    for (int i = 0; i < N; i++){
        dEdt[i] = 0;
        for (int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            dSqr = pow(x[i] - x[iP], 2)
                        + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            r = sqrt(dSqr);
            fij = 1 - 1 / m[iP]
                    * kernelSize / (omega[i] * DIM) * drho[i]
                        / (1 + kernelSize /(omega[i] * DIM) * dn[i]);
            //fij = 1;
            dEdt[i]  += m[iP]
                        * ((vx[i]-vx[iP])
                        * (x[i]-x[iP])
                        + (vy[i]-vy[iP])
                        * (y[i]-y[iP])
#if DIM == 3
                        + (vz[i]-vz[iP])*(z[i]-z[iP])
#endif
                        ) / r * P[i]/pow(rho[i],2)*fij*Kernel::dWdr(r, kernelSize)
                        ;
        }
    }
}

// For HLLC solver: Calculate norm vector between two particles:
// n_unit needs to be pre-allocated
void calcNunit(const int i, const int j, double* n_unit){
    double dSqr = pow(x[j] - x[i], 2)
                        + pow(y[j] - y[i], 2);
#if DIM == 3
    dSqr += pow(z[j] - z[i], 2);
#endif
    double r = sqrt(dSqr);
    n_unit[0] = (x[j] - x[i]) / r;
    n_unit[1] = (y[j] - y[i]) / r;
#if DIM == 3
    n_unit[2] = (z[j] - z[i]) / r;
#endif
}

#if PERIODIC_BOUNDARIES

// For HLLC solver: Calculate norm vector between two particles, one of which is a ghost:
// n_unit needs to be pre-allocated
void calcNunit(const Particles &ghostParticles, const int i, const int j, double* n_unit){
    double dSqr = pow(ghostParticles.x[j] - x[i], 2)
                        + pow(ghostParticles.y[j] - y[i], 2);
#if DIM == 3
    dSqr += pow(ghostParticles.z[j] - z[i], 2);
#endif
    double r = sqrt(dSqr);
    n_unit[0] = (ghostParticles.x[j] - x[i]) / r;
    n_unit[1] = (ghostParticles.y[j] - y[i]) / r;
#if DIM == 3
    n_unit[2] = (ghostParticles.z[j] - z[i]) / r;
#endif
}

// Computes the density via kernel smoothing w/ ghost particles
void Particles::compDensitySPH(const Particles &ghostParticles, const double &kernelSize){
    int iP;
    double dnst;
    double dSqr;
    double r;
    double ghostMass;
    for (int i = 0; i < N; ++i){
        dnst = 0;
        for(int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            dSqr = pow((x[i] - x[iP]), 2)
                        + pow((y[i] - y[iP]), 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            r = sqrt(dSqr);

            //Logger(DEBUG) << "density: i = " << i << ", j = " << iP << ", r = " << r;

            dnst += m[iP] * kernel(r, kernelSize);
        }
        // For normalization: Self interaction
        dnst += m[i] * kernel(0., kernelSize);

        // Ghosts
        //Logger(DEBUG) << "i: " << i << " noiGhosts: " << noiGhosts[i];
        //Logger(DEBUG) << "Random ghost mass: " << ghostParticles.m[0];
#if PERIODIC_BOUNDARIES
        for(int k = 0; k < noiGhosts[i]; ++k){
            iP = nnlGhosts[k+i*MAX_NUM_GHOST_INTERACTIONS];
            dSqr = pow(x[i] - ghostParticles.x[iP], 2)
                          + pow(y[i] - ghostParticles.y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - ghostParticles.z[iP], 2);
#endif
            r = sqrt(dSqr);
            // TODO: include mass in update ghost state
            dnst += m[ghostParticles.parent[iP]] * kernel(r, kernelSize);
            //Logger(DEBUG) << "k = " << k << " mass:  " << m[ghostParticles.parent[iP]] << " W: " << kernel(r,kernelSize);
        }
#endif
        // if ((i - 0) % 30 == 0){
        //     Logger(DEBUG) << "density from ghosts: i = " << i << ", dnst = " << dnst;
        // }
        rho[i] = dnst;
    }
}


void Particles::compAccSPH(const Particles &ghostParticles, const double &kernelSize){
    int iP;
    // Calculate acceleration for each particle, w/o Ghosts
    // c.f. eq 8 in Monaghan: SPH and its diverse Applications, Annu.Rev. Fluid mechanics, 2012
    double dSqr;
    double r;
    double PRhoHost;
    double PRhoTarget;
    double PIij;
    for (int i = 0; i < N; i++){
        ax[i] = 0;
        ay[i] = 0;
#if DIM == 3
        az[i] = 0;
#endif
        PRhoHost = P[i] / pow(rho[i], 2);
        for (int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];

            dSqr = pow(x[i] - x[iP], 2)
                        + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            r = sqrt(dSqr);
            PRhoTarget =  P[iP] / pow(rho[iP], 2);
#if ARTVISC
            PIij = compPIij(i, iP, kernelSize); // TODO: what is 1.5 and 3. ???
            ax[i] += m[iP] * (PRhoHost + PRhoTarget + PIij)  * (x[iP] - x[i])/r * Kernel::dWdr(r, kernelSize);
            ay[i] += m[iP] * (PRhoHost + PRhoTarget + PIij)  * (y[iP] - y[i])/r * Kernel::dWdr(r, kernelSize);

            //ax[i] -= m[j] * (PRhoHost + PRhoTarget) * (x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
            //ay[i] -= m[j] * (PRhoHost + PRhoTarget) * (y[i] - y[j])/r * Kernel::dWdr(r, kernelSize);
            // Logger(DEBUG) << i << " " << j << " " << m[j]*(PRhoHost+PRhoTarget)*(x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
#if DIM == 3
            az[i] -= m[j]*(PRhoHost+PRhoTarget+PIij)*(z[i] - z[j])/r * Kernel::dWdr(r, kernelSize);
#endif
#else
            ax[i] += m[iP] * (PRhoHost+PRhoTarget)  * (x[iP] - x[i])/r * Kernel::dWdr(r, kernelSize);
            ay[i] += m[iP] * (PRhoHost+PRhoTarget)  * (y[iP] - y[i])/r * Kernel::dWdr(r, kernelSize);

//ax[i] -= m[j] * (PRhoHost + PRhoTarget) * (x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
//ay[i] -= m[j] * (PRhoHost + PRhoTarget) * (y[i] - y[j])/r * Kernel::dWdr(r, kernelSize);
// Logger(DEBUG) << i << " " << j << " " << m[j]*(PRhoHost+PRhoTarget)*(x[i] - x[j])/r * Kernel::dWdr(r, kernelSize);
#if DIM == 3
            az[i] -= m[j]*(PRhoHost+PRhoTarget)*(z[i] - z[j])/r * Kernel::dWdr(r, kernelSize);
#endif
#endif
        }

        // Ghosts
        for (int k = 0; k < noiGhosts[i]; k++){
            iP = nnlGhosts[k+i*MAX_NUM_GHOST_INTERACTIONS];
            dSqr = pow(x[i] - ghostParticles.x[iP], 2)
                          + pow(y[i] - ghostParticles.y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - ghostParticles.z[iP], 2);
#endif
            r = sqrt(dSqr);
            PRhoTarget = ghostParticles.P[iP]
                / pow(ghostParticles.rho[iP], 2);
#if ArtVisc
            PIij = compPIij(ghostParticles, i, iP, kernelSize);
            ax[i] -= m[ghostParticles.parent[iP]]*(PRhoHost+PRhoTarget+PIij)
                * (x[i] - ghostParticles.x[iP])
                    / r * Kernel::dWdr(r, kernelSize);
            ay[i] -= m[ghostParticles.parent[iP]]*(PRhoHost+PRhoTarget+PIij)
                * (y[i] - ghostParticles.y[iP])
                    / r * Kernel::dWdr(r, kernelSize);
#if DIM == 3
            az[i] -= m[ghostParticles.parent[iP]]*(PRhoHost+PRhoTarget+PIij)
                * (z[i] - ghostParticles.z[iP]) / r * Kernel::dWdr(r, kernelSize);
#endif
#else
            ax[i] -= m[ghostParticles.parent[iP]]*(PRhoHost+PRhoTarget)
                * (x[i] - ghostParticles.x[iP])
                    / r * Kernel::dWdr(r, kernelSize);
            ay[i] -= m[ghostParticles.parent[iP]]*(PRhoHost+PRhoTarget)
                * (y[i] - ghostParticles.y[iP])
                    / r * Kernel::dWdr(r, kernelSize);
#if DIM == 3
            az[i] -= m[ghostParticles.parent[iP]]*(PRhoHost+PRhoTarget)
                * (z[i] - ghostParticles.z[iP]) / r * Kernel::dWdr(r, kernelSize);
#endif
#endif
        }
        // Logger(DEBUG) << "For i = " << i << ", ax and ay are " << ax[i] << " and " << ay[i];
    }
}

// Compute all omegas
void Particles::compOmegas(const Particles &ghostParticles, const double &kernelSize){
    for (int i = 0; i < N; i++){
        compOmega(i, ghostParticles, kernelSize);
        //std::cout << i << " " << omega[i] << std::endl;
    }
}

// Implementing equations F5 and F6 in Hopkins` GIZMO Paper
void Particles::calcdndrho(const Particles &ghostParticles, const double &kernelSize){
    int iP;
    double r, dSqr, tmp;
    for (int i = 0; i < N; i++){
        dn[i] = 0;
        drho[i] = 0;
        for (int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            dSqr = pow(x[i] - x[iP], 2)
                        + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            r = sqrt(dSqr);

            dn[i] -= 1 / kernelSize
                * (DIM *Kernel::cubicSpline(r, kernelSize) + r / kernelSize * Kernel::dWdh(r, kernelSize));
            drho[i] -= m[iP] / kernelSize
                * (DIM *Kernel::cubicSpline(r, kernelSize) + r / kernelSize * Kernel::dWdh(r, kernelSize));
        }
        for (int k = 0; k < noiGhosts[i]; k++){
            iP = nnl[k+i*MAX_NUM_GHOST_INTERACTIONS];
            dSqr = pow(x[i] - ghostParticles.x[iP], 2)
            + pow(y[i] - ghostParticles.y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - ghostParticles.z[iP], 2);
#endif
            r = sqrt(dSqr);
            dn[i] -= 1 / kernelSize
                * (DIM * Kernel::cubicSpline(r, kernelSize) + r/kernelSize * Kernel::dWdh(r, kernelSize));
            drho[i] -= m[ghostParticles.parent[iP]] / kernelSize
                * (DIM * Kernel::cubicSpline(r, kernelSize) + r/kernelSize * Kernel::dWdh(r, kernelSize));
        }
        //std::cout << "For i = " << i <<", drho is " << drho[i] << " dn " << dn[i] << " m " << m[i] << std::endl;
    }
}


// Implementing equation F3 in Hopkins` GIZMO Paper
void Particles::calcdE(const Particles &ghostParticles, const double &kernelSize){
    int iP;
    double dSqr, r, fij;
    for (int i = 0; i < N; i++){
        dEdt[i] = 0;
        for (int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            dSqr = pow(x[i] - x[iP], 2)
                        + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            r = sqrt(dSqr);
            fij = 1 - 1 / m[iP]
                    * kernelSize / (omega[i] * DIM) * drho[i]
                        / (1 + kernelSize /(omega[i] * DIM) * dn[i]);
            //fij = 1;
            dEdt[i]  += m[iP]
                        * ((vx[i]-vx[iP])
                        * (x[i]-x[iP])
                        + (vy[i]-vy[iP])
                        * (y[i]-y[iP])
#if DIM == 3
                        + (vz[i]-vz[iP])*(z[i]-z[iP])
#endif
                        ) / r * P[i]/pow(rho[i],2)*fij*Kernel::dWdr(r, kernelSize)
                        ;
        }
        // std::cout << "For i = " << i <<", fij is " << fij << std::endl;
        for (int k = 0; k < noiGhosts[i]; k++){
            iP = nnlGhosts[k+i*MAX_NUM_GHOST_INTERACTIONS];
            dSqr = pow(x[i] - ghostParticles.x[iP], 2)
            + pow(y[i] - ghostParticles.y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - ghostParticles.z[iP], 2);
#endif
            r = sqrt(dSqr);

            fij = 1 - 1 / m[ghostParticles.parent[iP]]
                    * kernelSize / omega[i] / DIM * drho[i] / (1 + kernelSize /omega[i] / DIM * dn[i]);
            //fij = 1;
            dEdt[i]  += m[ghostParticles.parent[iP]]
                        * ((vx[i]-ghostParticles.vx[iP])* (x[i]-ghostParticles.x[iP])
                        + (vy[i]-ghostParticles.vy[iP]) * (y[i]-ghostParticles.y[iP])
#if DIM == 3
                        + (vz[i]-vz[nnl[k+i*MAX_NUM_GHOST_INTERACTIONS]])*(z[i]-z[nnl[k+i*MAX_NUM_GHOST_INTERACTIONS]])
#endif
                        ) / r * P[i]/pow(rho[i],2)*fij*Kernel::dWdr(r, kernelSize)
                    ;
        }
        //std::cout << "dE for i=" << i << " is " << dEdt[i] << std::endl;
    }
}

#endif

#if ARTVISC
// Artificial Viscocity!
// Compute cs
void Particles::compCs(const double gamma){
    for (int i = 0; i < N; i++){
        cs[i] = sqrt((gamma - 1) * u[i]);
    }
}

// Compute Mu_ij
double Particles::compMuij(int i, int j, const double &kernelSize){
    double numerator;
    numerator = ((vx[i] - vx[j]) * (x[i] - x[j])
            + (vy[i] - vy[j]) * (y[i] - y[j]));
#if DIM == 3
    numerator += (vz[i] - vz[j]) * (z[i] - z[j]);
#endif

    double rSqr;
    rSqr = pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2);
#if DIM == 3
    rSqr += pow(z[i] - z[j], 2);
#endif

    return kernelSize * numerator / (sqrt(rSqr) + EPSMU * pow(kernelSize, 2));
}

// Compute PI_ij:
double Particles::compPIij(int i, int j, const double &kernelSize){
    double numerator;
    numerator = ((vx[i] - vx[j]) * (x[i] - x[j])
    + (vy[i] - vy[j]) * (y[i] - y[j]));
    #if DIM == 3
    numerator += (vz[i] - vz[j]) * (z[i] - z[j]);
    #endif

    if (numerator < 0){
        return 0;
    }
    else{
        double Muij = compMuij(i, j, kernelSize);

        double cij = (cs[i] + cs[j])/2;
        double rhoij = (rho[i] + rho[j]) / 2;

        return (- ALPHA_VISC * cij * Muij + BETA_VISC * pow(Muij, 2)) / rhoij;
    }
}

// Compute additional acceleration term for artificial viscocity
void Particles::compAccArtVisc(const double &kernelSize){
    double PIij, dSqr, r;
    // Interaction partner:
    int iP;
    for (int i = 0; i < N; i++){
        axArtVisc[i] = 0;
        ayArtVisc[i] = 0;
        #if DIM == 3
        azArtVisc[i] = 0;
        #endif
        for (int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];

            dSqr = pow(x[i] - x[iP], 2)
            + pow(y[i] - y[iP], 2);
            #if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
            #endif
            r = sqrt(dSqr);
            PIij = compPIij(i, iP, kernelSize);
            //Logger(DEBUG) << "PIij is " << PIij;
            axArtVisc[i] -= m[iP] * PIij * (x[i] - x[iP]) / r *  Kernel::dWdr(r, kernelSize);
            axArtVisc[i] -= m[iP] * PIij * (y[i] - y[iP]) / r *  Kernel::dWdr(r, kernelSize);
            #if DIM == 3
            azArtVisc[i] -= m[iP] * PIij * (z[i] - z[iP]) / r *  Kernel::dWdr(r, kernelSize);
            #endif
        }
    }
}

// Compute additional energy terms for artificial viscocity
void Particles::compUiArtVisc(const double &kernelSize){
    double PIij, dSqr, r;
    // Interaction partner:
    int iP;
    for (int i = 0; i < N; i++){
        dudtArtVisc[i] = 0;
        for (int j = 0; j < noi[i]; j++){
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];

            dSqr = pow(x[i] - x[iP], 2)
            + pow(y[i] - y[iP], 2);
            #if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
            #endif
            r = sqrt(dSqr);
            if (r <= 0){
                Logger(DEBUG) << "DANGER";
            }
            PIij = compPIij(i, iP, kernelSize);
            dudtArtVisc[i] += .5 * m[iP] * PIij *  Kernel::dWdr(r, kernelSize)
            * ((vx[i] - vx[iP]) * (x[i] - x[iP])
            + (vy[i] - vy[iP]) * (y[i] - y[iP])
            #if DIM == 3
            + (vz[i] - vy[iP]) * (z[i] - z[iP])
            #endif
        ) / r;
    }
}
}

#if PERIODIC_BOUNDARIES
// Compute Mu_ij
double Particles::compMuij(const Particles &ghostParticles, int i, int j, const double &kernelSize){
    double numerator;
    numerator = ((vx[i] - ghostParticles.vx[j]) * (x[i] - ghostParticles.x[j])
            + (vy[i] - ghostParticles.vy[j]) * (y[i] - ghostParticles.y[j]));
#if DIM == 3
    numerator += (vz[i] - ghostParticles.vz[j]) * (z[i] - ghostParticles.z[j]);
#endif

    double rSqr;
    rSqr = pow(x[i] - ghostParticles.x[j], 2) + pow(y[i] - ghostParticles.y[j], 2);
#if DIM == 3
    rSqr += pow(z[i] - ghostParticles.z[j], 2);
#endif
    return (kernelSize * numerator / (sqrt(rSqr) + 0.000025 * pow(kernelSize, 2)));
}


// Compute PI_ij:
double Particles::compPIij(const Particles &ghostParticles, int i, int j, const double &kernelSize){
    double numerator;
    numerator = ((vx[i] - ghostParticles.vx[j]) * (x[i] - ghostParticles.x[j])
            + (vy[i] - ghostParticles.vy[j]) * (y[i] - ghostParticles.y[j]));
#if DIM == 3
    numerator += (vz[i] - ghostParticles.vz[j]) * (z[i] - ghostParticles.z[j]);
#endif

    if (numerator < 0){
        return 0;
    }
    else{
        double Muij = compMuij(ghostParticles, i, j, kernelSize);

        double cij = (cs[i] + cs[ghostParticles.parent[j]])/2;
        double rhoij = (rho[i] + ghostParticles.rho[j]) / 2;

        return (- ALPHA_VISC * cij * Muij + BETA_VISC * pow(Muij, 2)) / rhoij;
    }
}

// Compute additional acceleration term for artificial viscocity w/ periodic boundaries
void Particles::compAccArtVisc(const Particles &ghostParticles, const double &kernelSize){
    double PIij, dSqr, r;
    // Interaction partner:
    int iP;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < noiGhosts[i]; j++){
            iP = nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS];

            dSqr = pow(x[i] - ghostParticles.x[iP], 2)
                        + pow(y[i] - ghostParticles.y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - ghostParticles.z[iP], 2);
#endif
            r = sqrt(dSqr);
            PIij = compPIij(ghostParticles, i, iP, kernelSize);
            //Logger(DEBUG) << "	> PIij ghost is " << PIij;
            axArtVisc[i] -= m[ghostParticles.parent[iP]] * PIij * (x[i] - ghostParticles.x[iP]) / r *  Kernel::dWdr(r, kernelSize);
            axArtVisc[i] -= m[ghostParticles.parent[iP]] * PIij * (y[i] - ghostParticles.y[iP]) / r *  Kernel::dWdr(r, kernelSize);
#if DIM == 3
            azArtVisc[i] -= m[ghostParticles.parent[iP]] * PIij * (z[i] - ghostParticles.z[iP]) / r *  Kernel::dWdr(r, kernelSize);
#endif
        }
    }
}



// Compute additional energy terms for artificial viscocity
void Particles::compUiArtVisc(const Particles &ghostParticles, const double &kernelSize){
    double PIij, dSqr, r;
    // Interaction partner:
    int iP;
    for (int i = 0; i < N; i++){
        for (int j = 0; j < noiGhosts[i]; j++){
            iP = nnlGhosts[j+i*MAX_NUM_GHOST_INTERACTIONS];

            dSqr = pow(x[i] - ghostParticles.x[iP], 2)
                        + pow(y[i] - ghostParticles.y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - ghostParticles.z[iP], 2);
#endif
            r = sqrt(dSqr);
            PIij = compPIij(ghostParticles, i, iP, kernelSize);
            dudtArtVisc[i] += .5 * m[ghostParticles.parent[iP]] * PIij *  Kernel::dWdr(r, kernelSize)
                * ((vx[i] - ghostParticles.vx[iP]) * (x[i] - ghostParticles.x[iP])
                    + (vy[i] - ghostParticles.vy[iP]) * (y[i] - ghostParticles.y[iP])
#if DIM == 3
                    + (vz[i] - ghostParticles.vy[iP]) * (z[i] - ghostParticles.z[iP])
#endif
                ) / r;
        }
    }
}

#endif

#endif
#endif
// Calculate acceleration for each particle, w/o Ghosts
// c.f. eq 8 in Monaghan: SPH and its diverse Applications, Annu.Rev. Fluid mechanics, 2012


// void Particles::printDensity(const double &gamma){
//     for (int i = 0; i < N; i += 25){
//         Logger(DEBUG) << "i=" << i << ":rho= " << rho[i] << ", u= " << u[i] << "< p= " << P[i] << "test=" << rho[i]*u[i]*(1 - gamma);
//     }
// }




// Todo: Remove
// Implementing equation F2 in Hopkins` GIZMO Paper
// Storing v_i*dP_i/dt in the variable unimportant, so that only a scalar is needed
// void Particles::calcdP(const Particles &ghostParticles, const double &kernelSize){
//     int iP;
//     double dSqr, r, fij, fji, summand;
//     for (int i = 0; i < N; i++){
//         unimportant[i] = 0;
//         for (int j = 0; j < noi[i]; j++){
//             dSqr = pow((x[i] - x[iP]), 2)
//                         + pow((y[i] - y[iP]), 2);
// #if DIM == 3
//             dSqr += pow(z[i] - z[iP], 2);
// #endif
//             r = sqrt(dSqr);
//
//             fij = 1 - 1 / m[iP]
//                     * kernelSize /omega[i] / DIM * drho[i]
//                     / (1 + kernelSize /omega[i] / DIM * dn[i]);
//
//             summand = m[i]*m[iP] * Kernel::dWdr(r, kernelSize)
//                     * (P[i]/pow(rho[i], 2)*fij
//                         +  P[iP]/pow(rho[iP], 2)*fji);
//
//
//             unimportant[i] -= summand * (vx[i]*(x[i]- x[iP]) + vy[i]*(y[i]-y[iP])
// #if DIM == 3
//                             + vz[i]*(z[i]- z[iP])
// #endif
//                             ) / r;
//         }
//         // Ghosts!
//         for (int k = 0; k < noiGhosts[i]; k++){
//             dSqr = pow((x[i] - x[nnl[k+i*MAX_NUM_INTERACTIONS]]), 2)
//                         + pow((y[i] - y[iP]), 2);
// #if DIM == 3
//             dSqr += pow(z[i] - z[iP], 2);
// #endif
//             r = sqrt(dSqr);
//
//             fij = 1 - 1 / m[ghostParticles.parent[iP]]
//                     * kernelSize /omega[i] / DIM * drho[i]
//                     / (1 + kernelSize /omega[i] / DIM * dn[i]);
//
//             fji = 1 - 1 / m[i]
//                     * kernelSize / omega[iP] / DIM * drho[iP]
//                     / ( 1 + kernelSize / omega[iP] / DIM * dn[iP]);
//
//             summand = m[i]*m[ghostParticles.parent[iP]] * Kernel::dWdr(r, kernelSize)
//                     * (P[i]/pow(rho[i], 2)*fij
//                         +  P[iP]/pow(rho[iP], 2)*fji);
//
//
//             unimportant[i] -= summand * (vx[i]*(x[i]- x[iP]) + vy[i]*(y[i]-y[iP])
// #if DIM == 3
//                             + vz[i]*(z[i]- z[iP])
// #endif
//                             ) / r;
//         }
//     }
// }
//
//


// Implementing equation F2 in Hopkins` GIZMO Paper
// Storing v_i*dP_i/dt in the variable unimportant, so that only a scalar is needed
// void Particles::calcdP(const double &kernelSize){
//     double fij, fji, summand;
//     for (int i = 0; i < N; i++){
//         unimportant[i] = 0;
//         for (int j = 0; j < noi[i]; j++){
//             dSqr = pow((x[i] - x[iP]), 2)
//                         + pow((y[i] - y[iP]), 2);
// #if DIM == 3
//             dSqr += pow(z[i] - z[iP], 2);
// #endif
//             r = sqrt(dSqr);
//
//             fij = 1 - 1 / m[iP]
//                     * kernelSize /omega[i] / DIM * drho[i] / (1 + kernelSize /omega[i] / DIM * dn[i]);
//
//             fji = 1 - 1 / m[i] * kernelSize / omega[iP]/DIM*drho[iP]
//                     / ( 1 + kernelSize / omega[iP] / DIM * dn[iP]);
//
//             summand = m[i]*m[nnList[j+i*MAX_NUM_INTERACTIONS]] * Kernel::dWdr(kernelSize)
//                     * (P[i]/pow(rho[i], 2)*fij
//                         +  P[iP]/pow(rho[iP], 2)*fji);
//
//
//             unimportant[i] += summand * (vx[i]*(x[i]- x[iP]) + vy[i]*(y[i]-y[iP]));
// #if DIM == 3
//             unimportant[i] += summand * vz[i]*(z[i]- z[iP]);
// #endif
//         }
//     }
// }



void Particles::compDensity(const double &kernelSize){
    for(int i=0; i<N; ++i){
//#if PERIODIC_BOUNDARIES
        compOmega(i, kernelSize);
        rho[i] = m[i]*omega[i];
        if(rho[i] <= 0.){
            Logger(WARN) << "Zero or negative density @" << i;
        }
        //if ((i - 0) % 30 == 0){
        //    Logger(DEBUG) << "density from ghosts: i = " << i << ", dnst = " << rho[i];
        //}
    }
}

void Particles::compOmega(int i, const double &kernelSize){
    double omg = 0.;
    int iP;
    for (int j=0; j<noi[i]; ++j){
        iP = nnl[j+i*MAX_NUM_INTERACTIONS];
        double dSqr = pow(x[i] - x[iP], 2)
                    + pow(y[i] - y[iP], 2);
#if DIM == 3
        dSqr += pow(z[i] - z[iP], 2);
#endif
        double r = sqrt(dSqr);
        omg += kernel(r, kernelSize);

        //Logger(DEBUG) << "x[" << i << "] = [" << x[i] << ", " <<  y[i] << "], x["
        //          << iP << "] = [" << x[iP] << ", "
        //          << y[iP] << "]";

    }
    omega[i] = omg + kernel(0., kernelSize); // add self interaction to normalization factor
}

void Particles::compPsijTilde(Helper &helper, const double &kernelSize){
    int iP;
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
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            double dSqr = pow(x[i] - x[iP], 2)
                          + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            double r = sqrt(dSqr);
            double psij_xi = kernel(r, kernelSize)/omega[i];

            xj[0] = x[iP];
            xj[1] = y[iP];
#if DIM==3
            xj[2] = z[iP];
#endif

            for (int alpha=0; alpha<DIM; ++alpha){
                for(int beta=0; beta<DIM; ++beta){
                    B[DIM*alpha+beta] += (xj[alpha] - xi[alpha])*(xj[beta] - xi[beta]) * psij_xi;
                }
            }
        }

        helper.inverseMatrix(B, DIM);

        for (int j=0; j<noi[i]; ++j) {
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            double dSqr = pow(x[i] - x[iP], 2)
                          + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            double r = sqrt(dSqr);
            double psij_xi = kernel(r, kernelSize) / omega[i];

            xj[0] = x[iP];
            xj[1] = y[iP];

#if DIM == 3
            xj[2] = z[iP];
#endif
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
        // if (i % 1 == 0){
        //     Logger(DEBUG) << "Paricle " << i << " has density " << rho[i] << " and u[i] " << u[i];
        // }
        P[i] = (gamma-1.)*rho[i]*u[i];
        // std::cout << P[i] << std::endl;
//#if DIM == 3
//        P[i] = (gamma-1.)*rho[i]*(u[i]+.5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]));
//#else
//        P[i] = (gamma-1.)*rho[i]*(u[i]+.5*(vx[i]*vx[i]+vy[i]*vy[i]));
//#endif
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
            int ij;
            for(ij=0; ij<noi[ji]; ++ij){
                if (nnl[ij+ji*MAX_NUM_INTERACTIONS] == i) break;
            }
            for (int alpha=0; alpha<DIM; ++alpha){
                Aij[i*MAX_NUM_INTERACTIONS+j][alpha] = 1./omega[i]*psijTilde_xi[i*MAX_NUM_INTERACTIONS+j][alpha]
                        - 1./omega[ji]*psijTilde_xi[ij+ji*MAX_NUM_INTERACTIONS][alpha];
            }

            //if(i < 10){
            //    Logger(DEBUG) << "A[" << i << " -> " << ji << "] = [" << Aij[i*MAX_NUM_INTERACTIONS+j][0] << ", "
            //                  << Aij[i*MAX_NUM_INTERACTIONS+j][1] << ", " << Aij[i*MAX_NUM_INTERACTIONS+j][2] << "], xi = ["
            //                  << x[i] << ", " << y[i] << ", " << z[i] << "], xj = " << x[ji] << ", " << y[ji] << ", " << z[ji] << "]";
            //}
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
    //Logger(DEBUG) << "#################### NEXT GRADIENT ####################";

    for (int i=0; i<N; ++i){
        double psiMaxNgb { std::numeric_limits<double>::min() };
        double psiMinNgb { std::numeric_limits<double>::max() };
        double psiMaxMid { std::numeric_limits<double>::min() };
        double psiMinMid { std::numeric_limits<double>::max() };

        for (int jn=0; jn<noi[i]; ++jn) {
            int j = nnl[i * MAX_NUM_INTERACTIONS + jn];

            //xij[0] = x[i] + kernelSize / 2. * (x[j] - x[i]);
            //xij[1] = y[i] + kernelSize / 2. * (y[j] - y[i]);
#if FIRST_ORDER_QUAD_POINT
            xij[0] = (x[i] + x[j])/2.;
            xij[1] = (y[i] + y[j])/2.;
#else
            xij[0] = x[i] + kernelSize / 4. * (x[j] - x[i]);
            xij[1] = y[i] + kernelSize / 4. * (y[j] - y[i]);
#endif

            xijxi[0] = xij[0] - x[i];
            xijxi[1] = xij[1] - y[i];
#if DIM == 3
#if FIRST_ORDER_QUAD_POINT
            xij[2] = (z[i] + z[j])/2.;
#else
            //xij[2] = z[i] + kernelSize/2. * (z[j] - z[i]);
            xij[2] = z[i] + kernelSize/4. * (z[j] - z[i]);
#endif
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

            //xij[0] = x[i] + kernelSize / 2. * (ghostParticles->x[j] - x[i]);
            //xij[1] = y[i] + kernelSize / 2. * (ghostParticles->y[j] - y[i]);

#if FIRST_ORDER_QUAD_POINT
            xij[0] = (x[i] + ghostParticles->x[j])/2.;
            xij[1] = (y[i] + ghostParticles->y[j])/2.;
#else
            xij[0] = x[i] + kernelSize / 4. * (ghostParticles->x[j] - x[i]);
            xij[1] = y[i] + kernelSize / 4. * (ghostParticles->y[j] - y[i]);
#endif

            xijxi[0] = xij[0] - x[i];
            xijxi[1] = xij[1] - y[i];
#if DIM == 3
#if FIRST_ORDER_QUAD_POINT
            xij[2] = (z[i] + ghostParticles->z[j])/2.;
#else
            //xij[2] = z[i] + kernelSize/2. * (ghostParticles->z[j] - z[i]);
            xij[2] = z[i] + kernelSize/4. * (ghostParticles->z[j] - z[i]);
#endif

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

        //if(i==6){
        //    Logger(DEBUG) << "alphaMin = " << alphaMin << ", alphaMax = " << alphaMax
        //                << ", psiMinNgb = " << psiMinNgb << ", psiMaxNgb = " << psiMaxNgb
        //                << ", psiMinNgb = " << psiMinMid << ", psiMaxNgb = " << psiMaxMid
        //                << ", f[i] = " << f[i];
        //}

        if (alphaMin <= alphaMax && BETA*alphaMin < 1.){
            //Logger(DEBUG) << "        > Limiting gradient with alphaMin@" << i << ", alpha = "
            //              << alphaMin << ", grad[0] = " << grad[i][0] << ", grad[1] = " << grad[i][1];
            grad[i][0] *= alphaMin;
            grad[i][1] *= alphaMin;
#if DIM==3
            grad[i][2] *= alphaMin;
#endif
        } else if (alphaMax <= alphaMin && BETA*alphaMax < 1.){
            //Logger(DEBUG) << "        > Limiting gradient with alphaMax@" << i << ", alpha = "
            //            << alphaMax << ", grad[0] = " << grad[i][0] << ", grad[1] = " << grad[i][1];
            grad[i][0] *= alphaMax;
            grad[i][1] *= alphaMax;
#if DIM==3
            grad[i][2] *= alphaMax;
#endif
        }
    }
}

double Particles::compGlobalTimestep(const double &gamma, const double &kernelSize){
    double dt_ = std::numeric_limits<double>::max();
    for (int i=0; i<N; ++i){

        double vSig = std::numeric_limits<double>::min();
        double ci = sqrt(gamma*P[i]/rho[i]); // soundspeed @i

        // searching for maximum signal speed
        for (int jn=0; jn<noi[i]; ++jn){
            int j = nnl[i*MAX_NUM_INTERACTIONS+jn];

            double cj = sqrt(gamma*P[j]/rho[j]); // soundspeed @j

            double xij[DIM], vij[DIM];

            xij[0] = x[i] - x[j];
            xij[1] = y[i] - y[j];

            vij[0] = vx[i] - vx[j];
            vij[1] = vy[i] - vy[j];

#if DIM == 3
            xij[2] = z[i] - z[j];
            vij[2] = vz[i] - vz[j];
#endif
            double vijxij = Helper::dotProduct(vij, xij)/sqrt(Helper::dotProduct(xij, xij));
            vijxij = vijxij < 0. ? vijxij : 0.;

            double vSig_i = ci+cj-vijxij;
            vSig = vSig_i > vSig ? vSig_i : vSig;

        }

        // TODO: Note: Kernel size is double the actual kernel size
        double dt = CFL*kernelSize/vSig;
        dt_ = dt < dt_ ? dt : dt_;
    }

    return dt_;
}


void Particles::compRiemannStatesLR(const double &dt, const double &kernelSize, const double &gamma){
    for (int i=0; i<N; ++i){
        double xij[DIM];
        //double vFrame[DIM];
        // helper vectors
        double xijxi[DIM], xjxi[DIM], xijxj[DIM];

        for (int jn=0; jn<noi[i]; ++jn){
            int j = nnl[i*MAX_NUM_INTERACTIONS+jn];

            xjxi[0] = x[j] - x[i];
            xjxi[1] = y[j] - y[i];

            //xij[0] = x[i] + kernelSize/2. * xjxi[0];
            //xij[1] = y[i] + kernelSize/2. * xjxi[1];

#if FIRST_ORDER_QUAD_POINT
            xij[0] = (x[i] + x[j])/2.;
            xij[1] = (y[i] + y[j])/2.;

            xijxj[0] = .5*(x[i] - x[j]);
            xijxj[1] = .5*(y[i] - y[j]);

            xijxi[0] = .5*(x[j] - x[i]);
            xijxi[1] = .5*(y[j] - y[i]);

#else
            xij[0] = x[i] + kernelSize/4. * xjxi[0];
            xij[1] = y[i] + kernelSize/4. * xjxi[1];

            xijxi[0] = xij[0] - x[i];
            xijxi[1] = xij[1] - y[i];

            xijxj[0] = xij[0] - x[j];
            xijxj[1] = xij[1] - y[j];
#endif

#if DIM==3
#if FIRST_ORDER_QUAD_POINT
            xij[2] = (z[i] + z[j])/2.;
            xijxj[2] = .5*(z[i] - z[j]);
            xijxi[2] = .5*(z[j] - z[i]);
#else
            xjxi[2] = z[j] - z[i];
            //xij[2] = z[i] + kernelSize/2. * xjxi[2];

            xij[2] = z[i] + kernelSize/4. * xjxi[2];

            xijxi[2] = xij[2] - z[i];
            xijxj[2] = xij[2] - z[j];
#endif
#endif

            int iW = i*MAX_NUM_INTERACTIONS+jn;
#if !MOVE_PARTICLES
            vFrame[iW][0] = 0.;
            vFrame[iW][1] = 0.;
#if DIM==3
            vFrame[iW][2] = 0.;
#endif
#else // MOVE_PARTICLES
#if FIRST_ORDER_QUAD_POINT
            vFrame[iW][0] = (vx[i] + vx[j])/2.;
            vFrame[iW][1] = (vy[i] + vy[j])/2.;
#if DIM==3
            vFrame[iW][2] = (vz[i] + vz[j])/2.;
#endif
#else
            double dotProd = Helper::dotProduct(xijxi, xjxi);
            double dSqr = Helper::dotProduct(xjxi, xjxi);

            vFrame[iW][0] = vx[i] + (vx[j]-vx[i]) * dotProd/dSqr;
            vFrame[iW][1] = vy[i] + (vy[j]-vy[i]) * dotProd/dSqr;
#if DIM==3
            vFrame[iW][2] = vz[i] + (vz[j]-vz[i]) * dotProd/dSqr;
#endif
#endif // FIRST_ORDER_QUAD_POINT
#endif // MOVE_PARTICLES
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

#if PAIRWISE_LIMITER
            double WijR_buf[DIM+2], WijL_buf[DIM+2];
            for(int nu=0; nu<DIM+2; ++nu){
                WijR_buf[nu] = WijR[iW][nu];
                WijL_buf[nu] = WijL[iW][nu];
            }
#endif

            //if (i == 46){// && jn == 28){
            //    Logger(DEBUG) << "        j = " << j
            //                  << ", rhoL = " << WijL[iW][0] << ", rhoR = " << WijR[iW][0]
            //                  << ", uL = " << WijL[iW][2] << ", uR = " << WijR[iW][2]
            //                  << ", PL = " << WijL[iW][1] << ", PR = " << WijR[iW][1];
            //}

            //if(i == 46){// && jn == 28){
            //    Logger(DEBUG) << "vFrame[iW] = [" << vFrame[iW][0]
            //                << ", " << vFrame[iW][1] << "]"
            //                << ", rhoGrad[i] = [" << rhoGrad[i][0]
            //                << ", " << rhoGrad[i][1] << "]";
            //    Logger(DEBUG) << "PGrad[i] = [" << PGrad[i][0] << ", " << PGrad[i][1]
            //                << "], PGrad[j] = [" << PGrad[j][0] << ", " << PGrad[j][1]
            //                << "], rho[i] = " << rho[i]
                //            << "], xijxi = [" << xijxi[0] << ", " << xijxi[1]
                //            << "], xijxj = [" << xijxj[0] << ", " << xijxj[1] << "]";
            //                << ", xj = [" << x[j] << ", " << y[j] << "] @" << j;
                //exit(5);
            //}
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

#if PAIRWISE_LIMITER
            double xijxi_abs = 0., xijxj_abs = 0., xjxi_abs = 0.;
            for(int alpha=0; alpha<DIM; ++alpha){
                xijxi_abs += xijxi[alpha]*xijxi[alpha];
                xijxj_abs += xijxj[alpha]*xijxj[alpha];
                xjxi_abs += xjxi[alpha]*xjxi[alpha];
            }
            xijxi_abs = sqrt(xijxi_abs);
            xijxj_abs = sqrt(xijxj_abs);
            xjxi_abs = sqrt(xjxi_abs);

            for(int nu=0; nu<DIM+2; ++nu){
                WijR[iW][nu] = pairwiseLimiter(WijR[iW][nu], WijR_buf[nu], WijL_buf[nu], xijxi_abs, xjxi_abs);
                WijL[iW][nu] = pairwiseLimiter(WijL[iW][nu], WijL_buf[nu], WijR_buf[nu], xijxj_abs, xjxi_abs);
            }
#endif

#if DEBUG_LVL
            if(WijR[iW][1] < 0. || WijL[iW][1] < 0.) {
                Logger(WARN) << "FACE RECONSTRUCTION > Negative pressure encountered@(i = " << i << ", j = " << j << ")";
                Logger(DEBUG) << "rhoL = " << WijL[iW][0] << ", rhoR = " << WijR[iW][0]
                              << ", uL = " << WijL[iW][2] << ", uR = " << WijR[iW][2]
                              << ", PL = " << WijL[iW][1] << ", PR = " << WijR[iW][1]
                              << ", PGrad[i] = [" << PGrad[i][0] << ", " << PGrad[i][1]
#if DIM == 3
                              << ", " << PGrad[i][2]
#endif
                              << "] , PGrad[j] = [" << PGrad[j][0] << ", " << PGrad[j][1]
#if DIM == 3
                              << ", " << PGrad[j][2]
#endif
                              << "]";
            }
#endif
            //if (i == 46){
            //    Logger(DEBUG) << "        j = " << j
            //                  << ", rhoL = " << WijL[iW][0] << ", rhoR = " << WijR[iW][0]
            //                  << ", uL = " << WijL[iW][2] << ", uR = " << WijR[iW][2]
            //                  << ", PL = " << WijL[iW][1] << ", PR = " << WijR[iW][1];
            //}

            // predict half a timestep
            double viDiv = vxGrad[i][0] + vyGrad[i][1];
            double vjDiv = vxGrad[j][0] + vyGrad[j][1];
#if DIM==3
            viDiv += vzGrad[i][2];
            vjDiv += vzGrad[j][2];
#endif
            // density
            WijR[iW][0] -= dt/2. * (rho[i] * viDiv + (vx[i]-vFrame[iW][0])*rhoGrad[i][0] + (vy[i]-vFrame[iW][1])*rhoGrad[i][1]);
            WijL[iW][0] -= dt/2. * (rho[j] * vjDiv + (vx[j]-vFrame[iW][0])*rhoGrad[j][0] + (vy[j]-vFrame[iW][1])*rhoGrad[j][1]);

            // energy
            WijR[iW][1] -= dt/2. * (gamma*P[i] * viDiv + (vx[i]-vFrame[iW][0])*PGrad[i][0] + (vy[i]-vFrame[iW][1])*PGrad[i][1]);
            WijL[iW][1] -= dt/2. * (gamma*P[j] * vjDiv + (vx[j]-vFrame[iW][0])*PGrad[j][0] + (vy[j]-vFrame[iW][1])*PGrad[j][1]);

            // velocities
            // TODO: center vL and vR and update vFrame (compare to GIZMO code hydro_core_meshless.h:178ff) ??
            WijR[iW][2] -= dt/2. * (PGrad[i][0]/rho[i] + (vx[i]-vFrame[iW][0])*vxGrad[i][0] + (vy[i]-vFrame[iW][1])*vxGrad[i][1]);
            WijL[iW][2] -= dt/2. * (PGrad[j][0]/rho[j] + (vx[j]-vFrame[iW][0])*vxGrad[j][0] + (vy[j]-vFrame[iW][1])*vxGrad[j][1]);
            WijR[iW][3] -= dt/2. * (PGrad[i][1]/rho[i] + (vx[i]-vFrame[iW][0])*vyGrad[i][0] + (vy[i]-vFrame[iW][1])*vyGrad[i][1]);
            WijL[iW][3] -= dt/2. * (PGrad[j][1]/rho[j] + (vx[j]-vFrame[iW][0])*vyGrad[j][0] + (vy[j]-vFrame[iW][1])*vyGrad[j][1]);
#if DIM==3
            // density
            WijR[iW][0] -= dt/2. * (vz[i]-vFrame[iW][2])*rhoGrad[i][2];
            WijL[iW][0] -= dt/2. * (vz[j]-vFrame[iW][2])*rhoGrad[j][2];

            // energy
            WijR[iW][1] -= dt/2. * (vz[i]-vFrame[iW][2])*PGrad[i][2];
            WijL[iW][1] -= dt/2. * (vz[j]-vFrame[iW][2])*PGrad[j][2];

#if DEBUG_LVL
            if(WijR[iW][1] < 0. || WijL[iW][1] < 0.) {
                Logger(WARN) << "TIME PREDICTION > Negative pressure encountered@(i = " << i << ", j = " << j << ")";

                double timePredP_i = dt / 2. * (gamma * P[i] * viDiv + (vx[i] - vFrame[iW][0]) * PGrad[i][0] +
                                                (vy[i] - vFrame[iW][1]) * PGrad[i][1]
#if DIM == 3
                                                + (vz[i] - vFrame[iW][2]) * PGrad[i][2]
#endif
                );
                Logger(DEBUG) << "Pressure timestep prediction term @i: " << timePredP_i;

                double timePredP_j = dt / 2. * (gamma * P[j] * viDiv + (vx[j] - vFrame[iW][0]) * PGrad[j][0] +
                                                (vy[j] - vFrame[iW][1]) * PGrad[j][1]
#if DIM == 3
                                                + (vz[j] - vFrame[iW][2]) * PGrad[j][2]
#endif
                );
                Logger(DEBUG) << "Pressure timestep prediction term @j: " << timePredP_j;
            }
#endif

            // velocities
            WijR[iW][2] -= dt/2. * (vz[i]-vFrame[iW][2])*vxGrad[i][2];
            WijL[iW][2] -= dt/2. * (vz[i]-vFrame[iW][2])*vxGrad[j][2];
            WijR[iW][3] -= dt/2. * (vz[i]-vFrame[iW][2])*vyGrad[i][2];
            WijL[iW][3] -= dt/2. * (vz[i]-vFrame[iW][2])*vyGrad[j][2];
            WijR[iW][4] -= dt/2. * (PGrad[i][2]/rho[i] + (vx[i]-vFrame[iW][0])*vzGrad[i][0] + (vy[i]-vFrame[iW][1])*vzGrad[i][1] + (vz[i]-vFrame[iW][2])*vzGrad[i][2]);
            WijL[iW][4] -= dt/2. * (PGrad[j][2]/rho[j] + (vx[j]-vFrame[iW][0])*vzGrad[j][0] + (vy[j]-vFrame[iW][1])*vzGrad[j][1] + (vz[j]-vFrame[iW][2])*vzGrad[j][2]);
#endif

            //if (i == 46){// && jn == 28){
            //    Logger(DEBUG) << "        j = " << j
            //                  << ", rhoL = " << WijL[iW][0] << ", rhoR = " << WijR[iW][0]
            //                  << ", uL = " << WijL[iW][2] << ", uR = " << WijR[iW][2]
            //                  << ", PL = " << WijL[iW][1] << ", PR = " << WijR[iW][1];
            //}

        }
    }
}

double Particles::pairwiseLimiter(double phi0, double phi_i, double phi_j, double xijxi_abs, double xjxi_abs) {
    double phi_ = phi_i;

    /// calculate helper values
    double phi_ij = phi_i + xijxi_abs / xjxi_abs * (phi_j - phi_i);
    double phiMin, phiMax;
    if (phi_i < phi_j) {
        phiMin = phi_i;
        phiMax = phi_j;
    } else {
        phiMin = phi_j;
        phiMax = phi_i;
    }
    double delta1 = PSI_1 * abs(phi_i - phi_j);
    double delta2 = PSI_2 * abs(phi_i - phi_j);
    double phiMinus, phiPlus;
    if ((phiMax + delta1 >= 0. && phiMax >= 0.) || (phiMax + delta1 < 0. && phiMax < 0.)) {
        phiPlus = phiMax + delta1;
    } else {
        phiPlus = phiMax / (1. + delta1 / abs(phiMax));
    }
    if ((phiMin - delta1 >= 0. && phiMin >= 0.) || (phiMin - delta1 < 0. && phiMin < 0.)) {
        phiMinus = phiMin - delta1;
    } else {
        phiMinus = phiMin / (1. + delta1 / abs(phiMin));
    }

    /// actually compute the effective face limited value
    if (phi_i < phi_j) {
        double minPhiD2;
        if (phi_ij + delta2 < phi0){
            minPhiD2 = phi_ij + delta2;
        } else {
            minPhiD2 = phi0;
        }

        phi_ = phiMinus > minPhiD2 ? phiMinus : minPhiD2;

    } else if (phi_i > phi_j){
        double maxPhiD2;
        if (phi_ij - delta2 > phi0){
            maxPhiD2 = phi_ij - delta2;
        } else {
            maxPhiD2 = phi0;
        }

        phi_ = phiPlus < maxPhiD2 ? phiPlus : maxPhiD2;

    }
    return phi_;
}

void Particles::solveRiemannProblems(const double &gamma, const Particles &ghostParticles){
#if USE_HLLC
    double n_unit[DIM];
#endif
    for (int i=0; i<N; ++i){

        //if (!(i % (N/VERBOSITY_PARTICLES))){
        //    Logger(DEBUG) << "        > i = " << i;
        //}
        //Logger(DEBUG) << "        > i = " << i << ", V = " << 1./omega[i];

        for (int j=0; j<noi[i]; ++j){
            int ii = i*MAX_NUM_INTERACTIONS+j; // interaction index
            //if (i == 9 && j == 11){
                //Logger(DEBUG) << "i = " << i << ", ii = " << ii << ", j = " << nnl[ii];
                //Logger(DEBUG) << "xi = [" << x[i] << ", " << y[i] << "], xj = ["
                //              << x[nnl[ii]] << ", " << y[nnl[ii]] << "] , vi = ["
                //              << vx[i] << ", " << vy[i] << "], vj = ["
                //              << vx[nnl[ii]] << ", " << vy[nnl[ii]] << "]";
                //Logger(DEBUG) << "vFrame = [" << vFrame[ii][0] << ", " << vFrame[ii][1] << "]";
                //Logger(DEBUG) << "rhoL = " << WijL[ii][0] << ", rhoR = " << WijR[ii][0]
                //              << ", uL = " << WijL[ii][2] << ", uR = " << WijR[ii][2]
                //              << ", PL = " << WijL[ii][1] << ", PR = " << WijR[ii][1]
                //              << ", Aij = [" << Aij[ii][0] << ", " << Aij[ii][1] << "]";
            //}

            //Logger(DEBUG) << "*WR = " << WijR[ii] << ", *WL = " << WijL[ii];
            //Logger(DEBUG) << "i = " << i << ", j = " << nnl[ii] << ", PR = " << WijR[ii][1] << ", PL = " << WijL[ii][1];

            if(WijR[ii][1] < 0. || WijL[ii][1] < 0.){
                Logger(WARN) << "Negative pressure encountered@(i = " << i << ", j = " << j << ") Very bad :( !!";
                Logger(DEBUG) << "    > rhoL = " << WijL[ii][0] << ", rhoR = " << WijR[ii][0]
                              << ", uL = " << WijL[ii][2] << ", uR = " << WijR[ii][2]
                              << ", PL = " << WijL[ii][1] << ", PR = " << WijR[ii][1]
                              << ", Aij = [" << Aij[ii][0] << ", " << Aij[ii][1]
#if DIM ==3
                              << ", " << Aij[ii][2]
#endif
                              << "]";

#if DEBUG_LVL > 1
                Logger(DEBUG) << "Aborting for debugging.";
                exit(6);
#endif
                //if(WijR[ii][1] < 0.){
                //    WijR[ii][1] = PRESSURE_FLOOR;
                //}
                //if(WijL[ii][1] < 0.){
                //    WijL[ii][1] = PRESSURE_FLOOR;
                //}
            }

            bool compute = true;
            int iij;
#if ENFORCE_FLUX_SYM
            //int ii = i*MAX_NUM_INTERACTIONS+j; // interaction index i->j
            int ji = nnl[ii]; // index i of particle j
            if(ji<i) {
                compute = false;
                // search neighbor i in nnl[] of j
                int ij;
                for (ij = 0; ij < noi[ji]; ++ij) {
                    if (nnl[ij + ji * MAX_NUM_INTERACTIONS] == i) break;
                }
                iij = ji * MAX_NUM_INTERACTIONS + ij; // interaction index j->i
            }
#endif
            if (compute){
#if USE_HLLC
                calcNunit(i, ii, n_unit);
                Riemann solver { WijL[ii], WijR[ii], vFrame[ii], Aij[ii], n_unit, i};
                solver.HLLC(Fij[ii], gamma);
#else
                Riemann solver { WijL[ii], WijR[ii], vFrame[ii], Aij[ii] , i };
                solver.exact(Fij[ii], gamma);
#endif
            } else {
                for(int d=0; d<DIM+2; ++d){
                    Fij[ii][d] = -Fij[iij][d];
                }
            }

            //if(i == 6){//&& j==11){
            //    Logger(DEBUG) << "Fluxes = [" << Fij[ii][0] << ", " << Fij[ii][1] << ", " << Fij[ii][2] << ", " << Fij[ii][3] << "]";
            //}
        }

#if PERIODIC_BOUNDARIES
        for (int j=0; j<noiGhosts[i]; ++j){
            int ii = i*MAX_NUM_GHOST_INTERACTIONS+j; // interaction index
            //Logger(DEBUG) << "i = " << i << ", ii = " << ii << ", jGhost = " << nnlGhosts[ii];
            //Logger(DEBUG) << "xi = [" << x[i] << ", " << y[i] << "], xj = ["
            //              << ghostParticles.x[nnlGhosts[ii]] << ", " << ghostParticles.y[nnlGhosts[ii]] << "]";
            //Logger(DEBUG) << "vFrameGhosts = [" << vFrameGhosts[ii][0] << ", " << vFrameGhosts[ii][1] << "]";

            if(WijRGhosts[ii][1] < 0. || WijLGhosts[ii][1] < 0.){
                Logger(WARN) << "Negative pressure encountered@(i = " << i << ", jGhost = " << j << ") Very bad :( !!";
                Logger(DEBUG) << "rhoL = " << WijLGhosts[ii][0] << ", rhoR = " << WijRGhosts[ii][0]
                              << ", uL = " << WijLGhosts[ii][2] << ", uR = " << WijRGhosts[ii][2]
                              << ", PL = " << WijLGhosts[ii][1] << ", PR = " << WijRGhosts[ii][1];
#if DEBUG_LVL
                Logger(DEBUG) << "Aborting for debugging.";
                exit(6);
#endif
            }

            bool compute = true;
            int iij;
#if ENFORCE_FLUX_SYM
            //int ii = i*MAX_NUM_GHOST_INTERACTIONS+j; // interaction index i->j
            int ji = nnlGhosts[ii]; // index i of particle j
            if (ghostParticles.parent[ji]<i){
                compute = false;
                int ij;
                // search neighbor i in nnlGhosts[] of j
                for (ij=0; ij<noiGhosts[ghostParticles.parent[ji]]; ++ij){
                    if(ghostParticles.parent[nnlGhosts[ij+ghostParticles.parent[ji]*MAX_NUM_GHOST_INTERACTIONS]] == i) break;
                }
                iij = ij+ghostParticles.parent[ji]*MAX_NUM_GHOST_INTERACTIONS;
            }
#endif

            if (compute) {
#if USE_HLLC
                calcNunit(i, ii, n_unit);
                Riemann solver{WijLGhosts[ii], WijRGhosts[ii], vFrameGhosts[ii], AijGhosts[ii], n_unit, i};
                solver.HLLC(FijGhosts[ii], gamma);
#else
                Riemann solver{WijLGhosts[ii], WijRGhosts[ii], vFrameGhosts[ii], AijGhosts[ii], i};
                solver.exact(FijGhosts[ii], gamma);
#endif
            } else {
                for(int d=0; d<DIM+2; ++d){
                    FijGhosts[ii][d] = -FijGhosts[iij][d];
                }
            }
        }
#endif
    }
}

void Particles::collectFluxes(Helper &helper, const Particles &ghostParticles){
    for (int i=0; i<N; ++i){
        mF[i] = 0.;

        vF[i][0] = 0.;
        vF[i][1] = 0.;
#if DIM == 3
        vF[i][2] = 0.;
#endif
        eF[i] = 0.;

        //Logger(DEBUG) << "      > i = " << i;

        for(int j=0; j<noi[i]; ++j){
            int ii = j+i*MAX_NUM_INTERACTIONS;
            //double AijNorm = sqrt(Helper::dotProduct(Aij[ii], Aij[ii]));

            /// MASS FLUXES
            //mF[i] += AijNorm*Fij[ii][0];
            mF[i] += Fij[ii][0];

            //Logger(DEBUG) << "xi = [" << x[i] << ", " << y[i] << "]"
            //          << ", xj = [" << x[nnl[ii]] << ", " << y[nnl[ii]] << "]"
            //          << ", mF[ii] = " << Fij[ii][0] << ", AijNorm = " << AijNorm;

            /// VELOCITY FLUXES
            // add de-boosted velocities
            //vF[i][0] += AijNorm * (Fij[ii][2] + Fij[ii][0]*vFrame[ii][0]);
            //vF[i][1] += AijNorm * (Fij[ii][3] + Fij[ii][0]*vFrame[ii][1]);

            vF[i][0] += Fij[ii][2];
            vF[i][1] += Fij[ii][3];
#if DIM==3
            vF[i][2] += Fij[ii][4];
#endif

            /// ENERGY FLUXES
            // allocate buffer for energy update
            //double Fv[DIM];
            //Fv[0] = Fij[ii][2];
            //Fv[1] = Fij[ii][3];
            //eF[i] += AijNorm*(Fij[ii][1] + .5*Helper::dotProduct(vFrame[ii], vFrame[ii])*Fij[ii][0]
            //                  + Helper::dotProduct(vFrame[ii], Fv));

            eF[i] += Fij[ii][1];

            //if (i == 46){
            //    Logger(DEBUG) << "  > j = " << nnl[ii] << ", AijNorm = " << AijNorm
            //              << ", vFrame = [" << vFrame[ii][0] << ", " << vFrame[ii][1] << "]";
            //    Logger(DEBUG) << "  > Fm = " << mF[i] << ", Fv = [" << vF[i][0] << ", " << vF[i][1]
            //                  << "], Fe = " << eF[i];
            //}

        }

#if PERIODIC_BOUNDARIES
        for(int j=0; j<noiGhosts[i]; ++j){
            int ii = j+i*MAX_NUM_GHOST_INTERACTIONS;
            //double AijNorm = sqrt(Helper::dotProduct(AijGhosts[ii], AijGhosts[ii]));

            /// MASS FLUXES
            //mF[i] += AijNorm*FijGhosts[ii][0];
            mF[i] += FijGhosts[ii][0];

            //Logger(DEBUG) << "xi = [" << x[i] << ", " << y[i] << "]"
            //              << ", xjGhost = [" << ghostParticles.x[nnlGhosts[ii]] << ", " << ghostParticles.y[nnlGhosts[ii]] << "]"
            //              << ", mF[ii] = " << FijGhosts[ii][0] << ", AijNorm = " << AijNorm;

            /// VELOCITY FLUXES
            // add de-boosted velocities
            //vF[i][0] += AijNorm * (FijGhosts[ii][2] + FijGhosts[ii][0]*vFrameGhosts[ii][0]);
            //vF[i][1] += AijNorm * (FijGhosts[ii][3] + FijGhosts[ii][0]*vFrameGhosts[ii][1]);

            vF[i][0] += FijGhosts[ii][2];
            vF[i][1] += FijGhosts[ii][3];


            // TODO: implement z-component to work for 3D

            /// ENERGY FLUXES
            // allocate buffer for energy update
            //double Fv[DIM];
            //Fv[0] = FijGhosts[ii][2];
            //Fv[1] = FijGhosts[ii][3];
            //eF[i] += AijNorm*(FijGhosts[ii][1] + .5*Helper::dotProduct(vFrameGhosts[ii], vFrameGhosts[ii])*FijGhosts[ii][0]
            //                  + Helper::dotProduct(vFrameGhosts[ii], Fv));

            eF[i] += FijGhosts[ii][1];

            //if (i == 46){
            //    Logger(DEBUG) << "  > jGhost = " << nnlGhosts[ii] << ", AijNorm = " << AijNorm
            //                  << ", vFrame = [" << vFrameGhosts[ii][0] << ", " << vFrameGhosts[ii][1] << "]";
            //    Logger(DEBUG) << "  > Fm = " << mF[i] << ", Fv = [" << vF[i][0] << ", " << vF[i][1]
            //                  << "], Fe = " << eF[i];
            //}
        }
#endif
    }
}

void Particles::updateStateAndPosition(const double &dt, const Domain &domain){

    for(int i=0; i<N; ++i){

        // store velocity for position update
        double vxi = vx[i], vyi = vy[i];
#if DIM==3
        double vzi = vz[i];
#endif

        // create state vector without mass
        double Q[DIM+1];
#if DIM==3
        Q[0] = m[i]*(u[i] + .5*(vxi*vxi+vyi*vyi+vzi*vzi));
#else
        Q[0] = m[i]*(u[i] + .5*(vxi*vxi+vyi*vyi));
#endif
        Q[1] = m[i]*vxi;
        Q[2] = m[i]*vyi;
#if DIM==3
        Q[3] = m[i]*vzi;
#endif
        //if(i == 46){
        //    Logger(DEBUG) << "i = " << i << ", mass flux = " << mF[i]
        //                  << ", momentum flux = [" << vF[i][0] << ", " << vF[i][1] << "], energy flux = " << eF[i];
        //}

        // UPDATE MASS
#if !MESHLESS_FINITE_MASS
        m[i] -= dt*mF[i];
#endif
        if (m[i] <= 0.){
            Logger(ERROR) << "Negative mass. m[" << i << "] =" << m[i] << ", mF = " << mF[i];
            //m[i] = MASS_FLOOR;
        }

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
        //Logger(DEBUG) << "New velocity v = [" << vx[i] << ", " << vy[i] << ", " << vz[i] << "]";
#endif
        // UPDATE INTERNAL ENERGY
        Q[0] -= dt*eF[i];
#if DIM==3
        u[i] = Q[0]/m[i]-.5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);

        //Logger(DEBUG) << "Total energy Q_E = " << Q[0] << ", kinetic energy E_kin = "
        //          << .5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);

#else
        u[i] = Q[0]/m[i]-.5*(vx[i]*vx[i]+vy[i]*vy[i]);
#endif
        if (u[i] <= 0.){
            Logger(ERROR) << "Negative internal energy. u[" << i << "] =" << u[i] << ", eF = " << eF[i];
            //u[i] = ENERGY_FLOOR;
        }

#if MOVE_PARTICLES
        // MOVE PARTICLES
        //x[i] += .5*(vx[i]+vxi)*dt;
        //y[i] += .5*(vy[i]+vyi)*dt;
        // stay consistent with choice of effective frame
        x[i] += vxi*dt;
        y[i] += vyi*dt;

#if DIM==3
        //z[i] += .5*(vz[i]+vzi)*dt;
        z[i] += vzi*dt;
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
#endif // PERIODIC_BOUNDARIES
#endif // MOVE_PARTICLES
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
            ghostParticles.parent[iGhost] = i;
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
            ghostParticles.parent[iGhost] = i;
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
            ghostParticles.parent[iGhost] = i;
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
            ghostParticles.vy[ghostMap[i]] = vy[i/(DIM+1)];
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

/*void Particles::updateGhostPsijTilde(Particles &ghostParticles){
    for (int i=0; i<N*(DIM+1); ++i){
        if (ghostMap[i] >= 0){
            for (int j=0; j<noiGhosts[i/(DIM+1)]; ++j){
                for (int alpha=0; alpha<DIM; ++alpha){
                    ghostParticles.psijTilde_xiGhosts[j+ghostMap[i]*MAX_NUM_GHOST_INTERACTIONS][alpha]
                                    = psijTilde_xiGhosts[j+i/(DIM+1)*MAX_NUM_GHOST_INTERACTIONS][alpha];
                }
            }
        }
    }
}*/

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
    int iP;
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
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            double dSqr = pow(x[i] - x[iP], 2)
                          + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            double r = sqrt(dSqr);
            //double psij_xi = kernel(r, kernelSize)/omega[nnl[j + i * MAX_NUM_INTERACTIONS]];
            double psij_xi = kernel(r, kernelSize)/omega[i];

            xj[0] = x[iP];
            xj[1] = y[iP];
#if DIM==3
            xj[2] = z[iP];
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
            iP = nnl[j+i*MAX_NUM_INTERACTIONS];
            double dSqr = pow(x[i] - x[iP], 2)
                          + pow(y[i] - y[iP], 2);
#if DIM == 3
            dSqr += pow(z[i] - z[iP], 2);
#endif
            double r = sqrt(dSqr);
            //double psij_xi = kernel(r, kernelSize) / omega[nnl[j + i * MAX_NUM_INTERACTIONS]];
            double psij_xi = kernel(r, kernelSize) / omega[i];

            xj[0] = x[iP];
            xj[1] = y[iP];

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

                //if(i == 6){
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

        //Logger(DEBUG) << "        > grad[" << i << "] = [" << grad[i][0] << ", "
        //          << grad[i][1] << "]";

        //for (int alpha=0; alpha < DIM; ++alpha ){
        //    if(i == 86) Logger(DEBUG) << "grad[" << i << "][" << alpha << "] = " << grad[i][alpha];
        //}
    }
}

void Particles::compEffectiveFace(const Particles &ghostParticles){
    for (int i=0; i<N; ++i){

        //Logger(DEBUG) << "      > i = " << i << ", x = [" << x[i] << ", " << y[i] <<
        //          "], omega = " << omega[i];

        for (int j=0; j<noiGhosts[i]; ++j){
            int ii = i*MAX_NUM_GHOST_INTERACTIONS+j;
            int ji = nnlGhosts[ii]; // index i of particle j
            int ij;
            for (ij=0; ij<noiGhosts[ghostParticles.parent[ji]]; ++ij){
                if(ghostParticles.parent[nnlGhosts[ij+ghostParticles.parent[ji]*MAX_NUM_GHOST_INTERACTIONS]] == i) break;
            }

            //Logger(DEBUG) << "j = " << j << ", ii = " << ii << ", ji = " << ji
            //          << ", ij = " << ij;
            //Logger(DEBUG) << "omegaGhost[ji] = " << ghostParticles.omega[ji]
            //          << ", psijTilde_xiGhost[ii] = [" << psijTilde_xiGhosts[ii][0] << ", " << psijTilde_xiGhosts[ii][1]
            //          << "], psijTilde_xiGhost[ji] = [" << psijTilde_xiGhosts[ij+ghostParticles.parent[ji]*MAX_NUM_GHOST_INTERACTIONS][0]
            //          << ", " << psijTilde_xiGhosts[ij+ghostParticles.parent[ji]*MAX_NUM_GHOST_INTERACTIONS][1] << "]";

            for (int alpha=0; alpha<DIM; ++alpha){
                AijGhosts[ii][alpha] = 1./omega[i]*psijTilde_xiGhosts[ii][alpha] - 1./ghostParticles.omega[ji]
                                       * psijTilde_xiGhosts[ij+ghostParticles.parent[ji]*MAX_NUM_GHOST_INTERACTIONS][alpha];
            }
        }
    }
}

void Particles::compRiemannStatesLR(const double &dt, const double &kernelSize, const double &gamma,
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

            //xij[0] = x[i] + kernelSize/2. * xjxi[0];
            //xij[1] = y[i] + kernelSize/2. * xjxi[1];


#if FIRST_ORDER_QUAD_POINT
            xij[0] = (x[i] + ghostParticles.x[j])/2.;
            xij[1] = (y[i] + ghostParticles.y[j])/2.;

            xijxj[0] = .5*(x[i] - ghostParticles.x[j]);
            xijxj[1] = .5*(y[i] - ghostParticles.y[j]);

            xijxi[0] = .5*(ghostParticles.x[j] - x[i]);
            xijxi[1] = .5*(ghostParticles.y[j] - y[i]);

#else
            xij[0] = x[i] + kernelSize/4. * xjxi[0];
            xij[1] = y[i] + kernelSize/4. * xjxi[1];
            xijxi[0] = xij[0] - x[i];
            xijxi[1] = xij[1] - y[i];

            xijxj[0] = xij[0] - ghostParticles.x[j];
            xijxj[1] = xij[1] - ghostParticles.y[j];
#endif

#if DIM==3
            xjxi[2] = ghostParticles.z[j] - z[i];
            //xij[2] = z[i] + kernelSize/2. * xjxi[2];
            xij[2] = z[i] + kernelSize/4. * xjxi[2];
            xijxi[2] = xij[2] - z[i];
            xijxj[2] = xij[2] - ghostParticles.z[j];
            // TODO: add FIRST_ORDER_QUAD_POINT
#endif

            int iW = i*MAX_NUM_GHOST_INTERACTIONS+jn;
#if !MOVE_PARTICLES
            vFrameGhosts[iW][0] = 0.;
            vFrameGhosts[iW][1] = 0.;
#if DIM==3
            vFrameGhosts[iW][2] = 0.;
#endif
#else // !MOVE_PARTICLES
#if FIRST_ORDER_QUAD_POINT
            vFrameGhosts[iW][0] = (vx[i] + ghostParticles.vx[j])/2.;
            vFrameGhosts[iW][1] = (vy[i] + ghostParticles.vy[j])/2.;
#if DIM==3
            vFrameGhosts[iW][2] = (vz[i] + ghostParticles.vz[j])/2.;
#endif
#else
            double dotProd = Helper::dotProduct(xijxi, xjxi);
            double dSqr = Helper::dotProduct(xjxi, xjxi);

            vFrameGhosts[iW][0] = vx[i] + (ghostParticles.vx[j]-vx[i]) * dotProd/dSqr;
            vFrameGhosts[iW][1] = vy[i] + (ghostParticles.vy[j]-vy[i]) * dotProd/dSqr;
#if DIM==3
            vFrameGhosts[iW][2] = vz[i] + (ghostParticles.vz[j]-vz[i]) * dotProd/dSqr;
#endif
#endif
#endif // !MOVE_PARTICLES
            /*if(i == 394){
                Logger(DEBUG) << "vFrame[0] = " << vFrameGhosts[iW][0]
                              << ", vFrame[1] = " << vFrameGhosts[iW][1]
                              << ", rhoGrad[i][0] = " << rhoGrad[i][0]
                              << ", rhoGrad[i][1] = " << rhoGrad[i][1]
                              << ", rho[i] = " << rho[i]
                              << ", xijxi = [" << xijxi[0] << ", " << xijxi[1]
                              << "], xij = [" << xij[0] << ", " << xij[1]
                              << "], xj = [" << ghostParticles.x[j] << ", " << ghostParticles.y[j] << "] @" << j;
                //exit(5);
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

            //if (i == 23 && iW == 28){
            //    Logger(DEBUG) << "        j = " << j
            //                  << ", rhoL = " << WijLGhosts[iW][0] << ", rhoR = " << WijRGhosts[iW][0]
            //                  << ", uL = " << WijLGhosts[iW][2] << ", uR = " << WijRGhosts[iW][2]
            //                  << ", PL = " << WijLGhosts[iW][1] << ", PR = " << WijRGhosts[iW][1];
            //}

            // predict half a timestep
            double viDiv = vxGrad[i][0] + vyGrad[i][1];
            double vjDiv = ghostParticles.vxGrad[j][0] + ghostParticles.vyGrad[j][1];
#if DIM==3
            viDiv += vzGrad[i][2];
            vjDiv += ghostParticles.vzGrad[j][2];
#endif
            WijRGhosts[iW][0] -= dt/2. * (rho[i] * viDiv + (vx[i]-vFrameGhosts[iW][0])*rhoGrad[i][0] + (vy[i]-vFrameGhosts[iW][1])*rhoGrad[i][1]);
            WijLGhosts[iW][0] -= dt/2. * (ghostParticles.rho[j] * vjDiv + (ghostParticles.vx[j]-vFrameGhosts[iW][0])*ghostParticles.rhoGrad[j][0] + (ghostParticles.vy[j]-vFrameGhosts[iW][1])*ghostParticles.rhoGrad[j][1]);
            WijRGhosts[iW][1] -= dt/2. * (gamma*P[i] * viDiv + (vx[i]-vFrameGhosts[iW][0])*PGrad[i][0] + (vy[i]-vFrameGhosts[iW][1])*PGrad[i][1]);
            WijLGhosts[iW][1] -= dt/2. * (gamma*ghostParticles.P[j] * vjDiv + (ghostParticles.vx[j]-vFrameGhosts[iW][0])*ghostParticles.PGrad[j][0] + (ghostParticles.vy[j]-vFrameGhosts[iW][1])*ghostParticles.PGrad[j][1]);
            WijRGhosts[iW][2] -= dt/2. * (PGrad[i][0]/rho[i] + (vx[i] - vFrameGhosts[iW][0])*vxGrad[i][0] + (vy[i] - vFrameGhosts[iW][1])*vxGrad[i][1]);
            WijLGhosts[iW][2] -= dt/2. * (ghostParticles.PGrad[j][0]/ghostParticles.rho[j] + (ghostParticles.vx[j]-vFrameGhosts[iW][0])*ghostParticles.vxGrad[j][0] + (ghostParticles.vy[j]-vFrameGhosts[iW][1])*ghostParticles.vxGrad[j][1]);
            WijRGhosts[iW][3] -= dt/2. * (PGrad[i][1]/rho[i] + (vx[i] - vFrameGhosts[iW][0])*vyGrad[i][0] + (vy[i] - vFrameGhosts[iW][1])*vyGrad[i][1]);
            WijLGhosts[iW][3] -= dt/2. * (ghostParticles.PGrad[j][1]/ghostParticles.rho[j] + (ghostParticles.vx[j]-vFrameGhosts[iW][0])*ghostParticles.vyGrad[j][0] + (ghostParticles.vy[j]-vFrameGhosts[iW][1])*ghostParticles.vyGrad[j][1]);
#if DIM==3
            // TODO: update below for 3D
            WijRGhosts[iW][4] -= dt/2. * PGrad[i][2]/rho[i];
            WijLGhosts[iW][4] -= dt/2. * ghostParticles.PGrad[j][2]/ghostParticles.rho[j];
#endif

           //if (i == 46) // && iW == 28){
           //     Logger(DEBUG) << "        j = " << j
           //                   << ", rhoL = " << WijLGhosts[iW][0] << ", rhoR = " << WijRGhosts[iW][0]
           //                   << ", uL = " << WijLGhosts[iW][2] << ", uR = " << WijRGhosts[iW][2]
           //                   << ", PL = " << WijLGhosts[iW][1] << ", PR = " << WijRGhosts[iW][1];
            //}

        }
    }
}

/// debugging function to printout the number of (ghost-)interactions
void Particles::printNoi(){
    for (int i=0; i<N; ++i) {
        Logger(DEBUG) << "          > i = " << i << ", x = [" << x[i] << ", " << y[i] <<
                  "], : noi = " << noi[i] <<
                  " + " << noiGhosts[i] << " = " << noi[i] + noiGhosts[i];
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

        std::vector<int> nnList(noiTot);

        std::vector<std::vector<double>> nnlPrtcls {};
        std::vector<std::vector<double>> nnlAij {};
        std::vector<std::vector<double>> nnl_vFrame {};
        std::vector<std::vector<double>> nnlWijR {};
        std::vector<std::vector<double>> nnlWijL {};

        std::vector<size_t> dataSpaceDims(2);
        dataSpaceDims[0] = std::size_t(noiTot);
        dataSpaceDims[1] = DIM;

        for (int j=0; j<noi[i]; ++j){

            int ii = j+i*MAX_NUM_INTERACTIONS;

            nnList[j] = nnl[ii];

            nnlPrtcls.push_back(std::vector<double>(DIM));
            nnlPrtcls[j][0] = x[nnl[ii]];
            nnlPrtcls[j][1] = y[nnl[ii]];
#if DIM == 3
            nnlPrtcls[j][2] = z[nnl[ii]];
#endif
            nnlAij.push_back(std::vector<double>(DIM));
            nnlAij[j][0] = Aij[ii][0];
            nnlAij[j][1] = Aij[ii][1];
#if DIM == 3
            nnlAij[j][2] = Aij[ii][2];
#endif

            nnl_vFrame.push_back(std::vector<double>(DIM));
            nnl_vFrame[j][0] = vFrame[ii][0];
            nnl_vFrame[j][1] = vFrame[ii][1];
#if DIM == 3
            nnl_vFrame[j][2] = vFrame[ii][2];
#endif
        }

        for (int j=0; j<noiGhosts[i]; ++j){

            int ii = j+i*MAX_NUM_GHOST_INTERACTIONS;

            nnList[j+noi[i]] = nnlGhosts[ii];

            nnlPrtcls.push_back(std::vector<double>(DIM));
            nnlPrtcls[j+noi[i]][0] = ghostParticles.x[nnlGhosts[ii]];
            nnlPrtcls[j+noi[i]][1] = ghostParticles.y[nnlGhosts[ii]];
#if DIM == 3
            nnlPrtcls[j+noi[i]][2] = ghostParticles.z[nnlGhosts[ii]];
#endif
            nnlAij.push_back(std::vector<double>(DIM));
            nnlAij[j+noi[i]][0] = AijGhosts[ii][0];
            nnlAij[j+noi[i]][1] = AijGhosts[ii][1];
#if DIM == 3
            nnlAij[j+noi[i]][2] = AijGhosts[ii][2];
#endif

            nnl_vFrame.push_back(std::vector<double>(DIM));
            nnl_vFrame[j+noi[i]][0] = vFrameGhosts[ii][0];
            nnl_vFrame[j+noi[i]][1] = vFrameGhosts[ii][1];
#if DIM == 3
            nnl_vFrame[j+noi[i]][2] = vFrameGhosts[ii][2];
#endif

        }

        HighFive::DataSet nnListDataSet = h5File.createDataSet<int>("/nnl" +std::to_string(i),
                                                                    HighFive::DataSpace(noiTot));
        nnListDataSet.write(nnList);

        HighFive::DataSet nnlDataSet = h5File.createDataSet<double>("/nnlPrtcls" + std::to_string(i),
                                                                    HighFive::DataSpace(dataSpaceDims));
        nnlDataSet.write(nnlPrtcls);

        HighFive::DataSet AijDataSet = h5File.createDataSet<double>("/Aij" + std::to_string(i),
                                                                    HighFive::DataSpace(dataSpaceDims));
        AijDataSet.write(nnlAij);

        HighFive::DataSet vFrameDataSet = h5File.createDataSet<double>("/vFrame" + std::to_string(i),
                                                                    HighFive::DataSpace(dataSpaceDims));
        vFrameDataSet.write(nnl_vFrame);

    }
}
#endif

// TODO: remove below
//void Particles::move(const double &dt, Domain &domain){
//
//    for(int i=0; i<N; ++i) {
//
//        //Logger(DEBUG) << "v@" << i <<  " = [" << vx[i] << ", " << vy[i] << "]"
//        //            << ", x[n] = [" << x[i] << ", " << y[i] << "]";
//
//        x[i] = x[i] + vx[i] * dt;
//        y[i] = y[i] + vy[i] * dt;
//#if DIM == 3
//        z[i] = z[i] +vz[i] * dt;
//#endif
//#if PERIODIC_BOUNDARIES
//        if (x[i] < domain.bounds.minX) {
//            x[i] = domain.bounds.maxX - (domain.bounds.minX - x[i]);
//        } else if (domain.bounds.maxX <= x[i]) {
//            x[i] = domain.bounds.minX + (x[i] - domain.bounds.maxX);
//        }
//        if (y[i] < domain.bounds.minY) {
//            y[i] = domain.bounds.maxY - (domain.bounds.minY - y[i]);
//        } else if (domain.bounds.maxY <= y[i]) {
//            y[i] = domain.bounds.minY + (y[i] - domain.bounds.maxY);
//        }
//#if DIM ==3
//        if (z[i] < domain.bounds.minZ) {
//            z[i] = domain.bounds.maxZ - (domain.bounds.minZ - z[i]);
//        } else if (domain.bounds.maxZ <= z[i]) {
//            z[i] = domain.bounds.minZ + (z[i] - domain.bounds.maxZ);
//        }
//#endif
//#endif
//        //Logger(DEBUG) << "                               x[n+1] = ["
//        //          << x[i] << ", " << y[i] << "]";
//    }
//}

/// Sanity check functions
double Particles::sumVolume(){
    double V = 0.;
    for (int i=0; i<N; ++i){
        V += 1./omega[i];
    }
    return V;
}

double Particles::sumMass(){
    double M = 0.;
    for (int i=0; i<N; ++i){
        if (std::isnan(m[i])){
            Logger(WARN) << "!! m[" << i <<"] " << "is nan. !!";
        }
        M += m[i];
    }
    return M;
}

double Particles::sumEnergy(){
    double E = 0.;
    for (int i=0; i<N; ++i){
#if DIM == 2
        E += m[i]*(u[i] + .5*(vx[i]*vx[i]+vy[i]*vy[i]));
#else
        E += m[i]*(u[i] + .5*(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]));
#endif
    }
    //std::cout << "E is " <<  E << std::endl;
    return E;
}

double Particles::sumMomentumX(){
    double momX = 0.;
    for (int i=0; i<N; ++i){
        momX += m[i]*vx[i];
    }
    return momX;
}

double Particles::sumMomentumY(){
    double momY = 0.;
    for (int i=0; i<N; ++i){
        momY += m[i]*vy[i];
    }
    return momY;
}

#if DIM == 3
double Particles::sumMomentumZ(){
    double momZ = 0.;
    for (int i=0; i<N; ++i){
        momZ += m[i]*vz[i];
    }
    return momZ;
}
#endif

void Particles::checkFluxSymmetry(Particles *ghostParticles){
    for (int i=0; i<N; ++i){
        for (int j=0; j<noi[i]; ++j){
            int ii = i*MAX_NUM_INTERACTIONS+j; // interaction index i->j

            int ji = nnl[ii]; // index i of particle j
            // search neighbor i in nnl[] of j
            int ij;
            for(ij=0; ij<noi[ji]; ++ij){
                if (nnl[ij+ji*MAX_NUM_INTERACTIONS] == i) break;
            }
            int iij = ji*MAX_NUM_INTERACTIONS+ij; // interaction index j->i

            bool notSym = false;
            if (Fij[iij][0] + Fij[ii][0] > FLUX_SYM_TOL){
                Logger(WARN) << "  > Mass fluxes are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << Fij[ii][0] << ", Fji = " << Fij[iij][0];
                notSym = true;
            }
            if (Fij[iij][1] + Fij[ii][1] > FLUX_SYM_TOL){
                Logger(WARN) << "  > Energy fluxes are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << Fij[ii][1] << ", Fji = " << Fij[iij][1];
                notSym = true;
            }
            if (Fij[iij][2] + Fij[ii][2] > FLUX_SYM_TOL){
                Logger(WARN) << "  > Momentum fluxes (x) are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << Fij[ii][2] << ", Fji = " << Fij[iij][2];
                notSym = true;
            }
            if (Fij[iij][3] + Fij[ii][3] > FLUX_SYM_TOL){
                Logger(WARN) << "  > Momentum fluxes (y) are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << Fij[ii][3] << ", Fji = " << Fij[iij][3];
                notSym = true;
            }
#if DIM == 3
            if (Fij[iij][4] + Fij[ii][4] > FLUX_SYM_TOL){
                Logger(WARN) << "  > Momentum fluxes (z) are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << Fij[ii][4] << ", Fji = " << Fij[iij][4];
                notSym = true;
            }
#endif
            if (notSym){
                Logger(INFO) << "  > Aij = [" << Aij[ii][0] << ", " << Aij[ii][1] << "], Aji = ["
                             << Aij[iij][0] << ", " << Aij[iij][1] << "]";
            }
        }
#if PERIODIC_BOUNDARIES
        for (int j=0; j<noiGhosts[i]; ++j) {
            int ii = i * MAX_NUM_GHOST_INTERACTIONS + j; // interaction index i->j

            int ji = nnlGhosts[ii]; // index i of particle j
            // search neighbor i in nnl[] of j
            int ij;
            // search neighbor i in nnlGhosts[] of j
            for (ij = 0; ij < noiGhosts[ghostParticles->parent[ji]]; ++ij) {
                if (ghostParticles->parent[nnlGhosts[ij + ghostParticles->parent[ji] * MAX_NUM_GHOST_INTERACTIONS]]
                    == i) break;
            }
            int iij = ij + ghostParticles->parent[ji] * MAX_NUM_GHOST_INTERACTIONS; // interaction index j->i

            bool notSym = false;
            if (FijGhosts[iij][0] + FijGhosts[ii][0] > FLUX_SYM_TOL) {
                Logger(WARN) << "  > Ghosts mass fluxes are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << FijGhosts[ii][0] << ", Fji = " << FijGhosts[iij][0];
                notSym = true;
            }
            if (FijGhosts[iij][1] + FijGhosts[ii][1] > FLUX_SYM_TOL) {
                Logger(WARN) << "  > Ghosts energy fluxes are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << FijGhosts[ii][1] << ", Fji = " << FijGhosts[iij][1];
                notSym = true;
            }
            if (FijGhosts[iij][2] + FijGhosts[ii][2] > FLUX_SYM_TOL) {
                Logger(WARN) << "  > Ghosts momentum fluxes (x) are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << FijGhosts[ii][2] << ", Fji = " << FijGhosts[iij][2];
                notSym = true;
            }
            if (FijGhosts[iij][3] + FijGhosts[ii][3] > FLUX_SYM_TOL) {
                Logger(WARN) << "  > Ghosts momentum fluxes (y) are NOT symmetric for " << i << " -> " << j
                             << " => Fij = " << FijGhosts[ii][3] << ", Fji = " << FijGhosts[iij][3];
                notSym = true;
            }
            if (notSym) {
                Logger(INFO) << "  > AijGhost = [" << AijGhosts[ii][0] << ", " << AijGhosts[ii][1] << "], Aji = ["
                             << AijGhosts[iij][0] << ", " << AijGhosts[iij][1] << "]";
            }
        }
#endif
    }
}

void Particles::dump2file(std::string filename, double simTime){
    // open output file
    HighFive::File h5File { filename, HighFive::File::ReadWrite |
                                      HighFive::File::Create |
                                      HighFive::File::Truncate };

    // dimensions for datasets containing vectors
    std::vector<size_t> dataSpaceDims(2);
    dataSpaceDims[0] = std::size_t(N); // number of particles
    dataSpaceDims[1] = DIM;

    // create datasets
    // TODO: Create a h5 object holding all meta data
    HighFive::DataSet timeDataSet = h5File.createDataSet<double>("/time", HighFive::DataSpace(1));
    HighFive::DataSet totalMassDataSet = h5File.createDataSet<double>("/totalMass", HighFive::DataSpace(1));
    HighFive::DataSet energyDataSet = h5File.createDataSet<double>("/energy", HighFive::DataSpace(1));
    HighFive::DataSet xMomDataSet = h5File.createDataSet<double>("/xMomentum", HighFive::DataSpace(1));
    HighFive::DataSet yMomDataSet = h5File.createDataSet<double>("/yMomentum", HighFive::DataSpace(1));
#if DIM == 3
    HighFive::DataSet zMomDataSet = h5File.createDataSet<double>("/zMomentum", HighFive::DataSpace(1));
#endif
    HighFive::DataSet rhoDataSet = h5File.createDataSet<double>("/rho", HighFive::DataSpace(N));
    HighFive::DataSet mDataSet = h5File.createDataSet<double>("/m", HighFive::DataSpace(N));
    HighFive::DataSet uDataSet = h5File.createDataSet<double>("/u", HighFive::DataSpace(N));
    HighFive::DataSet posDataSet = h5File.createDataSet<double>("/x", HighFive::DataSpace(dataSpaceDims));
    HighFive::DataSet velDataSet = h5File.createDataSet<double>("/v", HighFive::DataSpace(dataSpaceDims));
    HighFive::DataSet rhoGradDataSet = h5File.createDataSet<double>("/rhoGrad", HighFive::DataSpace(dataSpaceDims));
    HighFive::DataSet PDataSet = h5File.createDataSet<double>("/P", HighFive::DataSpace(N));
    HighFive::DataSet noiDataSet = h5File.createDataSet<int>("/noi", HighFive::DataSpace(N));


    // containers for particle data
    std::vector<double> timeVec({ simTime });
    std::vector<double> totalMassVec({ sumMass() });
    std::vector<double> energyVec({ sumEnergy() });
    std::vector<double> xMomVec({ sumMomentumX() });
    std::vector<double> yMomVec({ sumMomentumY() });
#if DIM ==3
    std::vector<double> zMomVec({ sumMomentumZ() });
#endif
    std::vector<double> rhoVec(rho, rho+N);
    std::vector<double> mVec(m, m+N);
    std::vector<double> uVec(u, u+N);
    std::vector<double> PVec(P, P+N);
    std::vector<int> noiVec(noi, noi+N);

    std::vector<std::vector<double>> posVec(N);
    std::vector<std::vector<double>> velVec(N);

    std::vector<std::vector<double>> rhoGradVec(N);

    // fill containers with data
    std::vector<double> posBuf(DIM);
    std::vector<double> velBuf(DIM);
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

        // velocity
        velBuf[0] = vx[i];
        velBuf[1] = vy[i];
#if DIM == 3
        velBuf[2] = vz[i];
#endif
        velVec[i] = velBuf;

        // density gradient
        rhoGradBuf[0] = rhoGrad[i][0];
        rhoGradBuf[1] = rhoGrad[i][1];
#if DIM == 3
        rhoGradBuf[2] = rhoGrad[i][2];
#endif
        rhoGradVec[i] = rhoGradBuf;
    }
    // write data
    timeDataSet.write(timeVec); // dummy vec containing one element
    totalMassDataSet.write(totalMassVec);
    energyDataSet.write(energyVec);
    xMomDataSet.write(xMomVec);
    yMomDataSet.write(yMomVec);
#if DIM == 3
    zMomDataSet.write(zMomVec);
#endif
    rhoDataSet.write(rhoVec);
    mDataSet.write(mVec);
    uDataSet.write(uVec);
    PDataSet.write(PVec);
    noiDataSet.write(noi);
    posDataSet.write(posVec);
    velDataSet.write(velVec);
    rhoGradDataSet.write(rhoGradVec);
}
