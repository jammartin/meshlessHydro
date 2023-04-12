//
// Created by Johannes Martin on 21.09.22.
//


#include "../include/SPH.h"


SPH::SPH(Configuration config, Particles *particles,
                               Domain::Cell bounds) : config { config }, particles { particles },
                                                      domain(bounds)
#if PERIODIC_BOUNDARIES
                                                      , ghostParticles(DIM*particles->N, true) // TODO: memory optimization
#endif
                                                      {

    Logger(INFO) << "    > Creating grid ... ";
    domain.createGrid(config.kernelSize);
    Logger(INFO) << "    > ... got " << domain.numGridCells << " cells";
    Logger(INFO) << "    > Running SPH with " << particles->N << " Particles";
}

void SPH::run(){
#if RUNSPH
    double t = 0;
    int step = 0;
    do {
        Logger(INFO) << "  > TIME: " << t << ", STEP: " << step;
        Logger(INFO) << "    > Assigning particles ...";
        particles->assignParticlesAndCells(domain);
        Logger(INFO) << "    > ... done.";
#if PERIODIC_BOUNDARIES
        Logger(INFO) << "    > Creating ghost particles ...";
        // Logger(DEBUG) << "      > Creating ghost grid";
        // domain.createGhostGrid();
        Logger(DEBUG) << "      > Creating ghost particles ... ";
        particles->createGhostParticles(domain, ghostParticles, config.kernelSize);
        Logger(DEBUG) << "      > ... found " << ghostParticles.N << " ghosts";
        Logger(INFO) << "    > ... done.";

#endif
        Logger(INFO) << "    > Nearest neighbor search";
        particles->gridNNS(domain, config.kernelSize);
#if PERIODIC_BOUNDARIES
        Logger(DEBUG) << "      > Ghosts NNS";
        particles->ghostNNS(domain, ghostParticles, config.kernelSize);
        //particles->printNoi();
#endif
        Logger(INFO) << "    > Computing density";
#if PERIODIC_BOUNDARIES
        particles->compDensitySPH(ghostParticles, config.kernelSize);
#else
        particles->compDensitySPH(config.kernelSize);
#endif
#if PERIODIC_BOUNDARIES
        particles->updateGhostState(ghostParticles);
#endif

        if (step == 0){
            Logger(INFO) << "   > Setting u so that P = 2.5";
            particles->setInternalEnergy(2.5, config.gamma);
#if PERIODIC_BOUNDARIES
            particles->updateGhostState(ghostParticles);
#endif
        }
        Logger(INFO) << "    > Computing pressure";
        particles->compPressure(config.gamma);
        // //particles->printDensity(config.gamma);
        // Logger(DEBUG) << "      SANITY CHECK > M_tot = " << particles->sumMass();
        // Logger(DEBUG) << "      SANITY CHECK > E_tot = " << particles->sumEnergy();
        // Logger(DEBUG) << "      SANITY CHECK > px_tot = " << particles->sumMomentumX();
        // Logger(DEBUG) << "      SANITY CHECK > py_tot = " << particles->sumMomentumY();

        // Logger(DEBUG) << "  > Some densities";
        // std::cout << particles->rho[2] << " " << particles->rho[27] << " " << particles->rho[4] << " " << particles->rho[25];
#if PERIODIC_BOUNDARIES
        particles->updateGhostState(ghostParticles);
#endif

#if SLOPE_LIMITING
        // TODO: Check slope limiter
        Logger(DEBUG) << "      > Limiting slopes";
        particles->slopeLimiter(config.kernelSize, &ghostParticles);
        Logger(DEBUG) << "      > Update limited ghost gradients";
        particles->updateGhostGradients(ghostParticles);
#endif


#if ARTVISC
Logger(INFO) << "   > Computing speed of sound";
particles->compCs(config.gamma);
        // Logger(INFO) << "   > Computing artificial viscocity for acceleration";
        // particles->compAccArtVisc(config.kernelSize);
// #if PERIODIC_BOUNDARIES
//         particles->compAccArtVisc(ghostParticles, config.kernelSize);
// #endif
#endif

#if PERIODIC_BOUNDARIES
Logger(DEBUG) << "	> Computing acceleration";
particles->compAccSPH(ghostParticles, config.kernelSize);
#else
Logger(DEBUG) << "	> Computing acceleration";
particles->compAccSPH(config.kernelSize);
#endif

#if PERIODIC_BOUNDARIES
        particles->updateGhostState(ghostParticles);
#endif

#if SLOPE_LIMITING
        // TODO: check how to properly limit gradiens
        Logger(DEBUG) << "      > Limiting slopes";
        particles->slopeLimiter(config.kernelSize);
#endif

        Logger(INFO) << "    > Euler integration";

        Logger(DEBUG) << "Old Energy:";
        particles->sumEnergy();

        particles->eulerSPH(config.timeStep, domain);

        Logger(DEBUG) << "New Energy";
        particles->sumEnergy();

         // ENERGY
         Logger(DEBUG) << " 	> Computing omega w/o ghost particcles:";
         particles->compOmegas(config.kernelSize); // Works

 #if PERIODIC_BOUNDARIES
        Logger(DEBUG) << " 	> Computing omega with ghostParticles:";
        particles->compOmegas(ghostParticles, config.kernelSize);

        Logger(DEBUG) << "	> Computing dn/dh, drho/dh";
        particles->calcdndrho(ghostParticles, config.kernelSize);

         //Logger(DEBUG) << "	> Computing dP/dt";
         //particles->calcdP(ghostParticles, config.kernelSize);

        Logger(DEBUG) << "      > Computing dE/dt";

        particles->calcdE(ghostParticles, config.kernelSize);
 #else
         particles->calcdndrho(config.kernelSize);
         //particles->calcdP(config.kernelSize);
         particles->calcdE(config.kernelSize);
 #endif

        particles->sumEnergy();
#if ARTVISC
        Logger(DEBUG) << "  > Computing artificial viscocity for energy";
        particles->compUiArtVisc(config.kernelSize);
#if PERIODIC_BOUNDARIES
        particles->compUiArtVisc(ghostParticles, config.kernelSize);
#endif
#endif

        Logger(DEBUG) << "Updating u[i]. Printing energy before and after";
        particles->compuis(config.timeStep, config.kernelSize);
        particles->sumEnergy();

// // #if PERIODIC_BOUNDARIES
// //         particles->updateGhostState(ghostParticles);
// // #endif

        if (step % config.h5DumpInterval == 0) {
            std::stringstream stepss;
            Logger(INFO) << "   > Dump particle distribution";
            stepss << std::setw(6) << std::setfill('0') << step;
            Logger(INFO) << "      > Dump particles to file";
            particles->dump2file(config.outDir + "/" + stepss.str() + std::string(".h5"));

#if DEBUG_LVL > 1
#if PERIODIC_BOUNDARIES
            Logger(INFO) << "      > Dump ghosts to file";
            ghostParticles.dump2file(config.outDir + "/" + stepss.str() + std::string("Ghosts.h5"));
            Logger(INFO) << "      > Dump NNL to file";
            particles->dumpNNL(config.outDir + "/" + stepss.str() + std::string("NNL.h5"), ghostParticles);
#endif
#endif
        }

        if (t>=config.timeEnd){
            Logger(INFO) << "    > t = " << t << " -> FINISHED!";
            break;
        }
        //Logger(DEBUG) << "      SANITY CHECK > V_tot = " << particles->sumVolume();
        // Logger(DEBUG) << "      SANITY CHECK > M_tot = " << particles->sumMass();
        // Logger(DEBUG) << "      SANITY CHECK > E_tot = " << particles->sumEnergy();
        // Logger(DEBUG) << "      SANITY CHECK > px_tot = " << particles->sumMomentumX();
        // Logger(DEBUG) << "      SANITY CHECK > py_tot = " << particles->sumMomentumY();
        //




        //Logger(INFO) << "    > Moving particles";
        //particles->move(config.timeStep, domain);

        t += config.timeStep;
        ++step;
        // if (step == 10){
	    //  	break;
        // }
        //break;
        //Logger(DEBUG) << "    > t = " << t << ", step =  " << step
        //          << ", t_end = " << config.timeEnd;

        // DEBUGGING
        // TODO: remove
        //Logger(DEBUG) << "      SANITY CHECK > M_tot = " << particles->sumMass();

        //stepss = std::stringstream();;
        //Logger(INFO) << "   > Dump particle distribution";
        //stepss << std::setw(6) << std::setfill('0') << step;
        //Logger(INFO) << "      > Dump particles to file";
        //particles->dump2file(config.outDir + "/" + stepss.str() + std::string(".h5"));
        // END DEBUGGING
    } while(t<config.timeEnd+config.timeStep);
#endif
}

SPH::~SPH(){}
