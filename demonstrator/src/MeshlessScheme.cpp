//
// Created by Johannes Martin on 21.09.22.
//

#include "../include/MeshlessScheme.h"


MeshlessScheme::MeshlessScheme(Configuration config, Particles *particles,
                               Domain::Cell bounds) : config { config }, particles { particles },
                                                      domain(bounds)
#if PERIODIC_BOUNDARIES
                                                      , ghostParticles(DIM*particles->N, true) // TODO: memory optimization
#endif
                                                      {

    Logger(INFO) << "    > Creating grid ... ";
    domain.createGrid(config.kernelSize);
    Logger(INFO) << "    > ... got " << domain.numGridCells << " cells";
}

void MeshlessScheme::run(){
    double t = 0;
    int step = 0;

#if ADAPTIVE_TIMESTEP
    int numDumpTimes = (int)(config.timeEnd/config.timeStep)/config.h5DumpInterval+1;
    Logger(DEBUG) << "      > Times for file dump: " << numDumpTimes;
    double dumpTimes[numDumpTimes];
    for(int iDump=0; iDump<numDumpTimes; ++iDump){
        dumpTimes[iDump] = iDump*config.timeStep*config.h5DumpInterval;
        Logger(DEBUG) << "        dumpTimes[" << iDump << "] = " << dumpTimes[iDump];
    }
    int dumpStep = 0;
    bool dump = true;
    bool dumpNext = false;
#endif // ADAPTIVE_TIMESTEP


    do {
        Logger(INFO) << "  > TIME: " << t << ", STEP: " << step;
#if !PERIODIC_BOUNDARIES
        Logger(INFO) << "    > Computing domain limits ...";
        double domainLimits[DIM*2];
        particles->getDomainLimits(domainLimits);
        Domain::Cell boundingBox { domainLimits };
        domain.bounds = boundingBox;
        domain.printout();
        Logger(DEBUG) << "      > ... creating grid ...";
        domain.createGrid(config.kernelSize);
        Logger(INFO) << "    > ... done.";
#endif
        Logger(INFO) << "    > Assigning particles ...";
        particles->assignParticlesAndCells(domain);
        Logger(INFO) << "    > ... done.";
#if PERIODIC_BOUNDARIES
        Logger(INFO) << "    > Creating ghost particles ...";
        //Logger(DEBUG) << "      > Creating ghost grid";
        //domain.createGhostGrid();
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
        particles->compDensity(config.kernelSize);
#if PERIODIC_BOUNDARIES
        particles->compDensity(ghostParticles, config.kernelSize);
#endif

        Logger(INFO) << "    > Computing pressure";
        particles->compPressure(config.gamma);
        //particles->printDensity(config.gamma);

        Logger(DEBUG) << "      SANITY CHECK > V_tot = " << particles->sumVolume();
        Logger(DEBUG) << "      SANITY CHECK > M_tot = " << particles->sumMass();
        Logger(DEBUG) << "      SANITY CHECK > E_tot = " << particles->sumEnergy();
        Logger(DEBUG) << "      SANITY CHECK > px_tot = " << particles->sumMomentumX();
        Logger(DEBUG) << "      SANITY CHECK > py_tot = " << particles->sumMomentumY();
#if DIM == 3
        Logger(DEBUG) << "      SANITY CHECK > pz_tot = " << particles->sumMomentumZ();
#endif

#if ADAPTIVE_TIMESTEP
        Logger(INFO) << "    > Selecting global timestep ... ";
        timeStep = particles->compGlobalTimestep(config.gamma, config.kernelSize);
        //Logger(INFO) << "Time  > dt = " << timeStep << " selected.";
        if(dumpStep >= numDumpTimes){
            Logger(ERROR) << "Simulation did not abort after reaching timeEnd. Exiting.";
            exit(9);
        } else if(t+timeStep>=dumpTimes[dumpStep+1]){
            dumpNext = true;
            timeStep = dumpTimes[dumpStep+1]-t;
        }
        Logger(INFO) << "Time  > dt = " << timeStep << " selected.";
#else // ADAPTIVE_TIMESTEP
        timeStep = config.timeStep;
#endif // ADAPTIVE_TIMESTEP

        Logger(INFO) << "    > Computing gradients";
#if PERIODIC_BOUNDARIES
        particles->updateGhostState(ghostParticles);
        particles->compPsijTilde(helper, ghostParticles, config.kernelSize);
        //Logger(DEBUG) << "      > Update ghost psij_xiTilde";
        //particles->updateGhostPsijTilde(ghostParticles);

        particles->gradient(particles->rho, particles->rhoGrad, ghostParticles.rho, ghostParticles);
        particles->gradient(particles->vx, particles->vxGrad, ghostParticles.vx, ghostParticles);
        particles->gradient(particles->vy, particles->vyGrad, ghostParticles.vy, ghostParticles);
#if DIM == 3
        particles->gradient(particles->vz, particles->vzGrad, ghostParticles.vz, ghostParticles);
#endif
        particles->gradient(particles->P, particles->PGrad, ghostParticles.P, ghostParticles);
        Logger(DEBUG) << "      > Update ghost gradients";
        particles->updateGhostGradients(ghostParticles);

#if SLOPE_LIMITING
        // TODO: Check slope limiter
        Logger(DEBUG) << "      > Limiting slopes";
        particles->slopeLimiter(config.kernelSize, &ghostParticles);
        Logger(DEBUG) << "      > Update limited ghost gradients";
        particles->updateGhostGradients(ghostParticles);
#endif // SLOPE_LIMITING
#else // PERIODIC_BOUNDARIES
        particles->compPsijTilde(helper, config.kernelSize);
        particles->gradient(particles->rho, particles->rhoGrad);
        particles->gradient(particles->vx, particles->vxGrad);
        particles->gradient(particles->vy, particles->vyGrad);
#if DIM == 3
        particles->gradient(particles->vz, particles->vzGrad);
#endif
        particles->gradient(particles->P, particles->PGrad);
#if SLOPE_LIMITING
        // TODO: check how to properly limit gradiens
        Logger(DEBUG) << "      > Limiting slopes";
        particles->slopeLimiter(config.kernelSize);
#endif // SLOPE_LIMITING
#endif
        Logger(INFO) << "    > Preparing Riemann solver";
        Logger(DEBUG) << "      > Computing effective faces";
        particles->compEffectiveFace();
#if PERIODIC_BOUNDARIES
        particles->compEffectiveFace(ghostParticles);
#endif // PERIODIC_BOUNDARIES
        Logger(DEBUG) << "      > Computing fluxes";
        particles->compRiemannStatesLR(timeStep, config.kernelSize, config.gamma);

#if PERIODIC_BOUNDARIES
        Logger(DEBUG) << "      > Computing ghost fluxes";
        particles->compRiemannStatesLR(timeStep, config.kernelSize, config.gamma,
                                     ghostParticles);
        //Logger(DEBUG) << "Aborting for debugging.";
        //exit(6);

#endif// PERIODIC_BOUNDARIES

#if ADAPTIVE_TIMESTEP
        if (dump){
            dump = false;

#else // ADAPTIVE_TIMESTEP
        if (step % config.h5DumpInterval == 0) {
#endif // ADAPTIVE_TIMESTEP

            std::stringstream stepss;
            Logger(INFO) << "   > Dump particle distribution";

            stepss << std::setw(6) << std::setfill('0')
#if ADAPTIVE_TIMESTEP
            << dumpStep;
#else
            << step;
#endif // ADAPTIVE_TIMESTEP

            Logger(INFO) << "      > Dump particles to file";
            particles->dump2file(config.outDir + "/" + stepss.str() + std::string(".h5"), t);

            ++dumpStep;

#if DEBUG_LVL > 1
#if PERIODIC_BOUNDARIES
            Logger(INFO) << "      > Dump ghosts to file";
            ghostParticles.dump2file(config.outDir + "/" + stepss.str() + std::string("Ghosts.h5"), t);
            Logger(INFO) << "      > Dump NNL to file";
            particles->dumpNNL(config.outDir + "/" + stepss.str() + std::string("NNL.h5"), ghostParticles);
#endif // PERIODIC_BOUNDARIES
#endif // DEBUG_LVL
        }
        if (t>=config.timeEnd){
            Logger(INFO) << "    > t = " << t << " -> FINISHED!";
            break;
        }


        Logger(INFO) << "    > Solving Riemann problems";
#if PERIODIC_BOUNDARIES
        particles->solveRiemannProblems(config.gamma, ghostParticles);
#else
        Particles ghostParticles { 0, true }; // DUMMY
        particles->solveRiemannProblems(config.gamma, ghostParticles);
#endif

#if DEBUG_LVL
        Logger(DEBUG) << "    > Checking flux symmetry";
#if PERIODIC_BOUNDARIES
        particles->checkFluxSymmetry(&ghostParticles);
#else
        particles->checkFluxSymmetry();
#endif
#endif

        Logger(INFO) << "    > Collecting fluxes";

        particles->collectFluxes(helper, ghostParticles);

        Logger(INFO) << "    > Updating state";
        particles->updateStateAndPosition(timeStep, domain);

        //Logger(INFO) << "    > Moving particles";
        //particles->move(config.timeStep, domain);

        t += timeStep;
        ++step;

#if ADAPTIVE_TIMESTEP
        if (dumpNext){
            dump = true;
            dumpNext = false;
        }
#endif // ADAPTIVE_TIMESTEP

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

    } while(t<config.timeEnd+timeStep);
}

MeshlessScheme::~MeshlessScheme(){}
