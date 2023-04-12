//
// Created by Johannes Martin on 17.09.21.
//

// header only libraries
#include <cxxopts.hpp>

#include "../include/Logger.h"
#include "../include/ConfigParser.h"
#include "../include/MeshlessScheme.h"
#include "../include/SPH.h"


structlog LOGCFG = {};

int main(int argc, char *argv[]){

    cxxopts::Options cmdLineOptions { "mlh",
                                      "Demonstrator for the meshless hydrodynamic simulation methods MFV and MFM. Also does SPH with smoothed gradient" };
    cmdLineOptions.add_options()
            ("c,config", "Path to config file", cxxopts::value<std::string>()->default_value("config.info"))
            ("v,verbose", "More printouts for debugging")
            ("s,silent", "Suppress normal printouts")
            ("h,help", "Show this help");

    auto cmdLineOpts = cmdLineOptions.parse(argc, argv);

    if (cmdLineOpts.count("help")) {
        std::cout << cmdLineOptions.help() << std::endl;
        exit(0);
    }

    ConfigParser confP { cmdLineOpts["config"].as<std::string>() };

    // initialize Logger
    LOGCFG.headers = true;
    LOGCFG.level = cmdLineOpts.count("verbose") ? DEBUG : INFO;

    if (cmdLineOpts.count("silent")){
        if(cmdLineOpts.count("verbose")){
            throw std::invalid_argument("Command line options -s and -v are incompatible");
        } else {
            LOGCFG.level = WARN;
        }
    }


    Logger(INFO) << "Reading configuration ... ";
#if RUNSPH
    SPH::Configuration config;
#else
    MeshlessScheme::Configuration config;
#endif

    config.initFile = confP.getVal<std::string>("initFile");
    Logger(INFO) << "    > Initial distribution: " << config.initFile;
    config.outDir = confP.getVal<std::string>("outDir");
    Logger(INFO) << "    > Output directory: " << config.outDir;
    config.timeStep = confP.getVal<double>("timeStep");
    Logger(INFO) << "    > Time step: " << config.timeStep;
    config.timeEnd = confP.getVal<double>("timeEnd");
    Logger(INFO) << "    > End of simulation: " << config.timeEnd;
    config.h5DumpInterval = confP.getVal<int>("h5DumpInterval");
    Logger(INFO) << "    > Dump data to h5 file every " << config.h5DumpInterval << " steps";
    config.kernelSize = confP.getVal<double>("kernelSize");
    Logger(INFO) << "    > Using global kernel size h = " << config.kernelSize;
    config.gamma = confP.getVal<double>("gamma");
    Logger(INFO) << "    > Adiabatic index for ideal gas EOS gamma = " << config.gamma;
#if PERIODIC_BOUNDARIES
    auto periodicBoxLimits = confP.getObj("periodicBoxLimits");
    config.periodicBoxLimits[0] = periodicBoxLimits.getVal<double>("lowerX");
    config.periodicBoxLimits[DIM] = periodicBoxLimits.getVal<double>("upperX");
    config.periodicBoxLimits[1] = periodicBoxLimits.getVal<double>("lowerY");
    config.periodicBoxLimits[DIM+1] = periodicBoxLimits.getVal<double>("upperY");
#if DIM == 3
    config.periodicBoxLimits[2] = periodicBoxLimits.getVal<double>("lowerZ");
    config.periodicBoxLimits[DIM+2] = periodicBoxLimits.getVal<double>("upperZ");
#endif
    std::string periodicBoxStr = "[";
    for (int i=0; i<2*DIM; i++){
        periodicBoxStr.append(std::to_string(config.periodicBoxLimits[i]));
        if(i<2*DIM-1) periodicBoxStr.append(", ");
    }
    Logger(INFO) << "    > Periodic boundaries within box: " << periodicBoxStr << "]";
#endif

    Logger(INFO) << "    > Reading initial distribution ...";

    InitialDistribution initDist { config.initFile };
    Particles particles { initDist.getNumberOfParticles() };
    initDist.getAllParticles(particles);

    Logger(INFO) << "    > N = " << particles.N;
    Logger(INFO) << "... done. Initializing simulation ...";

#if PERIODIC_BOUNDARIES
    double *domainLimits = config.periodicBoxLimits;
#else
    double domainLimits[DIM*2];
    particles.getDomainLimits(domainLimits);
#endif
    Domain::Cell boundingBox { domainLimits };


#if RUNSPH
    SPH algorithm {config, &particles, boundingBox};
#else

    MeshlessScheme algorithm { config, &particles, boundingBox };
#endif

    Logger(INFO) << "... done.";

    Logger(INFO) << "Starting time integration ...";
    algorithm.run();
    Logger(INFO) << "... done.";
    return 0;
}
