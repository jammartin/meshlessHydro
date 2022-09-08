//
// Created by Johannes Martin on 01.06.21.
//

#ifndef MESHLESSHYDRO_H5PROFILER_H
#define MESHLESSHYDRO_H5PROFILER_H

#include <string>
#include <unordered_map>
#include <vector>

#include <mpi.h>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

#include "Logger.h"

namespace ProfilerIds {
    const char* const N { "/general/numberOfParticles" };
    const char* const timePos { "/time/compPosition" };
    const char* const timeMove { "/time/moveParticles" };
    const char* const timePseudo { "/time/compPseudoParticles" };
    const char* const timeForce { "/time/compForce" };
    const char* const timeVel { "/time/compVelocity" };
    // only for parallel mode
    const char* const timeCommonCoarse { "/time/commonCoarsePseudoParticles" };
    const char* const timeForceExchange { "/time/forceParticleExchange" };
    const char* const forceRcv { "/force/receiveLength" };
    const char* const timeLb { "/time/loadBalancing" };
    const char* const lbRcv { "/loadBalancing/receiveLength" };

}

// Singleton
class H5Profiler {
public:
    // constructor is only called ONCE due to the singleton pattern
    // arguments have to be passed when first calling getInstance()
    static H5Profiler& getInstance(const std::string& outfile = "", int _myRank = -1, int _numProcs = -1){
        static H5Profiler instance_(outfile, _myRank, _numProcs);
        return instance_;
    }

    // deleting methods we don't want for the singleton design pattern (c++11)
    H5Profiler(H5Profiler const&) = delete;
    void operator=(H5Profiler const&) = delete;

    // setter
    void setStep(const int &_step) { step = _step; }
    void disableWrite() { disabled = true; }
    void enableWrite() { disabled = false; }

    void createTimeDataSet(const std::string& path, int steps);
    void time(const std::string &path);
    void timePause(const std::string &path);
    void time2file(const std::string &path, bool onlyWrite=false);

    // TEMPLATE FUNCTIONS

    // track single values
    template<typename T>
    void createValueDataSet(const std::string& path, int steps){
        dataSets[path] = h5file.createDataSet<T>(path, HighFive::DataSpace({std::size_t(steps),
                                                                                 std::size_t(numProcs)}));
    }

    template<typename T>
    void value2file(const std::string& path, T value){
        if (!disabled) dataSets[path].select({std::size_t(step), std::size_t(myRank)}, {1, 1}).write(value);
    }

    // track vector values
    template<typename T>
    void createVectorDataSet(const std::string& path, int steps, std::size_t size){
        dataSets[path] = h5file.createDataSet<T>(path, HighFive::DataSpace({std::size_t(steps),
                                                                            std::size_t(numProcs), size}));
        vectorSizes[path] = size;
    }

    template<typename T>
    void vector2file(const std::string& path, std::vector<T> data){
        if (!disabled) dataSets[path].select({std::size_t(step), std::size_t(myRank), 0},
                                             {1, 1, vectorSizes[path]}).write(data);
    }

private:
    H5Profiler(const std::string& outfile, int _myRank, int _numProcs);

    // basic containers for meta information
    HighFive::File h5file;
    std::unordered_map<std::string, HighFive::DataSet> dataSets {};
    int myRank;
    int numProcs;
    int step;
    bool disabled;

    // timing variables
    std::unordered_map<std::string, double> timeStart {};
    std::unordered_map<std::string, double> timeElapsed {};

    // vector sizes
    std::unordered_map<std::string, std::size_t> vectorSizes {};

};


#endif //MESHLESSHYDRO_H5PROFILER_H
