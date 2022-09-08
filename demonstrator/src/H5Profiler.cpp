//
// Created by Johannes Martin on 01.06.21.
//

#include "../include/H5Profiler.h"

H5Profiler::H5Profiler(const std::string& outfile, int _myRank, int _numProcs) :
    h5file { outfile,HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Truncate,
             HighFive::MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL) },
    myRank { _myRank }, numProcs { _numProcs }, step { 0 }, disabled { false }
{
    Logger(DEBUG) << "H5Profiler instance created (Singleton)";
}

// TIMING FUNCTIONS
void H5Profiler::createTimeDataSet(const std::string& path, int steps){
    dataSets[path] = h5file.createDataSet<double>(path, HighFive::DataSpace({std::size_t(steps),
                                                                                 std::size_t(numProcs)}));
    timeElapsed[path] = 0.;
}

void H5Profiler::time(const std::string &path){
    timeStart[path] = MPI_Wtime();
}

void H5Profiler::timePause(const std::string& path){
    timeElapsed[path] += MPI_Wtime() - timeStart[path];
}

void H5Profiler::time2file(const std::string& path, bool onlyWrite){
    if (!onlyWrite) timePause(path);
    if (!disabled) dataSets[path].select({std::size_t(step), std::size_t(myRank)}, {1, 1}).write(timeElapsed[path]);
    timeElapsed[path] = 0.;
}



