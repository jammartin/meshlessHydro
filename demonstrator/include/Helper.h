//
// Created by Johannes Martin on 26.09.22.
//

#ifndef MESHLESSHYDRO_HELPER_H
#define MESHLESSHYDRO_HELPER_H

#include "parameter.h"

/// Matrix inversion taken from: https://stackoverflow.com/a/3525136/6208997

extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}

class Helper {

public:
    void inverseMatrix(double *A, int N);

private:
    double WORK[DIM*DIM];
    int IPIV[DIM];
};

#endif //MESHLESSHYDRO_HELPER_H
