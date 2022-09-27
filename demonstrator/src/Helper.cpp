//
// Created by Johannes Martin on 26.09.22.
//

#include "../include/Helper.h"

void Helper::inverseMatrix(double *A, int N){
    //int *IPIV = new int[N];
    int LWORK = N*N;
    //double *WORK = new double[LWORK];
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    //delete[] IPIV;
    //delete[] WORK;
}

double Helper::dotProduct(double *a, double *b){
    double res = 0.;
    for (int alpha=0; alpha<DIM; ++alpha){
        res += a[alpha] * b[alpha];
    }
    return res;
}