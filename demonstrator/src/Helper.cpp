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

void Helper::crossProduct(double *a, double *b, double *crossProduct){
#if DIM == 3
    crossProduct[0] = a[1]*b[2] - a[2]*b[1];
    crossProduct[1] = a[2]*b[0] - a[0]*b[2];
    crossProduct[2] = a[0]*b[1] - a[1]*b[0];
#else
    Logger(ERROR) << "Cross product not defined for 2D. - Aborting.";
    exit(9);
#endif
}

void Helper::rotationMatrix2D(double *a, double *b, double *Lambda){
    Lambda[0] = a[0]*b[0] + a[1]*b[1];
    Lambda[1] = -(a[0]*b[1] - a[1]*b[0]);
    //Lambda[2] = a[0]*b[1] - a[1]*b[0];
    Lambda[2] = -Lambda[1];
    Lambda[3] = Lambda[0];
}

#if DIM==3
void Helper::rotationMatrix3D(double *a, double *b, double *Lambda){
    double v[DIM];
    crossProduct(a, b, v);
    double cosAB = dotProduct(a, b); // a and b MUST be normed

    // TODO: cosAB == -1 fails !!!

    double n = 1./(1. + cosAB);

    Lambda[0] = 1.-n*(v[2]*v[2]+v[1]*v[1]);
    Lambda[1] = -v[2]+n*v[0]*v[1];
    Lambda[2] = v[1]+n*v[0]*v[2];
    Lambda[3] = v[2]+n*v[0]*v[1];
    Lambda[4] = 1.-n*(v[2]*v[2]+v[0]*v[0]);
    Lambda[5] = -v[0]+n*v[1]*v[2];
    Lambda[6] = -v[1]+n*v[0]*v[2];
    Lambda[7] = v[0]+n*v[1]*v[2];
    Lambda[8] = 1.-n*(v[1]*v[1]+v[0]*v[0]);

    // check for aligned vectors
    //double cosAB = ab/(sqrt(dotProduct(a, a))*sqrt(dotProduct(b, b)));
    //if (abs(cosAB) < 1. + ROT_3D_ALIGN_TOL && abs(cosAB) > 1. - ROT_3D_ALIGN_TOL){
    //    Logger(WARN) << "Very small angle between a and b. - Check Lambda!!";
    //    Logger(DEBUG) << "a = [" << a[0] << ", " << a[1] << ", " << a[2] << "]";
    //    Logger(DEBUG) << "b = [" << b[0] << ", " << b[1] << ", " << b[2] << "]";
    //    Logger(DEBUG) << "Lambda = [" << Lambda[0] << ", " << Lambda[1] << ", " << Lambda[2];
    //    Logger(DEBUG) << "          " << Lambda[3] << ", " << Lambda[4] << ", " << Lambda[5];
    //    Logger(DEBUG) << "          " << Lambda[6] << ", " << Lambda[7] << ", " << Lambda[8] << "]";
    //}
}
#endif
