//include "mpi.h"
#include <stdlib.h>
#include <stdio.h>


void fillXVector(double* x, const int N) {
    for (int i = 0; i < N; ++i) {
        x[i] = 0;
    }
}

void fillMatrix(double* A, const int N) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i * N + j] = 1.0;
        }
        A[i * N + i] = 2.0;
    }
}

void fillBVector(double* b, const int N) {
    for (int i = 0; i < N; ++i) {
        b[i] = N + 1;
    }
}

void f(double* x, const double* A, double* Axn, const double* b, const int N, const double t) {
    for (int i = 0; i < N; ++i) {
        double xn = 0;
        for (int j = 0; j < N; ++j) {
            xn += A[i * N + j] * x[j];
        }
        Axn[i] = xn - b[i];
    }
    for (int i = 0; i < N; ++i) {
        x[i] -= t*(Axn[i]);
    }
}

double getSquareLen(const double* vector, const int N) {
    double len = 0;
    for (int i = 0; i < N; ++i) {
        len += vector[i] * vector[i];
    }

    return len;
}


int main(int argc, char* argv[]) {

    const double epsilon = 1e-4;
    const double t = 1e-6;
    const int maxIters = 10000;
    const int N = 1440;

//    MPI_Init(&argc,&argv);

//    double start = MPI_Wtime();

    double* x = new double[N];
    double* b = new double[N];
    double* A = new double[N * N];
    double* Axn = new double[N];
    double lenAxnMinusB;

    fillXVector(x, N);
    fillMatrix(A, N);
    fillBVector(b, N);

    const double squareEpsilonMulSquareLenB = getSquareLen(b, N) * epsilon * epsilon * 3;
    double lenAxnMinusB1 = __DBL_MAX__;
    double lenAxnMinusB2 = __DBL_MAX__;
    double lenAxnMinusB3 = __DBL_MAX__;
    

    int i = 0;
    do {
        f(x, A, Axn, b, N, t);
        lenAxnMinusB1 = lenAxnMinusB2;
        lenAxnMinusB2 = lenAxnMinusB3;
        lenAxnMinusB3 = getSquareLen(Axn, N);
        lenAxnMinusB = lenAxnMinusB1 + lenAxnMinusB2 + lenAxnMinusB3;
    } while (lenAxnMinusB >= squareEpsilonMulSquareLenB && i++ < maxIters);

//    printf("time : %lf\niters : %d\n", MPI_Wtime() - start, i);

    delete[] x;
    delete[] b;
    delete[] A;
    delete[] Axn;

    return 0;
}
