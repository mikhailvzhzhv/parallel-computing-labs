#include <omp.h>
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
    #pragma omp for schedule(runtime)
    for (int i = 0; i < N; ++i) {
        double xn = 0;
        for (int j = 0; j < N; ++j) {
            xn += A[i * N + j] * x[j];
        }
        Axn[i] = xn - b[i];
    }

    #pragma omp for schedule(runtime)
    for (int i = 0; i < N; ++i) {
        x[i] -= t * Axn[i];
    }
}

double getSquareLen(const double* vector, const int N) {
    double len = 0;
    #pragma omp parallel for schedule(runtime) reduction(+:len)
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

    double* x = new double[N];
    double* b = new double[N];
    double* A = new double[N * N];
    double* Axn = new double[N];

    double lenAxnMinusB;

    fillXVector(x, N);
    fillMatrix(A, N);
    fillBVector(b, N);

    const double squareEpsilonMulSquareLenB = getSquareLen(b, N) * epsilon * epsilon;

    int i = 0;
    int done = 0;
    int exitFlag = 0;
    
    double start = omp_get_wtime();

    #pragma omp parallel shared(done, i, exitFlag, lenAxnMinusB)
    {
        do {
            f(x, A, Axn, b, N, t);
            #pragma omp for schedule(runtime) reduction(+:lenAxnMinusB)
            for (int i = 0; i < N; ++i) {
                lenAxnMinusB += Axn[i] * Axn[i];
            }

            #pragma omp single
            {
                if (lenAxnMinusB < squareEpsilonMulSquareLenB) {
                    ++exitFlag;
                } else {
                    exitFlag = 0;
                }
                done = (exitFlag == 3) || (i++ >= maxIters);
                lenAxnMinusB = 0;
            }
        } while (!done);
    }

    printf("time : %lf\niters : %d\n", omp_get_wtime() - start, i);

    delete[] x;
    delete[] b;
    delete[] A;
    delete[] Axn;

    return 0;
}
