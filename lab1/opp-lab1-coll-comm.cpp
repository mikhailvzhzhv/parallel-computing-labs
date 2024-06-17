#include "mpi.h"
#include <stdlib.h>
#include <stdio.h>

const int N = 384*250;

void calculatePartSum(int* a, int* b, int* partSum, int size) {
    for (int i = 0; i < N / size; ++i) {
        for (int j = 0; j < N; ++j) {
            *partSum += a[i] * b[j];
        }
    }
}

int main(int argc, char* argv[]) {

    int size, rank;
    MPI_Init(&argc,&argv);

    double start = MPI_Wtime();

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int* a = (int*) malloc(N * sizeof(int));
    int* b = (int*) malloc(N * sizeof(int));
    int sum = 0;
    int partSum = 0;

    if (rank == 0) {
        for (int i = 0; i < N; ++i) {
            a[i] = 1 + rand() % 10;
            b[i] = 1 + rand() % 10;
        }
    }

    MPI_Bcast(b, N, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(a, N / size, MPI_INT, a, N, MPI_INT, 0, MPI_COMM_WORLD);

    calculatePartSum(a, b, &partSum, size);

    MPI_Reduce(&partSum, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("sum = %d; time = %lf\n", sum, MPI_Wtime() - start);
    }

    MPI_Finalize();

    free(a);
    free(b);

    return 0;
}
