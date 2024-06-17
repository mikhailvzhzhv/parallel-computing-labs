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

        for (int i = 1; i < size; ++i) {
            MPI_Send(a + i * N / size, N / size, MPI_INT, i, 123, MPI_COMM_WORLD);
            MPI_Send(b, N, MPI_INT, i, 123, MPI_COMM_WORLD);
        }
        calculatePartSum(a, b, &partSum, size);
        sum += partSum;

        for (int i = 1; i < size; ++i) {
            MPI_Recv(&partSum, 1, MPI_INT, i, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            sum += partSum;
        }
        printf("sum = %d; time = %lf\n", sum, MPI_Wtime() - start);
    }
    else {
        MPI_Recv(a + rank * N / size, N / size, MPI_INT, 0, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(b, N, MPI_INT, 0, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        calculatePartSum(a + rank * N / size, b, &partSum, size);
        MPI_Send(&partSum, 1, MPI_INT, 0, 123, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    free(a);
    free(b);

    return 0;
}
