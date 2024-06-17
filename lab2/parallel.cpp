#include <mpi.h>
#include <iostream>
#include <string>

void fillAMatrix(double* matrix, int rows, int columns);
void fillBMatrix(double* matrix, int rows, int columns);
void mulMatrix(double* A, double* B, double* C, int N1, int N2, int N3);

int main(int argc, char *argv[]) {

    if (argc < 3) {
        std::cerr << "too few args\n";
        return -1;
    }

    MPI_Init(&argc, &argv);

    double start = MPI_Wtime();
    
    int dims[2];
    dims[0] = std::stoi(argv[1]);
    dims[1] = std::stoi(argv[2]);

    int periods[2] = {0,0};
    int coords[2];
    int reorder = 0;
    int size, rank, sizey, sizex, ranky, rankx;

    MPI_Comm comm_2d;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    sizey = dims[0];
    sizex = dims[1];

    const int N1 = 10 * sizey;
    const int N2 = 14;
    const int N3 = 10 * sizex;

    double* A;
    double* B;
    double* C;

    double* partA = new double[N1 / sizey * N2];
    double* partB = new double[N3 / sizex * N2];
    double* partC = new double[N1 / sizey * N3 / sizex];
    
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &comm_2d);
    MPI_Cart_get(comm_2d, 2, dims, periods, coords);
    MPI_Cart_rank(comm_2d, coords, &rank);
    ranky = coords[0];
    rankx = coords[1];

    if (rankx == 0 && ranky == 0) {
        A = new double[N1 * N2];
        B = new double[N2 * N3];
        C = new double[N1 * N3];

        fillAMatrix(A, N1, N2);
        fillBMatrix(B, N2, N3);
    }

    MPI_Datatype col_vector, col_type, row_type, col_vector_resized;
    MPI_Type_contiguous(N3 / sizex * N2, MPI_DOUBLE, &col_type);
    MPI_Type_commit(&col_type);
    MPI_Type_contiguous(N1 / sizey * N2, MPI_DOUBLE, &row_type);
    MPI_Type_commit(&row_type);
    MPI_Type_vector(N2, N3 / sizex, N3, MPI_DOUBLE, &col_vector);
    MPI_Type_commit(&col_vector);
    MPI_Type_create_resized(col_vector, 0, N3 / sizex * sizeof(double), &col_vector_resized);
    MPI_Type_commit(&col_vector_resized);
    MPI_Comm rowComm, colComm;
    int subdims[2];

    subdims[0] = 0;
    subdims[1] = 1;
    MPI_Cart_sub(comm_2d, subdims, &rowComm);

    subdims[0] = 1;
    subdims[1] = 0;
    MPI_Cart_sub(comm_2d, subdims, &colComm);

    if (rankx == 0) {
        MPI_Scatter(A, 1, row_type, partA, 1, row_type, 0, colComm);
    }
    if (ranky == 0) {
        MPI_Scatter(B, 1, col_vector_resized, partB, 1, col_type, 0, rowComm);
    }

    MPI_Bcast(partA, 1, row_type, 0, rowComm);
    MPI_Bcast(partB, 1, col_type, 0, colComm);

    mulMatrix(partA, partB, partC, N1 / sizey, N2, N3 / sizex);

    MPI_Datatype c_vector, resized_vector_c;
    MPI_Type_vector(N1 / sizey, N3 / sizex, N3, MPI_DOUBLE, &c_vector);
    MPI_Type_commit(&c_vector);
    MPI_Type_create_resized(c_vector, 0, N3 / sizex * sizeof(double), &resized_vector_c);
    MPI_Type_commit(&resized_vector_c);

    int recvcounts[sizex * sizey];
    int displs[sizex * sizey];
    for (int i = 0; i < sizey; ++i) {
        for (int j = 0; j < sizex; ++j) {
            displs[i * sizex + j] = i * N1 / sizey * sizex + j;
            recvcounts[i * sizex + j] = 1;
        }
    }

    MPI_Gatherv(partC, N1 / sizey * N3 / sizex, MPI_DOUBLE, C, recvcounts, displs, resized_vector_c, 0, comm_2d);

    if (rank == 0) {
	std::cout << "time: "  << MPI_Wtime() - start << "\n";
    }

    if (rankx == 0 && ranky == 0) {
        for (int i = 0; i < N1 * N3; ++i) {
            std::cout << C[i] << " ";
        }
    }

    if (rankx == 0 && ranky == 0) {
        delete[] A;
        delete[] B;
        delete[] C;
    }
    delete[] partA;
    delete[] partB;
    delete[] partC;

    MPI_Type_free(&col_vector);
    MPI_Type_free(&col_type);
    MPI_Type_free(&row_type);
    MPI_Type_free(&col_vector_resized);
    MPI_Type_free(&c_vector);
    MPI_Type_free(&resized_vector_c);

    MPI_Comm_free(&comm_2d);
    MPI_Comm_free(&rowComm);
    MPI_Comm_free(&colComm);

    MPI_Finalize();

    return 0;
}

void fillAMatrix(double* matrix, int rows, int columns) {
    for (int i = 0; i < rows * columns; i++) {
        matrix[i] = 1;
    }
}

void fillBMatrix(double* matrix, int rows, int columns) {
    for (int i = 0; i < rows * columns; i++) {
        matrix[i] = 1;
    }
}

void fillMatrixZeros(double* matrix, int rows, int columns) {
    for (int i = 0; i < rows * columns; i++) {
        matrix[i] = 0;
    }
}

void mulMatrix(double* A, double* B, double* C, int N1, int N2, int N3) {
    fillMatrixZeros(C, N1, N3);
    for (int i = 0; i < N1; ++i) {
        for (int j = 0; j < N2; ++j) {
            for (int k = 0; k < N3; ++k) {
                C[i * N3 + k] += A[i * N2 + j] * B[j * N3 + k];
            }
        }
    }
}
