#include <mpi.h>
#include <iostream>
#include <string>

void fillAMatrix(double* matrix, const int rows, const int columns);
void fillBMatrix(double* matrix, const int rows, const int columns);
void mulMatrix(double* A, double* B, double* C, const int rowsA, const int common_side, const int colB);

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
    int size, rank, size_y, size_x, rank_y, rank_x;
    const int dims_size = 2;
    const int root = 0;

    MPI_Comm comm_2d;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    size_y = dims[0];
    size_x = dims[1];

    const int rows_A = 24 * 170;
    const int common_side = 1200;
    const int cols_B = 24 * 130;

    const int rows_partA = rows_A / size_y;
    const int cols_partB = cols_B / size_x;

    double* A;
    double* B;
    double* C;

    double* partA = new double[rows_partA * common_side];
    double* partB = new double[cols_partB * common_side];
    double* partC = new double[rows_partA * cols_partB];
    
    MPI_Cart_create(MPI_COMM_WORLD, dims_size, dims, periods, reorder, &comm_2d);
    MPI_Cart_get(comm_2d, dims_size, dims, periods, coords);
    MPI_Cart_rank(comm_2d, coords, &rank);
    rank_y = coords[0];
    rank_x = coords[1];

    if (rank_x == 0 && rank_y == 0) {
        A = new double[rows_A * common_side];
        B = new double[common_side * cols_B];
        C = new double[rows_A * cols_B];

        fillAMatrix(A, rows_A, common_side);
        fillBMatrix(B, common_side, cols_B);
    }

    MPI_Datatype vector_col, type_col, type_row, vector_col_resized;
    MPI_Type_contiguous(cols_partB * common_side, MPI_DOUBLE, &type_col);
    MPI_Type_commit(&type_col);
    MPI_Type_contiguous(rows_partA * common_side, MPI_DOUBLE, &type_row);
    MPI_Type_commit(&type_row);
    MPI_Type_vector(common_side, cols_partB, cols_B, MPI_DOUBLE, &vector_col);
    MPI_Type_commit(&vector_col);
    MPI_Type_create_resized(vector_col, 0, cols_partB * sizeof(double), &vector_col_resized);
    MPI_Type_commit(&vector_col_resized);
    MPI_Comm comm_row, comm_col;

    int subdims[2];

    subdims[0] = 0;
    subdims[1] = 1;
    MPI_Cart_sub(comm_2d, subdims, &comm_row);

    subdims[0] = 1;
    subdims[1] = 0;
    MPI_Cart_sub(comm_2d, subdims, &comm_col);

    if (rank_x == 0) {
        MPI_Scatter(A, 1, type_row, partA, 1, type_row, root, comm_col);
    }
    if (rank_y == 0) {
        MPI_Scatter(B, 1, vector_col_resized, partB, 1, type_col, root, comm_row);
    }

    MPI_Bcast(partA, 1, type_row, root, comm_row);
    MPI_Bcast(partB, 1, type_col, root, comm_col);

    mulMatrix(partA, partB, partC, rows_partA, common_side, cols_partB);

    MPI_Datatype vector_c, vector_c_resized;
    MPI_Type_vector(rows_partA, cols_partB, cols_B, MPI_DOUBLE, &vector_c);
    MPI_Type_commit(&vector_c);
    MPI_Type_create_resized(vector_c, 0, cols_partB * sizeof(double), &vector_c_resized);
    MPI_Type_commit(&vector_c_resized);

    int recvcounts[size_x * size_y];
    int displs[size_x * size_y];
    for (int i = 0; i < size_y; ++i) {
        for (int j = 0; j < size_x; ++j) {
            displs[i * size_x + j] = i * rows_partA * size_x + j;
            recvcounts[i * size_x + j] = 1;
        }
    }

    MPI_Gatherv(partC, rows_partA * cols_partB, MPI_DOUBLE, 
                C, recvcounts, displs, vector_c_resized, root, comm_2d);

    if (rank_x == 0 && rank_y == 0) {
        std::cout << "time: " << MPI_Wtime() - start << "\n";
    }

    if (rank_x == 0 && rank_y == 0) {
        for (int i = 0; i < rows_A * cols_B; ++i) {
            std::cout << C[i] << " ";
        }
    }

    if (rank_x == 0 && rank_y == 0) {
        delete[] A;
        delete[] B;
        delete[] C;
    }
    delete[] partA;
    delete[] partB;
    delete[] partC;

    MPI_Type_free(&vector_col);
    MPI_Type_free(&type_col);
    MPI_Type_free(&type_row);
    MPI_Type_free(&vector_col_resized);
    MPI_Type_free(&vector_c);
    MPI_Type_free(&vector_c_resized);

    MPI_Comm_free(&comm_2d);
    MPI_Comm_free(&comm_row);
    MPI_Comm_free(&comm_col);

    MPI_Finalize();

    return 0;
}

void fillAMatrix(double* matrix, const int rows, const int columns) {
    for (int i = 0; i < rows * columns; i++) {
        matrix[i] = 1;
    }
}

void fillBMatrix(double* matrix, const int rows, const int columns) {
    for (int i = 0; i < rows * columns; i++) {
        matrix[i] = 1;
    }
}

void fillMatrixZeros(double* matrix, const int rows, const int columns) {
    for (int i = 0; i < rows * columns; i++) {
        matrix[i] = 0;
    }
}

void mulMatrix(double* A, double* B, double* C, const int rowsA, const int common_side, const int colB) {
    fillMatrixZeros(C, rowsA, colB);
    for (int i = 0; i < rowsA; ++i) {
        for (int j = 0; j < common_side; ++j) {
            for (int k = 0; k < colB; ++k) {
                C[i * colB + k] += A[i * common_side + j] * B[j * colB + k];
            }
        }
    }
}