#include <mpi.h>
#include <iostream>


void fill_array_zeros(char* arr, const int size);
void create_slider(char* arr, const int cols);
void check_neighbours(char* arr_origin, char* arr_copy, const int cols, const int rows);
void check_neighbours_line(char* arr_origin, char* arr_copy, char* arr_of_cells_pos, const int cols, const int offset, int mode);
void update_world(char* arr_origin, char* arr_copy, const int cols, const int rows);
void print_world(char* arr, const int cols, const int rows);
void print_board(const int cols);
void calculate_vector_stop(char* vector_stop, char* arr, char** arrays, const int len, const int iter);
char compare_arrays(char* arr1, char* arr2, const int len);
void check_vector_stop(char* vector, const int rows, const int cols, char* flag);
void swap_arrays(char* arr1, char* arr2, const int len);


int main(int argc, char* argv[]) {

    MPI_Init(&argc, &argv);

    if (argc < 4) {
        std::cerr << "Too few args\n";
        MPI_Finalize();
        return 1;
    }

    int rank, size;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double start;
    if (rank == 0) {
        start = MPI_Wtime();
    }

    const int cols = std::stoi(argv[1]);;
    const int rows = std::stoi(argv[2]);;
    const int max_gen = std::stoi(argv[3]);;
    int rows_per_rank = rows / size;
    if (rank < rows % size) {
        ++rows_per_rank;
    }
    const int len_arr = rows_per_rank * cols;
    const int len_vector = max_gen + size - max_gen % size;
    const int send_by_proc = len_vector / size;

    char* arr_of_cells = new char[len_arr];
    char* ptr_on_arr_of_cells = arr_of_cells;
    char** arr_worlds = new char*[max_gen];
    char* arr_of_cells_top = new char[cols];
    char* arr_of_cells_bot = new char[cols];
    char* vector_stop = new char[len_vector];
    char* vector_stop_buf = new char[len_vector];
    fill_array_zeros(arr_of_cells, len_arr);
    fill_array_zeros(vector_stop, len_vector);
    
    MPI_Request req_send_bot, req_send_top, req_recv_bot, req_recv_top;
    MPI_Request req_vector;
    MPI_Status status;

    if (rank == 0) {
        create_slider(arr_of_cells, cols);
    }

    const int mode_bot = -1;
    const int mode_top = 1;
    char flag = 0;
    char stop = 0;

    int gen;
    for (gen = 0; gen < max_gen && !stop; ++gen) {
        #ifdef TEST
        if (rank == 0) {
            print_board(cols);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        for (int r = 0; r < size; ++r) {
            if (rank == r) {
                print_world(arr_of_cells, cols, rows_per_rank);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) {
            print_board(cols);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        #endif

        MPI_Isend(arr_of_cells, cols, MPI_CHAR, (rank + size - 1) % size, 111, MPI_COMM_WORLD, &req_send_top);
        MPI_Isend(arr_of_cells + cols * (rows_per_rank - 1), cols, MPI_CHAR, (rank + 1) % size, 222, MPI_COMM_WORLD, &req_send_bot);

        MPI_Irecv(arr_of_cells_top, cols, MPI_CHAR, (rank + size - 1) % size, 222, MPI_COMM_WORLD, &req_recv_bot);
        MPI_Irecv(arr_of_cells_bot, cols, MPI_CHAR, (rank + 1) % size, 111, MPI_COMM_WORLD, &req_recv_top);

        calculate_vector_stop(vector_stop, arr_of_cells, arr_worlds, len_arr, gen);
        MPI_Ialltoall(vector_stop, send_by_proc, MPI_CHAR, vector_stop_buf, send_by_proc, MPI_CHAR, MPI_COMM_WORLD, &req_vector);
        
        char* arr_of_cells_copy = new char[len_arr];
        arr_worlds[gen] = arr_of_cells_copy;
        fill_array_zeros(arr_of_cells_copy, len_arr);
        check_neighbours(arr_of_cells, arr_of_cells_copy, cols, rows_per_rank);

        MPI_Wait(&req_send_top, &status);
        MPI_Wait(&req_recv_bot, &status);
        check_neighbours_line(arr_of_cells, arr_of_cells_copy, arr_of_cells_top, cols, 0, mode_top);

        MPI_Wait(&req_send_bot, &status);
        MPI_Wait(&req_recv_top, &status);
        check_neighbours_line(arr_of_cells, arr_of_cells_copy, arr_of_cells_bot, cols, cols * (rows_per_rank - 1), mode_bot);

        MPI_Wait(&req_vector, &status);
        check_vector_stop(vector_stop_buf, size, send_by_proc, &flag);

        MPI_Allreduce(&flag, &stop, 1, MPI_CHAR, MPI_LOR, MPI_COMM_WORLD);
        
        update_world(arr_of_cells, arr_of_cells_copy, cols, rows_per_rank);
        swap_arrays(arr_of_cells, arr_of_cells_copy, len_arr);
    }

    if (rank == 0) {
        std::cout << "time: " << MPI_Wtime() - start << "\n";
        std::cout << "iters: " << gen << "\n";
    }
    for (int i = 0; i < gen; ++i) {
        delete[] arr_worlds[i];
    }
    delete[] arr_worlds;
    delete[] ptr_on_arr_of_cells;
    delete[] vector_stop;
    delete[] vector_stop_buf;
    delete[] arr_of_cells_top;
    delete[] arr_of_cells_bot;

    MPI_Finalize();

    return 0;
}


void fill_array_zeros(char* arr, const int size) {
    for (int i = 0; i < size; ++i) {
        arr[i] = 0;
    }
}

void create_slider(char* arr, const int cols) {
    arr[1] = 1;
    arr[cols + 2] = 1;
    for (int i = 0; i < 3; ++i) { 
        arr[cols * 2 + i] = 1;
    }
}

void check_neighbours(char* arr_origin, char* arr_copy, const int cols, const int rows) {
    int idx;
    int idx_top;
    int idx_bot;
    for (int row = 1; row < rows - 1; ++row) {
        for (int col = 0; col < cols; ++col) {
            idx = cols * row + col;
            idx_top = cols * (row - 1) + col;
            idx_bot = cols * (row + 1) + col;
            if (arr_origin[(idx + 1) % cols + cols * row]) {
                arr_copy[idx] += 1;
            }
            if (arr_origin[(idx - 1) % cols + cols * row]) {
                arr_copy[idx] += 1;
            }
            if (arr_origin[idx_top]) {
                arr_copy[idx] += 1;
            }
            if (arr_origin[(idx_top + 1) % cols + cols * (row - 1)]) {
                arr_copy[idx] += 1;
            }
            if (arr_origin[(idx_top + cols - 1) % cols + cols * (row - 1)]) {
                arr_copy[idx] += 1;
            }
            if (arr_origin[idx_bot]) {
                arr_copy[idx] += 1;
            }
            if (arr_origin[(idx_bot + 1) % cols + cols * (row + 1)]) {
                arr_copy[idx] += 1;
            }
            if (arr_origin[(idx_bot - 1) % cols + cols * (row + 1)]) {
                arr_copy[idx] += 1;
            }
        }
    }
}

void check_neighbours_line(char* arr_origin, char* arr_copy, char* arr_of_cells_pos, const int cols, const int offset, int mode) {
    int idx;
    for(int col = 0; col < cols; ++col) {
        idx = offset + col;
        if (arr_origin[(idx + 1) % cols + offset]) {
            arr_copy[idx] += 1;
        }
        if (arr_origin[(idx + cols - 1) % cols + offset]) {
            arr_copy[idx] += 1;
        }
        if (arr_of_cells_pos[col]) {
            arr_copy[idx] += 1;
        }
        if (arr_of_cells_pos[(col + 1) % cols]) {
            arr_copy[idx] += 1;
        }
        if (arr_of_cells_pos[(col + cols - 1) % cols]) {
            arr_copy[idx] += 1;
        }
        if (arr_origin[idx + cols * mode]) {
            arr_copy[idx] += 1;
        }
        if (arr_origin[(col + cols + 1) % cols + cols * mode + offset]) {
            arr_copy[idx] += 1;
        }
        if (arr_origin[(col + cols - 1) % cols + cols * mode + offset]) {
            arr_copy[idx] += 1;
        }
    }
}

void update_world(char* arr_origin, char* arr_copy, const int cols, const int rows) {
    for (int i = 0; i < rows * cols; ++i) {
        if (arr_origin[i] && (arr_copy[i] == 2 || arr_copy[i] == 3)) {
            arr_copy[i] = 1;
        } else if (!arr_origin[i] && arr_copy[i] == 3) {
            arr_copy[i] = 1;
        } else {
            arr_copy[i] = 0;
        }
    }
}

void print_board(const int cols) {
    std::cout << " ";
    for (int i = 0; i < cols; ++i) {
        std::cout << "-";
    }
    std::cout << "\n";
}

void print_world(char* arr, const int cols, const int rows) {
    for (int row = 0; row < rows; ++row) {
        std::cout << "|";
        for (int col = 0; col < cols; ++col) {
            if (arr[row * cols + col]) {
                std::cout << "+";
            } else {
                std::cout << " ";
            }
        }
        std::cout << "|" << "\n";
    }
}

void calculate_vector_stop(char* vector_stop, char* arr, char** arrays, const int len, const int iter) {
    for (int i = 0; i < iter; ++i) {
        vector_stop[i] = compare_arrays(arr, arrays[i], len);
    }
}

char compare_arrays(char* arr1, char* arr2, const int len) {
    for (int i = 0; i < len; ++i) {
        if (arr1[i] != arr2[i]) {
            return 0;
        }
    }
    return 1;
}

void check_vector_stop(char* vector, const int rows, const int cols, char* flag) {
    int sum;
    for (int col = 0; col < cols; ++col) {
        sum = 0;
        for (int row = 0; row < rows; ++row) {
            if (vector[row * cols + col] != 1) {
               break;
            }
            ++sum;
        }
        if (sum == rows) {
            *flag = 1;
            return;
        }
    }
}

void swap_arrays(char* arr1, char* arr2, const int len) {
    char sw;
    for (int i = 0; i < len; ++i) {
        sw = arr1[i];
        arr1[i] = arr2[i];
        arr2[i] = sw;
    }
}
