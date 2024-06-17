#include <mpi.h>
#include <pthread.h>
#include <math.h>
#include <stdio.h>

int iterCountMax = 3;
int L = 3000;
int *tasks;
int taskNum;

double globalRes = 0;
int tasksCount = 48000;
int tasksPerProcess;

int size, rank;
int root = 0;

pthread_t thread_exec;
pthread_t thread_support;
pthread_mutex_t mutex;
pthread_attr_t attr;

void fillTasks(int *tasks, int count, int iterCount);
void initThreads();
void *executeTasks(void *);
void *sendTask(void *);

int main(int argc, char* argv[]) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    pthread_mutex_init(&mutex, nullptr);
    
    tasksPerProcess = tasksCount / size;
    tasks = new int[tasksPerProcess];

    double startTime = MPI_Wtime();
    initThreads();
    double allTime = MPI_Wtime() - startTime;
    pthread_mutex_destroy(&mutex);
    pthread_attr_destroy(&attr);
    if(rank == root) {
        printf("Time:%f\n", allTime);
    }

    delete[] tasks;

    MPI_Finalize();

    return 0;
}

void initThreads() {

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    pthread_create(&thread_exec, &attr, executeTasks, nullptr);
    pthread_create(&thread_support, &attr, sendTask, nullptr);

    pthread_join(thread_exec, nullptr);
    pthread_join(thread_support, nullptr);
}

void fillTasks(int *tasks, int count, int iterCount) {
    for(int i = 0; i < count; i++) {
        tasks[i] = abs(50 - i % 100) * abs(rank - (iterCount % size)) * L;
    }
}

void *executeTasks(void *args) {
    int taskDone;
    double start, time;
    double m, n;
    int weight;
    int req, res;
    double GlobalCommonRes;
    
    for (int iterCount = 0; iterCount < iterCountMax; ++iterCount) {
        fillTasks(tasks, tasksPerProcess, iterCount);
        taskDone = 0;
        taskNum = 0;
        start = MPI_Wtime();

        while(taskNum < tasksPerProcess) {
            pthread_mutex_lock(&mutex);
            weight = tasks[taskNum];
            taskNum++;
            pthread_mutex_unlock(&mutex);

            for(int i = 0; i < weight; i++) {
                globalRes += sqrt(i);
            }
            taskDone++;
        }

        for(int r = 1; r < size; r++){
            req = 1;

            while(true) {
                MPI_Send(&req, 1, MPI_INT, (rank + r) % size, 0, MPI_COMM_WORLD);
                MPI_Recv(&res, 1, MPI_INT, (rank + r) % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if(res == -1) {
                    break;
                }
                int* tasksToRecv = new int[res];
                MPI_Recv(tasksToRecv, res, MPI_INT, (rank + r) % size, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for(int i = 0; i < res; i++) {
                    weight = tasksToRecv[i];
                    for(int k = 0; k < weight; k++) {
                        globalRes += sqrt(k);
                    }
                    taskDone++;
                }
                delete tasksToRecv;
            }
        }

        time = MPI_Wtime() - start;
        printf("rank = %d | globalRes = %.3lf | taskDone = %d | time = %lf\n", rank, globalRes, taskDone, time);

        MPI_Reduce(&time, &m, 1, MPI_DOUBLE, MPI_MAX, root, MPI_COMM_WORLD);
        MPI_Reduce(&time, &n, 1, MPI_DOUBLE, MPI_MIN, root, MPI_COMM_WORLD);

        if(rank == root) {
            double delta = m - n;
            printf("disbalance = %.3lf\n", delta);
            printf("part of disbalance = %.2lf\n", delta / m * 100);
        }

        MPI_Allreduce(&globalRes, &GlobalCommonRes, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if(rank == root) {
            printf("GlobalCommonRes  = %.3f\n", GlobalCommonRes);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    req = 0;
    MPI_Send(&req, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);

    pthread_exit(nullptr);
}

void *sendTask(void *args) {
    MPI_Status status;
    int req;
    int res;

    while(true) {
        MPI_Recv(&req, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        if(req == 0) {
            break;
        }

        pthread_mutex_lock(&mutex);
        int *taskToSend;
        if(taskNum < tasksPerProcess / 100 * 99) {
            res = tasksPerProcess / 200;
            taskToSend = new int[res];
            for(int i = 0; i < res; i++){
                taskToSend[i] = tasks[taskNum];
                taskNum++;
            }
        } else {
            res = -1;
        }
        pthread_mutex_unlock(&mutex);

        MPI_Send(&res, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
        if(res != -1) {
            MPI_Send(taskToSend, res, MPI_INT, status.MPI_SOURCE, 2, MPI_COMM_WORLD);
            delete[] taskToSend;
        }
    }

    pthread_exit(nullptr);
}
