#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <vector>
#include <math.h>

#include "fio.h"

int64_t size,k;
int thread_num, threadid;
pthread_barrier_t barrier;
double *pbd, *pfc ;
int *lb,*ub;
int32_t T ;

void  *worker(void *arg) {
    int64_t i, j, w ;
    int32_t t=0 ,cir;
    int myID = __sync_fetch_and_add(&threadid, 1);
    int64_t loc_size = (size) / thread_num;
    int64_t rest = (size) % thread_num;

    if (myID < rest) {
        lb[myID] = loc_size * myID + myID;
        ub[myID] = lb[myID] + loc_size + 1;
    } else {
        lb[myID] = loc_size * myID + rest;
        ub[myID] = lb[myID] + loc_size;
    }

    int num = (thread_num + 1) / 2;
    if (myID < num && thread_num % 2 == 0) {
        cir = num + 1;//w  Cycles
    } else {
        cir = num;
    }
    printf("cir[%d]=%ld",myID,cir);
    pthread_barrier_wait(&barrier);

    while (t < T) {
        for (w = myID; w < myID + cir; w++) {
            if (w < thread_num) {
                for (i = lb[myID]; i < ub[myID]; i++) {
                    int64_t bi = (i << 2);
                    int64_t fi = bi - i;
                    for (j = lb[w]; j < ub[w]; j++) {
                        int64_t bj = (j << 2);
                        int64_t fj = bj - j;

                        double dx = pbd[bi + 1] - pbd[bj + 1];
                        double dy = pbd[bi + 2] - pbd[bj + 2];
                        double dz = pbd[bi + 3] - pbd[bj + 3];
                        double sq = dx * dx + dy * dy + dz * dz;
                        double dist = sqrt(sq);
                        double fac = G * pbd[bi] * pbd[bj] / (dist * sq);
                        double fx = fac * dx;
                        double fy = fac * dy;
                        double fz = fac * dz;

                        pfc[fi] -= fx;
                        pfc[fi + 1] -= fy;
                        pfc[fi + 2] -= fz;
                        pfc[fj] += fx;
                        pfc[fj + 1] += fy;
                        pfc[fj + 2] += fz;
                    }
                }
            }
            if (w >= thread_num) {
                for (i = lb[w - thread_num]; i < ub[w - thread_num]; i++) {
                    int64_t bi = (i << 2);
                    int64_t fi = bi - i;
                    for (j = lb[myID]; j < ub[myID]; j++) {
                        int64_t bj = (j << 2);
                        int64_t fj = bj - j;

                        double dx = pbd[bi + 1] - pbd[bj + 1];
                        double dy = pbd[bi + 2] - pbd[bj + 2];
                        double dz = pbd[bi + 3] - pbd[bj + 3];
                        double sq = dx * dx + dy * dy + dz * dz;
                        double dist = sqrt(sq);
                        double fac = G * pbd[bi] * pbd[bj] / (dist * sq);
                        double fx = fac * dx;
                        double fy = fac * dy;
                        double fz = fac * dz;

                        pfc[fi] -= fx;
                        pfc[fi + 1] -= fy;
                        pfc[fi + 2] -= fz;
                        pfc[fj] += fx;
                        pfc[fj + 1] += fy;
                        pfc[fj + 2] += fz;
                    }
                }
            }
            pthread_barrier_wait(&barrier);
        }
        pthread_barrier_wait(&barrier);


        /*for (i = lb[myID]; i < ub[myID]; i++) {
            int64_t bi = (i << 2);
            int64_t fi = bi - i;
            pbd[bi + 1] = pbd[bi + 1] + pfc[fi] / pbd[bi];

            pfc[fi] = 0;
            pbd[bi + 2] = pbd[bi + 2] + pfc[fi + 1] / pbd[bi];
            pfc[fi + 1] = 0;
            pbd[bi + 3] = pbd[bi + 3] + pfc[fi + 2] / pbd[bi];
            pfc[fi + 2] = 0;
        }*/
        pthread_barrier_wait(&barrier);
        t++;
        printf("t=%d",t);
    }
    return (void *) 0;
}

int main(int argc, char** argv ) {
    int64_t i, j;

    thread_num = atoi(argv[1]);
    size = atoll(argv[2]);
    T = atoi(argv[3]);
    printf("T=%d",T);
    size = ((int64_t)1<<size);
    pbd = (double*)malloc(4*size*sizeof(double));

    lb = (int*)malloc(thread_num*sizeof(int));
    ub = (int*)malloc(thread_num*sizeof(int));

    pthread_barrier_init(&barrier, NULL, thread_num);
    threadid = 0;

    input_data(pbd, 4*size*sizeof(double));
    pfc = (double*)malloc(3*size*sizeof(double));
    memset(pfc, 0, 3*size*sizeof(double));

    pthread_t *threads = new pthread_t[thread_num];
    for (int i = 0; i < thread_num; i++) pthread_create(&(threads[i]), NULL, worker, &thread_num);
    for (int i = 0; i < thread_num; i++) pthread_join(threads[i], NULL);

    output_data(pbd, 4*size*sizeof(double));
    free(pbd);
    free(pfc);
    pthread_barrier_destroy(&barrier);
    return EXIT_SUCCESS;
}
