//
// Created by hongyao on 2018/11/1.
//
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>
#include "fio.h"

//Global variable
int64_t n,m,size,size_center;
int thread_num,threadid;
pthread_barrier_t barrier;
float *pa, *pb;

//********************************//
//**************子****************//
//**************线****************//
//**************程****************//
//********************************//
void  *worker(void *arg) {
    int64_t i,j,lb,ub;
    int myID = __sync_fetch_and_add(&threadid,1);
    int64_t loc_size = (n-2)/thread_num;
    int64_t rest = (n-2)%thread_num;
    //**********给线程分配计算资源*************
    if (myID < rest) {
        lb = loc_size * myID + myID;
        ub = lb + loc_size + 1;
    } else {
        lb = loc_size * myID + rest;
        ub = lb + loc_size;
    }
    if (pthread_barrier_wait(&barrier) == PTHREAD_BARRIER_SERIAL_THREAD) {
        memcpy(pb, pa, sizeof(float) * m);
        memcpy(&pb[m*(n-1)], &pa[m*(n-1)], sizeof(float) * m);
    }

    for (i=lb+1;i<ub+1;i++){
        pb[m*i]=pa[m*i];
        pb[m*(i+1)-1]=pa[m*(i+1)-1];
        for (j=i*m+1;j<i*m+m-1;j++){
            pb[j] = (pa[j] + pa[j - 1] + pa[j + 1] + pa[j-m] +pa[j+m]) / 5.0;
        }
    }
}

//********************************//
//**************主****************//
//**************线****************//
//**************程****************//
//********************************//
int main(int argc, char** argv ) {
    thread_num = atoi(argv[1]);
    n = atoll(argv[2]);
    m = atoll(argv[3]);
    size = ((int64_t) 1 << (n + m));
    size_center = size-2*m-2*n+4;
    threadid=0;
    n = 1L<<n;
    m = 1L<<m;
    if (posix_memalign((void **) &pa, getpagesize(), size * sizeof(float))) {
        perror("posix_memalign");
        return EXIT_SUCCESS;
    }
    if (posix_memalign((void **) &pb, getpagesize(), size * sizeof(float))) {
        perror("posix_memalign");
        return EXIT_SUCCESS;
    }
    pthread_barrier_init(&barrier,NULL,thread_num);    //设置栅樟
    input_data(pa, size * sizeof(float));
    pthread_t *threads = new pthread_t[thread_num];
    for (int i = 0; i < thread_num; i++) pthread_create(&(threads[i]), NULL, worker, &thread_num);
    for (int i = 0; i < thread_num; i++) pthread_join(threads[i], NULL);
    output_data(pb, size*sizeof(float));
    pthread_barrier_destroy(&barrier);      //删除栅樟*/
    free(pa);
    free(pb);
    return EXIT_SUCCESS;
}