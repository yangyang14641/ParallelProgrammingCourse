//
// Created by hongyao on 2018/10/15.
//

#include <pthread.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <vector>
#include <math.h>
#include "fio.h"

int64_t n,n_num,sizeofA,sizeofB,single_num;
int32_t thread_num,threadid;
float *A,*B,alfa;
pthread_barrier_t barrier;

void  *worker(void *arg){
    int64_t i;
    int myID = __sync_fetch_and_add(&threadid, 1);
    for (i=myID*single_num;i<(myID+1)*single_num;i++) *(B+i)=(*(A+i))*alfa+*(B+i);
    return (void*)0;
}

int main(int argc, char *argv[]){
    int i;
    thread_num=atoi(argv[1]);        //      设置P
    n=atoll(argv[2]);                        //      设置N
    n_num=((int64_t)1<<n);           //      设置2^N大小的数组
    single_num=n_num/thread_num;      //单个线程的分配数量
    alfa=atof(argv[3]);
    threadid=0;
    sizeofA = (sizeof(float)) *n_num;       //A数组空间大小
    sizeofB = (sizeof(float)) *n_num;       //B数组空间大小
    A= (float *)malloc(sizeof(float)*n_num);   //开辟A数组空间
    B= (float *)malloc(sizeof(float)*n_num);     //开辟B数组空间
    input_data(A, sizeofA, B, sizeofB);        //输入AB。
    printf("hello: thread is running!\n");
    // printf("A[2]=%f\n;B[2]=%f\nalfa=%f", A[2], B[2],alfa);  //检查
    pthread_barrier_init(&barrier,NULL,thread_num);    //设置栅樟
    pthread_t *threads = new pthread_t[thread_num];
    for(int i=0; i<thread_num; i++)  pthread_create(&(threads[i]),NULL,worker,&thread_num);
    for(int i=0; i<thread_num; i++)  pthread_join(threads[i],NULL);
    output_data(B,sizeofB ); //输出B
    //printf("B[2]=%f\n", B[2]);
    //printf("B[1]=%f\n", *(B+1));
    delete[ ] threads;
    pthread_barrier_destroy(&barrier);      //删除栅樟
    return EXIT_SUCCESS;
}
