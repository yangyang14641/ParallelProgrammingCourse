#include <pthread.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <vector>
#include <math.h>
#include "fio.h"

int64_t n,n_num,m,m_num,sizeofA,sizeofB,sizeofC;
int32_t thread_num,threadid;
int64_t *A,*B,*C,*pp,*qq,*qqa;
pthread_barrier_t barrier;
struct SGROUP {int64_t pianduan;int64_t geshu;int64_t locA;}; //pianduan为该片段，geshu为该片段的个数，locA为排序后该片段的首位在A中的位置
struct SGROUP *groupA;
//比较函数指针
int myCompar(const void *arg1,const void *arg2){
    int64_t *pa=(int64_t*)arg1,*pb=(int64_t*)arg2;
    return *pa>*pb;
}
//***********************//
//**********子***********//
//**********线***********//
//**********程***********//
//***********************//
void  *worker(void *arg) {
    int64_t i,j,lb,ub;
    int myID = __sync_fetch_and_add(&threadid, 1);
//    printf("xianchengyfhdshd s");
    int64_t loc_size = (m_num /2)/thread_num;
    int64_t rest = (m_num /2)%thread_num;

    //**********给线程分配计算资源*************
    if (myID < rest) {
        lb = loc_size * myID + myID;
        ub = lb + loc_size + 1;
    } else {
        lb = loc_size * myID + rest;
        ub = lb + loc_size;
    }

    //**********将B数组两位两位地保存在pp中***********
    for (i =lb; i < ub; i++) {
        for (j = 0; j < m_num /2; j++) {
            if (B[2*j] == i) {
                pp[2*i]=B[2*j];
                pp[2*i+1]=B[2*j+1];
            }
        }
    }
    pthread_barrier_wait(&barrier);
//**************开辟空间得到qqa（排序前累积）以及qq（排序后累积）***************
    if (pthread_barrier_wait(&barrier) == PTHREAD_BARRIER_SERIAL_THREAD) {
        for (i = 1; i < m_num / 2; i++) {
            qqa[i] = qqa[i - 1] + B[2*i - 1];
            qq[i] = qq[i - 1] + pp[2*i - 1];
        }
    }
    pthread_barrier_wait(&barrier);
//***************保存结构体的前两个数据**************
    for (i = lb; i <ub; i++) {
        groupA[i].pianduan = B[2*i];
        groupA[i].geshu = B[2*i+1];
    }
    pthread_barrier_wait(&barrier);
    //**************保存结构体的第三个数据**************
    for (i = lb; i <ub; i++) {
        groupA[i].locA = qq[groupA[i].pianduan];
    }
    pthread_barrier_wait(&barrier);
    //***************先对每个片段排序，再拷贝排序前A的片段到C的相应的位置****************
     for (i = lb; i < ub; i++) {
         qsort(&A[qqa[i]], groupA[i].geshu, sizeof(int64_t), myCompar);
         memcpy(&C[groupA[i].locA], &A[qqa[i]], sizeof(int64_t) * groupA[i].geshu);
     }

    return (void *) 0;
}

int main(int argc, char *argv[ ] ){
    int i;
    thread_num=atoi(argv[1]);        //      设置P
    n=atoll(argv[2]);                        //      设置N
    n_num=((int64_t) 1<<n);           //       设置2^N大小的数组
    m=atoll(argv[3]);                        //      设置M
    m_num=((int64_t) 1<<(m+1));           //      设置2*2^M大小的数组
    threadid=0;
    pp= (int64_t *)malloc(sizeof(int64_t)*n_num);   //将B重新排序后的数组
    qq= (int64_t *)malloc(sizeof(int64_t)*m_num/2);   //A排序后B中个数的累加数组
    qqa= (int64_t *)malloc(sizeof(int64_t)*m_num/2);   //A排序前B中个数的累加数组
    qq[0]=0;
    qqa[0]=0;
    sizeofA = (sizeof(int64_t)) *n_num;              //A数组空间大小
    sizeofB = (sizeof(int64_t)) *m_num;              //B数组空间大小
    sizeofC = (sizeof(int64_t)) *n_num;              //C数组空间大小
    A= (int64_t *)malloc(sizeof(int64_t)*n_num);     //开辟A数组空间
    B= (int64_t *)malloc(sizeof(int64_t)*m_num);     //开辟B数组空间
    C= (int64_t *)malloc(sizeof(int64_t)*n_num);     //开辟输出数组C数的组空间
    groupA=(struct SGROUP*)malloc((m_num/2)*sizeof(SGROUP));
    input_data(A, sizeofA, B, sizeofB);        //导入A、B。
    pthread_barrier_init(&barrier,NULL,thread_num);    //设置栅樟
    pthread_t *threads = new pthread_t[thread_num];
    for(int i=0; i<thread_num; i++)  pthread_create(&(threads[i]),NULL,worker,&thread_num);
    for(int i=0; i<thread_num; i++)  pthread_join(threads[i],NULL);
    output_data(C,sizeofC ); //输出C
    delete[ ] threads;
    pthread_barrier_destroy(&barrier);      //删除栅樟*/
    return EXIT_SUCCESS;
}








//
// Created by 洪垚 on 2018/10/29.
//

