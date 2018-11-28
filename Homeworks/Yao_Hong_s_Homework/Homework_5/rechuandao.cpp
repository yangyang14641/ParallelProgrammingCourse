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
int32_t *pinit;
int64_t n,m,size,k;
int thread_num,threadid;
pthread_barrier_t barrier;
float *pa, *pb ,*lowT;
float eps;
bool FLAG=true;
float loweresT = 0;
//********************************//
//**************子****************//
//**************线****************//
//**************程****************//
//********************************//
void  *worker(void *arg) {
    int64_t lb, ub, lc, uc, ld, ud, le, ue;
    int64_t pre, cur, nxt;
    int64_t i, j;
   // bool flag = true;
    float loweresT_p;
    int myID = __sync_fetch_and_add(&threadid, 1);
    int64_t loc_size = (n - 2) / thread_num;
    int64_t loc_size_lie = (m - 2) / thread_num;
    int64_t rest = (n - 2) % thread_num;
    int64_t rest_lie = (m - 2) % thread_num;
    int64_t loc_size_cpn = n / thread_num;
    int64_t loc_size_cpn_lie = m / thread_num;
    int64_t rest_cpn = (n) % thread_num;
    int64_t rest_cpn_lie = m % thread_num;
    //**********给线程分配计算资源*************
    if (myID < rest) {
        lb = loc_size * myID + myID;
        ub = lb + loc_size + 1;
    }
    else {
        lb = loc_size * myID + rest;
        ub = lb + loc_size;
    }

    if (myID < rest_lie) {
        lc = loc_size_lie * myID + myID;
        uc = lc + loc_size_lie + 1;
    }
    else {
    lc = loc_size_lie * myID + rest_lie;
    uc = lc + loc_size_lie;
    }

    if (myID < rest_cpn) {
        ld = loc_size_cpn * myID + myID;
        ud = ld + loc_size_cpn + 1;
    }
    else {
        ld = loc_size_cpn * myID + rest_cpn;
        ud = ld + loc_size_cpn;
    }

    if (myID < rest_cpn_lie) {
        le = loc_size_cpn_lie * myID + myID;
        ue = le + loc_size_cpn_lie + 1;
    }
    else {
        le = loc_size_cpn_lie * myID + rest_cpn_lie;
        ue = le + loc_size_cpn_lie;
    }
    //**********开始迭代*************
    while (FLAG) {
        nxt = m;
        for (j = lc + 1; j < uc + 1; j++) {
            pb[j] = 0.325 * pa[j] + 0.225 * (pa[nxt + j] + pa[j + 1] + pa[j - 1]);//第一行
        }
        cur = (lb + 1) * m;
        pre = lb * m;
        nxt = cur + m;
        for (i = lb + 1; i < ub + 1; i++) {
            pb[cur] = 0.325 * pa[cur] + 0.225 * (pa[pre] + pa[cur + 1] + pa[nxt]);
            for (j = 1; j < m - 1; j++)
                pb[cur + j] =
                        0.1 * pa[cur + j] + 0.225 * (pa[pre + j] + pa[nxt + j] + pa[cur + j + 1] + pa[cur + j - 1]);
            pb[cur + j] = 0.325 * pa[cur + j] + 0.225 * (pa[pre + j] + pa[cur + j - 1] + pa[nxt + j]);
            pre = cur;
            cur = nxt;
            nxt += m;
        }
        cur = (n-1)*m;
        pre = cur-m;
        for (j = lc + 1; j < uc + 1; j++) {
            pb[cur + j] = 0.325 * pa[cur + j] + 0.225 * (pa[pre + j] + pa[cur + j + 1] + pa[cur + j - 1]);//最后一行
        }
        pthread_barrier_wait(&barrier);
        if (myID==0) {
            pb[0] = 0.1 * pa[0] + 0.45 * (pa[1] + pa[m]);//四个角点
            pb[m - 1] = 0.1 * pa[m - 1] + 0.45 * (pa[m - 2] + pa[m + m - 1]);
            cur = (n - 1) * m;
            pb[cur] = 0.1 * pa[cur] + 0.45 * (pa[cur + 1] + pa[cur - m]);
            pb[cur + m - 1] = 0.1 * pa[cur + m - 1] + 0.45 * (pa[cur + m - 2] + pa[cur - 1]);
            for (i = 0; i < k; i++)     {
                pb[pinit[i * 3] * m + pinit[i * 3 + 1]] = pa[pinit[i * 3] * m + pinit[i * 3 + 1]];
            }                        //定热源
            FLAG = false;
        }
        pthread_barrier_wait(&barrier);
        //**********判断是否继续迭代*************
        cur = ld * m;
        for (i = ld; i < ud; i++) {
            for (j = 0; j < m; j++) {
                if (fabs(pb[cur + j] - pa[cur + j]) >= eps) {
                    FLAG = true;
                    break;
                }
            }
            if (FLAG) {
                break;
            }
            cur += m;
        }
     	pthread_barrier_wait(&barrier);
        if (myID==0) {
     		float *temp = pa;
     			pa = pb;
     			pb = temp;
        }
        pthread_barrier_wait(&barrier);
    }
    //**********开始寻找局部最低温*************
    cur = ld * m;
    loweresT_p = loweresT;
    for (i = ld; i < ud; i++) {
        for (j = 0; j < m; j++)
            if (loweresT_p > pa[cur + j]) {
                loweresT_p = pa[cur + j];
            }
        cur += m;
    }
    lowT[myID] = loweresT_p;
    return (void *) 0;
    }

//********************************//
//**************主****************//
//**************线****************//
//**************程****************//
//********************************//
int main(int argc, char** argv ) {
    int64_t i, j;
    thread_num = atoi(argv[1]); //p
    n = atoll(argv[2]);     //n
    m = atoll(argv[3]);     //m
    k = atoi(argv[4]);      //k
    eps = atof(argv[5]);      //eps

    size = ((int64_t) 1 << (n + m));
    pa = (float *) malloc(size * sizeof(float));
    pb = (float *) malloc(size * sizeof(float));
    lowT = (float *) malloc(thread_num * sizeof(float));
    n = 1L << n;
    m = 1L << m;
    pinit = (int32_t *) malloc(3 * k * sizeof(int32_t));  //热源位置保存
    pthread_barrier_init(&barrier, NULL, thread_num);    //初始化栅樟
    threadid = 0;

    input_data(pinit, 3 * k * sizeof(int32_t));
    for (int i = 0; i < k; i++) {
        float *fp = (float *) &pinit[i * 3 + 2];
        pa[pinit[i * 3] * m + pinit[i * 3 + 1]] = *fp;
        if (loweresT < *fp) loweresT = *fp;
    }
    pthread_t *threads = new pthread_t[thread_num];
    for (int i = 0; i < thread_num; i++) pthread_create(&(threads[i]), NULL, worker, &thread_num);
    for (int i = 0; i < thread_num; i++) pthread_join(threads[i], NULL);
    for (j = 0; j < thread_num; j++)
    for (j = 0; j < thread_num; j++)
        if (loweresT > lowT[j]) {
            loweresT = lowT[j];
        }                         //找极值*/
    output_data(&loweresT, sizeof(float));
    pthread_barrier_destroy(&barrier);      //删除栅樟*/
    free(pa);
    free(pb);
    free(pinit);
    return EXIT_SUCCESS;
}