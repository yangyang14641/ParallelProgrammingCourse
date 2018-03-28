#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <pthread.h>               

#define NANO           10000000000           
#define Max_Thread_Num 8                           

static long num_steps = 1000000;           
double step, sum = 0.0;                       
double serial_PI() {                                                           // serial computing
       int i;                                 
       double x;                                
       
       sum = 0.0;
       step = 1.0/(double) num_steps;
       for (i=1;i<= num_steps; i++) {
         x = (i-0.5)*step;
         sum = sum + 4.0/(1.0+x*x);
        }
        
        return step * sum;
}

int totalThread;
pthread_t        threads[Max_Thread_Num];
pthread_cond_t   cond;
pthread_mutex_t  mtx;

void *syn_worker(void *arg) {                                                  // Parallel computing using
     int    myID, lbound, ubound, i, loc_size; 
     double loc_sum, x;
     
     pthread_mutex_lock(&mtx);                                                 // { mutex lock
     myID = totalThread;
     totalThread++;
     if (totalThread<(*(int*)arg)) pthread_cond_wait(&cond, &mtx);
     else pthread_cond_broadcast(&cond);
     pthread_mutex_unlock(&mtx);                                               // mutex unlock }
     
     loc_size = ( num_steps + totalThread - 1 ) / totalThread;
     lbound = 1 + loc_size * myID;
     ubound = lbound + loc_size - 1;
     if (ubound>num_steps) ubound = num_steps;
     for (i=lbound;i<= ubound; i++) {
         x = (i-0.5)*step;
         pthread_mutex_lock(&mtx);                                             // { mutex lock
         sum = sum + 4.0/(1.0+x*x);                                            // mutex variable is sum
         pthread_mutex_unlock(&mtx);                                           // mutex unlock }
     }
     
     return (void*)0;
}

double pthread_PI_syn() {
       int i, thread_num=sysconf(_SC_NPROCESSORS_ONLN); 
       double x;
       
       sum = 0.0;
       step = 1.0/(double) num_steps;
       totalThread = 0;
       pthread_cond_init(&cond, NULL);
       pthread_mutex_init(&mtx, NULL);
       if ( thread_num>Max_Thread_Num ) thread_num = Max_Thread_Num;
       for(i=0; i<thread_num; i++) 
           pthread_create(&(threads[i]), NULL, syn_worker, &thread_num);     // create threads
       
       for(i=0; i<thread_num; i++) pthread_join(threads[i], NULL);
       pthread_cond_destroy(&cond);                                          // cond  destory
       pthread_mutex_destroy(&mtx);                                          // mutex destory
       return step * sum;
}

void *asyn_worker(void *arg) {
     int    myID, lbound, ubound, i, loc_size; 
     double loc_sum, x;
     
     pthread_mutex_lock(&mtx);                                               //{ mutex lock
     myID = totalThread;
     totalThread++;
     if (totalThread<(*(int*)arg)) pthread_cond_wait(&cond, &mtx);
     else pthread_cond_broadcast(&cond);
     pthread_mutex_unlock(&mtx);                                             // mutex unlock }
     
     loc_size = ( num_steps + totalThread - 1 ) / totalThread;
     lbound = 1 + loc_size * myID;
     ubound = lbound + loc_size - 1;
     if (ubound>num_steps) ubound = num_steps;
     loc_sum = 0;                                                            // initial local private varibale sum 
     for (i=lbound;i<= ubound; i++) {
         x = (i-0.5)*step;
         loc_sum = loc_sum + 4.0/(1.0+x*x);                                  // local sum   ! that's why program different
     }
     
     pthread_mutex_lock(&mtx);
     sum = sum + loc_sum;                                                    // critical area operation
     pthread_mutex_unlock(&mtx);         
     
    return (void*)0;
}

double pthread_PI_asyn() {
       int i, thread_num=sysconf(_SC_NPROCESSORS_ONLN); 
       double x;
       
       sum = 0.0;
       step = 1.0/(double) num_steps;
       totalThread = 0;
       pthread_cond_init(&cond, NULL);
       pthread_mutex_init(&mtx, NULL);
       if ( thread_num>Max_Thread_Num ) thread_num = Max_Thread_Num;
       for(i=0; i<thread_num; i++) 
           pthread_create(&(threads[i]), NULL, asyn_worker, &thread_num);
       
       for(i=0; i<thread_num; i++) pthread_join(threads[i], NULL);
       pthread_cond_destroy(&cond);
       pthread_mutex_destroy(&mtx);
       return step * sum;
}


int main (int argc, char* argv[]) { 
     double pi;
     struct timespec ts,te;
     double serial_cost, syn_cost, asyn_cost;
     
     
     
     clock_gettime(CLOCK_REALTIME, &ts);
     pi = serial_PI();
     clock_gettime(CLOCK_REALTIME, &te);
     serial_cost = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;
     printf("serial: PI=%20.18f  cost=%-15.10f\n", pi, serial_cost);
     
     clock_gettime(CLOCK_REALTIME, &ts);
     pi = pthread_PI_syn();
     clock_gettime(CLOCK_REALTIME, &te);
     syn_cost = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;
     printf("syn:    PI=%20.18f  cost=%-15.10f    speedup=%f\n", 
                                         pi, syn_cost, serial_cost/syn_cost);

     clock_gettime(CLOCK_REALTIME, &ts);
     pi = pthread_PI_asyn();
     clock_gettime(CLOCK_REALTIME, &te);
     asyn_cost = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;
     printf("asyn:   PI=%20.18f  cost=%-15.10f    speedup=%f\n", 
                                         pi, asyn_cost, serial_cost/asyn_cost);
     
     return EXIT_SUCCESS;
}
