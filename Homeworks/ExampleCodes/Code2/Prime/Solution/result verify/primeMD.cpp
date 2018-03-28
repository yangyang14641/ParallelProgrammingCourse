// head files
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h> 
//#include <math.h>
#include <stdlib.h>
#include <pthread.h>



// define parameters

/*
#define NANO           1000000000
#define Max_Thread_Num 256                                                  // define using how many threads
#define MAXIMUM        0x7fffffffffffffff
#define BLOCK_SIZE     65536


// global vars
long int n = 30000000;                                                      // how many prime number
*/


#define NANO           1000000000
#define Max_Thread_Num 4
#define MAXIMUM        0x7fffffffffffffff
#define BLOCK_SIZE     10


long int n = 200;                                                           // how many prime number


long int *vecPrime, nPrime, maximum; 


// worker's status vars struct
struct WorkerStatus {
       pthread_t  id;
       long int   *vecPrime;
       int        nPrime;
       int        maximum;
}  threads[Max_Thread_Num];


// Thread manager vars
int totalThread;
pthread_cond_t   cMaster, cWorker;
pthread_mutex_t  mtx;
long int lbound, ubound, task_size;




void serial_prime(long int arg) {
    long int  i, j, k, *temp;
    
    maximum = BLOCK_SIZE;
    vecPrime = (long int*)malloc(maximum*sizeof(long int));
    vecPrime[0] = 2;
    vecPrime[1] = 3;
    j = 2;
    nPrime = 1;
    lbound = 5;
    while (lbound<arg) {
         
         ubound = vecPrime[nPrime]*vecPrime[nPrime];                    
         if ( ubound<0 || ubound>arg ) ubound = arg;
         
         for ( i=lbound; i<ubound; i+=2 ) {
            while ( i>vecPrime[j-1]*vecPrime[j-1] ) j++;                        // define the search area              
            for ( k=1; k<j; k++ ) if ( i%vecPrime[k] == 0 ) break;
            if ( k<j ) continue;                                                // This is not a prime number
            nPrime++;                                                           // we find a prime because it can't be division by...
            if ( nPrime==maximum ) {                                            // alloc memory for new prime numbers 
               maximum += BLOCK_SIZE;
               temp = vecPrime;
               vecPrime = (long int*) malloc(maximum*sizeof(long int));
               memcpy(vecPrime, temp, (maximum - BLOCK_SIZE) * sizeof(long int));
               free(temp);
            }
            vecPrime[nPrime]=i;                                                  // store prime number
         }
         
         lbound = ubound + 2;
    }
    nPrime++;
    
    printf("\nPrime serial:\n");
    printf("\nAfter sort:\n");

    for (i = 0; i < nPrime; i++)
    printf("No. %ld prime Number: %ld\n",i,*(vecPrime+i)); 
    
    return;
} 

int cmpLongInt(const void *p1, const void *p2) {
    long int val1 = *((long int *)p1), val2 = *((long int *)p2);
    if ( val1<val2 ) return -1;
    if ( val1==val2 ) return 0;
    return 1;
}
   
void *mtx_worker(void *arg) {
     long int i, j, k, *temp;
     long int loc_lbound, loc_ubound;
     WorkerStatus  *pMyStatus=threads;
     
     while ( pMyStatus->id!=pthread_self() ) pMyStatus++; 
     pMyStatus->maximum = BLOCK_SIZE;
     pMyStatus->vecPrime = (long int*) malloc(pMyStatus->maximum*sizeof(long int));
     
     pthread_mutex_lock(&mtx);     
     totalThread++;          
     pthread_cond_broadcast(&cWorker); 
     j = 2;    
     while(lbound<MAXIMUM) {
          pthread_cond_wait(&cMaster, &mtx);
          if ( lbound==MAXIMUM ) continue;
          pthread_mutex_unlock(&mtx);
          
          pMyStatus->nPrime = 0;
          while ( true ) {
          	pthread_mutex_lock(&mtx);
                loc_lbound = lbound;
          	lbound += task_size;
//          	loc_lbound = lbound;                 // There is a problem in original sorce code
          	pthread_mutex_unlock(&mtx);
          	if ( loc_lbound >= ubound ) break;
          	loc_ubound = loc_lbound + task_size;
          	if ( loc_ubound > ubound ) loc_ubound = ubound;
          	
          	for ( i=loc_lbound; i<loc_ubound; i+=2 ) {
          	    while ( i>vecPrime[j-1]*vecPrime[j-1] ) j++;
          	    for ( k=1; k<j; k++ ) if ( i%vecPrime[k] == 0 ) break;
          	    if ( k<j ) continue;
          	    pMyStatus->vecPrime[pMyStatus->nPrime]=i;
          	    pMyStatus->nPrime++;
          	    if ( pMyStatus->nPrime == pMyStatus->maximum ) {
          	       temp = pMyStatus->vecPrime; 
          	       pMyStatus->maximum += BLOCK_SIZE;
          	       pMyStatus->vecPrime = (long int*) malloc(pMyStatus->maximum*sizeof(long int));
          	       memcpy(pMyStatus->vecPrime, temp, pMyStatus->nPrime*sizeof(long int));
          	       free(temp);
          	    }          	    
          	}
          }         
          
          pthread_mutex_lock(&mtx);
          totalThread++;
          pthread_cond_broadcast(&cWorker);
     }
     pthread_mutex_unlock(&mtx);         
     free(pMyStatus->vecPrime);
    return (void*)0;
}

void pthread_prime_mtx(long int arg) {
     int thread_num=sysconf(_SC_NPROCESSORS_ONLN);
     long int  i, j, k, *temp;    

     maximum = BLOCK_SIZE;
     vecPrime = (long int*)malloc(maximum*sizeof(long int));
     vecPrime[0] = 2;
     vecPrime[1] = 3;
     j = 2;
     nPrime = 2;
     lbound = 5;


     totalThread = 0;
     pthread_mutex_init(&mtx, NULL);
     pthread_cond_init(&cWorker, NULL);
     pthread_cond_init(&cMaster, NULL);
     if ( thread_num>Max_Thread_Num ) thread_num = Max_Thread_Num;
     for(i=0; i<thread_num; i++) pthread_create(&(threads[i].id), NULL, mtx_worker, NULL);          // create worker threads
     pthread_mutex_lock(&mtx);
     while ( totalThread!=thread_num ) pthread_cond_wait(&cWorker, &mtx);       // barrier master threads wait for worker threads
     
     while (lbound<arg) {
         ubound = vecPrime[nPrime-1]*vecPrime[nPrime-1];                      //  task division for worker threads
         if ( ubound<0 || ubound>arg ) ubound = arg;
         task_size = (ubound-lbound)/(10*thread_num);
         if ( task_size<10 ) task_size = (ubound-lbound+thread_num-1)/thread_num;
         if ( task_size%2==1 ) task_size++;
         
         totalThread = 0;
         pthread_cond_broadcast(&cMaster);   // broadcast of master thread state to worker threads tell them task division finished       
         while ( totalThread!=thread_num ) pthread_cond_wait(&cWorker, &mtx);  // barried the master thread wait for worker threads
       
         // this part is to sort after every thread finished 
         for (i=0; i<thread_num; i++) {
             if ( threads[i].nPrime==0 ) continue;              // means this thread didn't find any prime
             if ( nPrime+threads[i].nPrime>maximum ) {
             	temp = vecPrime;
             	while( nPrime+threads[i].nPrime>maximum ) maximum += BLOCK_SIZE;
             	vecPrime = (long int*) malloc(maximum*sizeof(long int));
             	memcpy(vecPrime, temp, nPrime*sizeof(long int));
             	free(temp);
             }
             memcpy(vecPrime+nPrime, threads[i].vecPrime, threads[i].nPrime*sizeof(long int)); 
             nPrime += threads[i].nPrime;       
         }
         
//         qsort(vecPrime, nPrime, sizeof(long int), cmpLongInt);   // quickly sort
         lbound = ubound + 2;
     }
     
      qsort(vecPrime, nPrime, sizeof(long int), cmpLongInt);    // quickly sort

      
     printf("\nPrime mutex:\n");
     printf("\nAfter sort:\n");
     for (i = 0; i < nPrime; i++)
     printf("No. %ld prime Number: %ld\n",i,*(vecPrime+i)); 

     lbound = MAXIMUM;
     pthread_cond_broadcast(&cMaster);
     pthread_mutex_unlock(&mtx);
     
     for(i=0; i<thread_num; i++) pthread_join(threads[i].id, NULL);
     pthread_mutex_destroy(&mtx);
     pthread_cond_destroy(&cWorker);
     pthread_cond_destroy(&cMaster);
     return;
}

void *atomic_worker(void *arg) {
     long int i, j, k, *temp;
     long int loc_lbound, loc_ubound;
     WorkerStatus  *pMyStatus=threads;
     
     while ( pMyStatus->id!=pthread_self() ) pMyStatus++; 
     pMyStatus->maximum = BLOCK_SIZE;
     pMyStatus->vecPrime = (long int*) malloc(pMyStatus->maximum*sizeof(long int));
     
     pthread_mutex_lock(&mtx);     
     totalThread++;          
     pthread_cond_broadcast(&cWorker); 
     j = 2;    
     while(lbound<MAXIMUM) {
          pthread_cond_wait(&cMaster, &mtx);
          if ( lbound==MAXIMUM ) continue;
          pthread_mutex_unlock(&mtx);
          
          pMyStatus->nPrime = 0;
          while ( true ) {
          	loc_lbound = __sync_fetch_and_add(&lbound, task_size);
          	if ( loc_lbound >= ubound ) break;
          	loc_ubound = loc_lbound + task_size;
          	if ( loc_ubound > ubound ) loc_ubound = ubound;
          	
          	for ( i=loc_lbound; i<loc_ubound; i+=2 ) {
          	    while ( i>vecPrime[j-1]*vecPrime[j-1] ) j++;
          	    for ( k=1; k<j; k++ ) if ( i%vecPrime[k] == 0 ) break;
          	    if ( k<j ) continue;
          	    pMyStatus->vecPrime[pMyStatus->nPrime]=i;
          	    pMyStatus->nPrime++;
          	    if ( pMyStatus->nPrime == pMyStatus->maximum ) {
          	       temp = pMyStatus->vecPrime; 
          	       pMyStatus->maximum += BLOCK_SIZE;
          	       pMyStatus->vecPrime = (long int*) malloc(pMyStatus->maximum*sizeof(long int));
          	       memcpy(pMyStatus->vecPrime, temp, pMyStatus->nPrime*sizeof(long int));
          	       free(temp);
          	    }          	    
          	}
          }         
          
          pthread_mutex_lock(&mtx);
          totalThread++;
          pthread_cond_broadcast(&cWorker);
     }
     pthread_mutex_unlock(&mtx);         
     free(pMyStatus->vecPrime);
    return (void*)0;
}

void pthread_prime_atomic(long int arg) {
     int thread_num=sysconf(_SC_NPROCESSORS_ONLN);
     long int  i, j, k, *temp;     

     maximum = BLOCK_SIZE;
     vecPrime = (long int*)malloc(maximum*sizeof(long int));
     vecPrime[0] = 2;
     vecPrime[1] = 3;
     j = 2;
     nPrime = 2;
     lbound = 5;

     totalThread = 0;
     pthread_mutex_init(&mtx, NULL);
     pthread_cond_init(&cWorker, NULL);
     pthread_cond_init(&cMaster, NULL);
     if ( thread_num>Max_Thread_Num ) thread_num = Max_Thread_Num;
     for(i=0; i<thread_num; i++) pthread_create(&(threads[i].id), NULL, atomic_worker, NULL);
     pthread_mutex_lock(&mtx);
     while ( totalThread!=thread_num ) pthread_cond_wait(&cWorker, &mtx);    
     
     while (lbound<arg) {
         ubound = vecPrime[nPrime-1]*vecPrime[nPrime-1];
         if ( ubound<0 || ubound>arg ) ubound = arg;
         task_size = (ubound-lbound)/(10*thread_num);
         if ( task_size<10 ) task_size = (ubound-lbound+thread_num-1)/thread_num;
         if ( task_size%2==1 ) task_size++;
         
         totalThread = 0;
         pthread_cond_broadcast(&cMaster);         
         while ( totalThread!=thread_num ) pthread_cond_wait(&cWorker, &mtx);
         for (i=0; i<thread_num; i++) {
             if ( threads[i].nPrime==0 ) continue;
             if ( nPrime+threads[i].nPrime>maximum ) {
             	temp = vecPrime;
             	while( nPrime+threads[i].nPrime>maximum ) maximum += BLOCK_SIZE;
             	vecPrime = (long int*) malloc(maximum*sizeof(long int));
             	memcpy(vecPrime, temp, nPrime*sizeof(long int));
             	free(temp);
             }
             memcpy(vecPrime+nPrime, threads[i].vecPrime, threads[i].nPrime*sizeof(long int)); 
             nPrime += threads[i].nPrime;       
         }
         
//         qsort(vecPrime, nPrime, sizeof(long int), cmpLongInt);
         lbound = ubound + 2;
     }

     qsort(vecPrime, nPrime, sizeof(long int), cmpLongInt);

     printf("\nPrime atomic:\n");
     printf("\nAfter sort:\n");

     for (i = 0; i < nPrime; i++)
     printf("No. %ld prime Number: %ld\n",i,*(vecPrime+i)); 

     lbound = MAXIMUM;
     pthread_cond_broadcast(&cMaster);
     pthread_mutex_unlock(&mtx);
     
     for(i=0; i<thread_num; i++) pthread_join(threads[i].id, NULL);
     pthread_mutex_destroy(&mtx);
     pthread_cond_destroy(&cWorker);
     pthread_cond_destroy(&cMaster);
     return;
}

void *dup_worker(void *arg) {
     long int i, j, k, *temp, loc_vecPrime[BLOCK_SIZE];
     long int loc_lbound, loc_ubound, loc_nPrime;
     WorkerStatus  *pMyStatus=threads;
     
     while ( pMyStatus->id!=pthread_self() ) pMyStatus++;
     
     pthread_mutex_lock(&mtx);     
     totalThread++;          
     pthread_cond_broadcast(&cWorker); 
     j = 2;    
     while(lbound<MAXIMUM) {
          pthread_cond_wait(&cMaster, &mtx);
          if ( lbound==MAXIMUM ) continue;
          pthread_mutex_unlock(&mtx);
          
          loc_nPrime = 0;
          while ( true ) {
          	loc_lbound = __sync_fetch_and_add(&lbound, task_size);
          	if ( loc_lbound >= ubound ) break;
          	loc_ubound = loc_lbound + task_size;
          	if ( loc_ubound > ubound ) loc_ubound = ubound;
          	
          	for ( i=loc_lbound; i<loc_ubound; i+=2 ) {
          	    while ( i>pMyStatus->vecPrime[j-1]*pMyStatus->vecPrime[j-1] ) j++;
          	    for ( k=1; k<j; k++ ) if ( i%pMyStatus->vecPrime[k] == 0 ) break;
          	    if ( k<j ) continue;
          	    loc_vecPrime[loc_nPrime]=i;
          	    loc_nPrime++;
          	    if ( loc_nPrime == BLOCK_SIZE ) {
          	       pthread_mutex_lock(&mtx);
          	       temp = vecPrime;          	       
          	       maximum += BLOCK_SIZE;
          	       vecPrime = (long int*) malloc(maximum*sizeof(long int));
          	       memcpy(vecPrime, temp, nPrime*sizeof(long int));
          	       memcpy(vecPrime+nPrime, loc_vecPrime, BLOCK_SIZE*sizeof(long int));
          	       nPrime += BLOCK_SIZE;
          	       pthread_mutex_unlock(&mtx);
          	       free(temp);
          	       loc_nPrime = 0;
          	    }          	    
          	}
          } 
          pthread_mutex_lock(&mtx);        
          if ( loc_nPrime>0 ) {
             if ( nPrime + loc_nPrime > maximum ) {
                temp = vecPrime;
                maximum += BLOCK_SIZE;
                vecPrime = (long int*) malloc(maximum*sizeof(long int));
                memcpy(vecPrime, temp, nPrime*sizeof(long int));
             }
             memcpy(vecPrime+nPrime, loc_vecPrime, loc_nPrime*sizeof(long int));
             nPrime += loc_nPrime;
          }
          totalThread++;
          pthread_cond_broadcast(&cWorker);
     }
     pthread_mutex_unlock(&mtx);         
     return (void*)0;
}

void pthread_prime_dup(long int arg) {
     int thread_num=sysconf(_SC_NPROCESSORS_ONLN);
     long int  i, j, k, *temp;     

     maximum = BLOCK_SIZE;
     vecPrime = (long int*)malloc(maximum*sizeof(long int));
     vecPrime[0] = 2;
     vecPrime[1] = 3;
     j = 2;
     nPrime = 2;
     lbound = 5;

     totalThread = 0;
     pthread_mutex_init(&mtx, NULL);
     pthread_cond_init(&cWorker, NULL);
     pthread_cond_init(&cMaster, NULL);
     if ( thread_num>Max_Thread_Num ) thread_num = Max_Thread_Num;
     for(i=0; i<thread_num; i++) {
        threads[i].vecPrime = NULL;
        pthread_create(&(threads[i].id), NULL, dup_worker, NULL);
     }
     pthread_mutex_lock(&mtx);                                                         
     while ( totalThread!=thread_num ) pthread_cond_wait(&cWorker, &mtx);    
     
     while (lbound<arg) {
         ubound = vecPrime[nPrime-1]*vecPrime[nPrime-1];
         if ( ubound<0 || ubound>arg ) ubound = arg;
         task_size = (ubound-lbound)/(10*thread_num);
         if ( task_size<10 ) task_size = (ubound-lbound+thread_num-1)/thread_num;
         if ( task_size%2==1 ) task_size++; 
         
         totalThread = 0;
         for (i=0; i<thread_num; i++) {
             if ( threads[i].vecPrime!=NULL ) free(threads[i].vecPrime);
             threads[i].nPrime = nPrime;
             threads[i].vecPrime = (long int*) malloc(nPrime*sizeof(long int));
             memcpy(threads[i].vecPrime, vecPrime, nPrime*sizeof(long int));
         }
         pthread_cond_broadcast(&cMaster);
         
         while ( totalThread!=thread_num ) pthread_cond_wait(&cWorker, &mtx);
         
//         qsort(vecPrime, nPrime, sizeof(long int), cmpLongInt);
         lbound = ubound + 2;
     }

     qsort(vecPrime, nPrime, sizeof(long int), cmpLongInt);
     
     printf("\nPrime dup:\n");
     printf("\nAfter sort:\n");

     for (i = 0; i < nPrime; i++)
     printf("No. %ld prime Number: %ld\n",i,*(vecPrime+i)); 
     
     lbound = MAXIMUM;
     pthread_cond_broadcast(&cMaster);
     pthread_mutex_unlock(&mtx);
     
     for(i=0; i<thread_num; i++) {
        pthread_join(threads[i].id, NULL);
        if ( threads[i].vecPrime != NULL ) free(threads[i].vecPrime);
     }



     pthread_mutex_destroy(&mtx);
     pthread_cond_destroy(&cWorker);
     pthread_cond_destroy(&cMaster);
     return;
}

int main(int argc, char* argv[] ){
    struct timespec ts,te;
    double serial_cost, mtx_cost, atomic_cost, dup_cost;
    
    if (n<=0) n=MAXIMUM;
    
    clock_gettime(CLOCK_REALTIME, &ts);
    serial_prime(n);
    free(vecPrime);
    clock_gettime(CLOCK_REALTIME, &te);
    serial_cost = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;
    printf("serial: found %ld primes  cost = %15.10f \n", nPrime, serial_cost);

    clock_gettime(CLOCK_REALTIME, &ts);
    pthread_prime_mtx(n);
    free(vecPrime);
    clock_gettime(CLOCK_REALTIME, &te);
    mtx_cost = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;
    printf("mtx   : found %ld primes  cost = %15.10f    speedup = %f \n", nPrime, mtx_cost, serial_cost/mtx_cost);
    
    clock_gettime(CLOCK_REALTIME, &ts);
    pthread_prime_atomic(n);
    free(vecPrime);
    clock_gettime(CLOCK_REALTIME, &te);
    atomic_cost = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;
    printf("atomic: found %ld primes  cost = %15.10f    speedup = %f \n", nPrime, atomic_cost, serial_cost/atomic_cost);

    clock_gettime(CLOCK_REALTIME, &ts);
    pthread_prime_dup(n);
    free(vecPrime);
    clock_gettime(CLOCK_REALTIME, &te);
    dup_cost = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;
    printf("dup   : found %ld primes  cost = %15.10f    speedup = %f \n", nPrime, dup_cost, serial_cost/dup_cost);
    
    return EXIT_SUCCESS;
}

