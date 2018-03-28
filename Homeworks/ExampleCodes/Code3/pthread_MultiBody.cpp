// Head files
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h> 
#include <unistd.h>
#include <pthread.h>

// Defines
#define NANO           1000000000
#define Max_Thread_Num 4
#define REAL           double

// Global vars
int BodyNum=0;
int TimeSteps=0;
REAL *body;
REAL *force;
struct STATUS {
       int lbound;
       int ubound;
       int *task;
       REAL *force;
       pthread_mutex_t  mtx;
} *status;   

// Pthread vars
int totalThread;
pthread_cond_t   cond;
pthread_mutex_t  mtx;

// Function declines 
double serial(); 

double pthread_reduce();
void *reduce_worker(void *arg);
int  freeWorker;  

double pthread_divide();
void *divide_worker(void *arg);
                                                             



//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
// The main function                                                                   
int main(int argc, char** argv ) {     
  int i;
  REAL ser_time, red_time, div_time;
  char *pStr;
  
  for ( i=1; i<argc; i++ ) {                                 // program parameters input in command line
	  pStr=strstr(argv[i], "-s=");
	  if ( pStr!=NULL) sscanf(pStr, "-s=%d", &BodyNum);
	  pStr=strstr(argv[i], "-t=");
	  if ( pStr!=NULL) sscanf(pStr, "-t=%d", &TimeSteps);

  }
  if ( BodyNum*TimeSteps==0) {
	  printf("usage: -s=number-of-bodies -t=number-of-steps\n");
	  return 0;
  }
  
  printf("Problem Scale: number of bodies: %d, number of time step: %d\n", BodyNum, TimeSteps);

  body = (REAL*)malloc(4*BodyNum*sizeof(REAL));
  
  ser_time = serial();
  printf("serial multibody: %f (sec)\n", ser_time);
  red_time = pthread_reduce();
  printf("prallel multibody pthread_reduce : %f (sec), speedup = %f\n", red_time, ser_time/red_time); 
  div_time = pthread_divide();
  printf("parallel multibody pthread_divide : %f (sec), speedup = %f\n", div_time, ser_time/div_time);

  free(body);  
  
}
//-------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------
// function pthread_divide
double pthread_divide() {
       pthread_t threads[Max_Thread_Num];
       int i, j, loc_size, task_num;
       int thread_num=sysconf(_SC_NPROCESSORS_ONLN); 
       struct timespec ts,te;
       double result;         
       
       /*  Initialize mass and positions in array p to make a test case   */
       for ( i=0; i<BodyNum; i++)    {
          body[4*i] = 10.05 + i;
          body[4*i+1] = 30.0*i;
          body[4*i+2] = 20.0*i;
          body[4*i+3] = 10.0*i;
       }
       
       clock_gettime(CLOCK_REALTIME, &ts);
       force = (REAL*)malloc(3*BodyNum*sizeof(REAL));
       for ( i=0; i<3*BodyNum; i++) force[i] = 0;
       totalThread = 0;
       freeWorker = 0;
       pthread_cond_init(&cond, NULL);                                       // pthread cond initial
       pthread_mutex_init(&mtx, NULL);                                       // pthread mutex initial
       if ( thread_num>Max_Thread_Num ) thread_num = Max_Thread_Num;
       
       // Computing task division
       loc_size = ( BodyNum + thread_num - 1 ) / thread_num;                  // computing local problem size
       status = (STATUS*)malloc(thread_num*sizeof(STATUS));                   // malloc sturcture for each thread
       task_num = thread_num/2+1;                                             // computing number of task (K^2)/2 / K + 1(diag)
       for(i=0; i<thread_num; i++) {                                           
          status[i].lbound = i*loc_size;                                      // each thread's low boundary
          status[i].ubound = status[i].lbound + loc_size;                     // each thread's up boundary
          if ( status[i].ubound>BodyNum ) status[i].ubound = BodyNum;         // up boundary modify for last thread
          status[i].force = (REAL*)malloc(3*loc_size*sizeof(REAL));           // allocate memory for each thread's particles' force
          status[i].task = (int*)malloc(task_num*sizeof(int));                // allocate memory for each thread's task 
          for(j=0; j<task_num; j++)               
             status[i].task[j] = (i+j)%thread_num;                            // determine subproblem task size
          pthread_mutex_init(&(status[i].mtx), NULL);                         // pthread mutex init
       }
       status[thread_num-1].ubound = BodyNum;

       // Create worker threads
       for(i=0; i<thread_num; i++) pthread_create(&(threads[i]), NULL, divide_worker, &thread_num);

       // Waiting untill Join Worker threads
       for(i=0; i<thread_num; i++) pthread_join(threads[i], NULL);
       
       // Destory pthread vars
       pthread_cond_destroy(&cond);
       pthread_mutex_destroy(&mtx);
       for(i=0; i<thread_num; i++) {
          free(status[i].task);
          free(status[i].force); 
          pthread_mutex_destroy(&(status[i].mtx));
       }
       
       // free memory 
       free(force);
       free(status);
       
       // time up
       clock_gettime(CLOCK_REALTIME, &te);
       result = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;
       
       // data outputs
       FILE *fResult=fopen("result_div_nbody.txt", "w");
       char str[50];
       for (i=0; i<BodyNum; i++)   { 	  
	  sprintf(str, "(%10.4f %10.4f %10.4f %10.4f)\n", body[4*i], body[4*i+1], body[4*i+2], body[4*i+3]);
	  fwrite(str, sizeof(char), strlen(str), fResult);
       }
       fclose(fResult);
       return result;

}
//-------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------
// function execute_task
void execute_task(int n) {
     REAL fac, fx, fy, fz;
     REAL dx, dy, dz, sq, dist; 
     int  i, j, bi,bj,fi,fj;
     
     for ( i=status[n].lbound; i<status[n].ubound; i++ ) {
         bi = 4*i;              
         fi = 3*i;              
         for ( j=i+1; j<status[n].ubound; j++ )  {
             bj = 4*j;
             fj = 3*j;
             /*Calculate force between particle i and j according to Newton's Law*/
             dx = body[bi+1] - body[bj+1];
             dy = body[bi+2] - body[bj+2];
             dz = body[bi+3] - body[bj+3];
             sq = dx*dx + dy*dy + dz*dz;
             dist = sqrt(sq);
             fac = body[bi] * body[bj] / ( dist * sq );
             fx = fac * dx;
             fy = fac * dy;
             fz = fac * dz;
             /*Add in force and opposite force to particle i and j */
             force[fi] -= fx;
             force[fi+1] -= fy;
             force[fi+2] -= fz;
             force[fj] += fx;
             force[fj+1] += fy;
             force[fj+2] += fz;
         }
     }
     return; 
}
//-------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------
// function execute_task
void execute_task(int n, int m) {
     REAL fac, fx, fy, fz;
     REAL dx, dy, dz, sq, dist; 
     int  i, j, bi,bj,fi,fj;
     
     for ( i=0; i<3*(status[m].ubound-status[m].lbound); i++ ) status[n].force[i] = 0;
     for ( i=status[n].lbound; i<status[n].ubound; i++ ) {
         bi = 4*i;              
         fi = 3*i;              
         for ( j=status[m].lbound; j<status[m].ubound; j++ )  {
             bj = 4*j;
             fj = 3*(j-status[m].lbound);
             /*Calculate force between particle i and j according to Newton's Law*/
             dx = body[bi+1] - body[bj+1];
             dy = body[bi+2] - body[bj+2];
             dz = body[bi+3] - body[bj+3];
             sq = dx*dx + dy*dy + dz*dz;
             dist = sqrt(sq);
             fac = body[bi] * body[bj] / ( dist * sq );
             fx = fac * dx;
             fy = fac * dy;
             fz = fac * dz;
             /*Add in force and opposite force to particle i and j */
             force[fi] -= fx;
             force[fi+1] -= fy;
             force[fi+2] -= fz;
             status[n].force[fj] += fx;
             status[n].force[fj+1] += fy;
             status[n].force[fj+2] += fz;
         }
     }
     return; 
}
//-------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------
// function execute_semi_task
void execute_semi_task(int n, int m) {
     REAL fac, fx, fy, fz;
     REAL dx, dy, dz, sq, dist; 
     int  i, j, bi,bj,fi,fj;
     
     for ( i=status[n].lbound; i<status[n].ubound; i++ ) {
         bi = 4*i;              
         fi = 3*i;              
         for ( j=status[m].lbound; j<status[m].ubound; j++ )  {
             bj = 4*j;
             fj = 3*(j-status[m].lbound);
             /*Calculate force between particle i and j according to Newton's Law*/
             dx = body[bi+1] - body[bj+1];
             dy = body[bi+2] - body[bj+2];
             dz = body[bi+3] - body[bj+3];
             sq = dx*dx + dy*dy + dz*dz;
             dist = sqrt(sq);
             fac = body[bi] * body[bj] / ( dist * sq );
             fx = fac * dx;
             fy = fac * dy;
             fz = fac * dz;
             /*Add in force and opposite force to particle i */
             force[fi] -= fx;
             force[fi+1] -= fy;
             force[fi+2] -= fz;
         }
     }
     return; 
}


//-------------------------------------------------------------------------------------------------
// function divide_worker
void *divide_worker(void *arg) {
     int t, i, j, bi,bj,fi,fj;
     
     int    myID, thread_num; 
     int    task_num, segID, lbound;
     
     pthread_mutex_lock(&mtx);
     myID = totalThread;
     totalThread++;
     pthread_mutex_unlock(&mtx);
     thread_num = (int)(*(int*)arg);
       
       // time step loop
       t = 0;
       while ( t<TimeSteps){
           /*  Loop over points calculating force between each pair.*/
           pthread_mutex_lock(&(status[myID].mtx));
           execute_task(myID);
           if (thread_num%2==0) {
              task_num = thread_num/2;
              execute_semi_task(myID, status[myID].task[task_num]);
           }
           else task_num = thread_num/2+1;
           
           for( i=1;i<task_num; i++) {
               execute_task(myID, status[myID].task[i]);
               pthread_mutex_unlock(&(status[myID].mtx));
               segID = status[myID].task[i];
               pthread_mutex_lock(&(status[segID].mtx));               
               lbound = 3*status[segID].lbound;
               for(j=0; j<3*(status[segID].ubound-status[segID].lbound); j++) force[j+lbound] += status[myID].force[j];
               pthread_mutex_unlock(&(status[segID].mtx));
               pthread_mutex_lock(&(status[myID].mtx));
           }
           
           // computing new position
           for ( i=status[myID].lbound; i<status[myID].ubound; i++ ){ 
              bi = 4*i;              
              fi = 3*i;
              body[bi+1] = body[bi+1] + force[fi] / body[bi];
              force[fi] = 0;
              body[bi+2] = body[bi+2] + force[fi+1] / body[bi];
              force[fi+1] = 0;
              body[bi+3] = body[bi+3] + force[fi+2] / body[bi];
              force[fi+2] = 0;
           }
           pthread_mutex_unlock(&(status[myID].mtx));
           pthread_mutex_lock(&mtx);
           freeWorker++;
           if ( freeWorker<thread_num ) pthread_cond_wait(&cond, &mtx);
           else {
             freeWorker = 0;
             pthread_cond_broadcast(&cond);
           }
           pthread_mutex_unlock(&mtx);

           t++;
       }
     return (void*)0;
}


//-------------------------------------------------------------------------------------------------
// function pthread_reduce
double pthread_reduce() {
       pthread_t threads[Max_Thread_Num];
       int i, thread_num=sysconf(_SC_NPROCESSORS_ONLN); 
       struct timespec ts,te;
       double result;         
       
       /*  Initialize mass and positions in array p to make a test case   */
       for ( i=0; i<BodyNum; i++)    {
          body[4*i] = 10.05 + i;
          body[4*i+1] = 30.0*i;
          body[4*i+2] = 20.0*i;
          body[4*i+3] = 10.0*i;
       }
       
       clock_gettime(CLOCK_REALTIME, &ts);
       force = (REAL*)malloc(3*BodyNum*sizeof(REAL));                                // alocation memory for force array
       for ( i=0; i<3*BodyNum; i++) force[i] = 0;
       totalThread = 0;
       freeWorker = 0;
       pthread_cond_init(&cond, NULL);                                               // initial pthread condiction var
       pthread_mutex_init(&mtx, NULL);                                               // initial pthread mutex var
       if ( thread_num>Max_Thread_Num ) thread_num = Max_Thread_Num;                 // maximum threads number 

       // Create pthread worker threads
       for(i=0; i<thread_num; i++) pthread_create(&(threads[i]), NULL, reduce_worker, &thread_num);
       
       // Waiting and join pthread worker threads
       for(i=0; i<thread_num; i++) pthread_join(threads[i], NULL);

       // destory pthread's vars
       pthread_cond_destroy(&cond);
       pthread_mutex_destroy(&mtx);
       
       // free memory
       free(force);
       
       // time up
       clock_gettime(CLOCK_REALTIME, &te);
       result = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;
       
       FILE *fResult=fopen("result_red_nbody.txt", "w");
       char str[50];
       for (i=0; i<BodyNum; i++)   { 	  
	  sprintf(str, "(%10.4f %10.4f %10.4f %10.4f)\n", body[4*i], body[4*i+1], body[4*i+2], body[4*i+3]);
	  fwrite(str, sizeof(char), strlen(str), fResult);
       }
       fclose(fResult);
       return result;
}


//-------------------------------------------------------------------------------------------------
// function reduce_worker
void *reduce_worker(void *arg) {
     REAL *locForce;
     REAL fac, fx, fy, fz;
     REAL dx, dy, dz, sq, dist; 
     int t, i, j, bi,bj,fi,fj;
     
     int    myID, thread_num, lbound, ubound, loc_size; 
     
     pthread_mutex_lock(&mtx);
     myID = totalThread;
     totalThread++;
     pthread_mutex_unlock(&mtx);
     thread_num = (int)(*(int*)arg);
     
     loc_size = ( BodyNum + thread_num - 1 ) / thread_num;                // thread's local problem size
     lbound = loc_size * myID;                                            // thread's local problem low boundary
     ubound = lbound + loc_size;                                          // thread's local problem up boundary
     if (myID==thread_num - 1) ubound = BodyNum;                          // up boundary modify

       locForce = (REAL*)malloc(3*BodyNum*sizeof(REAL));                  // allocation memory for threads' local force 
       for ( i=0; i<3*BodyNum; i++) locForce[i] = 0;
       
       // begin time loop
       t = 0;
       while ( t<TimeSteps){

           // Computing local force 
           /*  Loop over points calculating force between each pair. */
           for ( i=myID; i<BodyNum; i+=thread_num ) {
              bi = 4*i;              
              fi = 3*i;              
              for ( j=i+1; j<BodyNum; j++ )  {
                 bj = 4*j;
                 fj = 3*j;
                 /*Calculate force between particle i and j according to Newton's second Law */
                 dx = body[bi+1] - body[bj+1];
                 dy = body[bi+2] - body[bj+2];
                 dz = body[bi+3] - body[bj+3];
                 sq = dx*dx + dy*dy + dz*dz;
                 dist = sqrt(sq);
                 fac = body[bi] * body[bj] / ( dist * sq );
                 fx = fac * dx;
                 fy = fac * dy;
                 fz = fac * dz;
                 /*Add in force and opposite force to particle i and j */
                 locForce[fi] -= fx;
                 locForce[fi+1] -= fy;
                 locForce[fi+2] -= fz;
                 locForce[fj] += fx;
                 locForce[fj+1] += fy;
                 locForce[fj+2] += fz;
              }
           }
           
           // sum of local force 
           pthread_mutex_lock(&mtx);
           for ( i=0; i<3*BodyNum; i++ ) {
              force[i] += locForce[i];
              locForce[i] = 0;
           }
           freeWorker++;
           if ( freeWorker<thread_num ) pthread_cond_wait(&cond, &mtx);
           else {
             freeWorker = 0;
             pthread_cond_broadcast(&cond);
           }
           pthread_mutex_unlock(&mtx);
           
           // particles' position update
           for ( i=lbound; i<ubound; i++ ){ 
              bi = 4*i;              
              fi = 3*i;
              body[bi+1] = body[bi+1] + force[fi] / body[bi];
              force[fi] = 0;
              body[bi+2] = body[bi+2] + force[fi+1] / body[bi];
              force[fi+1] = 0;
              body[bi+3] = body[bi+3] + force[fi+2] / body[bi];
              force[fi+2] = 0;
           }
           
           pthread_mutex_lock(&mtx);
           freeWorker++;
           if ( freeWorker<thread_num ) pthread_cond_wait(&cond, &mtx);
           else {
             freeWorker = 0;
             pthread_cond_broadcast(&cond);
           }
           pthread_mutex_unlock(&mtx);

           t++;
       }
     free(locForce); 
     return (void*)0;
}


//-------------------------------------------------------------------------------------------------
// function serial
double serial() {
       REAL fac, fx, fy, fz;
       REAL dx, dy, dz, sq, dist; 
       int t, i, j, bi,bj,fi,fj;
       struct timespec ts,te;
       double result;
        
       /*  Initialize mass and positions in array p to make a test case   */
       for ( i=0; i<BodyNum; i++)    {
          body[4*i] = 10.05 + i;
          body[4*i+1] = 30.0*i;
          body[4*i+2] = 20.0*i;
          body[4*i+3] = 10.0*i;
       }
              
       clock_gettime(CLOCK_REALTIME, &ts);
          
       force = (REAL*)malloc(3*BodyNum*sizeof(REAL));
       for ( i=0; i<3*BodyNum; i++)  force[i] = 0;
       
       t = 0;
       while ( t<TimeSteps){
           /*  Loop over points calculating force between each pair.*/
           for ( i=0; i<BodyNum; i++ ) {
              bi = 4*i;              
              fi = 3*i;              
              for ( j=i+1; j<BodyNum; j++ )  {
                 bj = 4*j;
                 fj = 3*j;
                 /*Calculate force between particle i and j according to Newton's Law*/
                 dx = body[bi+1] - body[bj+1];
                 dy = body[bi+2] - body[bj+2];
                 dz = body[bi+3] - body[bj+3];
                 sq = dx*dx + dy*dy + dz*dz;
                 dist = sqrt(sq);
                 fac = body[bi] * body[bj] / ( dist * sq );
                 fx = fac * dx;
                 fy = fac * dy;
                 fz = fac * dz;
                 /*Add in force and opposite force to particle i and j */
                 force[fi] -= fx;
                 force[fi+1] -= fy;
                 force[fi+2] -= fz;
                 force[fj] += fx;
                 force[fj+1] += fy;
                 force[fj+2] += fz;
              }
           }
           for ( i=0; i<BodyNum; i++ ){ 
              bi = 4*i;              
              fi = 3*i;
              body[bi+1] = body[bi+1] + force[fi] / body[bi];
              force[fi] = 0;
              body[bi+2] = body[bi+2] + force[fi+1] / body[bi];
              force[fi+1] = 0;
              body[bi+3] = body[bi+3] + force[fi+2] / body[bi];
              force[fi+2] = 0;
           }
           t++;
       }
       free(force);
       
       clock_gettime(CLOCK_REALTIME, &te);
       result = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;
       
       FILE *fResult=fopen("result_ser_nbody.txt", "w");
       char str[50];
       for (i=0; i<BodyNum; i++)   { 	  
	  sprintf(str, "(%10.4f %10.4f %10.4f %10.4f)\n", body[4*i], body[4*i+1], body[4*i+2], body[4*i+3]);
	  fwrite(str, sizeof(char), strlen(str), fResult);
       }
       fclose(fResult);
       return result;
}
