/****************************************************************************
 * DESCRIPTION:  
 *   Serial HEAT2D Example - C Version
 *   This example is based on a simplified 
 *   two-dimensional heat equation domain decomposition.  The initial 
 *   temperature is computed to be high in the middle of the domain and 
 *   zero at the boundaries.  The boundaries are held at zero throughout 
 *   the simulation.  During the time-stepping, an array containing two 
 *   domains is used; these domains alternate between old data and new data.
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#define NANO           1000000000
#define Max_Thread_Num 256

#define NXPROB 110
#define NYPROB 110
struct Parms
{ 
  float cx;
  float cy;
  int nts;
} parms = {0.1, 0.1, 5000};

void   inidat(int nx, int ny, float *u1);
void   prtdat(int nx, int ny, float *u1, char fnam[]);

double serial(int nx, int ny, float *u);
void   update_serial(int nx, int ny, float *u1, float *u2);

void   *worker(void *arg);
double ptread_heat2D(int nx, int ny, float *u);

int main(int argc, char* argv[]) {
  float u[2][NXPROB][NYPROB];
  int ix, iy, iz, it;
  
  printf("\n Problem Scale: NXPROB = %d, NYPROB = %d， nts = %d\n\n", NXPROB, NYPROB, parms.nts);

  /************************************************************************
  **  Initialize grid.
  *************************************************************************/
  inidat(NXPROB, NYPROB, &u[0][0][0]);
  //prtdat(NXPROB, NYPROB, &u[0][0][0], "initial.dat");
  for (ix = 0; ix <= NXPROB-1; ix++) {
    u[1][ix][0] = u[0][ix][0];
    u[1][ix][NYPROB-1] = u[0][ix][NYPROB-1];
  }
  for (iy = 0; iy <= NYPROB-1; iy++) {
    u[1][0][iy] = u[0][0][iy];
    u[1][NXPROB-1][iy] = u[0][NXPROB-1][iy];
  }

  /***********************************************************************
  **  Iterate over all timesteps.
  ************************************************************************/
  /*
  iz = 0;
  for (it = 1; it <= parms.nts; it++) {
    update_serial(NXPROB, NYPROB, &u[iz][0][0], &u[1-iz][0][0]);
    iz = 1 - iz;
  }

  prtdat(NXPROB, NYPROB, &u[iz][0][0], "final.dat");
  */
  
  printf("serial: %f\n", serial(NXPROB, NYPROB, &u[0][0][0]));
  printf("pthread: %f\n", ptread_heat2D(NXPROB, NYPROB, &u[0][0][0]));
}

struct PROBLEM {
      int nx, ny;
      float *u1,*u2;
      int threadnum;
} problem;
int totalThread;
pthread_cond_t   cond;
pthread_mutex_t  mtx;

double ptread_heat2D(int nx, int ny, float *u) {
     pthread_t threads[Max_Thread_Num];
     int i, thread_num=sysconf(_SC_NPROCESSORS_ONLN);
     struct timespec ts,te;
     double time_cost;
     float  *temp;
     
     clock_gettime(CLOCK_REALTIME, &ts);
     
     pthread_cond_init(&cond, NULL);
     pthread_mutex_init(&mtx, NULL);
     if ( thread_num>Max_Thread_Num ) thread_num = Max_Thread_Num;
     problem.nx = nx;
     problem.ny = ny;
     problem.u1 = u;
     problem.u2 = u + nx*ny;
     problem.threadnum = 0;
     totalThread = 0;
     for(i=0; i<thread_num; i++) pthread_create(&(threads[i]), NULL, worker, &thread_num);     
       
     for(i=0; i<thread_num; i++) pthread_join(threads[i], NULL);
     pthread_cond_destroy(&cond);
     pthread_mutex_destroy(&mtx);
     
     clock_gettime(CLOCK_REALTIME, &te);
     time_cost = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;

     //prtdat(NXPROB, NYPROB, problem.u1, "final.dat");
     
     return time_cost;

}

void *worker(void *arg) {
     int    ix, iy, ny, it;
     int    myID, thread_num, lbound, ubound, loc_size; 
     float *u1, *u2, *temp;
     
     thread_num = (int)(*(int*)arg);
     
     pthread_mutex_lock(&mtx);
     myID = totalThread;
     totalThread++;
     pthread_mutex_unlock(&mtx);     
     
     ny = problem.ny;
     loc_size = ( problem.nx + thread_num - 3 ) / thread_num;
     lbound = 1 + loc_size * myID;
     ubound = lbound + loc_size - 1;
     if (myID==thread_num - 1) ubound = problem.nx-2;
     
     for (it = 1; it <= parms.nts; it++) {
         u1 = problem.u1;
         u2 = problem.u2;
         for (ix = lbound; ix <= ubound; ix++) {
             for (iy = 1; iy <= problem.ny-2; iy++) {
                *(u2+ix*ny+iy) = *(u1+ix*ny+iy)  + 
	 	       parms.cx * (*(u1+(ix+1)*ny+iy) + *(u1+(ix-1)*ny+iy) - 
	 			   2.0 * *(u1+ix*ny+iy)                      ) +
	 	       parms.cy * (*(u1+ix*ny+iy+1) + *(u1+ix*ny+iy-1) - 
	 			   2.0 * *(u1+ix*ny+iy)                  );
            }
         }
         
         pthread_mutex_lock(&mtx);
         problem.threadnum++;
         if ( problem.threadnum<thread_num ) pthread_cond_wait(&cond, &mtx);
         else {
           temp = problem.u1;
           problem.u1 = problem.u2;
           problem.u2 = temp;
           problem.threadnum = 0;
           pthread_cond_broadcast(&cond);
         }
         pthread_mutex_unlock(&mtx);
     }                   
     
    return (void*)0;
}

double serial(int nx, int ny, float *u){
     int it, iz;
     struct timespec ts,te;
     double time_cost;
     float  *u1, *u2, *temp;
     
     clock_gettime(CLOCK_REALTIME, &ts);
     u1 = u; 
     u2 = u + nx*ny;
     for (it = 1; it <= parms.nts; it++) {
       update_serial(NXPROB, NYPROB, u1, u2);
       temp = u1;
       u1 = u2;
       u2 = temp;
     }
     clock_gettime(CLOCK_REALTIME, &te);
     time_cost = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;

     //prtdat(NXPROB, NYPROB, u1, "final.dat");
     
     return time_cost;
}
/****************************************************************************
 *  subroutine update
 ****************************************************************************/
void update_serial(int nx, int ny, float *u1, float *u2) {
  int ix, iy;

  for (ix = 1; ix <= nx-2; ix++)
  {
    for (iy = 1; iy <= ny-2; iy++)
    {
      *(u2+ix*ny+iy) = *(u1+ix*ny+iy)  + 
		       parms.cx * (*(u1+(ix+1)*ny+iy) + *(u1+(ix-1)*ny+iy) - 
				   2.0 * *(u1+ix*ny+iy)                      ) +
		       parms.cy * (*(u1+ix*ny+iy+1) + *(u1+ix*ny+iy-1) - 
				   2.0 * *(u1+ix*ny+iy)                  );
    }
  }
}

/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny, float *u1) {
  int ix, iy;

  for (ix = 0; ix <= nx-1; ix++)
  {
    for (iy = 0; iy <= ny-1; iy++)
    {
      /* u1[ix][iy] = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1)); */
      *(u1+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
    }
  }
}

/**************************************************************************
 * subroutine prtdat
 **************************************************************************/
void prtdat(int nx, int ny, float *u1, char fnam[]) {
  int ix, iy;
  FILE *fp;

  fp = fopen(fnam, "w");
  for (iy = ny-1; iy >= 0; iy--)
  {
    for (ix = 0; ix <= nx-1; ix++)
    {
      fprintf(fp, "%8.3f", *(u1+ix*ny+iy));
      if (ix != nx-1)
      {
        fprintf(fp, " ");
      }
      else
      {
        fprintf(fp, "\n");
      }
    }
  }
  fclose(fp);
  printf("Wrote file: %s\n",fnam);
}
