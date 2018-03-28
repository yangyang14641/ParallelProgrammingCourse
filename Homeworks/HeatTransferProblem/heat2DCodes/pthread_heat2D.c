/****************************************************************************
 * FILE: parallel_heat2D.c
 * DESCRIPTION:  
 *   Serial HEAT2D Example - C Version
 *   This example is based on a simplified 
 *   two-dimensional heat equation domain decomposition.  The initial 
 *   temperature is computed to be high in the middle of the domain and 
 *   zero at the boundaries.  The boundaries are held at zero throughout 
 *   the simulation.  During the time-stepping, an array containing two 
 *   domains is used; these domains alternate between old data and new data.
 ****************************************************************************/
/*
Parallel code was developed by Yang Yang at Peking University in Nov 1 2016.

PDE: $\nabla^{2} u = 0$
Numerical Scheme: five point centered differential scheme

                              (i,j+1)
                                 |
                                 |
                     (i-1,j)---(i,j)---(i,j+1)
                                 |
                                 |
                              (i,j-1)

Solution method jacobi iteration
PDE: $\frac{\partial u}{\partial t} = \alpha \nabla^2 u$
Numerical Scheme: $\frac{u_{i,j}^{n+1} - u_{i,j}^{n}}{\Delta t} = \alpha (\frac{u_{i+1,j}^{n} + u_{i-1,j}^{n} 
- 2u_{i,j}^{n}}{\Delta x^2} + \frac{u_{i,j+1}^{n} + u_{i,j-1}^{n} - u_{i,j}^{n}}{\Delta y^2})$
Stability condition: The Fourier number \frac{\alpha \Delta t}{\Delta x^2} less than \frac{1}{4}
*/

/*
Blocks(Threads) division:

                 Y     -------------------
                 ^     |     |     |     |
                 |     | [6] | [7] | [8] |
                 |     |     |     |     |
                 |     -------------------
                 |     |     |     |     |
                 |     | [3] | [4] | [5] |
                 |     |     |     |     |
                 |     -------------------
                 |     |     |     |     |
                 |     | [0] | [1] | [2] |
                 |     |     |     |     |
                 |     -------------------
                 |----------------------------> X

This pthread program using threads Master-Slaves sturcture. Master thread is manager while other threads are computing workers.
*/

//-----------------------------------------------------------------------------------------------------------------

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h> 
//#include <math.h>
#include <stdlib.h>
#include <pthread.h>


#define NXPROB 120                       // grid num in x direction
#define NYPROB 120                       // grid num in y direction
#define NANO   1000000000              

//#define Max_Thread_Num  6              // maximun threads used in this program
#define BLOCKS_IN_X 3                    // Domain blocks in x direction
#define BLOCKS_IN_Y 2                    // Domain blocks in y diriction
#define Max_Thread_Num  BLOCKS_IN_X*BLOCKS_IN_Y  // maximun threads used in this program


// Sturcture for mesh parameters
struct Parms
{ 
  double cx;                             // $cx = \frac{1}{\Delta x^2}$
  double cy;                             // $cy = \frac{1}{\Delta y^2}$
  long int nts;
} parms = {0.1, 0.1, 5000000};




/****************************************************************************
 *  subroutine update
 ****************************************************************************/
void update(long int nx, long int ny, double *u1, double *u2)   /*double u1[nx][ny], u2[nx][ny];*/
{
  long int ix, iy;

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
void inidat(long int nx, long int ny, double *u1)     /*double u1[nx][ny];*/
{
  long int ix, iy;

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
void prtdat(long int nx, long int ny, double *u1, char *fnam)  /*double u1[nx][ny];*/
{
  long int ix, iy;
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



//-------------------------------------------------------------------------
// Serial heat 2D
//-------------------------------------------------------------------------

void serial_heat2D()
{ 
  double u[2][NXPROB][NYPROB];
  long int ix, iy, iz, it;
  double field_sum;
//  void inidat(), prtdat(), update();

  /************************************************************************
  **  Initialize grid.
  *************************************************************************/
  inidat(NXPROB, NYPROB, &u[0][0][0]);
  prtdat(NXPROB, NYPROB, &u[0][0][0], "initial_serial.dat");


  for (ix = 0; ix <= NXPROB-1; ix++)
  {
    u[1][ix][0] = u[0][ix][0];
    u[1][ix][NYPROB-1] = u[0][ix][NYPROB-1];
  }
  for (iy = 0; iy <= NYPROB-1; iy++)
  {
    u[1][0][iy] = u[0][0][iy];
    u[1][NXPROB-1][iy] = u[0][NXPROB-1][iy];
  }

  /***********************************************************************
  **  Iterate over all timesteps.
  ************************************************************************/
  iz = 0;
  for (it = 1; it <= parms.nts; it++)
  {
    update(NXPROB, NYPROB, &u[iz][0][0], &u[1-iz][0][0]);
    iz = 1 - iz;
  }


  /***********************************************************************
  **  Output final result.
  ************************************************************************/

  prtdat(NXPROB, NYPROB, &u[iz][0][0], "final_serial.dat");


  field_sum = 0.;
  for(ix=0; ix<NXPROB; ix++)
      for(iy=0; iy<NYPROB; iy++)
           field_sum += u[0][ix][iy];
   printf("Sum of whole field is %f\n",field_sum);

}












//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//     Parallel heat 2D 
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

// worker's status vars structure   (local)
struct WorkerStatus 
{
       pthread_t  id;                 // Pthread_t id
                          
       long int tid;                  // thread id
       long int nx_id;                // thread id in x direction for computing domain
       long int ny_id;                // thread id in y direction for computing domain
       
       long int nx_local;             // local(thread or block) index boundary in x direction
       long int ny_local;             // local(thread or block) index boundary in y direction
       
       long int offset_x;             // coordinate offset in x direction        
       long int offset_y;             // coordinate offset in y direction

       long int neighbor_up_id;       // up neighbor's thread id
       long int neighbor_bottom_id;   // bottom neighbor's thread id 
       long int neighbor_left_id;     // left neighbor's thread id
       long int neighbor_right_id;    // right neighbor's thread id
           
}  threads[Max_Thread_Num];


// problem's structure (global)
struct PROBLEM 
{
      long int nx, ny;
      double *u1, *u2;
} problem;


// Pthread global varibles
int totalThread;                                     // Thread counter
pthread_cond_t   cBarried;                           // Pthread Condition Varibles
pthread_mutex_t  mtx;                                // Pthread Mutex Varibles




// Mutex Worker Threads
//--------------------------------------------------------------------------
void *mutex_worker(void *mydata)
{
  double *u1, *u2, *temp;                                             // double u1[nx][ny], u2[nx][ny];
  long int nx_local, ny_local;                                        // local nx and local ny
  long int nx, ny;                                                    // problem's nx and ny
  long int ix_local, iy_local;                                        // local index
  long int ix, iy;                                                    // problem's index
  long int offset_x, offset_y;                                        // offset in x direction and y direction
  long int t_step;                                                    // Time step
  int thread_num;                                                     // total worker thread number                         
  struct WorkerStatus  *pMyStatus = threads;                          // this form is in C style
  
  thread_num = (int)(mydata);                                         // type conversion   
  while ( pMyStatus->id!=pthread_self() ) pMyStatus++;                // jump to my thread status sturcture 
  
  printf("Thread tid: %ld created successful.\n",pMyStatus->tid);
  
  nx_local = pMyStatus->nx_local;
  ny_local = pMyStatus->ny_local;
  offset_x = pMyStatus->offset_x;
  offset_y = pMyStatus->offset_y;

  u1 = problem.u1;
  u2 = problem.u2;
  nx = problem.nx;
  ny = problem.ny;
 
  // Barried and Sychronization of worker threads
  pthread_mutex_lock(&mtx);
  totalThread++ ;
  if(totalThread == thread_num)                                       // thread counter ++
  {
    totalThread = 0;
    pthread_cond_broadcast(&cBarried);
  } 
  else
  {
    while(pthread_cond_wait(&cBarried, &mtx) != 0);
  }
  pthread_mutex_unlock(&mtx);

  // time resolution step
  printf("thread id: %ld, stating computing step...\n",pMyStatus->tid);  
  t_step = 1;
  while(t_step < parms.nts)
  {    

    // part I of super computing step: computing
    //----------------------------------------------------------------------------
  
    //-----------------------------------------------------------------------------------------------------   
    // computing inner area 
    for (ix_local = 0; ix_local < nx_local; ix_local++)
    {
      ix = offset_x + ix_local;                     // local position to global position (Cartesian)
      for (iy_local = 0; iy_local < ny_local; iy_local++)
      {  
         // local global coordinate transformation
         iy = offset_y + iy_local;                  // local position to global position (Cartesian)
         // *u2 is write only whie *u1 is read only.
         *(u2+ix*ny+iy) = *(u1+ix*ny+iy)  + 
		           parms.cx * (*(u1+(ix+1)*ny+iy) + *(u1+(ix-1)*ny+iy) - 
				       2.0 * *(u1+ix*ny+iy)                      ) +
		           parms.cy * (*(u1+ix*ny+iy+1) + *(u1+ix*ny+iy-1) - 
				       2.0 * *(u1+ix*ny+iy)                  );
      }
      // local global coordinate transformation
    }
    //printf("thread id: %ld, Here\n", pMyStatus->tid);
    //------------------------------------------------------------------------------------------------------
    // Barried and Sychronization of worker threads
    pthread_mutex_lock(&mtx);
    totalThread++ ;
    if(totalThread == thread_num)                                       // thread counter ++
    {
      totalThread = 0;
      // next time step
      temp = u1;
      u1 = u2;
      u2 = temp;
      pthread_cond_broadcast(&cBarried);
    } 
    else
    {
      while(pthread_cond_wait(&cBarried, &mtx) != 0);
    }
    pthread_mutex_unlock(&mtx);
    t_step++;
    //if(pMyStatus->tid == 1)
    //printf("thread id: %ld, timesteps: %ld, totalThread: %d\n",pMyStatus->tid,t_step,totalThread);

  }
  printf("thread id: %ld, finished computing step...\n",pMyStatus->tid);  


  return (void*)0;
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------






// Mutex parallel Master Thread
//--------------------------------------------------------------------------
void parallel_heat2D()
{
  int thread_num = sysconf(_SC_NPROCESSORS_ONLN);          // return how many threads can use in this mechine 
  long int nBlockX, nBlockY;                               // number of blocks in computing domain  
  long int offset_x, offset_y;                             // offset in x direction and y direction

  double u[2][NXPROB][NYPROB];                             // Temperature
  long int ix, iy, iz, it;                                 // index vars
  long int i,j,k,m,n;                                      // index vars
  long x_temp,y_temp;

  double field_sum;

  /************************************************************************
  **  Initialize grid.
  *************************************************************************/
  inidat(NXPROB, NYPROB, &u[0][0][0]);                              // field initial 
  prtdat(NXPROB, NYPROB, &u[0][0][0], "initial_parallel.dat");      // output initial field
  printf("\n");
  // Boundary condition
  for (ix = 0; ix <= NXPROB-1; ix++)
  {
    u[1][ix][0] = u[0][ix][0];
    u[1][ix][NYPROB-1] = u[0][ix][NYPROB-1];
  }
  for (iy = 0; iy <= NYPROB-1; iy++)
  {
    u[1][0][iy] = u[0][0][iy];
    u[1][NXPROB-1][iy] = u[0][NXPROB-1][iy];
  }
  

  /************************************************************************
  **  Task division and initial for parallel computing.
  *************************************************************************/
  
  // determine how many threads we use
  if ( thread_num>Max_Thread_Num )                // do as what we expected  
     {  
        thread_num = Max_Thread_Num;              // determine number of threads
        nBlockX = BLOCKS_IN_X;
        nBlockY = BLOCKS_IN_Y;
     }
  else                                           // do as simple
     {
        nBlockX = thread_num;
        nBlockY = 1; 
     } 
      

  // domain blocks define
  offset_x = 1;
  offset_y = 1;
  for (iy = 0; iy < BLOCKS_IN_Y; iy++)
  { 
    for (ix = 0; ix < BLOCKS_IN_X; ix++)
     {  
        // Thread index
        i = ix + iy*BLOCKS_IN_X;                 
        
        // Thread number
        threads[i].tid = i;                      
        
        // Thread(Block) position
        threads[i].nx_id = ix;                   
        threads[i].ny_id = iy;                   
        
        
        // define blocks' size
        if (ix < BLOCKS_IN_X - 1 && iy < BLOCKS_IN_Y - 1)             // Inner blocks (task's amount is same)
          {
            threads[i].nx_local = (long int) (NXPROB-2)/BLOCKS_IN_X;
            threads[i].ny_local = (long int) (NYPROB-2)/BLOCKS_IN_Y;
          }
        else if(ix < BLOCKS_IN_X - 1)                                 // boundary blocks in Y direction
          {
            threads[i].nx_local = (long int) (NXPROB-2)/BLOCKS_IN_X;
            threads[i].ny_local = (NYPROB-2) - (BLOCKS_IN_Y-1)*((long int) (NYPROB-2)/BLOCKS_IN_Y);   
          }
        else if(iy < BLOCKS_IN_Y - 1)                                 // boundary blocks in X direction
          {
            threads[i].nx_local = (NXPROB-2) - (BLOCKS_IN_X-1)*((long int) (NXPROB-2)/BLOCKS_IN_X);
            threads[i].ny_local = (long int) (NYPROB-2)/BLOCKS_IN_Y;
          }
        else                                                          // boundary blocks both in X and Y direction
          {
            threads[i].nx_local = (NXPROB-2) - (BLOCKS_IN_X-1)*((long int) (NXPROB-2)/BLOCKS_IN_X);
            threads[i].ny_local = (NYPROB-2) - (BLOCKS_IN_Y-1)*((long int) (NYPROB-2)/BLOCKS_IN_Y);   
          }
        
        printf("Thread id: %ld, nx_local = %ld, ny_local = %ld \n",i,threads[i].nx_local, threads[i].ny_local); 
        
        // define offset in x 
        threads[i].offset_x = offset_x;
        threads[i].offset_y = offset_y;
        offset_x += threads[i].nx_local;
        printf("Thread id: %ld, offset_x = %ld, offset_y = %ld \n",i,threads[i].offset_x, threads[i].offset_y);

        // define blocks' neighbor blocks
        if (ix == 0 && iy == 0)                                          // left bottom conner block                                      
          {
            threads[i].neighbor_up_id      = ix + (iy+1)*BLOCKS_IN_X;
            threads[i].neighbor_bottom_id  = -1;                         // -1 means no neighbor in this direction
            threads[i].neighbor_right_id   = (ix+1) + iy*BLOCKS_IN_X;
            threads[i].neighbor_left_id    = -1;
           }
        else if(ix == BLOCKS_IN_X - 1 && iy == 0)                        // right bottom conner block
           {
            threads[i].neighbor_up_id      = ix + (iy+1)*BLOCKS_IN_X;
            threads[i].neighbor_bottom_id  = -1;
            threads[i].neighbor_right_id   = -1;
            threads[i].neighbor_left_id    = (ix-1) + iy*BLOCKS_IN_X;                 
           }
        else if(ix == 0 && iy == BLOCKS_IN_Y - 1)                        // left up conner block
           {
            threads[i].neighbor_up_id      = -1;
            threads[i].neighbor_bottom_id  = ix + (iy-1)*BLOCKS_IN_X;
            threads[i].neighbor_right_id   = (ix+1) + iy*BLOCKS_IN_X;
            threads[i].neighbor_left_id    = -1;                        
           }
        else if(ix == BLOCKS_IN_X - 1 && iy == BLOCKS_IN_Y - 1)          // right up conner block
           {
            threads[i].neighbor_up_id      = -1;
            threads[i].neighbor_bottom_id  = ix + (iy-1)*BLOCKS_IN_X;
            threads[i].neighbor_right_id   = -1;
            threads[i].neighbor_left_id    = (ix-1) + iy*BLOCKS_IN_X;
           }
        else if(iy == 0)                                                 // bottom edge blocks
           {
            threads[i].neighbor_up_id      = ix + (iy+1)*BLOCKS_IN_X;
            threads[i].neighbor_bottom_id  = -1;
            threads[i].neighbor_right_id   = (ix+1) + iy*BLOCKS_IN_X;
            threads[i].neighbor_left_id    = (ix-1) + iy*BLOCKS_IN_X;       
           }
        else if(iy == BLOCKS_IN_Y - 1)                                   // up edge blocks
           {
            threads[i].neighbor_up_id      = -1;
            threads[i].neighbor_bottom_id  = ix + (iy-1)*BLOCKS_IN_X;
            threads[i].neighbor_right_id   = (ix+1) + iy*BLOCKS_IN_X;
            threads[i].neighbor_left_id    = (ix-1) + iy*BLOCKS_IN_X;                 
           }
        else if(ix == 0)                                                 // left edge blocks
           {
            threads[i].neighbor_up_id      = ix + (iy+1)*BLOCKS_IN_X;
            threads[i].neighbor_bottom_id  = ix + (iy-1)*BLOCKS_IN_X;
            threads[i].neighbor_right_id   = (ix+1) + iy*BLOCKS_IN_X;
            threads[i].neighbor_left_id    = -1;
           }
        else if(ix == BLOCKS_IN_X - 1)                                  // right edge blocks
           {
            threads[i].neighbor_up_id      = ix + (iy+1)*BLOCKS_IN_X;
            threads[i].neighbor_bottom_id  = ix + (iy-1)*BLOCKS_IN_X;
            threads[i].neighbor_right_id   = -1;
            threads[i].neighbor_left_id    = (ix-1) + iy*BLOCKS_IN_X;                 
           }
        else                                                            // Inner blocks     
           {
            threads[i].neighbor_up_id      = ix + (iy+1)*BLOCKS_IN_X;
            threads[i].neighbor_bottom_id  = ix + (iy-1)*BLOCKS_IN_X;
            threads[i].neighbor_right_id   = (ix+1) + iy*BLOCKS_IN_X;
            threads[i].neighbor_left_id    = (ix-1) + iy*BLOCKS_IN_X;                   
           }
        printf("Thread id: %ld, up_id = %ld, bottom_id = %ld, ", i, threads[i].neighbor_up_id, threads[i].neighbor_bottom_id);
        printf("left_id = %ld, right_id = %ld\n", threads[i].neighbor_left_id, threads[i].neighbor_right_id);
        printf("\n");
      }
      offset_y += threads[i].ny_local;
      offset_x = 1;
  }


  // Pthread mutex and condition vars initial
  pthread_mutex_init(&mtx, NULL);
  pthread_cond_init(&cBarried, NULL);
  
  // threads counter reset
  totalThread = 0;     // initial threads counter 
  
  // define problem sturcture
  problem.nx = NXPROB;
  problem.ny = NYPROB;
  problem.u1 = &u[0][0][0];
  problem.u2 = &u[1][0][0];


  // threads create (fork)
  for(i=0; i<thread_num; i++) 
      pthread_create(&(threads[i].id), NULL, mutex_worker, (void *)thread_num);
  
  // threads join (join)
  for(i=0; i<thread_num; i++) pthread_join(threads[i].id, NULL);

  // Pthread mutex and condition vars destory
  pthread_mutex_destroy(&mtx);
  pthread_cond_destroy(&cBarried);

  /***********************************************************************
  **  Output final result.
  ************************************************************************/

  prtdat(NXPROB, NYPROB, &u[1][0][0], "final_parallel.dat");



  field_sum = 0.;
  for(ix=0; ix<NXPROB; ix++)
      for(iy=0; iy<NYPROB; iy++)
           field_sum += u[0][ix][iy];
   printf("Sum of whole field is %f\n",field_sum);


}




/**************************************************************************
 * Program main
 **************************************************************************/
int main(int argc, char* argv[] )
{
    struct timespec ts,te;
    double serial_cost,mtx_cost;           //serial_cost, mtx_cost
    


    clock_gettime(CLOCK_REALTIME, &ts);
    parallel_heat2D();
    clock_gettime(CLOCK_REALTIME, &te);
    mtx_cost = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;
    printf("mtx parallel:  cost = %15.10lf \n",  mtx_cost);

    clock_gettime(CLOCK_REALTIME, &ts);
    serial_heat2D();
    clock_gettime(CLOCK_REALTIME, &te);
    serial_cost = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;
    printf("serial:  cost = %15.10lf \n",  serial_cost);

    printf("Speed Up = %lf\n",serial_cost/mtx_cost);

    return EXIT_SUCCESS;
}







