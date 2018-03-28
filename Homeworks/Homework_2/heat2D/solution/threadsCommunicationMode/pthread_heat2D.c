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


#define NXPROB 60                       // grid num in x direction
#define NYPROB 60                       // grid num in y direction
#define NANO   1000000              

//#define Max_Thread_Num  6             // maximun threads used in this program
#define BLOCKS_IN_X 3                   // Domain blocks in x direction
#define BLOCKS_IN_Y 2                   // Domain blocks in y diriction
#define Max_Thread_Num  BLOCKS_IN_X*BLOCKS_IN_Y  // maximun threads used in this program


// Sturcture for mesh parameters
struct Parms
{ 
  double cx;                             // $cx = \frac{1}{\Delta x^2}$
  double cy;                             // $cy = \frac{1}{\Delta y^2}$
  long int nts;
} parms = {0.1, 0.1, 500000};



// worker's status vars structure
struct WorkerStatus 
{
       pthread_t  id;                 // Pthread_t id
                          
       long int tid;                  // thread id
       long int nx_id;                // thread id in x direction for computing domain
       long int ny_id;                // thread id in y direction for computing domain
       
       long int nx_local;             // local(thread or block) index boundary in x direction
       long int ny_local;             // local(thread or block) index boundary in y direction
       
       long int neighbor_up_id;       // up neighbor's thread id
       long int neighbor_bottom_id;   // bottom neighbor's thread id 
       long int neighbor_left_id;     // left neighbor's thread id
       long int neighbor_right_id;    // right neighbor's thread id
       
       double *up_data;               // up neighbor's boundary data
       double *bottom_data;           // bottom neighbor's boundary data
       double *left_data;             // left neighbor's boundary data
       double *right_data;            // right beighbor's boundary data
   
       double *u1;                    // thread's temperature field vars
       double *u2;                    // thread's temperature field vars
           
}  threads[Max_Thread_Num];




// Pthread global varibles
int totalThread,threadCounter;                       // Thread counter
pthread_cond_t   cMaster, cWorker, cBarried;         // Pthread Condition Varibles
pthread_mutex_t  mtx, mtxB;                          // Pthread Mutex Varibles






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




// Mutex Worker Threads
//--------------------------------------------------------------------------
void *mutex_worker(void *mydata)
{
  long int nx, ny;
  double *u1, *u2;                                                    /*double u1[nx][ny], u2[nx][ny];*/
  long int ix, iy;
  long int id_up,id_bottom,id_right,id_left; 
  long int t_step;                                                    // Time step
  int thread_num;                                                     // total worker thread number                         
  struct WorkerStatus  *pMyStatus = threads;                          // this form is in C style
  
  thread_num = (int)(mydata);                                        // type conversion   
  while ( pMyStatus->id!=pthread_self() ) pMyStatus++;                // jump to my thread status sturcture 
  
  printf("Thread tid: %ld created successful.\n",pMyStatus->tid);
  
  nx = pMyStatus->nx_local;
  ny = pMyStatus->ny_local;

  u1 = pMyStatus->u1;
  u2 = pMyStatus->u2;
  
  pthread_mutex_lock(&mtx);     
  totalThread++;                         
  pthread_cond_broadcast(&cWorker);        // let master thread wait for worker threads  
  pthread_cond_wait(&cMaster, &mtx);
  pthread_mutex_unlock(&mtx);  
  
    

  t_step = 1;
  while(t_step < parms.nts)
  {    

    // part I of super computing step: computing
    //----------------------------------------------------------------------------
    
    // computing boundary area           
      // four coner points 
         if(pMyStatus->neighbor_left_id != -1 && pMyStatus->neighbor_bottom_id != -1)       // left bottom
         {
           ix = 0; iy = 0;
           *(u2+iy*nx+ix) = *(u1+iy*nx+ix)  + 
	    	             parms.cx * (*(u1+iy*nx+(ix+1)) + *(pMyStatus->left_data+iy) - 
	 		       	         2.0 * *(u1+iy*nx+ix)                      ) +
	 	             parms.cy * (*(u1+iy*nx+iy+1) + *(pMyStatus->bottom_data+ix) - 
	 			         2.0 * *(u1+iy*nx+ix)                  );
         }
          
         if(pMyStatus->neighbor_right_id != -1 && pMyStatus->neighbor_bottom_id != -1)      // right bottom
         {
           ix = nx-1; iy = 0;
           *(u2+iy*nx+ix) = *(u1+iy*nx+ix)  + 
	    	             parms.cx * (*(pMyStatus->right_data+iy) + *(u1+(ix-1)+iy*nx) - 
			      	         2.0 * *(u1+ix+iy*nx)                      ) +
		             parms.cy * (*(u1+ix+(iy+1)*nx) + *(pMyStatus->bottom_data+ix) - 
				         2.0 * *(u1+ix+iy*nx)                  );
         }

         if(pMyStatus->neighbor_left_id != -1 && pMyStatus->neighbor_up_id != -1)           // left up
         {
           ix = 0; iy = ny-1;
           *(u2+ix+iy*nx) = *(u1+ix+iy*nx)  + 
	    	             parms.cx * (*(u1+(ix+1)+iy*nx) + *(pMyStatus->left_data+iy) - 
			   	          2.0 * *(u1+ix+iy*nx)                      ) +
		             parms.cy * (*(pMyStatus->up_data+ix) + *(u1+ix+(iy-1)*nx) - 
				          2.0 * *(u1+ix+iy*nx)                  );
         }

         if(pMyStatus->neighbor_right_id != -1 && pMyStatus->neighbor_up_id != -1)          // right up
         {
           ix = nx-1; iy = ny-1;
           *(u2+ix+iy*nx) = *(u1+ix+iy*nx)  + 
	    	             parms.cx * (*(pMyStatus->right_data+iy) + *(u1+(ix-1)+iy*nx) - 
			   	          2.0 * *(u1+ix+iy*nx)                      ) +
		             parms.cy * (*(pMyStatus->up_data+ix) + *(u1+ix+(iy-1)*nx) - 
				          2.0 * *(u1+ix+iy*nx)                  );
         }

       //----------------------------------------------------------------------------------------------
       // four edge points
       if(pMyStatus->neighbor_bottom_id != -1)     // bottom edge
       {
          iy = 0; 
          for (ix = 1; ix < nx-1; ix++)
            {
               *(u2+iy*nx+ix) = *(u1+iy*nx+ix)  + 
	    	              parms.cx * (*(u1+iy*nx+(ix+1)) + *(u1+iy*nx+(ix-1)) - 
	  	  	       	          2.0 * *(u1+iy*nx+ix)                      ) +
	  	              parms.cy * (*(u1+iy*nx+iy+1) + *(pMyStatus->bottom_data+ix) - 
	  			          2.0 * *(u1+iy*nx+ix)                  );
            }            
       }
       
       if(pMyStatus->neighbor_up_id != -1)         // up edge
       {
           iy = ny-1;
           for (ix = 1; ix < nx-1; ix++)
             {
                  *(u2+ix+iy*nx) = *(u1+ix+iy*nx)  + 
	    	             parms.cx * (*(u1+(ix+1)+iy*nx) + *(u1+(ix-1)+iy*nx) - 
			   	          2.0 * *(u1+ix+iy*nx)                      ) +
		             parms.cy * (*(pMyStatus->up_data+ix) + *(u1+ix+(iy-1)*nx) - 
				          2.0 * *(u1+ix+iy*nx)                  );
             }
       }
       
       if(pMyStatus->neighbor_left_id != -1)       // left edge
       {
          ix = 0;
          for (iy = 1; iy < ny-1; iy++)
            {     
                  *(u2+ix+iy*nx) = *(u1+ix+iy*nx)  + 
	    	             parms.cx * (*(u1+(ix+1)+iy*nx) + *(pMyStatus->left_data+iy) - 
			   	          2.0 * *(u1+ix+iy*nx)                      ) +
		             parms.cy * (*(u1+ix+(iy+1)*nx) + *(u1+ix+(iy-1)*nx) - 
				          2.0 * *(u1+ix+iy*nx)                  );
            }
       }
       
       if(pMyStatus->neighbor_right_id != -1)      // right edge
       {
           ix = nx-1;
           for (iy = 1; iy < ny-1; iy++)
              {
                  *(u2+ix+iy*nx) = *(u1+ix+iy*nx)  + 
	    	             parms.cx * (*(pMyStatus->right_data+iy) + *(u1+(ix-1)+iy*nx) - 
			   	          2.0 * *(u1+ix+iy*nx)                      ) +
		             parms.cy * (*(u1+ix+(iy+1)*nx) + *(u1+ix+(iy-1)*nx) - 
				          2.0 * *(u1+ix+iy*nx)                  );
              }            
       } 
  
    //-----------------------------------------------------------------------------------------------------   
    // computing inner area 
    for (ix = 1; ix < nx -1; ix++)
    {
      for (iy = 1; iy < nx -1; iy++)
      {
         *(u2+iy*nx+ix) = *(u1+iy*nx+ix)  + 
	    	          parms.cx * (*(u1+(ix+1)+iy*nx) + *(u1+(ix-1)+iy*nx) - 
			   	      2.0 * *(u1+ix+iy*nx)                      ) +
		          parms.cy * (*(u1+ix+(iy+1)*nx) + *(u1+ix+(iy-1)*nx) - 
				      2.0 * *(u1+ix+iy*nx)                  );
      }
    }

    //------------------------------------------------------------------------------------------------------
    // Barried
    pthread_mutex_lock(&mtxB);
    threadCounter++ ;
    printf("thread id: %ld, timesteps: %ld, threadCounter: %d\n",pMyStatus->tid,t_step,threadCounter);
    if(threadCounter == thread_num)                                       // thread counter ++
    {
      threadCounter = 0;
      pthread_cond_broadcast(&cBarried);
    } 
    else
    {
      while(pthread_cond_wait(&cBarried, &mtxB) != 0);
    }
    pthread_mutex_unlock(&mtxB);

    // part II of super computing step: communication
    //----------------------------------------------------------------------------
    id_up = pMyStatus->neighbor_up_id;
    id_bottom = pMyStatus->neighbor_bottom_id;

    // up boundary
    if(id_bottom != -1)
        for (ix = 0; ix < nx; ix++)      
            *(pMyStatus->bottom_data + ix) = *(threads[id_bottom].u2 + ((threads[id_bottom].ny_local - 1)*threads[id_bottom].nx_local) + ix);
    // bottom boundary
    if(id_up != -1)
        for (ix = 0; ix < nx; ix++)  
            *(pMyStatus->up_data + ix) = *(threads[id_up].u2 + (0*threads[id_up].nx_local) + ix);

    
    id_right = pMyStatus->neighbor_right_id;
    id_left  = pMyStatus->neighbor_left_id; 
    // left boundary
    if (id_right != -1)
        for (iy = 0; iy < ny; iy++)  
            *(pMyStatus->right_data + iy) = *(threads[id_right].u2 + (iy*threads[id_right].nx_local) + 0);
    // right boundary      
    if (id_left != -1)
        for (iy = 0; iy < ny; iy++)         
            *(pMyStatus->left_data + iy) = *(threads[id_left].u2 + (iy*threads[id_left].nx_local) + (threads[id_left].nx_local -1)); 
    
    //printf("thread id = %ld, time step = %ld boundary data communication done...\n", pMyStatus->tid, t_step);        

    // Barried
    pthread_mutex_lock(&mtxB);
    threadCounter++ ;
    printf("thread id: %ld, timesteps: %ld, threadCounter: %d\n",pMyStatus->tid,t_step,threadCounter);
    if(threadCounter == thread_num)                                       // thread counter ++
    {
      threadCounter = 0;
      pthread_cond_broadcast(&cBarried);
    } 
    else
    {
      while(pthread_cond_wait(&cBarried, &mtxB) != 0);
    }
    pthread_mutex_unlock(&mtxB);
    
    // next time step
    memcpy(u1,u2,(nx*ny) * sizeof(double));
    t_step++;

  }

  pthread_mutex_lock(&mtx);     
  totalThread++;                          
  pthread_cond_broadcast(&cWorker);        // let master thread wait for worker threads
  pthread_cond_wait(&cMaster, &mtx);
  pthread_mutex_unlock(&mtx);

  // wait for master thread read data from worker threads
  pthread_mutex_lock(&mtx);
  pthread_cond_wait(&cMaster, &mtx);
  pthread_mutex_unlock(&mtx);         

  // free block filed vars' memory
  free(u1);
  free(u2);
  free(pMyStatus->bottom_data);
  free(pMyStatus->up_data);
  free(pMyStatus->right_data);
  free(pMyStatus->left_data);
  
  return (void*)0;
}
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------






// Nutex Master Thread
//--------------------------------------------------------------------------
void parallel_heat2D()
{
  int thread_num = sysconf(_SC_NPROCESSORS_ONLN);          // return how many threads can use in this mechine 
  long int nBlockX, nBlockY;                               // number of blocks in computing domain  
  
  double u[2][NXPROB][NYPROB];                             // Temperature
  long int ix, iy, iz, it;                                 // index vars
  long int i,j,k,m,n;                                      // index vars
  long x_temp,y_temp;

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
  for (iy = 0; iy < BLOCKS_IN_Y; iy++) 
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
            threads[i].nx_local = (long int) NXPROB/BLOCKS_IN_X;
            threads[i].ny_local = (long int) NYPROB/BLOCKS_IN_Y;
          }
        else if(ix < BLOCKS_IN_X - 1)                                 // boundary blocks in Y direction
          {
            threads[i].nx_local = (long int) NXPROB/BLOCKS_IN_X;
            threads[i].ny_local = NYPROB - (BLOCKS_IN_Y-1)*((long int) NYPROB/BLOCKS_IN_Y);   
          }
        else if(iy < BLOCKS_IN_Y - 1)                                 // boundary blocks in X direction
          {
            threads[i].nx_local = NXPROB - (BLOCKS_IN_X-1)*((long int) NXPROB/BLOCKS_IN_X);
            threads[i].ny_local = (long int) NYPROB/BLOCKS_IN_Y;
          }
        else                                                          // boundary blocks both in X and Y direction
          {
            threads[i].nx_local = NXPROB - (BLOCKS_IN_X-1)*((long int) NXPROB/BLOCKS_IN_X);
            threads[i].ny_local = NYPROB - (BLOCKS_IN_Y-1)*((long int) NYPROB/BLOCKS_IN_Y);   
          }
        
        printf("Thread id: %ld, nx_local = %ld, ny_local = %ld \n",i,threads[i].nx_local,threads[i].ny_local); 
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
        // malloc memory for boundary data
        threads[i].up_data      = (double*)malloc( threads[i].nx_local * sizeof(double) );
        threads[i].bottom_data  = (double*)malloc( threads[i].nx_local * sizeof(double) );
        threads[i].right_data   = (double*)malloc( threads[i].ny_local * sizeof(double) );
        threads[i].left_data    = (double*)malloc( threads[i].ny_local * sizeof(double) );
        printf("memory allocation for boundary data finished...\n"); 
        // malloc memory for field vars
        threads[i].u1 = (double*)malloc( (threads[i].nx_local*threads[i].ny_local) * sizeof(double) );
        threads[i].u2 = (double*)malloc( (threads[i].nx_local*threads[i].ny_local) * sizeof(double) );
        printf("memory allocation for thread field vars finished...\n\n");
        
      }


  // Initial data for each blocks
  x_temp = 0;      // block initial mesh node shift in x coordinate
  y_temp = 0;      // block initial mesh node shift in y coordinate
  for (iy = 0; iy < BLOCKS_IN_Y; iy++)
    { 
      x_temp = 0;  
      for (ix = 0; ix < BLOCKS_IN_X; ix++)
      {  
        // Thread index
        i = ix + iy*BLOCKS_IN_X;               
        // initial data for field vars
        for (j = 0; j< threads[i].ny_local; j++)
          for (k = 0; k< threads[i].nx_local; k++)
             {
                *(threads[i].u1 + j*threads[i].nx_local + k) = u[0][x_temp+k][y_temp+j]; 
                *(threads[i].u2 + j*threads[i].nx_local + k) = u[1][x_temp+k][y_temp+j];
             }
        x_temp += threads[i].nx_local;             // global-local coordinate transformation         
       }
       y_temp += threads[i].ny_local;              // global-local coordinate transformation
     }

  
  // initial data for boundary                     // (core dump happened here error!)
  for (iy = 0; iy < BLOCKS_IN_Y; iy++) 
    for (ix = 0; ix < BLOCKS_IN_X; ix++)
     { 
        // Thread index
        i = ix + iy*BLOCKS_IN_X;
         
        m = threads[i].neighbor_up_id;
        if (m != -1) 
          for (k = 0; k < threads[i].nx_local; k++) 
               threads[i].up_data[k]       = *(threads[m].u1 + 0*threads[m].nx_local + k);

        m = threads[i].neighbor_bottom_id;
        if (m != -1)
          for (k = 0; k < threads[i].nx_local; k++)       
               threads[i].bottom_data[k]   = *(threads[m].u1 + (threads[m].ny_local - 1)*threads[m].nx_local + k);

        m = threads[i].neighbor_right_id;
        if (m != -1)
          for (j = 0; j < threads[i].ny_local; j++)     
               threads[i].left_data[j]     = *(threads[m].u1 + j*threads[m].nx_local + 0);
      
        m = threads[i].neighbor_left_id;
        if (m != -1)
          for (j = 0; j < threads[i].ny_local; j++)     
               threads[i].right_data[j]    = *(threads[m].u1 + j*threads[m].nx_local + (threads[i].nx_local -1));
     }


  // Pthread mutex and condition vars initial
  pthread_mutex_init(&mtx, NULL);
  pthread_mutex_init(&mtxB, NULL);
  pthread_cond_init(&cWorker, NULL);
  pthread_cond_init(&cMaster, NULL);
  pthread_cond_init(&cBarried, NULL);

  totalThread = 0;     // initial threads counter
  threadCounter = 0;   // initial threads counter 

  // threads create
  for(i=0; i<thread_num; i++) 
      pthread_create(&(threads[i].id), NULL, mutex_worker, (void *)thread_num);
  
  // block the master thread and wait condition send form workers
  pthread_mutex_lock(&mtx);
  while ( totalThread!=thread_num ) pthread_cond_wait(&cWorker, &mtx);
  totalThread = 0;
  pthread_cond_broadcast(&cMaster);     
  pthread_mutex_unlock(&mtx);
  printf("Master thread: All worker threads are ready for computing step....\n");
    
  // Parallel computing 
  /***********************************************************************
  **  Iterate over all timesteps.
  ************************************************************************/
  // block the master thread and wait condition send form workers
  pthread_mutex_lock(&mtx);
  while ( totalThread!=thread_num ) pthread_cond_wait(&cWorker, &mtx);    
  pthread_cond_broadcast(&cMaster);
  pthread_mutex_unlock(&mtx);
  printf("Master thread: All worker threads finished computing step....\n");



  // Assemble the final whole temperature field
  x_temp = 0;      // block initial mesh node shift in x coordinate
  y_temp = 0;      // block initial mesh node shift in y coordinate
  for (iy = 0; iy < BLOCKS_IN_Y; iy++)
    { 
      x_temp = 0;  
      for (ix = 0; ix < BLOCKS_IN_X; ix++)
      {  
        // Thread index
        i = ix + iy*BLOCKS_IN_X;               
        // initial data for field vars
        for (j = 0; j < threads[i].ny_local; j++)
          for (k = 0; k < threads[i].nx_local; k++)
             {
                u[0][x_temp+k][y_temp+j] = *(threads[i].u1 + j*threads[i].nx_local + k); 
                u[1][x_temp+k][y_temp+j] = *(threads[i].u2 + j*threads[i].nx_local + k);
             }
        x_temp += threads[i].nx_local;         
       }
       y_temp += threads[i].ny_local; 
     }

  // tell worker threads to destory themselves
  pthread_cond_broadcast(&cMaster);
  

  // Pthread mutex and condition vars destory
  pthread_mutex_destroy(&mtx);
  pthread_mutex_destroy(&mtxB);
  pthread_cond_destroy(&cWorker);
  pthread_cond_destroy(&cMaster);
  pthread_cond_destroy(&cBarried);

  /***********************************************************************
  **  Output final result.
  ************************************************************************/

  prtdat(NXPROB, NYPROB, &u[1][0][0], "final_serial.dat");

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
    printf("mtx parallel:  cost = %15.10f \n",  mtx_cost);

    clock_gettime(CLOCK_REALTIME, &ts);
    serial_heat2D();
    clock_gettime(CLOCK_REALTIME, &te);
    serial_cost = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;
    printf("serial:  cost = %15.10f \n",  serial_cost);

    printf("Speed Up = frac{serial_cost}{mtx_cost}: %ld\n",serial_cost/mtx_cost);

    return EXIT_SUCCESS;
}







