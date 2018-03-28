
/****************************************************************************
 * FILE: mpi_heat2D.c
 * OTHER FILES: draw_heat.c  
 * DESCRIPTIONS:  
 *   HEAT2D Example - Parallelized C Version
 *   This example is based on a simplified two-dimensional heat 
 *   equation domain decomposition.  The initial temperature is computed to be 
 *   high in the middle of the domain and zero at the boundaries.  The 
 *   boundaries are held at zero throughout the simulation.  During the 
 *   time-stepping, an array containing two domains is used; these domains 
 *   alternate between old data and new data.
 *
 *   In this parallelized version, the grid is decomposed by the master
 *   process and then distributed by rows to the worker processes.  At each 
 *   time step, worker processes must exchange border data with neighbors, 
 *   because a grid point's current temperature depends upon it's previous
 *   time step value plus the values of the neighboring grid points.  Upon
 *   completion of all time steps, the worker processes return their results
 *   to the master process.
 *
 *   Two data files are produced: an initial data set and a final data set.
 *   An X graphic of these two states displays after all calculations have
 *   completed.
 * AUTHOR: Blaise Barney - adapted from D. Turner's serial C version. Converted
 *   to MPI: George L. Gusciora (1/95)
 * LAST REVISED: 06/12/13 Blaise Barney
 ****************************************************************************/
/*  This domain split in two dimension version was developed by Yang Yang at Peking University  
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

    LAST REVISED: Novenber 2016
*/


/*   Computing domian split (Process Mapping)

     -----------------------------
     |    i,j     ||    i,j      |
     |   (0,1)    ||   (1,1)     |
     |  Block[2]  ||  Block[3]   |
     |  PID = 3   ||  PID = 4    |
     |            ||             |
     |---------------------------|
     |---------------------------|
     |    i,j     ||    i,j      |
     |   (0,0)    ||   (1,0)     |
     |  Block[0]  ||  Block[1]   |
     |  PID = 1   ||  PID = 2    |
     |            ||             |
     -----------------------------
     index: PID = i + j*BLOCKS_IN_X + 1
     Master PID = 0
*/



#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
// extern void draw_heat(int nx, int ny);       /* X routine to create graph */


#define NXPROB      120                  /* x dimension of problem grid */
#define NYPROB      120                  /* y dimension of problem grid */
#define STEPS       5000000              /* number of time steps */

#define MAXWORKER   8                   /* maximum number of worker tasks in user's machine */
#define MINWORKER   1                   /* minimum number of worker tasks in user's machine */
#define BLOCKS_IN_X    2                /* Process blocks in x direction */
#define BLOCKS_IN_Y    2                /* Process blocks in y direction */

#define BEGIN       1                   /* message tag */
#define LTAG        2                   /* message tag */
#define RTAG        3                   /* message tag */
#define UTAG        4                   /* message tag */
#define BTAG        5                   /* message tag */
#define DONE        4                   /* message tag */

#define NONE        0                   /* indicates no neighbor */
#define MASTER      0                   /* taskid of first process */


struct Parms { 
  double cx;
  double cy;
} parms = {0.1, 0.1};


void inidat(), prtdat(), update();

int main (int argc, char *argv[])
{

double  u[2][NXPROB][NYPROB];                          /* array for grid */
int	taskid,                                            /* this task's unique id */
	numworkers,                                        /* number of worker processes */
	numtasks,                                          /* number of tasks */
	averow,avecol,rows,cols,
    offsetrow,offsetcol,extrarow,extracol,             /* for sending rows and col of data */
	dest, source,                                      /* to - from for message send-receive */
	left,right,up,bottom,                              /* neighbor tasks */
	msgtype,                                           /* for message types */
	rc,rowstart,colstart,rowend,colend,                /* misc */
	i,j,k,m,n,ix,iy,iz,it;                             /* loop variables */
double *tempdatasend,*tempdatarecieve,
       *temprowdatasend,*temprowdatarecieve;            /* temporary row data */    
double field_sum;
double tstart,tend;                                    /* time counter */
MPI_Status status;


/* First, find out my taskid and how many tasks are running */
   MPI_Init(&argc,&argv);
   tstart = MPI_Wtime();
   MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
   MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
   
   numworkers = BLOCKS_IN_X*BLOCKS_IN_Y;
  
   if (taskid == MASTER) 
   {
      printf("\nProbelm Scale: NXPROB = %d, NYPROB = %d, STEPS = %d\n\n",NXPROB,NYPROB,STEPS);
      /************************* master code *******************************/
      /* Check the number of process entered by user is legal or illegal*/
      printf(" Number of Process entered by user is : %d\n", numtasks);
      printf(" Number of Worker Process seted in Program is: %d\n", numworkers);
      if (numworkers != numtasks-1)
      {  
         printf("numtask is not equals to numwokers, Please enter right number of processes. Error!\n");
         printf("Program Exit...\n");
         MPI_Abort(MPI_COMM_WORLD, rc);
         exit(1);
      }
      /* Check if numworkers is within range - quit if not */
      if ((numworkers > MAXWORKER) || (numworkers < MINWORKER)) 
      {
         printf("ERROR: the number of tasks must be between %d and %d.\n",
                 MINWORKER+1,MAXWORKER+1);
         printf("Quitting...\n");
         MPI_Abort(MPI_COMM_WORLD, rc);
         exit(1);
      }
      printf ("Starting mpi_heat2D with %d worker tasks.\n", numworkers);

      /* Initialize grid */
      printf("Grid size: X= %d  Y= %d  Time steps= %d\n",NXPROB,NYPROB,STEPS);
      printf("Grid Domain split numbers in X: %d, and split numbers in Y: %d\n",BLOCKS_IN_X,BLOCKS_IN_Y);
      printf("Initializing grid and writing initial.dat file...\n");
      inidat(NXPROB, NYPROB, u);
      prtdat(NXPROB, NYPROB, u, "initial.dat");

      /* Distribute work to workers.  Must first figure out how many rows and colums */
      /* send and what to do with extra rows.  */
      averow = NXPROB/BLOCKS_IN_X;
      avecol = NYPROB/BLOCKS_IN_Y;
      extrarow = NXPROB%BLOCKS_IN_X;
      extracol = NYPROB%BLOCKS_IN_Y;
      
      offsetrow = 0;
      offsetcol = 0;

      for (j=0; j<BLOCKS_IN_Y; j++)
      {  
         offsetrow = 0;
         for (i=0; i<BLOCKS_IN_X; i++)
         {
            rows = (i < extrarow) ? averow+1 : averow;
            cols = (j < extracol) ? avecol+1 : avecol;  
            /* Tell each worker who its neighbors are, since they must exchange */
            /* data with each other. */  
            if (i == 0) 
                left = NONE;
            else
                left = (i-1)+ j*BLOCKS_IN_X + 1;
            
            if (i == BLOCKS_IN_X -1)
                right = NONE;
            else
                right = (i+1) + j*BLOCKS_IN_X + 1;

            if (j == 0)
                bottom = NONE;
            else
                bottom = i + (j-1)*BLOCKS_IN_X + 1;
            
            if (j == BLOCKS_IN_Y -1)
                up = NONE;
            else
                up = i + (j+1)*BLOCKS_IN_X + 1;

            /*  Now send startup information to each worker  */
            dest = i + j*BLOCKS_IN_X + 1;         // Destination Process
            MPI_Send(&offsetrow, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&offsetcol, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&cols, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&left, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&right, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&up, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
            MPI_Send(&bottom, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);   
            /* Prepare and send computing domain data to each worker*/
            tempdatasend = (double *)malloc(rows*cols*sizeof(double));
            for(m=0; m<rows; m++)
            {
               for(n=0; n<cols; n++)
               {
                    k = m*cols + n; 
                    *(tempdatasend+k) = u[0][offsetrow+m][offsetcol+n];
               }
            }              
            MPI_Send(tempdatasend, rows*cols, MPI_DOUBLE, dest, BEGIN, 
                  MPI_COMM_WORLD);
            free(tempdatasend);      


            printf("Sent to task %d: rows= %d cols= %d offsetrow= %d offsetcol= %d\n",dest,rows,cols,offsetrow,offsetcol);
            printf("left= %d right= %d bottom= %d up= %d\n",left,right,bottom,up);
            offsetrow = offsetrow + rows;
         }
         offsetcol = offsetcol + cols;
      }

      /* Now wait for results from all worker tasks */
      for (j=0; j<BLOCKS_IN_Y; j++)
      { 
         for (i=0; i<BLOCKS_IN_X; i++)
         {
           source = i + j*BLOCKS_IN_X + 1;             /* Source Process */
           msgtype = DONE;
           MPI_Recv(&offsetrow, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, 
                  &status);
           MPI_Recv(&offsetcol, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, 
                  &status);
           MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
           MPI_Recv(&cols, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
           /* Receive and reshape data from each worker process*/
           tempdatarecieve = (double *)malloc(rows*cols*sizeof(double));
           MPI_Recv(tempdatarecieve, rows*cols, MPI_DOUBLE, source, msgtype, 
                  MPI_COMM_WORLD, &status);
           for(m=0; m<rows; m++)
            {
               for(n=0; n<cols; n++)
               {
                    k = m*cols + n; 
                    u[0][offsetrow+m][offsetcol+n] = *(tempdatarecieve+k);
               }
            }          
           free(tempdatarecieve);

         }           
                   
      }

      /* Write final output, call X graph and finalize MPI */
      printf("Writing final.dat file and generating graph...\n");
      prtdat(NXPROB, NYPROB, &u[0][0][0], "final.dat");
      tend = MPI_Wtime();
      printf("Click on MORE button to view initial/final states.\n");
      printf("Click on EXIT button to quit program.\n");
//      draw_heat(NXPROB,NYPROB);
      field_sum = 0.;
      for(i=0;i<NXPROB;i++)
         for(j=0;j<NYPROB;j++)
              field_sum += u[0][i][j];
      printf("Sum of whole field is %f\n",field_sum);
      printf("Total time use is %22.16E second\n",tend-tstart);
      MPI_Finalize();
   }   /* End of master code */



   /************************* workers code **********************************/
   if (taskid != MASTER) 
   {
      /* Initialize everything - including the borders - to zero */
      for (iz = 0; iz < 2; iz++)
         for (ix = 0; ix < NXPROB; ix++) 
            for (iy = 0; iy < NYPROB; iy++) 
                 u[iz][ix][iy] = 0.;

      /* Receive my offset, rows, neighbors and grid partition from master */
      source = MASTER;
      msgtype = BEGIN;
      MPI_Recv(&offsetrow, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&offsetcol, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&cols, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&left, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&right, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&up, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
      MPI_Recv(&bottom, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
      /* Receive and reshape data from each worker process*/
      tempdatarecieve = (double *)malloc(rows*cols*sizeof(double));
      MPI_Recv(tempdatarecieve, rows*cols, MPI_DOUBLE, source, msgtype, 
             MPI_COMM_WORLD, &status);
      for(m=0; m<rows; m++)
       {
          for(n=0; n<cols; n++)
          {
               k = m*cols + n; 
               u[0][offsetrow+m][offsetcol+n] = *(tempdatarecieve+k);
          }
       }          
      free(tempdatarecieve);     
     
      /* Determine border elements.  Need to consider first and last columns. */
      /* Obviously, row 0 can't exchange with row 0-1.  Likewise, the last */
      /* row can't exchange with last+1.  */
      rowstart = offsetrow;
      colstart = offsetcol;
      rowend = offsetrow + rows-1;
      colend = offsetcol + cols-1;
      /* Boundary situation */
      if (offsetrow == 0) 
         rowstart = 1;
      if ((offsetrow+rows) == NXPROB) 
         rowend--;
      if (offsetcol == 0) 
         colstart = 1;
      if ((offsetcol+cols) == NXPROB) 
         colend--;         
      printf("task= %d  rowstart= %d  rowend= %d  colstart= %d  colend= %d\n",taskid,rowstart,rowend,colstart,colend);

      /* Begin doing STEPS iterations.  Must communicate border rows and colums with */
      /* neighbors.  If I have the first or last grid row or colum, then I only need */
      /*  to  communicate with one neighbor in row or colum */
      temprowdatasend = (double*)malloc(rows*sizeof(double));               /* allocate continue memory space for send rowdata */
      temprowdatarecieve = (double *)malloc(rows*sizeof(double));           /* allocate continue memory space for recieve rowdata */
      printf("Task %d received work. Beginning time steps...\n",taskid);
      iz = 0;
      for (it = 1; it <= STEPS; it++)
      {  
         /* Send boundary data by using message passing interface */ 
         //printf("it = %d\n",it);
         if (left != NONE)
         {
            /* Send column data to left neighbor */  
            MPI_Send(&u[iz][offsetrow][offsetcol], cols, MPI_DOUBLE, left,
                     RTAG, MPI_COMM_WORLD);
         }
         if (right != NONE)
         {
            /* Send column data to right neighbor */
            MPI_Send(&u[iz][offsetrow+rows-1][offsetcol], cols, MPI_DOUBLE, right,
                      LTAG, MPI_COMM_WORLD);
         }
         /* Watch out for up and bottom the memory is not continue */
         if (up != NONE)
         {  
            /* Send data to up neighbor */ 
            for(k=0; k<rows; k++)
                *(temprowdatasend+k) =  u[iz][offsetrow+k][offsetcol+cols-1];
            MPI_Send(temprowdatasend, rows, MPI_DOUBLE, up,
                      BTAG, MPI_COMM_WORLD);       
         }
         if (bottom != NONE)
         {
            /* Send data to bottom neighbor */ 
            for(k=0; k<rows; k++)
                *(temprowdatasend+k) =  u[iz][offsetrow+k][offsetcol];
            MPI_Send(temprowdatasend, rows, MPI_DOUBLE, bottom,
                      UTAG, MPI_COMM_WORLD);
         }

         /* Recieve Boundary data by using message passing interface */
         if (left != NONE)
         {  //printf("taskid : %d  I'm alive! left != NONE Before Reciev\n",taskid); 
            /* Recieve column data from left neighbor */
            source = left;
            msgtype = LTAG;           
            MPI_Recv(&u[iz][offsetrow-1][offsetcol], cols, MPI_DOUBLE, source,
                      msgtype, MPI_COMM_WORLD, &status);
            //printf("taskid : %d  I'm alive! left != NONE After Recieve\n",taskid);          
         }
         if (right != NONE)
         {  //printf("taskid : %d  I'm alive! right != NONE Before Recieve\n",taskid); 
            /* Recieve column data from right neighbor */
            source = right;
            msgtype = RTAG;
            MPI_Recv(&u[iz][offsetrow+rows][offsetcol], cols, MPI_DOUBLE, source, msgtype,
                      MPI_COMM_WORLD, &status);
            //printf("taskid : %d  I'm alive! right != NONE After Recieve\n",taskid);          
         }
         /* Watch out for up and bottom the memory is not continue */
         if (up != NONE)
         {  //printf("taskid : %d  I'm alive! up != NONE Before Recieve\n",taskid); 
            /* Receive data from up neoghbor */
            source = up;
            msgtype = UTAG;
            
            MPI_Recv(temprowdatarecieve, rows, MPI_DOUBLE, source, msgtype,
                      MPI_COMM_WORLD, &status);
            for(k=0; k<rows; k++)
                u[iz][offsetrow+k][offsetcol+cols] = *(temprowdatarecieve+k);   
            //printf("taskid : %d  I'm alive! up != NONE After Recieve\n",taskid);       
         }
         if (bottom != NONE)
         {  //printf("taskid : %d  I'm alive! bottom != NONE Before Recieve\n",taskid); 
            /* Recieve data from bottom neighbor */
            source = bottom;
            msgtype = BTAG;
            
            MPI_Recv(temprowdatarecieve, rows, MPI_DOUBLE, source, msgtype,
                      MPI_COMM_WORLD, &status);
            for(k=0; k<rows; k++)
                u[iz][offsetrow+k][offsetcol-1] = *(temprowdatarecieve+k); 
            //printf("taskid : %d  I'm alive! bottom != NONE After Recieve\n",taskid);    
         }

         /* Now call update to update the value of grid points */
         update(rowstart,rowend,colstart,colend,NYPROB,&u[iz][0][0],&u[1-iz][0][0]);
         iz = 1 - iz;
      }
      free(temprowdatasend);
      free(temprowdatarecieve);
      /* Finally, send my portion of final results back to master */
      MPI_Send(&offsetrow, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
      MPI_Send(&offsetcol, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);      
      MPI_Send(&rows, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
      MPI_Send(&cols, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
     
      /* Prepare and send computing domain data to the master process */
      tempdatasend = (double *)malloc(rows*cols*sizeof(double));
      for(m=0; m<rows; m++)
      {
         for(n=0; n<cols; n++)
         {
              k = m*cols + n; 
              *(tempdatasend+k) = u[iz][offsetrow+m][offsetcol+n];
         }
      }              
      MPI_Send(tempdatasend, rows*cols, MPI_DOUBLE, MASTER, DONE, 
            MPI_COMM_WORLD);
      free(tempdatasend);      
      //printf("rows*cols = %d\n",rows*cols);
      MPI_Finalize();
   }
}





/**************************************************************************
 *  subroutine update
 ****************************************************************************/
void update(int rowstart, int rowend, int colstart, int colend, int ny, double *u1, double *u2)
{
   int ix, iy;
   for (ix = rowstart; ix <= rowend; ix++) 
      for (iy = colstart; iy <= colend; iy++) 
          *(u2+ix*ny+iy) = *(u1+ix*ny+iy)  + 
                           parms.cx * (*(u1+(ix+1)*ny+iy) +
                           *(u1+(ix-1)*ny+iy) - 
                           2.0 * *(u1+ix*ny+iy)) +
                           parms.cy * (*(u1+ix*ny+iy+1) +
                           *(u1+ix*ny+iy-1) - 
                           2.0 * *(u1+ix*ny+iy));
}




/*****************************************************************************
 *  subroutine inidat
 *****************************************************************************/
void inidat(int nx, int ny, double *u) 
{
   int ix, iy;

   for (ix = 0; ix <= nx-1; ix++) 
      for (iy = 0; iy <= ny-1; iy++)
          *(u+ix*ny+iy) = (double)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
}





/**************************************************************************
 * subroutine prtdat
 **************************************************************************/

/*
void prtdat(int nx, int ny, float *u1, char *fnam) {
int ix, iy;
FILE *fp;

fp = fopen(fnam, "w");
for (iy = ny-1; iy >= 0; iy--) {
  for (ix = 0; ix <= nx-1; ix++) {
    fprintf(fp, "%8.1f", *(u1+ix*ny+iy));
    if (ix != nx-1) 
      fprintf(fp, " ");
    else
      fprintf(fp, "\n");
    }
  }
fclose(fp);
}
*/

void prtdat(int nx, int ny, double *u1, char *fnam) 
{
   int ix, iy;
   FILE *fp;

   fp = fopen(fnam, "w");
   fprintf(fp,"VARIABLES = X, Y, T\n");
   fprintf(fp,"ZONE T="" I=%d, J=%d, F=POINT",nx,ny);
   for (iy = 0; iy < ny; iy++) 
   {
      for (ix = 0; ix < nx; ix++) 
      {
         fprintf(fp, "%d\t,%d\t,%E24.15\n", ix, iy, *(u1+ix*ny+iy));
         fprintf(fp, "\n");
      }
   }
   fclose(fp);
}