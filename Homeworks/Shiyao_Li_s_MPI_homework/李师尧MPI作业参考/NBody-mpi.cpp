#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h> 
#include <unistd.h>
#include <pthread.h>
#include<mpi.h>
#define NANO           1000000000
#define Max_Thread_Num 256
#define REAL           double

int BodyNum=0;
int TimeSteps=0;
REAL *body;
REAL *force;
MPI_Status status;
int myid, size;
int totalThread;
int bcastint[10];

pthread_cond_t   cond;
pthread_mutex_t  mtx;

double serial(); 
double mpi_nbody_blocking();
double mpi_nbody_noblocking();

double pthread_reduce();
void *reduce_worker(void *arg);
int  freeWorker;  

double pthread_divide();
void *divide_worker(void *arg);
// struct STATUS {
       // int lbound;
       // int ubound;
       // int *task;
       // REAL *force;
       // pthread_mutex_t  mtx;
// } *status;                                                                
                                                                   
int main(int argc, char** argv ) {     
	int i;
	REAL ser_time, red_time, div_time,mpi_time_blocking,mpi_time_noblocking;
	char *pStr;
	
   
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &myid );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
	if (myid==0) {
		for ( i=1; i<argc; i++ ) {
		  pStr=strstr(argv[i], "-s=");
		  if ( pStr!=NULL) sscanf(pStr, "-s=%d", &BodyNum);
		  pStr=strstr(argv[i], "-t=");
		  if ( pStr!=NULL) sscanf(pStr, "-t=%d", &TimeSteps);

		}
		if ( BodyNum*TimeSteps==0) {
		  printf("usage: -s=number-of-bodies -t=number-of-steps\n");
		  return 0;
		}
		bcastint[0]=BodyNum;
		bcastint[1]=TimeSteps;
	}
	MPI_Bcast(bcastint,2,MPI_INT,0,MPI_COMM_WORLD);
	BodyNum=bcastint[0];
	TimeSteps=bcastint[1];



	if (myid==0) {
		ser_time = serial();
		printf("serial: %f\n", ser_time);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	mpi_time_blocking=mpi_nbody_blocking();
	if (myid==0) printf("mpi-blocking: %f  speedup=%f\n",mpi_time_blocking,ser_time/mpi_time_blocking);
	mpi_time_noblocking=mpi_nbody_noblocking();
	if (myid==0) printf("mpi-noblocking: %f  speedup=%f\n",mpi_time_noblocking,ser_time/mpi_time_noblocking);
	free(body); 
	free(force);  
	MPI_Finalize();

}

double serial() {
       REAL fac, fx, fy, fz;
       REAL dx, dy, dz, sq, dist; 
       int t, i, j, bi,bj,fi,fj;
       struct timespec ts,te;
       double result;
        
		body = (REAL*)malloc(4*BodyNum*sizeof(REAL));
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

double mpi_nbody_blocking() {
       REAL fac, fx, fy, fz;
       REAL dx, dy, dz, sq, dist; 
       int t, i, j, bi,bj,fi,fj;
       struct timespec ts,te;
       double result;
        int locbodynum,loca,locb,locnum;
		REAL *bodysend,*forcesend;
		if (myid<size-1) {
			locbodynum=int(BodyNum/size);
			loca=0+myid*locbodynum;
			locb=loca+locbodynum-1;
		}
		else{
			locbodynum=BodyNum-int(BodyNum/size)*(size-1);
			loca=0+(size-1)*int(BodyNum/size);
			locb=BodyNum-1;			
		}			
		//printf("%d: %d %d %d \n",myid,loca, locb, locbodynum);
		body=(REAL*)malloc(4*BodyNum*sizeof(REAL));		
		force=(REAL*)malloc(3*BodyNum*sizeof(REAL));
       /*  Initialize mass and positions in array p to make a test case   */
       for ( i=0; i<BodyNum; i++)    {
          bi=4*i;
		  body[bi] = 10.05 + i;
          body[bi+1] = 30.0*i;
          body[bi+2] = 20.0*i;
          body[bi+3] = 10.0*i;
       }
              
       clock_gettime(CLOCK_REALTIME, &ts);
          
       for ( i=0; i<3*locbodynum; i++)  force[loca+i] = 0;
       
       t = 0;
	   int j1;
       while ( t<TimeSteps){
           /*  Loop over points calculating force between each pair.*/
           for ( i=0; i<locbodynum; i++ ) {
              bi = loca+4*i;              
              fi = loca+3*i;       
			  if (2*i<=BodyNum) 
				  locnum=int(BodyNum/2);
			  else
				  locnum=BodyNum-int(BodyNum/2);
              for ( j1=0; j1<locnum; j1++ )  {
				 j=(j1+i+1)%BodyNum;
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
		   
			for (i=0; i<size; i++){
				if (i==myid) continue;			   
				MPI_Send(force,3*BodyNum,MPI_DOUBLE,i,myid,MPI_COMM_WORLD);
			}
			for (i=0; i<size; i++) {
				if (i==myid) continue;	
				REAL *temp;
				temp=(REAL*)malloc(3*BodyNum*sizeof(REAL));
				MPI_Recv(temp,3*BodyNum,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
				for (j=loca; j<=locb; j++) {
					fj=3*j;
					force[fj]+=temp[fj];
					force[fj+1]+=temp[fj+1];
					force[fj+2]+=temp[fj+2];
				}
					
				delete temp;
			}
			
			// REAL *temp;
			// temp=(REAL*)malloc(3*size*BodyNum*sizeof(REAL));
			// for (i=0; i<size; i++) {
				// if (i==myid) continue;					
				// MPI_Recv(temp+3*i*BodyNum,3*BodyNum,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
			// }
			// for (i=0; i<size; i++) {
				// if (i==myid) continue;	
				// for (j=loca; j<=locb; j++) {
					// fj=3*j;
					// fi=fj+size*i;
					// force[fj]+=temp[fi];
					// force[fj+1]+=temp[fi+1];
					// force[fj+2]+=temp[fi+2];
				// }
			// }
			// delete temp;
			
           for ( i=0; i<locbodynum; i++ ){ 
              bi = 4*i;              
              fi = 3*i;
              body[bi+1] = body[bi+1] + force[fi] / body[bi];
              force[fi] = 0;
              body[bi+2] = body[bi+2] + force[fi+1] / body[bi];
              force[fi+1] = 0;
              body[bi+3] = body[bi+3] + force[fi+2] / body[bi];
              force[fi+2] = 0;
           }
		   
		   for (i=0; i<size; i++){
				if (i==myid) continue;			   
				MPI_Send(body+4*loca,4*locbodynum,MPI_DOUBLE,i,myid,MPI_COMM_WORLD);
			}
			for (i=0; i<size; i++) {
				if (i==myid) continue;
				MPI_Recv(body+4*locbodynum*i,4*locbodynum,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&status);
			}
			
		   
           t++;
       }
       free(force);
       
       clock_gettime(CLOCK_REALTIME, &te);
       result = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;
       
       FILE *fResult=fopen("result_mpi_nbody.txt", "w");
       char str[50];
       for (i=0; i<BodyNum; i++)   { 	  
		  sprintf(str, "(%10.4f %10.4f %10.4f %10.4f)\n", body[4*i], body[4*i+1], body[4*i+2], body[4*i+3]);
		  fwrite(str, sizeof(char), strlen(str), fResult);
       }
       fclose(fResult);
       return result;
}

double mpi_nbody_noblocking() {
       REAL fac, fx, fy, fz;
       REAL dx, dy, dz, sq, dist; 
       int t, i, j, bi,bj,fi,fj;
       struct timespec ts,te;
       double result;
	   MPI_Request *request;
	   MPI_Status  *state;
        int locbodynum,loca,locb,locnum;
		REAL *bodysend,*forcesend;
		if (myid<size-1) {
			locbodynum=int(BodyNum/size);
			loca=0+myid*locbodynum;
			locb=loca+locbodynum-1;
		}
		else{
			locbodynum=BodyNum-int(BodyNum/size)*(size-1);
			loca=0+(size-1)*int(BodyNum/size);
			locb=BodyNum-1;			
		}			
		//printf("%d: %d %d %d \n",myid,loca, locb, locbodynum);
		body=(REAL*)malloc(4*BodyNum*sizeof(REAL));		
		force=(REAL*)malloc(3*BodyNum*sizeof(REAL));
       /*  Initialize mass and positions in array p to make a test case   */
       for ( i=0; i<BodyNum; i++)    {
          bi=4*i;
		  body[bi] = 10.05 + i;
          body[bi+1] = 30.0*i;
          body[bi+2] = 20.0*i;
          body[bi+3] = 10.0*i;
       }
              
       clock_gettime(CLOCK_REALTIME, &ts);
          
       for ( i=0; i<3*locbodynum; i++)  force[loca+i] = 0;
       
       t = 0;
	   int j1;
	   int req;
	   request=new MPI_Request[size*2];
	   state=new MPI_Status[size*2];
       while ( t<TimeSteps){
           /*  Loop over points calculating force between each pair.*/
           for ( i=0; i<locbodynum; i++ ) {
              bi = loca+4*i;              
              fi = loca+3*i;   
			  if (2*i<=BodyNum) 
				  locnum=int(BodyNum/2);
			  else
				  locnum=BodyNum-int(BodyNum/2);			  
              for ( j1=0; j1<locnum; j1++ )  {
				 j=(j1+i+1)%BodyNum;
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
		   
			for (i=0; i<size; i++){
				if (i==myid) continue;			   
				MPI_Isend(force,3*BodyNum,MPI_DOUBLE,i,myid,MPI_COMM_WORLD,&request[i]);
			}
			REAL *temp;
			temp=(REAL*)malloc(3*size*BodyNum*sizeof(REAL));
			req=0;
			for (i=0; i<size; i++) {
				if (i==myid) continue;					
				MPI_Irecv(temp+3*i*BodyNum,3*BodyNum,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&request[size+req]);
				req++;
			}
			MPI_Waitall(size-1,request+size,state+size);
			//MPI_Barrier(MPI_COMM_WORLD);
			for (i=0; i<size; i++) {
				if (i==myid) continue;	
				for (j=loca; j<=locb; j++) {
					fj=3*j;
					fi=fj+size*i;
					force[fj]+=temp[fi];
					force[fj+1]+=temp[fi+1];
					force[fj+2]+=temp[fi+2];
				}
			}
			delete temp;
			
		   
		   
           for ( i=0; i<locbodynum; i++ ){ 
              bi = 4*i;              
              fi = 3*i;
              body[bi+1] = body[bi+1] + force[fi] / body[bi];
              force[fi] = 0;
              body[bi+2] = body[bi+2] + force[fi+1] / body[bi];
              force[fi+1] = 0;
              body[bi+3] = body[bi+3] + force[fi+2] / body[bi];
              force[fi+2] = 0;
           }
		   
		   for (i=0; i<size; i++){
				if (i==myid) continue;			   
				MPI_Isend(body+4*loca,4*locbodynum,MPI_DOUBLE,i,myid,MPI_COMM_WORLD,&request[i]);
			}
			req=0;
			for (i=0; i<size; i++) {
				if (i==myid) continue;
				MPI_Irecv(body+4*locbodynum*i,4*locbodynum,MPI_DOUBLE,i,i,MPI_COMM_WORLD,&request[size+req]);
				req++;
			}
			MPI_Waitall(size-1,request+size,state+size);
			//MPI_Barrier(MPI_COMM_WORLD);
		   
           t++;
       }
	   
       
       
       clock_gettime(CLOCK_REALTIME, &te);
       result = te.tv_sec - ts.tv_sec + (double)(te.tv_nsec-ts.tv_nsec)/NANO;
       
       FILE *fResult=fopen("result_mpi_nbody.txt", "w");
       char str[50];
       for (i=0; i<BodyNum; i++)   { 	  
		  sprintf(str, "(%10.4f %10.4f %10.4f %10.4f)\n", body[4*i], body[4*i+1], body[4*i+2], body[4*i+3]);
		  fwrite(str, sizeof(char), strlen(str), fResult);
       }
       fclose(fResult);
       return result;
}