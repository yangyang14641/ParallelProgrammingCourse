#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

int BodyNum=0;
int TimeSteps=0;                                                                   
                                                                   
int main(int argc, char** argv )
{
     
  int n, t, i, j;
  double *pBody;//´æ´¢Á£×ÓµÄ»ù±¾ÐÅÏ¢£¬Ã¿¸öÁ£×ÓÕ¼ÓÃ4¸öÁ¬ÐøµÄ¸¡µãÊý£ºmass¡¢x¡¢y¡¢z
  double *pForce;//´æ´¢Á£×ÓµÄÊÜÁ¦£¬Ã¿¸öÁ£×ÓÕ¼ÓÃ3¸öÁ¬ÐøµÄ¸¡µãÊý£ºFx¡¢Fy¡¢Fz
  double fac, fx, fy, fz;
  double dx, dy, dz, sq, dist; 
  clock_t c_start, c_end;
  double run_time;
  char *pStr;
  FILE *fResult;

  for ( i=1; i<argc; i++ ) {
	  pStr=strstr(argv[i], "-s=");
	  if ( pStr!=NULL) sscanf(pStr, "-s=%d", &BodyNum);
	  pStr=strstr(argv[i], "-t=");
	  if ( pStr!=NULL) sscanf(pStr, "-t=%d", &TimeSteps);

  }

  if ( BodyNum*TimeSteps==0) {
	  printf("usage: ser_nbody -s=number-of-bodies -t=number-of-steps\n");
	  return 0;
  }
  c_start = clock(); 

  pForce = new double[3*BodyNum];
  pBody = new double[4*BodyNum];

  /*  Initialize mass and positions in array p to make a test case
      Initialize force to 0
  */
  for ( i=0; i<BodyNum; i++)
    {
      *(pBody+4*i) = 10.05 + i;
      *(pBody+4*i+1) = 30.0*i;
      *(pBody+4*i+2) = 20.0*i;
      *(pBody+4*i+3) = 10.0*i;
      *(pForce+3*i) = 0;
      *(pForce+3*i+1) = 0;
      *(pForce+3*i+2) = 0;
    }

  t = 0;
  while ( t<TimeSteps){
    /*  Loop over points calculating force between each pair.*/

    for ( i=0; i<BodyNum; i++ )
      for ( j=i+1; j<BodyNum; j++ )
        {/*Calculate force between particle i and j according to Newton's Law*/
          dx = *(pBody+4*i+1) - *(pBody+4*j+1);
          dy = *(pBody+4*i+2) - *(pBody+4*j+2);
          dz = *(pBody+4*i+3) - *(pBody+4*j+3);
          sq = dx*dx + dy*dy + dz*dz;
          dist = sqrt(sq);
          fac = (*(pBody+4*i)) * (*(pBody+4*j)) / ( dist * sq );
          fx = fac * dx;
          fy = fac * dy;
          fz = fac * dz;

          /*Add in force and opposite force to particle i and j */
          *(pForce+3*i) = *(pForce+3*i) - fx;
          *(pForce+3*i+1) = *(pForce+3*i+1) - fy;
          *(pForce+3*i+2) = *(pForce+3*i+2) - fz;
          *(pForce+3*j) = *(pForce+3*j) + fx;
          *(pForce+3*j+1) = *(pForce+3*j+1) + fy;
          *(pForce+3*j+2) = *(pForce+3*j+2) + fz;
        }
    for ( i=0; i<BodyNum; i++ ){ 
	  *(pBody+4*i+1) = *(pBody+4*i+1) + (*(pForce+3*i)) / (*(pBody+4*i));
	  *(pForce+3*i) = 0;
	  *(pBody+4*i+2) = *(pBody+4*i+2) + (*(pForce+3*i+1)) / (*(pBody+4*i));
	  *(pForce+3*i+1) = 0;
	  *(pBody+4*i+3) = *(pBody+4*i+3) + (*(pForce+3*i+2)) / (*(pBody+4*i));
	  *(pForce+3*i+2) = 0;
    }
    t++;
  }


  fResult=fopen("result_ser_nbody.txt", "w");

  char result[50];
  for (i=0; i<BodyNum; i++)   { 	  
	  sprintf(result, "(%10.4f %10.4f %10.4f %10.4f)\n", *(pBody+4*i), *(pBody+4*i+1), *(pBody+4*i+2), *(pBody+4*i+3));
	  fwrite(result, sizeof(char), strlen(result), fResult);
   }
   fclose(fResult);

  delete[] pForce;
  delete[] pBody;  
  
  c_end =clock();
  run_time = (double)( c_end - c_start) / CLOCKS_PER_SEC; 
  printf("runtime is : %f\n", run_time);
}
