// program name: cudaThreadGrid.cu
// this program is designed for showing thread grid example.
// author: Shane Cook (Nvidia .Inc) 
// modified by Yang Yang @ Peking University July 2017
// 
//
// built in variables:
// gridDim.x -- number of thread blocks in X dim of thread grid 
// gridDim.y -- number of thread blocks in Y dim of thread grid
//
// blockDim.x -- number of threads in X dim of thread block
// blockDim.y -- number of threads in Y dim of thread block
//
// threadIdx.x -- thread index in X dim of thread block
// threadIdx.y -- thread index in Y dim of thread block
//
//
// Sketch diagram for thread grid for an array mapping:
//  o------> X
//  |
//  |
//  V Y
// ---------------------------------------------------------------------------------------------------- ---  ---     ---
// | array element 0  || array element 1  || array elemnt 2   || array element 3  || array element 4  |  ^    ^       ^
// |     X = 0        ||      X = 1       ||     X = 2        ||      X = 3       ||      X = 4       |  |    |       |
// |     Y = 0        ||      Y = 0       ||     Y = 0        ||      Y = 0       ||      Y = 0       |  |    V       |
// ---------------------------------------------------------------------------------------------------|  |   ---      |
// | array element 5  || array element 6  || array element 7  || array element 8  || array element 9  |  |blockDim.y  |
// |     X = 0        ||      X = 1       ||     X = 2        ||      X = 3       ||      X = 4       |  |            V
// |     Y = 1        ||      Y = 1       ||     Y = 1        ||      Y = 1       ||      Y = 1       |  |           ---
// ----------------------------------------------------------------------------------------------------  |         threadIdx.y
// | array element 10 || array element 11 || array element 12 || array element 13 || array element 14 |  |
// |     X = 0        ||      X = 1       ||     X = 2        ||      X = 3       ||      X = 4       |  |
// |     Y = 1        ||      Y = 1       ||     Y = 1        ||      Y = 1       ||      Y = 1       |  v
// ---------------------------------------------------------------------------------------------------- --- gridDim.y
//|<--------------------------------------------(gridDim.x)------------------------------------------>|
//|<---(blockDim.x)-->|
//|<------>| threadIdx.x

/*--------------------------------------------------------------------------------------------------------*/
// head files
#include <stdio.h>
#include <stdlib.h>
//#include <conio.h>




/* Cuda Kernel function: waht is my id */
__global__ void what_is_my_id_2d_A(unsigned int * const block_x,
                                   unsigned int * const block_y,
                                   unsigned int * const thread,
                                   unsigned int * const calc_thread,
                                   unsigned int * const x_thread,
                                   unsigned int * const y_thread,
                                   unsigned int * const grid_dimx,
                                   unsigned int * const grid_dimy,
                                   unsigned int * const block_dimx,
                                   unsigned int * const block_dimy)
{
    /* Thread absolute id and id in X dim and Y dim */
    const unsigned int idx        = (blockIdx.x * blockDim.x) + threadIdx.x; 
    const unsigned int idy        = (blockIdx.y * blockDim.y) + threadIdx.y;
    const unsigned int thread_idx = ((gridDim.x * blockDim.x) * idy) + idx;

    block_x[thread_idx] = blockIdx.x;
    block_y[thread_idx] = blockIdx.y;
    thread[thread_idx] = threadIdx.x;
    calc_thread[thread_idx] = thread_idx;
    x_thread[thread_idx] = idx;
    y_thread[thread_idx] = idy;
    grid_dimx[thread_idx] = gridDim.x;
    grid_dimy[thread_idx] = gridDim.y;
    block_dimx[thread_idx] = blockDim.x;
    block_dimy[thread_idx] = blockDim.y;
}





/* Macro definition */
#define ARRAY_SIZE_X 32
#define ARRAY_SIZE_Y 16
#define ARRAY_SIZE_IN_BYTES ((ARRAY_SIZE_X) * (ARRAY_SIZE_Y) * (sizeof(unsigned int)))

/* Declare statically four arrays of ARRAY_SIZE each */
unsigned int cpu_block_x[ARRAY_SIZE_Y][ARRAY_SIZE_X];
unsigned int cpu_block_y[ARRAY_SIZE_Y][ARRAY_SIZE_X];
unsigned int cpu_thread[ARRAY_SIZE_Y][ARRAY_SIZE_X];
unsigned int cpu_warp[ARRAY_SIZE_Y][ARRAY_SIZE_X];
unsigned int cpu_calc_thread[ARRAY_SIZE_Y][ARRAY_SIZE_X];
unsigned int cpu_xthread[ARRAY_SIZE_Y][ARRAY_SIZE_X];
unsigned int cpu_ythread[ARRAY_SIZE_Y][ARRAY_SIZE_X];
unsigned int cpu_grid_dimx[ARRAY_SIZE_Y][ARRAY_SIZE_X];
unsigned int cpu_grid_dimy[ARRAY_SIZE_Y][ARRAY_SIZE_X];
unsigned int cpu_block_dimx[ARRAY_SIZE_Y][ARRAY_SIZE_X];
unsigned int cpu_block_dimy[ARRAY_SIZE_Y][ARRAY_SIZE_X];




/* The main function */
int main(void)
{
    /* Total thread count = 32 * 4 = 128 */
    const dim3 threads_rect(32, 4);  /* 32 * 4 */
    const dim3 blocks_rect(1, 4);  
    
    /* Total thread count = 16 * 8 = 128 */
    const dim3 threads_square(16, 8);
    const dim3 blocks_square(2, 2);

    /* program pause wait for a getchar() in C++ */
    char ch;

    /* Decalre pointers for GPU based params */
    unsigned int * gpu_block_x;
    unsigned int * gpu_block_y;
    unsigned int * gpu_thread;
    unsigned int * gpu_warp;
    unsigned int * gpu_calc_thread;
    unsigned int * gpu_xthread;
    unsigned int * gpu_ythread;
    unsigned int * gpu_grid_dimx;
    unsigned int * gpu_grid_dimy;
    unsigned int * gpu_block_dimx;
    unsigned int * gpu_block_dimy;

    /* Allocate four arrays on the GPU */
    cudaMalloc((void **)&gpu_block_x, ARRAY_SIZE_IN_BYTES);         // Why here type is (void **)?
    cudaMalloc((void **)&gpu_block_y, ARRAY_SIZE_IN_BYTES);
    cudaMalloc((void **)&gpu_thread, ARRAY_SIZE_IN_BYTES);
    cudaMalloc((void **)&gpu_calc_thread, ARRAY_SIZE_IN_BYTES);
    cudaMalloc((void **)&gpu_xthread, ARRAY_SIZE_IN_BYTES);
    cudaMalloc((void **)&gpu_ythread, ARRAY_SIZE_IN_BYTES);
    cudaMalloc((void **)&gpu_grid_dimx, ARRAY_SIZE_IN_BYTES);
    cudaMalloc((void **)&gpu_grid_dimy, ARRAY_SIZE_IN_BYTES);
    cudaMalloc((void **)&gpu_block_dimx, ARRAY_SIZE_IN_BYTES);
    cudaMalloc((void **)&gpu_block_dimy, ARRAY_SIZE_IN_BYTES);

    /* Execute our cuda kernel */
    for (int kernel = 0; kernel < 2 ; kernel ++)
    {
      switch (kernel)
      { 
          case 0:
          {
            /* Excute our kernel function */    
            what_is_my_id_2d_A<<<blocks_rect, threads_rect>>>(gpu_block_x, gpu_block_y, gpu_thread, gpu_calc_thread, 
            gpu_xthread, gpu_ythread, gpu_grid_dimx, gpu_grid_dimy, gpu_block_dimx, gpu_block_dimy);
          } break;

          case 1:
          {
            what_is_my_id_2d_A<<<blocks_square, threads_square>>>(gpu_block_x, gpu_block_y, gpu_thread, gpu_calc_thread, 
            gpu_xthread, gpu_ythread, gpu_grid_dimx, gpu_grid_dimy, gpu_block_dimx, gpu_block_dimy);
          } break;

          default: exit(1); break;
      }

      /* Copy back the gpu results to the CPU, from display RAM to RAM in physical */
      cudaMemcpy(cpu_block_x, gpu_block_x, ARRAY_SIZE_IN_BYTES, cudaMemcpyDeviceToHost);         
      cudaMemcpy(cpu_block_y, gpu_block_y, ARRAY_SIZE_IN_BYTES, cudaMemcpyDeviceToHost);
      cudaMemcpy(cpu_thread, gpu_thread, ARRAY_SIZE_IN_BYTES, cudaMemcpyDeviceToHost);
      cudaMemcpy(cpu_calc_thread, gpu_calc_thread, ARRAY_SIZE_IN_BYTES, cudaMemcpyDeviceToHost);
      cudaMemcpy(cpu_xthread, gpu_xthread, ARRAY_SIZE_IN_BYTES, cudaMemcpyDeviceToHost);
      cudaMemcpy(cpu_ythread, gpu_ythread, ARRAY_SIZE_IN_BYTES, cudaMemcpyDeviceToHost);
      cudaMemcpy(cpu_grid_dimx, gpu_grid_dimx, ARRAY_SIZE_IN_BYTES, cudaMemcpyDeviceToHost);
      cudaMemcpy(cpu_grid_dimy, gpu_grid_dimy, ARRAY_SIZE_IN_BYTES, cudaMemcpyDeviceToHost);
      cudaMemcpy(cpu_block_dimx, gpu_block_dimx, ARRAY_SIZE_IN_BYTES, cudaMemcpyDeviceToHost);
      cudaMemcpy(cpu_block_dimy, gpu_block_dimy, ARRAY_SIZE_IN_BYTES, cudaMemcpyDeviceToHost);

      printf("\nKernel %d\n", kernel);
       
      /* Iterate through the arrays and print */
      for (int y = 0; y < ARRAY_SIZE_Y; y++)
      {
        for (int x = 0; x < ARRAY_SIZE_Y; x++)
        {
           printf("CT: %2u BKX: %1u BKY: %1u TID: %2u YTID: %2u XTID: %2u GDX: %1u GDY %1u BDX %1u BDY %1u\n",
                  cpu_calc_thread[y][x], cpu_block_x[y][x], cpu_block_y[y][x], cpu_thread[y][x], cpu_ythread[y][x],
                  cpu_xthread[y][x], cpu_grid_dimx[y][x], cpu_grid_dimx[y][x], cpu_block_dimx[y][x], cpu_block_dimy[y][x]);

           /* program pause and wait for a keyboard input */
           ch = getchar(); 
        }
      }
      /* waiting for any key so we can see the console window */
      printf("Press any key to continue\n");
      ch = getchar();

   }

      /* Free the arrays on the GPU as now we're done with them */
      cudaFree(gpu_block_x);
      cudaFree(gpu_block_y);
      cudaFree(gpu_thread);
      cudaFree(gpu_calc_thread);
      cudaFree(gpu_xthread);
      cudaFree(gpu_ythread);
      cudaFree(gpu_grid_dimx);
      cudaFree(gpu_grid_dimy);
      cudaFree(gpu_block_dimx);
      cudaFree(gpu_block_dimy);

      /* To avoid program exit automatically */
      ch = getchar();
}