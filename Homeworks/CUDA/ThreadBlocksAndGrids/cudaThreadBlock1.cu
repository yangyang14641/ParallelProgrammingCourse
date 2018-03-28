// program name: cudaThreadBlock.cu
// this program is designed for showing thread block example.
// author: Shane Cook (Nvidia .Inc) 
// modified by Yang Yang @ Peking University July 2017
// 
//
// Sketch diagram for thread blocks:
// -----------------------------------------------------------------------------------
// | thread block 0    || thread block 0   || thread block 1   || thread block 1    ||
// | thread bundle 0   || thread bundle 1  || thread bundle 0  || thread bundle 1   ||
// |  (thread 0~31)    || (thread 32~63)   || (thread 64~95)   || (thread 96~127)   ||
// -----------------------------------------------------------------------------------
//
// Adress space:
// -----------------------------------------------------------------------------------
// |                   ||                  ||                  ||                   ||
// |     adress        ||      adress      ||      adress      ||      adress       ||
// |     (0~31)        ||      (32~63)     ||      (64~95)     ||      (96~127)     ||
// -----------------------------------------------------------------------------------


// head files
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>




/* Cuda Kernel function: waht is my id */
__global__ void what_is_my_id(unsigned int * const block,
                              unsigned int * const thread,
                              unsigned int * const warp,
                              unsigned int * const calc_thread)
{
    /* Thread id equals to block index * block size + thread offset into the block */
    const unsigned int thread_idx = (blockIdx.x * blockDim.x) + threadIdx.x;

    block[thread_idx] = blockIdx.x;
    thread[thread_idx] = threadIdx.x;

    /* Calculate warp using built in variable warpSize */
    warp[thread_idx] = threadIdx.x / warpSize;

    calc_thread[thread_idx] = thread_idx;
}





#define ARRAY_SIZE 128
#define ARRAY_SIZE_IN_BYTES (sizeof(unsigned int) * ARRAY_SIZE)

/* Declare statically four arrays of ARRAY_SIZE each */

unsigned int cpu_block[ARRAY_SIZE];
unsigned int cpu_thread[ARRAY_SIZE];
unsigned int cpu_warp[ARRAY_SIZE];
unsigned int cpu_calc_thread[ARRAY_SIZE];

int main(void)
{
    /* Total thread count = 2 * 64 = 128 */
    const unsigned int num_blocks = 2;
    const unsigned int num_threads = 64;
    char ch;

    /* Decalre pointers for GPU based params */
    unsigned int * gpu_block;
    unsigned int * gpu_thread;
    unsigned int * gpu_warp;
    unsigned int * gpu_calc_thread;

    /* Declare loop counter for use later */
    unsigned int i;

    /* Allocate four arrays on the GPU */
    cudaMalloc((void **)&gpu_block, ARRAY_SIZE_IN_BYTES);         // Why here type is (void **)?
    cudaMalloc((void **)&gpu_thread, ARRAY_SIZE_IN_BYTES);
    cudaMalloc((void **)&gpu_warp, ARRAY_SIZE_IN_BYTES);
    cudaMalloc((void **)&gpu_calc_thread, ARRAY_SIZE_IN_BYTES);

    /* Execute our cuda kernel */
    what_is_my_id<<<num_blocks, num_threads>>>(gpu_block, gpu_thread, gpu_warp, gpu_calc_thread);

    /* Copy back the gpu results to the CPU, from display RAM to RAM in physical*/
    cudaMemcpy(cpu_block, gpu_block, ARRAY_SIZE_IN_BYTES, cudaMemcpyDeviceToHost);
    cudaMemcpy(cpu_thread, gpu_thread, ARRAY_SIZE_IN_BYTES, cudaMemcpyDeviceToHost);
    cudaMemcpy(cpu_warp, gpu_warp, ARRAY_SIZE_IN_BYTES, cudaMemcpyDeviceToHost);
    cudaMemcpy(cpu_calc_thread, gpu_calc_thread, ARRAY_SIZE_IN_BYTES, cudaMemcpyDeviceToHost);

    /* Free the arrays on the GPU as now we're done with them */
    cudaFree(gpu_block);
    cudaFree(gpu_thread);
    cudaFree(gpu_warp);
    cudaFree(gpu_calc_thread);

    /* Iterate through the arrays and print */
    for (i = 0; i < ARRAY_SIZE; i++)
    {  
       if(i>0)
       { 
         if(cpu_warp[i] == 1 && cpu_warp[i-1] == 0) printf("\n");
         if(cpu_block[i] == 1 && cpu_block[i-1] ==0) printf("\n\n");
       } 
       printf("Calculated Thread: %3u - Block: %2u - Warp %2u - Thread %3u\n",
               cpu_calc_thread[i], cpu_block[i], cpu_warp[i], cpu_thread[i]);       
    }

    /* To avoid program exit automatically */
    ch = getchar();
}