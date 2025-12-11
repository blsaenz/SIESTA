#include <cuda.h>
#include "math_functions.h"

// example of decaling constants/globals in CUDA device memory
// static __device__ __constant__ __TOptionData d_OptionData[MAX_OPTIONS];

// example of setting constant/global in CUDA device memory, i think
// cutilSafeCall( cudaMemcpyToSymbol(device_constant,host_constant, size) );

// CUDA Kernel: sia2_edd_GPU_kernel
// =======================================================================
// USAGE:
// sia2_edd_GPU_kernel<<<blocks, thread_per_block, grids>>>sia2_edd_solution_GPU(...)
//
// where blocks and threads/block are some division of the number of solutions to
// cells that need to be computed . Blocks can contain a maximum of 512/1024 threads,
// depending on the GPU.  Block numbers may be capped at 65535, so divide into blocks
// and threads accordingly to divy up total cells when calling kernel.  Standard is
// to use 256 threads/block

__global__ void sia2_edd_gpu_kernel_absorp (
    const int        *g_i32_in, 
    const float      *g_f32_in, 
    float            *g_f32_scr, 
    float            *g_f32_out, 
    int              i32_in_p_th,
    int              f32_in_p_th,
    int              f32_scr_p_th,
    int              f32_out_p_th,
    int              wl,
    int              zz)
{
      int i,j,k,ki,nilyr,nslyr,klev,wl_temp;
      float sd_or_bv,smalg,det,k2st;


      i = blockDim.x * blockIdx.x * i32_in_p_th + threadIdx.x;
      j = blockDim.x * blockIdx.x * f32_in_p_th + threadIdx.x;
      ki = blockDim.x * blockIdx.x * f32_scr_p_th + threadIdx.x + 7*(zz+1)*blockDim.x;  // write to 8th scratch slot - 7 contiguous needed in edd calc
      //  DEGUG - write to output instead of scrtch
      //ki = blockDim.x * blockIdx.x * f32_out_p_th + threadIdx.x;
      nilyr = g_i32_in[i];

      if (nilyr > 0) 
      {

      i += blockDim.x;
      nslyr = g_i32_in[i];
      klev = nilyr + nslyr;

      for (k=0;k<klev;k++)
      {

          sd_or_bv = g_f32_in[j];    
          j += blockDim.x;

          if (k < nslyr) 
          {  
              // find snow absorption 
              if (wl < wavl) 
              { 
                  k2st = aice[wl]*sd_or_bv*inv_iced + a_factor;    // account for ice part only	
              } else {
                  k2st = a_ice_ir*sd_or_bv*inv_iced + a_factor;  // near-IR absorption
              }
              j += 2*blockDim.x;
              wl_temp = wl;
             
          } else {
          
              // find ice and associated algae, detrital aborption
              if (wl < wavl) 
              { 
                  // Fresh Ice PAR abosrption ...
                  k2st = aice[wl]*(c1 - sd_or_bv) + sd_or_bv*awater[wl];
                  wl_temp = wl;
                  
              } else {
              
                  // Fresh Ice NIR abosrption ...
                  k2st = a_ice_ir;
                  wl_temp = wl-1;

              }

              // phyoplankton absorption                  
              smalg = g_f32_in[j];
              j += blockDim.x;
              k2st += aph[wl_temp]*max(c1,smalg*inv_c_chl*sd_or_bv);

              // detrital absorption                  
              det = g_f32_in[j];
              j += blockDim.x;
              k2st += det/smalg*inv_ad_denom*exp(p_08*float(wl_temp));

          }

          //if (threadIdx.x == 0) {
          //    printf("%i - k2st: %g %g %g %g %g %g\n",k,k2st,smalg,det/2.551,aice[wl_temp],awater[wl_temp]);
          //}
          
          // write out absorp to scratch memory
          g_f32_scr[ki] = k2st;
          ki += blockDim.x;          
          
          //g_f32_out[ki] = k2st;
          //ki += blockDim.x;
          //g_f32_out[ki] = smalg;
          //ki += blockDim.x;

      }

      }
}