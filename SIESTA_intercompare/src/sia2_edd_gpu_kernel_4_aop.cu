
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

__global__ void sia2_edd_gpu_kernel_aop (
    const int        *g_i32_in, 
    float            *g_f32_scr, 
    float            *g_f32_out, 
    int              i32_in_p_th,
    int              f32_scr_p_th,
    int              f32_out_p_th,
    int              zz)  // total layers
{
      int i,ki,kii,k,
      //kiii,
      nilyr,nslyr,klevp;

      // layer edd inputs
      float rdir_up,rdif_a_up,rdif_b_up,tdir_up,tdif_a_up, 
          tdif_b_up,trnlay_up;
      
      // layer aop outputs (4 out of 6, anyway)
      float trndir_up,trndir,trntdr_up,trntdr,trndif_up,trndif,
          rdndif_up,rdndif;

      // bookkeeping vars
      float refkm1,tdrrdir,tdndif;

      // find CUDA indices to cell data
      i = blockDim.x * blockIdx.x * i32_in_p_th + threadIdx.x;
      ki = blockDim.x * blockIdx.x * f32_scr_p_th + threadIdx.x;
      kii = blockDim.x * blockIdx.x * f32_scr_p_th + threadIdx.x + 7*(zz+1)*blockDim.x;
      //kiii = blockDim.x * blockIdx.x * f32_out_p_th + threadIdx.x;

      nilyr = g_i32_in[i];          

      if (nilyr > 0) 
      {

      i = i + blockDim.x;
      nslyr = g_i32_in[i];    // nslyr in
      klevp = nslyr + nilyr + 1; // , &  // number of radiation layers - 1

      // initialize
      trndir_up =   c1;
      trntdr_up =   c1;
      trndif_up =   c1;
      rdndif_up =   c0;

      //write out and aops layers
      g_f32_scr[kii] = trndir_up;
      kii += blockDim.x;
      g_f32_scr[kii] = trntdr_up;
      kii += blockDim.x;
      g_f32_scr[kii] = trndif_up;
      kii += blockDim.x;
      g_f32_scr[kii] = rdndif_up;
      kii += blockDim.x;

      for (k=1;k<klevp;k++) 
      { 

          
          // load layer above edd solution parameters;
          rdir_up = g_f32_scr[ki];
          ki = ki + blockDim.x;
          rdif_a_up = g_f32_scr[ki];
          ki = ki + blockDim.x;
          rdif_b_up = g_f32_scr[ki];
          ki = ki + blockDim.x;
          tdir_up = g_f32_scr[ki];
          ki = ki + blockDim.x;
          tdif_a_up = g_f32_scr[ki];
          ki = ki + blockDim.x;
          tdif_b_up = g_f32_scr[ki];
          ki = ki + blockDim.x;
          trnlay_up = g_f32_scr[ki];
          ki = ki + blockDim.x;

          // Calculate the solar beam transmission, total transmission, and
          // reflectivity for diffuse radiation from below at interface k, 
          // the top of the current layer k:
          //
          //              layers       interface
          //         
          //       ---------------------  kk-1 
          //                kk-1
          //       ---------------------  kk
          //                 kk
          //       ---------------------  

          trndir = trndir_up*trnlay_up;
          refkm1        = c1/(c1 - rdndif_up*rdif_a_up);
          tdrrdir       = trndir_up*rdir_up;
          tdndif        = trntdr_up - trndir_up;
          trntdr = trndir_up*tdir_up + 
           (tdndif + tdrrdir*rdndif_up)*refkm1*tdif_a_up;
          rdndif = rdif_b_up + 
            (tdif_b_up*rdndif_up*refkm1*tdif_a_up);
          trndif = trndif_up*refkm1*tdif_a_up;
          
           //DEBUG
          //if (k < 69) {
          //g_f32_out[kiii] = trntdr_up;
          //kiii += blockDim.x;
          //g_f32_out[kiii] = trndif_up;
          //kiii += blockDim.x;
          //}

 
          //write out and aops layers
          g_f32_scr[kii] = trndir;
          trndir_up = trndir;
          kii += blockDim.x;
          g_f32_scr[kii] = trntdr;
          trntdr_up = trntdr;
          kii += blockDim.x;
          g_f32_scr[kii] = trndif;
          trndif_up = trndif;
          kii += blockDim.x;
          g_f32_scr[kii] = rdndif;
          rdndif_up = rdndif;
          kii += blockDim.x;

      }
      }
}

 
 