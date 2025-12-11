
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

__global__ void sia2_edd_gpu_kernel_pur (
    const int        *g_i32_in, 
    const float      *g_f32_in, 
    float            *g_f32_out, 
    const float      *g_f32_scr, 
    int              i32_in_p_th,
    int              f32_in_p_th,
    int              f32_out_p_th,
    int              f32_scr_p_th,
    int              zz,
    int              wl)  // total layers
{
      int i,j,ki,kii,k,kiii,nilyr,nslyr,klev;

      // layer edd inputs
      float rdir,rdif_a,rdif_b,tdir,tdif_a, 
          tdif_b,trnlay;
      
      // layer aop outputs
      float trndir,trntdr,trndif,rdndif,rupdir_dn,rupdir_up,
          rupdif_dn,rupdif_up;

      // layer surface-relative transmissions
      float fdirup_dn,fdirup_up,fdirdn_dn,fdirdn_up,fdifup_dn,fdifup_up,
          fdifdn_dn,fdifdn_up;

      // bookkeeping vars
      float ed0_dir,ed0_dif,absorp,layer_in,layer_out,par_top,par_bot,
          par_mid,refk;

      // find CUDA indices to cell data
      i = blockDim.x * blockIdx.x * i32_in_p_th + threadIdx.x;

      nilyr = g_i32_in[i];          

      if (nilyr > 0) 
      {

      i = i + blockDim.x;
      nslyr = g_i32_in[i];    // nslyr in
      klev = nslyr + nilyr; // , &  // number of radiation layers - 1

      j = blockDim.x * blockIdx.x * f32_in_p_th + threadIdx.x + (zz*6 + 1 + wl*2)*blockDim.x;
      ki = blockDim.x * blockIdx.x * f32_scr_p_th + threadIdx.x + 
          (zz+1)*7*blockDim.x + ((klev+1)*4-1)*blockDim.x;
      kiii = blockDim.x * blockIdx.x * f32_scr_p_th + threadIdx.x +
          (klev*7-1)*blockDim.x;

      ed0_dir = g_f32_in[j];
      j = j + blockDim.x;      
      ed0_dif = g_f32_in[j];

      // zero output vectors

      if (wl == 0) 
      {
          // set kii for forward iteration
          kii = blockDim.x * blockIdx.x * f32_out_p_th + threadIdx.x;
          for (k=0;k<zz;k++)
          {          
              g_f32_out[kii] = c0;
              kii += blockDim.x;
              g_f32_out[kii] = c0;
              kii += blockDim.x;
          }
      }

      // reset kii for output
      kii = blockDim.x * blockIdx.x * f32_out_p_th + threadIdx.x + (klev*2-1)*blockDim.x;

      // initialize bottom interface for aops not already calculated
      rupdir_up =   c0;
      rupdif_up =   c0;

      // initial load of aops - remembering to load backwards from input
      rdndif =   g_f32_scr[ki];
      ki -= blockDim.x;
      trndif =   g_f32_scr[ki];
      ki -= blockDim.x;
      trntdr =   g_f32_scr[ki];
      ki -= blockDim.x;
      trndir =   g_f32_scr[ki];
      ki -= blockDim.x;

      // find top-relative interface transmission
      refk          = c1/(c1 - rdndif*rupdif_up);
      fdirup_up = (trndir*rupdir_up + 
                       (trntdr-trndir)   
                       *rupdif_up)*refk;
      fdirdn_up = trndir + (trntdr  
                      - trndir + trndir
                      *rupdir_up*rdndif)*refk;
      fdifup_up = trndif*rupdif_up*refk;
      fdifdn_up = trndif*refk;

      for (k=klev-1;k>=0;k--) 
      { 

          // transfer upper interface to lower interface, as we move up through layers
          fdirup_dn = fdirup_up;
          fdirdn_dn = fdirdn_up;
          fdifup_dn = fdifup_up;
          fdifdn_dn = fdifdn_up;
          rupdir_dn = rupdir_up;
          rupdif_dn = rupdif_up;

          // load layer edd solution parameters - backwards load!
          trnlay = g_f32_scr[kiii];
          kiii -= blockDim.x;
          tdif_b = g_f32_scr[kiii];
          kiii -= blockDim.x;
          tdif_a = g_f32_scr[kiii];
          kiii -= blockDim.x;
          tdir = g_f32_scr[kiii];
          kiii -= blockDim.x;
          rdif_b = g_f32_scr[kiii];
          kiii -= blockDim.x;
          rdif_a = g_f32_scr[kiii];
          kiii -= blockDim.x;
          rdir = g_f32_scr[kiii];
          kiii -= blockDim.x;

          // calculate remaining 2/6 aop
          refk = c1/(c1 - rdif_b*rupdif_dn);
          rupdir_up = rdir + 
            ( trnlay*rupdir_dn + 
             (tdir-trnlay)*rupdif_dn ) * 
              refk*tdif_b;
          rupdif_up = rdif_a + 
              tdif_a*rupdif_dn* 
              refk*tdif_b;

          // load remaining pre-calcualted 4/6 aops
          rdndif =   g_f32_scr[ki];
          ki -= blockDim.x;
          trndif =   g_f32_scr[ki];
          ki -= blockDim.x;
          trntdr =   g_f32_scr[ki];
          ki -= blockDim.x;
          trndir =   g_f32_scr[ki];
          ki -= blockDim.x;
          
          // find top-relative interface transmission
          refk          = c1/(c1 - rdndif*rupdif_up);
          fdirup_up = (trndir*rupdir_up + 
                           (trntdr-trndir)   
                           *rupdif_up)*refk;
          fdirdn_up = trndir + (trntdr  
                          - trndir + trndir
                          *rupdir_up*rdndif)*refk;
          fdifup_up = trndif*rupdif_up*refk;
          fdifdn_up = trndif*refk;


          // find total PAR in/out or layer
            layer_in = 
                (fdirdn_up + fdirup_dn) * ed0_dir +
                (fdifdn_up + fdifup_dn) * ed0_dif;
            
            layer_out = 
                (fdirup_up + fdirdn_dn) * ed0_dir +
                (fdifup_up + fdifdn_dn) * ed0_dif;


            // find absorption
            absorp = abs(layer_in - layer_out);  

            if (wl < wavl) 
            {

                absorp = absorp*quanta_to_watts[wl]*c10;

                // find and store PUR
                par_top = 
                    (fdirdn_dn+fdirup_dn) * ed0_dir +
                    (fdifdn_dn+fdifup_dn) * ed0_dif;
                par_bot = 
                    (fdirdn_up+fdirup_up) * ed0_dir +
                    (fdifdn_up+fdifup_up) * ed0_dif;
    
                // calc mean layer light using simple average
                par_mid = (par_top+par_bot)*p5;
    
                //calc mean layer light using exponential curve
                //if (par_bot .eq. 0.) then                                
                //    par_mid = (par_top+par_bot)*c_5
                //else  
                //    tmp_k = log(par_bot/par_top)/ice(ic,mi)%th(m_row)  
                //    par_mid = (par_bot - par_top)/tmp_k/ice(ic,mi)%th(m_row)
                //endif
                
                // store PAR
                g_f32_out[kii] = g_f32_out[kii] + par_mid*aph[wl]/aph_max*c10;
                //g_f32_out[kii] = trndir;

            }
            
            // store absoption
            kii -= blockDim.x;
            g_f32_out[kii] = g_f32_out[kii] + absorp;
            //g_f32_out[kii] = rdndif;
            kii -= blockDim.x;
                        
      }
      }
}