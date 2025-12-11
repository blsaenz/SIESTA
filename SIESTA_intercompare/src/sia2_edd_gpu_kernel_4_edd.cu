
#include <cuda.h>
#include "math_functions.h"

// example of decaling constants/globals in CUDA device memory
// static __device__ __constant__ __TOptionData d_OptionData[MAX_OPTIONS];
/**
const float c0 = 0.0f;
const float c1 = 1.0f;
const float c1p5 = 1.5f;
const float c2 = 2.0f;
const float c3 = 3.0f;
const float c4 = 4.0f;

const float p01 = 0.01f;
const float p001 = 0.002f;
const float p5 = 0.5f;
const float p75 = 0.75f;

// refractive index for sea ice, water; pre-computed, band-independent,
// diffuse fresnel reflectivities
const float refindx = 1.310f;  // refractive index of sea ice (used for water also)
const float cp063   = 0.063f;  // diffuse fresnel reflectivity from above
const float cp455   = 0.455f;  // diffuse fresnel reflectivity from below

__device__ float gauspt[8];  // gaussian angles (radians)
__device__ float gauswt[8];  // gaussian weights      
**/

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

__global__ void sia2_edd_gpu_kernel_edd (
    const int        *g_i32_in, 
    const float      *g_f32_in, 
    float            *g_f32_scr, 
    float            *g_f32_out, 
    int              i32_in_p_th,
    int              f32_in_p_th,
    int              f32_scr_p_th,
    int              f32_out_p_th,
    int              zz,
    float            lda_m)  // snow depth distribution multiplier
{

      int //(kind=int_kind) :: & 
         kfrsnl;      // radiation interface index for fresnel layer
 
      float //(kind=dbl_kind), dimension (0:klev,icells_DE) :: &
//         trndir_tot,
         rdir,  // layer reflectivity to direct radiation
         rdif_a,  // layer reflectivity to diffuse radiation from above
         rdif_b,  // layer reflectivity to diffuse radiation from below
         tdir,  // layer transmission to direct radiation (solar beam + diffuse)
         tdif_a,  // layer transmission to diffuse radiation from above
         tdif_b,  // layer transmission to diffuse radiation from below
         trnlay;      // solar beam transm for layer (direct beam only)
 
      int //(kind=int_kind) :: & 
         i,
         j,
         k,           // level index
         ki,
         kii,
//         kiii,
         klev,
         srftyp,
         nilyr,
         nslyr;
 
      float //(kind=dbl_kind) :: &
         tautot   , // layer optical depth
         wtot     , // layer single scattering albedo
         gtot     , // layer asymmetry parameter
         ts       , // layer scaled extinction optical depth
         ws       , // layer scaled single scattering albedo
         gs       , // layer scaled asymmetry parameter
         rintfc   ; // reflection (multiple) at an interface
 
      // perpendicular and parallel relative to plane of incidence and scattering
      float //(kind=dbl_kind) :: &
         R1       ,  // perpendicular polarization reflection amplitude
         R2       ,  // parallel polarization reflection amplitude
         T1       ,  // perpendicular polarization transmission amplitude
         T2       ,  // parallel polarization transmission amplitude
         Rf_dir_a ,  // fresnel reflection to direct radiation
         Tf_dir_a ,  // fresnel transmission to direct radiation
         Rf_dif_a ,  // fresnel reflection to diff radiation from above
         Rf_dif_b ,  // fresnel reflection to diff radiation from below
         Tf_dif_a ,  // fresnel transmission to diff radiation from above
         Tf_dif_b;   // fresnel transmission to diff radiation from below
  
      float //(kind=dbl_kind) :: &
         mu0      ,  // cosine solar zenith angle incident
         mu0n;         // cosine solar zenith angle in medium
 
      float //(kind=dbl_kind) :: &
         lm       ,  // temporary for el
         mu       ,  // cosine solar zenith for either snow or water
         ne;         // temporary for n
 
 
      float //(kind=dbl_kind) :: &
         alp      ,  // temporary for alpha
         gam      ,  // temporary for gamma
         ue       ,  // temporary for u
//         arg      ,  // exponential argument
         extins   ,  // extinction
         amg      ,  // alp - gam
         apg;          // alp + gam
 
      int //(kind=int_kind) :: &
         ng;           // gaussian integration index
 
      float //(kind=dbl_kind) :: &
         gwt      ,  // gaussian weight
         swt      ,  // sum of weights
         trn      ,  // layer transmission
         rdr      ,  // rdir for gaussian integration
         tdr      ,  // tdir for gaussian integration
         smr      ,  // accumulator for rdif gaussian integration
         smt      ,  // accumulator for tdif gaussian integration
         r2st     ,
         k2st;   

      // find CUDA indices to cell data
      i = blockDim.x * blockIdx.x * i32_in_p_th + threadIdx.x;
      j = blockDim.x * blockIdx.x * f32_in_p_th + threadIdx.x + (zz*3)*blockDim.x;
      ki = blockDim.x * blockIdx.x * f32_scr_p_th + threadIdx.x + 7*(zz+1)*blockDim.x;
      kii = blockDim.x * blockIdx.x * f32_scr_p_th + threadIdx.x;
//      kiii = blockDim.x * blockIdx.x * f32_out_p_th + threadIdx.x;

      // Copy slow global memory into local mem for processing
      nilyr = g_i32_in[i];

      if (nilyr > 0) 
      {

      i = i + blockDim.x;
      nslyr = g_i32_in[i];    // nslyr in
      i = i + blockDim.x;
      srftyp = g_i32_in[i];    // srftyp in
      klev = nslyr + nilyr; // , &  // number of radiation layers - 1
     
       // compute level of fresnel refraction
       if( srftyp < 2 )
       {
         // if snow over sea ice or bare sea ice, fresnel level is
         // at base of sea ice SSL (and top of the sea ice DL); the
         // snow SSL counts for one, then the number of snow layers,
         // then the sea ice SSL which also counts for one:
         kfrsnl = nslyr + 2;
       } else {
         // if ponded sea ice, fresnel level is the top of the pond 
         kfrsnl = 0;
       }

   
        // mu0 is cosine solar zenith angle above the fresnel level; make 
        // sure mu0 is large enough for stable and meaningful radiation
        // solution: .01 is like sun just touching horizon with its lower edge

        mu0  = fmax(g_f32_in[j],p01);   // coszen in
        j = j + blockDim.x;

         // mu0n is cosine solar zenith angle used to compute the layer
         // Delta-Eddington solution; it is initially computed to be the
         // value below the fresnel level, i.e. the cosine solar zenith 
         // angle below the fresnel level for the refracted solar beam:

         mu0n = sqrt(c1-((c1-mu0*mu0)/(refindx*refindx)));
 
        // proceed down one layer at a time; if the total transmission (trndir_tot) to
        // the interface just above a given layer is less than trmin, then no
        // Delta-Eddington computation for that layer is done.

 
        // init these vars to begin calc of trndir_tot, which tells calc
        // to stop of transmission is 1/1000 of surface
//            trndir_tot = c1;
//            trnlay = c1;
        
        // begin main level loop
        for (k=0;k<klev;k++) 
        { 
  
//              trndir_tot = trndir_tot*trnlay;
   
          // compute next layer Delta-eddington solution only if total transmission
          // of radiation to the interface just above the layer exceeds trmin.
//              if (trndir_tot > p001 ) 
//              {
   
           // bring in iops, thickness
           gtot    = g_f32_in[j];      // g in
           j = j + blockDim.x;
           r2st    = g_f32_in[j];   // scattering(r2st) in               
           j = j + blockDim.x;
           tautot = g_f32_in[j];    // dth in
           j = j + blockDim.x;

           k2st  = g_f32_scr[ki];  // tau in
           ki += blockDim.x;

           if (k < nslyr) 
           {
               tautot  = tautot * lda_m;
           }
           
           // find w0, tau
           wtot = r2st/(r2st+k2st);
           tautot = tautot*(r2st+k2st);
           

           ts   = (c1 - wtot*gtot*gtot)*tautot;
           ws   = (c1 - gtot*gtot)*wtot/(c1 - wtot*gtot*gtot);
           gs   = (gtot - gtot*gtot)/(c1 - gtot*gtot);
   
           lm = sqrt(c3*(c1-ws)*(c1 - ws*gs));
           //ue   = u_macro(ws,gs,lm);
           ue = c1p5*(c1 - ws*gs)/lm;
   
           // if level k is above fresnel level and the cell is non-pond, use the
           // non-refracted beam instead
   
           if( srftyp < 2 && k < kfrsnl ) 
           {
              mu0n = mu0;
           }

           //extins = max(exp_min, exp(-lm*ts/dx_exp))
           extins = exp(-lm*ts);
           //ne = n_macro(ue,extins);
           ne = ((ue+c1)*(ue+c1)/extins ) - ((ue-c1)*(ue-c1)*extins);
   
           // first calculation of rdif, tdif using Delta-Eddington formulas
   
           rdif_a = (ue+c1)*(ue-c1)*(c1/extins - extins)/ne;
           tdif_a = c4*ue/ne;
   
           // evaluate rdir,tdir for direct beam
           //trnlay(k) = max(exp_min, exp(-ts/(mu0n*dx_exp)))
           trnlay = exp(-ts/mu0n);
           //alp = alpha(ws,mu0n,gs,lm); 
           alp = p75*ws*mu0n*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu0n*mu0n));
           //gam = gamma(ws,mu0n,gs,lm);
           gam = p5*ws*((c1 + c3*gs*(c1-ws)*mu0n*mu0n) / (c1-lm*lm*mu0n*mu0n));
           apg = alp + gam;
           amg = alp - gam;
           rdir = amg*(tdif_a*trnlay - c1) + apg*rdif_a;
           tdir = apg*tdif_a + (amg*rdif_a - (apg-c1))*trnlay;
   
           // recalculate rdif,tdif using direct angular integration over rdir,tdir,
           // since Delta-Eddington rdif formula is not well-behaved (it is usually
           // biased low and can even be negative); use ngmax angles and gaussian
           // integration for most accuracy:
           swt = c0;
           smr = c0;
           smt = c0;
           for (ng=0;ng<8;ng++) 
           {

             mu  = gauspt[ng];
             gwt = gauswt[ng];
             rdr = (c1 - lm*lm*mu*mu);

             if (abs(rdr) > p001) 
             {
             
             swt = swt + mu*gwt;
             //trn = max(exp_min, exp(-ts/(mu*dx_exp)))
             trn = exp(-ts/mu);
             
             //alp = alpha(ws,mu,gs,lm);
             alp = p75*ws*mu*((c1 + gs*(c1-ws))/rdr);
             //gam = gamma(ws,mu,gs,lm);
             gam = p5*ws*((c1 + c3*gs*(c1-ws)*mu*mu) /rdr);
             apg = alp + gam;
             amg = alp - gam;
             rdr = amg*(tdif_a*trn-c1) + apg*rdif_a;
             tdr = apg*tdif_a + (amg*rdif_a-(apg-c1))*trn;
             smr = smr + mu*rdr*gwt;
             smt = smt + mu*tdr*gwt;
             }
           }      // ng
           rdif_a = smr/swt;
           tdif_a = smt/swt;
   
           // homogeneous layer
           rdif_b = rdif_a;
           tdif_b = tdif_a;
   
           // add fresnel layer to top of desired layer if either 
           // air or snow overlies ice; we ignore refraction in ice 
           // if a melt pond overlies it:
   
           if( k == kfrsnl )
           {
             // compute fresnel reflection and transmission amplitudes
             // for two polarizations: 1=perpendicular and 2=parallel to
             // the plane containing incident, reflected and refracted rays.
             R1 = (mu0 - refindx*mu0n) / (mu0 + refindx*mu0n);
             R2 = (refindx*mu0 - mu0n) / (refindx*mu0 + mu0n);
             T1 = c2*mu0 / (mu0 + refindx*mu0n);
             T2 = c2*mu0 / (refindx*mu0 + mu0n);
   
             // unpolarized light for direct beam
             Rf_dir_a = p5 * (R1*R1 + R2*R2);
             Tf_dir_a = p5 * (T1*T1 + T2*T2)*refindx*mu0n/mu0;
   
             // precalculated diffuse reflectivities and transmissivities
             // for incident radiation above and below fresnel layer, using
             // the direct albedos and accounting for complete internal
             // reflection from below; precalculated because high order
             // number of gaussian points (~256) is required for convergence:
   
             // above
             Rf_dif_a = cp063;
             Tf_dif_a = c1 - Rf_dif_a;
             // below
             Rf_dif_b = cp455;
             Tf_dif_b = c1 - Rf_dif_b;
   
             // the k = kfrsnl layer properties are updated to combined 
             // the fresnel (refractive) layer, always taken to be above
             // the present layer k (i.e. be the top interface):
   
             rintfc   = c1 / (c1-Rf_dif_b*rdif_a);
             tdir   = Tf_dir_a*tdir + 
                                  Tf_dir_a*rdir * 
                                  Rf_dif_b*rintfc*tdif_a;
             rdir   = Rf_dir_a + 
                                  Tf_dir_a*rdir * 
                                  rintfc*Tf_dif_b;
             rdif_a = Rf_dif_a + 
                                  Tf_dif_a*rdif_a * 
                                  rintfc*Tf_dif_b;
             rdif_b = rdif_b + 
                                  tdif_b*Rf_dif_b * 
                                  rintfc*tdif_a;
             tdif_a = Tf_dif_a*rintfc*tdif_a;
             tdif_b = tdif_b*rintfc*Tf_dif_b;
   
             // update trnlay to include fresnel transmission
             trnlay = Tf_dir_a*trnlay;
           }      // k = kfrsnl
   
//             } else {
           
//               rdir = c0;  // layer reflectivity to direct radiation
//               rdif_a = c0; // layer reflectivity to diffuse radiation from above
//               rdif_b = c0;  // layer reflectivity to diffuse radiation from below
//               tdir = c0;  // layer transmission to direct radiation (solar beam + diffuse)
//               tdif_a = c0;  // layer transmission to diffuse radiation from above
//               tdif_b = c0;  // layer transmission to diffuse radiation from below
//               trnlay = c0;      // solar beam transm for layer (direct beam only)                             
   
//             }  // trntdr(k) > p001 

       g_f32_scr[kii] = rdir;  // layer reflectivity to direct radiation
       kii += blockDim.x;
       g_f32_scr[kii] = rdif_a; // layer reflectivity to diffuse radiation from above
       kii += blockDim.x;
       g_f32_scr[kii] = rdif_b;  // layer reflectivity to diffuse radiation from below
       kii += blockDim.x;
       g_f32_scr[kii] = tdir;  // layer transmission to direct radiation (solar beam + diffuse)
       kii += blockDim.x;
       g_f32_scr[kii] = tdif_a;  // layer transmission to diffuse radiation from above
       kii += blockDim.x;
       g_f32_scr[kii] = tdif_b;  // layer transmission to diffuse radiation from below
       kii += blockDim.x;
       g_f32_scr[kii] = trnlay;      // solar beam transm for layer (direct beam only)                             
       kii += blockDim.x;

       //debug
        //  if (k < 69) {
        //  g_f32_out[kiii] = rdir;
        //  kiii += blockDim.x;
        //  g_f32_out[kiii] = rdif_a;
        //  kiii += blockDim.x;
        //  }


       }       // k   end main level loop
 
       }
} // end sia2_edd_gpu_kernel_edd
 