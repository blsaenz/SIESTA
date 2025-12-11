
#include <cuda.h>
#include "math_functions.h"

// example of decaling constants/globals in CUDA device memory
// static __device__ __constant__ __TOptionData d_OptionData[MAX_OPTIONS];

// example of setting constant/global in CUDA device memory, i think
// cutilSafeCall( cudaMemcpyToSymbol(device_constant,host_constant, size) );

#define z1 69  // define the max num layers

// Delta-Eddington solution expressions
// ======================================================================
/**
#define alpha(w,uu,gg,e) p75*w*uu*((c1 + gg*(c1-w))/(c1 - e*e*uu*uu))
#define gamma(w,uu,gg,e) p5*w*((c1 + c3*gg*(c1-w)*uu*uu) / (c1-e*e*uu*uu))
#define n_macro(uu,et)   ((uu+c1)*(uu+c1)/et ) - ((uu-c1)*(uu-c1)*et)
#define u_macro(w,gg,e)  c1p5*(c1 - w*gg)/e
#define el(w,gg)         sqrt(c3*(c1-w)*(c1 - w*gg))
#define taus(w,f,t)      (c1 - w*f)*t
#define omgs(w,f)        (c1 - f)*w/(c1 - w*f)
#define asys(gg,f)       (gg - f)/(c1 - f)
**/

// Constants
// =======================================================================
//const float nc1 = -1.0;
const float c0 = 0.0f;
const float c1 = 1.0f;
const float c1p5 = 1.5f;
const float c2 = 2.0f;
const float c3 = 3.0f;
const float c4 = 4.0f;
//const float c5 = 5.0f;
//const float c10 = 10.0f;

const float p01 = 0.01f;
const float p5 = 0.5f;
const float p75 = 0.75f;
//const float p99 = 0.99f;    

// refractive index for sea ice, water; pre-computed, band-independent,
// diffuse fresnel reflectivities
const float refindx = 1.310f;  // refractive index of sea ice (used for water also)
const float cp063   = 0.063f;  // diffuse fresnel reflectivity from above
const float cp455   = 0.455f;  // diffuse fresnel reflectivity from below

//const float gauspt[] =  // gaussian angles (radians)
//             {.9894009,  .9445750, 
//             .8656312,  .7554044, 
//             .6178762,  .4580168, 
//             .2816036,  .0950125};
//const float gauswt[] =  // gaussian weights      
//             {.0271525,  .0622535, 
//             .0951585,  .1246290, 
//             .1495960,  .1691565, 
//             .1826034,  .1894506};
//const float logm9[] = 
//             {0.102, 0.272, 0.427, 0.532, 0.721,
//              0.952, 1.31, 1.74, 3.31};
//const float logm6[] = 
//             {0.145, 0.385, 0.585, 0.860, 1.345, 2.70};



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

__global__ void sia2_edd_gpu_kernel (
    const float      *g_coszen, 
    const int        *g_srftyp, 
    const float      *g_tau, 
    const float      *g_w0,
    const float      *g_g,
    const int        *g_nslyr,
    const int        *g_nilyr,
    float            *g_trndir,
    float            *g_trntdr,
    float            *g_trndif,
    float            *g_rupdir,
    float            *g_rupdif,
    float            *g_rdndif,
    int              kk, // the light category number
    int              lda_n)  // the total number of light categories
{

      int //(kind=int_kind) :: & 
         kfrsnl;      // radiation interface index for fresnel layer
 
      float //(kind=dbl_kind), dimension (0:klev,icells_DE) :: &
         tau[z1], 
         w0[z1],
         g[z1],
         trndir[z1],
         trntdr[z1],
         trndif[z1],
         rupdir[z1],
         rupdif[z1],
         rdndif[z1],
         rdir[z1]    ,  // layer reflectivity to direct radiation
         rdif_a[z1]  ,  // layer reflectivity to diffuse radiation from above
         rdif_b[z1]  ,  // layer reflectivity to diffuse radiation from below
         tdir[z1]    ,  // layer transmission to direct radiation (solar beam + diffuse)
         tdif_a[z1]  ,  // layer transmission to diffuse radiation from above
         tdif_b[z1]  ,  // layer transmission to diffuse radiation from below
         trnlay[z1];      // solar beam transm for layer (direct beam only)
 
      int //(kind=int_kind) :: & 
         i,
         k,           // level index
         k1,
         ki,
         klev,
         klevp,
         srftyp,
         nslyr,
         nilyr;
 
      float //(kind=dbl_kind), parameter :: &
         trmin = 0.001f;   // minimum total transmission allowed

      float //(kind=dbl_kind) :: &
         tautot   , // layer optical depth
         wtot     , // layer single scattering albedo
         gtot     , // layer asymmetry parameter
         ftot     , // layer forward scattering fraction
         ts       , // layer scaled extinction optical depth
         ws       , // layer scaled single scattering albedo
         gs       , // layer scaled asymmetry parameter
         rintfc   , // reflection (multiple) at an interface
         refkp1   , // interface multiple scattering for k+1
         refkm1   , // interface multiple scattering for k-1
         tdrrdir  , // direct tran times layer direct ref 
         tdndif   , // total down diffuse = tot tran - direct tran
         coszen;    
 
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
//         alpha    ,  // term in direct reflectivity and transmissivity
//         gamma    ,  // term in direct reflectivity and transmissivity
//         el       ,  // term in alpha,gamma,n,u
//         taus     ,  // scaled extinction optical depth
//         omgs     ,  // scaled single particle scattering albedo
//         asys     ,  // scaled asymmetry parameter
//         u        ,  // term in diffuse reflectivity and transmissivity
//         n        ,  // term in diffuse reflectivity and transmissivity
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
 
      int  //(kind=int_kind), parameter :: &
         ngmax = 8;    // number of gaussian angles in hemisphere
 
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
         lda_m;          

      float //(kind=dbl_kind) :: &
         albodr = c0;  // direct reflectivity of ocean - zero under ice
      float //(kind=dbl_kind) :: &
         albodf = c0;  // diffuse reflectivity of ocean - zero under ice

      // constant arrays that won't compile as such, so they are placed here
      float gauspt[] =  // gaussian angles (radians)
             {.9894009f,  .9445750f, 
             .8656312f,  .7554044f, 
             .6178762f,  .4580168f, 
             .2816036f,  .0950125f};
      float gauswt[] =  // gaussian weights      
             {.0271525f,  .0622535f, 
             .0951585f,  .1246290f, 
             .1495960f,  .1691565f, 
             .1826034f,  .1894506f};
      float logm9[] = 
             {0.102f, 0.272f, 0.427f, 0.532f, 0.721f,
              0.952f, 1.31f, 1.74f, 3.31f};
      float logm6[] = {0.145f, 0.385f, 0.585f, 0.860f, 1.345f, 2.70f};


      // find CUDA indices to cell data
      i = blockDim.x * blockIdx.x + threadIdx.x;

      if (g_nilyr[i] != c0) {
    
    
          k1 = i*z1;
          
          // Copy slow global memory into local mem for processing
          coszen = g_coszen[i];
          srftyp = g_srftyp[i];
          nslyr = g_nslyr[i];
          nilyr = g_nilyr[i];
          klev = nslyr + nilyr; // , &  // number of radiation layers - 1
          klevp = klev  + 1; //              // number of radiation interfaces - 1
          for (k=0;k<klevp;k++) 
          { 
             ki = k + k1;
             g[k] = g_g[ki];
             tau[k] = g_tau[ki];
             w0[k] = g_w0[ki];
          }         
     
          // revise tau based on snow/light distribution category
          if (lda_n == 6) 
          {
              lda_m = logm6[kk-1];
          } else if (lda_n == 9)  {
              lda_m = logm9[kk-1];
          } else {
              lda_m = c1;
          }          
     
          for (k=0;k<klevp;k++) 
          { 
    
             // initialize all layer apparent optical properties to 0
             rdir[k]   = c0;
             rdif_a[k] = c0;
             rdif_b[k] = c0;
             tdir[k]   = c0;
             tdif_a[k] = c0;
             tdif_b[k] = c0;
             trnlay[k] = c0;
    
    
             // initialize all output to 0
             trndir[k] = c0;
             trntdr[k] = c0;
             trndif[k] = c0;
             rupdir[k] = c0;
             rupdif[k] = c0;
             rdndif[k] = c0;
          }

          // initialize top interface of top layer 
          trndir[0] =   c1;
          trntdr[0] =   c1;
          trndif[0] =   c1;
          rdndif[0] =   c0;
    
    
           // compute level of fresnel refraction
           if( srftyp < 2 )
           {
             // if snow over sea ice or bare sea ice, fresnel level is
             // at base of sea ice SSL (and top of the sea ice DL); the
             // snow SSL counts for one, then the number of snow layers,
             // then the sea ice SSL which also counts for one:
             kfrsnl = nslyr + 1;
           } else {
             // if ponded sea ice, fresnel level is the top of the pond 
             kfrsnl = 0;
           }
       
     
          // proceed down one layer at a time; if the total transmission to
          // the interface just above a given layer is less than trmin, then no
          // Delta-Eddington computation for that layer is done.
    
            // begin main level loop
            for (k=0;k<klev;k++) 
            { 
       
            // initialize current layer properties to zero; only if total
            // transmission to the top interface of the current layer exceeds the
            // minimum, will these values be computed below:
            if ( k > 0 ) 
            {
                // Calculate the solar beam transmission, total transmission, and
                // reflectivity for diffuse radiation from below at interface k, 
                // the top of the current layer k:
                //
                //              layers       interface
                //         
                //       ---------------------  k-1 
                //                k-1
                //       ---------------------  k
                //                 k
                //       ---------------------  
       
                  trndir[k] = trndir[k-1]*trnlay[k-1];
                  refkm1        = c1/(c1 - rdndif[k-1]*rdif_a[k-1]);
                  tdrrdir       = trndir[k-1]*rdir[k-1];
                  tdndif        = trntdr[k-1] - trndir[k-1];
                  trntdr[k] = trndir[k-1]*tdir[k-1] + 
                   (tdndif + tdrrdir*rdndif[k-1])*refkm1*tdif_a[k-1];
                  rdndif[k] = rdif_b[k-1] + 
                    (tdif_b[k-1]*rdndif[k-1]*refkm1*tdif_a[k-1]);
                  trndif[k] = trndif[k-1]*refkm1*tdif_a[k-1];
            }       // k > 0  
       
            // compute next layer Delta-eddington solution only if total transmission
            // of radiation to the interface just above the layer exceeds trmin.
              if (trntdr[k] > trmin ) 
              {
       
            // calculation over layers with penetrating radiation
               
               if (k < nslyr) 
               {
                   tautot  = tau[k] * lda_m;
               } else {
                   tautot  = tau[k];
               }
               wtot    = w0[k];
               gtot    = g[k];
               ftot    = gtot*gtot;

               ts   = (c1 - wtot*ftot)*tautot;
               ws   = (c1 - ftot)*wtot/(c1 - wtot*ftot);
               gs   = (gtot - ftot)/(c1 - ftot);
       
               //ts   = taus(wtot,ftot,tautot);
               //ws   = omgs(wtot,ftot);
               //gs   = asys(gtot,ftot);
               //lm   = el(ws,gs);
               lm = sqrt(c3*(c1-ws)*(c1 - ws*gs));
               //ue   = u_macro(ws,gs,lm);
               ue = c1p5*(c1 - ws*gs)/lm;
       
               // mu0 is cosine solar zenith angle above the fresnel level; make 
               // sure mu0 is large enough for stable and meaningful radiation
               // solution: .01 is like sun just touching horizon with its lower edge
       
               mu0  = fmax(coszen,p01);
       
               // mu0n is cosine solar zenith angle used to compute the layer
               // Delta-Eddington solution; it is initially computed to be the
               // value below the fresnel level, i.e. the cosine solar zenith 
               // angle below the fresnel level for the refracted solar beam:
       
               mu0n = sqrt(c1-((c1-mu0*mu0)/(refindx*refindx)));
       
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
       
               rdif_a[k] = (ue+c1)*(ue-c1)*(c1/extins - extins)/ne;
               tdif_a[k] = c4*ue/ne;
       
               // evaluate rdir,tdir for direct beam
               //trnlay(k) = max(exp_min, exp(-ts/(mu0n*dx_exp)))
               trnlay[k] = exp(-ts/mu0n);
               //alp = alpha(ws,mu0n,gs,lm); 
               alp = p75*ws*mu0n*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu0n*mu0n));
               //gam = gamma(ws,mu0n,gs,lm);
               gam = p5*ws*((c1 + c3*gs*(c1-ws)*mu0n*mu0n) / (c1-lm*lm*mu0n*mu0n));
               apg = alp + gam;
               amg = alp - gam;
               rdir[k] = amg*(tdif_a[k]*trnlay[k] - c1) + apg*rdif_a[k];
               tdir[k] = apg*tdif_a[k] + (amg*rdif_a[k] - (apg-c1))*trnlay[k];
       
               // recalculate rdif,tdif using direct angular integration over rdir,tdir,
               // since Delta-Eddington rdif formula is not well-behaved (it is usually
               // biased low and can even be negative); use ngmax angles and gaussian
               // integration for most accuracy:
               swt = c0;
               smr = c0;
               smt = c0;
               for (ng=0;ng<ngmax;ng++) 
               {
                 mu  = gauspt[ng];
                 gwt = gauswt[ng];
                 swt = swt + mu*gwt;
                 //trn = max(exp_min, exp(-ts/(mu*dx_exp)))
                 trn = exp(-ts/mu);
                 //alp = alpha(ws,mu,gs,lm);
                 alp = p75*ws*mu*((c1 + gs*(c1-ws))/(c1 - lm*lm*mu*mu));
                 //gam = gamma(ws,mu,gs,lm);
                 gam = p5*ws*((c1 + c3*gs*(c1-ws)*mu*mu) / (c1-lm*lm*mu*mu));
                 apg = alp + gam;
                 amg = alp - gam;
                 rdr = amg*(tdif_a[k]*trn-c1) + apg*rdif_a[k];
                 tdr = apg*tdif_a[k] + (amg*rdif_a[k]-(apg-c1))*trn;
                 smr = smr + mu*rdr*gwt;
                 smt = smt + mu*tdr*gwt;
               }      // ng
               rdif_a[k] = smr/swt;
               tdif_a[k] = smt/swt;
       
               // homogeneous layer
               rdif_b[k] = rdif_a[k];
               tdif_b[k] = tdif_a[k];
       
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
       
                 rintfc   = c1 / (c1-Rf_dif_b*rdif_a[kfrsnl]);
                 tdir[kfrsnl]   = Tf_dir_a*tdir[kfrsnl] + 
                                      Tf_dir_a*rdir[kfrsnl] * 
                                      Rf_dif_b*rintfc*tdif_a[kfrsnl];
                 rdir[kfrsnl]   = Rf_dir_a + 
                                      Tf_dir_a*rdir[kfrsnl] * 
                                      rintfc*Tf_dif_b;
                 rdif_a[kfrsnl] = Rf_dif_a + 
                                      Tf_dif_a*rdif_a[kfrsnl] * 
                                      rintfc*Tf_dif_b;
                 rdif_b[kfrsnl] = rdif_b[kfrsnl] + 
                                      tdif_b[kfrsnl]*Rf_dif_b * 
                                      rintfc*tdif_a[kfrsnl];
                 tdif_a[kfrsnl] = Tf_dif_a*rintfc*tdif_a[kfrsnl];
                 tdif_b[kfrsnl] = tdif_b[kfrsnl]*rintfc*Tf_dif_b;
       
                 // update trnlay to include fresnel transmission
                 trnlay[kfrsnl] = Tf_dir_a*trnlay[kfrsnl];
               }      // k = kfrsnl
       
               } // trntdr(k) > trmin
       
            }       // k   end main level loop
       
          // compute total direct beam transmission, total transmission, and
          // reflectivity for diffuse radiation (from below) for all layers
          // above the underlying ocean; note that we ignore refraction between 
          // sea ice and underlying ocean:
          //
          //       For k = klevp
          //
          //              layers       interface
          //
          //       ---------------------  k-1 
          //                k-1
          //       ---------------------  k
          //               ocean    
       
            k = klevp;
            trndir[k] = trndir[k-1]*trnlay[k-1];
            refkm1        = c1/(c1 - rdndif[k-1]*rdif_a[k-1]);
            tdrrdir       = trndir[k-1]*rdir[k-1];
            tdndif        = trntdr[k-1] - trndir[k-1];
            trntdr[k] = trndir[k-1]*tdir[k-1] + 
              (tdndif + tdrrdir*rdndif[k-1])*refkm1*tdif_a[k-1];
            rdndif[k] = rdif_b[k-1] + 
              (tdif_b[k-1]*rdndif[k-1]*refkm1*tdif_a[k-1]);
            trndif[k] = trndif[k-1]*refkm1*tdif_a[k-1];
       
          // compute reflectivity to direct and diffuse radiation for layers 
          // below by adding succesive layers starting from the underlying 
          // ocean and working upwards:
          //
          //              layers       interface
          //
          //       ---------------------  k
          //                 k
          //       ---------------------  k+1
          //                k+1
          //       ---------------------
       
            rupdir[klevp] = albodr;
            rupdif[klevp] = albodf;
    
          for (k=klev;k>=0;k--)
          {
              // interface scattering
              refkp1 = c1/( c1 - rdif_b[k]*rupdif[k+1]);
              // dir from top layer plus exp tran ref from lower layer, interface
              // scattered and tran thru top layer from below, plus diff tran ref
              // from lower layer with interface scattering tran thru top from below
              rupdir[k] = rdir[k] + 
                            ( trnlay[k]*rupdir[k+1] + 
                             (tdir[k]-trnlay[k])*rupdif[k+1] ) * 
                              refkp1*tdif_b[k];
              // dif from top layer from above, plus dif tran upwards reflected and
              // interface scattered which tran top from below
              rupdif[k] = rdif_a[k] + 
                              tdif_a[k]*rupdif[k+1]* 
                              refkp1*tdif_b[k];
          }       // k

          for (k=0;k<klevp;k++) 
          { 
             ki = k + k1;
             g_trndir[ki] = trndir[k];
             g_trntdr[ki] = trntdr[k];
             g_trndif[ki] = trndif[k];
             g_rupdir[ki] = rupdir[k];
             g_rupdif[ki] = rupdif[k];
             g_rdndif[ki] = rdndif[k];
          }         

      } // end of 'is this a valid cell' test (g != 0.0)

} // end sia2_edd_gpu_kernel
 