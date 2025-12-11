// sia2 delta eddinton solution CUDA kernel
// ---------------------------------------------------------------------
// Ben Saenz 
// July 2010

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cuda.h>
#include <cutil.h>
#include <cutil_inline.h>


// Typedefs
//========================================================================
typedef unsigned long devPtr;


// Constants
//========================================================================
const int wavl = 31;

const float c0 = 0.0f;
const float c1 = 1.0f;
const float c1p5 = 1.5f;
const float c2 = 2.0f;
const float c3 = 3.0f;
const float c4 = 4.0f;
const float c10 = 10.0f;
const float p01 = 0.01f;
const float p001 = 0.002f;
const float p5 = 0.5f;
const float p75 = 0.75f;
const float p_08 = -0.08f;
const float refindx = 1.310f;  // refractive index of sea ice (used for water also)
const float cp063   = 0.063f;  // diffuse fresnel reflectivity from above
const float cp455   = 0.455f;  // diffuse fresnel reflectivity from below
const float aph_max = 0.0161f;

__device__  float inv_c_chl;
__device__  float inv_iced;
__device__  float a_ice_ir;
__device__  float inv_ad_denom;
__device__  float a_factor;
__device__  float scr_block;

__device__ float gauspt[8];  // gaussian angles (radians)
__device__ float gauswt[8];  // gaussian weights      
__device__ float aice[31];  // spectral ice absorption
__device__ float awater[31];  // spectral water absorption
__device__ float aph[31];  // spectral phytoplanton absorption
__device__ float quanta_to_watts[31];


// Kernel function include statements
//========================================================================

#include "sia2_edd_gpu_kernel_4_absorp.cu"
#include "sia2_edd_gpu_kernel_4_edd.cu"
#include "sia2_edd_gpu_kernel_4_aop.cu"
#include "sia2_edd_gpu_kernel_4_pur.cu"


// Function Declarations (as needed)
//========================================================================

extern "C" void sia2_edd_gpu_init_malloc_
    (devPtr      *g_i32_in,    
     devPtr      *g_f32_in,
     devPtr      *g_f32_out,
     devPtr      *g_f32_scr,
     int         *i32_in_size,
     int         *f32_in_size,
     int         *f32_out_size,
     int         *f32_scr_size,
     float       *h_inv_c_chl,    
     float       *h_inv_iced,     
     float       *h_a_ice_ir,     
     float       *h_inv_ad_denom,     
     float       *h_a_factor     
    );

extern "C" void sia2_edd_gpu_copy_step_mem_
    (int         *c_i32_in,
     float       *c_f32_in,         
     devPtr      *g_i32_in,    
     devPtr      *g_f32_in,
     int         *i32_in_size,
     int         *f32_in_size
    );

extern "C" void sia2_edd_gpu_solution_                                 
    (devPtr      *g_i32_in,    
     devPtr      *g_f32_in,
     devPtr      *g_f32_scr,
     devPtr      *g_f32_out,
     int         *i32_in_p_th,
     int         *f32_in_p_th,
     int         *f32_scr_p_th,
     int         *f32_out_p_th,
     int         *th_PerBlock,
     int         *tot_threads,
     int         *zz,
     float       *lda_m
    );

extern "C" void sia2_edd_gpu_return_
    (devPtr      *g_f32_out,
     float       *c_f32_out,
     int         *f32_out_p_th,
     int         *th_PerBlock,
     int         *tot_threads);

extern "C" void sia2_edd_gpu_free_
    (devPtr      *g_i32_in,    
     devPtr      *g_f32_in,
     devPtr      *g_f32_out,
     devPtr      *g_f32_scr); 


// Main 
//========================================================================

/** Compile me with:
/usr/local/cuda/bin/nvcc -g -o cudatest -m64 -deviceemu -O0 -arch compute_10 -code compute_10 sia2_edd_gpu.cu -I"/Developer/GPU Computing/C/common/inc" -L/usr/local/cuda/lib -lcuda -lcudart -L"/Developer/GPU Computing/C/lib" -lcutil_x86_64 -L/usr/lib -lstdc++
**/
/**
int main () 
{
     float      *c_coszen = NULL;    
     int        *c_srftyp = NULL;
     float      *c_tau = NULL;
     float      *c_w0 = NULL;
     float      *c_g = NULL;         
     int        *c_nslyr = NULL; 
     int        *c_nilyr = NULL;
     float      *c_trndir = NULL;
     float      *c_trntdr = NULL;
     float      *c_trndif = NULL;
     float      *c_rupdir = NULL;
     float      *c_rupdif = NULL;
     float      *c_rdndif = NULL;     

     devPtr      g_coszen = 0;    
     devPtr        g_srftyp = 0;
     devPtr      g_tau = 0;
     devPtr      g_w0 = 0;
     devPtr      g_g = 0;         
     devPtr        g_nslyr = 0; 
     devPtr        g_nilyr = 0;
     devPtr      g_trndir = 0;
     devPtr      g_trntdr = 0;
     devPtr      g_trndif = 0;
     devPtr      g_rupdir = 0;
     devPtr      g_rupdif = 0;
     devPtr      g_rdndif = 0;     
     int         nblocks = max_cells; // the max number of cells for which we need to compute edd

    unsigned int big_size = (nblocks) * z1 * sizeof(float);
    unsigned int small_size = (nblocks) * sizeof(float);
    unsigned int small_size_int = (nblocks) * sizeof(int);
    
    int ida_n = 9;
    int kk = 5;
    
    int compute_this_many_cells = max_cells;
    
    int i,j;
    
     c_coszen = (float *)malloc(small_size);    
     c_srftyp = (int *)malloc(small_size_int);
     c_tau = (float *)malloc(big_size);
     c_w0 = (float *)malloc(big_size);
     c_g = (float *)malloc(big_size);         
     c_nslyr = (int *)malloc(small_size_int); 
     c_nilyr = (int *)malloc(small_size_int);
     c_trndir = (float *)malloc(big_size);
     c_trntdr = (float *)malloc(big_size);
     c_trndif = (float *)malloc(big_size);
     c_rupdir = (float *)malloc(big_size);
     c_rupdif = (float *)malloc(big_size);
     c_rdndif = (float *)malloc(big_size);         

     for (i=0;i<nblocks;i++) {
         for (j=0;j>z1;j++) {
             c_g[i*j+j] = 0.94f;
             c_w0[i*j+j] = 0.996f;
             c_tau[i*j+j] = 8.66f;
         }
         c_coszen[i] = -0.2f;
         c_srftyp[i] = 1;
         c_nslyr[i] = 26; 
         c_nilyr[i] = 30;
     }
     
     sia2_edd_gpu_init_malloc_(&g_coszen,&g_srftyp,&g_tau, 
                  &g_w0,&g_g,&g_nslyr,&g_nilyr,&g_trndir,&g_trntdr, 
                  &g_trndif, &g_rupdir, &g_rupdif, &g_rdndif, &nblocks);

              printf("Calling CUDA Delta-Eddington solution ...\n");

              sia2_edd_gpu_solution_(&g_coszen,&g_srftyp,&g_tau,&g_w0,&g_g,        
                  &g_nslyr,&g_nilyr,c_trndir,c_trntdr,c_trndif,c_rupdir,c_rupdif, 
                  c_rdndif,&g_trndir,&g_trntdr,&g_trndif,&g_rupdir,&g_rupdif,&g_rdndif,    
                  &compute_this_many_cells,&kk,&ida_n,c_g);

              printf("Done: CUDA Delta-Eddington solution\n");
                  
     sia2_edd_gpu_free_(&g_coszen,&g_srftyp,&g_tau, 
                  &g_w0,&g_g,&g_nslyr,&g_nilyr,&g_trndir,&g_trntdr, 
                  &g_trndif, &g_rupdir, &g_rupdif, &g_rdndif);

     
                  
}
**/

// Function Code 
//========================================================================

void sia2_edd_gpu_init_malloc_
    (devPtr      *g_i32_in,    
     devPtr      *g_f32_in,
     devPtr      *g_f32_out,
     devPtr      *g_f32_scr,
     int         *i32_in_size,
     int         *f32_in_size,
     int         *f32_out_size,
     int         *f32_scr_size,
     float       *h_inv_c_chl,    
     float       *h_inv_iced,     
     float       *h_a_ice_ir,     
     float       *h_inv_ad_denom,     
     float       *h_a_factor     
    )
{

	CUdevice dev;
	int major = 0, minor = 0;
	int deviceCount = 0;
	char deviceName[256];
    int nGpuArchCoresPerSM[] = { -1, 8, 32 };

    float h_gauspt[] =  // gaussian angles (radians)
         {.9894009f,  .9445750f, 
         .8656312f,  .7554044f, 
         .6178762f,  .4580168f, 
         .2816036f,  .0950125f};
    float h_gauswt[] =  // gaussian weights      
         {.0271525f,  .0622535f, 
         .0951585f,  .1246290f, 
         .1495960f,  .1691565f, 
         .1826034f,  .1894506f};

    float ai[] = {
        0.0885236842105263f,0.0772037974683544f,0.0679112500000000f,
        0.0608333333333333f,0.0537764705882353f,0.0448395348837209f,
        0.0419500000000000f,0.0417887640449438f,0.0431516483516484f,
        0.0456212765957447f,0.0482808510638298f,0.0525010309278351f,
        0.0548303030303030f,0.0603360000000000f,0.0677450980392157f,
        0.0710714285714286f,0.0739858490566038f,0.0788036697247706f,
        0.0889272727272727f,0.103463716814159f,0.121070175438597f,
        0.143871794871795f,0.174333333333333f,0.207666666666667f,
        0.240735537190083f,0.276822580645161f,0.315078740157480f,
        0.353062500000000f,0.386656250000000f,0.444000000000000f,
        0.520338345864662f};
        
    float aw[] = {
        0.006630f,0.004730f,0.004540f,0.004950f,0.006350f,0.009220f,
        0.009790f,0.01060f,0.01270f,0.01500f,0.02040f,0.03250f,
        0.04090f,0.04340f,0.04740f,0.05650f,0.06190f,0.06950f,
        0.08960f,0.1351f,0.2224f,0.2644f,0.2755f,0.2916f,0.3108f,
        0.3400f,0.4100f,0.4390f,0.4650f,0.5160f,0.6240f};

  float h_aph[] = {
        0.0135f,0.0142f,0.015f,0.0158f,0.016f,0.0149f,0.0139f,0.0131f,
        0.0124f,0.0116f,0.0107f,0.0097f,0.0086f,0.0075f,0.0065f,0.0054f,
        0.0045f,0.0036f,0.003f,0.0028f,0.0027f,0.0031f,0.0035f,0.0039f,
        0.0036f,0.0034f,0.0062f,0.0097f,0.0097f,0.0061f,0.0005f};

  float qtw[] = {
        0.273600f, 0.267270f, 0.261080f, 0.255030f, 0.249120f, 0.243350f, 
        0.237720f, 0.232230f, 0.226880f, 0.221670f, 0.216600f, 0.211670f, 
        0.206880f, 0.202230f, 0.197720f, 0.193350f, 0.189120f, 0.185030f, 
        0.181080f, 0.177270f, 0.173600f, 0.170070f, 0.166680f, 0.163430f, 
        0.160320f, 0.157350f, 0.154520f, 0.151830f, 0.149280f, 0.146870f, 
        0.144600f};

    
	// note your project will need to link with cuda.lib files on windows
	printf("CUDA Device Query (Driver API) statically linked version \n");

	CUresult err = cuInit(0);
    CU_SAFE_CALL_NO_SYNC(cuDeviceGetCount(&deviceCount));
	// This function call returns 0 if there are no CUDA capable devices.
	if (deviceCount == 0) {
        printf("There is no device supporting CUDA\n");
	}
    for (dev = 0; dev < deviceCount; ++dev) {

        // get number of SMs on this GPU
		CU_SAFE_CALL_NO_SYNC( cuDeviceComputeCapability(&major, &minor, dev) );

        if (dev == 0) {
			// This function call returns 9999 for both major & minor fields, if no CUDA capable devices are present
            if (major == 9999 && minor == 9999)
                printf("There is no device supporting CUDA.\n");
            else if (deviceCount == 1)
                printf("There is 1 device supporting CUDA\n");
            else
                printf("There are %d devices supporting CUDA\n", deviceCount);
        }
		CU_SAFE_CALL_NO_SYNC( cuDeviceGetName(deviceName, 256, dev) );
        printf("\nDevice %d: \"%s\"\n", dev, deviceName);
//		size_t totalGlobalMem;
//		CU_SAFE_CALL_NO_SYNC( cuDeviceTotalMem(&totalGlobalMem, dev) );
//        printf("  Total amount of global memory:                 %u bytes\n", (unsigned long long) totalGlobalMem);
 	    int multiProcessorCount;
		CU_SAFE_CALL_NO_SYNC( cuDeviceGetAttribute( &multiProcessorCount, CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, dev ) );
        printf("  Number of multiprocessors:                     %d\n", multiProcessorCount);
        printf("  Number of cores:                               %d\n", nGpuArchCoresPerSM[major] * multiProcessorCount);
 	    int maxThreadsPerBlock;
		CU_SAFE_CALL_NO_SYNC( cuDeviceGetAttribute( &maxThreadsPerBlock, CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK, dev ) );
		printf("  Maximum number of threads per block:           %d\n",	maxThreadsPerBlock);
 	    int blockDim[3];
		CU_SAFE_CALL_NO_SYNC( cuDeviceGetAttribute( &blockDim[0], CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X, dev ) );
		CU_SAFE_CALL_NO_SYNC( cuDeviceGetAttribute( &blockDim[1], CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Y, dev ) );
		CU_SAFE_CALL_NO_SYNC( cuDeviceGetAttribute( &blockDim[2], CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_Z, dev ) );
        printf("  Maximum sizes of each dimension of a block:    %d x %d x %d\n", blockDim[0], blockDim[1], blockDim[2]);
 	    int gridDim[3];
		CU_SAFE_CALL_NO_SYNC( cuDeviceGetAttribute( &gridDim[0], CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_X, dev ) );
		CU_SAFE_CALL_NO_SYNC( cuDeviceGetAttribute( &gridDim[1], CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Y, dev ) );
		CU_SAFE_CALL_NO_SYNC( cuDeviceGetAttribute( &gridDim[2], CU_DEVICE_ATTRIBUTE_MAX_GRID_DIM_Z, dev ) );
        printf("  Maximum sizes of each dimension of a grid:     %d x %d x %d\n", gridDim[0], gridDim[1], gridDim[2]);
  	    int memPitch;
		CU_SAFE_CALL_NO_SYNC( cuDeviceGetAttribute( &memPitch, CU_DEVICE_ATTRIBUTE_MAX_PITCH, dev ) );
        printf("  Maximum memory pitch:                          %u bytes\n", memPitch);
  	    int textureAlign;
		CU_SAFE_CALL_NO_SYNC( cuDeviceGetAttribute( &textureAlign, CU_DEVICE_ATTRIBUTE_TEXTURE_ALIGNMENT, dev ) );
        printf("  Texture alignment:                             %u bytes\n", textureAlign);
  	    int clockRate;
		CU_SAFE_CALL_NO_SYNC( cuDeviceGetAttribute( &clockRate, CU_DEVICE_ATTRIBUTE_CLOCK_RATE, dev ) );
        printf("  Clock rate:                                    %.2f GHz\n", clockRate * 1e-6f);
    }

    //cudaSetDevice(cutGetMaxGflopsDeviceId());    // this call required cuda_inline_runtime.h
    cudaSetDevice(0);  // GTX 560 shows up as device 1, GTX 470 as device 0

//       unsigned int freem;
//       unsigned int totm;
//       CU_SAFE_CALL_NO_SYNC( cuMemGetInfo( &freem, &totm) );
//       printf("  Mem:                  %u free                  %u total\n", freem,totm);
    int sof = sizeof(float);


    unsigned int i32s = *i32_in_size * sizeof(int);
    unsigned int f32s = *f32_in_size * sof;
    unsigned int f32s_o = *f32_out_size * sof;
    unsigned int f32s_scr = *f32_scr_size * sof;

    // Allocate intput vectors in device memory
    printf("Allocating g_i32_in on GPU...");
    CUDA_SAFE_CALL( cudaMalloc((void**)g_i32_in, i32s) );
    printf("%u\n",*g_i32_in);
    printf("Allocating g_f32_in on GPU...");
    CUDA_SAFE_CALL( cudaMalloc((void**)g_f32_in, f32s) );
    printf("%u\n",*g_f32_in);
    printf("Allocating g_f32_out on GPU...");
    CUDA_SAFE_CALL( cudaMalloc((void**)g_f32_out, f32s_o) );
    printf("%u\n",*g_f32_out);
    printf("Allocating g_f32_scr on GPU...");
    CUDA_SAFE_CALL( cudaMalloc((void**)g_f32_scr, f32s_scr) );
    printf("%u\n",*g_f32_scr);

    // copy over constant arrays used by kernel
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(gauspt, h_gauspt, 8*sof,0,cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(gauswt, h_gauswt, 8*sof,0,cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(aice,ai,31*sof,0,cudaMemcpyHostToDevice) );  // spectral ice absorption
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(awater,aw,31*sof,0,cudaMemcpyHostToDevice) );  // spectral water absorption
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(aph,h_aph,31*sof,0,cudaMemcpyHostToDevice) );  // spectral phytoplanton absorption
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(quanta_to_watts,qtw,31*sof,0,cudaMemcpyHostToDevice) );  // spectral phytoplanton absorption

    // copy over constant floats used by kernels
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(inv_c_chl,h_inv_c_chl,sof,0,cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(inv_iced,h_inv_iced,sof,0,cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(a_ice_ir,h_a_ice_ir,sof,0,cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(inv_ad_denom,h_inv_ad_denom,sof,0,cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(a_factor,h_a_factor,sof,0,cudaMemcpyHostToDevice) );

}

void sia2_edd_gpu_copy_step_mem_
    (int         *c_i32_in,
     float       *c_f32_in,         
     devPtr      *g_i32_in,    
     devPtr      *g_f32_in,
     int         *wr_i32_size,
     int         *wr_f32_size
    )
{
    unsigned int i32s = *wr_i32_size * sizeof(int);
    unsigned int f32s = *wr_f32_size * sizeof(float);

    
    // copy over memory
    printf("CUDA Threads IN: %i INTs, %i FLOATS\n",*wr_i32_size,*wr_f32_size);
    CUDA_SAFE_CALL(cudaMemcpy((void*)*g_i32_in, c_i32_in, i32s, cudaMemcpyHostToDevice) );
    CUDA_SAFE_CALL(cudaMemcpy((void*)*g_f32_in, c_f32_in, f32s, cudaMemcpyHostToDevice) );

}

void sia2_edd_gpu_solution_                                 
    (devPtr      *g_i32_in,    
     devPtr      *g_f32_in,
     devPtr      *g_f32_scr,
     devPtr      *g_f32_out,
     int         *i32_in_p_th,
     int         *f32_in_p_th,
     int         *f32_scr_p_th,
     int         *f32_out_p_th,
     int         *th_PerBlock,
     int         *tot_threads,
     int         *zz,
     float       *lda_m
    )
{
    int wl;

    // find current of CUDA vectors
    int threadsPerBlock = *th_PerBlock;
    int blocksPerGrid =
        (*tot_threads + threadsPerBlock - 1) / threadsPerBlock;

    for (wl=0;wl<=wavl;wl++) 
    {

        // GPU calc layer absorption
        sia2_edd_gpu_kernel_absorp<<<blocksPerGrid, threadsPerBlock>>>
        (
            (const int*)      *g_i32_in, 
            (const float*)    *g_f32_in, 
            (float *)         *g_f32_scr, 
            (float *)         *g_f32_out, 
            *i32_in_p_th,
            *f32_in_p_th,
            *f32_scr_p_th,
            *f32_out_p_th,
            wl,
            *zz
         );
        
        cutilSafeCall( cudaThreadSynchronize() );

        // GPU calc delta-Eddington parameters
        sia2_edd_gpu_kernel_edd<<<blocksPerGrid, threadsPerBlock>>>
        (
            (const int*)       *g_i32_in, 
            (const float*)    *g_f32_in, 
            (float *)         *g_f32_scr, 
            (float *)         *g_f32_out, 
            *i32_in_p_th,
            *f32_in_p_th,
            *f32_scr_p_th,
            *f32_out_p_th,
            *zz,
            *lda_m
         );
        
        cutilSafeCall( cudaThreadSynchronize() );

        // GPU calc 4 of 6 (downstream) layer AOPs
        sia2_edd_gpu_kernel_aop<<<blocksPerGrid, threadsPerBlock>>>
        (
            (const int*)      *g_i32_in, 
            (float *)         *g_f32_scr, 
            (float *)         *g_f32_out, 
            *i32_in_p_th,
            *f32_scr_p_th,
            *f32_out_p_th,
            *zz
         );
        
        cutilSafeCall( cudaThreadSynchronize() );

        // GPU calc 2 of 6 (upstream) AOPs, along with final outputs PUR and Watts absorbed
        sia2_edd_gpu_kernel_pur<<<blocksPerGrid, threadsPerBlock>>>
        (
            (const int*)      *g_i32_in, 
            (const float*)    *g_f32_in, 
            (float *)         *g_f32_out, 
            (const float *)   *g_f32_scr, 
            *i32_in_p_th,
            *f32_in_p_th,
            *f32_out_p_th,
            *f32_scr_p_th,
            *zz,
            wl
         );
        
        cutilSafeCall( cudaThreadSynchronize() );

     }       
}

void sia2_edd_gpu_return_
    (devPtr      *g_f32_out,
     float       *c_f32_out,
     int         *f32_out_p_th,
     int         *th_PerBlock,
     int         *tot_threads)     // # floats
{
    // find current of CUDA vectors
    int threadsPerBlock = *th_PerBlock;
    int blocksPerGrid =
        (*tot_threads + threadsPerBlock - 1) / threadsPerBlock;
    unsigned int f32s_o = *f32_out_p_th * blocksPerGrid * 
        threadsPerBlock * sizeof(float);

    CUDA_SAFE_CALL(cudaMemcpy(c_f32_out, (void *)*g_f32_out, f32s_o, cudaMemcpyDeviceToHost) );
    printf("CUDA Threads OUT: %i BigArraySize: %u\n",*tot_threads,f32s_o);
}
   

void sia2_edd_gpu_free_
    (devPtr      *g_i32_in,    
     devPtr      *g_f32_in,
     devPtr      *g_f32_out,
     devPtr      *g_f32_scr)     
{
    CUDA_SAFE_CALL(cudaFree((void*)*g_i32_in) );
    CUDA_SAFE_CALL(cudaFree((void*)*g_f32_in) );
    CUDA_SAFE_CALL(cudaFree((void*)*g_f32_out) );
    CUDA_SAFE_CALL(cudaFree((void*)*g_f32_scr) );
}

