! program to test CUDA C linking to ifort fortran

      program cudaifort
 
      use IFCORE
      implicit none
      
! Declarations
! ======================================================================

      integer :: &
          g_coszen, &    
          g_tau, &
          g_w0, &
          g_g, &         
          g_trndir, &
          g_trntdr, &
          g_trndif, &
          g_rupdir, &
          g_rupdif, &
          g_rdndif, &     
          g_nslyr, &
          g_nilyr, &
          g_srftyp
          
      integer :: &
          tcells, &
          fpe_old, &
          fpe_new
      
      tcells = 189959

      fpe_new = 0
      fpe_old = FOR_SET_FPE(fpe_new)      
      print *,'Set FPE Trap from: ',fpe_old,' to: ', fpe_new

      call sia2_edd_gpu_init_malloc(g_coszen,g_srftyp,g_tau, &
      g_w0,g_g,g_nslyr,g_nilyr,g_trndir,g_trntdr, &
      g_trndif, g_rupdir, g_rupdif, g_rdndif, tcells)

      fpe_new = 983054
      fpe_old = FOR_SET_FPE(fpe_new)      
      print *,'Set FPE Trap from: ',fpe_old,' to: ', fpe_new

      call sia2_edd_gpu_free(g_coszen,g_srftyp,g_tau, &
      g_w0,g_g,g_nslyr,g_nilyr,g_trndir,g_trntdr, & 
      g_trndif, g_rupdir, g_rupdif, g_rdndif)

      end program cudaifort