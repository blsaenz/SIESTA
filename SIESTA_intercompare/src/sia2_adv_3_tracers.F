! file: SIA2_advect_tracers.F
! Sea Ice Algae Model 2 - Saenz & Arrigo
! Version beta
! ======================================================================

     SUBROUTINE sia2_advect_tracers(f,m,ice,adv)

          use sia2_globals

          implicit none

          ! Function Arguments
          ! ------------------------------------------------------------
          type (forcing_type) :: f(tcells)
          type (meta_type) :: m(tcells)
          type (ice_type) :: ice(sda_n,ida_n,tcells)
          type (adv_type) :: adv(tcells)

          ! shared internal variables
          !----------------------------------------------------
          ! none

          ! private internal variables
          ! ------------------------------------------------------------
          integer :: ii,jj,iii,jjj,mi,mi1,ic,sc,i_ease,j_ease,int_z,sk_1,sk_z, &
              int_z_new,z_snow
          double precision :: tmp1,testvar,h_scale,a_scale,v_scale, &
              new_layer_bot,new_layer_top,old_layer_top,old_layer_bot, &
              interim_top,interim_bot,z1,z2,dz,dz_total,tracer1
          double precision, dimension(z_max) :: debug_z,tracer,th_old
          double precision, pointer :: temp3d(:,:)

#ifdef omp_on
!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ii,jj,iii,jjj,mi,mi1,ic,sc,i_ease,j_ease,int_z,sk_1,sk_z, &
!$OMP int_z_new,z_snow, &
!$OMP tmp1,testvar,h_scale,a_scale,v_scale, &
!$OMP new_layer_bot,new_layer_top,old_layer_top,old_layer_bot, &
!$OMP interim_top,interim_bot,z1,z2,dz,dz_total,tracer1, &
!$OMP debug_z,tracer,th_old, &
!$OMP temp3d)

!$OMP SECTIONS
!$OMP SECTION
#endif

          ! advect snow density/mass
          ! ----------------------------------------------------

#define s_tracer_var 1
#include "sia2_adv_3_3d_snow.inc.f90"

#ifdef omp_on
!$OMP SECTION
#endif

          ! advect snow heat
          ! ----------------------------------------------------

#undef s_tracer_var
#define s_tracer_var 2
#include "sia2_adv_3_3d_snow.inc.f90"



#ifdef omp_on
!$OMP SECTION
#endif

          ! advect melted/compacted snow indicator
          ! ----------------------------------------------------

#undef s_tracer_var
#define s_tracer_var 3
#include "sia2_adv_3_3d_snow.inc.f90"



#ifdef omp_on
!$OMP SECTION
#endif
          !print *,'Advecting S_layer...'
          ! advect S_layer
          ! ----------------------------------------------------
#define tracer_var 1
#include "sia2_adv_3_3d.inc.f90"


#ifdef omp_on
!$OMP SECTION
#endif


          !print *,'heat before: ',heat(147,81,:)

          ! advect heat
          ! ----------------------------------------------------
#undef tracer_var
#define tracer_var 3
#include "sia2_adv_3_3d.inc.f90"

          !print *,'heat after: ',heat(147,81,:)
#ifdef omp_on
!$OMP SECTION
#endif


          ! advect ice_d - this approximates conservation of mass
          ! ----------------------------------------------------
#undef tracer_var
#define tracer_var 2
#include "sia2_adv_3_3d.inc.f90"

#ifdef omp_on
!$OMP SECTION
#endif

          ! advect smalg
          ! ----------------------------------------------------
#undef tracer_var
#define tracer_var 4
#include "sia2_adv_3_3d.inc.f90"

#ifdef omp_on
!$OMP SECTION
#endif

          ! advect NO3
          ! ----------------------------------------------------
#undef tracer_var
#define tracer_var 5
#include "sia2_adv_3_3d.inc.f90"

#ifdef omp_on
!$OMP SECTION
#endif

          ! advect NH4
          ! ----------------------------------------------------
#undef tracer_var
#define tracer_var 6
#include "sia2_adv_3_3d.inc.f90"

#ifdef omp_on
!$OMP SECTION
#endif

          ! advect PO4
          ! ----------------------------------------------------
#undef tracer_var
#define tracer_var 7
#include "sia2_adv_3_3d.inc.f90"

#ifdef omp_on
!$OMP SECTION
#endif

          ! advect SiOH4
          ! ----------------------------------------------------
#undef tracer_var
#define tracer_var 8
#include "sia2_adv_3_3d.inc.f90"

#ifdef omp_on
!$OMP SECTION
#endif

          ! advect poc
          ! ----------------------------------------------------
#undef tracer_var
#define tracer_var 9
#include "sia2_adv_3_3d.inc.f90"

#ifdef omp_on
!$OMP SECTION
#endif

          ! advect age (2d tracer)
          ! ----------------------------------------------------
#define tracer1_var 1
#include "sia2_adv_3_2d.inc.f90"


#ifdef omp_on
!$OMP SECTION
#endif

          ! advect ridged fraction (2d tracer)
          ! ----------------------------------------------------
#undef tracer1_var
#define tracer1_var 2
#include "sia2_adv_3_2d.inc.f90"

          ! advect snow_dist fraction (2d tracer)
          ! ----------------------------------------------------
#undef tracer1_var
#define tracer1_var 3
#include "sia2_adv_3_2d.inc.f90"

          ! advect snow_rd fraction (2d tracer)
          ! ----------------------------------------------------
#undef tracer1_var
#define tracer1_var 4
#include "sia2_adv_3_2d.inc.f90"
          ! advect ml layer no3 (2d tracer)
          ! ----------------------------------------------------
#undef tracer1_var
#define tracer1_var 5
#include "sia2_adv_3_2d.inc.f90"
          ! advect ml layer nh4 (2d tracer)
          ! ----------------------------------------------------
#undef tracer1_var
#define tracer1_var 6
#include "sia2_adv_3_2d.inc.f90"
          ! advect ml layer po4 (2d tracer)
          ! ----------------------------------------------------
#undef tracer1_var
#define tracer1_var 7
#include "sia2_adv_3_2d.inc.f90"
          ! advect ml layer sioh3 (2d tracer)
          ! ----------------------------------------------------
#undef tracer1_var
#define tracer1_var 8
#include "sia2_adv_3_2d.inc.f90"



#ifdef omp_on
!$OMP END SECTIONS
!$OMP END PARALLEL
#endif

      end SUBROUTINE sia2_advect_tracers