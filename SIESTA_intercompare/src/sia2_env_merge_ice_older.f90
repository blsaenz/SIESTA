! SUBROUTINE: sia2_env_merge_ice
! ======================================================================
! merges ice conservatively from two categories

      subroutine sia2_env_merge_ice(sc_in,ic_in,sc_out,ic_out,mi,ice,pur,f,pl)

          use sia2_globals
          implicit none

          ! FUNCTION ARGUMENTS
		  !----------------------------------------------------
          integer :: &
              sc_in, &         ! snow distribution index of merge-in category
              ic_in, &         ! ice distribution index of merge-in category
              sc_out, &     ! snow distribution index of merge-out category
              ic_out, &     ! ice distribution index of merge-out category
              mi            ! model index of current ice cell

          type (ice_type) :: ice(sda_n,ida_n,tcells)
          integer (kind=2) :: pur(z_max,dt_per_day,lda_n,sda_n,ida_n,tcells)
          type (forcing_type) :: f(tcells)
          type (platelet_type) :: pl

          ! PRIVATE VARIABLES
		  !----------------------------------------------------
          !indices
          integer :: i,j,ii,jj,kk,jjj,iii,int_z,sk_1,sk_z,sc,ic,i_ease,j_ease

          ! sia2_env_brine_calc
          double precision :: t_melt,t_mean,s_mean,d_mean,bs_mean,bd_mean,bv_mean,heat_mean
          double precision :: alpha0,alpha1,alpha2,alpha3,tmp1,tmp2,tmp3,tmp4

          ! sia2_env_temp_from_heat
          ! double precision :: d_mean,heat_mean,t_mean,s_mean
          double precision ::  dz_mean,aq,bq,cq

          ! sia2_env_layer_thickness
          logical :: keep_growing,maxed_depth,ignore_h_max
          integer :: int_z_new,max_d_exceeded,z_snow_new
          double precision ::  r_depth,r_depth_last,t_depth,target_depth,layer_divisor,testvar
          double precision, dimension(z_max) :: th_new,id_new,debug_z,debug_start, &
              d_snow_new,heat_snow_new,th_snow_new

          ! sia2_env_adj_ice_boundaries
          ! integer :: ii,jj
          logical :: alg_fixed,remapping
          double precision :: interim_top,interim_bot,old_layer_top,new_layer_bot
          double precision :: old_layer_bot,z1,z2,dz,heat_total,dz_total,new_layer_top
          double precision :: flooded,t_flood,s_flood,d_flood,bs_flood,bd_flood
          double precision :: bv_flood,heat_flood,c_gl,c_gl_sal,s_gl,f_depth,snowh1m
          double precision :: dz_sk,sk_th_new,percent_converge,dhdt_cm_s,smalg_lost
		  double precision :: s_pond_new,h_pond_new,d_pond_new,smalg_pond_new
		  double precision :: poc_pond_new,no3_pond_new,nh4_pond_new,po4_pond_new
          double precision :: sioh4_pond_new,th_pond_new,snow_melt_now,melt_flood
          double precision :: mf_depth,melt_drain,c_gl_copy
          double precision :: sh_prev_sum,sh_offset_sum,age_new,ridged_new
          double precision, dimension(z_max) :: t_new,s_new,d_new,bs_new,bd_new
          double precision, dimension(z_max) :: bv_new,smalg_new,poc_new,smalg_pre_map
          double precision, dimension(z_max) :: melt_new,th_old,heat_new,d_new_snow
          double precision, dimension(z_max) :: heat_new_snow,th_new_snow
          double precision, dimension(z_max+1) :: NO3_new,NH4_new,PO4_new,SiOH4_new

          ! sia2_env_adj_snow_boundaries
          integer :: z_snow,z_last
          double precision :: airtemp1c,snowfall_d,c_snow

          ! redist local vars
          double precision :: af_total,vf_total

          ! temp erase me below
          double precision, dimension(z_max) :: t_prev,t_last,t_next,ki

          ! SUBROUTINE CODE
          ! ------------------------------------------------------------

          ! assign internal index to input category indices, so that #included
          ! code will work
          sc = sc_in
          ic = ic_in

          ! initializations
          snowh1m = 0.
          melt_drain = 0.
          max_d_exceeded = 0
          ignore_h_max = .false.

          !print *, 'Re-distribution from category',ic_out,'to category',ic

          ! redisribution of internal ice
          ! ------------------------------------------------------------
          int_z = max(ice(sc_out,ic_out,mi)%z-z_sk,ice(sc,ic,mi)%z-z_sk)
          th_new = 0.  ! vector assignment
          id_new = 0.  ! vector assignment

          af_total = ice(sc_out,ic_out,mi)%af + ice(sc,ic,mi)%af ! total area fraction we are combining

          do ii=1,int_z

              ! calculate fractional mutlipliers for both categories
              if (ii .le. ice(sc_out,ic_out,mi)%z-z_sk) then
                  z1 = ice(sc_out,ic_out,mi)%th(ii)*ice(sc_out,ic_out,mi)%af
              else
                  z1 = 0.
              endif

              z2 = 0.
              if (ice(sc,ic,mi)%af .gt. 0.) then
                  if (ii .le. ice(sc,ic,mi)%z-z_sk) then
                      z2 = ice(sc,ic,mi)%th(ii)*ice(sc,ic,mi)%af
                  endif
              endif

              ! find new bounds
              th_new(ii) = (z1+z2)/af_total
              if (ii .eq. 1) then
                  id_new(ii) = th_new(ii)
              else
                  id_new(ii) = id_new(ii-1) + th_new(ii)
              endif
              vf_total = af_total*th_new(ii)

              ! find new temp, salinity
              heat_mean = (ice(sc_out,ic_out,mi)%heat(ii)*z1 + ice(sc,ic,mi)%heat(ii)*z2) &
                       /vf_total
              d_mean = (ice(sc_out,ic_out,mi)%d(ii)*z1 + ice(sc,ic,mi)%d(ii)*z2) &
                       /vf_total
              s_mean = (ice(sc_out,ic_out,mi)%s(ii)*z1 + ice(sc,ic,mi)%s(ii)*z2) &
                       /vf_total
              dz_mean = 1.

#include "sia2_env_temp_from_heat.inc.f90"

              ! check to make sure we don't get warmer than the ocean
              !if (ice(sc_out,ic_out,mi)%id(ii) .le. ice(sc_out,ic_out,mi)%fbh) then
              !	  t_mean = min(t_mean,S_mean*mu)
              !else
              !	  t_mean = min(t_mean,f(mi)%t)
              !endif

              t_new(ii) = t_mean
              s_new(ii) = s_mean

              ! record tracers
              smalg_new(ii) = (ice(sc_out,ic_out,mi)%smalg(ii)*z1 &
                  + ice(sc,ic,mi)%smalg(ii)*z2)/vf_total
              poc_new(ii) = (ice(sc_out,ic_out,mi)%poc(ii)*z1 &
                  + ice(sc,ic,mi)%poc(ii)*z2)/vf_total
              no3_new(ii) = (ice(sc_out,ic_out,mi)%no3(ii)*z1 &
                  + ice(sc,ic,mi)%no3(ii)*z2)/vf_total
              nh4_new(ii) = (ice(sc_out,ic_out,mi)%nh4(ii)*z1 &
                  + ice(sc,ic,mi)%nh4(ii)*z2)/vf_total
              po4_new(ii) = (ice(sc_out,ic_out,mi)%po4(ii)*z1 &
                  + ice(sc,ic,mi)%po4(ii)*z2)/vf_total
              sioh4_new(ii) = (ice(sc_out,ic_out,mi)%sioh4(ii)*z1 &
                  + ice(sc,ic,mi)%sioh4(ii)*z2)/vf_total

          enddo

          ! redistribute skeletal layers
          ! ------------------------------------------------------------
          do ii=1,z_sk

              ! find corresponding skeletal layers
              jj = ice(sc_out,ic_out,mi)%z-z_sk+ii ! ic skeletal layer
              if (ice(sc,ic,mi)%af .gt. 0.) then
                  jjj = ice(sc,ic,mi)%z-z_sk+ii ! ic_n skeletal layer
              else
                  jjj=1 ! no ice in category ic
              endif
              kk = int_z+ii ! new skeletal layer

              ! calculate fractional mutlipliers for both categories
              z1 = ice(sc_out,ic_out,mi)%th(jj)*ice(sc_out,ic_out,mi)%af
              z2 = ice(sc,ic,mi)%th(jjj)*ice(sc,ic,mi)%af

              ! find new bounds
              th_new(kk) = (z1+z2)/af_total
              id_new(kk) = id_new(kk-1) + th_new(kk)
              vf_total = af_total*th_new(kk)

              ! record tracers
              smalg_new(kk) = (ice(sc_out,ic_out,mi)%smalg(jj)*z1 &
                  + ice(sc,ic,mi)%smalg(jjj)*z2)/vf_total
              poc_new(kk) = (ice(sc_out,ic_out,mi)%poc(jj)*z1 &
                  + ice(sc,ic,mi)%poc(jjj)*z2)/vf_total
              no3_new(kk) = (ice(sc_out,ic_out,mi)%no3(jj)*z1 &
                  + ice(sc,ic,mi)%no3(jjj)*z2)/vf_total
              nh4_new(kk) = (ice(sc_out,ic_out,mi)%nh4(jj)*z1 &
                  + ice(sc,ic,mi)%nh4(jjj)*z2)/vf_total
              po4_new(kk) = (ice(sc_out,ic_out,mi)%po4(jj)*z1 &
                  + ice(sc,ic,mi)%po4(jjj)*z2)/vf_total
              sioh4_new(kk) = (ice(sc_out,ic_out,mi)%sioh4(jj)*z1 &
                  + ice(sc,ic,mi)%sioh4(jjj)*z2)/vf_total

              ! keep skeletal physics constant
              s_new(kk) = f(mi)%s*c_5
              t_new(kk) = min(f(mi)%t,s_new(kk)*mu)

          enddo

          ! update surface temperature (not heat conservative yet)
          ! ------------------------------------------------------------
          if (ice(sc_out,ic_out,mi)%snow%z .eq. 0) then
              if (ice(sc,ic,mi)%snow%z .eq. 0) then
                  ! average ts if no snow - not heat conservative!!
                  ice(sc,ic,mi)%snow%ts = (ice(sc_out,ic_out,mi)%snow%ts*ice(sc_out,ic_out,mi)%af &
                  + ice(sc,ic,mi)%snow%ts*ice(sc,ic,mi)%af)/af_total
              else
                  ! set ts equal to ts of layer w/ snow
                  ice(sc,ic,mi)%snow%ts = ice(sc_out,ic_out,mi)%snow%ts
              endif
          elseif (ice(sc,ic,mi)%snow%z .gt. 0) then
              ! both have snow - take areal average (not heat conservative)!!
              ice(sc,ic,mi)%snow%ts = (ice(sc_out,ic_out,mi)%snow%ts*ice(sc_out,ic_out,mi)%af &
              + ice(sc,ic,mi)%snow%ts*ice(sc,ic,mi)%af)/af_total
          endif

          ! update ice physics, carry-over parameters, boundaries
          ! ------------------------------------------------------------
          ice(sc,ic,mi)%z = int_z + z_sk
          do ii=1,ice(sc,ic,mi)%z

              ! physics
              t_mean = t_new(ii)
              s_mean = s_new(ii)

#include "sia2_env_brine_calc.inc.f90"

              ice(sc,ic,mi)%t(ii) = t_mean
              ice(sc,ic,mi)%s(ii) = s_mean
              ice(sc,ic,mi)%d(ii) = d_mean
              ice(sc,ic,mi)%bs(ii) = bs_mean
              ice(sc,ic,mi)%bd(ii) = bd_mean
              ice(sc,ic,mi)%bv(ii) = bv_mean
              ice(sc,ic,mi)%heat(ii) = heat_mean

              ! update tracers
              ice(sc,ic,mi)%no3(ii) = NO3_new(ii)
              ice(sc,ic,mi)%nh4(ii) = NH4_new(ii)
              ice(sc,ic,mi)%po4(ii) = PO4_new(ii)
              ice(sc,ic,mi)%sioh4(ii) = SiOH4_new(ii)
              ice(sc,ic,mi)%poc(ii) = poc_new(ii)
              ice(sc,ic,mi)%smalg(ii) = smalg_new(ii)

              ! transfer carry-over paramaters
              if (ii .le. ice(sc_out,ic_out,mi)%z) then
                  z1 = ice(sc_out,ic_out,mi)%th(ii)*ice(sc_out,ic_out,mi)%af
              else
                  z1 = 0.
              endif
              if (ii .le. ice(sc,ic,mi)%z) then
                  z2 = ice(sc,ic,mi)%th(ii)*ice(sc,ic,mi)%af
              else
                  z2 = 0.
              endif
              vf_total = af_total*th_new(ii)

              ! vars with units x/pixel
              ice(sc,ic,mi)%prod(ii) = ice(sc,ic,mi)%prod(ii) + ice(sc_out,ic_out,mi)%prod(ii)

              ! vars with units x/m^2
              ice(sc,ic,mi)%fbv(ii) = (ice(sc,ic,mi)%fbv(ii)*ice(sc,ic,mi)%af &
                  + ice(sc_out,ic_out,mi)%fbv(ii)*ice(sc_out,ic_out,mi)%af)/af_total

              ! vars with units of x/m^3 and descriptive stats
!                      ice(sc,ic,mi)%dhdt_conv(ii) = (ice(sc,ic,mi)%dhdt_conv(ii)*z2 &
!                          + ice(sc_out,ic_out,mi)%dhdt_conv(ii)*z1)/vf_total
!                      ice(sc,ic,mi)%f0(ii) = (ice(sc,ic,mi)%f0(ii)*z2 &
!                          + ice(sc_out,ic_out,mi)%f0(ii)*z1)/vf_total
!                      do jj=1,dt_per_day
!                          do kk=1,lda_n
!                              pur(ii,jj,kk,ic,mi) = nint((pur(ii,jj,kk,ic,mi)*z2 &
!                                  + pur(ii,jj,kk,ic_out,mi)*z1)/vf_total)
!                          enddo
!                      enddo

              ! instantaneous/snapshot parameters that we will ignore
              !ice(sc,ic,mi)%bced(ii) =
              !ice(sc,ic,mi)%dsdt(ii) = 0.
              !ice(sc,ic,mi)%Ik1(ii) = 0.
              !ice(sc,ic,mi)%llim(ii) = 0.
              !ice(sc,ic,mi)%nlim(ii) = 0.
              !ice(sc,ic,mi)%plim(ii) = 0.
              !ice(sc,ic,mi)%silim(ii) = 0.
              !ice(sc,ic,mi)%slim(ii) = 0.
              !ice(sc,ic,mi)%dsdt3(ii) = 0.
              !ice(sc,ic,mi)%tgrad(ii) = 0.


              ! dimensions
              ice(sc,ic,mi)%th(ii) = th_new(ii)
              ice(sc,ic,mi)%id(ii) = id_new(ii)

          enddo

          ! update these values for remapping
          sk_z = int_z + z_sk
          sk_1 = int_z + 1

          r_depth = ice(sc,ic,mi)%id(int_z)

#include "sia2_env_ice_grid.inc.f90"

          ! start new layer top at old layer top - don't exclude any ice
          new_layer_top = -1.
          ! no growth/melt here
          s_gl = 0.
          c_gl = 0.
          ! no flooding here
          flooded = 0.
          melt_flood = 0.
          dhdt_cm_s = 0.   ! required in algal migration
          maxed_depth = .false.

#include "sia2_env_ice_remap.inc.f90"


          ! transfer pond paramaters
          ! --------------------------------------------------------

          if (use_ponds .eq. 1) then
              tmp1 = ice(sc,ic,mi)%af*ice(sc,ic,mi)%pond%th
              tmp2 = ice(sc_out,ic_out,mi)%af*ice(sc_out,ic_out,mi)%pond%th
              vf_total = tmp1 + tmp2

              if (vf_total .gt. 0.d0) then

                  ice(sc,ic,mi)%pond%s = (ice(sc,ic,mi)%pond%s*tmp1 &
                      + ice(sc_out,ic_out,mi)%pond%s*tmp2)/vf_total
                  ice(sc,ic,mi)%pond%d = (ice(sc,ic,mi)%pond%d*tmp1 &
                      + ice(sc_out,ic_out,mi)%pond%d*tmp2)/vf_total
                  ice(sc,ic,mi)%pond%heat = (ice(sc,ic,mi)%pond%heat*tmp1 &
                      + ice(sc_out,ic_out,mi)%pond%heat*tmp2)/vf_total
                  ice(sc,ic,mi)%pond%smalg = (ice(sc,ic,mi)%pond%smalg*tmp1 &
                      + ice(sc_out,ic_out,mi)%pond%smalg*tmp2)/vf_total     ! brine-based algal concentration (mgC/m^3)
                  ice(sc,ic,mi)%pond%poc = (ice(sc,ic,mi)%pond%poc*tmp1 &
                      + ice(sc_out,ic_out,mi)%pond%poc*tmp2)/vf_total     ! particulate organic carbon (detritus) (mgC/m^3)
                  ice(sc,ic,mi)%pond%no3 = (ice(sc,ic,mi)%pond%no3*tmp1 &
                      + ice(sc_out,ic_out,mi)%pond%no3*tmp2)/vf_total        ! ice layer brine NO3 concentration (µMol)
                  ice(sc,ic,mi)%pond%nh4 = (ice(sc,ic,mi)%pond%nh4*tmp1 &
                      + ice(sc_out,ic_out,mi)%pond%nh4*tmp2)/vf_total      ! ice layer brine NH4 concentration (µMol)
                  ice(sc,ic,mi)%pond%po4 = (ice(sc,ic,mi)%pond%po4*tmp1 &
                      + ice(sc_out,ic_out,mi)%pond%po4*tmp2)/vf_total     ! ice layer brine PO4 concentration (µMol)
                  ice(sc,ic,mi)%pond%sioh4 = (ice(sc,ic,mi)%pond%sioh4*tmp1 &
                      + ice(sc_out,ic_out,mi)%pond%sioh4*tmp2)/vf_total    ! ice layer brine SiOH4 concentration (µMol)

                  ice(sc,ic,mi)%pond%th = (tmp1+tmp2)/af_total
                  ice(sc,ic,mi)%pond%perc = (ice(sc,ic,mi)%af*ice(sc,ic,mi)%pond%perc + &
                    ice(sc_out,ic_out,mi)%af*ice(sc_out,ic_out,mi)%pond%perc)/af_total

              endif
          endif




          ! redistribute snow
          ! ------------------------------------------------------------

          ! calculate fractional mutlipliers for both categories
          z1 = ice(sc_out,ic_out,mi)%af/(ice(sc_out,ic_out,mi)%af + ice(sc,ic,mi)%af)
          z2 = ice(sc,ic,mi)%af/(ice(sc_out,ic_out,mi)%af + ice(sc,ic,mi)%af)

          ice(sc,ic,mi)%sh_prev = ice(sc_out,ic_out,mi)%sh_prev*z1 + ice(sc,ic,mi)%sh_prev*z2
          ice(sc,ic,mi)%sh_offset = ice(sc_out,ic_out,mi)%sh_offset*z1 + ice(sc,ic,mi)%sh_offset*z2


          jjj = max(ice(sc_out,ic_out,mi)%snow%z,ice(sc,ic,mi)%snow%z)
          ice(sc,ic,mi)%snow%depth = 0.
          do ii=1,jjj

              ! calculate fractional mutlipliers for both categories
              z1 = ice(sc_out,ic_out,mi)%snow%th(ii)*ice(sc_out,ic_out,mi)%af
              z2 = ice(sc,ic,mi)%snow%th(ii)*ice(sc,ic,mi)%af

              ! find new bounds
              th_new(ii) = (z1+z2)/af_total
              if (ii .eq. 1) then
                  id_new(ii) = th_new(ii)
              else
                  id_new(ii) = id_new(ii-1) + th_new(ii)
              endif
              vf_total = af_total*th_new(ii)

              ! find new temp
              heat_total = (ice(sc_out,ic_out,mi)%snow%heat(ii)*z1 + &
                  ice(sc,ic,mi)%snow%heat(ii)*z2)/vf_total
              d_mean = (ice(sc_out,ic_out,mi)%snow%d(ii)*z1 + &
                  ice(sc,ic,mi)%snow%d(ii)*z2)/vf_total
              dz_total = 1.

              ! Solve quadradic for new snow temp
              heat_total = heat_total/d_mean
              t_mean = (0.2309-sqrt(0.05331481+0.0136* &
                  (heat_snow0-heat_total)))/(-0.0068)
              t_mean = min(t_mean,kelvin0)

#include "sia2_env_snow_calc.inc.f90"

              ! update all snow parameters
              ice(sc,ic,mi)%snow%th(ii) = th_new(ii)
              ice(sc,ic,mi)%snow%t(ii) = t_mean - kelvin0
              ice(sc,ic,mi)%snow%d(ii) = d_mean
              ice(sc,ic,mi)%snow%heat(ii) = heat_mean
              ice(sc,ic,mi)%snow%melt(ii) = (ice(sc_out,ic_out,mi)%snow%melt(ii)*z1 + &
                  ice(sc,ic,mi)%snow%melt(ii)*z2)/vf_total

              ! find new total snow depth
              ice(sc,ic,mi)%snow%depth = ice(sc,ic,mi)%snow%depth + th_new(ii)

          enddo

          ! record final number snow layers
          ice(sc,ic,mi)%snow%z = jjj

          ! setup for snow regridding
          r_depth = ice(sc,ic,mi)%snow%depth
          ice(sc,ic,mi)%sh_prev = r_depth

#include "sia2_env_snow_grid.inc.f90"

          flooded = 0.
          melt_drain = 0.

#include "sia2_env_snow_remap.inc.f90"

          tmp1 = ice(sc,ic,mi)%af + ice(sc_out,ic_out,mi)%af

          ! update age and ridged tracers
          ! --------------------------------------------------------
          ice(sc,ic,mi)%age = (ice(sc,ic,mi)%age*ice(sc,ic,mi)%af + &
            ice(sc_out,ic_out,mi)%age*ice(sc_out,ic_out,mi)%af)/tmp1
          ice(sc,ic,mi)%ridged = (ice(sc,ic,mi)%ridged*ice(sc,ic,mi)%af + &
            ice(sc_out,ic_out,mi)%ridged*ice(sc_out,ic_out,mi)%af)/tmp1

					! update ml layer tracers
          ! --------------------------------------------------------
					ice(sc,ic,mi)%no3(ml_z) = &
						(ice(sc_out,ic_out,mi)%no3(ml_z)*ice(sc_out,ic_out,mi)%af + &
						ice(sc,ic,mi)%nh4(ml_z)*ice(sc,ic,mi)%af)/tmp1
					ice(sc,ic,mi)%nh4(ml_z) = &
						(ice(sc_out,ic_out,mi)%nh4(ml_z)*ice(sc_out,ic_out,mi)%af + &
						ice(sc,ic,mi)%nh4(ml_z)*ice(sc,ic,mi)%af)/tmp1
					ice(sc,ic,mi)%po4(ml_z) = &
						(ice(sc_out,ic_out,mi)%po4(ml_z)*ice(sc_out,ic_out,mi)%af + &
						ice(sc,ic,mi)%po4(ml_z)*ice(sc,ic,mi)%af)/tmp1
					ice(sc,ic,mi)%sioh4(ml_z) = &
						(ice(sc_out,ic_out,mi)%sioh4(ml_z)*ice(sc_out,ic_out,mi)%af + &
						ice(sc,ic,mi)%sioh4(ml_z)*ice(sc,ic,mi)%af)/tmp1
!					ice(sc,ic,mi)%poc(ml_z) = &
!						(ice(sc_out,ic_out,mi)%poc(ml_z)*ice(sc_out,ic_out,mi)%af + &
!						ice(sc,ic,mi)%poc(ml_z)*ice(sc,ic,mi)%af)/tmp1

          ! update af & nullify old category
          ! --------------------------------------------------------
          ice(sc,ic,mi)%af = ice(sc,ic,mi)%af + ice(sc_out,ic_out,mi)%af
          ic = ic_out
          sc = sc_out

#include "sia2_env_null_params.f90"
           pur(:,:,:,sc,ic,mi) = pur_0


      end subroutine sia2_env_merge_ice
