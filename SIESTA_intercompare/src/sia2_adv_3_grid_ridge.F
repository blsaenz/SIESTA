! file: SIA2_advect.F90
! Sea Ice Algae Model 2 - Saenz & Arrigo
! Version beta
! ======================================================================

      SUBROUTINE sia2_adv_grid_ridge(pj,ease_f,f,m,ice,adv)

#ifdef omp_on
          use omp_lib
#endif
          use sia2_globals

          implicit none

          ! Function Arguments
          ! ------------------------------------------------------------
          type (proj_type) :: pj
          type (ease_f_type) :: ease_f
          type (forcing_type) :: f(tcells)
          type (meta_type) :: m(tcells)
          type (ice_type) :: ice(sda_n,ida_n,tcells)
          type (adv_type) :: adv(tcells)

          ! platelet crap
          type (platelet_type) :: pl

          ! shared internal variables
          !----------------------------------------------------
          integer :: max_d_exceeded,chunk_size,convg_maxed
          double precision :: apn_denom
          type(ice_type),pointer :: ice_ridged(:,:,:),ice_temp(:),ice_temp2

          ! private internal variables
          ! ------------------------------------------------------------
          integer :: i,j,ii,jj,kk,jjj,iii,i_ease,j_ease,chunk,num_threads, &
              int_z,sk_1,sk_1_rdg,sk_z,ic,sc,bad_int_cv,bad_ext_cv, &
              update_p_comp, refine_dist
          ! u=horizontal(x in EASE grid) ice velocity, v=vertical(y in EASE grid) ice zelocity
          integer :: mi,mi1,sc1
          integer, dimension(8) :: mi_adv,i_adv,j_adv
          integer, dimension(0:8) :: p_order
          integer, dimension(3) :: out_case
          double precision :: u,v,flag,lost,sklost,adv_in_m3,adv_in_m2
          double precision :: z_th_min_multiplier,target_h,h_in,max_2d
          double precision :: adv_corr_static,smalg_lost,adv_sum,tmp1,tmp2,tmp3
          double precision :: ice_p, a_convg1, p_comp, mmntm, r_adv_lambda, hm
          double precision :: a_rdg_tot, f_ridge_sum,a_leftover,rdg_f,h_scale
          double precision :: a_scale,a_melt,af_total,id_insert,id_sw_gap
          double precision :: a_rdg_rdg,a_scale_sum,void_f
          double precision , dimension(0:8) :: p_arr
          double precision,  dimension(ida_n) :: apn, ff,id_ic,af_ic

          ! sia2_env_brine_calc
          double precision :: t_melt,t_mean,s_mean,d_mean,bs_mean,bd_mean,bv_mean,heat_mean
          double precision :: alpha0,alpha1,alpha2,alpha3

          ! sia2_env_temp_from_heat
          ! double precision :: d_mean,heat_mean,t_mean,s_mean
          double precision ::  dz_mean,aq,bq,cq

          ! sia2_env_layer_thickness
          logical :: keep_growing,maxed_depth,ignore_h_max
          integer :: int_z_new,z_snow_new
          double precision ::  r_depth,r_depth_last,t_depth,target_depth,layer_divisor,testvar
          double precision, dimension(z_max) :: th_new,id_new,debug_z,debug_start

          ! sia2_env_adj_ice_boundaries
          ! integer :: ii,jj
          logical :: alg_fixed,remapping
          double precision :: interim_top,interim_bot,old_layer_top,new_layer_bot
          double precision :: old_layer_bot,z1,z2,dz,heat_total,dz_total,new_layer_top
          double precision :: flooded,t_flood,s_flood,d_flood,bs_flood,bd_flood
          double precision :: bv_flood,heat_flood,c_gl,c_gl_sal,s_gl,f_depth,snowh1m
          double precision :: dz_sk,sk_th_new,percent_converge,dhdt_cm_s,sk_gl
		  double precision :: s_pond_new,h_pond_new,d_pond_new,smalg_pond_new
		  double precision :: poc_pond_new,no3_pond_new,nh4_pond_new,po4_pond_new
          double precision :: sioh4_pond_new,th_pond_new,snow_melt_now,melt_flood
          double precision :: mf_depth,melt_drain,c_gl_copy
          double precision :: sh_prev_sum,sh_offset_sum,age_new,ridged_new,snow_dist_new
          double precision :: snow_rd_new
          double precision, dimension(z_max) :: t_new,s_new,d_new,bs_new,bd_new
          double precision, dimension(z_max) :: bv_new,smalg_new,poc_new,smalg_pre_map
          double precision, dimension(z_max) :: melt_new,th_old,heat_new,d_new_snow
          double precision, dimension(z_max) :: heat_new_snow,th_new_snow,d_snow_new,th_snow_new,heat_snow_new
          double precision, dimension(z_max+1) :: NO3_new,NH4_new,PO4_new,SiOH4_new

          ! ridging temporary ice grid cells
          integer :: thread_n

          ! temp erase me below
          double precision, dimension(z_max) :: t_prev,t_last,t_next,ki


          ! sia2_env_adj_snow_boundaries
          integer :: z_snow,z_last
          double precision :: airtemp1c=0.d0,snowfall_d,c_snow

          double precision, dimension(z_max) :: th_debug



          ! var to see if ice gets too deep anywhere for current physical model
          max_d_exceeded = 0
          convg_maxed = 0
          maxed_depth = .false.
          ignore_h_max = .false.

          ! chunk size calc
          chunk_size = 16

          ! pre-calculate apn (participation calc) denominator
          apn_denom = 1./(1. - exp(-1.*r_a_star))

          ! allocate temporary structures
					thread_n = 1
#ifdef omp_on
          thread_n = omp_get_max_threads()
#endif
          allocate(ice_ridged(sda_n,ida_n,thread_n),ice_temp(thread_n),ice_temp2)


!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i,j,ii,jj,kk,jjj,iii,i_ease,j_ease,chunk,num_threads, &
!$OMP int_z,sk_1,sk_1_rdg, &
!$OMP sk_z,ic,sc,bad_int_cv,bad_ext_cv, update_p_comp, refine_dist, &
!$OMP mi,mi1,sc1, &
!$OMP mi_adv,i_adv,j_adv, &
!$OMP p_order, &
!$OMP out_case, &
!$OMP u,v,flag,lost,sklost,adv_in_m3,adv_in_m2, &
!$OMP z_th_min_multiplier,target_h,h_in,max_2d, &
!$OMP adv_corr_static,smalg_lost,adv_sum,tmp1,tmp2,tmp3, &
!$OMP ice_p, a_convg1, p_comp, mmntm, r_adv_lambda, hm, &
!$OMP a_rdg_tot, f_ridge_sum,a_leftover,rdg_f,h_scale, &
!$OMP a_scale,a_melt,af_total,id_insert,id_sw_gap, &
!$OMP a_rdg_rdg,a_scale_sum,void_f, &
!$OMP p_arr, &
!$OMP apn, ff,id_ic,af_ic, &
!$OMP t_melt,t_mean,s_mean,d_mean,bs_mean,bd_mean,bv_mean,heat_mean, &
!$OMP alpha0,alpha1,alpha2,alpha3, &
!$OMP dz_mean,aq,bq,cq, &
!$OMP keep_growing,maxed_depth,ignore_h_max, &
!$OMP int_z_new,z_snow_new, &
!$OMP r_depth,r_depth_last,t_depth,target_depth,layer_divisor,testvar, &
!$OMP th_new,id_new,debug_z,debug_start, &
!$OMP alg_fixed,remapping, &
!$OMP interim_top,interim_bot,old_layer_top,new_layer_bot, &
!$OMP old_layer_bot,z1,z2,dz,heat_total,dz_total,new_layer_top, &
!$OMP flooded,t_flood,s_flood,d_flood,bs_flood,bd_flood, &
!$OMP bv_flood,heat_flood,c_gl,c_gl_sal,s_gl,f_depth,snowh1m, &
!$OMP dz_sk,sk_th_new,percent_converge,dhdt_cm_s,sk_gl, &
!$OMP s_pond_new,h_pond_new,d_pond_new,smalg_pond_new, &
!$OMP poc_pond_new,no3_pond_new,nh4_pond_new,po4_pond_new, &
!$OMP sioh4_pond_new,th_pond_new,snow_melt_now,melt_flood, &
!$OMP mf_depth,melt_drain,c_gl_copy, &
!$OMP sh_prev_sum,sh_offset_sum,age_new,ridged_new,snow_dist_new, &
!$OMP snow_rd_new, &
!$OMP t_new,s_new,d_new,bs_new,bd_new, &
!$OMP bv_new,smalg_new,poc_new,smalg_pre_map, &
!$OMP melt_new,th_old,heat_new,d_new_snow, &
!$OMP heat_new_snow,th_new_snow,d_snow_new,th_snow_new,heat_snow_new, &
!$OMP NO3_new,NH4_new,PO4_new,SiOH4_new, &
!$OMP thread_n, &
!$OMP t_prev,t_last,t_next,ki, &
!$OMP z_snow,z_last, &
!$OMP airtemp1c,snowfall_d,c_snow, &
!$OMP th_debug)
!$OMP DO SCHEDULE(DYNAMIC,chunk_size)

		  do mi=1,tcells
		  if (m(mi)%status .eq. 1) then
!              do sc=1,sda_n
!              do ic=1,ida_n
!                  if (ice(sc,ic,mi)%af .gt. 0.) then
!                       tmp1=0.
!                       do ii=1,ice(sc,ic,mi)%z
!                       if (ice(sc,ic,mi)%smalg(ii) .gt. 1.e6) then
!                           tmp1 = 1.
!                       endif
!                       enddo
!                       if (tmp1 .eq. 1.) then
!                           print *,'Before Advection Gridding'
!                           print *, mi,'--',ic,' ------------------------------------'
!                           print *,ice(sc,ic,mi)%smalg
!                           testvar = -1
!                           testvar = sqrt(testvar)
!                       endif
!                   endif
!              enddo
!              enddo

					thread_n = 1
#ifdef omp_on
					thread_n = omp_get_thread_num() ! find thread index to temp ice structures
					thread_n = thread_n+1
#endif
          do ic=1,ida_n
          do sc=1,sda_n
          if (adv(mi)%ice(sc,ic)%af .ge. af_min) then


              ! find full EASE domain reference for current cell
              i_ease = m(mi)%grid_h
              j_ease = m(mi)%grid_v
              i = pj%mdh(i_ease,j_ease)
              j = pj%mdv(i_ease,j_ease)
              sk_z = adv(mi)%ice(sc,ic)%z
              int_z = sk_z - z_sk
              sk_1 = int_z + 1

              ! UPDATE SNOW VARS
              ! ----------------------------------------------------

              ! update total snow depth
              ice(sc,ic,mi)%sh_prev = adv(mi)%snow(sc,ic)%sh_prev

              ! update gridded snow
              if (adv(mi)%snow(sc,ic)%z .gt. 0) then

                  ! update snow grid
                  ice(sc,ic,mi)%snow%z = adv(mi)%snow(sc,ic)%z
                  ice(sc,ic,mi)%snow%depth = adv(mi)%snow(sc,ic)%depth
                  ice(sc,ic,mi)%sh_prev = adv(mi)%snow(sc,ic)%sh_prev
                  ice(sc,ic,mi)%snow%th = adv(mi)%snow(sc,ic)%th_new

                  ! update snow physics
                  do ii=1,ice(sc,ic,mi)%snow%z
                       heat_total = ice(sc,ic,mi)%snow%heat(ii)/ice(sc,ic,mi)%snow%d(ii)
                       t_mean = (0.2309-sqrt(0.05331481d0+0.0136d0* &
                           (heat_snow0-heat_total)))/(-0.0068d0)
                       ice(sc,ic,mi)%snow%t(ii) = min(t_mean - kelvin0,d0_)
                  enddo

              endif

              ! UPDATE ICE VARS
              ! ----------------------------------------------------

              ! debugging variables
              adv(mi)%ice(sc,ic)%th_old = ice(sc,ic,mi)%th
              adv(mi)%ice(sc,ic)%th_debug = adv(mi)%ice(sc,ic)%th_new

              ! update ice grid - skeletal layer shold never be above
              ! sk_th_max, since only area is advecting ...
              ice(sc,ic,mi)%z = adv(mi)%ice(sc,ic)%z
              ice(sc,ic,mi)%th = adv(mi)%ice(sc,ic)%th_new
              ice(sc,ic,mi)%id = adv(mi)%ice(sc,ic)%id_new
              int_z = ice(sc,ic,mi)%z - z_sk

             ! update ice physics
              do ii=1,ice(sc,ic,mi)%z

                  if(ii .le. int_z) then
                      dz_mean = 1. ! division by dz_mean is included in calc below
                      d_mean = ice(sc,ic,mi)%d(ii)
                      s_mean = ice(sc,ic,mi)%s(ii)
                      heat_mean = ice(sc,ic,mi)%heat(ii)

#include "sia2_env_temp_from_heat.inc.f90"

                  else
                      ! skeletal paramaters are constant
                      s_mean = f(mi)%s*c_5
                      t_mean = min(f(mi)%t,s_mean*mu)
                  endif

#include "sia2_env_brine_calc.inc.f90"

                  ice(sc,ic,mi)%t(ii) = t_mean
                  ice(sc,ic,mi)%s(ii) = s_mean
                  ice(sc,ic,mi)%d(ii) = d_mean
                  ice(sc,ic,mi)%bs(ii) = bs_mean
                  ice(sc,ic,mi)%bd(ii) = bd_mean
                  ice(sc,ic,mi)%bv(ii) = bv_mean
                  ice(sc,ic,mi)%heat(ii) = heat_mean

              enddo

              ! Update areal fraction before convergence/divergence calcs
              ! ---------------------------------------------------------
              ice(sc,ic,mi)%af = adv(mi)%ice(sc,ic)%af


          elseif (adv(mi)%ice(sc,ic)%af .gt. 0.d0 .or. ice(sc,ic,mi)%z .gt. 0) then
              ! category has too little ice and must be elimintated.

              adv(mi)%ice(sc,ic)%af = 0.d0

              ! zero layer
#include "sia2_env_null_params.f90"

          endif ! end if does category have (enough) ice? test
          enddo ! end of snow category loop
          enddo ! end of ice category loop


!              do sc=1,sda_n
!              do ic=1,ida_n
!                  if (ice(sc,ic,mi)%af .gt. 0.) then
!                       tmp1=0.
!                       do ii=1,ice(sc,ic,mi)%z
!                       if (ice(sc,ic,mi)%smalg(ii) .gt. 1.e6) then
!                           tmp1 = 1.
!                       endif
!                       enddo
!                       if (tmp1 .eq. 1.) then
!                           print *,'After Advection Gridding'
!                           print *, mi,'--',ic,' ------------------------------------'
!                           print *,ice(sc,ic,mi)%smalg
!                           testvar = -1
!                           testvar = sqrt(testvar)
!                       endif
!                   endif
!              enddo
!              enddo


          ! Convergence calculations (divergence is already accounted for above)
          ! --------------------------------------------------------------------
          adv(mi)%af_new = 0.


! 6. Melt or drop ice, rather than converge if convergence was not allowed
! =================================================================================

          if (adv(mi)%a_drop .gt. 0. .or. adv(mi)%a_melt .gt. 0.) then

              ! first melt/drop the smallest category, since these are
              ! most often created erroneously, and are also the most
              ! likely to actually melt
              ! ---------------------------------------------------------
              ic = 1
              do iii=1,2

                  if (iii .eq. 1 .and. adv(mi)%a_melt .gt. 0.d0) then
                      af_total =  adv(mi)%a_melt
                  elseif (iii .eq. 2 .and. adv(mi)%a_drop .gt. 0.d0) then
                      af_total = adv(mi)%a_drop
                  else
                      af_total = 0.d0
                  endif

                  ! find category 1 area
                  tmp1 = 0.d0
                  do sc=1,sda_n
                      tmp1 = tmp1 + ice(sc,ic,mi)%af
                  enddo

                  if (af_total .gt. 0.d0 .and. tmp1 .gt. 0.d0) then

                      tmp1 = tmp1*cell_area ! convert to m^2
                      tmp2 = tmp1 - af_total
                      if (tmp2 .lt. af_min) then
                          tmp2 = 0.d0
                          af_total = af_total -  tmp1
                      else
                          af_total = 0.d0
                      endif

                      tmp1 = tmp2/tmp1 ! ratio of new categoy area to old category area

                      do sc=1,sda_n
                        if (ice(sc,ic,mi)%z > 0) then
                          ! record lost ice and algae
                          tmp3 = ice(sc,ic,mi)%af*(1.d0 - tmp1)

                          ! record ice volume as advective loss domain loss
                          int_z = ice(sc,ic,mi)%z - z_sk
                          m(mi)%adv_loss = m(mi)%adv_loss + &
                            ice(sc,ic,mi)%id(ice(sc,ic,mi)%z)*tmp3*cell_area*1.e6
                          ! record biomass lost from model domain if melting
                          if (iii .eq. 1) then
                              do ii=1,ice(sc,ic,mi)%z
                                 if (ice(sc,ic,mi)%smalg(ii) .gt. 0.) then
                                    bv_mean = ice(sc,ic,mi)%bv(ii)*c_001
                                    smalg_lost = (ice(sc,ic,mi)%smalg(ii)*bv_mean - min_alg)
                                    smalg_lost = max(0.,smalg_lost)
                                    ! converting to g/pixel - (cell_area*1.e6)*1g/1000mg = 1000. scale factor
                                    m(mi)%bm_lost = m(mi)%bm_lost + smalg_lost* &
                                      ice(sc,ic,mi)%th(ii)*cell_area*1.e3*tmp3 ! g/pixel
                                  endif
                              enddo
                           endif

                          ice(sc,ic,mi)%af = max(0.d0,ice(sc,ic,mi)%af*tmp1)
                        endif
                      enddo
                  endif
              enddo

              ! drop equal percentage from each remaining category if
              ! there is anything left to drop
              ! ---------------------------------------------------------

              if (af_total .gt. 0.d0) then

                  ! find new area
                  tmp1 = 0.d0
                  do ic=2,ida_n
                  do sc=1,sda_n
                      tmp1 = tmp1 + ice(sc,ic,mi)%af
                  enddo
                  enddo

                  tmp2 = af_total/(tmp1*cell_area)  ! new area multiplier (af_total is a_drop)

                  do ic=2,ida_n
                  do sc=1,sda_n
                      if (ice(sc,ic,mi)%z .gt. 0) then
                          tmp3 = tmp2*ice(sc,ic,mi)%af

                          ! record lost ice and algae
                          ! record ice volume as advective loss domain loss
                          int_z = ice(sc,ic,mi)%z - z_sk
                          m(mi)%adv_loss = m(mi)%adv_loss + &
                            ice(sc,ic,mi)%id(ice(sc,ic,mi)%z)*tmp3*cell_area*1.e6

                          ice(sc,ic,mi)%af = max(0.d0,ice(sc,ic,mi)%af - tmp3)

                      endif

                  enddo
                  enddo

              endif

              ! reduce total convergence area by area melted
              adv(mi)%a_convg = adv(mi)%a_convg - adv(mi)%a_drop - adv(mi)%a_melt

              do ic=1,ida_n
              do sc=1,sda_n
                  if (ice(sc,ic,mi)%af .lt. af_min) then
                       adv(mi)%ice(sc,ic)%af = 0.d0
                       ! zero layer
#include "sia2_env_null_params.f90"

                  endif
              enddo
              enddo


          endif ! end if a_drop



          if (adv(mi)%a_convg .gt. 0.) then

!          if (mi .eq. 22505) then
!              print *,'af: ',ice(:,:,mi)%af
!              print *,'aconvg: ',adv(mi)%a_convg
!              print *,'id: ',ice(2,4,mi)%id(40)
!              print *,'s: ',ice(2,4,mi)%s
!              print *,'t: ',ice(2,4,mi)%t
!              print *,'snow: ',ice(2,4,mi)%snow%depth
!              print *,'snow: ',ice(2,4,mi)%snow%th
!          endif

! 7. Calculate ridging participation and extent
! =================================================================================
#include "sia2_adv_3_allocation.inc.f90"


! 8. Perform ridging redistribution
! =================================================================================
#include "sia2_adv_3_ridge.inc.f90"

              !print *, ice(:,mi)%snow%depth

!           if (mi .eq. 22505) then
!              print *,'--------------------------------------------------'
!              print *,'af after: ',ice(:,:,mi)%af
!              print *,'id: ',ice(2,4,mi)%id(40)
!              print *,'s: ',ice(2,4,mi)%s
!              print *,'t: ',ice(2,4,mi)%t
!              print *,'snow: ',ice(2,4,mi)%snow%depth
!              print *,'snow: ',ice(2,4,mi)%snow%th
!         endif

          else

! 9. Setup for new ice growth (happens later in sia2_env_redist.f)
! =================================================================================

				  if (f(mi)%ws .eq. c9999) then
				      tmp1 = sqrt(f(mi)%u10**2 + f(mi)%v10**2)
				  else
				      tmp1 = f(mi)%ws
				  endif

                  if(((f(mi)%at-kelvin0)-f(mi)%t)*tmp1 .le. f_limit) then
                      adv(mi)%af_new = abs(adv(mi)%a_convg/cell_area)
                  endif

          endif
          ! end of convergence test

          if (mi .eq. 20512) then
              tmp1 = 0.
          endif

!              do sc=1,sda_n
!              do ic=1,ida_n
!                  if (ice(sc,ic,mi)%af .gt. 0.) then
!                       tmp1=0.
!                       do ii=1,ice(sc,ic,mi)%z
!                       if (ice(sc,ic,mi)%smalg(ii) .gt. 1.e6) then
!                           tmp1 = 1.
!                       endif
!                       enddo
!                       if (tmp1 .eq. 1.) then
!                           print *,'After Ridging Gridding'
!                           print *, mi,'--',ic,' ------------------------------------'
!                           print *,ice(sc,ic,mi)%smalg
!                           testvar = -1
!                           testvar = sqrt(testvar)
!                       endif
!                   endif
!              enddo
!              enddo



          ! return tracers to brine concentrations
          do ic=1,ida_n
          do sc=1,sda_n
              if (ice(sc,ic,mi)%af .gt. 0.d0) then
                  do ii=1,ice(sc,ic,mi)%z
                      if (ice(sc,ic,mi)%bv(ii) .eq. 0.d0) then
                          testvar = -1
                      endif
                      bv_mean = 1./(ice(sc,ic,mi)%bv(ii)*c_001)
                      ice(sc,ic,mi)%NO3(ii) = ice(sc,ic,mi)%NO3(ii)*bv_mean
                      ice(sc,ic,mi)%NH4(ii) = ice(sc,ic,mi)%NH4(ii)*bv_mean
                      ice(sc,ic,mi)%PO4(ii) = ice(sc,ic,mi)%PO4(ii)*bv_mean
                      ice(sc,ic,mi)%SiOH4(ii) = ice(sc,ic,mi)%SiOH4(ii)*bv_mean
                      ice(sc,ic,mi)%smalg(ii) = ice(sc,ic,mi)%smalg(ii)*bv_mean
                      ice(sc,ic,mi)%poc(ii) = ice(sc,ic,mi)%poc(ii)*bv_mean
                  enddo
              endif
          enddo
          enddo


          endif ! end status test
          enddo ! end tcell loop


!$OMP END DO
!$OMP END PARALLEL

          if (max_d_exceeded .gt. 0) then
              print *, 'Maximum depth exceeded in ',max_d_exceeded,' cells.'
          endif
          if (convg_maxed .gt. 0) then
              print *, 'Convergence routine maxed out in ',convg_maxed,' cells.'
          endif

          deallocate(ice_ridged,ice_temp,ice_temp2)

      end SUBROUTINE sia2_adv_grid_ridge
