
      subroutine sia2_env_create_ice(m,ice,f,v_f,pl,mi)

          use sia2_globals
          implicit none

          ! function arguments
          !----------------------------------------------------
          type (meta_type) :: m(tcells)
          type (ice_type) :: ice(sda_n,ida_n,tcells)
          type (forcing_type) :: f(tcells)
          type (boundary_f_type) :: v_f
          type (platelet_type) :: pl
          integer :: ic,sc,mi

          ! private variables
          !----------------------------------------------------
          logical :: &
          keep_growing,maxed_depth,ignore_h_max


          integer :: &
          z_last,ii,kk,jj,z_snow,int_z_new, sk_1,sk_z,int_z,max_d_exceeded, &
          z_snow_new

          double precision :: &
          t_mean,s_mean,d_mean,bs_mean,bd_mean,bv_mean,heat_mean, &
          alpha0,alpha1,alpha2,alpha3,t_melt,layer_divisor, &
          t_depth,snowh1m,f_snow,tmp1,airtemp1c,t_interface,r_depth, &
          flooded,old_layer_top,new_layer_bot,old_layer_bot,new_layer_top, &
          dz_total,heat_total,interim_top,interim_bot,z1,z2,dz,snowfall_d, &
          testvar,top_f,bot_f,T_flood,S_flood,fresh_tmp,dsdt_tmp, &
          d_flood,bs_flood,bd_flood,bv_flood,heat_flood,melt_drain

          double precision, dimension(z_max) :: &
          t_new,th_new,d_new,id_new,melt_new,th_snow_new,d_snow_new,heat_snow_new

          ! start create ice
          !----------------------------------------------------
          ignore_h_max = .false.

          airtemp1c = f(mi)%at - kelvin0

          ! update status to that of an active cell
          m(mi)%status = 1
          m(mi)%pxl_h_offset = cell_side/2.
          m(mi)%pxl_v_offset = cell_side/2.


          if (start .and. (override_ic .ne. 1)) then
              if (boundary_file .eq. 1) then
                  if (init_ice_th > 0.) then
                      t_depth = init_ice_th
                  else
                      if(f(mi)%ic_interp .gt. 0.9) then
                        t_depth = 1.28           ! max ice height, for icecon >= 80%

                      else
                        t_depth = h_min + (1.28-h_min)*(f(mi)%ic_interp-0.1)/0.8 ! varied height between min mand max height for 0.1 < icecon < 0.8
                      endif
                  endif
! ----------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! init lines for Grossi/McMurdo 1982 validation only
! ----------------------------------------------------------------------------------
!                      t_depth = 1.95
! ----------------------------------------------------------------------------------
! init lines for Grossi/McMurdo 1982 validation only
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! ----------------------------------------------------------------------------------
              else
                  if (init_ice_th > 0.) then
                      t_depth = init_ice_th
                  else
                      if(f(mi)%ic_interp .gt. 0.9) then
                        t_depth = 1.28           ! max ice height, for icecon >= 80%

                      else
                        t_depth = h_min + (1.28-h_min)*(f(mi)%ic_interp-0.1)/0.8 ! varied height between min mand max height for 0.1 < icecon < 0.8
                      endif
                  endif
              endif
          else
              if ((override_ic .eq. 1) .and. (init_ice_th > 0.)) then
                  t_depth = init_ice_th
              else
                  t_depth = h_min !0.21d0 <- why did we need to set this to 21? not sure
              endif
          endif

          ! find new ice category
          ic = 1
          sc = 1
          if (ida_n .gt. 1) then
              do while (ic_h_max(ic) .lt. t_depth)
                  ic = ic + 1
              enddo
          endif

          ice(sc,ic,mi)%af = f(mi)%ic_interp ! set default areal fraction
          m(mi)%a_new = m(mi)%a_new + ice(sc,ic,mi)%af*cell_area  ! record new area

          ice(sc,ic,mi)%pur = 0.d0
          ice(sc,ic,mi)%ed_w = 0.d0
          ice(sc,ic,mi)%snow%ed_w = 0.d0

          ! heat test -----------------------------------------------
          ! t_depth = 0.32
          ! heat test -----------------------------------------------

          ! initially, assume no infiltration, but zero freeboard
          ice(sc,ic,mi)%fbh = 0.

          if ((snow_model .eq. 1) .and. &
          ((start .eqv. .FALSE.) .or. (boundary_file .eq. 1) .or. (override_ic .eq. 1))) then
              if ((override_ic .eq. 1) .or. (init_snow_th .gt. 0.)) then
                  ice(sc,ic,mi)%sh_prev = init_snow_th
              else
                  ice(sc,ic,mi)%sh_prev = 0.d0
              endif

! ----------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! testing for SHEBA snow
! ----------------------------------------------------------------------------------
!              ice(sc,ic,mi)%sh_prev = f(mi)%sh_interp*0.01
! ----------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! testing for SHEBA snow
! ----------------------------------------------------------------------------------

          else
            if (init_snow_th .gt. 0.) then
              ice(sc,ic,mi)%sh_prev = init_snow_th
            else
              if (f(mi)%sh_interp .eq. c9999) then
                  ice(sc,ic,mi)%sh_prev = 0.d0
              else
                  ice(sc,ic,mi)%sh_prev = f(mi)%sh_interp*0.01
              endif
            endif
          endif

          ! test
          !ice(sc,ic,mi)%sh_prev = 0.4

          snowh1m = ice(sc,ic,mi)%sh_prev

              ! test
              !snowh1m = 0.

          ! fractional snow cover
          if (h_snowpatch .eq. 0.) then
              f_snow = 1.
          else
              f_snow = snowh1m/(snowh1m + h_snowpatch)
          endif
          ! find mean ice surface temp from mean snowdepth
          !t_interface = airtemp1c + snowh1m*(f(mi)%t - airtemp1c) &
          !          /(snowh1m + t_depth)
          tmp1 = ksnow*2.2/(2.2*snowh1m + ksnow*t_depth) ! gamma in Arrigo et al. 1993 appendix
          t_interface = (0.2456008*airtemp1c + tmp1*f(mi)%t)/(0.2456008+tmp1)
          t_interface = min(f(mi)%t,t_interface)


! ---------------------------------------------------------------------------------
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! init lines for McMurdo 1982/Grossi forcing validation only !!!!!!!!!
! comments these out for all other runs
! ----------------------------------------------------------------------------------


          !t_interface = -8.0

! ---------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! ----------------------------------------------------------------------------------



!          if (airtemp1c .gt. 0.) then
              snowfall_d = den_s_wet*1.e6
!          else
!              snowfall_d = den_s_dry*1.e6
!          endif

          r_depth = snowh1m
          ice(sc,ic,mi)%snow%z = 0 ! tells grid/remap routines to create new snow
          ! tells remap to use linear interpolation of temperature in ice/snow pack
          if (iif .gt. 0) then
              ice(sc,ic,mi)%snow%ts = f(mi)%t
          else
              ice(sc,ic,mi)%snow%ts = t_interface
          endif

#include "sia2_env_snow_grid.inc.f90"

          flooded = 0.d0
          melt_drain = 0.d0

#include "sia2_env_snow_remap.inc.f90"

             ! set melt status for snow cover to 1 - this is old wet snow
             do ii=1,ice(sc,ic,mi)%snow%z
                 ice(sc,ic,mi)%snow%d(ii) = den_s_wet*1.e6
                 ice(sc,ic,mi)%snow%melt(ii) = 1.
             enddo

          r_depth = t_depth

#include "sia2_env_ice_grid.inc.f90"

          sk_1 = int_z_new + 1
          sk_z = int_z_new + z_sk
          int_z = int_z_new
          ice(sc,ic,mi)%z = sk_z

          ! setup flooding vars in case there is initial flooding specified
          top_f = fl_ratio
          bot_f = 1.-top_f ! remaining vol fraction taken up by seawater

          T_flood = f(mi)%t
          S_flood = f(mi)%s*bot_f
          S_mean = S_flood
          T_mean = f(mi)%t

#include "sia2_env_brine_calc.inc.f90"

          d_flood = d_mean
          bs_flood = bs_mean
          bd_flood = bd_mean
          bv_flood = bv_mean
          heat_flood = heat_mean

          ! setup individual internal ice layers
          ! -----------------------------------------------------------
          do ii=1,sk_z

              ! initial ice does not have fatty brine channels
              ice(sc,ic,mi)%bced(ii) = 0

              ! update thickness & ice depth
              ice(sc,ic,mi)%th(ii) = th_new(ii)
              ice(sc,ic,mi)%id(ii) = id_new(ii)

               if (ii .lt. sk_1) then

                  ! new ice layer salinity
               if (init_ice_type .gt. 0) then
                      if (init_ice_type .eq. 2) then
                          ! interpolate default new-ice salinity curve to get initial salinity
                          call sia2_interpolate_salinity2( &
                          (ice(sc,ic,mi)%id(ii)-ice(sc,ic,mi)%th(ii)/2)/t_depth,ice(sc,ic,mi)%s(ii))
                      elseif (init_ice_type .eq. 1) then
                          ! interpolate default new-ice salinity curve to get initial salinity
                          call sia2_interpolate_salinity1( &
                          (ice(sc,ic,mi)%id(ii)-ice(sc,ic,mi)%th(ii)/2)/t_depth,ice(sc,ic,mi)%s(ii))
                      else
                          ! constant ice salinity
                          ice(sc,ic,mi)%s(ii) = s_const
                      endif
               else
                  if (boundary_file .eq. 1) then
                      if (iis .eq. 3) then
                          ! multi-year Norwegian 2014 intercompare option
                          call sia2_interp_value_from_data(norwegian_sal, &
                                                      norwegian_dth_s, &
                                                      15, &
                                                      (ice(sc,ic,mi)%id(ii)-ice(sc,ic,mi)%th(ii)/2), &
                                                      ice(sc,ic,mi)%s(ii))

                      elseif (iis .eq. 2) then
                          ! interpolate default new-ice salinity curve to get initial salinity
                          call sia2_interpolate_salinity2( &
                          (ice(sc,ic,mi)%id(ii)-ice(sc,ic,mi)%th(ii)/2)/t_depth,ice(sc,ic,mi)%s(ii))
                      elseif (iis .eq. 1) then
                          ! interpolate default new-ice salinity curve to get initial salinity
                          call sia2_interpolate_salinity1( &
                          (ice(sc,ic,mi)%id(ii)-ice(sc,ic,mi)%th(ii)/2)/t_depth,ice(sc,ic,mi)%s(ii))
                      else
                          ! constant ice salinity
                          ice(sc,ic,mi)%s(ii) = s_const
                      endif
                  elseif (start .or. iis .eq. 2) then
                      ! interpolate default new-ice salinity curve to get initial salinity
                      call sia2_interpolate_salinity2( &
                      (ice(sc,ic,mi)%id(ii)-ice(sc,ic,mi)%th(ii)/2)/t_depth,ice(sc,ic,mi)%s(ii))
                  elseif (iis .eq. 1) then
                      ! interpolate default new-ice salinity curve to get initial salinity
                      call sia2_interpolate_salinity1( &
                      (ice(sc,ic,mi)%id(ii)-ice(sc,ic,mi)%th(ii)/2)/t_depth,ice(sc,ic,mi)%s(ii))
                  else
                      ! constant ice salinity
                      ice(sc,ic,mi)%s(ii) = s_const
                  endif
               endif

                  ! new ice layer temperature
                  if (iif .gt. 0.) then
                       ice(sc,ic,mi)%t(ii) = min(f(mi)%t,mu*ice(sc,ic,mi)%s(ii))
                  else
                      if (iit .eq. 1 .and. .not. start) then
!                      if (iit .eq. 1) then
                          ! isothermal at Tw
                         ice(sc,ic,mi)%t(ii) = min(f(mi)%t,mu*ice(sc,ic,mi)%s(ii))
                      else
                          ! interpolate linearly between t_interface and water temp
                          ice(sc,ic,mi)%t(ii) = t_interface - &
                              (t_interface-min(f(mi)%t,mu*f(mi)%s)) * &
                              (ice(sc,ic,mi)%id(ii) - ice(sc,ic,mi)%th(ii)/2)/t_depth
                      endif
                  endif


              else

                  ice(sc,ic,mi)%s(ii) = f(mi)%s*c_5
                  ice(sc,ic,mi)%t(ii) = min(f(mi)%t,mu*f(mi)%s)

              endif

              ! set default brine paramaters
              T_mean = ice(sc,ic,mi)%t(ii)
              S_mean = ice(sc,ic,mi)%s(ii)

#include "sia2_env_brine_calc.inc.f90"

              ! update brine-based physics
              ice(sc,ic,mi)%d(ii) = d_mean
              ice(sc,ic,mi)%bs(ii) = bs_mean
              ice(sc,ic,mi)%bd(ii) = bd_mean
              ice(sc,ic,mi)%bv(ii) = bv_mean
              ice(sc,ic,mi)%heat(ii) = heat_mean

              ! give minimum algae concentration
              ice(sc,ic,mi)%smalg(ii) = alg_wc/(ice(sc,ic,mi)%bv(ii)*c_001)  !concentrate algae into brine volume


              ! start with zero POC
              ice(sc,ic,mi)%poc(ii) = 0.

              if (start .and. boundary_file .eq. 0) then
                  tmp1 = f(mi)%s/ice(sc,ic,mi)%bs(ii)
              elseif (iin .eq. 1) then
                  ! scale nutrients linearly from 0->full from top->bottom, mimicking fall surface depletetion
                  tmp1 = ((ice(sc,ic,mi)%id(ii)-ice(sc,ic,mi)%th(ii)/2)/(t_depth+sk_th_min*z_sk))
              elseif (iin .eq. 2) then
                  tmp1 = ((ice(sc,ic,mi)%id(ii)-ice(sc,ic,mi)%th(ii)/2)/(t_depth+sk_th_min*z_sk))**2
              elseif (iin .eq. 3) then
                  tmp1 = n_f
              elseif (iin .eq. 4) then
                  ! same as seawater
                  tmp1 = f(mi)%s/ice(sc,ic,mi)%bs(ii)
              else
                  tmp1 = 1.
              endif


              ! intial nutrients
              ice(sc,ic,mi)%no3(ii)= f(mi)%no3*ice(sc,ic,mi)%bs(ii)/f(mi)%s*tmp1  ! scale conservatively with salinity
              ice(sc,ic,mi)%nh4(ii)= f(mi)%nh4*ice(sc,ic,mi)%bs(ii)/f(mi)%s*tmp1  ! scale conservatively with salinity
              ice(sc,ic,mi)%po4(ii)= f(mi)%po4*ice(sc,ic,mi)%bs(ii)/f(mi)%s*tmp1  ! scale conservatively with salinity
              ice(sc,ic,mi)%sioh4(ii)= f(mi)%sioh4*ice(sc,ic,mi)%bs(ii)/f(mi)%s*tmp1  ! scale conservatively with salinity

              if (iis .eq. 3) then
                  ! multi-year Norwegian 2014 intercompare option
                  call sia2_interp_value_from_data(norwegian_N, &
                                              norwegian_dth, &
                                              16, &
                                              (ice(sc,ic,mi)%id(ii)-ice(sc,ic,mi)%th(ii)/2), &
                                              ice(sc,ic,mi)%no3(ii))
                  ice(sc,ic,mi)%no3(ii)= ice(sc,ic,mi)%no3(ii)*ice(sc,ic,mi)%bs(ii)/f(mi)%s
                  call sia2_interp_value_from_data(norwegian_Si, &
                                              norwegian_dth, &
                                              16, &
                                              (ice(sc,ic,mi)%id(ii)-ice(sc,ic,mi)%th(ii)/2), &
                                              ice(sc,ic,mi)%sioh4(ii))
                  ice(sc,ic,mi)%sioh4(ii)= ice(sc,ic,mi)%sioh4(ii)*ice(sc,ic,mi)%bs(ii)/f(mi)%s
                  call sia2_interp_value_from_data(norwegian_smalg, &
                                              norwegian_dth, &
                                              16, &
                                              (ice(sc,ic,mi)%id(ii)-ice(sc,ic,mi)%th(ii)/2), &
                                              ice(sc,ic,mi)%smalg(ii))
                  ice(sc,ic,mi)%smalg(ii)= ice(sc,ic,mi)%smalg(ii)*ice(sc,ic,mi)%bs(ii)/f(mi)%s
                  ice(sc,ic,mi)%smalg(ii)=ice(sc,ic,mi)%smalg(ii)*c_n*gC_mC ! convert microM N to gC
                endif


           enddo


           if (use_pl .eq. 1) then
               pl%t = f(mi)%t
               pl%s = 0.8*f(mi)%s
               pl%bs = f(mi)%s
               pl%bd = 1000000. + pl%bs*800.
               pl%bv = 800.
               do ii=1,pl_max
                   pl%no3(ii) = f(mi)%no3
                   pl%nh4(ii) = f(mi)%nh4
                   pl%po4(ii) = f(mi)%po4
                   pl%sioh4(ii) = f(mi)%sioh4
                   pl%smalg(ii) = min_alg/(pl%bv*c_001)
                   pl%poc(ii) = 0.
               enddo
           endif

           ! set molecular sub-layer nutrient concentrations to seawater values
           ! -----------------------------------------------------------
           ice(sc,ic,mi)%poc(ml_z) = f(mi)%poc
           ice(sc,ic,mi)%no3(ml_z) = f(mi)%no3
           ice(sc,ic,mi)%nh4(ml_z) = f(mi)%nh4
           ice(sc,ic,mi)%po4(ml_z) = f(mi)%po4
           ice(sc,ic,mi)%sioh4(ml_z) = f(mi)%sioh4

           !print *,'Forcing nutrients: ',f(mi)%no3,f(mi)%po4,f(mi)%sioh4

           ! record new ice volume as new to model domain
           m(mi)%a_new = m(mi)%a_new + ice(sc,ic,mi)%id(int_z)* &
             ice(sc,ic,mi)%af*cell_area*1.e6

                        ! total agua and salt accounting
                        if (.not. start) then
                                call sia2_ice_mass(ice(sc,ic,mi),fresh_tmp,dsdt_tmp)
                                ! record fresh water flux (kg)
                                m(mi)%h2o_flux = m(mi)%h2o_flux - &
                                    fresh_tmp*ice(sc,ic,mi)%af*cell_area*1.e6
                                ! record salt flux (kg)
                                m(mi)%salt_flux = m(mi)%salt_flux - &
                                    dsdt_tmp*ice(sc,ic,mi)%af*cell_area*1.e6
                        endif


           ! set age and ridged vars
           ! -----------------------------------------------------------
           if (boundary_file .eq. 1) then
               ice(sc,ic,mi)%age = 0.d0
               ice(sc,ic,mi)%ridged = 0.d0
           else
               if (start) then
                   ice(sc,ic,mi)%age = 0.9d0
                   ice(sc,ic,mi)%ridged = 0.2d0
               else
                   ice(sc,ic,mi)%age = 0.0d0
                   ice(sc,ic,mi)%ridged = 0.0d0
               endif
           endif
           ice(sc,ic,mi)%snow_dist = dble(sd_num)*0.5
           ice(sc,ic,mi)%snow_rd = 0.d0





      end subroutine sia2_env_create_ice