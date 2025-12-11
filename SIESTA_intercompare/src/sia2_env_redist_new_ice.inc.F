! ======================================================================
! START: Add new ice
! ======================================================================


      if (adv(mi)%af_new .gt. af_min) then

          ! record new area
          m(mi)%a_new = m(mi)%a_new + adv(mi)%af_new*cell_area

          ic = 1  ! we are dealing with the 1st ice category
          sc = 1  ! putting new ice in non-flooded category


 					! zero salt/h2o accounting vars
 					dsdt_net = 0.d0
 					fresh_net = 0.d0

          if (ice(sc,ic,mi)%z .gt. 0) then
              sk_z = ice(sc,ic,mi)%z
              int_z = sk_z - z_sk
              sk_1 = int_z + 1

							! record dsdt_net, fresh_net
							call sia2_ice_mass(ice(sc,ic,mi),fresh_net,dsdt_net)
							fresh_net = fresh_net*ice(sc,ic,mi)%af ! scale by area here since it is increasing
							dsdt_net = dsdt_net*ice(sc,ic,mi)%af ! scale by area here since it is increasing

          else
              int_z = z_int_min
              sk_1 = int_z + 1
              sk_z = int_z + z_sk
              ice(sc,ic,mi)%sh_prev = d0_
          endif
          af_total = adv(mi)%af_new + ice(sc,ic,mi)%af ! total area fraction we are combining
          vf_total = adv(mi)%af_new*dble(z_int_min)*z_th_min + ice(sc,ic,mi)%af*ice(sc,ic,mi)%id(int_z)
          r_depth = vf_total/af_total

          if (i_ease .eq. 279 .and. j_ease .eq. 290) then
              testvar = 1
          endif

          do ii=1,sk_z

              if (ii .le. int_z) then

                  ! calculate fractional ice volumes
                  z1 = adv(mi)%af_new*dble(z_int_min)*z_th_min/vf_total
                  z2 = ice(sc,ic,mi)%af*ice(sc,ic,mi)%id(int_z)/vf_total

                  ! update th,id for combined internal ice
                  if (ice(sc,ic,mi)%z .gt. 0) then
                      ice(sc,ic,mi)%th(ii) = ice(sc,ic,mi)%th(ii)*r_depth/ice(sc,ic,mi)%id(int_z)
                  else
                      ice(sc,ic,mi)%th(ii) = z_th_min
                  endif
                  if (ii .eq. 1) then
                      ice(sc,ic,mi)%id(ii) = ice(sc,ic,mi)%th(ii)
                  else
                      ice(sc,ic,mi)%id(ii) = ice(sc,ic,mi)%id(ii-1) + ice(sc,ic,mi)%th(ii)
                  endif

                  ! find physics of new ice
                  dz = ice(sc,ic,mi)%id(ii) - ice(sc,ic,mi)%th(ii)/2.  ! depth

                  if (iis .eq. 2) then
                      ! interpolate default new-ice salinity curve to get initial salinity
                      call sia2_interpolate_salinity2(dz/r_depth,s_mean)
                  elseif (iis .eq. 1) then
                      ! interpolate default new-ice salinity curve to get initial salinity
                      call sia2_interpolate_salinity1(dz/r_depth,s_mean)
                  else
                      ! constant ice salinity
                      s_mean = s_const
                  endif

                  if (iit .eq. 1) then
                      ! isothermal at Tw
                      t_mean = min(f(mi)%t,s_mean*mu)
                  else

                      ! interpolate linearly between t_interface and water temp
                      t_mean = min(airtemp1c,f(mi)%t,f(mi)%s*mu)
                      t_mean = t_mean - (t_mean-min(f(mi)%t,f(mi)%s*mu)) * dz/r_depth
                  endif


#include "sia2_env_brine_calc.inc.f90"

									! record new ice salt and h2o flux



                  ! find new temp, salinity
                  heat_mean = heat_mean*z1 + ice(sc,ic,mi)%heat(ii)*z2
                  d_mean = d_mean*z1 + ice(sc,ic,mi)%d(ii)*z2
                  s_mean = s_mean*z1 + ice(sc,ic,mi)%s(ii)*z2
                  dz_mean = 1.

#include "sia2_env_temp_from_heat.inc.f90"

              else

                  ! calculate fractional ice volumes - different from internal ice, since
                  ! this is the skeletal ice volume is considred differently
                  tmp1 = adv(mi)%af_new*z_th_min + ice(sc,ic,mi)%af*ice(sc,ic,mi)%th(ii)
                  z1 = adv(mi)%af_new*z_th_min/tmp1
                  z2 = ice(sc,ic,mi)%af*ice(sc,ic,mi)%th(ii)/tmp1

                  ! update th,id for combined skeletal ice - try to conserve mass
                  ice(sc,ic,mi)%th(ii) = tmp1/af_total
                  ice(sc,ic,mi)%id(ii) = ice(sc,ic,mi)%id(ii-1) + ice(sc,ic,mi)%th(ii)

                  ! find physics of new ice skeletal layer
                  t_mean = f(mi)%s*mu
                  s_mean = f(mi)%s*c_5

              endif

#include "sia2_env_brine_calc.inc.f90"

              ! record ice physics
              ice(sc,ic,mi)%t(ii) = min(t_mean,min_t) ! ?? why min
              ice(sc,ic,mi)%s(ii) = s_mean
              ice(sc,ic,mi)%d(ii) = d_mean
              ice(sc,ic,mi)%bs(ii) = bs_mean
              ice(sc,ic,mi)%bd(ii) = bd_mean
              ice(sc,ic,mi)%bv(ii) = bv_mean
              ice(sc,ic,mi)%heat(ii) = heat_mean


              ! record tracers
              bv_mean = bv_mean*c_001
              ice(sc,ic,mi)%smalg(ii) = alg_wc*bv_mean*z1 + ice(sc,ic,mi)%smalg(ii)*z2
              ice(sc,ic,mi)%poc(ii) = ice(sc,ic,mi)%poc(ii)*z2  ! assume no POC in seawater
              ice(sc,ic,mi)%no3(ii) = f(mi)%no3*bv_mean*z1 + ice(sc,ic,mi)%no3(ii)*z2
              ice(sc,ic,mi)%nh4(ii) = f(mi)%nh4*bv_mean*z1 + ice(sc,ic,mi)%nh4(ii)*z2
              ice(sc,ic,mi)%po4(ii) = f(mi)%po4*bv_mean*z1 + ice(sc,ic,mi)%po4(ii)*z2
              ice(sc,ic,mi)%sioh4(ii) = f(mi)%sioh4*bv_mean*z1 + ice(sc,ic,mi)%sioh4(ii)*z2

          enddo

          r_depth = ice(sc,ic,mi)%id(int_z)
          ice(sc,ic,mi)%z = sk_z

#include "sia2_env_ice_grid.inc.f90"

          if (maxed_depth) then
             print *,'mi: ',mi,'maxed depth in sia2_env_redist_new_ice.inc.F - weird, should be smallest category'
             maxed_depth = .false.
          endif

          ! start new layer top at old layer top - don't exclude any ice
          new_layer_top = -1.
          ! no growth/melt here
          s_gl = 0.
          c_gl = 0.
          ! no flooding here
          flooded = 0.
          melt_flood = 0.
          dhdt_cm_s = 0.   ! required in algal migration

#include "sia2_env_ice_remap.inc.f90"

          ! scale snow by area, reset ts
          dz = ice(sc,ic,mi)%af/af_total
          ice(sc,ic,mi)%sh_prev = ice(sc,ic,mi)%sh_prev*dz

          if (ice(sc,ic,mi)%snow%z .gt. 0) then

              do ii=1,ice(sc,ic,mi)%snow%z
                  ice(sc,ic,mi)%snow%th(ii) = ice(sc,ic,mi)%snow%th(ii)*dz
              enddo
              ice(sc,ic,mi)%snow%depth = ice(sc,ic,mi)%snow%depth*dz
              r_depth = ice(sc,ic,mi)%sh_prev

#include "sia2_env_snow_grid.inc.f90"

              flooded = 0.d0
              melt_drain = 0.d0

#include "sia2_env_snow_remap.inc.f90"

          else

              ! average ts over area if no snow
              t_mean = min(airtemp1c,f(mi)%t,mu*ice(sc,ic,mi)%s(1))
              ice(sc,ic,mi)%snow%ts = (t_mean*adv(mi)%af_new &
                  + ice(sc,ic,mi)%snow%ts*ice(sc,ic,mi)%af)/af_total

          endif

          ! update pond, age, ridged, snowdist tracers
          dz = ice(sc,ic,mi)%af/af_total
          ice(sc,ic,mi)%pond%th = ice(sc,ic,mi)%pond%th*dz
          ice(sc,ic,mi)%pond%perc = ice(sc,ic,mi)%pond%perc*dz
          ice(sc,ic,mi)%age = ice(sc,ic,mi)%age*dz
          ice(sc,ic,mi)%ridged = ice(sc,ic,mi)%ridged*dz
          ice(sc,ic,mi)%snow_dist =  &
              (dble(sd_num)*0.5d0*adv(mi)%af_new + &
               ice(sc,ic,mi)%snow_dist*ice(sc,ic,mi)%af) / af_total
          ice(sc,ic,mi)%snow_rd = 0.d0

					! update ml layer tracers
					ice(sc,ic,mi)%no3(ml_z) = (f(mi)%no3*adv(mi)%af_new + &
						ice(sc,ic,mi)%no3(ml_z)*ice(sc,ic,mi)%af)/af_total
					ice(sc,ic,mi)%nh4(ml_z) = (f(mi)%nh4*adv(mi)%af_new + &
						ice(sc,ic,mi)%nh4(ml_z)*ice(sc,ic,mi)%af)/af_total
					ice(sc,ic,mi)%po4(ml_z) = (f(mi)%po4*adv(mi)%af_new + &
						ice(sc,ic,mi)%po4(ml_z)*ice(sc,ic,mi)%af)/af_total
					ice(sc,ic,mi)%sioh4(ml_z) = (f(mi)%sioh4*adv(mi)%af_new + &
						ice(sc,ic,mi)%sioh4(ml_z)*ice(sc,ic,mi)%af)/af_total

          ! update areal fraction to include newly created ice
          ice(sc,ic,mi)%af = af_total

					! find ending fresh and salt mass
					call sia2_ice_mass(ice(sc,ic,mi),fresh_tmp,dsdt_tmp)

					! record fresh water flux (kg)
					fresh_net = (fresh_tmp*ice(sc,ic,mi)%af-fresh_net)
					m(mi)%h2o_flux = m(mi)%h2o_flux - fresh_net*cell_area*1.e6

					! record salt flux (kg)
					dsdt_net = (dsdt_tmp*ice(sc,ic,mi)%af-dsdt_net)
					m(mi)%salt_flux = m(mi)%salt_flux - dsdt_net*cell_area*1.e6

      endif

! ======================================================================
! END: Add new ice
! ======================================================================
