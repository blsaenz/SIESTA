
		   ! Adjust Ice Layers
		   ! -----------------------------------------------------------
		   if (((abs(c_gl) + abs(s_gl) + abs(crunch_tot)) .ge. gl_max) &
		     .or. (flooded .ge. fl_max) &
		     .or. (iii .eq. dtt-1)) then

			   ! find new thickness/icedepth structure with th_new and id_new
			   r_depth = ice(sc,ic,mi)%id(int_z) + c_gl + s_gl
			   if (r_depth .lt. h_min) then
				   r_depth = ice(sc,ic,mi)%id(int_z) + c_gl
				   s_gl = 0.d0
				   if (r_depth .lt. h_min) then
					   r_depth = -1.
				   endif
               endif

	           ! grow/melt ice column resulting from heat flux
			   ! check to make sure we don't grow below minimum
			   if (r_depth .ge. h_min) then

				   ! record bm_lost and freshwater input for melting
				   if (c_gl .lt. 0.d0) then
					   dz_mean = abs(c_gl)
					   ii = ice(sc,ic,mi)%z
					   do while((dz_mean .gt. 0.d0) .and. (ii .gt. 0))
						   smalg_lost = ice(sc,ic,mi)%smalg(ii) - min_alg
						   smalg_lost = max(0.,smalg_lost)
						   if (dz_mean .gt. ice(sc,ic,mi)%th(ii)) then
							  ! converting to g/pixel - (cell_area*1.e6)*1g/1000mg = 1000. scale factor
							  m(mi)%bm_lost =  m(mi)%bm_lost + smalg_lost* &
							  ice(sc,ic,mi)%th(ii)*cell_area*1.e3*ice(sc,ic,mi)%af   ! g/pixel
							  dz_mean = dz_mean - ice(sc,ic,mi)%th(ii)
							  ii = ii-1
						   else
							  ! converting to g/pixel - (cell_area*1.e6)*1g/1000mg = 1000. scale factor
							  m(mi)%bm_lost =  m(mi)%bm_lost + smalg_lost* &
							  dz_mean*cell_area*1.e3*ice(sc,ic,mi)%af   ! g/pixel
							  dz_mean = 0.d0
						   endif
					   enddo

                       ! for salt accounting, because the pesky-constant salt skeletal layer has a minimum, we take salt and freshwater fluxes
                       ! from interior layers
					   dz_mean = abs(c_gl)
					   ii = ice(sc,ic,mi)%z - z_sk
					   do while((dz_mean .gt. 0.d0) .and. (ii .gt. 0))
						   if (dz_mean .gt. ice(sc,ic,mi)%th(ii)) then
							  ! converting to kg/m^2 freshwater input 1000kg/m3 * m * fraction_water by weight
							  !m(mi)%fresh_out = m(mi)%fresh_out + ice(sc,ic,mi)%th(ii)*ice(sc,ic,mi)%af*1.d3*(1.d0-ice(sc,ic,mi)%s(ii)*c_001)
							  ! converting to kg/m^2 salt input 1000kg/m3 * m * fraction_salt by weight
							  !m(mi)%dsdt_out = m(mi)%dsdt_out + ice(sc,ic,mi)%th(ii)*ice(sc,ic,mi)%af*ice(sc,ic,mi)%s(ii)

							  dz_mean = dz_mean - ice(sc,ic,mi)%th(ii)
							  ii = ii-1
						   else
							  ! converting to kg/m^2 freshwater input 1000kg/m3 * m * fraction_water by weight
							  !m(mi)%fresh_out = m(mi)%fresh_out + dz_mean*ice(sc,ic,mi)%af*1.d3*(1.d0-ice(sc,ic,mi)%s(ii)*c_001)
							  ! converting to kg/m^2 salt input 1000kg/m3 * m * fraction_salt by weight
							  !m(mi)%dsdt_out = m(mi)%dsdt_out + dz_mean*ice(sc,ic,mi)%af*ice(sc,ic,mi)%s(ii)
							  dz_mean = 0.d0
						   endif
					   enddo


				   endif


                   ! record freshwater input for surface melting - done in meltwater percolation routine


                   ! Determine surface flooding - only one type of surface flooding
                   ! may occur at once (seawater flooding or meltwater flooding)
                   ! and seawater flooding takes precedence

	               melt_flood = 0.d0
	               melt_drain = 0.d0
                   z_snow = ice(sc,ic,mi)%snow%z

!				   if ((flooded .ge. fl_max .and. ice(sc,ic,mi)%bv(1) .ge. vb_crit) &
!				   .or. (flooded .ge. fl_max*2.)) then
				   if (flooded .ge. fl_max) then

					  did_flood = .true.
                      flood_1x = .false.

                      tmp1 = max(IceD/2.d0,ice(sc,ic,mi)%snow%d(1))

                      f_depth = flooded*ice(sc,ic,mi)%snow%d(1)/tmp1

					  r_depth = r_depth + f_depth

					  ! record ice volume change due to snow-ice gain
					  m(mi)%snow_growth = m(mi)%snow_growth + f_depth* &
						  ice(sc,ic,mi)%af*cell_area*1.e6

					  ! Conservation of heat: mix snow and seawater
					  ! assuming mass/volume changes are negligable inside a layer initially...
					  !top_f = snow_mass/snowh1m/IceD  ! volume fraction of fresh ice in snow

					  top_f = tmp1/IceD ! snow/water volume ratio
					  bot_f = 1.-top_f ! water/snow volume ratio

					  s_mean = max(min_sal,f(mi)%s*bot_f)
					  d_mean = f(mi)%d*bot_f + IceD*top_f

					 ! snow heat is assumed to be just latent heat of melting + 0.02*t*d
					  heat_mean = cw*f(mi)%t*f(mi)%d*bot_f - &
							(ice(sc,ic,mi)%snow%t(z_snow)*0.02 + Lf)*IceD*top_f
					  dz_mean = 1.

#include "sia2_env_temp_from_heat.inc.f90"

                      t_mean = min(s_mean*mu*1.01,t_mean)
                      !t_mean = ice(sc,ic,mi)%t(1)


!					  S_mean = S_flood
!					  T_mean = f(mi)%t

#include "sia2_env_brine_calc.inc.f90"

					  T_flood = t_mean
					  S_flood = s_mean
					  d_flood = d_mean
					  bs_flood = bs_mean
					  bd_flood = bd_mean
					  bv_flood = bv_mean
					  heat_flood = heat_mean

                  elseif (use_drain .eq. 1) then

                       jj=0
                       z_snow_new = 0
                       do ii=1,ice(sc,ic,mi)%z
                          if (ice(sc,ic,mi)%bv(ii) .gt. bv_conv .and. &
                              ice(sc,ic,mi)%id(ii) .le. ice(sc,ic,mi)%fbh) then
                              jj = ii
                          endif
                       enddo

                       ! drain ice layers above and eq. to jj of brine, turn it into snow
                       if (jj .gt. 0) then
                           did_drain = .true.
                           th_snow_new = 0.d0
                           d_snow_new = 0.d0
                           heat_snow_new = 0.d0
                           z_snow_new = jj
                           do ii=1,jj
                               melt_drain = melt_drain + ice(sc,ic,mi)%th(ii)
                               th_snow_new(ii) = ice(sc,ic,mi)%th(ii)
                               d_snow_new(ii) = (1.-ice(sc,ic,mi)%bv(jj)*c_001)*IceD
                               T_mean = ice(sc,ic,mi)%t(ii) + kelvin0
                               heat_snow_new(ii) = d_snow_new(ii)*(heat_snow0 - 0.2309*T_mean - &
                                     0.0034*T_mean**2)
                           enddo


                           ! update s_gl to drop layer from ice (and move to snow)
                           s_gl = s_gl - melt_drain

                           r_depth = r_depth - melt_drain

                       endif ! end of is there draining?
                  ! superimposed ice formation
!                  elseif (ice(sc,ic,mi)%pond%th .gt. 0. .and. z_snow .gt. 0 .and. &
!                      ice(sc,ic,mi)%bv(1) .lt. bv_conv) then
!
!					  ! make snow-ice
!					  did_flood = .true.
!					  flood_1x = .false.
!
!					  top_f = ice(sc,ic,mi)%snow%d(z_snow)/IceD ! snow/water volume ratio
!					  bot_f = 1.-top_f ! water/snow volume ratio
!
!					  ! find height of new snow ice based on flooding of snow
!
!					  melt_flood = ice(sc,ic,mi)%pond%th/bot_f
!					  melt_flood = min(melt_flood,ice(sc,ic,mi)%snow%depth)
!					  dz = melt_flood*bot_f
!
!					  r_depth = r_depth + melt_flood
!
!					  ! record ice volume change due to snow-ice gain
!					  m(mi)%snow_growth = m(mi)%snow_growth + melt_flood* &
!						  ice(sc,ic,mi)%af*cell_area*1.e6
!
!					  S_mean = max(min_sal,ice(sc,ic,mi)%pond%s*bot_f)
!					  d_mean = ice(sc,ic,mi)%pond%d*bot_f + IceD*top_f
!					 ! snow heat is assumed to be just latent heat of melting + 0.02*t*d
!					  heat_mean = cw*ice(sc,ic,mi)%pond%t*ice(sc,ic,mi)%pond%d*bot_f - &
!							(ice(sc,ic,mi)%snow%t(z_snow)*0.02 + Lf)*IceD*top_f
!					  dz_mean = 1.
!
!#include "sia2_env_temp_from_heat.inc.f90"
!
!                      t_mean = min(s_mean*mu*1.01,t_mean)
!
!#include "sia2_env_brine_calc.inc.f90"
!
!					  T_flood = t_mean
!					  S_flood = s_mean
!					  d_flood = d_mean
!					  bs_flood = bs_mean
!					  bd_flood = bd_mean
!					  bv_flood = bv_mean
!					  heat_flood = heat_mean
!
!					  ! store pond tracers for use in flooding routine later
!					  smalg_pond_new = ice(sc,ic,mi)%pond%smalg
!					  poc_pond_new = ice(sc,ic,mi)%pond%poc
!					  no3_pond_new = ice(sc,ic,mi)%pond%no3
!					  nh4_pond_new = ice(sc,ic,mi)%pond%nh4
!					  po4_pond_new = ice(sc,ic,mi)%pond%po4
!					  sioh4_pond_new = ice(sc,ic,mi)%pond%sioh4
!
!                     ice(sc,ic,mi)%pond%th = ice(sc,ic,mi)%pond%th - dz
!                     if (ice(sc,ic,mi)%pond%th .le. 0.) then
!						  ice(sc,ic,mi)%pond%th = 0.d0
!						  ice(sc,ic,mi)%pond%t = 0.d0
!						  ice(sc,ic,mi)%pond%s = 0.d0
!						  ice(sc,ic,mi)%pond%d = 0.d0
!						  ice(sc,ic,mi)%pond%heat = 0.d0
!						  ice(sc,ic,mi)%pond%smalg = 0.d0     ! brine-based algal concentration (mgC/m^3)
!						  ice(sc,ic,mi)%pond%poc = 0.d0     ! particulate organic carbon (detritus) (mgC/m^3)
!						  ice(sc,ic,mi)%pond%no3 = 0.d0        ! ice layer brine NO3 concentration (µMol)
!						  ice(sc,ic,mi)%pond%nh4 = 0.d0      ! ice layer brine NH4 concentration (µMol)
!						  ice(sc,ic,mi)%pond%po4 = 0.d0     ! ice layer brine PO4 concentration (µMol)
!						  ice(sc,ic,mi)%pond%sioh4 = 0.d0    ! ice layer brine SiOH4 concentration (µMol)
!                     endif

                   endif

#include "sia2_env_ice_grid.inc.f90"


                   snow_melt_now = 0.d0
                   if (use_ponds .eq. 1 .and. s_gl .lt. 0.d0) then

#include "sia2_env_pond_update.inc.f90"

                   endif

                   ! if surface melt, use this assigment get rid of surface ice
			       new_layer_top = abs(s_gl)

#include "sia2_env_ice_remap.inc.f90"

                   if (melt_drain .gt. 0.d0) then

                      ! pond_percolate requires th_new to be assinged to brine volumes!
                      do ii=1,sk_z
                          ! collapse layers to bv thickness for pond drainage
                          tmp1 = ice(sc,ic,mi)%bv(ii)*c_001
                          th_new(ii) = ice(sc,ic,mi)%th(ii)*tmp1
                      enddo

#include "sia2_env_pond_percolate3.inc.f90"

                      ice(sc,ic,mi)%pond%th = 0.d0
                      ice(sc,ic,mi)%pond%t = 0.d0
                      ice(sc,ic,mi)%pond%s = 0.d0
                      ice(sc,ic,mi)%pond%d = 0.d0
                      ice(sc,ic,mi)%pond%heat = 0.d0
                      ice(sc,ic,mi)%pond%smalg = 0.d0     ! brine-based algal concentration (mgC/m^3)
                      ice(sc,ic,mi)%pond%poc = 0.d0     ! particulate organic carbon (detritus) (mgC/m^3)
                      ice(sc,ic,mi)%pond%no3 = 0.d0        ! ice layer brine NO3 concentration (µMol)
                      ice(sc,ic,mi)%pond%nh4 = 0.d0      ! ice layer brine NH4 concentration (µMol)
                      ice(sc,ic,mi)%pond%po4 = 0.d0     ! ice layer brine PO4 concentration (µMol)
                      ice(sc,ic,mi)%pond%sioh4 = 0.d0    ! ice layer brine SiOH4 concentration (µMol)

                   endif

			  ! end of check to make sure we don't grow below minimum
			  endif

			  ! record total ice growth
              if (c_gl .gt. 0.d0) then
				   m(mi)%cong_growth = m(mi)%cong_growth + &
					 c_gl*ice(sc,ic,mi)%af*cell_area*1.e6
			  else
				   m(mi)%melt_loss = m(mi)%melt_loss + &
					 c_gl*ice(sc,ic,mi)%af*cell_area*1.e6
			  endif
              if (s_gl .gt. 0.d0) then
				   m(mi)%cong_growth = m(mi)%cong_growth + &
					 s_gl*ice(sc,ic,mi)%af*cell_area*1.e6
			  else
				   m(mi)%melt_loss = m(mi)%melt_loss + &
					 s_gl*ice(sc,ic,mi)%af*cell_area*1.e6
			  endif

			  ! record growth/melt totals for timestep
			  c_gl_tot = c_gl_tot + c_gl ! congelation growth/melt
			  s_gl_tot = s_gl_tot + s_gl ! surface growth/melt

			  ! reset current recorders to zero
			  c_gl = 0.
			  s_gl = 0.
			  crunch_tot = 0.

			  ! Adjust Skeletal/bottom layer growth/melt
			  ! -----------------------------------------------------------

			  if (sk_gl .ne. 0.d0) then
				   ! grow/melt skeletal layer
				   dz_sk = sk_gl/DBLE(z_sk)

				   if (dz_sk .lt. 0.d0) then
                      do ii=sk_1,sk_z
						  bv_mean = ice(sc,ic,mi)%bv(ii)*c_001
						  smalg_lost = ice(sc,ic,mi)%smalg(ii) - min_alg
						  smalg_lost = max(0.,smalg_lost)

						  ! biomass lost to ocean - converting to g/pixel - (cell_area*1.e6)*1g/1000mg = 1000. scale factor
						  m(mi)%bm_lost = m(mi)%bm_lost + smalg_lost*1.e3* &
						  abs(dz_sk)*cell_area*ice(sc,ic,mi)%af   ! g/pixel

                          ! freshwater and salt accounting
                          ! converting to kg/m^2 freshwater input 1000kg/m3 * m * fraction_water by weight
                          !m(mi)%fresh_out = m(mi)%fresh_out + abs(dz_sk)*ice(sc,ic,mi)%af*1.d3*(1.d0-ice(sc,ic,mi)%s(ii)*c_001)
                          ! converting to kg/m^2 salt input 1000kg/m3 * m * fraction_salt by weight
                          !m(mi)%dsdt_out = m(mi)%dsdt_out + abs(dz_sk)*ice(sc,ic,mi)%af*ice(sc,ic,mi)%s(ii)

					  enddo
				  endif

#include "sia2_env_update_skeletal.inc.f90"

                   ! reset sk_gl
				   sk_gl = 0.d0
			   endif
		  endif


		  ! Adjust Snow Layers
		  ! -----------------------------------------------------------
          ! find if projected snowdepth is above minimum
		  valid_snowh = (snowh1m + snow_gl) .ge. snow_min

		  ! add fresh water to melt water pond
		  if ((use_ponds .eq. 1) .and. (snow_melt .gt. 0.) .and. &
		  (z_snow .gt. 0)) then

			   snow_melt_now = snow_melt
			   new_layer_top = 0.

#include "sia2_env_pond_update.inc.f90"

			   snow_melt = 0.

		  endif

          if (dt_step .eq. 211 .and. iii .eq. 5) then
          testvar = -1
          endif

		  if (((z_snow .eq. 0) .and. valid_snowh) .or. &
		  ((z_snow .gt. 0) .and. .not. valid_snowh) .or. &
		  (abs(snow_gl) .gt. gl_max) .or. did_flood .or. did_drain .or. (iii .eq. dtt-1)) then
			  ! calc new total snow height, accounting for snowfall, melt, surface flooding
			  snowh1m = max(snowh1m + snow_gl, 0.)

			  if (did_flood) then
                  if (boundary_file .eq. 1) then
                      ice(sc,ic,mi)%sh_offset = ice(sc,ic,mi)%sh_offset + flooded
                  endif
				  snowh1m = max(snowh1m - flooded, 0.)

              elseif (did_drain) then

                  ! draining surface ice layer and turning to snow
			      if(melt_drain .gt. 0.d0) then
                      do ii = 1,z_snow_new
                          ! make sure that drained layer still has some ice left in it
                          if (d_snow_new(ii) .gt. 0.d0) then
                              snowh1m = snowh1m + th_snow_new(ii) ! add drained ice as snow
                              if (boundary_file .eq. 1) then
                                  ice(sc,ic,mi)%sh_offset = ice(sc,ic,mi)%sh_offset + th_snow_new(ii)
                              endif
                          endif
                      enddo
                  endif

			  endif

!              if (snow_gl .lt. 0.) then
!                  ice(sc,ic,mi)%snow%d(z_snow) = ice(sc,ic,mi)%snow%d(z_snow)*ice(sc,ic,mi)%snow%th(z_snow)/ &
!                      (ice(sc,ic,mi)%snow%th(z_snow)-snow_gl)
!              endif

              r_depth = snowh1m

#include "sia2_env_snow_grid.inc.f90"

#include "sia2_env_snow_remap.inc.f90"

              snow_gl = 0.
			  if (did_flood) then
				  did_flood = .false.
                  flooded = 0.
              endif
              if (did_drain) then
				  did_drain = .false.
                  melt_drain = 0.d0
              endif
		  endif
