! sia_adv_ridge_70.inc.F
! ---------------------------------------------------------------------
! This file does the actual ridging between ice categories, but reduces
! the areal ridging allocatation by 30% and subsitutes in 30% seawater
!
! The routine randomly chooses whether to put the 30% height seawater
! near the top, or near the middle of ridged ice

	! create copy of grid cell
			do ic=1,ida_n
			do sc=1,sda_n
				ice_ridged(sc,ic,thread_n) = ice(sc,ic,mi)
			enddo
			enddo
			ice_temp2 = ice(1,1,mi)

			if(ff(1) .lt. 0.5) then
			testvar = 1
			endif


	! if (ida_n .gt. 1) then

			! perform inter-category ridging up to largest ice category
			do ic=1,ida_n
			do sc=1,sda_n

					s_new = 0.d0
					d_new = 0.d0
					heat_new = 0.d0
					smalg_new = 0.d0
					poc_new = 0.d0
					no3_new = 0.d0
					nh4_new = 0.d0
					po4_new = 0.d0
					sioh4_new = 0.d0
					age_new = 0.d0
					ridged_new = 0.d0
					snow_dist_new = 0.d0
					snow_rd_new = 0.d0

					z_snow = 0
					d_new_snow = 0.d0
					heat_new_snow = 0.d0
					melt_new = 0.d0
					th_new_snow = 0.d0 ! have to track snow mass height inputs, since it doesn't 'ridge' like ice
					! find density of snowfall - may be needed if, through advection,
					! the snowdepth crosses the min threshold and is initialized
					if ((f(mi)%at-kelvin0) .gt. -2.) then
							snowfall_d = den_s_wet*1.e6
					else
							snowfall_d = den_s_dry*1.e6
					endif

					s_pond_new = 0.d0
					h_pond_new = 0.d0
					d_pond_new = 0.d0
					smalg_pond_new = 0.d0
					poc_pond_new = 0.d0
					no3_pond_new = 0.d0
					nh4_pond_new = 0.d0
					po4_pond_new = 0.d0
					sioh4_pond_new = 0.d0
					th_pond_new = 0.d0

					sh_prev_sum = 0.d0
					sh_offset_sum = 0.d0

					! find ice area in category subtracting ridge-participating ice
					adv(mi)%ice(sc,ic)%af = cell_area*ice(sc,ic,mi)%af
					if ((adv(mi)%a_rdg(ic) .gt. 0.d0) .and. (ice(sc,ic,mi)%z .gt. 0)) then
							adv(mi)%ice(sc,ic)%af = adv(mi)%ice(sc,ic)%af - &
								adv(mi)%a_rdg(ic)*ice(sc,ic,mi)%af/af_ic(ic)
							r_depth = ice(sc,ic,mi)%id(ice(sc,ic,mi)%z-z_sk)
					else
							r_depth = ic_h_med(ic)
					endif
					a_scale_sum = 0.d0

					! determine height scaling of grid used for ridging category

					r_depth = r_depth*adv(mi)%ice(sc,ic)%af  ! store area-weighted ice thickness sum
					hm = adv(mi)%ice(sc,ic)%af ! store sum of areas
					!if (ic .gt. 1) then
					do jjj=1,ida_n
						if (sda_n .eq. 1 .or. sc .eq. 2) then
							do sc1 = 1,sda_n
								if (adv(mi)%rf(ic,jjj) .gt. 0.d0 .and. ice(sc1,jjj,mi)%af .gt. 0.d0) then
									tmp1 = adv(mi)%a_rdg(jjj) * adv(mi)%rf(ic,jjj) * &
										ice(sc1,jjj,mi)%af/af_ic(jjj) ! area ridging in (before ridging)
									if (jjj .eq. ic .and. ic .eq. ida_n) then
										! self-ridging -- need to take into account ff (which was probably
										! set arbitrarily to accomodate ridging) to find new category height
										! ---------------------------------------------------------
										h_scale = ice(sc1,jjj,mi)%id(ice(sc1,jjj,mi)%z-z_sk) / ff(jjj)
										tmp1 = tmp1 * ff(jjj) ! area after ridging
    											! adjust cv_void_f
    											void_f = cv_void_f
    											if (id_ic(jjj) .lt. cvf_switch) then
    												void_f = void_f*cvf_thin
    											endif
										tmp2 = h_scale/(1.d0-void_f) ! self-ridging height
									else
										! adding to total areal height in ridged cell, knowing that
										! ice from other categories is ridged in @ ic_h_med height
										! ----------------------------------------------------
    											! adjust cv_void_f
    											void_f = cv_void_f
    											if (id_ic(jjj) .lt. cvf_switch) then
    												void_f = void_f*cvf_thin
    											endif
										h_scale = ice(sc1,jjj,mi)%id(ice(sc1,jjj,mi)%z-z_sk) / &
											ic_h_med(ic)*(1.d0-void_f)
										tmp1 = tmp1*h_scale  ! area after ridging
										tmp2 = ic_h_med(ic) ! ridging height
									endif
									r_depth = r_depth + tmp1*tmp2
									hm = hm + tmp1
								endif
							enddo
						endif
					enddo
					!endif
					if (hm .gt. 0.d0) then
						r_depth = r_depth/hm ! divide to find new ice thickness
						hm = r_depth ! assign to hm for calcs below
					else
						r_depth = ic_h_med(ic)
						hm = ic_h_med(ic)
					endif

					! We may ridge higher than h_max temporarily
					if (ic .eq. ida_n) ignore_h_max = .true.

#include "sia2_env_ice_grid.inc.f90"

          ignore_h_max = .false.

          sk_1_rdg = int_z_new + 1

          ! record tracers that remain in layer (are not a part of ice that ridges)
          if (ice(sc,ic,mi)%af .gt. 0.d0 .and. apn(ic) .ge. 0.d0) then

				  h_scale = ice(sc,ic,mi)%id(ice(sc,ic,mi)%z-z_sk)/hm
				  a_scale = adv(mi)%ice(sc,ic)%af*h_scale
				  a_scale_sum = a_scale_sum + a_scale
				  th_old = ice(sc,ic,mi)%th/h_scale ! vector multiply = streching layer, so account for ridging below
				  sk_z = ice(sc,ic,mi)%z
				  int_z = sk_z - z_sk
				  sk_1 = int_z+1

				  new_layer_bot = 0.

				  ! add incoming internal ice layers
				  do ii=1,int_z_new

					  ! re-adjust new layer bottom for next snow layer temp adjustment
					  new_layer_top = new_layer_bot
					  new_layer_bot = new_layer_bot + th_new(ii) ! th_new comes from snow_grid calc

					  ! find 1st old layer that contains part of new layer, going down
					  jj = 0
					  old_layer_bot = 0.
					  do while((new_layer_top .ge. old_layer_bot) .and. (jj .lt. sk_z))
						   jj=jj+1
						   ! old thickness...
						   old_layer_bot = old_layer_bot + th_old(jj)
					   enddo
					   old_layer_top = old_layer_bot - th_old(jj)
					   ! now jj is OLD layer where NEW layer ii starts...

					   ! find total heat/salt from multiple layers/partial layers that make up
					   ! the new layer
					   do while ((old_layer_top .lt. new_layer_bot) .and. (jj .le. sk_z))

						   ! ----- NEW LAYER GEOMETRIES ------------------------------
						   interim_top = max(old_layer_top,new_layer_top)
						   interim_bot = min(old_layer_bot,new_layer_bot)

						   z1 = interim_top - old_layer_top
						   z2 = interim_bot - old_layer_top
						   dz = z2-z1

						   ! record tracers
						   tmp1 = dz*a_scale
						   s_new(ii) = s_new(ii) + ice(sc,ic,mi)%s(jj)*tmp1
						   d_new(ii) = d_new(ii) + ice(sc,ic,mi)%d(jj)*tmp1
						   heat_new(ii) = heat_new(ii) + ice(sc,ic,mi)%heat(jj)*tmp1
						   no3_new(ii) = no3_new(ii) + ice(sc,ic,mi)%no3(jj)*tmp1
						   nh4_new(ii) = nh4_new(ii) + ice(sc,ic,mi)%nh4(jj)*tmp1
						   po4_new(ii) = po4_new(ii) + ice(sc,ic,mi)%po4(jj)*tmp1
						   sioh4_new(ii) = sioh4_new(ii) + ice(sc,ic,mi)%sioh4(jj)*tmp1
						   smalg_new(ii) = smalg_new(ii) + ice(sc,ic,mi)%smalg(jj)*tmp1
						   poc_new(ii) = poc_new(ii) + ice(sc,ic,mi)%poc(jj)*tmp1
						   ! ----- SETUP VARIABLES FOR NEXT PARTIAL LAYER ------------
						   ! keeping track of emerging new layer thickness, for
						   dz_total = dz_total + dz

						   jj=jj+1
						   if (jj .le. sk_z) then
							   ! find boundaries of next old layer
							   old_layer_top = old_layer_bot
							   old_layer_bot = old_layer_bot + th_old(jj)
						   endif

					   enddo  ! end of layer ii mapping

				  enddo ! end of ii layer loop for internal ice

				  ! add incoming skeletal ice layers
				  do ii=1,z_sk

					  ! find corresponding skeletal layers
					  jj = ice(sc,ic,mi)%z-z_sk+ii ! ic skeletal layer
					  kk = sk_1_rdg - 1 + ii

					  tmp1 = adv(mi)%ice(sc,ic)%af*ice(sc,ic,mi)%th(jj)

					  th_new(kk) = th_new(kk) + tmp1
					  s_new(kk) = s_new(kk) + ice(sc,ic,mi)%s(jj)*tmp1
					  d_new(kk) = d_new(kk) + ice(sc,ic,mi)%d(jj)*tmp1
					  heat_new(kk) = heat_new(kk) + ice(sc,ic,mi)%heat(jj)*tmp1
					  smalg_new(kk) = smalg_new(kk) + ice(sc,ic,mi)%smalg(jj)*tmp1
					  poc_new(kk) = poc_new(kk) + ice(sc,ic,mi)%poc(jj)*tmp1
					  no3_new(kk) = no3_new(kk) + ice(sc,ic,mi)%no3(jj)*tmp1
					  nh4_new(kk) = nh4_new(kk) + ice(sc,ic,mi)%nh4(jj)*tmp1
					  po4_new(kk) = po4_new(kk) + ice(sc,ic,mi)%po4(jj)*tmp1
					  sioh4_new(kk) = sioh4_new(kk) + ice(sc,ic,mi)%sioh4(jj)*tmp1

				  enddo

				  ! account for snow during ridging
				  tmp1 = adv(mi)%ice(sc,ic)%af

				  if (ice(sc,ic,mi)%snow%z .gt. 0) then

					  z_snow = max(z_snow,ice(sc,ic,mi)%snow%z)
					  do ii=1, ice(sc,ic,mi)%snow%z
						  tmp2 = tmp1*ice(sc,ic,mi)%snow%th(ii)
						  th_new_snow(ii) = th_new_snow(ii) + tmp2
						  d_new_snow(ii) = d_new_snow(ii) + ice(sc,ic,mi)%snow%d(ii)*tmp2
						  heat_new_snow(ii) = heat_new_snow(ii) + ice(sc,ic,mi)%snow%heat(ii)*tmp2
						  melt_new(ii) = melt_new(ii) + ice(sc,ic,mi)%snow%melt(ii)*tmp2
					  enddo

				  else

					  sh_prev_sum = sh_prev_sum + ice(sc,ic,mi)%sh_prev*tmp1

				  endif  ! end of snow input check

				  ice(sc,ic,mi)%snow%z = z_snow
				  sh_offset_sum = sh_offset_sum + ice(sc,ic,mi)%sh_offset*tmp1

				  ! update ml layer tracers
				  no3_new(ml_z) = no3_new(ml_z) + &
				  	ice(sc,ic,mi)%no3(ml_z)*adv(mi)%ice(sc,ic)%af
				  nh4_new(ml_z) = nh4_new(ml_z) + &
				  	ice(sc,ic,mi)%nh4(ml_z)*adv(mi)%ice(sc,ic)%af
				  po4_new(ml_z) = po4_new(ml_z) + &
				  	ice(sc,ic,mi)%po4(ml_z)*adv(mi)%ice(sc,ic)%af
				  sioh4_new(ml_z) = sioh4_new(ml_z) + &
				  	ice(sc,ic,mi)%sioh4(ml_z)*adv(mi)%ice(sc,ic)%af

				  ! update age and ridged tracers from non-ridging ice
				  age_new = age_new + ice(sc,ic,mi)%age*adv(mi)%ice(sc,ic)%af
				  ridged_new = ridged_new + ice(sc,ic,mi)%ridged*adv(mi)%ice(sc,ic)%af
				  snow_dist_new = snow_dist_new + ice(sc,ic,mi)%snow_dist*adv(mi)%ice(sc,ic)%af
				  snow_rd_new = snow_rd_new + ice(sc,ic,mi)%snow_rd*adv(mi)%ice(sc,ic)%af



              endif

              ! perhaps change this restriction later, if ridging is allowed into the first category
              !if (ic .gt. 1) then  ! commented out, so that we can have voids with ic==1

                  do jjj=1,ida_n
                      if (sda_n .eq. 1 .or. sc .eq. 2) then

                      do sc1 = 1,sda_n

                      if (adv(mi)%rf(ic,jjj) .gt. 0.d0 .and. ice(sc1,jjj,mi)%af .gt. 0.d0) then

    											! adjust cv_void_f
    											void_f = cv_void_f
    											if (id_ic(jjj) .lt. cvf_switch) then
    												void_f = void_f*cvf_thin
    											endif

                          ! -----------------------------------------------
                          ! modified h_scale, a_scale, th_old and to accomodate reduced ridging effectiveness
                          ! and gap layers
                          tmp1 = adv(mi)%a_rdg(jjj)* &
                          	ice(sc1,jjj,mi)%af/af_ic(jjj)* &
                          	adv(mi)%rf(ic,jjj) ! pre-ridging area input
                          h_scale = &
                          	ice(sc1,jjj,mi)%id(ice(sc1,jjj,mi)%z-z_sk) / &
                          		(hm * (1.d0-void_f))
                          a_scale = tmp1 * h_scale
                          a_scale_sum = a_scale_sum + a_scale

                          th_old = ice(sc1,jjj,mi)%th/h_scale ! vector multiply = streching layer, so account for ridging below

													! determine ridged area post-ridging
													if(ic .eq. jjj .and. ic .eq. ida_n) then
														a_rdg_rdg = tmp1*ff(jjj) ! area*thickness/ridging_thickness
													else
                          	a_rdg_rdg = tmp1 * &
                          		ice(sc1,jjj,mi)%id(ice(sc1,jjj,mi)%z-z_sk) / &
                          		(ic_h_med(ic)*(1.d0-void_f))
                          endif

                          ! -----------------------------------------------
                          sk_z = ice(sc1,jjj,mi)%z
                          int_z = sk_z - z_sk
                          sk_1 = int_z+1

                          ! add incoming area
                          adv(mi)%ice(sc,ic)%af = adv(mi)%ice(sc,ic)%af + a_rdg_rdg

                          new_layer_bot = 0.

                          if (void_f .gt. 0.d0) then

                              ! --------------------------------------------
                              ! gen random number that deterimines where seawater pockets will exist
                              ! vertically in ice pack
                              ! --------------------------------------------
                              call random_number(id_insert)
                              if(id_insert .lt. 0.666666) then
                                  if(id_insert .lt. 0.333333) then
                                      id_insert = 0.25
                                  else
                                      id_insert = 0.5
                                  endif
                              else
                                  id_insert = 0.75
                              endif
                              id_insert = id_insert*hm*(1.d0-void_f) ! height of insertion
                              id_sw_gap = hm*void_f ! thickness of seawater gap insertion

                              ! record void volume
                          		m(mi)%void_growth = m(mi)%void_growth + id_sw_gap*a_rdg_rdg*1.e6

                              ! --------------------------------------------

															! find snow for possible snow in gap
															tmp2 = 0.d0 ! density*thickness
															tmp1 = 0.d0 ! temp*thickness*density  !  -- taken out - seemed not really needed since a hack anyway
															if((snow_in_gaps .eq. 1) .and. &
																(ice(sc1,jjj,mi)%snow%z) .gt. 0) then
																do jj=1,ice(sc1,jjj,mi)%snow%z
																	tmp2 = tmp2 + ice(sc1,jjj,mi)%snow%th(jj) * &
																		ice(sc1,jjj,mi)%snow%d(jj)
																	tmp1 = tmp1 + ice(sc1,jjj,mi)%snow%th(jj) * &
																		ice(sc1,jjj,mi)%snow%t(jj) * &
																		ice(sc1,jjj,mi)%snow%d(jj)
																enddo
																if (tmp2 .gt. 0.d0) then ! protect against error in snow%z - should not be needed!
																	tmp1 = tmp1/tmp2 ! mean snow temp
                                endif
																tmp2 = tmp2/IceD ! pure ice contained in snow (m^3)
																tmp2 = min(0.5d0,tmp2/id_sw_gap) ! % ice in gaps
															endif

															! find ice/seawater mix properties in gaps
                              heat_mean = (-1.d0*Lf + tmp1*c0)*tmp2*IceD + &  ! ice heat*density
																	(1.d0-tmp2)*cw*f(mi)%t*f(mi)%d
															d_mean = tmp2*IceD + f(mi)%d*(1.d0-tmp2)
															dz_mean = 1.d0
															s_mean = f(mi)%s*(1.d0-tmp2) ! % ice in gaps

															! note - sia2_env_brine_calc.inc.f90 uses tmp1 & tmp2
#include "sia2_env_temp_from_heat.inc.f90"


#include "sia2_env_brine_calc.inc.f90"

															bv_mean = c_001*bv_mean ! brine volume,fraction

                          else

                              id_sw_gap = 0.d0

                         endif


                          ! add incoming internal ice layers
                          do ii=1,int_z_new

                              ! re-adjust new layer bottom for next snow layer temp adjustment
                              new_layer_top = new_layer_bot
                              new_layer_bot = new_layer_bot + th_new(ii) ! th_new comes from snow_grid calc
                              dz_total = 0.d0

                               ! ----------------------------------------
                               ! Add seawater layers here!
                               ! ----------------------------------------

                               if (id_sw_gap .gt. 0.d0 .and. new_layer_bot .ge. id_insert) then
                                   id_sw_gap = id_sw_gap - th_new(ii)
                                   if (id_sw_gap .lt. 0.d0) then
                                       dz = id_sw_gap + th_new(ii)
                                       id_sw_gap = 0.d0
                                   else
                                       dz = th_new(ii)
                                   endif

                                   if (dz .lt. 0.d0) then
                                       dz = sqrt(dz)
                                   endif

                                   tmp1 = dz*a_scale
                                   s_new(ii) = s_new(ii) + s_mean*tmp1
                                   d_new(ii) = d_new(ii) + d_mean*tmp1
                                   heat_new(ii) =  heat_new(ii) + heat_mean*tmp1
																	 tmp1 = tmp1*bv_mean
                                   NO3_new(ii) = NO3_new(ii) + f(mi)%no3*tmp1
                                   NH4_new(ii) = NH4_new(ii) + f(mi)%nh4*tmp1
                                   PO4_new(ii) = PO4_new(ii) + f(mi)%po4*tmp1
                                   SiOH4_new(ii) = SiOH4_new(ii) + f(mi)%sioh4*tmp1
                                   poc_new(ii) = poc_new(ii) + f(mi)%poc*tmp1
                                   smalg_new(ii) = smalg_new(ii) + alg_wc*tmp1

                                   dz_total = dz
                               endif

                               if (dz_total .lt. th_new(ii)) then

                                   ! if layer is already partialy filled with seawater,
                                   ! change "new_layer_top" by subtracting the size of
                                   ! the gap, so that it maps correctly to the old ice
                                   ! in the remapping routine below
                                   if (dz_total .gt. 0.d0) then
                                       new_layer_top = new_layer_top - hm*void_f
                                       new_layer_bot = new_layer_top + th_new(ii) - dz_total
                                   endif

                               ! --------------------------------------------

                                  ! find 1st old layer that contains part of new layer, going down
                                  jj = 0
                                  old_layer_bot = 0.
                                  do while((new_layer_top .ge. old_layer_bot) .and. (jj .lt. sk_z))
                                       jj=jj+1
                                       ! old thickness...
                                       old_layer_bot = old_layer_bot + th_old(jj)
                                   enddo
                                   old_layer_top = old_layer_bot - th_old(jj)
                                   ! now jj is OLD layer where NEW layer ii starts...

                                   ! find total heat/salt from multiple layers/partial layers that make up
                                   ! the new layer
                                   do while ((old_layer_top .lt. new_layer_bot) .and. (jj .le. sk_z))

                                       ! ----- NEW LAYER GEOMETRIES ------------------------------
                                       interim_top = max(old_layer_top,new_layer_top)
                                       interim_bot = min(old_layer_bot,new_layer_bot)

                                       z1 = interim_top - old_layer_top
                                       z2 = interim_bot - old_layer_top
                                       dz = z2-z1

                                    if (dz .lt. 0.d0) then
                                       dz = sqrt(dz)
                                   endif


                                       ! record tracers
                                       tmp1 = dz*a_scale
                                       s_new(ii) = s_new(ii) + ice(sc1,jjj,mi)%s(jj)*tmp1
                                       d_new(ii) = d_new(ii) + ice(sc1,jjj,mi)%d(jj)*tmp1
                                       heat_new(ii) = heat_new(ii) + ice(sc1,jjj,mi)%heat(jj)*tmp1
                                       no3_new(ii) = no3_new(ii) + ice(sc1,jjj,mi)%no3(jj)*tmp1
                                       nh4_new(ii) = nh4_new(ii) + ice(sc1,jjj,mi)%nh4(jj)*tmp1
                                       po4_new(ii) = po4_new(ii) + ice(sc1,jjj,mi)%po4(jj)*tmp1
                                       sioh4_new(ii) = sioh4_new(ii) + ice(sc1,jjj,mi)%sioh4(jj)*tmp1
                                       smalg_new(ii) = smalg_new(ii) + ice(sc1,jjj,mi)%smalg(jj)*tmp1
                                       poc_new(ii) = poc_new(ii) + ice(sc1,jjj,mi)%poc(jj)*tmp1
                                       ! ----- SETUP VARIABLES FOR NEXT PARTIAL LAYER ------------
                                       ! keeping track of emerging new layer thickness, for
                                       dz_total = dz_total + dz

                                       jj=jj+1
                                       if (jj .le. sk_z) then
                                           ! find boundaries of next old layer
                                           old_layer_top = old_layer_bot
                                           old_layer_bot = old_layer_bot + th_old(jj)
                                       endif

                                   enddo  ! end of layer ii mapping

                               endif ! end of "is this layer full of seawater already?"

                          enddo ! end of ii layer loop for internal ice

                          ! add incoming skeletal ice layers
                   				! note this is a lossy process - not accounting for losses
                   				! during ridging currently!!!
                         do ii=1,z_sk

                              ! find corresponding skeletal layers
                              jj = ice(sc1,jjj,mi)%z-z_sk+ii ! ic skeletal layer
                              kk = sk_1_rdg - 1 + ii

                              !tmp1 = adv(mi)%a_rdg(jjj)*ice(sc1,jjj,mi)%af/af_ic(jjj)
                              !tmp1 = adv(mi)%rf(ic,jjj)*tmp1*ice(sc1,jjj,mi)%th(jj)
                              tmp1 = a_rdg_rdg*ice(sc1,jjj,mi)%th(jj)

                              th_new(kk) = th_new(kk) + tmp1
                              s_new(kk) = s_new(kk) + ice(sc1,jjj,mi)%s(jj)*tmp1
                              d_new(kk) = d_new(kk) + ice(sc1,jjj,mi)%d(jj)*tmp1
                              heat_new(kk) = heat_new(kk) + ice(sc1,jjj,mi)%heat(jj)*tmp1
                              smalg_new(kk) = smalg_new(kk) + ice(sc1,jjj,mi)%smalg(jj)*tmp1
                              poc_new(kk) = poc_new(kk) + ice(sc1,jjj,mi)%poc(jj)*tmp1
                              no3_new(kk) = no3_new(kk) + ice(sc1,jjj,mi)%no3(jj)*tmp1
                              nh4_new(kk) = nh4_new(kk) + ice(sc1,jjj,mi)%nh4(jj)*tmp1
                              po4_new(kk) = po4_new(kk) + ice(sc1,jjj,mi)%po4(jj)*tmp1
                              sioh4_new(kk) = sioh4_new(kk) + ice(sc1,jjj,mi)%sioh4(jj)*tmp1

                          enddo

                          ! account for snow during ridging
                          tmp1 = adv(mi)%a_rdg(jjj)*ice(sc1,jjj,mi)%af/af_ic(jjj)
                          tmp1 = adv(mi)%rf(ic,jjj)*tmp1 ! ridging area - before ridging reduces it
                          if (snow_ridging .eq. 0) then
															! record fresh water flux (kg)
															call sia2_snow_mass(ice(sc,ic,mi),tmp2)
															m(mi)%h2o_flux = m(mi)%h2o_flux + &
																tmp2*(tmp1-a_rdg_rdg)*1.e6
															! take out increase in depth via ridging here (i.e. lose mass!)
                              tmp1 = a_rdg_rdg
                          endif

                          if (ice(sc1,jjj,mi)%snow%z .gt. 0) then

                              ! instead of 'ridging' snow mixed between categories as a mass only -
                              ! the mass of snow 'ridged in' on ridging ice is spread evenly across
                              ! the category.  Then the snow distribtion adjustment should compensate
                              ! for the snow depth variation, especially since it blows around a lot
                              ! anyway

                              z_snow = max(z_snow,ice(sc1,jjj,mi)%snow%z)
                              do ii=1, ice(sc1,jjj,mi)%snow%z
                                  tmp2 = tmp1*ice(sc1,jjj,mi)%snow%th(ii)
                                  th_new_snow(ii) = th_new_snow(ii) + tmp2
                                  d_new_snow(ii) = d_new_snow(ii) + ice(sc1,jjj,mi)%snow%d(ii)*tmp2
                                  heat_new_snow(ii) = heat_new_snow(ii) + ice(sc1,jjj,mi)%snow%heat(ii)*tmp2
                                  melt_new(ii) = melt_new(ii) + ice(sc1,jjj,mi)%snow%melt(ii)*tmp2
                              enddo

                          else

                              sh_prev_sum = sh_prev_sum + ice(sc1,jjj,mi)%sh_prev*tmp1

                          endif  ! end of snow input check

                          ice(sc,ic,mi)%snow%z = z_snow ! update category snow layers
                          sh_offset_sum = sh_offset_sum + ice(sc1,jjj,mi)%sh_offset*tmp1

													! add ml layer tracers
													no3_new(ml_z) = no3_new(ml_z) + &
														ice(sc1,jjj,mi)%no3(ml_z)*a_rdg_rdg
													nh4_new(ml_z) = nh4_new(ml_z) + &
														ice(sc1,jjj,mi)%nh4(ml_z)*a_rdg_rdg
													po4_new(ml_z) = po4_new(ml_z) + &
														ice(sc1,jjj,mi)%po4(ml_z)*a_rdg_rdg
													sioh4_new(ml_z) = sioh4_new(ml_z) + &
														ice(sc1,jjj,mi)%sioh4(ml_z)*a_rdg_rdg

                          ! update age and ridged tracers from ridging ice
                          age_new = age_new + a_rdg_rdg*ice(sc1,jjj,mi)%age
                          ridged_new = ridged_new + a_rdg_rdg
                          snow_dist_new = snow_dist_new + a_rdg_rdg*ice(sc1,jjj,mi)%snow_dist
                          snow_rd_new = snow_rd_new + a_rdg_rdg*ice(sc1,jjj,mi)%snow_rd

                      endif ! end of 'does ice category (jjj) ridge into this ice category(ic)?'

                      enddo ! iterate snow/flooding category
                      endif ! test to determine whether ridged ice should be added to this category

                  enddo ! end of iteration over ice category inputs

              !endif ! end of ic > 1 test

              if (adv(mi)%ice(sc,ic)%af/cell_area .ge. af_min) then

 							if (m(mi)%grid_h .eq. 298 .and. m(mi)%grid_v .eq. 288) then
                  testvar = -1
              endif


                  ! update final tracer values by dividing by area
                  int_z = int_z_new
                  sk_z = int_z + z_sk
                  sk_1 = int_z+1

									! temporary sanity check
									if (a_scale_sum/cell_area .gt. 1.00001 .or. &
									adv(mi)%ice(sc,ic)%af/cell_area .gt. 1.00001) then
										testvar = -1
										testvar=sqrt(testvar)
									endif
!									testvar = ice(sc,ic,mi)%af*cell_area ! total area before ridging
!									tmp1 = adv(mi)%a_convg
!									h_scale = testvar/(testvar - tmp1)
!									tmp2 = ice(sc,ic,mi)%id(ice(sc,ic,mi)%z-z_sk)
!									h_scale = tmp2*h_scale
!
!									tmp3 = adv(mi)%a_rdg(1)*ff(1) * & ! ridged area after ridging
!									 	tmp2/ff(1)/(1.d0-void_f)* void_f ! height injected seawater
!									tmp3 = tmp3/(testvar - tmp1) ! injected water spread across new area
!
!									if (abs((h_scale+tmp3) - hm) .gt. 1.d-10) then
!										testvar=-1
!										testvar=sqrt(testvar)
!									endif
!


                  do ii=1,int_z

                      ! update tracers
                      tmp1 = a_scale_sum*th_new(ii)
                      ice(sc,ic,mi)%no3(ii) = NO3_new(ii)/tmp1
                      ice(sc,ic,mi)%nh4(ii) = NH4_new(ii)/tmp1
                      ice(sc,ic,mi)%po4(ii) = PO4_new(ii)/tmp1
                      ice(sc,ic,mi)%sioh4(ii) = SiOH4_new(ii)/tmp1
                      ice(sc,ic,mi)%poc(ii) = poc_new(ii)/tmp1
                      ice(sc,ic,mi)%smalg(ii) = smalg_new(ii)/tmp1

                      ! ice physics
                      if (ii .le. int_z) then

                          ! ----- SOLVE FOR T_new using quadratic eqn -------------------
                          dz_mean = 1.
                          d_new(ii) = d_new(ii)/tmp1
                          d_mean = d_new(ii)
                          s_new(ii) = s_new(ii)/tmp1
                          S_mean = s_new(ii)
                          if (s_mean .lt. 0.09999d0) then
                          	testvar = -1
                          	testvar=sqrt(testvar)
                          endif
                          heat_mean = heat_new(ii)/tmp1 ! division by dz is included in calc below

#include "sia2_env_temp_from_heat.inc.f90"

                          ! check to make sure we don't get warmer than the ocean
                          !if (ice(sc,ic,mi)%id(ii) .le. ice(sc,ic,mi)%fbh) then
                          !	  t_mean = min(t_mean,S_mean*mu)
                          !else
                          !	  t_mean = min(t_mean,f(mi)%t)
                          !endif

                          t_new(ii) = t_mean
                      else
                          ! skeletal paramaters are constant
                          s_new(ii) = f(mi)%s*c_5
                          t_new(ii) = min(f(mi)%t,f(mi)%s*mu)
                      endif


                      ice(sc,ic,mi)%t(ii) = t_new(ii)
                      ice(sc,ic,mi)%s(ii) = s_new(ii)
                      T_mean = ice(sc,ic,mi)%t(ii)
                      S_mean = max(min_sal,ice(sc,ic,mi)%s(ii))

#include "sia2_env_brine_calc.inc.f90"

                      ice(sc,ic,mi)%d(ii) = d_mean
                      ice(sc,ic,mi)%bs(ii) = bs_mean
                      ice(sc,ic,mi)%bd(ii) = bd_mean
                      ice(sc,ic,mi)%bv(ii) = bv_mean
                      ice(sc,ic,mi)%heat(ii) = heat_mean

                      ! update th and id (in case ridging has intialized a new category)
                      ice(sc,ic,mi)%th(ii) = th_new(ii)
                      if (ii .eq. 1) then
                          ice(sc,ic,mi)%id(ii) = th_new(ii)
                      else


                          ice(sc,ic,mi)%id(ii) = ice(sc,ic,mi)%id(ii-1) + th_new(ii)

!  								if(ii .eq. 40) then
!  								if (abs((h_scale+tmp3) - ice(sc,ic,mi)%id(ii)) .gt. 1.d-10) then
!          		print *,'mi i j ',mi,m(mi)%grid_h,m(mi)%grid_v
!							print *,'hm ',hm
!							print *,'h_scale ',h_scale
!							print *,'ff ',ff(1)
!							print *,'id_new(int_z_new) ',id_new(int_z_new)
!							print *,'int_z_new ',int_z_new
!							print *,'ice%z ',ice(1,1,mi)%z
!							print *,'id_new - ',id_new
!							print *,'ice_id ',ice(1,1,mi)%id(ice(1,1,mi)%z-z_sk)
!							print *,'th_new - ',th_new
!							print *,'id_old ',ice_ridged(1,1,thread_n)%id(ice_ridged(1,1,thread_n)%z-z_sk)
!							print *,'adv--af ',adv(mi)%ice(1,1)%af
!							print *,'ice--af ',ice(1,1,mi)%af
!							print *,'a_rdg ',adv(mi)%a_rdg(1)
!							print *,'a_convg ',adv(mi)%a_convg
!							print *,'ts ',ice(1,1,mi)%ts
!							print *,'t1 ',ice(1,1,mi)%t(1)
!							print *,'s1 ',ice(1,1,mi)%s(1)
!							print *,'snowd ',ice(1,1,mi)%snow%depth
!							print *,'snow z ',ice(1,1,mi)%snow%z
!							print *,'old snowd ',ice_ridged(1,1,thread_n)%snow%depth
!							print *,'old snow z ',ice_ridged(1,1,thread_n)%snow%z
!										pause
!										testvar=-1
!										testvar=sqrt(testvar)
!									endif
!									endif

                      endif

                  enddo

                  ! update skeletal layers
                  do jj=sk_1,sk_z

                      ice(sc,ic,mi)%th(jj) = th_new(jj)/adv(mi)%ice(sc,ic)%af
                      ice(sc,ic,mi)%smalg(jj) = smalg_new(jj)/th_new(jj)
                      ice(sc,ic,mi)%poc(jj) = poc_new(jj)/th_new(jj)
                      ice(sc,ic,mi)%no3(jj) = no3_new(jj)/th_new(jj)
                      ice(sc,ic,mi)%nh4(jj) = nh4_new(jj)/th_new(jj)
                      ice(sc,ic,mi)%po4(jj) = po4_new(jj)/th_new(jj)
                      ice(sc,ic,mi)%sioh4(jj) = sioh4_new(jj)/th_new(jj)

                  enddo
                  dz_sk = 0.

#include "sia2_env_update_skeletal.inc.f90"

                  ! update TS if it doesn't exist yet/category is new
                  ! --------------------------------------------------------
                  if (ice(sc,ic,mi)%z .eq. 0) then
                      if (z_snow .gt. 0) then
                          ice(sc,ic,mi)%snow%ts = ice(sc,ic,mi)%snow%t(z_snow)
                      else
                          ice(sc,ic,mi)%snow%ts = ice(sc,ic,mi)%t(1)
                      endif
                  endif
                  ice(sc,ic,mi)%z  = sk_z

                  ! Update snow in category
                  ! --------------------------------------------------------

                  ice(sc,ic,mi)%sh_offset = sh_offset_sum/adv(mi)%ice(sc,ic)%af
                  r_depth = sh_prev_sum/adv(mi)%ice(sc,ic)%af

                  if (z_snow .gt. 0) then
                      do ii=1,z_snow

                          if (th_new_snow(ii) .eq. 0. .or. d_new_snow(ii) .eq. 0.) then
                              testvar = -1
                          endif

                          ice(sc,ic,mi)%snow%th(ii) = th_new_snow(ii)/adv(mi)%ice(sc,ic)%af
                          r_depth = r_depth + ice(sc,ic,mi)%snow%th(ii)
                          ice(sc,ic,mi)%snow%d(ii) = d_new_snow(ii)/th_new_snow(ii)
                          heat_total = heat_new_snow(ii)/th_new_snow(ii)/ice(sc,ic,mi)%snow%d(ii)
                          ice(sc,ic,mi)%snow%t(ii) = (0.2309-sqrt(0.05331481+0.0136* &
                             (heat_snow0-heat_total)))/(-0.0068)
                          ice(sc,ic,mi)%snow%t(ii) = min(ice(sc,ic,mi)%snow%t(ii) - kelvin0,d0_)
                          ice(sc,ic,mi)%snow%melt(ii) = melt_new(ii)/th_new_snow(ii)
                          t_mean = ice(sc,ic,mi)%snow%t(ii) + kelvin0
                          d_mean = ice(sc,ic,mi)%snow%d(ii)
#include "sia2_env_snow_calc.inc.f90"
                          ice(sc,ic,mi)%snow%heat(ii) = heat_mean

                      enddo
                  endif

                  ice(sc,ic,mi)%sh_prev = r_depth

#include "sia2_env_snow_grid.inc.f90"

                  flooded = 0.d0
                  melt_drain = 0.d0

#include "sia2_env_snow_remap.inc.f90"


								! Move PUR, pond, other tracers??
								! ------------------------------------------------------------------

								ice(sc,ic,mi)%no3(ml_z) = no3_new(ml_z)/adv(mi)%ice(sc,ic)%af
								ice(sc,ic,mi)%nh4(ml_z) = nh4_new(ml_z)/adv(mi)%ice(sc,ic)%af
								ice(sc,ic,mi)%po4(ml_z) = po4_new(ml_z)/adv(mi)%ice(sc,ic)%af
								ice(sc,ic,mi)%sioh4(ml_z) = sioh4_new(ml_z)/adv(mi)%ice(sc,ic)%af

								ice(sc,ic,mi)%age = age_new/adv(mi)%ice(sc,ic)%af
								ice(sc,ic,mi)%ridged = ridged_new/adv(mi)%ice(sc,ic)%af
								ice(sc,ic,mi)%snow_dist = snow_dist_new/adv(mi)%ice(sc,ic)%af
								ice(sc,ic,mi)%snow_rd = snow_rd_new/adv(mi)%ice(sc,ic)%af

              else

                  ! null category if it drops below minimum area
                  adv(mi)%ice(sc,ic)%af = 0.

              endif ! end of test to see if category has any ice to update

              ! copy and move ice category to new ridged grid cell
              ! ---------------------------------------------------------
!							if (m(mi)%grid_h .eq. 298 .and. m(mi)%grid_v .eq. 288) then
!          		print *,'mi i j ',mi,m(mi)%grid_h,m(mi)%grid_v
!								print *,'hm ',hm
!							print *,'h_scale ',h_scale
!							print *,'ff ',ff(1)
!							print *,'id_new ',id_new(40)
!							print *,'ice_id ',ice(1,1,mi)%id(ice(1,1,mi)%z-z_sk)
!							print *,'id_old ',ice_ridged(1,1,thread_n)%id(ice_ridged(1,1,thread_n)%z-z_sk)
!							print *,'adv--af ',adv(mi)%ice(1,1)%af
!							print *,'ice--af ',ice(1,1,mi)%af
!							print *,'a_rdg ',adv(mi)%a_rdg(1)
!							print *,'a_convg ',adv(mi)%a_convg
!							print *,'ts ',ice(1,1,mi)%ts
!							print *,'t1 ',ice(1,1,mi)%t(1)
!							print *,'s1 ',ice(1,1,mi)%s(1)
!							print *,'snowd ',ice(1,1,mi)%snow%depth
!							print *,'snow z ',ice(1,1,mi)%snow%z
!							print *,'old snowd ',ice_ridged(1,1,thread_n)%snow%depth
!							print *,'old snow z ',ice_ridged(1,1,thread_n)%snow%z
!							endif



              ice_temp(thread_n) = ice(sc,ic,mi) ! ridged ice copy
              ice(sc,ic,mi) = ice_ridged(sc,ic,thread_n) ! copy current (non-ridged) ice back to memory
            	ice_ridged(sc,ic,thread_n) = ice_temp(thread_n) ! update ridged ice grid cell w/ ridged ice

          enddo ! end if sda_n tracer ridging loop
          enddo ! end if ida_n tracer ridging loop

      !endif ! end of 'are we using ice categories' test


			! move ridged grid cell into place of current grid cell
			! ----------------------------------------------------------------
			do ic=1,ida_n
			do sc=1,sda_n
				ice(sc,ic,mi) = ice_ridged(sc,ic,thread_n)
			enddo
			enddo

      ! special case of largest ice category - or when running without ice categories (ida_n == 1)
			! largest ice category may grow to thick - must reduce thickness
			! and account for non-conservative mass/bio losses
			! ------------------------------------------------------------
      ic = ida_n
      do sc=1,sda_n
      if ((adv(mi)%ice(sc,ic)%af .gt. 0.d0)  .and. (ice(sc,ic,mi)%z > 0)) then
        if(ice(sc,ic,mi)%id(ice(sc,ic,mi)%z-z_sk) .gt. h_max) then

          ! update these values for remapping below
          sk_z = ice(sc,ic,mi)%z
          int_z = sk_z - z_sk
          sk_1 = int_z + 1

          h_scale = h_max/ice(sc,ic,mi)%id(int_z)

					do ii=1,int_z
						dz = ice(sc,ic,mi)%th(ii)*h_scale

            ! report biomass disappeared to model domain
            smalg_lost = ice(sc,ic,mi)%smalg(ii) - min_alg
					  ! converting to g/pixel - (cell_area*1.e6)*1g/1000mg = 1000. scale factor
					  m(mi)%bm_lost = m(mi)%bm_lost + smalg_lost* &
					  dz*1.e3*adv(mi)%ice(sc,ic)%af   ! g/pixel

					  ! put salt and h20 back in ocean, I guess?
						call sia2_ice_mass(ice(sc,ic,mi),tmp1,tmp2)
						! record fresh water flux (kg)
						m(mi)%h2o_flux = m(mi)%h2o_flux + &
							tmp1*adv(mi)%ice(sc,ic)%af*1.e6*h_scale
						! record salt flux (kg)
						m(mi)%salt_flux = m(mi)%salt_flux + &
							tmp2*adv(mi)%ice(sc,ic)%af*1.e6*h_scale

					enddo

          ice(sc,ic,mi)%th = ice(sc,ic,mi)%th*h_scale
          ice(sc,ic,mi)%id = ice(sc,ic,mi)%id*h_scale

          ! find th_new, id_new and int_z_new in terms of proper height allocation for grid
          r_depth = h_max

#include "sia2_env_ice_grid.inc.f90"

          ! start new layer top at old layer top - don't exclude any ice
          new_layer_top = -1.
          ! no growth/melt here
          s_gl = 0.
          c_gl = 0.
          sk_gl = 0.d0
          ! no flooding here
          flooded = 0.
          melt_flood = 0.
          dhdt_cm_s = 0.   ! required in algal migration


#include "sia2_env_ice_remap.inc.f90"

        endif
      endif  ! end of test to see if the largest category is involved w/ ridging
      enddo ! end interation over snow categories


!
!      ! special case of largest ice category - or when running without ice categories (ida_n == 1)
!      ! perform height scaling to account for ridging of this category
!      ! ---------------------------------------------------------------------------
!			ice_temp(thread_n) = ice(1,1,mi)
!			ice(1,1,mi) = ice_temp2
!
!      ic = ida_n
!      do sc=1,sda_n
!      if (adv(mi)%ice(sc,ic)%af .gt. 0.) then
!
!          ! areal reduction in largest ice category - we only ridge larger, so all ridging in this category stays...
!          if (ic .gt. 1) then
!              tmp2 = 0.d0
!              do jj=1,sda_n
!                  tmp2 = tmp2 + adv(mi)%ice(jj,ic)%af
!              enddo
!              tmp2 = adv(mi)%ice(sc,ic)%af/tmp2  ! snow category percentage of total ice category
!              tmp1 = tmp2*adv(mi)%a_rdg(ic)*(1. - ff(ic)) ! area to remove
!          else
!              tmp2 = ice(sc,ic,mi)%af/af_ic(ic)  ! fractional of ice catagory area that is this snow category
!              adv(mi)%ice(sc,ic)%af = ice(sc,ic,mi)%af*cell_area ! this is normally taken care of in ridging, but not ridging w/o ice categories...
!              tmp1 = adv(mi)%a_convg*tmp2
!          endif
!          ! find height scaling of category
!          h_scale = adv(mi)%ice(sc,ic)%af/(adv(mi)%ice(sc,ic)%af - tmp1) ! find amount to scale height by
!          adv(mi)%ice(sc,ic)%af = adv(mi)%ice(sc,ic)%af - tmp1 ! shrink ridging ice area
!          ice(sc,ic,mi)%th = ice(sc,ic,mi)%th*h_scale
!          ice(sc,ic,mi)%id = ice(sc,ic,mi)%id*h_scale
!
!          ! find th_new, id_new and int_z_new in terms of proper height allocation for grid
!          r_depth = ice(sc,ic,mi)%id(ice(sc,ic,mi)%z-z_sk)
!
!          ! update these values for remapping below
!          sk_z = ice(sc,ic,mi)%z
!          int_z = sk_z - z_sk
!          sk_1 = int_z + 1
!
!#include "sia2_env_ice_grid.inc.f90"
!
!          ! start new layer top at old layer top - don't exclude any ice
!          new_layer_top = -1.
!          ! no growth/melt here
!          s_gl = 0.
!          c_gl = 0.
!          ! no flooding here
!          flooded = 0.
!          melt_flood = 0.
!          dhdt_cm_s = 0.   ! required in algal migration
!
!          ! record bm_lost for exceeding maxmimum depth
!          if (maxed_depth) then
!              dz_mean = r_depth
!              ii = sk_z
!              do while(dz_mean .lt. 0.)
!                  smalg_lost = ice(sc,ic,mi)%smalg(ii) - min_alg
!                  smalg_lost = max(0.,smalg_lost)
!                  if (dz_mean .gt. ice(sc,ic,mi)%th(ii)) then
!                     ! converting to g/pixel
!                     m(mi)%bm_lost = m(mi)%bm_lost + smalg_lost* &
!                     ice(sc,ic,mi)%th(ii)*cell_area*1.e3*ice(sc,ic,mi)%af   ! g/pixel
!                     dz_mean = dz_mean + ice(sc,ic,mi)%th(ii)
!                     ii = ii-1
!                  else
!                     ! converting to g/pixel - (cell_area*1.e6)*1g/1000mg = 1000. scale factor
!                     m(mi)%bm_lost = m(mi)%bm_lost + smalg_lost* &
!                     dz_mean*cell_area*1.e3*ice(sc,ic,mi)%af   ! g/pixel
!                     dz_mean = 0.
!                  endif
!              enddo
!             ! short-circuit built-in protection from exceeding max depth - always need
!             ! to adjust boundaries during advection
!             maxed_depth = .false.
!         endif
!
!#include "sia2_env_ice_remap.inc.f90"
!
!          ! scale snow in largest category to conserve mass
!          if (snow_ridging .eq. 1) then
!              if (ice(sc,ic,mi)%snow%z .gt. 0) then
!
!                  ice(sc,ic,mi)%snow%depth = ice(sc,ic,mi)%snow%depth*h_scale
!                  ice(sc,ic,mi)%snow%th = ice(sc,ic,mi)%snow%th*h_scale
!
!                  r_depth = ice(sc,ic,mi)%snow%depth
!                  ice(sc,ic,mi)%sh_prev = r_depth
!
!#include "sia2_env_snow_grid.inc.f90"
!
!                  flooded = 0.d0
!                  melt_drain = 0.d0
!
!#include "sia2_env_snow_remap.inc.f90"
!
!              endif
!           endif
!
!      endif  ! end of test to see if the largest category is involved w/ ridging
!      enddo ! end interation over snow categories
!
!			if(abs(ice(1,1,mi)%id(ice(1,1,mi)%z-z_sk) - &
!				ice_temp(thread_n)%id(ice_temp(thread_n)%z-z_sk)) .gt. &
!				1.d-3) then
!				testvar=-1
!				testvar=sqrt(testvar)
!			endif
!
!			ice(1,1,mi) = ice_temp(thread_n)

      ! update areal fractions - remember that right now the adv%af values
      ! are the total category area, so...
      tmp1 = 0.
      do ic=1,ida_n
      do sc=1,sda_n
          if(adv(mi)%ice(sc,ic)%af/cell_area .ge. af_min) then
              tmp1 = tmp1 + adv(mi)%ice(sc,ic)%af/cell_area
          endif
      enddo
      enddo
      if (tmp1 .gt. 1.01) then
          print *,'Holy Crap! an ice fraction greater than one!', mi
      endif

      tmp1 = 0.
      do ic=1,ida_n
      do sc=1,sda_n
          if(adv(mi)%ice(sc,ic)%af/cell_area .lt. af_min) then

              ice(sc,ic,mi)%af = 0.

              ! zero layer
#include "sia2_env_null_params.f90"

          else
              ice(sc,ic,mi)%af = adv(mi)%ice(sc,ic)%af/cell_area
              tmp1 = tmp1 + ice(sc,ic,mi)%af
          endif
      enddo
      enddo

!							if (m(mi)%grid_h .eq. 298 .and. m(mi)%grid_v .eq. 288) then
!							print *,'------------------------------------------------'
!							print *,'hm ',hm
!							print *,'h_scale ',h_scale
!							print *,'ff ',ff(1)
!							print *,'id_new ',id_new(40)
!							print *,'ice_id ',ice(1,1,mi)%id(ice(1,1,mi)%z-z_sk)
!							print *,'id_old ',ice_ridged(1,1,thread_n)%id(ice_ridged(1,1,thread_n)%z-z_sk)
!							print *,'adv--af ',adv(mi)%ice(1,1)%af
!							print *,'ice--af ',ice(1,1,mi)%af
!							print *,'a_rdg ',adv(mi)%a_rdg(1)
!							print *,'a_convg ',adv(mi)%a_convg
!							print *,'ts ',ice(1,1,mi)%ts
!							print *,'t1 ',ice(1,1,mi)%t(1)
!								print *,'s1 ',ice(1,1,mi)%s(1)
!						print *,'snowd ',ice(1,1,mi)%snow%depth
!							print *,'snow z ',ice(1,1,mi)%snow%z
!							print *,'old snowd ',ice_ridged(1,1,thread_n)%snow%depth
!							print *,'old snow z ',ice_ridged(1,1,thread_n)%snow%z
!							endif



      ! correct machine precision error on areal fractions
!      tmp2 = f(mi)%ic_interp_next - tmp1
!      if (abs(tmp2) .lt. 1.e-8) then
!          ice(1,mi)%af = ice(1,mi)%af + tmp2
!      else
!          testvar = -1
!      endif


