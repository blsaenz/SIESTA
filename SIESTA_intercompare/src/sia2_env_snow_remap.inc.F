

		   ! Adjust boundaries of snow layers and calculate new layer temperatures
		   ! using snow heat capacity to conserve total heat.
		   ! ------------------------------------------------------------------

	       z_snow = ice(sc,ic,mi)%snow%z
	       keep_growing = .TRUE.  ! set this variable to trigger rest of layer-tracking vars

           if (z_snow .gt. 0) then

			   ! setup vars to track new and old layer heights
			   old_layer_top = d0_
			   new_layer_bot = d0_
			   old_layer_bot = d0_

			   ! allocate new layer tracers
			   do ii=1,z_snow

				   ! These vars keep track of new ice physics while layers are mixing
				   dz_total = d0_
				   d_new(ii) = d0_
				   heat_total = d0_
				   melt_new(ii) = d0_

				   if (melt_drain .gt. d0_) then

                       ! re-adjust new layer bottom for next snow layer temp adjustment
                       new_layer_top = new_layer_bot
                       new_layer_bot = new_layer_bot + th_new(ii) ! th_new comes from snow_grid calc

                       ! find 1st old layer that contains part of new layer, going down
                       jj = z_snow_new+1
                       old_layer_bot = d0_
                       do while((new_layer_top .ge. old_layer_bot) .and. (jj .gt. 1))
                           jj=jj-1
                           ! old thickness...
                           old_layer_bot = old_layer_bot + th_snow_new(jj)
                       enddo

                       ! add heat to new layer from old layer
                       if ((jj .le. z_snow_new) .and. (old_layer_bot .gt. new_layer_top)) then
                           old_layer_top = old_layer_bot - th_snow_new(jj)
                           ! jj is OLD layer where NEW layer ii starts...

                           ! add heat to new layer from old layer
                           do while ((old_layer_top .lt. new_layer_bot) .and. (jj .gt. 0))
                               ! make sure that drained layer still has some ice left in it
                               if (d_snow_new(jj) .gt. 0.d0) then
                                   ! ----- NEW LAYER GEOMETRIES ------------------------------

                                   interim_top = max(old_layer_top,new_layer_top)
                                   interim_bot = min(old_layer_bot,new_layer_bot)

                                   z1 = interim_top - old_layer_top
                                   z2 = interim_bot - old_layer_top
                                   dz = z2-z1

                                   ! record density*thickness, heat*thickness
                                   d_new(ii) = d_new(ii) + d_snow_new(jj)*dz
                                   heat_total = heat_total + heat_snow_new(jj)*dz
                                   melt_new(ii) = melt_new(ii) + dz

                                   ! ----- SETUP VARIABLES FOR NEXT PARTIAL LAYER ------------
                                   ! keeping track of emerging new layer thickness, for
                                   dz_total = dz_total + dz

                               endif
                               ! find boundaries of next old layer
                               jj=jj-1
                               if (jj .gt. 0) then
                                   old_layer_top = old_layer_bot
                                   old_layer_bot = old_layer_bot + th_snow_new(jj)
                               endif
                           enddo

                       endif

                   endif


				   if (dz_total .lt. th_new(ii)) then

                       if (keep_growing) then
                           ! setup vars to track new and old layer heights
                           old_layer_top = d0_
                           if (flooded .gt. d0_) then
                               new_layer_bot = flooded ! start new layers higher up if some snow flooded
                           else
                               new_layer_bot = d0_
                           endif
                           old_layer_bot = d0_
                           keep_growing = .FALSE.
                           melt_drain = d0_
                       endif

					   ! add newly-fallen snow to layer
					   if (z_last .eq. 0) then

						   ! fill layer with new snow
						   dz = th_new(ii) - dz_total

						   ! new snow density and 'melt' status
						   if (airtemp1c .lt. des_s_switch) then
							   d_mean = den_s_dry*1.e6
						   else
							   d_mean = den_s_wet*1.e6
							   melt_new(ii) = melt_new(ii) + dz
						   endif

						   ! find median depth of new snow sub-layer, store in tmp1
						   tmp1 = d0_
						   do jj=1,ii
							   tmp1 = tmp1 + th_new(jj)
						   enddo
						   tmp1 = tmp1 - dz_total - th_new(ii)*c_5

						   ! find temp of snow, interpolated between surface temp and airtemp
						   t_mean = ice(sc,ic,mi)%snow%ts + tmp1*(min(airtemp1c,d0_) - &
							   ice(sc,ic,mi)%snow%ts)/ice(sc,ic,mi)%snow%depth

                           t_mean = t_mean + kelvin0

						   ! find heat mean
						   heat_mean = d_mean*(heat_snow0 - 0.2309*T_mean - &
							   0.0034*T_mean**2)

						   ! record tracers (melt already done above)
						   d_new(ii) = d_new(ii) + d_mean*dz
						   heat_total = heat_total + heat_mean*dz
						   dz_total = th_new(ii)

					   ! map old snow to new snow layer
					   else

						   ! re-adjust new layer bottom for next snow layer temp adjustment
						   new_layer_top = new_layer_bot
						   new_layer_bot = new_layer_bot + th_new(ii) ! th_new comes from snow_grid calc

						   ! find 1st old layer that contains part of new layer, going down
						   jj = 0
						   old_layer_bot = d0_
						   do while((new_layer_top .ge. old_layer_bot) .and. (jj .lt. z_last))
							   jj=jj+1
							   ! old thickness...
							   old_layer_bot = old_layer_bot + ice(sc,ic,mi)%snow%th(jj)
						   enddo

						   ! add heat to new layer from old layer
						   if ((jj .le. z_last) .and. (old_layer_bot .gt. new_layer_top)) then
							   old_layer_top = old_layer_bot - ice(sc,ic,mi)%snow%th(jj)
							   ! jj is OLD layer where NEW layer ii starts...

							   ! add heat to new layer from old layer
							   do while ((old_layer_top .lt. new_layer_bot) .and. (jj .le. z_last))
								   ! ----- NEW LAYER GEOMETRIES ------------------------------

								   interim_top = max(old_layer_top,new_layer_top)
								   interim_bot = min(old_layer_bot,new_layer_bot)

								   z1 = interim_top - old_layer_top
								   z2 = interim_bot - old_layer_top
								   dz = z2-z1

								   ! record density*thickness, heat*thickness
								   d_new(ii) = d_new(ii) + ice(sc,ic,mi)%snow%d(jj)*dz
								   heat_total = heat_total + ice(sc,ic,mi)%snow%heat(jj)*dz
								   melt_new(ii) = melt_new(ii) + ice(sc,ic,mi)%snow%melt(jj)*dz

								   ! ----- SETUP VARIABLES FOR NEXT PARTIAL LAYER ------------
								   ! keeping track of emerging new layer thickness, for
								   dz_total = dz_total + dz

								   ! find boundaries of next old layer
								   jj=jj+1
								   old_layer_top = old_layer_bot
								   old_layer_bot = old_layer_bot + ice(sc,ic,mi)%snow%th(jj)
							   enddo

							   ! add heat from 'new snow' to layer (using surface temp)
							   if (dz_total .lt. th_new(ii)) then
								   ! record heat*thickness
								   T_mean = ice(sc,ic,mi)%snow%ts+kelvin0
								   d_new(ii) = d_new(ii) + snowfall_d*(th_new(ii) - dz_total)
								   heat_mean = snowfall_d*(heat_snow0 - 0.2309*T_mean - &
									 0.0034*T_mean**2)
								   heat_total = heat_total + heat_mean*(th_new(ii) - dz_total)
								   if (airtemp1c .gt. -1.0) then
									   melt_new(ii) = melt_new(ii) + (th_new(ii) - dz_total)
								   endif
								   dz_total = th_new(ii)

							   endif
						   else
							   ! add heat from 'new snow' to layer (from surface temp)
							   ! record heat*thickness
							   T_mean = ice(sc,ic,mi)%snow%ts+kelvin0
							   d_new(ii) = snowfall_d*th_new(ii)
							   heat_mean = snowfall_d*(heat_snow0 - 0.2309*T_mean - &
								 0.0034*T_mean**2)
							   heat_total = heat_total + heat_mean*th_new(ii)
							   if (airtemp1c .gt. -1.0) then
								   melt_new(ii) = melt_new(ii) + th_new(ii)
							   endif
							   dz_total = th_new(ii)
						   endif

					   endif  ! end of z_last .eq. 0 test

				   endif ! end of dz_total .lt. th_new test

				   ! Solve quadradic for new snow temp
				   d_new(ii) = d_new(ii)/dz_total
				   heat_total = heat_total/dz_total/d_new(ii)
				   t_new(ii) = (0.2309-sqrt(0.05331481+0.0136* &
					   (heat_snow0-heat_total)))/(-0.0068)
				   t_new(ii) = min(t_new(ii) - kelvin0,d0_)
				   melt_new(ii) = melt_new(ii)/dz_total
			   enddo

			   ! update snow thickness & temp
			   do ii=1,z_snow
	               ice(sc,ic,mi)%snow%th(ii) = th_new(ii)
	               ice(sc,ic,mi)%snow%t(ii) = t_new(ii)
	               T_mean = t_new(ii) + kelvin0
	               d_mean = d_new(ii)

#include "sia2_env_snow_calc.inc.f90"

				  ice(sc,ic,mi)%snow%d(ii) = d_mean
				  ice(sc,ic,mi)%snow%heat(ii) = heat_mean
				  ice(sc,ic,mi)%snow%melt(ii) = melt_new(ii)
			   enddo
			   ! zero out unused array members
               if (z_snow .lt. z_snow_max) then
				   do ii=z_snow+1,z_snow_max
					   ice(sc,ic,mi)%snow%th(ii) = d0_
					   ice(sc,ic,mi)%snow%t(ii) = d0_
					   ice(sc,ic,mi)%snow%d(ii) = d0_
					   ice(sc,ic,mi)%snow%heat(ii) = d0_
					   ice(sc,ic,mi)%snow%melt(ii) = d0_
				   enddo
               endif

		   else
               ! no snow anymore - need to do anything? right now, not zeroing old snow values
               if (z_last .gt. 0) then
                   ice(sc,ic,mi)%snow%ts = ice(sc,ic,mi)%t(1) ! move surface temp to prevent numerical heat issues - needed?
               endif
           endif


