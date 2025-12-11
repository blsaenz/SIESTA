      ! allocate temp3d array
      allocate(temp3d(z_max,tcells))

      ! advect into temporary holding variable (temp3d)
      do ic=1,ida_n
      do sc=1,sda_n
      do mi=1,tcells
      if (m(mi)%status .ge. 1) then

          ! test to see whether this category is involved in advection at all
          if (adv(mi)%ice(sc,ic)%af .gt. af_min) then

              ! initialize temporary holding array
              temp3d(:,mi) = d0_  ! vector scalar assignment

              ! find full EASE domain reference for current cell
              i_ease = m(mi)%grid_h
              j_ease = m(mi)%grid_v
              int_z_new = adv(mi)%ice(sc,ic)%z-z_sk

              do jjj=1,9

                  if (jjj .eq. 9) then
                      mi1 = mi
                      a_scale = (1.-adv(mi)%out)*ice(sc,ic,mi)%af
                  elseif (adv(mi)%in(jjj) .gt. d0_) then
                      mi1 = m(mi)%mia(jjj)
                      a_scale = adv(mi)%in(jjj)*ice(sc,ic,mi1)%af
                  else
                      a_scale = 0.
                  endif

                  if (a_scale .gt. 0.) then

					  h_scale = ice(sc,ic,mi1)%id(ice(sc,ic,mi1)%z-z_sk)/adv(mi)%ice(sc,ic)%depth
					  v_scale = a_scale*h_scale
                      th_old = ice(sc,ic,mi1)%th/h_scale ! vector multiply = streching layer, so account for ridging below
                      sk_z = ice(sc,ic,mi1)%z
                      int_z = sk_z - z_sk
                      sk_1 = int_z+1
                      tracer = &
#if tracer_var==1
                          ice(sc,ic,mi1)%s
#elif tracer_var==2
                          ice(sc,ic,mi1)%d
#elif tracer_var==3
                          ice(sc,ic,mi1)%heat
#elif tracer_var==4
                          ice(sc,ic,mi1)%smalg
#elif tracer_var==5
                          ice(sc,ic,mi1)%NO3(1:z_max)
#elif tracer_var==6
                          ice(sc,ic,mi1)%NH4(1:z_max)
#elif tracer_var==7
                          ice(sc,ic,mi1)%PO4(1:z_max)
#elif tracer_var==8
                          ice(sc,ic,mi1)%SiOH4(1:z_max)
#elif tracer_var==9
                          ice(sc,ic,mi1)%poc(1:z_max)
#endif

                      ! remap tracers to internal ice layers
                      ! ------------------------------------------------------
                      new_layer_bot = 0.

                      do ii=1,int_z_new

                          ! re-adjust new layer bottom for next snow layer temp adjustment
                          new_layer_top = new_layer_bot
                          new_layer_bot = new_layer_bot + adv(mi)%ice(sc,ic)%th_new(ii)

                          ! find 1st old layer that contains part of new layer, going down
                          jj = 0
                          old_layer_bot = 0.
                          do while((new_layer_top .ge. old_layer_bot) .and. (jj .lt. int_z))
                               jj=jj+1
                               ! old thickness...
                               old_layer_bot = old_layer_bot + th_old(jj)
                           enddo
                           old_layer_top = old_layer_bot - th_old(jj)
                           ! now jj is OLD layer where NEW layer ii starts...

                           do while ((old_layer_top .lt. new_layer_bot) .and. (jj .le. int_z))

                               ! ----- NEW LAYER GEOMETRIES ------------------------------
                               interim_top = max(old_layer_top,new_layer_top)
                               interim_bot = min(old_layer_bot,new_layer_bot)

                               z1 = interim_top - old_layer_top
                               z2 = interim_bot - old_layer_top
                               dz = z2-z1

                               ! record tracer
                               temp3d(ii,mi) = temp3d(ii,mi) + tracer(jj)*dz*v_scale

                               ! ----- SETUP VARIABLES FOR NEXT PARTIAL LAYER ------------
                               ! keeping track of emerging new layer thickness, for
                               dz_total = dz_total + dz

                               jj=jj+1
                               if (jj .le. int_z) then
                                   ! find boundaries of next old layer
                                   old_layer_top = old_layer_bot
                                   old_layer_bot = old_layer_bot + th_old(jj)
                               endif

                           enddo  ! end of layer ii mapping

                      enddo ! end of ii layer loop for internal ice

                      ! add incoming skeletal ice layers
                      do ii=1,z_sk
						  iii = ice(sc,ic,mi1)%z-z_sk+ii
						  jj = adv(mi)%ice(sc,ic)%z-z_sk+ii
                          if (jj .eq. -1) then
                              testvar = 1
                          endif
						  temp3d(jj,mi) = temp3d(jj,mi) + tracer(iii)*a_scale*ice(sc,ic,mi1)%th(iii)
                      enddo

                  endif ! end of v_scale test (is there and volume advecting in from category?)

              enddo ! end of 8 neighbor cells + local cell loop

              ! find final tracer concentrations by dividng by new advective
              ! area and thickness, and assing back into local cell
			  do ii=1,adv(mi)%ice(sc,ic)%z
				  temp3d(ii,mi) = temp3d(ii,mi)/adv(mi)%ice(sc,ic)%th_new(ii)/adv(mi)%ice(sc,ic)%af
			  enddo

          endif ! end advection-active category test


	  endif ! end status test
	  enddo ! end mi loop

	  ! copy back into original variable
	  do mi=1,tcells
	      if (m(mi)%status .ge. 1) then
	          if (adv(mi)%ice(sc,ic)%af .gt. 0.) then
#if tracer_var==1
				  ice(sc,ic,mi)%s &
#elif tracer_var==2
				  ice(sc,ic,mi)%d &
#elif tracer_var==3
				  ice(sc,ic,mi)%heat &
#elif tracer_var==4
				  ice(sc,ic,mi)%smalg &
#elif tracer_var==5
				  ice(sc,ic,mi)%NO3(1:z_max)  &
#elif tracer_var==6
				  ice(sc,ic,mi)%NH4(1:z_max)  &
#elif tracer_var==7
				  ice(sc,ic,mi)%PO4(1:z_max)  &
#elif tracer_var==8
				  ice(sc,ic,mi)%SiOH4(1:z_max)  &
#elif tracer_var==9
				  ice(sc,ic,mi)%poc(1:z_max)  &
#endif
            		  = temp3d(:,mi)

                  endif ! end advection-active category test
              endif ! end status test
          enddo ! end mi loop
          enddo ! end snow category loop
          enddo ! end ice category loop

          ! deallocate temp3d array
          deallocate(temp3d)
