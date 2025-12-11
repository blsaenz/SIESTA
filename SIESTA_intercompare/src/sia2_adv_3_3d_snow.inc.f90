      ! allocate temp3d array
      allocate(temp3d(z_max,tcells))

      ! advect into temporary holding variable (temp3d)
      do ic=1,ida_n
      do sc=1,sda_n
      do mi=1,tcells
      if (m(mi)%status .ge. 1) then

          ! test to see whether this category is involved in advection at all
          if (adv(mi)%snow(sc,ic)%z .gt. 0) then

              ! initialize temporary holding array
              temp3d(:,mi) = d0_  ! vector scalar assignment

              ! find full EASE domain reference for current cell
              i_ease = m(mi)%grid_h
              j_ease = m(mi)%grid_v

              if(mi .eq. 20513) then
                 testvar = -1
              endif
              if(mi .eq. 20035) then
                 testvar = -1
              endif

              do jjj=1,9

                  a_scale = d0_
                  if (jjj .eq. 9) then
                      mi1 = mi
                      if (ice(sc,ic,mi1)%snow%z .gt. 0) then
                          a_scale = (1.-adv(mi)%out)*ice(sc,ic,mi)%af
                      endif
                  elseif (adv(mi)%in(jjj) .gt. d0_) then
                      mi1 = m(mi)%mia(jjj)
                      if (ice(sc,ic,mi1)%snow%z .gt. 0) then
                          a_scale = adv(mi)%in(jjj)*ice(sc,ic,mi1)%af
                      endif
                  endif

                  if (a_scale .gt. d0_) then

!                      if (adv(mi)%snow(sc,ic)%sh_prev .gt. adv(mi)%snow(sc,ic)%depth) then
                          h_scale = ice(sc,ic,mi1)%snow%depth/adv(mi)%snow(sc,ic)%sh_prev
!                       else
!                          h_scale = ice(sc,ic,mi1)%snow%depth/adv(mi)%snow(sc,ic)%depth
!                      endif
                      th_old = ice(sc,ic,mi1)%snow%th/h_scale
                      v_scale = a_scale*h_scale

                      tracer = &
#if s_tracer_var==1
                          ice(sc,ic,mi1)%snow%d
#elif s_tracer_var==2
                          ice(sc,ic,mi1)%snow%heat
#elif s_tracer_var==3
                          ice(sc,ic,mi1)%snow%melt
#endif

                      ! remap tracers to internal ice layers
                      ! ------------------------------------------------------
                      new_layer_bot = d0_

                      do ii=1,adv(mi)%snow(sc,ic)%z

                          ! re-adjust new layer bottom for next snow layer temp adjustment
                          new_layer_top = new_layer_bot
                          new_layer_bot = new_layer_bot + adv(mi)%snow(sc,ic)%th_new(ii) ! th_new comes from snow_grid calc

                          ! find 1st old layer that contains part of new layer, going down
                          jj = 0
                          old_layer_bot = d0_
                          do while((new_layer_top .ge. old_layer_bot) .and. (jj .lt. z_snow_max))
                               jj=jj+1
                               ! old thickness...
                               old_layer_bot = old_layer_bot + th_old(jj)
                           enddo
                           old_layer_top = old_layer_bot - th_old(jj)
                           ! now jj is OLD layer where NEW layer ii starts...

                           do while ((old_layer_top .lt. new_layer_bot) .and. (jj .le. z_snow_max))

                               ! ----- NEW LAYER GEOMETRIES ------------------------------
                               interim_top = max(old_layer_top,new_layer_top)
                               interim_bot = min(old_layer_bot,new_layer_bot)

                               z1 = interim_top - old_layer_top
                               z2 = interim_bot - old_layer_top
                               dz = z2-z1

                               !if (jj .eq. ice(sc,ic,mi)%snow%z .and. ii .eq. adv(mi)%snow(sc,ic)%z) then
                               !    if ((dz_total + dz) .lt. adv(mi)%snow(sc,ic)%th_new(ii)) then
                               !        dz = dz + (adv(mi)%snow(sc,ic)%th_new(ii)-(dz_total + dz))
                               !    endif
                               !endif

                               ! record tracer
                               temp3d(ii,mi) = temp3d(ii,mi) + tracer(jj)*dz*v_scale

                               ! ----- SETUP VARIABLES FOR NEXT PARTIAL LAYER ------------
                               ! keeping track of emerging new layer thickness, for
                               dz_total = dz_total + dz

                               jj=jj+1
                               if (jj .le. z_snow_max) then
                                   ! find boundaries of next old layer
                                   old_layer_top = old_layer_bot
                                   old_layer_bot = old_layer_bot + th_old(jj)
                               endif

                           enddo  ! end of layer ii mapping

                      enddo ! end of ii snow layer loop

                  endif ! end of v_scale test (is there and volume advecting in from category?)

              enddo ! end of 8 neighbor cells + local cell loop

              ! find final tracer concentrations by dividing by new advective
              ! area and thickness, and passing back into local cell
              tmp1 = 1.d0/adv(mi)%ice(sc,ic)%af
!              if (adv(mi)%snow(sc,ic)%sh_prev .gt. adv(mi)%snow(sc,ic)%depth) then
                  tmp1 = tmp1*adv(mi)%snow(sc,ic)%sh_prev/adv(mi)%snow(sc,ic)%depth
!              endif
              do ii=1,adv(mi)%snow(sc,ic)%z
                  temp3d(ii,mi) = temp3d(ii,mi)*tmp1/adv(mi)%snow(sc,ic)%th_new(ii)
              enddo

          endif ! end advection-active category test


      endif ! end status test
      enddo ! end mi loop

      ! copy back into original variable
          do mi=1,tcells
              if (m(mi)%status .ge. 1) then
                  if (adv(mi)%ice(sc,ic)%af .gt. d0_ .and. adv(mi)%snow(sc,ic)%z .gt. 0) then

                  do ii=1,adv(mi)%snow(sc,ic)%z
#if s_tracer_var==1
!                              if (abs(400000 - temp3d(ii,mi)) .gt. 0.1) then
!                                  print *,mi,ic
!                                  print *,ice(sc,ic,mi)%snow%d
!                                  print *,temp3d(:,mi)
!                                  testvar = -1
!                              endif
                              ice(sc,ic,mi)%snow%d(ii) &
#elif s_tracer_var==2
                              ice(sc,ic,mi)%snow%heat(ii) &
#elif s_tracer_var==3
                              ice(sc,ic,mi)%snow%melt(ii) &
#endif
                      = temp3d(ii,mi)
                  enddo

              endif ! end advection-active category test
          endif ! end status test
      enddo ! end mi loop
      enddo ! end snow category loop
      enddo ! end ice category loop

      ! deallocate temp3d array
      deallocate(temp3d)

