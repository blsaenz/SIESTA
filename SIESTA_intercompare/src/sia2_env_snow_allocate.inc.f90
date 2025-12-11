! ======================================================================
! START: Allocate snowfall/prescribed melt
! ======================================================================


!      f(mi)%sh_interp = 15.
!      f(mi)%sh_interp_next = 15.

          snow_model_temp = snow_model

          tmp1 = 0.


          ! From probability of snow blowing plots (Li and Pomery 1997)
          if (f(mi)%ws .eq. c9999) then
              v10 = sqrt(f(mi)%u10**2+f(mi)%v10**2)
          else
              v10 = f(mi)%ws
          endif
          !print *,'V10: ',v10

          v10 = min(7.0,max(0.d0,v10-8.d0))/7.d0

          do ic=1,ida_n
              sc=1
              if (ice(sc,ic,mi)%af .gt. 0.d0 .and. ice(sc,ic,mi)%snow%z .gt. 0.d0) then

                  if (ice(sc,ic,mi)%snow_dist .lt. dble(sd_num) .and. &
                      ice(sc,ic,mi)%snow%melt(1) .lt. 0.5d0) then

                      ice(sc,ic,mi)%snow_dist = ice(sc,ic,mi)%snow_dist + &
                          v10*dt/snow_rd_lim * dble(sd_num)/2.d0

                      ! add to snow cummulaitve redistribution
!                      ice(sc,ic,mi)%snow_rd = ice(sc,ic,mi)%snow_rd + v10*dt

!                      if (ice(sc,ic,mi)%snow_rd .gt. snow_rd_lim) then
                          ! reset distribution
!                          ice(sc,ic,mi)%snow_dist = dble(sd_num)
!                          ice(sc,ic,mi)%snow_rd = 0.d0
!                      endif

!                  else

!                      ice(sc,ic,mi)%snow_rd = 0.d0

                  endif

                  ice(sc,ic,mi)%snow_dist = &
                      min(ice(sc,ic,mi)%snow_dist,dble(sd_num))

              endif
          enddo


          r_depth = 0.
          af_total = 0.
          do ic=1,ida_n
          do sc=1,sda_n
              if (ice(sc,ic,mi)%af .gt. 0.d0) then

                  if (sda_n .gt. 1 .and. sc .eq. 1 .and. ice(sc,ic,mi)%snow_dist .lt. dble(sd_num)) then
                      tmp2 = 0.d0
                      tmp1 = ice(sc,ic,mi)%snow_dist
                      do jj=1,sd_num
                          tmp2 = tmp2 + min(1.d0,tmp1)*sda_multiplier(jj)
                          tmp1 = max(0.d0,tmp1-1.d0)
                      enddo
                      if (tmp2 .gt. 0.d0) then
                          sh_mean = ice(sc,ic,mi)%sh_prev* &
                             ice(sc,ic,mi)%snow_dist/tmp2
                      else
                          sh_mean = 0.d0
                      endif
                  else
                      sh_mean = ice(sc,ic,mi)%snow%depth
                  endif


                  if (snow_model .eq. 0) then
                      sh_diff(sc,ic) = f(mi)%sh_interp*c_01 - (sh_mean+ice(sc,ic,mi)%sh_offset)
                      ! above was modified to become below for draining surface layers - doen't work for flooding though?
                      !sh_diff(sc,ic) = f(mi)%sh_interp*c_01+ice(sc,ic,mi)%sh_offset - (ice(sc,ic,mi)%sh_prev))
                  elseif (snow_model .eq. 1) then
                      sh_diff(sc,ic) = f(mi)%pr*dt_s/snowfall_d*1000. ! kg/m^2/s * s * m^3/kg = m

                      m(mi)%pr_clim = m(mi)%pr_clim + sh_diff(sc,ic)*ice(sc,ic,mi)%af

                      if (f(mi)%sh_interp .ne. c9999 .and. f(mi)%sh_interp_next .ne. c9999 .or. &
!                      f(mi)%sh_interp*c_01 .gt. max_snow_h .or. tmp2 .gt. -1.d0) then
                      f(mi)%sh_interp*c_01 .le. max_snow_h) then
                          m(mi)%pr_ssmi = m(mi)%pr_ssmi + &
                            (f(mi)%sh_interp_next*c_01 - sh_mean)*ice(sc,ic,mi)%af
                      endif

                  elseif (snow_model .ge. 2) then
                      tmp1 = f(mi)%pr*dt_s/snowfall_d*1000.
                      tmp2 = ice(sc,ic,mi)%snow%ts
                      do ii=1,ice(sc,ic,mi)%snow%z
                         tmp2 = max(tmp2,ice(sc,ic,mi)%snow%t(ii))
                      enddo
                      if (f(mi)%sh_interp .eq. c9999 .or. f(mi)%sh_interp_next .eq. c9999 .or. &
!                      f(mi)%sh_interp*c_01 .gt. max_snow_h .or. tmp2 .gt. -1.d0) then
                      f(mi)%sh_interp*c_01 .gt. max_snow_h) then
                          ! use ncep precip
                          sh_diff(sc,ic) = tmp1
                          snow_model_temp = 1
                      else
                          ! use
                          tmp2 = f(mi)%sh_interp_next*c_01 - sh_mean
                          ! sh_increment = max(0.,tmp2) ! limit to snow growth only
                          sh_diff(sc,ic) = tmp2

                          m(mi)%pr_clim = m(mi)%pr_clim + tmp1*ice(sc,ic,mi)%af
                          m(mi)%pr_ssmi = m(mi)%pr_ssmi + tmp2*ice(sc,ic,mi)%af
                      endif
                  endif

                  ! sum differences to find mean snow difference r_depth
                  r_depth = r_depth + sh_diff(sc,ic)*ice(sc,ic,mi)%af
                  af_total = af_total + ice(sc,ic,mi)%af

              else
                  sh_diff(sc,ic) = 0.
                  sh_increment(sc,ic) = 0.
              endif

          enddo
          enddo

          ! find mean snow difference
          r_depth = r_depth/af_total

          !print *,'Snow: ',f(mi)%pr,r_depth

 !         if (ida_n .eq. 1) then
              ! new snow is simply the difference if only one ice category
!              sh_increment(1,1) = sh_diff(1,1)/dt_s
!          else

              ! init housekeep vars for loop below - sh_diff is used differently here,
              ! sort of as a valid stamp for r_depth allocation
              sh_diff = 1. ! vector assignment
              keep_growing = .true.

              ! find individual category adjustments
              do while (keep_growing)
                  tmp2 = r_depth ! find current snow depth change target
                  keep_growing = .false.
                  do ic=1,ida_n
                  do sc=1,sda_n
                      if (ice(sc,ic,mi)%af .gt. 0. .and. sh_diff(sc,ic) .gt. 0.) then

                          if (snow_model .eq. 0) then
                              tmp1 = ice(sc,ic,mi)%sh_prev+ice(sc,ic,mi)%sh_offset
                          else
                              tmp1 = ice(sc,ic,mi)%sh_prev
                          endif

                          if ((tmp2 + tmp1) .lt. 0.)  then
                              ! snowdepth reduction exceeds snowdepth - reallocate to other categories if they exist
                              kk = 0
                              if (af_total .ne. ice(sc,ic,mi)%af) then
                                  if (sc .lt. sda_n) then
                                      do jj=sc+1,sda_n
                                          if (ice(jj,ic,mi)%sh_prev .gt. 0.) then
                                              kk = 1
                                          endif
                                      enddo
                                  endif
                                  if (ic .lt. ida_n) then
                                  do ii=ic+1,ida_n
                                  do jj=1,sda_n
                                      if (ice(jj,ii,mi)%sh_prev .gt. 0.) then
                                          kk = 1
                                      endif
                                  enddo
                                  enddo
                                  endif
                              endif
                              if (kk .eq. 1) then
                                  r_depth = r_depth + abs(tmp1 + tmp2)*ice(sc,ic,mi)%af/(af_total - ice(sc,ic,mi)%af)
                                  keep_growing = .true.
                              endif
                              sh_diff(sc,ic) = 0.
                              sh_increment(sc,ic) = -1.*tmp1/dt_s
                          else
                              ! everything is fine - use current snowdepth change
                              sh_diff(sc,ic) = tmp2
                              sh_increment(sc,ic) = tmp2/dt_s
                          endif
                      endif
                  enddo  ! end of snow cateogry iteration
                  enddo  ! end of ice cateogry iteration
              enddo ! end while loop

!          endif ! end of multiple ice category check

          print *,'sh_diff: ',sh_diff
          print *,'sh_increment: ',sh_increment











! ======================================================================
! END: Allocate snowfall/prescribed melt
! ======================================================================

