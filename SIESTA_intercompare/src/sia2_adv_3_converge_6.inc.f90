
! CASE 1 - accept all convergence/divergence
! ---------------------------------------------------------------------------------

      v_mod = 0  ! init velocity modification counter

!$OMP DO SCHEDULE(DYNAMIC,chunk_size)

      do mi=1,tcells

          if (m(mi)%status .ge. 1) then

              ! initialize
              adv(mi)%a_melt = d0_
              adv(mi)%a_allowed = d0_

              ! total area and mean ice depth for each ICE catagory across SNOW categories
              if (mi .eq. 20737) then
                   af_ic = 0.d0
              endif
              af_ic = 0.d0
              id_ic = 0.d0
              do ic=1,ida_n
                  do sc=1,sda_n
                      if (ice(sc,ic,mi)%af .gt. 0.d0) then
                      af_ic(ic) = af_ic(ic) + ice(sc,ic,mi)%af
                      id_ic(ic) = id_ic(ic) + ice(sc,ic,mi)%id(ice(sc,ic,mi)%z-z_sk)*ice(sc,ic,mi)%af
                      endif
                  enddo
                  if (af_ic(ic) .gt. 0.d0) then
                     id_ic(ic) = id_ic(ic)/af_ic(ic)
                  endif
              enddo

              if (adv(mi)%a_convg .gt. d0_) then

                  ! melt ice, instead of converging, if mean temp is
                  ! greater than -2.5
                  ! --------------------------------------------------------
                  ic=1
                  do sc=1,sda_n
                      if (ice(sc,ic,mi)%af .gt. d0_ .and. adv(mi)%ice(sc,ic)%t_mean .ge. -2.0d0) then

                          ! find new ice category area, after advection (tmp1)
                          tmp1 = (1.-adv(mi)%out)*ice(sc,ic,mi)%af
                          do jj=1,8
                              mi1 = m(mi)%mia(jj)
                              if(mi1 .gt. 0) then
                                tmp1 = tmp1 + adv(mi)%in(jj)*ice(sc,ic,mi1)%af
                              endif
                          enddo
                          tmp1 = tmp1*cell_area

                          ! set category up for melting insead of advection
                          tmp2 = max(d0_,adv(mi)%a_convg - adv(mi)%a_melt)  ! a_convg not yet melted
                          adv(mi)%a_melt = adv(mi)%a_melt + min(tmp2,tmp1)
                      endif
                  enddo

                  ! find input convergence percentage allowed, based on ice strength (depth in this case)
                  tmp2 = min(1.d0,max(d0_,((adv(mi)%id_mean-2.d0)/3.d0)))  ! linear decrease over 2-5m

                  ! find max allowed areal convergence
                  adv(mi)%a_allowed = (adv(mi)%a_convg - adv(mi)%a_melt) * (1.d0-tmp2)
                  tmp3 = adv(mi)%a_convg - adv(mi)%a_melt - adv(mi)%a_allowed  ! area not allowed

                  if (tmp3 .gt. 0.d0) then

                      ! prevent un-allowed convergence by retaining ice area, if less than 100%
                      tmp1 = adv(mi)%a_new + tmp3  ! new ice area w/ un-allowed convergence added back in
                      adv(mi)%a_convg = adv(mi)%a_convg - tmp3 + max(0.d0,tmp1-cell_area)

                     ! if convergence is still happening after correction by ice area, then
                     tmp3 = adv(mi)%a_convg - adv(mi)%a_melt - adv(mi)%a_allowed  ! area not allowed
                     if (tmp3 .gt. 0.d0) then

                          ! drop it, if we can't figure out what to do with it
                          adv(mi)%a_drop = adv(mi)%a_drop + tmp3

                     endif

                  endif

              endif ! end of convergence check

          endif ! end of status check

      enddo ! end of tcell interation

!$OMP END DO
!$OMP END PARALLEL

    !print *,'v_mod: ',v_mod

    ! see if convergence calc passes test
    if (v_mod .le. 1) then
        do_adv_calc = .false.
    endif

! ---------------------------------------------------------------------------------
! CASE 1 - End
