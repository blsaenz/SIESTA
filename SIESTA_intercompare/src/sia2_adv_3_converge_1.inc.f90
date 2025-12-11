
! CASE 1 - accept all convergence/divergence
! ---------------------------------------------------------------------------------

      v_mod = 0  ! init velocity modification counter

!$OMP DO SCHEDULE(DYNAMIC,chunk_size)

      do mi=1,tcells

          if (m(mi)%status .eq. 1) then

              if (iii .eq. 0) then              

             
                  ! start with zero ice vectors, and work up from there
                  adv(mi)%v_scale = d0_
                  adv(mi)%v_incr = v_incr

                  if (m(mi)%status .ge. 1 .and. adv(mi)%a_convg .gt. 0.) then
 
                      ! melt ice, instead of converging, if surface temp is 
                      ! greater than -2.5
                      ! --------------------------------------------------------
                      do ic=1,ida_n
                          if (ice(ic,mi)%af .gt. 0 .and. ice(ic,mi)%t(1) .gt. -2.5d0) then
        
                              ! find new ice category area, after advection (tmp1)
                              tmp1 = (1.-adv(mi)%out)*ice(ic,mi)%af 
                              do jj=1,8
                                  mi1 = m(mi)%mia(jj)                      
                                  tmp1 = tmp1 + adv(mi)%in(jj)*ice(ic,mi1)%af
                              enddo
                              tmp1 = tmp1*cell_area  
                              
                              ! set category up for melting insead of advection
                              tmp2 = max(d0_,adv(mi)%a_convg - adv(mi)%a_melt)  ! a_convg not yet melted via smaller category
                              adv(mi)%a_melt = adv(mi)%a_melt + min(tmp2,tmp1)
                          endif
                      enddo
    
                      ! find advective input that might can account for convergence
                      tmp1 = min(adv(mi)%a_convg - adv(mi)%a_melt,adv(mi)%a_in)
                      
                      ! find input convergence allowed based on ice stregnth
                      tmp2 = min(1.d0,adv(mi)%p/108000.d0)  ! 4m * 270000 Nm is max
                      
                      ! sum allowed input convergence and internel convergence to
                      ! find max allowed areal convergence
                      adv(mi)%a_allowed = adv(mi)%a_convg + tmp2*(1.d0-tmp1)

                  endif

              else


!                 if (adv(mi)%a_allowed .gt. d0_) then

                     if (adv(mi)%a_convg .gt. adv(mi)%a_allowed) then 
    
                          do jj=1,8
                              if (adv(mi)%in(jj) .gt. d0_) then
                                  mi1 = m(mi)%mia(jj)
                                  adv(mi1)%v_incr = -1.d0*v_incr
                               endif
                           enddo

                      endif

 !                 endif
               
              endif

 
          endif ! end of status check

      enddo ! end of tcell interation

!$OMP END DO
!$OMP END PARALLEL

    !print *,'v_mod: ',v_mod

    ! see if convergence calc passes test
    if (iii .eq. v_incr_count) then
        do_adv_calc = .false.
    endif

! ---------------------------------------------------------------------------------
! CASE 1 - End
