
! CASE 1 - accept all convergence/divergence
! ---------------------------------------------------------------------------------

      v_mod = 0  ! init velocity modification counter

!$OMP DO SCHEDULE(DYNAMIC,chunk_size)

      do mi=1,tcells

          if (m(mi)%status .ge. 1) then 
              
              ! initialize
              adv(mi)%a_melt = d0_
              adv(mi)%a_allowed = d0_          
          
              if (adv(mi)%a_convg .gt. 0.) then
              
                  ! melt ice, instead of converging, if mean temp is 
                  ! greater than -2.5
                  ! --------------------------------------------------------
                  do ic=1,ida_n
                      if (ice(ic,mi)%af .gt. 0 .and. adv(mi)%ice(ic)%t_mean .ge. -2.5d0) then
    
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
                  
                  ! find input convergence percentage allowed, based on ice strength (depth in this case)
                  tmp2 = min(1.d0,max(0.d0,((adv(mi)%id_mean-2.d0)/2.d0)))  ! linear decrease over 1-4m 
                  
                  ! sum allowed input convergence and internel convergence to
                  ! find max allowed areal convergence
                   adv(mi)%a_allowed = adv(mi)%a_convg + tmp1*(tmp2 - 1.d0)
    
                  !
                  ! adv(mi)%a_allowed = (adv(mi)%a_convg - adv(mi)%a_melt) * (1.d0-tmp2)
    
    
                  if ((adv(mi)%a_convg - adv(mi)%a_melt) .gt. adv(mi)%a_allowed) then
                  
                     ! incrementally reduce input ice velocities, if below max reduction
                      do jj=1,8
                          if (adv(mi)%in(jj) .gt. 0.d0) then
                               mi1 = m(mi)%mia(jj)
        
                              if (adv(mi1)%v_scale .gt. d0_) then
                                  
                                  adv(mi1)%v_scale = max(0.d0,adv(mi1)%v_scale - v_incr)
                                  v_mod = v_mod + 1                          
        
                              endif                        
    
                          endif
                                             
                      enddo
    
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
