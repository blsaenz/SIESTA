
! CASE 1 - accept all convergence/divergence
! ---------------------------------------------------------------------------------

      v_mod = 0  ! init velocity modification counter

!$OMP DO SCHEDULE(DYNAMIC,chunk_size)

          do mi=1,tcells
             if (m(mi)%status .ge. 1 .and. adv(mi)%a_convg .gt. 0.) then
                                    
                  ! reduce incoming u & v using the "momentum"/strength ratio: if strength
                  ! is greater than "momentum," where "momentum" = momentum*constant,
                  ! and the constant is a tuned parameter
                  
                  ! momentum = kg * m / s                  
                  ! ice strength = N * m = kg * m / s^2
                  
                  ! so we could turn the momentum into a force by applying it over a 
                  ! period of seconds, but that is just a scalar factor anyway,
                  ! so we will just use one scalar factor

                  ! find compsite ice height that is effected by convergence                  
                   


                  do jj=1,8
                      if (adv(mi)%in(jj) .gt. 0.) then
                          mi1 = m(mi)%mia(jj)
                          u = f(mi1)%ivu_interp/1.e6*dt_s*adv(mi1)%v_scale ! km/timestep
                          v = f(mi1)%ivv_interp/1.e6*dt_s*adv(mi1)%v_scale! km/timestep
                          ! km/timestep * g/m^3 * 1kg/1000g * 1000m/1km * 1timestep/?s = kg*m/s
                          tmp1 = sqrt(u**2 + v**2)*adv(mi)%mass/dt_s ! momentum (kg*m/s)

                          tmp2 = min(1.d0,adv(mi)%p/4.)  ! find ice height multiplier - at 10m we are completely resiting convergence

                          ! reduce speed of incoming ice, count velocity mod
!                          if (tmp2 .gt. tmp1*p_factor) then
                              tmp3 = min(adv(mi)%a_in,adv(mi)%a_convg)/adv(mi)%a_in  ! percentage of convergence due to input
                              tmp3 = tmp3*adv(mi)%in1a(jj)/adv(mi)%a_in      ! percentage of input resulting in convergence
                              ! assigning tmp3 to vscale to reduce velocity should result in prevention of convergence now
                              ! now allow some convergence based on ice strength, and modify tmp3 by ice stength based tmp2                             
!                              tmp2 = tmp1*p_factor / adv(mi)%p
!                              tmp2 = (1.d0-tmp2*tmp3*0.33d0)  ! add 0.33 factor to allow proper convergence 
                              tmp2 = (1.-tmp3)*tmp2  ! 
                              tmp2 = max(0.d0,(1.-tmp2))
                              if (tmp2 .lt. adv(mi1)%v_scale) then
                                  adv(mi1)%v_scale = tmp2
                                  v_mod = v_mod + 1
                              endif
!                          endif
                      endif
                  enddo

              endif ! end of status check
          enddo ! end of tcell interation

!$OMP END DO
!$OMP END PARALLEL

    !print *,'v_mod: ',v_mod

    ! see if convergence calc passes test
    if (v_mod .lt. 1) then
        do_adv_calc = .false.
    endif

! ---------------------------------------------------------------------------------
! CASE 1 - End
