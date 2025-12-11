      ! allocate temp3d array
      allocate(temp3d(1,tcells))

      ! advect into temporary holding variable (temp3d)
      do ic=1,ida_n
      do sc=1,sda_n
      do mi=1,tcells
      if (m(mi)%status .ge. 1) then

          ! test to see whether this category is involved in advection at all
          if (adv(mi)%ice(sc,ic)%af .gt. af_min) then

              ! initialize temporary holding array
              temp3d(1,mi) = d0_

              ! find full EASE domain reference for current cell
              i_ease = m(mi)%grid_h
              j_ease = m(mi)%grid_v

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

                      temp3d(1,mi) = temp3d(1,mi) + a_scale* &
#if tracer1_var==1
                          ice(sc,ic,mi1)%age
#elif tracer1_var==2
                          ice(sc,ic,mi1)%ridged
#elif tracer1_var==3
                          ice(sc,ic,mi1)%snow_dist
#elif tracer1_var==4
                          ice(sc,ic,mi1)%snow_rd
#elif tracer1_var==5
                          ice(sc,ic,mi1)%no3(ml_z)
#elif tracer1_var==6
                          ice(sc,ic,mi1)%nh4(ml_z)
#elif tracer1_var==7
                          ice(sc,ic,mi1)%sioh4(ml_z)
#elif tracer1_var==8
                          ice(sc,ic,mi1)%poc(ml_z)
#endif
                  endif

              enddo ! end of 8 neighbor cells + local cell loop

              ! find final tracer concentrations by dividng by new advective area
			  temp3d(1,mi) = temp3d(1,mi)/adv(mi)%ice(sc,ic)%af

          endif ! end advection-active category test


	  endif ! end status test
	  enddo ! end mi loop

	  ! copy back into original variable
	  do mi=1,tcells
	      if (m(mi)%status .ge. 1) then
	          if (adv(mi)%ice(sc,ic)%af .gt. 0.) then
#if tracer1_var==1
                  ice(sc,ic,mi)%age &
#elif tracer1_var==2
                  ice(sc,ic,mi)%ridged &
#elif tracer1_var==3
                  ice(sc,ic,mi)%snow_dist &
#elif tracer1_var==4
                  ice(sc,ic,mi)%snow_rd &
#elif tracer1_var==5
                  ice(sc,ic,mi)%no3(ml_z) &
#elif tracer1_var==6
                  ice(sc,ic,mi)%nh4(ml_z) &
#elif tracer1_var==7
                  ice(sc,ic,mi)%sioh4(ml_z) &
#elif tracer1_var==8
                  ice(sc,ic,mi)%poc(ml_z) &
#endif
            		  = temp3d(1,mi)

                  endif ! end advection-active category test
              endif ! end status test
          enddo ! end mi loop
          enddo ! end snow category loop
          enddo ! end ice category loop

          ! deallocate temp3d array
          deallocate(temp3d)
