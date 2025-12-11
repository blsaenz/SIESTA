! file: SIA2_advect_setup.F90
! Sea Ice Algae Model 2 - Saenz & Arrigo
! Version beta
! ======================================================================

     SUBROUTINE sia2_advect_setup(pj,ease_f,f,m,ice,adv)

          use sia2_globals

          implicit none

          ! Function Arguments
          ! ------------------------------------------------------------
          type (proj_type) :: pj
          type (ease_f_type) :: ease_f
          type (forcing_type) :: f(tcells)
          type (meta_type) :: m(tcells)
          type (ice_type) :: ice(sda_n,ida_n,tcells)
          type (adv_type) :: adv(tcells)

          ! shared internal variables
          !----------------------------------------------------
          integer :: chunk_size,max_d_exceeded,v_mod,v_incr_count
          double precision :: v_incr

          ! private internal variables
          ! ------------------------------------------------------------
          logical :: keep_growing,maxed_depth,do_adv_calc,ignore_h_max
          integer :: ii,jj,iii,ic,sc,mi,mi1,i_ease,j_ease,int_z,sk_1,sk_z, &
             int_z_new,z_last,z_snow_new
          integer, dimension(3) :: out_case
          integer, dimension(8) :: mi_adv,i_adv,j_adv
          double precision :: u,v,flag,lost,sklost,smalg_lost,tmp1,tmp2, &
              tmp3,tmp4,bv_mean,r_depth,layer_divisor
          double precision, dimension(z_max) :: th_new,id_new
          double precision,  dimension(ida_n) :: id_ic,af_ic

           ! chunk size calc
          chunk_size = 16

          ! v_scale increment size (fractional)
          v_incr = 0.1d0
          v_incr_count = int(1.d0/v_incr)

          ! init over-max-depth counter
          max_d_exceeded = 0
          ignore_h_max = .false.

          do_adv_calc = .true.
          iii = 0

!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(keep_growing,maxed_depth,ignore_h_max, &
!$OMP ii,jj,iii,ic,sc,mi,mi1,i_ease,j_ease,int_z,sk_1,sk_z, &
!$OMP int_z_new,z_last,z_snow_new, &
!$OMP out_case, &
!$OMP mi_adv,i_adv,j_adv, &
!$OMP u,v,flag,lost,sklost,smalg_lost,tmp1,tmp2, &
!$OMP tmp3,tmp4,bv_mean,r_depth,layer_divisor, &
!$OMP th_new,id_new)

!$OMP DO SCHEDULE(DYNAMIC,chunk_size)
          do mi=1,tcells
              if (m(mi)%status .ge. 1) then
                  adv(mi)%v_scale = 1.d0  ! init velocity scale
                  adv(mi)%in1 = d0_ ! vector assignment
                  adv(mi)%a_melt = 0.d0  ! init velocity scale
                  adv(mi)%a_drop = 0.d0  ! init velocity scale
                  ! find t_mean
                  do ic=1,ida_n
                      adv(mi)%ice(:,ic)%t_mean = d0_
                      tmp1 = d0_
                      do sc=1,sda_n
                          if (ice(sc,ic,mi)%af .gt. d0_) then
                              do ii=1,ice(sc,ic,mi)%z-z_sk
                                  adv(mi)%ice(sc,ic)%t_mean = adv(mi)%ice(sc,ic)%t_mean &
                                    + ice(sc,ic,mi)%t(ii)*ice(sc,ic,mi)%th(ii)*ice(sc,ic,mi)%af
                                  tmp1 = tmp1 + ice(sc,ic,mi)%th(ii)*ice(sc,ic,mi)%af
                              enddo

							  ! set mass-balace paramaters to bulk volume-based concentrations to
							  ! simplify advection
							  ! --------------------------------------------------------
                              do ii=1,ice(sc,ic,mi)%z
								  bv_mean = ice(sc,ic,mi)%bv(ii)*c_001
								  ice(sc,ic,mi)%NO3(ii) = ice(sc,ic,mi)%NO3(ii)*bv_mean
								  ice(sc,ic,mi)%NH4(ii) = ice(sc,ic,mi)%NH4(ii)*bv_mean
								  ice(sc,ic,mi)%PO4(ii) = ice(sc,ic,mi)%PO4(ii)*bv_mean
								  ice(sc,ic,mi)%SiOH4(ii) = ice(sc,ic,mi)%SiOH4(ii)*bv_mean
								  ice(sc,ic,mi)%smalg(ii) = ice(sc,ic,mi)%smalg(ii)*bv_mean
								  ice(sc,ic,mi)%poc(ii) = ice(sc,ic,mi)%poc(ii)*bv_mean
                              enddo

                          endif
                      enddo
                      if (tmp1 .gt. d0_) then
                          adv(mi)%ice(1,ic)%t_mean = adv(mi)%ice(1,ic)%t_mean/tmp1
                      endif
                      if (sda_n .gt. 1) then
                          do sc=2,sda_n
                              adv(mi)%ice(sc,ic)%t_mean = adv(mi)%ice(1,ic)%t_mean
                          enddo
                      endif
                  enddo
               endif
          enddo
!$OMP END DO
!$OMP END PARALLEL

          do while (do_adv_calc)


!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(keep_growing,maxed_depth, &
!$OMP ii,jj,iii,ic,mi,mi1,i_ease,j_ease,int_z,sk_1,sk_z, &
!$OMP int_z_new,z_last,z_snow_new, &
!$OMP out_case, &
!$OMP mi_adv,i_adv,j_adv, &
!$OMP u,v,flag,lost,sklost,smalg_lost,tmp1,tmp2, &
!$OMP tmp3,tmp4,bv_mean,r_depth,layer_divisor, &
!$OMP th_new,id_new, &
!$OMP id_ic,af_ic)

!$OMP DO SCHEDULE(DYNAMIC,chunk_size)

! 1.
! =================================================================================


          do mi=1,tcells
              if (m(mi)%status .ge. 1) then

                  ! find full EASE domain reference for current cell
                  i_ease = m(mi)%grid_h
                  j_ease = m(mi)%grid_v

                  ! clear advection fractions, temporary vars
                  adv(mi)%in = d0_ ! vector assignment
                  adv(mi)%in1a = d0_ ! vector assignment
                  adv(mi)%out = d0_
                  adv(mi)%v_mod = 0

				  ! determine indexes to adjacent cells that current cell (mi)
				  ! may be advecting with
                  do jj=1,8
                      i_adv(jj) = m(mi)%ia(jj)
                      j_adv(jj) = m(mi)%ja(jj)
                      mi_adv(jj) = m(mi)%mia(jj)

                      ! find inputs into current cell
					  if (mi_adv(jj) .ne. 0) then
						  if (m(mi_adv(jj))%status .ge. 1) then

                              v = f(mi_adv(jj))%ivu_interp/1.e6* &
                                  dt_s*adv(mi_adv(jj))%v_scale ! km/timestep
                              u = f(mi_adv(jj))%ivv_interp/1.e6* &
                                  dt_s*adv(mi_adv(jj))%v_scale! km/timestep

							  flag = ease_f%ice_vec(i_adv(jj),j_adv(jj),3)
							  select case (jj)
								  case (1)
									  if (v .lt. 0.) then
										  adv(mi)%in(1) = (cell_side-abs(u))*abs(v)/cell_area ! (fraction of input cell)/timestep
									  endif
								  case (2)
									  if ((u .lt. 0.) .and. (v .lt. 0.)) then
										  adv(mi)%in(2) = abs(u)*abs(v)/cell_area ! (fraction of input cell)/timestep
									  endif
								  case (3)
									  if (u .lt. 0.) then
										  adv(mi)%in(3) = (cell_side-abs(v))*abs(u)/cell_area ! (fraction of input cell)/timestep
									  endif
								  case (4)
									  if ((u .lt. 0.) .and. (v .gt. 0.)) then
										  adv(mi)%in(4) = abs(u)*v/cell_area ! (fraction of input cell)/timestep
									  endif
								  case (5)
									  if (v .gt. 0.) then
										  adv(mi)%in(5) = (cell_side-abs(u))*v/cell_area ! (fraction of input cell)/timestep
									  endif
								  case (6)
									  if ((u .gt. 0.) .and. (v .gt. 0.)) then
										  adv(mi)%in(6) = u*v/cell_area ! (fraction of input cell)/timestep
									  endif
								  case (7)
									  if (u .gt. 0.) then
										  adv(mi)%in(7) = (cell_side-abs(v))*u/cell_area ! (fraction of input cell)/timestep
									  endif
								  case (8)
									  if ((u .gt. 0.) .and. (v .lt. 0.)) then
										  adv(mi)%in(8) = u*abs(v)/cell_area ! (fraction of input cell)/timestep
									  endif
							  end select

                              ! record original inputs
                              if (iii .eq. 0) then
                                  adv(mi)%in1(jj) = adv(mi)%in(jj)
                              endif

						  endif
                      endif
                  enddo  ! end of jj=1,8 do loop

! 2.
! =================================================================================
! all this checking is necessary because of two things: 1) if the ice vector points
! into land, the biomass is conserved in the ice and not lost in the landward
! direction. 2) biomass is lost if the non-land receiving cell does not have ice
! (status .eq. 0), and the lost biomass needs to be recorded

                  ! case: self (output)
                  ! OOO    xxO     Oxx                            812
                  ! OXO -> xXO     OXx                            7 3
                  ! OOO    OOO, or OOO, etc... where case below = 654
                  if (ease_f%ice_vec(i_ease,j_ease,3) .gt. d0_) then

                      ! find scaled u,v
                      v = f(mi)%ivu_interp/1.e6*dt_s*adv(mi)%v_scale ! km/timestep
                      u = f(mi)%ivv_interp/1.e6*dt_s*adv(mi)%v_scale! km/timestep

                      lost = d0_  ! init temporary biomass loss variable to zero
                      sklost = d0_  ! init temporary biomass loss variable to zero

                      if ((u .gt. d0_) .and. (v .lt. d0_)) then
                          !               000
                          ! cases 3,4,5:  0 3
                          !               054
                           out_case(1:3) = (/3, 4, 5/)
                           !out_case(1) = 4
                           !out_case(1) = 5
                      elseif ((u .lt. d0_) .and. (v .lt. d0_)) then
                          !               000
                          ! cases 5,6,7:  7 0
                          !               650
                           out_case(1:3) = (/5, 6, 7/)
                      elseif ((u .lt. d0_) .and. (v .gt. d0_)) then
                          !               810
                          ! cases 1,7,8:  7 0
                          !               000
                           out_case(1:3) = (/7, 8, 1/)
                      elseif ((u .gt. d0_) .and. (v .gt. d0_)) then
                          !               012
                          ! cases 1,2,3:  0 3
                          !               000
                           out_case(1:3) = (/1, 2, 3/)
                      elseif (u .eq. d0_) then
                          !               010
                          ! cases 5 or 1: 0 0
                          !               050
                          if (v .lt. 0) then
                              out_case(1:3) = (/5, 0, 0/)
                          elseif (v .ge. d0_) then
                              out_case(1:3) = (/1, 0, 0/)
                          endif
                      elseif (v .eq. d0_) then
                          !               000
                          ! cases 3 or 7: 7 3
                          !               000
                          if (u .ge. 0) then
                              out_case(1:3) = (/3, 0, 0/)
                          elseif (u .lt. 0) then
                              out_case(1:3) = (/7, 0, 0/)
                          endif
                      endif

                      do jj=1,3
                          select case (out_case(jj))
                              case (1)
                                  tmp1 = (cell_side-abs(u))*v
                              case (2)
                                  tmp1 = u*v
                              case (3)
                                  tmp1 = (cell_side-abs(v))*u
                              case (4)
                                  tmp1 = u*abs(v)
                              case (5)
                                  tmp1 = (cell_side-abs(u))*abs(v)
                              case (6)
                                  tmp1 = abs(u)*abs(v)
                              case (7)
                                  tmp1 = (cell_side-abs(v))*abs(u)
                              case (8)
                                  tmp1 = abs(u)*v
                              case default  ! zero case
                                  tmp1 = d0_
                          end select

                          if (tmp1 .gt. d0_) then
                              if ((i_adv(jj) .gt. mdh1) .and. (j_adv(jj) .gt. mdv1) .and. &
                              (i_adv(jj) .le. mdh2) .and. (j_adv(jj) .le. mdv2)) then
                                  if(pj%mask(i_adv(jj),j_adv(jj)) .eq. 0) then
                                      tmp1 = 0. ! land mask - can't advect out into land
                                  elseif (mi_adv(out_case(jj)) .eq. 0) then
                                      lost = lost + tmp1 ! cell is not a part of the model domain
                                  elseif(m(mi_adv(out_case(jj)))%status .eq. 0) then
                                      lost = lost + tmp1 ! cell is in model domain, but not active
                                  endif
                              else
                                  tmp1 = 0.  ! out of model domain
                              endif
                          endif

                          ! record the areal fraction leaving the current cell
                          adv(mi)%out = adv(mi)%out + tmp1 ! (fraction of cell)/timestep

                      enddo

                      adv(mi)%out = adv(mi)%out/cell_area  ! (fraction of cell)/timestep

                  endif ! end of case: self


                  ! determine if any advection is occuring, to protect against divide by zero later
                  adv(mi)%adv_sum = adv(mi)%in(1)+adv(mi)%in(2)+adv(mi)%in(3) &
                      +adv(mi)%in(4)+adv(mi)%in(5)+adv(mi)%in(6) &
                      +adv(mi)%in(7)+adv(mi)%in(8)+adv(mi)%out

                  ! sanity check for advection percentages - must always be positive
                  do ii=1,8
                     if (adv(mi)%in(ii) .lt. 0.) then
                        print *,'Bad NEGATIVE Advection -',i_ease,j_ease,'In:',adv(mi)%in, &
                        'Out:',adv(mi)%out
                        call exit(0)
                     endif
                  enddo
                  if (adv(mi)%out .lt. 0.) then
                      print *,'Bad NEGATIVE Advection -',i_ease,j_ease,'In:',adv(mi)%in, &
                      'Out:',adv(mi)%out
                      call exit(0)
                  endif

                  ! record area statistics for divergence/convergence accounting
                  ! use sum of current categories so we don't end up
                  ! drifting off the measured concentrations
                  adv(mi)%af = 0.
                  do ic=1,ida_n
                  do sc=1,sda_n
                      adv(mi)%af = adv(mi)%af + ice(sc,ic,mi)%af
                  enddo
                  enddo
                  adv(mi)%af_open =  1.d0 - adv(mi)%af ! find current open water fraction
                  adv(mi)%a = adv(mi)%af*cell_area  ! find current total ice fraction

                  adv(mi)%a_new = f(mi)%ic_interp_next*cell_area
                  adv(mi)%a_out = adv(mi)%af*adv(mi)%out*cell_area
                  adv(mi)%a_open = cell_area - adv(mi)%a
                  adv(mi)%a_in = 0.d0

                  do ii=1,8
                      if (adv(mi)%in(ii) .gt. 0.) then
                          tmp1 = 0.
                          do ic=1,ida_n
                          do sc=1,sda_n
                              tmp1 = tmp1 + ice(sc,ic,mi_adv(ii))%af
                          enddo
                          enddo
                          adv(mi)%in1a(ii) = adv(mi)%in(ii)*cell_area*tmp1
                          adv(mi)%a_in = adv(mi)%a_in + adv(mi)%in1a(ii)
                      endif
                  enddo

                  ! find composite ice stength
                  tmp1 = d0_   ! mean height var
                  tmp2 = d0_   ! mean area var
                  do ic=1,ida_n
                  do sc=1,sda_n
                      if (ice(sc,ic,mi)%af .gt. 0.) then
                          int_z = ice(sc,ic,mi)%z - z_sk
                          tmp1 = tmp1 + ice(sc,ic,mi)%id(int_z)*ice(sc,ic,mi)%af
                          tmp2 = tmp2 + ice(sc,ic,mi)%af
                      endif
                  enddo
                  enddo
                  tmp1 = tmp1/tmp2
                  adv(mi)%p = 27000.*tmp1*exp(-20.*(1.-tmp2))
                  adv(mi)%id_mean = tmp1  ! recording mean ice height

                  ! find mean ice/snow mass
                   adv(mi)%mass = 0.
                   do ic=1,ida_n
                   do sc=1,sda_n
                       do ii = 1,ice(sc,ic,mi)%z
                           adv(mi)%mass = adv(mi)%mass + ice(sc,ic,mi)%d(ii)*ice(sc,ic,mi)%th(ii)*ice(sc,ic,mi)%af
                       enddo
                       ! add snow mass
                       do ii = 1,ice(sc,ic,mi)%snow%z
                           adv(mi)%mass = adv(mi)%mass + ice(sc,ic,mi)%snow%d(ii) &
                             *ice(sc,ic,mi)%snow%th(ii)*ice(sc,ic,mi)%af
                       enddo
                  enddo
                  enddo
                  adv(mi)%mass = adv(mi)%mass/1000.d0 ! converted to kg/m^3

                  ! find default convergence area
                  adv(mi)%a_convg = (adv(mi)%a - adv(mi)%a_out + adv(mi)%a_in) - adv(mi)%a_new
                  adv(mi)%f_convg = adv(mi)%a_convg/cell_area


              endif ! end of if (status .ge. 1) at begining of do loop
          enddo

!$OMP END DO


! BEGIN CONVERGENCE OPTIMIZATION
! =================================================================================

! Examine convergence for 8 adjacent squares, reduce by min strength, fix incoming vectors
! and parition conv according to strength vs strength.
!
! Smoothing kernel over ice velocities - what type?
! run new conv. numbers
!
! somehow examine convergence vs. strength and then decide whether to repeat optimization and
! smoothing process

#include "sia2_adv_3_converge_6.inc.f90"
! un-comment 2 lines below if not using any convergence correction
!          do_adv_calc = .false.
!!$OMP END PARALLEL

! =================================================================================
! END CONVERGENCE OPTIMIZATION

          iii = iii + 1  ! increment advection caclution counter

      enddo ! end of "do_adv_calc" while loop

      print *,'Convergence iterations: ',iii

! 3.
! =================================================================================

!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(keep_growing,maxed_depth, &
!$OMP ii,jj,iii,ic,mi,mi1,i_ease,j_ease,int_z,sk_1,sk_z, &
!$OMP int_z_new,z_last,z_snow_new, &
!$OMP out_case, &
!$OMP mi_adv,i_adv,j_adv, &
!$OMP u,v,flag,lost,sklost,smalg_lost,tmp1,tmp2, &
!$OMP tmp3,bv_mean,r_depth,layer_divisor, &
!$OMP th_new,id_new)

!$OMP DO SCHEDULE(DYNAMIC,chunk_size)

          do mi=1,tcells
             if (m(mi)%status .ge. 1) then

                  ! record convergence area
                  m(mi)%a_convg = m(mi)%a_convg + max(0.d0,adv(mi)%a_convg - adv(mi)%a_melt)

                  ! record area dropped during adevection
                  m(mi)%a_drop = m(mi)%a_drop + adv(mi)%a_drop

                  do ic=1,ida_n
                  do sc=1,sda_n

                      ! initialize ice thickness vars
                      adv(mi)%ice(sc,ic)%af = 0.
                      adv(mi)%ice(sc,ic)%th_new = 0.
                      adv(mi)%ice(sc,ic)%id_new = 0.
                      adv(mi)%ice(sc,ic)%z = 0
                      adv(mi)%ice(sc,ic)%depth = 0.
                      adv(mi)%snow(sc,ic)%z = 0
                      adv(mi)%snow(sc,ic)%th_new = 0.  ! vector assigment
                      adv(mi)%snow(sc,ic)%sh_prev = 0.
                      adv(mi)%snow(sc,ic)%depth = 0.


                      if (ice(sc,ic,mi)%af .gt. 0.d0) then
                          ! add local sea ice area
                          adv(mi)%ice(sc,ic)%af = (1.-adv(mi)%out)*ice(sc,ic,mi)%af

                          ! add local ice volume
                          adv(mi)%ice(sc,ic)%depth = &
                            adv(mi)%ice(sc,ic)%af*ice(sc,ic,mi)%id(ice(sc,ic,mi)%z-z_sk)

                          ! add to local skeletal volume
                          adv(mi)%ice(sc,ic)%th_new(sk_1_adv) = &
                            adv(mi)%ice(sc,ic)%af*ice(sc,ic,mi)%th(ice(sc,ic,mi)%z-z_sk+1)

                          ! add local snow volume
                          adv(mi)%snow(sc,ic)%depth = &
                            adv(mi)%ice(sc,ic)%af*ice(sc,ic,mi)%snow%depth
                          adv(mi)%snow(sc,ic)%sh_prev =  &
                            adv(mi)%ice(sc,ic)%af*ice(sc,ic,mi)%sh_prev

                      endif

                      do ii=1,8
                          mi1 = m(mi)%mia(ii)
                          if (adv(mi)%in(ii) .gt. 0.) then

                              ! add advected sea ice area
                              adv(mi)%ice(sc,ic)%af = adv(mi)%ice(sc,ic)%af + &
                                adv(mi)%in(ii)*ice(sc,ic,mi1)%af

                              if (ice(sc,ic,mi1)%af .gt. 0.) then

                                  ! add advected ice volume
                                  adv(mi)%ice(sc,ic)%depth = adv(mi)%ice(sc,ic)%depth + &
                                    adv(mi)%in(ii)*ice(sc,ic,mi1)%af* &
                                    ice(sc,ic,mi1)%id(ice(sc,ic,mi1)%z-z_sk)

                                  ! add advected skeletal volume
                                  adv(mi)%ice(sc,ic)%th_new(sk_1_adv) = &
                                    adv(mi)%ice(sc,ic)%th_new(sk_1_adv) + &
                                    adv(mi)%in(ii)*ice(sc,ic,mi1)%af* &
                                    ice(sc,ic,mi1)%th(ice(sc,ic,mi1)%z-z_sk+1)

                                  ! add advected snow volume
                                  adv(mi)%snow(sc,ic)%depth = &
                                    adv(mi)%snow(sc,ic)%depth + &
                                    adv(mi)%in(ii)*ice(sc,ic,mi1)%af* &
                                    ice(sc,ic,mi1)%snow%depth
                                  adv(mi)%snow(sc,ic)%sh_prev = &
                                    adv(mi)%snow(sc,ic)%sh_prev + &
                                    adv(mi)%in(ii)*ice(sc,ic,mi1)%af* &
                                    ice(sc,ic,mi1)%sh_prev

                              endif
                          endif
                      enddo

                      ! test to see whether this category is involved in advection at all
                      if (adv(mi)%ice(sc,ic)%af .gt. af_min) then

                           ! find advective ice, skeletal, and snow depths
                           adv(mi)%ice(sc,ic)%depth = adv(mi)%ice(sc,ic)%depth/adv(mi)%ice(sc,ic)%af
                           adv(mi)%ice(sc,ic)%th_new(sk_1_adv) = &
                             adv(mi)%ice(sc,ic)%th_new(sk_1_adv)/adv(mi)%ice(sc,ic)%af
                           adv(mi)%snow(sc,ic)%depth = adv(mi)%snow(sc,ic)%depth/adv(mi)%ice(sc,ic)%af
                           adv(mi)%snow(sc,ic)%sh_prev = adv(mi)%snow(sc,ic)%sh_prev/adv(mi)%ice(sc,ic)%af

                           ! find advective ice grid for remaping of tracers
                           r_depth = adv(mi)%ice(sc,ic)%depth

#include "sia2_env_ice_grid.inc.f90"

                           ! we should account for mass/biomass lost if we max out depth here
                           ! when transfering tracers, i think - but if everything coming in
                           ! is less than h_max, then everything resulting here will be
                           ! under h_max too (ridging is later...) - so we should never
                           ! max depth here

                           ! update advective ice grid
                           adv(mi)%ice(sc,ic)%z = int_z_new + z_sk
                           adv(mi)%ice(sc,ic)%th_new(1:int_z_new) = th_new(1:int_z_new)
                           adv(mi)%ice(sc,ic)%id_new(1:int_z_new) = id_new(1:int_z_new)

                           ! transfer skeletal heights to correct grid layer
                           do ii=int_z_new+1,int_z_new+z_sk
                               adv(mi)%ice(sc,ic)%th_new(ii) = adv(mi)%ice(sc,ic)%th_new(sk_1_adv)
                               adv(mi)%ice(sc,ic)%id_new(ii) = adv(mi)%ice(sc,ic)%id_new(ii-1) + &
                                   adv(mi)%ice(sc,ic)%th_new(ii)
                           enddo

                           ! find advective snow grid for remaping of tracers
                           ! but don't bother if below snow_min
!                           if (adv(mi)%snow(sc,ic)%sh_prev .gt. adv(mi)%snow(sc,ic)%depth) then
                               r_depth = adv(mi)%snow(sc,ic)%sh_prev
!                           else
!                               r_depth = adv(mi)%snow(sc,ic)%depth
!                           endif

                           ! only advect snow if there is any to advect in -
                           ! by ding this we avoid the situation where snow
                           ! layers are created during avection. Snow creation
                           ! is handled in the environment routine
                           if (adv(mi)%snow(sc,ic)%depth .gt. snow_min .and. &
                               adv(mi)%snow(sc,ic)%depth .gt. z_th_min) then

#include "sia2_adv_3_snow_grid.inc.f90"

                               ! update advective snow grid
                               adv(mi)%snow(sc,ic)%z = z_snow_new
                               adv(mi)%snow(sc,ic)%th_new(1:z_snow_new) = th_new(1:z_snow_new)

                           else

                               adv(mi)%snow(sc,ic)%z = 0

                           endif

                          do ii=1,ice(sc,ic,mi)%z

                              ! determine model domain & advective ice volume gain/loss
                              ! -----------------------------------------------------------
                              ! record ice volume lost from model domain
                              m(mi)%md_loss = m(mi)%md_loss + ice(sc,ic,mi)%th(ii)*lost* &
                                cell_area*1.e6*ice(sc,ic,mi)%af ! m^3

                              ! record biomass (bm_lost) lost from model domain
                              smalg_lost = ice(sc,ic,mi)%smalg(ii) - min_alg ! mg/m^3
                              smalg_lost = max(0.,smalg_lost)
                              m(mi)%bm_lost = m(mi)%bm_lost + lost*smalg_lost* &
                                1.e3*ice(sc,ic,mi)%th(ii)*cell_area*ice(sc,ic,mi)%af  ! g/pixel

                          enddo


                      else





                      endif ! end of advection-active category test

                  enddo ! end of snow category do loop
                  enddo ! end of ice category do loop
              endif ! end of if (status .ge. 1) at begining of do loop
          enddo

!$OMP END DO
!$OMP END PARALLEL

      end SUBROUTINE sia2_advect_setup

