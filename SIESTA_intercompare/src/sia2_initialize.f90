! file: SIA2_initialize.f90
! Sea Ice Algae Model 2 - Saenz & Arrigo
! Version beta
! ======================================================================
! ======================================================================


! ======================================================================
! Subroutine: sia2_initialize.F90
! Purpose: zeros out memory for 1st step
! ======================================================================

      SUBROUTINE sia2_initialize(pj,m,f,ice,pur,adv,mp_f,ease_f,expt)

          use sia2_globals
          implicit none

          ! function arguments
		  !----------------------------------------------------
          type (proj_type) :: pj
          type (meta_type) :: m(tcells)
          type (forcing_type) :: f(tcells)
          type (ice_type) :: ice(sda_n,ida_n,tcells)
          integer (kind=2) :: pur(z_max,dt_per_day,lda_n,sda_n,ida_n,tcells)
          type (adv_type) :: adv(tcells)
          type (mp_f_type) :: mp_f
          type (ease_f_type) :: ease_f
          double precision :: expt(2,exp_bins)

          ! internal variables
		  !----------------------------------------------------
		  integer i,j,ii,jj,kk,iii,jjj,mi,i_ease,j_ease
		  double precision :: exp_arg2,exp_arg1,exp2,exp1


          ! zero out everything
          ! ---------------------------------------------------
          print *,'Initializing memory to zero...'

          ! zero PUR 24-hours storage array
          pur = pur_0

          ! zero out ease projection forcing vars
          ease_f%icecon = 0.d0
          ease_f%icecon_next = 0.d0
          ease_f%icecon_past = 0.d0
          ease_f%sh = 0.d0
          ease_f%sh_next = 0.d0
          ease_f%ice_vec = 0.d0
          ease_f%icevec_past = 0.d0
          ease_f%sh_past = 0.d0

          do j=1,grid_v
              do i=1,grid_h

                  ! zero out modeled variables
                  if (pj%mi(i,j) .gt. 0) then
                      ii = pj%mi(i,j)

                      ! advection init
                      adv(ii)%af_new = 0.d0

                      ! record modeled cell meta data (grid location, etc.)
                      m(ii)%grid_h = i
                      m(ii)%grid_v = j
                      m(ii)%lat = pj%lat(i,j)
                      m(ii)%lon = pj%lon(i,j)
                      m(ii)%x_mp = pj%x_mp(i,j)
                      m(ii)%y_mp = pj%y_mp(i,j)
                      m(ii)%x_dg = pj%x_dg(i,j)
                      m(ii)%y_dg = pj%y_dg(i,j)
                      m(ii)%x_ec = pj%x_ec(i,j)
                      m(ii)%y_ec = pj%y_ec(i,j)
                      m(ii)%x_eci = pj%x_eci(i,j)
                      m(ii)%y_eci = pj%y_eci(i,j)

                      ! zero single vars
                      m(ii)%status = 0
					  m(ii)%melt_loss = 0.d0
					  m(ii)%adv_loss = 0.d0
					  m(ii)%md_loss = 0.d0
					  m(ii)%cong_growth = 0.d0
					  m(ii)%snow_growth = 0.d0
					  m(ii)%void_growth = 0.d0
					  m(ii)%adv_gain = 0.d0
					  m(ii)%md_gain = 0.d0
					  m(ii)%a_convg = 0.d0
					  m(ii)%a_new = 0.d0
					  m(ii)%a_drop = 0.d0
					  m(ii)%pr_clim = 0.d0
					  m(ii)%pr_ssmi = 0.d0
					  m(ii)%prod_sum_int = 0.d0
					  m(ii)%prod_sum_bot = 0.d0
					  m(ii)%bm_lost = 0.d0
					  m(ii)%tlim = 0.d0
					  m(ii)%nlim = 0.d0
					  m(ii)%silim = 0.d0
					  m(ii)%plim = 0.d0
					  m(ii)%llim = 0.d0
					  m(ii)%slim = 0.d0
                      m(ii)%salt_flux = 0.d0
                      m(ii)%h2o_flux = 0.d0

                      m(ii)%p_wgt_int = 0.d0
                      m(ii)%p_wgt_af_int = 0.d0
                      m(ii)%p_wgt_bot = 0.d0
                      m(ii)%p_wgt_af_bot = 0.d0

                      m(ii)%pwid_sum_int = 0.d0
                      m(ii)%pwid_sum_af_int = 0.d0
                      m(ii)%pwsd_sum_int = 0.d0
                      m(ii)%pwsd_sum_af_int = 0.d0
                      m(ii)%pwt_sum_int = 0.d0
                      m(ii)%pwt_sum_af_int = 0.d0
                      m(ii)%pwid_sum_bot = 0.d0
                      m(ii)%pwid_sum_af_bot = 0.d0
                      m(ii)%pwsd_sum_bot = 0.d0
                      m(ii)%pwsd_sum_af_bot = 0.d0
                      m(ii)%pwt_sum_bot = 0.d0
                      m(ii)%pwt_sum_af_bot = 0.d0

					  f(ii)%sh_interp = 0.d0
					  f(ii)%sh_interp_next = 0.d0
					  f(ii)%sh_interp_last = 0.d0
					  f(ii)%ic_interp = 0.d0
					  f(ii)%ic_interp_next = 0.d0
					  f(ii)%ivu_interp = c9999
					  f(ii)%ivv_interp = c9999
					  f(ii)%ivu_interp_next = c9999
					  f(ii)%ivv_interp_next = c9999
                      f(ii)%at = 0.d0
                      f(ii)%p = 0.d0
                      f(ii)%h = 0.d0
                      f(ii)%u10 = 0.d0
                      f(ii)%v10 = 0.d0
                      f(ii)%fc = 0.d0
                      f(ii)%pr = 0.d0
                      f(ii)%t = 0.d0
                      f(ii)%s = 0.d0
                      f(ii)%d = 0.d0
                      f(ii)%no3 = 0.d0
                      f(ii)%po4 = 0.d0
                      f(ii)%nh4 = 0.d0
                      f(ii)%sioh4 = 0.d0

                      do kk=1,sda_n
                      do jjj=1,ida_n

						  ice(kk,jjj,ii)%Ed0 = 0.d0
						  ice(kk,jjj,ii)%Ed1 = 0.d0
						  ice(kk,jjj,ii)%PAR_bot = 0.d0
						  ice(kk,jjj,ii)%PAR_bot_pl = 0.d0
						  ice(kk,jjj,ii)%af = 0.d0
						  ice(kk,jjj,ii)%pur = 0.d0
						  ice(kk,jjj,ii)%ed_w = 0.d0
						  ice(kk,jjj,ii)%snow%ed_w = 0.d0

                          ice(kk,jjj,ii)%z = 0
						  ice(kk,jjj,ii)%Ts = 0.d0
						  ice(kk,jjj,ii)%sh_prev = 0.d0
						  ice(kk,jjj,ii)%fbh = 0.d0
						  ice(kk,jjj,ii)%sh_offset = 0.d0

						  ice(kk,jjj,ii)%snow%z = 0
						  ice(kk,jjj,ii)%snow%depth = 0.d0
						  ice(kk,jjj,ii)%snow%ts = 0.d0
                          ice(kk,jjj,ii)%age = 0.d0
                          ice(kk,jjj,ii)%ridged = 0.d0
                          ice(kk,jjj,ii)%snow_dist = 0.d0
                          ice(kk,jjj,ii)%snow_rd = 0.d0

                          ! zero melt pond vars
                          ice(kk,jjj,ii)%pond%t = 0.d0
                          ice(kk,jjj,ii)%pond%s = 0.d0
                          ice(kk,jjj,ii)%pond%d = 0.d0
                          ice(kk,jjj,ii)%pond%heat = 0.d0
                          ice(kk,jjj,ii)%pond%th = 0.d0
                          ice(kk,jjj,ii)%pond%perc = 0.d0
                          ice(kk,jjj,ii)%pond%smalg = 0.d0     ! brine-based algal concentration (mgC/m^3)
                          ice(kk,jjj,ii)%pond%poc = 0.d0     ! particulate organic carbon (detritus) (mgC/m^3)
                          ice(kk,jjj,ii)%pond%no3 = 0.d0        ! ice layer brine NO3 concentration (然ol)
                          ice(kk,jjj,ii)%pond%nh4 = 0.d0      ! ice layer brine NH4 concentration (然ol)
                          ice(kk,jjj,ii)%pond%po4 = 0.d0     ! ice layer brine PO4 concentration (然ol)
                          ice(kk,jjj,ii)%pond%sioh4 = 0.d0    ! ice layer brine SiOH4 concentration (然ol)

						  ! zero 3d params
						  do jj=1,z_max
							  ice(kk,jjj,ii)%smalg(jj) = 0.d0
							  ice(kk,jjj,ii)%prod(jj) = 0.d0
							  ice(kk,jjj,ii)%Ik1(jj) = 0.d0
							  ice(kk,jjj,ii)%poc(jj) = 0.d0
							  ice(kk,jjj,ii)%th(jj) = 0.d0
							  ice(kk,jjj,ii)%id(jj) = 0.d0
							  ice(kk,jjj,ii)%t(jj) = 0.d0
							  ice(kk,jjj,ii)%s(jj) = 0.d0
							  ice(kk,jjj,ii)%d(jj) = 0.d0
							  ice(kk,jjj,ii)%bs(jj) = 0.d0
							  ice(kk,jjj,ii)%bd(jj) = 0.d0
							  ice(kk,jjj,ii)%bv(jj) = 0.d0
							  ice(kk,jjj,ii)%heat(jj) = 0.d0
                              ice(kk,jjj,ii)%llim(jj) = 0.d0
                              ice(kk,jjj,ii)%nlim(jj) = 0.d0
                              ice(kk,jjj,ii)%plim(jj) = 0.d0
                              ice(kk,jjj,ii)%silim(jj) = 0.d0
                              ice(kk,jjj,ii)%slim(jj) = 0.d0
                              ice(kk,jjj,ii)%bced(jj) = 0
                              ice(kk,jjj,ii)%drained(jj) = 0
                              ice(kk,jjj,ii)%dsdt(jj) = 0.d0
                              ice(kk,jjj,ii)%fbv(jj) = 0.d0
                              ice(kk,jjj,ii)%dhdt_conv(jj) = 0.d0
                              ice(kk,jjj,ii)%f0(jj) = 0.d0
                              ice(kk,jjj,ii)%dsdt3(jj) = 0.d0
                              ice(kk,jjj,ii)%tgrad(jj) = 0.d0
							  ice(kk,jjj,ii)%snow%t(jj) = 0.d0
							  ice(kk,jjj,ii)%snow%th(jj) = 0.d0
							  ice(kk,jjj,ii)%snow%d(jj) = 0.d0
							  ice(kk,jjj,ii)%snow%heat(jj) = 0.d0
							  ice(kk,jjj,ii)%snow%melt(jj) = 0.d0
						  enddo
						  do jj=1,z_max+1
							  ice(kk,jjj,ii)%no3(jj) = 0.d0
							  ice(kk,jjj,ii)%nh4(jj) = 0.d0
							  ice(kk,jjj,ii)%po4(jj) = 0.d0
							  ice(kk,jjj,ii)%sioh4(jj) = 0.d0
							  ice(kk,jjj,ii)%no3_mean(jj) = 0.d0
							  ice(kk,jjj,ii)%nh4_mean(jj) = 0.d0
							  ice(kk,jjj,ii)%po4_mean(jj) = 0.d0
							  ice(kk,jjj,ii)%sioh4_mean(jj) = 0.d0
						  enddo
                      enddo ! end of sda loop
                      enddo ! end of ida loop
                  endif
              enddo
          enddo

          ! record neighbor cell references
          do mi=1,tcells

			  ! find full EASE domain reference for current cell
			  i_ease = m(mi)%grid_h
			  j_ease = m(mi)%grid_v
			  i = pj%mdh(i_ease,j_ease)
			  j = pj%mdv(i_ease,j_ease)

			  do jj=1,8

				  ! determine indexes to adjacent cells that current cell (ic,mi)
				  ! may be advecting with

				  select case (jj)
					  ! case 1: upper (input)
					  ! OXO
					  ! O O
					  ! OOO
					  case (1)
						  m(mi)%ia(1) = i_ease
						  m(mi)%ja(1) = j_ease-1

					  ! case 2: upper right (input)
					  ! OOX
					  ! O O
					  ! OOO
					  case (2)
						  m(mi)%ia(2) = i_ease+1
						  m(mi)%ja(2) = j_ease-1
					  ! case 3: right (input)
					  ! OOO
					  ! O X
					  ! OOO
					  case (3)
						  m(mi)%ia(3) = i_ease+1
						  m(mi)%ja(3) = j_ease
					  ! case 4: lower right (input)
					  ! OOO
					  ! O O
					  ! OOX
					  case (4)
						  m(mi)%ia(4) = i_ease+1
						  m(mi)%ja(4) = j_ease+1
					  ! case 5: lower (input)
					  ! OOO
					  ! O O
					  ! OXO
					  case (5)
						  m(mi)%ia(5) = i_ease
						  m(mi)%ja(5) = j_ease+1
					  ! case 6: lower right (input)
					  ! OOO
					  ! O O
					  ! XOO
					  case (6)
						  m(mi)%ia(6) = i_ease-1
						  m(mi)%ja(6) = j_ease+1
					  ! case 7: left (input)
					  ! OOO
					  ! X O
					  ! OOO
					  case (7)
						  m(mi)%ia(7) = i_ease-1
						  m(mi)%ja(7) = j_ease
					  ! case 8: upper left (input)
					  ! XOO
					  ! OOO
					  ! OOO
					  case (8)
						  m(mi)%ia(8) = i_ease-1
						  m(mi)%ja(8) = j_ease-1
				  end select

				  if ((m(mi)%ia(jj) .gt. mdh1) .and. (m(mi)%ja(jj) .gt. mdv1) .and. &
				  (m(mi)%ia(jj) .le. mdh2) .and. (m(mi)%ja(jj) .le. mdv2) .and. &
				  (pj%mi(m(mi)%ia(jj),m(mi)%ja(jj)) .ne. 0)) then
					  ! record neighbor cell
					  m(mi)%mia(jj) = pj%mi(m(mi)%ia(jj),m(mi)%ja(jj))
                  else
                      ! there is no neighbor cell
                      m(mi)%mia(jj) = 0
					  m(mi)%ia(jj) = 0
					  m(mi)%ja(jj) = 0
                  endif
              enddo  ! end of 8 neighbors
          enddo  ! end of tcells loop

          ! zero out mercator projection forcing vars - just in case
          do j=1,mp_y
              do i=1,mp_x
                  mp_f%ncep%at(i,j) = 0.d0
                  mp_f%ncep%p(i,j) = 0.d0
                  mp_f%ncep%h(i,j) = 0.d0
                  mp_f%ncep%fc(i,j) = 0.d0
                  mp_f%ncep%u10(i,j) = 0.d0
                  mp_f%ncep%v10(i,j) = 0.d0
                  mp_f%ncep%pr(i,j) = 0.d0
                  mp_f%ncep_next%at(i,j) = 0.d0
                  mp_f%ncep_next%p(i,j) = 0.d0
                  mp_f%ncep_next%h(i,j) = 0.d0
                  mp_f%ncep_next%fc(i,j) = 0.d0
                  mp_f%ncep_next%u10(i,j) = 0.d0
                  mp_f%ncep_next%v10(i,j) = 0.d0
                  mp_f%ncep_next%pr(i,j) = 0.d0
                  mp_f%ncep_interp%at(i,j) = 0.d0
                  mp_f%ncep_interp%p(i,j) = 0.d0
                  mp_f%ncep_interp%h(i,j) = 0.d0
                  mp_f%ncep_interp%fc(i,j) = 0.d0
                  mp_f%ncep_interp%u10(i,j) = 0.d0
                  mp_f%ncep_interp%v10(i,j) = 0.d0
                  mp_f%ncep_interp%pr(i,j) = 0.d0
                  do jj = 1,wavl
                      mp_f%Edir(1)%Ed(i,j,jj) = 0.d0
                      mp_f%Edir(2)%Ed(i,j,jj) = 0.d0
                      mp_f%Edif(1)%Ed(i,j,jj) = 0.d0
                      mp_f%Edif(2)%Ed(i,j,jj) = 0.d0
                  enddo
              enddo
          enddo

          print *,'Pre-calculating exp lookup table'

          do i=1,exp_bins
              if (i .eq. 1) then
                  exp_arg1 = exp_a_max
                  exp1 = exp_max
              else
                  exp_arg1 = expls*dble(i-1)
                  exp1 = exp(exp_arg1)
              endif
              exp_arg2 = expls*dble(i)
              exp2 = exp(exp_arg2)
              expt(1,i) = (exp2-exp1)/(exp_arg2-exp_arg1)  ! slope
              expt(2,i) = exp2 - exp_arg2*expt(1,i)        ! intercept
          enddo



          print *,'Done Initializing memory.'

      end SUBROUTINE sia2_initialize

