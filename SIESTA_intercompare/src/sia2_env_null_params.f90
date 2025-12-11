! ==============================================================================
! START NULL PARAMS
! ==============================================================================

		  ice(sc,ic,mi)%af = 0.d0

		  ice(sc,ic,mi)%Ed0 = 0.d0
          ice(sc,ic,mi)%Ed1 = 0.d0
		  ice(sc,ic,mi)%PAR_bot = 0.d0
		  ice(sc,ic,mi)%PAR_bot_pl = 0.d0
          ice(sc,ic,mi)%pur = 0.d0
          ice(sc,ic,mi)%ed_w = 0.d0
          ice(sc,ic,mi)%snow%ed_w = 0.d0

		  ice(sc,ic,mi)%z = 0
		  ice(sc,ic,mi)%no3(z_max+1) = 0.d0
		  ice(sc,ic,mi)%NH4(z_max+1) = 0.d0
		  ice(sc,ic,mi)%PO4(z_max+1) = 0.d0
		  ice(sc,ic,mi)%sioh4(z_max+1) = 0.d0
          ice(sc,ic,mi)%sh_prev = 0.d0
		  ice(sc,ic,mi)%sh_offset = 0.d0


		  ice(sc,ic,mi)%snow%ts = 0.d0
          ice(sc,ic,mi)%age = 0.d0
          ice(sc,ic,mi)%ridged = 0.d0
          ice(sc,ic,mi)%snow_dist = 0.d0
          ice(sc,ic,mi)%snow_rd = 0.d0

		  ice(sc,ic,mi)%snow%depth = 0.d0
		  ice(sc,ic,mi)%snow%z = 0

		  ! zero melt pond vars
		  ice(sc,ic,mi)%pond%t = 0.d0
		  ice(sc,ic,mi)%pond%s = 0.d0
		  ice(sc,ic,mi)%pond%d = 0.d0
		  ice(sc,ic,mi)%pond%heat = 0.d0
		  ice(sc,ic,mi)%pond%th = 0.d0
		  ice(sc,ic,mi)%pond%perc = 0.d0
		  ice(sc,ic,mi)%pond%smalg = 0.d0     ! brine-based algal concentration (mgC/m^3)
		  ice(sc,ic,mi)%pond%poc = 0.d0     ! particulate organic carbon (detritus) (mgC/m^3)
		  ice(sc,ic,mi)%pond%no3 = 0.d0        ! ice layer brine NO3 concentration (µMol)
		  ice(sc,ic,mi)%pond%nh4 = 0.d0      ! ice layer brine NH4 concentration (µMol)
		  ice(sc,ic,mi)%pond%po4 = 0.d0     ! ice layer brine PO4 concentration (µMol)
		  ice(sc,ic,mi)%pond%sioh4 = 0.d0    ! ice layer brine SiOH4 concentration (µMol)

		  do ii=1,z_max

			  ! null 3d params
			  ice(sc,ic,mi)%smalg(ii) = 0.d0
			  ice(sc,ic,mi)%prod(ii) = 0.d0
			  ice(sc,ic,mi)%poc(ii) = 0.d0
			  ice(sc,ic,mi)%Ik1(ii) = 0.d0
			  ice(sc,ic,mi)%gmax(ii) = 0.d0
			  ice(sc,ic,mi)%no3(ii) = 0.d0
			  ice(sc,ic,mi)%nh4(ii) = 0.d0
			  ice(sc,ic,mi)%po4(ii) = 0.d0
			  ice(sc,ic,mi)%sioh4(ii) = 0.d0
			  ice(sc,ic,mi)%no3_mean(ii) = 0.d0
			  ice(sc,ic,mi)%nh4_mean(ii) = 0.d0
			  ice(sc,ic,mi)%po4_mean(ii) = 0.d0
			  ice(sc,ic,mi)%sioh4_mean(ii) = 0.d0
			  ice(sc,ic,mi)%th(ii) = 0.d0
			  ice(sc,ic,mi)%bv(ii) = 0.d0
			  ice(sc,ic,mi)%bd(ii) = 0.d0
			  ice(sc,ic,mi)%bs(ii) = 0.d0
			  ice(sc,ic,mi)%id(ii) = 0.d0
			  ice(sc,ic,mi)%d(ii) = 0.d0
			  ice(sc,ic,mi)%s(ii) = 0.d0
			  ice(sc,ic,mi)%t(ii) = 0.d0
			  ice(sc,ic,mi)%heat(ii) = 0.d0
			  ice(sc,ic,mi)%llim(ii) = 0.d0
			  ice(sc,ic,mi)%nlim(ii) = 0.d0
			  ice(sc,ic,mi)%plim(ii) = 0.d0
			  ice(sc,ic,mi)%silim(ii) = 0.d0
			  ice(sc,ic,mi)%slim(ii) = 0.d0
              ice(sc,ic,mi)%bced(ii) = 0
			  ice(sc,ic,mi)%drained(ii) = 0
			  ice(sc,ic,mi)%dsdt(ii) = 0.d0
			  ice(sc,ic,mi)%fbv(ii) = 0.d0
			  ice(sc,ic,mi)%dhdt_conv(ii) = 0.d0
			  ice(sc,ic,mi)%f0(ii) = 0.d0
			  ice(sc,ic,mi)%dsdt3(ii) = 0.d0
			  ice(sc,ic,mi)%tgrad(ii) = 0.d0
			  ice(sc,ic,mi)%snow%t(ii) = 0.d0
			  ice(sc,ic,mi)%snow%th(ii) = 0.d0
			  ice(sc,ic,mi)%snow%d(ii) = 0.d0
			  ice(sc,ic,mi)%snow%heat(ii) = 0.d0
			  ice(sc,ic,mi)%snow%melt(ii) = 0.d0
!              do iii = 1,lda_n
!				  do kk=1,dt_per_day
!					  PUR(ii,kk,iii,sc,ic,mi) = 0.d0
!				  enddo
!              enddo
		   enddo

!           pur(:,:,:,sc,ic,mi) = pur_0

! ==============================================================================
! END NULL PARAMS
! ==============================================================================
