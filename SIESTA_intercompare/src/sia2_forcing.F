! file: SIA2_forcing.f90
! Sea Ice Algae Model 2 - Saenz & Arrigo
! Version beta
! ======================================================================
! ======================================================================

!=======================================================================
! Subroutine: sia2_forcing_load
! Purpose: is called at every timestep, and loads appropriate forcing
! variables given timing
!=======================================================================

      SUBROUTINE sia2_forcing_load(pj,mp_f,dg_f,ec_f,eci_f,ease_f,v_f)

        use sia2_globals
        implicit none

        ! function arguments
    !----------------------------------------------------
        type (proj_type) :: pj
        type (mp_f_type) :: mp_f
        type (dg_f_type) :: dg_f
        type (ec_f_type) :: ec_f
        type (eci_f_type) :: eci_f
        type (ease_f_type) :: ease_f
        type (boundary_f_type) :: v_f

        ! internal variables
    !----------------------------------------------------
        logical exists,tmp_leap
        integer sizeof,i,j,kk,cur_date80,offset,index,j_mdh,i_mdv,mi
        integer ncep_next_yr,today,next_yr,tomorrow,month,next_month, &
        	this_5_day,next_5_day,next_5_year,icevec_dim
        !integer*2 temp_row(1441)
        !integer*2 temp_row_sm(721)
        !integer*2 temp_row_vec(963)
        !character(len=77) :: iceconfilename
        !character(len=79) :: icetempfilename
        !character(len=83) :: snowdepthfilename
        !character(len=64) :: iceconfilename_ssmi
        !character(len=60) :: icevecfilename

        double precision :: last_icecon,ncep_next_hr,smalg_lost,bv_mean,next_hr


        !---------------------------------------------------------------------
        ! calculate timing variables
        !---------------------------------------------------------------------

        call sia2_days_1980(cur_hour,cur_year,cur_date80)
        today = int(cur_hour/24.)+1
        call sia2_forcing_next_day(cur_year,cur_hour,next_yr,next_hr,leap_year)
        tomorrow = int(next_hr/24.)+1
        call sia2_getmonth(cur_hour,leap_year,month)
        next_month = month+1
        if (next_month .gt. 12) then
            next_month = 1
        endif
        this_5_day = int(cur_hour/120.)+1
        next_5_day = this_5_day + 1
        next_5_year = cur_year
        if (next_5_day .gt. 73) then
        	next_5_day = 1
        	next_5_year = cur_year + 1
        endif

        !---------------------------------------------------------------------
        ! Load hourly physics
        !---------------------------------------------------------------------
        if (int(cur_hour)+1 .gt. last_hour) then

            if (start .or. restart .eq. 1) then

        print *, cur_year,int(cur_hour)+1,'loading Edif...'
        !  Open Gregg & Carter diffuse surface spectral PAR, cloud corrected
        !-----------------------------------------------------------------
        call gc_read(cur_year, int(cur_hour)+1, mp_f%Edif(f_index), 2)

        print *, cur_year,int(cur_hour)+1,'loading Edir...'
        !  Open Gregg & Carter direct surface spectral PAR, cloud corrected
        !-----------------------------------------------------------------
        call gc_read(cur_year, int(cur_hour)+1, mp_f%Edir(f_index), 1)

            endif

            ! Load hourly interpolated station-specific forcing data
            if (boundary_file .eq. 1) then
               print *, cur_year,int(cur_hour)+1,'loading boundary data ...'
               call sia2_forcing_load_v_f(v_f,cur_year,int(cur_hour)+1)
            endif

        endif


        if (atm_f .eq. 2) then
          if (int(cur_hour/3.)+1 .gt. last3hour) then

            if (start .or. restart .eq. 1 .or. (int(cur_hour/12.)+1 .gt. last12hour)) then
              eci_f%ecmwf%pr = 0.  ! zero previous accumulation if we have past the 12hour accumulation limit
            else
              eci_f%ecmwf%pr = eci_f%ecmwf_next%pr   ! save old forcing for accumulation subtraction
            endif

            print *, cur_year,int(cur_hour/3.)+1,'loading ECMWF Interim 3-hour climatology...'
            ! load following forcing for interpolation
            call ecmwf3_int_read(cur_year, int(cur_hour/3.)+1, eci_f%ecmwf_next)

          endif
        endif

        !---------------------------------------------------------------------
        ! Load 4xdaily physics (every 6 hours), if using them
        !---------------------------------------------------------------------
        if (int(cur_hour/6.)+1 .gt. last6hour) then

          ! find indexes for NEXT timestep, to read the next forcing data
          ! used for interpolation between between 24 hour forcing...
          ! --------------------------------------------------------
          call sia2_forcing_next_ncep(cur_year,cur_hour,6.d0,ncep_next_yr, &
            ncep_next_hr,leap_year)
          !print *,'NCEP forcing: ',cur_year,cur_hour,int(cur_hour/6)+1,ncep_next_yr,ncep_next_hr

          if (atm_f .eq. 0) then

            print *, cur_year,int(cur_hour/6)+1,'loading NCEP DOE Re. II 6-hour climatology...'
            ! NCEP Surface Air Temp
            ! --------------------------------------------------------
            ! load current airtemp for 1st step...
            if (start .or. restart .eq. 1) then
              call ncep6_read(cur_year, int(cur_hour/6)+1, mp_f%ncep_next%at, ncep_at)   ! static temp for 1st 6 hours...
            endif
            !print *, cur_year,int(cur_hour/6)+1,'loading air temp climatology...'
            ! save old airtemp for interpolation
            mp_f%ncep%at = mp_f%ncep_next%at
            ! load following airtemp for interpolation
            call ncep6_read(ncep_next_yr,int(ncep_next_hr/6)+1, mp_f%ncep_next%at, ncep_at)  ! what happens at year change?

            ! NCEP 2m Air Pressure
            ! --------------------------------------------------------
            ! load current air pressure for 1st step...
            if (start .or. restart .eq. 1) then
              call ncep6_read(cur_year, int(cur_hour/6)+1, mp_f%ncep_next%p, ncep_p)   ! static temp for 1st 6 hours...
            endif
            !print *, cur_year,int(cur_hour/6)+1,'loading air pressure climatology...'
            ! save old airtemp for interpolation
            mp_f%ncep%p = mp_f%ncep_next%p
            ! load following airtemp for interpolation
            call ncep6_read(ncep_next_yr,int(ncep_next_hr/6)+1, mp_f%ncep_next%p, ncep_p)  ! what happens at year change?

            ! NCEP 2m Air Pressure
            ! --------------------------------------------------------
            ! load current humidity for 1st step...
            if (start .or. restart .eq. 1) then
              call ncep6_read(cur_year, int(cur_hour/6)+1, mp_f%ncep_next%h, ncep_h)   ! static temp for 1st 6 hours...
            endif
            !print *, cur_year,int(cur_hour/6)+1,'loading humidity climatology...'
            ! save old airtemp for interpolation
            mp_f%ncep%h = mp_f%ncep_next%h
            ! load following airtemp for interpolation
            call ncep6_read(ncep_next_yr,int(ncep_next_hr/6)+1, mp_f%ncep_next%h, ncep_h)  ! what happens at year change?

            ! NCEP surface specific humidity
            ! --------------------------------------------------------
            ! load current fractional cloud cover for 1st step...
            if (start .or. restart .eq. 1) then
              call ncep6_read(cur_year, int(cur_hour/6)+1, mp_f%ncep_next%fc, ncep_fc)   ! static temp for 1st 6 hours...
            endif
            !print *, cur_year,int(cur_hour/6)+1,'loading cloud cover climatology...'
            ! save old airtemp for interpolation
            mp_f%ncep%fc = mp_f%ncep_next%fc
            ! load following airtemp for interpolation
            call ncep6_read(ncep_next_yr,int(ncep_next_hr/6)+1, mp_f%ncep_next%fc, ncep_fc)  ! what happens at year change?

            ! NCEP 10m u-direction wind speed
            ! --------------------------------------------------------
            ! load current V10-u wind speed for 1st step...
            if (start .or. restart .eq. 1) then
              call ncep6_read(cur_year, int(cur_hour/6)+1, mp_f%ncep_next%u10, ncep_u10)   ! static temp for 1st 6 hours...
            endif
            !print *, cur_year,int(cur_hour/6)+1,'loading V10-u wind speed climatology...'
            ! save old airtemp for interpolation
            mp_f%ncep%u10 = mp_f%ncep_next%u10
            ! load following airtemp for interpolation
            call ncep6_read(ncep_next_yr,int(ncep_next_hr/6)+1, mp_f%ncep_next%u10, ncep_u10)  ! what happens at year change?

            ! NCEP 10m v-direction wind speed
            ! --------------------------------------------------------
            ! load current V10-v wind speed for 1st step...
            if (start .or. restart .eq. 1) then
              call ncep6_read(cur_year, int(cur_hour/6)+1, mp_f%ncep_next%v10, ncep_v10)   ! static temp for 1st 6 hours...
            endif
            !print *, cur_year,int(cur_hour/6)+1,'loading V10-v wind speed climatology...'
            ! save old airtemp for interpolation
            mp_f%ncep%v10 = mp_f%ncep_next%v10
            ! load following airtemp for interpolation
            call ncep6_read(ncep_next_yr,int(ncep_next_hr/6)+1, mp_f%ncep_next%v10, ncep_v10)  ! what happens at year change?

            ! NCEP precipitation rate
            ! --------------------------------------------------------
            ! load current precip rate for 1st step...
            if (start .or. restart .eq. 1) then
              call ncep6_read(cur_year, int(cur_hour/6)+1, mp_f%ncep_next%pr, ncep_pr)   ! static temp for 1st 6 hours...
            endif
            !print *, cur_year,int(cur_hour/6)+1,'loading precip. rate climatology...'
            ! save old airtemp for interpolation
            mp_f%ncep%pr = mp_f%ncep_next%pr
            ! load following airtemp for interpolation
            call ncep6_read(ncep_next_yr,int(ncep_next_hr/6)+1, mp_f%ncep_next%pr, ncep_pr)  ! what happens at year change?

          elseif (atm_f .eq. 1) then

            ! Load ECMWF
            ! --------------------------------------------------------
            ! load ECMWF for 1st step...
            if (start .or. restart .eq. 1) then
              call ecmwf6_read(cur_year, int(cur_hour/6)+1, ec_f%ecmwf_next)   ! load precip rate in m for next 6 hours
            endif
            !print *, cur_year,int(cur_hour/6)+1,'loading ECMWF climatology...'
            ! save old airtemp for interpolation
            ec_f%ecmwf = ec_f%ecmwf_next
            ! load following airtemp for interpolation
            call ecmwf6_read(ncep_next_yr,int(ncep_next_hr/6)+1, ec_f%ecmwf_next)

          elseif (atm_f .eq. 2) then

            ! Load ECMWF Interim 6 hourly
            ! --------------------------------------------------------
            ! load ECMWF Interim for 1st step...
            if (start .or. restart .eq. 1) then
              call ecmwf6_int_read(cur_year, int(cur_hour/6)+1, eci_f%ecmwf_next)   ! load precip rate in m for next 6 hours
            endif
            !print *, cur_year,int(cur_hour/6)+1,'loading ECMWF Interim 6-hour climatology...'

            ! save old forcing for interpolation
            eci_f%ecmwf%at = eci_f%ecmwf_next%at          ! 2m air temperature (kelvin)
            eci_f%ecmwf%p = eci_f%ecmwf_next%p            ! surface air pressure (Pa)
            eci_f%ecmwf%dpt = eci_f%ecmwf_next%dpt        ! dew point temperature (kelvin)
            eci_f%ecmwf%fc = eci_f%ecmwf_next%fc          ! fraction cloud cover (%)
            eci_f%ecmwf%u10 = eci_f%ecmwf_next%u10        ! 10m horizontal wind speed vector (m/s)
            eci_f%ecmwf%v10 = eci_f%ecmwf_next%v10            ! vertical wind speed vector (m/s)

            ! load following forcing for interpolation
            call ecmwf6_int_read(ncep_next_yr,int(ncep_next_hr/6)+1, eci_f%ecmwf_next)

          endif  ! end of ECMWF forcing load

        endif


          !---------------------------------------------------------------------
          ! Load 2xdaily physics (every 12 hours), if using them
          !---------------------------------------------------------------------
!          if (int(cur_hour/12.)+1 .gt. last12hour .and. atm_f .eq. 2) then

			  ! find indexes for NEXT timestep, to read the next forcing data
			  ! used for interpolation between between 24 hour forcing...
			  ! --------------------------------------------------------
!			  call sia2_forcing_next_ncep(cur_year,cur_hour,12.d0,ncep_next_yr, &
!				  ncep_next_hr,leap_year)

			  ! Load ECMWF Interim 12 hourly
			  ! --------------------------------------------------------
			  ! load ECMWF Interim for 1st step...
!			  if (start .or. restart .eq. 1) then
!				  call ecmwf12_int_read(cur_year, int(cur_hour/12.)+1, eci_f%ecmwf_next)   ! load precip rate in m for next 6 hours
!			  endif
!			  print *, cur_year,int(cur_hour/12)+1,'loading ECMWF Interim 12-hour climatology...'

			  ! save old forcing for interpolation
!			  eci_f%ecmwf%pr = eci_f%ecmwf_next%pr              ! total precipitation (m)
!             eci_f%ecmwf%ssr = eci_f%ecmwf_next%ssr            ! surface solar radiation (W/m^2)

			  ! load following forcing for interpolation
!			  call ecmwf12_int_read(ncep_next_yr,int(ncep_next_hr/12.)+1, eci_f%ecmwf_next)

!         endif


        !---------------------------------------------------------------------
        ! Load daily physics
        !---------------------------------------------------------------------

        ! pre-load snow/icecon 1 timestep before day change, so that interpolation is
        if (start .or. restart .eq. 1) then

            ! pre-load and init snow
          call sia2_forcing_ssmisnow(ease_f)

! ----------------------------------------------------------------------
! with smoothing -------------------------------------------------------
!
!          ! pre-load and init ice vectors
!  	       call sia2_forcing_icevec(ease_f)
!
!          ! pre-load and init icecon
!  	       call sia2_forcing_icecon(ease_f)
!
! with smoothing -------------------------------------------------------
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! without smoothing -----------------------------------------------------
          call SSMI_icecon_read(cur_year, today)   ! reading into global SSMI_grid_int2

          do j = 1,721
            do i=1,721
              ease_f%icecon_next(i,j) = dfloat(SSMI_grid_int2(i,j))/100.  ! make icecon a fraction
              if(ease_f%icecon_next(i,j) .ge. 0.2d0 .and. &
              ease_f%icecon_next(i,j) .le. 1.d0) then
                if(icecon_mod .eq. 1) then
                  ease_f%icecon_next(i,j) = &
                    (ease_f%icecon_next(i,j)-0.2d0)*8.d0/7.d0 + 0.2d0
                  ease_f%icecon_next(i,j) = min(ease_f%icecon_next(i,j),1.d0)
                elseif(icecon_mod .eq. 2) then
                  ease_f%icecon_next(i,j) = &
                    (ease_f%icecon_next(i,j)-0.2d0)*8.d0/6.d0 + 0.2d0
                  ease_f%icecon_next(i,j) = min(ease_f%icecon_next(i,j),1.d0)
                endif
              endif
            enddo
          enddo

          call sia2_forcing_next_day(cur_year,cur_hour,next_yr,next_hr,leap_year)
          call SSMI_icecon_read(next_yr, tomorrow)   ! reading into global SSMI_grid_int2
          do j = 1,721
            do i=1,721
              ease_f%icecon(i,j) = ease_f%icecon_next(i,j)
              ease_f%icecon_next(i,j) = dfloat(SSMI_grid_int2(i,j))/100.  ! make icecon a fraction
              if(ease_f%icecon_next(i,j) .ge. 0.2d0 .and. &
              ease_f%icecon_next(i,j) .le. 1.d0) then
                if(icecon_mod .eq. 1) then
                  ease_f%icecon_next(i,j) = &
                    (ease_f%icecon_next(i,j)-0.2d0)*8.d0/7.d0 + 0.2d0
                  ease_f%icecon_next(i,j) = min(ease_f%icecon_next(i,j),1.d0)
                elseif(icecon_mod .eq. 2) then
                  ease_f%icecon_next(i,j) = &
                    (ease_f%icecon_next(i,j)-0.2d0)*8.d0/6.d0 + 0.2d0
                  ease_f%icecon_next(i,j) = min(ease_f%icecon_next(i,j),1.d0)
                endif
              endif
            enddo
          enddo

          ! pre-load and ice vectors
          call SSMI_icevec_read(cur_year, today)   ! reading into global icevec_grid_int2
          if (hemi .eq. 0) then
              icevec_dim = 321
              offset = 200 ! 25km ice vectors files are only 321x321
          else
              icevec_dim = 361
              offset = 180
          endif
          do kk=1,3
              do j=1,icevec_dim
                  do i=1,icevec_dim
                      ease_f%ice_vec(offset+i,offset+j,kk) = &
                        dfloat(icevec_grid_int2(i,j,kk))
                  enddo
              enddo
          enddo
          call SSMI_icevec_read(next_yr, tomorrow)   ! reading into global icevec_grid_int2
          do kk=1,3
              do j=1,icevec_dim
                  do i=1,icevec_dim
                      ease_f%ice_vec_next(offset+i,offset+j,kk) = &
                        dfloat(icevec_grid_int2(i,j,kk))
                  enddo
              enddo
          enddo
! without smoothing -----------------------------------------------------
! ----------------------------------------------------------------------


        elseif (int((cur_hour+dt)/24)+1 .gt. today) then

          !  Open 25km snow depth
          !---------------------------------------------------------
          print *, cur_year,tomorrow,'loading SMMI snow depths from netcdf...'
          call sia2_forcing_ssmisnow(ease_f)

! ----------------------------------------------------------------------
! with smoothing -------------------------------------------------------
!
!          ! open/load ice vectors
!	       call sia2_forcing_icevec(ease_f)
!
!          ! open/load ice concentrations
!	       call sia2_forcing_icecon(ease_f)
!
! with smoothing -------------------------------------------------------
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! without smoothing -----------------------------------------------------

          !  Open 25km ice concentration
          !-----------------------------------------------------------------
          ! load current icecon from icecon_next, and read in icecon_next
          print *, next_yr, tomorrow,'loading next SMMI icecon from netcdf...'
          call SSMI_icecon_read(next_yr, tomorrow)   ! reading into global SSMI_grid_int2
          ease_f%icecon = ease_f%icecon_next
          ease_f%icecon_next = dble(SSMI_grid_int2)/100.
          do j=1,721
            do i=1,721
              if(ease_f%icecon_next(i,j) .ge. 0.2d0 .and. &
              ease_f%icecon_next(i,j) .le. 1.d0) then
                if(icecon_mod .eq. 1) then
                  ease_f%icecon_next(i,j) = &
                    (ease_f%icecon_next(i,j)-0.2d0)*8.d0/7.d0 + 0.2d0
                  ease_f%icecon_next(i,j) = min(ease_f%icecon_next(i,j),1.d0)
                elseif(icecon_mod .eq. 2) then
                  ease_f%icecon_next(i,j) = &
                    (ease_f%icecon_next(i,j)-0.2d0)*8.d0/6.d0 + 0.2d0
                  ease_f%icecon_next(i,j) = min(ease_f%icecon_next(i,j),1.d0)
                endif
              endif
            enddo
          enddo

                !  Open 25km ice motion vectors
                !-----------------------------------------------------------------
          print *, next_yr,tomorrow,'loading current SMMI ice vectors from netcdf...'
          call SSMI_icevec_read(next_yr, tomorrow)   ! reading into global icevec_grid_int2
          ease_f%ice_vec = ease_f%ice_vec_next
          offset = 200 ! 25km ice vectors files are only 321x321
          do kk=1,3
              do j=1,321
                  do i=1,321
                      ease_f%ice_vec_next(offset+i,offset+j,kk) = &
                      dfloat(icevec_grid_int2(i,j,kk))
                  enddo
              enddo
          enddo

! without smoothing -----------------------------------------------------
! ----------------------------------------------------------------------


        endif

        ! ------------------------------------------------------------
        ! 5 DAY PHYSICS
        ! ------------------------------------------------------------
        if (this_5_day .ne. last_5_day) then

					if (use_sose .eq. 1 .and. cur_year >= 2005) then

            print *, cur_year,month,'loading SOSE fields...'

						if (start .or. restart .eq. 1 .or. &
						(cur_year == 2005 .and. this_5_day < 5)) then
  						call sia2_sose_read(this_5_day,cur_year,sose_level,sose_f%sose_next)
  					endif
  					sose_f%sose = sose_f%sose_next
  					call sia2_sose_read(next_5_day,next_5_year,sose_level,sose_f%sose_next)

					endif

				endif

        ! ------------------------------------------------------------
        ! MONTHLY PHYSICS
        ! ------------------------------------------------------------

        if (month .gt. dble(last_month)) then

            print *, cur_year,month,'loading netcdf World Ocean Atlas climatology...'

            ! load water temp climatology
            if (start .or. restart .eq. 1) then
                ! load current stuff
                call sia2_woa_read(month,woa_depth,dg_f%woa_next%t,woa_t)
            endif
            dg_f%woa%t = dg_f%woa_next%t
            call sia2_woa_read(next_month,woa_depth,dg_f%woa_next%t,woa_t)

            ! load salinity climatology
            if (start .or. restart .eq. 1) then
                ! load current stuff
                call sia2_woa_read(month,woa_depth,dg_f%woa_next%s,woa_s)
            endif
            dg_f%woa%s = dg_f%woa_next%s
            call sia2_woa_read(next_month,woa_depth,dg_f%woa_next%s,woa_s)

            ! load nitrogen climatology
            if (start .or. restart .eq. 1) then
                ! load current stuff
                call sia2_woa_read(month,woa_depth,dg_f%woa_next%n,woa_n)
            endif
            dg_f%woa%n = dg_f%woa_next%n
            call sia2_woa_read(next_month,woa_depth,dg_f%woa_next%n,woa_n)

            ! load phosphorus climatology
            if (start .or. restart .eq. 1) then
                ! load current stuff
                call sia2_woa_read(month,woa_depth,dg_f%woa_next%p,woa_p)
            endif
            dg_f%woa%p = dg_f%woa_next%p
            call sia2_woa_read(next_month,woa_depth,dg_f%woa_next%p,woa_p)

            ! load silica climatology
            if (start .or. restart .eq. 1) then
                ! load current stuff
                call sia2_woa_read(month,woa_depth,dg_f%woa_next%si,woa_si)
            endif
            dg_f%woa%si = dg_f%woa_next%si
            call sia2_woa_read(next_month,woa_depth,dg_f%woa_next%si,woa_si)

        endif


      end SUBROUTINE sia2_forcing_load


!=======================================================================
! Subroutine: sia2_forcing_ssmisnow
! Purpose: loads snow depths over ice, keeping seven
! days in memory, for use later inconstructing a 7-day moving average
! (six days behind + current day)
!=======================================================================
     SUBROUTINE sia2_forcing_ssmisnow(ease_f)

		  use sia2_globals
		  implicit none

          ! function arguments
          ! ------------------------------------------------------------
          type (ease_f_type) :: ease_f

          ! local variables
          ! ------------------------------------------------------------
          integer :: next_yr
          double precision :: next_hr

          ! pre-load past snowdepths for use in first time step
          ! ------------------------------------------------------------
		  if (start .or. restart .eq. 1) then
              ! load 7 days worth of back snow depths
			  call sia2_forcing_ssmisnow1(ease_f,1,cur_year,cur_hour-144.)
			  call sia2_forcing_ssmisnow1(ease_f,2,cur_year,cur_hour-120.)
			  call sia2_forcing_ssmisnow1(ease_f,3,cur_year,cur_hour-96.)
			  call sia2_forcing_ssmisnow1(ease_f,4,cur_year,cur_hour-72.)
			  call sia2_forcing_ssmisnow1(ease_f,5,cur_year,cur_hour-48.)
			  call sia2_forcing_ssmisnow1(ease_f,6,cur_year,cur_hour-24.)
			  call sia2_forcing_ssmisnow1(ease_f,7,cur_year,cur_hour)
              ! create first mean snowdepth
              call sia2_forcing_update_snow(.true.,ease_f)
              ! find next day's snow & load
              call sia2_forcing_next_day(cur_year,cur_hour,next_yr,next_hr,leap_year)
    		  call sia2_forcing_ssmisnow1(ease_f,1,next_yr,next_hr)
              ! find mean snow for next day
              call sia2_forcing_update_snow(.false.,ease_f)
              snowd_index = 1
		  else

              ! find yr/hr which will load next timestep's snow
              call sia2_forcing_next_day(cur_year,cur_hour+dt,next_yr,next_hr,leap_year)

              ! increment index of snowpast array
			  if (snowd_index .eq. 7) then
				  snowd_index = 1
			  else
				  snowd_index = snowd_index+1
			  endif

			  ! load next timestep's snow
			  call sia2_forcing_ssmisnow1(ease_f,snowd_index,next_yr,next_hr)

			  ! interpolate over 7 days, and move snowh_next -> snowh
			  call sia2_forcing_update_snow(.false.,ease_f)
          endif

     end SUBROUTINE sia2_forcing_ssmisnow



!=======================================================================
! Subroutine: sia2_forcing_ssmisnow1
! Purpose: loads one day of snow depths
!=======================================================================
      SUBROUTINE sia2_forcing_ssmisnow1(ease_f,sindex,year,hour)
              use sia2_globals
              implicit none

          ! function arguments
          ! ------------------------------------------------------------
          integer :: year,sindex
          type (ease_f_type) :: ease_f
          double precision :: hour

		  !  Open 25km snow depth
		  !----------------------
		  !  Snow depth in cm &
		  !  251=bad brightness temp
		  !  252=no sea ice
		  !  values above 60cm should be capped, and values < zero set to zero
		  !  use a 5-day average, mark all pixels with derivative of
		  !  greater than +or-5 cm/day as bad!  Interpolate?
		  !-----------------------------------------------------------------
		  call SSMI_snow_read(year, int(hour/24)+1)   ! reading into global SSMI_grid_int2
		  ease_f%sh_past(:,:,sindex) = dfloat(SSMI_grid_int2)

      end SUBROUTINE sia2_forcing_ssmisnow1

!=======================================================================
! Subroutine: sia2_forcing_icecon1
! Purpose: loads one day of sea ice concentrations
!=======================================================================
      SUBROUTINE sia2_forcing_icecon1(ease_f,icindex,year,hour)
              use sia2_globals
              implicit none

          ! function arguments
          ! ------------------------------------------------------------
          integer :: year,icindex,i,j
          type (ease_f_type) :: ease_f
          double precision :: hour

        !  Open 25km snow ice concentrations
        !-----------------------------------------------------------------
        call SSMI_icecon_read(year, int(hour/24)+1)   ! reading into global SSMI_grid_int2
        ease_f%icecon_past(:,:,icindex) = dfloat(SSMI_grid_int2)/100.d0

        ! implement icecon mod if used
          do j=1,721
            do i=1,721
              if(ease_f%icecon_past(i,j,icindex) .ge. 0.2d0 .and. &
              ease_f%icecon_past(i,j,icindex) .le. 1.d0) then
                if(icecon_mod .eq. 1) then
                  ease_f%icecon_past(i,j,icindex) = &
                    (ease_f%icecon_past(i,j,icindex)-0.2d0)*8.d0/7.d0 + 0.2d0
                  ease_f%icecon_past(i,j,icindex) = min(ease_f%icecon_past(i,j,icindex),1.d0)
                elseif(icecon_mod .eq. 2) then
                  ease_f%icecon_past(i,j,icindex) = &
                    (ease_f%icecon_past(i,j,icindex)-0.2d0)*8.d0/6.d0 + 0.2d0
                  ease_f%icecon_past(i,j,icindex) = min(ease_f%icecon_past(i,j,icindex),1.d0)
                endif
              endif
            enddo
          enddo

      end SUBROUTINE sia2_forcing_icecon1


!=======================================================================
! Subroutine: sia2_forcing_update_snow
! Purpose: at each time step to estimate/interpolate snow depth based on a 7 day
! mean of smmi calculated snow depth
!=======================================================================
      SUBROUTINE sia2_forcing_update_snow(init,ease_f)
          use sia2_globals
          implicit none


          ! function arguments
          ! ------------------------------------------------------------
          logical :: init
          type (ease_f_type) :: ease_f

          ! local variables
          ! ------------------------------------------------------------
          double precision :: sum,new_snowh,variance
          integer :: i,j,ii,count

          print *,'Interpolating snowdepth over last 7 days...'

          do j= 1,721
              do i = 1,721

				  sum = 0.
				  count = 0
				  variance = 0.
				  do ii = 1,7
                      if (snow_model .lt. 2) then

						  if ((ease_f%sh_past(i,j,ii) .lt. 251.) .and. &
						  (ease_f%sh_past(i,j,ii) .ge. 0.)) then
							  sum = sum + ease_f%sh_past(i,j,ii)
							  count = count+1;
						  endif
                      else
						  if ((ease_f%sh_past(i,j,ii) .lt. 60.) .and. &
						  (ease_f%sh_past(i,j,ii) .ge. 0.)) then
							  sum = sum + ease_f%sh_past(i,j,ii)
							  count = count+1;
						  endif
					      if ((ease_f%sh_past(i,j,ii) .lt. 251.) .and. &
						  (ease_f%sh_past(i,j,ii) .ge. 0.) .and. (ii .gt. 1)) then
							  variance = variance + abs(ease_f%sh_past(i,j,ii) - &
                                  ease_f%sh_past(i,j,ii-1))
                          endif

                      endif
				  enddo
                  if (snow_model .lt. 2) then
					  if (count .gt. 0) then
						  new_snowh = sum/dfloat(count)
						  if (init) then
							  ! intial snow assignment
							  ease_f%sh_next(i,j) = new_snowh
						  endif
						  ease_f%sh(i,j) = ease_f%sh_next(i,j)
						  ease_f%sh_next(i,j) = new_snowh
					  endif
                  else
					  if (count .gt. 0) then
					      variance = variance/dble(count)
                          if (variance .lt. 20.) then
							  new_snowh = sum/dfloat(count)
							  if (init) then
								  ! intial snow assignment to prepare for future steps
								  ease_f%sh_next(i,j) = new_snowh
							  endif
							  ease_f%sh(i,j) = ease_f%sh_next(i,j)
							  ease_f%sh_next(i,j) = new_snowh
                          else
                              ease_f%sh_next(i,j) = 9999.
                          endif
                      else
                          ease_f%sh_next(i,j) = 9999.
                      endif
                  endif

              enddo
          enddo

      end SUBROUTINE sia2_forcing_update_snow


!=======================================================================
! Subroutine: sia2_forcing_update_icecon
! Purpose: at each time step to estimate/interpolate ice concenctration
! based on a 5 day mean of ice concectrations
!=======================================================================
      SUBROUTINE sia2_forcing_update_icecon(init,ease_f)
        use sia2_globals
        implicit none


        ! function arguments
        ! ------------------------------------------------------------
        logical :: init
        type (ease_f_type) :: ease_f

        ! local variables
        ! ------------------------------------------------------------
        double precision :: sum,new_snowh,variance
        integer :: i,j,ii,count

        print *,'Smoothing ice concentratation over last 5 days...'

        do j= 1,721
          do i = 1,721

            sum = 0.
            count = 0
            variance = 0.
!            do ii = 1,5
!            do ii = 1,3
            do ii = 1,2

              if ((ease_f%icecon_past(i,j,ii) .le. 1.) .and. &
              (ease_f%icecon_past(i,j,ii) .ge. 0.)) then
                sum = sum + ease_f%icecon_past(i,j,ii)
                count = count+1;
              endif

            enddo

            if (count .gt. 0) then
                if (init) then
                   ease_f%icecon_next(i,j) = sum/dble(count)
                endif
                ease_f%icecon(i,j) = ease_f%icecon_next(i,j)
                ease_f%icecon_next(i,j) = sum/dble(count)
            else
                if (init) then
                   ease_f%icecon_next(i,j) = 0.d0
                endif
                ease_f%icecon(i,j) = ease_f%icecon_next(i,j)
                ease_f%icecon_next(i,j) = 0.d0
            endif

          enddo
        enddo

      end SUBROUTINE sia2_forcing_update_icecon


!=======================================================================
! Subroutine: sia2_forcing_icevec
! Purpose: loads ssmi sea ice vectors, keeping 5
! days in memory, for use later inconstructing a 5-day moving average
! (4 days behind + current day)
!=======================================================================

     SUBROUTINE sia2_forcing_icevec(ease_f)

		  use sia2_globals
		  implicit none

      ! function arguments
      ! ------------------------------------------------------------
      type (ease_f_type) :: ease_f

      ! local variables
      ! ------------------------------------------------------------
      integer :: next_yr
      double precision :: next_hr

      print *,'Smoothing ice vectors over last 5 days...'

      ! pre-load past snowdepths for use in first time step
      ! ------------------------------------------------------------
		  if (start .or. restart .eq. 1) then
              ! load 5 days worth of back ice vectors
!			  call sia2_forcing_icevec1(ease_f,1,cur_year,cur_hour-96.)
!			  call sia2_forcing_icevec1(ease_f,2,cur_year,cur_hour-72.)
!			  call sia2_forcing_icevec1(ease_f,3,cur_year,cur_hour-48.)
!			  call sia2_forcing_icevec1(ease_f,4,cur_year,cur_hour-24.)
!			  call sia2_forcing_icevec1(ease_f,5,cur_year,cur_hour)

!			  call sia2_forcing_icevec1(ease_f,1,cur_year,cur_hour-48.)
!			  call sia2_forcing_icevec1(ease_f,2,cur_year,cur_hour-24.)
!			  call sia2_forcing_icevec1(ease_f,3,cur_year,cur_hour)

			  call sia2_forcing_icevec1(ease_f,1,cur_year,cur_hour-24.)
			  call sia2_forcing_icevec1(ease_f,2,cur_year,cur_hour)

        ! create first mean snowdepth
        call sia2_forcing_update_icevec(.true.,ease_f)
        ! find next day's snow & load
        call sia2_forcing_next_day(cur_year,cur_hour,next_yr,next_hr,leap_year)
        call sia2_forcing_icevec1(ease_f,1,next_yr,next_hr)
        ! find mean snow for next day
        call sia2_forcing_update_icevec(.false.,ease_f)
        icevec_index = 1

		  else

        ! find yr/hr which will load next timestep's snow
        call sia2_forcing_next_day(cur_year,cur_hour+dt,next_yr,next_hr,leap_year)

        ! increment index of snowpast array
!			  if (icevec_index .eq. 5) then
!			  if (icevec_index .eq. 3) then
			  if (icevec_index .eq. 2) then
				  icevec_index = 1
			  else
				  icevec_index = icevec_index+1
			  endif

			  ! load next timestep's snow
			  call sia2_forcing_icevec1(ease_f,icevec_index,next_yr,next_hr)

			  ! interpolate over 7 days, and move snowh_next -> snowh
			  call sia2_forcing_update_icevec(.false.,ease_f)

      endif

     end SUBROUTINE sia2_forcing_icevec

!=======================================================================
! Subroutine: sia2_forcing_icecon
! Purpose: loads ssmi ice concectrations, keeping 5
! days in memory, for use later inconstructing a 5-day moving average
! (4 days behind + current day)
!=======================================================================

     SUBROUTINE sia2_forcing_icecon(ease_f)

		  use sia2_globals
		  implicit none

          ! function arguments
          ! ------------------------------------------------------------
          type (ease_f_type) :: ease_f

          ! local variables
          ! ------------------------------------------------------------
          integer :: next_yr
          double precision :: next_hr

          ! pre-load past snowdepths for use in first time step
          ! ------------------------------------------------------------
		  if (start .or. restart .eq. 1) then
              ! load 5 days worth of back ice vectors
!			  call sia2_forcing_icecon1(ease_f,1,cur_year,cur_hour-96.)
!			  call sia2_forcing_icecon1(ease_f,2,cur_year,cur_hour-72.)
!			  call sia2_forcing_icecon1(ease_f,3,cur_year,cur_hour-48.)
!			  call sia2_forcing_icecon1(ease_f,4,cur_year,cur_hour-24.)
!			  call sia2_forcing_icecon1(ease_f,5,cur_year,cur_hour)

!			  call sia2_forcing_icecon1(ease_f,1,cur_year,cur_hour-48.)
!			  call sia2_forcing_icecon1(ease_f,2,cur_year,cur_hour-24.)
!			  call sia2_forcing_icecon1(ease_f,3,cur_year,cur_hour)

			  call sia2_forcing_icecon1(ease_f,1,cur_year,cur_hour-24.)
			  call sia2_forcing_icecon1(ease_f,2,cur_year,cur_hour)

        ! create first mean snowdepth
        call sia2_forcing_update_icecon(.true.,ease_f)
        ! find next day's snow & load
        call sia2_forcing_next_day(cur_year,cur_hour,next_yr,next_hr,leap_year)
        call sia2_forcing_icecon1(ease_f,1,next_yr,next_hr)
        ! find mean snow for next day
        call sia2_forcing_update_icecon(.false.,ease_f)
        icecon_index = 1

		  else

        ! find yr/hr which will load next timestep's snow
        call sia2_forcing_next_day(cur_year,cur_hour+dt,next_yr,next_hr,leap_year)

        ! increment index of snowpast array
!			  if (icecon_index .eq. 5) then
!			  if (icecon_index .eq. 3) then
			  if (icecon_index .eq. 2) then
				  icecon_index = 1
			  else
				  icecon_index = icecon_index+1
			  endif

			  ! load next timestep's snow
			  call sia2_forcing_icecon1(ease_f,icecon_index,next_yr,next_hr)

			  ! interpolate over 7 days, and move snowh_next -> snowh
			  call sia2_forcing_update_icecon(.false.,ease_f)

      endif

     end SUBROUTINE sia2_forcing_icecon


!=======================================================================
! Subroutine: sia2_forcing_icevec1
! Purpose: loads one day of ice vectors
!=======================================================================
      SUBROUTINE sia2_forcing_icevec1(ease_f,ivindex,year,hour)
              use sia2_globals
              implicit none

          ! function arguments
          ! ------------------------------------------------------------
          integer :: year,ivindex,i,j,kk,offset,icevec_dim
          type (ease_f_type) :: ease_f
          double precision :: hour

		  !  Open 25km snow depth
		  !----------------------
		  !  Snow depth in cm &
		  !  251=bad brightness temp
		  !  252=no sea ice
		  !  values above 60cm should be capped, and values < zero set to zero
		  !  use a 5-day average, mark all pixels with derivative of
		  !  greater than +or-5 cm/day as bad!  Interpolate?
		  !-----------------------------------------------------------------
      if (hemi .eq. 0) then
          icevec_dim = 321
          offset = 200 ! 25km ice vectors files are only 321x321
      else
          icevec_dim = 361
          offset = 180
      endif

		  call SSMI_icevec_read(year, int(hour/24)+1)   ! reading into global SSMI_grid_int2
		  do kk=1,3
			  do j=1,icevec_dim
				  do i=1,icevec_dim
					  ease_f%icevec_past(offset+i,offset+j,kk,ivindex) = &
						dfloat(icevec_grid_int2(i,j,kk))
				  enddo
			  enddo
		  enddo

      end SUBROUTINE sia2_forcing_icevec1

!=======================================================================
! Subroutine: sia2_forcing_update_icevec
! Purpose: at each time step to estimate/interpolate ice vectors
! on a 5 day mean
!=======================================================================
      SUBROUTINE sia2_forcing_update_icevec(init,ease_f)
          use sia2_globals
          implicit none


          ! function arguments
          ! ------------------------------------------------------------
          logical :: init
          type (ease_f_type) :: ease_f

          ! local variables
          ! ------------------------------------------------------------
          double precision :: sum_u,sum_v,new_u,new_v
          integer :: i,j,ii,count,offset

          print *,'Smoothing ice vectors over last 5 days...'

          offset = 200
          do j= 1+offset,321+offset
              do i = 1+offset,321+offset

				  sum_u = 0.
				  sum_v = 0.
				  count = 0
!				  do ii = 1,5
!				  do ii = 1,3
				  do ii = 1,2
                      if (ease_f%icevec_past(i,j,3,ii) .gt. 0.) then
                          sum_u = sum_u + ease_f%icevec_past(i,j,1,ii)
                          sum_v = sum_v + ease_f%icevec_past(i,j,2,ii)
                          count = count + 1
                      endif
                  enddo

                  if (count .gt. 0) then
                      if (init) then
                         ease_f%ice_vec_next(i,j,1) = sum_u/dble(count)
                         ease_f%ice_vec_next(i,j,2) = sum_v/dble(count)
                         ease_f%ice_vec_next(i,j,3) = 1.
                      endif
                      ease_f%ice_vec(i,j,:) = ease_f%ice_vec_next(i,j,:)
                      ease_f%ice_vec_next(i,j,1) = sum_u/dble(count)
                      ease_f%ice_vec_next(i,j,2) = sum_v/dble(count)
                      ease_f%ice_vec_next(i,j,3) = 1.
                  else
                      if (init) then
                         ease_f%ice_vec_next(i,j,1) = 0.d0
                         ease_f%ice_vec_next(i,j,2) = 0.d0
                         ease_f%ice_vec_next(i,j,3) = 0.d0
                      endif
                      ease_f%ice_vec(i,j,:) = ease_f%ice_vec_next(i,j,:)
                      ease_f%ice_vec_next(i,j,1) = 0.d0
                      ease_f%ice_vec_next(i,j,2) = 0.d0
                      ease_f%ice_vec_next(i,j,3) = 0.d0
                  endif

              enddo
          enddo

      end SUBROUTINE sia2_forcing_update_icevec


!=======================================================================
! Subroutine: sia2_forcing_next_ncep
! Purpose: finds the year and hour of the next available ncep forcing data
!=======================================================================
      SUBROUTINE sia2_forcing_next_ncep(year,hour,freq,next_year,next_hour,leap)

		  implicit none

		  integer :: year,next_year
          double precision :: hour,next_hour,freq
		  logical :: leap,eoy,fake_leap

		  fake_leap = leap
		  next_year = year
		  next_hour = hour + freq

		  call sia2_end_of_year(next_year,next_hour,eoy,fake_leap)

		  if (eoy) then
			  call sia2_increment_year(next_year,next_hour,fake_leap)
              print *,'Year incremented!! EOY = ',eoy
		  endif

      end SUBROUTINE sia2_forcing_next_ncep


!=======================================================================
! Subroutine: sia2_forcing_next_day
! Purpose: finds the year and hour of the next available daily forcing data
!=======================================================================
      SUBROUTINE sia2_forcing_next_day(year,hour,next_year,next_hour,leap)

		  implicit none

		  integer :: year,next_year
          double precision :: hour,next_hour
		  logical :: leap,eoy,fake_leap

		  fake_leap = leap
		  next_year = year
		  next_hour = hour + 24.

		  call sia2_end_of_year(next_year,next_hour,eoy,fake_leap)

		  if (eoy) then
			  call sia2_increment_year(next_year,next_hour,fake_leap)
		  endif

      end SUBROUTINE sia2_forcing_next_day


!=======================================================================
! Subroutine: sia2_swapendian_int2
! Purpose: Some data is in little endian, int2 format.  This utilty should
! swap it to big endian, which this program was written to operate in
!=======================================================================
      SUBROUTINE sia2_swapendian_int2(int_in)

          implicit none

          integer*2 :: int_in
          integer*1 :: bytes_in(2),bytes_out(2)

          bytes_in = transfer(int_in,bytes_in)
          bytes_out(1) = bytes_in(2)
          bytes_out(2) = bytes_in(1)
          int_in = transfer(bytes_out,int_in)

      end SUBROUTINE sia2_swapendian_int2

!=======================================================================
! Subroutine: sia2_swapendian_int8
! Purpose: Some data is in little endian, int8 format.  This utilty should
! swap it to big endian, which this program was written to operate in
!=======================================================================
      SUBROUTINE sia2_swapendian_int8(int_in)

          implicit none

          integer*8 :: int_in
          integer*1 :: bytes_in(8),bytes_out(8)

          bytes_in = transfer(int_in,bytes_in)
          bytes_out(1) = bytes_in(8)
          bytes_out(2) = bytes_in(7)
          bytes_out(3) = bytes_in(6)
          bytes_out(4) = bytes_in(5)
          bytes_out(5) = bytes_in(4)
          bytes_out(6) = bytes_in(3)
          bytes_out(7) = bytes_in(2)
          bytes_out(8) = bytes_in(1)

          int_in = transfer(bytes_out,int_in)

      end SUBROUTINE sia2_swapendian_int8

