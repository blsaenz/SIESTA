! SUBROUTINE: sia2_interpolate_forcing
! at each time step to estimate/interpolate air and/or surface temperature
! between forcing data
! ======================================================================

      SUBROUTINE sia2_interpolate_forcing(pj,m,f,ease_f,mp_f,dg_f,ec_f,eci_f, &
      v_f, stations_mdh,stations_mdv)

          use sia2_globals
          implicit none

          ! Subroutine Arguments
          ! ------------------------------------------------------------
          type (proj_type) :: pj
          type (meta_type) :: m(tcells)
          type (forcing_type) :: f(tcells)
          type (ease_f_type) :: ease_f
          type (mp_f_type) :: mp_f
          type (dg_f_type) :: dg_f
          type (ec_f_type) :: ec_f
          type (eci_f_type) :: eci_f
          type (boundary_f_type) :: v_f
          integer, dimension(n_stations) :: stations_mdh, stations_mdv

          ! Shared Internal Variables
          ! ------------------------------------------------------------
          integer :: chunk_size
          double precision :: timing,timing2,s_timing,month,today

          ! Private Internal Variables
          ! ------------------------------------------------------------
          integer i,j,mi,i_mp,j_mp,i_dg,j_dg,i_ec,j_ec,i_eci,j_eci,ii,mo, &
              is,js,fw_high
          real :: d_real
          double precision :: u,v,flag,tmp1,tmp2,p_h2o,p_h2o_sat,row_air

          ! snow depth/ice concentration is being interpolated for the next time step, adding dt
          ! to cur_hour
          if (start) then
              s_timing = dt
              last_day = int(cur_hour/24)+1
              last3hour = int(cur_hour/3)+1
              last6hour = int(cur_hour/6)+1
              last12hour = int(cur_hour/12)+1
              last_year = cur_year
              last_hour = int(cur_hour)+1
          else
              s_timing = ((cur_hour+dt)/24+1)-last_day
          endif

          if (s_timing .ge. 1.) then
              s_timing = s_timing-1.
          endif

          today = int(cur_hour/24.)+1

          chunk_size = 16

          fw_high = 0

!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i,j,mi,i_mp,j_mp,i_dg,j_dg,i_ec,j_ec,i_eci,j_eci,ii,mo, &
!$OMP d_real,month,u,v,flag,tmp1,tmp2,p_h2o,p_h2o_sat,row_air,is,js)

          ! Interpolate atmospheric forcing over time
          ! ------------------------------------------------------------

          if (atm_f .eq. 0) then

              ! find how far between NCEP forcing data the current time step is
              if (ncep_f .eq. 6) then
                  timing = (cur_hour/6+1)-last6hour
              else
                  timing = (cur_hour/24+1)-last_day
              endif
              if (timing .ge. 1. .or. start) then
                  timing = 0.
              endif

!$OMP DO SCHEDULE (STATIC,chunk_size)

              ! Interpolate NCEP
              do j = 1,mp_y
                  do i = 1,mp_x
                      mp_f%ncep_interp%at(i,j) = mp_f%ncep%at(i,j) + &
                          (mp_f%ncep_next%at(i,j)-mp_f%ncep%at(i,j))*timing
                      mp_f%ncep_interp%p(i,j) = mp_f%ncep%p(i,j) + &
                          (mp_f%ncep_next%p(i,j)-mp_f%ncep%p(i,j))*timing
                      mp_f%ncep_interp%h(i,j) = mp_f%ncep%h(i,j) + &
                          (mp_f%ncep_next%h(i,j)-mp_f%ncep%h(i,j))*timing
                      mp_f%ncep_interp%fc(i,j) = mp_f%ncep%fc(i,j) + &
                          (mp_f%ncep_next%fc(i,j)-mp_f%ncep%fc(i,j))*timing
                      mp_f%ncep_interp%u10(i,j) = mp_f%ncep%u10(i,j) + &
                          (mp_f%ncep_next%u10(i,j)-mp_f%ncep%u10(i,j))*timing
                      mp_f%ncep_interp%v10(i,j) = mp_f%ncep%v10(i,j) + &
                          (mp_f%ncep_next%v10(i,j)-mp_f%ncep%v10(i,j))*timing
                      mp_f%ncep_interp%pr(i,j) = mp_f%ncep%pr(i,j) + &
                          (mp_f%ncep_next%pr(i,j)-mp_f%ncep%pr(i,j))*timing
                  enddo
              enddo

!$OMP END DO

          elseif (atm_f .eq. 1) then

              timing = (cur_hour/6+1)-last6hour
              if (timing .ge. 1. .or. start) then
                  timing = 0.
              endif

!$OMP DO SCHEDULE (STATIC,chunk_size)

              ! Interpolate ECMWF
              do j = 1,ec_y
                  do i = 1,ec_x

                      ec_f%ecmwf_interp%at(i,j) = ec_f%ecmwf%at(i,j) + &
                          (ec_f%ecmwf_next%at(i,j)-ec_f%ecmwf%at(i,j))*timing
                      ec_f%ecmwf_interp%p(i,j) = ec_f%ecmwf%p(i,j) + &
                          (ec_f%ecmwf_next%p(i,j)-ec_f%ecmwf%p(i,j))*timing
                      ec_f%ecmwf_interp%dpt(i,j) = ec_f%ecmwf%dpt(i,j) + &
                          (ec_f%ecmwf_next%dpt(i,j)-ec_f%ecmwf%dpt(i,j))*timing
                      ec_f%ecmwf_interp%fc(i,j) = ec_f%ecmwf%fc(i,j) + &
                          (ec_f%ecmwf_next%fc(i,j)-ec_f%ecmwf%fc(i,j))*timing
                      ec_f%ecmwf_interp%u10(i,j) = ec_f%ecmwf%u10(i,j) + &
                          (ec_f%ecmwf_next%u10(i,j)-ec_f%ecmwf%u10(i,j))*timing
                      ec_f%ecmwf_interp%v10(i,j) = ec_f%ecmwf%v10(i,j) + &
                          (ec_f%ecmwf_next%v10(i,j)-ec_f%ecmwf%v10(i,j))*timing
                      ec_f%ecmwf_interp%pr(i,j) = ec_f%ecmwf%pr(i,j) + &
                          (ec_f%ecmwf_next%pr(i,j)-ec_f%ecmwf%pr(i,j))*timing
                      ec_f%ecmwf_interp%ssr(i,j) = ec_f%ecmwf%ssr(i,j) + &
                          (ec_f%ecmwf_next%ssr(i,j)-ec_f%ecmwf%ssr(i,j))*timing
                  enddo
              enddo

!$OMP END DO

          elseif (atm_f .eq. 2) then

              timing = (cur_hour/6+1)-last6hour
              if (timing .ge. 1. .or. start) then
                  timing = 0.
              endif

              ! if a 12-hour period has passed, we start accumulation
              ! anew, so that the loaded pr_next is the total AND incremental precip
              ! and we don't need to subtract the last (pr) value
              if(int(cur_hour/12.)+1 .gt. last12hour .or. start) then
                  timing2 = 0.
              else
                  timing2 = 1.
              endif

!$OMP DO SCHEDULE (STATIC,chunk_size)

              ! Interpolate ECMWF
              do j = 1,eci_y
                  do i = 1,eci_x

                      ! 6 hour forcing variables
                      eci_f%ecmwf_interp%at(i,j) = eci_f%ecmwf%at(i,j) + &
                          (eci_f%ecmwf_next%at(i,j)-eci_f%ecmwf%at(i,j))*timing
                      eci_f%ecmwf_interp%p(i,j) = eci_f%ecmwf%p(i,j) + &
                          (eci_f%ecmwf_next%p(i,j)-eci_f%ecmwf%p(i,j))*timing
                      eci_f%ecmwf_interp%dpt(i,j) = eci_f%ecmwf%dpt(i,j) + &
                          (eci_f%ecmwf_next%dpt(i,j)-eci_f%ecmwf%dpt(i,j))*timing
                      eci_f%ecmwf_interp%fc(i,j) = eci_f%ecmwf%fc(i,j) + &
                          (eci_f%ecmwf_next%fc(i,j)-eci_f%ecmwf%fc(i,j))*timing
                      eci_f%ecmwf_interp%u10(i,j) = eci_f%ecmwf%u10(i,j) + &
                          (eci_f%ecmwf_next%u10(i,j)-eci_f%ecmwf%u10(i,j))*timing
                      eci_f%ecmwf_interp%v10(i,j) = eci_f%ecmwf%v10(i,j) + &
                          (eci_f%ecmwf_next%v10(i,j)-eci_f%ecmwf%v10(i,j))*timing

                      ! 3 hour forcing variables (accumulation at 12 internvals)
                      ! don't interpolate precip - it is a focecast total sum, not a rate
                       eci_f%ecmwf_interp%pr(i,j) = eci_f%ecmwf_next%pr(i,j) - eci_f%ecmwf%pr(i,j)

                  enddo
              enddo

!$OMP END DO

          endif


          ! Interpolate WOA 1 degree data over time
          ! ------------------------------------------------------------
          if (start .or. restart .eq. 1) then
              call sia2_patch_dg_grid(pj,dg_f,v_f,m)  ! patch holes in WOA data
          endif


          call sia2_getmonth_dble(cur_hour,leap_year,month)
          timing = month-dble(last_month)  ! % of month elapsed

!$OMP DO SCHEDULE (STATIC,chunk_size)
          do j = 1,dg_y
              do i = 1,dg_x
                  dg_f%woa_interp%t(i,j) = dg_f%woa%t(i,j) + &
                      (dg_f%woa_next%t(i,j)-dg_f%woa%t(i,j))*real(timing)
                  dg_f%woa_interp%s(i,j) = dg_f%woa%s(i,j) + &
                      (dg_f%woa_next%s(i,j)-dg_f%woa%s(i,j))*real(timing)
                  dg_f%woa_interp%n(i,j) = dg_f%woa%n(i,j) + &
                      (dg_f%woa_next%n(i,j)-dg_f%woa%n(i,j))*real(timing)
                  dg_f%woa_interp%p(i,j) = dg_f%woa%p(i,j) + &
                      (dg_f%woa_next%p(i,j)-dg_f%woa%p(i,j))*real(timing)
                  dg_f%woa_interp%si(i,j) = dg_f%woa%si(i,j) + &
                      (dg_f%woa_next%si(i,j)-dg_f%woa%si(i,j))*real(timing)
              enddo
          enddo

!$OMP END DO

          ! Interpolate SOSE data over time
          ! ------------------------------------------------------------
          if (use_sose .eq. 1) then
                  timing = (cur_hour/120.+1)-dble(last_5_day)  ! % of 5-day period elapsed
                    if (timing .ge. 1. .or. start) then
                        timing = 0.
                    endif
!$OMP DO SCHEDULE (STATIC,chunk_size)
                    do j = 1,grids_v
                        do i = 1,grids_h
                            sose_f%sose_interp%t(i,j) = sose_f%sose%t(i,j) + &
                                (sose_f%sose_next%t(i,j)-sose_f%sose%t(i,j))*timing
                            sose_f%sose_interp%s(i,j) = sose_f%sose%s(i,j) + &
                                (sose_f%sose_next%s(i,j)-sose_f%sose%s(i,j))*timing
                            sose_f%sose_interp%u(i,j) = sose_f%sose%u(i,j) + &
                                (sose_f%sose_next%u(i,j)-sose_f%sose%u(i,j))*timing
                            sose_f%sose_interp%v(i,j) = sose_f%sose%v(i,j) + &
                                (sose_f%sose_next%v(i,j)-sose_f%sose%v(i,j))*timing
                        enddo
                    enddo
!$OMP END DO
          endif


          ! iterate over model grid, interpolate daily vars and assign to forcing
          ! ---------------------------------------------------------------------

!$OMP DO SCHEDULE(DYNAMIC,chunk_size)
          do mi = 1,tcells

              ! Update position if performing a 1D/boundary_file model run
              ! -------------------------------------------------------------
              if (boundary_file .eq. 1 .and. v_f%lat .ne. c9999 .and. v_f%lon .ne. c9999) then

                  ! update model domain indexes, to account for the possibility that our station has moved
                  ! could be crap, effects writing of data
                  call sia2_ease25_from_ll(pj,v_f%lat,v_f%lon,v_f%mdh,v_f%mdv)
                  ! update station write-out vars
                  stations_mdh(mi) = v_f%mdh
                  stations_mdv(mi) = v_f%mdv
                  ! update modeled var projection reference
                  m(mi)%grid_h = v_f%mdh
                  m(mi)%grid_v = v_f%mdv
                  ! update projection references
                  pj%mi(v_f%mdh,v_f%mdv) = 1
                  mdh1 = v_f%mdh
                  mdv1 = v_f%mdv
                  mdh2=mdh1
                  mdv2=mdv1

                  m(mi)%lon = v_f%lon
                  m(mi)%lat = v_f%lat
                  m(mi)%x_mp = v_f%x_mp
                  m(mi)%y_mp = v_f%y_mp
                  m(mi)%x_ec = v_f%x_ec
                  m(mi)%y_ec = v_f%y_ec
                  m(mi)%x_eci = v_f%x_eci
                  m(mi)%y_eci = v_f%y_eci
                  m(mi)%x_dg = v_f%x_dg
                  m(mi)%y_dg = v_f%y_dg

              else  ! non-boundary_file, normal case

                  ! if only one cell, advection moves single cell about
                  if ((adv_on .eq. 1) .and. ((tcells .eq. 1) .or. (stn_only .eq. 1))) then
                      !if (start) then
                      !    m(mi)%pxl_h_offset = cell_side/2.
                      !    m(mi)%pxl_v_offset = cell_side/2.
                      !endif
                      u = ease_f%ice_vec(m(mi)%grid_h,m(mi)%grid_v,1)/1.e6*dt_s ! km/timestep
                      v = ease_f%ice_vec(m(mi)%grid_h,m(mi)%grid_v,2)/1.e6*dt_s ! km/timestep
                      flag = ease_f%ice_vec(m(mi)%grid_h,m(mi)%grid_v,3)

                      m(mi)%pxl_h_offset = m(mi)%pxl_h_offset + u
                      if (m(mi)%pxl_h_offset .gt. cell_side) then
                          m(mi)%pxl_h_offset = m(mi)%pxl_h_offset - cell_side
                          m(mi)%grid_h = m(mi)%grid_h + 1
                      elseif (m(mi)%pxl_h_offset .lt. 0.) then
                          m(mi)%pxl_h_offset = m(mi)%pxl_h_offset + cell_side
                          m(mi)%grid_h = m(mi)%grid_h - 1
                      endif

                      m(mi)%pxl_v_offset = m(mi)%pxl_v_offset + v
                      if (m(mi)%pxl_v_offset .gt. cell_side) then
                          m(mi)%pxl_v_offset = m(mi)%pxl_v_offset - cell_side
                          m(mi)%grid_v = m(mi)%grid_v + 1
                      elseif (m(mi)%pxl_h_offset .lt. 0.) then
                          m(mi)%pxl_v_offset = m(mi)%pxl_v_offset + cell_side
                          m(mi)%grid_v = m(mi)%grid_v - 1
                      endif

                      m(mi)%lon = pj%lon(m(mi)%grid_h,m(mi)%grid_v)
                      m(mi)%lat = pj%lat(m(mi)%grid_h,m(mi)%grid_v)
                      m(mi)%x_mp = pj%x_mp(m(mi)%grid_h,m(mi)%grid_v)
                      m(mi)%y_mp = pj%y_mp(m(mi)%grid_h,m(mi)%grid_v)
                      m(mi)%x_ec = pj%x_ec(m(mi)%grid_h,m(mi)%grid_v)
                      m(mi)%y_ec = pj%y_ec(m(mi)%grid_h,m(mi)%grid_v)
                      m(mi)%x_dg = pj%x_dg(m(mi)%grid_h,m(mi)%grid_v)
                      m(mi)%y_dg = pj%y_dg(m(mi)%grid_h,m(mi)%grid_v)
                  endif
              endif

  !print *,'mdh,mdv:',m(1)%grid_h,m(1)%grid_v
  !print *,'mdh,mdv:',v_f%mdh,v_f%mdv
  !print *,'lat,lon:',v_f%lat,v_f%lon

              ! find projection indexes for current grid cell
              ! --------------------------------------------------------
              i = m(mi)%grid_h  ! ease projection x
              j = m(mi)%grid_v  ! ease projection y

              is = i-200  ! ease projection x, 321x321
              js = j-200  ! ease projection y, 321x321


              i_mp = m(mi)%x_mp  ! mercator/ncep projection x
              j_mp = m(mi)%y_mp  ! mercator/ncep projection y

              i_ec = m(mi)%x_ec  ! ECMWF 2.5 degree projection x
              j_ec = m(mi)%y_ec  ! ECMWF 2.5 degree projection y

              i_eci = m(mi)%x_eci  ! ECMWF 1.5 projection x
              j_eci = m(mi)%y_eci  ! ECMWF 1.5 projection y

              i_dg = m(mi)%x_dg   ! 1 degree grid projection x
              j_dg = m(mi)%y_dg    ! 1 degree grid projection y


              ! Assign/Correct Atmospheric Forcing
              ! --------------------------------------------------------

              ! Air Temp Forcing
              if (boundary_file .eq. 1 .and. v_f%at .ne. c9999) then
                  f(mi)%at = v_f%at + kelvin0
                  !f(mi)%at = mp_f%ncep_interp%at(i_mp,j_mp)  ! turn on to bypass boundary_file temp
              elseif (atm_f .eq. 0) then
                  f(mi)%at = mp_f%ncep_interp%at(i_mp,j_mp)
              elseif (atm_f .eq. 1) then
                  f(mi)%at = ec_f%ecmwf_interp%at(i_ec,j_ec)
              elseif (atm_f .eq. 2) then
                  f(mi)%at = eci_f%ecmwf_interp%at(i_eci,j_eci)
              endif
              if (f(mi)%at .lt. 218.15) then
                  f(mi)%at = 218.15 ! limit air temp to -55, b/c tubulent latent heat flux craps out if airtemp goes too low ...
              endif

              ! Air Pressure Forcing
              v_f%p = c9999
              if (boundary_file .eq. 1 .and. v_f%p .ne. c9999) then
                  f(mi)%p = v_f%p*100.
              elseif (atm_f .eq. 0) then
                  f(mi)%p = mp_f%ncep_interp%p(i_mp,j_mp)
              elseif (atm_f .eq. 1) then
                  f(mi)%p = ec_f%ecmwf_interp%p(i_ec,j_ec)
              elseif (atm_f .eq. 2) then
                  f(mi)%p = eci_f%ecmwf_interp%p(i_eci,j_eci)
              endif

              ! 10m Wind Speed Forcing
              v_f%ws = c9999
              if (boundary_file .eq. 1 .and. v_f%ws .ne. c9999) then
                  f(mi)%ws = v_f%ws
                  f(mi)%u10 = c9999
                  f(mi)%v10 = c9999
              elseif (atm_f .eq. 0) then
                  f(mi)%ws = c9999
                  f(mi)%u10 = mp_f%ncep_interp%u10(i_mp,j_mp)
                  f(mi)%v10 = mp_f%ncep_interp%v10(i_mp,j_mp)
              elseif (atm_f .eq. 1) then
                  f(mi)%ws = c9999
                  f(mi)%u10 = ec_f%ecmwf_interp%u10(i_ec,j_ec)
                  f(mi)%v10 = ec_f%ecmwf_interp%v10(i_ec,j_ec)
              elseif (atm_f .eq. 2) then
                  f(mi)%ws = c9999
                  f(mi)%u10 = eci_f%ecmwf_interp%u10(i_eci,j_eci)
                  f(mi)%v10 = eci_f%ecmwf_interp%v10(i_eci,j_eci)
              endif

              ! Humidity Forcing
              v_f%rhum = c9999
                ! water vapor saturation pressure
                p_h2o_sat = 6.1078d0*10**((7.5d0*f(mi)%at-2048.625d0)/ &
                (f(mi)%at-35.85d0)) ! mb (1 mb = 100 Pa = 1 hPa)
              if (boundary_file .eq. 1 .and. ((v_f%rhum .ne. c9999) .or. &
                  (v_f%shum .ne. c9999))) then
                  if (v_f%shum .eq. c9999) then
                      f(mi)%rhum = v_f%rhum/100.
                      if (f(mi)%rhum .gt. 1) then
                          f(mi)%rhum = 1.
                      endif
                      p_h2o = p_h2o_sat*f(mi)%rhum
                      f(mi)%h = 0.622d0*p_h2o/(f(mi)%p*0.01d0 - 1.622d0*p_h2o)
                  else
                      f(mi)%h = v_f%shum ! specific humidity
                      row_air = &
                          f(mi)%p/(f(mi)%at*(R_air - f(mi)%h*R_air + f(mi)%h*R_h2o))
                      f(mi)%rhum = row_air*f(mi)%h*R_h2o*f(mi)%at/(p_h2o_sat*100.)
                  endif
              elseif (atm_f .eq. 0) then
                  f(mi)%h = mp_f%ncep_interp%h(i_mp,j_mp) ! specific humidity
                  row_air = &
                      f(mi)%p/(f(mi)%at*(R_air - f(mi)%h*R_air + f(mi)%h*R_h2o))
                  f(mi)%rhum = row_air*f(mi)%h*R_h2o*f(mi)%at/(p_h2o_sat*100.)
              elseif (atm_f .ge. 1) then
                  ! saturation vapor pressure at air temp
                  f(mi)%h = 6.1078*10**((7.5*f(mi)%at-2048.625)/(f(mi)%at-35.85))  ! mb
                  ! saturation vapor pressure at dewpoint
                  if (atm_f .eq. 1) then
                      f(mi)%rhum = ec_f%ecmwf_interp%dpt(i_ec,j_ec) ! K
                  elseif (atm_f .eq. 2) then
                      f(mi)%rhum = eci_f%ecmwf_interp%dpt(i_eci,j_eci) ! K
                  endif
                  f(mi)%rhum = 6.1078*10**((7.5*f(mi)%rhum-2048.625)/(f(mi)%rhum-35.85))  ! mb
                  ! relative humidity = dewpoint saturation pressure / airtemp saturation pressure
                  f(mi)%rhum = f(mi)%rhum/f(mi)%h ! fraction
                  p_h2o = p_h2o_sat*f(mi)%rhum
                  f(mi)%h = 0.622d0*p_h2o/(f(mi)%p*0.01d0 - 1.622d0*p_h2o)
              endif

              ! Short Wave irradiance forcing
              if (boundary_file .eq. 1) then
                  f(mi)%swd = v_f%swd
                  !f(mi)%sens_heat_flux = v_f%sens_heat_flux
                  !f(mi)%late_heat_flux = v_f%late_heat_flux

                  !f(mi)%swd = ec_f%ecmwf_interp%ssr(i_ec,j_ec)/21960.  ! turn on to bypass boundary_file downweling short wave
                  !f(mi)%swd = c9999
              else
                  f(mi)%swd = c9999
                  !f(mi)%sens_heat_flux = c9999
                  !f(mi)%late_heat_flux = c9999
              endif

              ! Long Wave irradiance forcing
              if (boundary_file .eq. 1) then
                  f(mi)%lwd = v_f%lwd
                  f(mi)%lwd = c9999
              else
                  f(mi)%lwd = c9999
              endif

              ! Cloud Cover Forcing
              v_f%fc = c9999
              if (boundary_file .eq. 1 .and. v_f%fc .ne. c9999) then
                  f(mi)%fc = v_f%fc*100. ! (%)
                  !f(mi)%fc = mp_f%ncep_interp%fc(i_mp,j_mp) ! fractional cloud cover (%)
              else if (atm_f .eq. 0) then
                  f(mi)%fc = mp_f%ncep_interp%fc(i_mp,j_mp) ! fractional cloud cover (%)
              elseif (atm_f .eq. 1) then
                  f(mi)%fc = ec_f%ecmwf_interp%fc(i_ec,j_ec)*100.
              elseif (atm_f .eq. 2) then
                  f(mi)%fc = eci_f%ecmwf_interp%fc(i_eci,j_eci)*100.
              endif
              if (f(mi)%fc .lt. 0.) then
                  f(mi)%fc = 0. ! fix fractional cloud cover coming from climatology files that drop below 0
              endif

              ! Precipitation Forcing
              if (boundary_file .eq. 1 .and. v_f%pr .ne. c9999) then
                  f(mi)%pr = v_f%pr               ! (kg/m^2/s)
!                  f(mi)%pr = ec_f%ecmwf_interp%pr(i_ec,j_ec) ! (m)
                  print *,'pr:',v_f%pr
!                  f(mi)%pr = f(mi)%pr*0.0463 ! m * 1/21600s * 1000 kg/m^3 = kg/m^3/s
!              elseif (boundary_file .eq. 1 .and. v_f%pr_snow .ne. c9999) then
!                  f(mi)%pr = v_f%pr_snow                  ! (kg/m^2/s)
                  !f(mi)%pr = ec_f%ecmwf_interp%pr(i_ec,j_ec) ! (m)
                  !f(mi)%pr = f(mi)%pr*0.0463 ! m * 1/21600s * 1000 kg/m^3 = kg/m^3/s
              elseif (atm_f .eq. 0) then
                  f(mi)%pr = mp_f%ncep_interp%pr(i_mp,j_mp) ! precipitation (kg/m^2/s)
              elseif (atm_f .eq. 1) then
                  f(mi)%pr = ec_f%ecmwf_interp%pr(i_ec,j_ec) ! (m)
                  f(mi)%pr = f(mi)%pr*0.0463 ! m * 1/21600s * 1000 kg/m^3 = kg/m^2/s
              elseif (atm_f .eq. 2) then
                  f(mi)%pr = eci_f%ecmwf_interp%pr(i_eci,j_eci) ! (m)
                  f(mi)%pr = f(mi)%pr*0.0925926 ! m * 1/10800s * 1000 kg/m^3 = kg/m^2/s
              endif
              if (f(mi)%pr .lt. 0.) then
                  f(mi)%pr = 0. ! fix precipitation coming from climatology files that drop below 0
              endif

              ! Assign/Correct Ocean Forcing
              ! --------------------------------------------------------

              ! nutrient climatology
              if (woa .eq. 1) then
                  f(mi)%no3 = dble(dg_f%woa_interp%n(i_dg,j_dg))
                  f(mi)%nh4 = 0.  ! assuming no nh4 in open water
                  f(mi)%po4 = dble(dg_f%woa_interp%p(i_dg,j_dg))
                  f(mi)%sioh4 = dble(dg_f%woa_interp%si(i_dg,j_dg))
                  f(mi)%poc = sw_poc
              else
                  f(mi)%no3 = sw_no3
                  f(mi)%nh4 = sw_nh4
                  f(mi)%po4 = sw_po4
                  f(mi)%sioh4 = sw_sio4
                  f(mi)%poc = sw_poc
              endif


              ! Temperature/Salinity forcing
              if (boundary_file .eq. 1 .and. v_f%t .ne. c9999) then
                  f(mi)%t = v_f%t
              !elseif (woa .eq. 1) then
              !    f(mi)%t = dble(dg_f%woa_interp%t(i_dg,j_dg))
              else
                  f(mi)%t = Tw
              endif

              if (boundary_file .eq. 1 .and. v_f%s .ne. c9999) then
                  f(mi)%s = v_f%s
              !elseif (woa .eq. 1) then
              !    f(mi)%s = dble(dg_f%woa_interp%s(i_dg,j_dg))
              else
                  f(mi)%s = Sw
              endif

              if (f(mi)%t .lt. mu*f(mi)%s) then
                  f(mi)%t = f(mi)%s*mu   ! correct temperature to prevent freezing
              endif

              if (boundary_file .eq. 1) then
                  if (v_f%no3 .ne. c9999) then
                    f(mi)%no3 = v_f%no3
                  endif
                  if (v_f%nh4 .ne. c9999) then
                      f(mi)%nh4 = v_f%nh4
                  endif
                  if (v_f%po4 .ne. c9999) then
                      f(mi)%po4 = v_f%po4
                  endif
                  if (v_f%sioh4 .ne. c9999) then
                      f(mi)%sioh4 = v_f%sioh4
                  endif
              endif

              !call sia2_sw_density(real(f(mi)%t),real(f(mi)%s),0.1,d_real)
              !f(mi)%d = dble(d_real)
              f(mi)%d = f(mi)%s*800.d0 + 1.d6


              ! interpolate ice concentration and snow depth
              ! ------------------------------------------------------------

              if (ease_f%icecon(i,j) .ge. fi .and. ease_f%icecon(i,j) .le. 1.d0) then
                  if (ease_f%icecon_next(i,j) .ge. fi .and. ease_f%icecon_next(i,j) .le. 1.d0) then

                      ! remember last values
                      f(mi)%sh_interp_last = f(mi)%sh_interp

                      ! assign next->current for 1st forcing
                      if (f(mi)%ic_interp_next .eq. 0.) then
                          f(mi)%sh_interp_next = ease_f%sh(i,j)
                          f(mi)%ic_interp_next = ease_f%icecon(i,j)
                      endif

                      ! carry over next sh, sd values to current values
                      f(mi)%sh_interp = f(mi)%sh_interp_next
                      f(mi)%ic_interp = f(mi)%ic_interp_next

                      ! determine next sh, sd values
                      if (ease_f%sh(i,j) .eq. c9999 .or. ease_f%sh_next(i,j) .eq. c9999) then
                          f(mi)%sh_interp_next = c9999
                      else
                          f(mi)%sh_interp_next = ease_f%sh(i,j) + &
                              (ease_f%sh_next(i,j)-ease_f%sh(i,j))*s_timing
                      endif
                      f(mi)%ic_interp_next = ease_f%icecon(i,j) + &
                          (ease_f%icecon_next(i,j)-ease_f%icecon(i,j))*s_timing

                  else
                      ! next icecon/snow is out of range/not valid.  Use current values

                      ! remember last values
                      f(mi)%sh_interp_last = f(mi)%sh_interp

                      ! assign next->current for 1st forcing
                      if (f(mi)%ic_interp_next .eq. 0.) then
                          f(mi)%sh_interp_next = ease_f%sh(i,j)
                          f(mi)%ic_interp_next = ease_f%icecon(i,j)
                      endif

                      ! carry over next sh, sd values to current values
                      f(mi)%sh_interp = f(mi)%sh_interp_next
                      f(mi)%ic_interp = f(mi)%ic_interp_next

                      ! carry over next sh, sd values to current values
                      f(mi)%sh_interp_next = ease_f%sh(i,j)
                      f(mi)%ic_interp_next = ease_f%icecon(i,j)
                  endif
              else
                  f(mi)%sh_interp_next = 0.
                  f(mi)%ic_interp_next = 0.
                  f(mi)%sh_interp = 0.
                  f(mi)%ic_interp = 0.
              endif

              if (ease_f%ice_vec(i,j,3) .gt. d0_) then

                  if (f(mi)%ivu_interp_next .eq. c9999) then
                      ! next value is crap, use current value
                      f(mi)%ivu_interp = ease_f%ice_vec(i,j,1)
                      f(mi)%ivv_interp = ease_f%ice_vec(i,j,2)
                  else
                      ! carry over next sh, sd values to current values
                      f(mi)%ivu_interp = f(mi)%ivu_interp_next
                      f(mi)%ivv_interp = f(mi)%ivv_interp_next
                  endif

                  if (ease_f%ice_vec_next(i,j,3) .gt. d0_) then

                      ! linearly interpolate to find next ice vectors
                      f(mi)%ivu_interp_next = ease_f%ice_vec(i,j,1) + &
                          (ease_f%ice_vec_next(i,j,1)-ease_f%ice_vec(i,j,1))*s_timing
                      f(mi)%ivv_interp_next = ease_f%ice_vec(i,j,2) + &
                          (ease_f%ice_vec_next(i,j,2)-ease_f%ice_vec(i,j,2))*s_timing

                  else

                      ! next ice_vec is out of range/not valid, so use current until
                       f(mi)%ivu_interp_next = ease_f%ice_vec(i,j,1)
                      f(mi)%ivv_interp_next = ease_f%ice_vec(i,j,2)

                  endif
              else

                  if (f(mi)%ivu_interp_next .eq. c9999) then
                      f(mi)%ivu_interp = d0_
                      f(mi)%ivv_interp = d0_
                  else
                      f(mi)%ivu_interp = f(mi)%ivu_interp_next
                      f(mi)%ivv_interp = f(mi)%ivv_interp_next
                  endif

                  f(mi)%ivu_interp_next = d0_
                  f(mi)%ivv_interp_next = d0_

              endif


              if (boundary_file .eq. 1) then

                  if (v_f%sh .ne. c9999) then
                       f(mi)%sh_interp = v_f%sh  ! snow height from boundary_file
                  else
                      if (f(mi)%sh_interp .eq. c9999) then
                          print *,'Invalid snow depth in boundary_file mode.'
                          print *,'Check boundary_file data and snow_model compatibility.'
                          call exit(0)
                      endif
                  endif
                   f(mi)%ic_interp = 1.
                  f(mi)%ic_interp_next = 1.

              elseif(override_ic .eq. 1) then

                   f(mi)%ic_interp = 1.
                  f(mi)%ic_interp_next = 1.

              endif


                    if (use_sose .eq. 1 .and. m(mi)%status .eq. 1 .and. cur_year >= 2005) then

                        f(mi)%t = sose_f%sose_interp%t(is,js)
                        f(mi)%s = sose_f%sose_interp%s(is,js)
                        if (f(mi)%t .lt. mu*f(mi)%s) then
                                f(mi)%t = f(mi)%s*mu   ! correct temperature to prevent freezing
                        endif
                        f(mi)%d = f(mi)%s*800.d0 + 1.d6

                        call sia2_env_basal_heat(f(mi)%t,f(mi)%s,f(mi)%d, &
                                f(mi)%ivu_interp*1.d-3,f(mi)%ivv_interp*1.d-3, &
                                sose_f%sose_interp%u(is,js),sose_f%sose_interp%v(is,js), &
                                f(mi)%fw)

                        !print *, 'Calculated OHF: ',f(mi)%fw,' Grid: ',i,j

                        ! make 2 W/m^2 the minimum ocean heat flux to ice
                        f(mi)%fw = max(f(mi)%fw,2.d0)

                        ! make 60 W/m^2 the maximum ocean heat flux to ice
                        if (f(mi)%fw .gt. 60.d0) then
                            fw_high = fw_high + 1
                            f(mi)%fw = min(f(mi)%fw,60.d0)
                        endif

                    else

              ! Ocean Heat Flux Forcing
              if (boundary_file .eq. 1 .and. v_f%fw .ne. c9999) then
                  f(mi)%fw = v_f%fw
              elseif (ohf .eq. 1) then
                  if (ohf .eq. 1) then
                      f(mi)%fw = f(mi)%d*cw*0.006*0.01*(f(mi)%t - mu*f(mi)%s) ! from McFee et al. 1999
                  elseif (ohf .eq. 2) then

                      if (today .lt. 16) then   ! Jand = 20 W/m^2
                          f(mi)%fw = 16. + (today+16)/31*4.   ! Jan 1-15
                      elseif (today .lt. 45) then   ! 20
                          f(mi)%fw = 20.
                      elseif (today .lt. 76) then   ! 16
                          f(mi)%fw = 20. - (today-45)/31*4.   ! Jan 1-15
                      elseif (today .lt. 106) then   ! 12
                          f(mi)%fw = 16. - (today-76)/30*4.   ! Jan 1-15
                      elseif (today .lt. 137) then   ! 9
                          f(mi)%fw = 12. - (today-106)/31*3.   ! Jan 1-15
                      elseif (today .lt. 167) then   ! 6
                          f(mi)%fw = 9. - (today-137)/30*3.   ! Jan 1-15
                      elseif (today .lt. 198) then   ! 4
                          f(mi)%fw = 6. - (today-167)/30*2.   ! Jan 1-15
                      elseif (today .lt. 229) then   ! 4
                          f(mi)%fw = 4.   ! Jan 1-15
                      elseif (today .lt. 259) then   ! 2
                          f(mi)%fw = 4. - (today-229)/30*2.   ! Jan 1-15
                      elseif (today .lt. 290) then   ! 2
                          f(mi)%fw = 2.   ! Jan 1-15
                      elseif (today .lt. 320) then   ! 5
                          f(mi)%fw = 2. + (today-290)/30*3.   ! Oct 16-Nov14
                      elseif (today .lt. 351) then   ! Dec = 14 W/m^2
                          f(mi)%fw = 5. + (today-320)/31*9.   ! Nov 15-Dec154
                      else
                          f(mi)%fw = 14. + (today-351)/31*6.   ! Dec 15-
                      endif

                  else
                      f(mi)%fw = Fw
                  endif
              endif
            endif


                  !f(mi)%fw = 5.+20.*abs(1.83+f(mi)%t)/(1.83-1.60)
                  !print *,'Fw: ',f(mi)%fw

              !print *,'T: ',f(mi)%t
              !print *,'Tf: ',mu*f(mi)%s
              !print *,'Ocean Heat Flux: ',f(mi)%d*cw*0.006*0.01*(f(mi)%t - mu*f(mi)%s)



              ! correct bogus snow depths back to max_snow_h
!                  if (f(mi)%sh_interp/100. .gt. max_snow_h) then
!                      f(mi)%sh_interp = max_snow_h*100.
!                  endif

              call sia2_getmonth(cur_hour,leap_year,mo)
              select case (mo)
                  case (1)
                      tmp1 = 0.4444  ! airtemp skew seasonal adjustment
                      tmp2 = 0.6486  ! snowdepth skew seasonal adjustment
                  case (2)
                      tmp1 = 0.4753  ! airtemp skew seasonal adjustment
                      tmp2 = 0.6682  ! snowdepth skew seasonal adjustment
                  case (3)
                      tmp1 = 0.5679  ! airtemp skew seasonal adjustment
                      tmp2 = 0.7267  ! snowdepth skew seasonal adjustment
                  case (4)
                      tmp1 = 0.7222  ! airtemp skew seasonal adjustment
                      tmp2 = 0.8243  ! snowdepth skew seasonal adjustment
                  case (5)
                      tmp1 = 0.9383  ! airtemp skew seasonal adjustment
                      tmp2 = 0.9610  ! snowdepth skew seasonal adjustment
                  case (6)
                      tmp1 = 1.2160  ! airtemp skew seasonal adjustment
                      tmp2 = 1.1366  ! snowdepth skew seasonal adjustment
                  case (7)
                      tmp1 = 1.5556  ! airtemp skew seasonal adjustment
                      tmp2 = 1.3514  ! snowdepth skew seasonal adjustment
                  case (8)
                      tmp1 = 1.2160  ! airtemp skew seasonal adjustment
                      tmp2 = 1.1366  ! snowdepth skew seasonal adjustment
                  case (9)
                      tmp1 = 0.9383  ! airtemp skew seasonal adjustment
                      tmp2 = 0.9610  ! snowdepth skew seasonal adjustment
                  case (10)
                      tmp1 = 0.7222  ! airtemp skew seasonal adjustment
                      tmp2 = 0.8243  ! snowdepth skew seasonal adjustment
                  case (11)
                      tmp1 = 0.5679  ! airtemp skew seasonal adjustment
                      tmp2 = 0.7267  ! snowdepth skew seasonal adjustment
                  case (12)
                      tmp1 = 0.4753  ! airtemp skew seasonal adjustment
                      tmp2 = 0.6682  ! snowdepth skew seasonal adjustment
              end select

              ! skew snow depth
              if (monte_carlo .eq. 1 .and. boundary_file .eq. 1) then

                  ii = int(mi/5)  ! get fives
                  ii = mi-ii*5    ! remainder

                  if (ii .eq. 1) then
                      f(mi)%at = -3.d0 + kelvin0
                  elseif (ii .eq. 2) then
                       f(mi)%at = -6.d0 + kelvin0
                  elseif (ii .eq. 3) then
                       f(mi)%at = -12.d0 + kelvin0
                  elseif (ii .eq. 4) then
                       f(mi)%at = -20.d0 + kelvin0
                  elseif (ii .eq. 0) then
                       f(mi)%at = -30.d0 + kelvin0
                  endif

                  if (mi .le. 5) then
                      f(mi)%sh_interp = 0.d0
                  elseif (mi .gt. 5 .and. mi .le. 10) then
                      f(mi)%sh_interp = 0.d0
                  elseif (mi .gt. 10 .and. mi .le. 15) then
                      f(mi)%sh_interp = 10.d0
                  elseif (mi .gt. 15 .and. mi .le. 20) then
                      f(mi)%sh_interp = 25.d0
                  elseif (mi .gt. 20 .and. mi .le. 25) then
                      f(mi)%sh_interp = 25.d0
                      f(mi)%fw = 15
                  elseif (mi .gt. 25 .and. mi .le. 30) then
                      f(mi)%sh_interp = 10.d0
                  elseif (mi .gt. 30 .and. mi .le. 35) then
                      f(mi)%sh_interp = 30.d0
                  elseif (mi .gt. 35 .and. mi .le. 40) then
                      f(mi)%sh_interp = 30.d0
                      f(mi)%fw = 15
                  endif

              else
                  if (snow_skew .ne. 0.) then ! if statement protects normal boundary_file case from tmp2 monthly adjustment
                      f(mi)%sh_interp = f(mi)%sh_interp*(snow_skew+1.d0)!*tmp2
                      f(mi)%pr = f(mi)%pr + f(mi)%pr*(snow_skew+1.d0)!*tmp2
                  endif
                  f(mi)%fw = f(mi)%fw + ohf_skew ! add skew to ocean heat flux
                  f(mi)%at = f(mi)%at + tmp1*at_skew ! add skew to air temp (0 is default so if doesn't matter if tmp1 is set)
              endif

              !f(mi)%fw = 2.


              !if (dt_step .eq. 0) then
              !f(mi)%at = kelvin0 + 3.
              !else
              !f(mi)%at = kelvin0 - 10.
              !endif

              ! print for debugging
              !print *,'SH Next: ',f(mi)%sh_interp_next
              !print *,'SH:      ',f(mi)%sh_interp


              ! heat test -----------------------------------------------
              !snowh1m = 0.0
              !if (dt_step .lt. 40) then
              !     airtemp1c = -5.
              !elseif (dt_step .lt. 50) then
              !     airtemp1c = -5. - (dt_step-40)*2.
              !elseif (dt_step .lt. 100) then
              !     airtemp1c = -24.
              !elseif (dt_step .lt. 110) then
              !     airtemp1c = -24. + (dt_step - 100)/2.
              !elseif (dt_step .lt. 160) then
              !     airtemp1c = -19.5
              !elseif (dt_step .lt. 250) then
              !     airtemp1c = 2.
              !elseif (dt_step .lt. 550) then
              !     airtemp1c = -30.
              !elseif (dt_step .eq. 550) then
              !     call exit(0)
              !endif

              !if (dt_step .lt. 400) then
              !     airtemp1c = -40.
              !else
              !     airtemp1c = 3.
              !endif

              !if (dt_step .eq. 1000) then
              !    call exit(0)
              !endif

              !if (dt_step .ge. 1000 .and. dt_step .le. 8000) then
              !    f(mi)%at = kelvin0+4.
              !endif

              !if (dt_step .ge. 4000) then
              !    f(mi)%sh_interp_next = 0.
              !endif

              ! mp_f%ncep_interp%at(i_mp,j_mp) = kelvin0 + airtemp1c
              ! heat test -----------------------------------------------

          enddo

!$OMP END DO
!$OMP END PARALLEL

                print *, 'Number of FW > 60: ',fw_high

      end SUBROUTINE sia2_interpolate_forcing


!=======================================================================
! Subroutine: sia2_patch_dg_grid
! Purpose: Around ocean boundaries, there may the cells where WOA
! forcing data is either bad or does not exist.  This subroutine will
! re-map these cells to WOA cells nearby with good data.
!=======================================================================
      SUBROUTINE sia2_patch_dg_grid(pj,dg_f,v_f,m)

          use sia2_globals
          implicit none

          ! function arguments
          !----------------------------------------------------
          type (proj_type) :: pj
          type (dg_f_type) :: dg_f
          type (boundary_f_type) :: v_f
          type (meta_type) :: m(tcells)

          ! local variables
          ! ------------------------------------------------------------
          integer :: mi,x,y,away

          do mi=1,tcells

              if (boundary_file .eq. 1) then
                  x = v_f%x_dg
                  y = v_f%y_dg
              else
                  x = m(mi)%x_dg
                  y = m(mi)%y_dg
              endif

              away = 1

              if (dg_f%woa%t(x,y) .gt. 50.) then

                  do while (dg_f%woa%t(x,y) .gt. 50.)  ! bad data is represented by a large number (1.e20)

                      ! search around cell for valid data, focusing
                      ! on northern neighbors first, since this is geared to southern hemisphere
                      if (dg_f%woa%t(x,y+away) .lt. 50.) then
                          x = x
                          y = y+away
                      else
                          if (dg_f%woa%t(x+away,y) .lt. 50.) then
                              x = x+away
                              y = y
                          else
                              if (dg_f%woa%t(x-away,y) .lt. 50.) then
                                  x = x-away
                                  y = y
                              else
                                  if (dg_f%woa%t(x+away,y+away) .lt. 50.) then
                                      x = x+away
                                      y = y+away
                                  else
                                      if (dg_f%woa%t(x-away,y+away) .lt. 50.) then
                                          x = x-away
                                          y = y+away
                                      else
                                          if (dg_f%woa%t(x+away,y-away) .lt. 50.) then
                                              x = x+away
                                              y = y-away
                                          else
                                              if (dg_f%woa%t(x-away,y-away) .lt. 50.) then
                                                  x = x-away
                                                  y = y-away
                                              else
                                                  if (dg_f%woa%t(x,y-away) .lt. 50.) then
                                                      x = x
                                                      y = y-away
                                                  endif
                                              endif
                                          endif
                                      endif
                                  endif
                              endif
                          endif
                      endif

                      away = away + 1

                  enddo

                  ! record new degree grid indexes
                  pj%x_dg(m(mi)%grid_h,m(mi)%grid_v) = x
                  pj%y_dg(m(mi)%grid_h,m(mi)%grid_v) = y
                  m(mi)%x_dg = x
                  m(mi)%y_dg = y

               endif

           enddo

       end SUBROUTINE sia2_patch_dg_grid


!=======================================================================
! Subroutine: sia2_sw_density
! Purpose: Uses Brydon et al. 1999 for seawater equation of state, returns
! density from potential temp,salinity,pres.
! t = -2 - 30 degC
! s = 30 - 38 %
! p = 0-50 MPa (surface pressure is ~100 kPa = 0.0 MPa
!=======================================================================
      SUBROUTINE sia2_sw_density(t,s,p,d)

          use sia2_globals
          implicit none

          ! function arguments
          !----------------------------------------------------
          real :: t,s,p,d

          ! local variables
          ! ------------------------------------------------------------
          real, dimension(7) :: C

          C = d_a + d_b*p + d_g*(p**2)
          d = 1.e6 + 1.e3*(C(1) + C(2)*t + C(3)*s + C(4)*t**2 + C(5)*s*t &
              + C(6)*t**3 + C(7)*s*t**2)

      end SUBROUTINE sia2_sw_density



!=======================================================================
! Subroutine: sia2_env_basal_heat
! Purpose: find the basal heat applied to the bottom of ice
!=======================================================================

    SUBROUTINE sia2_env_basal_heat(tmix,smix,dmix,ui,vi,uo,vo,fbot)

    ! Globals and USE Statements
    !----------------------------------------------------

        use sia2_globals
        implicit none

    ! Function Arguments
    !----------------------------------------------------
        double precision, intent (in) :: &
            tmix,                &    ! mixed layer temperature (degC)
            smix,                &    ! mixed layer salinity (psu)
            dmix,                &    ! mixed layer density (kg/m^3)
            ui,                    &    ! u-direction ice velocity (m/s)
            vi,                    &    ! v-direction ice velocity (m/s)
            uo,                    &    ! u-direction ocean velocity (m/s)
            vo                        ! v-direction ocean velocity (m/s)
        double precision, intent (out) :: &
            fbot                    ! basal heat flux to ice (positive = upward (ice melt))

    ! Internal Parameters
    !----------------------------------------------------
        double precision :: &
            ustar_min     = 1.d-3,      & ! minimum friction velocity
            ch_io             = 0.006d0,    & ! heat transfer coefficient (must be dimensionless, otherwise units don't seem to work...)
            cw_io             = 0.0055d0        ! ice-ocean drag coefficient

    ! Internal Variables
    !----------------------------------------------------
        double precision :: &
            tauw,                & ! ice-ocean stress
            ustar,            & ! friction velocity
            tf                        ! mixed layer freezing temperature

        tauw = cw_io*dmix*(sqrt( (ui+uo)**2 + (vi+vo)**2)) ! ice-ocean stress

        ustar = sqrt(tauw/dmix)

        ustar = max(ustar,ustar_min)

        tf = mu*smix

        fbot = dmix*cw*ch_io*ustar*(tmix - tf)

    end SUBROUTINE sia2_env_basal_heat


