

      SUBROUTINE sia2_edd_gpu_calc(mp_f,f,m,ice)

          use sia2_gpu
          use sia2_globals
          use omp_lib
          implicit none

          ! subroutine rrguments (shared)
          ! ------------------------------------------------------------
          type (mp_f_type) :: mp_f
          type (forcing_type) :: f(tcells)
          type (meta_type) :: m(tcells)
          type (ice_type) :: ice(sda_n,ida_n,tcells)

          ! shared internal variables
		  !----------------------------------------------------
          logical :: tcells_not_scanned,gpu_not_full
          integer :: i,chunk_size,wr_i32_size,gpu_i,mi_master,&
              wr_f32_size,read_f32_size,tcell_blocks,tcell_block_size,zz
          real :: logm
          double precision :: j_date,rad_date,solar_dec,rad_time,lda_f
          character (len=8) :: dte
          character (len=10) :: tme

          ! private internal variables
		  !----------------------------------------------------
          integer :: j,k,ii,jj,kk,mi,ic,sc,ki,k1,kz,mi_wavl,gpu_cell, &
              m_row,srftyp,i_mp,j_mp,int_z,sk_z,thr,mi_c,ic_c,sc_c, &
              sk_1,ice_1,ice_z,z_snow,nm_i,block_i,thread_i,tcbi, &
              dummy_i,gpu_reserve

          double precision :: ed_w_dir,ed_w_dif,wv_mul,bv_mean, &
              w_temp,zenith_angle,coszen,tmp1,tmp2,tmp3, &
              airtemp1c

          double precision, dimension(z_max+z_snow_max) :: &
              dth,     &
              r2st,    &
              g2st,    &
              k2st,    &
              am_sum

          ! code on, dude
		  !----------------------------------------------------

          chunk_size = 64 ! set OPENMP chunk size
          zz = z_max + z_snow_max

          lda_f = lda_d  ! init light/snow distribution divider

          mi_master = 0 ! initalize mi for counting below
          mi_c = 0 ! initalize mi for counting below
          ic_c = 100 ! inintalize current ice category (high)

          mi = 0  ! to catch bugs
          ic = 0  ! to catch bugs

          ! convert julian date to radians
          j_date = cur_hour/24+1
          rad_date = (j_date/365)*2.*pi

          ! Calculate solar declination (between -23deg27' and 23deg27') for use in
          ! Solar elevation equation
          solar_dec = 0.39637 - 22.9133*cos(rad_date) + 4.02543*sin(rad_date)  &
              - 0.3872*cos(2*rad_date) + 0.052*sin(2*rad_date)

          ! Convert degrees to radians
          solar_dec = pi*solar_dec/180.

          ! Convert time to radians - these equations seem to give midnight as noon,
          ! so I subtracted pi to the time to correct -
          rad_time = (j_date - int(j_date))*2.*pi - pi

          tcells_not_scanned = .TRUE.

          ! find interation chunk for OPENMP division
          tcell_block_size = gpu_tpb/ida_n
          tcell_blocks = tcells/tcell_block_size+1

          call DATE_AND_TIME(date=dte, time=tme)
          print *,'T before setup: ',tme

          do while (tcells_not_scanned)

! ----------------------------------------------------------------------
! filling GPU arrays - decided if cell is valid and prep counters, etc
! ----------------------------------------------------------------------
          gpu_i = -1   ! starts with an increment
          gpu_cell = 0
          gpu_not_full = .TRUE.

!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(j,k,ii,jj,kk,mi,ic,sc,ki,k1,kz,mi_wavl,gpu_cell, &
!$OMP m_row,srftyp,i_mp,j_mp,int_z,sk_z,thr,mi_c,ic_c,sc_c, &
!$OMP sk_1,ice_1,ice_z,z_snow,nm_i,block_i,thread_i,tcbi, &
!$OMP dummy_i,gpu_reserve, &
!$OMP ed_w_dir,ed_w_dif,wv_mul,bv_mean, &
!$OMP w_temp,zenith_angle,coszen,tmp1,tmp2,tmp3, &
!$OMP airtemp1c, &
!$OMP dth,     &
!$OMP r2st,    &
!$OMP g2st,    &
!$OMP k2st,    &
!$OMP am_sum)

!$OMP DO SCHEDULE(STATIC,chunk_size)

          do tcbi=1,tcell_blocks

          !print *,'Shit: ',mi_master,tcells,gpu_i,gpu_not_full,tcells_not_scanned
          gpu_reserve = 1
#ifdef omp_on
          gpu_reserve = omp_get_num_threads()
#endif
          gpu_reserve = gpu_reserve*ida_n*sda_n

          if (gpu_not_full .and. tcells_not_scanned) then

          do dummy_i=1,tcell_block_size

          if (gpu_not_full .and. tcells_not_scanned) then


!$OMP CRITICAL

              ! increment cell
              mi_master = mi_master + 1
              ! reserve gpu threads for each ice category in cell
              if ((gpu_i + gpu_reserve) .ge. gpu_max_threads) then
                  gpu_not_full = .FALSE.
                  mi_master = mi_master - 1
              elseif (mi_master .gt. tcells) then
                  tcells_not_scanned = .FALSE.
              else
                  mi_c = mi_master
              endif
              !kk = omp_get_thread_num()
              !print *,'mi_c: ',mi_c,'gpu_i:',gpu_i

!$OMP END CRITICAL

              if (gpu_not_full .and. tcells_not_scanned .and. m(mi_c)%status .gt. 0) then
                  i_mp = m(mi_c)%x_mp  ! mercator/ncep projection x
                  j_mp = m(mi_c)%y_mp  ! mercator/ncep projection y

                  ! init absorption, pur arrays
                  do j=1,ida_n
                  do k=1,sda_n
                      ice(k,j,mi_c)%pur = 0.d0
                      ice(k,j,mi_c)%PAR_bot = 0.d0
                      ! these below redundantly calculated for every valid category - could improve eff. some
                      ice(k,j,mi_c)%ed_w = 0.d0
                      ice(k,j,mi_c)%snow%ed_w = 0.d0
                      ice(k,j,mi_c)%ed0_nir_dir = 0.d0
                      ice(k,j,mi_c)%ed0_nir_dif = 0.d0
                      ice(k,j,mi_c)%Ed0 = 0.
                  enddo
                  enddo

                  ! Find Total shortwave irradiance in Watts/m^2
                  do ii=1,wavl

                      ! wavelength multiplier - 400 and 700 nm are valued at ~half the others, to make a total of 301
                      if (ii .eq. 1 .or. ii .eq. wavl) then
                          wv_mul = 5.5d0
                      else
                          wv_mul = 10.d0
                      endif

                      ! µEin -> watts conversion for wavelength
                      w_temp = quanta_to_watts(ii)*wv_mul  ! Watts / (µEin/m/s)

                      ! surface irradiance in Watts
                      ed_w_dir = mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii)*w_temp
                      ed_w_dif = mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii)*w_temp

                      ! find Ed0, NIR watts
                      ice(1,1,mi_c)%ed0_nir_dir = ice(1,1,mi_c)%ed0_nir_dir + ed_w_dir
                      ice(1,1,mi_c)%ed0_nir_dif = ice(1,1,mi_c)%ed0_nir_dif + ed_w_dif
                      ice(1,1,mi_c)%Ed0 = ice(1,1,mi_c)%Ed0 + ed_w_dir + ed_w_dif

                  enddo

                  ! copy Eds to other categories
                  do j=1,ida_n
                  do k=1,sda_n
                      ice(k,j,mi_c)%ed0_nir_dir = ice(1,1,mi_c)%ed0_nir_dir
                      ice(k,j,mi_c)%ed0_nir_dif = ice(1,1,mi_c)%ed0_nir_dif
                      ice(k,j,mi_c)%Ed0 = ice(1,1,mi_c)%Ed0
                  enddo
                  enddo

                  if (ice(1,1,mi_c)%Ed0 .gt. 0.d0) then

                      airtemp1c = f(mi_c)%at - kelvin0
                      zenith_angle = acos(sin(solar_dec)*sin(m(mi_c)%lat/360*6.2832) &
                      + cos(solar_dec)*cos(m(mi_c)%lat/360*6.2832)*cos(rad_time) )   ! radians
                      if (zenith_angle .lt. 0.d0) then
                          zenith_angle = 0.d0  ! prevent negative light values
                      endif
                      coszen = cos(zenith_angle)

                      do ic_c=1,ida_n
                      do sc_c=1, sda_n

                          if (ice(sc_c,ic_c,mi_c)%z .gt. 0) then

!$OMP CRITICAL

                              ! increment active gpu cell and gpu memory counters - after since C-style indexes
                              gpu_i = gpu_i + 1
                              gpu_cell = gpu_i
                              !print *, 'grabbing gpu_i ... ',gpu_i, '(',mi_c, ic_c,')'
!$OMP END CRITICAL

                              block_i = gpu_cell/gpu_tpb ! C-style reference, starting at 0
                              thread_i = gpu_cell - block_i*gpu_tpb ! C-style reference, starting at 0
                              k1 = block_i*f32_in_p_th*gpu_tpb + thread_i + 1 ! start of block&thread data


                              sk_z = ice(sc_c,ic_c,mi_c)%z
                              int_z = sk_z - z_sk
                              sk_1 = int_z + 1
                              z_snow = ice(sc_c,ic_c,mi_c)%snow%z
                              ice_z = sk_z + z_snow
                              ice_1 = z_snow + 1

                              ! 1nm wavelength index
                              nm_i = ii*10-9

! ----------------------------------------------------------------------
! calc scattering function(g2st) and scattering(r2st)
! ----------------------------------------------------------------------

                              ! estimation particle absorption from chlorophyll and detrital concentrations in ice
                              ! also estimate scattering in ice due to brine volume, which is constant over wavelength
                              am_sum = 0.
                              do jj=1,sk_z
                                  tmp3 = ice(sc_c,ic_c,mi_c)%bv(jj)*c_001
                                  tmp1 = max(1.d0,ice(sc_c,ic_c,mi_c)%smalg(jj)/c_chl*tmp3)  ! minimum 1mg/m^3 chla absorption
                                  do ii=1,wavl
                                      ! sum algal spectral absorption, used to find detrital absorption, ad in cuda kernel
                                      am_sum(jj) = am_sum(jj)+aph(ii*10-9)*tmp1*10.   ! x10 b/c 10 nm wavelength bins
                                  enddo

                                  ! sea ic_ce scattering
                                  if (ice(sc_c,ic_c,mi_c)%t(jj) .le. -22.) then
                                      ! hydrohalite crystal scattering regime
                                      r2st(jj+z_snow) = 2400.
                                      g2st(jj+z_snow) = 0.94
                                  else
                    !                  if (z_snow .eq. 0 .and. jj .le. 5) then
                    !                  if (jj .le. 5) then
                    !                      if (jj .eq. 1) then
                                             ! surface scattering layer
                    !                          r2st(jj+z_snow) = 900.
                    !                          g2st(jj+z_snow) = 0.94
                    !                      else ! drained layers
                    !                          r2st(jj+z_snow) = 100.
                    !                          g2st(jj+z_snow) = 0.94
                    !                      endif
                    !                  else

                                          tmp1 = max(2.,abs(ice(sc_c,ic_c,mi_c)%t(jj)))-2.
                    !                      r2st(jj+z_snow) = (18.*ice(ic_c,mi_c)%s(jj))*sqrt(tmp1)+15;   ! power law increase
                                          tmp2 = 9.6*ice(sc_c,ic_c,mi_c)%s(jj)
                                          r2st(jj+z_snow) = min(tmp2*8.,tmp2*tmp1)+15;   ! linear low-temp increase
                                          g2st(jj+z_snow) = 0.98

                    !                  endif


                    !                      if (ice(ic_c,mi_c)%t(jj) .le. -10.) then
                                              ! mirabilite crystal scattering regime
                    !                          r2st(jj+z_snow) = 300.-100.*(ice(ic_c,mi_c)%t(jj)+10.)/12.
                    !                          r2st(jj+z_snow) = 200.
                    !                          g2st(jj+z_snow) = 0.98
                    !                      else
                                              ! brine pocket scattering regime
                    !                          tmp2 = min(-1.8,ice(ic_c,mi_c)%t(jj))
                                              !r2st(jj+z_snow) = -971.1/(tmp2)**2+309.7
                    !                          r2st(jj+z_snow) = -636.5/abs(tmp2)+363.6240
                                              !r2st(jj+z_snow) = -307./abs(tmp2)+181.
                                              !r2st(jj+z_snow) = max(r2st(jj+z_snow),300.)
                    !                          g2st(jj+z_snow) = 0.98
                    !                      endif
                     !                 endif
                                  endif


                                  ! algal scattering - from Babin et al. 2003, scattering = 1.0m^2/gDryMass* mgC/m^3 * 1gC/1000mgC * 1gDryMass/0.19gC (Sicko-Goad et al. 1984)
                                  r2st(jj+z_snow) = r2st(jj+z_snow) + &
                                    (ice(sc_c,ic_c,mi_c)%smalg(jj) + ice(sc_c,ic_c,mi_c)%poc(jj)) &
                                    *tmp3*5.2632d-3

                                  ! record ice layer depth
                                  dth(jj+z_snow) = ice(sc_c,ic_c,mi_c)%th(jj)

                              enddo

                              ! deal with snow wavelength independent inherent optical properties
                              if (ice(sc_c,ic_c,mi_c)%snow%z .gt. 0) then

                                  ! surface type over ice = snow
                                  srftyp = 1

                                  ii = 1
                                  do jj=z_snow,1,-1

                                      ! find snow layer scattering (r)
                                      if (ice(sc_c,ic_c,mi_c)%snow%d(jj) .gt. 0.6) then
                                          if ((rs_switch .eq. 0 .and. airtemp1c .lt. 0.) .or. &
                                          (rs_switch .eq. 1 .and. ice(sc_c,ic_c,mi_c)%snow%t(jj) .le. -1.) .or. &
                                          (rs_switch .eq. 2 .and. ice(sc_c,ic_c,mi_c)%snow%melt(jj) .lt. 0.5)) then
                                              r2st(ii) = 3000
                                              g2st(ii) = 0.89
                                          else
                                              r2st(ii) = 900
                                              g2st(ii) = 0.94
                                              !r2st(ii) = 2000
                                              !g2st(ii) = 0.89
                                          endif
                                      else
                                          ! surface scattering layer
                                          r2st(ii) = 900
                                          g2st(ii) = 0.94
                                      endif

                                      ! record effective snow layer attenuation depth
                                      dth(ii) = ice(sc_c,ic_c,mi_c)%snow%th(jj)*snow_fudge

                                      ii = ii + 1
                                  enddo

                              else

                                  ! surface type over ice = air
                                  srftyp = 1

                              endif

        ! ----------------------------------------------------------------------
        ! store in input GPU arrays
        ! ----------------------------------------------------------------------
                              ki = k1
                              do jj=1,ice_z

                                  !ki = k1 + (3*(jj-1)+1)*gpu_tpb
                                  if (jj .le. z_snow) then

                                      ii = z_snow - jj + 1
                                      c_f32_in(ki) = real(ice(sc_c,ic_c,mi_c)%snow%d(ii))
                                      ki = ki + 3*gpu_tpb

                                  else
                                      ii = jj - z_snow
                                      c_f32_in(ki) = real(ice(sc_c,ic_c,mi_c)%bv(ii)*c_001)
                                      ki = ki + gpu_tpb
                                      c_f32_in(ki) = real(ice(sc_c,ic_c,mi_c)%smalg(ii))
                                      ki = ki + gpu_tpb
                                      c_f32_in(ki) = real(ice(sc_c,ic_c,mi_c)%poc(ii)*am_sum(ii))
                                      ki = ki + gpu_tpb

                                  endif

                              enddo

                              ki = k1 + (3*zz)*gpu_tpb
                              c_f32_in(ki) = real(coszen)
                              ki = ki + gpu_tpb
                              do jj=1,ice_z
                                  c_f32_in(ki) = real(g2st(jj))
                                  ki = ki + gpu_tpb
                                  c_f32_in(ki) = real(r2st(jj))
                                  ki = ki + gpu_tpb
                                  c_f32_in(ki) = real(dth(jj))
                                  ki = ki + gpu_tpb
                              enddo


                              ki = k1 +  (6*zz + 1)*gpu_tpb
                              do ii=1,wavl+1

                                  if (ii .le. wavl) then

                                      c_f32_in(ki) = real(mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii))
                                      ki = ki + gpu_tpb
                                      c_f32_in(ki) = real(mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii))
                                      ki = ki + gpu_tpb

                                  else

                                      c_f32_in(ki) = real(ice(sc_c,ic_c,mi_c)%ed0_nir_dir*(par_to_swd-1.d0))
                                      ki = ki + gpu_tpb
                                      c_f32_in(ki) = real(ice(sc_c,ic_c,mi_c)%ed0_nir_dif*(par_to_swd-1.d0))
                                      ki = ki + gpu_tpb

                                  endif

                              enddo

                              ! record GPU column calc metadata
                              mi_gpu_index(gpu_cell+1) = mi_c*100 + ic_c*10 + sc_c  ! record value that indicated both mi & ic (cell & category)

                              ki = block_i*i32_in_p_th*gpu_tpb + thread_i + 1 ! start i32_in block&thread data
                              c_i32_in(ki) = sk_z
                              ki = ki +  gpu_tpb ! start i32_in block&thread data
                              c_i32_in(ki) = z_snow
                              ki = ki +  gpu_tpb ! start i32_in block&thread data
                              c_i32_in(ki) = srftyp

                          endif ! end - is ice category active (z > 0)?

                      enddo ! end of snow category loop
                      enddo ! end of ice category loop

                  endif ! end - is there any light?

              endif ! end if stats/ice category active, or 3rd gpu full, or done with cells check

          endif ! end of 2nd gpu full, or done with cells check?

          enddo ! end of iteration over this mi block

          endif ! end of 1st gpu full, or done with cells check

          enddo ! mi block iteration

!$OMP END DO

!$OMP END PARALLEL

          ! Fill remaining threads in current gpu thread block so they don't
          ! execute
          ! ----------------------------------------------------------------
          block_i = gpu_i/gpu_tpb ! C-style reference, starting at 0
          thread_i = gpu_i - block_i*gpu_tpb ! C-style reference, starting at 0
          do ii = (thread_i+1),(gpu_tpb-1)
              k1 = block_i*i32_in_p_th*gpu_tpb + ii + 1 ! start of block&thread data
             ! set number of ice layers to 0, so gpu knows to skip calc
              c_i32_in(k1) = 0    ! nilyr = 0
          enddo

          ! if for some reason gpu_i exactly equals tcells, prevent extra iteration
          if (mi_master .eq. tcells) then
              tcells_not_scanned = .FALSE.
          endif

! ----------------------------------------------------------------------
! do GPU calc(s)
! ----------------------------------------------------------------------

          call DATE_AND_TIME(date=dte, time=tme)
          print *,'T before copy: ',tme

          ! convert zero-start counter to total number of active threads
          gpu_i = gpu_i + 1

          if (gpu_i .gt. 0) then

          wr_i32_size = (gpu_i + gpu_tpb - 1) / gpu_tpb  ! num valid data blocks
          wr_f32_size = wr_i32_size*f32_in_p_th*gpu_tpb
          wr_i32_size = wr_i32_size*i32_in_p_th*gpu_tpb

          call sia2_edd_gpu_copy_step_mem( &
              c_i32_in,    &
              c_f32_in,    &
              g_i32_in,    &
              g_f32_in,    &
              wr_i32_size, &
              wr_f32_size)

          call DATE_AND_TIME(date=dte, time=tme)
          print *,'T after copy: ',tme

          do i=1,lda_n

              logm = find_logm(i,lda_n)
              print *,'Calling CUDA Delta-Eddington solution ...', logm


              call sia2_edd_gpu_solution( &
                   g_i32_in,    &
                   g_f32_in,    &
                   g_f32_scr,   &
                   g_f32_out,   &
                   i32_in_p_th, &
                   f32_in_p_th, &
                   f32_scr_p_th,&
                   f32_out_p_th,&
                   gpu_tpb,     &
                   gpu_i,   &
                   zz,          &
                   logm)

              call sia2_edd_gpu_return( &
                   g_f32_out,     &
                   c_f32_out,     &
                   f32_out_p_th,  &
                   gpu_tpb,       &
                   gpu_i)


!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(j,k,ii,jj,kk,mi,ic,ki,k1,kz,mi_wavl,gpu_cell, &
!$OMP m_row,srftyp,i_mp,j_mp,int_z,sk_z,thr,mi_c,ic_c, &
!$OMP sk_1,ice_1,ice_z,z_snow,nm_i,block_i,thread_i,tcbi, &
!$OMP dummy_i,gpu_reserve, &
!$OMP ed_w_dir,ed_w_dif,wv_mul,bv_mean, &
!$OMP w_temp,zenith_angle,coszen,tmp1,tmp2,tmp3, &
!$OMP airtemp1c, &
!$OMP dth,     &
!$OMP r2st,    &
!$OMP g2st,    &
!$OMP k2st,    &
!$OMP am_sum)

!$OMP DO SCHEDULE(DYNAMIC,chunk_size)

              do j=1,gpu_i

                  ! determine memory indices, bookkeeping
                  ! -----------------------------------------------------
                  !thr = omp_get_thread_num()

                  mi = mi_gpu_index(j)/100
                  ic = (mi_gpu_index(j) - mi*100)/10
                  sc = mi_gpu_index(j) - mi*100 - ic*10

                  sk_z = ice(sc,ic,mi)%z
                  int_z = sk_z - z_sk
                  sk_1 = int_z + 1
                  z_snow = ice(sc,ic,mi)%snow%z
                  ice_z = sk_z + z_snow
                  ice_1 = z_snow + 1

                  ! find CUDA c_f32_out array mapping
                  block_i = (j-1)/gpu_tpb ! C-style reference, starting at 0
                  thread_i = (j-1) - block_i*gpu_tpb! C-style reference, starting at 0
                  ki = block_i*f32_out_p_th*gpu_tpb + thread_i + 1 ! start of block&thread data

                  ! map back apparent optical properties
                  ! -----------------------------------------------------
                  do k=1,ice_z
                      if (k .le. z_snow) then
                          ice(sc,ic,mi)%snow%ed_w(k) = ice(sc,ic,mi)%snow%ed_w(k) + &
                            dble(c_f32_out(ki))*lda_f
                          ki = ki + gpu_tpb  ! skip array index to ignore PUR in snow - not used
                      else
                          m_row = k - z_snow  ! find correct ice layer index
                          ice(sc,ic,mi)%ed_w(m_row) =ice(sc,ic,mi)%ed_w(m_row) + &
                            dble(c_f32_out(ki))*lda_f ! store absorbed watts
                          ki = ki + gpu_tpb
                          ice(sc,ic,mi)%pur(m_row,i) = dble(c_f32_out(ki))  ! store PUR
                      endif
                      ki = ki + gpu_tpb
                  enddo

                  !if (mi .eq. 10) then
                  !    print *,'PUR: ',ice(sc,ic,mi)%pur(:,i)
                  !    print *,'Ed_W: ',ice(sc,ic,mi)%ed_w
                  !    print *,'Ed_W_snow: ',ice(sc,ic,mi)%snow%ed_w
                  !endif

              enddo ! end of recording gpu run(s) data

 !$OMP END DO
 !$OMP END PARALLEL

              print *,'Delta-Eddington solution completed and stored. ',gpu_i,' total cells.'

          enddo  ! end of snow/light distribution loop
          endif

          enddo  ! end of all_cell_loaded test

          call DATE_AND_TIME(date=dte, time=tme)
          print *,'T after: ',tme

          !print *,'pur: ',ice(1,1)%pur
          !print *,'ed_w: ',ice(1,1)%ed_w
          !print *,'ed_w_snow: ',ice(1,1)%snow%ed_w


      end SUBROUTINE sia2_edd_gpu_calc

