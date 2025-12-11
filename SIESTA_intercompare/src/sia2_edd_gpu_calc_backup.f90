

      SUBROUTINE sia2_edd_gpu_calc(mp_f,f,m,ice,pur)      

          use sia2_gpu
          use sia2_globals
          use OMP_LIB
          implicit none
          
          ! subroutine rrguments (shared)
          ! ------------------------------------------------------------
          type (mp_f_type) :: mp_f
          type (forcing_type) :: f(tcells)
          type (meta_type) :: m(tcells)
          type (ice_type) :: ice(ida_n,tcells)
          integer (kind=2) :: pur(z_max,dt_per_day,lda_n,ida_n,tcells)

          ! shared internal variables
		  !----------------------------------------------------          
          logical :: all_cells_loaded
          integer :: i,chunk_size,irr_z,mi_c,ic_c,wr_i32_size, &
              wr_f32_size,read_f32_size
          real :: logm
          double precision :: j_date,rad_date,solar_dec,rad_time,lda_f

          ! private internal variables
		  !----------------------------------------------------          
          integer :: j,k,ii,jj,kk,mi,ic,ki,k1,kz,mi_wavl, &
              gpu_i,m_row,srftyp,i_mp,j_mp,int_z,sk_z,thr, &
              sk_1,ice_1,ice_z,z_snow,nm_i,block_i,thread_i,zz
          double precision :: ed_w_dir,ed_w_dif,wv_mul,bv_mean, &
              w_temp,zenith_angle,coszen,am_sum,tmp1,tmp2,tmp3, &
              airtemp1c

          double precision, dimension(z_max+z_snow_max) :: &
              dth,     &
              r2st,    &
              g2st,    &
              k2st

          ! code on, dude
		  !----------------------------------------------------          

          chunk_size = 64 ! set OPENMP chunk size
          zz = z_max + z_snow_max
          irr_z = zz + 1  

          lda_f = lda_d  ! init light/snow distribution divider
                  
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

          all_cells_loaded = .FALSE.

          do while (.not. all_cells_loaded)

! ----------------------------------------------------------------------
! filling GPU arrays - decided if cell is valid and prep counters, etc
! ----------------------------------------------------------------------
          gpu_i = 0
          gpu_cells = 0 ! init var to keep track of how many gpu calcs are valid
  
 
          do while (gpu_i .lt. gpu_max_threads)
          if (.not. all_cells_loaded) then

!!$OMP CRITICAL
 
              ic_c = ic_c + 1
              r2st = 0.d0
              g2st = 0.d0
              if (ic_c .gt. ida_n) then
                  ic_c = 1    
                  if (mi_c .eq. tcells) then
                      all_cells_loaded = .TRUE.
                  else    
                      mi_c = mi_c + 1
                      airtemp1c = f(mi_c)%at - kelvin0
                      i_mp = m(mi_c)%x_mp  ! mercator/ncep projection x
                      j_mp = m(mi_c)%y_mp  ! mercator/ncep projection y
                      zenith_angle = acos(sin(solar_dec)*sin(m(mi_c)%lat/360*6.2832) &
                      + cos(solar_dec)*cos(m(mi_c)%lat/360*6.2832)*cos(rad_time) )   ! radians	  
                      if (zenith_angle .lt. 0.d0) then
                          zenith_angle = 0.d0  ! prevent negative light values
                      endif
                      coszen = cos(zenith_angle)

                      !print *,'zenith: ',zenith_angle,'coszen:',coszen

                      ! init absorption, pur arrays -- can this be done somewhere else?
                      do j=1,ida_n
                          ice(j,mi_c)%pur = 0.d0
                          ice(j,mi_c)%PAR_bot = 0.d0
                          ! these below redundantly calculated for every valid category - could improve eff. some
                          ice(j,mi_c)%ed_w = 0.d0
                          ice(j,mi_c)%snow%ed_w = 0.d0
                          ice(j,mi_c)%ed0_nir_dir = 0.d0
                          ice(j,mi_c)%ed0_nir_dif = 0.d0
                          ice(j,mi_c)%Ed0 = 0.
                      enddo
                  endif
              endif

!!$OMP END CRITICAL

              if (.not. all_cells_loaded .and. &
                  m(mi_c)%status .eq. 1 .and. ice(ic_c,mi_c)%z .gt. 0) then
                  
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
                      ice(ic_c,mi_c)%ed0_nir_dir = ice(ic_c,mi_c)%ed0_nir_dir + ed_w_dir
                      ice(ic_c,mi_c)%ed0_nir_dif = ice(ic_c,mi_c)%ed0_nir_dif + ed_w_dif
                      ice(ic_c,mi_c)%Ed0 = ice(ic_c,mi_c)%Ed0 + ed_w_dir + ed_w_dif
                  
                  enddo
                       
                  ! final test to see if this wavelength radiation transfer
                  ! need to be calculated
                  if (ice(ic_c,mi_c)%Ed0 .gt. 0.d0) then

                      block_i = gpu_cells/gpu_tpb ! C-style reference, starting at 0
                      thread_i = gpu_cells - block_i*gpu_tpb ! C-style reference, starting at 0
                      k1 = block_i*f32_in_p_th*gpu_tpb + thread_i + 1 ! start of block&thread data

                      ! increment active gpu cell and gpu memory counters - after since C-style indexes
                      ! start at 0
                      gpu_cells = gpu_cells + 1
                      gpu_i = gpu_i + 1
 
                      sk_z = ice(ic_c,mi_c)%z
                      int_z = sk_z - z_sk
                      sk_1 = int_z + 1
                      z_snow = ice(ic_c,mi_c)%snow%z
                      ice_z = sk_z + z_snow
                      ice_1 = z_snow + 1

                      ! 1nm wavelength index
                      nm_i = ii*10-9

! ----------------------------------------------------------------------
! calc scattering function(g2st) and scattering(r2st)
! ----------------------------------------------------------------------

                      ! estimation particle absorption from chlorophyll and detrital concentrations in ice
                      ! also estimate scattering in ice due to brine volume, which is constant over wavelength
                      do jj=1,sk_z

                          am_sum = c0
                          tmp3 = ice(ic_c,mi_c)%bv(jj)*c_001
                          tmp1 = max(1.d0,ice(ic_c,mi_c)%smalg(jj)/c_chl*tmp3)  ! minimum 1mg/m^3 chla absorption
                          do ii=1,wavl          
                              ! sum algal spectral absorption, used to find detrital absorption, ad in cuda kernel
                              am_sum = am_sum+aph(ii*10-9)*tmp1*10.   ! x10 b/c 10 nm wavelength bins                  
                          enddo
            
                          ! sea ic_ce scattering
                          if (ice(ic_c,mi_c)%t(jj) .le. -22.) then
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
            
                                  tmp1 = max(2.,abs(ice(ic_c,mi_c)%t(jj)))-2.
            !                      r2st(jj+z_snow) = (18.*ice(ic_c,mi_c)%s(jj))*sqrt(tmp1)+15;   ! power law increase
                                  tmp2 = 9.6*ice(ic_c,mi_c)%s(jj)
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
                            (ice(ic_c,mi_c)%smalg(jj) + ice(ic_c,mi_c)%poc(jj)) &
                            *tmp3*5.2632d-3
            
                          ! record ice layer depth
                          dth(jj+z_snow) = ice(ic_c,mi_c)%th(jj)
                          
                      enddo

                      ! deal with snow wavelength independent inherent optical properties
                      if (ice(ic_c,mi_c)%snow%z .gt. 0) then
            
                          ! surface type over ice = snow
                          srftyp = 1
            
                          ii = 1
                          do jj=z_snow,1,-1
                           
                              ! find snow layer scattering (r)
                              if (ice(ic_c,mi_c)%snow%d(jj) .gt. 0.6) then
                                  if ((rs_switch .eq. 0 .and. airtemp1c .lt. 0.) .or. &
                                  (rs_switch .eq. 1 .and. ice(ic_c,mi_c)%snow%t(jj) .le. -1.) .or. &
                                  (rs_switch .eq. 2 .and. ice(ic_c,mi_c)%snow%melt(jj) .lt. 0.5)) then
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
                              dth(ii) = ice(ic_c,mi_c)%snow%th(jj)*snow_fudge
            
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
                      c_f32_in(ki) = real(am_sum)
                      ki = ki + gpu_tpb
                      do jj=1,ice_z			  

                          !ki = k1 + (3*(jj-1)+1)*gpu_tpb
                          if (jj .le. z_snow) then

                              ii = z_snow - jj + 1                      
                              c_f32_in(ki) = real(ice(ic_c,mi_c)%snow%d(ii))
                              ki = ki + 3*gpu_tpb

                          else
                              ii = jj - z_snow                          
                              c_f32_in(ki) = real(ice(ic_c,mi_c)%bv(ii)*c_001)
                              ki = ki + gpu_tpb
                              c_f32_in(ki) = real(ice(ic_c,mi_c)%smalg(ii))
                              ki = ki + gpu_tpb
                              c_f32_in(ki) = real(ice(ic_c,mi_c)%poc(ii))
                              ki = ki + gpu_tpb

                          endif

                      enddo

                      ki = k1 + (3*zz+1)*gpu_tpb
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

                                            
                      ki = k1 +  (6*zz + 2)*gpu_tpb
                      do ii=1,wavl+1
                          
                          if (ii .le. wavl) then

                              c_f32_in(ki) = real(mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii))
                              ki = ki + gpu_tpb                          
                              c_f32_in(ki) = real(mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii))
                              ki = ki + gpu_tpb                          

                          else

                              c_f32_in(ki) = real(ice(ic_c,mi_c)%ed0_nir_dir*(par_to_swd-1.d0))
                              ki = ki + gpu_tpb                          
                              c_f32_in(ki) = real(ice(ic_c,mi_c)%ed0_nir_dir*(par_to_swd-1.d0))
                              ki = ki + gpu_tpb                          

                          endif

                      enddo                              
                          
                      ! record GPU column calc metadata
                      mi_gpu_index(gpu_i) = mi_c*10 + ic_c  ! record value that indicated both mi & ic (cell & category)

                      ki = block_i*i32_in_p_th*gpu_tpb + thread_i + 1 ! start i32_in block&thread data 
                      c_i32_in(ki) = sk_z
                      ki = ki +  gpu_tpb ! start i32_in block&thread data 
                      c_i32_in(ki) = z_snow
                      ki = ki +  gpu_tpb ! start i32_in block&thread data 
                      c_i32_in(ki) = srftyp

                  endif ! end - is there any light?
  
              endif ! end if stats/ice category active check
    
          else

              block_i = gpu_i/gpu_tpb ! C-style reference, starting at 0
              thread_i = gpu_i - block_i*gpu_tpb! C-style reference, starting at 0
              k1 = block_i*i32_in_p_th*gpu_tpb + thread_i + 1 ! start of block&thread data

              gpu_i = gpu_i + 1 ! increment k since all cells have been scanned
              
             ! set number of ice laters to 0, so gpu knows to skip calc
              c_i32_in(k1) = 0    ! nilyr = 0
          
          endif ! end of have all cells be loaded?
          enddo ! end of loading GPU arrays for GPU calc(s)
          
! ----------------------------------------------------------------------
! do GPU calc(s)
! ----------------------------------------------------------------------

  
          if (gpu_cells .gt. 0) then

          wr_i32_size = (gpu_cells + gpu_tpb - 1) / gpu_tpb  ! num valid data blocks
          wr_f32_size = wr_i32_size*f32_in_p_th*gpu_tpb
          wr_i32_size = wr_i32_size*i32_in_p_th*gpu_tpb

          call sia2_edd_gpu_copy_step_mem( &
              c_i32_in,    &
              c_f32_in,    &       
              g_i32_in,    &
              g_f32_in,    &
              wr_i32_size, &
              wr_f32_size)
  
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
                   gpu_cells,   &
                   zz,          &
                   logm)

              call sia2_edd_gpu_return( &
                   g_f32_out,     &
                   c_f32_out,     &
                   f32_out_p_th,  &
                   gpu_tpb,       &
                   gpu_cells)


!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(j,k,ii,jj,kk,mi,ic,ki,k1,kz,mi_wavl, &
!$OMP gpu_i,m_row,srftyp,i_mp,j_mp,int_z,sk_z,thr, &
!$OMP sk_1,ice_1,ice_z,z_snow,nm_i,block_i,thread_i,zz, &
!$OMP ed_w_dir,ed_w_dif,wv_mul,bv_mean, &
!$OMP w_temp,zenith_angle,coszen,am_sum,tmp1,tmp2,tmp3, &
!$OMP airtemp1c, &
!$OMP dth,     &
!$OMP r2st,    &
!$OMP g2st,    &
!$OMP k2st)

!$OMP DO SCHEDULE(DYNAMIC,chunk_size)

              do j=1,gpu_cells

                  ! determine memory indices, bookkeeping
                  ! -----------------------------------------------------
                  !thr = omp_get_thread_num()

                  mi = mi_gpu_index(j)/10
                  ic = mi_gpu_index(j) - mi*10
          
                  sk_z = ice(ic,mi)%z
                  int_z = sk_z - z_sk
                  sk_1 = int_z + 1
                  z_snow = ice(ic,mi)%snow%z
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
                          ice(ic,mi)%snow%ed_w(k) = ice(ic,mi)%snow%ed_w(k) + &
                            dble(c_f32_out(ki))*lda_f
                          ki = ki + gpu_tpb  ! skip array index to ignore PUR in snow - not used 
                      else
                          m_row = k - z_snow  ! find correct ice layer index
                          ice(ic,mi)%ed_w(m_row) =ice(ic,mi)%ed_w(m_row) + &
                            dble(c_f32_out(ki))*lda_f ! store absorbed watts
                          ki = ki + gpu_tpb
                          ice(ic,mi)%pur(m_row,i) = dble(c_f32_out(ki))  ! store PUR
                      endif
                      ki = ki + gpu_tpb
                  enddo                    

                  !if (mi .eq. 27) then
                  !    print *,'PUR: ',ice(ic,mi)%pur(:,i)
                  !    print *,'Ed_W: ',ice(ic,mi)%ed_w
                  !endif

              enddo ! end of recording gpu run(s) data
             
 !$OMP END DO         
 !$OMP END PARALLEL          
          
              print *,'Delta-Eddington solution completed and stored. ',gpu_cells,' total cells.'

          enddo  ! end of snow/light distribution loop
          endif
          
          enddo  ! end of all_cell_loaded test
 
          !print *,'pur: ',ice(1,1)%pur
          !print *,'ed_w: ',ice(1,1)%ed_w
          !print *,'ed_w_snow: ',ice(1,1)%snow%ed_w
   
   
      end SUBROUTINE sia2_edd_gpu_calc
      
