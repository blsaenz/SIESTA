

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
          integer :: i,chunk_size,irr_z,mi_c,ic_c,wavl_c,wr_i32_size, &
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
              absorp,refk,par_top,par_bot,par_mid,layer_in,layer_out, &
              tmp_k,airtemp1c, &
              rintfc   , & ! reflection (multiple) at an interface
              refkp1   , & ! interface multiple scattering for k+1
              refkm1   , & ! interface multiple scattering for k-1
              tdrrdir  , & ! direct tran times layer direct ref 
              tdndif       ! total down diffuse = tot tran - direct tran

          double precision, dimension(z_max+z_snow_max) :: &
              dth,     &
              r2st,    &
              g2st,    &
              k2st,    &
              rdir,    & ! layer reflectivity to direct radiation
              rdif_a,  & ! layer reflectivity to diffuse radiation from above
              rdif_b,  & ! layer reflectivity to diffuse radiation from below
              tdir,    & ! layer transmission to direct radiation (solar beam + diffuse)
              tdif_a,  & ! layer transmission to diffuse radiation from above
              tdif_b,  & ! layer transmission to diffuse radiation from below
              trnlay,  & ! solar beam transm for layer (direct beam only)
              fdirup,  & ! up   flux at model interface due to direct beam at top surface
              fdirdn,  & ! down flux at model interface due to direct beam at top surface
              fdifup,  & ! up   flux at model interface due to diffuse beam at top surface
              fdifdn     ! down flux at model interface due to diffuse beam at top surface

          double precision, dimension(z_max+z_snow_max+1) :: &
              trndir,  &
              trntdr,  &
              trndif,  &
              rdndif,   &
              rupdir,  &
              rupdif

          double precision, dimension(z_max,wavl) :: &
              am,ad 

          ! code on, dude
		  !----------------------------------------------------          

          chunk_size = 64 ! set OPENMP chunk size
          irr_z = z_max + z_snow_max + 1  

          lda_f = lda_d  ! init light/snow distribution divider
                  
          mi_c = 0 ! initalize mi for counting below
          ic_c = 100 ! inintalize current ice category (high)
          wavl_c = 1000 ! initialize current wavelength (high)

          mi = 0  ! to catch bugs
          ic = 0  ! to catch bugs
          
          zz = z_max + z_snow_max

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
 
              ! increment cell (mi), ice category (ic), and wavlength (mi_wavl)
              ! and update cell-specific parameters and indices
              wavl_c = wavl_c + 1  
              if (wavl_c .gt. (wavl+1)) then
                  wavl_c = 1
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
              endif                      

!!$OMP END CRITICAL

              if (.not. all_cells_loaded .and. &
                  m(mi_c)%status .eq. 1 .and. ice(ic_c,mi_c)%z .gt. 0) then
                  
                  ! Find Total shortwave irradiance in Watts/m^2
                  if (wavl_c .le. wavl) then

                      ! wavelength multiplier - 400 and 700 nm are valued at ~half the others, to make a total of 301
                      if (wavl_c .eq. 1 .or. wavl_c .eq. wavl) then
                          wv_mul = 5.5d0
                      else
                          wv_mul = 10.d0
                      endif
                      
                      ! µEin -> watts conversion for wavelength
                      w_temp = (7.e-7*lambda(wavl_c)**2 - 0.0012*lambda(wavl_c) + &
                               0.6416)*wv_mul  ! Watts / (µEin/m/s)            
                      
                      ! surface irradiance in Watts
                      ed_w_dir = mp_f%Edir(f_index)%Ed(i_mp,j_mp,wavl_c)*w_temp                      
                      ed_w_dif = mp_f%Edif(f_index)%Ed(i_mp,j_mp,wavl_c)*w_temp

                      ! find Ed0, NIR watts
                      ice(ic_c,mi_c)%ed0_nir_dir = ice(ic_c,mi_c)%ed0_nir_dir + ed_w_dir
                      ice(ic_c,mi_c)%ed0_nir_dif = ice(ic_c,mi_c)%ed0_nir_dif + ed_w_dif
                      ice(ic_c,mi_c)%Ed0 = ice(ic_c,mi_c)%Ed0 + ed_w_dir + ed_w_dif
                  
                  endif
                       
                  ! final test to see if this wavelength radiation transfer
                  ! need to be calculated
                  if ((ed_w_dir+ed_w_dif) .gt. 0.d0 .or. &
                      ((ice(ic_c,mi_c)%ed0_nir_dir+ice(ic_c,mi_c)%ed0_nir_dif) .gt. 0. .and. &
                       wavl_c .eq. wavl+1)) then

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
                      nm_i = wavl_c*10-9

! ----------------------------------------------------------------------
! calc absporption and scattering, if it hasn't been done yet
! ----------------------------------------------------------------------
                      if (r2st(1) .eq. 0.d0) then 

                          ! estimation particle absorption from chlorophyll and detrital concentrations in ice
                          ! also estimate scattering in ice due to brine volume, which is constant over wavelength
                          do jj=1,sk_z
                              am_sum = 0.d0
                              tmp3 = ice(ic_c,mi_c)%bv(jj)*c_001
                              tmp1 = max(1.0,ice(ic_c,mi_c)%smalg(jj)/c_chl*tmp3)  ! minimum 1mg/m^3 chla absorption
                              do ii=1,wavl          
                                  ! algal spectral absorption
                                  am(jj,ii) = aph(ii*10-9)*tmp1
                                  am_sum = am_sum+am(jj,ii)*10.   ! x10 b/c 10 nm wavelength bins                  
                              enddo
                              ad(jj,1) = ice(ic_c,mi_c)%poc(jj)/ice(ic_c,mi_c)%smalg(jj)* &
                                  am_sum/ad_denom ! using fractional amount of poc compared to live algae              
                              do ii=2,wavl
                                  ad(jj,ii) = ad(jj,1)*exp(-0.008*((dble(ii)-1.)*10.))              
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

                      endif ! end of has scattering/absorption been calculated yet for this category (am_sum == 0)?

! ----------------------------------------------------------------------
! calc apparent optical properties, store in GPU arrays
! ----------------------------------------------------------------------

                      k2st = 0.d0
                      do jj=1,ice_z			  
        
                          ! case for snow absorption
                          if (jj .lt. ice_1) then
                              ii = z_snow - jj + 1
    
                              if (wavl_c .le. wavl) then	
                                  ! find layer absorption (k)
                                  k2st(jj) = aice(nm_i)*ice(ic_c,mi_c)%snow%d(ii)/IceD + &
                                      a_factor    ! account for ice part only	
                              else                          
                                  ! NIR absorption
                                  k2st(jj) = a_ice_ir*ice(ic_c,mi_c)%snow%d(ii)/IceD + &
                                      a_factor! account for ice part only	
                              endif
        
                          ! case for ice absortion
                          else
        
                              m_row = jj-z_snow
        
                              if (wavl_c .le. wavl) then		 
                 
                                  bv_mean = ice(ic_c,mi_c)%bv(m_row)*c_001
                
                                  ! Fresh Ice abosrption ...
                                  !k2st(jj) = 0.92 - 3.522e-3*lambda(wavl_c) + 3.632e-6*lambda(wavl_c)**2 ! pure ice - from arrigo
                                  k2st(jj) = (aice(nm_i)*(1. - bv_mean) + bv_mean* &
                                  (awater(nm_i) - awater_sc(nm_i)*ice(ic_c,mi_c)%s(m_row) - &
                                  (0. - 22.)*awater_tc(nm_i)))                 
            
                                  ! phyoplankton & detrital absorption
                                  k2st(jj) = k2st(jj) + (am(m_row,wavl_c)+ad(m_row,wavl_c))
        
                              else
                              
                                  ! NIR absorption of ice + particle/chla absorption at 700nm
                                  k2st(jj) = a_ice_ir + (am(m_row,wavl)+ad(m_row,wavl))
        
                              endif
        
                          endif
                          
                          ! record AOPs and metadata needed for 
                          ! delta-Eddington attenuation GPU calc
                          ki = k1 + ((jj-1)*3+1)*gpu_tpb  ! memory index calc
                          c_f32_in(ki) = real(g2st(jj))
                          ki = ki + gpu_tpb                       
                          tmp1 = r2st(jj)/(r2st(jj)+k2st(jj))
                          c_f32_in(ki) = real(r2st(jj)/(r2st(jj)+k2st(jj)))
                          if (abs(dth(jj)) .lt. 0.0001 .or. abs(dth(jj)) .gt. 10.) then
                              print *,'t: ',ice(ic,mi)%t
                              print *,'s: ',ice(ic,mi)%s
                              print *,'bv: ',ice(ic,mi)%bv
                              print *,'th: ',ice(ic,mi)%bv
                              print *,'smalg: ',ice(ic,mi)%smalg
                              print *,'absorp: ',k2st
                              print *,'scat: ',r2st
                              print *,'status: ',m(mi_c)%status
                              print *,'ice_z: ',ice_z
                              print *,'sk_z: ',sk_z
                              print *,'wavl_c: ',wavl_c
                              print *,'jj: ',jj
                              bv_mean = -1
                              tmp1 = sqrt(bv_mean)
                          endif
                         ki = ki + gpu_tpb                       
                          c_f32_in(ki) = real((r2st(jj)+k2st(jj))*dth(jj))

                          !print '(a15,i6,i4,i4,i3,i6)','gi bi thi jj ki ',gpu_i,block_i,thread_i,jj,ki

                      enddo

                      ! record GPU column calc metadata
                      mi_gpu_index(gpu_i) = mi_c*1000 + 10*wavl_c + ic_c  ! record value that indicated both mi & ic (cell & category)
                      
                      c_f32_in(k1) = real(coszen)

                      k1 = block_i*i32_in_p_th*gpu_tpb + thread_i + 1 ! start i32_in block&thread data 
                      c_i32_in(k1) = sk_z
                      k1 = k1 +  gpu_tpb ! start i32_in block&thread data 
                      c_i32_in(k1) = z_snow
                      k1 = k1 +  gpu_tpb ! start i32_in block&thread data 
                      c_i32_in(k1) = srftyp
                                            
                  endif ! end of 'is this wavl-ic-mi combo valid for calculation?'
  
                  !if (wavl_c .eq. 31 .and. mi_c .eq. 12) then
                  !  print *,'BEFORE_____________________________________'
                  !  print *,'g2st: ',g2st
                  !  print *,'r2st: ',r2st
                  !  print *,'k2st: ',k2st
                  !endif
  
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

!                  if (ii .eq. 32 .and. mi .eq. 27) then
!                      print *,'pur','--',ice(1,27)%pur
!                      print *,'ed_w','--',ice(1,27)%ed_w
!                      print *,'ed_w_snow','--',ice(1,27)%snow%ed_w
!                  endif

  
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
                   g_f32_out,   &
                   c_f32_out,   &      
                   i32_in_p_th, &
                   f32_in_p_th, &
                   f32_out_p_th,&
                   gpu_tpb,     &
                   gpu_cells,   &
                   logm)


!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(j,k,ii,jj,kk,mi,ic,ki,k1,kz,mi_wavl, &
!$OMP gpu_i,m_row,srftyp,i_mp,j_mp,int_z,sk_z,thr, &
!$OMP sk_1,ice_1,ice_z,z_snow,nm_i,block_i,thread_i,zz, &
!$OMP ed_w_dir,ed_w_dif,wv_mul,bv_mean, &
!$OMP w_temp,zenith_angle,coszen,am_sum,tmp1,tmp2,tmp3, &
!$OMP absorp,refk,par_top,par_bot,par_mid,layer_in,layer_out, &
!$OMP tmp_k,airtemp1c, &
!$OMP rintfc, & 
!$OMP refkp1, & 
!$OMP refkm1, & 
!$OMP tdrrdir, & 
!$OMP tdndif,   &  
!$OMP dth,     &
!$OMP r2st,    &
!$OMP g2st,    &
!$OMP k2st,    &
!$OMP rdir,    & 
!$OMP rdif_a,  & 
!$OMP rdif_b,  & 
!$OMP tdir,    & 
!$OMP tdif_a,  & 
!$OMP tdif_b,  & 
!$OMP trnlay,  & 
!$OMP fdirup,  & 
!$OMP fdirdn,  & 
!$OMP fdifup,  &
!$OMP fdifdn,   &
!$OMP trndir,  &
!$OMP trntdr,  &
!$OMP trndif,  &
!$OMP rdndif,   &
!$OMP rupdir,  &
!$OMP rupdif, &
!$OMP am,ad)

!$OMP DO SCHEDULE(DYNAMIC,chunk_size)

              do j=1,gpu_cells

                  ! determine memory indices, bookkeeping
                  ! -----------------------------------------------------
                  thr = omp_get_thread_num()

                  mi = mi_gpu_index(j)/1000
                  ii = (mi_gpu_index(j) - mi*1000)/10
                  ic = mi_gpu_index(j) - mi*1000 - ii*10

                  i_mp = m(mi)%x_mp  ! mercator/ncep projection x
                  j_mp = m(mi)%y_mp  ! mercator/ncep projection y

                  ! wavelength multiplier - 400 and 700 nm are valued at ~half the others, to make a total of 301
                  if (ii .eq. 1 .or. ii .eq. wavl) then
                      wv_mul = 5.5d0
                  else
                      wv_mul = 10.d0
                  endif

                  ! 1nm wavelength index
                  nm_i = ii*10-9
          
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
                      rdir(k) = dble(c_f32_out(ki))
                      ki = ki + gpu_tpb
                      rdif_a(k) = dble(c_f32_out(ki))
                      ki = ki + gpu_tpb
                      rdif_b(k) = dble(c_f32_out(ki))
                      ki = ki + gpu_tpb
                      tdir(k) = dble(c_f32_out(ki))
                      ki = ki + gpu_tpb
                      tdif_a(k) = dble(c_f32_out(ki))
                      ki = ki + gpu_tpb
                      tdif_b(k) = dble(c_f32_out(ki))
                      ki = ki + gpu_tpb
                      trnlay(k) = dble(c_f32_out(ki))
                      ki = ki + gpu_tpb

                         kk = block_i*f32_in_p_th*gpu_tpb + thread_i + 1 + ((k-1)*3+1)*gpu_tpb  ! memory index calc
                          g2st(k) = c_f32_in(kk)
                          kk = kk + gpu_tpb                       
                          r2st(k) = c_f32_in(kk)
                          kk = kk + gpu_tpb                       
                          k2st(k) = c_f32_in(kk)

                  enddo

                  !if (ii .eq. 31 .and. mi .eq. 12) then
                  !  print *,'After_____________________________________'
                  !  print *,'g: ',g2st
                  !  print *,'w0: ',r2st
                  !  print *,'tau: ',k2st
                  !endif

                  !if (ii .eq. 27) then
                  !if (block_i .eq. 3 .and. thread_i .eq. 93) then
                  !    print *,'rdir',mi,ii,'--',rdir(1:ice_z)
                  !    print *,'rdif_a',mi,ii,'--',rdif_a(1:ice_z)
                  !    print *,'rdif_b',mi,ii,'--',rdif_b(1:ice_z)
                  !    print *,'tdir',mi,ii,'--',tdir(1:ice_z)
                  !    print *,'tdif_a',mi,ii,'--',tdif_a(1:ice_z)
                  !    print *,'tdif_b',mi,ii,'--',tdif_b(1:ice_z)
                  !    print *,'trnlay',mi,ii,'--',trnlay(1:ice_z)

                  !    print *,'g',mi,ii,'--',g2st(1:ice_z)
                  !    print *,'w0',mi,ii,'--',r2st(1:ice_z)
                  !    print *,'tau',mi,ii,'--',k2st(1:ice_z)
                  !endif

                  ! initialize top interface of top layer 
                   trndir(1) =   1.d0
                   trntdr(1) =   1.d0
                   trndif(1) =   1.d0
                   rdndif(1) =   0.d0
                  
                  ! find transmission & reflection
                  ! -----------------------------------------------------
                  do k=2,ice_z+1
                      trndir(k) = trndir(k-1)*trnlay(k-1)
                      refkm1        = 1.d0/(1.d0 - rdndif(k-1)*rdif_a(k-1))
                      tdrrdir       = trndir(k-1)*rdir(k-1)
                      tdndif        = trntdr(k-1) - trndir(k-1)
                      trntdr(k) = trndir(k-1)*tdir(k-1) + &
                       (tdndif + tdrrdir*rdndif(k-1))*refkm1*tdif_a(k-1)
                      rdndif(k) = rdif_b(k-1) + &
                        (tdif_b(k-1)*rdndif(k-1)*refkm1*tdif_a(k-1))
                      trndif(k) = trndif(k-1)*refkm1*tdif_a(k-1)
                  enddo 
                  rupdir(ice_z+1) = 0.d0
                  rupdif(ice_z+1) = 0.d0
                  do k=ice_z,1,-1
                      refkp1        = 1.d0/( 1.d0 - rdif_b(k)*rupdif(k+1))
                      rupdir(k) = rdir(k) + &
                        ( trnlay(k)*rupdir(k+1) + &
                         (tdir(k)-trnlay(k))*rupdif(k+1) ) * &
                          refkp1*tdif_b(k)
                      rupdif(k) = rdif_a(k) + &
                          tdif_a(k)*rupdif(k+1)* &
                          refkp1*tdif_b(k)
                  enddo  

!                if (ii .eq. 27) then
!                print *,mi,'trndir: ',trndir
!                endif
 
                do jj=1,ice_z+1
                   
                    ! interface scattering
                    refk          = 1.d0/(1.d0 - rdndif(jj)*rupdif(jj))
                    ! dir tran ref from below times interface scattering, plus diff
                    ! tran and ref from below times interface scattering
                    fdirup(jj) = (trndir(jj)*rupdir(jj) + &
                                    (trntdr(jj)-trndir(jj)) &  
                                    *rupdif(jj))*refk
                    ! dir tran plus total diff trans times interface scattering plus
                    ! dir tran with up dir ref and down dif ref times interface scattering 
                    fdirdn(jj) = trndir(jj) + (trntdr(jj)  &
                                  - trndir(jj) + trndir(jj)  & 
                                  *rupdir(jj)*rdndif(jj))*refk
                    ! diffuse tran ref from below times interface scattering
                    fdifup(jj) = trndif(jj)*rupdif(jj)*refk
                    ! diffuse tran times interface scattering
                    fdifdn(jj) = trndif(jj)*refk
        
                    if (jj .gt. 1) then

                        if (ii .le. wavl) then
                        
                            ! find layer absorption
                            layer_in = &
                                (fdirdn(jj-1) + fdirup(jj)) * &
                                mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii) + &
                                (fdifdn(jj-1) + fdifup(jj)) * &
                                mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii)
                            
                            layer_out = &
                                (fdirup(jj-1) + fdirdn(jj)) * &
                                mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii) + &
                                (fdifup(jj-1) + fdifdn(jj)) * &
                                mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii)
                        
                        else

                            ! find NIR layer absorption
                            layer_in = &
                                ((fdirdn(jj-1) + fdirup(jj)) * &
                                ice(ic,mi)%ed0_nir_dir + &
                                (fdifdn(jj-1) + fdifup(jj)) * &
                                ice(ic,mi)%ed0_nir_dif) * &
                                (par_to_swd-1.d0)
                            layer_out = &
                                ((fdirup(jj-1) + fdirdn(jj)) * &
                                ice(ic,mi)%ed0_nir_dir + &
                                (fdifup(jj-1) + fdifdn(jj)) * &
                                ice(ic,mi)%ed0_nir_dif) * &
                                (par_to_swd-1.d0)

                        endif

                        absorp = abs(layer_in - layer_out)  ! units depend on whether this is NIR or a PAR wavelenth
                        absorp = max(0.,absorp)

                        if ((jj-1) .lt. ice_1) then
    
                            ! index to snow layer -- ed_w_snow is stored top of snow in 1st array position
                            m_row = jj-1

                            ! record absorption for wavelength in snow
                            tmp1 = absorp*lda_f  ! Watts/m^2 or µEin/m^2/s
                            
                            if (ii .le. wavl) then
                                tmp1 = tmp1*quanta_to_watts(ii)*wv_mul  ! Watts/m^2
                            endif
    
                            ice(ic,mi)%snow%ed_w(m_row) = &
                                ice(ic,mi)%snow%ed_w(m_row) + tmp1   ! Watts/m^2 

                        else
    
                            ! index to ice layer
                            m_row = jj-1-z_snow
                       
                            if (ii .le. wavl) then
    
                                ! record absorption for wavelength in ice
                                tmp1 = absorp*quanta_to_watts(ii)*wv_mul  ! Watts/m^2                                
                                !if (mi .eq. 27) then
                                !    print *,'thread: ',thr,' wavl: ',ii,' m_row: ',m_row
                                !endif
                                ice(ic,mi)%ed_w(m_row) = ice(ic,mi)%ed_w(m_row) + tmp1*lda_f ! Watts/m^2	
        
                                ! find total irradiance at top and bottom of layer
                                par_top = &
                                  ((fdirdn(jj-1)+fdirup(jj-1)) &
                                  *mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii) + &
                                   (fdifdn(jj-1)+fdifup(jj-1)) &
                                  *mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii))
                                par_bot = &
                                  ((fdirdn(jj)+fdirup(jj)) &
                                  *mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii) + &
                                   (fdifdn(jj)+fdifup(jj)) &
                                  *mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii))

                                ! calc mean layer light using simple average
                                !par_mid = (par_top+par_bot)*0.5
                
                                !calc mean layer light using exponential curve
                                if (par_bot .eq. 0.) then                                
                                    par_mid = (par_top+par_bot)*c_5
                                else  
                                    tmp_k = log(par_bot/par_top)/ice(ic,mi)%th(m_row)  
                                    par_mid = (par_bot - par_top)/tmp_k/ice(ic,mi)%th(m_row)
                                endif

                                !print *,'par_mid ',m_row,' -- ',par_mid
                                      
                                ! find wavelength contribution to PUR
                                ice(ic,mi)%pur(m_row,i) = ice(ic,mi)%pur(m_row,i) &
                                  + par_mid*aph(nm_i)/aph_max*wv_mul  ! times wv_mul, b/c 10nm wavelength bins

                                !print *,'par_mid: ',par_mid,aph(nm_i),wv_mul

                              ! record bottom PAR
                              if (jj .eq. ice_z+1) then
                                  ice(ic,mi)%PAR_bot = ice(ic,mi)%PAR_bot + par_bot*wv_mul*lda_f
                              endif
     
                            else
                          
                                ! record NIR absorption for heat conduction in ice
                                ice(ic,mi)%ed_w(m_row) = ice(ic,mi)%ed_w(m_row) + &
                                   absorp*lda_f 
                                !if (mi .eq. 27) then
                                !    print *,'thread: ',thr,' wavl: ',ii,' m_row: ',m_row
                                !endif

                            endif
                         
                        endif
 
                    endif 
                                
                enddo  ! end of layer loop


                  !if (ii .eq. 32 .and. mi .eq. 27) then
                  !    print *,'pur',mi,'--',ice(ic,mi)%pur(:,i)
                  !    print *,'ed_w',mi,'--',ice(ic,mi)%ed_w
                  !   print *,'ed_w_snow',mi,'--',ice(ic,mi)%snow%ed_w
                  !endif

                  !if (ii .eq. 32 .and. mi .eq. 12) then
                  !    print *,'i -- ',ice(ic,mi)%ed_w
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
      
