! ==============================================================================
! START CALC LIGHT
! ==============================================================================

! uneeded vars ?
! ed0_temp, Ed,




          ! find denominator for ad estimation from chlorophyll
          ! moved to outside model domain loop - doesn't change
!          ad_denom = 0.
!          do ii=1,wavl
!              ad_denom = ad_denom + exp(-0.008*((ii-1)*10))*10. ! times 10, b/c 10nm wavelength bins
!          enddo

          dth = 0.d0  ! vector assignment

		  ! get number of optical layers
		  z_snow = ice(sc,ic,mi)%snow%z
		  ice_1 = z_snow+1
		  ice_z = z_snow + ice(sc,ic,mi)%z

          Ed_IR(1) = 0.d0  ! zero IR band - it will be added to incrementally later in Watts
          Ed0_nir_dir = 0.d0
          Ed0_nir_dif = 0.d0
          PUR_dble(:,kk) = 0.d0
		  ice(sc,ic,mi)%PAR_bot = 0.d0

          albodr = 0.d0
          albodf = 0.d0

		  ! estimation particle absorption from chlorophyll and detrital concentrations in ice
		  ! also estimate scattering in ice due to brine volume, which is constant over wavelength
		  do jj=1,sk_z
			  am_sum = 0.d0
			  tmp3 = ice(sc,ic,mi)%bv(jj)*c_001
		      tmp1 = max(1.0,ice(sc,ic,mi)%smalg(jj)/c_chl*tmp3)  ! minimum 1mg/m^3 chla absorption
			  do ii=1,wavl
				  ! algal spectral absorption
				  am(jj,ii) = aph(ii*10-9)*tmp1
				  am_sum = am_sum+am(jj,ii)*10.d0   ! x10 b/c 10 nm wavelength bins
			  enddo
			  ad(jj,1) = ice(sc,ic,mi)%poc(jj)/ice(sc,ic,mi)%smalg(jj)*am_sum/ad_denom ! using fractional amount of poc compared to live algae
			  do ii=2,wavl
				  ad(jj,ii) = ad(jj,1)*exp(-0.008d0*((dble(ii)-1.d0)*10.d0))
			  enddo

              ! sea ice scattering
			  if (ice(sc,ic,mi)%t(jj) .le. -22.) then
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

                      tmp1 = max(2.d0,abs(ice(sc,ic,mi)%t(jj)))-2.d0
!                      r2st(jj+z_snow) = (18.*ice(sc,ic,mi)%s(jj))*sqrt(tmp1)+15;   ! power law increase
                      tmp2 = 9.6d0*ice(sc,ic,mi)%s(jj)
                      r2st(jj+z_snow) = min(tmp2*8.,tmp2*tmp1)+15.d0;   ! linear low-temp increase
				      g2st(jj+z_snow) = 0.98d0

!                  endif


!                      if (ice(sc,ic,mi)%t(jj) .le. -10.) then
                          ! mirabilite crystal scattering regime
!                          r2st(jj+z_snow) = 300.-100.*(ice(sc,ic,mi)%t(jj)+10.)/12.
!                          r2st(jj+z_snow) = 200.
!                          g2st(jj+z_snow) = 0.98
!                      else
                          ! brine pocket scattering regime
!                          tmp2 = min(-1.8,ice(sc,ic,mi)%t(jj))
                          !r2st(jj+z_snow) = -971.1/(tmp2)**2+309.7
!                          r2st(jj+z_snow) = -636.5/abs(tmp2)+363.6240
                          !r2st(jj+z_snow) = -307./abs(tmp2)+181.
                          !r2st(jj+z_snow) = max(r2st(jj+z_snow),300.)
!                          g2st(jj+z_snow) = 0.98
!                      endif
 !                 endif
              endif

              ! algal scattering - from Babin et al. 2003, scattering = 1.0m^2/gDryMass* mgC/m^3 * 1gC/1000mgC * 1gDryMass/0.19gC (Sicko-Goad et al. 1984)
              r2st(jj+z_snow) = r2st(jj+z_snow) + (ice(sc,ic,mi)%smalg(jj) + ice(sc,ic,mi)%poc(jj)) &
                  *tmp3*5.2632d-3

              ! record ice layer depth
              dth(jj+z_snow) = ice(sc,ic,mi)%th(jj)

		  enddo

          ! deal with snow wavelength independent inherent optical properties
          if (ice(sc,ic,mi)%snow%z .gt. 0) then

              ! surface type over ice = snow
              srftyp = 1

              ii = 1
              do jj=z_snow,1,-1

				  ! find snow layer scattering (r)
                  if (ice(sc,ic,mi)%snow%d(jj) .lt. 6.d5) then
                      if ((rs_switch .eq. 0 .and. airtemp1c .lt. 0.d0) .or. &
                      (rs_switch .eq. 1 .and. ice(sc,ic,mi)%snow%t(jj) .le. -1.) .or. &
                      (rs_switch .eq. 2 .and. ice(sc,ic,mi)%snow%melt(jj) .lt. 0.5d0)) then
                          r2st(ii) = 3000.d0
                          g2st(ii) = 0.89d0
                      else
                          r2st(ii) = 900.d0
                          g2st(ii) = 0.94d0
                          !r2st(ii) = 2000
                          !g2st(ii) = 0.89
                      endif
                  else
                      ! surface scattering layer
                      r2st(ii) = 900.d0
                      g2st(ii) = 0.94
                  endif

				  ! record effective snow layer attenuation depth
				  dth(ii) = ice(sc,ic,mi)%snow%th(jj)*snow_fudge*lda_multiplier(kk)

                  ii = ii + 1
              enddo

          else

              ! surface type over ice = air
              srftyp = 1

          endif

          ! calculate column optical attenuation - main spectral loop
          do ii=1,wavl+1

			  ! wavelength index
			  nm_i = ii*10-9

			  ! wavelength multiplier - 400 and 700 nm are valued at ~half the others, to make a total of 301
			  if (ii .eq. 1 .or. ii .eq. wavl) then
			      wv_mul = 5.5d0
			  else
			      wv_mul = 10.d0
			  endif

              if (ii .le. wavl) then

				  ! find total downwelling irradiance
				  Ed(1,ii) = mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii) + &
					  mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii)

				  ! find total irradiance in watts from PAR (400-700nm), for use later in finding incident near-IR
				  Ed0_nir_dir = Ed0_nir_dir + &
				    mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii)*quanta_to_watts(ii)*wv_mul  ! Watts/m^2
				  Ed0_nir_dif = Ed0_nir_dif + &
				    mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii)*quanta_to_watts(ii)*wv_mul  ! Watts/m^2

                  ! record surface downwelling for integration below
                  Ed0_temp = Ed(1,ii)

              else

				  ! convert PAR watts to NIR watts
				  Ed0_nir_dir = Ed0_nir_dir*(par_to_swd-1.)
				  Ed0_nir_dif = Ed0_nir_dif*(par_to_swd-1.)

                  ! record surface downwelling for integration below
                  Ed0_temp = Ed0_nir_dir + Ed0_nir_dif

              endif

			  do jj=1,ice_z

                  ! case for snow absorption
                  if (jj .lt. ice_1) then
                      jjj = z_snow - jj + 1

                      if (ii .le. wavl) then
						  ! find layer absorption (k)
						  k2st(jj) = aice(nm_i)*ice(sc,ic,mi)%snow%d(jjj)/IceD + &
						      a_factor    ! account for ice part only
                      else
                          ! NIR absorption
                          k2st(jj) = a_ice_ir*ice(sc,ic,mi)%snow%d(jjj)/IceD + &
                              a_factor! account for ice part only
                      endif

                  ! case for ice absortion
                  else

		              m_row = jj-z_snow

                      if (ii .le. wavl) then

		                  bv_mean = ice(sc,ic,mi)%bv(m_row)*c_001

						  ! Fresh Ice abosrption ...
						  !k2st(jj) = 0.92 - 3.522e-3*lambda(ii) + 3.632e-6*lambda(ii)**2 ! pure ice - from arrigo
						  k2st(jj) = (aice(nm_i)*(1. - bv_mean) + bv_mean* &
						  (awater(nm_i) - awater_sc(nm_i)*ice(sc,ic,mi)%s(m_row) - &
						  (0. - 22.)*awater_tc(nm_i)))

						  ! phyoplankton & detrital absorption
						  k2st(jj) = k2st(jj) + (am(m_row,ii)+ad(m_row,ii))

                      else

                          ! NIR absorption of ice + particle/chla absorption at 700nm
                          k2st(jj) = a_ice_ir + (am(m_row,wavl)+ad(m_row,wavl))

                      endif

                  endif

                  ! calculate omega & tau IOPs needed for delta-Eddington attenuation
                  w2st(jj) = r2st(jj)/(r2st(jj)+k2st(jj))
                  tau2st(jj) = (r2st(jj)+k2st(jj))*dth(jj)

              enddo

              ! moved to outside loop
              !albodr = 0.
              !albodf = 0.

              call sia2_env_edd_solution(coszen, srftyp, tau2st, w2st, g2st, expt, &
              albodr, albodf, trndir, trntdr, trndif, rupdir, rupdif, &
              rdndif, z_snow, ice(sc,ic,mi)%z)

                ! the interface reflectivities and transmissivities required
                ! to evaluate interface fluxes are returned from solution_dEdd;
                ! now compute up and down fluxes for each interface, using the
                ! combined layer properties at each interface:
                !
                !              layers       interface
                !
                !       ---------------------  k
                !                 k
                !       ---------------------

                 ! if (ii .eq. 31 .and. mi .eq. 12 .and. kk .eq. 1) then
                 !   print *,'BEFORE_____________________________________'
                 !   print *,'g2st: ',g2st
                 !   print *,'r2st: ',r2st
                 !   print *,'k2st: ',k2st
                 !   print *,'After_____________________________________'
                 !   print *,'g: ',g2st
                 !   print *,'w0: ',w2st
                 !   print *,'tau: ',tau2st
                 ! endif

                do jj=1,ice_z+1
                    ! interface scattering
                    refk          = 1./(1. - rdndif(jj)*rupdif(jj))
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
                                (fdirdn(jj-1) + fdirup(jj)) * &
                                Ed0_nir_dir + &
                                (fdifdn(jj-1) + fdifup(jj)) * &
                                Ed0_nir_dif
                            layer_out = &
                                (fdirup(jj-1) + fdirdn(jj)) * &
                                Ed0_nir_dir + &
                                (fdifup(jj-1) + fdifdn(jj)) * &
                                Ed0_nir_dif

                        endif

                        absorp = abs(layer_in - layer_out)  ! units depend on whether this is NIR or a PAR wavelenth
                        absorp = max(0.,absorp)

                        if ((jj-1) .lt. ice_1) then

                            ! index to snow layer -- ed_w_snow is stored with top of snow in 1st array position
                            m_row = jj-1

                            ! record absorption for wavelength in snow
                            tmp1 = absorp*lda_f  ! Watts/m^2 or µEin/m^2/s

                            if (ii .le. wavl) then
                                tmp1 = tmp1*quanta_to_watts(ii)*wv_mul  ! Watts/m^2
                            endif

                            Ed_W_snow(m_row) = Ed_W_snow(m_row) + tmp1   ! Watts/m^2

                        else

                            ! index to ice layer
                            m_row = jj-1-z_snow

                            if (ii .le. wavl) then

                                ! record absorption for wavelength in ice
                                tmp1 = absorp*quanta_to_watts(ii)*wv_mul  ! Watts/m^2
                                Ed_W_mean(m_row) = Ed_W_mean(m_row) + tmp1*lda_f ! Watts/m^2

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

                                ! find downwelling irradiance at top and bottom of layer
                                !par_top = &
                                !  ((fdirdn(jj-1)) &
                                !  *mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii) + &
                                !   (fdifdn(jj-1)) &
                                !  *mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii))
                                !par_bot = &
                                !  ((fdirdn(jj)) &
                                !  *mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii) + &
                                !   (fdifdn(jj)) &
                                !  *mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii))

                                ! calc mean layer light using simple average
                                !par_mid = (par_top+par_bot)*0.5

                                ! calc mean layer light using exponential curve
                                if (par_bot .eq. 0.d0) then
                                    par_mid = (par_top+par_bot)*c_5
                                else
                                    tmp_k = log(par_bot/par_top)/dth(jj-1)
                                    par_mid = (par_bot - par_top)/tmp_k/dth(jj-1)
                                endif

                                ! record Ed at middle of the layer, just in case other parts of
                                ! the model use it?
                                Ed(m_row+1,ii) = par_mid

                                ! find wavelength contribution to PUR
                                PUR_dble(m_row,kk) = PUR_dble(m_row,kk) &
                                  + par_mid*aph(nm_i)/aph_max*wv_mul  ! times wv_mul, b/c 10nm wavelength bins

                              ! record bottom PAR
                              if (jj .eq. ice_z+1) then
                                  ice(sc,ic,mi)%PAR_bot = ice(sc,ic,mi)%PAR_bot + par_bot*wv_mul*lda_f
                              endif

                            else

                                ! record NIR absorption for heat conduction in ice
                                Ed_W_mean(m_row) = Ed_W_mean(m_row) + &
                                   absorp*lda_f

                            endif

                        endif

                    endif

                enddo  ! end of layer loop

            enddo  ! end of spectral loop

! ==============================================================================
! END CALC LIGHT
! ==============================================================================
