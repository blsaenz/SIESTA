! ==============================================================================
! START CALC LIGHT
! req'd variables:
! kk = sda number
! many others....
! ==============================================================================

          ! find denominator for ad estimation from chlorophyll
          ! moved to outside model domain loop - doesn't change
!          ad_denom = 0.
!          do ii=1,wavl
!              ad_denom = ad_denom + exp(-0.008*((ii-1)*10))*10. ! times 10, b/c 10nm wavelength bins
!          enddo

          Ed_IR(1) = 0.

          ! if there is snow, specular reflection doesn't really happen
          if (ice(ic,mi)%snow%z .gt. 0) then
              tmp1 = rsky0
          else
              tmp1 = rsun0
          endif

          ! surface specular reflection and near-IR calculation
          do ii=1,wavl

              ! find total downwelling irradiance right under ice/snow surface (i.e. after reflection)
              Ed(1,ii) = (1.-tmp1)*mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii) + &
                  (1.-rsky0)*mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii)

              ! correct for clouds if using clear-sky par
              Ed(1,ii) = Ed(1,ii)*par_cf_f

              ! find total irradiance in watts from PAR (400-700nm), for use later in finding incident near-IR
              Ed_IR(1) = Ed_IR(1) + Ed(1,ii)*quanta_to_watts(ii)*10.  ! Watts/m^2

          enddo

          ! convert PAR watts to near-IR watts
          Ed_IR(1) = Ed_IR(1)*(par_to_swd-1.)

          ! attenuate PAR and near-IR through snow layers
          if (ice(ic,mi)%snow%z .gt. 0) then
              do jj=ice(ic,mi)%snow%z,1,-1

                  ! find snowdepth for attenuation
                  th1 = ice(ic,mi)%snow%th(jj)*snow_fudge*lda_multiplier(kk)
                  ! find snow-ice equivalent depth for absorption
                  th2 = th1*ice(ic,mi)%snow%d(jj)/IceD

                  ! deal with absorbed and then full attenuation of PAR
                  do ii=1,wavl

                      nm_i = ii*10-9

					  ! snow spectral extinction coeffcient
					  if ((rs_switch .eq. 0 .and. airtemp1c .lt. 0.) .or. &
					  (rs_switch .eq. 1 .and. ice(ic,mi)%snow%t(jj) .le. -1.) .or. &
					  (rs_switch .eq. 2 .and. ice(ic,mi)%snow%melt(jj) .lt. 0.5)) then
					      Kds = Kds_dry(ii)
					  else
					      Kds = Kds_wet(ii)
					  endif

                      ! tmp3 = irradiance at lambda after absorption
                      tmp3 = Ed(1,ii)*exp(-1.*aice(nm_i)*th2) ! µEin/m^2/s
                      ! difference between original value and tmp3 = absorbed, convert to watts and add to Ed_W_Snow
			          Ed_W_snow(jj) = Ed_W_snow(jj) + (Ed(1,ii) - tmp3)*lda_f* &
			              quanta_to_watts(ii)*10.  ! Watts/m^2
                      ! full attenuation of PAR wavelength through snow layer
                      Ed(1,ii) = Ed(1,ii) * exp(-1.*Kds*th1)
                  enddo

			      ! absorb near ir
			      tmp3 = Ed_IR(1)*exp(-1.*a_ice_ir_cos*th2)  ! tmp3 = near-IR irradiance (Watts/m^2)
			      Ed_W_snow(jj) = Ed_W_snow(jj) + (Ed_IR(1) - tmp3)*lda_f  ! record absorbed near-IR in layer

			      ! full attenuation of near-IR for next timestep
			      ! substituting using <Kds from 700nm> for <near-IR backscattering> - probably results
			      ! in a bit of over-estimation of near-IR attenuation, under-estimation of
			      ! absorption
			      Ed_IR(1) = Ed_IR(1)*exp(-1.*(a_ice_ir_cos*th2 + Kds*th1))

              enddo ! end of snow layer loop

          !else
          ! if snow is below consideration, pass all PAR and near-IR to ice routines below
			  !if (h_snowpatch .eq. 0.) then
			  !    tmp2 = 1.
			  !else
                  !tmp2 = tmp1/(tmp1 + h_snowpatch)
			  !endif
    	      !Ed_IR(1) = Ed_IR(1)*(par_to_swd-1.)*(1.-tmp2)
          endif

          ! estimation particle absorption from chlorophyll and detrital concentrations in ice
          do jj=1,sk_z
              am_sum = 0.

              do ii=1,wavl
                  nm_i = ii*10-9
                  am(jj,ii) = aph(nm_i)*ice(ic,mi)%smalg(jj)/c_chl*(ice(ic,mi)%bv(jj)*c_001)  ! algae in bulk ice conc.
                  am_sum = am_sum+am(jj,ii)*10.   ! x10 b/c 10 nm wavelength bins
              enddo
              !ad(jj,1) = delta*am_sum/ad_denom
              ad(jj,1) = ice(ic,mi)%poc(jj)/ice(ic,mi)%smalg(jj)*am_sum/ad_denom ! using fractional amount of poc comprated to live algae
              do ii=2,wavl
                  ad(jj,ii) = ad(jj,1)*exp(-0.008*((dble(ii)-1.)*10.))
              enddo
          enddo

          ! perform 1-way attentuation with total attentuation k = kdice + kdp, and
          ! determine absorption component for heat flux
          PUR_dble(:,kk) = 0.
          do jj=1,sk_z
              do ii=1,wavl
				  ! calc ice spectral attenuation
				  nm_i = ii*10-9

				  ! Fresh Ice Kd ...
				  !tmp3 = 0.92 - 3.522e-3*lambda(ii) + 3.632e-6*lambda(ii)**2 ! pure ice - from arrigo
                  bv_mean = ice(ic,mi)%bv(jj)*c_001
                  tmp3 = (aice(nm_i)*(1. - bv_mean) + bv_mean* &
                  (awater(nm_i) - awater_sc(nm_i)*ice(ic,mi)%s(jj) - (0. - 22.)*awater_tc(nm_i)))/0.656

				  ! Salty Ice Kd ...
				  if (ice(ic,mi)%t(jj) .lt. -21.) then
					  ! Equation from Arrigo 1991 paper - in meters, not nm?
					  tmp2 = 65.114 - 0.23551*lambda(ii) + 0.00023078*lambda(ii)**2
					  tmp2 = tmp2*ice(ic,mi)%s(jj)/10.    ! scale linearly by bulk salinity, assuming Perovich 1979 used 10ppt ice for measurments here...
					  ! Equation from Arrigo BASIC code - seems to be in nm, not meters
					  !tmp2 = 1.2226-(0.005585*lambda)+(8.237E-06*lambda**2) &
					  !-(3.6E-09*lambda**3)
				  else
					  ex(jj,ii) = 4.31-0.021295*lambda(ii)+3.657e-5*lambda(ii)**2-1.9819e-8*lambda(ii)**3
!                      if (ice(ic,mi)%id(jj) .le. ice(ic,mi)%fbh) then
                          ! bubbly porous surface ice
!                          bv_mean_real = 80.
!					      tmp2 = dble(bv_mean_real**ex(jj,ii))
!                      else
                          ! internal ice
						  if (ice(ic,mi)%bv(jj) .gt. 80.) then
							  bv_mean_real = 80.
						  else
							  bv_mean_real = real(ice(ic,mi)%bv(jj))
						  endif
						  ! Equation from
						  ! arrigo BASIC code: Kdice = 0.027*brine_v^ex, as illustated in Fig. 1,
						  ! Arrigo 1991
						  tmp2 = dble(0.27*bv_mean_real**ex(jj,ii))
!					  endif
				  endif
				  ! use max Kd of fresh/salt ice equations
				  Kdice(jj,ii) = max(tmp3,tmp2);
				  ! calc particle spectral absorption
				  Kdp(jj,ii) = (am(jj,ii)+ad(jj,ii))/0.656    ! mu=0.656
				  ! attenuate downwards
				  Ed(jj+1,ii) = Ed(jj,ii)*exp(-1.*(Kdp(jj,ii)+Kdice(jj,ii)) &
				  *ice(ic,mi)%th(jj))
                  ! record bottom PAR
                  if (jj .eq. sk_z) then
                  	  ice(ic,mi)%PAR_bot = ice(ic,mi)%PAR_bot*10.*lda_f
                  endif
				  ! find Ed absorbed, convert to watts, and store in Ed_watts for layer
                  tmp1 = Ed(jj,ii)*exp(-1.*(Kdp(jj,ii) + tmp3)*ice(ic,mi)%th(jj)) ! absorption only - no scattering
   				  Ed_W_mean(jj) = Ed_W_mean(jj) + (Ed(jj,ii) - tmp1)* &
					quanta_to_watts(ii)*lda_f*10.  ! Watts/m^2
				  ! integrate mean Ed across spectrum to find PUR (add Ed*10 to total for layer)
				  PUR_dble(jj,kk) = PUR_dble(jj,kk) &
				      + 0.5*(Ed(jj+1,ii)+Ed(jj,ii))* &
					  aph(nm_i)/aph_max*10. ! times 10, b/c 10nm wavelength bins
               enddo ! end spectral do
               ! calculation Ed_W absorbed for heat flux later
               ! -------------------------------------------------------

			   ! find near IR absorption - consider scattering from 700nm
			   kdice_ir = a_ice_ir_cos + Kdice(jj,wavl) ! total near-IR ice attentuation = absorption + scattering at 700nm

			   ! absorb near ir
			   tmp2 = Ed_IR(jj)*exp(-1.*a_ice_ir_cos*ice(ic,mi)%th(jj))
			   Ed_W_IR = (Ed_IR(jj) - tmp2)

			   ! total attenuation of near-IR for next timestep
			   Ed_IR(jj+1) = Ed_IR(jj)*exp(-1.*kdice_ir*ice(ic,mi)%th(jj))

			   ! update
			   Ed_W_mean(jj) = Ed_W_mean(jj) +  Ed_W_IR*lda_f

          enddo ! end ice layer do

! ==============================================================================
! END CALC LIGHT
! ==============================================================================
