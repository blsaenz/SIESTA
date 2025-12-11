
!      if (ice(ic,mi)%id(int_z) .ge. 0.2) then

          bv_total = 0.
          bv_flux = 0.
          no3_total = 0.
          nh4_total = 0.
          po4_total = 0.
          sioh4_total = 0.
          poc_total = 0.
          heat_total = 0.
          b_flux_max = 0.
          jjj = 0
          heat_stop = int_z + 1

          debug_z = ice(ic,mi)%no3(1:z_max)

          if (dt_step .eq. 3 .and. iii .eq. (dtt - 1)) then
              testvar = -1
          endif


          do ii=1,int_z

              dsdt_cw(ii) = 0.
              dsdt1 = 0.
              dsdt2 = 0.
              dsdt3 = 0.
              dcheat(ii) = 0.
              b_flux(ii) = 0.
              h_flux(ii) = 0.
		      dbvdt = 0.
		      F0 = 0.
              bv_mean = ice(ic,mi)%bv(ii)  ! ppt
              dz = ice(ic,mi)%bv(ii)*c_001*ice(ic,mi)%th(ii) ! layer brine volume (m^3)

              ! gravity drainage: only if above critical brine volume,
              ! and density (salinity) of brine is higher than of surrounding seawater
              if ((bv_mean .ge. vb_crit) .and. ((ii .gt. vb_open) &
!                .and. (ice(ic,mi)%bs(ii) .gt. f(mi)%s)) &
!                .or. ((ice(ic,mi)%bced(ii) .eq. 1) .and. (bv_mean .ge. vb_crit)) &
                )) then

                  if (ice(ic,mi)%bs(ii) .gt. f(mi)%s) then

                       ! temp gradient used in gravity drainage eqn. is calculated
                      ! as the average between gradients above and below, then
                      ! averaged across this and last timestep

                      ! Cox-Weeks gravity drainage
                      ! -----------------------------------------------
                      if (ii .eq. 1) then
!                          if (z_snow .gt. 0) then
!                              T_grad = (ice(ic,mi)%snow%t(z_snow) - ice(ic,mi)%t(ii+1)) / &
!                                 (c_5*(ice(ic,mi)%th(ii+1) + ice(ic,mi)%th(ii) + ice(ic,mi)%snow%t(z_snow)))
!                          else
                              T_grad = (ice(ic,mi)%t(ii) - ice(ic,mi)%t(ii+1)) / &
                                  (c_5*(ice(ic,mi)%th(ii+1) + ice(ic,mi)%th(ii)))
!                          endif

                      elseif (ii .eq. int_z) then
                           T_grad = (ice(ic,mi)%t(ii-1) - f(mi)%t) / &
                               (c_5*ice(ic,mi)%th(ii-1) + ice(ic,mi)%th(ii))
                      else
                           T_grad = (ice(ic,mi)%t(ii-1) - ice(ic,mi)%t(ii+1)) / &
                               (c_5*(ice(ic,mi)%th(ii-1) + ice(ic,mi)%th(ii+1)) + ice(ic,mi)%th(ii))
                      endif

                      ice(ic,mi)%tgrad(ii) = T_grad

                      !dsdt1 = -1.*T_grad*c_01* & ! 0.01 factor is to convert from 1/m to 1/cm
                      !(1.68E-5 - 3.37E-7*bv_mean)*dtt_s ! dsdt in ppt/s, so multiplying by sec., the minus sign gives a negative dsdt

                      dsdt1 = 4.2e-6*dtt_s*T_grad*(max(bv_mean*c_001 - 0.05,1.d-5))**(1.2) ! Petrich et al. desal rate

                      ! record instantaneous dsdt rate (psu/s)
                      ice(ic,mi)%dsdt(ii) = -1.*T_grad*c_01* & ! 0.01 factor is to convert from 1/m to 1/cm
                         (1.68E-5 - 3.37E-7*bv_mean)


                      ! estimate conductive heat flux into layer
                      if (ii .eq. 1) then
                          if (z_snow .gt. 1) then
                             tmp1 = (ice(ic,mi)%t(1) - ice(ic,mi)%snow%t(1))/dth(ice_1)
                          else
                             tmp1 = (ice(ic,mi)%t(1) - ice(ic,mi)%snow%ts)/dth(ice_1)
                          endif
                          tmp2 = (ice(ic,mi)%t(1) - ice(ic,mi)%t(2))/dth(ice_1+1)
                      else
                          tmp1 = (ice(ic,mi)%t(ii) - ice(ic,mi)%t(ii-1))/dth(z_snow+ii)
                          tmp2 = (ice(ic,mi)%t(ii) - ice(ic,mi)%t(ii+1))/dth(z_snow+ii+1)
                      endif

                      F0 = tmp1*ki(z_snow+ii) + tmp2*ki(z_snow+ii+1) - ed_w_mean(ii)


                      ! Flux-based desal
                      ! -----------------------------------------------
!                      if (iii .ge. 1) then

                          ! F0 calc from temp gradient - moved outside loop for testing...
                          ! --------------------------------------------
                          !ki_down = (d_new(ii)/IceD)*(c0-0.011*ice(ic,mi)%t(ii) + &
                          !  0.09*s_new(ii)/ice(ic,mi)%t(ii) - (d_new(ii)-IceD)/1.e6)
                          !F0 = -1.*T_grad*ki_down

                          ! F0 = heat_last(ii)*ice(ic,mi)%th(ii)/dtt_s

                          !dhdt = 5.356e-12*F0**3 - 6.069e-10*F0**2 + 4.815e-7*F0  ! (cm/s) dhdt vs. F0 regression at T=-2 freezing temp, cox weeks keff

                          !dhdt = 9.759e-12*F0**3 - 6.529e-10*F0**2 + 4.463e-7*F0  ! (cm/s) dhdt vs. F0 regression at T=-2 freezing temp, petrich keff
                          dhdt = 9.873e-12*F0**3 - 6.621e-10*F0**2 + 4.474e-7*F0  ! (cm/s) dhdt vs. F0 regression at T=-1.8 freezing temp, petrich keff
                          !dhdt = 9.988e-12*F0**3 - 6.714e-10*F0**2 + 4.484e-7*F0  ! (cm/s) dhdt vs. F0 regression at T=-1.6 freezing temp, petrich keff


!                      else
                          ! use dhdt from last timestep if this is 1st sub-dy timestep and there is no heat_diff/t-diff yet
!                          dhdt = ice(ic,mi)%dhdt_conv(ii)
!                      endif

                      ! record instaneous flux-based desal vars
                      ice(ic,mi)%dhdt_conv(ii) = dhdt
                      ice(ic,mi)%f0(ii) = F0

                      ! assume F0 is used to freeze brine, with stable salinity type desalintion
                      ! calculated using keff

                      if (dhdt .gt. 0.) then
                          ! keff from cox/weeks
                          !if (dhdt .gt. 3.6e-5) then
                          !    keff=0.26/(0.26+0.74*exp(-7243.0*dhdt))
                          !else
                          !    if (dhdt .lt. 2.0e-6) then
                          !       keff=0.12
                          !    else
                          !       keff=0.8925+0.0568*log(dhdt)
                          !    endif
                          !endif

                          ! keff from petrich
                          keff = 0.19d0*(dhdt*7.4074d4)**(0.46)
                          !keff=max(keff,0.12)
                          !print *,'keffs: ',keff,tmp4


                          ! record convective desal rate
                          ice(ic,mi)%dsdt3(ii) = -1.*dhdt*c_01*ice(ic,mi)%bs(ii)*(1.-keff)/ice(ic,mi)%th(ii)

                      endif

                      ! decide how to operate - if salinity is not stable,
                      ! hold temperature constant and do and convective desalination
                      ! otherwise use cox-weeks not-too-convective desal
                      if (bv_mean .ge. bv_conv .and. dhdt .gt. 0.) then

                          if (s_new(ii) .gt. keff*ice(ic,mi)%bs(ii)) then

                              ! record start if convective brine movement, or add to it
                              jjj = jjj+1

                              ! add heat back to layer for next time step to
                              ! keep temp relativey constant
                              !dcheat(ii) = heat_last(ii)*ice(ic,mi)%th(ii)/dtt_s
                              dcheat(ii) = F0
                              ! m/cm * cm/s * s * m/m * g/m^3 * J/g / m = J

                              !dcheat(ii) = c_01*dhdt*keff*IceD/ice(ic,mi)%th(ii)
                              !print *,ii,'F0: ',F0,'dcheat: ',dcheat(ii)

                              if (heat_stop .gt. int_z) then
                                  heat_stop = ii
                              endif

                              ! do agressive desal
                              dsdt3 = -1.*dtt_s*dhdt*c_01*ice(ic,mi)%bs(ii)*(1.-keff)/ice(ic,mi)%th(ii)

                              ! scale dsdt according to bv^3
                              tmp1 = (ice(ic,mi)%bv(ii)-bv_conv)   ! difference
                              tmp2 = (400.-bv_conv)       ! scale
                              dsdt3 = min(dsdt1,dsdt3)   ! find bigger one
                                                         ! assuming bv = 0.5 is wide open
!                              if (tmp1 .lt. tmp2) then
!                                  dsdt3 = dsdt1 + (dsdt3-dsdt1)*(tmp1/tmp2)**3
!                              endif

                              ! take larger amount of two desal values (more negative)
                              dsdt1 = min(dsdt1,dsdt3)

                          else
                              ! if convection has started above somewhere already, mix with below layers
                              ! that are above bv_conv below
                              if (jjj .gt. 0 .and. bv_mean .ge. bv_conv) then
                                  jjj = jjj+1
                              endif
                          endif

                      endif

                      ! prevent increasing salinity
                      dsdt1 = min(0.,dsdt1)

                      ! bound salinity at sal_min, don't gravity drain upward...
                      if ((s_new(ii) + dsdt1) .lt. min_sal) then
                          dsdt1 = min_sal - s_new(ii)
                      endif

                      ! record desal for use later during environment re-calculation
                      dsdt_cw(ii) = dsdt_cw(ii) + dsdt1

                  else


                      ! no gravity desal
                      ice(ic,mi)%dsdt(ii) = 0.
                      ice(ic,mi)%dhdt_conv(ii) = 0.
                      ice(ic,mi)%dsdt3(ii) = 0.
                      ice(ic,mi)%f0(ii) = 0.
                      ice(ic,mi)%tgrad(ii) = 0.

                      ! if convection has started above somewhere already, mix with below layers
                      ! that are above bv_conv below
                      if (jjj .gt. 0 .and. bv_mean .ge. bv_conv) then
                          jjj = jjj+1
                      endif

                  endif  ! end of bs gradient check for gravity drainage

				  ! brine expulsion - only do ice volume adjustment if ice is warming!
				  ! only do brine explusion if gravity drainage is not operating?
!	               if (t_change(ii) .gt. 0.) then
!                  if (abs(bs_diff(ii)) .ge. 1.d-10) then

                     ! simplified eqn below
!	                  dsdt2 = s_new(ii)*(  &
!	                     (bd_last(ii))* &
!	                     (bs_last(ii))**(1.-1.e6/IceD)* &
!	                     exp(800./IceD*(bs_diff(ii))) - 1.)

!	                  dsdt2 = s_new(ii)*(  &
!	                     (bd_last(ii))* &
!	                     (bs_last(ii))**(0.8911)* &
!	                     exp((8.7061e-04)*(bs_diff(ii))) - 1.d0)

                      ! bound salinity at sal_min
!                      if ((s_new(ii) + dsdt1 + dsdt2) .lt. min_sal) then
!                        dsdt2 = min_sal - s_new(ii) - dsdt1
!                      endif
!	                  dsdt_cw(ii) = dsdt_cw(ii) + dsdt2

!	              endif
!	              endif

              endif ! end of check for open brine channels

              ! determine brine/seawater fluxes from dsdt
			  if (dsdt1 .lt. 0.) then
                  bv_mean = ice(ic,mi)%bv(ii)*c_001

                  ! assuming that this desal happens at a constant temperature,
                  ! we ignore the volume change in ice, if any - since we don't know
                  ! the convective heat flux yet, the amount of ice grown
                  ! is independent of the desal process

                  ! limit dilution to a minimum 2 psu different - lower can cause
                  ! excessive fluid movement
                  tmp1 = f(mi)%s - ice(ic,mi)%bs(ii)
!                  tmp2 = ice(ic,mi)%bs(ii+1) - ice(ic,mi)%bs(ii)
                  tmp2 = tmp1
                  tmp1 = min(-2.,tmp1)
                  tmp2 = min(-2.,tmp2)

                  ! calculate the percent of brine volume fluxed assuming dilution
                  ! at constant temperature
!                  if(ice(ic,mi)%bced(ii) .eq. 1) then
                  if (ice(ic,mi)%bv(ii) .ge. bv_conv) then
                      ! dilute the original bs with seawater - assume direct connection
                      ! with seawater
                      dbvdt = ((ice(ic,mi)%bs(ii)*bv_mean + dsdt1)/bv_mean - &
                          ice(ic,mi)%bs(ii))/(tmp1)
                  else
                      ! dilute the original bs to bs_dilute with brine below
                      dbvdt = ((ice(ic,mi)%bs(ii)*bv_mean + dsdt1)/bv_mean - &
                          ice(ic,mi)%bs(ii))/(tmp2)
                  endif

                  ! convert to m^3/timestep
                  dbvdt = dbvdt*ice(ic,mi)%th(ii)*bv_mean

				  ! old dbvdt calc
!				  S_mean = s_new(ii) + dsdt1
!			      bs_dilute = 1000.*ice(ic,mi)%d(ii)*S_mean/ &
!					  (ice(ic,mi)%bd(ii)*ice(ic,mi)%bv(ii))
!				  bs_dilute = ice(ic,mi)%bs(ii)*(s_mean/s_new(ii))   ! effectively the same thing as above
!                 if (ice(ic,mi)%bv(ii) .ge. bv_conv) then
                      ! dilute the original bs to bs_dilute with seawater - assume direct connection
                      ! with seawater
!                      dbvdt = 0.918*ice(ic,mi)%bv(ii)*ice(ic,mi)%th(ii)*c_001* &
!					      (bs_dilute - ice(ic,mi)%bs(ii))/(f(mi)%s - bs_dilute)
                      ! replacement dilution
!                   else
                      ! dilute the original bs to bs_dilute with brine below
!                      dbvdt = 0.918*ice(ic,mi)%bv(ii)*ice(ic,mi)%th(ii)*c_001* &
!                        (bs_dilute - ice(ic,mi)%bs(ii))/(ice(ic,mi)%bs(ii+1) - bs_dilute) ! m^3/m^2/s
!                      dbvdt = 0.918*ice(ic,mi)%bv(ii)*ice(ic,mi)%th(ii)*c_001* &
!                        (bs_dilute - ice(ic,mi)%bs(ii))/(f(mi)%s - bs_dilute) ! m^3/m^2/s
!                  endif



		           !fvs = -9.652e-9 + 0.387*dhdt - 1029.1*dhdt**2  ! cm^3/cm^2/s - from arrigo basic bodel
		           !dbvdt = c_01*(ice(ic,mi)%bv(ii)*c_001)*fvs                       ! m^3/m^2/s







                  ! scale volume flux
                  dbvdt = dbvdt*dbvdt_scale

                  if (ice(ic,mi)%bv(ii) .ge. bv_conv .and. dhdt .gt. 0.) then

                      ! avaiable convective heat
                      tmp1 = dbvdt*cw*(f(mi)%t-ice(ic,mi)%t(ii))*f(mi)%d

                      ! heat capacity of layer --> seawater temp
                      T_melt = mu*s_new(ii)
			          tmp2 = -1.*d_new(ii)*(c0*(f(mi)%t - ice(ic,mi)%t(ii)) &
				         - Lf*mu*s_new(ii)*(1./f(mi)%t - 1./ice(ic,mi)%t(ii))) &
				         * ice(ic,mi)%th(ii) ! j/m^3 - not j/m^2

				      !print *,'F0: ',F0,'heat brine: ',tmp1, 'layer_heat: ',tmp2

!				      dcheat(ii) = min(tmp1,tmp2)
!				      heat_total = heat_total + tmp1

				  else

                      ! avaiable convective heat
                      tmp1 = dbvdt*cw*(ice(ic,mi)%t(ii+1)-ice(ic,mi)%t(ii))*f(mi)%d

                      ! heat capacity of layer --> seawater temp
                      T_melt = mu*s_new(ii)
			          tmp2 = -1.*d_new(ii)*(c0*(ice(ic,mi)%t(ii+1) - ice(ic,mi)%t(ii)) &
				         - Lf*mu*s_new(ii)*(1./ice(ic,mi)%t(ii+1) - 1./ice(ic,mi)%t(ii))) &
				         * ice(ic,mi)%th(ii) ! j/m^3 - not j/m^2

				      !print *,'F0, slow: ',F0,'heat brine: ',tmp1, 'layer_heat: ',tmp2

!				      dcheat(ii) = min(tmp1,tmp2)
!				      heat_total = heat_total + tmp1


				  endif



			  elseif (dsdt2 .gt. 0.) then

!                  if(ice(ic,mi)%bced(ii) .eq. 1) then
                  if (ice(ic,mi)%bv(ii) .ge. bv_conv) then
                      ! add volume of seawater needed to reach new salinity
                      dbvdt = dsdt2*ice(ic,mi)%th(ii)/f(mi)%s
                  else
                      ! add volume of brine needed to reach new salinity
                      dbvdt = dsdt2*ice(ic,mi)%th(ii)/ice(ic,mi)%bs(ii+1)
                  endif

			  endif

			  ! protect against negative b_flux - happens when bs_dilute
			  ! is below the seawater salinity
			  dbvdt = max(0.,dbvdt)

              ! test to see if desal is happening too fast to account for heat/nutrient transfer

              if (ice(ic,mi)%bv(ii) .ge. bv_conv) then
                  ! stop nutrient flux from above, if there is a layer above
                  b_flux(ii) = 0.
                  h_flux(ii) = dbvdt ! record assumed horizontal/direct ocean flux for ki adjustment

                  ! test to see if desal is happening too fast to account for full nutrient transfer
                  ! if (dbvdt .gt. dz) then
                  !    print *,'Barf! desal convection too fast: ',ii,dbvdt,dz
                  !    if (dtt_s .eq. dtt_s_3) then
                  !        integr_valid = .false.
                  !    endif
                  !endif

                  ! add to volume fluxed in convective layer
                  bv_flux = bv_flux + dbvdt

                  ! add to convective layer brine volume & tracer totals
                  bv_total = bv_total + dz
                  no3_total = no3_total + dz*no3_new(ii)
                  nh4_total = nh4_total + dz*nh4_new(ii)
                  po4_total = po4_total + dz*po4_new(ii)
                  sioh4_total = sioh4_total + dz*sioh4_new(ii)
                  poc_total = poc_total + dz*poc_new(ii)

              endif

              if ((ice(ic,mi)%bv(ii) .lt. bv_conv) .or. (ii .eq. int_z)) then

                  ! test to see if final ice layer is also 'open' convecting
                  ! with the ocean.  if so, we need to include it in the
                  ! covection routine below by modifying mo
                  if ((ii .eq. int_z) .and. ice(ic,mi)%bv(ii) .ge. bv_conv) then
                      mo = 1
                  else
                      mo = 0

                      ! take care of current layer, which uses the standard brine volume
                      ! flux rountine to mix nutrients - heat is supposedly taken
                      ! care of in this case already by the ki of sea ice eqn.
                      if (ii .eq. 1) then
                          b_flux(ii) = dbvdt
                      else
                          b_flux(ii) = b_flux(ii-1) + dbvdt
                      endif

                      ! ignoring convective heat through brine flux - is very small

                      ! don't modify
                      h_flux(ii) = 0.

                      tmp1 = dz*0.3
                      do while (b_flux(ii)/dble(bfrf) .ge. tmp1)
                          bfrf = bfrf*2
                          print *,ii,' - Internal Volume fluxed (',b_flux(ii)/dble(bfrf),') greater than 1/3 BV (',tmp1,').'
                          !print *,ii,' - increasing flux resolution: ',bfrf
                      enddo

                  endif

                  ! solve for new nutrients and heat flux for above 'open'
                  ! convecting layers, now that we are not convecting anymore
                  ! (or we have reached the last ice layer above water)
                  if (jjj .gt. 0) then

                      bv_flux = min(bv_total,bv_flux)
                      tmp1 = 1./(bv_total)
                      tmp2 = (bv_total-bv_flux)*tmp1
                      no3_total = (no3_total*tmp2 + bv_flux*f(mi)%no3)*tmp1
                      nh4_total = (nh4_total*tmp2 + bv_flux*f(mi)%nh4)*tmp1
                      po4_total = (po4_total*tmp2 + bv_flux*f(mi)%po4)*tmp1
                      sioh4_total = (sioh4_total*tmp2 + bv_flux*f(mi)%sioh4)*tmp1
                      poc_total = (poc_total*tmp2 + bv_flux*f(mi)%poc)*tmp1

                      do jj=ii-1+mo,ii-jjj+mo,-1

                          ! update salinity and density (common to all layers)

                          if (jj .eq. 0) then
                              testvar = -1.
                          endif

                          ! change convective nutrients, using bv_mean to account for
                          ! brine volume change that happens later
!                          if (dsdt_cw(jj) .ne. 0.) then
!                              s_mean = s_new(jj) + dsdt_cw(jj)
!                              bs_mean = bs_new(jj)
                              ! calculate new brine density (c=800 g m-3 ppt-1)
!                              bd_mean = 1000000. + bs_mean*800.       ! g/m^3
                              ! new ice density
!                              d_mean = (1.-bb_f)*IceD*bd_mean*bs_mean/(bd_mean*bs_mean - s_mean* &
!                                  (bd_mean - IceD))
                              ! new brine volume again after desalination...
!                              bv_mean = 1000.*(d_mean*s_mean)/(bd_mean*bs_mean)  ! ppt

!                              tmp2 = bv_mean/bv_new(jj)
!                          else
                              tmp2 = 1.
!                          endif
                          dno3(jj) = dno3(jj) + (no3_total*tmp2 - no3_new(jj))
                          dnh4(jj) = dnh4(jj) + (nh4_total*tmp2 - nh4_new(jj))
                          dpo4(jj) = dpo4(jj) + (po4_total*tmp2 - po4_new(jj))
                          dsioh4(jj) = dsioh4(jj) + (sioh4_total*tmp2 - sioh4_new(jj))
                          dpoc(jj) = dpoc(jj) + (poc_total*tmp2 - poc_new(jj))

                          ! heat capacity of layer --> seawater temp
!                          T_melt = mu*s_new(ii)
!                          tmp2 = -1.*d_new(ii)*(c0*(f(mi)%t - ice(ic,mi)%t(ii)) &
!                             - Lf*mu*s_new(ii)*(1./f(mi)%t - 1./ice(ic,mi)%t(ii))) &
!                             * ice(ic,mi)%th(ii) ! j/m^3 - not j/m^2

                          ! available heat into layer
!                          tmp3 = heat_total*ice(ic,mi)%bv(ii)*c_001* &
!                              ice(ic,mi)%th(ii)*tmp1

!                          dcheat(jj) = min(tmp2,tmp3)*ice(ic,mi)%bv(ii)*c_001* &
!                              ice(ic,mi)%th(ii)*tmp1

                      enddo

                      ! reset convecting brine volumes & tracers to 0
                      bv_flux = 0.
                      bv_total = 0.
                      no3_total = 0.
                      nh4_total = 0.
                      po4_total = 0.
                      sioh4_total = 0.
                      poc_total = 0.
                      heat_total = 0.
                      jjj = 0

                  endif

              endif  ! end of convection/bottom of ice check

          enddo    ! end of desal+flux volume calculation

          if (integr_valid) then

              do ii=1,int_z

                  ! record vars
                  ice(ic,mi)%dsdt(ii) = ice(ic,mi)%dsdt(ii) + dsdt_cw(ii)
!                  ice(ic,mi)%fbv(ii) = ice(ic,mi)%fbv(ii) + h_flux(ii)

              enddo

              ! pass on brine convection through skeletal layers
              do ii=sk_1,sk_z
                  b_flux(ii) = b_flux(int_z)
                  tmp1 = ice(ic,mi)%bv(ii)*c_001*ice(ic,mi)%th(ii)*0.3
                  do while (b_flux(ii)/dble(bfrf) .ge. tmp1)
                      bfrf = bfrf*2
                      print *,ii,' - Skeletal Volume fluxed (',b_flux(ii)/dble(bfrf),') greater than 1/3 BV (',tmp1,').'
                      !print *,ii,' - increasing flux resolution: ',bfrf
                  enddo
              enddo

          endif ! end of valid integration check

!      endif  ! end of minimum thickness for desal check

!          print *,'bced: ',ice(ic,mi)%bced
!          print *,'dsdt: ',dsdt
!          print *,'dheat: ',dheat
!          print *,'no3: ',ice(ic,mi)%no3
!          print *,'dcheat: ',dcheat(1:15)
!          print *,'light heat: ',Ed_W_mean
!          print *,'t: ',ice(ic,mi)%t
!          print *,'===================================================='

