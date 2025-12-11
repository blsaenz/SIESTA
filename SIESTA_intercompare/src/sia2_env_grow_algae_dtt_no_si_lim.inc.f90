              if (alg_dtt .eq. 1) then
                  do ii=1,sk_z
                      ! return to brine conc.
                      bv_mean = ice(sc,ic,mi)%bv(ii)*c_001
                      ice(sc,ic,mi)%no3(ii) = ice(sc,ic,mi)%no3(ii)/bv_mean
                      ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii)/bv_mean
                      ice(sc,ic,mi)%po4(ii) = ice(sc,ic,mi)%po4(ii)/bv_mean
                      ice(sc,ic,mi)%sioh4(ii) = ice(sc,ic,mi)%sioh4(ii)/bv_mean
                      ice(sc,ic,mi)%smalg(ii) = ice(sc,ic,mi)%smalg(ii)/bv_mean
                      ice(sc,ic,mi)%poc(ii) = ice(sc,ic,mi)%poc(ii)/bv_mean
                  enddo


! ==============================================================================
! START GROW ALGAE
! ==============================================================================

      ! test to see if there is any light - if not, don't calc algal growth, just die
      ! ----------------------------------------------------------------
      if (ice(sc,ic,mi)%Ed0 .gt. 0.) then

          ! init some vars for keeping track of things though lda looping
          do ii=1,sk_z
              smalg_new(ii) = 0.
              poc_new(ii) = 0.
              dN_all(ii) = 0.
              !NO3_new(ii) = 0.
              dNH4(ii) = 0.
              dPO4(ii) = 0.
              dSiOH4(ii) = 0.
          enddo

          bot_id = ice(sc,ic,mi)%id(ice(sc,ic,mi)%z) - bot_th

          ! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
          ! ------------------------------------------------------------
          ! START LIGHT ADJUSTMENT (LDA) LOOP
          ! ------------------------------------------------------------
          ! ------------------------------------------------------------
          do kk=1,lda_n
          ! ------------------------------------------------------------
          ! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
              ! what areal fraction does this lda category represent?
              lda_f = lda_d

              ! loop thorugh ice layers, calc limitations, ik', growth in each
              do ii = 1,sk_z

                  ! Calc Ikprime
                  !--------------------------------------------------------------
!                  PUR_mean24 = 0.
!                  do jj=1,dt_per_day
!                      PUR_mean24 = PUR_mean24 + PUR(ii,jj,kk,sc,ic,mi)
!                  enddo
!                  PUR_mean24 = PUR_mean24/dble(dt_per_day)
!
!                  ice(sc,ic,mi)%Ik1(ii) = Ek_max/(1. + A*exp(-1.*B*PUR_mean24))
!                  !tmp1 = explt(-1.*B*PUR_mean24,use_expt,exp_bins,exp_mul,expt)
!                  !ice(sc,ic,mi)%Ik1(ii) = Ek_max/(1. + A*tmp1)

                  ! Calculate nutrient, light, and salinity limitations
                  ! -------------------------------------------------------------
                  NO3_lim = ice(sc,ic,mi)%no3(ii) / (Ks_NO3 + ice(sc,ic,mi)%no3(ii))
                  NH4_lim = ice(sc,ic,mi)%nh4(ii) / (Ks_NH4 + ice(sc,ic,mi)%nh4(ii))
                  ! find total N limitation term, as the last limiting form of N
                  NO3_lim = max(NO3_lim,NH4_lim)
                  PO4_lim = ice(sc,ic,mi)%po4(ii) / (Ks_PO4 + ice(sc,ic,mi)%po4(ii))
                  SiO4_lim = ice(sc,ic,mi)%sioh4(ii) / (Ks_SiO4 + ice(sc,ic,mi)%sioh4(ii))

                  SiO4_lim = 1.d0

!                  if (ice(sc,ic,mi)%Ik1(ii) .le. 0.) then
!                      light_lim = 0.
!                  else
!                      light_lim = 1. - exp(-1.*PUR(ii,pur_clock,kk,sc,ic,mi)/ice(sc,ic,mi)%Ik1(ii))
!                      !light_lim = 1. - explt(-1.*PUR(ii,pur_clock,kk,sc,ic,mi)/ice(sc,ic,mi)%Ik1(ii), &
!                      !    use_expt,exp_bins,exp_mul,expt)
!
!                  endif
                  !tmp1 = 0.
                  !do jj=1,dt_per_day
                  !    tmp1 = tmp1 + 1. - exp(-1.*PUR(ii,jj,kk,sc,ic,mi)/ice(sc,ic,mi)%Ik1(ii))
                  !enddo
                  !tmp1 = tmp1/dble(dt_per_day)
                  !print *,ii,': light_lim=',light_lim,' tmp1=',tmp1

                  ! find salinity limitation: Arrigo (199?)
                  if (ice(sc,ic,mi)%bs(ii) .gt. 103) then
                      sal_lim = 0.01
                  else
                      sal_lim = a0 + a1*ice(sc,ic,mi)%bs(ii)  &
                          + a2*ice(sc,ic,mi)%bs(ii)**2 &
                          + a3*ice(sc,ic,mi)%bs(ii)**3
                  endif

                  ! figure out if ice resolution/layers has grown, and if so invent an ll
                  if (ll(ii,kk) .eq. -1.) then
                      ! assume the difference between the last two ll's is the same difference
                      ! between new last two ll's - an underestamation most likely.
                      ll(ii,kk) = ll(ii-1,kk) - abs(ll(ii-2,kk) - ll(ii-1,kk))
                  endif

                  ! record growth limitations
                  ice(sc,ic,mi)%llim(ii) = ll(ii,kk)
                  ice(sc,ic,mi)%nlim(ii) = NO3_lim
                  ice(sc,ic,mi)%plim(ii) = PO4_lim
                  ice(sc,ic,mi)%silim(ii) = SiO4_lim
                  ice(sc,ic,mi)%slim(ii) = sal_lim

                  ! Determine which factor (nutrient or light) is limiting
                  ! -------------------------------------------------------------
                  min_lim = min(NO3_lim, PO4_lim, SiO4_lim, ll(ii,kk))
!                  if (ii.eq. 41) then
!                      print *,'Light Lim: ',light_lim
!                      print *,'Min Lim: ',min_lim
!                  endif

                  ! Calculate growth rate for this layer = Gmax*limitation - Respiration
                  ! -------------------------------------------------------------
                  Gmax = G0*exp(rg*ice(sc,ic,mi)%t(ii))*min_lim*sal_lim !&
                  !tmp1 = explt(rg*ice(sc,ic,mi)%t(ii),use_expt,exp_bins,exp_mul,expt)
                  !Gmax = G0*tmp1*min_lim*sal_lim !&

                  ! Grow
                  ! -------------------------------------------------------------
                  newp = ice(sc,ic,mi)%smalg(ii)*exp(Gmax*dtt_days) - ice(sc,ic,mi)%smalg(ii)

                  if (newp .gt. 0.) then
                      ! Calc new growth, production, nut. requirements, use up nutrients,
                      ! don't go negative
                      ! -------------------------------------------------------------
                      newp = newp*mC_gC  ! mMolC/m^3

                      ! find potential production from available nutrients
                      ! -------------------------------------------------------------
                      newpN = (ice(sc,ic,mi)%no3(ii)+ice(sc,ic,mi)%nh4(ii))*c_n  ! mMolN/m^3 * mMolC/mMolN = mMolC/m^3
                      newpP = ice(sc,ic,mi)%po4(ii)*c_p
                      newpSi = ice(sc,ic,mi)%sioh4(ii)*c_si
                      newp = min(newp,newpN,newpP,newpSi)

                      ! use up N
                      dN_all_temp = -1.*newp*n_c
                      ! use up DIP, Silica
                      dSiO4_temp = -1.*newp*si_c ! mMolSi/m^3
                      dPO4_temp = -1.*newp*p_c ! mMolP/m^3

                      newp = newp*gC_mC ! mg/m^3 [within brine]
                  else
                      dN_all_temp = 0.
                      dSiO4_temp = 0.
                      dPO4_temp = 0.
                  endif
                  dNH4_temp = 0.  ! default is zero - no remin

                  ! add to standing biomass
                  ! -------------------------------------------------------------
                  smalg_temp = ice(sc,ic,mi)%smalg(ii) + newp ! new standing biomass before grazing/death
                  death = smalg_temp - smalg_temp*exp(-1.*xp*dtt_days) ! death/grazing calucalted using new biomass ?
                  !tmp1 = explt(-1.*xp*dtt_days,use_expt,exp_bins,exp_mul,expt)
                  !death = smalg_temp - smalg_temp*tmp1 ! death/grazing calucalted using new biomass ?
                  smalg_temp = smalg_temp - death ! die
                  if (newp .lt. 0.) then
                     death = death - newp ! if respiration is causing net death, include in death term for remineralization
                  endif

                  !if (ii .eq. 41) then
                  !print *,'smalg_z initial: ',bio(sc,ic,mi)%smalg(ii)
                  !print *,'gmax: ',gmax*dt_days
                  !print *,'resp: ',R0*exp(rr*ice(sc,ic,mi)%t(ii))*dt_days
                  !print *,'growth: ',growth
                  !print *,'new_p: ',newp*gC_mC
                  !print *,'smalg new: ',smalg_temp
                  !endif

                  ! Remineralize
                  ! -------------------------------------------------------------
                  poc_temp = ice(sc,ic,mi)%poc(ii)
                  tmp1 = poc_temp*remin*dtt_days
                  dNH4_temp = tmp1*mC_gC*n_c*remin_f  ! all happens in same volume, so don't need to take into account thickness or brine volume...
                  dPO4_temp = dPO4_temp + tmp1*mC_gC*p_c*remin_f  ! all happens in same volume, so don't need to take into account thickness or brine volume...
                  poc_temp = poc_temp - tmp1

                  ! make sure algal concentration does not fall below minimum
                  min_alg_brine = min_alg/(ice(sc,ic,mi)%bv(ii)*c_001)
                  if (smalg_temp .lt. min_alg_brine) then
                      smalg_temp = min_alg_brine  !algae concentrated into brine volume
                  else
                      ! can add to poc, if not at algal minimum. Have to check because if we added
                      ! from death at algal minimum, poc would be created out of nothing
                      poc_temp = poc_temp + death
                  endif

                  ! record production
                  ! -------------------------------------------------------------
                  totalp = newp
                  if (totalp .gt. 0.) then
                      ! converting to g/pixel - (brine_v/1000)*(cell_area*1.e6)*1g/1000mg = brine_v*cell_area
                      production_temp = totalp*ice(sc,ic,mi)%bv(ii)*ice(sc,ic,mi)%th(ii)* &
                          ice(sc,ic,mi)%af*cell_area
                      ! dividing by number of snow depths to arrive at average at end of lda loop
                      ice(sc,ic,mi)%prod(ii) = ice(sc,ic,mi)%prod(ii) + production_temp*lda_f

                      ! record bottom or internal production
                      if (ice(sc,ic,mi)%id(ii) .gt. bot_id) then
                          if ((ice(sc,ic,mi)%id(ii)-ice(sc,ic,mi)%th(ii)) .ge. bot_id) then
                              m(mi)%prod_sum_bot = m(mi)%prod_sum_bot + production_temp
                              p_tot_bot = p_tot_bot + totalp*ice(sc,ic,mi)%bv(ii)*ice(sc,ic,mi)%th(ii)
                          else
                              bot_f = (ice(sc,ic,mi)%id(ii)-bot_id)/ice(sc,ic,mi)%th(ii) ! bottom ice fraction of layer
                              m(mi)%prod_sum_bot = m(mi)%prod_sum_bot + production_temp*bot_f
                              m(mi)%prod_sum_int = m(mi)%prod_sum_int + production_temp*(1.d0-bot_f)
                              p_tot_bot = p_tot_bot + totalp*ice(sc,ic,mi)%bv(ii)*ice(sc,ic,mi)%th(ii)*bot_f
                              p_tot_int = p_tot_int + totalp*ice(sc,ic,mi)%bv(ii)*ice(sc,ic,mi)%th(ii)*(1.d0-bot_f)
                          endif
                      else
                          m(mi)%prod_sum_int = m(mi)%prod_sum_int + production_temp
                          p_tot_int = p_tot_int + totalp*ice(sc,ic,mi)%bv(ii)*ice(sc,ic,mi)%th(ii)
                      endif

                      if (monte_carlo .eq. 1) then
                          production_temp = totalp*ice(sc,ic,mi)%bv(ii)*c_001*ice(sc,ic,mi)%th(ii)*c_001 ! (g)
                          mc_prod(mi) =  mc_prod(mi) + production_temp*lda_f
                      endif

                  endif


                  ! record new vars for updating main vars at end of lda loop,
                  ! dividing by the number of different snow depths to end up with an
                  ! average value
                  ! -------------------------------------------------------------

                  dSiO4_temp = 0.d0

                  smalg_new(ii) = smalg_new(ii)+smalg_temp*lda_f
                  poc_new(ii) = poc_new(ii) + poc_temp*lda_f
                  dN_all(ii) = dN_all(ii)+dN_all_temp*lda_f
                  !NO3_new(ii) = NO3_new(ii)+NO3_temp*lda_f  ! NO3 is now paritioned below, after lda loop
                  dNH4(ii) = dNH4(ii)+dNH4_temp*lda_f
                  dPO4(ii) = dPO4(ii)+dPO4_temp*lda_f
                  dSiOH4(ii) = dSiOH4(ii)+dSiO4_temp*lda_f

              enddo  ! end of ice layer loop


          ! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
          ! ------------------------------------------------------------
          ! END LIGHT ADJUSTMENT LOOP
          enddo

          ! update vars that were adjusted during the lda loop
          ice(sc,ic,mi)%smalg = smalg_new
          ice(sc,ic,mi)%poc(1:sk_z) = poc_new(1:sk_z)

          ! allocate new N concentrations to pools
          do ii=1,sk_z
              ice(sc,ic,mi)%PO4(ii) = ice(sc,ic,mi)%PO4(ii) + dPO4(ii)
              ice(sc,ic,mi)%SiOH4(ii) = ice(sc,ic,mi)%SiOH4(ii) + dSiOH4(ii)

              ! START - preferential uptake of nh4
              !if (ice(sc,ic,mi)%nh4(ii) .gt. 0.) then
              !    ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii) + dN_all(ii)
              !    if (ice(sc,ic,mi)%nh4(ii) .lt. 0.) then
              !        ice(sc,ic,mi)%no3(ii) = ice(sc,ic,mi)%no3(ii) + ice(sc,ic,mi)%nh4(ii)
              !        ice(sc,ic,mi)%nh4(ii) = 0.
              !    endif
              !else
              !    ice(sc,ic,mi)%no3(ii) = ice(sc,ic,mi)%no3(ii) + dN_all(ii) ! mMolN/m^3
              !endif
              !ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii) + dNH4(ii)
              ! END - preferential uptake of nh4

              ! START - proportional uptake of nh4
              tmp2 = (ice(sc,ic,mi)%no3(ii)+ice(sc,ic,mi)%nh4(ii))
              if (tmp2 .gt. 0.) then
                  tmp1 = ice(sc,ic,mi)%nh4(ii)/tmp2
              else
                  tmp1 = 0.
              endif
              ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii) + dN_all(ii)*tmp1
              ice(sc,ic,mi)%no3(ii) = ice(sc,ic,mi)%no3(ii) + dN_all(ii)*(1.-tmp1)
              ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii) + dNH4(ii)
              ! END - proportional uptake of nh4

              ! START - equal uptake of nh4
              !ice(sc,ic,mi)%no3(ii) = ice(sc,ic,mi)%no3(ii) + dN_all(ii)*c_5
              !ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii) + dN_all(ii)*c_5 ! mMolN/m^3
              !if (ice(sc,ic,mi)%no3(ii) .lt. 0.) then
              !      ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii) - ice(sc,ic,mi)%no3(ii)
              !      ice(sc,ic,mi)%no3(ii) = 0.
              !elseif (ice(sc,ic,mi)%nh4(ii) .lt. 0.) then
              !      ice(sc,ic,mi)%no3(ii) = ice(sc,ic,mi)%no3(ii) - ice(sc,ic,mi)%nh4(ii)
              !      ice(sc,ic,mi)%nh4(ii) = 0.
              !endif
              !ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii) + dNH4(ii)
              ! END - equal uptake of nh4

          enddo

      else ! case where there is no light...

          do ii = 1,sk_z

              if (ice(sc,ic,mi)%th(ii) .gt. 0.) then

                  death  = -1.*xp

                  ! convert death rate to mgC/m^3
                  death = ice(sc,ic,mi)%smalg(ii) - ice(sc,ic,mi)%smalg(ii)*exp(death*dtt_days) ! death/grazing
                  !tmp1 = explt(death*dtt_days,use_expt,exp_bins,exp_mul,expt)
                  !death = ice(sc,ic,mi)%smalg(ii) - ice(sc,ic,mi)%smalg(ii)*tmp1 ! death/grazing

                  ice(sc,ic,mi)%smalg(ii) = ice(sc,ic,mi)%smalg(ii) - death  ! mg/m^3 [within brine]

                  ! need bm_lost??

                  min_alg_brine = min_alg/(ice(sc,ic,mi)%bv(ii)*c_001)
                  if (ice(sc,ic,mi)%smalg(ii) .lt. min_alg_brine) then
                      ice(sc,ic,mi)%smalg(ii) = min_alg_brine  !algae concentrated into brine volume
                  else
                      ice(sc,ic,mi)%poc(ii) = ice(sc,ic,mi)%poc(ii) + death
                  endif

                  ! Remineralize
                  ! -------------------------------------------------------------
                  tmp1 = ice(sc,ic,mi)%poc(ii)*remin*dtt_days
                  ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii) + tmp1*mC_gC*n_c*remin_f  ! all happens in same volume, so don't need to take into account thickness or brine volume...
                  ice(sc,ic,mi)%po4(ii) = ice(sc,ic,mi)%po4(ii) + tmp1*mC_gC*p_c*remin_f  ! all happens in same volume, so don't need to take into account thickness or brine volume...
                  ice(sc,ic,mi)%poc(ii) = ice(sc,ic,mi)%poc(ii) - tmp1

              endif

              ice(sc,ic,mi)%Ik1(ii) = 0.
              ice(sc,ic,mi)%llim(ii) = 0.
              NO3_lim = ice(sc,ic,mi)%no3(ii) / (Ks_NO3 + ice(sc,ic,mi)%no3(ii))
              NH4_lim = ice(sc,ic,mi)%nh4(ii) / (Ks_NH4 + ice(sc,ic,mi)%nh4(ii))
              ice(sc,ic,mi)%nlim(ii) = max(NO3_lim,NH4_lim)
              ice(sc,ic,mi)%plim(ii) = ice(sc,ic,mi)%po4(ii) / (Ks_PO4 + ice(sc,ic,mi)%po4(ii))
              ice(sc,ic,mi)%silim(ii) = ice(sc,ic,mi)%sioh4(ii) / (Ks_SiO4 + ice(sc,ic,mi)%sioh4(ii))

              ice(sc,ic,mi)%silim(ii) = 1.d0

              if (ice(sc,ic,mi)%bs(ii) .gt. 103) then
                  sal_lim = 0.01
              else
                  sal_lim = a0 + a1*ice(sc,ic,mi)%bs(ii)  &
                      + a2*ice(sc,ic,mi)%bs(ii)**2 &
                      + a3*ice(sc,ic,mi)%bs(ii)**3
              endif
              ice(sc,ic,mi)%slim(ii) = sal_lim

          enddo

      endif    ! end of light/no light test

! ==============================================================================
! END GROW ALGAE DTT
! ==============================================================================


      ! START PLATELET LAYER CRAP
      ! ---------------------------------------------------------------------------------------
      ! ---------------------------------------------------------------------------------------

      if (use_pl .eq. 1) then
          if (pl%pur(1) .gt. 0.) then
              do ii=1,pl_max

                  ! use Ik' from last skeletal layer - I'm not saving 24h par values for platelet layer
                    if (ice(sc,ic,mi)%Ik1(sk_z) .le. 0.) then
                      light_lim = 0.
                  else
                      light_lim = 1. - exp(-1.*pl%pur(ii)/ice(sc,ic,mi)%Ik1(sk_z))
                  endif

                  NO3_lim = pl%no3(ii) / (Ks_NO3 + pl%no3(ii))
                  NH4_lim = pl%nh4(ii) / (Ks_NH4 + pl%nh4(ii))
                  ! find total N limitation term, as the last limiting form of N
                  NO3_lim = max(NO3_lim,NH4_lim)
                  PO4_lim = pl%po4(ii) / (Ks_PO4 + pl%po4(ii))
                  SiO4_lim = pl%sioh4(ii) / (Ks_SiO4 + pl%sioh4(ii))

                  ! find salinity limitation: Arrigo (199?)
                  if (pl%bs .gt. 103) then
                      sal_lim = 0.01
                  else
                      sal_lim = a0 + a1*pl%bs  &
                          + a2*pl%bs**2 &
                          + a3*pl%bs**3
                  endif
                  ! Determine which factor (nutrient or light) is limiting
                  ! -------------------------------------------------------------
                  min_lim = min(NO3_lim, PO4_lim, SiO4_lim, light_lim)

                  ! growth rate
                  ! -------------------------------------------------------------
                  Gmax = G0*exp(rg*pl%t)*min_lim*sal_lim !&

                  ! Grow
                  ! -------------------------------------------------------------
                  newp = pl%smalg(ii)*exp((Gmax - resp)*dtt_days) - pl%smalg(ii)

                  if (growth .gt. 0.) then
                      ! Calc new growth, production, nut. requirements, use up nutrients,
                      ! don't go negative
                      ! -------------------------------------------------------------
                      newp = newp*mC_gC  ! mMolC/m^3

                      ! find potential production from available nutrients
                      ! -------------------------------------------------------------
                      newpN = (pl%no3(ii)+pl%nh4(ii))*c_n  ! mMolN/m^3 * mMolC/mMolN = mMolC/m^3
                      newpP = pl%po4(ii)*c_p
                      newpSi = pl%sioh4(ii)*c_si
                      newp = min(newp,newpN,newpP,newpSi)

                      ! use up N
                      dN_all_temp = -1.*newp*n_c
                      ! use up DIP, Silica
                      dSiO4_temp = -1.*newp*si_c ! mMolSi/m^3
                      dPO4_temp = -1.*newp*p_c ! mMolP/m^3

                      newp = newp*gC_mC ! mg/m^3 [within brine]
                  else
                      dN_all_temp = 0.
                      dSiO4_temp = 0.
                      dPO4_temp = 0.
                  endif
                  dNH4_temp = 0.  ! default is zero - no remin

                  ! add to standing biomass, detritus
                  ! -------------------------------------------------------------
                  smalg_temp = pl%smalg(ii) + newp ! new standing biomass before grazing/death
                  death = smalg_temp - smalg_temp*exp(-1.*xp*dtt_days) ! death/grazing calucalted using new biomass ?
                  smalg_temp = smalg_temp - death ! die
                  if (newp .lt. 0.) then
                     death = death - newp ! if respiration is causing net death, include in death term for remineralization
                  endif

                  ! Remineralize
                  ! -------------------------------------------------------------
                  poc_temp = pl%poc(ii)
                  tmp1 = poc_temp*remin*dtt_days
                  dNH4_temp = tmp1*mC_gC*n_c*remin_f  ! all happens in same volume, so don't need to take into account thickness or brine volume...
                  dPO4_temp = dPO4_temp + tmp1*mC_gC*p_c*remin_f  ! all happens in same volume, so don't need to take into account thickness or brine volume...
                  poc_temp = poc_temp - tmp1

                  ! make sure algal concentration does not fall below minimum
                  min_alg_brine = min_alg/(pl%bv*c_001)
                  if (smalg_temp .lt. min_alg_brine) then
                      smalg_temp = min_alg_brine  !algae concentrated into brine volume
                  else
                      ! can add to poc, if not at algal minimum. Have to check because if we added
                      ! from death at algal minimum, poc would be created out of nothing
                      poc_temp = poc_temp + death
                  endif

                  ! record new vars for updating main vars at end of lda loop,
                  ! dividing by the number of different snow depths to end up with an
                  ! average value
                  ! -------------------------------------------------------------
                  pl%smalg(ii) = smalg_temp
                  pl%poc(ii) = poc_temp
                  pl%PO4(ii) = pl%PO4(ii) + dPO4_temp
                  pl%SiOH4(ii) = pl%SiOH4(ii) + dSiO4_temp
                  ! preferential update of nh4
                  if (pl%nh4(ii) .gt. 0.) then
                      pl%nh4(ii) = pl%nh4(ii) + dN_all_temp
                      if (pl%nh4(ii) .lt. 0.) then
                          pl%no3(ii) = pl%no3(ii) + pl%nh4(ii)
                          pl%nh4(ii) = 0.
                      endif
                  else
                      pl%no3(ii) = pl%no3(ii) + dN_all_temp ! mMolN/m^3
                  endif
                  pl%nh4(ii) = pl%nh4(ii) + dNH4_temp

             enddo
         else

             ! no light in platelet layer
              do ii = 1,pl_max

                  death  = -1.*xp

                  ! convert death rate to mgC/m^3
                  death = pl%smalg(ii) - pl%smalg(ii)*exp(death*dtt_days) ! death/grazing
                  pl%smalg(ii) = pl%smalg(ii) - death  ! mg/m^3 [within brine]

                  ! need bm_lost??

                  min_alg_brine = min_alg/(pl%bv*c_001)
                  if (pl%smalg(ii) .lt. min_alg_brine) then
                      pl%smalg(ii) = min_alg_brine  !algae concentrated into brine volume
                  else
                      pl%poc(ii) = pl%poc(ii) + death
                  endif

                  ! Remineralize
                  ! -------------------------------------------------------------
                  tmp1 = pl%poc(ii)*remin*dtt_days
                  pl%nh4(ii) = pl%nh4(ii) + tmp1*mC_gC*n_c*remin_f  ! all happens in same volume, so don't need to take into account thickness or brine volume...
                  pl%po4(ii) = pl%po4(ii) + tmp1*mC_gC*p_c*remin_f  ! all happens in same volume, so don't need to take into account thickness or brine volume...
                  pl%poc(ii) = pl%poc(ii) - tmp1

              enddo

          endif
      endif

      ! ---------------------------------------------------------------------------------------
      ! ---------------------------------------------------------------------------------------
      ! END PLATELET LAYER CRAP



                  do ii=1,sk_z
                      ! return to bulk conc.
                      bv_mean = ice(sc,ic,mi)%bv(ii)*c_001
                      ice(sc,ic,mi)%no3(ii) = ice(sc,ic,mi)%no3(ii)*bv_mean
                      ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii)*bv_mean
                      ice(sc,ic,mi)%po4(ii) = ice(sc,ic,mi)%po4(ii)*bv_mean
                      ice(sc,ic,mi)%sioh4(ii) = ice(sc,ic,mi)%sioh4(ii)*bv_mean
                      ice(sc,ic,mi)%smalg(ii) = ice(sc,ic,mi)%smalg(ii)*bv_mean
                      ice(sc,ic,mi)%poc(ii) = ice(sc,ic,mi)%poc(ii)*bv_mean
                  enddo
              endif

