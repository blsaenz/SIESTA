
      ! areal participation calculation
      ! ======================================================================

      ! zero out advection areas
      adv(mi)%rf = 0.d0 !vector assignment
      ff = 0.d0  !vector assignment
      af_ic = 0.d0
      id_ic = 0.d0
      do ic=1,ida_n
          do sc=1,sda_n
              if (ice(sc,ic,mi)%af .gt. 0.d0) then
              af_ic(ic) = af_ic(ic) + ice(sc,ic,mi)%af
              id_ic(ic) = id_ic(ic) + ice(sc,ic,mi)%id(ice(sc,ic,mi)%z-z_sk)*ice(sc,ic,mi)%af
              endif
          enddo
          if (af_ic(ic) .gt. 0.d0) then
             id_ic(ic) = id_ic(ic)/af_ic(ic)
          endif
      enddo

      if (ida_n .gt. 1) then

          tmp1 = 0.  ! temporary var for tracking total ice area in categories
          do ic=1,ida_n-1

              if (af_ic(ic) .gt. 0.) then

                  ! distribute ridged area of input ice category into other ice categories
                  r_depth = id_ic(ic)
                  hm = 2.0d0*r_depth
                  r_adv_lambda = 1./(sqrt(r_depth)*adv_mu)
                  tmp2 = 0.  ! store total fraction, since redistribution is not precise

                  ! find depth of starting category
                  if (af_ic(1) .gt. 0.d0) then
                      r_depth = id_ic(1)
                  else
                      r_depth = ic_h_med(1) ! maybe this should be zero, or h_min, but we don't model zero ice?
                  endif

                  do ii=2,ida_n
                      ! figure out depths
                      r_depth_last = r_depth
                      if (af_ic(ii) .gt. 0.d0) then
                          r_depth = id_ic(ii)
                      else
                          r_depth = ic_h_med(ii)
                      endif
                      ! below is the fraction of rigding ice of category ic that ridges
                      ! into the ice category 'ii'
                      if (r_depth .ge. hm) then
                          adv(mi)%rf(ii,ic) = &
                              (exp(-1.*r_adv_lambda*(r_depth_last - hm)) - &
                               exp(-1.*r_adv_lambda*(r_depth - hm)))
                          tmp2 = tmp2 + adv(mi)%rf(ii,ic)
                      endif
                      if (ii .eq. ida_n .and. tmp2 .lt. 1.) then
                          ! add difference into largest category
                          adv(mi)%rf(ida_n,ic) = adv(mi)%rf(ida_n,ic) + (1.-tmp2)
                          tmp2 = 1.d0
                      endif
                  enddo

                  ! go back and correct rf as needed, and create ff for later optimization
                  do ii=1,ida_n
                      ! multiply by participation, correct imprecise distribution equation
                      adv(mi)%rf(ii,ic) = adv(mi)%rf(ii,ic)/tmp2
                      ! add to constant parts up as ff variable, for category areal optimization
!                      if (af_ic(ii) .gt. 0) then
!                          r_depth = id_ic(ii)
!                      else
                          r_depth = ic_h_med(ii)
!                      endif

    											void_f = cv_void_f
    											if (id_ic(ic) .lt. cvf_switch) then
    												void_f = void_f*cvf_thin
    											endif

                      ff(ic) = ff(ic) +  adv(mi)%rf(ii,ic)* &
                          id_ic(ic)/(r_depth*(1.d0-void_f))
                  enddo

              endif
          enddo

          ! ------------------------------------------------------------
          ! ------------------------------------------------------------
          ! reduce overall ridging effectiveness to account for void spaces
          ! this should make the code below allocate more total area for ridging
          ! to accomodate the reduced compaction

          !if (cv_void_f .gt. 0.d0) then
          !    ff = ff/(1.d0-cv_void_f)  ! vector operation
          !endif

          ! ------------------------------------------------------------
          ! ------------------------------------------------------------
          ! force the ridging fraction in the last category to 0.5, and
          ! add to totals (does not participate in void fraction shit)
          if (af_ic(ida_n) .gt. 0.) then
              adv(mi)%rf(ida_n,ida_n) = 1.d0
              ff(ida_n) = 0.5
          endif

          ! find areal convergence allocation
          ! use ff (the mutliplier which gives the area of ridged ice from a category after ridging)
          ! to optimize the total ridging area to produce the correct final area
          ! that accounts for the "observed" convergence
          ! ----------------------------------------------------------

          refine_dist = 1
          adv(mi)%a_rdg = 0.d0 ! vector assigment
          do jj=1,ida_n

              if (refine_dist .eq. 1) then
                  a_convg1 = adv(mi)%a_convg
                  if (jj .eq. 1) then
                      ! find initial convergence allocation between ice categories
                      ! ----------------------------------------------------------
                      do ic=1,ida_n
                          ! find fractional ice (open water area not used here)
                          tmp3 = af_ic(ic)/adv(mi)%af
                          tmp2 = tmp1 + tmp3
                          ! find participation of input ice category in ridging
                          apn(ic) = (exp((-1.)*tmp1*r_a_star) - &
                                     exp((-1.)*tmp2*r_a_star)) &
                                     *apn_denom
                          tmp1 = tmp2  ! recycle cumulatve fraction for next iteration

                          ! sum to find total areal convergence fraction
                          f_ridge_sum = f_ridge_sum + apn(ic)*ff(ic)

                      enddo
                  else
                      ! figure out remaining convergence area we need to account for
                      ! MATLAB: aconv_n = aconv - (sum(a_final(1:i)) - sum(a_final(1:i).*fff1(1:i)));
                      ! ----------------------------------------------------------
                      do  ii=1,jj-1
                          a_convg1 = a_convg1 - adv(mi)%a_rdg(ii)*ff(ii)
                      enddo
                      ! normalize remainaing apn values to a total of 1
                      tmp1 = 0.d0    ! apn sum
                      do ii=jj,ida_n
                          tmp1 = tmp1 + apn(ii)
                      enddo
                      do ii=jj,ida_n
                          apn(ii) = apn(ii)/tmp1
                      enddo
                  endif

                  ! find running total ridging area
                  tmp2 = 0.    ! ff*apn sum = mean ridging efficiency
                  do ii=jj,ida_n
                      tmp2 = tmp2 + apn(ii)*ff(ii)
                  enddo
                  a_rdg_tot = a_convg1/(1.d0-tmp2)

                  ! set refine_dist to 0
                  refine_dist = 0
              endif

              ! test to make sure ridging area does not exceed total ice
              ! category area, and correct if it does
              if (apn(jj)*a_rdg_tot .gt. cell_area*af_ic(jj)) then
                  if (jj .lt. ida_n) then
                      ! find next category that has ice/jj assigned
                      ii = jj+1
                      do while ((ff(ii) .eq. 0.) .and. (ii .lt. ida_n))
                          ii = ii+1
                      enddo
                      if (ii .eq. ida_n .and. ff(ii) .eq. 0.) then
                          print *,'Ridging scheme cannot accommodate this much ridging - cell ',mi
                          adv(mi)%a_rdg(jj) = cell_area*af_ic(jj)
                      else
                          tmp1 = cell_area*af_ic(jj)/a_rdg_tot*ff(jj)/ff(ii)
                          apn(ii) = apn(ii) + tmp1 ! add leftover psuedo-allocation to next larger ice category
                          apn(jj) = -1 ! mark this allocation as bad, just in case we try to use it again
                          adv(mi)%a_rdg(jj) = cell_area*af_ic(jj) ! select entire category for ridging
                          refine_dist = 1 ! need to refine the remaining distribution now, since this category was maxed
                      endif
                  else
                      a_convg1 = adv(mi)%a_convg
                      do ii=1,ida_n-1
                          a_convg1 = a_convg1 - adv(mi)%a_rdg(ii)
                      enddo
                      ff(ida_n) = 1.d0-a_convg1/(cell_area*af_ic(jj))
                      ! also need to fix biggest ice category ridging fraction somewhere else?
                      ! or maybe I can just use ff later? -- yes that's what is happening now
                      adv(mi)%a_rdg(jj) = cell_area*af_ic(jj)
                      convg_maxed = convg_maxed + 1
                      apn(jj) = -1 ! mark this allocation as bad, just in case we try to use it again
                  endif
              else
                  adv(mi)%a_rdg(jj) = a_rdg_tot*apn(jj)
              endif

          enddo ! end jj loop over ice category dist refinement

      else  ! only one 'ice category' (ida=1)

				! set ridging fraction to 1 - meaning that all ridging occurs from this category, into this category
				adv(mi)%rf = 1.d0

 				! set ff (ridging efficiency), here meaning ridged ice is 2x non-ridged
        ff(1) = 0.5 ! -- note cv_voids_f (voids fraction during ridging) is factored in during ridging

        apn(1) = 1.d0 ! mark allocation

				adv(mi)%a_rdg(1) = adv(mi)%a_convg/(1.d0-ff(1))

        if(adv(mi)%a_rdg(1) .gt. cell_area*af_ic(1)) then
					ff(1) = 1.d0-adv(mi)%a_convg/(cell_area*af_ic(1))
					! also need to fix biggest ice category ridging fraction somewhere else?
					! or maybe I can just use ff later? -- yes that's what is happening now
					adv(mi)%a_rdg(1) = cell_area*af_ic(1)
					convg_maxed = convg_maxed + 1
	        apn(1) = -1.d0 ! mark allocation as bad, just in case we try to use it again
				endif

      endif  ! end of 'are we using ice categories' test
