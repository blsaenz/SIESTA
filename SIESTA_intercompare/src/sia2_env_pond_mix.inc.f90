! Melt pond height - Notz 2009 using Darcy's Law
! ----------------------------------------------------------------------
! z(t+1) = z(t)exp(g*d*perm/(µh)
! where
! z = pond height (m)
! g = gravity (m/s^2)
! d = liquid density (kg/m^3)
! perm = minimum ice permeability
! µ = fluid viscosity = 1.79  (kg/m/s)
! h = ice thickness (m)

! viscosity of seawater (sal=35ppt) equation: visc = 0.0015*T^2 - 0.0608*T + 1.8284 m^2/s*10e6


!Mix pelt pond water with brine in uppermost layers that are above bv_conv

          if (ice(ic,mi)%bv(1) .gt. bv_conv .and. ice(ic,mi)%pond%th .gt. 0.) then

              ! convert pond tracers to total mass
              th_pond_new = ice(ic,mi)%pond%th
              s_pond_new = ice(ic,mi)%pond%s*ice(ic,mi)%pond%th
              h_pond_new = ice(ic,mi)%pond%heat*ice(ic,mi)%pond%th
              poc_pond_new = ice(ic,mi)%pond%poc*ice(ic,mi)%pond%th
              no3_pond_new = ice(ic,mi)%pond%no3*ice(ic,mi)%pond%th
              nh4_pond_new = ice(ic,mi)%pond%nh4*ice(ic,mi)%pond%th
              po4_pond_new = ice(ic,mi)%pond%po4*ice(ic,mi)%pond%th
              sioh4_pond_new = ice(ic,mi)%pond%sioh4*ice(ic,mi)%pond%th
              bv_total = ice(ic,mi)%pond%th

              jj = 1
              do while (ice(ic,mi)%bv(jj) .gt. bv_conv)
                  bv_mean = ice(ic,mi)%th(jj)*ice(ic,mi)%bv(jj)*c_001
                  s_pond_new = s_pond_new + ice(ic,mi)%bs(jj)*bv_mean
                  h_pond_new = h_pond_new + ice(ic,mi)%t(jj)*cw*ice(ic,mi)%bd(jj)*bv_mean
                  poc_pond_new = poc_pond_new + ice(ic,mi)%poc(jj)*bv_mean
                  no3_pond_new = no3_pond_new + ice(ic,mi)%no3(jj)*bv_mean
                  nh4_pond_new = nh4_pond_new + ice(ic,mi)%nh4(jj)*bv_mean
                  po4_pond_new = po4_pond_new + ice(ic,mi)%po4(jj)*bv_mean
                  sioh4_pond_new = sioh4_pond_new + ice(ic,mi)%sioh4(jj)*bv_mean
                  bv_total = bv_total + bv_mean
                  jj = jj+1
              enddo

              ! save for update if heat integration is valid
              s_pond_new = s_pond_new/bv_total
              h_pond_new = h_pond_new/bv_total
              poc_pond_new = poc_pond_new/bv_total
              no3_pond_new = no3_pond_new/bv_total
              nh4_pond_new = nh4_pond_new/bv_total
              po4_pond_new = po4_pond_new/bv_total
              sioh4_pond_new = sioh4_pond_new/bv_total

              ! save tracer changes for layers convecting with melt layer
              do ii=1,jj-1
                  bv_mean = ice(ic,mi)%bv(ii)*c_001
                  dsdt(ii) = dsdt(ii) + (s_pond_new*bv_mean - ice(ic,mi)%s(ii))
                  dpoc(ii) = dpoc(ii) + (poc_pond_new - ice(ic,mi)%poc(ii))
                  dno3(ii) = dno3(ii) + (no3_pond_new - ice(ic,mi)%no3(ii))
                  dnh4(ii) = dnh4(ii) + (nh4_pond_new - ice(ic,mi)%nh4(ii))
                  dpo4(ii) = dpo4(ii) + (po4_pond_new - ice(ic,mi)%po4(ii))
                  dsioh4(ii) = dsioh4(ii) + (sioh4_pond_new - ice(ic,mi)%sioh4(ii))
              enddo

              ! set pond to drain - assuming draining instantly if bv_conv is
              ! reached in top ice layer
              th_pond_new = 0.

          endif

! 2. Drain pond water laterally only according to Darcy's law using the minimum
!    permeability above freeboard.  Don't move into ice that much.
! ----------------------------------------------------------------------


          ! currently we are draining instantly
          ! if there is no snow, assume it just goes away by
          ! pushing some of the mixed brine out somewhere
          if (ice(ic,mi)%snow%z .eq. 0) then
              th_pond_new = 0.
          endif

          ! start of a lateral drainage scheme - never tested


!          tmp1 = 0.
!          jj = 1
!          min_perm = 1.e10.
!          bv_open = 100
!          do while (tmp1 .lt. ice(ic,mi)%fbh)
!              min_perm = min(min_perm,perm(jj))
!              if(ice(ic,mi)%bv(jj) .gt. vb_crit .and. bv_open .gt. 0) bv_open = jj
!              tmp1 = ice(ic,mi)%id(jj)
!          enddo

          ! determine total drainage volume (height) using freeboard height
!          dbvdt = (ice(ic,mi)%pond%th) * &
!              exp( dtt_s*grav*ice(ic,mi)%pond%d*c_001*min_perm / (visc*ice(ic,mi)%id(sk_z)) )  ! m/s

          ! can't drain more than is currently in pond
!          dbvdt = min(dbvdt,ice(ic,mi)%pond%th) ! m = m^3 (b/c layer is 100% brine)

