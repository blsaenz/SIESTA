
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




!      if (ice(ic,mi)%pond%th .gt. 0. .and. vb_open .eq. 1) then
      if (ice(ic,mi)%pond%th .gt. 0.) then

          ! determine total drainage volume (height) using freeboard height
          dbvdt = ice(ic,mi)%pond%th - (ice(ic,mi)%pond%th + ice(ic,mi)%fbh) * &
              exp( -1.*dtt_s*grav*ice(ic,mi)%pond%d*min_perm / (visc*20.) )  ! m/s

          ! units:    m * exp( s * m/s^2 * kg/m^3 * m^2 * 1/m * ms/kg) = m


          ! determine total drainage volume (height) using freeboard height
!          dbvdt = (ice(ic,mi)%pond%th) - (ice(ic,mi)%pond%th) * &
!              exp( -1.*dtt_s*grav*ice(ic,mi)%pond%d*min_perm / (visc*ice(ic,mi)%id(sk_z)) )  ! m/s

          ! can't drain more than is currently in pond
          dbvdt = min(dbvdt,ice(ic,mi)%pond%th) ! m = m^3 (b/c layer is 100% brine)

          ! scale by horizontal trasport ratio
          bv_total = dbvdt*(1. - hbd_f)

          remapping = .FALSE.

          do ii = 1,sk_z

              dz_total = 0.
              s_new(ii) = 0.
              d_new(ii) = 0.
              heat_total = 0.
              smalg_new(ii) = 0.
              poc_new(ii) = 0.
              no3_new(ii) = 0.
              nh4_new(ii) = 0.
              po4_new(ii) = 0.
              sioh4_new(ii) = 0.

              if (bv_total .gt. 0.) then

                  ! find depth of drainage through brine channels of current layer
                  bv_total = bv_total - th_new(ii)
                  if (bv_total .lt. 0.) then
                      dz = bv_total + th_new(ii)
                      bv_total = 0.
                  else
                      dz = th_new(ii)
                  endif

                  ! add to holding vars for tracers - in total conentration
                  s_new(ii) = s_new(ii) + ice(ic,mi)%pond%s*dz
                  heat_total = heat_total + ice(ic,mi)%pond%heat*dz  ! density already included
                  smalg_new(ii) = smalg_new(ii) + ice(ic,mi)%pond%smalg*dz
                  poc_new(ii) = poc_new(ii) + ice(ic,mi)%pond%poc*dz
                  no3_new(ii) = no3_new(ii) + ice(ic,mi)%pond%no3*dz
                  nh4_new(ii) = nh4_new(ii) + ice(ic,mi)%pond%nh4*dz
                  po4_new(ii) = po4_new(ii) + ice(ic,mi)%pond%po4*dz
                  sioh4_new(ii) = sioh4_new(ii) + ice(ic,mi)%pond%sioh4*dz

                  dz_total = dz
              endif

              if (dz_total .lt. th_new(ii)) then

                   if (.NOT. remapping) then
                       remapping = .TRUE.
                       new_layer_top = 0.
                   else
                       new_layer_top = new_layer_bot
                   endif
                   new_layer_bot = new_layer_top + th_new(ii) - dz_total

                   ! find 1st old layer that contains part of new layer, going down
                   jj = 0.
                   old_layer_bot = 0.
                   do while((new_layer_top .ge. old_layer_bot) .and. (jj .lt. sk_z))
                       jj=jj+1
                       ! old thickness...
                       old_layer_bot = old_layer_bot + th_new(jj)
                   enddo
                   old_layer_top = old_layer_bot - th_new(jj)
                   ! now jj is OLD layer where NEW layer ii starts...

                   ! find total heat/salt from multiple layers/partial layers that make up
                   ! the new layer
                   do while ((old_layer_top .lt. new_layer_bot) .and. (jj .le. sk_z))

					   ! ----- NEW LAYER GEOMETRIES ------------------------------
					   interim_top = max(old_layer_top,new_layer_top)
					   interim_bot = min(old_layer_bot,new_layer_bot)

					   z1 = interim_top - old_layer_top
					   z2 = interim_bot - old_layer_top
					   dz = z2-z1

					   ! record tracers
					   s_new(ii) = s_new(ii) + ice(ic,mi)%bs(jj)*dz
					   heat_total = heat_total + ice(ic,mi)%t(jj)*ice(ic,mi)%bd(jj)*cw*dz
					   no3_new(ii) = no3_new(ii) + ice(ic,mi)%no3(jj)*dz
					   nh4_new(ii) = nh4_new(ii) + ice(ic,mi)%nh4(jj)*dz
					   po4_new(ii) = po4_new(ii) + ice(ic,mi)%po4(jj)*dz
					   sioh4_new(ii) = sioh4_new(ii) + ice(ic,mi)%sioh4(jj)*dz
					   smalg_new(ii) = smalg_new(ii) + ice(ic,mi)%smalg(jj)*dz
					   poc_new(ii) = poc_new(ii) + ice(ic,mi)%poc(jj)*dz
					   ! ----- SETUP VARIABLES FOR NEXT PARTIAL LAYER ------------
					   ! keeping track of emerging new layer thickness, for
					   dz_total = dz_total + dz

                       jj=jj+1
                       if (jj .le. sk_z) then
                           ! find boundaries of next old layer
                           old_layer_top = old_layer_bot
                           old_layer_bot = old_layer_bot + th_new(jj)
                       endif

                   enddo

              endif

              if (ii .lt. sk_1) then

                  ! save salinity changes for finsing mean during heat flux, and
                  ! later updating
                  bv_mean = ice(ic,mi)%bv(ii)*c_001
                  bs_mean = s_new(ii)/th_new(ii)
                  dsdt1 = bs_mean*bv_mean - ice(ic,mi)%s(ii)
                  dsdt(ii) = dsdt(ii) + dsdt1

                  ! add to heat change from brine movement
                  bd_mean = 1000000. + bs_mean*800.       ! find new brine density to calculate true heat change
                  heat_mean = ice(ic,mi)%bd(ii)*ice(ic,mi)%t(ii)*cw ! old brine heat
                  heat_total = heat_total/th_new(ii)  ! new brine heat
                  dcheat(ii) = dcheat(ii) + (heat_total-heat_mean)*bv_mean/bd_mean

              endif

              ! update state variables
              dsmalg(ii) = dsmalg(ii) + (smalg_new(ii)/th_new(ii) - ice(ic,mi)%smalg(ii))
              dpoc(ii) = dpoc(ii) + (poc_new(ii)/th_new(ii) - ice(ic,mi)%poc(ii))
              dno3(ii) = dno3(ii) + (no3_new(ii)/th_new(ii) - ice(ic,mi)%no3(ii))
              dnh4(ii) = dnh4(ii) + (nh4_new(ii)/th_new(ii) - ice(ic,mi)%nh4(ii))
              dpo4(ii) = dpo4(ii) + (po4_new(ii)/th_new(ii) - ice(ic,mi)%po4(ii))
              dsioh4(ii) = dsioh4(ii) + (sioh4_new(ii)/th_new(ii) - ice(ic,mi)%sioh4(ii))

          enddo

          ! update pond depth
          ice(ic,mi)%pond%th = ice(ic,mi)%pond%th - dbvdt
          if (ice(ic,mi)%pond%th .le. 0.) then
              ice(ic,mi)%pond%t = 0.
              ice(ic,mi)%pond%s = 0.
              ice(ic,mi)%pond%d = 0.
              ice(ic,mi)%pond%heat = 0.
              ice(ic,mi)%pond%th = 0.
              ice(ic,mi)%pond%smalg = 0.     ! brine-based algal concentration (mgC/m^3)
              ice(ic,mi)%pond%poc = 0.     ! particulate organic carbon (detritus) (mgC/m^3)
              ice(ic,mi)%pond%no3 = 0.        ! ice layer brine NO3 concentration (µMol)
              ice(ic,mi)%pond%nh4 = 0.      ! ice layer brine NH4 concentration (µMol)
              ice(ic,mi)%pond%po4 = 0.     ! ice layer brine PO4 concentration (µMol)
              ice(ic,mi)%pond%sioh4 = 0.    ! ice layer brine SiOH4 concentration (µMol)
          endif


      endif  ! end of use_pond .eq. 1 check

