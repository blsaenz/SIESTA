
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




!      if (ice(sc,ic,mi)%pond%th .gt. 0. .and. vb_open .eq. 1) then
!      if (ice(sc,ic,mi)%pond%th .gt. 0.d0) then

          ! get drainage volume
          bv_total = th_pond_new*pond_f_perc

           ! record out fluxes of water and salt
           dbvdt = bv_total
           ii = sk_z
           do while(dbvdt .gt. 0.d0)
               if (dbvdt .gt. th_new(ii)) then
                  ! converting to kg/m^2 freshwater input 1000kg/m3 * m * fraction_water by weight
                  !m(mi)%fresh_out = m(mi)%fresh_out + &
                    !th_new(ii)*ice(sc,ic,mi)%af*1.d3*(1.d0-ice(sc,ic,mi)%s(ii)*c_001)
                  ! converting to kg/m^2 salt input 1000kg/m3 * m * fraction_salt by weight
                  !m(mi)%dsdt_out = m(mi)%dsdt_out + &
                   ! th_new(ii)*ice(sc,ic,mi)%af*ice(sc,ic,mi)%s(ii)

                  dbvdt = dbvdt - th_new(ii)
                  ii = ii-1
               else
                  ! converting to kg/m^2 freshwater input 1000kg/m3 * m * fraction_water by weight
                  !m(mi)%fresh_out = m(mi)%fresh_out + &
                  !  dbvdt*ice(sc,ic,mi)%af*1.d3*(1.d0-ice(sc,ic,mi)%s(ii)*c_001)
                  ! converting to kg/m^2 salt input 1000kg/m3 * m * fraction_salt by weight
                  !m(mi)%dsdt_out = m(mi)%dsdt_out + &
                  !  dbvdt*ice(sc,ic,mi)%af*ice(sc,ic,mi)%s(ii)
                  dbvdt = 0.d0
               endif
           enddo
           ! record pond fraction going srtaight to ocean
           ! converting to kg/m^2 freshwater input 1000kg/m3 * m * fraction_water by weight
           !m(mi)%fresh_out = m(mi)%fresh_out + &
           !  th_pond_new*(1.d0-pond_f_perc)* &
           !  ice(sc,ic,mi)%af*1.d3*(1.d0-ice(sc,ic,mi)%pond%s*c_001)
           ! converting to kg/m^2 salt input 1000kg/m3 * m * fraction_salt by weight
           !m(mi)%dsdt_out = m(mi)%dsdt_out + &
           !  th_pond_new*(1.d0-pond_f_perc)*ice(sc,ic,mi)%af*ice(sc,ic,mi)%pond%s


          dbvdt = bv_total ! reset dbvdt for different salinity accounting below - need indenpendent copies for salt and other tracers

          do ii = 1,sk_z

              ! percolation of salinity - different from other tracers b/c of minimum salinity
              ! ------------------------------------------------------------------------------

              if (ii .eq. 1) then
                  remapping = .FALSE.
              else
                  id_mean = new_layer_bot  ! store new_layer_bot for use with tracers below also
              endif
              dz_total = 0.d0
              s_new(ii) = 0.d0


              ! add min salinity to bv_total - otherwise flushing volume
              ! can be lost in min brine volume and neevr flush anything
              bv_mean = min(1.d0,min_sal/ice(sc,ic,mi)%s(ii)) ! percentage of brine vol due to min sal
              bv_mean = (1.d0 - bv_mean)*th_new(ii) ! effective brine volume for salinity purposes

              if (dbvdt .gt. 0.d0) then

                  ! find depth of drainage through brine channels of current layer
                  dbvdt = dbvdt - bv_mean
                  if (dbvdt .lt. 0.) then
                      dz = dbvdt + bv_mean
                      dbvdt = 0.
                  else
                      dz = bv_mean
                  endif

                  ! add to holding vars for salinity
                  s_new(ii) = s_new(ii) + ice(sc,ic,mi)%pond%s*dz

                  dz_total = dz
              endif

              if (dz_total .lt. th_new(ii)) then

                   ! scale dz_total back to th_new size someehow? no, fill rest of th_new

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

					   ! record remapped salinity
					   s_new(ii) = s_new(ii) + ice(sc,ic,mi)%bs(jj)*dz

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

                  ! save salinity changes for finishing mean during heat flux, and
                  ! later updating
                  bs_mean = s_new(ii)/th_new(ii)
                  dsdt1 = ice(sc,ic,mi)%s(ii)*bs_mean/ice(sc,ic,mi)%bs(ii) - ice(sc,ic,mi)%s(ii)
                  dsdt(ii) = dsdt(ii) + dsdt1

                  ! add to heat change from brine movement
  !                bd_mean = 1000000. + bs_mean*800.       ! find new brine density to calculate true heat change
  !                heat_mean = ice(sc,ic,mi)%bd(ii)*ice(sc,ic,mi)%t(ii)*cw ! old brine heat
  !                heat_total = heat_total/th_new(ii)  ! new brine heat
  !                dcheat(ii) = dcheat(ii) + (heat_total-heat_mean)*bv_mean/bd_mean

              endif


              ! percolation of tracers
              ! ------------------------------------------------------------------------------

              if (ii .eq. 1) then
                  remapping = .FALSE.
              else
                  new_layer_bot = id_mean
              endif
              dz_total = 0.d0
              d_new(ii) = 0.d0
              heat_total = 0.d0
              smalg_new(ii) = 0.d0
              poc_new(ii) = 0.d0
              no3_new(ii) = 0.d0
              nh4_new(ii) = 0.d0
              po4_new(ii) = 0.d0
              sioh4_new(ii) = 0.d0

              if (bv_total .gt. 0.d0) then

                  ! find depth of drainage through brine channels of current layer
                  bv_total = bv_total - th_new(ii)
                  if (bv_total .lt. 0.) then
                      dz = bv_total + th_new(ii)
                      bv_total = 0.
                  else
                      dz = th_new(ii)
                  endif

                  ! add to holding vars for tracers - in total conentration
                  heat_total = heat_total + ice(sc,ic,mi)%pond%heat*dz  ! density already included
!                  smalg_new(ii) = smalg_new(ii) + ice(sc,ic,mi)%pond%smalg*dz
                  poc_new(ii) = poc_new(ii) + ice(sc,ic,mi)%pond%poc*dz
                  no3_new(ii) = no3_new(ii) + ice(sc,ic,mi)%pond%no3*dz
                  nh4_new(ii) = nh4_new(ii) + ice(sc,ic,mi)%pond%nh4*dz
                  po4_new(ii) = po4_new(ii) + ice(sc,ic,mi)%pond%po4*dz
                  sioh4_new(ii) = sioh4_new(ii) + ice(sc,ic,mi)%pond%sioh4*dz

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
					   heat_total = heat_total + ice(sc,ic,mi)%t(jj)*ice(sc,ic,mi)%bd(jj)*cw*dz
					   no3_new(ii) = no3_new(ii) + ice(sc,ic,mi)%no3(jj)*dz
					   nh4_new(ii) = nh4_new(ii) + ice(sc,ic,mi)%nh4(jj)*dz
					   po4_new(ii) = po4_new(ii) + ice(sc,ic,mi)%po4(jj)*dz
					   sioh4_new(ii) = sioh4_new(ii) + ice(sc,ic,mi)%sioh4(jj)*dz
!					   smalg_new(ii) = smalg_new(ii) + ice(sc,ic,mi)%smalg(jj)*dz
					   poc_new(ii) = poc_new(ii) + ice(sc,ic,mi)%poc(jj)*dz
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

              ! update state variables
!              dsmalg(ii) = dsmalg(ii) + (smalg_new(ii)/th_new(ii) - ice(sc,ic,mi)%smalg(ii))
              dpoc(ii) = dpoc(ii) + (poc_new(ii)/th_new(ii) - ice(sc,ic,mi)%poc(ii))
              dno3(ii) = dno3(ii) + (no3_new(ii)/th_new(ii) - ice(sc,ic,mi)%no3(ii))
              dnh4(ii) = dnh4(ii) + (nh4_new(ii)/th_new(ii) - ice(sc,ic,mi)%nh4(ii))
              dpo4(ii) = dpo4(ii) + (po4_new(ii)/th_new(ii) - ice(sc,ic,mi)%po4(ii))
              dsioh4(ii) = dsioh4(ii) + (sioh4_new(ii)/th_new(ii) - ice(sc,ic,mi)%sioh4(ii))

          enddo

!      endif  ! end of use_pond .eq. 1 check

