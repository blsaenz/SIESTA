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

              keep_growing = .true.
              m_order = 2;

              ! center diagonal - always 1 for implicit volume mxing
              DC(1) = 1.
              DC(2) = 1.

              ! linear eqn constants
              s_new(1) = ice(ic,mi)%pond%s
              NO3_new(1) = ice(ic,mi)%pond%no3
              NH4_new(1) = ice(ic,mi)%pond%nh4
              PO4_new(1) = ice(ic,mi)%pond%po4
              SiOH4_new(1) = ice(ic,mi)%pond%sioh4
              poc_new(1) = ice(ic,mi)%pond%poc
              dheat(1) = ice(ic,mi)%pond%heat
              s_new(2) = ice(ic,mi)%bs(1)
              NO3_new(2) = ice(ic,mi)%NO3(1)
              NH4_new(2) = ice(ic,mi)%NH4(1)
              PO4_new(2) = ice(ic,mi)%PO4(1)
              SiOH4_new(2) = ice(ic,mi)%SiOH4(1)
              poc_new(2) = ice(ic,mi)%poc(1)
              dheat(2) = ice(ic,mi)%t(1)*ice(ic,mi)%bd(1)*cw

              ! volume flux calculation & timestep, and drain amount
              ! of pond fluxed
              flux_down = ice(ic,mi)%pond%th*0.33
              tmp1 = z_th_min*bv_conv*c_001*0.33
              if (flux_down .gt. tmp1) then
                  flux_down = tmp1
                  th_pond_new = th_pond_new - tmp1
              else
                  th_pond_new = 0.
              endif
              jjj = 3

              ! calc matrix coefficients
              Fd = flux_down/(ice(ic,mi)%pond%th)
              Fu = flux_down/(ice(ic,mi)%th(1)*ice(ic,mi)%bv(1)*c_001)

              ! upper diagonal
              DU(1) = 	-1.*Fd ! flux from current layer DOWN
              ! lower diagonal
              DL(1) =  -1.*Fu ! flux from layer below UP
              ! append to center diagonal - only for implict euler
              DC(1) = DC(1) - DU(1)
              DC(2) = DC(2) - DL(1)

              do ii=2,sk_z

                  if (ice(ic,mi)%bv(ii) .gt. bv_conv .and. keep_growing) then

                      m_order = m_order + 1

                      ! Set center diagonal default
                      DC(m_order) = 1.

                      ! linear eqn constants
                      s_new(m_order) = ice(ic,mi)%bs(ii)
                      NO3_new(m_order) = ice(ic,mi)%NO3(ii)
                      NH4_new(m_order) = ice(ic,mi)%NH4(ii)
                      PO4_new(m_order) = ice(ic,mi)%PO4(ii)
                      SiOH4_new(m_order) = ice(ic,mi)%SiOH4(ii)
                      poc_new(m_order) = ice(ic,mi)%poc(ii)
                      dheat(m_order) = ice(ic,mi)%t(ii)*ice(ic,mi)%bd(ii)*cw

                      ! calc matrix coefficients
                      Fd = flux_down/(ice(ic,mi)%th(ii-1)*ice(ic,mi)%bv(ii-1)*c_001)
                      Fu = flux_down/(ice(ic,mi)%th(ii)*ice(ic,mi)%bv(ii)*c_001)

                      ! upper diagonal - flux from above layer DOWN
                      DU(ii) = 	-1.*Fd
                      ! lower diagonal - flux from current layer UP
                      DL(ii) =  -1.*Fu
                      ! subtract what's coming from above from above center diagonal
                      DC(ii) = DC(ii) - DU(ii)
                      ! subtract what's leaving here (to go up) from this center diagonal
                      DC(m_order) = DC(m_order) - DL(ii)

                  else

                      keep_growing = .false.

                  endif

              enddo

              do ii=1,jj
                  ! flux 'em
                  ! calc flux NO3
                  DC_calc = DC
                  DU_calc = DU
                  DL_calc = DL
                  call DGTSV(m_order,1,DL_calc,DC_calc,DU_calc,NO3_new,m_order,info)
                  ! calc flux NH4
                  DC_calc = DC
                  DU_calc = DU
                  DL_calc = DL
                  call DGTSV(m_order,1,DL_calc,DC_calc,DU_calc,NH4_new,m_order,info)
                  ! calc flux PO4
                  DC_calc = DC
                  DU_calc = DU
                  DL_calc = DL
                  call DGTSV(m_order,1,DL_calc,DC_calc,DU_calc,PO4_new,m_order,info)
                  ! calc flux SiO4
                  DC_calc = DC
                  DU_calc = DU
                  DL_calc = DL
                  call DGTSV(m_order,1,DL_calc,DC_calc,DU_calc,SiOH4_new,m_order,info)
                  ! calc flux poc
                  DC_calc = DC
                  DU_calc = DU
                  DL_calc = DL
                  call DGTSV(m_order,1,DL_calc,DC_calc,DU_calc,poc_new,m_order,info)
                  ! calc salt flux
                  DC_calc = DC
                  DU_calc = DU
                  DL_calc = DL
                  call DGTSV(m_order,1,DL_calc,DC_calc,DU_calc,s_new,m_order,info)
                  ! calc salt flux
                  DC_calc = DC
                  DU_calc = DU
                  DL_calc = DL
                  call DGTSV(m_order,1,DL_calc,DC_calc,DU_calc,dheat,m_order,info)
              enddo

              ! save pond info
              s_pond_new = s_new(1)
              h_pond_new = dheat(1)
              poc_pond_new = poc_new(1)
              no3_pond_new = no3_new(1)
              nh4_pond_new = NH4_new(1)
              po4_pond_new = PO4_new(1)
              sioh4_pond_new = SiOH4_new(1)

              ! record tracer changes for later consideration if integration is valid
              do ii=1,m_order-1

                  jj = ii + 1 ! flux result index
                  bv_mean = ice(ic,mi)%bv(ii)*c_001
                  dsdt(ii) = dsdt(ii) + (s_new(jj) - ice(ic,mi)%bs(ii))*bv_mean
                  dcheat(ii) = dcheat(ii) + (dheat(jj) - ice(ic,mi)%t(ii)*ice(ic,mi)%bd(ii)*cw)
                  dpoc(ii) = dpoc(ii) + (poc_new(jj) - ice(ic,mi)%poc(ii))
                  dno3(ii) = dno3(ii) + (no3_new(jj) - ice(ic,mi)%no3(ii))
                  dnh4(ii) = dnh4(ii) + (NH4_new(jj) - ice(ic,mi)%nh4(ii))
                  dpo4(ii) = dpo4(ii) + (PO4_new(jj) - ice(ic,mi)%po4(ii))
                  dsioh4(ii) = dsioh4(ii) + (SiOH4_new(jj) - ice(ic,mi)%sioh4(ii))

              enddo


          endif
