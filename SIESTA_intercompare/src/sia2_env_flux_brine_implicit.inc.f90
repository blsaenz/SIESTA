! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! START Backwards/Implicit Euler Nutrient Flux
! ----------------------------------------------------------------------
!
! Derivation of nutrient flux:
! ----------------------------
!    V = volume of brine (m^3/m^2) = brine_v*thickness = m^3/m^3*m
!    deltaV = brine flux volume (m^3/m^2) = flux*dt = m^3/m^2/s*s
!    N(z) = nutrient concentration @ level z (mMol/m^3)
!    N(z+1) = nutrient concentration @ level z+1 (mMol/m^3)
!    N(z)' = N(z) + dN(z)/dt

!    ((V - deltaVu - deltaVd)*N(z) + deltaVu*N(z+1) + deltaVd*N(z-1))/V = N(z)'
!    (V - deltaVu - deltaVd)*N(z) + deltaVu*N(z+1) + deltaVd*N(z-1) = V*N(z)'
!    N(z)*V - N(z)*deltaVu - N(z)*deltaVd + N(z+1)*deltaVu + N(z-1)*deltaVd = V*N(z)'
!    V*N(z)' =  N(z)*(V - deltaVu - deltaVd) + N(z+1)*(deltaVu) + N(z-1)*(deltaVd)
!    N(z)' =  N(z)*(1 - deltaVu/V - deltaVd/V) + N(z+1)*(deltaVu/V) + N(z-1)*(deltaVd/V)
!    rearanging into discretization of dN(z)/dt:
!    N(z)' - N(z) =  deltaVu/V*(N(z+1) - N(z)) + deltaVd/V*(N(z-1) - N(z))
!    making it implicit:
!    N(z)' - N(z) =  deltaVu/V*(N(z+1)' - N(z)') + deltaVd/V*(N(z-1)' - N(z)')
!    rearanging into tridiagonal format:
!    N(z)' - N(z) =  deltaVu/V*(N(z+1)' - N(z)') + deltaVd/V*(N(z-1)' - N(z)')
!    N(z) = N(z)'*(1 + deltaVu/V + deltaVd/V) - deltaVu/V*N(z+1)' - deltaVd/V*N(z-1)'
!    where:
!    deltaVu/V = flux_up*dt/(brine_v*thickness), deltaVd/V = flux_down*dt/(brine_v*thickness)
!
! Matrix layout & solution:
! --------------
!
! SUBROUTINE DGTSV( N, NRHS, DL, D, DU, B, LDB, INFO )
! tri-diangonal matrix solution to D*X = B, where
!
! D is matrix with diagonals DC (center) DU (upper) and DL(lower)
! B is the solution
! X is the result, which is returned in the variable B
! NRHS = dimension of solution, result, which for us = 1
! LDB = leading dimension of B array, which for us  = N
! INFO = and integer
!
! IMPUT VARS:
! ice(ic,mi)%bv - brine vol (ppt)
! ice(ic,mi)%NO3 - (mmol/m^3 in brine)
! ice(ic,mi)%NH4 - (mmol/m^3 in brine)
! ice(ic,mi)%PO4 - (mmol/m^3 in brine)
! ice(ic,mi)%SiOH4 - (mmol/m^3 in brine)
! flux_d - diffusive flux (m^3/m^2)
! flux_ml - molecular flux (m^3/m^2)
! flux_sk - skeletal flux, based on dhdt (m^3/m^2)
! dt_flux - time step length of flux (sec)
! hv - molecular diffusion layer (m)
! sk_z - number of bottom (skeletal) later
! vb_crit - brine volume above which flux can happen (ppt - constant)
! m_order - size of flux matrix
!
! OUTPUT VARS:
! ice(ic,mi)%NO3
! ice(ic,mi)%NH4
! ice(ic,mi)%PO4
! ice(ic,mi)%SiOH4
!
! INTERNAL VARS:
!
! integer :: ii,info
! double precision :: flux_down,th_down,bv_down,Fd,Fu,Fdiff
! double precision, dimension(z_max+1) :: dNO3,dNH4,dPO4,dSiOH4,DC,DC_calc
! double precision, dimension(z_max) :: DL,DL_calc,DU,DU_calc,r_vec
!
! ----------------------------------------------------------------------


       debug_z = ice(ic,mi)%no3(1:z_max)

	  ! return to brine conc.
	  do ii=1,sk_z
		  tmp1 = 1000./ice(ic,mi)%bv(ii)
		  ice(ic,mi)%no3(ii) = ice(ic,mi)%no3(ii)*tmp1
		  ice(ic,mi)%nh4(ii) = ice(ic,mi)%nh4(ii)*tmp1
		  ice(ic,mi)%po4(ii) = ice(ic,mi)%po4(ii)*tmp1
		  ice(ic,mi)%sioh4(ii) = ice(ic,mi)%sioh4(ii)*tmp1
		  ice(ic,mi)%smalg(ii) = ice(ic,mi)%smalg(ii)*tmp1
		  ice(ic,mi)%poc(ii) = ice(ic,mi)%poc(ii)*tmp1
	  enddo

	  ! boundary layer physics to determine supply of nutrients
	  ! to bottom of ice

	  ! ut = friction velocity = 8e-3 - 4e-3 m/s, depending on tide
	  ! v = kinematic viscosity of seawater = 1.83e-6 m^2/s
	  ! hv = molecular sub-layer thickness = v/ut (m)
	  !
	  ! Fn = flux of nutrients across "mixed layer" = dN/dz*K, where
	  ! dz = mixed layer = 2.5m
	  ! dN = N(top mixed layer/bottom moelecular sub-layer) - N(seawter)
	  ! K = eddy diffusion coeff. = 5e-4 to 1.5e-3, depending on tide
	  !
	  ! Fn_skel = flux of nutrients across molecular sub-layer = dN_o/hv*D, where
	  ! dN_o = N(skeletal layer) - N(bottom molecular sub-layer)
	  ! D = molecular diffusion coeff. of Nutrients = 1e-9 m^2/s


	  flux_ml = 1.e-3 ! midway between tidal extremes
      if (use_mdiff .eq. 1) then
	      flux_d = 4.2e-10 ! m^2/s for sioh4 at -2 deg C
	  else
	      flux_d = 0.
	  endif
	  !hv = 3.05e-4 ! (m) calculated with ut = 6e-3, mid-way between 8e-3 - 4e-3 m/s
	  hv = 1.e-3 ! set arbitrarily - diffusion is happening really quick!
	  fvb = 1.0e-6   ! cm^3/cm^s/s - Arrigo et al. 1993, from Reeburgh 1984
      !fvb = 0.

      ! viscosity of seawater (sal=35ppt) equation: visc = 0.0015*T^2 - 0.0608*T + 1.8284 m^2/s*10e6

	  ! if growing, skeletal transport can be huge, so we increase the
	  ! number of timesteps so things down't blow up
	  if (dhdt_cm_s .gt. 0) then
		  ! calculate skeletal transport from total dhdt
		  !fvs = -9.652e-9 + 4.49e-6*dhdt_cm_s - 1.39e-7*dhdt_cm_s**2  ! cm^3/cm^2/s - apparently bad
		  fvs = -9.652e-9 + 0.387*dhdt_cm_s - 1029.1*dhdt_cm_s**2  ! cm^3/cm^2/s - from arrigo basic bodel
		  flux_sk = c_01*(ice(ic,mi)%bv(sk_1)*c_001)*fvs                       ! m^3/m^2/s
	  else
		  flux_sk = 0.
	  endif

	  ! set size of matrix for flux solution
	  m_order = sk_z+1

	  ! initial calcs
      min_perm = 0.

	  do ii=1,sk_z+1
          if (ii .eq. sk_z+1) then
              jj = ml_z
              dheat(ii) = f(mi)%t*cw*f(mi)%d
          else
              jj = ii
          endif

          p_flux(ii) = 0.

	      ! Set center diagonal default
	      DC(ii) = 1.

		  ! linear eqn constants
		  dNO3(ii) = ice(ic,mi)%NO3(jj)
		  dNH4(ii) = ice(ic,mi)%NH4(jj)
		  dPO4(ii) = ice(ic,mi)%PO4(jj)
		  dSiOH4(ii) = ice(ic,mi)%SiOH4(jj)
		  dpoc(ii) = ice(ic,mi)%poc(jj)
!		  dheat(ii) = ice(ic,mi)%t(jj)*ice(ic,mi)%bd(jj)

	  enddo

      h_flux = dheat ! vector assignment

      ! these are calculated independent of layer
      Fdiff = dt_flux*flux_d/(hv**2)
      Feddy = dt_flux*flux_ml/6.25  ! 6.25 = 2.5m mixed layer, squared
      !Feddy = dt_flux*2.3148e-6/(0.6**2) ! simulate 0.6m platelet layer


      do ii = 1,sk_z+1
          flux_down = d0_

		  if (ii .eq. sk_z) then
			  ! this layer fluxes to the molecular layer below (z+1) using both
			  ! volume and diffusive flux

		      ! calc matrix coefficients, check to make sure brine_v is high enough
			  if (ice(ic,mi)%bv(ii) .ge. vb_crit) then
                  ! find largest flux between brine convection, skeletal convection, and gravity drainage
				  flux_down = (1.e-5)*ice(ic,mi)%bv(ii)*fvb
                  ! find max flux
				  flux_down = max(flux_sk,flux_down)*dt_flux
				  flux_down = max(flux_down,b_flux(ii)) ! no time in b_flux - already taken care of it desal
				  Fd = flux_down/(ice(ic,mi)%th(ii)*ice(ic,mi)%bv(ii)*c_001)
			  else
                  Fd = 0.
			  endif

			  ! use the largest flux between diffusion and convection
			  if (Fd .gt. Fdiff) then
				  Fu = flux_down/2.5 ! flux sk with with entire 2.5m mixed layer
			  else
				  Fd = Fdiff
				  Fu = Fdiff
			  endif

			  ! linear eqn constant - only for explicit euler
			  !dNO3(ii) = dNO3(ii) - ice(ic,mi)%NO3(ii)*Fd
			  !dNH4(ii) = dNH4(ii) - ice(ic,mi)%NH4(ii)*Fd
			  !dPO4(ii) = dPO4(ii) - ice(ic,mi)%PO4(ii)*Fd
			  !dSiOH4(ii) = dSiOH4(ii) - ice(ic,mi)%SiOH4(ii)*Fd
			  !dpoc(ii) = dpoc(ii) - ice(ic,mi)%poc(ii)*Fd

              !if (boundary_file .eq. 1 .and. use_pl .eq. 1) then
				  !dNO3(ii+1) = pl%NO3(1) - pl%NO3(1)*Fu
				  !dNH4(ii+1) = pl%NH4(1) - pl%NH4(1)*Fu
				  !dPO4(ii+1) = pl%PO4(1) - pl%PO4(1)*Fu
				  !dSiOH4(ii+1) = pl%SiOH4(1) - pl%SiOH4(1)*Fu
				  !dpoc(ii+1) = pl%poc(1) - pl%poc(1)*Fu
              !else
				  !dNO3(ii+1) = dNO3(ii+1) - ice(ic,mi)%NO3(ml_z)*Fu
				  !dNH4(ii+1) = dNH4(ii+1) - ice(ic,mi)%NH4(ml_z)*Fu
				  !dPO4(ii+1) = dPO4(ii+1) - ice(ic,mi)%PO4(ml_z)*Fu
				  !dSiOH4(ii+1) = dSiOH4(ii+1) - ice(ic,mi)%SiOH4(ml_z)*Fu
 				  !dpoc(ii+1) = dpoc(ii+1) - ice(ic,mi)%poc(ml_z)*Fu
             !endif

              ! upper diagonal
			  DU(ii) = 	-1.*Fd ! flux from current layer DOWN
			  ! lower diagonal
			  DL(ii) =  -1.*Fu ! flux from layer below UP
              ! append to center diagonal - only implicit euler
              DC(ii) = DC(ii) - DU(ii)
              DC(ii+1) = DC(ii+1) - DL(ii)

		  elseif (ii .eq. sk_z+1) then
			  ! this molecular layer fluxes with the skeletal layer diffusively,
			  ! and using eddy diffusion through a 2.5m mixed layer to the bulk
			  ! seawater values

			  ! modify linear eqn constant, since seawater values below the mixed
			  ! layer are assumed to be constant

              ! linear eqn constant - only for explicit euler
              if (use_pl .eq. 1) then
                  tmp1 = 1. ! do nothing, placeholder for lazy coding
              else
				  ! linear eqn constant - only for explicit euler
				  !dNO3(ii) = dNO3(ii) + Feddy*(f(mi)%no3 - ice(ic,mi)%NO3(ml_z))
				  !dNH4(ii) = dNH4(ii) + Feddy*(f(mi)%nh4 - ice(ic,mi)%NH4(ml_z))
				  !dPO4(ii) = dPO4(ii) + Feddy*(f(mi)%po4 - ice(ic,mi)%PO4(ml_z))
				  !dSiOH4(ii) = dSiOH4(ii) + Feddy*(f(mi)%sioh4 - ice(ic,mi)%SiOH4(ml_z))
				  !dpoc(ii) = dpoc(ii) + Feddy*(f(mi)%poc - ice(ic,mi)%poc(ml_z))

				  ! linear eqn constant - only for implicit euler
				  dNO3(ii) = dNO3(ii) + Feddy*f(mi)%no3
				  dNH4(ii) = dNH4(ii) + Feddy*f(mi)%nh4
				  dPO4(ii) = dPO4(ii) + Feddy*f(mi)%po4
				  dSiOH4(ii) = dSiOH4(ii) + Feddy*f(mi)%sioh4
				  dpoc(ii) = dpoc(ii) + Feddy*f(mi)%poc
              endif

		  else

              ! check to make sure brine vol. is high enough to flux
			  if (ice(ic,mi)%bv(ii) .ge. vb_crit .and. ice(ic,mi)%bv(ii+1) .ge. vb_crit) then

                  ! find inter-layer geomentric mean brine volume, and flux distance
				  th_down = ice(ic,mi)%th(ii) + ice(ic,mi)%th(ii+1)
				  bv_down = ice(ic,mi)%bv(ii)*ice(ic,mi)%bv(ii+1)/ &
				    (ice(ic,mi)%th(ii)*ice(ic,mi)%bv(ii+1) + &
				     ice(ic,mi)%th(ii+1)*ice(ic,mi)%bv(ii))*th_down
				  th_down = th_down*c_5

				  ! determine standard brine flux
				  flux_down = (1.e-5)*bv_down*fvb  ! convert fvb to m^3/m^2/s, times brine fraction

                  ! if a skeletal layer, use skeletal flux if greater than brine flux
				  if (ii .ge. sk_1) then
					  flux_down = max(flux_sk,flux_down)
				  endif

                  ! if gravity drainage is faster, use that instead
				  flux_down = flux_down*dt_flux
				  flux_down = max(flux_down,b_flux(ii)) ! no time in b_flux - already taken care of it desal

				  ! calc matrix coefficients
				  Fd = flux_down/(ice(ic,mi)%th(ii)*ice(ic,mi)%bv(ii)*c_001)
                  Fu = flux_down/(ice(ic,mi)%th(ii+1)*ice(ic,mi)%bv(ii+1)*c_001)

				  ! linear eqn constant - only for explicit euler
				  !dNO3(ii) = dNO3(ii) - ice(ic,mi)%NO3(ii)*Fd
				  !dNH4(ii) = dNH4(ii) - ice(ic,mi)%NH4(ii)*Fd
				  !dPO4(ii) = dPO4(ii) - ice(ic,mi)%PO4(ii)*Fd
				  !dSiOH4(ii) = dSiOH4(ii) - ice(ic,mi)%SiOH4(ii)*Fd
				  !dpoc(ii) = dpoc(ii) - ice(ic,mi)%poc(ii)*Fd

				  !dNO3(ii+1) = dNO3(ii+1) - ice(ic,mi)%NO3(ii+1)*Fu
				  !dNH4(ii+1) = dNH4(ii+1) - ice(ic,mi)%NH4(ii+1)*Fu
				  !dPO4(ii+1) = dPO4(ii+1) - ice(ic,mi)%PO4(ii+1)*Fu
				  !dSiOH4(ii+1) = dSiOH4(ii+1) - ice(ic,mi)%SiOH4(ii+1)*Fu
				  !dpoc(ii+1) = dpoc(ii+1) - ice(ic,mi)%poc(ii+1)*Fu

				  ! upper diagonal
				  DU(ii) = 	-1.*Fd ! flux from current layer DOWN
				  ! lower diagonal
				  DL(ii) =  -1.*Fu ! flux from layer below UP
				  ! append to center diagonal - only for implict euler
				  DC(ii) = DC(ii) - DU(ii)
				  DC(ii+1) = DC(ii+1) - DL(ii)

			  else
				  ! lower layer diagonal
				  DL(ii) = 0.
				  DU(ii) = 0.
			  endif

		  endif

          if (ii .ne. sk_z+1) then
              ice(ic,mi)%fbv(ii) = ice(ic,mi)%fbv(ii) + flux_down
          endif

      enddo

      ! calc flux NO3
      DC_calc = DC
      DU_calc = DU
      DL_calc = DL
      call DGTSV(m_order,1,DL_calc,DC_calc,DU_calc,dNO3,m_order,info)
      ! calc flux NH4
      DC_calc = DC
      DU_calc = DU
      DL_calc = DL
      call DGTSV(m_order,1,DL_calc,DC_calc,DU_calc,dNH4,m_order,info)
      ! calc flux PO4
      DC_calc = DC
      DU_calc = DU
      DL_calc = DL
      call DGTSV(m_order,1,DL_calc,DC_calc,DU_calc,dPO4,m_order,info)
      ! calc flux SiO4
      DC_calc = DC
      DU_calc = DU
      DL_calc = DL
      call DGTSV(m_order,1,DL_calc,DC_calc,DU_calc,dSiOH4,m_order,info)
      ! calc flux poc
      DC_calc = DC
      DU_calc = DU
      DL_calc = DL
      call DGTSV(m_order,1,DL_calc,DC_calc,DU_calc,dpoc,m_order,info)
      ! calc flux heat
!      DC_calc = DC
!      DU_calc = DU
!      DL_calc = DL
!      call DGTSV(m_order,1,DL_calc,DC_calc,DU_calc,h_flux,m_order,info)

      ! update concentrations
	  do ii=1,sk_z+1
          if (ii .eq. sk_z+1) then
              jj = ml_z
          else
              jj = ii
!              dcheat(ii) = dcheat(ii) + (h_flux(ii) - dheat(ii))
          endif

		  ice(ic,mi)%NO3(jj) = dNO3(ii)
		  ice(ic,mi)%NH4(jj) = dNH4(ii)
		  ice(ic,mi)%PO4(jj) = dPO4(ii)
		  ice(ic,mi)%SiOH4(jj) = dSiOH4(ii)
		  ice(ic,mi)%poc(jj) = dpoc(ii)

	      ! test flux
		  if (ice(ic,mi)%NO3(jj) .lt. 0.) then
			  ice(ic,mi)%NO3(jj) = 0.
		  endif

		  if (ice(ic,mi)%NH4(jj) .lt. 0.) then
			  !print *, 'NH4 less than zero',ii,ice(ic,mi)%NH4(ii),i,j
			  ice(ic,mi)%NH4(jj) = 0.
		  endif

		  if (ice(ic,mi)%PO4(jj) .lt. 0.) then
			  ice(ic,mi)%PO4(jj) = 0.
		  endif

		  if (ice(ic,mi)%SiOH4(jj) .lt. 0.) then
			  ice(ic,mi)%SiOH4(jj) = 0.
		  endif

		  if (ice(ic,mi)%poc(jj) .lt. 0.) then
			  ice(ic,mi)%poc(jj) = 0.
		  endif

	  enddo

      ! platelet crap
      ! -------------------------------------------------------------------
      ! bypasses mixed layer...
     if (use_pl .eq. 1) then
          pl%no3(1) = dno3(sk_z+1)
          pl%nh4(1) = dnh4(sk_z+1)
          pl%po4(1) = dpo4(sk_z+1)
          pl%sioh4(1) = dsioh4(sk_z+1)
          pl%poc(1) = dpoc(sk_z+1)
          Feddy = c_001*dt_flux*2.1412e-6/(pl_th**2) ! simulate platelet layer
          do jj = 1,100
			  do ii=1,pl_max
				  if (ii .eq. 1) then
					 dno3(ii) = Feddy*(pl%no3(ii+1) - pl%no3(ii))
					 dnh4(ii) = Feddy*(pl%nh4(ii+1) - pl%nh4(ii))
					 dpo4(ii) = Feddy*(pl%po4(ii+1) - pl%po4(ii))
					 dsioh4(ii) = Feddy*(pl%sioh4(ii+1) - pl%sioh4(ii))
					 dpoc(ii) = Feddy*(pl%poc(ii+1) - pl%poc(ii))
				  elseif (ii .eq. pl_max) then
					 dno3(ii) = Feddy*(f(mi)%no3 - pl%no3(ii))
					 dnh4(ii) = Feddy*(f(mi)%nh4 - pl%nh4(ii))
					 dpo4(ii) = Feddy*(f(mi)%po4 - pl%po4(ii))
					 dsioh4(ii) = Feddy*(f(mi)%sioh4 - pl%sioh4(ii))
					 dpoc(ii) = Feddy*(f(mi)%poc - pl%poc(ii))
				  else
					 dno3(ii) = Feddy*(pl%no3(ii+1) - 2*pl%no3(ii) + pl%no3(ii-1))
					 dnh4(ii) = Feddy*(pl%nh4(ii+1) - 2*pl%nh4(ii) + pl%nh4(ii-1))
					 dpo4(ii) = Feddy*(pl%po4(ii+1) - 2*pl%po4(ii) + pl%po4(ii-1))
					 dsioh4(ii) = Feddy*(pl%sioh4(ii+1) - 2*pl%sioh4(ii) + pl%sioh4(ii-1))
					 dpoc(ii) = Feddy*(pl%poc(ii+1) - 2*pl%poc(ii) + pl%poc(ii-1))
				  endif
			  enddo

			  do ii=1,pl_max
				 pl%no3(ii) = pl%no3(ii) + dno3(ii)
				 pl%nh4(ii) = pl%nh4(ii) + dnh4(ii)
				 pl%po4(ii) = pl%po4(ii) + dpo4(ii)
				 pl%sioh4(ii) = pl%sioh4(ii) + dsioh4(ii)
				 pl%poc(ii) = pl%poc(ii) + dpoc(ii)
			  enddo
          enddo
      endif

	  !return to absolute volume concentrations
	  do ii=1,sk_z
		  bv_mean = ice(ic,mi)%bv(ii)*c_001
		  ice(ic,mi)%no3(ii) = ice(ic,mi)%no3(ii)*bv_mean
		  ice(ic,mi)%nh4(ii) = ice(ic,mi)%nh4(ii)*bv_mean
		  ice(ic,mi)%po4(ii) = ice(ic,mi)%po4(ii)*bv_mean
		  ice(ic,mi)%sioh4(ii) = ice(ic,mi)%sioh4(ii)*bv_mean
		  ice(ic,mi)%smalg(ii) = ice(ic,mi)%smalg(ii)*bv_mean
		  ice(ic,mi)%poc(ii) = ice(ic,mi)%poc(ii)*bv_mean
	  enddo


! ----------------------------------------------------------------------
! END Backwards/Implicit Euler Nutrient Flux
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
