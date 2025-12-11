! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! START: brine calc
! ----------------------------------------------------------------------
! File: SIA2_env_brine_calc.inc.f90
! Purpose: inlude-able code to calculate brine-related paramaters from
! ice temp and salinity
!
! Required:
! T_mean - ice temp (degC)
! S_mean - bulk ice salinity (psu)
!
! Results:
! bs_mean - brine salinity (psu)
! bd_mean - brine density (g/m^3)
! d_mean - bulk ice density (g/m^3)
! bv_mean - brine vol (ppt)
! heat_mean - heat required to warm ice from T_ref to T_mean
! ----------------------------------------------------------------------

		  T_melt = mu*s_mean

          if (T_mean .lt. T_melt) then

			  ! calc new brine salinity
			  alpha0 = 1.d0
			  if (T_mean .gt. -3.058d0) then
				  ! use linear freezing point calc for brine salinity above temp=-3.058
				  bs_mean = T_mean/mu
			  else
				  ! use polynomial brine sal estimation for temps colder than -3.058
				  if (T_mean .le. -3.058d0 .and. T_mean .ge. -22.9d0) then
					  alpha0 = -3.9921d0
					  alpha1 = -22.7d0
					  alpha2 = -1.0015d0
					  alpha3 = -0.019956d0
				  elseif (T_mean .lt. -22.9d0 .and. T_mean .ge. -44.d0) then
					  alpha0 = 206.24d0
					  alpha1 = -1.8907d0
					  alpha2 = -0.060868d0
					  alpha3 = -0.0010247d0
				  elseif (T_mean .lt. -44.d0 .and. T_mean .ge. -54.d0) then
					   alpha0 = -4442.1d0
					   alpha1 = -277.86d0
					   alpha2 = -5.501d0
					   alpha3 = -0.03669d0
				  elseif (T_mean .lt. -54.d0) then
					   alpha0 = 0.d0
				  endif

				  if (alpha0 .eq. 0.d0) then
					  ! fix upper brine sal limit at 300 psu
					  bs_mean = 300.d0
				  else
					  ! use selected polynomial to calc brine sal
					  bs_mean = alpha0 + alpha1*T_mean + &
					  alpha2*T_mean**2 + alpha3*T_mean**3     ! ppt
				  endif
			  endif

			  ! error check brine sal
			  if (bs_mean .lt. 0.d0) then
				  bs_mean = 0.d0
			  endif

			  ! calculate new brine density (c=800 g m-3 ppt-1)
			  bd_mean = 1.d6 + bs_mean*800.d0       ! g/m^3

			  ! new ice density
			  d_mean = (1.d0-bb_f)*IceD*bd_mean*bs_mean/(bd_mean*bs_mean - S_mean* &
				  (bd_mean - IceD))

			  ! new brine volume again after desalination...
			  bv_mean = 1.d3*(d_mean*S_mean)/(bd_mean*bs_mean)  ! ppt

              ! d(bv)/dt
              !-(IceD*s*((1e6 + (800*t)/mu)/mu - (800*s)/mu + (800*t)/mu^2))/(s*(1e6 - IceD + (800*t)/mu) - (t*(1e6 + (800*t)/mu))/mu)^2
              !-(rho_i*s*((rho_w + (d_s*t)/mu)/mu - (d_s*s)/mu + (d_s*t)/mu^2))/(s*(rho_w - rho_i + (d_s*t)/mu) - (t*(rho_w + (d_s*t)/mu))/mu)^2

              ! d(enthalpy)/dt
              ! -(di*mu*s*(db^2*t^2 + 2*db*dw*mu*t - di*s*db*mu^2 + dw^2*mu^2))/(db*t^2 + di*mu^2*s - dw*mu^2*s + dw*mu*t - db*mu*s*t)^2
              ! -(rho_i*mu*c*(d_s^2*t^2 + 2*d_s*rho_w*mu*t - rho_i*c*d_s*mu^2 + rho_w^2*mu^2))/(d_s*t^2 + rho_i*mu^2*c - rho_w*mu^2*c + rho_w*mu*t - d_s*mu*c*t)^2
              ! -(rho_i*mu_Tm*cg*(d_s^2*t^2 + 2*d_s*rho_w*mu_Tm*t - rho_i*cg*d_s*mu_Tm^2 + rho_w^2*mu_Tm^2))/(d_s*t^2 + rho_i*mu_Tm^2*cg - rho_w*mu_Tm^2*cg + rho_w*mu_Tm*t - d_s*mu_Tm*cg*t)^2

			  heat_mean = -1.d0*d_mean*(c0*(T_melt - T_mean) &
				  + Lf*(1. - T_melt/T_mean) - cw*T_melt) ! j/m^3 - not j/m^2

          else

              bs_mean = s_mean
              bd_mean = 1.d6 + bs_mean*800.d0       ! g/m^3
              d_mean = bd_mean
              bv_mean = 1.d3
              heat_mean = cw*T_mean*d_mean

          endif




! ----------------------------------------------------------------------
! End: brine calc
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
