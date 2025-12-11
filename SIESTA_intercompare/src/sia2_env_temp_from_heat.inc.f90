! ---- SIA2_env_temp_from_heat.inc.f90 ---------------------------------
! this inline routine solves for a temp based on other layer parameters,
! including heat referenced to a eference temp.  It requires these variables to be set:
! heat_mean = the total heat in the layer in J/m^2
! d_mean = layer density (g/m^3)
! S_mean = layer salinity (psu)
! dz_mean = layer thickness (m)
! ----------------------------------------------------------------------

      tmp1 = heat_mean/dz_mean/d_mean
      tmp2 = s_mean*mu
      if (tmp1 .lt. tmp2*cw) then

		  aq = c0  ! quadratic a term
		  bq = (cw-c0)*tmp2 - Lf - tmp1 ! quadratic b term
		  cq = Lf*tmp2  ! quadradic c term

		  ! solve heat eqn. backwards to find T_mean
		  ! the "negative" solutions seem to always be right
		  T_mean = (-1.d0*bq - sqrt(bq**2 - 4.d0*aq*cq))/(2.d0*aq)

      else

          T_mean = tmp1/cw  ! T = Q/cw

      endif
