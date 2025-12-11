! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! START: snow calc
! ----------------------------------------------------------------------
! File: SIA2_env_snow_calc.inc.f90
! Purpose: re-calculates internal snow characterstics based on changed temeprature
!
! Required:
! T_mean - new snow temp (degC)
! d_mean - current snow density
!
! Results:
! d_mean - snow density (g/m^3)
! heat_mean - heat required to warm snow to melting
! ----------------------------------------------------------------------

	   heat_mean = d_mean*(heat_snow0 - 0.2309*T_mean - &
		 0.0034*T_mean**2)

! ----------------------------------------------------------------------
! End: snow calc
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

