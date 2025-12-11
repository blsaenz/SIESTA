

!=======================================================================
!BOP
!
! %DESCRIPTION:
!
! Atmospheric boundary interface (stability based flux calculations)
!
! !REVISION HISTORY:
!  SVN:$Id: ice_atmo.F90 140 2008-07-25 20:15:53Z eclare $
!
! author: Elizabeth C. Hunke, LANL
!
! 2003: Vectorized by Clifford Chen (Fujitsu) and William Lipscomb
! 2004: Block structure added by William Lipscomb
! 2006: Converted to free source form (F90) by Elizabeth Hunke
!=======================================================================
!BOP
!
! !IROUTINE: atmo_boundary_layer - compute coefficients for atm-ice fluxes,
!                                  stress and Tref/Qref
!
! !INTERFACE:
!
subroutine sia2_env_atmo(Tsf,potT,wind,zlvl,Qa,rhoa,lhcoef,shcoef)

      implicit none
      save

! %DESCRIPTION:
!
! Compute coefficients for atm/ice fluxes, stress, and reference
! temperature and humidity. NOTE: \\
! (1) all fluxes are positive downward,  \\
! (2) here, tstar = (WT)/U*, and qstar = (WQ)/U*,  \\
! (3) wind speeds should all be above a minimum speed (eg. 1.0 m/s). \\
!
! ASSUME:
!  The saturation humidity of air at T(K): qsat(T)  (kg/m**3)
!
! Code originally based on CSM1
!
! %REVISION HISTORY: same as module
!
! %USES:
!
! %INPUT/OUTPUT PARAMETERS:
!

       double precision, intent(in) :: &
          Tsf      , & ! surface temperature of ice or ocean (degC)
          potT     , & ! air potential temperature  (K)
          wind     , & ! wind speed (m/s)
          zlvl     , & ! atm level height (m)
          Qa       , & ! specific humidity (kg/kg)
          rhoa         ! air density (kg/m^3)

       double precision, intent(out) :: &
          shcoef   , & ! transfer coefficient for sensible heat
          lhcoef       ! transfer coefficient for latent heat
! !
! !EOP
! !
!
       integer :: &
          k

       double precision :: &
          TsfK  , & ! surface temperature in Kelvin (K)
          xqq   , & ! temporary variable
          psimh , & ! stability function at zlvl   (momentum)
          psimhs, & ! stable profile
          ssq   , & ! sat surface humidity     (kg/kg)
          qqq   , & ! for qsat, dqsfcdt
          TTT   , & ! for qsat, dqsfcdt
          qsat  , & ! the saturation humidity of air (kg/m^3)
          delt  , & ! potential T difference   (K)
          delq  , & ! humidity difference      (kg/kg)
          Lheat     ! Lvap or Lsub, depending on surface type

       double precision :: &
          ustar , & ! ustar (m/s)
          tstar , & ! tstar
          qstar , & ! qstar
          rdn   , & ! sqrt of neutral exchange coefficient (momentum)
          rhn   , & ! sqrt of neutral exchange coefficient (heat)
          ren   , & ! sqrt of neutral exchange coefficient (water)
          rd    , & ! sqrt of exchange coefficient (momentum)
          re    , & ! sqrt of exchange coefficient (water)
          rh    , & ! sqrt of exchange coefficient (heat)
          vmag  , & ! surface wind magnitude   (m/s)
          alz   , & ! ln(zlvl  /z10)
          thva  , & ! virtual temperature      (K)
          cp    , & ! specific heat of moist air
          hol   , & ! H (at zlvl  ) over L
          stable, & ! stability factor
          psixh     ! stability function at zlvl   (heat and water)

         ! from CICE 4 ice_constants.F90
      double precision, parameter :: &
         c0 = 0.d0, &
         c1 = 1.d0, &
         c2 = 2.d0, &
         c8 = 8.d0, &
         c10 = 10.d0, &
         c16 = 16.d0, &
         p5 = 0.5d0, &
         iceruf   = 0.0005d0   ,&! ice surface roughness (m)
         zref   = 10.d0   ,&! reference height for stability (m)
         qqqice  = 11637800.d0   ,&! for qsat over ice
         TTTice  = 5897.80d0      ,&! for qsat over ice
         Lsub      = 2.835d6 ,&! latent heat, sublimation freshwater (J/kg)
         Lvap      = 2.501d6 ,&! latent heat, vaporization freshwater (J/kg)
         zvir      = 0.606d0   ,&! rh2o/rair - 1.0
         R_air = 285.05d0, & ! J/kg/K
         R_h2o = 461.495d0, &
         vonkar    = 0.4d0     ,&! von Karman constant
         cp_wv     = 1.81d3  ,&! specific heat of water vapor (J/kg/K)
         cp_air    = 1005.0d0  ,&! specific heat of air (J/kg/K)
         ! (Briegleb JGR 97 11475-11485  July 1992)
         Tffresh   = 273.15d0  ,&! freezing temp of fresh ice (K)
         gravit    = 9.80616d0 , &
         cpvir = cp_wv/cp_air-c1 , & ! defined as cp_wv/cp_air - 1.
         umin  = c1, &  !             , &    ! minimum wind speed (m/s)
         pih = 0.5d0*3.14159265358979323846d0

       ! local functions
       double precision :: &
          xd    , & ! dummy argument
          psimhu, & ! unstable part of psimh
          psixhu    ! unstable part of psimx

      !------------------------------------------------------------
      ! Define functions
      !------------------------------------------------------------

       psimhu(xd)  = log((c1+xd*(c2+xd))*(c1+xd*xd)/c8) &
                   - c2*atan(xd) + pih
 !ech                  - c2*atan(xd) + 1.571_dbl_kind

       psixhu(xd)  =  c2 * log((c1 + xd*xd)/c2)

      !------------------------------------------------------------
      ! Initialize
      !------------------------------------------------------------

!          print *,'------------------------------------------------------------'
!          print *,'Tsf: ',Tsf      ! surface temperature of ice or ocean (degC)
!          print *,'potT: ', potT    ! air potential temperature  (K)
!          print *,'wind: ', wind    ! wind speed (m/s)
!          print *,'zlvl: ', zlvl   ! atm level height (m)
!          print *,'Qa: ', Qa      ! specific humidity (kg/kg)
!          print *,'rhoa: ', rhoa         ! air density (kg/m^3)



         shcoef = c0
         lhcoef = c0

      !------------------------------------------------------------
      ! define some more needed variables
      !------------------------------------------------------------

         qqq  = qqqice          ! for qsat
         TTT  = TTTice          ! for qsat
         if (Tsf .ge. c0) then
             Lheat = Lvap           ! water to vapor
         else
             Lheat = Lsub           ! ice to vapor
         endif
         vmag = max(umin, wind)
         rdn  = vonkar/log(zref/iceruf) ! neutral coefficient


         TsfK       = Tsf + Tffresh     ! surface temp (K)
         qsat       = qqq * exp(-TTT/TsfK)   ! saturation humidity (kg/m^3)
         ssq        = qsat / rhoa       ! sat surf hum (kg/kg)

         thva   = potT * (c1 + zvir * Qa) ! virtual pot temp (K)
         delt  = potT - TsfK       ! pot temp diff (K)
         delq  = Qa - ssq          ! spec hum dif (kg/kg)
         alz    = log(zlvl/zref)
         cp     = cp_air*(c1 + cpvir*ssq)

      !------------------------------------------------------------
      ! first estimate of Z/L and ustar, tstar and qstar
      !------------------------------------------------------------

         ! neutral coefficients, z/L = 0.0
         rhn = rdn
         ren = rdn

         ! ustar,tstar,qstar
         ustar = rdn * vmag
         tstar = rhn * delt
         qstar = ren * delq

      !------------------------------------------------------------
      ! iterate to converge on Z/L, ustar, tstar and qstar
      !------------------------------------------------------------

      do k=1,5

        ! compute stability & evaluate all stability functions
            hol = vonkar * gravit * zlvl &
                   * (tstar/thva &
                    + qstar/(c1/zvir+Qa)) &
                   / ustar**2

            hol    = sign( min(abs(hol),c10), hol )
            !if (hol >= 0)
            !    hol = min(abs(hol),c10);
            !else
            !    hol = -1.0*min(abs(hol),c10);
            !end

            stable = p5 + sign(p5 , hol)
            !if (hol >= 0)
            !    stable = p5 + p5;
            !else
            !    stable = p5 - p5;
            !end

            xqq    = max(sqrt(abs(c1 - c16*hol)) , c1)
            xqq    = sqrt(xqq)

            ! Jordan et al 1999
            psimhs = -(0.7d0*hol &
                   + 0.75d0*(hol-14.3d0) &
                   * exp(-0.35d0*hol) + 10.7d0)
            psimh  = psimhs*stable &
                    + (c1 - stable)*psimhu(xqq)
            psixh  = psimhs*stable &
                    + (c1 - stable)*psixhu(xqq)

        ! shift all coeffs to measurement height and stability
            rd = rdn / (c1+rdn/vonkar*(alz-psimh))
            rh = rhn / (c1+rhn/vonkar*(alz-psixh))
            re = ren / (c1+ren/vonkar*(alz-psixh))

        ! update ustar, tstar, qstar using updated, shifted coeffs
            ustar = rd * vmag
            tstar = rh * delt
            qstar = re * delq

      enddo                     ! end iteration

      !------------------------------------------------------------
      ! coefficients for turbulent flux calculation
      !------------------------------------------------------------
      ! add windless coefficient for sensible heat flux
      ! as in Jordan et al (JGR, 1999)
      !------------------------------------------------------------

         shcoef = rhoa * ustar * cp * rh + c1
         lhcoef = rhoa * ustar * Lheat  * re


      end subroutine sia2_env_atmo



!=======================================================================
!BOP
!
! !IROUTINE: atmo_boundary_const - compute coeeficients for atm-ice fluxes
!
!
! !INTERFACE:
!
      subroutine sia2_env_atmo_const(Tsf,wind,rhoa,lhcoef,shcoef)

! !DESCRIPTION:
!
! Compute coefficients for atm/ice fluxes, stress
! NOTE: \\
! (1) all fluxes are positive downward,  \\
! (2) reference temperature and humidity are NOT computed
!
! !REVISION HISTORY: same as module
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:

      double precision, intent(in) :: &
         Tsf      , & ! surface temperature in Kelvin (K)
         wind     , & ! wind speed (m/s)
         rhoa         ! air density (kg/m^3)

      double precision, intent(out):: &
         shcoef   , & ! transfer coefficient for sensible heat
         lhcoef       ! transfer coefficient for latent heat

! internal vars

      double precision :: &
         Lheat   ! Lvap or Lsub, depending on surface type

      ! from CICE 4 ice_constants.F90
      double precision, parameter :: &
         cp_air    = 1005.0d0  ,& ! specific heat of air (J/kg/K)
         Lsub      = 2.835d6 ,&   ! latent heat, sublimation freshwater (J/kg)
         Lvap      = 2.501d6      ! latent heat, vaporization freshwater (J/kg)

      !------------------------------------------------------------
      ! Initialize
      !------------------------------------------------------------

         shcoef = 0.d0
         lhcoef = 0.d0

      !------------------------------------------------------------
      ! define variables that depend on surface type
      !------------------------------------------------------------
         if (Tsf .ge. 0.d0) then
             Lheat = Lvap           ! water to vapor
         else
             Lheat = Lsub           ! ice to vapor
         endif

      !------------------------------------------------------------
      ! coefficients for turbulent flux calculation
      !------------------------------------------------------------

         shcoef = (1.20d-3)*cp_air*rhoa*wind
         lhcoef = (1.50d-3)*Lheat*rhoa*wind

      end subroutine sia2_env_atmo_const


