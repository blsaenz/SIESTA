

!=======================================================================
!BOP
!
! !IROUTINE: solution_dEdd - evaluate solution for Delta-Edddington solar
!
! !INTERFACE:
!
      subroutine sia2_env_edd_solution (coszen, srftyp, tau, w0, g, expt, &
             albodr, albodf, trndir, trntdr, trndif, rupdir, rupdif, &
             rdndif, nslyr, nilyr)
!
! !DESCRIPTION:
!
! Given input vertical profiles of optical properties, evaluate the
! monochromatic Delta-Eddington solution.
!
! !REVISION HISTORY:
!
! author:  Bruce P. Briegleb, NCAR
! updated: Dec 2009, Benjamin Saenz, Stanford
!
! !USES:
!
      use sia2_globals
      implicit none

! !INPUT/OUTPUT PARAMETERS:

      double precision,  &
         intent(in) :: &
         coszen      ! cosine solar zenith angle

      double precision,  &
         intent(in) :: &
         srftyp      ! surface type over ice: (0=air, 1=snow, 2=pond)

      integer, &
         intent(in) :: &
         nslyr,   &      ! number snow layers
         nilyr          ! number ice layers

      double precision, dimension(0:2*z_max), &
         intent(in) :: &
         tau     , & ! layer extinction optical depth
         w0      , & ! layer single scattering albedo
         g           ! layer asymmetry parameter

      double precision, dimension(2,exp_bins), &
         intent(in) :: &
         expt        ! exp() lookup table

      double precision,  &
         intent(in) :: &
         albodr  , & ! ocean albedo to direct rad
         albodf      ! ocean albedo to diffuse rad

      ! following arrays are defined at model interfaces; 0 is the top of the
      ! layer above the sea ice; klevp is the sea ice/ocean interface.
      double precision, dimension (0:2*z_max), &
         intent(out) :: &
         trndir  , & ! solar beam down transmission from top
         trntdr  , & ! total transmission to direct beam for layers above
         trndif  , & ! diffuse transmission to diffuse beam for layers above
         rupdir  , & ! reflectivity to direct radiation for layers below
         rupdif  , & ! reflectivity to diffuse radiation for layers below
         rdndif      ! reflectivity to diffuse radiation for layers above

!  INTERFACES

      INTERFACE
          FUNCTION explt(n,use_expt,exp_bins,exp_mul,expt)
              double precision :: n         ! exp argument
              integer :: use_expt           ! switch
              integer :: exp_bins           ! number of indexes to lookup table
              double precision :: exp_mul   ! muliplier to convert n to appropriate index
              double precision, dimension(2,exp_bins) :: expt
              double precision :: explt
              integer :: i
          end FUNCTION explt
      end INTERFACE


!
!EOP
!-----------------------------------------------------------------------
!
! Delta-Eddington solution for snow/air/pond over sea ice
!
! Generic solution for a snow/air/pond input column of klev+1 layers,
! with srftyp determining at what interface fresnel refraction occurs.
!
! Computes layer reflectivities and transmissivities, from the top down
! to the lowest interface using the Delta-Eddington solutions for each
! layer; combines layers from top down to lowest interface, and from the
! lowest interface (underlying ocean) up to the top of the column.
!
! Note that layer diffuse reflectivity and transmissivity are computed
! by integrating the direct over several gaussian angles. This is
! because the diffuse reflectivity expression sometimes is negative,
! but the direct reflectivity is always well-behaved. We assume isotropic
! radiation in the upward and downward hemispheres for this integration.
!
! Assumes monochromatic (spectrally uniform) properties across a band
! for the input optical parameters.
!
! If total transmission of the direct beam to the interface above a particular
! layer is less than trmin, then no further Delta-Eddington solutions are
! evaluated for layers below.
!
! The following describes how refraction is handled in the calculation.
!
! First, we assume that radiation is refracted when entering either
! sea ice at the base of the surface scattering layer, or water (i.e. melt
! pond); we assume that radiation does not refract when entering snow, nor
! upon entering sea ice from a melt pond, nor upon entering the underlying
! ocean from sea ice.
!
! To handle refraction, we define a "fresnel" layer, which physically
! is of neglible thickness and is non-absorbing, which can be combined to
! any sea ice layer or top of melt pond. The fresnel layer accounts for
! refraction of direct beam and associated reflection and transmission for
! solar radiation. A fresnel layer is combined with the top of a melt pond
! or to the surface scattering layer of sea ice if no melt pond lies over it.
!
! Some caution must be exercised for the fresnel layer, because any layer
! to which it is combined is no longer a homogeneous layer, as are all other
! individual layers. For all other layers for example, the direct and diffuse
! reflectivities/transmissivities (R/T) are the same for radiation above or
! below the layer. This is the meaning of homogeneous! But for the fresnel
! layer this is not so. Thus, the R/T for this layer must be distinguished
! for radiation above from that from radiation below. For generality, we
! treat all layers to be combined as inhomogeneous.
!
!-----------------------------------------------------------------------

! Local
      integer :: &
         klev, &        ! radiation layer count
         klevp          ! radiation interface count

      integer :: &
         kfrsnl      ! radiation interface index for fresnel layer

      ! following variables are defined for each layer; 0 refers to the top
      ! layer. In general we must distinguish directions above and below in
      ! the diffuse reflectivity and transmissivity, as layers are not assumed
      ! to be homogeneous (apart from the single layer Delta-Edd solutions);
      ! the direct is always from above.
      double precision, dimension (0:2*z_max) :: &
         rdir    , & ! layer reflectivity to direct radiation
         rdif_a  , & ! layer reflectivity to diffuse radiation from above
         rdif_b  , & ! layer reflectivity to diffuse radiation from below
         tdir    , & ! layer transmission to direct radiation (solar beam + diffuse)
         tdif_a  , & ! layer transmission to diffuse radiation from above
         tdif_b  , & ! layer transmission to diffuse radiation from below
         trnlay      ! solar beam transm for layer (direct beam only)

      integer :: &
         kk, testvar           ! level index

      double precision, parameter :: &
         trmin = 0.001d0   ! minimum total transmission allowed
      ! total transmission is that due to the direct beam; i.e. it includes
      ! both the directly transmitted solar beam and the diffuse downwards
      ! transmitted radiation resulting from scattering out of the direct beam
      double precision :: &
         tautot   , & ! layer optical depth
         wtot     , & ! layer single scattering albedo
         gtot     , & ! layer asymmetry parameter
         ftot     , & ! layer forward scattering fraction
         ts       , & ! layer scaled extinction optical depth
         ws       , & ! layer scaled single scattering albedo
         gs       , & ! layer scaled asymmetry parameter
         rintfc   , & ! reflection (multiple) at an interface
         refkp1   , & ! interface multiple scattering for k+1
         refkm1   , & ! interface multiple scattering for k-1
         tdrrdir  , & ! direct tran times layer direct ref
         tdndif       ! total down diffuse = tot tran - direct tran

      ! perpendicular and parallel relative to plane of incidence and scattering
      double precision :: &
         R1       , & ! perpendicular polarization reflection amplitude
         R2       , & ! parallel polarization reflection amplitude
         T1       , & ! perpendicular polarization transmission amplitude
         T2       , & ! parallel polarization transmission amplitude
         Rf_dir_a , & ! fresnel reflection to direct radiation
         Tf_dir_a , & ! fresnel transmission to direct radiation
         Rf_dif_a , & ! fresnel reflection to diff radiation from above
         Rf_dif_b , & ! fresnel reflection to diff radiation from below
         Tf_dif_a , & ! fresnel transmission to diff radiation from above
         Tf_dif_b     ! fresnel transmission to diff radiation from below

      ! refractive index for sea ice, water; pre-computed, band-independent,
      ! diffuse fresnel reflectivities
      double precision, parameter :: &
         refindx = 1.310d0  , & ! refractive index of sea ice (used for water also)
         cp063   = 0.063d0  , & ! diffuse fresnel reflectivity from above
         cp455   = 0.455d0  , & ! diffuse fresnel reflectivity from below
         p75     = 0.75d0   , &
         p5      = 0.5d0    , &
         p01     = 0.01d0   , &
         c1      = 1.0d0    , &
         c2      = 2.0d0    , &
         c3      = 3.0d0    , &
         c4      = 4.0d0

      double precision :: &
         mu0      , & ! cosine solar zenith angle incident
         mu0n         ! cosine solar zenith angle in medium

      double precision :: &
         alpha    , & ! term in direct reflectivity and transmissivity
         gamma    , & ! term in direct reflectivity and transmissivity
         el       , & ! term in alpha,gamma,n,u
         taus     , & ! scaled extinction optical depth
         omgs     , & ! scaled single particle scattering albedo
         asys     , & ! scaled asymmetry parameter
         u        , & ! term in diffuse reflectivity and transmissivity
         n        , & ! term in diffuse reflectivity and transmissivity
         lm       , & ! temporary for el
         muu       , & ! cosine solar zenith for either snow or water
         ne           ! temporary for n

      double precision :: &
         w        , & ! dummy argument for statement function
         uu       , & ! dummy argument for statement function
         gg       , & ! dummy argument for statement function
         e        , & ! dummy argument for statement function
         f        , & ! dummy argument for statement function
         t        , & ! dummy argument for statement function
         et           ! dummy argument for statement function

      double precision :: &
         alp      , & ! temporary for alpha
         gam      , & ! temporary for gamma
         ue       , & ! temporary for u
         arg      , & ! exponential argument
         extins   , & ! extinction
         amg      , & ! alp - gam
         apg          ! alp + gam

      integer, parameter :: &
         ngmax = 8    ! number of gaussian angles in hemisphere

      double precision, dimension (ngmax) :: &
         gauspt   , & ! gaussian angles (radians)
         gauswt       ! gaussian weights

      data gauspt/ &
                 .9894009d0,  .9445750d0, &
                 .8656312d0,  .7554044d0, &
                 .6178762d0,  .4580168d0, &
                 .2816036d0,  .0950125d0  /
      data gauswt/ &
                 .0271525d0,  .0622535d0, &
                 .0951585d0,  .1246290d0, &
                 .1495960d0,  .1691565d0, &
                 .1826034d0,  .1894506d0  /

      integer :: &
         ng           ! gaussian integration index

      double precision :: &
         gwt      , & ! gaussian weight
         swt      , & ! sum of weights
         trn      , & ! layer transmission
         rdr      , & ! rdir for gaussian integration
         tdr      , & ! tdir for gaussian integration
         smr      , & ! accumulator for rdif gaussian integration
         smt          ! accumulator for tdif gaussian integration

      ! Delta-Eddington solution expressions
      alpha(w,uu,gg,e) = p75*w*uu*((c1 + gg*(c1-w))/(c1 - e*e*uu*uu))
      gamma(w,uu,gg,e) = p5*w*((c1 + c3*gg*(c1-w)*uu*uu) &
                        / (c1-e*e*uu*uu))
      n(uu,et)         = ((uu+c1)*(uu+c1)/et ) - ((uu-c1)*(uu-c1)*et)
      u(w,gg,e)        = 1.5*(c1 - w*gg)/e
      el(w,gg)         = sqrt(c3*(c1-w)*(c1 - w*gg))
      taus(w,f,t)      = (c1 - w*f)*t
      omgs(w,f)        = (c1 - f)*w/(c1 - w*f)
      asys(gg,f)       = (gg - f)/(c1 - f)

!-----------------------------------------------------------------------

      ! these numbers are counts from 0, not the total number of layers...
      klev = nslyr + nilyr - 1  ! radiation layer count (0-klev, with klev+1 = total layers)
      klevp = klev + 1      ! radiation interface count(0-klevp, with klevp+1 = total interfaces)

      ! initialize all output to 0
      do kk = 0, klevp
         trndir(kk) = 0.d0
         trntdr(kk) = 0.d0
         trndif(kk) = 0.d0
         rupdir(kk) = 0.d0
         rupdif(kk) = 0.d0
         rdndif(kk) = 0.d0
      enddo

      ! initialize all layer apparent optical properties to 0
      do kk = 0, klev
         rdir(kk)   = 0.d0
         rdif_a(kk) = 0.d0
         rdif_b(kk) = 0.d0
         tdir(kk)   = 0.d0
         tdif_a(kk) = 0.d0
         tdif_b(kk) = 0.d0
         trnlay(kk) = 0.d0
      enddo

      ! initialize top interface of top layer
        trndir(0) =   c1
        trntdr(0) =   c1
        trndif(0) =   c1
        rdndif(0) =   0.d0

        ! compute level of fresnel refraction
       if( srftyp < 2 ) then
         ! if snow over sea ice or bare sea ice, fresnel level is
         ! at base of sea ice SSL (and top of the sea ice DL); the
         ! snow SSL counts for one, then the number of snow layers,
         ! then the sea ice SSL which also counts for one:
         kfrsnl = nslyr + 2
       else
         ! if ponded sea ice, fresnel level is the top of the pond
         kfrsnl = 0
       endif


      ! proceed down one layer at a time; if the total transmission to
      ! the interface just above a given layer is less than trmin, then no
      ! Delta-Eddington computation for that layer is done.

        ! begin main level loop
        do kk=0,klev

        ! initialize current layer properties to zero; only if total
        ! transmission to the top interface of the current layer exceeds the
        ! minimum, will these values be computed below:
        if ( kk .gt. 0 ) then
            ! Calculate the solar beam transmission, total transmission, and
            ! reflectivity for diffuse radiation from below at interface k,
            ! the top of the current layer k:
            !
            !              layers       interface
            !
            !       ---------------------  kk-1
            !                kk-1
            !       ---------------------  kk
            !                 kk
            !       ---------------------

              trndir(kk) = trndir(kk-1)*trnlay(kk-1)
              refkm1        = c1/(c1 - rdndif(kk-1)*rdif_a(kk-1))
              tdrrdir       = trndir(kk-1)*rdir(kk-1)
              tdndif        = trntdr(kk-1) - trndir(kk-1)
              trntdr(kk) = trndir(kk-1)*tdir(kk-1) + &
               (tdndif + tdrrdir*rdndif(kk-1))*refkm1*tdif_a(kk-1)
              rdndif(kk) = rdif_b(kk-1) + &
                (tdif_b(kk-1)*rdndif(kk-1)*refkm1*tdif_a(kk-1))
              trndif(kk) = trndif(kk-1)*refkm1*tdif_a(kk-1)
        endif       ! k > 0

        ! compute next layer Delta-eddington solution only if total transmission
        ! of radiation to the interface just above the layer exceeds trmin.
          if (trntdr(kk) > trmin ) then

        ! calculation over layers with penetrating radiation

           tautot  = tau(kk)
           wtot    = w0(kk)
           gtot    = g(kk)
           ftot    = gtot*gtot

           ts   = taus(wtot,ftot,tautot)
           ws   = omgs(wtot,ftot)
           gs   = asys(gtot,ftot)
           if (c3*(c1-ws)*(c1 - ws*gs) .lt. 0.) then
               testvar = -1
           endif

           lm   = el(ws,gs)
           ue   = u(ws,gs,lm)

           ! mu0 is cosine solar zenith angle above the fresnel level; make
           ! sure mu0 is large enough for stable and meaningful radiation
           ! solution: .01 is like sun just touching horizon with its lower edge

           mu0  = max(coszen,p01)

           ! mu0n is cosine solar zenith angle used to compute the layer
           ! Delta-Eddington solution; it is initially computed to be the
           ! value below the fresnel level, i.e. the cosine solar zenith
           ! angle below the fresnel level for the refracted solar beam:

           mu0n = sqrt(c1-((c1-mu0*mu0)/(refindx*refindx)))

           ! if level k is above fresnel level and the cell is non-pond, use the
           ! non-refracted beam instead

           if( srftyp .lt. 2 .and. kk .lt. kfrsnl ) mu0n = mu0

           !extins = max(exp_min, exp(-lm*ts/dx_exp))
           extins = exp(-lm*ts)
           !extins = explt(-lm*ts,use_expt,exp_bins,exp_mul,expt)
           ne = n(ue,extins)

           ! first calculation of rdif, tdif using Delta-Eddington formulas

           rdif_a(kk) = (ue+c1)*(ue-c1)*(c1/extins - extins)/ne
           tdif_a(kk) = c4*ue/ne

           ! evaluate rdir,tdir for direct beam
           !trnlay(k,ij) = max(exp_min, exp(-ts/(mu0n*dx_exp)))
           trnlay(kk) = exp(-ts/mu0n)
           !trnlay(kk) = explt(-ts/mu0n,use_expt,exp_bins,exp_mul,expt)

           alp = alpha(ws,mu0n,gs,lm)
           gam = gamma(ws,mu0n,gs,lm)
           apg = alp + gam
           amg = alp - gam
           rdir(kk) = amg*(tdif_a(kk)*trnlay(kk) - c1) + &
                       apg*rdif_a(kk)
           tdir(kk) = apg*tdif_a(kk) + &
                       (amg*rdif_a(kk) - (apg-c1))*trnlay(kk)

           ! recalculate rdif,tdif using direct angular integration over rdir,tdir,
           ! since Delta-Eddington rdif formula is not well-behaved (it is usually
           ! biased low and can even be negative); use ngmax angles and gaussian
           ! integration for most accuracy:
           swt = 0.
           smr = 0.
           smt = 0.
           do ng=1,ngmax
             muu  = gauspt(ng)
             gwt = gauswt(ng)
             swt = swt + muu*gwt
             !trn = max(exp_min, exp(-ts/(muu*dx_exp)))
             trn = exp(-ts/muu)
             !trn = explt(-ts/muu,use_expt,exp_bins,exp_mul,expt)
             alp = alpha(ws,muu,gs,lm)
             gam = gamma(ws,muu,gs,lm)
             apg = alp + gam
             amg = alp - gam
             rdr = amg*(tdif_a(kk)*trn-c1) + &
                   apg*rdif_a(kk)
             tdr = apg*tdif_a(kk) + &
                   (amg*rdif_a(kk)-(apg-c1))*trn
             smr = smr + muu*rdr*gwt
             smt = smt + muu*tdr*gwt
           enddo      ! ng
           rdif_a(kk) = smr/swt
           tdif_a(kk) = smt/swt

           ! homogeneous layer
           rdif_b(kk) = rdif_a(kk)
           tdif_b(kk) = tdif_a(kk)

           ! add fresnel layer to top of desired layer if either
           ! air or snow overlies ice; we ignore refraction in ice
           ! if a melt pond overlies it:

           if( kk == kfrsnl ) then
             ! compute fresnel reflection and transmission amplitudes
             ! for two polarizations: 1=perpendicular and 2=parallel to
             ! the plane containing incident, reflected and refracted rays.
             R1 = (mu0 - refindx*mu0n) / &
                  (mu0 + refindx*mu0n)
             R2 = (refindx*mu0 - mu0n) / &
                  (refindx*mu0 + mu0n)
             T1 = c2*mu0 / &
                  (mu0 + refindx*mu0n)
             T2 = c2*mu0 / &
                  (refindx*mu0 + mu0n)

             ! unpolarized light for direct beam
             Rf_dir_a = p5 * (R1*R1 + R2*R2)
             Tf_dir_a = p5 * (T1*T1 + T2*T2)*refindx*mu0n/mu0

             ! precalculated diffuse reflectivities and transmissivities
             ! for incident radiation above and below fresnel layer, using
             ! the direct albedos and accounting for complete internal
             ! reflection from below; precalculated because high order
             ! number of gaussian points (~256) is required for convergence:

             ! above
             Rf_dif_a = cp063
             Tf_dif_a = c1 - Rf_dif_a
             ! below
             Rf_dif_b = cp455
             Tf_dif_b = c1 - Rf_dif_b

             ! the kk = kfrsnl layer properties are updated to combined
             ! the fresnel (refractive) layer, always taken to be above
             ! the present layer kk (i.e. be the top interface):

             rintfc   = c1 / (c1-Rf_dif_b*rdif_a(kfrsnl))
             tdir(kfrsnl)   = Tf_dir_a*tdir(kfrsnl) + &
                                  Tf_dir_a*rdir(kfrsnl) * &
                                  Rf_dif_b*rintfc*tdif_a(kfrsnl)
             rdir(kfrsnl)   = Rf_dir_a + &
                                  Tf_dir_a*rdir(kfrsnl) * &
                                  rintfc*Tf_dif_b
             rdif_a(kfrsnl) = Rf_dif_a + &
                                  Tf_dif_a*rdif_a(kfrsnl) * &
                                  rintfc*Tf_dif_b
             rdif_b(kfrsnl) = rdif_b(kfrsnl) + &
                                  tdif_b(kfrsnl)*Rf_dif_b * &
                                  rintfc*tdif_a(kfrsnl)
             tdif_a(kfrsnl) = Tf_dif_a*rintfc*tdif_a(kfrsnl)
             tdif_b(kfrsnl) = tdif_b(kfrsnl)*rintfc*Tf_dif_b

             ! update trnlay to include fresnel transmission
             trnlay(kfrsnl) = Tf_dir_a*trnlay(kfrsnl)
           endif      ! kk = kfrsnl

          endif ! trntdr(kk,ij) > trmin

        enddo       ! kk   end main level loop

      ! compute total direct beam transmission, total transmission, and
      ! reflectivity for diffuse radiation (from below) for all layers
      ! above the underlying ocean; note that we ignore refraction between
      ! sea ice and underlying ocean:
      !
      !       For kk = klevp
      !
      !              layers       interface
      !
      !       ---------------------  kk-1
      !                kk-1
      !       ---------------------  kk
      !       \\\\\\\ ocean \\\\\\\

        kk = klevp
        trndir(kk) = trndir(kk-1)*trnlay(kk-1)
        refkm1        = c1/(c1 - rdndif(kk-1)*rdif_a(kk-1))
        tdrrdir       = trndir(kk-1)*rdir(kk-1)
        tdndif        = trntdr(kk-1) - trndir(kk-1)
        trntdr(kk) = trndir(kk-1)*tdir(kk-1) + &
          (tdndif + tdrrdir*rdndif(kk-1))*refkm1*tdif_a(kk-1)
        rdndif(kk) = rdif_b(kk-1) + &
          (tdif_b(kk-1)*rdndif(kk-1)*refkm1*tdif_a(kk-1))
        trndif(kk) = trndif(kk-1)*refkm1*tdif_a(kk-1)

      ! compute reflectivity to direct and diffuse radiation for layers
      ! below by adding succesive layers starting from the underlying
      ! ocean and working upwards:
      !
      !              layers       interface
      !
      !       ---------------------  kk
      !                 kk
      !       ---------------------  kk+1
      !                kk+1
      !       ---------------------

        rupdir(klevp) = albodr
        rupdif(klevp) = albodf

      do kk=klev,0,-1
          ! interface scattering
          refkp1        = c1/( c1 - rdif_b(kk)*rupdif(kk+1))
          ! dir from top layer plus exp tran ref from lower layer, interface
          ! scattered and tran thru top layer from below, plus diff tran ref
          ! from lower layer with interface scattering tran thru top from below
          rupdir(kk) = rdir(kk) + &
                        ( trnlay(kk)*rupdir(kk+1) + &
                         (tdir(kk)-trnlay(kk))*rupdif(kk+1) ) * &
                          refkp1*tdif_b(kk)
          ! dif from top layer from above, plus dif tran upwards reflected and
          ! interface scattered which tran top from below
          rupdif(kk) = rdif_a(kk) + &
                          tdif_a(kk)*rupdif(kk+1)* &
                          refkp1*tdif_b(kk)
      enddo       ! kk

      end subroutine sia2_env_edd_solution

