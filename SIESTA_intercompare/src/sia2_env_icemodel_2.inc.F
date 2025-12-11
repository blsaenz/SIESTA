! ==============================================================================
! START GROW ICE - ICE MODEL 2
! --------------------------------------------------------------------
      ! set up growth layer recorders to record ice/snow growth for time step
      ! and intialize workins vars and arrays
      ! ----------------------------------------------------------------

      ! vector initializations
      t_prev = 0.d0 ! vector assignment -- used by desal rountine
      heat_last = 0.d0 ! vector assignment -- used by desal rountine
      bs_last = 1.d0 ! vector assignment -- used by desal rountine
      bd_last = 1.d0 ! vector assignment -- used by desal rountine
      bv_last = 0.d0 ! vector assignment -- used by desal rountine
      bs_diff = 0.d0 ! vector assignment -- used by desal rountine
      ice(sc,ic,mi)%prod = 0.d0 ! vector assignment
      p_tot_int = 0.d0
      p_tot_bot = 0.d0
      heat_diff = d0_
      t_next = 0.d0
      t_last = 0.d0
      DC = 0.d0
      DL = 0.d0
      DU = 0.d0
      ice(sc,ic,mi)%dsdt = d0_
      ice(sc,ic,mi)%fbv = d0_
      !if (desal .lt. 2) then
      !      t_prev = ice(sc,ic,mi)%t
      !      bs_last = ice(sc,ic,mi)%bs
    !      bd_last = ice(sc,ic,mi)%bd
    !      bv_last = ice(sc,ic,mi)%bv
     ! endif

      !print *,'af: ',sc,ic,ice(sc,ic,mi)%af
      !print *,'snow_dist: ',sc,ic,ice(sc,ic,mi)%snow_dist
      !print *,'snow_rd: ',sc,ic,ice(sc,ic,mi)%snow_rd
      !print *,'snow_h: ',sc,ic,ice(sc,ic,mi)%sh_prev

      ! scalar initializations
      c_gl = 0.d0   ! growth layer congelation height to 0
      c_gl_tot = 0.d0
      c_gl_sal = 0.d0   ! growth layer salinity to 0
      sk_gl = 0.d0   ! growth layer skeletal height to 0
      s_gl = 0.d0   ! growth layer surface (only for sublimation - should never be positive)
      s_gl_tot = 0.d0
      snow_gl = 0.d0
      snow_melt = 0.d0
      melt_drain = 0.d0
      crunch_tot = 0.d0
      fl_gl = 0.d0
      snow_gl = 0.d0
      flooded = 0.d0
      dhdt_cm_s = 0.d0   ! required in brine flux, and algal migration
      b_flux_max_tot = 0.d0
      bfrf = 1 ! brine flux reduction factor - higher number will cause finer resolution brine flux
      flood_1x = .true.
      did_flood = .false.       ! init var that tells snow routine that flooding happened
      did_drain = .false.       ! init var that tells snow routine that flooding happened

      ! initial bookeeping operations
      ! --------------------------------------------------------------------
      ! find end of bottom internal layer, top & bottom skeletal layers
      sk_z = ice(sc,ic,mi)%z
      int_z = sk_z - z_sk
      sk_1 = int_z + 1

      ! Calc solar zenith angle, which varies sinusoidally over the day at a given
      ! latitude at a given date during the year
      zenith_angle = acos(sin(solar_dec)*sin(m(mi)%lat/360*6.2832) &
      + cos(solar_dec)*cos(m(mi)%lat/360*6.2832)*cos(rad_time) )   ! radians

      ! If zenith angle is below zero, then make it zero so we don't calculate
      ! negative light values
      if (zenith_angle .lt. 0.) then
          zenith_angle = 0.d0
      endif
      coszen = cos(zenith_angle)

      ! find spectral reflection at air/snow or ice interface
      if (zenith_angle .le. 0.8) then
          rsun0=0.018736*10**(0.15316*zenith_angle)
      else
          rsun0=0.00036779*10**(2.1505*zenith_angle)
      endif

      if (airtemp1c .gt. c0) then
          rsky0 = 0.d0
      else
          rsky0 = 0.05d0
      endif

      ! Set up for heat transfer by determining new surface flux equation
      ! --------------------------------------------------------------------
      ! water vapor saturation pressure
      p_h2o_sat = 6.1078d0*10**((7.5d0*f(mi)%at-2048.625d0)/ &
      (f(mi)%at-35.85d0))  ! mb (1 mb = 100 Pa = 1 hPa)

       if (f(mi)%rhum .eq. c9999) then  ! calculate rhum from specific humidity
           ! find air density

           !                     total air pressure
           ! air density = -------------------------------- , sh = specific humidity
           !               airtemp(R_air-sh*R_air+sh*R_h2o)
           row_air = &
               f(mi)%p/(f(mi)%at*(R_air - f(mi)%h*R_air + f(mi)%h*R_h2o))
           rhum = row_air*f(mi)%h*R_h2o*f(mi)%at &
               /(p_h2o_sat*100.)
           p_h2o = p_h2o_sat*rhum
           shum = f(mi)%h
       else
           rhum = f(mi)%rhum
           p_h2o = p_h2o_sat*rhum
           row_air = (f(mi)%p*R_h2o + 100.*rhum*p_h2o_sat*(R_air - R_h2o)) / &
               (f(mi)%at*R_air*R_h2o)
           shum = 0.622d0*p_h2o/(f(mi)%p*0.01d0 - 1.622d0*p_h2o)
       endif

!      rhum = f(mi)%rhum
!      shum = f(mi)%h
!      row_air = (f(mi)%p*R_h2o + 100.*f(mi)%rhum*p_h2o_sat*(R_air - R_h2o)) / &
!              (f(mi)%at*R_air*R_h2o)
      !print *,'rhum, shum, at: ',rhum,shum,f(mi)%at-kelvin0


      !potenial T = T(kelvin)*(P0/P)^(R/cp)
      potT = f(mi)%at*(100000.d0/f(mi)%p)**(R_air/cp)

      ! find longwave emissivity basced on atmospheric cloud cover;  (König-Langlo & Augstein 1994)
      epsa = 0.765 + 0.22*(f(mi)%fc*0.01)**3 ! longwave emissivity of air

      ! find one static part of turbulent heat flux equation
      if (f(mi)%ws .eq. c9999) then
          v10 = sqrt(f(mi)%u10**2+f(mi)%v10**2)
      else
          v10 = f(mi)%ws
      endif
      Fe_sub1 = 0.622*(Lv*1000.)*row_air*Ce*v10/(f(mi)%p*c_01)  ! J/kg/K * kg/m^3 * (dimensionless) * m/s * 1/mbar = J/m^2/mbar/K

      ! find second static part of turbulent heat flux equation
      Fe_sub2 = Fe_sub1*(rhum*(fe_A*f(mi)%at**4 &
          + fe_B*f(mi)%at**3 &
          + fe_C*f(mi)%at**2 &
          + fe_D*f(mi)%at + fe_E) - fe_E)

      ! surface longwave in
      if (f(mi)%lwd .eq. c9999) then
          FL = epsa*steph_boltz*f(mi)%at**4
      else
          Fl = f(mi)%lwd
      endif

      ! static part of sensible latent heat flux (Fs)
      Fs_sub = row_air*cp*Ch*v10*f(mi)%at

      ! assign c_gl_bv to that of last layer for 1st sub-dt step - it will be
      ! always one step behind due to the interdependence for keff and dhdt...
      c_gl_bv = ice(sc,ic,mi)%bv(int_z)*0.001

      ! setup cloud correction factor
      if (par_cf .eq. 1) then
          par_cf_f = (1. - 0.53*sqrt(f(mi)%fc*c_01))
      else
          par_cf_f = 1.
      endif

      ! get current snow depth, and assign to snowh1m
      snowh1m = ice(sc,ic,mi)%sh_prev
      if (h_snowpatch .eq. 0.) then
          f_snow = 1.
      else
          f_snow = snowh1m/(snowh1m + h_snowpatch)
      endif

      ! initialize shortwave storage vectors
      ice(sc,ic,mi)%PAR_bot = 0.  ! vector assignment

      lda_f = lda_d  ! init light/snow distribution divider

      ! find Ed0 is it has not been found yet
#ifdef GPU
      if (use_gpu .eq. 1) then
          ! recall ed_w arrays
          Ed_W_mean = ice(sc,ic,mi)%ed_w ! vector assignment
          Ed_W_snow = ice(sc,ic,mi)%snow%ed_w  ! vector assignment
      else
#endif
          Ed_W_mean = 0. ! vector assignment
          Ed_W_snow = 0.  ! vector assignment
          ice(sc,ic,mi)%Ed0 = 0.d0
          do ii = 1,wavl
              ! Find Total shortwave irradiance in Watts/m^2
              ice(sc,ic,mi)%Ed0 = ice(sc,ic,mi)%Ed0 + &
              (mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii) + mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii)) * &
                  (7.e-7*lambda(ii)**2 - 0.0012*lambda(ii) + 0.6416)*par_cf_f*10.  ! Watts/m^2
          enddo

          ! rectify spectral watts with shortwave down?
          print *,'Ed0-1',ice(sc,ic,mi)%Ed0
          print *,'swd-1',f(mi)%swd
          if (boundary_file .eq. 1 .and. f(mi)%swd .ne. c9999) then

            ! case where swd shows sunlight, but gc says no
            if (ice(sc,ic,mi)%Ed0 .eq. 0.d0 .and. f(mi)%swd .gt. 0.d0) then
                ! find last valid spectral profile to use with swd
                mp_f%Edir(f_index)%Ed(i_mp,j_mp,:) = Edir_last_valid
                mp_f%Edif(f_index)%Ed(i_mp,j_mp,:) = Edif_last_valid
                do ii = 1,wavl
                    ! Find Total shortwave irradiance in Watts/m^2
                    ice(sc,ic,mi)%Ed0 = ice(sc,ic,mi)%Ed0 + &
                    (mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii) + mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii)) * &
                        (7.e-7*lambda(ii)**2 - 0.0012*lambda(ii) + 0.6416)*par_cf_f*10.  ! Watts/m^2
                enddo
            endif

            if (ice(sc,ic,mi)%Ed0 > 0.d0) then
              tmp1 = f(mi)%swd / (ice(sc,ic,mi)%Ed0 * par_to_swd) ! factor to modify gc_par by
              do ii = 1,wavl
                mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii) = mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii)*tmp1
                mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii) = mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii)*tmp1
              enddo
              ice(sc,ic,mi)%Ed0 = f(mi)%swd / par_to_swd
              !f(mi)%swd = c9999

              ! only one copy over model domain - only works with boundary_file
              Edir_last_valid = mp_f%Edir(f_index)%Ed(i_mp,j_mp,:)
              Edif_last_valid = mp_f%Edif(f_index)%Ed(i_mp,j_mp,:)


            endif
          endif
          print *,'Ed0-2',ice(sc,ic,mi)%Ed0
          print *,'swd-2',f(mi)%swd
#ifdef GPU
      endif
#endif

      ! calculate light up front - once per time step
      ! ----------------------------------------------------------------
      ! sub-dt timestep is too costly b/c of power operator

      !print *,'zenith: ',zenith_angle,'coszen:',coszen
      !print *,ice(sc,ic,mi)%Ed0

      do kk = 1,lda_n

          if (ice(sc,ic,mi)%Ed0 .gt. 0.) then


#ifdef GPU

              if (use_gpu .eq. 1) then

                  ! light is pre-calculated using GPU
                  PUR_dble(:,kk) = ice(sc,ic,mi)%pur(:,kk)

              else
#endif

!#include "sia2_env_light.inc.f90"
!#include "sia2_env_light_2-stream_bulk.inc.f90"
#include "sia2_env_edd.inc.f90"

#include "sia2_env_light_platelet.inc.f90"

#ifdef GPU
              endif
#endif

              !ed_w_snow = ed_w_snow*1.415
              !ed_w_mean = ed_w_mean*1.415
              !PUR_dble = PUR_dble*1.415

              do ii = 1,sk_z
                  ! convert to int and save to daily PUR array
                  ! ----------------------------------------------------
!                  tmp1 = sqrt(PUR_dble(ii,kk)*one_over_pur_c)   ! for byte size pur array
!                  pur(ii,pur_clock,kk,sc,ic,mi) = min(nint(tmp1),255)-127   ! for byte size pur array

                  pur(ii,pur_clock,kk,sc,ic,mi) = min(nint(PUR_dble(ii,kk)*one_over_pur_c_2),65535)-32767   ! for int(2) size pur array

                  ! calc light limitation for this layer (ii) and light distribution (kk)
                  ! ----------------------------------------------------
                  PUR_mean24 = 0.
                  do jj=1,dt_per_day
                      if (jj .eq. pur_clock) then
                          PUR_mean24 = PUR_mean24 + PUR_dble(ii,kk)
                      else
                          PUR_mean24 = PUR_mean24 + &
!                            pur_c*(dble(pur(ii,jj,kk,sc,ic,mi))+127.d0)**2  ! for byte size pur array
                            pur_c_2*dble(pur(ii,jj,kk,sc,ic,mi) + 32767)     ! for int(2) size pur array
                      endif
                  enddo
                  PUR_mean24 = PUR_mean24/dble(dt_per_day)
                  Ik1(ii,kk) = Ek_max/(1. + A*exp(-1.*B*PUR_mean24))
                  !tmp1 = explt(-1.*B*PUR_mean24,use_expt,exp_bins,exp_mul,expt)
                  !Ik1(ii,kk) = Ek_max/(1. + A*tmp1)

                  if (Ik1(ii,kk) .le. 0.) then
                      ll(ii,kk) = 0.
                  else
                      ll(ii,kk) = 1. - exp(-1.*PUR_dble(ii,kk)/Ik1(ii,kk))
                      !ll(ii,kk) = 1. - explt(-1.*PUR(ii,pur_clock,kk,sc,ic,mi)/Ik1(ii,kk), &
                      !    use_expt,exp_bins,exp_mul,expt)

                  endif
              enddo
              if (sk_z .lt. z_max) then
                  do ii=sk_z+1,z_max
                      ll(ii,kk) = -1.
                  enddo
              endif

          else

              ! no light - set light vars to zero
              PUR_dble(:,kk) = 0.
!              pur(:,pur_clock,kk,sc,ic,mi) = -127
              pur(:,pur_clock,kk,sc,ic,mi) = pur_0

              ll(:,kk) = 0.
!              do ii = 1,sk_z
!                  PUR_dble(ii,kk) = 0.
!                  pur(ii,pur_clock,kk,sc,ic,mi) = -127
!                 ll(ii,kk) = 0.
!              enddo
              if (use_pl .eq. 1) then
                  do ii = 1,pl_max
                      pl%pur(ii) = 0.
                  enddo
              endif
              ice(sc,ic,mi)%PAR_bot = 0.
              ice(sc,ic,mi)%PAR_bot_pl = 0.
          endif

      enddo

!      if (i_ease .eq. 281 .and. j_ease .eq. 306) then
!      if (mi .eq. 27) then
!          print *,'pur -- ',PUR_dble(:,1)
!          print *,'edw_i -- ',ed_w_mean
!          print *,'edw_s -- ',ed_w_snow
!          if (use_gpu .eq. 1) then
!              print *,'ed_nir_dir: ',ice(sc,ic,mi)%ed0_nir_dir
!              print *,'ed_nir_dif: ',ice(sc,ic,mi)%ed0_nir_dif
!          print *,'i-again -- ',ice(sc,ic,mi)%ed_w
!          print *,'i-again -- ',ice(sc,ic,mi)%snow%ed_w
!          else
!              print *,'ed_nir_dir: ',Ed0_nir_dir
!              print *,'ed_nir_dif: ',Ed0_nir_dif
 !         endif
!      endif

      !DEBUG - REMOVE ME
      ! save over par_bot with total absorption
!      ice(sc,ic,mi)%par_bot = 0.
!     do ii=1,sk_z
!         ice(sc,ic,mi)%par_bot = ice(sc,ic,mi)%par_bot + ed_w_mean(ii)
!         if (ii .le. z_snow) then
!             ice(sc,ic,mi)%par_bot = ice(sc,ic,mi)%par_bot + ed_w_snow(ii)
!         endif
!     enddo


      ! print *,"Pond Thickness: ",ice(sc,ic,mi)%pond%th

!#include "sia2_env_kevin_light.inc.f90"

      ! convert tracers to bulk volume concentrations,zero other stuff
      do ii=1,sk_z
          tmp1 = ice(sc,ic,mi)%bv(ii)*c_001
          ice(sc,ic,mi)%no3(ii) = ice(sc,ic,mi)%no3(ii)*tmp1
          ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii)*tmp1
          ice(sc,ic,mi)%po4(ii) = ice(sc,ic,mi)%po4(ii)*tmp1
          ice(sc,ic,mi)%sioh4(ii) = ice(sc,ic,mi)%sioh4(ii)*tmp1
          ice(sc,ic,mi)%smalg(ii) = ice(sc,ic,mi)%smalg(ii)*tmp1
          ice(sc,ic,mi)%poc(ii) = ice(sc,ic,mi)%poc(ii)*tmp1
      enddo

      ! record dsdt_net, fresh_net masses
        call sia2_ice_mass(ice(sc,ic,mi),fresh_net,dsdt_net)

      ! ----------------------------------------------------------------
      ! BEGIN SUB-DT PHYSICS TIMESTEPS
      ! ----------------------------------------------------------------

      iii = 0
      dtt = dtt_1
      dtt_s = dtt_s_1

           if (dt_step .eq. 125 .and. mi .eq. 29361 .and. ic .eq. 2) then
               testvar = -1
               !call exit (0)
           endif

        do while (iii .lt. dtt)

          ! check to make sure flooding didn't nullify entire snow category mid-step
          ! af compare to 0.d0 not reliable on windows!
          if ((ice(sc,ic,mi)%z > 0) .AND. (ice(sc,ic,mi)%af .gt. 0.d0)) then


          if (dt_step .eq. 211 .and. iii .eq. 0) then
          testvar = -1
          endif

          if (ice(1,1,mi)%snow%ts .lt. -100.) then
              testvar = -1
              testvar = sqrt(testvar)
          endif

          !print *, "dtt step: ",iii

          ! sub-dt loop initializations
          ! ------------------------------------------------------------

          ! vector inits
          heat_temp = d0_ ! zero mean heat storage vector
          dsdt = d0_
          dno3 = d0_
          dnh4 = d0_
          dpo4 = d0_
          dsioh4 = d0_
          dpoc = d0_
          dsmalg = d0_
          dcheat = d0_      ! req'd by desal routine, and later in heat flux
          b_flux = d0_      ! req'd by desal routines, and later in heat flux
          h_flux = d0_      ! req'd by desal routines, and later in heat flux
          s_new = ice(sc,ic,mi)%s ! req'd by desal routines, and later in heat flux
          d_new = ice(sc,ic,mi)%d  ! req'd by desal routines, and later in heat flux

          !scalar inits
          heat_stop = int_z  ! assume entire ice column is involved in heat-flux calc
          integr_valid = .true.   ! default is valid integration
!          print *,'Pond th: ', ice(sc,ic,mi)%pond%th
          th_pond_new = ice(sc,ic,mi)%pond%th
          s_pond_new = ice(sc,ic,mi)%pond%s
          h_pond_new = ice(sc,ic,mi)%pond%heat
          poc_pond_new = ice(sc,ic,mi)%pond%poc
          no3_pond_new = ice(sc,ic,mi)%pond%no3
          nh4_pond_new = ice(sc,ic,mi)%pond%nh4
          po4_pond_new = ice(sc,ic,mi)%pond%po4
          sioh4_pond_new = ice(sc,ic,mi)%pond%sioh4

          ! find fraction of ice that has snow cover, compute areal multitpliers
          if (h_snowpatch .eq. 0.) then
              f_snow = 1.
          else
              f_snow = snowh1m/(snowh1m + h_snowpatch)
          endif

          ! find snow/ice thermal conductivities and inter-layer thickesses
#include "sia2_env_ki.inc.f90"

          if (desal .eq. 7) then
#include "sia2_env_desalinate_7.inc.f90"
          endif

          ! save temp and heat changes - required for next desal cycle, and so that
          ! further updates below don't barf when layer numbers change
          t_prev = ice(sc,ic,mi)%t
          t_snow_prev = ice(sc,ic,mi)%snow%t
          heat_last = ice(sc,ic,mi)%heat
          bs_last = ice(sc,ic,mi)%bs
          bd_last = ice(sc,ic,mi)%bd
          bs_diff = ice(sc,ic,mi)%bs

          ! use mean Ed_W for use in heat conduction
          Ed_W_temp = Ed_W_mean

#include "sia2_env_flux_heat.inc.f90"

          if (z_snow .eq. 0) then
              T_melt = ice(sc,ic,mi)%s(1)*mu
              ts_snow = min(Ts_next,0.)
              !ts_snow = Ts_next
          else
              ! snow can't be above 0degC
              ts_snow = min(Ts_next,0.)
          endif
          Fc_snow = Fc_top

          ! record new snow temp/flux
          if (z_snow .gt. 0) then
              t_snow_temp(1:z_snow) = T_next(1:z_snow)
          endif

          ! add to mean ice heat
          do ii=1,int_z
              !T_mean = min(T_next(ice_1+ii-1),ice(sc,ic,mi)%s(ii)*mu)
              T_mean = T_next(ice_1+ii-1)
              S_mean = ice(sc,ic,mi)%s(ii)

#include "sia2_env_brine_calc.inc.f90"

              heat_temp(ii) = heat_mean
          enddo

          if (integr_valid .or. (dtt_s .eq. dtt_s_3)) then

              ! figure out timestep in days
              dtt_days = dtt_s/86400.

              ! RECORD NEW ICE TEMPERATURES, UPDATE ICE PHYSCS
              ! ---------------------------------------------------------------

              ! record changes for desal rountine, before changing salinity
              ! to determine temp/brine changes purely due to heat conduction
              ! (i.e. excluding convection as a souce of convection!)
              if (desal .ge. 2) then
                  do ii=1,int_z
                      !ice(sc,ic,mi)%t(ii) = min(T_next(ice_1+ii-1),ice(sc,ic,mi)%s(ii)*mu)
                      t_mean = T_next(ice_1+ii-1)
                      ! check to make sure we don't get too warm (4degC)
                      s_mean = ice(sc,ic,mi)%s(ii)
                      t_mean = min(t_mean,min_t)

#include "sia2_env_brine_calc.inc.f90"

                      ! do calcs used by desal routine in next sub-dt timestep
                      t_change(ii) = t_prev(ii) - t_mean
                      heat_last(ii) = heat_last(ii) - heat_mean
                      bs_last(ii) = bs_mean/bs_last(ii)
                      bd_last(ii) = bd_mean/bd_last(ii)
                      bs_diff(ii) = bs_diff(ii) - bs_mean
                  enddo
              endif


              ! update ice
              ! --------------------------------------------------------
              if (use_ponds .eq. 1) then
                  if (th_pond_new .gt. 0.d0) then
                      ice(sc,ic,mi)%pond%th = th_pond_new
                      ice(sc,ic,mi)%pond%perc = min(th_pond_new,ice(sc,ic,mi)%pond%perc)
                      ice(sc,ic,mi)%pond%s = s_pond_new
                      ice(sc,ic,mi)%pond%d = 1000000.+800.*ice(sc,ic,mi)%pond%s  ! find density of melt pond water
                      ice(sc,ic,mi)%pond%heat = h_pond_new
                      ice(sc,ic,mi)%pond%t = ice(sc,ic,mi)%pond%heat/ice(sc,ic,mi)%pond%d/cw
                      ice(sc,ic,mi)%pond%poc = poc_pond_new
                      ice(sc,ic,mi)%pond%no3 = no3_pond_new
                      ice(sc,ic,mi)%pond%nh4 = nh4_pond_new
                      ice(sc,ic,mi)%pond%po4 = po4_pond_new
                      ice(sc,ic,mi)%pond%sioh4 = sioh4_pond_new
                  else
                      ice(sc,ic,mi)%pond%th = 0.d0
                      ice(sc,ic,mi)%pond%perc = 0.d0
                      ice(sc,ic,mi)%pond%t = 0.d0
                      ice(sc,ic,mi)%pond%s = 0.d0
                      ice(sc,ic,mi)%pond%d = 0.d0
                      ice(sc,ic,mi)%pond%heat = 0.d0
                      ice(sc,ic,mi)%pond%smalg = 0.d0     ! brine-based algal concentration (mgC/m^3)
                      ice(sc,ic,mi)%pond%poc = 0.d0     ! particulate organic carbon (detritus) (mgC/m^3)
                      ice(sc,ic,mi)%pond%no3 = 0.d0        ! ice layer brine NO3 concentration (µMol)
                      ice(sc,ic,mi)%pond%nh4 = 0.d0      ! ice layer brine NH4 concentration (µMol)
                      ice(sc,ic,mi)%pond%po4 = 0.d0     ! ice layer brine PO4 concentration (µMol)
                      ice(sc,ic,mi)%pond%sioh4 = 0.d0    ! ice layer brine SiOH4 concentration (µMol)
                  endif
               endif



              do ii=1,int_z

                  ! update modified tracer concentrations - dtracer is in brine concentration
                  bv_mean=ice(sc,ic,mi)%bv(ii)*c_001
                  ice(sc,ic,mi)%smalg(ii) = ice(sc,ic,mi)%smalg(ii) + dsmalg(ii)*bv_mean
                  if (ice(sc,ic,mi)%smalg(ii) .lt. 1.e-9) ice(sc,ic,mi)%smalg(ii) = 0.
                  ice(sc,ic,mi)%poc(ii) = ice(sc,ic,mi)%poc(ii) + dpoc(ii)*bv_mean
                  if (ice(sc,ic,mi)%poc(ii) .lt. 1.e-9) ice(sc,ic,mi)%poc(ii) = 0.
                  ice(sc,ic,mi)%no3(ii) = ice(sc,ic,mi)%no3(ii) + dno3(ii)*bv_mean
                  if (ice(sc,ic,mi)%no3(ii) .lt. 1.e-9) ice(sc,ic,mi)%no3(ii) = 0.
                  ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii) + dnh4(ii)*bv_mean
                  if (ice(sc,ic,mi)%nh4(ii) .lt. 1.e-9) ice(sc,ic,mi)%nh4(ii) = 0.
                  ice(sc,ic,mi)%po4(ii) = ice(sc,ic,mi)%po4(ii) + dpo4(ii)*bv_mean
                  if (ice(sc,ic,mi)%po4(ii) .lt. 1.e-9) ice(sc,ic,mi)%po4(ii) = 0.
                  ice(sc,ic,mi)%sioh4(ii) = ice(sc,ic,mi)%sioh4(ii) + dsioh4(ii)*bv_mean
                  if (ice(sc,ic,mi)%sioh4(ii) .lt. 1.e-9) ice(sc,ic,mi)%sioh4(ii) = 0.

                  !ice(sc,ic,mi)%t(ii) = min(T_next(ice_1+ii-1),ice(sc,ic,mi)%s(ii)*mu)
!                  if (ii .lt. heat_stop) then
                      ice(sc,ic,mi)%t(ii) = T_next(ice_1+ii-1)
!                  endif

                  if (desal .gt. 1) then
                       ice(sc,ic,mi)%s(ii) = s_new(ii)
                      !ice(sc,ic,mi)%s(ii) = ice(sc,ic,mi)%s(ii) + dsdt(ii)
                  endif

                  ! total agua and salt accounting
                  ! convert bulk dsdt (psu or ppt) to kg/m^2: g/kg * fractional_area * m(thickness) * 1e3 kg/m^3 * 1kg/1000g
                  !if (ice(sc,ic,mi)%dsdt(ii) .lt. 0.d0) then
                      !m(mi)%dsdt_out = m(mi)%dsdt_out - ice(sc,ic,mi)%dsdt(ii)*ice(sc,ic,mi)%af*ice(sc,ic,mi)%th(ii) ! kg/m^2
                  !endif

                  ! check to make sure we don't get too warm (4degC)
                  ice(sc,ic,mi)%t(ii) = min(ice(sc,ic,mi)%t(ii),min_t)
                  !if (ice(sc,ic,mi)%id(ii) .le. ice(sc,ic,mi)%fbh) then
                  !      ice(sc,ic,mi)%t(ii) = min(ice(sc,ic,mi)%t(ii),ice(sc,ic,mi)%s(ii)*(-0.054))
                  !else
                  !      ice(sc,ic,mi)%t(ii) = min(ice(sc,ic,mi)%t(ii),f(mi)%t)
                  !endif

                  T_mean = ice(sc,ic,mi)%t(ii)
                  S_mean = ice(sc,ic,mi)%s(ii)

#include "sia2_env_brine_calc.inc.f90"

                  ! update brine params for next sub-dt growth step - which are needed? just d?
                  ice(sc,ic,mi)%d(ii) = d_mean
                  ice(sc,ic,mi)%bs(ii) = bs_mean
                  ice(sc,ic,mi)%bd(ii) = bd_mean
                  ice(sc,ic,mi)%bv(ii) = bv_mean
                  ice(sc,ic,mi)%heat(ii) = heat_mean
              enddo

              ! RECORD NEW SNOW TEMPERATURES, UPDATE SNOW PHYSICS
              ! -----------------------------------------------------------
              z_snow = ice(sc,ic,mi)%snow%z
              ! update surface temperatures
              ice(sc,ic,mi)%snow%ts = ts_snow
              if (z_snow .gt. 0) then
                  jj=1
                  do ii = z_snow,1,-1
                      ! relcalc snow paramaters
                      T_mean = min(t_snow_temp(ii),0.)

                      T_mean = T_mean + kelvin0
                      d_mean = ice(sc,ic,mi)%snow%d(jj)

#include "sia2_env_snow_calc.inc.f90"

                      ice(sc,ic,mi)%snow%t(jj) = T_mean - kelvin0
                      ice(sc,ic,mi)%snow%d(jj) = d_mean
                      ice(sc,ic,mi)%snow%heat(jj) = heat_mean
                      if (ice(sc,ic,mi)%snow%t(jj) .ge. -2.) then
                          ice(sc,ic,mi)%snow%melt(jj) = 1.
                      endif

                      jj=jj+1
                  enddo
!                      if (ice(sc,ic,mi)%snow%t(z_snow) .ge. -1.) then
!                          ice(sc,ic,mi)%snow%melt(1:z_snow) = 1.
!                      endif

              endif

              ! flux before grid change, since desal-induced flux is calcuated
              ! above is for the current number of layers, and grid calcs
              ! may change the number of layers
              if (ice(sc,ic,mi)%id(int_z) .ge. 0.2) then
                  dt_flux = dtt_s/dble(bfrf)
                  do jjj=1,bfrf

#include "sia2_env_flux_brine.inc.f90"

                  enddo
              endif

              ! record pre-growth nutrient concentrations available to algae
              ice(sc,ic,mi)%no3_mean = ice(sc,ic,mi)%no3
              ice(sc,ic,mi)%nh4_mean = ice(sc,ic,mi)%nh4
              ice(sc,ic,mi)%po4_mean = ice(sc,ic,mi)%po4
              ice(sc,ic,mi)%sioh4_mean = ice(sc,ic,mi)%sioh4



#include "sia2_env_grid_calcs.inc.f90"

              ! protect against nullification of current category during flooding
              ! af compare to 0.d0 not reliable on windows!
              if ((ice(sc,ic,mi)%z > 0) .AND. (ice(sc,ic,mi)%af .gt. 0.d0)) then

#include "sia2_env_grid_remap.inc.f90"

#include "sia2_env_grow_algae_dtt.inc.f90"

                  ! calculate mean nutrient delivery available to algae
                  ice(sc,ic,mi)%no3_mean = c_5*(ice(sc,ic,mi)%no3_mean + ice(sc,ic,mi)%no3)/dtt_s
                  ice(sc,ic,mi)%nh4_mean = c_5*(ice(sc,ic,mi)%nh4_mean + ice(sc,ic,mi)%nh4)/dtt_s
                  ice(sc,ic,mi)%po4_mean = c_5*(ice(sc,ic,mi)%po4_mean + ice(sc,ic,mi)%po4)/dtt_s
                  ice(sc,ic,mi)%sioh4_mean = c_5*(ice(sc,ic,mi)%sioh4_mean + ice(sc,ic,mi)%sioh4)/dtt_s

                  ! increment sub-dt time step
                  iii = iii+1

              else  ! end of mid-step nullifaction-due-to-flooding check

                  iii = dtt

              endif

        else ! end of adaptive timestep check --> if (integr_valid) ...

              ! reduce timestep
              if (dtt_s .eq. dtt_s_1) then
                    dtt = dtt_2
                    dtt_s = dtt_s_2
                    iii = iii*dtt_2/dtt_1
                    sh_increment = sh_increment*dtt_s_2/dtt_s_1
              elseif (dtt_s .eq. dtt_s_2) then
                    dtt = dtt_3
                    dtt_s = dtt_s_3
                    iii = iii*dtt_3/dtt_2
                    sh_increment = sh_increment*dtt_s_3/dtt_s_2
              endif

              ! reset inter-dtt timestep variables
              t_prev = 0. ! vector assignment -- used by desal rountine
              heat_last = 0. ! vector assignment -- used by desal rountine
              bs_last = 1. ! vector assignment -- used by desal rountine
              bd_last = 1. ! vector assignment -- used by desal rountine
              bv_last = 0. ! vector assignment -- used by desal rountine
              bs_diff = 0. ! vector assignment -- used by desal rountine

          endif

          else

            iii = dtt

          endif  ! end of whole-step nullifaction-due-to-flooding check

    ! ----------------------------------------------------------------
      ! END SUB-DT PHYSICS TIMESTEPS
      ! ----------------------------------------------------------------
      enddo

      ! record timestep used
      if (dtt .eq. dtt_1) then
          n_dtt_1 = n_dtt_1 + 1
      elseif (dtt .eq. dtt_2) then
          n_dtt_2 = n_dtt_2 + 1
      elseif (dtt .eq. dtt_3) then
          n_dtt_3 = n_dtt_3 + 1
      endif

            ! find sub-step ending fresh and salt mass
            call sia2_ice_mass(ice(sc,ic,mi),fresh_tmp,dsdt_tmp)

            ! record fresh water flux (kg)
            fresh_net = (fresh_tmp-fresh_net)*ice(sc,ic,mi)%af*cell_area*1.e6
            m(mi)%h2o_flux = m(mi)%h2o_flux - fresh_net

            ! record salt flux (kg)
            dsdt_net = (dsdt_tmp-dsdt_net)*ice(sc,ic,mi)%af*cell_area*1.e6
            m(mi)%salt_flux = m(mi)%salt_flux - dsdt_net

      tmp1 = ice(sc,ic,mi)%id(ice(sc,ic,mi)%z)
      tmp2 = tmp1-bot_th

      ! record productivity weights for habitat characterization
      tmp3 = p_tot_int + p_tot_bot ! total production - m^-2
      if (tmp1 .gt. 0.) then
          tmp4 = tmp3/tmp1 ! total productivity - m^-3
            else
                tmp4 = 0.d0
            endif
!            tmp3 = p_tot_int
!        tmp4 = p_tot_bot

      m(mi)%p_wgt_int = m(mi)%p_wgt_int + tmp3
      m(mi)%p_wgt_af_int = m(mi)%p_wgt_af_int + tmp3*ice(sc,ic,mi)%af
      m(mi)%p_wgt_bot = m(mi)%p_wgt_bot + tmp4
      m(mi)%p_wgt_af_bot = m(mi)%p_wgt_af_bot + tmp4*ice(sc,ic,mi)%af

      ! record productivity-weighted ice and snow thickness and ice surf temp
      m(mi)%pwid_sum_int = m(mi)%pwid_sum_int + tmp3*tmp1
      m(mi)%pwid_sum_af_int = m(mi)%pwid_sum_af_int + tmp3*tmp1*ice(sc,ic,mi)%af
      m(mi)%pwid_sum_bot = m(mi)%pwid_sum_bot + tmp4*tmp1
      m(mi)%pwid_sum_af_bot = m(mi)%pwid_sum_af_bot + tmp4*tmp1*ice(sc,ic,mi)%af
      tmp1 = ice(sc,ic,mi)%snow%depth
      m(mi)%pwsd_sum_int = m(mi)%pwsd_sum_int + tmp3*tmp1
      m(mi)%pwsd_sum_af_int = m(mi)%pwsd_sum_af_int + tmp3*tmp1*ice(sc,ic,mi)%af
      m(mi)%pwsd_sum_bot = m(mi)%pwsd_sum_bot + tmp4*tmp1
      m(mi)%pwsd_sum_af_bot = m(mi)%pwsd_sum_af_bot + tmp4*tmp1*ice(sc,ic,mi)%af
      tmp1 = ice(sc,ic,mi)%t(1)
      m(mi)%pwt_sum_int = m(mi)%pwt_sum_int + tmp3*tmp1
      m(mi)%pwt_sum_af_int = m(mi)%pwt_sum_af_int + tmp3*tmp1*ice(sc,ic,mi)%af
      m(mi)%pwt_sum_bot = m(mi)%pwt_sum_bot + tmp4*tmp1
      m(mi)%pwt_sum_af_bot = m(mi)%pwt_sum_af_bot + tmp4*tmp1*ice(sc,ic,mi)%af

      if (ice(sc,ic,mi)%af .gt. 0.d0) then

          ! record last snow depth for sub-dt interpolation
          ice(sc,ic,mi)%sh_prev = snowh1m
          !print *,'Actual:    ',snowh1m

          ! record precipitation
          ice(sc,ic,mi)%pr_sum = ice(sc,ic,mi)%pr_sum + f(mi)%pr*dt_s


          do ii=1,sk_z
              ! return to brine conc.
              bv_mean = ice(sc,ic,mi)%bv(ii)*c_001
              ice(sc,ic,mi)%no3(ii) = ice(sc,ic,mi)%no3(ii)/bv_mean
              ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii)/bv_mean
              ice(sc,ic,mi)%po4(ii) = ice(sc,ic,mi)%po4(ii)/bv_mean
              ice(sc,ic,mi)%sioh4(ii) = ice(sc,ic,mi)%sioh4(ii)/bv_mean
              ice(sc,ic,mi)%smalg(ii) = ice(sc,ic,mi)%smalg(ii)/bv_mean
              ice(sc,ic,mi)%poc(ii) = ice(sc,ic,mi)%poc(ii)/bv_mean
          enddo

          ! get older
          ice(sc,ic,mi)%age = ice(sc,ic,mi)%age + dt_years

      endif


       if (boundary_file .eq. 1 .and. monte_carlo .eq. 0) then
           print *,"total_flooded = ",total_flooded
       endif



! ==============================================================================
! END GROW ICE - ICE MODEL 2
! ==============================================================================
