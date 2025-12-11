! file: read_parameters.F90
! Sea Ice Algae Model 2 - Saenz & Arrigo
! Version beta
! ======================================================================
! ======================================================================

! ======================================================================
! Subroutine: sia2_read_parameters
! Purpose: Loads constants from file
! ======================================================================


      SUBROUTINE sia2_read_parameters()
          use sia2_globals
          implicit none

          integer :: i


! Load Constants
! ======================================================================

      open(unit=21, file='constants.txt', form='FORMATTED', access='sequential')
      read(21,*) restart    ! 1) use retart files to resume simulation 2) start new simulation
      read(21,*) stn_only   ! only compute stations, ignoring model domain (mdh1,mdv1, etc.) below
      read(21,*) adv_on        ! advect ice according to ice motion product
      read(21,*) use_gpu    ! perform delta-eddinton light calculations using GPU card
      read(21,*) write_f    ! write frequency (hours)
      read(21,*) write_disk ! actually write to disk - may not want to if debugging (1=yes, 0=no)
      read(21,*) wr_stations! write out station netcdf files - costly disk access, maybe want to supress
      read(21,*) desal        ! desalination type (1 = makut, non-mechanistic, 2 = convective)
      read(21,*) bv_conv    ! critical brine volume at which desal switches from direct convection to relative conv.
      read(21,*) dhdt_desal   ! 1 = allow desal to be calculated using the bottom dhdt equations. 0 = cox/weeks only
      read(21,*) mix_each   ! 1 = desal each layer independently (don't mix bv_conv layers), 0 = mix bv_conv layers
      read(21,*) conv_flux  ! 1 = trasfer heat to convecting brine, 0 = no extra heat
      read(21,*) dbvdt_scale ! scaling factor for desal dilution
      read(21,*) hemi        ! hemisphere - 0 = Southern/Antarctic, 1 = Northern, Arctic
      read(21,*) mdh1        ! model domain - hozirontal upper left (RossSeaBox = 372/744, Full = 202/404)
      read(21,*) mdv1        ! model domain - vertical upper left (RossSeaBox = 424/848, Full = 194/388)
      read(21,*) mdh2        ! model domain - hozirontal lower right (RossSeaBox = 374/748, Full = 518/1036)
      read(21,*) mdv2        ! model domain - vertical lower right (RossSeaBox = 426/752, Full = 526/1052)
      read(21,*) atm_f        !atmospheric forcing source (0=NCEP,1=ECMWF)
      read(21,*) atmo       ! atmospheric boundary calculation type (0=maykut/kottmeier 2003, 1=CICE 4)
      read(21,*) woa        ! toggle for use of specified forcing (0) or world ocean atlas forcing (1)
      read(21,*) ohf        ! ocean heat flux is specificied (0) or derived from world ocean atlas (1)
      read(21,*) use_sose ! toggle for use of SOSE output for ocean fields
      read(21,*) sose_level ! SOSE level for sea surface temp and salinity (1-10)
      read(21,*) boundary_file ! toggle for modeling of (1) extra boundary data, (0) or specified model domain only
      read(21,*) monte_carlo    ! toggle for running purmutations of boundary_file forcing based on modified snowfall,temperature
      read(21,*) override_ic    ! override the built-in dependence on satellite ice conectration, assuming ice the whole run
      read(21,*) icecon_mod    ! modify satellite ice concentration (0=no modification, 1=90%->100%, 2=80%->100%)
      read(21,*) grid_model     ! ice model type (1 = variable layer (no advection), 2 = constant layers)
      read(21,*) snow_model ! snow model type
      read(21,*) snow_ridging ! 1=conserve snow mass dring ridging,0=don't ridge up snow (w/ mass loss)
      read(21,*) alg_dtt    ! algal model calculation frequency (0=once per time step, 1=same frequency as sub-dt ice physics)
      read(21,*) use_ponds  ! don't disappear surface melted ice - retain pond and use Darcy's law for brine drainage
      read(21,*) use_drain  ! 1=drain surface ice if above bv_conv and turn into snow, 0=don't
      read(21,*) hbd_f      ! horizontal brine drainage fraction - fraction of surfance melt that drains away through cracks, etc and not thorugh vertical brine network
      read(21,*) init_ice_th    ! 0. = use satellite-based ice function for initial ice thickness, other positive number = force init ice thickness to this number (m)
      read(21,*) init_ice_type  ! 0 = use satellite-based ice function for inital ice sal profile, 1 = use 1st-year ice profile, 2 = use multi-year ice profile
      read(21,*) init_snow_th   ! 0. = use satellite-based snow initial snow thickness, other positive number = force init snow thickness to this number (m)
      read(21,*) iit        ! intial ice temp = water (0 = use linear gradient w/ airtemp as inital ice temp, 1 = use water temp)
      read(21,*) iis        ! initial ice salinity (0 = use s_cont constant salinity 1 = use standard 1st year ice "C" curve, 2 = multiyear ice curve)
      read(21,*) iin        ! initial ice nutrients (0 = full nutrients upon ice creation, 1 = linear depleted nutrients)
      read(21,*) iif        ! initial snow flooding (m) - used only in boundary_file mode
      read(21,*) s_const    ! salinity constant to use if iis paramater is set to 0
      read(21,*) n_f        ! nutrient fraction - used to scale inital forcing nutirent concentrations if iin set to 2
      read(21,*) ohf_skew   ! ocean heat flux skew.  This number of added to the ocean heat flux (but it is always kept above 0)
      read(21,*) at_skew    ! air temperature flux skew.  This number of added to the air temp
      read(21,*) snow_skew  ! snow depth skew - this number is multiplied by the snow depth
      read(21,*) alg_mig    ! algae migration (0=algae stay put verticaly, 1=algae move with their respective layers while ice grows (quasi-movement))
      read(21,*) bins         ! num of wavelength bins
      read(21,*) dt            ! length of time step (days)
      read(21,*) dt_sub_1   ! fraction of main timestep for ice physics
      read(21,*) dt_sub_2   ! fraction of main timestep for fast-changing ice physics
      read(21,*) dt_sub_3   ! fraction of main timestep for very-changing ice physics
      read(21,*) dt_sub_1_tol   ! minimum temperature tolerance for use of dt_sub_1 (deg C)
      read(21,*) dt_sub_2_tol   ! minimum temperature tolerance for use of dt_sub_2 (deg C)
      read(21,*) Ek_max        ! Spectral photoaclimation parameter (microEin*m-2*s-1)
      read(21,*) A             ! parameter of light utilization equation (dimensionless)
      read(21,*) B             ! parameter of light utilization equation (dimensionless)
      read(21,*) begin_j_date ! start date (decimal julian day)
      read(21,*) begin_year ! start year (year)
      read(21,*) end_j_date ! start date (decimal julian day)
      read(21,*) end_year   ! start year (year)
      read(21,*) xp         ! phytoplankton death rate (dimensionless)
      read(21,*) remin      ! rate of poc/detritus remineralization (1/day)
      read(21,*) remin_f    ! fraction poc remineralization to available N,P
      read(21,*) c_chl         ! Carbon:Chlorophyl ratio (grams)
      read(21,*) c_n         ! Carbon:Nitrogen Ratio (dimensionless)
      read(21,*) c_si         ! Carbon:Silicon Ratio (dimensionless)
      read(21,*) c_p         ! Carbon:Phosphorus Ratio (dimensionless)
      read(21,*) Ks_NH4        ! Half-saturation rate constant (microMolar)
      read(21,*) Ks_NO3     ! Half-saturation rate constant (microMolar)
      read(21,*) Ks_SiO4     ! Half-saturation rate constant (microMolar)
      read(21,*) Ks_PO4        ! Half-saturation rate constant (microMolar)
      read(21,*) Fw         ! Oceanic heat flux from water (cal/m^2/sec)
      read(21,*) Sw          ! salinty of seawater (ppt)
      read(21,*) Sd         ! density of seawater  (g/m^3)
      read(21,*) IceD          ! density of pure ice (g/m^3)
      read(21,*) rg            ! growth rate constant rg (degC^-1) - from Epply et al 1972
      read(21,*) G0         ! growth rate @ zero dec C (d^-1) - from Epply at al 1972
      read(21,*) z_sk        ! number of modeled layers in skeletal layer
      read(21,*) sk_th_min    ! minimum height of iindividual skeletal layer
      read(21,*) sk_th_max    ! maximum height of iindividual skeletal layer
      read(21,*) z_int_min    ! defines minimum number of ice model layers
      read(21,*) z_th_min    ! minimum height of iindividual layer
      read(21,*) z_th_max    ! maximum height of iindividual layer
      read(21,*) z_th_crit  ! critical height of individual layer, below which a smaller physics timestep mut be used (m)
      read(21,*) max_snow_h ! maximum height of snow over sea ice (m) - used b/c performance of snow data degrades above 60cm
      read(21,*) Tw            ! temperature of freezing of seawater
      read(21,*) ksnow        ! standard snow conductivity (Cal/m/s/degC)
      read(21,*) f_limit    ! delta_T * wind_speed limit, above which no ice grows in model (Km/s)
      read(21,*) f_sk        ! fraction of skeletal layer open to convection
      read(21,*) fb            ! fraction of brine tube layer that is brine tubes
      read(21,*) fi            ! minimum fractional sea ice coverage neccessary for algal growth
      read(21,*) sw_NH4        ! seawater NH4 concentration, used if climatologies not used (然olar)
      read(21,*) sw_NO3        ! seawater N03 concentration, used if climatologies not used (然olar)
      read(21,*) sw_SiO4    ! seawater SiO2 concentration, used if climatologies not used (然olar)
      read(21,*) sw_PO4        ! seawater PO4 concentration, used if climatologies not used (然olar)
      read(21,*) sw_poc     ! seawater detritus concentration (g/m^3)
      read(21,*) vb_crit    ! brine volume above which gravity drainage occurs (ppt)
      read(21,*) vb_max        ! brine volume above which layer thickness will be reduced/melted (ppt)
      read(21,*) conv_max   ! maximum convective flux for use in desalination (cm^3/cm^2/s)
      read(21,*) dhdt_i        ! interval for ice growth calculation (hours)
      read(21,*) ncep_f     ! ncep/doe II forcing frequency (6 = every 6 hours, 24 = daily)
      read(21,*) min_alg    ! minimum microalgal concentration (mg/m^3)
      read(21,*) alg_wc     ! water column algal concentration - used when freezing/flooding/creating new ice (mgC/m^3)
      read(21,*) min_sal    ! minimum bulk ice salinity (ppt)
      read(21,*) nr_tol     ! Newton-Raphson method tolerance for surface temp (T0)
      read(21,*) max_it        ! Newton-Raphson maximum iterations
      read(21,*) alb_s_dry     ! standard albedo dry snow (%)
      read(21,*) alb_s_wet     ! standard albedo wet snow (%)
      read(21,*) alb_i_dry     ! standard albedo of sea ice (%)
      read(21,*) alb_i_wet     ! standard albedo of sea ice (%)
      read(21,*) den_s_dry  ! density dry snow (g/cm^3)
      read(21,*) den_s_wet  ! density wet snow (g/cm^3)
      read(21,*) des_s_switch ! temperature above which wet snow falls, below which dry snow falls
      read(21,*) eps_snow   ! long wave emissivity of snow (%)
      read(21,*) eps_ice    ! long wave emissivity of ice (%)
      read(21,*) alg_mig_crit !maximum ice growth rate under which algae maintain position (cm/day)
      read(21,*) gl_max_f      ! fraction of minimum layer height that is the trigger for boundary adjustment
      read(21,*) fl_max_f      ! fraction of minimum layer height that is the trigger for boundary adjustment for flooding
      read(21,*) T_ref        ! temp from which all heat intergrals are calculated (degC) (hopefully the lower than the lowest ice temp)
      read(21,*) temp_tol    ! temp difference tolerance for heat methods - if the difference in temp get smaller than this, take corrective action
      read(21,*) pr_on        ! on/off toggle for causing precipitation to go directly into snow ice formation, if appropriate
      read(21,*) bb_f         ! bubble fraction in sea ice - causes density to be decreased by this amount
      read(21,*) fl_crit    ! minimum freeboard required to flood ice surface (m)
      read(21,*) fl_ratio    ! ratio of ice/water in flooded snow
      read(21,*) woa_depth    ! level of World Ocean Altas data to use as forcing (1=surface,2=10m,3=20m,4=30m,5=50m,6=75m,7=100m)
      read(21,*) h_snowpatch ! contant that determines the bare ice/snow covered ice ratio for a given snowdepth
      read(21,*) par_to_swd ! energy conversion factor between PAR (300-700nm) irradiance to total shortwave irradiance (250-4000nm) (energy/energy units)
      read(21,*) a_ice_ir   ! weighted-mean ice absorption coefficient for near-IR (700-4000nm)
      read(21,*) da_f       ! distribution adjustment factor - effectively changes the s.d. of the log-normal distribution - value between 1 and 0
      read(21,*) snow_min   ! minimum modeled snow depth (m) - below this value, the model assumes no snow for thermodynamic/light/albedo putposes
      read(21,*) gc_offset  ! offset in hours for indexing gc_par_ice file (240*24=5760 for arrigo 1982 sim)
      read(21,*) par_cf     ! use climatology could fraction to correct par (0=off (cf already included in par), 1=yes)
      read(21,*) use_mdiff  ! use molecular diffusion to supply nutrients to ice bottom (0=off 1=on)
      read(21,*) ts_is_at   ! short-circuit atmospheric heat transfer by fixing surface temp to air temp (0=off 1=on)
      read(21,*) kevin_light ! use Kevin's light model instead of gregg&carder 1990 (0=off 1=on)
      read(21,*) use_pl     ! use a platelet layer - only works if validaiton is on (0=off 1=on)
      read(21,*) snow_fudge ! snow depth fudge factor for multiplicatively adjecting the snow depth for light purposes
      read(21,*) r_s_dry    ! scattering coefficient (1/m) for dry snow
      read(21,*) r_s_wet    !  scattering coefficient (1/m) for wet snow
      read(21,*) r_i_23     ! scattering coefficient (1/m) for very cold ice
      read(21,*) r_i_8      ! scattering coefficient (1/m) for cold ice (zero if not used)
      read(21,*) min_scat   ! minimum scattering coefficient (1/m)
      read(21,*) rs_switch  ! basis for determining snow scattering (0 = airtemp>0, 1 = layer_temp>-1, 2 = melt variable)
      read(21,*) a_factor   ! wavelength independent snow absorption for delta eddington light routine
      read(21,*) use_expt   ! 1=perform exp lookup table approximation,0=system exp()
      read(21,*) exp_bins   ! number of exp lookup table bins
      read(21,*) exp_a_min  ! minimum exponential function argument
      read(21,*) no_flooding ! prevents flooding if non-zero
      read(21,*) flood_brine ! (1) pushed brine upward or (2) floods directly with seawater
      read(21,*) p_factor   ! ice strength scalar, used to compare ice strength to ice momentum during advection refinement
      read(21,*) snow_rd_lim ! number that determines when snow depth distribution is re-distributed
      read(21,*) cv_void_f   ! convergence void fraction
      read(21,*) cvf_switch  ! ice height below which cv_void_f is affected by cvf_thin (m)
      read(21,*) cvf_thin    ! multiplier of cv_void_f when ice is less than cvf_switch (fraction)
      read(21,*) snow_in_gaps ! determines whether snow is inserted into ridging gaps, or just seawater
      close(21)

      open(unit=21, file='aph.txt', form='FORMATTED', access='sequential')
      aph_max = 0.
      do i=1,301
          read(21,*) aph(i)
          if (aph(i) .gt. aph_max) then
              aph_max = aph(i)
          endif
      enddo
      close(21)

      open(unit=21, file='awater.txt', form='FORMATTED', access='sequential')
      do i=1,301
          read(21,*) awater(i)
      enddo
      close(21)

      open(unit=21, file='awater_tc.txt', form='FORMATTED', access='sequential')
      do i=1,301
          read(21,*) awater_tc(i)
      enddo
      close(21)

      open(unit=21, file='awater_sc.txt', form='FORMATTED', access='sequential')
      do i=1,301
          read(21,*) awater_sc(i)
      enddo
      close(21)

      open(unit=21, file='awater_v.txt', form='FORMATTED', access='sequential')
      do i=1,301
          read(21,*) awater_v(i)
      enddo
      close(21)

      open(unit=21, file='aice.txt', form='FORMATTED', access='sequential')
      do i=1,301
          read(21,*) aice(i)
      enddo
      close(21)

      open(unit=21, file='aozone.txt', form='FORMATTED', access='sequential')
      do i=1,301
          read(21,*) aozone(i)
      enddo
      close(21)

      open(unit=21, file='surfacelight.txt', form='FORMATTED', access='sequential')
      do i=1,301
          read(21,*) surfacelight(i)
          surfacelight(i)=surfacelight(i)*dble(i+399)*.008355 ! convert from w/m^2 to 湲in/m^2/s
      enddo
      close(21)

      end SUBROUTINE sia2_read_parameters
