! ======================================================================
! Global Paramaters, Types, & Variables
! ======================================================================

      module sia2_globals

! ======================================================================
! Parameters
! ======================================================================

          ! ------------------------------------------------------------
          ! sizes of arrays, and array indexes
          ! ------------------------------------------------------------
          integer, parameter :: ic_n = 1
          integer, parameter :: sc_n = 1          ! currently, either 1 or 2
          integer, parameter :: z_max = 42
          integer, parameter :: z_snow_max = 26
          integer, parameter :: la = 1
          integer, parameter :: sd_num = 9      ! resolution of lognrnmal snow depth distrbution

          integer, parameter :: pl_max = 30
          double precision, parameter :: pl_th = 0.02d0
          double precision, parameter :: h_max = 10.d0  ! 10m ice max
          double precision, parameter :: bot_th = 0.2d0  ! thickness of bottom section considered bottom ice (m from bottom)

          double precision, parameter :: pond_f_perc = 0.3d0 ! fraction of ponded/melt water that pushed down through brine network

          integer, parameter :: mp_x=192
          integer, parameter :: mp_y=94

          integer, parameter :: ec_x=144
          integer, parameter :: ec_y=73

          integer, parameter :: eci_x=240
          integer, parameter :: eci_y=121

          integer, parameter :: dg_x=360
          integer, parameter :: dg_y=180

          integer, parameter :: wavl=31
          integer, parameter :: logn_n = 9

          ! ------------------------------------------------------------
          ! input data resolution based on southern hemisphere EASE grid
          ! ------------------------------------------------------------
          integer, parameter :: grid_v=721
          integer, parameter :: grid_h=721
          integer, parameter :: grids_v=321
          integer, parameter :: grids_h=321
          integer, parameter :: h1_ncep_i=202
          integer, parameter :: h2_ncep_i=518
          integer, parameter :: v1_ncep_i=194
          integer, parameter :: v2_ncep_i=526

          ! ------------------------------------------------------------
          ! input data resolution based on northern hemisphere EASE grid
          ! ------------------------------------------------------------
          integer, parameter :: gridn_v=361
          integer, parameter :: gridn_h=361

          ! ------------------------------------------------------------
          ! Enumerations
          ! ------------------------------------------------------------
          integer, parameter :: ncep_at=1
          integer, parameter :: ncep_p=2
          integer, parameter :: ncep_h=3
          integer, parameter :: ncep_fc=4
          integer, parameter :: ncep_u10=5
          integer, parameter :: ncep_v10=6
          integer, parameter :: ncep_pr=7

          integer, parameter :: woa_t=1
          integer, parameter :: woa_s=2
          integer, parameter :: woa_n=3
          integer, parameter :: woa_p=4
          integer, parameter :: woa_si=5

          ! ------------------------------------------------------------
          ! Empirical Constants
          ! ------------------------------------------------------------
          double precision, parameter :: pi=3.141592d0
          double precision, parameter :: steph_boltz=5.66d-8 !W/m^2/K^4
          double precision, parameter :: gC_mC = 12.01d0 !gramsC/molesC
          double precision, parameter :: mC_gC = 1.0d0/12.01d0 !molesC/gramC
          double precision, parameter :: kelvin0 = 273.16d0 !Kelvin-Celsius offset
          double precision, parameter :: mu = -0.054d0 ! linear freezing point coefficient for sea water
          double precision, parameter :: cell_side = 25.067525d0 ! km
          double precision, parameter :: cell_area = cell_side**2 ! km^2
          double precision, parameter :: atan_max = 1.557407724654902d0 ! = tan(1), used for atan function mulitplier
          double precision, parameter :: grav = 9.80616d0 ! gravity, m/s^2
          double precision, parameter :: Lv = 2501.0d0 ! kJ/kg latent heat of vaporization of water (2260)
          double precision, parameter :: Lf = 334.0d0 ! kJ/kg latent heat of fusion of water (334)
          double precision, parameter :: cw = 3.96d0 ! specific heat of seawater around 0 deg C - J/g/k
          double precision, parameter :: c0 = 2.113d0   ! specific heat of pure ice - J/g/k
          double precision, parameter :: R_air = 287.058d0  ! Universal gas constant for dry air - J/kg/K; pressure must be in Pascals for: density = pressure/RT
          double precision, parameter :: R_h2o = 461.495d0 ! Universal gas constant for water vapor - J/kg/K; pressure must be in Pascals for: density = pressure/RT

          ! ------------------------------------------------------------
          ! Empirical Parameterizations
          ! ------------------------------------------------------------
          double precision, parameter :: Ce = 2.1d-3      ! coefficient of turbulent latent heat transfer (water vapor?)
          double precision, parameter :: Ch = 2.0d-3      ! coefficient of turbulent sensible heat transfer
          double precision, parameter :: cp = 1005.0d0      ! J/kg/K
          double precision, parameter :: eps0_s = 0.97d0        ! long-wave emissivity of snow surface -
          double precision, parameter :: eps0_i = 0.99d0        ! long-wave emissivity of ice surface -
          double precision, parameter :: fe_A = 2.7798d-6   ! hPa/K^4, turbutent heat flux coefficient
          double precision, parameter :: fe_B = -2.6913d-3  ! hPa/K^3, turbutent heat flux coefficient
          double precision, parameter :: fe_C = 0.97920d0     ! hPa/K^2, turbutent heat flux coefficient
          double precision, parameter :: fe_D = -158.64d0     ! hPa/K, turbutent heat flux coefficient
          double precision, parameter :: fe_E = 9653.2d0      ! hPa, turbutent heat flux coefficient
          double precision, parameter :: qq1 = 1.16378d7      ! q1 constant in turbulent latentent heat flux calc, CICE v. 4
          double precision, parameter :: qq2 = 5897.8d0       ! q2 constant in turbulent latentent heat flux calc, CICE v. 4

          double precision, parameter :: ss0 = 8.0d0     ! 1st year ice salinilty standard curve value (top)
          double precision, parameter :: ss1 = 6.3d0     ! 1st year ice salinilty standard curve value (mid)
          double precision, parameter :: ss2 = 5.6d0     ! 1st year ice salinilty standard curve value (mid)
          double precision, parameter :: ss3 = 5.3d0    ! 1st year ice salinilty standard curve value (mid)
          double precision, parameter :: ss4 = 5.2d0     ! 1st year ice salinilty standard curve value (mid)
          double precision, parameter :: ss5 = 5.1d0     ! 1st year ice salinilty standard curve value (mid)
          double precision, parameter :: ss6 = 4.9d0     ! 1st year ice salinilty standard curve value (mid)
          double precision, parameter :: ss7 = 4.8d0     ! 1st year ice salinilty standard curve value (mid)
          double precision, parameter :: ss8 = 4.8d0     ! 1st year ice salinilty standard curve value (mid)
          double precision, parameter :: ss9 = 6.2d0     ! 1st year ice salinilty standard curve value (bottom)

          double precision, parameter :: ssm0 = 0.1d0     ! multi-year salinilty standard curve value (top)
          double precision, parameter :: ssm1 = 0.2d0     ! multi-year salinilty standard curve value (mid)
          double precision, parameter :: ssm2 = 0.2d0     ! multi-year salinilty standard curve value (mid)
          double precision, parameter :: ssm3 = 0.6d0     ! multi-year salinilty standard curve value (mid)
          double precision, parameter :: ssm4 = 1.9d0     ! multi-year salinilty standard curve value (mid)
          double precision, parameter :: ssm5 = 3.1d0     ! multi-year salinilty standard curve value (mid)
          double precision, parameter :: ssm6 = 3.3d0     ! multi-year salinilty standard curve value (mid)
          double precision, parameter :: ssm7 = 3.4d0     ! multi-year salinilty standard curve value (mid)
          double precision, parameter :: ssm8 = 3.9d0     ! multi-year salinilty standard curve value (mid)
          double precision, parameter :: ssm9 = 6.2d0     ! multi-year salinilty standard curve value (bottom)

          double precision, parameter :: a0 = -0.12073d0    ! salinity limiting poly coefficient 0 (dimensionless)
          double precision, parameter :: a1 = 0.07097d0     ! salinity limiting poly coefficient 1 (dimensionless)
          double precision, parameter :: a2 = -0.00133d0    ! salinity limiting poly coefficient 2 (dimensionless)
          double precision, parameter :: a3 = 6.3427d-06  ! salinity limiting poly coefficient 3 (dimensionless)

          double precision, parameter :: a_star = 0.05d0  ! mean value of areal fraction participating in ridging
          double precision, parameter :: r_a_star = 20.d0  ! 1/a_star
          double precision, parameter :: adv_mu = 4.d0  !  used to calc lambda constant in ice category redistribution function

          double precision, parameter :: ki_min = 0.563d0  !  minimum sea ice thermal conductivity (W/m/K)

          ! seawater equation of state coefficient alpha
          real, parameter :: d_a(7) = (/ -1.36471e-1, 4.68181e-2, 8.07004e-1, -7.45353e-3, -2.94418e-3, 3.43570e-5, 3.48658e-5/)
          ! seawater equation of state coefficient beta
          real, parameter :: d_b(7) = (/ 5.06423e-1, -3.57109e-3, -8.76148e-4, 5.25243e-5, 1.57976e-5, -3.46686e-7, -1.68764e-7/)
          ! seawater equation of state coefficient gamma
          real, parameter :: d_g(7) = (/ -5.52640e-4, 4.88584e-6, 9.96027e-7, -7.25139e-8, -3.98736e-9, 4.00631e-10, 8.26368e-11/)

          double precision, parameter :: logn_multiplier9(9) = (/0.102d0, 0.272d0, 0.427d0, 0.532d0, 0.721d0, &
          0.952d0, 1.31d0, 1.74d0, 3.31d0/)

          double precision, parameter :: logn_multiplier6(6) = (/0.145d0, 0.385d0, 0.585d0, 0.860d0, 1.345d0, 2.70d0/)
          double precision, parameter :: logn_multiplier5(5) = (/1.82d-1, 4.08d-1, 6.96d-1, 1.14d0, 2.59d0/)

          double precision, target :: gauss_bins_3(4) = (/0.d0, 1.3717d-01, 3.0777d-01, 1.d0/)
          double precision, target :: gauss_bins_4(5) = (/0.d0, 1.0153d-01, 2.1483d-01, 3.6633d-01, 1.d0/)
          double precision, target :: gauss_bins_5(6) = (/0.d0, 8.0522d-02, 1.6677d-01, 2.6766d-01, 4.0738d-01, 1.d0/)
          double precision, target :: gauss_bins_6(7) = (/0.d0, 6.7155d-02, 1.3749d-01, 2.1515d-01, 3.0840d-01, 4.4080d-01, 1.d0/)
          double precision, target :: gauss_bins_7(8) = (/0.d0, 5.7288d-02, 1.1649d-01, &
                                      1.8014d-01, 2.5207d-01, 3.3991d-01, 4.6626d-01, 1.d0/)
          double precision, target :: gauss_bins_8(9) = (/0.d0, 5.0286d-02, 1.0185d-01, &
                                      1.5595d-01, 2.1515d-01, 2.8294d-01, 3.6696d-01, 4.8950d-01, 1.d0/)
          double precision, target :: gauss_bins_9(10) = (/0.d0, 4.4558d-02, 9.0070d-02, &
                                      1.3749d-01, 1.8810d-01, 2.4411d-01, 3.0872d-01, 3.8956d-01, 5.0859d-01, 1.d0/)
          double precision, target :: gauss_bins_10(11) = (/0.d0, 4.0102d-02, 8.0840d-02, &
                                      1.2285d-01, 1.6709d-01, 2.1483d-01, 2.6798d-01, 3.3004d-01, &
                                      4.0802d-01, 5.2355d-01, 1.d0/)

          double precision, parameter :: c9999 = 9999.0d0
          double precision, parameter :: c_5 = 0.5d0
          double precision, parameter :: cn_5 = -0.5d0
          double precision, parameter :: c_1 = 0.1d0
          double precision, parameter :: c_01 = 0.01d0
          double precision, parameter :: c_001 = 0.001d0
          double precision, parameter :: c1_5 = 3.0d0/2.0d0
          double precision, parameter :: c_333 = 1.0d0/3.0d0
          double precision, parameter :: c_111 = 1.0d0/9.0d0

          double precision, parameter :: d0_ = 0.0d0

          double precision, parameter :: pur_c = 0.0246059207996924d0
          double precision, parameter :: one_over_pur_c = 40.6406250000000d0
          double precision, parameter :: pur_c_2 = 0.0244140625000000d0
          double precision, parameter :: one_over_pur_c_2 = 40.96d0
!          integer (kind=1), parameter :: pur_0 = -127
          integer (kind=2), parameter :: pur_0 = -32767
          double precision, parameter :: af_min = 5.d-9

          character*20, save :: version_string = 'beta - initial      '
          double precision, allocatable :: ida_multiplier(:)
          double precision, allocatable :: lda_multiplier(:)
          double precision, allocatable :: sda_multiplier(:)


          double precision, parameter :: norwegian_sal(15) = (/1.1d0,0.1d0,0.2d0,2.5d0,2.6d0,2.2d0,3.3d0,4.8d0, &
                                         3.4d0,3.3d0,3.1d0,4.0d0,4.0d0,3.8d0,4.9d0/)
          double precision, parameter :: norwegian_dth_s(15) = (/0.d0,0.05d0,0.15d0,0.25d0,0.35d0,0.45d0,0.55d0, &
                                         0.65d0,0.75d0,0.85d0,0.95d0,1.05d0,1.15d0,1.25d0,1.35d0/)

          double precision, parameter :: norwegian_N(16) = (/1.2133d0,1.6352d0,0.7084d0,1.1116d0,0.580d0,0.400d0, &
          0.422d0,0.436d0,0.574d0,1.792d0,0.690d0,0.400d0,0.400d0,0.4352d0,3.1164d0,16.80d0/)
          double precision, parameter :: norwegian_Si(16) = (/0.70d0,0.70d0,0.70d0,0.70d0,0.70d0,0.7d0,0.70d0, &
          0.70d0,0.70d0,0.70d0,0.70d0,0.70d0,0.70d0,0.70d0,0.70d0,0.70d0/)
          double precision, parameter :: norwegian_smalg(16) = (/0.166238571428571d0,0.305519047619048d0, &
          0.472243333333333d0,0.797851428571429d0,1.12345952380952d0,1.07686000000000d0,0.716543333333333d0, &
          0.383952380952381d0,0.470720476190476d0,0.557488095238095d0,0.505974761904762d0,0.354899047619048d0, &
          0.214843333333333d0,0.173122857142857d0,0.147333333333333d0,0.171904761904762d0/)
          double precision, parameter :: norwegian_dth(16) = (/0.d0,0.05d0,0.15d0,0.25d0,0.35d0,0.45d0,0.55d0, &
                                         0.65d0,0.75d0,0.85d0,0.95d0,1.05d0,1.15d0,1.25d0,1.35d0,1.4d0/)


! ======================================================================
! Global Types
! ======================================================================

          ! ------------------------------------------------------------
          ! Static array types (forcing & projection data)
          ! ------------------------------------------------------------

          ! Projection type - holds projection conversion information for various projections
          type, public :: proj_type
              integer, dimension(grid_h,grid_v) :: &
                  mdh, &          ! model domain horizontal grid index
                  mdv, &          ! model domain vertical grid index
                  x_mp, &         ! NCEP/mercator grid horizontal(longitude) index
                  y_mp, &         ! NCEP/mercator grid vertical(latitude) index
                  x_ec, &         ! ECMWF 2.5 grid horizontal(longitude) index
                  y_ec, &         ! ECMWF 2.5 grid vertical(latitude) index
                  x_eci, &        ! ECMWF 1.5 grid horizontal(longitude) index
                  y_eci, &        ! ECMWF 1.5 grid vertical(latitude) index
                  x_dg, &         ! 1 degree grid horizontal(longitude) index
                  y_dg, &         ! 1 degree grid vertical(latitude) index
                  mi, &           ! grid index to modeled variables
                  mask            ! ocean mask - 1 = ocean, 0 = excluded from analysis (land or coastline)
              double precision, dimension(grid_h,grid_v) :: &
                  lat, &          ! central latitude of model grid cell
                  lon             ! central longitude of model grid cell
          end type proj_type

          ! NCEP type - holds an entire set of NCEP forcing information
          type, public :: ncep_type
              ! air temp, 2m air pressure, specific humidity, cloud fraction,
              ! u direction wind speed, v direciton wind speed, precipitation rate
              double precision, dimension(mp_x,mp_y) :: &
                  at, &           ! 2m air temperature (kelvin)
                  p, &            ! surface air pressure (mb)
                  h, &            ! specific humidity ()
                  fc, &           ! fraction cloud cover (%)
                  u10, &          ! 10m horizontal wind speed vector (m/s)
                  v10, &          ! 10m vertical wind speed vector (m/s)
                  pr              ! precipitation rate (kg/m^2/s)
          end type ncep_type

          ! ed_type - holds a set of global spectral radiation based on the
          ! NCEP/DOE grid
          type, public :: ed_type
              double precision, dimension(mp_x,mp_y,wavl) :: &
                  Ed              ! 湲in/m^2/s
          end type ed_type

          ! mp_f_type - Mercator projection forcing type, holds all NCEP forcing data
          ! using the ~2 degree NCEP grid
          type, public :: mp_f_type
              type (ncep_type) :: &
                  ncep            ! current NCEP forcing dataset
              type (ncep_type) :: &
                  ncep_next       ! next available NCEP forcing dataset
              type (ncep_type) :: &
                  ncep_interp     ! NCEP forcing interpolated for current time step
              type (ed_type), dimension(2) :: &
                  Edir, &         ! Direct irradience (湲in/m^2/s)
                  Edif            ! Diffuse irradience (湲in/m^2/s)
          end type mp_f_type

          ! ECMWF type - holds an entire set of ECMWF forcing information
          type, public :: ecmwf_type
              ! air temp, 2m air pressure, specific humidity, cloud fraction,
              ! u direction wind speed, v direciton wind speed, precipitation rate
              double precision, dimension(ec_x,ec_y) :: &
                  at, &           ! 2m air temperature (kelvin)
                  p, &            ! surface air pressure (Pa)
                  dpt, &          ! dew point temperature (kelvin)
                  fc, &           ! fraction cloud cover (%)
                  u10, &          ! 10m horizontal wind speed vector (m/s)
                  v10, &          ! 10m vertical wind speed vector (m/s)
                  pr, &           ! total precipitation (m)
                  ssr             ! surface solar radiation (W/m^2)
          end type ecmwf_type

          ! ECMWF Interim type - holds an entire set of ECMWF forcing information
          type, public :: ecmwf_int_type
              ! air temp, 2m air pressure, specific humidity, cloud fraction,
              ! u direction wind speed, v direciton wind speed, precipitation rate
              double precision, dimension(eci_x,eci_y) :: &
                  at, &           ! 2m air temperature (kelvin)
                  p, &            ! surface air pressure (Pa)
                  dpt, &          ! dew point temperature (kelvin)
                  fc, &           ! fraction cloud cover (%)
                  u10, &          ! 10m horizontal wind speed vector (m/s)
                  v10, &          ! 10m vertical wind speed vector (m/s)
                  pr, &           ! total precipitation (m)
                  ssr             ! surface solar radiation (W/m^2)
          end type ecmwf_int_type

          ! ec_f_type - ECMWF projection forcing type, holds all ECMWF forcing data
          ! using the ECMWF grid
          type, public :: ec_f_type
              type (ecmwf_type) :: &
                  ecmwf           ! current NCEP forcing dataset
              type (ecmwf_type) :: &
                  ecmwf_next      ! next available NCEP forcing dataset
              type (ecmwf_type) :: &
                  ecmwf_interp    ! NCEP forcing interpolated for current time step
          end type ec_f_type

          ! eci_f_type - ECMWF Interim projection forcing type, holds all ECMWF Interim forcing data
          ! using the ECMWF grid
          type, public :: eci_f_type
              type (ecmwf_int_type) :: &
                  ecmwf           ! current NCEP forcing dataset
              type (ecmwf_int_type) :: &
                  ecmwf_next      ! next available NCEP forcing dataset
              type (ecmwf_int_type) :: &
                  ecmwf_interp    ! NCEP forcing interpolated for current time step
          end type eci_f_type

          ! woa_type - holds 1 degrees World Ocean Atlas data for a specific depth,
          ! which is determined by setting the woa_depth constant
          type, public :: woa_type
              real, dimension(dg_x,dg_y) :: &
                  t, &            ! World Ocean Atlas water temp (dec C)
                  s, &            ! World Ocean Atlas salinity (psu)
                  n, &            ! World Ocean Atlas nitrate concentration (然ol)
                  p, &            ! World Ocean Atlas phosphate concentration (然ol)
                  si              ! World Ocean Atlas silica concentration (然ol)
          end type woa_type

          ! degree grid forcing type - holds forcing datasets available in a
          ! 1 degree grid
          type, public :: dg_f_type
              type (woa_type) :: &
                  woa             ! current WOA forcing dataset
              type (woa_type) :: &
                  woa_next        ! next available WOA forcing dataset
              type (woa_type) :: &
                  woa_interp      ! WOA forcing interpolated for current time step
          end type dg_f_type

          ! sose_type - holds ease-gridded SOSE output
          type, public :: sose_type
              double precision, dimension(grids_h,grids_v) :: &
                  t, &            ! temperature (degC)
                  s, &            ! salinity (psu)
                  u, &            ! u-velocity (m/s)
                  v             	! v-velocity (m/s)
          end type sose_type

          ! sose_f - holds SOSE forcing data
          type, public :: sose_f_type
              type (sose_type) :: &
                  sose             ! current sose forcing dataset
              type (sose_type) :: &
                  sose_next        ! next available sose forcing dataset
              type (sose_type) :: &
                  sose_interp      ! sose forcing interpolated for current time step
          end type sose_f_type


          ! EASE projection forcing type - holds forcing variables that are
          ! projected using ease grid
          type, public :: ease_f_type
              double precision, dimension(grid_h,grid_v) :: &
                  icecon, &       ! ssm/i ice concentration for current timestep (fraction)
                  icecon_next     ! pre-loaded ssm/i ice concentration for next step (fraction)
              double precision, dimension(grid_h,grid_v) :: &
                  sh, &           ! ssm/i snow depth for current timestep (cm)
                  sh_next         ! pre-loaded ssm/i snow depth for next step (cm)
              double precision, dimension(grid_h,grid_v,7) :: &
                  sh_past         ! array of past snow depths for validaiton of snow depth algorithm
              double precision, dimension(grid_h,grid_v,5) :: &
                  icecon_past         ! array of past snow depths for validaiton of snow depth algorithm
              double precision, dimension(grid_h,grid_v,3) :: &
                  ice_vec, &       ! ease grid ice vectors (u,v,flag)
                  ice_vec_next     ! pre-loaded ease grid ice vectors for next step(u,v,flag)
              double precision, dimension(grid_h,grid_v,3,5) :: &
                  icevec_past
          end type ease_f_type

          type, public :: boundary_f_type ! boundary_file station forcing type
              double precision :: &
                  at, &
                  sh, &
                  ih, &
                  p, &
                  lat, &
                  lon, &
                  iceh1, &
                  ws, &
                  rhum, &
                  fc, &
                  shum, &
                  sens_heat_flux, &
                  late_heat_flux, &
                  pr_snow, &
                  pr
              double precision :: &
                  fw, &
                  t, &
                  s, &
                  no3, &
                  nh4, &
                  po4, &
                  sioh4, &
                  swd, &
                  lwd
              integer :: &
                  time, &
                  time_offset, &
                  x_mp, &
                  y_mp, &
                  x_ec, &
                  y_ec, &
                  x_eci, &
                  y_eci, &
                  x_dg, &
                  y_dg, &
                  length, &
                  begin_year, &
                  mdh, &
                  mdv
          end type boundary_f_type

          ! ------------------------------------------------------------
          ! Dynamically allocated array types (modeled data)
          ! ------------------------------------------------------------

          ! Interpolated/corrected forcing data is stored for one modeled grid point
          ! is stored in this type
          type, public :: forcing_type
              double precision :: &
                  at, &
                  p, &
                  h, &
                  fc, &
                  u10, &
                  v10, &
                  pr, &
                  ws, &
                  rhum, &
                  swd, &
                  lwd, &
                  !sens_heat_flux, &
                  !late_heat_flux, &
                  fw, &
                  t, &
                  s, &
                  d, &
                  no3, &
                  nh4, &
                  po4, &
                  sioh4, &
                  poc, &
                  sh_interp, &
                  sd_interp, &
                  sh_new, &
                  sh_interp_next, &
                  sd_interp_next, &
                  ic_interp, &
                  ic_interp_next, &
                  sh_interp_last, &
                  sd_interp_last, &
                  ivu_interp, &
                  ivv_interp, &
                  ivu_interp_next, &
                  ivv_interp_next
          end type forcing_type

          type, public :: meta_type
              integer :: &
                  status, &
                  grid_h, &
                  grid_v, &
                  x_mp, &
                  y_mp, &
                  x_ec, &
                  y_ec, &
                  x_eci, &
                  y_eci, &
                  x_dg, &
                  y_dg
              double precision :: &
                  melt_loss, &
                  adv_loss, &
                  md_loss, &
                  cong_growth, &
                  snow_growth, &
                  void_growth, &
                  adv_gain, &
                  md_gain, &
                  a_convg, &
                  a_new, &
                  a_drop, &
                  lat, &
                  lon, &
                  pr_clim, &
                  pr_ssmi, &
                  pxl_h_offset, &
                  pxl_v_offset, &
                  prod_sum_int, &
                  prod_sum_bot, &
                  bm_lost, &
                  tlim, &
                  llim, &         ! light limitation fraction (dimensionless)
                  nlim, &         ! nitrogen limitation fraction (dimensionless)
                  plim, &         ! phosphorus limitation fraction (dimensionless)
                  silim, &        ! silica limitation fraction (dimensionless)
                  slim, &        ! salinity limitation fraction (dimensionless)
                  salt_flux, &
                  h2o_flux, &
                  p_wgt_int, &
                  p_wgt_af_int, &
                  p_wgt_bot, &
                  p_wgt_af_bot, &
                  pwid_sum_int, &
                  pwid_sum_af_int, &
                  pwsd_sum_int, &
                  pwsd_sum_af_int, &
                  pwt_sum_int, &
                  pwt_sum_af_int, &
                  pwid_sum_bot, &
                  pwid_sum_af_bot, &
                  pwsd_sum_bot, &
                  pwsd_sum_af_bot, &
                  pwt_sum_bot, &
                  pwt_sum_af_bot

              integer, dimension(8) :: &
                  ia, &           ! i_ease index of 8 surround cells
                  ja, &           ! j_ease index of 8 surrounding cells
                  mia             ! mi index of 8 surroundind cells
          end type meta_type

          type, public :: snow_type
              integer :: &
                  z
              double precision :: &
                  depth, &
                  ts
              double precision, dimension(z_max) :: &
                  t, &
                  d, &
                  heat, &
                  th, &
                  melt, &
                  ed_w
          end type snow_type

          type, public :: pond_type
              double precision :: &
                  t, &
                  s, &
                  d, &
                  heat, &
                  th, &
                  perc, &
                  smalg, &        ! brine-based algal concentration (mgC/m^3)
                  poc, &          ! particulate organic carbon (detritus) (mgC/m^3)
                  no3, &          ! ice layer brine NO3 concentration (然ol)
                  nh4, &          ! ice layer brine NH4 concentration (然ol)
                  po4, &          ! ice layer brine PO4 concentration (然ol)
                  sioh4           ! ice layer brine SiOH4 concentration (然ol)
          end type pond_type


          type, public :: ice_type
		      integer :: &
		          z               ! number of active layers
		      integer, dimension(z_max) :: &
		          bced, &         ! brine channeled - indicated whether layer has good connection with underlying seawater
                  drained         ! layer is drained of brine, affecting optical properties and thermal conductivity and capactitance
              double precision :: &
                  ts, &           ! surface temperature of top layer (ice OR snow) (degC)
                  fbh, &          ! freeboard height of ice pack (m)
                  pr_sum, &       ! precipitation sum - record total precipitation in cell (kg/m^2)
                  c_gl_bv, &      ! congelation growth layer brine volume - used to determine density new ice growth
                  sh_offset, &    ! offset for calculating thue snow height following flooding (m) - used with boundary_file data
                  sh_prev, &      ! snow height from previous time step (m) - used to calculate incremental snow height change during timestep
                  af, &           ! area fraction - records fraction of grid cell that this ice type represents
                  Ed0, &
                  Ed1, &
                  PAR_bot, &
                  PAR_bot_pl, &
                  ed0_nir_dir, &
                  ed0_nir_dif, &
                  age, &
                  ridged, &
                  snow_dist, &
                  snow_rd

              double precision, dimension(z_max) :: &
                  t, &            ! ice layer temperature (degC)
                  s, &            ! ice layer bulk salinity (psu)
                  d, &            ! ice layer density (g/m^2)
                  bs, &           ! ice layer brine salinity (psu)
                  bd, &           ! ice layer brine density (g/m^2)
                  bv, &           ! ice layer brine volume (ppt)
                  llim, &         ! light limitation fraction (dimensionless)
                  nlim, &         ! nitrogen limitation fraction (dimensionless)
                  plim, &         ! phosphorus limitation fraction (dimensionless)
                  silim, &        ! silica limitation fraction (dimensionless)
                  slim, &         ! salinity limitation fraction (dimensionless)
                  heat, &         ! ice layer weirdo-enthalpy (J/g)
                  th, &           ! ice layer thickness (m)
                  id, &           ! ice layer bottom depth (m)
                  smalg, &        ! brine-based algal concentration (mgC/m^3)
                  prod, &         ! time-step productivity (mgC/grid cell)
                  Ik1, &          ! time-step Ik-prime (湲in/m^2/s)
                  gmax, &         ! time-step maximum growth rate
                  dsdt, &         ! desal rate (psu/timestep)
                  fbv, &          ! brine volume flux m^3/m^2/timestep
                  dhdt_conv, &       ! storage for timestep-derived dhdt from last timestep, since desal occurs in reponse to heat flux
                  f0, &
                  dsdt3, &
                  tgrad, &
                  ed_w
              double precision, dimension(z_max+1) :: &
                  poc, &          ! particulate organic carbon (detritus) (mgC/m^3)
                  no3, &          ! ice layer brine NO3 concentration (然ol) - instantaneous concentration
                  nh4, &          ! ice layer brine NH4 concentration (然ol) - instantaneous concentration
                  po4, &          ! ice layer brine PO4 concentration (然ol) - instantaneous concentration
                  sioh4, &        ! ice layer brine SiOH4 concentration (然ol) - instantaneous concentration
                  no3_mean, &          ! ice layer brine NO3 concentration (然ol) - conc. available to algae, over time
                  nh4_mean, &          ! ice layer brine NH4 concentration (然ol) - conc. available to algae, over time
                  po4_mean, &          ! ice layer brine PO4 concentration (然ol) - conc. available to algae, over time
                  sioh4_mean           ! ice layer brine SiOH4 concentration (然ol) - conc. available to algae, over time
              type (snow_type) :: &
                  snow         ! snow_type - an ice grid cell may have more than one snow_type
              type (pond_type) :: &
                  pond            ! pond_type - melt pond variables
              double precision, dimension(z_max,la) :: &
                  pur             ! Photosynthetically useable radiation by depth and wavelength (湲in/m^s/s/wavelength)
          end type ice_type


          type, public :: platelet_type
              double precision :: &
                  t, &
                  s, &
                  d, &
                  bs, &
                  bd, &
                  bv
              double precision, dimension(pl_max) :: &
                  no3, &
                  nh4, &
                  po4, &
                  sioh4, &
                  Ik1, &
                  smalg, &
                  poc, &
                  pur
          end type platelet_type

          type, public :: snow_adv_type
              integer :: &
                  z
              double precision :: &
                  sh_prev, &
                  depth           ! new advected snow depth, before ridging
              double precision, dimension(z_max) :: &
                  th_new
          end type snow_adv_type

          type, public :: ice_adv_type
              integer :: &
                  z               ! new number of ice layers resulting from advection
              double precision :: &
                  af, &           ! areal fraction of advected ice in this category, not yet ridged, could be > 1!
                  depth, &           ! new advected internal ice depth depth (before ridging)
                  t_mean
              double precision, dimension(z_max) :: &
                  th_new, &       ! new ice layer thicknesses resulting from advection (m)
                  id_new, &       ! new ice layer depths resulting from advection (m)
                  th_corr, &      ! corrected ice layer thicknesses resulting from boundaries during advection (m)
                  th_debug, &     ! used to record particular ice thicknesses during advection for debugging purposes (m)
                  th_snow_new, &  ! new snow layer thicknesses resulting from advection (m)
                  th_old          ! previous ice layer thicknesses (before advection) (m)
          end type ice_adv_type

          type, public :: adv_type
              integer :: &
                  v_mod
              double precision :: &
                  adv_sum, &      ! sum of in and out advection fractions
                  out, &          ! areal percentage of cell advected out
                  out1, &         ! convergence-optimzed areal percentage of cell advected out
                  icecon_corr, &  ! fractional correction of icecon due to next ice concentration forcing (dimensionless)
                  a, &            ! ice-covered area of cell (m^2)
                  af, &            ! ice-covered fraction of cell (m^2)
                  a_new, &        ! ice-covered area of cell after advection(m^2)
                  a_out, &        ! ice-covered area advected out of cell(m^2)
                  a_in, &         ! ice-covered area advected into cell(m^2)
                  p, &            ! ice strength (N)
                  cvm0, &         ! convergence minimization parameter, original
                  cvm1, &         ! convergence minimization parameter, currently optimized
                  v_scale, &      ! velocity scalar, used to change velocity using shared memory
                  mass, &         ! mean ice pack mass (kg/m^3)
                  a_open, &       ! pre-advection area open water
                  af_open, &      ! pre-advection fraction of open water
                  af_new, &       ! area fraction of new ice, used to ad new ice during redistribution of ice categories
				  a_convg, &      ! area of input converging
				  f_convg, &      ! fraction of input converging
                  a_melt, &
                  a_allowed, &
                  id_mean, &
                  a_drop
			  double precision, dimension(8)  :: &
                  in, &           ! convergence-optimzed areal percentage inputs
                  in1, &          ! areal percentage of bordering cell advected in - before advection
                  in1a            ! area of advection in
!			  double precision, dimension(0:8)  :: &
!				  a_convg, &      ! area of input converging
!				  f_convg         ! fraction of input converging
			  double precision, dimension(ic_n)  :: &
			      a_rdg           ! final total area that goes ends up ridging according to rf (km^2)
			  double precision, dimension(ic_n,ic_n)  :: &
                  rf              ! ridging factor - areal reduction scalar for ridged ice
              type (ice_adv_type), dimension(sc_n,ic_n) :: &
                  ice            ! snow_adv_type - holds advection variables for snow layers
              type (snow_adv_type), dimension(sc_n,ic_n) :: &
                  snow            ! snow_adv_type - holds advection variables for snow layers
         end type adv_type

          type, public :: stn_write_type
              logical :: &
                  valid           ! turns on/off writing for station
              integer :: &
                  step            ! index to the current writeout step
              character (LEN=80) :: &
                  fname           ! string holds output filename for station
              double precision :: &
                  time            ! write out model time
              double precision, dimension(18) :: &
                  s               ! array to hold single-value station variables
              double precision, dimension(30,z_max) :: &
                  z               ! array to hold z-dimension station variables
              double precision, dimension(5,z_max+1) :: &
                  z1              ! array to hold z1-dimension station variables
              double precision, dimension(3,wavl) :: &
                  wavl            ! array to hold wavl dimension station variables
              double precision, dimension(3) :: &
                  icevec          ! array to station icevec data
          end type stn_write_type


! ======================================================================
! Global variables & data for referencing in all subroutines
! ======================================================================

          logical, save :: &
              leap_year,start,do_write,z_odd
          integer, save :: &
              z, bins, begin_year, begin_j_date, cur_year, stn_only, &
              end_j_date, end_year, steps, last_day, last_hour, last_year, woa, &
              grid_model, mdh1,mdh2,mdv1,mdv2, boundary_file, ohf,use_pl, &
              z_il, z_fb, z_bt, z_sk, mdh, mdv, ncep_f,last6hour,last12hour,write_step, &
              snowd_index, hour24,tcells,woa_depth,ts_is_at, sda_n, &
              pur_clock, dt_per_day, n_stations, adv_on,skdump,alg_dtt, &
              write_disk,dtt_1,dtt_2,dtt_3,z_int_max,alg_mig,iin, icevec_index,&
              ncep_n, max_it,dt_step,snow_loaded,desal,kevin_light, &
              z_int_min,i_temp,pr_on, snow_model, n_dtt_1,n_dtt_2,n_dtt_3, &
              last_month,atm_f,iit,iis,ida_n,ida_norm_n,lda_n,lda_norm_n, &
              f_index,f_index_next,gc_offset,par_cf,use_mdiff,ml_z, &
              sk_z_adv,int_z_adv,sk_1_adv,use_ponds,monte_carlo,rs_switch, &
              use_drain, dhdt_desal, mix_each, conv_flux,exp_bins,use_expt, &
              restart, snow_ridging,no_flooding,last3hour,use_gpu,flood_brine, &
              atmo,override_ic,wr_stations,icecon_mod,icecon_index,snow_in_gaps, &
              use_sose, sose_level, last_5_day, hemi, init_ice_type
          double precision, save :: &
              dt, Ek, A, B, xp, c_chl, c_n, cur_hour, &
              c_si, c_p, Ek_max, Ks_NH4, Ks_NO3, Ks_SiO4, Ks_PO4, fw, &
              Sw, IceD, rg, G0, fb_h,pur_stepper, &
              il_h, sk_h, Tw, bt_h, ksnow,gl_max, sw_poc, &
              f_limit, f_sk, fb, fi, sw_NO3, sw_NH4, write_f, aph_max, &
              sw_SiO4, sw_PO4, vb_crit, dhdt_i, next_write_hour, &
              min_alg, dt_s, max_snow_h,iif, &
              min_sal, remin_f, pur_divisor,bv_conv, &
              nr_tol,s_const,b_flux_max_tot,z_th_crit,h_crit, &
              remin, dt_years, snow_rd_lim,&
              n_c,p_c,si_c,chl_c,dt_days,cv_void_f, &
              alb_s_dry,alb_s_wet,alb_i_dry,alb_i_wet,den_s_dry, &
              den_s_wet,Sd,eps_snow,eps_ice,alg_mig_crit,n_f, &
              h_min,max_h,gl_max_f,Lf_h,T_ref,temp_tol,d_temp,sk_q, &
              sk_bs,sk_bd,sk_d,sk_bv,sk_s,sk_h_max,sk_heat,sk_h_min, &
              conv_max,vb_max,z_th_min,z_th_max,sk_th_min,sk_th_max, &
              dt_sub_1,dt_sub_2,dt_sub_3,dt_sub_1_tol,dt_sub_2_tol, &
              dtt_s_1,dtt_s_2,dtt_s_3,heat_snow0,bb_f,fl_max,fl_max_f, &
              fl_crit,fl_ratio,cur_month,atan_c_i,h_snowpatch, &
              par_to_swd,ida_d,da_f,lda_d,snow_min,alg_wc, &
              ohf_skew,at_skew,snow_skew,snow_fudge,ad_denom,hbd_f, &
              min_t,des_s_switch,a_ice_ir,a_ice_ir_cos,r_s_dry,r_s_wet, &
              r_i_23,r_i_8,min_scat,a_factor, dbvdt_scale,exp_max,exp_min, &
              exp_mul,exp_a_min,exp_a_max, expls, total_flooded, p_factor, &
              cvf_thin,cvf_switch, init_ice_th, init_snow_th
          double precision, save, dimension (301) :: &
              aice, &
              aph, &
              awater, &
              awater_tc, &
              awater_sc, &
              awater_v, &
              aozone, &
              rdsnow, &
              rwsnow, &
              surfacelight
          character (LEN=14) :: &
              datestr
          integer(2), save, pointer :: &
              SSMI_grid_int2(:,:), &
              icevec_grid_int2(:,:,:), &
              Ed_int2(:,:,:), &
              mp_grid_int2(:,:), &
              ec_grid_int2(:,:), &
              eci_grid_int2(:,:)
          double precision, save, pointer :: &
              kds_wet(:), &
              kds_dry(:), &
              lambda(:), &
              quanta_to_watts(:)
					type (sose_f_type), save, pointer :: &
							sose_f
          double precision, save, dimension (130) :: &
              mc_prod
          double precision, save, dimension (ic_n) :: &
              ic_h_max, &
              ic_h_med, &
              ic_h_min

     end module sia2_globals

