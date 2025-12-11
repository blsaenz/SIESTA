! file: seaice2.f90
! Sea Ice Algae Model 2 - Saenz & Arrigo
! Version beta
! ======================================================================


! openmp environment variables
! ----------------------------------------------------------------------
! OMP_SCHEDULE 	Sets the runtime schedule type and chunk size. 	Default: STATIC ( no chunk size specified)
! OMP_NUM_THREADS 	Sets the number of threads to use during execution. Default: Number of processors
! OMP_DYNAMIC 	Enables (TRUE) or disables (FALSE) the dynamic adjustment of the number of threads. Default: FALSE
! OMP_NESTED 	Enables (TRUE) or disables (FALSE) nested parallelism.  Default: FALSE

      program seaice2

      use sia2_globals
#ifdef GPU
      use sia2_gpu
#endif
#ifdef IFORT
      use IFCORE
#endif
!      use OMP_LIB
      implicit none


! Declarations
! ======================================================================

      ! PROGRAM, THREAD, and MEMORY MANAGEMENT VARIABLES
      logical :: eoy,exists
      integer :: i,ii,j,jj,start_step,fpe_old,fpe_new
      integer :: now(3)
      double precision :: h_layer,r_depth  ! temp var used in some max/min layer calcs

      ! PROJECTION, GRIDDING, and STATION VARIABLES
      type (proj_type), pointer :: pj
      type (meta_type), allocatable :: m(:)
      integer, allocatable :: stations_mdh(:), stations_mdv(:)

      ! STATIC FORCING DATA VARIABLES
      type (mp_f_type), pointer :: mp_f
      type (dg_f_type), pointer :: dg_f
      type (ease_f_type), pointer :: ease_f
      type (ec_f_type), pointer :: ec_f
      type (eci_f_type), pointer :: eci_f
      type (boundary_f_type), pointer :: v_f

      ! MODELED GRID VARIABLES
      type (forcing_type), allocatable :: f(:)
      type (ice_type), allocatable :: ice(:,:,:)
      type (ice_type), allocatable :: slush(:,:)
      type (adv_type), allocatable :: adv(:)
      type (platelet_type) :: pl
	  integer (kind=2), allocatable :: PUR(:,:,:,:,:,:)

      ! PRINTING/OUTPUT VARIABLES
      double precision, allocatable :: sum_var(:,:,:),temp3d(:,:)
      type (stn_write_type), allocatable :: stn_output(:,:,:)
      character (len=80) :: out_fname
      character (len=80) :: stat_string
      character (len=6) :: station_format
      character (len=2) :: stat_hr1
      character (len=3) :: stat_day
      character (len=4) :: stat_hr2,stat_year
      character (len=6) :: stat_step
      character (len=2) :: cn_stations

      ! LOOKUP TABLE & MATH VARIABLES
      double precision, allocatable :: expt(:,:)
      double precision, pointer :: gauss_bins(:)

! Read Parameters
! ======================================================================
	  ! find date & time string for output filenames
	  call sia2_do_datestr(datestr)

      call sia2_read_parameters()

! Model Run Initializations - getting everything ready for step 1
! ======================================================================

      ! determine seawater temperature from salinity,
      ! adding a little bit to avoid numerical with freezing (i.e. should be frozen given temp...)
      if(Tw .lt. (Sw*mu + 0.0001)) then
          Tw = Sw*mu + 0.0001
      endif
      min_t = 4. !min_sal*mu*1.01

      ! figure out skeletal paramaters, used to reduce future calculations
      sk_s = Sw*0.5   ! psu
	  sk_bs = Tw/mu ! psu
	  sk_bd = 1000000. + sk_bs*800.  ! g/m^3
	  sk_d = (1.-bb_f)*IceD*sk_bd*sk_bs/(sk_bd*sk_bs - sk_s*(sk_bd - IceD)) ! g/m^3
	  sk_bv = 1000.*(sk_d*sk_s)/(sk_bd*sk_bs)  ! ppt
	  sk_heat = -1.*sk_d*(c0*(sk_s*mu - Tw) &
		+ Lf*(1. - sk_s*mu/Tw) - cw*sk_s*mu) ! j/m^3 - not j/m^2
      sk_q = sk_d*(c0*(sk_s*mu - Tw) &
		+ Lf*(1. - sk_s*mu/Tw))

	  sk_h_max = DBLE(z_sk)*sk_th_max
      sk_h_min = DBLE(z_sk)*sk_th_min! (m) - minimum total skeletal thickness

      heat_snow0 = 0.2309*kelvin0 + 0.0034*kelvin0**2

      ! determine min, max, critical, growth layer heighta for ice
      ! gridding routines
      ! --------------------------------------------------------------
      z_int_max = z_max - z_sk
      ml_z = z_max+1 ! index to mixed layer values - always at the end of ice array
	  ! skeletal layer references to bottom of of arrays,used in advection routines
	  sk_z_adv = z_max
	  int_z_adv = z_max-z_sk
	  sk_1_adv = int_z_adv+1
      ! minimum congelation ice height
      !if (grid_model .lt. 2) then
      !    z_int_min = z_int_max ! non-changeable number of layers
      !endif
      !z_int_min = 10
	  h_min = z_th_min*z_int_min ! (m) skeletal layers don't count for height calcs
      h_crit = z_th_min*z_int_max ! ice height at which the expanding layer gridmodel is used (grid_model = 2 only)

      ! find h_max (max ice height) by accordian-expanding inner layers to full 1.3x between-layer height

!      h_max = 2.*z_th_min
!      h_layer = z_th_min
!	  do ii=2,int(z_int_max/2)
!	      h_layer = h_layer*1.3
!	      h_max = h_max + h_layer*2.
!	  enddo


	  gl_max = z_th_min*gl_max_f ! (m) largest step that boundary adjustment may take
	  fl_max = z_th_min*fl_max_f ! (m) largest step that boundary adjustment may take for flooding
! ----------------------
! for testing
!	  fl_max = 0.001 ! (m) largest step that boundary adjustment may take for flooding
! ----------------------

      total_flooded = d0_

      atan_c_i = atan_max/c_5   ! 0.5 m is the cutoff for the atan function that describes albedo

      ! find inverse nutrient ratios, timing vars, and other often things
      ! for the purpose of calculation efficiency
      ! ---------------------------------------------------------------
      n_c = 1./c_n
      p_c = 1./c_p
      si_c = 1./c_si
      chl_c = 1./c_chl

	  if (la .eq. 1) then
		  lda_n =1
		  allocate(lda_multiplier(lda_n))
		  lda_multiplier(1) = 1.
		  lda_d = 1.
	  elseif (la .eq. 9) then
		  lda_n = logn_n
		  allocate(lda_multiplier(lda_n))
		  lda_multiplier = logn_multiplier9
		  lda_norm_n = 6
		  lda_d = 1./dble(logn_n)
	  elseif (la .eq. 6) then
		  lda_n = 6
		  allocate(lda_multiplier(lda_n))
		  lda_multiplier = logn_multiplier6
		  lda_norm_n = 4
		  lda_d = 1./dble(logn_n)
	  elseif (la .eq. 5) then
		  lda_n = 5
		  allocate(lda_multiplier(lda_n))
		  lda_multiplier = logn_multiplier5
		  lda_norm_n = 3
		  lda_d = 1./dble(logn_n)
	  else
          print *,'Light category resolution (la) must be 1, 5, 6, or 9'
          call exit(0)
      endif


      ! find mean cosine near-IR absorption in ice
      a_ice_ir_cos = a_ice_ir/0.656
      ! pre-calculate snow spectral attenuation curves
      allocate(kds_wet(wavl))
      allocate(kds_dry(wavl))
      allocate(lambda(wavl))
      allocate(quanta_to_watts(wavl))
	  ad_denom = d0_
	  do ii=1,wavl
	      ! record wavelength of spectral bin
	      lambda(ii) = dble(ii)*10.+390.  ! nm
	      quanta_to_watts(ii) = 7.e-7*lambda(ii)**2 - 0.0012*lambda(ii) + 0.6416
          ! find denominator for ad (detrital absorption) estimation from chlorophyll
		  ad_denom = ad_denom + exp(-0.008*dble((ii-1)*10))*10. ! times 10, b/c 10nm wavelength bins
		  Kds_dry(ii) = 104. - 0.38*lambda(ii) + 0.0004248*lambda(ii)**2    ! grenfell & maykut 1977
		  Kds_wet(ii) = 72.36 - 0.27*lambda(ii) + 0.0002883*lambda(ii)**2  ! grenfell & maykut 1977
	  enddo

      ! determine if using effective snow catagories for better flooding estimation
      ! ---------------------------------------------------------------
      if (sc_n .gt. 1) then
          sda_n = 2
          if (sd_num .eq. 9) then
              allocate(sda_multiplier(sd_num))
              sda_multiplier = logn_multiplier9
          elseif(sd_num .eq. 6) then
              allocate(sda_multiplier(sd_num))
              sda_multiplier = logn_multiplier9
          elseif(sd_num .eq. 5) then
              allocate(sda_multiplier(sd_num))
              sda_multiplier = logn_multiplier9
          else
              print *,'Snow category distribution resolution (sd_num) out of range: must be 5 or 6 or 9'
              call exit(0)
          endif
      else
          if (sc_n .ne. 1) then
              print *,'Snow category distribution resolution (sd_num) must be 1 for sd_c==1'
              call exit(0)
          endif
          sda_n = 1
          allocate(sda_multiplier(1))
          sda_multiplier(1) = 1.d0
      endif

      ! determine if using an ice distribution model (ida) and set globals
      ! and ice category boundaries
      ! ---------------------------------------------------------------
      if (ic_n .gt. 0) then

		  ic_h_min(1) = h_min

          if (ic_n .eq. 1) then
              ida_n = ic_n
			  ic_h_max(1) = h_max
          elseif (ic_n .ge. 3 .and. ic_n .le. 10) then
              ida_n = ic_n
			  select case (ida_n)
				  case (3)
                      gauss_bins => gauss_bins_3
				  case (4)
                      gauss_bins => gauss_bins_4
				  case (5)
                      gauss_bins => gauss_bins_5
				  case (6)
                      gauss_bins => gauss_bins_6
				  case (7)
                      gauss_bins => gauss_bins_7
				  case (8)
                      gauss_bins => gauss_bins_8
				  case (9)
                      gauss_bins => gauss_bins_9
				  case (10)
                      gauss_bins => gauss_bins_10
              end select
              do jj=2,ida_n
				  ic_h_max(jj-1) = h_min + (h_max-h_min)*gauss_bins(jj)
				  ic_h_min(jj) = ic_h_max(jj-1)
			  enddo
			  ic_h_max(ida_n) = h_max

              ! subverting ida=5 with custom categories
              if (ida_n .eq. 5) then
                  i=1
                  j=i+1
                  ic_h_max(i) = h_min + (h_max-h_min)*0.02040816326531
                  ic_h_min(j) = ic_h_max(i)
                  i=j
                  j=i+1
                  ic_h_max(i) = h_min + (h_max-h_min)*0.05102040816327
                  ic_h_min(j) = ic_h_max(i)
                  i=j
                  j=i+1
                  ic_h_max(i) = h_min + (h_max-h_min)*0.10204081632653
                  ic_h_min(j) = ic_h_max(i)
                  i=j
                  j=i+1
                  ic_h_max(i) = h_min + (h_max-h_min)*0.28571428571429
                  ic_h_min(j) = ic_h_max(i)
!                  i=j
!                  j=i+1
!                  ic_h_max(i) = h_min + (h_max-h_min)*0.41
!                  ic_h_min(j) = ic_h_max(i)
              endif

              do jj=1,ida_n
                  ic_h_med(jj) =  (ic_h_min(jj) + ic_h_max(jj))/2.
              enddo

          else
              print *,'Ice category parameter out of range: must be 1 or 3 to 9'
              call exit(0)
          endif
      else
          print *,'Ice category parameter out of range: must be 1 or 3 to 9'
          call exit(0)
      endif

      allocate(v_f)


          if (boundary_file .eq. 0) then

          ! can't use platelet layer unless this is a boundary_file run
          use_pl = 0
          iif = d0_

          monte_carlo = 0

		  ! read in locations of stations where detailed records are written out
		  ! ---------------------------------------------------------------
		  open(unit=21, file='stations.txt', form='FORMATTED', access='sequential')
		  read(21,*) n_stations  ! number of stations to read in from stations.txt
		  allocate(stations_mdh(n_stations))
		  allocate(stations_mdv(n_stations))
		  allocate(stn_output(sda_n,ida_n,n_stations))
		  write(cn_stations,'(i2)') n_stations
		  !station_format = '('//cn_stations//'i5)'
      station_format = '(1i5)'
		  !write(station_format,'(i1)') n_stations
		  !station_format = '('
		  !write(station_format(2),'(i1)') n_stations
		  !station_format(3:5) = 'i5)'
		  read(21,FMT=station_format) stations_mdh
		  read(21,FMT=station_format) stations_mdv
		  close(21)

		  ! figure out how large model domain is for array memory allocation
		  ! ---------------------------------------------------------------
		  mdh = mdh2-mdh1+1
		  mdv = mdv2-mdv1+1

      else
		  ! setup for boundary_file run
		  ! ---------------------------------------------------------------
          mdh = 1
          mdv = 1
          ! set these to some random, valid cell so that some functions don't barf
          mdh1 = 281
          mdh2 = 281
          mdv1 = 306
          mdv2 = 306

          ! figure out where boundary_file location is, so model domain is setup correctly
          call sia2_forcing_load_v_f(v_f,-1,-1)
          if (monte_carlo .eq. 1) then
              n_stations = 40
          else
              n_stations = 1
          endif
		  allocate(stations_mdh(n_stations))
		  allocate(stations_mdv(n_stations))
		  allocate(stn_output(sda_n,ida_n,n_stations))

      endif



      ! figure out how many model steps per day (dt_per_day) - needed compute 24-hour
      ! PUR mean
      ! ---------------------------------------------------------------
      dt_per_day = int(24./dt)

      if (dt_per_day .lt. 1) then
          dt_per_day = 1
      endif

      pur_divisor = 1./dble(dt_per_day)
      pur_stepper = 1. - pur_divisor

      print *,'Allocating memory...'

      ! Allocate Projection & Forcing Data Arrays
      ! ---------------------------------------------------------------
      allocate(pj)      ! allocate model projection info type
      allocate(ease_f)  ! allocate EASE projection forcing type
      allocate(mp_f)    ! allocate mercator projection (NCEP) forcing type
      allocate(dg_f)    ! allocate 1 degree grid projection forcing type
      allocate(ec_f)    ! allocate ECMWF projection forcing type
      allocate(eci_f)    ! allocate ECMWF projection forcing type

      ! read in projection and grid data
      ! ---------------------------------------------------------------
      ! read in mask, lat, lon, files and determine model domain
      call sia2_read_pj(pj,v_f,stations_mdh,stations_mdv)
      if (monte_carlo .eq. 1) then
          tcells = n_stations
      endif

      ! Allocate temporary data arrays
      ! ---------------------------------------------------------------
      allocate(SSMI_grid_int2(grid_h,grid_v))
      if (hemi .eq. 0) then
        allocate(icevec_grid_int2(grids_h,grids_v,3))
      else
        allocate(icevec_grid_int2(gridn_h,gridn_v,3))
      endif
      allocate(Ed_int2(mp_x,mp_y,wavl))
      allocate(mp_grid_int2(mp_x,mp_y))
      allocate(ec_grid_int2(ec_x,ec_y))
      allocate(eci_grid_int2(eci_x,eci_y))
      allocate(temp3d(tcells,z_max))

      ! Create Calculated Variable Grids
      ! ---------------------------------------------------------------

      ! MODELED GRID VARIABLES
      allocate(f(tcells))
      allocate(m(tcells))
      allocate(ice(sda_n,ida_n,tcells))
      allocate(pur(z_max,dt_per_day,lda_n,sda_n,ida_n,tcells))
      ! ADVECTION VARIABLES
      allocate(adv(tcells))

      ! PRINTING/OUTPUT VARIABLES
      allocate(sum_var(mdh,mdv,105))

      ! Allocate SOSE field arrays
      ! ---------------------------------------------------------------
      if (use_sose .eq. 1) then
				allocate(sose_f)
			endif


#ifdef GPU

      ! Allocate GPU arrays
      ! ---------------------------------------------------------------
!      if (use_gpu .eq. 1 .and. tcells .ge. 40) then
      if (use_gpu .eq. 1) then
          ! allocate host arrays, and figure out gpu max threads & blocks
          call gpu_allocate(z_max+z_snow_max)

          ! initialize GPU arrays to zero
          g_i32_in = 0.
          g_f32_in = 0.
          g_f32_scr = 0.
          g_f32_out = 0.

          ! assign gpu calc constants
           h_inv_c_chl = real(1.d0/c_chl)       ! host inverse carbon/chla ratio
		   h_inv_iced = real(1.d0/IceD)         ! host inverse pure ice densiry
		   h_a_ice_ir = real(a_ice_ir)          ! host near-ir abosroption
		   h_inv_ad_denom = real(1.d0/ad_denom) ! host inverse ad_denom
		   h_a_factor = real(a_factor)          ! host snow extra absorption factor

      else
          use_gpu = 0
      endif

#endif


      ! Allocate lookup tables
      ! ---------------------------------------------------------------
	  allocate(expt(2,exp_bins))
	  expls = exp_a_min/dble(exp_bins)
	  exp_a_max = d0_
	  exp_max = 1.d0
	  exp_min = exp(exp_a_min)
	  exp_mul = dble(exp_bins)/exp_a_min

      print *,'Done allocating memory.'

      ! initialize memory
      ! ---------------------------------------------------------------
      call sia2_initialize(pj,m,f,ice,pur,adv,mp_f,ease_f,expt)

      ! Initialize timing variables
      !-----------------------------------------------------------------

      if (boundary_file .eq. 1) then
          cur_hour = v_f%time
          begin_year = v_f%begin_year
          begin_j_date = int((cur_hour)/24)+1
          !steps = int(v_f%length/dt)-1
      else
          cur_hour = begin_j_date*24-24
      endif
		  call sia2_num_steps(begin_year,begin_j_date,end_year,end_j_date,dt,steps)

      f_index = 1
      write_step = 1
      do_write = .false.
      pur_clock = 1
      dt_s = dt*3600.  ! hrs->seconds
      dt_days = dt/24.d0 ! days
      dt_years = dt/8760.d0 ! years
      dtt_1 = int(dt/dt_sub_1) ! number of sub-steps in ice growth
      dtt_2 = int(dt/dt_sub_2) ! number of sub-steps in ice growth
      dtt_3 = int(dt/dt_sub_3) ! number of sub-steps in ice growth
      dtt_s_1 = dt_sub_1*dt*3600 ! secs
      dtt_s_2 = dt_sub_2*dt*3600 ! secs
      dtt_s_3 = dt_sub_3*dt*3600 ! secs

      mc_prod = 0.


      ! check for restart files - is so, reload saved memory and update vars
      ! -----------------------------------------------------------------
      if (restart .eq. 1) then

				call sia2_restart(m,ice,adv,pur,out_fname,stn_output,start_step,1)
				call sia2_isleapyear(cur_year,leap_year)
				hour24 = cur_hour-int(cur_hour/24)*24
						start = .false.

      else

				cur_year = begin_year
				call sia2_isleapyear(cur_year,leap_year)
				hour24 = cur_hour-int(cur_hour/24)*24
				last_day = 0
				last_year = cur_year-1
				last3hour = 0
				last6hour = 0
				last_5_day = 0
				last12hour = 0
				last_hour = 0
				last_month = 0
				start_step = 0
				start = .true.

				! construct output filenames with run date and time
				call sia2_write_construct_fname(out_fname,'output              ',6,'nc ',2,datestr)
				print *, 'Data Output File: ',trim(out_fname)

      endif

! Model Loop
! ======================================================================

      do i=start_step,steps

	  ! write out restart file 8x per year
	  !-----------------------------------------------------------------
!          if ((cur_hour .eq. 0 .or. cur_hour .eq. 1080 .or. cur_hour .eq. 2160 .or. &
!              cur_hour .eq. 3240 .or. cur_hour .eq. 4320 .or. cur_hour .eq. 5400 .or. &
!              cur_hour .eq. 6480 .or. cur_hour .eq. 7560) .and. (restart .ne. 1)) then
!              call sia2_restart(m,ice,adv,pur,out_fname,stn_output,i,0)
!          endif

		  ! write out restart file every 5 days!
		  !-----------------------------------------------------------------
          if ((mod(cur_hour,120.) .eq. 0.) .and. &
              (restart .ne. 1) .and. &
              (.not. start)) then
!              call sia2_restart(m,ice,adv,pur,out_fname,stn_output,i,0)
          endif

		  ! write out restart file every day!
		  !-----------------------------------------------------------------
          !if (((int(cur_hour/24.)+1) .gt. last_day) .and. &
          !    (restart .ne. 1) .and. &
          !    (.not. start)) then
          !   call sia2_restart(m,ice,adv,pur,out_fname,stn_output,i,0)
          !endif

          dt_step=i
          f_index_next = f_index
          print *, '========================================================================='
          write(stat_year,'(i4)')cur_year
          write(stat_day,'(i3)')int(cur_hour/24)+1
          write(stat_hr1,'(i2)')hour24
          write(stat_hr2,'(i4)')nint(cur_hour)
          write(stat_step,'(i6)')i
          stat_string = 'Year: '//stat_year//'   Day: '//stat_day//'   Hour_24: '//stat_hr1
          stat_string = TRIM(stat_string)//'   Hrs since Jan1: '//stat_hr2//'   Step: '//stat_step
          print *,stat_string
          !print *,  'Year:',cur_year,'Day:', int(cur_hour/24)+1,'Hour_24:',hour24, &
          !      'Hrs since Jan1:',cur_hour, 'Step:',i
          print *, '-------------------------------------------------------------------------'

      ! run model steps
      !-----------------------------------------------------------------

          ! get data required for current timestep
          ! --------------------------------------------------------------
          print *,'Loading forcing data ...'
          call sia2_forcing_load(pj,mp_f,dg_f,ec_f,eci_f,ease_f,v_f)

          ! interpolate ncep and ssmi-snow forcing data for current timestep
          ! --------------------------------------------------------------
          print *,'Interpolating forcing data over timestep...'
          call sia2_interpolate_forcing(pj,m,f,ease_f,mp_f,dg_f,ec_f,eci_f, &
              v_f,stations_mdh,stations_mdv)

#ifdef GPU
          if (use_gpu .eq. 1) then
              if (i .eq. start_step) then
				  fpe_new = 0
				  fpe_old = FOR_SET_FPE(fpe_new)
				  print *,'Set FPE Trap from: ',fpe_old,' to: ', fpe_new
                  call sia2_edd_gpu_init_malloc(g_i32_in,g_f32_in,g_f32_scr, &
                    g_f32_out,i32_in_size,f32_in_size,f32_scr_size,f32_out_size, &
                    h_inv_c_chl,h_inv_iced,h_a_ice_ir,h_inv_ad_denom,h_a_factor)
                  if (fpe_old .ne. 0) then
					  fpe_new = 983054
					  fpe_old = FOR_SET_FPE(fpe_new)
    				  print *,'Set FPE Trap from: ',fpe_old,' to: ', fpe_new
                  endif
              endif
			  call itime(now)
			  print '(a25,a40,i3,a1,i2,a1,i2)',' Starting GPU Light Calc:',' ',now(1),':',now(2),':',now(3)
              call sia2_edd_gpu_calc(mp_f,f,m,ice)
          endif
#endif

          ! Environment & Growth Calc
          ! --------------------------------------------------------------
          call itime(now)
          print '(a27,a38,i3,a1,i2,a1,i2)',' Starting Environment Calc:',' ',now(1),':',now(2),':',now(3)
          ! this ugly function 1) creates/grows ice, 2) calculates the ice environment,
          ! and 3) calculates PUR and grows algae.  It is ugly and mostly inline, to
          ! allow fast parallel exectution.  Sorry in advance to anyone who has to read it.
          call sia2_environment_calc(pj,m,ice,pur,f,v_f,mp_f,ease_f,dg_f,adv, &
              sum_var,stn_output,out_fname,pl,expt)

          ! Advection
          ! --------------------------------------------------------------
!          if ((adv_on .eq. 1) .and. (tcells .gt. 1) .and. (stn_only .eq. 0)) then
          if ((adv_on .eq. 1) .and. (boundary_file .ne. 1) .and. (stn_only .eq. 0)) then
              call itime(now)
              print '(a39,a26,i3,a1,i2,a1,i2)',' End Environment Calc, Start Advection:',' ', &
              now(1),':',now(2),':',now(3)

 !             call sia2_advect_2(pj,ease_f,f,m,ice,adv)
              call sia2_advect_setup(pj,ease_f,f,m,ice,adv)
              call sia2_advect_tracers(f,m,ice,adv)
              call sia2_adv_grid_ridge(pj,ease_f,f,m,ice,adv)

              call itime(now)
              print '(a31,a34,i3,a1,i2,a1,i2)',' End Advection, Starting Write:',' ', &
              now(1),':',now(2),':',now(3)
          else
              call itime(now)
              print '(a53,a12,i3,a1,i2,a1,i2)',' End Environment Calc (No Advection), Starting Write:', &
              ' ',now(1),':',now(2),':',now(3)
          endif

          ! Write Out
          ! --------------------------------------------------------------
          call sia2_write_output(pj,m,ice,pur,adv,f,mp_f,ease_f,sum_var, &
          out_fname,stations_mdh,stations_mdv,stn_output,i+1)

          call itime(now)
          print '(a11,a54,i3,a1,i2,a1,i2)',' End Write:',' ',now(1),':',now(2),':',now(3)

      ! Increment timing variables
      !-----------------------------------------------------------------
          f_index = f_index_next
          last_day = int(cur_hour/24)+1
          last3hour = int(cur_hour/3)+1
          last6hour = int(cur_hour/6)+1
          last12hour = int(cur_hour/12)+1
          last_5_day = int(cur_hour/120)+1
          last_year = cur_year
          last_hour = int(cur_hour)+1
          call sia2_getmonth(cur_hour,leap_year,last_month)

          cur_hour = cur_hour + dt;
          call sia2_end_of_year(cur_year,cur_hour,eoy,leap_year)
          if(eoy) then
              next_write_hour = next_write_hour - cur_hour + dt
              call sia2_increment_year(cur_year,cur_hour,leap_year)
              last_day = 0
              last3hour = 0
              last6hour = 0
              last12hour= 0
              last_5_day = 0
              last_year = cur_year-1
              last_hour = 0
          endif

          hour24 = cur_hour-int(cur_hour/24)*24
          ! increment PUR clock, so that we're remebering 24-hour PUR in a different
          ! spot in th array - pur_clock just moves around for 1-24
          if (pur_clock .eq. dt_per_day) then
              pur_clock = 1
          else
              pur_clock = pur_clock + 1;
          endif

          if (start) then
              start = .false.
          endif

          if (restart .eq. 1) then
              restart = 0
          endif

      enddo

! Deallocate Grids
! ======================================================================

      deallocate(pj,ease_f,mp_f,ec_f,eci_f,f,m,ice,pur,adv,sum_var,temp3d, &
        stations_mdh,stations_mdv,v_f,SSMI_grid_int2,icevec_grid_int2, &
        Ed_int2,mp_grid_int2,lda_multiplier,stn_output, &
        kds_wet,kds_dry,lambda,expt)

      end program seaice2

