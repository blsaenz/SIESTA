! file: SIA2_environment.F90
! Sea Ice Algae Model 2 - Saenz & Arrigo
! Version beta
! ======================================================================

! SUBROUTINE: sia2_environment_calc
! ======================================================================


      SUBROUTINE sia2_environment_calc(pj,m,ice,pur,f,v_f,mp_f,ease_f,dg_f,adv, &
          sum_var,st_out,out_fname,pl,expt)

          use sia2_globals
!          use OMP_LIB
          implicit none

          ! function arguments
          !----------------------------------------------------
          type (proj_type) :: pj
          type (meta_type) :: m(tcells)
          type (ice_type) :: ice(sda_n,ida_n,tcells)
          integer (kind=2) :: pur(z_max,dt_per_day,lda_n,sda_n,ida_n,tcells)
          type (forcing_type) :: f(tcells)
          type (boundary_f_type) :: v_f
          type (mp_f_type) :: mp_f
          type (ease_f_type) :: ease_f
          type (dg_f_type) :: dg_f
          type (adv_type) :: adv(tcells)
          double precision, dimension(mdh,mdv,105) :: sum_var
          type (stn_write_type) :: st_out(sda_n,ida_n,n_stations)
          character (LEN=60) :: out_fname
          type (platelet_type) :: pl
          double precision :: expt(2,exp_bins)

          ! shared internal variables
          !----------------------------------------------------
          logical :: load_ed,find_lim
          integer :: max_d_exceeded,chunk_size
          double precision :: j_date,rad_date,solar_dec,rad_time ! constant per timestep over earth

          integer         nf_open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             ncid)
          external        nf_open
          integer         nf_close

!                          integer             ncid)
          external        nf_close

          ! private internal variables
          !----------------------------------------------------

          ! general vars
          integer :: mi,i,j,ii,iii,jj,jjj,kk,ic,sc,i_ease,j_ease,i_mp,j_mp,thread,slot
          integer :: i_dg,j_dg,now(3),sk_1,sk_z,int_z,ic_out
          double precision :: airtemp1c,snowfall,snowfall_d,t_depth,snowh1m,tmp1,tmp2,tmp3,tmp4

          ! sia2_env_brine_calc
          double precision :: t_mean,s_mean,d_mean,bs_mean,bd_mean,bv_mean,heat_mean
          double precision :: alpha0,alpha1,alpha2,alpha3

          ! sia2_env_temp_from_heat
          ! double precision :: d_mean,heat_mean,t_mean,s_mean
          double precision ::  dz_mean,aq,bq,cq

          ! sia2_env_layer_thickness
          logical :: keep_growing,maxed_depth,eoy,valid_snowh,ignore_h_max
          integer :: int_z_new,z_last,z_snow_new
          double precision ::  r_depth,target_depth,layer_divisor,sh_diff_mean
          double precision, dimension(sda_n,ida_n) :: sh_increment,sh_diff
          double precision, dimension(z_max) :: th_new,id_new

          ! sia2_env_adj_boundaries
          ! integer :: ii,jj
          double precision :: interim_top,interim_bot,old_layer_top,new_layer_bot
          double precision :: old_layer_bot,z1,z2,dz,heat_total,dz_total,new_layer_top
          double precision :: flooded,t_flood,s_flood,d_flood,bs_flood,bd_flood
          double precision :: bv_flood,heat_flood,t_interface,d_total,af_total,vf_total
          double precision :: s_pond_new,h_pond_new,d_pond_new,smalg_pond_new
          double precision :: poc_pond_new,no3_pond_new,nh4_pond_new,po4_pond_new
          double precision :: sioh4_pond_new,th_pond_new,af_tot,sh_mean
          double precision :: s_flood_new,heat_flood_new,d_flood_new,no3_flood_new
          double precision :: nh4_flood_new,po4_flood_new,sioh4_flood_new,poc_flood_new
          double precision, dimension(z_max) :: t_new,s_new,d_new,bs_new,bd_new,t_change
          double precision, dimension(z_max) :: bv_new,smalg_new,poc_new,smalg_pre_map
          double precision, dimension(z_max) :: melt_new
          double precision, dimension(z_max+1) :: NO3_new,NH4_new,PO4_new,SiOH4_new

          ! Heat Flux & Boundary Adjustment section
          !integer :: ii,jjj
          ! double precision :: d_mean
          logical :: integr_valid,did_flood,alg_fixed,remapping,flood_1x,did_drain
          integer :: info,m_order,z_snow,ice_1,ice_z,dtt,mo,heat_stop,snow_model_temp
          double precision :: sk_gl,c_gl,c_gl_sal,s_gl,sk_cur,smalg_lost,fc_int
          double precision :: p_h2o_sat,row_air,rhum,shum,epsa,v10,Fe_sub1,Fe_sub2,Fs_sub
          double precision :: FL,albedo,eps0,ks,top_a,top_b,Fr,F0_constants,bot_a,bot_b
          double precision :: F0_4,F0_3,F0_2,F0_1,dtt_s_snowh,ki_up,ki_down,Ts_kelvin
          double precision :: Ts_last,Tsi_last,Ts_next,Tsi_next,T_diff,c_snow,F0,dF0,Fc
          double precision :: mat_sd,mat_sb,mat_sc,mat_0a,mat_0b,mat_0c,rs,r_0,xl,xls
          double precision :: th1,th2,thz,th2z,Fc_top,Fc_bot,albedo_s,albedo_i,f_snow
          double precision :: mat_a,mat_b,mat_c,r1,dhdt_g_m2,dhdt,c_gl_bv,dhdt_cm_s
          double precision :: keff,dz_c,dz_sk,sk_th_new,f_depth,ts_prev,th_crunch
          double precision :: c_gl_tot,s_gl_tot,dtt_s,T_step,snow_gl,crunch_tot,fl_gl
          double precision :: t_melt,q_melt,dtt_days,snow_melt,snow_melt_now,melt_flood
          double precision :: mf_depth,melt_drain,ki1,heat_brine,c_gl_copy
          double precision :: T_up, T_dn, T_ii,id_mean, potT,shcoef,lhcoef,p_h2o
          double precision :: Fm0,F_little_L,Fe,F_s ! debug vars
          double precision, dimension(z_max+1) :: dNO3,dNH4,dPO4,dSiOH4,dpoc,dsmalg,dth_ed_w
          double precision, dimension(z_max+z_snow_max) :: dth,t_next,t_last,th
          double precision, dimension(z_max+z_snow_max) :: DU,DU_calc,DL,DL_calc
          double precision, dimension(z_max+z_snow_max) :: DC,DC_calc,r_vec,ki,ci,ki_alt
          double precision, dimension(z_max) :: heat_temp,dcheat,debug_z,debug_start
          double precision, dimension(z_max) :: t_snow_temp,d_snow_new,th_snow_new,heat_snow_new
          double precision :: ts_snow,Fc_snow


          ! Desalination Section
          integer :: vb_open,cv_open,bfrf
          double precision :: T_grad,bs_dilute,dsdt1,dsdt2,bulk_nc,dbvdt, &
              dsdt3,dsdt4,t_brine,d_brine,dsdt_net,fresh_net,dsdt_tmp,fresh_tmp
      double precision :: bs_dif,superimposed
          !double precision, dimension(z_max) :: s_new
          double precision, dimension(z_max+1) :: dsdt,dheat,dheat_last,dsdt_cw
          double precision, dimension(z_max) :: bs_last,bd_last,bv_last,t_prev,t_snow_prev
          double precision, dimension(z_max) :: heat_last,bs_diff,heat_diff

          ! Brine Flux Section
          !integer :: ii,m_order,info
          integer :: flux_steps,conv_start
          double precision :: flux_down,th_down,bv_down,Fd,Fu,Fdiff,Feddy,fvs
          double precision :: dt_flux,flux_ml,flux_d,hv,flux_sk,fvb,fs
          double precision :: b_flux_max,min_perm,bv_flux,bv_total,visc
          double precision :: no3_total,nh4_total,po4_total,sioh4_total,poc_total
          !double precision, dimension(z_max+1) :: dNO3,dNH4,dPO4,dSiOH4,DC,DC_calc
          !double precision, dimension(z_max) :: DL,DL_calc,DU,DU_calc,r_vec
          double precision, dimension(z_max+1) :: b_flux,porosity,perm,pres,p_flux,h_flux


          ! Surface Infiltration Section
          double precision :: fbh_new,top_f,bot_f,ice_mass,snow_mass,snow_dist_new,snow_rd_new

          ! vars required by sia2_env_light.inc.f90
          integer :: nm_i,m_row,a_i,b_i,c_i,d_i,srftyp
!          integer, dimension(8*z_max) :: ipiv
          real :: bv_mean_real
          real, dimension(z_max,wavl) :: ex
          double precision :: zenith_angle,E_sum,lda_f,par_cf_f,kdice_ir,Ed_W_IR
          double precision :: rsun0,rsky0,temp_brine_v,am_sum,Kds,testvar
          double precision :: Ed0_new,Ed0_temp,Ed1_new,Ed1_temp,Ed_temp,kds_z
          double precision :: mean_spec,absorp,ed_med,albodr,albodf,coszen
!          double precision, dimension(wavl) :: Ed_snow1_5,Ed_snow4,Ed_snow10
          double precision, dimension(z_max,wavl) :: Kdice,Kdp,am,ad !,ex
          double precision, dimension(z_max+1,wavl) :: Ed
          double precision, dimension(wavl) :: Edir_last_valid,Edif_last_valid
!          double precision, dimension(8*z_max) :: ABCD
!          double precision, dimension(10,8*z_max) :: AB_mat
          double precision, dimension(2*z_max) :: k2st,r2st,ex2st,g2st,w2st,tau2st, &
             trndir, trntdr, trndif, rupdir, rupdif, rdndif
          ! irradiance summary calculation vars
          double precision :: &
              refk,par_top,par_bot,par_mid,heat_wavl,layer_in,layer_out,tmp_k, &
              Ed0_nir_dir,Ed0_nir_dif,wv_mul
          double precision, dimension(2*z_max) ::    &
              fdirup,  & ! up   flux at model interface due to direct beam at top surface
              fdirdn,  & ! down flux at model interface due to direct beam at top surface
              fdifup,  & ! up   flux at model interface due to diffuse beam at top surface
              fdifdn     ! down flux at model interface due to diffuse beam at top surface
          double precision, dimension(z_max,lda_n) ::    &
              pur_dble,  &    ! temp var for storing current pur in double format
              Ik1,       &    ! calcuated Ik prime
              ll              ! light limitation term
          ! vars required by sia2_env_grow_algae.inc.f90
          double precision :: Gmax,NO3_lim,NH4_lim,SiO4_lim,PO4_lim,light_lim,sal_lim
          double precision :: min_lim, growth, PUR_mean24, death, newpN
          double precision :: newpP,newpSi,newp,min_alg_brine,totalp,p_tot_int,p_tot_bot
          double precision :: smalg_temp,dNO3_temp,dNH4_temp,dPO4_temp,dSiO4_temp
          double precision :: production_temp,dN_all_temp,resp,poc_temp,bot_id
          !double precision, dimension(z_max) :: smalg_new,NO3_new,NH4_new,PO4_new,SiOH4_new
          double precision, dimension(z_max) :: PUR_new, PUR_temp, dN_all, Ed_W_temp
          double precision, dimension(z_max) :: Ed_W_mean,Ed_W_snow
          double precision, dimension(z_max+1) :: Ed_IR

!          double precision :: zenith_deg,kl_M,kl_Tr,kl_alpha,kl_Bn,kl_Ta,kl_Tw
!          double precision :: kl_Moz,kl_Ozone,kl_T0,kl_au,kl_Tu,kl_Id,kl_ao,kl_PARID
!          double precision :: kl_PARD,kl_PARG,kl_Dr,kl_Da,kl_v,kl_w,kl_x,kl_y,kl_z
!          double precision :: kl_Dm,kl_d,kl_g,kl_AtmAlbedo,kl_cc,elevation_deg
!          double precision, dimension (wavl) :: kl_par

          type (ice_type) :: ice_temp  ! for use in splitting partial categories

          ! interfaces
          !----------------------------------------------------

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


          ! whole model domain calculations for current timestep
          ! ------------------------------------------------------------

          ! var to see if ice gets too deep anywhere for current physical model
          max_d_exceeded = 0
          maxed_depth = .false.
          ignore_h_max = .false.

          ! convert julian date to radians
          j_date = cur_hour/24+1
          rad_date = (j_date/365)*2.*pi

          ! Calculate solar declination (between -23deg27' and 23deg27') for use in
          ! Solar elevation equation
          solar_dec = 0.39637 - 22.9133*cos(rad_date) + 4.02543*sin(rad_date)  &
              - 0.3872*cos(2*rad_date) + 0.052*sin(2*rad_date)

          ! Convert degrees to radians
          solar_dec = pi*solar_dec/180.

          ! Convert time to radians - these equations seem to give midnight as noon,
          ! so I added pi to the time to correct -
          rad_time = (j_date - int(j_date))*2.*pi - pi

          ! set variable timestep counters to zero
          n_dtt_1 = 0
          n_dtt_2 = 0
          n_dtt_3 = 0

          ! chunk size calc
          chunk_size = 20

          ! initialize ed_loaded so it is not done by multiple threads
          load_ed = .true.

          ! determine whether to calc limitation variable for write out
          find_lim = .FALSE.
          if (hour24 .eq. 12) then
              find_lim = .TRUE.
          endif

          ! for use if boundary_file
          Edir_last_valid = 0.d0
          Edif_last_valid = 0.d0

!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(mi,i,j,ii,iii,jj,jjj,kk,ic,sc,i_ease,j_ease,i_mp,j_mp,thread,slot, &
!$OMP i_dg,j_dg,now,sk_1,sk_z,int_z,ic_out, &
!$OMP airtemp1c,snowfall,snowfall_d,t_depth,snowh1m,tmp1,tmp2,tmp3,tmp4, &
!$OMP t_mean,s_mean,d_mean,bs_mean,bd_mean,bv_mean,heat_mean, &
!$OMP alpha0,alpha1,alpha2,alpha3, &
!$OMP dz_mean,aq,bq,cq, &
!$OMP keep_growing,maxed_depth,eoy,valid_snowh,ignore_h_max, &
!$OMP int_z_new,z_last,z_snow_new, &
!$OMP r_depth,target_depth,layer_divisor,sh_diff_mean, &
!$OMP sh_increment,sh_diff, &
!$OMP th_new,id_new, &
!$OMP interim_top,interim_bot,old_layer_top,new_layer_bot, &
!$OMP old_layer_bot,z1,z2,dz,heat_total,dz_total,new_layer_top, &
!$OMP flooded,t_flood,s_flood,d_flood,bs_flood,bd_flood, &
!$OMP bv_flood,heat_flood,t_interface,d_total,af_total,vf_total, &
!$OMP s_pond_new,h_pond_new,d_pond_new,smalg_pond_new, &
!$OMP poc_pond_new,no3_pond_new,nh4_pond_new,po4_pond_new, &
!$OMP sioh4_pond_new,th_pond_new,af_tot,sh_mean, &
!$OMP s_flood_new,heat_flood_new,d_flood_new,no3_flood_new, &
!$OMP nh4_flood_new,po4_flood_new,sioh4_flood_new,poc_flood_new, &
!$OMP t_new,s_new,d_new,bs_new,bd_new,t_change, &
!$OMP bv_new,smalg_new,poc_new,smalg_pre_map, &
!$OMP melt_new, &
!$OMP NO3_new,NH4_new,PO4_new,SiOH4_new, &
!$OMP integr_valid,did_flood,alg_fixed,remapping,flood_1x,did_drain, &
!$OMP info,m_order,z_snow,ice_1,ice_z,dtt,mo,heat_stop,snow_model_temp, &
!$OMP sk_gl,c_gl,c_gl_sal,s_gl,sk_cur,smalg_lost,fc_int, &
!$OMP p_h2o_sat,row_air,rhum,shum,epsa,v10,Fe_sub1,Fe_sub2,Fs_sub, &
!$OMP FL,albedo,eps0,ks,top_a,top_b,Fr,F0_constants,bot_a,bot_b, &
!$OMP F0_4,F0_3,F0_2,F0_1,dtt_s_snowh,ki_up,ki_down,Ts_kelvin, &
!$OMP Ts_last,Tsi_last,Ts_next,Tsi_next,T_diff,c_snow,F0,dF0,Fc, &
!$OMP mat_sd,mat_sb,mat_sc,mat_0a,mat_0b,mat_0c,rs,r_0,xl,xls, &
!$OMP th1,th2,thz,th2z,Fc_top,Fc_bot,albedo_s,albedo_i,f_snow, &
!$OMP mat_a,mat_b,mat_c,r1,dhdt_g_m2,dhdt,c_gl_bv,dhdt_cm_s, &
!$OMP keff,dz_c,dz_sk,sk_th_new,f_depth,ts_prev,th_crunch, &
!$OMP c_gl_tot,s_gl_tot,dtt_s,T_step,snow_gl,crunch_tot,fl_gl, &
!$OMP t_melt,q_melt,dtt_days,snow_melt,snow_melt_now,melt_flood, &
!$OMP mf_depth,melt_drain,ki1,heat_brine,c_gl_copy, &
!$OMP T_up, T_dn, T_ii, id_mean,potT,shcoef,lhcoef,p_h2o, &
!$OMP Fm0,F_little_L,Fe,F_s, &
!$OMP dNO3,dNH4,dPO4,dSiOH4,dpoc,dsmalg,dth_ed_w, &
!$OMP dth,t_next,t_last,th, &
!$OMP DU,DU_calc,DL,DL_calc, &
!$OMP DC,DC_calc,r_vec,ki,ci,ki_alt, &
!$OMP heat_temp,dcheat,debug_z,debug_start, &
!$OMP t_snow_temp,d_snow_new,th_snow_new,heat_snow_new, &
!$OMP ts_snow,Fc_snow, &
!$OMP vb_open,cv_open,bfrf, &
!$OMP T_grad,bs_dilute,dsdt1,dsdt2,bulk_nc,dbvdt, &
!$OMP dsdt3,dsdt4,t_brine,d_brine,dsdt_net,fresh_net,dsdt_tmp, &
!$OMP bs_dif,superimposed,fresh_tmp, &
!$OMP dsdt,dheat,dheat_last,dsdt_cw, &
!$OMP bs_last,bd_last,bv_last,t_prev,t_snow_prev, &
!$OMP heat_last,bs_diff,heat_diff            , &
!$OMP flux_steps,conv_start, &
!$OMP flux_down,th_down,bv_down,Fd,Fu,Fdiff,Feddy,fvs, &
!$OMP dt_flux,flux_ml,flux_d,hv,flux_sk,fvb,fs, &
!$OMP b_flux_max,min_perm,bv_flux,bv_total,visc, &
!$OMP no3_total,nh4_total,po4_total,sioh4_total,poc_total, &
!$OMP b_flux,porosity,perm,pres,p_flux,h_flux, &
!$OMP fbh_new,top_f,bot_f,ice_mass,snow_mass,snow_dist_new,snow_rd_new, &
!$OMP nm_i,m_row,a_i,b_i,c_i,d_i,srftyp, &
!!$OMP ipiv, &
!$OMP bv_mean_real, &
!$OMP ex, &
!$OMP zenith_angle,E_sum,lda_f,par_cf_f,kdice_ir,Ed_W_IR, &
!$OMP rsun0,rsky0,temp_brine_v,am_sum,Kds,testvar, &
!$OMP Ed0_new,Ed0_temp,Ed1_new,Ed1_temp,Ed_temp,kds_z, &
!$OMP mean_spec,absorp,ed_med,albodr,albodf,coszen, &
!!$OMP Ed_snow1_5,Ed_snow4,Ed_snow10, &
!$OMP Kdice,Kdp,am,ad, &
!$OMP Ed, &
!!$OMP ABCD, &
!!$OMP AB_mat, &
!$OMP k2st,r2st,ex2st,g2st,w2st,tau2st, &
!$OMP trndir, trntdr, trndif, rupdir, rupdif, rdndif, &
!$OMP refk,par_top,par_bot,par_mid,heat_wavl,layer_in,layer_out,tmp_k, &
!$OMP Ed0_nir_dir,Ed0_nir_dif,wv_mul, &
!$OMP fdirup, &
!$OMP fdirdn, &
!$OMP fdifup, &
!$OMP fdifdn, &
!$OMP pur_dble,  &
!$OMP Ik1,       &
!$OMP ll,        &
!$OMP Gmax,NO3_lim,NH4_lim,SiO4_lim,PO4_lim,light_lim,sal_lim, &
!$OMP min_lim, growth, PUR_mean24, death, newpN, &
!$OMP newpP,newpSi,newp,min_alg_brine,totalp,p_tot_int,p_tot_bot, &
!$OMP smalg_temp,dNO3_temp,dNH4_temp,dPO4_temp,dSiO4_temp, &
!$OMP production_temp,dN_all_temp,resp,poc_temp,bot_id, &
!$OMP PUR_new, PUR_temp, dN_all, Ed_W_temp, &
!$OMP Ed_W_mean,Ed_W_snow, &
!$OMP Ed_IR, Edif_last_valid, Edir_last_valid, &
!$OMP ice_temp)
!!$OMP zenith_deg,kl_M,kl_Tr,kl_alpha,kl_Bn,kl_Ta,kl_Tw, &
!!$OMP kl_Moz,kl_Ozone,kl_T0,kl_au,kl_Tu,kl_Id,kl_ao,kl_PARID, &
!!$OMP kl_PARD,kl_PARG,kl_Dr,kl_Da,kl_v,kl_w,kl_x,kl_y,kl_z, &
!!$OMP kl_Dm,kl_d,kl_g,kl_AtmAlbedo,kl_cc,elevation_deg, &
!!$OMP kl_par, &


!$OMP MASTER


          if (write_disk .eq. 1) then

#include "sia2_env_write_to_disk.inc.f90"

          endif

          kk = int(cur_hour+dt)+1
          if (load_ed .and. (kk .gt. last_hour)) then

              load_ed = .false.

              if (f_index .eq. 1) then
                  f_index_next = 2
              else
                  f_index_next = 1
              endif

              if (leap_year .and. kk .gt. 8784) then
                  jj = cur_year+1
                  kk = 1
              elseif (kk .gt. 8760) then
                  jj = cur_year+1
                  kk = 1
              else
                  jj = cur_year
              endif

              print *, jj,kk,'loading Edif...'
              !  Open Gregg & Carter diffuse surface spectral PAR, cloud corrected
              !-----------------------------------------------------------------
              call gc_read(jj,kk,mp_f%Edif(f_index_next), 2)

              print *, jj,kk,'loading Edir...'
              !  Open Gregg & Carter direct surface spectral PAR, cloud corrected
              !-----------------------------------------------------------------
              call gc_read(jj,kk,mp_f%Edir(f_index_next), 1)

          endif
!$OMP END MASTER


!$OMP DO SCHEDULE(DYNAMIC,chunk_size)


          ! Run environment calculations for model domain only, as specified in the
          ! constants.txt file
          ! ----------------------------------------------------------------
              do mi=1,tcells

!                  if (boundary_file .eq. 1) then
                      ! fake EASE grid stuff
!                      i_ease = 373
!                      j_ease = 425
                      ! print *,'Cloud Fraction: ', f(mi)%fc
!                  else
                      ! get current model index (mi)
                      ! mi = mem_map(slot,thread)
                      ! find full EASE domain reference for current grid cell
                      i_ease = m(mi)%grid_h
                      j_ease = m(mi)%grid_v
!                  endif

                  ! find mercator projection indexes for current grid cell
                  i_mp = m(mi)%x_mp  ! mercator/ncep projection x
                  j_mp = m(mi)%y_mp  ! mercator/ncep projection y

                  ! find degree projection indexes for current grid cell
                  i_dg = m(mi)%x_dg  ! mercator/ncep projection x
                  j_dg = m(mi)%y_dg  ! mercator/ncep projection y

                  ! find air temp in degC
                  airtemp1c = f(mi)%at - kelvin0

                  ! find density of snowfall - perhaps this can get more complicated sometime
                  if (airtemp1c .gt. -2.) then
                      snowfall_d = den_s_wet*1.e6
                  else
                      snowfall_d = den_s_dry*1.e6
                  endif

                  if (m(mi)%status .eq. 0) then
                      ! check to make sure that minimum ice concentration is met...
                      ! -------------------------------------------------------------
                      if ((boundary_file .eq. 1) .or. (override_ic .eq. 1) .or. &
                          ((f(mi)%ic_interp .ge. fi) .and. (f(mi)%ic_interp .le. 1.) &
                          .and. (ease_f%ice_vec(i_ease,j_ease,3) .gt. 0.))) then
                          ! check to see if start, set heights, salinities, nutrients, bio
                          ! -------------------------------------------------------------

! ***********************************************************************
! ***********************************************************************
! INCLUDED CREATE ICE CODE HERE
! ***********************************************************************
! ***********************************************************************
!#include "sia2_env_create_ice.inc.f90"
                              call sia2_env_create_ice(m,ice,f,v_f,pl,mi)

                      endif ! end of ice concentration check

                  endif ! end of status=0/new ice check

                  ! compute cell if active
                  ! -------------------------------------------------------------
                  if ((boundary_file .eq. 1) .or. (override_ic .eq. 1) .or. &
                  ((m(mi)%status .ge. 1) .and. ((f(mi)%ic_interp .ge. fi) &
                  .and. (f(mi)%ic_interp .le. 1.)))) then

#include "sia2_env_redist.inc.f90"

#include "sia2_env_snow_allocate.inc.f90"

!              do sc=1,sda_n
!              do ic=1,ida_n
!                  if (ice(sc,ic,mi)%af .gt. 0.) then
!                       tmp1=0.
!                       do ii=1,ice(sc,ic,mi)%z
!                       if (ice(sc,ic,mi)%smalg(ii) .gt. 1.e6) then
!                           tmp1 = 1.
!                       endif
!                       enddo
!                       if (tmp1 .eq. 1.) then
!                           print *,'After redist, before ice model'
!                           print *, mi,'--',ic,' ------------------------------------'
!                           print *,ice(sc,ic,mi)%smalg
!                           testvar = -1
!                           testvar = sqrt(testvar)
!                       endif
!                   endif
!              enddo
!              enddo


                      do sc=1,sda_n
                      do ic=1,ida_n
              if ((ice(sc,ic,mi)%z > 0) .and. (ice(sc,ic,mi)%af .gt. 0.d0)) then
! ***********************************************************************
! ***********************************************************************
! INCLUDED PHYSICAL ICE MODEL CODE HERE
! ***********************************************************************
! ***********************************************************************
#include "sia2_env_icemodel_2.inc.f90"

                          endif

                      enddo
                      enddo

!              do sc=1,sda_n
!              do ic=1,ida_n
!                  if (ice(sc,ic,mi)%af .gt. 0.) then
!                       tmp1=0.
!                       do ii=1,ice(sc,ic,mi)%z
!                       if (ice(sc,ic,mi)%smalg(ii) .gt. 1.e6) then
!                           tmp1 = 1.
!                       endif
!                       enddo
!                       if (tmp1 .eq. 1.) then
!                           print *,'After ice model'
!                           print *, mi,'--',ic,' ------------------------------------'
!                           print *,ice(sc,ic,mi)%smalg
!                           testvar = -1
!                           testvar = sqrt(testvar)
!                       endif
!                   endif
!              enddo
!              enddo



                      ! caculate & encode limitation terms
                      if (find_lim) then
                          tmp1 = 0.d0
                          t_mean = 0.d0
                          NO3_lim = 0.d0
                          SiO4_lim = 0.d0
                          PO4_lim = 0.d0
                          light_lim = 0.d0
                          sal_lim = 0.d0
                          do ic=1,ida_n
                          do sc=1,sda_n
                              do ii=1,ice(sc,ic,mi)%z
                                  !tmp3 = ice(sc,ic,mi)%smalg(ii)*ice(sc,ic,mi)%bv(ii)* &
                                    !ice(sc,ic,mi)%th(ii)*cell_area*ice(sc,ic,mi)%af
                                    tmp3 = ice(sc,ic,mi)%th(ii)*ice(sc,ic,mi)%af

                                  tmp1 = tmp1 + tmp3
                                  t_mean = t_mean + ice(sc,ic,mi)%t(ii)*tmp3
                                  NO3_lim = NO3_lim + ice(sc,ic,mi)%nlim(ii)*tmp3
                                  PO4_lim = PO4_lim + ice(sc,ic,mi)%plim(ii)*tmp3
                                  SiO4_lim = SiO4_lim + ice(sc,ic,mi)%silim(ii)*tmp3
                                  light_lim = light_lim + ice(sc,ic,mi)%llim(ii)*tmp3
                                  sal_lim = sal_lim + ice(sc,ic,mi)%slim(ii)*tmp3
                              enddo
                          enddo
                          enddo
                          m(mi)%tlim = t_mean/tmp1
                          m(mi)%nlim = NO3_lim/tmp1
                          m(mi)%silim = SiO4_lim/tmp1
                          m(mi)%plim = PO4_lim/tmp1
                          m(mi)%llim = light_lim/tmp1
                          m(mi)%slim = sal_lim/tmp1
                      endif

                  else

                      if (m(mi)%status .gt. 0) then ! don't remember anything.

                          ! null 2d params
                          m(mi)%status = 0
                          adv(mi)%af_new = 0.

                          do ic=1,ida_n
                          do sc=1,sda_n

                              if (ice(sc,ic,mi)%af .gt. 0.) then

                                  ! record ice volume as model domain loss
                                  int_z = ice(sc,ic,mi)%z - z_sk
                                  m(mi)%md_loss = m(mi)%md_loss + &
                                    ice(sc,ic,mi)%id(ice(sc,ic,mi)%z)* &
                                    ice(sc,ic,mi)%af*cell_area*1.e6

                                  ! record biomass mass lost from ice / model domain
                                  do ii=1,ice(sc,ic,mi)%z
                                      if (ice(sc,ic,mi)%smalg(ii) .gt. 0.) then
                                           bv_mean = ice(sc,ic,mi)%bv(ii)*c_001
                                           smalg_lost = (ice(sc,ic,mi)%smalg(ii)*bv_mean - min_alg)
                                           smalg_lost = max(0.,smalg_lost)
                                          ! converting to g/pixel - (cell_area*1.e6)*1g/1000mg = 1000. scale factor
                                          m(mi)%bm_lost = m(mi)%bm_lost + smalg_lost* &
                                          ice(sc,ic,mi)%th(ii)*cell_area*1.e3*ice(sc,ic,mi)%af   ! g/pixel
                                      endif
                                enddo


                                    call sia2_ice_mass(ice(sc,ic,mi),fresh_tmp,dsdt_tmp)

                                    ! record fresh water flux (kg)
                                    m(mi)%h2o_flux = m(mi)%h2o_flux + &
                                        fresh_tmp*ice(sc,ic,mi)%af*cell_area*1.e6

                                    ! record salt flux (kg)
                                    m(mi)%salt_flux = m(mi)%salt_flux + &
                                        dsdt_tmp*ice(sc,ic,mi)%af*cell_area*1.e6



                              endif

#include "sia2_env_null_params.f90"
           pur(:,:,:,sc,ic,mi) = pur_0

                          enddo
                          enddo

                      endif  ! end status=1 test

                  endif ! end ice conentration test

          enddo ! end tcells loop

!$OMP END DO
!$OMP END PARALLEL

          if (max_d_exceeded .gt. 0) then
              print *, 'Maximum depth exceeded in ',max_d_exceeded,' cells.'
          endif
          print *, 'DTT 1,2,3: ',n_dtt_1,n_dtt_2,n_dtt_3



      end SUBROUTINE sia2_environment_calc



! SUBROUTINE: sia2_interpolate_salinity
! used to find intital conditions from a 10-value distribution
! for 1st year ice salinity
! ======================================================================

      PURE SUBROUTINE sia2_interpolate_salinity1(h_rel,value)
          use sia2_globals
          implicit none

          double precision,intent(in) :: h_rel
          double precision, intent(out) :: value

          !h_rel is the height of the desired measuremnt normalized to 1.
          !i.e., h_rel = interp_height/total_height

          ! interpolate between 10 points centered on integers -
          ! this means the below 0.05 an abosse 9.5 values are
          ! simply clipped to the end values
          !------------------------------------------------
          if (h_rel .lt. 0.05) then
              value=ss0
          elseif ((h_rel .ge. 0.05) .and. (h_rel .lt. 0.15)) then
              value = ss0+(ss1-ss0)*(h_rel-0.05)/0.1
          elseif ((h_rel .ge. 0.15) .and. (h_rel .lt. 0.25)) then
              value = ss1+(ss2-ss1)*(h_rel-0.15)/0.1
          elseif ((h_rel .ge. 0.25) .and. (h_rel .lt. 0.35)) then
              value = ss2+(ss3-ss2)*(h_rel-0.25)/0.1
          elseif ((h_rel .ge. 0.35) .and. (h_rel .lt. 0.45)) then
              value = ss3+(ss4-ss3)*(h_rel-0.35)/0.1
          elseif ((h_rel .ge. 0.45) .and. (h_rel .lt. 0.55)) then
              value = ss4+(ss5-ss4)*(h_rel-0.45)/0.1
          elseif ((h_rel .ge. 0.55) .and. (h_rel .lt. 0.65)) then
              value = ss5+(ss6-ss5)*(h_rel-0.55)/0.1
          elseif ((h_rel .ge. 0.65) .and. (h_rel .lt. 0.75)) then
              value = ss6+(ss7-ss6)*(h_rel-0.65)/0.1
          elseif ((h_rel .ge. 0.75) .and. (h_rel .lt. 0.85)) then
              value = ss7+(ss8-ss7)*(h_rel-0.75)/0.1
          elseif ((h_rel .ge. 0.85) .and. (h_rel .lt. 0.95)) then
              value = ss8+(ss9-ss8)*(h_rel-0.85)/0.1
          elseif ((h_rel .ge. 0.95) .and. (h_rel .le. 1.0)) then
              value = ss9
          endif

          !print *,'interp:', scale,h,h_rel,value


      end SUBROUTINE sia2_interpolate_salinity1

! SUBROUTINE: sia2_interpolate_salinity2
! used to find intital conditions from a 10-value distribution
! for multi-year ice salinity
! ======================================================================

      PURE SUBROUTINE sia2_interpolate_salinity2(h_rel,value)
          use sia2_globals
          implicit none

          double precision,intent(in) :: h_rel
          double precision, intent(out) :: value

          !h_rel is the height of the desired measuremnt normalized to 1.
          !i.e., h_rel = interp_height/total_height

          ! interpolate between 10 points centered on integers -
          ! this means the below 0.05 an abosse 9.5 values are
          ! simply clipped to the end values
          !------------------------------------------------
          if (h_rel .lt. 0.05) then
              value=ssm0
          elseif ((h_rel .ge. 0.05) .and. (h_rel .lt. 0.15)) then
              value = ssm0+(ssm1-ssm0)*(h_rel-0.05)/0.1
          elseif ((h_rel .ge. 0.15) .and. (h_rel .lt. 0.25)) then
              value = ssm1+(ssm2-ssm1)*(h_rel-0.15)/0.1
          elseif ((h_rel .ge. 0.25) .and. (h_rel .lt. 0.35)) then
              value = ssm2+(ssm3-ssm2)*(h_rel-0.25)/0.1
          elseif ((h_rel .ge. 0.35) .and. (h_rel .lt. 0.45)) then
              value = ssm3+(ssm4-ssm3)*(h_rel-0.35)/0.1
          elseif ((h_rel .ge. 0.45) .and. (h_rel .lt. 0.55)) then
              value = ssm4+(ssm5-ssm4)*(h_rel-0.45)/0.1
          elseif ((h_rel .ge. 0.55) .and. (h_rel .lt. 0.65)) then
              value = ssm5+(ssm6-ssm5)*(h_rel-0.55)/0.1
          elseif ((h_rel .ge. 0.65) .and. (h_rel .lt. 0.75)) then
              value = ssm6+(ssm7-ssm6)*(h_rel-0.65)/0.1
          elseif ((h_rel .ge. 0.75) .and. (h_rel .lt. 0.85)) then
              value = ssm7+(ssm8-ssm7)*(h_rel-0.75)/0.1
          elseif ((h_rel .ge. 0.85) .and. (h_rel .lt. 0.95)) then
              value = ssm8+(ssm9-ssm8)*(h_rel-0.85)/0.1
          elseif ((h_rel .ge. 0.95) .and. (h_rel .le. 1.0)) then
              value = ssm9
          endif

          !print *,'interp:', scale,h,h_rel,value


      end SUBROUTINE sia2_interpolate_salinity2


    SUBROUTINE sia2_interp_value_from_data(array,x,len,x1,value)

      implicit none

      double precision, intent(in) :: array(len)
      double precision, intent(in) :: x(len)
      double precision, intent(in) :: x1
      integer, intent(in) :: len
      double precision, intent(out) :: value
      integer :: i

      do i=1,len
        if (x1 <= x(i)) then
          exit
        endif
      enddo

      if ((i >= len) .or. (x1 >= x(len))) then
        value = array(len)  ! z beyond range
      elseif (i == 1) then
        value = array(1)    ! z before range
      else
        value = array(i-1)+(array(i)-array(i-1))*(x1-x(i-1))/(x(i)-x(i-1))
      endif

      return

    end SUBROUTINE sia2_interp_value_from_data



! SUBROUTINE: sia2_ice_mass
! used to find masses of h2o and salt for flux accounting purposes
! note: does not take ice area into account - masses are in kg/m^2
! ======================================================================
    subroutine sia2_ice_mass(ice,h20,salt)

        use sia2_globals

        ! Input
        type (ice_type), intent(in) :: &
            ice             ! ice structure

        ! out
        double precision, intent(out) :: &
            h20,          & ! total h20 mass (kg/m^2)
            salt            ! total salt mass (kg/m^2)

        ! Internal Variables
        double precision :: &
            mass            ! layer mass (kg/m^2)

        h20 = 0.d0
        salt = 0.d0
      do ii=1,ice%snow%z
          h20 = h20 + &
              ice%snow%d(ii)*ice%snow%th(ii)*1.d-3 ! kg h2o
      enddo
      do ii=1,ice%z
          mass = ice%th(ii)*ice%d(ii)*1.d-3
          h20 = h20 + mass*(1.d3 / (1.d3 + ice%s(ii)))  ! kg h2o- bubble fraction already accounted for in density
          salt = salt + mass* &
              (ice%s(ii) / (1.d3 + ice%s(ii))) ! kg salt
      enddo
      if (use_ponds .eq. 1 .and. ice%pond%th .gt. 0.d0) then
          mass = ice%pond%th*ice%pond%d*1.d-3
          h20 = h20 + mass*(1.d3 / (1.d3 + ice%pond%s)) ! kg h2o
          salt = salt + mass* &
              (ice%pond%s / (1.d3 + ice%pond%s)) ! kg salt
      endif

    end subroutine sia2_ice_mass


! SUBROUTINE: sia2_snow_mass
! used to find mass of h2o int snow for flux accounting purposes
! note: does not take ice area into account - masses are in kg/m^2
! ======================================================================
    subroutine sia2_snow_mass(ice,h20)

        use sia2_globals

        ! Input
        type (ice_type), intent(in) :: &
            ice             ! ice structure

        ! out
        double precision, intent(out) :: &
            h20            ! total h20 mass (kg/m^2)

        ! Internal Variables
        double precision :: &
            mass            ! layer mass (kg/m^2)

        h20 = 0.d0
      do ii=1,ice%snow%z
          h20 = h20 + &
              ice%snow%d(ii)*ice%snow%th(ii)*1.d-3 ! kg h2o
      enddo

    end subroutine sia2_snow_mass


! SUBROUTINE: siesta_remap_tracer
! ======================================================================
    subroutine siesta_remap_tracer(th_old,th_new,tracer,offset_down, &
                                   int_z,int_z_new,stretch,offset_top,tracer_new)

        use sia2_globals
        implicit none

        double precision, dimension(z_max),intent(in) :: th_old,th_new,tracer
        double precision, intent(in) :: offset_down,offset_top
        integer, intent(in) :: int_z,int_z_new
        logical, intent(in) :: stretch
        double precision, dimension(z_max), intent(inout) :: tracer_new

        double precision, dimension(z_max) :: tracer_mod, th_mod
        double precision :: dz_total,new_layer_top,new_layer_bot
        double precision :: old_layer_bot,old_layer_top,interim_top,interim_bot
        double precision :: dz,z1,z2,offset_total
        logical :: mapping_started
        integer :: ii,sk_1,sk_z,jj


        old_layer_top = 0.
        new_layer_bot = 0.
        old_layer_bot = 0.
        if (offset_top .ne. 0.d0) then
          new_layer_top = offset_top
        else
          new_layer_top = -1.
        endif

        sk_1 = int_z + 1
        sk_z = int_z+z_sk

        print *,'offset_down',offset_down
        print *,'offset_top',offset_top

        ! defaults
        th_mod = th_old  ! vector assign
        tracer_mod = tracer
        offset_total = 0.d0
        tracer_new = 0.d0  ! vector assign
        mapping_started = .FALSE.


        ! decide how to deal with changes at bottom of ice pack - offset_down
        if (offset_down .ne. 0.) then
          if (stretch .or. offset_down .lt. 0.d0) then
              dz = 0.d0
              do ii=1,sk_z
                dz = dz + th_old(ii)
              enddo
              z1 = (offset_down + dz) / dz ! fractional change
              th_mod = th_old*z1
              tracer_mod = tracer / z1
              offset_total = 0.d0

              print *,'stretching/shrinking'
              print *,'z1',z1
              print *,'dz',dz
              print *,'th_old',th_old
              print *,'th_mod',th_mod
              print *,'th_new',th_new
              print *,'tracer',tracer
              print *,'tracer',tracer_mod
          else
              offset_total = offset_down
          endif
        endif


        do ii=1,int_z_new+z_sk

            if (offset_total .gt. 0.d0) then
              if (offset_total .gt. th_new(ii)) then
                dz_total = th_new(ii)
                offset_total = offset_total - th_new(ii)
              else
                dz_total = offset_total
                offset_total = 0.d0
              endif
            else
              dz_total = 0.d0
            endif

            if (dz_total .lt. th_new(ii)) then

                ! find heights of new layer with regard to old thicknesses
                ! new_layer_top should be zero even if flooding has occured,
                ! because it is really a reference to the old ice pack.
                if (.NOT. mapping_started) then
                    mapping_started = .TRUE.
                    if (new_layer_top .lt. 0.d0) then
                        new_layer_top = 0.d0
                    endif
                    new_layer_bot = new_layer_top + th_new(ii) - dz_total
                else
                    new_layer_top = new_layer_bot
                    new_layer_bot = new_layer_bot + th_new(ii)
                endif

                ! find 1st old layer that contains part of new layer, going down
                jj = 0
                old_layer_bot = 0.d0
                do while(new_layer_top .ge. old_layer_bot .and. jj .le. sk_z)
                    jj=jj+1
                    old_layer_bot = old_layer_bot + th_mod(jj)
                enddo
                if (jj .eq. 0) then
                  print *,ii
                  print *,new_layer_top
                  print *,new_layer_bot
                  print *,old_layer_top
                  print *,th_new
                  print *,th_mod
                  print *,tracer
                  print *,tracer_new
                endif

                old_layer_top = old_layer_bot - th_mod(jj)
                ! now jj is OLD layer where NEW layer ii starts...

                ! find total from multiple layers/partial layers that make up
                do while ((old_layer_top .lt. new_layer_bot) .and. (jj .le. sk_z))

                    ! ----- NEW LAYER GEOMETRIES ------------------------------
                    interim_top = max(old_layer_top,new_layer_top)
                    interim_bot = min(old_layer_bot,new_layer_bot)
                    z1 = interim_top - old_layer_top
                    z2 = interim_bot - old_layer_top
                    dz = z2-z1

                    tracer_new(ii) = tracer_new(ii) + tracer(jj)*dz

                    !print *,ii
                    !print *,new_layer_top,new_layer_bot
                    !print *,old_layer_top,old_layer_bot
                    !print *,jj,dz

                    ! setup variable for next partial layer addition
                    dz_total = dz_total + dz

                    jj=jj+1
                    if (jj .le. sk_z) then
                        ! find boundaries of next old layer
                        old_layer_top = old_layer_bot
                        old_layer_bot = old_layer_bot + th_mod(jj)
                    endif

                enddo

                ! round out dz_total whether we are shrinking or growing,
                ! otherwise division below can create/take away tracers
                dz_total = th_new(ii)

            endif     ! end of "if dz_total .lt. th_new(ii)

            ! divide by thickness to return to m-3 from m-2 concentration
            tracer_new(ii) = tracer_new(ii) / th_new(ii)
            !print *,'ii,tracer',ii,tracer_new(ii),th_new(ii)

        enddo ! end of new grid portioning
        print *,'tracer_new', tracer_new

    end subroutine siesta_remap_tracer




! FUNCTION: explt
! calculates approximate exp function from linear approximation lookup table
! ======================================================================

      FUNCTION explt(n,use_expt,exp_bins,exp_mul,expt)

          double precision :: n         ! exp argument
          integer :: use_expt           ! switch
          integer :: exp_bins           ! number of indexes to lookup table
          double precision :: exp_mul   ! muliplier to convert n to appropriate index
          double precision, dimension(2,exp_bins) :: expt
          double precision :: explt
          integer :: i

          if (use_expt .eq. 1) then
              i=int(anint(n*exp_mul))       ! find lookup table index
              if (i .lt. 1 .or. i .gt. exp_bins) then
                  explt = exp(n)                  ! outside of lookup table bounds, use system exp()
              else
                  explt = n*expt(1,i) + expt(2,i) ! linear approximation for exp
              endif
          else
              explt = exp(n)
          endif

          !print *, n,abs(exp(n)-explt)

      end FUNCTION explt


#include "sia2_env_merge_ice.f90"


#include "sia2_env_atmo.f90"


