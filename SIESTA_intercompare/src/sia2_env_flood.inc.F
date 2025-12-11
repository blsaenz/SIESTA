
          ! Calculate Flooding
		  ! -----------------------------------------------------------

          ! find mean density
          ice_mass = 0.
          do ii = 1,sk_z
              ice_mass = ice_mass + ice(sc,ic,mi)%d(ii)*ice(sc,ic,mi)%th(ii)
          enddo
		  ! find snow mass
		  snow_mass = 0.
		  do ii = 1,ice(sc,ic,mi)%snow%z
		      snow_mass = snow_mass + ice(sc,ic,mi)%snow%d(ii) &
		        *ice(sc,ic,mi)%snow%th(ii)
		  enddo

          ! Calculate freeboard height of mean snow depth
          fbh_new = ice(sc,ic,mi)%id(sk_z) - (snow_mass + ice_mass)/f(mi)%d ! m


          ! determine if bine channels are open enough to flood
          tmp1 = 1.
          do ii=1,int_z
              if (ice(sc,ic,mi)%bv(ii) .lt. 70.) then
                  tmp1 = 0.
              endif
          enddo

          if (dt_step .eq. 372 .and. mi .eq. 41151 .and. ic .eq. 1) then
          testvar = -1
          endif


          if (tmp1 .gt. 0.d0 .and. flood_1x .and. (no_flooding .eq. 0) &
              .and. snow_mass .gt. 0.d0) then

          ! Regular ice flooding - calculate & execute
          ! transfer to flooded category
		  ! -----------------------------------------------------------
          if (sda_n .gt. 1 .and. sc .eq. 1) then

              ! find lognormal distribtion mean, for estimating true snow depth distribution
              if (ice(1,ic,mi)%snow_dist .lt. dble(sd_num)) then
                  tmp2 = 0.d0
                  tmp1 = ice(1,ic,mi)%snow_dist
                  do jj=1,sd_num
                      tmp2 = tmp2 + min(1.d0,tmp1)*sda_multiplier(jj)
                      tmp1 = max(0.d0,tmp1-1.d0)
                  enddo
                  if (tmp2 .gt. 0.d0) then
                      sh_mean = ice(1,ic,mi)%snow%depth* &
                         ice(1,ic,mi)%snow_dist/tmp2
                  else
                      sh_mean = 0.d0
                  endif
              else
                  sh_mean = ice(1,ic,mi)%snow%depth
              endif

              ! starting with tallest snow distribution available, work downward to
              ! find how much of the snow distribution is flooding,
              ! and calcualte how to parition snow
              af_total = 0.d0
              flooded = 0.d0
              tmp4 = 0.d0
              snow_dist_new = ice(1,ic,mi)%snow_dist

              do jj=1,sd_num

                  ! limit flooding to 0.05 dist increments, to prevent
                  ! numerical errors
                  tmp1 = ice(1,ic,mi)%snow_dist - dble(jj-1)

                  ! check for flooding if whole snow depth increment is present, or
                  ! flooding has already started in shorter categories
                  if (tmp1 .ge. 1.d0 .or. (flooded .gt. 0.d0 .and. tmp1 .ge. 0.05)) then

                      ! find new snow mass
                      tmp1 = snow_mass * sda_multiplier(jj) * sh_mean / ice(1,ic,mi)%snow%depth ! snow_mass was calc'ed for current mean snow depth - make into distribution mass

                      ! find freeboard for snow category
                      tmp2 = ice(sc,ic,mi)%id(sk_z) - (tmp1 + ice_mass)/f(mi)%d ! m

                      ! test to see if freeboard is underwater,
                      ! if so find flood depth and areal extent of snow distribution
                      if (tmp2 .lt. fl_crit)  then

                          ! areal extent of snow category
                          tmp3 = min(1.d0,ice(1,ic,mi)%snow_dist - dble(jj-1))  ! percentage of snow dist
                          snow_dist_new = snow_dist_new - tmp3
                          tmp3 = tmp3*ice(1,ic,mi)%af/ice(1,ic,mi)%snow_dist

                          ! record snow depth of flooded snow
                          tmp4 = tmp4 + sh_mean*sda_multiplier(jj)*tmp3

                          ! record area flooded
                          af_total = af_total + tmp3

                          ! record flooded depth
                          flooded = flooded + abs(tmp2)*tmp3

                      endif
                  endif
              enddo

              if (snow_dist_new .gt. dble(sd_num)) then
                  testvar = -1
              endif


              ! backup cat 1 ice, flood cat 1 ice, merge with cat 2 ice
              !(current flooded ice), copy cat 1 ice back from backup,
              ! re-adjust cat 1 ice area for loss to cat 2.
              if (flooded .gt. 0.d0) then

                  ! backup entire ice category
    		          ! -----------------------------------------------------------
                  ice_temp = ice(1,ic,mi)

                  ! grow snow layer to mean depth of flooded ice
                  ! ------------------------------------------------------
                  tmp4 = tmp4/af_total ! new snow height
                  tmp1 = tmp4/ice(1,ic,mi)%snow%depth ! height scale factor

                  ice(1,ic,mi)%snow%th = ice(1,ic,mi)%snow%th*tmp1 ! vector multiply
                  ice(1,ic,mi)%snow%depth = tmp4


                  ! flood category 1
		          ! -----------------------------------------------------------

                  ! find areal mean flood depth
                  flooded = flooded/af_total

                  ! start new layer top at old layer top - don't exclude any ice
                  new_layer_top = -1.
	              melt_flood = 0.d0
	              melt_drain = 0.d0
                  dhdt_cm_s = 0.d0   ! required in algal migration

                  did_flood = .true.
                  flood_1x = .false.

                  ! find new thickness/icedepth structure with th_new and id_new
                  r_depth = ice(1,ic,mi)%id(ice(1,ic,mi)%z-z_sk) + c_gl + s_gl

                  tmp1 = max(IceD/2.d0,ice(sc,ic,mi)%snow%d(1))
                  f_depth = flooded*ice(sc,ic,mi)%snow%d(1)/tmp1
                  r_depth = r_depth + f_depth

                  if (r_depth .lt. h_min) then
                      print *,'sc ic mi:',sc,ic,mi,'- Weird r_depth less than h_min during flooding'
                  endif

                  ! record ice volume change due to snow-ice gain
                  m(mi)%snow_growth = m(mi)%snow_growth + f_depth* &
                      af_total*cell_area*1.e6

                  ! Conservation of heat: mix snow and seawater
                  ! assuming mass/volume changes are negligable inside a layer initially...
                  !top_f = snow_mass/snowh1m/IceD  ! volume fraction of fresh ice in snow

                  top_f = ice(sc,ic,mi)%snow%d(1)/IceD ! snow/water volume ratio
                  bot_f = 1.-top_f ! water/snow volume ratio

                  s_mean = max(min_sal,f(mi)%s*bot_f)
                  d_mean = f(mi)%d*bot_f + IceD*top_f

                  ! snow heat is assumed to be just latent heat of melting + 0.02*t*d
                  heat_mean = cw*f(mi)%t*f(mi)%d*bot_f - &
                        (ice(sc,ic,mi)%snow%t(z_snow)*0.02 + Lf)*IceD*top_f
                  dz_mean = 1.

#include "sia2_env_temp_from_heat.inc.f90"

                  t_mean = min(s_mean*mu*1.01,t_mean)
                  !t_mean = ice(sc,ic,mi)%t(1)


!				  S_mean = S_flood
!			      T_mean = f(mi)%t

#include "sia2_env_brine_calc.inc.f90"

                  T_flood = t_mean
                  S_flood = s_mean
                  d_flood = d_mean
                  bs_flood = bs_mean
                  bd_flood = bd_mean
                  bv_flood = bv_mean
                  heat_flood = heat_mean

#include "sia2_env_ice_grid.inc.f90"

                  ! if surface melt, get rid of some ice
                  ! s_gl should never be positive, should rename to s_ml (melt layer)
			      new_layer_top = abs(s_gl)

#include "sia2_env_ice_remap.inc.f90"

                  ! record total ice growth
                  if (c_gl .gt. 0.d0) then
                       m(mi)%cong_growth = m(mi)%cong_growth + &
                         c_gl*af_total*cell_area*1.e6
                  else
                       m(mi)%melt_loss = m(mi)%melt_loss + &
                         c_gl*af_total*cell_area*1.e6
                  endif
                  if (s_gl .gt. 0.d0) then
                       m(mi)%cong_growth = m(mi)%cong_growth + &
                         s_gl*af_total*cell_area*1.e6
                  else
                       m(mi)%melt_loss = m(mi)%melt_loss + &
                         s_gl*af_total*cell_area*1.e6
                  endif

                  ! set category area of flooded area
                  ice(1,ic,mi)%af = af_total

                  ! flooding takes from snow layers
                  ! ---------------------------------------------------------
                  if (boundary_file .eq. 1) then
                      ice(sc,ic,mi)%sh_offset = ice(sc,ic,mi)%sh_offset + flooded
                  endif

                  r_depth = max(ice(1,ic,mi)%snow%depth - flooded, 0.) + snow_gl
                  ice(1,ic,mi)%sh_prev = r_depth

#include "sia2_env_snow_grid.inc.f90"

#include "sia2_env_snow_remap.inc.f90"

                  did_flood = .false.
                  flooded = 0.d0

                  ! perform salt and freshwater accounting
                  ! (fresh_net and dsdt_net should be set from beginning of timestep)
                  ! ---------------------------------------------------------
									call sia2_ice_mass(ice(sc,ic,mi),fresh_tmp,dsdt_tmp)
									! record fresh water flux (kg)
									fresh_tmp = (fresh_tmp-fresh_net)*af_total*cell_area*1.e6
									m(mi)%h2o_flux = m(mi)%h2o_flux - fresh_tmp
									! record salt flux (kg)
									dsdt_tmp = (dsdt_tmp-dsdt_net)*af_total*cell_area*1.e6
									m(mi)%salt_flux = m(mi)%salt_flux - dsdt_tmp


		          ! Merge newly flooded ice into category 2 (flooded category)
		          ! -----------------------------------------------------------
                  id_mean = ice(1,ic,mi)%id(ice(1,ic,mi)%z-z_sk) ! newly flooded ice depth
                  jjj = ic
                  do jj=ic,ida_n
                      if (id_mean .gt. ic_h_min(jj) .and. id_mean .le. ic_h_max(jj)) then
                          jjj = jj
                      endif
                  enddo
                  ! must convert category values that we are merging INTO to bulk formulations
                  if (ice(2,jjj,mi)%af .gt. 0.d0) then
					  do ii=1,ice(2,jjj,mi)%z
						  tmp1 = ice(2,jjj,mi)%bv(ii)*c_001
						  ice(2,jjj,mi)%no3(ii) = ice(2,jjj,mi)%no3(ii)*tmp1
						  ice(2,jjj,mi)%nh4(ii) = ice(2,jjj,mi)%nh4(ii)*tmp1
						  ice(2,jjj,mi)%po4(ii) = ice(2,jjj,mi)%po4(ii)*tmp1
						  ice(2,jjj,mi)%sioh4(ii) = ice(2,jjj,mi)%sioh4(ii)*tmp1
						  ice(2,jjj,mi)%smalg(ii) = ice(2,jjj,mi)%smalg(ii)*tmp1
						  ice(2,jjj,mi)%poc(ii) = ice(2,jjj,mi)%poc(ii)*tmp1
					  enddo
                  endif

                  call sia2_env_merge_ice(2,jjj,1,ic,mi,ice,pur,f,pl)

                  ! put flooded category back in brine concentrations format
                  do ii=1,ice(2,jjj,mi)%z
                      tmp1 = 1000.d0/ice(2,jjj,mi)%bv(ii)
                      ice(2,jjj,mi)%no3(ii) = ice(2,jjj,mi)%no3(ii)*tmp1
                      ice(2,jjj,mi)%nh4(ii) = ice(2,jjj,mi)%nh4(ii)*tmp1
                      ice(2,jjj,mi)%po4(ii) = ice(2,jjj,mi)%po4(ii)*tmp1
                      ice(2,jjj,mi)%sioh4(ii) = ice(2,jjj,mi)%sioh4(ii)*tmp1
                      ice(2,jjj,mi)%smalg(ii) = ice(2,jjj,mi)%smalg(ii)*tmp1
                      ice(2,jjj,mi)%poc(ii) = ice(2,jjj,mi)%poc(ii)*tmp1
                  enddo

                  ! re-instate solid ice
                  ! -----------------------------------------------------------
                  ice(1,ic,mi) = ice_temp

                  if ((ice(1,ic,mi)%af - af_total) .gt. af_min .and. snow_dist_new .gt. 0.d0) then

                      ice(1,ic,mi)%af = ice(1,ic,mi)%af - af_total
                      ice(1,ic,mi)%snow_dist = snow_dist_new
                      flooded = 0.d0
                      sk_z = ice(1,ic,mi)%z
                      int_z = sk_z - z_sk
                      sk_1 = int_z + 1

                      ! change snow layer to mean depth of flooded ice
                      ! ------------------------------------------------------
                      r_depth = 0.d0
                      tmp1 = ice(1,ic,mi)%snow_dist
                      do jj=1,sd_num
                          r_depth = r_depth + min(1.d0,tmp1)*sh_mean*sda_multiplier(jj)
                          tmp1 = max(0.d0,tmp1-1.d0)
                      enddo
                      r_depth = r_depth/ice(1,ic,mi)%snow_dist  ! new depth

                      tmp4 = r_depth/ice(1,ic,mi)%snow%depth    ! depth change multiplier

                      ! Calculate freeboard height of new snow depth
                      fbh_new = ice(1,ic,mi)%id(sk_z) - &
                          (snow_mass*tmp4 + ice_mass)/f(mi)%d ! m

                      ice(1,ic,mi)%snow%th = ice(1,ic,mi)%snow%th*tmp4 ! vector multiply
                      ice(1,ic,mi)%snow%depth = r_depth
                      ice(1,ic,mi)%sh_prev = r_depth

                      r_depth = ice(1,ic,mi)%snow%depth + snow_gl

                      ! remap snow back to correct resolution grid

#include "sia2_env_snow_grid.inc.f90"

#include "sia2_env_snow_remap.inc.f90"

                      ! reset snow_gl, since it has not been applied to both remainig
                      ! and outgoing ice (as opposed yo s_gl, c_gl, which were only applied
                      ! to outogin ice, and will be considered later in sia_env_grid_remap
                      snow_gl = 0.

                      if (ice(1,ic,mi)%snow%z .gt. 0.d0) then
                          snowh1m = ice(1,ic,mi)%snow%depth
                      else
                          snowh1m = 0.d0
                      endif

                  else

                      !print *,'Entire snow depth distribution flooded away ... is this for real?'
                      !print *,'sc ic mi: ',sc,ic,mi,snow_dist_new,ice(sc,ic,mi)%af,af_total
                      !print *,'is sh: ',ice(sc,ic,mi)%id(ice(sc,ic,mi)%z),sh_mean,ice(sc,ic,mi)%snow%depth
                      !print *,'sh_increment',sh_increment
                      !print *,'all snow depths: ',ice(:,:,mi)%snow%depth
                      !print *,'all afs: ',ice(:,:,mi)%af

#include "sia2_env_null_params.f90"

                  endif


          	  endif ! end if flooded

          else ! end if sc==1 && sda_n > 1

              ! determine if flooding in category 2 (flooded category), and
              ! re-flood if necessary
              ! -----------------------------------------------------------
              flooded = 0.d0

              ! Calculate freeboard height of mean snow depth
              fbh_new = ice(sc,ic,mi)%id(sk_z) - (snow_mass + ice_mass)/f(mi)%d ! m

              ! flood based on surface depression
              if (fbh_new .lt. fl_crit) then
                  flooded = abs(fbh_new)
                  fbh_new = fbh_new + flooded
                  if (boundary_file .eq. 1) then
                      total_flooded = total_flooded + flooded
                  endif
              endif

               ! flood based on wicking
!              if (ice(sc,ic,mi)%bv(1) .ge. bv_conv .and. fbh_new .lt. 0.055) then
!                  flooded = 0.055 - fbh_new
!                  ice(sc,ic,mi)%fbh = fbh_new  ! for now use this - don't adjust since it is wicking above fbh

!                  if (boundary_file .eq. 1) then
!                      total_flooded = total_flooded + flooded
!                 endif
!              endif

          endif ! end sc switch

          endif ! end of are all flooding pre-conditinons met?

          ! reset freeboard height for this category
          ice(sc,ic,mi)%fbh = fbh_new


