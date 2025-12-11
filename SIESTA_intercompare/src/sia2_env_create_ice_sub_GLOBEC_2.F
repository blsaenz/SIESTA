
      subroutine sia2_env_create_ice(m,ice,f,v_f,pl,mi)

          use sia2_globals
          implicit none

          ! function arguments
		  !----------------------------------------------------          
          type (meta_type) :: m(tcells)
          type (ice_type) :: ice(ida_n,tcells)
          type (forcing_type) :: f(tcells)
          type (validation_f_type) :: v_f
          type (platelet_type) :: pl
          integer :: ic,mi

          ! private variables
		  !----------------------------------------------------                    
          logical :: &
          keep_growing,maxed_depth

		  integer :: &
          z_last,ii,kk,jj,z_snow,int_z_new, sk_1,sk_z,int_z,max_d_exceeded

          double precision :: &
          t_mean,s_mean,d_mean,bs_mean,bd_mean,bv_mean,heat_mean, &
          alpha0,alpha1,alpha2,alpha3,t_melt,layer_divisor, &
          t_depth,snowh1m,f_snow,tmp1,airtemp1c,t_interface,r_depth, &
          flooded,old_layer_top,new_layer_bot,old_layer_bot,new_layer_top, &
          dz_total,heat_total,interim_top,interim_bot,z1,z2,dz,snowfall_d, &
          testvar,top_f,bot_f,T_flood,S_flood, &
          d_flood,bs_flood,bd_flood,bv_flood,heat_flood,melt_drain, &
          d_snow_new,heat_snow_new

          double precision, dimension(z_max) :: &
          t_new,th_new,d_new,id_new,melt_new

          ! start create ice
		  !----------------------------------------------------                    
          
		  airtemp1c = f(mi)%at - kelvin0

          ! update status to that of an active cell
          m(mi)%status = 1 
          m(mi)%pxl_h_offset = cell_side/2.
          m(mi)%pxl_v_offset = cell_side/2.


          if (start) then
              if (validation .eq. 1) then
                  t_depth = abs(v_f%iceh1/100.)
! ----------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! init lines for Grossi/McMurdo 1982 validation only
! ----------------------------------------------------------------------------------
!                      t_depth = 1.95
! ----------------------------------------------------------------------------------
! init lines for Grossi/McMurdo 1982 validation only
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! ----------------------------------------------------------------------------------

              else
				  if(f(mi)%ic_interp .gt. 0.9) then
					  t_depth = 1.28           ! max ice height, for icecon >= 80%

				  else
					  t_depth = h_min + (1.28-h_min)*(f(mi)%ic_interp-0.1)/0.8 ! varied height between min mand max height for 0.1 < icecon < 0.8
				  endif
              endif              
          else
              t_depth = h_min
          endif
          

          ! find new ice category
          ic = 1
          if (ida_n .gt. 1) then
              do while (ic_h_min(ic) .gt. t_depth) 
                  ic = ic + 1
              enddo
          endif

          ice(ic,mi)%af = f(mi)%ic_interp ! set default areal fraction
          m(mi)%a_new = m(mi)%a_new + ice(ic,mi)%af*cell_area  ! record new area

          ! heat test ----------------------------------------------- 
          ! t_depth = 0.32
          ! heat test ----------------------------------------------- 
		            
          ! initially, assume no infiltration, but zero freeboard
          ice(ic,mi)%fbh = 0.
                       
          if ((snow_model .eq. 1) .and. &
          ((start .eqv. .FALSE.) .or. (validation .eq. 1))) then
              ice(ic,mi)%sh_prev = 0.
  
! ----------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! testing for SHEBA snow
! ----------------------------------------------------------------------------------
!              ice(ic,mi)%sh_prev = f(mi)%sh_interp*0.01
! ----------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! testing for SHEBA snow
! ----------------------------------------------------------------------------------
 
              snowh1m = ice(ic,mi)%sh_prev
          else           
              if (f(mi)%sh_interp .eq. c9999) then
                  ice(ic,mi)%sh_prev = 0.
              else
                  ice(ic,mi)%sh_prev = f(mi)%sh_interp*0.01
              endif
              snowh1m = ice(ic,mi)%sh_prev
          endif

              ! test
              !snowh1m = 0.

          ! fractional snow cover
		  if (h_snowpatch .eq. 0.) then
		      f_snow = 1.
		  else
			  f_snow = snowh1m/(snowh1m + h_snowpatch)			  
		  endif
		  ! find mean ice surface temp from mean snowdepth
		  !t_interface = airtemp1c + snowh1m*(f(mi)%t - airtemp1c) &
		  !	      /(snowh1m + t_depth)
          tmp1 = ksnow*2.2/(2.2*snowh1m + ksnow*t_depth) ! gamma in Arrigo et al. 1993 appendix
          t_interface = (0.2456008*airtemp1c + tmp1*f(mi)%t)/(0.2456008+tmp1)
		  t_interface = min(f(mi)%t,t_interface)

! ---------------------------------------------------------------------------------
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! init lines for ISW Site B forcing validation only !!!!!!!!!
! comments these out for all other runs
! ----------------------------------------------------------------------------------


!		  t_interface = 2. !f(mi)%t

! ---------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! ----------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! init lines for ISW Site A forcing validation only !!!!!!!!!
! comments these out for all other runs
! ----------------------------------------------------------------------------------


!		  t_interface = f(mi)%t

! ---------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! ----------------------------------------------------------------------------------



! ---------------------------------------------------------------------------------
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! init lines for McMurdo 1982/Grossi forcing validation only !!!!!!!!!
! comments these out for all other runs
! ----------------------------------------------------------------------------------


!		  t_interface = -14.0

! ---------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! ----------------------------------------------------------------------------------



!		  if (airtemp1c .gt. 0.) then
		      snowfall_d = den_s_wet*1.e6
!		  else
!		      snowfall_d = den_s_dry*1.e6
!		  endif

          r_depth = snowh1m
          ice(ic,mi)%snow%z = 0 ! tells grid/remap routines to create new snow
          ! tells remap to use linear interpolation of temperature in ice/snow pack
          if (iif .gt. 0) then
              ice(ic,mi)%snow%ts = f(mi)%t  
          else
              ice(ic,mi)%snow%ts = t_interface  
          endif		  
      
#include "sia2_env_snow_grid.inc.f90"               

          flooded = 0.
          melt_drain = 0.

#include "sia2_env_snow_remap.inc.f90"               

             ! set melt status for snow cover to 1 - this is old wet snow
             do ii=1,ice(ic,mi)%snow%z
                 ice(ic,mi)%snow%d(ii) = den_s_wet*1.e6
                 ice(ic,mi)%snow%melt(ii) = 1.
             enddo
          
          r_depth = t_depth

#include "sia2_env_ice_grid.inc.f90"               
          
          sk_1 = int_z_new + 1
          sk_z = int_z_new + z_sk
          int_z = int_z_new
          ice(ic,mi)%z = sk_z

          ! setup flooding vars in case there is initial flooding specified
		  top_f = fl_ratio					  
		  bot_f = 1.-top_f ! remaining vol fraction taken up by seawater

		  T_flood = f(mi)%t
		  S_flood = f(mi)%s*bot_f
		  S_mean = S_flood
		  T_mean = f(mi)%t
		  
#include "sia2_env_brine_calc.inc.f90"

		  d_flood = d_mean
		  bs_flood = bs_mean
		  bd_flood = bd_mean
		  bv_flood = bv_mean
		  heat_flood = heat_mean

! ---------------------------------------------------------------------------------
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! init lines for ISW Site B forcing validation only !!!!!!!!!
! comments these out for all other runs
! ----------------------------------------------------------------------------------


!		  t_interface = f(mi)%t

! ---------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! ----------------------------------------------------------------------------------


          ! setup individual internal ice layers
          ! ----------------------------------------------------------- 
          do ii=1,sk_z

              ! initial ice does not have fatty brine channels
              ice(ic,mi)%bced(ii) = 0

			  ! update thickness & ice depth
			  ice(ic,mi)%th(ii) = th_new(ii)
			  ice(ic,mi)%id(ii) = id_new(ii)
 
               if (ii .lt. sk_1) then 
	
				  ! new ice layer salinity
                  if (validation .eq. 1) then
                      if (iis .eq. 2) then
						  ! interpolate default new-ice salinity curve to get initial salinity
						  call sia2_interpolate_salinity2( & 
							(ice(ic,mi)%id(ii)-ice(ic,mi)%th(ii)/2)/t_depth,ice(ic,mi)%s(ii))
                      elseif (iis .eq. 1) then
						  ! interpolate default new-ice salinity curve to get initial salinity
						  call sia2_interpolate_salinity1( & 
							(ice(ic,mi)%id(ii)-ice(ic,mi)%th(ii)/2)/t_depth,ice(ic,mi)%s(ii))
					  else
						  ! constant ice salinity
						  ice(ic,mi)%s(ii) = s_const
					  endif
				  elseif (start .or. iis .eq. 2) then
					  ! interpolate default new-ice salinity curve to get initial salinity
					  call sia2_interpolate_salinity2( & 
						(ice(ic,mi)%id(ii)-ice(ic,mi)%th(ii)/2)/t_depth,ice(ic,mi)%s(ii))
				  elseif (iis .eq. 1) then
					  ! interpolate default new-ice salinity curve to get initial salinity
					  call sia2_interpolate_salinity1( & 
						(ice(ic,mi)%id(ii)-ice(ic,mi)%th(ii)/2)/t_depth,ice(ic,mi)%s(ii))
				  else
					  ! constant ice salinity
					  ice(ic,mi)%s(ii) = s_const
                  endif
	
				  ! new ice layer temperature
                  if (iif .gt. 0.) then
 				      ice(ic,mi)%t(ii) = f(mi)%t                 
                  else
					  if (iit .eq. 1) then
						  ! isothermal at Tw
						 ice(ic,mi)%t(ii) = f(mi)%t
					  else
						  ! interpolate linearly between t_interface and water temp
						  ice(ic,mi)%t(ii) = t_interface - (t_interface-f(mi)%t) * (ice(ic,mi)%id(ii) - &
						  ice(ic,mi)%th(ii)/2)/t_depth
					  endif
                  endif

! ---------------------------------------------------------------------------------
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! init lines for ISW Site B forcing validation only !!!!!!!!!
! comments these out for all other runs
! ----------------------------------------------------------------------------------
                  
!                 if (ii .lt. 8) then
!                    ice(ic,mi)%t(ii) = -2. - (-2.-f(mi)%t) * (ice(ic,mi)%id(ii) - &
!					  ice(ic,mi)%th(ii)/2)/(0.07)
!                 endif

! ---------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! ----------------------------------------------------------------------------------
				  
              else
	
				  ice(ic,mi)%t(ii) = f(mi)%t
				  ice(ic,mi)%s(ii) = f(mi)%s*c_5

              endif              


! ---------------------------------------------------------------------------------
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! init lines for ISW Site B forcing validation only !!!!!!!!!
! comments these out for all other runs
! ----------------------------------------------------------------------------------

              ! use below for 122 layers
 !             if (ii .lt. 8) then
 !                 ! interpolate linearly between 6 and 16
 !                 ice(ic,mi)%s(ii) = 6. + 12. * (ice(ic,mi)%id(ii) - &
 !                 ice(ic,mi)%th(ii)/2)/0.08
                  !ice(ic,mi)%t(ii) = -2. - (-2.-f(mi)%t) * (ice(ic,mi)%id(ii) - &
                  !ice(ic,mi)%th(ii)/2)/0.08
 !             endif                  
 !             if (ii .ge. 8 .and. ii .le. 29) then
 !                 ice(ic,mi)%s(ii) = 18
 !                 ice(ic,mi)%bced(ii) = 1
 !             endif
 !             if (ii .ge. 30 .and. ii .le. 48) then
 !                 ice(ic,mi)%s(ii) = 8.
 !                 ice(ic,mi)%bced(ii) = 1
 !             endif


              ! use below for 42 layers
!              if (ii .lt. 12) then
                  ! interpolate linearly between 6 and 16
!                 ice(ic,mi)%s(ii) = 6. + 14. * (ice(ic,mi)%id(ii) - &
!                  ice(ic,mi)%th(ii)/2)/0.13
!              endif                  
!              if (ii .ge. 12 .and. ii .le. 17) then
!                  ice(ic,mi)%s(ii) = 20.
!                  ice(ic,mi)%bced(ii) = 1
!              endif
!              ice(ic,mi)%s(18) = 13.
! ---------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! ----------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! init lines for ISW Site A forcing validation only !!!!!!!!!
! comments these out for all other runs
! ----------------------------------------------------------------------------------
!              if (ii .lt. 19) then
!                  ice(ic,mi)%s(ii) = 16.
!                  ice(ic,mi)%bced(ii) = 1
!              endif                  
!              if (ii .ge. 19 .and. ii .le. 36) then
!                  ice(ic,mi)%s(ii) = 8.
!                  ice(ic,mi)%bced(ii) = 1
!              endif

! ---------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! ----------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! init lines for GLOBEC1 forcing only !!!!!!!!!
! comments these out for all other runs
! ----------------------------------------------------------------------------------
!              if (ii .lt. 20) then
!                  ice(ic,mi)%s(ii) = min(ice(ic,mi)%s(ii),4.)
!              endif                  
!              if (ii .ge. 19 .and. ii .le. 36) then
!                  ice(ic,mi)%s(ii) = 8.
!                  ice(ic,mi)%bced(ii) = 1
!              endif

! ---------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! ----------------------------------------------------------------------------------


			  ! set default brine paramaters
			  T_mean = ice(ic,mi)%t(ii)
			  S_mean = ice(ic,mi)%s(ii)

#include "sia2_env_brine_calc.inc.f90"
	
			  ! update brine-based physics
			  ice(ic,mi)%d(ii) = d_mean
			  ice(ic,mi)%bs(ii) = bs_mean
			  ice(ic,mi)%bd(ii) = bd_mean
			  ice(ic,mi)%bv(ii) = bv_mean
			  ice(ic,mi)%heat(ii) = heat_mean

              ! give minimum algae concentration
              ice(ic,mi)%smalg(ii) = alg_wc/(ice(ic,mi)%bv(ii)*c_001)  !concentrate algae into brine volume

! ---------------------------------------------------------------------------------
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! init lines for ISW Site B forcing validation only !!!!!!!!!
! comments these out for all other runs
! ----------------------------------------------------------------------------------


               ! use below for 122 layers
 !             if (ii .ge. 8 .and. ii .le. 48) then
 !                 ice(ic,mi)%smalg(ii) = 200./(ice(ic,mi)%bv(ii)*c_001)
 !             endif
                  
!                  ice(ic,mi)%smalg(ii) = 410.



               ! use below for 42 layers
!              if (ii .ge. 12 .and. ii .le. 17) then
!                  ice(ic,mi)%smalg(ii) = 300./(ice(ic,mi)%bv(ii)*c_001)
!              endif
!              ice(ic,mi)%smalg(18) = 150.

! ---------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! ----------------------------------------------------------------------------------



              ! start with zero POC
              ice(ic,mi)%poc(ii) = 0.

              if (start .and. validation .eq. 0) then
                  tmp1 = f(mi)%s/ice(ic,mi)%bs(ii)
              elseif (iin .eq. 1) then
                  ! scale nutrients linearly from 0->full from top->bottom, mimicking fall surface depletetion
                  tmp1 = ((ice(ic,mi)%id(ii)-ice(ic,mi)%th(ii)/2)/(t_depth+sk_th_min*z_sk))
              elseif (iin .eq. 2) then
                  tmp1 = ((ice(ic,mi)%id(ii)-ice(ic,mi)%th(ii)/2)/(t_depth+sk_th_min*z_sk))**2
              elseif (iin .eq. 3) then
                  tmp1 = n_f
              elseif (iin .eq. 4) then
                  ! same as seawater
                  tmp1 = f(mi)%s/ice(ic,mi)%bs(ii)
              else
                  tmp1 = 1.
              endif

! ---------------------------------------------------------------------------------
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! init lines for ISW Site B forcing validation only !!!!!!!!!
! comments these out for all other runs
! ----------------------------------------------------------------------------------

!              if (ii .lt. 49) then
!                  tmp1 = f(mi)%s/ice(ic,mi)%bs(ii)*0.2
!              else
!                  tmp1 = f(mi)%s/ice(ic,mi)%bs(ii)*0.3
!              endif



! ---------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! ----------------------------------------------------------------------------------


              ! intial nutrients 
              ice(ic,mi)%no3(ii)= f(mi)%no3*ice(ic,mi)%bs(ii)/f(mi)%s*tmp1  ! scale conservatively with salinity
              ice(ic,mi)%nh4(ii)= f(mi)%nh4*ice(ic,mi)%bs(ii)/f(mi)%s*tmp1  ! scale conservatively with salinity
              ice(ic,mi)%po4(ii)= f(mi)%po4*ice(ic,mi)%bs(ii)/f(mi)%s*tmp1  ! scale conservatively with salinity
              ice(ic,mi)%sioh4(ii)= f(mi)%sioh4*ice(ic,mi)%bs(ii)/f(mi)%s*tmp1  ! scale conservatively with salinity

           enddo

! ---------------------------------------------------------------------------------
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! init lines for ISW Site B forcing validation only !!!!!!!!!
! comments these out for all other runs
! ----------------------------------------------------------------------------------

             ! THIS WAS ONLY USED FOR TESTING AN HIGH_BRINE VOL FLUID TRANSFER METHOD
             ! do not use under normal circumstances
             
             ! ice(ic,mi)%no3(32:122)= f(mi)%no3
             ! ice(ic,mi)%nh4(32:122)= f(mi)%nh4
             ! ice(ic,mi)%po4(32:122)= f(mi)%po4
             ! ice(ic,mi)%sioh4(32:122)= f(mi)%sioh4

! ---------------------------------------------------------------------------------
! \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
! ----------------------------------------------------------------------------------


           if (use_pl .eq. 1) then
			   pl%t = f(mi)%t
			   pl%s = 0.8*f(mi)%s
			   pl%bs = f(mi)%s
			   pl%bd = 1000000. + pl%bs*800.
			   pl%bv = 800.
               do ii=1,pl_max
                   pl%no3(ii) = f(mi)%no3                   
                   pl%nh4(ii) = f(mi)%nh4                   
                   pl%po4(ii) = f(mi)%po4                   
                   pl%sioh4(ii) = f(mi)%sioh4                   
                   pl%smalg(ii) = min_alg/(pl%bv*c_001) 
                   pl%poc(ii) = 0.
               enddo
           endif



           ! set molecular sub-layer nutrient concentrations to seawater values
           ! ----------------------------------------------------------- 
           ice(ic,mi)%poc(ml_z) = f(mi)%poc
           ice(ic,mi)%no3(ml_z) = f(mi)%no3
           ice(ic,mi)%nh4(ml_z) = f(mi)%nh4
           ice(ic,mi)%po4(ml_z) = f(mi)%po4
           ice(ic,mi)%sioh4(ml_z) = f(mi)%sioh4

           !print *,'Forcing nutrients: ',f(mi)%no3,f(mi)%po4,f(mi)%sioh4

           ! record new ice volume as new to model domain
           m(mi)%md_gain = m(mi)%md_gain + ice(ic,mi)%id(int_z)* &
             ice(ic,mi)%af*cell_area*1.e6

!      if(ice(ic,mi)%s(1) .le. 0.12 .and. ((ice(ic,mi)%no3(1)*ice(ic,mi)%bv(1)*c_001) .gt. 10.)) then
!          testvar = -1
!          testvar = sqrt(testvar)
!      endif


      end subroutine sia2_env_create_ice