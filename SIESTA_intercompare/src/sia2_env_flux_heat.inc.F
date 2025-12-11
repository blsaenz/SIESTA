
          ! get number of snow layers
		  z_snow = ice(sc,ic,mi)%snow%z

         ! find albedo(s)
		  if (ice(sc,ic,mi)%snow%ts .gt. -1.) then
			  ! scale albedo by difference between wet and dry ice between 0 and -1 degC
!			  albedo_i = alb_i_wet + (alb_i_dry - alb_i_wet)*min(1.,abs(ice(sc,ic,mi)%snow%ts))
!			  albedo_s = alb_s_wet + (alb_s_dry - alb_s_wet)*min(1.,abs(ice(sc,ic,mi)%snow%ts))
        albedo_i = alb_i_wet
        albedo_s = alb_s_wet
		  else
			  albedo_i = alb_i_dry
			  albedo_s = alb_s_dry
		  endif

		  !if (ice(sc,ic,mi)%bv(1) .gt. 150.) then
		  !    albedo_i = 0.3
		  !endif

		  if (ice(sc,ic,mi)%id(int_z) .lt. c_5) then
			  ! scale using atan function toward open water albedo as thickness approaches 0
			  albedo_i = albedo_i - (albedo_i-0.06)* &        ! 0.06 = open water albedo
				  (1.-atan(atan_c_i*ice(sc,ic,mi)%id(int_z)))
		  endif

		  if (z_snow .eq. 0) then
              ! bare ice albedo
			  albedo = albedo_i
              ! set longwave emmissivity to that of ice
			  eps0 = eps_ice
          else
              ! set albedo to snow (weighted value)
			  albedo = f_snow*albedo_s + (1.-f_snow)*albedo_i
              ! set longwave emmissivity to that of snow
			  eps0 = eps_snow
		  endif

		  ! surface shortwave (needs albedo from above)
          if (f(mi)%swd .eq. c9999) then
			  ! Fr = e(sc,ic,mi)%Ed0*par_to_swd*(1.-albedo) !-- old crap - duplicated in Ed_W variable
			  Fr = 0.
          else
              Fr = f(mi)%swd*(1.-albedo)
			  ! Fr = e(sc,ic,mi)%Ed0*par_to_swd*(1.-albedo) ! used
              ! must zero internal shortwave heating, otherwise it would be duplicate
              Ed_W_temp = 0. ! vector assignment
              Ed_W_snow = 0. ! vector assignment
          endif
		  if (Fr .lt. 0.) then
			  Fr = 0.
	      endif

		  ! vars below are pre-calculated to reduce computation in iteration
		  F0_constants = Fr + FL + Fe_sub2 + Fs_sub
		  F0_4 = -1.*(eps0*steph_boltz+Fe_sub1*fe_A)
		  F0_3 = -1.*Fe_sub1*fe_B
		  F0_2 = -1.*Fe_sub1*fe_C
		  F0_1 = -1.*(Fe_sub1*fe_D + row_air*cp*Ch*v10)

          ! debug temperature problem
          tmp1 = 0.
          do ii=1,sk_z
             if ((ice(sc,ic,mi)%t(ii) .gt. min_t) .or. (ice(sc,ic,mi)%s(ii) .lt. min_sal) &
             .or. (s_new(ii) .lt. min_sal)) then
                 if ((ice(sc,ic,mi)%t(ii) .gt. min_t) .and. (tmp1 .eq. 0.)) then
					 print '(a15,i5,a8,i3,a8,i3)','Invalid T - mi=',mi,' grid_h=', &
					 	m(mi)%grid_h,' grid_v=',m(mi)%grid_v
					 !print *,'s=',ice(sc,ic,mi)%s
					 print '(i3,a4,i3,a6E14.7,a7,E14.7)',ii,' of ',sk_z,' -- t=', &
					 	ice(sc,ic,mi)%t(ii),' s_new=',s_new(ii)
                     ice(sc,ic,mi)%t(ii) = min_t
					 tmp1 = 1.
                 elseif((ice(sc,ic,mi)%s(ii) .lt. min_sal) .and. (tmp1 .eq. 0.)) then
					 !print '(a15,i5,a8,i3,a8,i3)','Invalid S - mi=',mi,' grid_h=', &
					 !	m(mi)%grid_h,' grid_v=',m(mi)%grid_v
					 !print '(i3,a4,i3,a6,E14.7,a3,E14.7,a7,E14.7)',ii,' of ',sk_z, &
					 !	' -- t=',ice(sc,ic,mi)%t(ii),' s=',ice(sc,ic,mi)%s(ii),' s_new=',s_new(ii)
					 !print *,'s=',ice(sc,ic,mi)%s
					 !print *,'t=',ice(sc,ic,mi)%t
					 !print *,'s_new=',s_new
                     ice(sc,ic,mi)%s(ii) = min_sal
					 tmp1 = 1.
                 elseif((s_new(ii) .lt. min_sal) .and. (tmp1 .eq. 0.)) then
					 !print '(a19,i5,a8,i3,a8,i3)','Invalid s_new - mi=',mi,' grid_h=',m(mi)%grid_h,' grid_v=',m(mi)%grid_v
					 !print *,'s=',ice(sc,ic,mi)%s
					 !print *,'t=',ice(sc,ic,mi)%t
					 !print *,'s_new=',s_new
					 !print '(i3,a4,i3,a6,E14.7,a7,E14.7)',ii,' of ',sk_z,' -- t=',ice(sc,ic,mi)%t(ii),' s_new=',s_new(ii)
                     s_new(ii) = min_sal
					 tmp1 = 1.
                 endif
             endif
          enddo

          ! prepare for heat flux convergence iteration loop
		  ! ----------------------------------------------------------------

		  ! remember current temps
		  Ts_last = ice(sc,ic,mi)%snow%ts
          if (z_snow .gt. 0) then
              jj=1
              do ii = z_snow,1,-1
                  T_last(jj) = ice(sc,ic,mi)%snow%t(ii)
                  jj=jj+1
              enddo
          endif
		  T_last(ice_1:ice_z) = ice(sc,ic,mi)%t(1:int_z)

		  ! seed new temps with old
          Ts_next = Ts_last
		  T_next = T_last

		  ! set iteration counter to zero
		  jjj = 0;

		  ! set iteration tolerance to something ridiculously large to enter the loop
		  T_diff = 1.e8

          if ((dt_step .eq. 6932) .or. (ice(sc,ic,mi)%t(1) .lt. (-40.))) then
              testvar = -1
          endif

          ! BEGIN CONVERGENCE
          ! ----------------------------------------------------------------
		  do while((abs(T_diff) .gt. nr_tol) .and. (jjj .lt. max_it))

			  ! convert surface temp to kelvin
			  Ts_kelvin = Ts_next + kelvin0

			  ! COMPUTE SURFACE/TOP LAYER MATRIX COEFFICIENTS
			  ! ----------------------------------------------------
              if (atmo .eq. 1) then

                  ! incoming longwave (FL) is already calculated at start of ice model

                  ! outgoing longwave
                  F_little_L = -1*eps0*steph_boltz*ts_kelvin**4 ! W/m^2/K^4 * K^4 = W/m^2

                  !call sia2_env_atmo(Ts_last,potT,v10,2.d0,shum,row_air,lhcoef,shcoef)
                  call sia2_env_atmo_const(Ts_next,v10,row_air,lhcoef,shcoef)

                  F_s = shcoef*(potT - Ts_kelvin)

                  Fe = lhcoef*(shum - (qq1/row_air)*exp(-qq2/Ts_kelvin))

                  F0 = Fr + FL + F_little_l + F_s + Fe

                  dF0 = -1.d0*( shcoef +  &  ! term from sensible turbulent flux (F_s)
                                4.d0*eps0*steph_boltz*Ts_kelvin**3 +  & ! term from outgoing longwave (F_little_l)
                                (lhcoef*qq1*qq2)/(row_air*exp(qq2/Ts_kelvin)*Ts_kelvin**2) &  ! term from turbulent latent flux (Fe)
                              )

              else

                  ! total surface heat flux minus conductive part
                  F0 = F0_constants + F0_4*Ts_kelvin**4 + F0_3*Ts_kelvin**3 + &
                    F0_2*Ts_kelvin**2 + F0_1*Ts_kelvin

                  ! slope of above equation
                  dF0 = 4.*F0_4*(Ts_kelvin**3) + 3.*F0_3*(Ts_kelvin**2) + &
                    2.*F0_2*Ts_kelvin + F0_1

                  ! debug vars ---------------
                  Fm0 = Fr + FL + Fe_sub2 - (eps0*steph_boltz+Fe_sub1*fe_A)*ts_kelvin**4 - &
                          Fe_sub1*fe_B*ts_kelvin**3 - Fe_sub1*fe_C*ts_kelvin**2 - (Fe_sub1*fe_D + &
                          row_air*cp*Ch*v10)*ts_kelvin + Fs_sub

                  F_little_L = -1*eps0*steph_boltz*ts_kelvin**4 ! W/m^2/K^4 * K^4 = W/m^2

                  Fe = Fe_sub1*(fe_A*(rhum* &
                      f(mi)%at**4-ts_kelvin**4)+fe_B*(rhum* &
                      f(mi)%at**3-ts_kelvin**3)+fe_C*(rhum* &
                      f(mi)%at**2-ts_kelvin**2)+fe_D*(rhum* &
                      f(mi)%at-ts_kelvin)+fe_E*(rhum-1)) ! W/m^2

                  F_s = row_air*cp*Ch*v10*(f(mi)%at-ts_kelvin)  ! kg/m^3 * J/kg/K * (dimensionless) * m/s * K = J/s/m^2 = W/m^2
                  ! debug vars ---------------

              endif

!              print *,'iii: ',iii,'------------------------------------------'
!              print *,'shcoef: ',shcoef
!              print *,'lhcoef: ',lhcoef
!              print *,'F_s: ',F_s
!              print *,'Fe: ',Fe
!              print *,'F_s: ',row_air*cp*Ch*v10*(f(mi)%at-ts_kelvin)
!              print *,'Fe: ',Fe_sub1*(fe_A*(rhum* &
!                      f(mi)%at**4-ts_kelvin**4)+fe_B*(rhum* &
!                      f(mi)%at**3-ts_kelvin**3)+fe_C*(rhum* &
!                      f(mi)%at**2-ts_kelvin**2)+fe_D*(rhum* &
!                      f(mi)%at-ts_kelvin)+fe_E*(rhum-1)) ! W/m^2
!
!              print *,'F0: ',F0
!              print *,'dF0: ',dF0
!              print *,'potT: ',potT
!              print *,'shum: ',shum
!              print *,'p: ',f(mi)%p
!              print *,'rhum: ',rhum
!              print *,'at: ',f(mi)%at


#include "sia2_env_flux_heat_solution7.inc.f90"

             ! increment interation counter
             jjj = jjj+1

			 ! END CONVERGENCE ITERATION
			 ! -------------------------------------------------------------
		  enddo

          ! test to see if the integration was valid, or whether it needs to
          ! be repeated with a smaller timestep
		  if ((T_step .gt. dt_sub_1_tol) .and. (dtt_s .eq. dtt_s_1)) then
			  integr_valid = .false.
		  elseif ((T_step .gt. dt_sub_2_tol) .and. (dtt_s .eq. dtt_s_2)) then
			  integr_valid = .false.
		  endif

