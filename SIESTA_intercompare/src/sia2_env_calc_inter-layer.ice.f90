          ! get number of snow layers
		  z_snow = ice(sc,ic,mi)%snow%z

          ! top ice layer conductivity
          s_mean = s_new(1) !(ice(sc,ic,mi)%s(1)+s_new(1))*c_5
          d_mean = d_new(1) !(ice(sc,ic,mi)%d(1)+d_new(1))*c_5
		  ki_down = (d_mean/IceD)*(2.11-0.011*ice(sc,ic,mi)%t(1) + &
			0.09*s_mean/ice(sc,ic,mi)%t(1) - (d_mean-IceD)/1.e6)
	      ki_down = max(ki_down,ki_min)

		  if (z_snow .eq. 0) then
 			  ! determines correct matrix order
			  th(1) = ice(sc,ic,mi)%th(1)
			  dth(1) = th(1)*c_5
			  ! indices for heat matrices
			  ice_1 = 1
			  ice_z = int_z
              ! conductivity
 			  ki(1) = ki_down
	      else
			  ! find inter-layer thickness
			  dth(1) = ice(sc,ic,mi)%snow%th(z_snow)*c_5
			  dth(z_snow+1) = (ice(sc,ic,mi)%snow%th(1) + ice(sc,ic,mi)%th(1))*c_5
			  ki(1) = ksnow
			  if (z_snow .gt. 1) then
				  jj = 2
				  do ii=z_snow,2,-1
                      ki(ii) = ks
					  dth(jj) = (ice(sc,ic,mi)%snow%th(ii) + ice(sc,ic,mi)%snow%th(ii-1))*c_5
					  jj = jj+1
				  enddo
			  endif
			  ice_1 = z_snow+1
			  ice_z = z_snow+int_z
			  ! equate fluxes to find inter-ice/snow ki
			  ki(ice_1) =  ks*ki_down/(ice(sc,ic,mi)%snow%th(1)*ki_down + &
				ice(sc,ic,mi)%th(1)*ks)*(ice(sc,ic,mi)%snow%th(1)+ice(sc,ic,mi)%th(1))
		  endif

		  ! find inter-layer thicknesses, conductivities
		  do ii=2,int_z
              ! mid-point to mid-point layer thickness
			  dth(ii+z_snow) = (ice(sc,ic,mi)%th(ii-1) + ice(sc,ic,mi)%th(ii))*c_5

              ! inter-layer conductivity
			  ki_up = ki_down
			  s_mean = s_new(ii) !(ice(sc,ic,mi)%s(ii)+s_new(ii))*c_5
			  d_mean = d_new(ii) !(ice(sc,ic,mi)%d(ii)+d_new(ii))*c_5
			  ki_down = (d_mean/IceD)*(2.11-0.011*ice(sc,ic,mi)%t(ii) + &
				0.09*s_mean/ice(sc,ic,mi)%t(ii) - (d_mean-IceD)/1.e6)
	          ki_down = max(ki_down,ki_min)

			  ! equate fluxes to find inter-layer ki
			  ki(ii+z_snow) = ki_up*ki_down/(ice(sc,ic,mi)%th(ii-1)*ki_down + &
				ice(sc,ic,mi)%th(ii)*ki_up)*(ice(sc,ic,mi)%th(ii-1)+ice(sc,ic,mi)%th(ii))


		  enddo

		  ! last thickness, conductivity to water
		  dth(ice_z+1) = ice(sc,ic,mi)%th(int_z)*c_5
		  ki(ice_z+1) = ki_down

