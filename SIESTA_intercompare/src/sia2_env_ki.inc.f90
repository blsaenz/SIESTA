
          ! get number of snow layers
		  z_snow = ice(sc,ic,mi)%snow%z

          ! top ice layer conductivity
          s_mean = s_new(1) !(ice(sc,ic,mi)%s(1)+s_new(1))*c_5
          d_mean = d_new(1) !(ice(sc,ic,mi)%d(1)+d_new(1))*c_5
          ki1 = (d_mean/IceD)*(2.11-0.011*ice(sc,ic,mi)%t(1) + &
            0.09*s_mean/ice(sc,ic,mi)%t(1) - (d_mean-IceD)/1.e6)
          ki1 = max(ki1,ki_min)

		  if (z_snow .eq. 0) then

			  ! indices for heat matrices
			  ice_1 = 1
			  ice_z = int_z

			  ! find intra-layer thickness
			  dth(1) = ice(sc,ic,mi)%th(1)*c_5

              ! top layer conductivity
 			  ki(1) = ki1
			  ki_down = ki1

	      else

			  ! indices for heat matrices
			  ice_1 = z_snow+1
			  ice_z = z_snow+int_z

			  ! find intra-layer thickness
			  dth(1) = ice(sc,ic,mi)%snow%th(z_snow)*c_5
			  dth(ice_1) = (ice(sc,ic,mi)%snow%th(1) + ice(sc,ic,mi)%th(1))*c_5

              ! top layer conductivity
			  d_mean = 1.d-6*ice(sc,ic,mi)%snow%d(z_snow)
			  ki_down = 0.138d0 - 1.01d0*d_mean + 3.233d0*d_mean**2 + 0.1495d0  ! sturm 1997 + 0.1495
              !ki_down = ksnow
			  ki(1) = ki_down

	          ! remaning snow intra-layer conductivities and thicknesses
	          if (z_snow .gt. 1) then
				  jj = 2
				  do ii=z_snow-1,1,-1
				      ! snow conductivity
				      ki_up = ki_down
                      d_mean = 1.d-6*ice(sc,ic,mi)%snow%d(ii)
                      ki_down = 0.138d0 - 1.01d0*d_mean + 3.233d0*d_mean**2 + 0.1495d0  ! sturm 1997 + 0.1495
                      !ki_down = ksnow
                      ! equate fluxes to find intra-layer ki
                      ki(jj) = ki_up*ki_down/ &
                        (ice(sc,ic,mi)%snow%th(ii+1)*ki_down + &
                        ice(sc,ic,mi)%snow%th(ii)*ki_up)* &
                        (ice(sc,ic,mi)%snow%th(ii+1)+ice(sc,ic,mi)%snow%th(ii))
                      ! intra-layer snow thickness
					  dth(jj) = (ice(sc,ic,mi)%snow%th(ii) + ice(sc,ic,mi)%snow%th(ii+1))*c_5
					  jj = jj+1
				  enddo
			  endif

			  ! equate fluxes to find intra-ice/snow ki
			  ki(ice_1) =  ki1*ki_down/(ice(sc,ic,mi)%snow%th(1)*ki1 + &
				ice(sc,ic,mi)%th(1)*ki_down)*(ice(sc,ic,mi)%snow%th(1)+ice(sc,ic,mi)%th(1))

              !

		  endif

		  ! find inter-layer ice thicknesses, conductivities
		  do ii=2,int_z

              ! mid-point to mid-point layer thickness
			  dth(ii+z_snow) = (ice(sc,ic,mi)%th(ii-1) + ice(sc,ic,mi)%th(ii))*c_5

              ! intra-layer conductivity
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

