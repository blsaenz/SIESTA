! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! START - FIND NEW SNOW LAYER THICKNESSES
! ----------------------------------------------------------------------
! File: SIA2_env_snow_thickness.ice.f90
! This code finds new layer heights the snow layers, and also adjusts
! the boundaries while conserving heat
!
! Requires:
! r_depth - the new depth to be assigned
! ice(sc,ic,mi)%snow%z - the old number of layers
! ice(sc,ic,mi)%snow%th - old layer heights
! ice(sc,ic,mi)%snow%t - old snow layer temps
!
! Results:
! ice(sc,ic,mi)%snow%z - new number of snow layers
! th_new - new layer heights
! ----------------------------------------------------------------------


		! Find new number of snow layers (m(mi)%z_snow and layer thicknesses
		! (th_new) based on current snow depth
		! ------------------------------------------------------------------
		 keep_growing = .true.
		 z_last = ice(sc,ic,mi)%snow%z
         ice(sc,ic,mi)%snow%depth = r_depth

		 if (r_depth .gt. snow_min .and. r_depth .gt. z_th_min) then

			 if (grid_model .eq. 0) then
				 ! determine snow layers required if less than z_snow_max
				 ! and allocate snow depth equally
				 if (r_depth .le. dble(z_snow_max)*z_th_min) then
					 ii = 1
					 do while ((ii .le. z_snow_max) .and. keep_growing)
						 if (r_depth .le. z_th_min*dble(ii)) then
							 ice(sc,ic,mi)%snow%z = ii
							 th_new = r_depth/dble(ii)
							 keep_growing = .false.
						 endif
						 ii = ii+1
					 enddo

				 ! allocate snow depth equally
				 else
					 ice(sc,ic,mi)%snow%z = z_snow_max
					 th_new = r_depth/dble(z_snow_max) ! vector/scalar opertion
				 endif

		     elseif (grid_model .ge. 1) then

				 ! determine snow layers required if less than z_snow_max
				 ! and allocate snow depth equally
				 if (r_depth .le. dble(z_snow_max)*z_th_min) then
					 ii = 1
					 do while ((ii .le. z_snow_max) .and. keep_growing)
						 if (r_depth .le. z_th_min*dble(ii)) then
							 ice(sc,ic,mi)%snow%z = ii
							 do jj=1,ii
								 th_new(jj) = r_depth/dble(ii)
							 enddo
							 keep_growing = .false.
						 endif
						 ii = ii+1
					 enddo

				 ! allocate snow depth to layers increasing in size up to
				 ! 130% of the previous lower layer
				 else
					  ice(sc,ic,mi)%snow%z = z_snow_max

					  ! remove minimum depth & find unit to be added to
					  r_depth = r_depth - z_snow_max*z_th_min

					  jj = z_snow_max/2
					  ! find divisor by which we will divide depth bins
					  ! geometric series sum (minus 2 b/c top and bottom are constant)
					  layer_divisor = 2.d0*(1.d0 - 1.3d0**jj)/(1.d0 - 1.3d0) - 2.d0

					  ! find new thicknesses
					  th_new(1) = z_th_min
					  th_new(z_snow_max) = z_th_min
					  do ii=2,z_snow_max/2
						  jj = ii-1
						  th_new(ii) = z_th_min + (1.3d0**jj)/layer_divisor*r_depth
						  th_new(z_snow_max-ii+1) = th_new(ii)
					  enddo
				 endif

             endif ! end of grid_model test

		 else ! end of r_depth > snow_min test
		    ice(sc,ic,mi)%snow%z = 0 ! no snow layers
		    ice(sc,ic,mi)%snow%depth = 0.d0

            ! vector assignments below
            ice(sc,ic,mi)%snow%t = 0.d0
            ice(sc,ic,mi)%snow%d = 0.d0
            ice(sc,ic,mi)%snow%heat = 0.d0
            ice(sc,ic,mi)%snow%th = 0.d0
            ice(sc,ic,mi)%snow%melt = 0.d0

		    ! save r_depth into sh_prev so mass is not wiped out
		    ice(sc,ic,mi)%sh_prev = r_depth

		 endif


! ----------------------------------------------------------------------
! END - FIND NEW SNOW LAYER THICKNESSES
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
