		! Find new number of snow layers (m(mi)%z_snow and layer thicknesses
		! (th_new) based on current snow depth
		! ------------------------------------------------------------------
		 keep_growing = .true.
		 z_last = ice(sc,ic,mi)%snow%z

         if (grid_model .eq. 0) then
             ! determine snow layers required if less than z_snow_max
             ! and allocate snow depth equally
             if (r_depth .le. dble(z_snow_max)*z_th_min) then
                 ii = 1
                 do while ((ii .le. z_snow_max) .and. keep_growing)
                     if (r_depth .le. z_th_min*dble(ii)) then
                         z_snow_new = ii
                         th_new = r_depth/dble(ii)
                         keep_growing = .false.
                     endif
                     ii = ii+1
                 enddo

             ! allocate snow depth equally
             else
                 z_snow_new = z_snow_max
                 th_new = r_depth/dble(z_snow_max) ! vector/scalar opertion
             endif

         elseif (grid_model .ge. 1) then

             ! determine snow layers required if less than z_snow_max
             ! and allocate snow depth equally
             if (r_depth .le. dble(z_snow_max)*z_th_min) then
                 ii = 1
                 do while ((ii .le. z_snow_max) .and. keep_growing)
                     if (r_depth .le. z_th_min*dble(ii)) then
                         z_snow_new = ii
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
                  z_snow_new = z_snow_max

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

