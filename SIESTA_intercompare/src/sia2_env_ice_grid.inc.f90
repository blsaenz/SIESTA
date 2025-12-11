! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
! START - FIND NEW ICE LAYER THICKNESSES
! ----------------------------------------------------------------------
! File: SIA2_env_layer_thickness.ice.f90
! This code finds new layer heights for a new total ice depth
!
! Requires:
! r_depth - the new depth to be assigned
!
! Results:
! th_new - contains new heights, up to th_z (all layers except skeletal layers)
! id_new - contains new ice depths for th_new
! int_z_new - new count of internal layers
! ----------------------------------------------------------------------

          keep_growing = .true.
          z_last = ice(sc,ic,mi)%z

          ! test for max_depth
          if ((r_depth .gt. h_max) .and. (.not. ignore_h_max)) then
              max_d_exceeded = max_d_exceeded + 1
              maxed_depth = .true.
              r_depth = h_max
              print *,'Max depth exceeded - sc ic mi:',sc,ic,mi
          else
              maxed_depth = .false.
          endif

          if (grid_model .eq. 0) then

              if (r_depth .lt. h_min) then
                  if ((r_depth+0.000001) .lt. h_min) then
                      print *,'Warning - minimum depth exceeded - unknown results may occur'
                      print *,'mi: ',mi,'in sia2_ice_grid line 26'
                  else
                      r_depth = h_min ! tollerate some floating point slop here, reset to h_min to avoid boundary creep
                  endif
              endif

              int_z_new = z_int_max

              ! divide up depth equally between layers
              th_new = r_depth/int_z_new

          elseif (grid_model .ge. 1) then

              if ((grid_model .eq. 2) .and. (r_depth .le. h_crit)) then
                  ! determine ice layers required if less than z_int_max
                  r_depth = r_depth - z_th_min*2.d0  ! top and bottom layers don't change thickness
                  ii = z_int_min-1
                  do while ((ii .le. z_int_max-2) .and. keep_growing)
                      if (r_depth .le. z_th_min*dble(ii)) then
                          ii = ii-1 ! take it back one layer so that we don't end up with layers less than z_th_min
                          int_z_new = ii+2 ! include top and bottom layers that are at z_th_min always...

                          ! update middle layers
                             do jj=2,ii+1
                              th_new(jj) = r_depth/dble(ii)
                          enddo
                          ! update top/bottom layers
                          th_new(1) = z_th_min
                          th_new(int_z_new) = z_th_min
                          keep_growing = .false.
                      endif
                      ii = ii+1
                  enddo
                  if (int_z_new .lt. z_int_max) then
                      do ii=int_z_new+1,z_int_max
                          th_new(ii) = 0.d0
                      enddo
                  endif
              else

                  if (r_depth .lt. h_min .and. (grid_model .eq. 1)) then
                      if ((r_depth+0.000001) .lt. h_min) then
                          print *,'Warning - minimum depth exceeded - unknown results may occur'
                          print *,'mi: ',mi,'in sia2_ice_grid line 26'
                      else
                          r_depth = h_min ! tollerate some floating point slop here, reset to h_min to avoid boundary creep
                      endif
                  endif

                  ! number of internal layers is z_int_max
                  int_z_new = z_int_max

                  ! remove minimum depth & find unit to be added to
                  r_depth = r_depth - dble(int_z_new)*z_th_min

                  jj = z_int_max/2
                  ! find divisor by which we will divide depth bins
                  ! geometric series sum (minus 2 b/c top and bottom are constant)
                  layer_divisor = 2.d0*(1.d0 - 1.3d0**jj)/(1.d0 - 1.3d0) - 2.d0

                  ! find new thicknesses
                  th_new(1) = z_th_min
                  th_new(int_z_new) = z_th_min
                  do ii=2,int_z_new/2
                      jj = ii-1
                      th_new(ii) = z_th_min + (1.3d0**jj)/layer_divisor*r_depth
                      th_new(int_z_new-ii+1) = th_new(ii)
                  enddo

              endif

          endif

           ! deal with skeletal thicknesses
           if (z_last .eq. 0) then
               ! new ice - create skeletal layer
               do ii=int_z_new+1,int_z_new+z_sk
                   th_new(ii) = sk_th_min
               enddo
           else
               ! re-grid existing skeletal layer to new position
               jj=z_last-z_sk+1
               do ii=int_z_new+1,int_z_new+z_sk
                   th_new(ii) = ice(sc,ic,mi)%th(jj)
                   jj=jj+1
               enddo
           endif

           ! update id_new
           ! -----------------------------------------------------------
           id_new(1) = th_new(1)
           do ii=2,int_z_new+z_sk
               id_new(ii) = id_new(ii-1) + th_new(ii)
           enddo

! ----------------------------------------------------------------------
! END - FIND NEW ICE LAYER THICKNESSES
! /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
