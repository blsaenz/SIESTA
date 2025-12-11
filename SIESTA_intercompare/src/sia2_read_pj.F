! file: SIA2_readpj.f90
! Sea Ice Algae Model 2 - Saenz & Arrigo
! Version beta
! ======================================================================
! ======================================================================


! ======================================================================
! Subroutine: sia2_readpj.F90
! Purpose: Reads in EASE projection grids, masks, and lat/lon data and prepares
! grid indexes
! ======================================================================


      SUBROUTINE sia2_read_pj(pj,v_f,stations_mdh,stations_mdv)

          use sia2_globals
          implicit none

          ! function arguments
      !----------------------------------------------------

          type (proj_type) :: pj
          type (boundary_f_type) :: v_f
          integer, dimension(n_stations) :: stations_mdh, stations_mdv

          ! internal variables
      !----------------------------------------------------
          integer :: i,j,ii,jj,sizeof,recl,i_md,j_md
          logical :: exists


      ! Load land mask from file
      !----------------------------------------------------
      if (hemi .eq. 0) then
        inquire(file='data/land_masks/EASE_25km_land_mask_721x721.txt',exist=exists)
        if (exists) then
          open(unit=22, file='data/land_masks/EASE_25km_land_mask_721x721.txt', form='FORMATTED')
          read(22,*) pj%mask
        else
          print *, 'err in sia2_read_pj: EASE 25km land mask grid not found'
          call exit(1010)
        endif
      else
        inquire(file='data/land_masks/EASE_25km_land_mask_721x721_NORTH.txt',exist=exists)
        if (exists) then
          open(unit=22, file='data/land_masks/EASE_25km_land_mask_721x721_NORTH.txt', form='FORMATTED')
          read(22,*) pj%mask
        else
          print *, 'err in sia2_read_pj: EASE 25km land mask grid not found'
          call exit(1010)
        endif
      endif
      ! Load latitude longitude grids for calculating solar elevation/PUR
      !----------------------------------------------------
      sizeof=8
      if (hemi .eq. 0) then
    !#ifdef Bigendian
    !      inquire(file='data/lat_lon_grids/EASE_25km_lat_invert.double',exist=exists)
    !      if (exists) then
    !        open(unit=23, file='data/lat_lon_grids/EASE_25km_lat_invert.double', form='UNFORMATTED',  &
    !             access='direct', recl=sizeof*720*721)
    !#else
          inquire(file='data/lat_lon_grids/EASE_25km_lat_invert.double.le',exist=exists)
          if (exists) then
            open(unit=23, file='data/lat_lon_grids/EASE_25km_lat_invert.double.le', form='UNFORMATTED',  &
                 access='direct', recl=sizeof*720*721)
    !#endif
            read(23,rec=1) pj%lat(:,1:720)
            ! Fill in blanks in projection grid left from above (don't ever use those cells anyway)
            do j= 1,grid_h
                pj%lat(j,grid_v) = -200.
            enddo
          else
            print *, 'err in sia2_read_pj: EASE 25km Latitude grid not found'
            call exit(1009)
          endif

    !#ifdef Bigendian
    !      inquire(file='data/lat_lon_grids/EASE_25km_lon_invert.double',exist=exists)
    !      if (exists) then
    !        open(unit=24, file='data/lat_lon_grids/EASE_25km_lon_invert.double', form='UNFORMATTED',  &
    !             access='direct', recl=sizeof*720*721)
    !#else
          inquire(file='data/lat_lon_grids/EASE_25km_lon_invert.double.le',exist=exists)
          if (exists) then
            open(unit=24, file='data/lat_lon_grids/EASE_25km_lon_invert.double.le', form='UNFORMATTED',  &
                 access='direct', recl=sizeof*720*721)
    !#endif
            read(24,rec=1) pj%lon(:,1:720)
            ! Fill in blanks in projection grid left from above (don't ever use those cells anyway)
            do j= 1,grid_h
                pj%lon(j,grid_v) = -200.
            enddo
          else
            print *, 'err in sia2_read_pj: EASE 25km Longitude grid not found'
            call exit(1008)
          endif

      else ! hemi==1, nothern hemisphere
          inquire(file='data/lat_lon_grids/EASE_25km_lat_invert_NORTH.double.le',exist=exists)
          if (exists) then
            open(unit=23, file='data/lat_lon_grids/EASE_25km_lat_invert_NORTH.double.le', form='UNFORMATTED',  &
                 access='direct', recl=sizeof*721*721)
            read(23,rec=1) pj%lat(:,1:721)
          else
            print *, 'err in sia2_read_pj: EASE 25km NORTH Latitude grid not found'
            call exit(1009)
          endif
          inquire(file='data/lat_lon_grids/EASE_25km_lon_invert_NORTH.double.le',exist=exists)
          if (exists) then
            open(unit=24, file='data/lat_lon_grids/EASE_25km_lon_invert_NORTH.double.le', form='UNFORMATTED',  &
                 access='direct', recl=sizeof*721*721)
            read(24,rec=1) pj%lon(:,1:721)
          else
            print *, 'err in sia2_read_pj: EASE 25km NORTH Longitude grid not found'
            call exit(1008)
          endif

      endif

      ! find index to NCEP/mercator-projection files
      !-----------------------------------------------------------------
      do j = 1,grid_v-1
        do i = 1,grid_h
                  if (pj%lon(i,j) .ne. -200. .and. pj%lat(i,j) .ne. -200.) then
            call gc_ll2ij(pj%lon(i,j),pj%lat(i,j),pj%x_mp(i,j),pj%y_mp(i,j))
            call ec_ll2ij(pj%lon(i,j),pj%lat(i,j),pj%x_ec(i,j),pj%y_ec(i,j))
            call ecint_ll2ij(pj%lon(i,j),pj%lat(i,j),pj%x_eci(i,j),pj%y_eci(i,j))
            call sia2_dgxy_from_ll(pj%lat(i,j),pj%lon(i,j),pj%x_dg(i,j),pj%y_dg(i,j))

                  else
                      pj%x_mp(i,j) = 1
                      pj%y_mp(i,j) = 1
                      pj%x_ec(i,j) = 1
                      pj%y_ec(i,j) = 1
                      pj%x_eci(i,j) = 1
                      pj%y_eci(i,j) = 1
                      pj%x_dg(i,j) = 1
                      pj%y_dg(i,j) = 1
                  endif
        enddo
       enddo

          close(22)
          close(23)
          close(24)

          ! find model domain indexes if using boundary_file
          if (boundary_file .eq. 1 .and. v_f%lat .ne. c9999 .and. v_f%lon .ne. c9999) then
              ! don't know where we are in the model domain, set these vars
              ! so sia2_ease25_from_ll knows to iterate over entire model domain
              v_f%mdh = 0
              v_f%mdv = 0
              call sia2_ease25_from_ll(pj,v_f%lat,v_f%lon,v_f%mdh,v_f%mdv)
        mdh1 = v_f%mdh
        mdv1 = v_f%mdv
        mdh2=mdh1
        mdv2=mdv1
          endif

          ! determine model domain
      !-----------------------------------------------------------------
          print *,'mdh:',mdh,'mdv:',mdv
          tcells = 0
          do j = 1,grid_v
              ! find model domain reference for current cell
              j_md = j - mdv1 + 1
              do i = 1,grid_h
                  ! find model domain reference for current cell
                  i_md = i - mdh1 + 1
                  if (stn_only .eq. 1) then
                      pj%mi(i,j) = 0
                      do ii=1,n_stations
                        if(i .eq. stations_mdh(ii) .and. j .eq. stations_mdv(ii)) then
                          tcells = tcells + 1
                          pj%mi(i,j) = tcells
                        endif
                      enddo
                  elseif (j_md .ge. 1 .and. j_md .le. mdv .and. i_md .ge. 1 &
                    .and. i_md .le. mdh) then
                      pj%mdh(i,j) = i_md
                      pj%mdv(i,j) = j_md
                     if (boundary_file .eq. 1) then
                      tcells = tcells + 1
                      ! map model memory to grid
                      pj%mi(i,j) = tcells
                    elseif (pj%mask(i,j) .eq. 1.) then !and. (pj%lat(i,j) .lt. -50.)) then
                      ! find total modeled variable memory required
                      tcells = tcells + 1
                      ! map model memory to grid
                      pj%mi(i,j) = tcells
                      !print *,'yup!','mdh:',mdh,'mdv:',mdv,i,j,i_md,j_md,tcells
                      !print *,'yup!',pj%lat(i,j),pj%lon(i,j)
                      !print *, pj%x_mp(i,j), pj%y_mp(i,j)
                      !print *, pj%x_ec(i,j), pj%y_ec(i,j)
                      !print *, pj%x_eci(i,j), pj%y_eci(i,j)
                      !print *, pj%x_dg(i,j), pj%y_dg(i,j)
                    else
                      pj%mi(i,j) = 0
                      !print *,'nope!','mdh:',mdh,'mdv:',mdv,i,j,i_md,j_md
                      !print *,'nope!',pj%lat(i,j),pj%lon(i,j)
                    endif
                   else
                      pj%mdh(i,j) = -200
                      pj%mdv(i,j) = -200
                   endif
              enddo
          enddo

      end SUBROUTINE sia2_read_pj


      SUBROUTINE sia2_ease25_from_ll(pj,lat,lon,ease_h,ease_v)

          use sia2_globals
          implicit none

          ! function arguments
      !----------------------------------------------------
          type (proj_type) :: pj
          double precision :: lat,lon
          integer :: ease_h,ease_v

          ! internal variables
      !----------------------------------------------------
          integer i,j,diff_i,diff_j
          double precision :: diff,min_diff


          if ((ease_h .eq. 0) .or. (ease_v .eq. 0)) then
             call sia2_ease25_from_ll_brute(pj,lat,lon,diff_i,diff_j)
          else
              diff_i = ease_h
              diff_j = ease_v
              ! check nearby cells to see if boundary_file station has moved
              min_diff = abs(lat-pj%lat(ease_h,ease_v))+abs(lon-pj%lon(ease_h,ease_v))
              if (min_diff < 1.0) then
                  call sia2_ease25_from_ll_brute(pj,lat,lon,diff_i,diff_j)
              else
                  do j=-1,1
                      do i=-1,1
                         diff =  abs(lat-pj%lat(ease_h+i,ease_v+j))+abs(lon-pj%lon(ease_h+i,ease_v+j))
                         if (diff .lt. min_diff) then
                             min_diff = diff
                             diff_i = ease_h+i
                             diff_j = ease_v+j
                          endif
                      enddo
                  enddo
              endif

         endif

         ease_h = diff_i
         ease_v = diff_j

      end SUBROUTINE sia2_ease25_from_ll

      SUBROUTINE sia2_ease25_from_ll_brute(pj,lat,lon,ease_h,ease_v)

          use sia2_globals
          implicit none

          ! function arguments
      !----------------------------------------------------
          type (proj_type) :: pj
          double precision :: lat,lon
          integer :: ease_h,ease_v

          ! internal variables
      !----------------------------------------------------
          integer i,j,diff_i,diff_j
          double precision :: diff,min_diff

         ! interate over whole damn thing
         min_diff = 1.e6
         do j=1,721
             do i=1,721
                 diff = abs(lat-pj%lat(i,j))+abs(lon-pj%lon(i,j))
                 if (diff .lt. min_diff) then
                     min_diff = diff
                     diff_i = i
                     diff_j = j
                  endif
              enddo
          enddo
         ease_h = diff_i
         ease_v = diff_j

      end SUBROUTINE sia2_ease25_from_ll_brute





      SUBROUTINE sia2_dgxy_from_ll(lat,lon,x_dg,y_dg)
          double precision :: lat,lon
          integer :: x_dg,y_dg

      x_dg = int(lon) + 1        ! 1 degree grid projection x (used by WOA data)
      y_dg = int(lat+90.) + 1    ! 1 degree grid projection y (used by WOA data)

      end SUBROUTINE sia2_dgxy_from_ll

