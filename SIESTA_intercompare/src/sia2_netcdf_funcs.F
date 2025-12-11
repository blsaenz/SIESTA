!=======================================================================
! Subroutine: sia2_forcing_load_stations
! Purpose: loads boundary data from the boundary_file netcdf page
!=======================================================================
     SUBROUTINE sia2_forcing_load_v_f(v_f,year,hour)

          use sia2_globals
          implicit none

          ! function arguments
          !----------------------------------------------------
          type (boundary_f_type) :: v_f
          integer :: year,hour
          double precision :: dtime

#include "netcdf.inc"
#include "gc_cdf.h"

          ! internal variables
          !----------------------------------------------------
          integer :: nc_index
          character (LEN=25) :: cdffile


          cdffile = 'data/Boundary/boundary.nc'
          rc = nf_open(cdffile,nf_nowrite,ncid)
          if (rc .ne. NF_NOERR) then
              print *, 'Error opening netcdf boundary file'
              call exit(0)
          endif
          rc = nf_inq_varid(ncid, 'time', varid)
          if (rc .ne. NF_NOERR) then
              print *, 'Error finding TIME variable in boundary file'
              call exit(0)
          endif

          if (hour .lt. 0) then
              nc_index = begin_j_date*24-24+1
              rc = nf_get_var1_double(ncid, varid, nc_index, dtime)
              v_f%time = NINT(dtime)
              if (rc .ne. NF_NOERR) then
                  print *, 'Error finding 1st time in boundary file'
                  call exit(0)
              endif
              v_f%time_offset = v_f%time
              rc = nf_inq_dimlen(ncid, varid, v_f%length);
              rc = nf_get_att_int(ncid, varid, 'begin_year',v_f%begin_year)
!             rc = nf_inq_varid(ncid, 'iceh', varid)
!             if (rc .eq. 0) then
!                 rc = nf_get_var1_double(ncid, varid, nc_index, v_f%iceh1)
!             else
!                 v_f%iceh1 = h_min
!             endif
          else
              nc_index = hour-v_f%time_offset + (year-v_f%begin_year)*8760  ! this will fuck up if begin_year is a leap year
              nc_index = hour
              rc = nf_get_var1_int(ncid, varid, nc_index, dtime)
              if (rc .ne. NF_NOERR) then
                  print *, 'Error finding time in boundary file'
                  call exit(0)
              endif
              v_f%time = NINT(dtime)
          endif

          rc = nf_inq_varid(ncid, 'at', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%at)
          else
              v_f%at = c9999;
          endif

          rc = nf_inq_varid(ncid, 'pressure', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%p)
          else
              v_f%p = c9999;
          endif

          rc = nf_inq_varid(ncid, 'ws', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%ws)
          else
              v_f%ws = c9999;
          endif

          rc = nf_inq_varid(ncid, 'cloud', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%fc)
          else
              v_f%fc = c9999;
          endif

          rc = nf_inq_varid(ncid, 'rh', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%rhum)
          else
              v_f%rhum = c9999;
          endif

          rc = nf_inq_varid(ncid, 'shortwave', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%swd)
          else
              v_f%swd = c9999;
          endif

          rc = nf_inq_varid(ncid, 'longwave_down', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%lwd)
          else
              v_f%lwd = c9999;
          endif

          rc = nf_inq_varid(ncid, 'ohf', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%fw)
          else
              v_f%fw = c9999;
          endif

          rc = nf_inq_varid(ncid, 'snowh', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%sh)
          else
              v_f%sh = c9999;
          endif

          rc = nf_inq_varid(ncid, 'snow', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%pr_snow) ! kg/m2/s
          else
              v_f%pr_snow = c9999;
          endif

          rc = nf_inq_varid(ncid, 'sh', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%shum)
          else
              v_f%shum = c9999;
          endif

          rc = nf_inq_varid(ncid, 'sensible', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%sens_heat_flux)
          else
              v_f%sens_heat_flux = c9999;
          endif

          rc = nf_inq_varid(ncid, 'latent', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%late_heat_flux)
          else
              v_f%late_heat_flux = c9999;
          endif

          rc = nf_inq_varid(ncid, 'lat', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%lat)
          else
              v_f%lat = c9999;
          endif

          rc = nf_inq_varid(ncid, 'lon', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%lon)
          else
              v_f%lon = c9999;
          endif

          rc = nf_inq_varid(ncid, 't', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%t)
          else
              v_f%t = c9999;
          endif

          rc = nf_inq_varid(ncid, 's', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%s)
          else
              v_f%s = c9999;
          endif

          rc = nf_inq_varid(ncid, 'no3', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%no3)
          else
              v_f%no3 = c9999;
          endif

          rc = nf_inq_varid(ncid, 'nh4', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%nh4)
          else
              v_f%nh4 = c9999;
          endif

          rc = nf_inq_varid(ncid, 'po4', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%po4)
          else
              v_f%po4 = c9999;
          endif

          rc = nf_inq_varid(ncid, 'sioh4', varid)
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%sioh4)
          else
              v_f%sioh4 = c9999;
          endif

          ! snow overrides pr
          rc = nf_inq_varid(ncid, 'snow', varid)
          v_f%pr_snow = c9999;
          if (rc .eq. 0) then
              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%pr)
          else
              v_f%pr = c9999;
          endif

          ! snow overrides pr
!          rc = nf_inq_varid(ncid, 'pr', varid)
!          if (rc .eq. 0) then
!              v_f%pr_snow = v_f%pr;
!              rc = nf_get_var1_double(ncid, varid, nc_index, v_f%pr)
!          else
!              v_f%pr = c9999;
!          endif

          rc = nf_close(ncid)

          if ((v_f%lon .ne. c9999) .and. (v_f%lat .ne. c9999)) then
              ! update NCEP indexes
              call gc_ll2ij(v_f%lon,v_f%lat,v_f%x_mp,v_f%y_mp)

              ! update ECMWF indexes
              call ec_ll2ij(v_f%lon,v_f%lat,v_f%x_ec,v_f%y_ec)

              ! update ECMWF indexes
              call ecint_ll2ij(v_f%lon,v_f%lat,v_f%x_eci,v_f%y_eci)

              ! update degree grid indexes
              call sia2_dgxy_from_ll(v_f%lat,v_f%lon,v_f%x_dg,v_f%y_dg)
          endif


      end SUBROUTINE sia2_forcing_load_v_f


      subroutine gc_read(year, hour, Ed, var)

! Reads Ed after reflection (Ed-0), Ed before reflection (Ed+0) and
! fraction of direct irradiance of total (other part is diffuse irrad)

! Input: year, hour (since 1-1-year 0:00), var = which variable to read
!     var=1: Edm0
!     var=2: Ed before reflection (Edp0)
!     var=3: Edir % of Total (before cloud)

! Output: global Ed* at 31 wavelength for that day and time

! Ed*0 = light at 31 wavelengths (400, 410, ... 700), units: uEin/m2/s/nm
! To convert to PAR: 5*wl(1), 5*wl(31), 10*wl(2:30) -> uEin/m2/s
! For Edir output is percentage
      use sia2_globals
      implicit none
#include "netcdf.inc"
#include "gc_cdf.h"

      integer year, hour, var, i,j,kk
      double precision :: Ed(mp_x,mp_y,wavl)
      character (len=4) :: varnm, cyear
      character (len=16) :: datadir
      character (len=28) :: cdffile


      datadir = 'data/gc_par_ice/'

      if (var .eq. 1) varnm='Edir'
      if (var .eq. 2) varnm='Edif'

      write(cyear,'(i4)')year

      cdffile = datadir//varnm//'_'//cyear//'.nc'

      rc = nf_open(cdffile,nf_nowrite,ncid)
      if (rc .ne. NF_NOERR) then
          print *, 'Error opening netcdf GC file'
          call exit(0)
      endif
      rc = nf_inq_varid(ncid, varnm, varid)
      rc = nf_get_att_real(ncid, varid, 'add_offset', offset)
      rc = nf_get_att_real(ncid, varid, 'scale_factor', sclf)

      corner4(4) = hour-gc_offset
      corner4(3) = 1
      corner4(2) = 1
      corner4(1) = 1
      edge4(4) = 1
      edge4(3) = 31
      edge4(2) = 94
      edge4(1) = 192

      rc = nf_get_vara_int2(ncid, varid, corner4, edge4, Ed_int2)  !Ed_int2 is global...
      if (rc .ne. NF_NOERR) then
          print *, 'Error ',rc,'opening reading GC PAR data from cdf file'
          call exit(0)
      endif



      rc = nf_close(ncid)

      do kk=1,wavl
         do j=1,mp_y
            do i=1,mp_x
               Ed(i,j,kk) = dfloat(Ed_int2(i,j,kk)) * sclf + offset
               if (Ed(i,j,kk) .lt. 0.0001) then
                    Ed(i,j,kk) = 0.
               endif
               !print *, Edi2(i,j,kk),'x', sclf,'+',offset,'=',Ed(i,j,kk)

               !Ed2(i,j) = Edi2(i,j,6)
            enddo
         enddo
      enddo


      !open(unit=31, file='gc_par_31.int2', form='UNFORMATTED', access='direct', recl=2*192*94)
      !write(31,rec=1) Ed2
      !close(31)
      !call exit(0)

      return
      end subroutine gc_read


! Routine above from Gert was copied below here (and modifed) to read in
! air temperature climatologies for caluclation of ice temp when using SSMI
!==========================================================================

      subroutine ncep_read(year, day, ncep_out, ncep_t)

      use sia2_globals
      implicit none
#include "netcdf.inc"
#include "gc_cdf.h"

      integer year, day, i, j, ncep_t
      !integer*2 air2(mp_x,mp_y)
      double precision, dimension(mp_x,mp_y) :: ncep_out
      character cyear*4
      character airdir*27, airvarnm*3, aircdffile*34
      character pdir*29, pvarnm*4, pcdffile*36
      character hdir*28, hvarnm*4, hcdffile*35
      character fcdir*30, fcvarnm*4, fccdffile*37
      character u10dir*29, u10varnm*4, u10cdffile*36
      character v10dir*27, v10varnm*4, v10cdffile*34
      character prdir*30, prvarnm*5, prcdffile*37

      write(cyear,'(i4)')year


      if (ncep_t .eq. ncep_at) then
          airdir = 'data/NCEP-DOE/air.2m.gauss.'
          airvarnm = 'air'
          aircdffile = airdir//cyear//'.nc'
          rc = nf_open(aircdffile,nf_nowrite,ncid)
          rc = nf_inq_varid(ncid, airvarnm, varid)
      elseif (ncep_t .eq. ncep_p) then
          pdir = 'data/NCEP-DOE/pres.sfc.gauss.'
          pvarnm = 'pres'
          pcdffile = pdir//cyear//'.nc'
          rc = nf_open(pcdffile,nf_nowrite,ncid)
          rc = nf_inq_varid(ncid, pvarnm, varid)
      elseif (ncep_t .eq. ncep_h) then
          hdir = 'data/NCEP-DOE/shum.2m.gauss.'
          hvarnm = 'shum'
          hcdffile = hdir//cyear//'.nc'
          rc = nf_open(hcdffile,nf_nowrite,ncid)
          rc = nf_inq_varid(ncid, hvarnm, varid)
      elseif (ncep_t .eq. ncep_fc) then
          fcdir ='data/NCEP-DOE/tcdc.eatm.gauss.'
          fcvarnm = 'tcdc'
          fccdffile = fcdir//cyear//'.nc'
          rc = nf_open(fccdffile,nf_nowrite,ncid)
          rc = nf_inq_varid(ncid, fcvarnm, varid)
      elseif (ncep_t .eq. ncep_u10) then
          u10dir ='data/NCEP-DOE/uwnd.10m.gauss.'
          u10varnm = 'uwnd'
          u10cdffile = u10dir//cyear//'.nc'
          rc = nf_open(u10cdffile,nf_nowrite,ncid)
          rc = nf_inq_varid(ncid, u10varnm, varid)
      elseif (ncep_t .eq. ncep_v10) then
          v10dir ='data/NCEP-DOE/vwnd.10m.gauss.'
          v10varnm = 'vwnd'
          v10cdffile = v10dir//cyear//'.nc'
          rc = nf_open(v10cdffile,nf_nowrite,ncid)
          rc = nf_inq_varid(ncid, v10varnm, varid)
      elseif (ncep_t .eq. ncep_pr) then
          prdir ='data/NCEP-DOE/prate.sfc.gauss.'
          prvarnm = 'prate'
          prcdffile = prdir//cyear//'.nc'
          rc = nf_open(prcdffile,nf_nowrite,ncid)
          rc = nf_inq_varid(ncid, prvarnm, varid)
      endif

      rc = nf_get_att_real(ncid, varid, 'add_offset', offset)
      rc = nf_get_att_real(ncid, varid, 'scale_factor', sclf)


      if ((ncep_t .eq. ncep_p) .or. (ncep_t .eq. ncep_fc) .or. (ncep_t .eq. ncep_pr)) then

          corner3(3) = day
          corner3(2) = 1
          corner3(1) = 1
          edge3(3) = 1
          edge3(2) = 94
          edge3(1) = 192

          rc = nf_get_vara_int2(ncid, varid, corner3, edge3, mp_grid_int2) ! mp_grid_int2 is a global...
      else

          corner4(4) = day
          corner4(3) = 1
          corner4(2) = 1
          corner4(1) = 1
          edge4(4) = 1
          edge4(3) = 1
          edge4(2) = 94
          edge4(1) = 192

          rc = nf_get_vara_int2(ncid, varid, corner4, edge4, mp_grid_int2) ! mp_grid_int2 is a global...
      endif

      rc = nf_close(ncid)


      do j=1,94
         do i=1,192
            ncep_out(i,j) = dfloat(mp_grid_int2(i,j))*sclf + offset
         enddo
      enddo

      return
      end subroutine ncep_read




! Routine above from Gert was copied below here (and modifed) to read in
! air temperature climatologies for caluclation of ice temp when using SSMI
!==========================================================================

      subroutine ncep6_read(year, day, ncep_out, ncep_t)

      use sia2_globals
      implicit none
#include "netcdf.inc"
#include "gc_cdf.h"

      integer year, day, i, j, ncep_t
      !integer*2 air2(mp_x,mp_y)
      double precision, dimension(mp_x,mp_y) :: ncep_out
       character cyear*4
      character airdir*35, airvarnm*3, aircdffile*42
      character pdir*37, pvarnm*4, pcdffile*44
      character hdir*36, hvarnm*4, hcdffile*43
      character fcdir*38, fcvarnm*4, fccdffile*45
      character u10dir*37, u10varnm*4, u10cdffile*44
      character v10dir*37, v10varnm*4, v10cdffile*44
      character prdir*38, prvarnm*5, prcdffile*45

      write(cyear,'(i4)')year


      if (ncep_t .eq. ncep_at) then
          airdir = 'data/NCEP-DOE/air.2m.gauss.4xdaily.'
          airvarnm = 'air'
          aircdffile = airdir//cyear//'.nc'
          rc = nf_open(aircdffile,nf_nowrite,ncid)
          rc = nf_inq_varid(ncid, airvarnm, varid)
      elseif (ncep_t .eq. ncep_p) then
          pdir = 'data/NCEP-DOE/pres.sfc.gauss.4xdaily.'
          pvarnm = 'pres'
          pcdffile = pdir//cyear//'.nc'
          rc = nf_open(pcdffile,nf_nowrite,ncid)
          rc = nf_inq_varid(ncid, pvarnm, varid)
      elseif (ncep_t .eq. ncep_h) then
          hdir = 'data/NCEP-DOE/shum.2m.gauss.4xdaily.'
          hvarnm = 'shum'
          hcdffile = hdir//cyear//'.nc'
          rc = nf_open(hcdffile,nf_nowrite,ncid)
          rc = nf_inq_varid(ncid, hvarnm, varid)
      elseif (ncep_t .eq. ncep_fc) then
          fcdir ='data/NCEP-DOE/tcdc.eatm.gauss.4xdaily.'
          fcvarnm = 'tcdc'
          fccdffile = fcdir//cyear//'.nc'
          rc = nf_open(fccdffile,nf_nowrite,ncid)
          rc = nf_inq_varid(ncid, fcvarnm, varid)
      elseif (ncep_t .eq. ncep_u10) then
          u10dir ='data/NCEP-DOE/uwnd.10m.gauss.4xdaily.'
          u10varnm = 'uwnd'
          u10cdffile = u10dir//cyear//'.nc'
          rc = nf_open(u10cdffile,nf_nowrite,ncid)
          rc = nf_inq_varid(ncid, u10varnm, varid)
      elseif (ncep_t .eq. ncep_v10) then
          v10dir ='data/NCEP-DOE/vwnd.10m.gauss.4xdaily.'
          v10varnm = 'vwnd'
          v10cdffile = v10dir//cyear//'.nc'
          rc = nf_open(v10cdffile,nf_nowrite,ncid)
          rc = nf_inq_varid(ncid, v10varnm, varid)
      elseif (ncep_t .eq. ncep_pr) then
          prdir ='data/NCEP-DOE/prate.sfc.gauss.4xdaily.'
          prvarnm = 'prate'
          prcdffile = prdir//cyear//'.nc'
          rc = nf_open(prcdffile,nf_nowrite,ncid)
          rc = nf_inq_varid(ncid, prvarnm, varid)
      endif

      rc = nf_get_att_real(ncid, varid, 'add_offset', offset)
      rc = nf_get_att_real(ncid, varid, 'scale_factor', sclf)

      if ((ncep_t .eq. ncep_p) .or. (ncep_t .eq. ncep_fc) .or. (ncep_t .eq. ncep_pr)) then

          corner3(3) = day
          corner3(2) = 1
          corner3(1) = 1
          edge3(3) = 1
          edge3(2) = 94
          edge3(1) = 192

          rc = nf_get_vara_int2(ncid, varid, corner3, edge3, mp_grid_int2) ! mp_grid_int2 is a global...
      else

          corner4(4) = day
          corner4(3) = 1
          corner4(2) = 1
          corner4(1) = 1
          edge4(4) = 1
          edge4(3) = 1
          edge4(2) = 94
          edge4(1) = 192

          rc = nf_get_vara_int2(ncid, varid, corner4, edge4, mp_grid_int2) ! mp_grid_int2 is a global...
      endif

      rc = nf_close(ncid)


      do j=1,94
         do i=1,192
            ncep_out(i,j) = dfloat(mp_grid_int2(i,j))*sclf + offset
         enddo
      enddo

      return
      end subroutine ncep6_read

! Routine above from Gert was copied below here (and modifed) to read in
! air temperature climatologies for caluclation of ice temp when using SSMI
!==========================================================================

      subroutine ecmwf6_int_read(year, time_slice, ecmwf_out)

      use sia2_globals
      implicit none
#include "netcdf.inc"
#include "gc_cdf.h"

      integer :: year, time_slice, i, j, ii
      double precision :: sclf_d,offset_d
      type (ecmwf_int_type) :: ecmwf_out
      character (len=4) :: cyear
      character (len=21) :: ecmwf_dir
      character (len=35) :: ecmwf_file
      character (len=3) :: avarnm='t2m'  ! 2m air temperature
      character (len=3) :: pvarnm='msl'  ! 2m air pressure
      character (len=3) :: hvarnm='d2m'  ! 2m dewpoint
      character (len=3) :: fcvarnm='tcc'  ! total cloud cover
      character (len=3) :: u10varnm='u10' ! 10m u windspeed
      character (len=3) :: v10varnm='v10' ! 10m v windspeed

      write(cyear,'(i4)')year

      ecmwf_dir = 'data/ECMWF_int/ECMWF_'
      ecmwf_file = ecmwf_dir//cyear//'_int_4x.nc'
      rc = nf_open(ecmwf_file,nf_nowrite,ncid)
      if (rc .ne. nf_noerr) then
          print *,'Error opening netcdf file: ',ecmwf_file
          print *,'NF_ERROR ',rc
      endif

      corner3(3) = time_slice
      corner3(2) = 1
      corner3(1) = 1
      edge3(3) = 1
      edge3(2) = eci_y
      edge3(1) = eci_x

      do ii=1,6

          select case (ii)
              case (1)
                  rc = nf_inq_varid(ncid, avarnm, varid)
              case (2)
                  rc = nf_inq_varid(ncid, pvarnm, varid)
              case (3)
                  rc = nf_inq_varid(ncid, hvarnm, varid)
              case (4)
                  rc = nf_inq_varid(ncid, fcvarnm, varid)
              case (5)
                  rc = nf_inq_varid(ncid, u10varnm, varid)
              case (6)
                  rc = nf_inq_varid(ncid, v10varnm, varid)
          end select

          if (rc .ne. nf_noerr) then
              print *,'Error finding netcdf var: ',ecmwf_file,ii
              print *,'NF_ERROR ',rc
          endif
          rc = nf_get_att_double(ncid, varid, 'add_offset', offset_d)
          rc = nf_get_att_double(ncid, varid, 'scale_factor', sclf_d)
          rc = nf_get_vara_int2(ncid, varid, corner3, edge3, eci_grid_int2) ! ec_grid_int2 is a global...
          if (rc .ne. nf_noerr) then
              print *,'Error finding netcdf var: ',ecmwf_file,ii
              print *,'NF_ERROR ',rc
          endif

          select case (ii)
              case (1)
                  do j=1,eci_y
                     do i=1,eci_x
                        ecmwf_out%at(i,j) = dble(eci_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
              case (2)
                  do j=1,eci_y
                     do i=1,eci_x
                        ecmwf_out%p(i,j) = dble(eci_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
              case (3)
                  do j=1,eci_y
                     do i=1,eci_x
                        ecmwf_out%dpt(i,j) = dble(eci_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
              case (4)
                  do j=1,eci_y
                     do i=1,eci_x
                        ecmwf_out%fc(i,j) = dble(eci_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
              case (5)
                  do j=1,eci_y
                     do i=1,eci_x
                        ecmwf_out%u10(i,j) = dble(eci_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
              case (6)
                  do j=1,eci_y
                     do i=1,eci_x
                        ecmwf_out%v10(i,j) = dble(eci_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
          end select

      enddo

      rc = nf_close(ncid)

      return
      end subroutine ecmwf6_int_read


      subroutine ecmwf12_int_read(year, time_slice, ecmwf_out)

      use sia2_globals
      implicit none
#include "netcdf.inc"
#include "gc_cdf.h"

      integer :: year, time_slice, i, j, ii
      double precision :: sclf_d,offset_d
      type (ecmwf_int_type) :: ecmwf_out
      character (len=4) :: cyear
      character (len=21) :: ecmwf_dir
      character (len=35) :: ecmwf_file
      character (len=2) :: prvarnm='tp'  ! total precipitation
      character (len=3) :: rvarnm='ssr'  ! surface solar radiation

      write(cyear,'(i4)')year

      ecmwf_dir = 'data/ECMWF_int/ECMWF_'
      ecmwf_file = ecmwf_dir//cyear//'_int_2x.nc'
      rc = nf_open(ecmwf_file,nf_nowrite,ncid)
      if (rc .ne. nf_noerr) then
          print *,'Error opening netcdf file: ',ecmwf_file
          print *,'NF_ERROR ',rc
      endif

      corner3(3) = time_slice
      corner3(2) = 1
      corner3(1) = 1
      edge3(3) = 1
      edge3(2) = eci_y
      edge3(1) = eci_x

      do ii=1,2

          select case (ii)
              case (1)
                  rc = nf_inq_varid(ncid, prvarnm, varid)
              case (2)
                  rc = nf_inq_varid(ncid, rvarnm, varid)
          end select

          if (rc .ne. nf_noerr) then
              print *,'Error finding netcdf var: ',ecmwf_file,ii
              print *,'NF_ERROR ',rc
          endif
          rc = nf_get_att_double(ncid, varid, 'add_offset', offset_d)
          rc = nf_get_att_double(ncid, varid, 'scale_factor', sclf_d)
          rc = nf_get_vara_int2(ncid, varid, corner3, edge3, eci_grid_int2) ! ec_grid_int2 is a global...
          if (rc .ne. nf_noerr) then
              print *,'Error finding netcdf var: ',ecmwf_file,ii
              print *,'NF_ERROR ',rc
          endif

          select case (ii)
              case (1)
                  do j=1,eci_y
                     do i=1,eci_x
                        ecmwf_out%pr(i,j) = dble(eci_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
              case (2)
                  do j=1,eci_y
                     do i=1,eci_x
                        ecmwf_out%ssr(i,j) = dble(eci_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
          end select

      enddo

      rc = nf_close(ncid)

      return
      end subroutine ecmwf12_int_read


      subroutine ecmwf3_int_read(year, time_slice, ecmwf_out)

      use sia2_globals
      implicit none
#include "netcdf.inc"
#include "gc_cdf.h"

      integer :: year, time_slice, i, j, ii
      double precision :: sclf_d,offset_d
      type (ecmwf_int_type) :: ecmwf_out
      character (len=4) :: cyear
      character (len=21) :: ecmwf_dir
      character (len=38) :: ecmwf_file
      character (len=2) :: prvarnm='tp'  ! total precipitation
!      character (len=3) :: rvarnm='ssr'  ! surface solar radiation

      write(cyear,'(i4)')year

      ecmwf_dir = 'data/ECMWF_int/ECMWF_'
      ecmwf_file = ecmwf_dir//cyear//'_int_3h_pr.nc'
      rc = nf_open(ecmwf_file,nf_nowrite,ncid)
      if (rc .ne. nf_noerr) then
          print *,'Error opening netcdf file: ',ecmwf_file
          print *,'NF_ERROR ',rc
      endif

      corner3(3) = time_slice
      corner3(2) = 1
      corner3(1) = 1
      edge3(3) = 1
      edge3(2) = eci_y
      edge3(1) = eci_x

!      do ii=1,2

!          select case (ii)
!              case (1)
                  rc = nf_inq_varid(ncid, prvarnm, varid)
!              case (2)
!                  rc = nf_inq_varid(ncid, rvarnm, varid)
!         end select

          if (rc .ne. nf_noerr) then
              print *,'Error finding netcdf var: ',ecmwf_file,ii
              print *,'NF_ERROR ',rc
          endif
          rc = nf_get_att_double(ncid, varid, 'add_offset', offset_d)
          rc = nf_get_att_double(ncid, varid, 'scale_factor', sclf_d)
          rc = nf_get_vara_int2(ncid, varid, corner3, edge3, eci_grid_int2) ! ec_grid_int2 is a global...
          if (rc .ne. nf_noerr) then
              print *,'Error finding netcdf var: ',ecmwf_file,ii
              print *,'NF_ERROR ',rc
          endif

!          select case (ii)
!              case (1)
                  do j=1,eci_y
                     do i=1,eci_x
                        ecmwf_out%pr(i,j) = dble(eci_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
!              case (2)
!                  do j=1,eci_y
!                     do i=1,eci_x
!                        ecmwf_out%ssr(i,j) = dble(eci_grid_int2(i,j))*sclf_d + offset_d
!                     enddo
!                  enddo
!         end select

!      enddo

      rc = nf_close(ncid)

      return
      end subroutine ecmwf3_int_read




      subroutine ecmwf6_read(year, time_slice, ecmwf_out)

      use sia2_globals
      implicit none
#include "netcdf.inc"
#include "gc_cdf.h"

      integer :: year, time_slice, i, j, ii
      double precision :: sclf_d,offset_d
      type (ecmwf_type) :: ecmwf_out
      character (len=4) :: cyear
      character (len=17) :: ecmwf_dir
      character (len=24) :: ecmwf_file
      character (len=3) :: avarnm='p2t'  ! 2m air temperature
      character (len=3) :: pvarnm='msl'  ! 2m air pressure
      character (len=3) :: hvarnm='p2d'  ! 2m dewpoint
      character (len=3) :: fcvarnm='tcc'  ! total cloud cover
      character (len=4) :: u10varnm='p10u' ! 10m u windspeed
      character (len=4) :: v10varnm='p10v' ! 10m v windspeed
      character (len=2) :: prvarnm='tp'  ! total precipitation
      character (len=3) :: rvarnm='ssr'  ! surface solar radiation

      write(cyear,'(i4)')year

      ecmwf_dir = 'data/ECMWF/ECMWF_'
      ecmwf_file = ecmwf_dir//cyear//'.nc'
      rc = nf_open(ecmwf_file,nf_nowrite,ncid)
      if (rc .ne. nf_noerr) then
          print *,'Error opening netcdf file: ',ecmwf_file
          print *,'NF_ERROR ',rc
      endif

      corner3(3) = time_slice
      corner3(2) = 1
      corner3(1) = 1
      edge3(3) = 1
      edge3(2) = ec_y
      edge3(1) = ec_x

      do ii=1,8

          select case (ii)
              case (1)
                  rc = nf_inq_varid(ncid, avarnm, varid)
                  if (rc .ne. NF_NOERR) then
                      rc = nf_inq_varid(ncid, '2t', varid)
                  endif
              case (2)
                  rc = nf_inq_varid(ncid, pvarnm, varid)
              case (3)
                  rc = nf_inq_varid(ncid, hvarnm, varid)
                  if (rc .ne. NF_NOERR) then
                      rc = nf_inq_varid(ncid, '2d', varid)
                  endif
              case (4)
                  rc = nf_inq_varid(ncid, fcvarnm, varid)
              case (5)
                  rc = nf_inq_varid(ncid, u10varnm, varid)
                  if (rc .ne. NF_NOERR) then
                      rc = nf_inq_varid(ncid, '10u', varid)
                  endif
              case (6)
                  rc = nf_inq_varid(ncid, v10varnm, varid)
                  if (rc .ne. NF_NOERR) then
                      rc = nf_inq_varid(ncid, '10v', varid)
                  endif
              case (7)
                  rc = nf_inq_varid(ncid, prvarnm, varid)
              case (8)
                  rc = nf_inq_varid(ncid, rvarnm, varid)
          end select

          if (rc .ne. nf_noerr) then
              print *,'Error finding netcdf var: ',ecmwf_file,ii
              print *,'NF_ERROR ',rc
          endif
          rc = nf_get_att_double(ncid, varid, 'add_offset', offset_d)
          rc = nf_get_att_double(ncid, varid, 'scale_factor', sclf_d)
          rc = nf_get_vara_int2(ncid, varid, corner3, edge3, ec_grid_int2) ! ec_grid_int2 is a global...
          if (rc .ne. nf_noerr) then
              print *,'Error finding netcdf var: ',ecmwf_file,ii
              print *,'NF_ERROR ',rc
          endif

          select case (ii)
              case (1)
                  do j=1,ec_y
                     do i=1,ec_x
                        ecmwf_out%at(i,j) = dble(ec_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
              case (2)
                  do j=1,ec_y
                     do i=1,ec_x
                        ecmwf_out%p(i,j) = dble(ec_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
              case (3)
                  do j=1,ec_y
                     do i=1,ec_x
                        ecmwf_out%dpt(i,j) = dble(ec_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
              case (4)
                  do j=1,ec_y
                     do i=1,ec_x
                        ecmwf_out%fc(i,j) = dble(ec_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
              case (5)
                  do j=1,ec_y
                     do i=1,ec_x
                        ecmwf_out%u10(i,j) = dble(ec_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
              case (6)
                  do j=1,ec_y
                     do i=1,ec_x
                        ecmwf_out%v10(i,j) = dble(ec_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
              case (7)
                  do j=1,ec_y
                     do i=1,ec_x
                        ecmwf_out%pr(i,j) = dble(ec_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
              case (8)
                  do j=1,ec_y
                     do i=1,ec_x
                        ecmwf_out%ssr(i,j) = dble(ec_grid_int2(i,j))*sclf_d + offset_d
                     enddo
                  enddo
          end select

      enddo

      rc = nf_close(ncid)

      return
      end subroutine ecmwf6_read



! Routine above from Gert was copied below here (and modifed) to read in
! SSMI icecon values stored in netcdf format
!==========================================================================

      subroutine SSMI_icecon_read(year, day)

      use sia2_globals
      implicit none
#include "netcdf.inc"
#include "gc_cdf.h"

      integer year, day, i, j
      character (LEN=37) :: datadir
      character (LEN=6) :: varnm
      character (LEN=44) :: cdffile
      character (LEN=4) :: cyear


      if (hemi .eq. 0) then
          datadir = 'data/SSMI/netcdf/EASE_SH_25km_icecon_'
      else
          datadir = 'data/SSMI/netcdf/EASE_NH_25km_icecon_'
      endif

      varnm = 'icecon'
      write(cyear,'(i4)')year

      cdffile = datadir//cyear//'.nc'

      rc = nf_open(cdffile,nf_nowrite,ncid)
      rc = nf_inq_varid(ncid, varnm, varid)

      corner3(3) = day
      corner3(2) = 1
      corner3(1) = 1
      edge3(3) = 1
      edge3(2) = 721
      edge3(1) = 721


      rc = nf_get_vara_int2(ncid, varid, corner3, edge3, SSMI_grid_int2)
      rc = nf_close(ncid)

      ! ----- taken out, b/c we need to account for biomass lost outside
      !do j=1,721
      !   do i=1,721
      !      ease_f%icecon(i,j) = dfloat(SSMI_grid_int2(i,j))/100.
      !   enddo
      !enddo

      return
      end


! Routine above from Gert was copied below here (and modifed) to read in
! SSMI snow depth values stored in netcdf format
!==========================================================================

      subroutine SSMI_snow_read(year, day)

      use sia2_globals
      implicit none
#include "netcdf.inc"
#include "gc_cdf.h"

      integer year, day, i, j
      character (LEN=40) :: datadir
      character (LEN=9) :: varnm
      character (LEN=47) :: cdffile
      character (LEN=4) :: cyear

      datadir = 'data/SSMI/netcdf/EASE_SH_25km_snowdepth_'
      varnm = 'snowdepth'
      write(cyear,'(i4)')year

      cdffile = datadir//cyear//'.nc'

      rc = nf_open(cdffile,nf_nowrite,ncid)
      rc = nf_inq_varid(ncid, varnm, varid)

      corner3(3) = day
      corner3(2) = 1
      corner3(1) = 1
      edge3(3) = 1
      edge3(2) = 721
      edge3(1) = 721


      rc = nf_get_vara_int2(ncid, varid, corner3, edge3, SSMI_grid_int2)
      rc = nf_close(ncid)

      return
      end



! Routine above from Gert was copied below here (and modifed) to read in
! SSMI ice vector values stored in netcdf format
!==========================================================================

      subroutine SSMI_icevec_read(year, day)

      use sia2_globals
      implicit none
#include "netcdf.inc"
#include "gc_cdf.h"

      integer year, day, i, j
      character (LEN=49) :: datadir
      character (LEN=1) :: varnm_u
      character (LEN=1) :: varnm_v
      character (LEN=4) :: varnm_flag
      character (LEN=56) :: cdffile
      character (LEN=4) :: cyear

      corner3(3) = day
      corner3(2) = 1
      corner3(1) = 1
      edge3(3) = 1
      if (hemi .eq. 0) then
          edge3(2) = 321
          edge3(1) = 321
          datadir = 'data/SSMI/netcdf/EASE_SH_25km_321x321_icevectors_'
      else
          edge3(2) = 361
          edge3(1) = 361
          datadir = 'data/SSMI/netcdf/EASE_NH_25km_361x361_icevectors_'
      endif

      varnm_u = 'u'
      varnm_v = 'v'
      varnm_flag = 'flag'
      write(cyear,'(i4)')year

      cdffile = datadir//cyear//'.nc'

      rc = nf_open(cdffile,nf_nowrite,ncid)
      if (rc .ne. 0) then
          print *,'icevec open failed',rc,cdffile
          call exit(0)
      endif

      rc = nf_inq_varid(ncid, varnm_u, varid)
      rc = nf_get_vara_int2(ncid, varid, corner3, edge3, icevec_grid_int2(:,:,1))

      rc = nf_inq_varid(ncid, varnm_v, varid)
      rc = nf_get_vara_int2(ncid, varid, corner3, edge3, icevec_grid_int2(:,:,2))

      rc = nf_inq_varid(ncid, varnm_flag, varid)
      rc = nf_get_vara_int2(ncid, varid, corner3, edge3, icevec_grid_int2(:,:,3))

      rc = nf_close(ncid)


      return
      end


! Routine above from Gert was copied below here (and modifed) to read in
! air temperature climatologies for caluclation of ice temp when using SSMI
!==========================================================================

      subroutine sia2_woa_read(mo,depth,woa_out,woa_var)

          use sia2_globals
          implicit none
#include "netcdf.inc"
#include "gc_cdf.h"

          integer :: mo, depth, woa_var
          real, dimension(dg_x,dg_y) :: woa_out
          character (LEN=4) :: cyear
          character (LEN=8) :: varnm
          character (LEN=25) :: cdffile

          if (woa_var .eq. woa_t) then
              varnm = 't0112an1'
              cdffile = 'data/WOA_2005/t0112an1.nc'
              rc = nf_open(cdffile,nf_nowrite,ncid)
              rc = nf_inq_varid(ncid, varnm, varid)
          elseif (woa_var .eq. woa_s) then
              varnm = 's0112an1'
              cdffile = 'data/WOA_2005/s0112an1.nc'
              rc = nf_open(cdffile,nf_nowrite,ncid)
              rc = nf_inq_varid(ncid, varnm, varid)
          elseif (woa_var .eq. woa_n) then
              varnm = 'n0112an1'
              cdffile = 'data/WOA_2005/n0112an1.nc'
              rc = nf_open(cdffile,nf_nowrite,ncid)
              rc = nf_inq_varid(ncid, varnm, varid)
          elseif (woa_var .eq. woa_p) then
              varnm = 'p0112an1'
              cdffile = 'data/WOA_2005/p0112an1.nc'
              rc = nf_open(cdffile,nf_nowrite,ncid)
              rc = nf_inq_varid(ncid, varnm, varid)
          elseif (woa_var .eq. woa_si) then
              varnm = 'i0112an1'
              cdffile = 'data/WOA_2005/i0112an1.nc'
              rc = nf_open(cdffile,nf_nowrite,ncid)
              rc = nf_inq_varid(ncid, varnm, varid)
          endif

          corner4(4) = mo
          corner4(3) = depth
          corner4(2) = 1
          corner4(1) = 1
          edge4(4) = 1
          edge4(3) = 1
          edge4(2) = dg_y
          edge4(1) = dg_x

          rc = nf_get_vara_real(ncid, varid, corner4, edge4, woa_out)

          rc = nf_close(ncid)

          return

      end subroutine sia2_woa_read




          subroutine sia2_sose_read(day5,yr,level,sose_out)

          use sia2_globals
          implicit none
#include "netcdf.inc"
#include "gc_cdf.h"

          integer :: day5,yr,level
          type (sose_type) :: sose_out
          character (LEN=20) :: datadir
          character (LEN=4) :: cyear
          character (LEN=1) :: varnm
          character (LEN=27) :: cdffile

      write(cyear,'(i4)')yr
      datadir = 'data/SOSE/EASE_SOSE_'
            cdffile = datadir//cyear//'.nc'

          corner4(4) = day5
          corner4(3) = level
          corner4(2) = 1
          corner4(1) = 1
          edge4(4) = 1
          edge4(3) = 1
          edge4(2) = grids_v
          edge4(1) = grids_h

          corner3(3) = day5
          corner3(2) = 1
          corner3(1) = 1
          edge3(3) = 1
          edge3(2) = grids_v
          edge3(1) = grids_h

            rc = nf_open(cdffile,nf_nowrite,ncid)

            varnm = 't'
            rc = nf_inq_varid(ncid, varnm, varid)
          rc = nf_get_vara_double(ncid, varid, corner4, edge4, sose_out%t)

            varnm = 's'
            rc = nf_inq_varid(ncid, varnm, varid)
          rc = nf_get_vara_double(ncid, varid, corner4, edge4, sose_out%s)

            varnm = 'u'
            rc = nf_inq_varid(ncid, varnm, varid)
          rc = nf_get_vara_double(ncid, varid, corner3, edge3, sose_out%u)

            varnm = 'v'
            rc = nf_inq_varid(ncid, varnm, varid)
          rc = nf_get_vara_double(ncid, varid, corner3, edge3, sose_out%v)

          rc = nf_close(ncid)

          return

      end subroutine sia2_sose_read



! Routine copied from Gert
!==========================================================================

      subroutine gc_ll2ij(lon, lat, il, jl)

! Converts input latitude/longitude to corresponding i,j value in Ed array
! latitude input range: -90-90
! longitude input range: 0-360

      implicit none
      real*8 lon, lat, latt, lont, p5lonr, lonr, latb(47), pm
      integer il, jl, first
      data latb /87.5975, 85.7031, 83.8020, 81.8991, 79.9954, 78.0914,  &
          76.1872, 74.2829, 72.3785, 70.4739, 68.5694, 66.6648,         &
          64.7602, 62.8556, 60.9510, 59.0462, 57.1416, 55.2369,         &
          53.3323, 51.4276, 49.5229, 47.6181, 45.7134, 43.8088,         &
          41.9040, 39.9993, 38.0946, 36.1898, 34.2852, 32.3805,         &
          30.4757, 28.5710, 26.6663, 24.7615, 22.8568, 20.9521,         &
          19.0473, 17.1426, 15.2379, 13.3331, 11.4284,  9.5237,         &
           7.6189,  5.7142,  3.8095,  1.9047,  0.0000/

! find latitude index from latb (lat bounds). (Gaussian grid -> cannot
! use center lat)
      latt = lat
      jl = 1

      if (latt .gt. 0) then
         do while (latt .lt. latb(jl))
            jl = jl+1
         enddo
      elseif (latt .lt. 0) then
         latt = -1.*latt
         do while (latt .le. latb(jl))  ! use le, bounds corners switched
            jl = jl+1
         enddo
         jl = 95-jl
      elseif (latt .eq. 0) then         ! use 0-1.9047S pixel
         jl = 48
      endif


! find longitude index, pixel size 1.875
      if(lon .lt. 0) then
          lon = lon+360
      endif

      lonr = 1.875
      p5lonr = 0.9375  ! half pixel size
      lont = lon + p5lonr
      il = int(lont/lonr) + 1
      if (il .eq. 193) il = 1  ! case of longitude gt 360-p5lonr


      return
      end



! Converts input latitude/longitude to corresponding i,j value in ECMWF ERA 40 grid
! this is a 2.5 degree grid
! latitude input range: -90-90
! longitude input range: 0-360

      subroutine ec_ll2ij(lon, lat, il, jl)

      implicit none

      integer :: il,jl
      double precision :: lon, lat, testvar


      if (lon .lt. 0.) then
          lon = 360 + lon ! convert negative longitudes to 0-360
      endif

      il = nint(lon/2.5)
      if (il .eq. 0) il = 144

      jl = 37 - nint(lat/2.5)

      if (jl .gt. 73 .or. il .gt. 144 .or. il .lt. 1 .or. jl .lt. 1) then
          testvar = -1
          testvar = sqrt(testvar)
      endif

      return
      end

! Converts input latitude/longitude to corresponding i,j value in ECMWF ERA INTERIM grid
! this is a 1.5 degree grid
! latitude input range: -90-90
! longitude input range: 0-360

      subroutine ecint_ll2ij(lon, lat, il, jl)

      implicit none

      integer :: il,jl
      double precision :: lon, lat, testvar


      if (lon .lt. 0.) then
          lon = 360 + lon ! convert negative longitudes to 0-360
      endif

      il = nint(lon/1.5)
      if (il .eq. 0) il = 240

      jl = 61 - nint(lat/1.5)

      if (jl .gt. 121 .or. il .gt. 240 .or. il .lt. 1 .or. jl .lt. 1) then
          testvar = -1
          testvar = sqrt(testvar)
      endif

      return
      end

