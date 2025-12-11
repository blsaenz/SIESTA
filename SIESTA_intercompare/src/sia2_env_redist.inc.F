
! ======================================================================
! START: Redistibute ice to different categories as it melts/grows out of the range
! of the current one
! ======================================================================

       if (dt_step .ge. 133 .and. mi .eq. 24928) then
           testvar = -1
       endif


       af_tot = 0.d0

	  ! convert tracers to bulk volume concentrations,zero other stuff
      do ic=1,ida_n
      do sc=1,sda_n

          if ((ice(sc,ic,mi)%af .gt. 0.d0) .and. (ice(sc,ic,mi)%z .gt. 0)) then
	          do ii=1,ice(sc,ic,mi)%z
                  tmp1 = ice(sc,ic,mi)%bv(ii)*c_001
                  ice(sc,ic,mi)%no3(ii) = ice(sc,ic,mi)%no3(ii)*tmp1
                  ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii)*tmp1
                  ice(sc,ic,mi)%po4(ii) = ice(sc,ic,mi)%po4(ii)*tmp1
                  ice(sc,ic,mi)%sioh4(ii) = ice(sc,ic,mi)%sioh4(ii)*tmp1
                  ice(sc,ic,mi)%smalg(ii) = ice(sc,ic,mi)%smalg(ii)*tmp1
                  ice(sc,ic,mi)%poc(ii) = ice(sc,ic,mi)%poc(ii)*tmp1
              enddo
              af_tot = af_tot + ice(sc,ic,mi)%af
          endif
	  enddo
	  enddo
      snowh1m = 0.d0
      melt_drain = 0.d0

! Add new ice to smallest category
! ----------------------------------------------------------------------
#include "sia2_env_redist_new_ice.inc.f90"


      do ic_out=1,ida_n

          ! merge categories based on depth
          ! ------------------------------------------------------------
          if ((ice(1,ic_out,mi)%af .gt. 0.d0) .and. (ice(1,ic_out,mi)%z .gt. 0)) then

              ic = 0
              id_mean = 0.d0
              af_total = 0.d0
              do sc=1,sda_n
                  if (ice(sc,ic_out,mi)%af .gt. 0.d0) then
                      id_mean = id_mean + &
                        ice(sc,ic_out,mi)%id(ice(sc,ic_out,mi)%z-z_sk)* &
                        ice(sc,ic_out,mi)%af
                      af_total = af_total + ice(sc,ic_out,mi)%af
                  endif
              enddo
              id_mean = id_mean/af_total

              if (id_mean .gt. ic_h_max(ic_out) .and. ic_out .lt. ida_n) then
                  !  Push current category into the next larger one
                  ic = ic_out + 1
              elseif (id_mean .lt. ic_h_min(ic_out) .and. ic_out .gt. 1) then
                  ! push current category into next smaller one, if available
                  ic = ic_out - 1
              endif

              if (ic .ne. 0) then

                  ! merge solid ice category subroutine
                  snow_dist_new = (ice(1,ic,mi)%snow_dist*ice(1,ic,mi)%af + &
                    ice(1,ic_out,mi)%snow_dist*ice(1,ic_out,mi)%af) / &
                    (ice(1,ic,mi)%af + ice(1,ic_out,mi)%af)
                  snow_rd_new = (ice(1,ic,mi)%snow_rd*ice(1,ic,mi)%af + &
                    ice(1,ic_out,mi)%snow_rd*ice(1,ic_out,mi)%af) / &
                    (ice(1,ic,mi)%af + ice(1,ic_out,mi)%af)
                  call sia2_env_merge_ice(1,ic,1,ic_out,mi,ice,pur,f,pl)
                  ice(1,ic,mi)%snow_dist = snow_dist_new
                  ice(1,ic,mi)%snow_rd = snow_rd_new
                  ! test to see if slush, to save cpu cycles
                  if (sda_n .gt. 1) then
                  if (ice(2,ic_out,mi)%af .gt. 0.d0) then
                      call sia2_env_merge_ice(2,ic,2,ic_out,mi,ice,pur,f,pl)
                  endif
                  endif

              endif

          endif

          ! merge flooded catagies back into standard categories
          ! if they aren't really flooded anymore
          ! ------------------------------------------------------------
          ! here we ignore the snow_dist var, and let it be reset as needed
          ! by wind re-distribution
          ! ------------------------------------------------------------
          if (sda_n .gt. 1) then
              if (ice(2,ic_out,mi)%af .gt. 0.) then

                  ! only look at top half of ice for slush desalination
                  ! jj = (ice(2,ic_out,mi)%z-z_sk)/2

                  ! look at all but last 5 layers
                  jj = ice(2,ic_out,mi)%z-z_sk-5

                  ic = 0
                  do ii=1,jj
                      if(ice(2,ic_out,mi)%bv(ii) .gt. bv_conv) then
                          ic = 1
                      endif
                  enddo

                  if (ic .eq. 0) then
                      ! merge no-longer-slush into solid ice category
                      call sia2_env_merge_ice(1,ic_out,2,ic_out,mi,ice,pur,f,pl)
                  endif
              endif
          endif
      enddo

	  ! convert tracers to brine concentrations
      do ic=1,ida_n
      do sc=1,sda_n
          if (ice(sc,ic,mi)%af .gt. 0.) then
	          do ii=1,ice(sc,ic,mi)%z
                  ! return to brine conc.
                  bv_mean = ice(sc,ic,mi)%bv(ii)*c_001
                  ice(sc,ic,mi)%no3(ii) = ice(sc,ic,mi)%no3(ii)/bv_mean
                  ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii)/bv_mean
                  ice(sc,ic,mi)%po4(ii) = ice(sc,ic,mi)%po4(ii)/bv_mean
                  ice(sc,ic,mi)%sioh4(ii) = ice(sc,ic,mi)%sioh4(ii)/bv_mean
                  ice(sc,ic,mi)%smalg(ii) = ice(sc,ic,mi)%smalg(ii)/bv_mean
                  ice(sc,ic,mi)%poc(ii) = ice(sc,ic,mi)%poc(ii)/bv_mean
              enddo
          endif
	  enddo
	  enddo


! ======================================================================
! END: Redistibute ice categories
! ======================================================================





