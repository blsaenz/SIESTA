     ! dz_sk is  the currency here - if this is set, actual growth/melt is happening
     ! otherwise, just advection/merge rectification
					debug_z = ice(sc,ic,mi)%no3(1:z_max)

					tmp1 = 1.
					tmp2 = ice(sc,ic,mi)%th(sk_z)
					if (ice(sc,ic,mi)%th(sk_z) .le. sk_th_min) then
					! skeletal layer is already below minimum (from advection)
						tmp1 = ice(sc,ic,mi)%th(sk_z)/sk_th_min	 ! dilute (can't create algae out of nothing)
						tmp2 = sk_th_min
						dz_sk = max(0.,dz_sk)
            print *,'sk_gl - diluting',sk_gl,dz_sk,ice(sc,ic,mi)%th(sk_z),dz_sk
					elseif ((ice(sc,ic,mi)%th(sk_z)+dz_sk) .lt. sk_th_min) then
					! dz_sk brings total under minimum - cut d_sk
						dz_sk = sk_th_min - ice(sc,ic,mi)%th(sk_z) !dz_sk takes it to sk_th_min
					endif
					if (ice(sc,ic,mi)%th(sk_z) .gt. sk_th_max) then
					! case where skeletal layer is already over the max
						tmp1 = ice(sc,ic,mi)%th(sk_z)/sk_th_max ! concentrate (don't lose algae)
						tmp2 = sk_th_max
						dz_sk = min(0.,dz_sk)
					elseif((ice(sc,ic,mi)%th(sk_z)+dz_sk) .gt. sk_th_max) then
					! case where dz_sk puts sk_th_new over the max
						 dz_sk = sk_th_max - ice(sc,ic,mi)%th(sk_z)
					endif
					sk_th_new = tmp2 + dz_sk


            ! algal migration - only melt - automatically retained/diluted above during growth
        tmp3 = dhdt_cm_s*86400.
        if ((alg_mig .gt. 0) .and. (dz_sk .lt. 0.)) then

            ! reconstructing sk_gl, since sk_gl is global and we may be updating skeletal
            ! more than once, and we don't want to grow/melt extra times
            z2 = dz_sk * DBLE(z_sk)

            tmp2 = abs(tmp3)
            if (alg_mig .eq. 1) then

              ! Algae can move at speed of ice bottom growth/melt up to alg_mig_crit speed
              if (tmp2 .gt. alg_mig_crit) then
                  tmp2 = z2 * alg_mig_crit/tmp2
              else
                  tmp2 = z2
              endif
            elseif (alg_mig .eq. 2) then
              tmp2 = abs(tmp3)
              if (tmp2 .gt. alg_mig_crit) then
                  tmp2 = 0.d0  !
              else
                  tmp2 = z2*(1.-tmp2/alg_mig_crit)
              endif
            endif
            dz = DBLE(z_sk)*ice(sc,ic,mi)%th(sk_1)
            z1 = (tmp2 + dz) / dz ! fractional change
            print *,'sk_gl - concentrating',z2,tmp3
        else
            z1 = 1.0d0
        endif
        print *,'tmp1,z1',tmp1,z1
        do ii=sk_1,sk_z
           smalg_pre_map(ii) = ice(sc,ic,mi)%smalg(ii) * tmp1 / z1
        enddo


					! maintain constant skeletal layer physics
					do ii=sk_1,sk_z

						! skeletal physics remain constant
						s_mean = f(mi)%s*c_5
						t_mean = min(f(mi)%t,f(mi)%s*mu)

#include "sia2_env_brine_calc.inc.f90"

						ice(sc,ic,mi)%t(ii) = t_mean
						ice(sc,ic,mi)%s(ii) = s_mean
						ice(sc,ic,mi)%d(ii) = d_mean
						ice(sc,ic,mi)%bs(ii) = bs_mean
						ice(sc,ic,mi)%bd(ii) = bd_mean
						ice(sc,ic,mi)%bv(ii) = bv_mean
						ice(sc,ic,mi)%heat(ii) = heat_mean

						! assuming tracers are in bulk format, so we don't have to
						! change values
						! skeletal layer may have somehow gotten above max (probably through advection)
						! don't loose any skeletal algae - concentrate in new depth

						ice(sc,ic,mi)%smalg(ii) = smalg_pre_map(ii)


!						if (ice(sc,ic,mi)%no3(ii) .gt. f(mi)%no3) then
!							print *,'Warning: no3 concentration higher than seawater in skeletal layer'
!							print *,'mi: ',mi,'layer: ',ii,'in sia2_env_update skeletal line 50'
!						endif

						! add seawater to fill in growing skeletal layer
						if (dz_sk .gt. 0.) then
							bv_mean = ice(sc,ic,mi)%bv(ii)*c_001
							ice(sc,ic,mi)%no3(ii) = (f(mi)%no3*dz_sk*bs_mean/f(mi)%s*bv_mean &
								+ ice(sc,ic,mi)%no3(ii)*ice(sc,ic,mi)%th(ii))/sk_th_new
							ice(sc,ic,mi)%nh4(ii) = (f(mi)%nh4*dz_sk*bs_mean/f(mi)%s*bv_mean &
								+ ice(sc,ic,mi)%nh4(ii)*ice(sc,ic,mi)%th(ii))/sk_th_new
							ice(sc,ic,mi)%po4(ii) = (f(mi)%po4*dz_sk*bs_mean/f(mi)%s*bv_mean &
								+ ice(sc,ic,mi)%po4(ii)*ice(sc,ic,mi)%th(ii))/sk_th_new
							ice(sc,ic,mi)%sioh4(ii) = (f(mi)%sioh4*dz_sk*bs_mean/f(mi)%s*bv_mean &
								+ ice(sc,ic,mi)%sioh4(ii)*ice(sc,ic,mi)%th(ii))/sk_th_new
							! particulates below are concentrated/sticky, not multiplied by bv_mean
							ice(sc,ic,mi)%smalg(ii) = (alg_wc*dz_sk*bs_mean/f(mi)%s &
								+ ice(sc,ic,mi)%smalg(ii)*ice(sc,ic,mi)%th(ii))/sk_th_new
							ice(sc,ic,mi)%poc(ii) = (f(mi)%poc*dz_sk*ice(sc,ic,mi)%bs(ii)/f(mi)%s &
								+ ice(sc,ic,mi)%poc(ii)*ice(sc,ic,mi)%th(ii))/sk_th_new
						endif

						ice(sc,ic,mi)%th(ii) = sk_th_new
						ice(sc,ic,mi)%id(ii) = ice(sc,ic,mi)%id(ii-1) + sk_th_new

!						if (ice(sc,ic,mi)%no3(ii) .gt. f(mi)%no3) then
!							print *,'Warning: no3 concentration higher than seawater in skeletal layer'
!							print *,'mi: ',mi,'layer: ',ii,'in sia2_env_update skeletal line 84'
!						endif

					 enddo
