

                   ! convert pond tracers to total mass
                   th_pond_new = ice(sc,ic,mi)%pond%th
                   s_pond_new = ice(sc,ic,mi)%pond%s*ice(sc,ic,mi)%pond%th
                   h_pond_new = ice(sc,ic,mi)%pond%heat*ice(sc,ic,mi)%pond%th
                   smalg_pond_new = ice(sc,ic,mi)%pond%smalg*ice(sc,ic,mi)%pond%th
                   poc_pond_new = ice(sc,ic,mi)%pond%poc*ice(sc,ic,mi)%pond%th
                   no3_pond_new = ice(sc,ic,mi)%pond%no3*ice(sc,ic,mi)%pond%th
                   nh4_pond_new = ice(sc,ic,mi)%pond%nh4*ice(sc,ic,mi)%pond%th
                   po4_pond_new = ice(sc,ic,mi)%pond%po4*ice(sc,ic,mi)%pond%th
                   sioh4_pond_new = ice(sc,ic,mi)%pond%sioh4*ice(sc,ic,mi)%pond%th

                   if (snow_melt_now .gt. 0.) then

                       ! assume that snow_gl is less than one layer, so only top layer can melt
                       ! should change this to be more flexible some time and allow
                       ! more that one snow layer to melt at once
                       z_snow = ice(sc,ic,mi)%snow%z
                       dz = snow_melt_now*ice(sc,ic,mi)%snow%d(z_snow)*1.e-6*0.918
                       th_pond_new = th_pond_new + dz
                       snow_melt_now = 0.

                       ! don't add any tracers or salt.  The heat added
                       ! is zero b/c the heat of fresh water at 0degC
                       ! is zero
                   else
                       if (melt_drain .gt. 0.d0) then
                           do jjj=1,z_snow_new
                               bv_mean = ice(sc,ic,mi)%bv(jjj)*c_001
                               dz = ice(sc,ic,mi)%th(jjj)*bv_mean
                               z1 = dz/bv_mean
                               th_pond_new = th_pond_new + dz
                               ice(sc,ic,mi)%pond%perc = ice(sc,ic,mi)%pond%perc + dz
                               s_pond_new = s_pond_new + ice(sc,ic,mi)%bs(jjj)*dz
                               h_pond_new = h_pond_new + ice(sc,ic,mi)%t(jjj)*cw*dz
                               no3_pond_new = no3_pond_new + ice(sc,ic,mi)%no3(jjj)*z1
                               nh4_pond_new = nh4_pond_new + ice(sc,ic,mi)%nh4(jjj)*z1
                               po4_pond_new = po4_pond_new + ice(sc,ic,mi)%po4(jjj)*z1
                               sioh4_pond_new = sioh4_pond_new + ice(sc,ic,mi)%sioh4(jjj)*z1
                               smalg_pond_new = smalg_pond_new + smalg_pre_map(jjj)*z1
                               poc_pond_new = poc_pond_new + ice(sc,ic,mi)%poc(jjj)*z1
                           enddo

                       endif
                       tmp1 = abs(s_gl) - melt_drain
                       if (tmp1 .gt. 0.d0) then
                           ! we assume here that jj (the layer where new ice top is currently found)
                           ! has already been detected.
                           jjj = 1
                           do while (tmp1 .gt. 0.d0)
                               dz = tmp1
                               if (dz .gt. ice(sc,ic,mi)%th(jjj)) then
                                   dz = ice(sc,ic,mi)%th(jjj)
                               endif
                               th_pond_new = th_pond_new + dz
                               ice(sc,ic,mi)%pond%perc = ice(sc,ic,mi)%pond%perc + dz
                               s_pond_new = s_pond_new + ice(sc,ic,mi)%s(jjj)*dz
                               tmp2 = 1000000.+800.*ice(sc,ic,mi)%s(jjj) ! find density of melted ice
                               h_pond_new = h_pond_new + ice(sc,ic,mi)%s(jjj)*mu*cw*tmp2*dz
                               no3_pond_new = no3_pond_new + ice(sc,ic,mi)%no3(jjj)*dz
                               nh4_pond_new = nh4_pond_new + ice(sc,ic,mi)%nh4(jjj)*dz
                               po4_pond_new = po4_pond_new + ice(sc,ic,mi)%po4(jjj)*dz
                               sioh4_pond_new = sioh4_pond_new + ice(sc,ic,mi)%sioh4(jjj)*dz
                               smalg_pond_new = smalg_pond_new + smalg_pre_map(jjj)*dz
                               poc_pond_new = poc_pond_new + ice(sc,ic,mi)%poc(jjj)*dz

                               tmp1 = tmp1 - dz
                               jjj = jjj+1
                           enddo
                       endif
                   endif
				   ! restore and integrate new tracers
				   ice(sc,ic,mi)%pond%th = th_pond_new
				   ice(sc,ic,mi)%pond%s = s_pond_new/th_pond_new
			       ice(sc,ic,mi)%pond%d = 1000000.+800.*ice(sc,ic,mi)%pond%s  ! find density of melt pond water
				   ice(sc,ic,mi)%pond%heat = h_pond_new/th_pond_new
				   ice(sc,ic,mi)%pond%t = ice(sc,ic,mi)%pond%heat/ice(sc,ic,mi)%pond%d/cw
                   ice(sc,ic,mi)%pond%smalg = smalg_pond_new/th_pond_new
                   ice(sc,ic,mi)%pond%poc = poc_pond_new/th_pond_new
                   ice(sc,ic,mi)%pond%no3 = no3_pond_new/th_pond_new
                   ice(sc,ic,mi)%pond%nh4 = nh4_pond_new/th_pond_new
                   ice(sc,ic,mi)%pond%po4 = po4_pond_new/th_pond_new
                   ice(sc,ic,mi)%pond%sioh4 = sioh4_pond_new/th_pond_new


                   !print *, ice(sc,ic,mi)%pond%th,ice(sc,ic,mi)%pond%s,ice(sc,ic,mi)%pond%t,ice(sc,ic,mi)%pond%no3
                   !if (ice(sc,ic,mi)%pond%no3 .gt. 5. .or. ice(sc,ic,mi)%bv(41) .lt. 10.) then
                   !    testvar = -1
                   !endif