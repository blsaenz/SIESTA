
          ! adjust snow height
          ! -----------------------------------------------------------
          tmp1 = sh_increment(sc,ic)*dtt_s
          snow_gl = snow_gl + tmp1


          ! FIND DHDT bottom
          ! -----------------------------------------------------------
          ! first we calculate dhdt without division by density - we don't know where
          ! to put the new growth yet
          dhdt_g_m2 = dtt_s*(Fc_bot - f(mi)%fw)
          dhdt_cm_s = 0.d0

          !print *,'Bottom Flux: ',Fc_bot - f(mi)%fw,dhdt_g_m2

          sk_cur = 0.d0
          do ii = sk_1,sk_z
              sk_cur = sk_cur + ice(sc,ic,mi)%th(ii)  ! find total current skeletal thickness
          enddo

          ! growing ice - allocate to either skeletal or internal ice
          if (dhdt_g_m2 .gt. 0.d0) then
              ! add to skeletal layer first, if it isn't fully formed (i.e. at sk_h_max)
              if ((sk_cur + sk_gl) .lt. sk_h_max) then
                  T_mean = f(mi)%s*mu
                  T_melt = T_melt*c_5  ! assume skeletal is 1/2 seawater salinity
                  S_mean = f(mi)%s*c_5
                    bs_mean = f(mi)%s
                    bd_mean = 1.d6 + bs_mean*800.d0
                    d_mean = (1.d0-bb_f)*IceD*bd_mean*bs_mean/(bd_mean*bs_mean - S_mean* &
                  (bd_mean - IceD))

                    q_melt = c0*(T_melt - T_mean) + Lf*(1. - T_melt/T_mean) ! heat w/out mixed layer temp factored in
                    q_melt = q_melt + cw*(max(c0,(f(mi)%t-T_melt))) ! add heat required to reach t_melt, if any
                    q_melt = d_mean*q_melt

                  if (q_melt .gt. 0.d0) then

                        dhdt = dhdt_g_m2/q_melt
                        if ((sk_cur + dhdt + sk_gl) .gt. sk_h_max) then
                            dhdt = sk_h_max - sk_cur - sk_gl
                            dhdt_g_m2 = dhdt_g_m2 - q_melt*dhdt ! remove energy used to grow skeletal layer
                        else
                            dhdt_g_m2 = 0.d0
                        endif
                        sk_gl = sk_gl + dhdt
                        dhdt_cm_s = dhdt_cm_s + 100.*dhdt/dtt_s

                    endif

                  ! salt and agua accounting
                  ! converting to kg/m^2 freshwater input 1000kg/m3 * m * fraction_water by weight
                  !m(mi)%fresh_in = m(mi)%fresh_in + dhdt*ice(sc,ic,mi)%af*1.d3*(1.d0-f(mi)%s*c_001)
                  ! converting to kg/m^2 salt input 1000kg/m3 * m * fraction_salt by weight
                  !m(mi)%dsdt_in = m(mi)%dsdt_in + dhdt*ice(sc,ic,mi)%af*f(mi)%s

                  ! converting to kg/m^2 salt input 1000kg/m3 * m * fraction_salt by weight
                  !m(mi)%dsdt_out = m(mi)%dsdt_out + dhdt*ice(sc,ic,mi)%af*(f(mi)%s*0.5d0)

              endif
              ! add to congelation ice
              if (dhdt_g_m2 .gt. 0.d0) then
          ! determine conductive heat flux at boundary (W/m^2)
          F0 = dhdt_g_m2/dtt_s
            ! use stable salinity regression to figure out dhdt
          dhdt_cm_s = 9.873e-12*F0**3 - 6.621e-10*F0**2 + 4.474e-7*F0  ! (cm/s) dhdt vs. F0 regression at T=-1.8 freezing temp, petrich keff
          ! back out stable salinity from dhdt_cm_s (petrich 2006 keff)
                  keff = 0.19d0*(dhdt_cm_s*7.4074d4)**(0.46)
                  keff=max(keff,0.12)

!          dhdt = dhdt_cm_s*c_01*dtt_s
!                    ! update congelation growth layer parameters
!                    c_gl_sal = (c_gl_sal*c_gl + f(mi)%s*keff*dhdt)/ &
!                        (c_gl + dhdt)
!                    c_gl = c_gl + dhdt

                    ! predicted new ice properties
                  T_mean = f(mi)%s*mu ! temp required to grow ice in seawater of salinity f(mi)%s
                    s_mean = f(mi)%s*keff
                    bs_mean = f(mi)%s
                    bd_mean = 1.d6 + bs_mean*800.d0
                    d_mean = (1.d0-bb_f)*IceD*bd_mean*bs_mean/(bd_mean*bs_mean - S_mean* &
                  (bd_mean - IceD))
                  T_melt = mu*s_mean

                    q_melt = c0*(T_melt - T_mean) + Lf*(1. - T_melt/T_mean) ! heat w/out mixed layer temp factored in
                    q_melt = q_melt + cw*(max(c0,(f(mi)%t-T_melt))) ! add heat required to reach t_melt, if any
                    q_melt = d_mean*q_melt

                    dhdt = (F0*dtt_s)/q_melt
                    c_gl_sal = (c_gl_sal*c_gl + s_mean*dhdt)/(c_gl + dhdt)
                    c_gl = c_gl + dhdt

              endif
          else ! melting ice
              ! if skeletal layer is present, melt it first
              if ((sk_cur + sk_gl) .gt. sk_h_min) then
                  t_melt = ice(sc,ic,mi)%s(sk_z)*mu
                  q_melt = ice(sc,ic,mi)%d(sk_z)*(c0*(T_melt - ice(sc,ic,mi)%t(sk_z)) &
                    + Lf*(1. - T_melt/ice(sc,ic,mi)%t(sk_z)))
                    if (q_melt .gt. 0.d0) then
                        dhdt = dhdt_g_m2/q_melt
                        if ((sk_cur + sk_gl + dhdt) .lt. sk_h_min) then
                            dhdt = sk_h_min - sk_cur - sk_gl
                            dhdt_g_m2 = dhdt_g_m2 - dhdt*q_melt
                        else
                            dhdt_g_m2 = 0.d0
                        endif
                        sk_gl = sk_gl + dhdt
                        print *,'sk_gl - ',sk_gl
                        dhdt_cm_s = dhdt_cm_s + 100.d0*dhdt/dtt_s
                        print *,'dhdt_cm_s - ',dhdt_cm_s
                    endif

              endif
              ! melt congelation ice
              tmp1 = ice(sc,ic,mi)%id(int_z)
              if ((dhdt_g_m2 .lt. 0.) .and. (tmp1 .gt. h_min)) then
                  t_melt = ice(sc,ic,mi)%s(int_z)*mu
                  t_mean = min(t_melt,ice(sc,ic,mi)%t(int_z)) ! can't melt ice that has already melted
                  q_melt = ice(sc,ic,mi)%d(int_z)*(c0*(T_melt - t_mean) &
                    + Lf*(1. - T_melt/t_mean))
                  if (q_melt .gt. 0.d0) then
                      dhdt = dhdt_g_m2/q_melt
                      if (abs(dhdt) .gt. gl_max) then
                          print *,'melting exceeded gl_max - dhdt: ',dhdt,dtt_s,q_melt
                      endif
!                      dhdt = max(-1*gl_max,dhdt)  ! prevent melting too fast --
                      if ((tmp1 + c_gl + dhdt) .lt. h_min) then
                              dhdt = h_min - tmp1 - c_gl
                              dhdt_g_m2 = dhdt_g_m2 - dhdt*q_melt
                      else
                              dhdt_g_m2 = 0.
                      endif
                      dhdt_cm_s = dhdt_cm_s + 100.*dhdt/dtt_s
                      c_gl = c_gl + dhdt
                      print *,'cgl - ',c_gl
                      print *,'dhdt_cm_s - ',dhdt_cm_s
                  else
                          dhdt = 0.
                      print *,'are we really here?'
                  endif
            endif
        endif

          ! subliming surface snow or ice
          ! ------------------------------------------------------------
          dhdt_g_m2 = dtt_s*(Fe)
          if (dhdt_g_m2 .lt. 0. .and. snow_model_temp .eq. 1) then
              if (ice(sc,ic,mi)%snow%depth .gt. 0.d0) then
                  if (ice(sc,ic,mi)%snow%z .gt. 0) then
                      q_melt = ice(sc,ic,mi)%snow%d(ice(sc,ic,mi)%snow%z)
                  else
                      q_melt = 4.d5
                  endif
                  !if (ice(sc,ic,mi)%snow%ts .lt. 0.d0) then
                      q_melt = q_melt*(Lv + Lf)  ! these were reversed before
                  !else
                  !    q_melt = q_melt*Lv ! these were reversed before
                  !endif

                  dhdt = dhdt_g_m2/q_melt ! snow sublimation height
                  if ((snowh1m + snow_gl + dhdt) .lt. 0.) then
                      dhdt_g_m2 = (ice(sc,ic,mi)%snow%depth+snow_gl+dhdt)*q_melt
                      snow_gl = -1.0*ice(sc,ic,mi)%snow%depth
                  else
                      dhdt_g_m2 = 0.d0
                      snow_gl = snow_gl + dhdt
                  endif
              endif
          endif

          if (ice(sc,ic,mi)%snow%depth .eq. 0.d0)  then
              tmp1 = ice(sc,ic,mi)%id(int_z)
              if ((dhdt_g_m2 .lt. 0.) .and. (tmp1 .gt. h_min)) then
                  t_melt = ice(sc,ic,mi)%s(1)*mu
                  t_mean = min(t_melt,ice(sc,ic,mi)%t(1)) ! can't sublime ice that has already melted
                  q_melt = (c0*(T_melt - t_mean) &
                    + Lf*(1. - T_melt/t_mean))
                  q_melt = ice(sc,ic,mi)%d(1)*(q_melt + Lv)   ! add latent heat of vaporization
                  if (q_melt .gt. 0) then
                      dhdt = dhdt_g_m2/q_melt
                      if (abs(dhdt) .gt. gl_max) then
                          print *,'subliming exceeded gl_max - dhdt: ',dhdt,dtt_s,q_melt,Fe
                      endif
                      dhdt = max(-1*gl_max,dhdt)  ! prevent melting too fast --
                      if ((tmp1 + c_gl + dhdt) .lt. h_min) then
                          dhdt = h_min - tmp1 - c_gl
                          dhdt_g_m2 = dhdt_g_m2 - dhdt*q_melt
                      else
                          dhdt_g_m2 = 0.
                      endif
                      s_gl = s_gl + dhdt
                  else
                      dhdt = 0.
                  endif
              endif
           endif

           ! FIND DHDT top
           ! -----------------------------------------------------------
           z_snow = ice(sc,ic,mi)%snow%z
           T_melt = ice(sc,ic,mi)%s(1)*mu

!                       print *,'ice(sc,ic,mi)%snow%ts',ice(sc,ic,mi)%snow%ts
!                       print *,'ts_snow',ts_snow
!                       print *,'ki(1)',ki(1)
!                       print *,'Ts_last',Ts_last
!                       print *,'Ts_next',Ts_last
!                       print *,'T_prev(1)',T_prev(1)
!                       print *,'T_last(1)',T_last(1)
!                       print *,'T_next(1)',T_next(1)


           if ((z_snow .eq. 0) .and. (ice(sc,ic,mi)%snow%ts .ge. T_melt)) then
               ! don't 'melt away ice to nothing' if more realistic brine movement is occuring
!                   if (desal .lt. 2) then
                   ! find total surface flux over bare ice
                    if (atmo .eq. 1) then
                        Fm0 = Fr + FL + F_little_l + F_s + Fe
                        Fc_top = ki(1)*(Ts_last-T_last(1)) &
                            /dth(1)
                    else
                     Ts_kelvin = ice(sc,ic,mi)%snow%ts + kelvin0
                     Fm0 = F0_constants + F0_4*Ts_kelvin**4 + F0_3*Ts_kelvin**3 + &
                        F0_2*Ts_kelvin**2 + F0_1*Ts_kelvin
                     Fc_top = 2.*ki(1)*(ice(sc,ic,mi)%snow%ts-(T_last(1)+ice(sc,ic,mi)%t(1))*c_5) &
                        /ice(sc,ic,mi)%th(1)
                    endif
                   F0 = Fm0 - Fc_top
                   if (F0 .gt. 0.) then
                       t_melt = ice(sc,ic,mi)%s(1)*mu
!                           if (t_melt .ge. ice(sc,ic,mi)%t(1)) then
!                               s_gl = s_gl - ice(sc,ic,mi)%th(1)
!                           else
                           q_melt = ice(sc,ic,mi)%d(1)*(c0*(T_melt - ice(sc,ic,mi)%t(1)) &
                             + Lf*(1. - T_melt/ice(sc,ic,mi)%t(1)))
                           if (q_melt .gt. 0) then ! protects against surface melt going backward, in case ice is at h_min
                               s_gl = s_gl - dtt_s*(F0)/q_melt
!                                   if ((ice(sc,ic,mi)%id(int_z) + c_gl + s_gl) .lt. h_min) then
!                                       s_gl = 0.  ! don't melt if we're going below h_min
!                                   endif
                           endif
!                           endif
                   endif
!                   endif
           elseif ((z_snow .gt. 0) .and. (ice(sc,ic,mi)%snow%ts .ge. 0.)) then
               if ((snow_model_temp .eq. 1) .or. (snow_model_temp .eq. 3)) then
                   ! find total surface flux
                    if (atmo .eq. 1) then
                        Fm0 = Fr + FL + F_little_l + F_s + Fe
                        Fc_top = ki(1)*(Ts_last-T_last(1)) &
                            /dth(1)
                    else
                       Ts_kelvin = ice(sc,ic,mi)%snow%ts + kelvin0
                       Fm0 = F0_constants + F0_4*Ts_kelvin**4 + F0_3*Ts_kelvin**3 + &
                            F0_2*Ts_kelvin**2 + F0_1*Ts_kelvin - Fc_snow
                       Fc_top = 2.*ki(1)*(ice(sc,ic,mi)%snow%ts - &
                       (t_snow_prev(z_snow) + ice(sc,ic,mi)%snow%t(z_snow))*c_5) &
                            /ice(sc,ic,mi)%snow%th(z_snow)
                    endif
                   F0 = Fm0 - Fc_top
                   ! if positive, F0 must melt snow
                   if (F0 .gt. 0) then
                       ! amount melted
                       tmp1 = dtt_s*(F0)/(Lf*ice(sc,ic,mi)%snow%d(z_snow))
                       if (use_ponds .eq. 1) then
                           snow_melt = snow_melt + tmp1
                       endif
                       snow_gl = snow_gl - tmp1
                   endif
!                       print *,'Fm0',Fm0
!                       print *,'Fc_top',Fc_top
!                       print *,'F0',F0
!                       print *,'tmp1 - melt thickness',tmp1
!                       print *,'snow_gl',snow_gl
               elseif ((snow_model_temp .eq. 0) .and. (sh_increment(sc,ic) .lt. 0.)) then
                   ! specificed snow has decreased under melt conditions - add to snowmelt
                   ! snow_gl is already taken care under snow routine above
                   if (use_ponds .eq. 1) then
                       snow_melt = snow_melt - sh_increment(sc,ic)
                   endif
               endif
           endif

           ! if snow model 2, then determine snow melt based on snow_gl and snow temp
           if (snow_model_temp .eq. 2 .and. sh_increment(sc,ic) .lt. 0.d0) then
               if (ice(sc,ic,mi)%snow%ts .ge. 0.d0 .and. use_ponds .eq. 1) then
                   snow_melt = snow_melt - sh_increment(sc,ic)
               endif
           endif


           if (dt_step .eq. 139 .and. mi .eq. 29361 .and. ic .eq. 2 .and. sc .eq. 2) then
               testvar = -1
               !call exit (0)
           endif


#include "sia2_env_flood.inc.f90"