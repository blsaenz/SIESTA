
      ! surface matrix constant
      rs = -1.*F0 + dF0*Ts_last

      ! pre-calc top layer ki
      ki_down = (ice(sc,ic,mi)%d(1)/IceD)*(2.11-0.011*ice(sc,ic,mi)%t(1) + &
        0.09*ice(sc,ic,mi)%s(1)/ice(sc,ic,mi)%t(1) - (ice(sc,ic,mi)%d(1)-IceD)/1.e6)
      ki_down = max(ki_down,ki_min)

      ! top layer is SNOW
      ! ----------------------------------------------------------------
      if (z_snow .gt. 0) then

          ! mean temperature
          tmp1 = (T_last(1)+T_next(1))*c_5
          tmp1 = min(tmp1,0.)

          ! snow heat capacity averaged between T_last and T_next, * snow density
          tmp2 = 0.185 + 0.689e-2*(tmp1+kelvin0)
          tmp2 = max(0.02,tmp2)
          ci(1) = ice(sc,ic,mi)%snow%d(z_snow)*tmp2
          xl = 2.*dtt_s/((dth(1)+dth(2))*ci(1))

          ! assign constant vector - no irradiance absorbed in snow
          !r1 = T_last(1) + xl*Ed_W_snow(1)
          tmp3 = dtt_s/(ci(1)*ice(sc,ic,mi)%snow%th(z_snow))
          r1 = T_last(1) + Ed_W_snow(1)*tmp3

          ! determine of snow is at the melting point.  If so, toggle mo variable
          ! for different surface heat balance equation
          Fc_top = ki(1)*(Ts_last-T_last(1))/dth(1)
          Fm0 = F0-Fc_top
          if (snow_model .eq. 1 .or. snow_model .eq. 3) then
              if ((Ts_last .ge. 0.) .and. (Fm0 .gt. 0.)) then
                  ! melting - no need to calc new surface temp
                  mo = 0
                  Ts_next = 0.
              else
                  mo = 1
              endif
          else
              mo = 1
          endif

      ! top layer is ICE
      ! ----------------------------------------------------------------
      else

          d_mean = (d_new(1) + ice(sc,ic,mi)%d(1))*c_5
!          s_mean = (ice(sc,ic,mi)%s(1) + s_new(1))*c_5
          s_mean = s_new(1)
          bv_mean = bv_new(1)*c_001
          t_melt = s_mean*mu

          tmp1 = min(T_last(1),T_melt)
          tmp2 = min(T_next(1),T_melt)

          !ci(1) = (1.d0-bb_f)*d_mean*( &             ! remove air fraction
          !    c0 - Lf*T_melt/(tmp1*tmp2))   ! approximated heat capacity

          ci(1) = (1.d0-bb_f) * ( &             ! remove air fraction
              bd_new(1)*bv_mean*cw +   &    ! heat capacity of pure ice
              IceD*(1.d0-bv_mean)*c0 - &    ! heat capacity of brine
              IceD*Lf*T_melt/(tmp1*tmp2))   ! latent heat of freezing based on temp change


          xl = 2.*dtt_s/((dth(1)+dth(2))*ci(1))
          ! assign constant vector - include irradiance absorbed in ice
          tmp3 = dtt_s/(ci(1)*ice(sc,ic,mi)%th(1))
          r1 = T_last(1) + tmp3*Ed_W_temp(1) + xl*dcheat(1)

          ! look at Fc beforehand to see whether melting is occuring
          Fc_top = ki(1)*(Ts_last-T_last(1))/dth(1)
          Fm0 = F0-Fc_top
!                  if (desal .lt. 2) then
              if ((Ts_last .ge. t_melt) .and. (Fm0 .gt. 0.)) then
                  ! melting - no need to calc new surface temp
                  mo = 0
                  Ts_next = t_melt
              else
                  mo = 1
              endif
!                  else
!                      mo = 1
!                  endif
      endif

      if (ts_is_at .gt. 0) then
          Ts_last = f(mi)%at - kelvin0
          Ts_next = Ts_last
          mo = 0
      endif

      if (mo .eq. 1) then

          ! surface matrix coeffs.
          mat_sc = ki(1)/dth(1)
          mat_sb = dF0 - ki(1)/dth(1)

          mat_a = -1.*xl*ki(1)/dth(1)
          mat_c = -1.*xl*ki(2)/dth(2)
          mat_b = 1. - mat_a - mat_c

          ! assign constant vector
          r_vec(1) = rs
          r_vec(2) = r1

          ! assign top row of matrix
          DU(1) = mat_sc
          DC(1) = mat_sb

          ! assign second row of matrix
          DU(2) = mat_c
          DC(2) = mat_b
          DL(1) = mat_a

     else
          ! assign constant vector
          r_vec(1) = r1 + Ts_last*xl*ki(1)/dth(1) ! <-- not needed when Ts_last = 0

          ! assign top row of matrix
          mat_a = -1.*xl*ki(1)/dth(1)
          mat_c = -1.*xl*ki(2)/dth(2)
          mat_b = 1. - mat_a - mat_c

          DU(1) = mat_c
          DC(1) = mat_b

     endif

     do ii=2,ice_z

         ! COMPUTE INTERNAL LAYER MATRIX COEFFICIENTS
         ! ----------------------------------------------------
         if (ii .lt. ice_1) then
             tmp1 = (T_last(ii)+T_next(ii))*c_5
             tmp1 = min(tmp1,0.)
             ! snow heat capacity averaged between T_last and T_next, *ice_d
             tmp2 = 0.185 + 0.689e-2*(tmp1+kelvin0)
             tmp2 = max(0.02,tmp2)
             ci(ii) = ice(sc,ic,mi)%snow%d(z_snow-ii+1)*tmp2
             xl = 2.*dtt_s/((dth(ii)+dth(ii+1))*ci(ii))
             ! assign constant vector w/ absorbed irradiance
             !r_vec(ii+mo) = T_last(ii) + xl*Ed_W_snow(ii)

             tmp3 = dtt_s/(ci(ii)*ice(sc,ic,mi)%snow%th(z_snow-ii+1))
             r_vec(ii+mo) = T_last(ii) + tmp3*Ed_W_snow(ii)


         else

              ! heat flux dependent desalination
              ! ------------------------------------------------------------
              jj = ii-z_snow
              d_mean = (d_new(jj) + ice(sc,ic,mi)%d(jj))*c_5
    !          s_mean = (ice(sc,ic,mi)%s(jj) + s_new(jj))*c_5
              s_mean = s_new(jj)
              t_melt = s_mean*mu
              bv_mean = bv_new(jj)*c_001

              tmp1 = min(T_last(ii),T_melt)
              tmp2 = min(T_next(ii),T_melt)

              !ci(ii) = (1.d0-bb_f)*d_mean*( &       ! remove air fraction
              !    c0 - Lf*T_melt/(tmp1*tmp2))   ! approximated heat capacity


              ci(ii) = (1.d0-bb_f) * ( &            ! remove air fraction
                  bd_new(jj)*bv_mean*cw +  &    ! heat capacity of pure ice
                  IceD*(1.d0-bv_mean)*c0 - &    ! heat capacity of brine
                  IceD*Lf*T_melt/(tmp1*tmp2))   ! latent heat of freezing based on temp change

              xl = 2.*dtt_s/((dth(ii)+dth(ii+1))*ci(ii))
              tmp3 = dtt_s/(ci(ii)*ice(sc,ic,mi)%th(jj))

              ! assign constant vector - include irradiance absorbed in ice
               r_vec(ii+mo) = T_last(ii) &
!                + tmp3*(dcheat(jj) + Ed_W_temp(jj))
                + tmp3*Ed_W_temp(jj) + xl*dcheat(jj)

         endif

         ! assign row of matrix
         DU(ii+mo) = -1.*xl*ki(ii+1)/dth(ii+1) ! replaced above to reduce computation
         DL(ii-1+mo) = -1.*xl*ki(ii)/dth(ii)
         DC(ii+mo) = 1. - DU(ii+mo) - DL(ii-1+mo)

         ! ammend matrix values for final ice layer, when bottom temp is fixed
         if (ii .eq. ice_z) then

             r_vec(ice_z+mo) = r_vec(ice_z+mo) - DU(ii+mo)*(f(mi)%s*mu)
             !DU(ii+mo) = d0_

             ! add heat absorbed in skeletal layers to bottom layer
             ! to maintain conservation of heat
             do jj=sk_1,sk_z
                 r_vec(ice_z+mo) = r_vec(ice_z+mo) + tmp3*(Ed_W_temp(jj))
             enddo

          endif

     enddo


     ! SOLVE MATRIX PROBLEM FOR IMPLICIT HEAT CONDUCTION
     ! ---------------------------------------------------------
     m_order = int_z+z_snow+mo
     DC_calc = DC
     DU_calc = DU
     DL_calc = DL

!             print *,'Iteration: ',jjj,z_snow,sk_z,ice_z
!             print *,'heat_stop: ',heat_stop,int_z
!             print *,'T_last: ',t_last
!             print *,'T_next: ',t_next
!             print *,'DU: ',du
!             print *,'DC: ',dc
!             print *,'DL: ',dl
!             print *,'r_vec: ',r_vec

     call DGTSV(m_order,1,DL_calc,DC_calc,DU_calc,r_vec,m_order,info)

 !            print *,'info: ',info
 !            print *,'r_vec: ',r_vec


     ! triage new temps
     if (mo .eq. 1) then
         T_diff = abs(r_vec(1) - Ts_next)
         T_step = abs(r_vec(1) - Ts_last)
         Ts_next = r_vec(1)
     else
         T_diff = 0.
         T_step = 0.
         !Ts_next = 0.
     endif
     do ii=1,ice_z
         T_diff=max(T_diff,abs(r_vec(ii+mo)-T_next(ii)))
         T_step=max(T_step,abs(r_vec(ii+mo)-T_last(ii)))
         T_next(ii) = r_vec(ii+mo)
     enddo

   ! conductive heat flux at boundaries
   Fc_bot = ki(ice_z+1)*(f(mi)%s*mu-T_next(ice_z))/dth(ice_z+1) ! mean flux over sub-dt timestep
   Fc_top = ki(1)*(Ts_next-T_next(1))/dth(1)  ! mean flux over sub-dt timestep

