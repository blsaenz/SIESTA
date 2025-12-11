      if (kevin_light .gt. 0) then
          kl_PARID = 0.
          kl_PARD = 0.
          kl_PARG = 0.
          zenith_deg = 180.*zenith_angle/pi
          elevation_deg = 90. - zenith_deg


          if (elevation_deg .gt. 0.) then

			  do ii=1,wavl

					lambda = dble(ii)*10.+390.  ! nm
					jj = ii*10-9                ! index to 1-nm absorption arrays

					!----Transmittance due to Rayleigh scattering
					kl_M=1./(cos(zenith_angle) + 0.15*(93.885-zenith_angle)**(-1.253))    ! relative air mass
					kl_Tr=exp(-.0088*kl_M*(1./(lambda*c_001)**4))    ! Transmittance due to

					!----Transmittance due to aerosol extinction
					kl_alpha = 0.6
					kl_Bn = 0.07
					kl_Ta=EXP(-1.*kl_Bn*kl_M*(lambda*c_001)**(-1.*kl_alpha))    !Transmittance due to

					!----Transmittance due to water vapor
					kl_w=1.42                !cm of precipitable water
					kl_Tw=exp(-0.2385*awater_v(jj)*kl_w*kl_M/(1+20.07*kl_w*awater_v(jj)*kl_M)**0.45)

					!----Transmittance due to ozone absorption
					kl_Moz=1.0035/(cos(zenith_angle)*cos(zenith_angle)+0.007)**0.5
					kl_Ozone=.344            !cm  default
					kl_T0=(exp(-aozone(jj)*kl_Ozone*kl_Moz))

					!----Transmittance due to gas absorption--this is 100% in the visible
					kl_au=0.
					kl_Tu=exp(-1.41*kl_au*kl_M/(1.+118.93*kl_au*kl_M)**0.45)

					kl_Id=surfacelight(jj)*kl_Tr*kl_Ta*kl_T0*kl_Tw*kl_Tu*COS(zenith_angle)  !direct light in units of W/m2/µm
					kl_PARID=kl_PARID+kl_Id*10.   !Id integrated from 400-700 nm

					!Calculation of diffuse radiation
					kl_Dr=surfacelight(jj)*COS(zenith_angle)*kl_T0*kl_Tu*kl_Tw*kl_Ta*(1-kl_Tr)*0.5
					kl_Da=surfacelight(jj)*COS(zenith_angle)*kl_T0*kl_Tu*kl_Tw*kl_Tr*(1-kl_Ta)*0.85*1.
					kl_v=(exp(aozone(jj)*kl_Ozone*1.9))
					kl_w=1.
					kl_x=1.
					kl_y=EXP(-.0088*1.9*(1./lambda**4))
					kl_z=EXP(-kl_Bn*1.9*(lambda*c_001)**(-1.*kl_alpha))
					kl_AtmAlbedo=kl_v*kl_w*kl_x*((1-kl_y)*kl_z*0.5+0.22*(1.-kl_z)*1.*kl_y)
					kl_Dm=(kl_Id+kl_Dr+kl_Da)*0.8*kl_AtmAlbedo/(1.-.8*kl_AtmAlbedo)
					kl_d=kl_Da+kl_Dr+kl_Dm                                 !in units of W/m2/µm
					kl_PARD=kl_PARD+kl_d*10.                  ! D integrated from 400-700 nm
					kl_g=kl_Id+kl_d
					kl_PARG=kl_PARID+kl_PARD
					kl_cc = 10.**(-0.099*f(mi)%fc*c_01)
					kl_par(ii)=kl_g*kl_cc

					!print *,lambda,'zenith_angle',zenith_angle,mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii),mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii)
					mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii) = kl_Id*kl_cc
					mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii) = kl_d*kl_cc
					!print *,mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii),mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii)
			  enddo
          else
			  do ii=1,wavl

					mp_f%Edir(f_index)%Ed(i_mp,j_mp,ii) = 0.
					mp_f%Edif(f_index)%Ed(i_mp,j_mp,ii) = 0.
			  enddo
		  endif

      endif ! end of kevin light
