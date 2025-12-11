

      ! test to see if we are doing any of this stuff
!      if ((ice(sc,ic,mi)%pond%th .gt. 0.) .or. (ice(sc,ic,mi)%id(int_z) .ge. 0.2)) then

          ! prep vars for routines
          vb_open = 0
          cv_open = 0
          visc = 1.79 ! g/m/s - need equation for viscosity as a func. of temp/salinity
          dsdt_cw = 0. ! vector assignment
          min_perm = c9999

          do ii=1,sk_z
              ! find lowermost layer open to brine convection
              if (ice(sc,ic,mi)%bv(ii) .lt. vb_crit .and. ii .le. int_z) then
                  vb_open = ii
                  perm(ii) = 0.
              else
                  !perm(ii) = (3.0e-8)*(ice(sc,ic,mi)%bv(ii)*c_001)**3
                  !perm(ii) = (3.0e-8)*((ice(sc,ic,mi)%bv(ii)-vb_crit)*c_001)**2
                  perm(ii) = (1.d-17)*(ice(sc,ic,mi)%bv(ii))**3.1
                  !perm(ii) = min(perm(ii),(3.0e-8)*((ice(sc,ic,mi)%bv(ii)-vb_crit)*c_001)**2)
              endif

              ! find minimum permeability
              min_perm = min(min_perm,perm(ii))

			  ! return tracers to brine conc.
              tmp1 = 1000./ice(sc,ic,mi)%bv(ii)
              ice(sc,ic,mi)%smalg(ii) = ice(sc,ic,mi)%smalg(ii)*tmp1
              ice(sc,ic,mi)%poc(ii) = ice(sc,ic,mi)%poc(ii)*tmp1
			  ice(sc,ic,mi)%no3(ii) = ice(sc,ic,mi)%no3(ii)*tmp1
			  ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii)*tmp1
			  ice(sc,ic,mi)%po4(ii) = ice(sc,ic,mi)%po4(ii)*tmp1
			  ice(sc,ic,mi)%sioh4(ii) = ice(sc,ic,mi)%sioh4(ii)*tmp1

              ! collapse layers to bv thickness for pond drainage
              tmp1 = ice(sc,ic,mi)%bv(ii)*c_001
              th_new(ii) = ice(sc,ic,mi)%th(ii)*tmp1
          enddo

#include "sia2_env_pond_mix3.inc.f90"

          ! pond draining take precedence of cox/week desal, so here
          ! we prodvide updated tracer/brine salinity/salinity vars for cox/weeks
          ! to operate with
          do ii=1,sk_z
              if ((ice(sc,ic,mi)%s(ii) + dsdt(ii)) .lt. min_sal) then
                  dsdt(ii) = min_sal - ice(sc,ic,mi)%s(ii)
              endif
              s_new(ii) = ice(sc,ic,mi)%s(ii) + dsdt(ii)

              ! find new parameters of convected brine
              bs_new(ii) = ice(sc,ic,mi)%bs(ii)*(s_new(ii)/ice(sc,ic,mi)%s(ii))
              bd_new(ii) = 1000000. + bs_new(ii)*800.       ! g/m^3
              bv_mean = ice(sc,ic,mi)%bv(ii)*c_001 ! brine volume doesn't change until after heat flux
              d_new(ii) = ice(sc,ic,mi)%d(ii) + bv_mean*(bd_new(ii) - ice(sc,ic,mi)%bd(ii))
              ! new brine volume again after desalination...
              d_mean = (1.-bb_f)*IceD*bd_new(ii)*bs_new(ii)/(bd_new(ii)*bs_new(ii) - &
                  s_new(ii)*(bd_new(ii) - IceD))
              bv_new(ii) = 1000.*(d_mean*s_new(ii))/(bd_new(ii)*bs_new(ii))  ! used to calc nutrient flux in cox/weeks

              smalg_new(ii) = ice(sc,ic,mi)%smalg(ii) + dsmalg(ii)
              poc_new(ii) = ice(sc,ic,mi)%poc(ii) + dpoc(ii)
              no3_new(ii) = ice(sc,ic,mi)%no3(ii) + dno3(ii)
              nh4_new(ii) = ice(sc,ic,mi)%nh4(ii) + dnh4(ii)
              po4_new(ii) = ice(sc,ic,mi)%po4(ii) + dpo4(ii)
              sioh4_new(ii) = ice(sc,ic,mi)%sioh4(ii) + dsioh4(ii)
          enddo

#include "sia2_env_desal_cw4.inc.f90"

          do ii=1,sk_z

              bv_mean = ice(sc,ic,mi)%bv(ii)*c_001

              ! assign s_new and d_new for consideration in heat flux
              if (ii .lt. sk_1) then
                  ! update salinity and density (common to all layers)
                  ! bound salinity at sal_min, don't gravity drain upward...
                  dsdt(ii) = dsdt(ii) + dsdt_cw(ii)
                  if ((ice(sc,ic,mi)%s(ii) + dsdt(ii)) .lt. min_sal) then
                      dsdt(ii) = min_sal - ice(sc,ic,mi)%s(ii)
                  endif
                  s_new(ii) = ice(sc,ic,mi)%s(ii) + dsdt(ii)
                  ! find new parameters of convected brine
                  bs_mean = ice(sc,ic,mi)%bs(ii)*(s_new(ii)/ice(sc,ic,mi)%s(ii))
                  bd_mean = 1000000. + bs_mean*800.       ! g/m^3
                  d_new(ii) = ice(sc,ic,mi)%d(ii) + bv_mean*(bd_mean - ice(sc,ic,mi)%bd(ii))
              endif

 			  ! bulk tracer conc.
              ice(sc,ic,mi)%smalg(ii) = ice(sc,ic,mi)%smalg(ii)*bv_mean
              ice(sc,ic,mi)%poc(ii) = ice(sc,ic,mi)%poc(ii)*bv_mean
			  ice(sc,ic,mi)%no3(ii) = ice(sc,ic,mi)%no3(ii)*bv_mean
			  ice(sc,ic,mi)%nh4(ii) = ice(sc,ic,mi)%nh4(ii)*bv_mean
			  ice(sc,ic,mi)%po4(ii) = ice(sc,ic,mi)%po4(ii)*bv_mean
			  ice(sc,ic,mi)%sioh4(ii) = ice(sc,ic,mi)%sioh4(ii)*bv_mean

              ! reduce b_flux by bfrf
              b_flux(ii) = b_flux(ii)/dble(bfrf)

          enddo
!    endif
