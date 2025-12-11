! Melt pond height - Notz 2009 using Darcy's Law
! ----------------------------------------------------------------------
! z(t+1) = z(t)exp(g*d*perm/(µh)
! where
! z = pond height (m)
! g = gravity (m/s^2)
! d = liquid density (kg/m^3)
! perm = minimum ice permeability
! µ = fluid viscosity = 1.79  (kg/m/s)
! h = ice thickness (m)

! viscosity of seawater (sal=35ppt) equation: visc = 0.0015*T^2 - 0.0608*T + 1.8284 m^2/s*10e6

          ! superimposed ice formation
          jj = 0
          superimposed = ice(sc,ic,mi)%pond%th - ice(sc,ic,mi)%pond%perc
          if (superimposed .gt. 0.d0 .and. z_snow .gt. 0) then

              !find layer that is less than 850000 g/m^3 density (850000 assumed by Nicolause et al. 2003)
              jj = 1
              do while (jj .le. ice(sc,ic,mi)%snow%z .and. superimposed .gt. 0.d0)
                  if (ice(sc,ic,mi)%snow%d(jj) .lt. 850000.d0) then

                      ! find new layer density from snow /water volume ratios
                      top_f = ice(sc,ic,mi)%snow%d(jj)/IceD*ice(sc,ic,mi)%snow%th(jj)  ! Volume Ice
                      bot_f = ice(sc,ic,mi)%pond%th*1.0893246187d0 ! Volume new ice
                      tmp1 = (top_f+bot_f)/ice(sc,ic,mi)%snow%th(jj) ! new ice fraction
                      tmp2 = 850000.d0/IceD ! max ice fraction in snow
                      if (tmp1 .gt. tmp2) then
                          superimposed = (tmp1 - tmp2)*superimposed*0.918 ! leftover melt
                          d_mean = 850000.d0 ! new snow densiy
                      else
                          superimposed = 0.d0
                          d_mean = tmp1*IceD ! new snow density
                      endif

                      ! find temp of snow snow layer with new melt water accreted
                      heat_mean = ice(sc,ic,mi)%snow%heat(jj)
                      tmp3 = ice(sc,ic,mi)%pond%th
                      tmp4 = ice(sc,ic,mi)%snow%t(jj)
                      heat_mean = ice(sc,ic,mi)%snow%heat(jj)*ice(sc,ic,mi)%snow%th(jj) - &
                        Lf*(ice(sc,ic,mi)%pond%th - th_pond_new)*1.d6 ! snow heat + latent heat of freezing of melt
                      heat_total = heat_mean/ice(sc,ic,mi)%snow%th(jj)/d_mean
                      t_mean = (0.2309-sqrt(0.05331481+0.0136* &
                        (heat_snow0-heat_total)))/(-0.0068)
                      t_mean = min(t_mean,kelvin0)

                      ! update snow layer physics

#include "sia2_env_snow_calc.inc.f90"

                      ice(sc,ic,mi)%snow%t(jj) = min(t_mean - kelvin0,d0_)
                      ice(sc,ic,mi)%snow%d(jj) = d_mean
                      ice(sc,ic,mi)%snow%heat(jj) = heat_mean
                      ice(sc,ic,mi)%snow%melt(jj) = 1.d0

                  endif

                  jj = jj+1

              enddo

              ! update pond
              th_pond_new = max(0.d0,ice(sc,ic,mi)%pond%perc + superimposed)


          endif

          ! remaining pond drains immediately
          if (th_pond_new .gt. 0.d0) then

#include "sia2_env_pond_percolate3.inc.f90"

              ! save for update if heat integration is valid
              th_pond_new = 0.d0

          endif

