      ! START PLATELET LAYER CRAP
      ! ---------------------------------------------------------------------------------------
      ! ---------------------------------------------------------------------------------------

          !print *,'use_pl',use_pl

          if (use_pl .eq. 1) then
			  do jj=1,pl_max
				  am_sum = 0.
				  do ii=1,wavl
					  lambda = dble(ii)*10.+390.  ! in nm, need to make into m?
					  am(jj,ii) = aph(ii*10-9)*pl%smalg(jj)/c_chl*(pl%bv*c_001)  ! algae in bulk ice conc.
					  am_sum = am_sum+am(jj,ii)*10.   ! x10 b/c 10 nm wavelength bins
				  enddo
			      ad(jj,1) = pl%poc(jj)/pl%smalg(jj)*am_sum/ad_denom ! using fractional amount of poc compared to live algae
				  do ii=2,wavl
					  lambda = dble(ii)*10.+390.  ! in nm, need to make into m?
					  ad(jj,ii) = ad(jj,1)*exp(-0.008*((dble(ii)-1.)*10.))
				  enddo
			  enddo

			  ! Copy ED from last skeletal layer
			  do ii=1,wavl
				  Ed(1,ii) = Ed(sk_z+1,ii)
			  enddo

              bv_mean_real = real(pl%bv)
              ice(sc,ic,mi)%PAR_bot_pl = 0.

              do jj=1,pl_max
                  pl%pur(jj) = 0
				  do ii=1,wavl
					  ! find wavelength in nm
					  lambda = (dble(ii)*10.+390.)  ! in nm, need to make into m?
					  ! calc ice spectral attenuation
                      Kdice(jj,ii) = 3.5233 - (1.83e-2)*lambda(ii) + (3.28e-5)*lambda(ii)**2 - (1.81e-8)*lambda(ii)**3
					  ! calc particle spectral attenuation
					  Kdp(jj,ii) = (am(jj,ii)+ad(jj,ii))/0.656    ! mu=0.656
					  ! attenuate downwards
					  Ed(jj+1,ii) = Ed(jj,ii)*exp(-1.*(Kdp(jj,ii)+Kdice(jj,ii))*pl_th)
                      ! sum to find PUR
					  pl%PUR(jj) = pl%PUR(jj) &
						  + 0.5*(Ed(jj+1,ii)+Ed(jj,ii))* &
						  aph(ii*10-9)/aph_max*10. ! times 10, b/c 10nm wavelength bins
                      if (jj .eq. pl_max) then
                          ice(sc,ic,mi)%PAR_bot_pl = ice(sc,ic,mi)%PAR_bot_pl + Ed(jj+1,ii)*lda_f*10.
                      endif

				   enddo ! end spectral do
			  enddo ! end platelet layer do

              !print *, 'PAR_bots: ',ice(sc,ic,mi)%PAR_bot,ice(sc,ic,mi)%PAR_bot_pl

          else

              ice(sc,ic,mi)%PAR_bot_pl = ice(sc,ic,mi)%PAR_bot

          endif

      ! ---------------------------------------------------------------------------------------
      ! ---------------------------------------------------------------------------------------
      ! END PLATELET LAYER CRAP
