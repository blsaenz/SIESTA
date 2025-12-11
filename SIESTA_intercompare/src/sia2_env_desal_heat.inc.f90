
		  jj = ii-z_snow
		  dsdt(jj) = d0_
		  if (ice(ic,mi)%id(int_z) .ge. 0.2 .and. (ice(ic,mi)%bv(jj) .ge. vb_crit) &
			  .and. (jj .gt. vb_open) .and. (ice(ic,mi)%bs(jj) .gt. f(mi)%s)) then
	
			  ! Cox-Weeks gravity drainage 
			  ! -----------------------------------------------
			  if (jj .eq. 1) then
				   T_grad = (ice(ic,mi)%t(1) - ice(ic,mi)%t(2)) / dth(2)
				  if (z_snow .gt. 0) then
					 tmp1 = (ice(ic,mi)%t(1) - ice(ic,mi)%snow%t(1))/dth(ice_1)
				  else
					 tmp1 = (ice(ic,mi)%t(1) - ice(ic,mi)%snow%ts)/dth(ice_1)
				  endif
			  else
				  if (jj .eq. int_z) then
				   T_grad = (ice(ic,mi)%t(jj-1) - f(mi)%t) / (dth(ii) + dth(ii+1))
				  else
				   T_grad = (ice(ic,mi)%t(jj-1) - ice(ic,mi)%t(jj+1)) / (dth(ii) + dth(ii+1))
				  endif
				  tmp1 = (ice(ic,mi)%t(jj) - ice(ic,mi)%t(jj-1))/dth(ii)
			  endif
		      tmp2 = (ice(ic,mi)%t(jj) - ice(ic,mi)%t(jj+1))/dth(ii+1)
	
			  ! dsdt rate (psu/s) 
			  dsdt(jj) = 4.2d-6*dtt_s*T_grad*(max(bv_mean*c_001 - 0.05,1.d-5))**(1.2) 
	
			  F0 = tmp1*ki(ii) + tmp2*ki(ii+1)                     
			 
			  dhdt = 9.873e-12*F0**3 - 6.621e-10*F0**2 + 4.474e-7*F0  ! (cm/s) dhdt vs. F0 regression at T=-1.8 freezing temp, petrich keff                     
	
			  if ((ice(ic,mi)%bv(jj) .ge. bv_conv) .and. (dhdt .gt. 0.) .and. &
				(ice(ic,mi)%s(jj) .gt. keff*ice(ic,mi)%bs(jj))) then
	
				  ! do agressive 'stable salinity' desal    
				  ! ----------------------------------------------------
				  keff = 0.19d0*(dhdt*7.4074d4)**(0.46)
				  dsdt3 = -1.*dtt_s*dhdt*c_01*ice(ic,mi)%bs(jj)*(1.-keff)/ice(ic,mi)%th(jj)
	
				  dsdt(jj) = max(dsdt(jj),dsdt3)
				  
			  endif
	
		  endif
