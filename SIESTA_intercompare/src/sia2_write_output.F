! file: write_output.F
! Sea Ice Algae Model 2 - Saenz & Arrigo
! Version beta
! ======================================================================
! ======================================================================

! ======================================================================
! Subroutine: sia2_write_output
! Purpose: writes out everything to NETCDF files as defined in the input constants
! ======================================================================

      SUBROUTINE sia2_write_output(pj,m,ice,pur,adv,f,mp_f,ease_f,sum_var, &
          out_fname,stations_mdh,stations_mdv,stn_output,step)

          use sia2_globals
          implicit none

          ! function arguments
      !----------------------------------------------------
          type (proj_type) :: pj
          type (meta_type) :: m(tcells)
          type (ice_type) :: ice(sda_n,ida_n,tcells)
          integer (kind=2) :: pur(z_max,dt_per_day,lda_n,sda_n,ida_n,tcells)
          type (adv_type) :: adv(tcells)
          type (forcing_type) :: f(tcells)
          type (mp_f_type) :: mp_f
          type (ease_f_type) :: ease_f
          double precision, dimension(mdh,mdv,105) :: sum_var
          character (LEN=80) :: out_fname
          integer, dimension(n_stations) :: stations_mdh, stations_mdv
          type (stn_write_type) :: stn_output(sda_n,ida_n,n_stations)
          integer :: step

      ! shared internal variables
      !----------------------------------------------------
          integer :: cur_date80,chunk_size
          double precision :: time

      ! private internal variables
      !----------------------------------------------------

          integer :: ii,jj,kk,ic,sc,mi,jjj,i_md,j_md,ref_i
          double precision :: bot_id,bot_f,af,af_tot,tmp1,tmp2, &
              bm_sum,bm_id,af_f,th_af_sum,s_sum,iv_sum,bv_sum, &
              top33,bot33,top_f



          chunk_size = 16

      ! find hours since 1980
      ! -----------------------------------------------------------
      call sia2_days_1980(d0_,cur_year,cur_date80)
      time = (cur_date80-1)*24 + cur_hour

      ! write out all data for specificed stations
      ! -----------------------------------------------------------
      do ii = 1,n_stations

        call sia2_write_station(pj,m,ice,pur,adv,f,mp_f,ease_f, &
            stations_mdh(ii),stations_mdv(ii),stn_output,time,step,ii)
      enddo

          if (stn_only .ne. 1) then

        if (start .or. restart .eq. 1) then
          ! initialize write counter to write out first step
                  if (start) then
              next_write_hour = cur_hour
                endif

          ! record lat/lon grid
          call sia2_write_cdf_double(mdh,'x',1,mdv,'y',1,0,' ',1,0,' ',1,2, &
          d0_,pj%lat(mdh1:mdh2,mdv1:mdv2),d0_,d0_,out_fname,'lats',4, &
          'degrees',7,'latitude position of EASE grid cell',35)

          call sia2_write_cdf_double(mdh,'x',1,mdv,'y',1,0,' ',1,0,' ',1,2, &
          d0_,pj%lon(mdh1:mdh2,mdv1:mdv2),d0_,d0_,out_fname,'lons',4, &
          'degrees',7,'longitude position of EASE grid cell',36)

        endif

        ! write out production and other data for entire model domain,
        ! but only as freqently as specifed in the write_f (write frequency)
        ! global variable
        ! -----------------------------------------------------------

        if (cur_hour .ge. next_write_hour) then

          ! write time to file
          ! -----------------------------------------------------------
          call sia2_write_cdf_double(write_step,'time',4,0,' ',1,0,' ',1,0,' ',1,1, &
          time,d0_,d0_,d0_,out_fname,'time',4,'hours',5,'Time',4)

          ! move desired vars to a grid for writing out to NETCDF output file
          ! -----------------------------------------------------------


!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ii,jj,kk,ic,sc,mi,jjj,i_md,j_md, &
!$OMP bot_id,bot_f,af,af_tot,tmp1,tmp2, &
!$OMP bm_sum,bm_id,af_f,th_af_sum,s_sum,iv_sum,bv_sum, &
!$OMP ref_i,top33,bot33,top_f)

!$OMP DO SCHEDULE (STATIC,chunk_size)
          do jj = mdv1,mdv2
            j_md = jj - mdv1 + 1
            do ii = mdh1,mdh2
              i_md = ii - mdh1 + 1
              mi = pj%mi(ii,jj)  ! find modeled varibles index from grid
              ! clean out temporary grids
              sum_var(i_md,j_md,:) = 0.d0

              if (mi .gt. 0) then
                if (m(mi)%status .ge. 1) then
                    af_tot = 0.d0
                    af_f = 0.d0
                    bm_sum = 0.d0
                    bm_id = 0.d0
                    th_af_sum = 0.d0
                    s_sum=0.d0
                    iv_sum=0.d0
                    bv_sum=0.d0

                  do ic=1,ida_n
                  do sc=1,sda_n
                                      af = ice(sc,ic,mi)%af
                                      if (af .gt. 0.d0) then
                                          af_tot = af_tot + af
                                          if(sc .gt. 1) then
                                              af_f = af_f + af
                                          endif
                                          sum_var(i_md,j_md,25) = sum_var(i_md,j_md,25) + ice(sc,ic,mi)%age*af
                                          sum_var(i_md,j_md,26) = sum_var(i_md,j_md,26) + ice(sc,ic,mi)%ridged*af
                      ! determine ice depth where "bottom ice" starts (? cm from congelation ice bottom)
                      bot_id = ice(sc,ic,mi)%id(ice(sc,ic,mi)%z) - bot_th
                      do kk = 1,ice(sc,ic,mi)%z

                        ! find mean layer depth
                        if (kk .eq. 1) then
                            tmp1 = ice(sc,ic,mi)%id(kk)*c_5
                        else
                            tmp1 = (ice(sc,ic,mi)%id(kk) - ice(sc,ic,mi)%id(kk-1))*c_5
                            tmp1 = tmp1 + ice(sc,ic,mi)%id(kk)
                        endif

                        ! recored biomass/area-weighted ice depth
                        tmp2 = ice(sc,ic,mi)%smalg(kk)*ice(sc,ic,mi)%bv(kk)*ice(sc,ic,mi)%th(kk)
                        bm_id = bm_id + tmp2*tmp1
                        bm_sum = bm_sum + tmp2

                        ! record mean bulk salinity
                        tmp2 = ice(sc,ic,mi)%th(kk)*af
                        th_af_sum = th_af_sum + tmp2
                        s_sum = s_sum + ice(sc,ic,mi)%s(kk)*tmp2

                        ! record mean bulk salinity
                        tmp2 = ice(sc,ic,mi)%th(kk)*af
                        th_af_sum = th_af_sum + tmp2
                        s_sum = s_sum + ice(sc,ic,mi)%s(kk)*tmp2

                        ! record total ice and brine volume
                        bv_sum = bv_sum + & ! brine volume
                          ice(sc,ic,mi)%bv(kk)*ice(sc,ic,mi)%th(kk)*c_001*cell_area*af
                        iv_sum = iv_sum + & ! ice volume
                          (1.d0-c_001*ice(sc,ic,mi)%bv(kk))*(1.d0-bb_f)* &
                          ice(sc,ic,mi)%th(kk)*cell_area*af

                        if (ice(sc,ic,mi)%id(kk) .gt. bot_id) then
                          if ((ice(sc,ic,mi)%id(kk)-ice(sc,ic,mi)%th(kk)) .ge. bot_id) then
                            ! layer is completely in "bottom ice"
                            sum_var(i_md,j_md,1) = sum_var(i_md,j_md,1) + ice(sc,ic,mi)%smalg(kk)* &
                            ice(sc,ic,mi)%bv(kk)*ice(sc,ic,mi)%th(kk)*cell_area*af  ! gC/pixel

                          else
                            ! layer spans bottom ice boundary - divide between internal
                            ! and bottom biomass, production vars
                            bot_f = (ice(sc,ic,mi)%id(kk)-bot_id)/ice(sc,ic,mi)%th(kk) ! bottom ice fraction of layer
                            ! converting to g/pixel - (brine_v/1000)*(cell_area*1.e6)*1g/1000mg = brine_v*cell_area
                            sum_var(i_md,j_md,1) = sum_var(i_md,j_md,1) + ice(sc,ic,mi)%smalg(kk)* &
                            ice(sc,ic,mi)%bv(kk)*ice(sc,ic,mi)%th(kk)*cell_area*bot_f*af  ! gC/pixel
                            sum_var(i_md,j_md,2) = sum_var(i_md,j_md,2) + ice(sc,ic,mi)%smalg(kk)* &
                            ice(sc,ic,mi)%bv(kk)*ice(sc,ic,mi)%th(kk)*cell_area*(1.-bot_f)*af  ! gC/pixel
                          endif
                        else ! layer is completely in "internal ice"
                          sum_var(i_md,j_md,2) = sum_var(i_md,j_md,2) + ice(sc,ic,mi)%smalg(kk)* &
                          ice(sc,ic,mi)%bv(kk)*ice(sc,ic,mi)%th(kk)*cell_area*af  ! gC/pixel
                        endif

                        ! new biomass stuff
                        ref_i = 54 + ((sc-1)*ic_n+(ic-1))*5+1
                        top33 = c_333*ice(sc,ic,mi)%id(ice(sc,ic,mi)%z)
                        !bot33 = ice(sc,ic,mi)%id(ice(sc,ic,mi)%z)-0.2
                        bot33 = ice(sc,ic,mi)%id(ice(sc,ic,mi)%z)-top33


                        if (ice(sc,ic,mi)%id(kk) .gt. bot33) then
                          if ((ice(sc,ic,mi)%id(kk)-ice(sc,ic,mi)%th(kk)) .ge. bot33) then
                            ! layer is completely in "bottom ice"
!                           sum_var(i_md,j_md,ref_i+2) = sum_var(i_md,j_md,ref_i+2) + ice(sc,ic,mi)%smalg(kk)* &
!                           ice(sc,ic,mi)%bv(kk)*ice(sc,ic,mi)%th(kk)*1.d-3  ! mgC/m^2
                            sum_var(i_md,j_md,ref_i+2) = sum_var(i_md,j_md,ref_i+2) + ice(sc,ic,mi)%s(kk)* &
                            ice(sc,ic,mi)%th(kk) ! psu

                          else
                            ! layer spans bottom ice boundary - divide between internal
                            ! and bottom biomass, production vars
                            bot_f = (ice(sc,ic,mi)%id(kk)-bot33)/ice(sc,ic,mi)%th(kk) ! bottom ice fraction of layer
                            ! converting to g/pixel - (brine_v/1000) = brine_v*m*1e-3 = mg C m-2
!                           sum_var(i_md,j_md,ref_i+2) = sum_var(i_md,j_md,ref_i+2) + ice(sc,ic,mi)%smalg(kk)* &
!                           ice(sc,ic,mi)%bv(kk)*ice(sc,ic,mi)%th(kk)*bot_f *1.d-3  ! mgC/m^2
!                           sum_var(i_md,j_md,ref_i+1) = sum_var(i_md,j_md,ref_i+1) + ice(sc,ic,mi)%smalg(kk)* &
!                           ice(sc,ic,mi)%bv(kk)*ice(sc,ic,mi)%th(kk)*(1.-bot_f)*1.d-3  ! mgC/m^2
                            sum_var(i_md,j_md,ref_i+2) = sum_var(i_md,j_md,ref_i+2) + ice(sc,ic,mi)%s(kk)* &
                            ice(sc,ic,mi)%th(kk)*bot_f  ! psu
                            sum_var(i_md,j_md,ref_i+1) = sum_var(i_md,j_md,ref_i+1) + ice(sc,ic,mi)%s(kk)* &
                            ice(sc,ic,mi)%th(kk)*(1.-bot_f)  !psu
                          endif
                        elseif (ice(sc,ic,mi)%id(kk) .lt. top33) then

                          if ((ice(sc,ic,mi)%id(kk) - ice(sc,ic,mi)%th(kk)) .lt. top33) then
                            ! layer totally in top 33%
!                           sum_var(i_md,j_md,ref_i) = sum_var(i_md,j_md,ref_i) + ice(sc,ic,mi)%smalg(kk)* &
!                           ice(sc,ic,mi)%bv(kk)*ice(sc,ic,mi)%th(kk)*1.d-3  ! mgC/m^2
                            sum_var(i_md,j_md,ref_i) = sum_var(i_md,j_md,ref_i) + ice(sc,ic,mi)%s(kk)* &
                            ice(sc,ic,mi)%th(kk)  ! psu
                          else
                            ! layer spans top ice boundary - divide between internal
                            ! and bottom biomass, production vars
                            top_f = (ice(sc,ic,mi)%id(kk)-top33)/ice(sc,ic,mi)%th(kk) ! top ice fraction of layer
                            ! converting to g/pixel - (brine_v/1000) = brine_v*m*1e-3 = mg C m-2
!                            sum_var(i_md,j_md,ref_i+1) = sum_var(i_md,j_md,ref_i+1) + ice(sc,ic,mi)%smalg(kk)* &
!                            ice(sc,ic,mi)%bv(kk)*ice(sc,ic,mi)%th(kk)*top_f *1.d-3  ! mgC/m^2
!                            sum_var(i_md,j_md,ref_i) = sum_var(i_md,j_md,ref_i) + ice(sc,ic,mi)%smalg(kk)* &
!                            ice(sc,ic,mi)%bv(kk)*ice(sc,ic,mi)%th(kk)*(1.-top_f)*1.d-3  ! mgC/m^2
                            sum_var(i_md,j_md,ref_i+1) = sum_var(i_md,j_md,ref_i+1) + ice(sc,ic,mi)%s(kk)* &
                            ice(sc,ic,mi)%th(kk)*top_f  ! psu
                            sum_var(i_md,j_md,ref_i) = sum_var(i_md,j_md,ref_i) + ice(sc,ic,mi)%s(kk)* &
                            ice(sc,ic,mi)%th(kk)*(1.-top_f)  ! psu

                          endif

                        else

                          ! layer is in middle 33%
!                          sum_var(i_md,j_md,ref_i+1) = sum_var(i_md,j_md,ref_i+1) + ice(sc,ic,mi)%smalg(kk)* &
!                          ice(sc,ic,mi)%bv(kk)*ice(sc,ic,mi)%th(kk)*1.d-3  ! mgC/m^2
                          sum_var(i_md,j_md,ref_i+1) = sum_var(i_md,j_md,ref_i+1) + ice(sc,ic,mi)%s(kk)* &
                          ice(sc,ic,mi)%th(kk)  ! psu

                        endif

                        ! record af, ice depth of category
                          ! layer is in middle 33%
                          sum_var(i_md,j_md,ref_i+4) = ice(sc,ic,mi)%id(ice(sc,ic,mi)%z)  ! m
                          ! layer is in middle 33%
                          sum_var(i_md,j_md,ref_i+3) = ice(sc,ic,mi)%af  ! fraction



                      enddo

                      ! grid non-layered vars for netcdf output
                      sum_var(i_md,j_md,7) = sum_var(i_md,j_md,7) + ice(sc,ic,mi)%id(ice(sc,ic,mi)%z)*af  ! total ice depth
                      sum_var(i_md,j_md,8) = sum_var(i_md,j_md,8) + ice(sc,ic,mi)%t(1)*af  ! top layer ice temp
                      sum_var(i_md,j_md,19) = sum_var(i_md,j_md,19) + ice(sc,ic,mi)%snow%ts*af  ! top layer ice temp
                      sum_var(i_md,j_md,10) = sum_var(i_md,j_md,10) + ice(sc,ic,mi)%sh_prev*af  ! time-step snow height

                                          sum_var(i_md,j_md,34) = sum_var(i_md,j_md,34) + af*bm_id/bm_sum


                                      endif ! end of valid ice category test

                  enddo ! end of sda sub-loop
                  enddo ! end of ida sub-loop

                                  ! divide by icecon to get mean values across categories
                                  sum_var(i_md,j_md,7) = sum_var(i_md,j_md,7)/af_tot
                                  sum_var(i_md,j_md,8) = sum_var(i_md,j_md,8)/af_tot
                                  sum_var(i_md,j_md,19) = sum_var(i_md,j_md,19)/af_tot
                                  sum_var(i_md,j_md,10) = sum_var(i_md,j_md,10)/af_tot

                                  sum_var(i_md,j_md,25) = sum_var(i_md,j_md,25)/af_tot  ! age
                                  sum_var(i_md,j_md,26) = sum_var(i_md,j_md,26)/af_tot  ! ridged

                                  sum_var(i_md,j_md,34) = sum_var(i_md,j_md,34)/af_tot  ! mean biomass location


                                  !  record pixel-wide stats
                                  sum_var(i_md,j_md,9) = f(mi)%ic_interp  ! ice concentration

                  sum_var(i_md,j_md,11) = m(mi)%pr_clim  !
                  sum_var(i_md,j_md,12) = m(mi)%melt_loss  !
                  sum_var(i_md,j_md,13) = m(mi)%adv_loss  !
                  sum_var(i_md,j_md,14) = m(mi)%md_loss  !

                                  sum_var(i_md,j_md,15) = m(mi)%a_convg  !
                  sum_var(i_md,j_md,16) = m(mi)%a_new  !
                  sum_var(i_md,j_md,17) = m(mi)%cong_growth  !
                  sum_var(i_md,j_md,18) = m(mi)%snow_growth  !

                  sum_var(i_md,j_md,20) = m(mi)%pr_ssmi

                  sum_var(i_md,j_md,21) = f(mi)%ivu_interp ! ice vec u
                  sum_var(i_md,j_md,22) = f(mi)%ivv_interp ! ice vec v
                  sum_var(i_md,j_md,23) = af_tot  ! ice vec correction
                  sum_var(i_md,j_md,24) = m(mi)%a_drop  ! ice dropped as a result of advection scheme

                  sum_var(i_md,j_md,23) = af_tot  !

                  sum_var(i_md,j_md,27) = m(mi)%tlim  ! mean biomass-weighted temp
                  sum_var(i_md,j_md,28) = m(mi)%nlim  ! biomass-weighed N limitation
                  sum_var(i_md,j_md,29) = m(mi)%plim  ! biomass-weighed P limitation
                  sum_var(i_md,j_md,30) = m(mi)%silim  ! biomass-weighed Si limitation
                  sum_var(i_md,j_md,31) = m(mi)%slim  ! biomass-weighed salinity limitation
                  sum_var(i_md,j_md,32) = m(mi)%llim  ! biomass-weighed light limitation
                  sum_var(i_md,j_md,33) = af_f/af_tot  ! fraction flooded

                  sum_var(i_md,j_md,49) = s_sum/th_af_sum   ! bulk ice salinity

                  sum_var(i_md,j_md,50) = bv_sum ! brine volume
                  sum_var(i_md,j_md,51) = iv_sum ! ice volume

                  sum_var(i_md,j_md,52) = m(mi)%void_growth ! ice volume


                  sum_var(i_md,j_md,105) = f(mi)%fw ! calculated ocean heat flux


                endif

                              ! record sums that may come from de-activated grid cells
                              sum_var(i_md,j_md,3) = m(mi)%prod_sum_int  ! internal ice production (g/C/pixel)
                              sum_var(i_md,j_md,4) = m(mi)%prod_sum_bot  ! bottom ice production (g/C/pixel)
                              sum_var(i_md,j_md,5) = m(mi)%bm_lost  ! biomass lost in water column (g/C/pixel)

                sum_var(i_md,j_md,35) =  m(mi)%salt_flux  ! dsdt
                sum_var(i_md,j_md,36) =  m(mi)%h2o_flux  ! equiv. freshwater input

                !if (m(mi)%p_wgt_int .gt. 0.d0) then
                  sum_var(i_md,j_md,37) =  m(mi)%pwid_sum_int!/m(mi)%p_wgt_int  ! interior prod-weighted ice depth
                  sum_var(i_md,j_md,39) =  m(mi)%pwsd_sum_int!/m(mi)%p_wgt_int  ! interior prod-weighted snow depth
                  sum_var(i_md,j_md,41) =  m(mi)%pwt_sum_int!/m(mi)%p_wgt_int  ! interior prod-weighted ice temp
                !endif
                !if (m(mi)%p_wgt_af_int .gt. 0.d0) then
                  sum_var(i_md,j_md,38) =  m(mi)%pwid_sum_af_int!/m(mi)%p_wgt_af_int  ! interior areal prod-weighted ice depth
                  sum_var(i_md,j_md,40) =  m(mi)%pwsd_sum_af_int!/m(mi)%p_wgt_af_int  ! interior areal prod-weighted snow depth
                  sum_var(i_md,j_md,42) =  m(mi)%pwt_sum_af_int!/m(mi)%p_wgt_af_int  ! interior areal prod-weighted ice temp
                !endif

                !if (m(mi)%p_wgt_bot .gt. 0.d0) then
                    sum_var(i_md,j_md,43) =  m(mi)%pwid_sum_bot!/m(mi)%p_wgt_bot  ! bottom prod-weighted ice depth
                    sum_var(i_md,j_md,45) =  m(mi)%pwsd_sum_bot!/m(mi)%p_wgt_bot  ! bottom prod-weighted snow depth
                    sum_var(i_md,j_md,47) =  m(mi)%pwt_sum_bot!/m(mi)%p_wgt_bot  ! bottom prod-weighted ice temp
                !endif
                !if (m(mi)%p_wgt_af_bot .gt. 0.d0) then
                    sum_var(i_md,j_md,44) =  m(mi)%pwid_sum_af_bot!/m(mi)%p_wgt_af_bot  ! bottom areal prod-weighted ice depth
                    sum_var(i_md,j_md,46) =  m(mi)%pwsd_sum_af_bot!/m(mi)%p_wgt_af_bot  ! bottom areal prod-weighted snow depth
                    sum_var(i_md,j_md,48) =  m(mi)%pwt_sum_af_bot!/m(mi)%p_wgt_af_bot  ! bottom areal prod-weighted ice temp
                !endif
                sum_var(i_md,j_md,53) = m(mi)%p_wgt_af_bot
                sum_var(i_md,j_md,54) = m(mi)%p_wgt_af_int

              endif
            enddo
          enddo


!$OMP END DO
!$OMP DO SCHEDULE (STATIC,chunk_size)

          ! clear cummulative vars
          do ii = 1,tcells
              m(ii)%prod_sum_int = 0.d0
              m(ii)%prod_sum_bot = 0.d0
              m(ii)%bm_lost = 0.d0
              m(ii)%melt_loss = 0.d0
              m(ii)%adv_loss = 0.d0
              m(ii)%md_loss = 0.d0
              m(ii)%cong_growth = 0.d0
              m(ii)%snow_growth = 0.d0
              m(ii)%void_growth = 0.d0
              m(ii)%adv_gain = 0.d0
              m(ii)%md_gain = 0.d0
              m(ii)%a_convg = 0.d0
              m(ii)%a_new = 0.d0
              m(ii)%a_drop = 0.d0
              m(ii)%pr_clim = 0.d0
              m(ii)%pr_ssmi = 0.d0
              m(ii)%tlim = 0.d0
              m(ii)%nlim = 0.d0
              m(ii)%silim = 0.d0
              m(ii)%plim = 0.d0
              m(ii)%llim = 0.d0
              m(ii)%slim = 0.d0
              m(ii)%slim = 0.d0

              m(ii)%salt_flux = 0.d0
              m(ii)%h2o_flux = 0.d0

              m(ii)%p_wgt_int = 0.d0
              m(ii)%p_wgt_af_int = 0.d0
              m(ii)%p_wgt_bot = 0.d0
              m(ii)%p_wgt_af_bot = 0.d0

              m(ii)%pwid_sum_int = 0.d0
              m(ii)%pwid_sum_af_int = 0.d0
              m(ii)%pwsd_sum_int = 0.d0
              m(ii)%pwsd_sum_af_int = 0.d0
              m(ii)%pwt_sum_int = 0.d0
              m(ii)%pwt_sum_af_int = 0.d0
              m(ii)%pwid_sum_bot = 0.d0
              m(ii)%pwid_sum_af_bot = 0.d0
              m(ii)%pwsd_sum_bot = 0.d0
              m(ii)%pwsd_sum_af_bot = 0.d0
              m(ii)%pwt_sum_bot = 0.d0
              m(ii)%pwt_sum_af_bot = 0.d0

          enddo

!$OMP END DO
!$OMP END PARALLEL


          ! increment counters
          ! -----------------------------------------------------------
                  ! tell threaded function to write to disk
                  do_write = .true.
                  ! index (number) of records written
                  write_step = write_step+1
          ! used to test when to write next
          next_write_hour = cur_hour + write_f

              endif ! end of write hour test
          endif ! end of stn_only test


      end SUBROUTINE sia2_write_output



! ======================================================================
! Subroutine: sia2_write_station
! Purpose: writes out all forcing and calculated variables for a
! particular station, as specified in the stations.txt file
! ======================================================================

      SUBROUTINE sia2_write_station(pj,m,ice,pur,adv,f,mp_f,ease_f, &
            sh,sv,stn_output,time,step,stn)

          use sia2_globals
          implicit none

          ! function arguments
      !----------------------------------------------------
          type (proj_type) :: pj
          type (meta_type) :: m(tcells)
          type (ice_type) :: ice(sda_n,ida_n,tcells)
          integer (kind=2) :: pur(z_max,dt_per_day,lda_n,sda_n,ida_n,tcells)
          type (adv_type) :: adv(ida_n,tcells)
          type (forcing_type) :: f(tcells)
          type (mp_f_type) :: mp_f
          type (ease_f_type) :: ease_f
          type (stn_write_type) :: stn_output(sda_n,ida_n,n_stations)
          integer :: sh,sv,step,stn
          double precision :: time,testvar

          ! internal variables
      !----------------------------------------------------
          integer :: mi,ic,lc,sc,ii,jj,write_h,write_v
          double precision :: id_temp,af_temp


          if (((sh .ge. mdh1) .and. (sh .le. mdh2) .and. (sv .ge. mdv1) .and. (sv .le. mdv2)) &
              .or. (boundary_file .eq. 1) .or. (stn_only .eq. 1) .or. (monte_carlo .eq. 1)) then

              do ic=1,ida_n
              do sc=1,sda_n

              ! turn on station for writing
              stn_output(sc,ic,stn)%valid = .true.

        ! find model domain reference to ease grid station coordinates
        mi = pj%mi(sh,sv)
              if (monte_carlo .eq. 1) then
                  mi = stn
              endif

        if (start) then

          call sia2_write_construct_fname(stn_output(sc,ic,stn)%fname,'station_--_------     ', &
              18,'nc ',2,datestr)

          ! find station output name
                  if(monte_carlo .eq. 1) then
                      write_h = 999
                      write_v = stn
          elseif (boundary_file .eq. 1) then
                      write_h = 999
                      write_v = 999
          else
                      write_h = sh
                      write_v = sv
          endif
          write(stn_output(sc,ic,stn)%fname(40:40),'(i1)') sc
          write(stn_output(sc,ic,stn)%fname(41:41),'(i1)') ic
          write(stn_output(sc,ic,stn)%fname(43:45),'(3i3)') write_h   ! output/station_XX_XXX-XXX...
          write(stn_output(sc,ic,stn)%fname(47:49),'(3i3)') write_v
          print *, 'Station Output File: ',trim(stn_output(sc,ic,stn)%fname)

                  if (wr_stations .ne. 0) then

                  ! write lat/lon to netcdf
          call sia2_write_cdf_double(1,'x',1,1,'y',1,step,'time',4,0,' ',1,3, &
          d0_,d0_,m(mi)%lat,d0_,stn_output(sc,ic,stn)%fname,'lats',4, &
          'degrees',7,'latitude position of EASE grid cell',35)

          call sia2_write_cdf_double(1,'x',1,1,'y',1,step,'time',4,0,' ',1,3, &
          d0_,d0_,m(mi)%lon,d0_,stn_output(sc,ic,stn)%fname,'lons',4, &
          'degrees',7,'longitude position of EASE grid cell',36)

          endif

        endif

              ! writeout timing
              stn_output(sc,ic,stn)%step = step ! index to write step
        stn_output(sc,ic,stn)%time = time ! actual time in hours since 1980

        ! record single (1d) vars
        stn_output(sc,ic,stn)%s(1) = ease_f%icecon(sh,sv)
        stn_output(sc,ic,stn)%s(2) = f(mi)%ic_interp
        stn_output(sc,ic,stn)%s(3) = f(mi)%sh_interp
        stn_output(sc,ic,stn)%s(4) = ice(sc,ic,mi)%sh_prev
        stn_output(sc,ic,stn)%s(7) = f(mi)%at
        stn_output(sc,ic,stn)%s(8) = f(mi)%pr
        stn_output(sc,ic,stn)%s(14) = m(mi)%lat
        stn_output(sc,ic,stn)%s(15) = m(mi)%lon
        stn_output(sc,ic,stn)%s(11) = mp_f%ncep_interp%at(m(mi)%x_mp,m(mi)%y_mp)
        stn_output(sc,ic,stn)%s(12) = dble(m(mi)%status)


              if (ice(sc,ic,mi)%z .gt. 0) then
                  stn_output(sc,ic,stn)%s(5) = ice(sc,ic,mi)%id(ice(sc,ic,mi)%z)
              else
                  stn_output(sc,ic,stn)%s(5) = 0.d0
        endif
        stn_output(sc,ic,stn)%s(6) = dble(ice(sc,ic,mi)%snow%z)
        stn_output(sc,ic,stn)%s(9) = ice(sc,ic,mi)%snow%ts
        stn_output(sc,ic,stn)%s(10) = ice(sc,ic,mi)%Ed0
        stn_output(sc,ic,stn)%s(13) = dble(ice(sc,ic,mi)%z)
        stn_output(sc,ic,stn)%s(16) = ice(sc,ic,mi)%PAR_bot
        stn_output(sc,ic,stn)%s(17) = ice(sc,ic,mi)%PAR_bot_pl
        stn_output(sc,ic,stn)%s(18) = ice(sc,ic,mi)%af

        ! record z-dimensional vars
              stn_output(sc,ic,stn)%z(1,:) = ice(sc,ic,mi)%snow%th
              stn_output(sc,ic,stn)%z(2,:) = ice(sc,ic,mi)%snow%t
              stn_output(sc,ic,stn)%z(3,:) = ice(sc,ic,mi)%snow%d
              stn_output(sc,ic,stn)%z(4,:) = ice(sc,ic,mi)%snow%heat
              stn_output(sc,ic,stn)%z(5,:) = ice(sc,ic,mi)%smalg
              stn_output(sc,ic,stn)%z(6,:) = ice(sc,ic,mi)%prod
              stn_output(sc,ic,stn)%z(7,:) = ice(sc,ic,mi)%heat
              stn_output(sc,ic,stn)%z(8,:) = ice(sc,ic,mi)%th
              stn_output(sc,ic,stn)%z(9,:) = ice(sc,ic,mi)%id
              stn_output(sc,ic,stn)%z(10,:) = ice(sc,ic,mi)%t
              stn_output(sc,ic,stn)%z(12,:) = ice(sc,ic,mi)%Ik1
              stn_output(sc,ic,stn)%z(13,:) = ice(sc,ic,mi)%bv
              stn_output(sc,ic,stn)%z(14,:) = ice(sc,ic,mi)%bd
              stn_output(sc,ic,stn)%z(15,:) = ice(sc,ic,mi)%bs
              stn_output(sc,ic,stn)%z(16,:) = ice(sc,ic,mi)%s
              stn_output(sc,ic,stn)%z(17,:) = ice(sc,ic,mi)%d
              stn_output(sc,ic,stn)%z(18,:) = ice(sc,ic,mi)%llim
              stn_output(sc,ic,stn)%z(19,:) = ice(sc,ic,mi)%nlim
              stn_output(sc,ic,stn)%z(20,:) = ice(sc,ic,mi)%plim
              stn_output(sc,ic,stn)%z(21,:) = ice(sc,ic,mi)%silim
              stn_output(sc,ic,stn)%z(22,:) = ice(sc,ic,mi)%slim
              stn_output(sc,ic,stn)%z(23,:) = ice(sc,ic,mi)%snow%melt
              stn_output(sc,ic,stn)%z(24,:) = ice(sc,ic,mi)%dsdt
              stn_output(sc,ic,stn)%z(25,:) = ice(sc,ic,mi)%fbv
              stn_output(sc,ic,stn)%z(26,:) = ice(sc,ic,mi)%dhdt_conv
              stn_output(sc,ic,stn)%z(27,:) = ice(sc,ic,mi)%f0
              stn_output(sc,ic,stn)%z(28,:) = ice(sc,ic,mi)%dsdt3
              stn_output(sc,ic,stn)%z(29,:) = ice(sc,ic,mi)%tgrad
              stn_output(sc,ic,stn)%z(30,:) = ice(sc,ic,mi)%gmax
              do ii=1,z_max
!                  stn_output(stn)%z(11,ii) = pur_c*(dble(pur(ii,pur_clock,lc,ic,mi))+127.d0)**2 ! for byte pur array
                  stn_output(sc,ic,stn)%z(11,ii) = 0.
                  do jj=1,lda_n
            stn_output(sc,ic,stn)%z(11,ii) = stn_output(sc,ic,stn)%z(11,ii) + &
              pur_c_2*dble(pur(ii,pur_clock,jj,sc,ic,mi) + 32767)*lda_d     ! for int(2) size pur array
                  enddo
!                  print *,ii,'--',stn_output(stn)%z(11,ii)
              enddo

        ! record z-dimensional vars
              stn_output(sc,ic,stn)%z1(1,:) = ice(sc,ic,mi)%no3
              stn_output(sc,ic,stn)%z1(2,:) = ice(sc,ic,mi)%nh4
              stn_output(sc,ic,stn)%z1(3,:) = ice(sc,ic,mi)%po4
              stn_output(sc,ic,stn)%z1(4,:) = ice(sc,ic,mi)%sioh4
              stn_output(sc,ic,stn)%z1(5,:) = ice(sc,ic,mi)%poc

        ! record wavl-dimension vars
              stn_output(sc,ic,stn)%wavl(1,:) = mp_f%Edif(f_index)%Ed(m(mi)%x_mp,m(mi)%y_mp,:)
              stn_output(sc,ic,stn)%wavl(2,:) = mp_f%Edir(f_index)%Ed(m(mi)%x_mp,m(mi)%y_mp,:)
              !stn_output(stn)%wavl(3,:) = ice(sc,ic,mi)%PAR_bot

          ! record icevecs
          stn_output(sc,ic,stn)%icevec(1) = f(mi)%ivu_interp
          stn_output(sc,ic,stn)%icevec(2) = f(mi)%ivv_interp

              enddo
              enddo

          elseif (start) then
              print *, 'Station location is outside of model domain - no station data written'
              stn_output(:,:,stn)%valid = .false.
          else
              stn_output(:,:,stn)%valid = .false.
          endif

      end SUBROUTINE sia2_write_station




! ======================================================================
! Subroutine: sia2_write_cdf_double
! Purpose: General routine to handle writing data out to a netcdf
! file.  If file does not exist, it creates a new ones and defines all possible
! dimensions for output files, so that anything can be added later to the same
! file.
! ======================================================================

       subroutine sia2_write_cdf_double(dimx,dimx_n,dimx_n_l,dimy,dimy_n, &
          dimy_n_l,dimz,dimz_n,dimz_n_l,dimt,dimt_n,dimt_n_l,ndims,value1, &
          values2,values3,values4,file_n,var_n,var_n_l,units,units_l,long_name, &
          long_name_l)

          use sia2_globals

#include "netcdf.inc"

          integer, dimension(1:4) :: dstart,dcount
          integer :: ncid,varid,err,var_n_l,units_l,long_name_l
          integer :: dimx,dimy,dimz,dimt,ndims,tmp
          integer :: dimx_n_l,dimy_n_l,dimz_n_l,dimt_n_l,i
!          character (LEN=20) :: dimx_n
!          character (LEN=20) :: dimy_n
!          character (LEN=20) :: dimz_n
!          character (LEN=20) :: dimt_n
!          character (LEN=80) :: file_n
!          character (LEN=20) :: var_n
!          character (LEN=60) :: units
!          character (LEN=60) :: long_name
          character (LEN=*) :: dimx_n
          character (LEN=*) :: dimy_n
          character (LEN=*) :: dimz_n
          character (LEN=*) :: dimt_n
          character (LEN=*) :: file_n
          character (LEN=*) :: var_n
          character (LEN=*) :: units
          character (LEN=*) :: long_name
          double precision :: value1
          double precision, dimension(dimx,dimy) :: values2
          double precision, dimension(dimx,dimy) :: values3
          double precision, dimension(dimx,dimy,dimz) :: values4

          ! attempt to open file - if not present, create it, with standard dimensions
          err = NF_OPEN(trim(file_n),NF_WRITE,ncid)
          if (err .ne. NF_NOERR) then
              ! create NETCDF file
              print *,'Creating NETCDF file: ',trim(file_n)
              call sia2_write_create_cdf(file_n,ncid,dimx,dimy)
              !print *,'Done creating file.'
          endif

          err = nf_set_fill(ncid,NF_NOFILL,tmp)

          err = nf_inq_varid(ncid,var_n(1:var_n_l),varid)
          if (err .ne. NF_NOERR) then
              ! create variable if it doesn't exist
              err = nf_redef(ncid)
              ! print *,'Creating NETCDF var: ',var_n(1:var_n_l)
              call sia2_write_create_double_var(ncid,varid,dimx,dimx_n,dimx_n_l, &
                  dimy,dimy_n,dimy_n_l,dimz,dimz_n,dimz_n_l,dimt,dimt_n,dimt_n_l, &
                  ndims,file_n,var_n,var_n_l,units,units_l, &
                  long_name,long_name_l)
              err = NF__ENDDEF(ncid,20000,4,0,4)
          endif

          ! define bounrdaries of data write
          dstart(4) = dimt
          dstart(3) = 1
          dstart(2) = 1
          dstart(1) = 1
          dcount(4) = 1
          dcount(3) = dimz
          dcount(2) = dimy
          dcount(1) = dimx
          if (ndims .eq. 3) then
              dstart(3) = dimz
              dcount(3) = 1
          elseif ((ndims .lt. 1) .or. (ndims .gt. 4)) then
              print *, 'Error: number of NETCDF dimensions must be between 1 & 4'
              call exit(999)
          endif

          ! write data to variable
          if (ndims .eq. 4) then
              err = nf_put_vara_double(ncid,varid,dstart(1:ndims),dcount(1:ndims), &
                  values4)
          elseif (ndims .eq. 3) then
              err = nf_put_vara_double(ncid,varid,dstart(1:ndims),dcount(1:ndims), &
                  values3)
          elseif (ndims .eq. 2) then
              err = nf_put_vara_double(ncid,varid,dstart(1:ndims),dcount(1:ndims), &
                  values2)
          elseif (ndims .eq. 1) then
              err = nf_put_var1_double(ncid,varid,dimx,value1)
          endif

          ! check to see if it wrote correctly
          if (err .ne. NF_NOERR) then
              print *,'Error writing to NETCDF file...'
              print *,'file: ',trim(file_n)
        print *,'var: ',var_n(1:var_n_l)
        print *,'nc error: ', err
              err = -1
              value1=sqrt(dble(err))
              call exit(999)
          endif

          err = nf_close(ncid)

      end subroutine sia2_write_cdf_double


! ======================================================================
! Subroutine: sia2_write_open_cdf_double
! Purpose: Same thing as sia2_write_cdf_double, except that it assumes
! you have already opened the file (ncid) and doesn't close it.
! Created in attempt to do faster writes.
! ======================================================================

       subroutine sia2_write_open_cdf_double(dimx,dimx_n,dimx_n_l,dimy,dimy_n, &
          dimy_n_l,dimz,dimz_n,dimz_n_l,dimt,dimt_n,dimt_n_l,ndims,value1, &
          values2,values3,values4,file_n,var_n,var_n_l,units,units_l,long_name, &
          long_name_l,ncid)

          use sia2_globals

#include "netcdf.inc"

          integer, dimension(1:4) :: dstart,dcount
          integer :: ncid,varid,err,var_n_l,units_l,long_name_l
          integer :: dimx,dimy,dimz,dimt,ndims,tmp
          integer :: dimx_n_l,dimy_n_l,dimz_n_l,dimt_n_l,i
!          character (LEN=20) :: dimx_n
!          character (LEN=20) :: dimy_n
!          character (LEN=20) :: dimz_n
!          character (LEN=20) :: dimt_n
!          character (LEN=80) :: file_n
!          character (LEN=20) :: var_n
!          character (LEN=60) :: units
!          character (LEN=60) :: long_name
          character (LEN=*) :: dimx_n
          character (LEN=*) :: dimy_n
          character (LEN=*) :: dimz_n
          character (LEN=*) :: dimt_n
          character (LEN=*) :: file_n
          character (LEN=*) :: var_n
          character (LEN=*) :: units
          character (LEN=*) :: long_name
          double precision :: value1
          double precision, dimension(dimx,dimy) :: values2
          double precision, dimension(dimx,dimy) :: values3
          double precision, dimension(dimx,dimy,dimz) :: values4

          err = nf_set_fill(ncid,NF_NOFILL,tmp)

          err = nf_inq_varid(ncid,var_n(1:var_n_l),varid)
          if (err .ne. NF_NOERR) then
              ! create variable if it doesn't exist
              err = nf_redef(ncid)
              !print *,'Creating NETCDF var: ',var_n(1:var_n_l)
              call sia2_write_create_double_var(ncid,varid,dimx,dimx_n,dimx_n_l, &
                  dimy,dimy_n,dimy_n_l,dimz,dimz_n,dimz_n_l,dimt,dimt_n,dimt_n_l, &
                  ndims,file_n,var_n,var_n_l,units,units_l, &
                  long_name,long_name_l)
              err = NF__ENDDEF(ncid,20000,4,0,4)
          endif

          ! define bounrdaries of data write
          dstart(4) = dimt
          dstart(3) = 1
          dstart(2) = 1
          dstart(1) = 1
          dcount(4) = 1
          dcount(3) = dimz
          dcount(2) = dimy
          dcount(1) = dimx
          if (ndims .eq. 3) then
              dstart(3) = dimz
              dcount(3) = 1
          elseif ((ndims .lt. 1) .or. (ndims .gt. 4)) then
              print *, 'Error: number of NETCDF dimensions must be between 1 & 4'
              call exit(999)
          endif

          ! write data to variable
          if (ndims .eq. 4) then
              err = nf_put_vara_double(ncid,varid,dstart(1:ndims),dcount(1:ndims), &
                  values4)
          elseif (ndims .eq. 3) then
              err = nf_put_vara_double(ncid,varid,dstart(1:ndims),dcount(1:ndims), &
                  values3)
          elseif (ndims .eq. 2) then
              err = nf_put_vara_double(ncid,varid,dstart(1:ndims),dcount(1:ndims), &
                  values2)
          elseif (ndims .eq. 1) then
              err = nf_put_var1_double(ncid,varid,dimx,value1)
          endif

          ! check to see if it wrote correctly
          if (err .ne. NF_NOERR) then
              print *,'Error writing to NETCDF file...'
              print *,'file: ',trim(file_n)
        print *,'var: ',var_n(1:var_n_l)
        print *,'nc error: ', err
              err = -1
              value1=sqrt(dble(err))
              call exit(999)
          endif

      end subroutine sia2_write_open_cdf_double



       subroutine sia2_nc_open(fname,ncid)

          implicit none
          integer         nf_open
!                         (character*(*)       path,
!                          integer             mode,
!                          integer             ncid)
          external        nf_open
!          character (LEN=80) :: fname
          character (LEN=*) :: fname
          integer :: ncid,err
          double precision :: testvar

          ! passing NF_WRITE (=1) to nf_open
          err = nf_open(trim(fname),1,ncid)
          if (err .ne. 0) then
              print *,'NETCDF file not found: ',trim(fname)
              print *,' NF_ERROR: ',err
              testvar = -1.
              testvar = sqrt(testvar)
              call exit(0)
          endif

      end subroutine sia2_nc_open

      subroutine sia2_nc_close(fname,ncid)

          implicit none
          integer         nf_close
!                         (integer             ncid)
          external        nf_close
          !character (LEN=80) :: fname
          character (LEN=*) :: fname
          integer :: ncid,err

          err = nf_close(ncid)
          if (err .ne. 0) then
              print *,'Failed to close NETCDF file: ',trim(fname)
              print *,' NF_ERROR: ',err
              call exit(0)
          endif

      end subroutine sia2_nc_close







! ======================================================================
! Subroutine: sia2_write_create_double_var
! Purpose: Creates an new variable in an already-open
! netcdf file (assumes all dimensions have already been defined)
! ======================================================================

      subroutine sia2_write_create_double_var(ncid,varid,dimx,dimx_n,dimx_n_l,dimy, &
          dimy_n,dimy_n_l,dimz,dimz_n,dimz_n_l,dimt,dimt_n,dimt_n_l,ndims, &
          file_n,var_n,var_n_l,units,units_l,long_name,long_name_l)

          use sia2_globals
          implicit none

#include "netcdf.inc"

          integer, dimension(1:4) :: all_dims
          integer :: ncid,varid,err,var_n_l,units_l,long_name_l
          integer :: dimx_n_l,dimy_n_l,dimz_n_l,dimt_n_l
          integer :: dimxid,dimyid,dimzid,dimtid
          integer :: dimx,dimy,dimz,dimt,ndims
!          character (LEN=20) :: dimx_n
!          character (LEN=20) :: dimy_n
!          character (LEN=20) :: dimz_n
!          character (LEN=20) :: dimt_n
!          character (LEN=80) :: file_n
!          character (LEN=20) :: var_n
!          character (LEN=60) :: units
!          character (LEN=60) :: long_name
          character (LEN=*) :: dimx_n
          character (LEN=*) :: dimy_n
          character (LEN=*) :: dimz_n
          character (LEN=*) :: dimt_n
          character (LEN=*) :: file_n
          character (LEN=*) :: var_n
          character (LEN=*) :: units
          character (LEN=*) :: long_name

          ! initialize vars to make ifort compiler happy
          dimxid = 0
          dimyid = 0
          dimzid = 0
          dimtid = 0

          ! find dimensions
          err = nf_inq_dimid(ncid,dimx_n(1:dimx_n_l),dimxid)
          if (err .ne. NF_NOERR) then
              print *,'Warning: cannot find NETCDF dimension: ',dimx_n(1:dimx_n_l)
          endif
          if (ndims .gt. 1) then
              err = nf_inq_dimid(ncid,dimy_n(1:dimy_n_l),dimyid)
              if (err .ne. NF_NOERR) then
                  print *,'Warning: cannot find NETCDF dimension: ',dimy_n(1:dimy_n_l)
              endif
              if (ndims .gt. 2) then
                  err = nf_inq_dimid(ncid,dimz_n(1:dimz_n_l),dimzid)
                  if (err .ne. NF_NOERR) then
                      print *,'Warning: cannot find NETCDF dimension: ',dimz_n(1:dimz_n_l)
                  endif
                  if (ndims .gt. 3) then
                      err = nf_inq_dimid(ncid,dimt_n(1:dimt_n_l),dimtid)
                      if (err .ne. NF_NOERR) then
                          print *,'Warning: cannot find NETCDF dimension: ',dimt_n(1:dimt_n_l)
                      endif
                  endif
              endif
          endif

          ! make dimensions variable from dims
          all_dims(1) = dimxid
          all_dims(2) = dimyid
          all_dims(3) = dimzid
          all_dims(4) = dimtid
!          all_dims(4) = dimxid
!          all_dims(3) = dimyid
!          all_dims(2) = dimzid
!          all_dims(1) = dimtid

          ! create variable variable
          err = nf_def_var(ncid,var_n(1:var_n_l),nf_double,ndims,all_dims(1:ndims),varid)
          if (err .ne. NF_NOERR) then
              ! create NETCDF file
        print *,'Error creating NETCDF var...'
        print *,'file: ',trim(file_n)
        print *,'var: ',var_n(1:var_n_l)
        print *,'nc error: ', err
              call exit(999)
          endif

          !err = nf_def_var_deflate(ncid, varid, 0, 1, 4)

          err = nf_put_att_real(ncid,varid,'_FillValue',nf_double,1,-9.91E36)
          err = nf_put_att_real(ncid,varid,'missing_value',nf_double,1,-9.91E36)
          err = nf_put_att_text(ncid,varid,'units', units_l, units)
          err = nf_put_att_text(ncid,varid,'long_name',long_name_l,long_name)

      end subroutine sia2_write_create_double_var

! ======================================================================
! Subroutine: sia2_write_create_cdf
! Purpose: Creates an new netcdf file call file_n, and assigns
! all dimensions possible from the current model run so that any vars
! can be added
! ======================================================================

      subroutine sia2_write_create_cdf(file_n,ncid,xlength,ylength)
          use sia2_globals
          implicit none

#include "netcdf.inc"

          integer :: err,ndim1,ncid,xdim,ydim,wavdim,zdim,tdim,varid,i,z1dim
          integer :: xout_dim,yout_dim,zout_dim,wavout_dim,tout_dim,z1out_dim
          integer :: ivdim,dim4,dim5,dim7,dim_dt_per_day,gridv_dim,gridh_dim
          integer :: pldim,plout_dim
          integer :: ivout_dim,out_dim4,out_dim5,out_dim7,out_dtpd,out_gridv,out_gridh
          integer :: xlength,ylength,bndsdim
          !character (LEN=80) :: file_n
          character (LEN=*) :: file_n



          !err=NF_CREATE(trim(file_n),IOR(NF_NOCLOBBER,NF_NETCDF4),ncid)
          err=NF_CREATE(trim(file_n),NF_NOCLOBBER,ncid)
          if (err .ne. NF_NOERR) then
              print *,'Error creating NETCDF file...'
              print *,'file: ',trim(file_n)
              print *,'nc error: ', err
              call exit(999)
          endif

          ! declare dimensions w/ length

          err = NF_DEF_DIM (ncid, 'x', xlength, xdim)
          err = NF_DEF_DIM (ncid, 'y', ylength, ydim)
          err = NF_DEF_DIM (ncid, 'z', z_max, zdim)
          err = NF_DEF_DIM (ncid, 'z1', z_max+1, z1dim)
          err = NF_DEF_DIM (ncid, 'bnds', 2, bndsdim)
          err = NF_DEF_DIM (ncid, 'ivd', 3, ivdim)
          err = NF_DEF_DIM (ncid, 'dim4', 4, dim4)
          err = NF_DEF_DIM (ncid, 'dim5', 5, dim5)
          err = NF_DEF_DIM (ncid, 'dim7', 7, dim7)
          err = NF_DEF_DIM (ncid, 'dt_per_day', dt_per_day, dim_dt_per_day)
          err = NF_DEF_DIM (ncid, 'grid_v', grid_v, gridv_dim)
          err = NF_DEF_DIM (ncid, 'grid_h', grid_h, gridh_dim)
          err = NF_DEF_DIM (ncid, 'wavl', wavl, wavdim)
          err = NF_DEF_DIM (ncid, 'time', NF_UNLIMITED, tdim)

          ! Create dimension characterstics
          ndim1 = 1

          if (use_pl .eq. 1) then
              err = NF_DEF_DIM (ncid, 'pl', pl_max, pldim)
        err = nf_def_var(ncid,'pl', nf_int,ndim1,pldim,plout_dim)
        err = nf_put_att_text(ncid,plout_dim, 'positive', 4, 'down')
        err = nf_put_att_text(ncid,plout_dim, 'units', 5, 'level')
        err = nf_put_att_text(ncid,plout_dim, 'long_name', 14, 'Platelet Level')
          endif

          ! setup X dimension as latitudinal model domain (subset of EASE grid)
          err = nf_def_var(ncid,'x', nf_int,ndim1,xdim,xout_dim)
          err = nf_put_att_text(ncid,xout_dim, 'units', 4, 'cell')
          err = nf_put_att_text(ncid,xout_dim, 'long_name', 20, 'EASE E/W cell number')

          ! setup Y dimension as longitudinal model domain (subset of EASE grid)
          err = nf_def_var(ncid,'y', nf_int,ndim1,ydim,yout_dim)
          err = nf_put_att_text(ncid,yout_dim, 'units', 4, 'cell')
          err = nf_put_att_text(ncid,yout_dim, 'long_name', 20, 'EASE N/S cell number')

          ! setup Z dimension as vertical model layers domain
          err = nf_def_var(ncid,'z', nf_int,ndim1,zdim,zout_dim)
          err = nf_put_att_text(ncid,zout_dim, 'positive', 4, 'down')
          err = nf_put_att_text(ncid,zout_dim, 'units', 5, 'level')
          err = nf_put_att_text(ncid,zout_dim, 'long_name', 9, 'Ice Level')

          ! setup Z1 dimension as vertical model layers domain
          err = nf_def_var(ncid,'z1', nf_int,ndim1,z1dim,z1out_dim)
          err = nf_put_att_text(ncid,z1out_dim, 'positive', 4, 'down')
          err = nf_put_att_text(ncid,z1out_dim, 'units', 5, 'level')
          err = nf_put_att_text(ncid,z1out_dim, 'long_name', 11, 'Ice Level+1')

          ! setup icevec dimension as vertical model layers domain
          err = nf_def_var(ncid,'ivd', nf_int,ndim1,ivdim,ivout_dim)
          err = nf_put_att_text(ncid,ivout_dim, 'units', 4, 'none')
          err = nf_put_att_text(ncid,ivout_dim, 'long_name', 14, 'storage size 3')

          err = nf_def_var(ncid,'dim4', nf_int,ndim1,dim4,out_dim4)
          err = nf_put_att_text(ncid,out_dim4, 'units', 4, 'none')
          err = nf_put_att_text(ncid,out_dim4, 'long_name', 14, 'storage size 4')

          err = nf_def_var(ncid,'dim5', nf_int,ndim1,dim5,out_dim5)
          err = nf_put_att_text(ncid,out_dim5, 'units', 4, 'none')
          err = nf_put_att_text(ncid,out_dim5, 'long_name', 14, 'storage size 5')

          err = nf_def_var(ncid,'dim7', nf_int,ndim1,dim7,out_dim7)
          err = nf_put_att_text(ncid,out_dim7, 'units', 4, 'none')
          err = nf_put_att_text(ncid,out_dim7, 'long_name', 14, 'storage size 7')

          ! setup dt_per_day dimension
          err = nf_def_var(ncid,'dt_per_day', nf_int,ndim1,dim_dt_per_day,out_dtpd)
          err = nf_put_att_text(ncid,out_dtpd, 'units', 4, 'none')
          err = nf_put_att_text(ncid,out_dtpd, 'long_name', 19, 'model steps per day')

          ! setup EASE grid_v dimension
          err = nf_def_var(ncid,'grid_v', nf_int,ndim1,gridv_dim,out_gridv)
          err = nf_put_att_text(ncid,out_gridv, 'units', 4, 'none')
          err = nf_put_att_text(ncid,out_gridv, 'long_name', 18, 'EASE grid vertical')

          ! setup EASE grid_h dimension
          err = nf_def_var(ncid,'grid_h', nf_int,ndim1,gridh_dim,out_gridh)
          err = nf_put_att_text(ncid,out_gridh, 'units', 4, 'none')
          err = nf_put_att_text(ncid,out_gridh, 'long_name', 20, 'EASE grid horizontal')

          ! setup wavl dimension as vertical model layers domain
          err = nf_def_var(ncid,'wavl', nf_int,ndim1,wavdim,wavout_dim)
          err = nf_put_att_text(ncid,wavout_dim, 'units', 2, 'nm')
          err = nf_put_att_text(ncid,wavout_dim, 'long_name', 10, 'Wavelength')

          ! setup time dimension
          err = nf_def_var(ncid,'time', nf_double,ndim1,tdim,tout_dim)
          err = nf_put_att_text(ncid,tout_dim, 'units', 4, 'hours')
          err = nf_put_att_text(ncid,tout_dim, 'long_name', 4, 'time')
          err = nf_put_att_text(ncid,tout_dim,'time_origin', 11, '1-Jan-1980')

          ! end netcdf definition
          err = NF__ENDDEF(ncid,20000,4,0,4)

          ! add values to x
          err = nf_inq_varid(ncid,'x',varid)
          do i=1,xlength
              err = nf_put_var1_int(ncid,varid,i,i)
          enddo

          ! add values to y
          err = nf_inq_varid(ncid,'y',varid)
          do i=1,ylength
              err = nf_put_var1_int(ncid,varid,i,i)
          enddo

          ! add values to z
          err = nf_inq_varid(ncid,'z',varid)
          do i=1,z_max
              err = nf_put_var1_int(ncid,varid,i,i)
          enddo

          ! add values to z1
          err = nf_inq_varid(ncid,'z1',varid)
          do i=1,z_max+1
              err = nf_put_var1_int(ncid,varid,i,i)
          enddo

          ! add values to ivd
          err = nf_inq_varid(ncid,'ivd',varid)
          do i=1,3
              err = nf_put_var1_int(ncid,varid,i,i)
          enddo

          ! add values to wavl
          err = nf_inq_varid(ncid,'wavl',varid)
          do i=1,wavl
              err = nf_put_var1_int(ncid,varid,i,400+10*(i-1))
          enddo

      end subroutine sia2_write_create_cdf


! ======================================================================
! Subroutine: sia2_blank_line
! Purpose: sets an input 60 character line to all blanks
! ======================================================================

      SUBROUTINE sia2_blank_line(line)
          character (LEN=80) :: line
          line(1:40) = '                                        '
          line(41:80) = '                                        '
      end SUBROUTINE sia2_blank_line


! ======================================================================
! Subroutine: sia2_do_datestr
! Purpose: makes custom runtime string with date & time for file output
! purposes, in the form DDMMMYYYY_HHMM
! ======================================================================

      SUBROUTINE sia2_do_datestr(date_time)

           ! date_and_time ( DATE, TIME, ZONE, VALUES )
           ! DATE is a character*8 variable that contains the date. CCYYMMDD
           ! TIME is a character*10 variable that contains the time. hhmmss.sss
           ! ZONE is a character*5 variable that contains the time zone. +hhmm
           ! VALUES is an integer array of dimensin 8, containing the year,
           !     the month, the day, the time difference with respect to UTC
           !     in minutes, the hour, the minutes, the seconds, and the
           !     milliseconds of the current date.


           ! fdate: Returns the current date and time as: Day Mon dd hh:mm:ss yyyy.

          character (LEN=30) :: fd_str
          character (LEN=8) :: d_str
          character (LEN=10) :: t_str
          character (LEN=5) :: z_str
          character (LEN=14) :: date_time
          integer, dimension(8) :: dt_values


          call fdate(fd_str)
          call date_and_time(d_str,t_str,z_str,dt_values)

          date_time(1:2) = fd_str(9:10)
          if (dt_values(3) .lt. 10) then
              date_time(1:1) = '0'
          endif
          date_time(3:5) = fd_str(5:7)
          date_time(6:9) = fd_str(21:24)
          date_time(10:10) = '_'
          date_time(11:12) = fd_str(12:13)
          date_time(13:14) = fd_str(15:16)

      end SUBROUTINE sia2_do_datestr




! ======================================================================
! Subroutine: sia2_write_construct_fname
! Purpose: makes output filename
! ======================================================================

      SUBROUTINE sia2_write_construct_fname(fname,fname_short, &
            fname_short_l,ftype,ftype_l,date_time)

          integer :: ftype_l,fname_short_l
          character (LEN=80) :: fname
          character (LEN=20) :: fname_short
          character (LEN=3) :: ftype
          character (LEN=14) :: date_time

          call sia2_blank_line(fname)

          fname(1:16) = "data/output/run_"
          fname(17:30) = date_time
          fname(31:31) = '/'
          fname(32:31+fname_short_l) = fname_short(1:fname_short_l)
          fname(32+fname_short_l:32+fname_short_l) = '_'
          fname(33+fname_short_l:46+fname_short_l) = date_time(1:14)
          fname(47+fname_short_l:47+fname_short_l) = '.'
          fname(48+fname_short_l:47+fname_short_l+ftype_l) = ftype(1:ftype_l)

      end SUBROUTINE sia2_write_construct_fname


! ======================================================================
! Subroutine: sia2_restart
! Purpose: writes/reads working arrays and data to/from restart files
! ======================================================================

      subroutine sia2_restart(m,ice,adv,pur,out_fname,stn_output,cur_step,rw)

          use sia2_globals
          implicit none

          ! included restart file vars
          ! ------------------------------------------------------------
          type (meta_type) :: m(tcells)
          type (ice_type) :: ice(sda_n,ida_n,tcells)
          type (adv_type) :: adv(tcells)
          integer (kind=2) :: pur(z_max,dt_per_day,lda_n,sda_n,ida_n,tcells)
          character (len=80) :: &
              out_fname
          type (stn_write_type) :: stn_output(sda_n,ida_n,n_stations)
          integer :: &
              cur_step, &
              rw

          ! local vars
          ! ------------------------------------------------------------
          logical :: &
              exists
          integer :: &
              sizeof, &
              snow_size, &
              pond_size, &
              adv_ice_size, &
              adv_snow_size, &
              ic, &
              mi, &
              i
          character (len=20) :: restart_fname
          character (len=80) :: restart_bin_fname
          integer :: int_size = 1
          integer :: double_size = 2

          if (rw .eq. 0) then
              print *,'Writing restart files... Step ',cur_step
          else
              print *,'Reading restart files...'
          endif

          snow_size = &                             ! find size of snow structure
            ! type size   ! array size  ! number vars &
              int_size *    1 *           1             + &
              double_size * 1 *           2             + &
              double_size * z_max *       6

          pond_size = &                             ! find size of pond structure
            ! type size   ! array size  ! number vars &
              double_size * 1 *           12

          ! write restart text file
          ! ------------------------------------------------------------
          restart_fname(1:20) = 'restart.params.txt  '       ! contruct filename
      call sia2_write_construct_fname(restart_bin_fname,'restart.params      ',14,'txt',3,datestr)
          if (rw .eq. 0) then
              do i=1,2
                  if (i .eq. 1) then
            open(unit=22, file=trim(restart_fname), form='FORMATTED', access='sequential')
          else
            open(unit=22, file=trim(restart_bin_fname), form='FORMATTED', access='sequential')
                  endif
          write(22,'(i6)') cur_step
          write(22,'(i6)') steps
          write(22,'(i4)') pur_clock
          write(22,'(i4)') cur_year
          write(22,'(d20.14)') cur_hour
          write(22,'(i3)') last_day
          write(22,'(i4)') last3hour
          write(22,'(i4)') last6hour
          write(22,'(i2)') last_5_day
          write(22,'(i4)') last_year
          write(22,'(i4)') last_hour
          write(22,'(d20.14)') next_write_hour
          write(22,'(i6)') write_step
          write(22,'(a80)') out_fname
          write(22,'(a14)') datestr
          close(22)
              enddo
          else
        inquire(file=trim(restart_fname),exist=exists)
              print *,'Reading restart params...'
        if (exists) then
          open(unit=24, file=trim(restart_fname),form='FORMATTED',access='sequential')
          read(24,'(i6)') cur_step
          read(24,'(i6)') steps
          read(24,'(i4)') pur_clock
          read(24,'(i4)') cur_year
          read(24,'(d20.14)') cur_hour
          read(24,'(i3)') last_day
          read(24,'(i4)') last3hour
          read(24,'(i4)') last6hour
          read(24,'(i2)') last_5_day
          read(24,'(i4)') last_year
          read(24,'(i4)') last_hour
          read(24,'(d20.14)') next_write_hour
          read(24,'(i6)') write_step
          read(24,'(a80)') out_fname
          read(24,'(a14)') datestr
          close(24)
        else
          print *, ' Can not find restart file: ',trim(restart_fname)
          call exit(0)
        endif
          endif ! end of read or write check

          ! write out m binary structure (can be humongous)
          ! ------------------------------------------------------------
      call sia2_write_construct_fname(restart_bin_fname,'restart_m           ',9,'bin',3,datestr)

          sizeof = &                                    ! find size of m structure
            ! type size     ! array size  ! number vars &
              int_size *    1 *           11           + &
              int_size *    8 *           3            + &
              double_size * 1 *           46

          sizeof = sizeof+2  ! need extra 2 bytes, at least for ifort

          if (rw .eq. 0) then

        open(unit=22, file=trim(restart_bin_fname), form='UNFORMATTED',  &
!          access='direct', recl=sizeof*tcells)
           access='direct', recl=sizeof)
              do mi=1,tcells
          write(22,rec=mi) m(mi)
          enddo
!       write(22,rec=1) m
        close(22)

          else

        inquire(file=trim(restart_bin_fname),exist=exists)
        if (exists) then
                  print *,'Reading meta...'
          open(unit=24, file=trim(restart_bin_fname), form='UNFORMATTED',  &
!              access='direct', recl=sizeof*tcells)
               access='direct', recl=sizeof)
          do mi=1,tcells
            read(24,rec=mi) m(mi)
          enddo
!         read(24,rec=1) m
          close(24)
        else
          print *, ' Can not find restart file: ',trim(restart_bin_fname)
          call exit(0)
        endif

          endif ! end of read or write check


          ! write out ice binary structure (can be humongous)
          ! ------------------------------------------------------------
      call sia2_write_construct_fname(restart_bin_fname,'restart_ice         ',11,'bin',3,datestr)

          sizeof = &                                    ! find size of ice structure
            ! type size     ! array size  ! number vars &
              int_size *    1 *           1             + &
              int_size *    z_max *       2             + &
              double_size * 1 *           17            + &
              double_size * z_max *       25            + &
              double_size * (z_max+1) *   9             + &
              double_size * z_max * la *  1             + &
              snow_size *   1 *           1             + &
              pond_size *   1 *           1

          sizeof = (sizeof+2)*ida_n*sda_n ! seems like structures need 2 extra 4-byte words, at least using ifort
          if (rw .eq. 0) then

        open(unit=22, file=trim(restart_bin_fname), form='UNFORMATTED',  &
           access='direct', recl=sizeof)
        i = 1
              do mi=1,tcells
!                  do ic=1,ida_n
                      if (m(mi)%status .eq. 1) then
                  write(22,rec=i) ice(:,:,mi)
            i = i+1
                      endif
!             enddo
          enddo
        close(22)

          else

        inquire(file=trim(restart_bin_fname),exist=exists)
        if (exists) then
                  print *,'Reading ice...'
          open(unit=24, file=trim(restart_bin_fname), form='UNFORMATTED',  &
             access='direct', recl=sizeof)
            i = 1
          do mi=1,tcells
!           do ic=1,ida_n
                      if (m(mi)%status .eq. 1) then
              read(24,rec=i) ice(:,:,mi)
                i = i+1
                      endif
!           enddo
          enddo
          close(24)
        else
          print *, ' Can not find restart file: ',trim(restart_bin_fname)
          call exit(0)
        endif

          endif ! end of read or write check


         adv_ice_size = &                             ! find size of adv_ice_type structure
            ! type size   ! array size  ! number vars &
              int_size *    1 *           1             + &
              double_size * 1 *           3             + &
              double_size * z_max *       6

          adv_snow_size = &                            ! find size of adv_snow_type structure
            ! type size   ! array size  ! number vars &
              int_size *    1 *           1             + &
              double_size * 1 *           2             + &
              double_size * z_max *       1

          ! write out adv binary structure (can be humongous)
          ! ------------------------------------------------------------
      call sia2_write_construct_fname(restart_bin_fname,'restart_adv         ',11,'bin',3,datestr)

          sizeof = &                                    ! find size of adv structure
            ! type size     ! array size  ! number vars &
              int_size *    1 *           1             + &
              double_size * 1 *           23            + &
              double_size * 8 *           3             + &
              double_size * ida_n *       1             + &
              double_size * ida_n*ida_n*  1             + &
             adv_ice_size * ida_n*sda_n*  1             + &
            adv_snow_size * ida_n*sda_n*  1

          sizeof = sizeof+8 + ida_n*sda_n*2 ! seems like structures need extra bytes?, at least using ifort
          if (rw .eq. 0) then

        open(unit=22, file=trim(restart_bin_fname), form='UNFORMATTED',  &
           access='direct', recl=sizeof)
        i = 1
              do mi=1,tcells
                      if (m(mi)%status .eq. 1) then
                  write(22,rec=i) adv(mi)
            i = i+1
                      endif
          enddo
        close(22)

          else

        inquire(file=trim(restart_bin_fname),exist=exists)
        if (exists) then
                  print *,'Reading adv...'
          open(unit=24, file=trim(restart_bin_fname), form='UNFORMATTED',  &
             access='direct', recl=sizeof)
            i = 1
          do mi=1,tcells
                      if (m(mi)%status .eq. 1) then
              read(24,rec=i) adv(mi)
                i = i+1
                      endif
          enddo
          close(24)
        else
          print *, ' Can not find restart file: ',trim(restart_bin_fname)
          call exit(0)
        endif

          endif ! end of read or write check


          ! write out PUR binary structure (can REALLY be humongous)
          ! ------------------------------------------------------------
      call sia2_write_construct_fname(restart_bin_fname,'restart_pur         ',11,'bin',3,datestr)

          sizeof = &                                    ! find size of pur array
            ! type size     ! array size  ! number vars &
!             1 * z_max*dt_per_day*lda_n*ida_n*tcells   / 2 + 1  ! for int*1 array
              1 * z_max*dt_per_day*lda_n*sda_n*ida_n / 2               ! for int*2 array

          if (rw .eq. 0) then

        open(unit=22, file=trim(restart_bin_fname), form='UNFORMATTED',  &
           access='direct', recl=sizeof)
        i = 1
        do mi=1,tcells
!         do ic=1,ida_n
                      if (m(mi)%status .eq. 1) then
              write(22,rec=i) pur(:,:,:,:,:,mi)
              i = i+1
                      endif
!         enddo
        enddo
        close(22)

          else

        inquire(file=trim(restart_bin_fname),exist=exists)
        if (exists) then
                  print *,'Reading pur...'
          open(unit=24, file=trim(restart_bin_fname), form='UNFORMATTED',  &
               access='direct', recl=sizeof)
            i = 1
          do mi=1,tcells
!           do ic=1,ida_n
                      if (m(mi)%status .eq. 1) then
              read(24,rec=i) pur(:,:,:,:,:,mi)
                i = i+1
                      endif
!           enddo
          enddo
          close(24)
        else
          print *, ' Can not find restart file: ',trim(restart_bin_fname)
          call exit(0)
        endif

          endif ! end of read or write check

          ! write out station binary structure (can be humongous)
          ! ------------------------------------------------------------
      call sia2_write_construct_fname(restart_bin_fname,'restart_station     ',15,'bin',3,datestr)

          sizeof = &                                    ! find size of station structure
            ! type size     ! array size  ! number vars &
              int_size *    1 *           1            + &  ! represents logical, which is same size as int
              int_size *    1 *           1            + &
              int_size *    20 *          1            + &  ! represents 80 byte (char) string
              double_size * 1 *           18           + &
              double_size * z_max *       30           + &
              double_size * (z_max+1) *   5            + &
              double_size * wavl *        3            + &
              double_size * 3 *           1

          sizeof = (sizeof+2)*ida_n*sda_n ! seems like structures need 2 extra 4-byte words, at least using ifort
          if (rw .eq. 0) then

       open(unit=22, file=trim(restart_bin_fname), form='UNFORMATTED',  &
          access='direct', recl=sizeof)
             do i=1,n_stations
                 write(22,rec=i) stn_output(:,:,i)
             enddo
       close(22)

          else

        inquire(file=trim(restart_bin_fname),exist=exists)
        if (exists) then
                  print *,'Reading stations...'
          open(unit=24, file=trim(restart_bin_fname), form='UNFORMATTED',  &
               access='direct', recl=sizeof)
                  do i=1,n_stations
              read(24,rec=i) stn_output(:,:,i)
          enddo
          close(24)
        else
          print *, ' Can not find restart file: ',trim(restart_bin_fname)
          call exit(0)
        endif

          endif ! end of read or write check



    end subroutine sia2_restart



