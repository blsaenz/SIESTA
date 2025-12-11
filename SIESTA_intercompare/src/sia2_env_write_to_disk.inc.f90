

      if (.not. (start)) then
          ! write model-domain array data to disk
          ! ------------------------------------------------------------
          if (do_write) then

              ! re-initialize output files just in case it need to be re-created
              ! (for instance if it was re-named to avoid gigantic file sizes)
              if (restart .eq. 1) then

				  ! record lat/lon grid
				  call sia2_write_cdf_double(mdh,'x',1,mdv,'y',1,0,' ',1,0,' ',1,2, &
				  d0_,pj%lat(mdh1:mdh2,mdv1:mdv2),d0_,d0_,out_fname,'lats',4, &
				  'degrees',7,'latitude position of EASE grid cell',35)

				  call sia2_write_cdf_double(mdh,'x',1,mdv,'y',1,0,' ',1,0,' ',1,2, &
				  d0_,pj%lon(mdh1:mdh2,mdv1:mdv2),d0_,d0_,out_fname,'lons',4, &
				  'degrees',7,'longitude position of EASE grid cell',36)

              endif

              ! open file for sequential writing of vars to netcdf output file
              call sia2_nc_open(out_fname,jj)

			  ! write biomass grids
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,2),d0_,out_fname,'bm_int',6, &
			  'gC/pixel',8,'Modeled Small Ice Algae Biomass - Internal layers',49,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,1),d0_,out_fname,'bm_bot',6, &
			  'gC/pixel',8,'Modeled Small Ice Algae Biomass - Bottom layers',47,jj)
			  ! write production grids
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,3),d0_,out_fname,'prod_int',8, &
			  'gC/pixel',8,'Modeled Small Ice Algae Production - Internal layers',52,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,4),d0_,out_fname,'prod_bot',8, &
			  'gC/pixel',8,'Modeled Small Ice Algae Production - Bottom layers',50,jj)
			  ! write bm_lost grids
!			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
!			  d0_,d0_,sum_var(:,:,6),d0_,out_fname,'bm_lost_int',11, &
!			  'gC/pixel',8,'Biomass lost to model domain',28,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,5),d0_,out_fname,'bm_lost',7, &
			  'gC/pixel',8,'Biomass lost to model domain',28,jj)
			  ! write ice depth grid
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,7),d0_,out_fname,'ice_depth',9, &
			  'm',1,'Modeled Ice Layer Depth-Down=positive',37,jj)
			  ! write ice depth grid
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,8),d0_,out_fname,'T1',2, &
			  'degC',4,'Top ice layer temp',18,jj)
			  ! write ice concentration grid
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,9),d0_,out_fname,'icecon',6, &
			  '%',1,'percent ice cover',17,jj)
			  ! write snow depth grid
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,10),d0_,out_fname,'snowh',5, &
			  'm',1,'snow depth',10,jj)

			  ! write precipitation comparison grids
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,11),d0_,out_fname,'pr_clim',7, &
			  'm',6,'Climatological Precipitation',28,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,20),d0_,out_fname,'pr_ssmi',7, &
			  'm',6,'SSMI Precipitation',18,jj)

			  ! write vol ice lost due to melting grid
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,12),d0_,out_fname,'melt_loss',9, &
			  'm^3',3,'Volume of sea ice lost due to bottom melt',41,jj)

			  ! write vol ice lost through advection to other cells
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,13),d0_,out_fname,'adv_loss',8, &
			  'm^3',3,'Volume of sea ice lost through advection to other cells',55,jj)

			  ! write vol ice lost to model boundaries
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,14),d0_,out_fname,'md_loss',7, &
			  'm^3',3,'Volume of sea ice lost to model boundaries',42,jj)


			  ! write ice strength
!			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
!			  d0_,d0_,sum_var(:,:,15),d0_,out_fname,'strength',8, &
!			  'N',3,'Sea ice strength',16,jj)

			  ! write ice area
!			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
!			  d0_,d0_,sum_var(:,:,16),d0_,out_fname,'a',1, &
!			  'm^2',3,'area ice area',12,jj)

			  ! write ice area - next time step
!			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
!			  d0_,d0_,sum_var(:,:,17),d0_,out_fname,'a_new',5, &
!			  'm^2',3,'new sea ice area',16,jj)

			  ! write ice out area
!			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
!			  d0_,d0_,sum_var(:,:,18),d0_,out_fname,'a_out',5, &
!			  'm^2',3,'ice area out',12,jj)

			  ! write ice in area
!			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
!			  d0_,d0_,sum_var(:,:,19),d0_,out_fname,'a_in',4, &
!			  'm^2',3,'ice area in',11,jj)

			  ! write surface temp
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,19),d0_,out_fname,'Ts',2, &
			  'degC',4,'Surface temp',12,jj)
			  ! write ice concentration grid

			  ! write vol ice added through congelation growth
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,17),d0_,out_fname,'cong_growth',11, &
			  'm^3',3,'Volume of sea ice added through congelation growth',50,jj)

			  ! write vol ice added thorugh snow ice formation
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,18),d0_,out_fname,'snow_growth',11, &
			  'm^3',3,'Volume of sea ice added thorugh snow ice formation',50,jj)

			  ! write convergece area
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,15),d0_,out_fname,'a_convg',7, &
			  'm^2',3,'convergence area',16,jj)

			  ! write new ice area
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,16),d0_,out_fname,'a_new',5, &
			  'm^2',3,'new ice created',15,jj)


			  ! write vol ice added thorugh advection to other cells
!			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
!			  d0_,d0_,sum_var(:,:,18),d0_,out_fname,'adv_gain',8, &
!			  'm^3',3,'Volume of sea ice added thorugh advection to other cells',56,jj)

			  ! write vol ice ice added due to model boundaries
!			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
!			  d0_,d0_,sum_var(:,:,19),d0_,out_fname,'md_gain',7, &
!			  'm^3',3,'Volume of sea ice added due to model boundaries',47,jj)

			  ! ice vectors
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,21),d0_,out_fname,'uvec',4, &
			  'mm/s',4,'u direction ice vector',22,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,22),d0_,out_fname,'vvec',4, &
			  'mm/s',4,'v direction ice vector',22,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,23),d0_,out_fname,'af_tot',6, &
			  'fraction',8,'actual ice concentration',24,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,24),d0_,out_fname,'a_drop',6, &
			  'km^2',84,'ice dropped from model grid b/c of convergence',67,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,25),d0_,out_fname,'age',3, &
			  'years',5,'mean pixel ice age',18,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,26),d0_,out_fname,'ridged',6, &
			  'fraction',8,'mean pixel ridged ice fraction',30,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,27),d0_,out_fname,'tmean',5, &
			  'degC',4,'biomass-weighted mean temp',26,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,28),d0_,out_fname,'nlim',4, &
			  'fraction',8,'N limitation',12,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,29),d0_,out_fname,'plim',4, &
			  'fraction',8,'P limitation',12,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,30),d0_,out_fname,'silim',5, &
			  'fraction',8,'Si limitation',13,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,31),d0_,out_fname,'slim',4, &
			  'fraction',8,'Salinity limitation',19,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,32),d0_,out_fname,'llim',4, &
			  'fraction',8,'light limitation',16,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,33),d0_,out_fname,'f_flooded',9, &
			  'fraction',8,'fraction of category flooded',28,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,34),d0_,out_fname,'bm_depth',8, &
			  'm',1,'mean biomass depth',18,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,35),d0_,out_fname,'salt_flux',9, &
			  'kg/m^2',6,'salt flux to ocean',18,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,36),d0_,out_fname,'h2o_flux',8, &
			  'kg/m^2',7,'freshwater flux to ocean',35,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,51),d0_,out_fname,'ice_vol',7, &
			  'm^3',3,'total ice volume',16,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,50),d0_,out_fname,'brine_vol',9, &
			  'm^3',3,'total brine volume',18,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,37),d0_,out_fname,'pwid_int',8, &
			  'm',1,'productivity-weighted ice depth-int',40,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,38),d0_,out_fname,'pwid_af_int',11, &
			  'm',1,'productivity-weighted ice mean depth-int',45,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,39),d0_,out_fname,'pwsd_int',8, &
			  'm',1,'productivity-weighted snow depth-int',41,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,40),d0_,out_fname,'pwsd_af_int',11, &
			  'm',1,'productivity-weighted snow mean depth-int',46,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,41),d0_,out_fname,'pwt_int',7, &
			  'degC',4,'productivity-weighted ice temp-int',39,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,42),d0_,out_fname,'pwt_af_int',10, &
			  'deg',4,'productivity-weighted mean ice temp-int',44,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,43),d0_,out_fname,'pwid_bot',8, &
			  'm',1,'productivity-weighted ice depth-bot',40,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,44),d0_,out_fname,'pwid_af_bot',11, &
			  'm',1,'productivity-weighted ice mean depth-bot',45,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,45),d0_,out_fname,'pwsd_bot',8, &
			  'm',1,'productivity-weighted snow depth-bot',41,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,46),d0_,out_fname,'pwsd_af_bot',11, &
			  'm',1,'productivity-weighted snow mean depth-bot',46,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,47),d0_,out_fname,'pwt_bot',7, &
			  'degC',4,'productivity-weighted ice temp-bot',39,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,48),d0_,out_fname,'pwt_af_bot',10, &
			  'deg',4,'productivity-weighted mean ice temp-bot',44,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,49),d0_,out_fname,'s',1, &
			  'psu',3,'mean bulk salinity',18,jj)

			  ! write vol ice added through void spaces during ridging
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,52),d0_,out_fname,'void_growth',11, &
			  'm^3',3,'Volume of sea ice added by voids during ridging',47,jj)


			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,53),d0_,out_fname,'pwt_bot_wgt',11, &
			  'm^2*gC/day',10,'bottom productivity-area mean weight',36,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,54),d0_,out_fname,'pwt_int_wgt',11, &
			  'm^2*gC/day',10,'interior productivity-area mean weight',38,jj)


			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,55),d0_,out_fname,'bm11t',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,56),d0_,out_fname,'bm11m',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,57),d0_,out_fname,'bm11b',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,58),d0_,out_fname,'af11',4, &
			  'fraction',8,'category area',13,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,59),d0_,out_fname,'id11',4, &
			  'm',1,'category ice depth',18,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,60),d0_,out_fname,'bm12t',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,61),d0_,out_fname,'bm12m',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,62),d0_,out_fname,'bm12b',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,63),d0_,out_fname,'af12',4, &
			  'fraction',8,'category area',13,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,64),d0_,out_fname,'id12',4, &
			  'm',1,'category ice depth',18,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,65),d0_,out_fname,'bm13t',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,66),d0_,out_fname,'bm13m',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,67),d0_,out_fname,'bm13b',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,68),d0_,out_fname,'af13',4, &
			  'fraction',8,'category area',13,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,69),d0_,out_fname,'id13',4, &
			  'm',1,'category ice depth',18,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,70),d0_,out_fname,'bm14t',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,71),d0_,out_fname,'bm14m',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,72),d0_,out_fname,'bm14b',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,73),d0_,out_fname,'af14',4, &
			  'fraction',8,'category area',13,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,74),d0_,out_fname,'id14',4, &
			  'm',1,'category ice depth',18,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,75),d0_,out_fname,'bm15t',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,76),d0_,out_fname,'bm15m',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,77),d0_,out_fname,'bm15b',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,78),d0_,out_fname,'af15',4, &
			  'fraction',8,'category area',13,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,79),d0_,out_fname,'id15',4, &
			  'm',1,'category ice depth',18,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,80),d0_,out_fname,'bm21t',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,81),d0_,out_fname,'bm21m',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,82),d0_,out_fname,'bm21b',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,83),d0_,out_fname,'af21',4, &
			  'fraction',8,'category area',13,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,84),d0_,out_fname,'id21',4, &
			  'm',1,'category ice depth',18,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,85),d0_,out_fname,'bm22t',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,86),d0_,out_fname,'bm22m',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,87),d0_,out_fname,'bm22b',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,88),d0_,out_fname,'af22',4, &
			  'fraction',8,'category area',13,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,89),d0_,out_fname,'id22',4, &
			  'm',1,'category ice depth',18,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,90),d0_,out_fname,'bm23t',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,91),d0_,out_fname,'bm23m',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,92),d0_,out_fname,'bm23b',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,93),d0_,out_fname,'af23',4, &
			  'fraction',8,'category area',13,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,94),d0_,out_fname,'id23',4, &
			  'm',1,'category ice depth',18,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,95),d0_,out_fname,'bm24t',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,96),d0_,out_fname,'bm24m',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,97),d0_,out_fname,'bm24b',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,98),d0_,out_fname,'af24',4, &
			  'fraction',8,'category area',13,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,99),d0_,out_fname,'id24',4, &
			  'm',1,'category ice depth',18,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,100),d0_,out_fname,'bm25t',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,101),d0_,out_fname,'bm25m',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,102),d0_,out_fname,'bm25b',5, &
			  'gC/m^2',6,'category biomass',16,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,103),d0_,out_fname,'af25',4, &
			  'fraction',8,'category area',13,jj)
			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,104),d0_,out_fname,'id25',4, &
			  'm',1,'category ice depth',18,jj)

			  call sia2_write_open_cdf_double(mdh,'x',1,mdv,'y',1,write_step,'time',4,0,' ',1,3, &
			  d0_,d0_,sum_var(:,:,105),d0_,out_fname,'ohf',3, &
			  'W m-2',1,'ocean heat flux',15,jj)



              ! close netcdf
              call sia2_nc_close(out_fname,jj)

              ! turn off write-to-disk
              do_write = .false.
          endif

          ! write stations data to disk
          ! ------------------------------------------------------------
          if (wr_stations .ne. 0) then
          do ii=1,n_stations
          do ic=1,ida_n
          do sc=1,sda_n

          if (st_out(sc,ic,ii)%valid) then

              ! open file for sequential writing of vars to netcdf station file
              call sia2_nc_open(st_out(sc,ic,ii)%fname,jj)


              ! Writeout timing

			  call sia2_write_open_cdf_double(st_out(sc,ic,ii)%step,'time',4,0,' ',1,0,' ',1,0,' ',1,1, &
			  st_out(sc,ic,ii)%time,d0_,d0_,d0_,st_out(sc,ic,ii)%fname,'time',4,'hours',5,'Time',4,jj)


              ! Writeout 1D (single) station vars

			  ! write out icecon
			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(1),d0_,st_out(sc,ic,ii)%fname,'ease_icecon',11, &
			  '%',1,'percent ice cover',17,jj)

			  ! write out icecon
			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(2),d0_,st_out(sc,ic,ii)%fname,'icecon',6, &
			  '%',1,'percent ice cover',17,jj)

			   ! write out ssmi snow depth
			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(3),d0_,st_out(sc,ic,ii)%fname,'snowh',5, &
			  'cm',2,'snow depth over sea ice',23,jj)

			   ! write out snow depth
			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(4),d0_,st_out(sc,ic,ii)%fname,'sh_prev',7, &
			  'm',1,'snow depth over sea ice',23,jj)

              ! find and write out mean ice depth
			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(5),d0_,st_out(sc,ic,ii)%fname,'id_mean',7, &
			  'm',1,'weighted mean ice depth of ice categories',41,jj)

			   ! write out snow layers
			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(6),d0_,st_out(sc,ic,ii)%fname,'z_snow',6, &
			  'none',4,'model snow layers',17,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(7),d0_,st_out(sc,ic,ii)%fname,'airtemp',7, &
			  'degK',4,'2m Air Temperature',18,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(8),d0_,st_out(sc,ic,ii)%fname,'precip',6, &
			  'kg/m^2',6,'Precipitation Rate',18,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(9),d0_,st_out(sc,ic,ii)%fname,'Ts',2, &
			  'degC',4,'Model Surface Temperature',25,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(10),d0_,st_out(sc,ic,ii)%fname,'Ed0',3, &
			  'uEin/m^2/s',10,'Surface ice/snow irradience',27,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(11),d0_,st_out(sc,ic,ii)%fname,'airtemp_interp',14, &
			  'degC',4,'Interpolated Air Temperature',28,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(12),d0_,st_out(sc,ic,ii)%fname,'status',6, &
			  'none',4,'Pixel Status Code',17,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(13),d0_,st_out(sc,ic,ii)%fname,'z_ice',5, &
			  'none',4,'model ice layers',16,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(14),d0_,st_out(sc,ic,ii)%fname,'lats',4, &
			  'degrees',7,'latitude position of EASE grid cell',35,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(15),d0_,st_out(sc,ic,ii)%fname,'lons',4, &
			  'degrees',7,'longitude position of EASE grid cell',36,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(16),d0_,st_out(sc,ic,ii)%fname,'PAR_bot',7, &
			  'uEin/m^2/s',10,'PAR at ice bottom',17,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(17),d0_,st_out(sc,ic,ii)%fname,'PAR_bot_pl',10, &
			  'uEin/m^2/s',10,'PAR at platelet bottom',22,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,st_out(sc,ic,ii)%step,'time',4,0,' ',1,3, &
			  d0_,d0_,st_out(sc,ic,ii)%s(18),d0_,st_out(sc,ic,ii)%fname,'af',2, &
			  'fraction',8,'pixel areal fraction',28,jj)

              ! Writeout z-dimensional station vars


			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(1,:),st_out(sc,ic,ii)%fname,'th_snow',7, &
			  'm',1,'Modeled Snow Layer Thickness',28,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(2,:),st_out(sc,ic,ii)%fname,'t_snow',6, &
			  'degC',4,'Modeled Snow Layer Temperature',33,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(3,:),st_out(sc,ic,ii)%fname,'d_snow',6, &
			  'g/m^3',5,'Modeled Snow Layer Density',29,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(4,:),st_out(sc,ic,ii)%fname,'heat_snow',9, &
			  '???',3,'Modeled Snow Layer Heat',29,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(5,:),st_out(sc,ic,ii)%fname,'smalg',5, &
			  'mgC/m^3',7,'Modeled Small Ice Algae Concentration',37,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(6,:),st_out(sc,ic,ii)%fname,'prod',4, &
			  'gC',2,'Ice Algae Production',20,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(7,:),st_out(sc,ic,ii)%fname,'heat',4, &
			  'J/g',3,'enthalpy required to heat layer to 0degC',40,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(8,:),st_out(sc,ic,ii)%fname,'thickness',9, &
			  'm',1,'Modeled Ice Layer Thickness',27,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(9,:),st_out(sc,ic,ii)%fname,'Ice_depth',9, &
			  'm',1,'Modeled Ice Layer Depth-Down=positive',37,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(10,:),st_out(sc,ic,ii)%fname,'T_layer',7, &
			  'degC',4,'Modeled Ice Layer Temperature',32,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(11,:),st_out(sc,ic,ii)%fname,'PUR',3, &
			  'uEin/m^2/s',10,'Photosynthetically Useable Radiation',36,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(12,:),st_out(sc,ic,ii)%fname,'Ik_prime',8, &
			  'uEin/m^2/s',10,'Spectral Photo Adaptation Parameter',35,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(13,:),st_out(sc,ic,ii)%fname,'brine_v',7, &
			  'ppt',3,'Brine Volume',12,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(14,:),st_out(sc,ic,ii)%fname,'brine_d',7, &
			  'g/m^3',5,'Brine Density',13,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(15,:),st_out(sc,ic,ii)%fname,'brine_sal',9, &
			  'ppt',3,'Brine Salinity',14,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(16,:),st_out(sc,ic,ii)%fname,'S_layer',7, &
			  'psu',3,'Sea Ice Salinity',16,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(17,:),st_out(sc,ic,ii)%fname,'ice_d',5, &
			  'g/m^3',3,'Sea Ice Density',15,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(18,:),st_out(sc,ic,ii)%fname,'llim',4, &
			  'fraction',8,'Light Limitation',15,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(19,:),st_out(sc,ic,ii)%fname,'nlim',4, &
			  'fraction',8,'N Limitation',12,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(20,:),st_out(sc,ic,ii)%fname,'plim',4, &
			  'fraction',8,'P Limitation',12,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(21,:),st_out(sc,ic,ii)%fname,'silim',5, &
			  'fraction',8,'Si Limitation',13,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(22,:),st_out(sc,ic,ii)%fname,'slim',4, &
			  'fraction',8,'Salinity Limitation',19,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(23,:),st_out(sc,ic,ii)%fname,'melt',4, &
			  'fraction',8,'melted, large grain snow fraction',33,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(24,:),st_out(sc,ic,ii)%fname,'dsdt',4, &
			  'ppt',3,'salinity lost',15,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(25,:),st_out(sc,ic,ii)%fname,'fbv',3, &
			  'm^3/m^2',7,'brine volume flux',17,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(26,:),st_out(sc,ic,ii)%fname,'dhdt_conv',9, &
			  'm',1,'estimated new internal ice',26,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(27,:),st_out(sc,ic,ii)%fname,'f0',2, &
			  'W/m^2',5,'estimated new ice flux',22,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(28,:),st_out(sc,ic,ii)%fname,'dsdt3',5, &
			  'ppt',3,'salinity lost',15,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(29,:),st_out(sc,ic,ii)%fname,'tgrad',5, &
			  'degC/m',6,'temperature gradient',20,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max,'z',1,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z(30,:),st_out(sc,ic,ii)%fname,'gmax',4, &
			  '1/day',5,'maximum daily growth rate',25,jj)

              ! Writeout z1-dimensional station vars


			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max+1,'z1',2,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z1(1,:),st_out(sc,ic,ii)%fname,'NO3',3, &
			  'mmol/m^3',8,'Modeled Nitrate Concentration',29,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max+1,'z1',2,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z1(2,:),st_out(sc,ic,ii)%fname,'NH4',3, &
			  'mmol/m^3',8,'Modeled Ammonia Concentration',29,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max+1,'z1',2,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z1(3,:),st_out(sc,ic,ii)%fname,'PO4',3, &
			  'mmol/m^3',8,'Modeled Phosphate Concentration',31,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max+1,'z1',2,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z1(4,:),st_out(sc,ic,ii)%fname,'SiOH4',5, &
			  'mmol/m^3',8,'Modeled Silicate Concentration',30,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,z_max+1,'z1',2,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%z1(5,:),st_out(sc,ic,ii)%fname,'poc',3, &
			  'mg/m^3',6,'Detrital Pool',13,jj)


              ! Writeout wavl-dimension station vars


			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,wavl,'wavl',4,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%wavl(1,:),st_out(sc,ic,ii)%fname,'Edif',4, &
			  'uEin/m^2/s/nm',13,'diffuse irradiance before reflection',36,jj)

			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,wavl,'wavl',4,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%wavl(2,:),st_out(sc,ic,ii)%fname,'Edir',4, &
			  'uEin/m^2/s/nm',13,'direct irradiance before reflection',35,jj)

!			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,wavl,'wavl',4,st_out(sc,ic,ii)%step,'time',4,4, &
!			  d0_,d0_,d0_,st_out(sc,ic,ii)%wavl(3,:),st_out(sc,ic,ii)%fname,'PAR_bot',7, &
!			  'uEin/m^2/s/nm',13,'PAR at ice bottom',17,jj)


              ! Writeout icevecs


			  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,3,'ivd',3,st_out(sc,ic,ii)%step,'time',4,4, &
			  d0_,d0_,d0_,st_out(sc,ic,ii)%icevec,st_out(sc,ic,ii)%fname,'ice_vectors',11, &
			  'cm/s',4,'u_magnitude,v_magnitude,validity',32,jj)


              ! Writeout platelet layers if used

			  if (use_pl .eq. 1) then

				  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,pl_max,'pl',2,st_out(sc,ic,ii)%step,'time',4,4, &
				  d0_,d0_,d0_,pl%no3,st_out(sc,ic,ii)%fname,'plno3',5, &
				  'mmol/m^3',8,'Modeled Nitrate Concentration',29,jj)

				  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,pl_max,'pl',2,st_out(sc,ic,ii)%step,'time',4,4, &
				  d0_,d0_,d0_,pl%nh4,st_out(sc,ic,ii)%fname,'plnh4',5, &
				  'mmol/m^3',8,'Modeled Ammonia Concentration',29,jj)

				  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,pl_max,'pl',2,st_out(sc,ic,ii)%step,'time',4,4, &
				  d0_,d0_,d0_,pl%po4,st_out(sc,ic,ii)%fname,'plpo4',5, &
				  'mmol/m^3',8,'Modeled Phosphate Concentration',31,jj)

				  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,pl_max,'pl',2,st_out(sc,ic,ii)%step,'time',4,4, &
				  d0_,d0_,d0_,pl%sioh4,st_out(sc,ic,ii)%fname,'plsioh4',7, &
				  'mmol/m^3',8,'Modeled Silicate Concentration',30,jj)

				  call sia2_write_open_cdf_double(1,'x',1,1,'y',1,pl_max,'pl',2,st_out(sc,ic,ii)%step,'time',4,4, &
				  d0_,d0_,d0_,pl%smalg,st_out(sc,ic,ii)%fname,'plsmalg',7, &
				  'mgC/m^3',7,'Modeled Small Ice Algae Concentration',37,jj)

			  endif

               ! close netcdf station
              call sia2_nc_close(st_out(sc,ic,ii)%fname,jj)

          endif
          enddo
          enddo
          enddo
          endif
    endif
