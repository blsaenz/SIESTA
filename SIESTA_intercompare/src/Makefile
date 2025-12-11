lib_path = ~/Dropbox/Projects/KPP_ECO_ICE/support/netcdf-3.6.3-gfortran_64_mac_10.11/lib
incl_path = ~/Dropbox/Projects/KPP_ECO_ICE/support/netcdf-3.6.3-gfortran_64_mac_10.11/include


opt1 = -g -O0 -fbounds-check -m64
opt2 = -g -O0 -fbounds-check -m64

optlink = -L/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/Current -lBLAS -lLAPACK $(lib_path)/libnetcdf.a 

ofile = sia2_globals.o  sia2_adv_3_setup.o sia2_adv_3_tracers.o sia2_adv_3_grid_ridge.o sia2_env_create_ice_sub.o sia2_env_edd_solution.o sia2_env.o sia2_netcdf_funcs.o sia2_interpolate_forcing.o sia2_forcing.o sia2_write_output.o sia2_initialize.o sia2_read_parameters.o sia2_read_pj.o sia2_timing.o seaice2.o

seaice2:		$(ofile)
	gfortran -o sia2temp $(opt1) $(ofile) $(optlink) 
seaice2.o:
	gfortran -c $(opt1) seaice2.f90

sia2_timing.o:
	gfortran -c $(opt1) sia2_timing.f90
sia2_read_pj.o:
	gfortran -c $(opt1) sia2_read_pj.f90
sia2_read_parameters.o:
	gfortran -c $(opt1) sia2_read_parameters.f90
sia2_initialize.o:
	gfortran -c $(opt1) sia2_initialize.f90
sia2_forcing.o:
	gfortran -c $(opt1) sia2_forcing.f90
sia2_interpolate_forcing.o:
	gfortran -c $(opt1) sia2_interpolate_forcing.f90
sia2_adv_3_setup.o:
	gfortran -c $(opt1) sia2_adv_3_setup.f90
sia2_adv_3_tracers.o:
	gfortran -c $(opt1) sia2_adv_3_tracers.f90
sia2_adv_3_grid_ridge.o:
	gfortran -c $(opt1) sia2_adv_3_grid_ridge.f90
sia2_env.o:
	gfortran -c $(opt1) sia2_env.f90
sia2_env_create_ice_sub.o:
	gfortran -c $(opt1) sia2_env_create_ice_sub.f90
sia2_env_edd_solution.o:
	gfortran -c $(opt2) sia2_env_edd_solution.f90
sia2_write_output.o:
	gfortran -c $(opt1) sia2_write_output.f90
sia2_netcdf_funcs.o:
	gfortran -c $(opt1) sia2_netcdf_funcs.f90
sia2_globals.o:
	gfortran -c $(opt1) sia2_globals.f90

pre:
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_env_update_skeletal.inc.F sia2_env_update_skeletal.inc.f90

	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_env_ice_remap.inc.F sia2_env_ice_remap.inc.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_env_snow_remap.inc.F sia2_env_snow_remap.inc.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_env_pond_mix3.inc.F sia2_env_pond_mix3.inc.f90

	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_env_redist_new_ice.inc.F sia2_env_redist_new_ice.inc.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_env_desalinate_7.inc.F sia2_env_desalinate_7.inc.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_env_flood.inc.F sia2_env_flood.inc.f90

	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_env_grid_remap.inc.F sia2_env_grid_remap.inc.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_env_grid_calcs.inc.F sia2_env_grid_calcs.inc.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_env_flux_heat.inc.F sia2_env_flux_heat.inc.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_env_redist.inc.F sia2_env_redist.inc.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_adv_3_ridge.inc.F sia2_adv_3_ridge.inc.f90

	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_env_create_ice_sub.F sia2_env_create_ice_sub.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path)  sia2_env_icemodel_2.inc.F sia2_env_icemodel_2.inc.f90

	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path)  sia2_edd_gpu_calc.F sia2_edd_gpu_calc.f90

	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_netcdf_funcs.F sia2_netcdf_funcs.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_write_output.F sia2_write_output.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_env.F sia2_env.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path)  sia2_adv_3_grid_ridge.F sia2_adv_3_grid_ridge.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_adv_3_setup.F sia2_adv_3_setup.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path)  sia2_adv_3_tracers.F sia2_adv_3_tracers.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path) sia2_forcing.F sia2_forcing.f90
	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path)  sia2_read_pj.F sia2_read_pj.f90



	/opt/local/bin/cpp-mp-5 -w -P -C -I$(incl_path)   seaice2.F seaice2.f90

clean:
	rm -f *.o *.mod seaice2.f90 sia2_env_flood.inc.f90 sia2_env_redist_new_ice.inc.f90 sia2_env_redist.inc.f90 sia2_adv_3_ridge.inc.f90 sia2_env_pond_mix3.inc.f90 sia2_env_create_ice_sub.f90 sia2_env_desalinate_7.inc.f90 sia2_netcdf_funcs.f90 sia2_write_output.f90 sia2temp sia2_env.f90 sia2_env_icemodel_2.inc.f90 sia2_adv_3_grid_ridge.f90 sia2_adv_3_tracers.f90 sia2_adv_3_setup.f90 sia2_env_create_ice.inc.f90 sia2_env_flux_heat.inc.f90 sia2_env_grid_remap.inc.f90 sia2_env_grid_calcs.inc.f90 sia2_env_ice_remap.inc.f90 sia2_env_snow_remap.inc.f90 sia2_env_desalinate_01.inc.f90 sia2_env_desalinate_2.inc.f90 sia2_env_desalinate_3.inc.f90 sia2_forcing.f90 sia2_read_pj.f90 sia2_env_update_skeletal.inc.f90 sia2_edd_gpu_calc.f90

