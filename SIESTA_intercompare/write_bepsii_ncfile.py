import os.path

import netCDF4
import numpy as np
import xarray
from netCDF4 import Dataset

def siesta_data_fix(var, z_ice):
    """
    Simulates the siesta_data_fix function from Matlab.
    Sets values to 0.0 for layers beyond the specified ice depth.

    Args:
        var (np.ndarray): The input data array (layers x time).
        z_ice (np.ndarray): The ice depth indices for each time step.

    Returns:
        np.ndarray: The fixed data array.
    """
    # Assuming var is (layers, steps) and z_ice is (steps,)
    steps, layers = var.shape
    for ii in range(steps):
        if z_ice[ii] < layers:
            # Set values beyond z_ice to 0.0
            var[ii,int(z_ice[ii]):] = 0.0
    return var

def sampledaily(var):
    """
    Simulates the sampledaily function from Matlab.
    Takes the first hourly sample for each day.

    Args:
        var (np.ndarray): The input data array.

    Returns:
        np.ndarray: The daily sampled data.
    """
    if var.ndim == 2:
        # Assuming var is (layers, time)
        time_steps, layers = var.shape
        days = time_steps // 24
        # Select the first sample of each 24-hour block
        md = var[::24,:]
    else:
        # Assuming var is (time,)
        time_steps = var.size
        days = time_steps // 24
        md = var[::24]
    return md

def meandaily(var):
    """
    Simulates the meandaily function from Matlab.
    Calculates the mean of each 24-hour period.

    Args:
        var (np.ndarray): The input data array.

    Returns:
        np.ndarray: The daily mean data.
    """
    if var.ndim == 2:
        time_steps ,layers = var.shape
        days = time_steps // 24
        # Reshape to (layers, days, 24) and sum over the last axis
        md = var.reshape(days, 24, layers).mean(axis=1)
    else:
        time_steps = var.size
        days = time_steps // 24
        # Reshape to (days, 24) and take the mean over the last axis
        md = var.reshape(days, 24).mean(axis=1)
    return md

def sumdaily(var):
    """
    Simulates the sumdaily function from Matlab.
    Calculates the sum of each 24-hour period.

    Args:
        var (np.ndarray): The input data array.

    Returns:
        np.ndarray: The daily sum data.
    """
    if var.ndim == 2:
        time_steps ,layers = var.shape
        days = time_steps // 24
        # Reshape to (layers, days, 24) and sum over the last axis
        md = var.reshape(days, 24, layers).sum(axis=1)
    else:
        time_steps = var.size
        days = time_steps // 24
        # Reshape to (days, 24) and sum over the last axis
        md = var.reshape(days, 24).sum(axis=1)
    return md

def snow_depth(station):
    """
    Simulates the snow_depth function from Matlab.
    Calculates the sum of snow thickness for each time step.

    Args:
        station (dict): The station data dictionary.

    Returns:
        np.ndarray: The daily snow depth.
    """
    th_snow = station['th_snow']
    # sum along the layers dimension, up to the snow depth index
    sd = np.array([np.sum(th_snow[ii,:int(station['z_ice'][ii])]) for ii in range(th_snow.shape[0])])
    return sd

def ice_depth(station):
    """
    Simulates the ice_depth function from Matlab.
    Calculates the sum of ice thickness for each time step.

    Args:
        station (dict): The station data dictionary.

    Returns:
        np.ndarray: The daily ice depth.
    """
    th = station['thickness']
    # sum along the layers dimension, up to the ice depth index
    id = np.array([np.sum(th[ii,:int(station['z_ice'][ii])]) for ii in range(th.shape[0])])
    return id

# Define variables and their attributes
vars2D = [
    ('z_u', 'upper coordinate of ice layers', 'm'),
    ('z_l', 'lower coordinate of ice layers', 'm'),
    ('z_i', 'mid-point of ice layers', 'm'),
    ('t_i', 'ice temperature', '°C'),
    ('s_i', 'ice salinity', 'g/kg'),
    ('bvf', 'brine volume fraction', '%'),
    ('sbr', 'brine salinity', 'g/kg'),
    ('PUR', 'Photosynthetically Usable Radiation', 'µmol photons/m2/s'),
    ('chla_alg', 'total ice algal chlorophyll', 'mg C/m3'),
    ('c_alg', 'total ice algal organic carbon', 'mg C/m3'),
    ('c_det', 'total detrital organic carbon', 'mg C/m3'),
    ('din', 'dissolved inorganic nitrogen (no3+no2+nh4) in sea ice', 'mmol/m3'),
    ('dip', 'dissolved phosphate in sea ice', 'mmol/m3'),
    ('dsi', 'dissolved silica in sea ice', 'mmol/m3'),
    ('lim_lig', 'mean daily light limitation factor', 'fractional'),
    ('lim_nut', 'mean daily ice nutrient limitation factor', 'fractional'),
    ('lim_sal', 'mean daily ice salinity inhibition factor', 'fractional'),
    ('lim_n', 'mean daily ice nitrogen limitation', 'fractional'),
    ('lim_p', 'mean daily ice phosphate limiation', 'fractional'),
    ('lim_si', 'mean daily ice silicate limitation', 'fractional')
]

vars1D = [
    ('doy', 'day of year', 'days'),
    ('PAR_ml', 'mean daily incident PAR at ocean surface', 'µmol photons/m2/s'),
    ('chla_ml', 'ocean mixed layer chl-a (fixed)', 'mg/m3'),
    ('c_alg_ml', 'ocean mixed layer algal carbon (fixed)', 'mgC/m3'),
    ('NCP_ice', 'net community production in sea ice', 'mgC/m2/day'),
    ('h_i', 'ice thickness', 'm'),
    ('h_s', 'snow depth', 'm'),
    ('Nutr1p', 'phosphate', 'mmol P m-2'),
    ('Nutr2n', 'nitrate', 'mmol N m-2'),
    ('Nutr3s', 'silicate', 'mmol Sim-2'),
    ('Algae1p', 'P in sea ice algae', 'mmol P m-2'),
    ('Algae1n', 'N in sea ice algae', 'mmol N m-2'),
    ('Algae1s', 'Si in sea ice algae', 'mmol Si m-2'),
    ('Algae1c', 'C in sea ice algae', 'mg C m-2'),
    ('Algae1l', 'Chl-a in sea ice algae', 'mg Chl m-2'),
    ('t_surface', 'ice/snow upper surface temperature', '°C')

]


def write_bepsii_ncfile(station_ncfile, doy_start, ncfile):
    """
    Translates the main Matlab function to Python.
    Creates and populates a NetCDF file with BEPSII data.

    Args:
        station (dict): A dictionary containing all the station data.
        doy_start (int): The starting day of year.
        ncfile (str): The name of the output NetCDF file.
    """

    nz = 42

    station_file = netCDF4.Dataset(station_ncfile)
    station_vars = ['thickness','Ice_depth','z_ice','T_layer','S_layer','brine_v','brine_sal','PUR',
                    'smalg','poc','NO3','NH4','PO4','SiOH4','prod','llim','nlim','silim','plim','slim',
                    'th_snow','z_snow','Ed0','t_snow','Ts']
    #for v in station_vars:
    #    print(station_file.variables[v])
    station = {v:station_file.variables[v][...].squeeze() for v in station_vars}


    # Pre-process data
    th = sampledaily(siesta_data_fix(station['thickness'], station['z_ice']))
    bepd = {}

    bepd['z_l'] = sampledaily(siesta_data_fix(station['Ice_depth'], station['z_ice']))
    bepd['z_u'] = np.zeros_like(bepd['z_l'])
    bepd['z_u'][:, 1:] = bepd['z_l'][:, :-1]
    bepd['z_u'][bepd['z_l']==0.0] = 0.0
    bepd['z_i'] = bepd['z_l'] - th / 2

    bepd['h_i'] = sampledaily(ice_depth(station))
    bepd['h_s'] = sampledaily(snow_depth(station))
    bepd['t_surface'] = meandaily(station['Ts'])
    bepd['t_i'] = sampledaily(siesta_data_fix(station['T_layer'], station['z_ice']))
    #bepd['t_snow'] = sampledaily(siesta_data_fix(station['t_snow'], station['z_snow']))
    bepd['s_i'] = sampledaily(siesta_data_fix(station['S_layer'], station['z_ice']))
    bepd['bvf'] = sampledaily(siesta_data_fix(station['brine_v'], station['z_ice'])) / 1000
    bepd['sbr'] = sampledaily(siesta_data_fix(station['brine_sal'], station['z_ice']))
    bepd['PUR'] = sampledaily(siesta_data_fix(station['PUR'], station['z_ice']))
    # bepd['PAR_bot'] = sampledaily(station['PAR_bot']) # Not in vars1D/2D

    smalg = siesta_data_fix(station['smalg'] * station['brine_v'] / 1000, station['z_ice'])
    chla = smalg / 35
    bepd['chla_alg'] = sampledaily(chla)
    bepd['c_alg'] = sampledaily(siesta_data_fix(station['smalg'] * station['brine_v'] / 1000, station['z_ice']))
    #bepd['c_det'] = sampledaily(siesta_data_fix(station['poc'][:, :, :42, :] * station['brine_v'] / 1000, station['z_ice']))
    bepd['c_det'] = sampledaily(siesta_data_fix(station['poc'][:,:42] * station['brine_v'] / 1000, station['z_ice']))

    #din = sampledaily(siesta_data_fix(station['NO3'][:, :, :42, :] * station['brine_v'] / 1000, station['z_ice']))
    #din += sampledaily(siesta_data_fix(station['NH4'][:, :, :42, :] * station['brine_v'] / 1000, station['z_ice']))
    din = sampledaily(siesta_data_fix(station['NO3'][:,:42] * station['brine_v'] / 1000, station['z_ice']))
    din += sampledaily(siesta_data_fix(station['NH4'][:,:42] * station['brine_v'] / 1000, station['z_ice']))
    bepd['din'] = din

    #bepd['dip'] = sampledaily(siesta_data_fix(station['PO4'][:, :, :42, :] * station['brine_v'] / 1000, station['z_ice']))
    #bepd['dsi'] = sampledaily(siesta_data_fix(station['SiOH4'][:, :, :42, :] * station['brine_v'] / 1000, station['z_ice']))
    bepd['dip'] = sampledaily(siesta_data_fix(station['PO4'][:,:42] * station['brine_v'] / 1000, station['z_ice']))
    bepd['dsi'] = sampledaily(siesta_data_fix(station['SiOH4'][:,:42] * station['brine_v'] / 1000, station['z_ice']))

    bepd['NCP_ice'] = sumdaily(siesta_data_fix(station['prod'] / 628380809.625625, station['z_ice'])).sum(axis=1) * 1000
    bepd['lim_lig'] = meandaily(siesta_data_fix(station['llim'], station['z_ice']))
    bepd['lim_n'] = meandaily(siesta_data_fix(station['nlim'], station['z_ice']))
    bepd['lim_p'] = meandaily(siesta_data_fix(station['plim'], station['z_ice']))
    bepd['lim_si'] = meandaily(siesta_data_fix(station['silim'], station['z_ice']))

    bepd['lim_nut'] = np.minimum(bepd['lim_n'], bepd['lim_p'])
    bepd['lim_nut'] = np.minimum(bepd['lim_nut'], bepd['lim_si'])

    bepd['lim_sal'] = meandaily(siesta_data_fix(station['slim'], station['z_ice']))
    bepd['PAR_ml'] = meandaily(station['Ed0'])

    # Nutrient and algae stocks
    #nutr1p_data = siesta_data_fix(station['PO4'][:, :, :42, :] * station['brine_v'] / 1000 * station['thickness'], station['z_ice'])
    nutr1p_data = siesta_data_fix(station['PO4'][:,:42] * station['brine_v'] / 1000 * station['thickness'], station['z_ice'])
    bepd['Nutr1p'] = sumdaily(nutr1p_data).sum(axis=1)

    #nutr2n_data = siesta_data_fix(station['NO3'][:, :, :42, :] * station['brine_v'] / 1000 * station['thickness'], station['z_ice'])
    nutr2n_data = siesta_data_fix(station['NO3'][:,:42] * station['brine_v'] / 1000 * station['thickness'], station['z_ice'])
    bepd['Nutr2n'] = sumdaily(nutr2n_data).sum(axis=1)

    #nutr3s_data = siesta_data_fix(station['SiOH4'][:, :, :42, :] * station['brine_v'] / 1000 * station['thickness'], station['z_ice'])
    nutr3s_data = siesta_data_fix(station['SiOH4'][:,:42] * station['brine_v'] / 1000 * station['thickness'], station['z_ice'])
    bepd['Nutr3s'] = sumdaily(nutr3s_data).sum(axis=1)

    bepd['c_alg'] = sampledaily(siesta_data_fix(station['smalg'] * station['brine_v'] / 1000, station['z_ice']))
    algae1c_data = siesta_data_fix(station['smalg'] * station['brine_v'] / 1000 * station['thickness'], station['z_ice'])
    bepd['Algae1c'] = sampledaily(algae1c_data).sum(axis=1)

    bepd['Algae1p'] = bepd['Algae1c'] / 106 / 12.01 # mmol / m2
    bepd['Algae1n'] = bepd['Algae1c'] / 7 / 12.01 # mmol / m2
    bepd['Algae1s'] = bepd['Algae1c'] / 4 / 12.01 # mmol / m2
    bepd['Algae1l'] = bepd['Algae1c'] / 35. # mg

    time_l = bepd['h_i'].size
    bepd['doy'] = np.arange(1, time_l + 1) + doy_start - 1

    # Fixed values
    bepd['chla_ml'] = np.full(time_l, 0.1)
    bepd['c_alg_ml'] = bepd['chla_ml'] * 35

    # Create NetCDF file
    with Dataset(ncfile, 'w', format='NETCDF4') as ncid:
        print(f"Creating dimensions...")
        tdim_id = ncid.createDimension('time', time_l)
        zdim_id = ncid.createDimension('z', nz)

        # Create 1D variables
        print("Creating 1D variables...")
        for name, long_name, units in vars1D:
            varid = ncid.createVariable(name, 'f8', ('time',))
            varid.long_name = long_name
            varid.units = units

        # Create 2D variables
        print("Creating 2D variables...")
        for name, long_name, units in vars2D:
            varid = ncid.createVariable(name, 'f8', ('time','z'))
            varid.long_name = long_name
            varid.units = units

        # Write data
        print("Writing 1D variables...")
        for name, _, _ in vars1D:
            print(f"  - Writing {name}...")
            ncid.variables[name][:] = bepd[name]

        print("Writing 2D variables...")
        for name, _, _ in vars2D:
            print(f"  - Writing {name}...")
            # Transpose to match (z, time) dimension order if needed
            ncid.variables[name][:] = bepd[name]

    print(f"Successfully created and wrote data to {ncfile}")


if __name__ == '__main__':

    # Call the main function
    # station_nc = "X:/data/SIESTA_intercompare/output/run_11Jun2019_0844_NICE_GC_SWD_Conform_3p5_wc_walb_0p508_algmig2_both_directions/station_11_999-999_11Jun2019_0844.nc"
    # cwd,_ = os.path.split(station_nc)
    # output_filename = os.path.join(cwd,'EXP2_NICE_Saenz_2025-08-07.nc')
    # write_bepsii_ncfile(station_nc, 113, output_filename)
    # exit()

    # f0 = "X:/data/SIESTA_intercompare/output/run_11Jun2019_0844_NICE_GC_SWD_Conform_3p5_wc_walb_0p508_algmig2_both_directions/EXP2_NICE_Saenz.nc"
    # f1 = "X:/data/SIESTA_intercompare/output/run_11Jun2019_0844_NICE_GC_SWD_Conform_3p5_wc_walb_0p508_algmig2_both_directions/EXP2_NICE_Saenz_2025-08-07.nc"
    #
    # import matplotlib.pyplot as plt
    # ds0 = xarray.open_dataset(f0)
    # ds1 = xarray.open_dataset(f1)
    #
    # try:
    #     for v,_,_ in vars1D:
    #         plt.plot(ds0[v],ds1[v])
    #         plt.title(v)
    #         plt.show()
    #     ds0.close()
    #     ds1.close()
    # except:
    #     ds0.close()
    #     ds1.close()
    #
    # exit()

    # Call the main function
    station_nc = "X:/data/SIESTA_intercompare/output/run_11Jun2019_0844_NICE_GC_SWD_Conform_3p5_wc_walb_0p508_algmig2_both_directions/station_11_999-999_11Jun2019_0844.nc"
    cwd,_ = os.path.split(station_nc)
    output_filename = os.path.join(cwd,'EXP2_NICE_Saenz_2025-08-07.nc')
    write_bepsii_ncfile(station_nc, 113, output_filename)

    station_nc = "X:/data/SIESTA_intercompare/output/run_11Jun2019_0931_RL_GC_SWD_Conform_3p5_wc_walb_0p508_algmig2_ksi4/station_11_999-999_11Jun2019_0931.nc"
    cwd,_ = os.path.split(station_nc)
    output_filename = os.path.join(cwd,'EXP1_Resolute_Saenz_2025-08-07.nc')
    write_bepsii_ncfile(station_nc, 113, output_filename)

    station_nc = "X:/data/SIESTA_intercompare/output/run_11Jun2019_0940_NICE_GC_SWD_Conform_3p5_wc_walb_0p508_algmig2_ksi4/station_11_999-999_11Jun2019_0940.nc"
    cwd,_ = os.path.split(station_nc)
    output_filename = os.path.join(cwd,'EXP1_NICE_Saenz_2025-08-07.nc')
    write_bepsii_ncfile(station_nc, 113, output_filename)

    station_nc = "X:/data/SIESTA_intercompare/output/run_12Jun2019_2056_RL_GC_SWD_Conform_3p5_wc_walb_0p3_algmig2_la9/station_11_999-999_12Jun2019_2056.nc"
    cwd,_ = os.path.split(station_nc)
    output_filename = os.path.join(cwd,'EXP2_Resolute_Saenz_2025-08-07.nc')
    write_bepsii_ncfile(station_nc, 113, output_filename)
