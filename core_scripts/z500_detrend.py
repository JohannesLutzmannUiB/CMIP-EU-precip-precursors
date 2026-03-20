import argparse
import xarray as xr
import numpy as np
import os
import subprocess
import glob

g = 9.80665

def parse_args():
    parser = argparse.ArgumentParser(description = "detrending of z500, removal of the stationary seasonal cycle based on ERA5")
    
    parser.add_argument("--model", required = True, 
                        help = 'name of model simulation to load as input, e.g. CNRM-CM6-1')
    parser.add_argument("--experiment", required=True, 
                        help = 'name of simulation experiment, e.g. historical, ssp370')
    parser.add_argument("--member", required=True, 
                        help = 'Which ensemble member, if any')

# Optional arguments
    parser.add_argument("--basedir", default="/Data/gfi/share/ModData/CMIP_EU_Precip_Precursors/raw", 
                        help = 'root directory')
    parser.add_argument("--varname", default="zg", 
                        help = 'name of the variable to detrend')
    parser.add_argument("--era5_cycle", default="/Data/gfi/share/ModData/CMIP_EU_Precip_Precursors/aux/cycles/ERA5.nc", 
                        help = 'file with the ERA5 seasonal cycle')
    parser.add_argument('--latmin', type=float, default=30, 
                        help='min lat')
    parser.add_argument('--latmax', type=float, default=85, 
                        help='max lat')

    return parser.parse_args()


# def main(args):
def main(
    model, experiment, member,
    basedir="/Data/gfi/share/ModData/CMIP_EU_Precip_Precursors/raw",
    varname="zg",
    era5_cycle=("/Data/gfi/share/ModData/CMIP_EU_Precip_Precursors/"
                +"aux/cycles/ERA5.nc"), latmin=30, latmax=85):
    
    input = f"{basedir}/{model}/z500/{experiment}/"
    # files=[file for file in os.listdir(input) if file.endswith('.nc')]
    files=[file for file in os.listdir(input) if member in file]
    if len(files) == 1:
        file = files[0]
    elif len(files) > 1:
        raise ValueError("Ambiguous Input: more than one matching file path detected")
    else:
        raise ValueError("No matching file found")
    
    input_pattern = f'{input}/{file}'
    date_range = file.split('_')[6].split(".")[0]
    print(date_range)
    

    # if not input_pattern :
    #     raise FileNotFoundError(f"No files found in : {input/input_pattern}")

    ds = xr.open_dataset(input_pattern)
    x = ds[varname]
    x.load()
    print(date_range[0:8], date_range[9:17])

    x = x.sel(time = slice(date_range[0:8], date_range[9:17]))

    x_detrended = detrend_seasonal_cycle(x, era5_cycle, latmin, latmax)
    if varname.split :
        x_detrended = x_detrended.rename(f'{varname}_detrend')
    else :
        x_detrended = x_detrended.rename(f'{varname}')

    output_pattern = f"{basedir}/{model}/z500_detrend/{experiment}/z500_detrend_day_{model}_{experiment}_{member}_gn_{date_range}.nc"
    print(output_pattern)
    os.makedirs(os.path.dirname(output_pattern), exist_ok=True)
    # subprocess.run(['chmod','-R','g+wrx',os.path.dirname(output_pattern)])   # add this only when dir above is created
    ds_out = x_detrended.to_dataset().assign_attrs(ds.attrs)
    ds_out.to_netcdf(output_pattern)




def detrend_seasonal_cycle(x, era5_cycle, latmin, latmax):
    # Subset latitudes
    x = x.sel(lat = slice(latmin, latmax), ) #add longitude extend here

    # Weighted mean lat/lon
    weights = np.cos(np.deg2rad(x.lat))
    xmean = x.weighted(weights).mean(('lat', 'lon'))

    # Group by (years, months)
    xmean = xmean.assign_coords(year = xmean['time.year'], month = xmean['time.month']).set_index(time = ('year', 'month'))
    monthly_x = xmean.groupby('time').mean().unstack('time')

    if 'plev' in monthly_x.coords:
        monthly_x = monthly_x.drop_vars('plev')

    # Gaussian smoothing
    y = monthly_x.to_dataset('month').to_dataframe()
    smooth_y = y.rolling(window = 31, min_periods = 1, win_type = 'gaussian', center = True).mean(std = 17)
    smooth_y = smooth_y.to_xarray().to_array('month')

    # Removing the smoothed trend
    stacked_x = x.assign_coords(year = x['time.year'], month = x['time.month']).set_index(time = ('year', 'month'))
    z = (smooth_y-smooth_y.isel(year=0)).stack(time = ('year', 'month')) # This works only if data goes back to early SSP
    x_detrended = stacked_x - z.sel(time = stacked_x.time)
    

    x_detrended = x_detrended.reset_index('time')  # supprimer le MultiIndex
    x_detrended = x_detrended.assign_coords(time = x.time)

    # Removing the stationary seasonal cycle and mean state (ERA5)
    # era_cycle = xr.open_dataset(era5_cycle).z500
    # era_cycle = era_cycle.sel(lat = slice(latmin, latmax + 1))
    # NH_mean_cycle = era_cycle.groupby('time.month').mean().weighted(np.cos(np.deg2rad(era_cycle.lat))).mean(('lat', 'lon'))

    # # Conversion of the unit (m²/s² → m)
    # NH_mean_cycle = NH_mean_cycle/g

    # # Substract the seasonal cycle
    # x_detrended = x_detrended.groupby('time.month') + NH_mean_cycle
    x_detrended = x_detrended.assign_attrs(x.attrs)    

    return x_detrended




if __name__ == "__main__":
    
    args = parse_args()
    main(args.model, args.experiment, args.member, basedir=args.basedir,
         varname=args.varname, era5_cycle=args.era5_cycle, latmin=args.latmin,
         latmax=args.latmax)
