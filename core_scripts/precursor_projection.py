from dask.distributed import Client, LocalCluster
import argparse
import xarray as xr
import os
import sys
import numpy as np
from domino.core import IndexGenerator
import itertools as it
from tqdm import tqdm
import subprocess
# import logging

###
#some debuggging code
import functools
import types
import sys

def log_function_call(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        print(f"\nCalling: {func.__name__} with args={args}, kwargs={kwargs}")
        return func(*args, **kwargs)
    return wrapper

def decorate_all_functions(module, decorator):
    for name in dir(module):
        attr = getattr(module, name)
        if isinstance(attr, types.FunctionType):
            setattr(module, name, decorator(attr))
###

def parse_args(arg_list=None):
    parser = argparse.ArgumentParser(description="Project model data onto precursor patterns.")

    # Required arguments
    parser.add_argument('--model', type=str, required=True,
                        help='name of model simulation to load as input.')
    parser.add_argument('--experiment', type=str, required=True,
                        help='name of simulation experiment, e.g. historical, ssp370')

    # Optional arguments
    parser.add_argument('--member', type=str, default='',
                        help='Which ensemble member, if any.')
    
    parser.add_argument('--variables', nargs='+', type=str, default=['z500_detrend','u850','v850'],
                    help='Which variables to project.')

    parser.add_argument('--overwrite',action='store_true',
                        help='Overwrite indices if they already exist.')
    
    parser.add_argument('--seasons', nargs='+', type=str, default=['DJF', 'MAM','JJA','SON'],
                    help='Which seasons to compute indices for.')
    
    parser.add_argument('--regions', nargs='+', type=int, default=None,
                help='Which regions to compute indices for. Default is to use all regions in maskfile.')
    
    parser.add_argument('--lags', nargs='+', type=int, default=[0],
                    help='Which lag times to use.')

    parser.add_argument('--inputdir',type=str,default='/Data/gfi/share/ModData/CMIP_EU_Precip_Precursors/raw/',
                        help='Directory in which to look for field data.')
    
    parser.add_argument('--auxdir',type=str,default='/Data/skd/projects/global/cmip6_precursors/aux/',
                    help='Directory in which to look for reference cycle, precursor patterns and precursor coefficients.')
    
    parser.add_argument('--savedir',type=str,default='/Data/skd/projects/global/cmip6_precursors/outputs/indices/',
                    help='Directory in which to save output.')
    
    parser.add_argument('--precursorid',type=str,default='standard',
                    help='Which precursor patterns to use.')
    
    parser.add_argument('--refid',type=str,default='ERA5',
                    help='Which reference dataset seasonal cycle and std params are taken from.')
    
    parser.add_argument('--debug',action='store_true',
                        help='Print all functions as they are called.')

    return parser.parse_args(arg_list)

def get_save_path(args):
    savedirs=[f'{args.savedir}{args.model}/{v.split("_detrend")[0]}/{args.experiment}/' for v in args.variables]
    for savedir in savedirs:
        os.makedirs(savedir,exist_ok=True)
        # only tun if dir was just made, otherwise we might be changing permissions of a dir that other processes are writing to.
        # subprocess.run(['chmod','-R','g+wrx',savedir])
    
    if args.member!='':
        savedirs=[dir+f'member_{args.member}_' for dir in savedirs]

    return savedirs

def load_input_field(args,v,v_in_file):
    indir=f'{args.inputdir}raw/{args.model}/{v}/{args.experiment}/'
    filenames=[file for file in os.listdir(indir) if file.endswith('.nc')]
    if args.member=='': #single model case
        if len(filenames)>1:
            raise(IOError(f'Expected single file in input directory {indir} when --member flag is absent. Found {len(filenames)}.'))

    else: #ensemble case
        filenames=[file for file in filenames if file.find(f'_{args.member}_')!=-1]
        if len(filenames)>1:
            raise(IOError(f'Expected to find only one file in input directory {indir} containing string _{args.member}_. Found {len(filenames)}.'))
    
    data=xr.open_dataset(indir+filenames[0],
        chunks=dict(time=-1,lat=-1,lon=-1))[v_in_file].chunk('auto')
    with open("./mylog.log", "a") as f:
        f.write("variable in the dataset right after the reading the input field for "+v_in_file+":\n")
        f.write(data.name+"\n")
    
    return data

def load_input_fields(args,name_dict):

<<<<<<< HEAD
    ds=xr.merge([load_input_field(args,v,name_dict[v]).rename(v) for v in args.variables], compat='override', join='outer')
=======
    # ds=xr.merge([load_input_field(args,v,name_dict[v]).rename(v.split("_detrend")[0]) for v in args.variables],compat='override')
    ds=xr.merge([load_input_field(args,v,name_dict[v]).rename(v) for v in args.variables],compat='override')
    with open("./mylog.log", "a") as f:
        f.write("all variable in the merged dataset right after the reading the input fields:\n")
        f.write(" ".join(ds.variables.keys())+"\n")
>>>>>>> d421045 (committing all changes to the CMIP precip precursor code made by Johannes Lutzmann up to 20/03/2026)
    return ds

def interp_to_1degree_grid_over_precursor_domain(ds,latmin=15,latmax=85,lonmin=-100,lonmax=60):
    """assumes dims are named lat and lon"""
    deg1_lats=np.arange(latmin,latmax+0.01)
    deg1_lons=np.arange(lonmin,lonmax+0.01)

    return ds.interp(lat=deg1_lats,lon=deg1_lons,method='linear') 



def apply_region_masking_and_average(mask,field,regions,lon='lon',lat='lat'):

    weights=np.cos(np.deg2rad(field[lat]))
    indices=[]

    for r in regions:
        M=((mask==r).astype(int).interp_like(field)>0)
        masked_field=field.where(M==1)
        collapsed_da=masked_field.mean(lon).weighted(weights).mean(lat)
        indices.append(collapsed_da.assign_coords(mask_id=r))

    return indices

def ensure_lon_180(da):
    try:
        lon=da.lon
    except:
        raise(ValueError('This script assumes the model data has a "lon" coordinate.'))
    
    #this does nothing to a da which already is in -180 to 180 degree longitude
    da=da.assign_coords(lon = (lon + 180) % 360 - 180)
    da=da.sortby("lon")
    return da

def load_precursor_patterns_and_params(args):
        
    pattern_path=f'{args.auxdir}/patterns/{args.precursorid}_'
    param_path=f'{args.auxdir}/std_params/{args.precursorid}_{args.refid}_'

    # with open("./mylog.log", "a") as f:
    #     f.write("variable names before the varname splitting:\n")
    #     f.write(" ".join(args.variables)+"\n")
    #if a variable ends _detrend, we load the pattern of the base variable.
    var_names=[v.lower().split('_detrend')[0] for v in args.variables]
<<<<<<< HEAD
=======
    # with open("./mylog.log", "a") as f:
    #     f.write("variable names after the varname splitting:\n")
    #     f.write(" ".join(var_names)+"\n")

>>>>>>> d421045 (committing all changes to the CMIP precip precursor code made by Johannes Lutzmann up to 20/03/2026)
    index_names=[f'{v}_lag{l}_index_val1' for v,l in it.product(var_names,args.lags)]
    pattern_all_seasons=[]
    param_all_seasons=[]

    for s in args.seasons:
        pattern_all_regions=[]
        param_all_regions=[]

        for r in args.regions:
            P=xr.open_dataset(pattern_path+f'{s}_region{r:02d}.nc')
            # with open("./mylog.log", "a") as f:
            #     f.write("variable names before the varname splitting:\n")
            #     f.write(" ".join(list(ds.variables.keys()))+"\n")
            
            
            P=P[var_names].sel(lag=args.lags,index_val=1).assign_coords(region_id=r,season=s)
            pattern_all_regions.append(P.load())

            sp=xr.open_dataset(param_path+f'{s}_region{r:02d}.nc')
            sp=sp[index_names].assign_coords(region_id=r,season=s)
            param_all_regions.append(sp.load())

        pattern_all_seasons.append(xr.concat(pattern_all_regions, 'region_id', coords='different', compat='equals', join='outer'))
        param_all_seasons.append(xr.concat(param_all_regions, 'region_id', coords='different', compat='equals', join='outer'))

    precursor_patterns=xr.concat(pattern_all_seasons, 'season', coords='different', compat='equals', join='outer')
    precursor_std_params=xr.concat(param_all_seasons, 'season', coords='different', compat='equals', join='outer')

    
    
    return precursor_patterns,precursor_std_params
    
def load_cycle(args):

    cycle_path=f'{args.auxdir}/cycles/{args.refid}.nc'

    #here we do treat detrended variables differently, so this is
    #different to the corresponding line in load_precursor_patterns_and_params
    #Renaming z500 to match variables in ERA5 cycle data.
    var_names=[v.lower().split('_detrend')[0] for v in args.variables]
    cycle=xr.open_dataset(cycle_path,chunks='auto')[var_names]
    return cycle.load()

def geopotential_height_to_geopotential_if_present(ds):

    for dv in ds.data_vars:
        da=ds[dv]
        if da.attrs['standard_name']=='geopotential_height' and da.attrs['units']=='m':
            da.attrs['standard_name']=='geopotential'
            da.attrs['units']=='m**2/s**2'
            ds[dv]=da*9.81
    return ds


def deseasonalise_field(targ_field,cycle):

    #we shift cycle from a time axis to a dayofyear axis
    cycle=cycle.assign_coords(dayofyear=cycle['time.dayofyear'])
    cycle=cycle.swap_dims(time='dayofyear')

    #Then resample the cycle's dayofyear climatology along the targ_field time axis
    expanded_cycle=cycle.sel(dayofyear=targ_field['time.dayofyear'].values)
    expanded_cycle=expanded_cycle.swap_dims(dayofyear='time').assign_coords(time=targ_field.time)
    
    #raises error if the cycle and targ field don't match up
    targ_field,expanded_cycle=xr.align(targ_field,expanded_cycle,join='exact')
    with open("./mylog.log", "a") as f:
        f.write("target field before deseasonalising in deseasonalise_field:\n")
        f.write(" ".join(list(targ_field.variables.keys()))+"\n")
        f.write("expanded_cycle before deseasonalising in deseasonalise_field:\n")
        f.write(" ".join(list(targ_field.variables.keys()))+"\n")

    #return deseasonalised field
    deseas_field=targ_field-expanded_cycle
    with open("./mylog.log", "a") as f:
        f.write("target field after deseasonalising in deseasonalise_field:\n")
        f.write(" ".join(list(deseas_field.variables.keys()))+"\n")
    return deseas_field

def drop_vars_for_existing_files_if_any(P,field,save_paths):

    files_exist=np.array([os.path.isfile(savepath) for savepath in save_paths])
    files_that_exist=np.array(save_paths)[files_exist]
    files_to_write=np.array(save_paths)[~files_exist]

    n_exist=sum(files_exist)
    if n_exist>0:
        print(f'The following files already exist, and --overwrite not set. Skipping.')
        print(files_that_exist)

        vars_to_compute=[p.split('/')[-3] for p in files_to_write]
        vars_to_compute=[v.lower().split('_detrend')[0] for v in vars_to_compute]

        P=P[vars_to_compute]
        field=field[vars_to_compute]

    return P,field,files_to_write

# Safely serialize args to be valid NetCDF attributes
def make_serializable_attrs(args):
    attr_dict = {}
    for key, value in vars(args).items():

        if isinstance(value, bool):
            attr_dict[key]=str(bool(value))

        elif isinstance(value, (str, int, float)):
            attr_dict[key] = value
        elif value is None:
            attr_dict[key] = "None"
        elif isinstance(value, (list, tuple)):
            attr_dict[key] = ",".join(str(v) for v in value)
        else:
            attr_dict[key] = str(value)  # fallback: best-effort string
    return attr_dict


def strip_char_dim_name_encoding(da):
    da = da.copy(deep=False)
    da.encoding = {key: value for key, value in da.encoding.items() if key != 'char_dim_name'}

    for coord_name in da.coords:
        coord = da.coords[coord_name]
        coord.encoding = {
            key: value for key, value in coord.encoding.items() if key != 'char_dim_name'
        }

    return da


def project_onto_precursor_indices_and_save(ds,patterns,params,args):
    
    #work out the full save directory
    outdirs=get_save_path(args)
    
    
    slices=[dict(lag=l,index_val=1) for l in args.lags]
    for r in tqdm(args.regions):
        for s in args.seasons:

            in_season=ds['time.season']==s
            field=ds.isel(time=in_season)

            P=patterns.sel(season=s,region_id=r)
            param=params.sel(season=s,region_id=r)
            
            save_paths=[dir+f'{s}_region{r}.nc' for dir in outdirs]
            # with open("./mylog.log", "a") as f:
            #     f.write("outputting the savepaths right after defining them:\n")
            #     f.write(" ".join(save_paths)+"\n")
            #     f.write("outputting the variables before the drop_vars and expand_dims operations:\n")
            #     f.write(" ".join(list(field.variables.keys()))+"\n")
                
            if not args.overwrite:
                P,field,save_paths=drop_vars_for_existing_files_if_any(P,field,save_paths)
                if len(P.data_vars)==0:
                    continue #skip to next region-season combination if all files already written.
            try:
                P=P.expand_dims('index_val').assign_coords(index_val=[1])
                P=P.expand_dims('lag').assign_coords(lag=args.lags)
            except Exception as e:
                pass
            
<<<<<<< HEAD
=======
            with open("./mylog.log", "a") as f:
                f.write("outputting the variables after the drop_vars and expand_dims operations:\n")
                f.write(" ".join(list(field.variables.keys()))+"\n")
                f.write("outputting the variables in the Patterns and Parameters at the same time:\n")
                f.write(" ".join(list(P.variables.keys()))+"\n")
                f.write(" ".join(list(param.variables.keys()))+"\n")
                
            if "z500_detrend" in list(field.variables.keys()):
                field = field.rename({"z500_detrend": "z500"})

>>>>>>> d421045 (committing all changes to the CMIP precip precursor code made by Johannes Lutzmann up to 20/03/2026)
            #generate the indices and standardise them
            IG=IndexGenerator()
            indices=IG.generate(P,field,slices=slices,
                ix_means=param.sel(param='mean'),ix_stds=param.sel(param='std'))
<<<<<<< HEAD
            #adds all metadata from index generation into each index
            #save
=======
            with open("./mylog.log", "a") as f:
                f.write("These are the variables that are left after IG.gnerate:\n")
                f.write(" ".join(list(indices.data_vars))+"\n")
            # adds all metadata from index generation into each index
            # save
            with open("./mylog.log", "a") as f:
                f.write("outputting all savepath and dv values:\n")
>>>>>>> d421045 (committing all changes to the CMIP precip precursor code made by Johannes Lutzmann up to 20/03/2026)
            for savepath,dv in zip(save_paths,list(indices.data_vars)):
                if 'z500_detrend' in savepath:
                    savepath = savepath.replace("z500_detrend", "z500")
                da=indices[dv]
                da.attrs = make_serializable_attrs(args)
<<<<<<< HEAD
                da = strip_char_dim_name_encoding(da)
=======
                sys.stdout.flush()
>>>>>>> d421045 (committing all changes to the CMIP precip precursor code made by Johannes Lutzmann up to 20/03/2026)
                da.to_netcdf(savepath)
                with open("./mylog.log", "a") as f:
                    f.write(savepath+" "+dv+"\n")
    
    return# save_paths, indices.data_vars

def run_precursor_projection(
    model, experiment, member="", variables=['z500_detrend','u850','v850'],
    overwrite=False, seasons=['DJF', 'MAM','JJA','SON'], regions=None,
    lags=[0], inputdir='/Data/gfi/share/ModData/CMIP_EU_Precip_Precursors/',
    auxdir='/Data/gfi/share/ModData/CMIP_EU_Precip_Precursors/aux/',
    savedir='/Data/skd/projects/global/cmip6_precursors/outputs/indices/', 
    precursorid='standard', refid='ERA5', debug=False):
    
    arg_list = (
        ["--model", model, "--experiment", experiment, "--member", member,
         "--variables"]+variables+["--seasons"]+seasons
        +["--lags"]+[str(l) for l in lags]
        +["--inputdir", inputdir, "--auxdir", auxdir, "--savedir", savedir,
          "--precursorid", precursorid, "--refid", refid])

    if overwrite:
        arg_list.append("--overwrite")
        
    if debug:
        arg_list.append("--debug")

    if not regions is None:
        arg_list+=["--regions"]+[str(r) for r in regions]
    
    print(arg_list)
    args = parse_args(arg_list)
    
    main(args)
    return
        

def main(args):
    #use multi-core for speed
    cluster = LocalCluster(n_workers=4, memory_limit='16GiB')
    client = Client(cluster)
    print('Access dask dashboard: ', client.dashboard_link)
    with open("./mylog.log", "w") as f:
        f.write("Starting a new log:\n")
    
    if args.debug:
        decorate_all_functions(sys.modules[__name__], log_function_call)

    #use a pre-established set of regions if none specified
    if args.regions is None:
        # region 2 was an uninhabited island so we dropped it.
        args.regions=[1,*np.arange(3,40)] 
    else:
        pass

    #load precursor patterns and standardisation params
    # which ensure precursor index has mean 0 and std 1
    patterns,params=load_precursor_patterns_and_params(args)

    #load reference seasonal cycle
    cycle=load_cycle(args)

    #load model field data and put it on a 1deg grid over the precursor domain
    targ_field=load_input_fields(args,var_name_dict)
    # Renaming z500_detrend into z500 to match cycle variable name
    for data_vars in targ_field.data_vars:
        if data_vars == 'z500_detrend':
            targ_field = targ_field.rename({data_vars: 'z500'})

    targ_field=ensure_lon_180(targ_field)
    targ_field=interp_to_1degree_grid_over_precursor_domain(targ_field)

    #Now we compute, everything was lazy up to now
    targ_field=targ_field.load()
    with open("./mylog.log", "a") as f:
        f.write("target field after loading:\n")
        f.write(" ".join(list(targ_field.variables.keys()))+"\n")
    
    targ_field=geopotential_height_to_geopotential_if_present(targ_field)
    with open("./mylog.log", "a") as f:
        f.write("target field after geopotential_height_to_geopotential_if_present:\n")
        f.write(" ".join(list(targ_field.variables.keys()))+"\n")
    #remove seasonal cycle from field data:
    targ_field=deseasonalise_field(targ_field,cycle)
<<<<<<< HEAD
    if 'z500_detrend' in targ_field:
        targ_field = targ_field.rename(z500_detrend='z500')
=======
    
    with open("./mylog.log", "a") as f:
        f.write("target field just before the last function call:\n")
        f.write(" ".join(list(targ_field.variables.keys()))+"\n")
        # f.write("patterns before the last function call:\n")
        # f.write(" ".join(list(patterns.variables.keys()))+"\n")
        # f.write("params field just before the last function call:\n")
        # f.write(" ".join(list(params.variables.keys()))+"\n")
        
>>>>>>> d421045 (committing all changes to the CMIP precip precursor code made by Johannes Lutzmann up to 20/03/2026)
    #Do the projection and save them to file
    # save_paths, data_vars = 
    
    project_onto_precursor_indices_and_save(targ_field,patterns,params,args)
    return #save_paths, data_vars

var_name_dict={
    'z500_detrend':'zg_detrend',
    'z500':'zg',
    'u850':'ua',
    'v850':'va'
}
if __name__=='__main__':

    # #use multi-core for speed
    # cluster = LocalCluster(n_workers=4, memory_limit='16GiB')
    # client = Client(cluster)
    # print('Access dask dashboard: ', client.dashboard_link)
    
    args = parse_args()
    main(args)
    # if args.debug:
    #     decorate_all_functions(sys.modules[__name__], log_function_call)

    # #use a pre-established set of regions if none specified
    # if args.regions is None:
    #     # region 2 was an uninhabited island so we dropped it.
    #     args.regions=[1,*np.arange(3,40)] 
    # else:
    #     pass

    # #load precursor patterns and standardisation params
    # # which ensure precursor index has mean 0 and std 1
    # patterns,params=load_precursor_patterns_and_params(args)

    # #load reference seasonal cycle
    # cycle=load_cycle(args)

    # #load model field data and put it on a 1deg grid over the precursor domain
    # targ_field=load_input_fields(args,var_name_dict)
    # targ_field=ensure_lon_180(targ_field)
    # targ_field=interp_to_1degree_grid_over_precursor_domain(targ_field)

    # #Now we compute, everything was lazy up to now
    # targ_field=targ_field.load()

    # targ_field=geopotential_height_to_geopotential_if_present(targ_field)

    # #remove seasonal cycle from field data:
    # targ_field=deseasonalise_field(targ_field,cycle)

    # #Do the projection and save them to file
    # project_onto_precursor_indices_and_save(targ_field,patterns,params,args)
