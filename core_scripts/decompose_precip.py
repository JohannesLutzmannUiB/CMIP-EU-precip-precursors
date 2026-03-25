import os
import numpy as np
import xarray as xr
import pandas as pd
from sklearn.decomposition import PCA
import argparse


def parse_args(arg_list=None):
    parser = argparse.ArgumentParser(
        description="Use precursor indices and precipitation indices to decompose model bias and trends."
    )

    # Required arguments
    parser.add_argument('--model', type=str, required=True,
                        help='Name of model to decompose.')

    parser.add_argument('--experiment', type=str, required=True,
                        help='Name of model experiment.')

    parser.add_argument('--time_period', type=int, nargs=2, required=True,
                        help='The time period to use, given as 2 years (inclusive), e.g. 1979 2015.')


    parser.add_argument('--historical_members', nargs='+', type=str, default=None,
                        help='Which ensemble member or members to use for the historical period. Defaults to assuming not an ensemble.')
    
    parser.add_argument('--future_members', nargs='+', type=str, default=None,
                        help='Which ensemble member or members to use for the future period. Defaults to assuming not an ensemble.')
    
    parser.add_argument('--auxdir', type=str, required=True,
                        help='Directory to store auxiliary artifacts like PCA patterns and bin edges.')

    # Data configuration
    parser.add_argument('--variables', nargs='+', type=str, default=['z500','u850','v850'],
                        help='List of precursor variable names.')
    
    parser.add_argument('--hazardvariable', type=str, default='pr',
                        help='Name of the hazard variable to decompose.')

    # Ensemble configuration
    parser.add_argument('--members', nargs='+', type=str, default='r1i1p1f1',
                        help='Which ensemble member or members to use. Defaults to r1i1p1f1. If multiple members are passed, output folder will be specified by "--ensname"')

    parser.add_argument('--ensname',type=str,default='ens',help='specifies name of member directory where the ensemble will be stored.')

    # Reference / precomputed parameters
    parser.add_argument('--ref', action='store_true',
                        help='If set, treat current run as reference (no experiment or member suffix).')
    parser.add_argument('--ref_pca', type=str, default=None,
                        help='Path to precomputed PCA patterns (optional).')
    parser.add_argument('--ref_bins', type=str, default=None,
                        help='Path to precomputed bin edges (optional).')
    parser.add_argument('--ref_thresh', type=str, default=None,
                        help='Path to precomputed categorical threshold (optional).')

    # Decomposition parameters
    parser.add_argument('--nprecursorbins', type=int, default=10,
                        help='Number of bins for precursor decomposition.')
    parser.add_argument('--eventthreshold', type=float, default=0.95,
                        help='Threshold for converting hazard to categorical event.')

    # Execution control
    parser.add_argument('--overwrite', action='store_true',
                        help='Overwrite existing output files.')

    # Filtering
    parser.add_argument('--regions', nargs='+', type=int, default=None,
                        help='Which regions to compute metrics for. Default is to use all regions.')
    parser.add_argument('--seasons', nargs='+', type=str, default=['DJF', 'MAM', 'JJA', 'SON'],
                        help='Seasons to process.')

    return parser.parse_args(arg_list)

def get_ref_data(args):
    bp=f'{args.inputdir}{args.referencemodel}/'
    # vars=np.atleast_1d(args.variables).append(args.hazardvariable)
    vars=list(np.atleast_1d(args.variables)) + [args.hazardvariable]
    ds=[]
    for v in vars:
        dir=f'{bp}{v}/{args.historical_experiment}/'
        a1=[]
        for r in args.regions:
            a2=[]
            for s in args.seasons:
                # da=xr.open_dataarray(f'{dir}+{args.season}_region{args.region_id}.nc')
                da=xr.open_dataarray(f'{dir}/{s}_region{r}.nc')
                a2.append(da.assign_coords(season=s))
            a1.append(xr.concat(
                a2, 'season', coords='different', compat='equals',
                join='outer').assign_coords(region_id=r))
        da=xr.concat(a1, 'region_id', coords='different', compat='equals', join='outer')
        # ds.append(da.rename(v))
        ds.append(da)
    ds=xr.merge(ds, compat='override', join='outer')

def get_decomp_output_path(args, suff='csv'):
    model=args.model

    if args.ref:
        path= f'{args.savedir}/decompositions/{model}/decomp.{suff}'
        return path
    
    exp=args.experiment
    members=np.atleast_1d(args.members)
    if len(members)>1:
        mem_lab=args.ensname
    else:
        mem_lab=members[0]
        
    path = f'{args.savedir}/decompositions/{model}/{exp}/{mem_lab}/decomp.{suff}'

    return path

def get_index_and_hazard_data(args):

    bp=f'{args.inputdir}{args.model}/'
    vars=list(np.atleast_1d(args.variables)) + [args.hazardvariable]
    ds=[]
    mems=np.atleast_1d(args.historical_members)
    not_ens= (len(mems)==1) and (mems[0] is None)
    
    for v in vars:
        vdir=f'{bp}{v}/{args.experiment}/'
        os.makedirs(vdir,exist_ok=True)
        a0=[]
        for mem in mems:
            if args.ref:
                dir = vdir
            else:
                dir = vdir + f'member_{mem}_'
            os.makedirs(dir,exist_ok=True)
            a1=[]
            for r in args.regions:
                a2=[]
                for s in args.seasons:
                    da=xr.open_dataarray(f'{dir}{s}_region{r}.nc')
                    a2.append(da.assign_coords(season=s))
                a1.append(xr.concat(a2, 'season', coords='different', compat='equals', join='outer').assign_coords(region_id=r))

            a1=xr.concat(a1, 'region_id', coords='different', compat='equals', join='outer')
            if not args.ref: a1=a1.assign_coords(member=mem)
            a0.append(a1)
        a0=xr.concat(a0, 'member', coords='different', compat='equals', join='outer')#.rename(v)
        ds.append(a0)

    ds=xr.merge(ds, compat='override', join='outer')
    y0,y1=args.time_period
    ds = ds.sel(time=slice(str(y0), str(y1)))
    return ds

def get_future_data(args):
    #Supports bias only analysis
    if args.future_experiment == "none":
        return None
        
    bp=f'{args.inputdir}{args.model}/'

    # vars=np.atleast_1d(args.variables).append(args.hazardvariable)
    vars=list(np.atleast_1d(args.variables)) + [args.hazardvariable]

    ds=[]
    mems=np.atleast_1d(args.future_members)
    not_ens= (len(mems)==1) and (mems[0] is None)
    
    for v in vars:
        vdir=f'{bp}{v}/{args.future_experiment}/'
        os.makedirs(vdir,exist_ok=True)
        # subprocess.run(["chmod", "-R", "g+rwx", vdir], check=True)
        a0=[]
        for mem in mems:
            if not_ens:
                dir=vdir
            else:
                dir=vdir+f'member_{mem}_'

            a1=[]
            for r in args.regions:
                a2=[]
                for s in args.seasons:
                    # da=xr.open_dataarray(f'{dir}_{args.season}_region{args.region_id}.nc')
                    da=xr.open_dataarray(f'{dir}{s}_region{r}.nc')
                    a2.append(da.assign_coords(season=s))
                a1.append(xr.concat(a2, 'season', coords='different', compat='equals', join='outer').assign_coords(region_id=r))

            a1=xr.concat(a1, 'region_id', coords='different', compat='equals', join='outer')
            if not not_ens: a1=a1.assign_coords(member=mem)
            a0.append(a1)

        a0=xr.concat(a0, 'member', coords='different', compat='equals', join='outer')#.rename(v)
        ds.append(a0)

    ds=xr.merge(ds, compat='override', join='outer')

    y0,y1=args.future_period
    # future_period=np.arange(y0,y1+1)
    # ds=ds.isel(np.isin(ds['time.year'],future_period))
    ds = ds.sel(time=slice(str(y0), str(y1)))
    return ds

def get_savepaths(args,s,r,suff='csv'):
    s1=f'{args.savedir}/{args.model}/'

    # check if the data is part of an ensemmble
    mems_hist = np.atleast_1d(args.historical_members)
    not_ens_hist = (len(mems_hist)==1) and (mems_hist[0] is None)

    mems_fut = np.atleast_1d(args.future_members)
    not_ens_fut = (len(mems_fut)==1) and (mems_fut[0] is None)

    not_ens = not_ens_hist or not_ens_fut    
    # Create a subdirectory for this specific cobination
    # of historical and future members
    if not not_ens:
        if len(mems_hist) > 3:
            mems_str_hist = ("historical_"+mems_hist[0]+"_through_"
                             +mems_hist[-1])
        else:
            mems_str_hist = "historical_"+"_".join(mems_hist)

        if len(mems_fut) > 3:
            mems_str_fut = (args.future_experiment+"_"+mems_fut[0]+"_through_"
                            +mems_fut[-1])
        else:
            mems_str_fut = (args.future_experiment+"_"+"_".join(mems_fut))

        
        s1 += mems_str_hist + "_" +  mems_str_fut + "/"
    
    s2=f'{s}_region{r}.{suff}'
    return s1+'decomp_'+s2, s1+'terms_'+s2



def run_decompose_precip(
    model, future_experiment, eventthreshold=0.95, nprecursorbins=10,
    seasons=['DJF', 'MAM','JJA','SON'], regions=None, historical_members=None,
    future_members=None, referencemodel='ERA5',
    historical_experiment='historical', hist_period=[1979,2014],
    future_period=[2060,2100], variables=['z500','u850','v850'],
    hazardvariable='pr', overwrite=False,
    inputdir='/Data/skd/projects/global/cmip6_precursors/outputs/indices/',
    auxdir='/Data/skd/projects/global/cmip6_precursors/aux/',
    savedir='/Data/skd/projects/global/cmip6_precursors/outputs/decompositions/'
):
    
    arg_list = (
        ["--model", model, "--future_experiment", future_experiment,
        '--eventthreshold', str(eventthreshold), '--nprecursorbins',
         str(nprecursorbins), "--seasons"]+seasons
        +['--referencemodel', referencemodel,
          '--historical_experiment', historical_experiment, "--variables"]
        +variables+['--hazardvariable', hazardvariable, '--inputdir', inputdir,
                    '--auxdir', auxdir, '--savedir', savedir])
    
    if overwrite:
        arg_list.append("--overwrite")

    if not regions is None:
        arg_list+=["--regions"]+[str(r) for r in regions]

    if not future_members is None:
        if type(future_members) in (list, tuple):
            arg_list+=["--future_members"]+[str(f) for f in future_members]
        elif type(future_members) is str:
            arg_list+=["--future_members", future_members]
    
    if not historical_members is None:
        if type(historical_members) in (list, tuple):
            arg_list+=(["--historical_members"]
                       +[str(h) for h in historical_members])
        elif type(historical_members) is str:
            arg_list+=["--historical_members", historical_members]
            

    if type(hist_period) is list:
        if len(hist_period)==2:
            arg_list+=["--hist_period"]+[str(h) for h in hist_period]

    if type(future_period) is list:
        if len(future_period)==2:
            arg_list+=["--future_period"]+[str(f) for f in future_period]

    args = parse_args(arg_list)
    main(args)
    return args
            


def decompose_hazard_odds_ratio(ref_ds,h_model_ds,
                                f_model_ds,h_var,s_var,bin_num=5,
                                p_dvs=None,make_h_var_cat=False,quantile=None):
    return _prep_and_decompose(binned_decomposition,return_decomp_as_dataarray,
                               ref_ds,h_model_ds,f_model_ds,h_var,s_var,bin_num,p_dvs,make_h_var_cat,quantile)

def _prep_and_decompose(decomp_func,output_func,ref_ds,
                        h_model_ds,f_model_ds,h_var,s_var,bin_num,p_dvs,make_h_var_cat,quantile):

    #needed if returned with full_output, but not used.
    ref_PCA_Solver=None
    quantile_thresh=None
    PCs=None
    EOFs=None
    #handle PCA projection (based on ref_ds variability) if needed.
    if s_var[:2]=='PC':

        if p_dvs is None:
            print('Using default p_dvs: Z500, U850 and V850')
            p_dvs=['z500_lag0_index_val1',
                   'u850_lag0_index_val1',
                   'v850_lag0_index_val1']

        pcN=int(s_var[2:])
        ref_ds,ref_PCA_Solver,PCs,EOFs=fit_PCA_to_ds(ref_ds,pcN,p_dvs)
        h_model_ds[s_var]=xr.DataArray(
            data=apply_PCA_to_ds(h_model_ds,ref_PCA_Solver,pcN,p_dvs),
            coords=dict(time=h_model_ds.time))
    
    f_model_ds[s_var]=xr.DataArray(
        data=apply_PCA_to_ds(f_model_ds,ref_PCA_Solver,pcN,p_dvs),
        coords=dict(time=f_model_ds.time))
    #handle categorical event def if needed.
    if make_h_var_cat:
        if quantile is None:
            raise(ValueError('Need quantile for exceedance calculation if make_h_var_cat is True'))
        
        quantile_thresh=ref_ds[h_var].quantile(quantile)

        ref_ds[h_var+'_cat']=(ref_ds[h_var]>quantile_thresh).astype(int)
        h_model_ds[h_var+'_cat']=(h_model_ds[h_var]>quantile_thresh).astype(int)
        f_model_ds[h_var+'_cat']=(f_model_ds[h_var]>quantile_thresh).astype(int)

        h_var=h_var+'_cat'
    
    #make sure we have binary event time series, or else I expect the computation will be wrong. If we change our mind we can delete this block.
    
    unique_h_vals=[np.unique(ds[h_var]) for ds in [ref_ds,h_model_ds,f_model_ds]]
    try:
        ds=ds.rename(tp='pr')
    except:
        ds.pr #fails if pr does not exist on ds.
    ds=ds.dropna('time',how='all')
    return ds

def compute_pca_pattern(ds,condition_var,p_dvs):

    if not condition_var.startswith('PC'):
        return None
    
    N = int(condition_var[2:])
    solver = PCA(n_components=N, whiten=True).fit(
        ds[p_dvs].squeeze().to_array('feature').T)
    PCA_pattern=solver.components_
    return PCA_pattern # a 2d numpy array

def handle_PCA_projection(ds,condition_var,p_dvs,PCA_pattern):

    if not condition_var.startswith('PC'):
        return ds
    
    X=ds[p_dvs].to_array('feature').T
    N = int(condition_var[2:])

    PC=(X.values @ PCA_pattern.T)[:,N-1][:,0]
    ds[condition_var]=xr.DataArray(data=PC,coords=dict(time=ds.time))
    return ds

def compute_bin_edges(ds, condition_var, bin_num):
    #our bins are based on the reference ds. We've modified this to include two unobserved bins, at the high and low end. 
    bins=np.array([-1000,*ds[condition_var].quantile(np.linspace(0,1,bin_num+1)),1000])
    bins[1]=bins[1]-0.01 # so min ref value is not included in unobserved first bin
    return bins #a list of floats

def compute_ref_thresh(ds, hazard_var, quantile):
    quantile_thresh=ds[hazard_var].quantile(quantile)
    return quantile_thresh #just a float

def handle_categorical_target(ds, h_var, is_cat, thresh):
    if not is_cat:
        return ds, h_var
    
    ds[h_var+'_cat']=(ds[h_var]>thresh).astype(int)
    h_var=h_var+'_cat'

    return ds, h_var

def apply_decomposition(ds, h_var, s_var, bins):
    Ph_s=ds.groupby_bins(s_var,bins=bins).mean()[h_var].fillna(0) #average value of hazard in bin. Is a probability for binary data. Bins with no hazard risk get a 0
    P_s=ds.groupby_bins(s_var,bins=bins).count()[s_var].fillna(0)/ds[s_var].time.size #occurence prob of synoptic bin. Bins that don't occur get a 0.

    return Ph_s, P_s

def decomp_to_df(Ph_s,P_s,model,season,region_id):
    
    data=dict(
        dyn= (P_s,'dynamical'),
        conv= (Ph_s,'conversion')
    )
    decomp_rows=[]
    for name, (vals, s1) in data.items():
        for b, v in enumerate(vals, 1):
            decomp_rows.append({"model": model, "season": season,
                        "region_id":region_id, "bin": b,
                        "source": s1, "value": np.atleast_1d(v.values)[0]})
                        
    return pd.DataFrame(decomp_rows)

def postprocess_decomp_dfs(df):
    return pd.concat(df,ignore_index=True)

def save_decomp_df(df,output_path):
    os.makedirs('/'.join(output_path.split('/')[:-1]),exist_ok=True)
    print('saving:')
    print(output_path)
    df.to_csv(output_path)
    return

def save_ensemble_documentation(output_path, members):
    if len(np.atleast_1d(members))<=1:
        return
    dirs=output_path.split('/')
    filename=dirs.pop(-1)

    output_dir='/'.join(dirs)
    doc_path=output_dir+f'members_included_in_{filename}'
    
    df=pd.DataFrame(data=list(members),columns=['member'])
    df.to_csv(doc_path)
    return

def save_pca(pca_dict,auxdir,model, experiment, p_dvs):
    filename=f'{model}_{experiment}_'
    for v in p_dvs:
        filename+=f'{v}_'
    filename+='PCA_patterns.npz'
    np.savez(auxdir+filename,**pca_dict)
    return
"""
def load_pca_pattern(auxdir,model,experiment,p_dvs):
    filename=f'{model}_{experiment}_'
    for v in p_dvs:
        filename+=f'{v}_'
    filename+='PCA_patterns.npz'
    return dict(np.load(auxdir+filename))

def load_bin_edges(auxdir, model, experiment, condition_var, binnum):
    filename=f'{model}_{experiment}_{condition_var}_n{binnum}_bin_edges.npz'
    return dict(np.load(auxdir+filename))

def load_cat_thresh(auxdir, model, experiment,hazardvar,eventthreshold):
    filename=f'{model}_{experiment}_{hazardvar}_p{eventthreshold}_thresholds.npz'
    return dict(np.load(auxdir+filename))

"""
def save_bin_edges(bin_dict, auxdir, model, experiment, condition_var, binnum):
    filename=f'{model}_{experiment}_{condition_var}_n{binnum}_bin_edges.npz'
    np.savez(auxdir+filename,**bin_dict)
    return

def load_param_data(path):
    return dict(np.load(path))

def save_cat_thresh(thresh_dict, auxdir, model, experiment,hazardvar,eventthreshold):

    filename=f'{model}_{experiment}_{hazardvar}_p{eventthreshold}_thresholds.npz'
    np.savez(auxdir+filename,**thresh_dict)
    return


def validate_sel_data(ds,s,r):

    
    nan_vals=ds.isnull().sum().to_array('feature')
    for dv,c in zip(nan_vals.feature.values,nan_vals.values):
        if c!=0: 
            N=ds[dv].size
            print(f'Warning: {c} nan values out of {N} for variable {dv} in data {s} region {r}.')
    return 0

def main(args):

    decomp_output_path=get_decomp_output_path(args)
    if os.path.isfile(decomp_output_path) and not args.overwrite:
        print(f'Decomp. file {decomp_output_path} exists, and --overwrite flag not set. Exiting.')
        exit()


    condition_var = args.variables[0] if len(args.variables) == 1 else 'PC1'
    p_dvs = [v.split('_detrend')[0] + '_lag0_index_val1' for v in args.variables]
    make_h_var_cat = True

    if args.regions is None:
        #region 2 was an uninhabited island so we dropped it.
        args.regions=[1,*np.arange(3,40)] 
    else:
        pass

    ref_data=get_ref_data(args)
    hist_data=get_hist_data(args)
    future_data=get_future_data(args)
    hazard_var=args.hazardvariable
    make_h_var_cat=True
    p=args.eventthreshold
    bin_num=args.nprecursorbins
    model=args.model
    all_summed_terms=[]
    for s in args.seasons[:]:
        for r in args.regions[:]:
            try:
                decomp_path,term_path=get_savepaths(args,s,r)
                rd=ref_data.sel(season=s,region_id=r).rename(tp='pr')
                hd=hist_data.sel(season=s,region_id=r)
                if future_data is None:
                    fd=hd
                else:
                    fd=future_data.sel(season=s,region_id=r)
                hd = xr.concat([hd.sel(member=k, drop=True) for k in hd.member],
                               dim='time', coords='different', compat='equals',
                               join='outer')
                fd = xr.concat([fd.sel(member=k, drop=True) for k in fd.member],
                               dim='time', coords='different', compat='equals',
                               join='outer')
                decomposed_hazard=decompose_hazard_odds_ratio(rd.dropna('time'),
                                                              hd.dropna('time'),
                                                              fd.dropna('time'),                                             
                                                            hazard_var,condition_var,
                                                            make_h_var_cat=make_h_var_cat,
                                                            p_dvs=p_dvs,
                                                            quantile=p,bin_num=bin_num)
                #format the decomposed quantities
                decomposed_df=decomp_to_pd_df(decomposed_hazard,model,s,r)
                #compute and format decomposition terms.
                #We don't currently save this.
                terms_df=decomp_to_term_pd_df(decomposed_hazard,model,s,r)
                
                #sum the terms over all bins, which we do save:
                summed_terms_df=terms_df.groupby(
                    ["model","season","region_id", "source", "term"], 
                    as_index=False
                )["value"].sum()
                
                os.makedirs('/'.join(decomp_path.split('/')[:-1]), exist_ok=True)
                os.makedirs('/'.join(term_path.split('/')[:-1]), exist_ok=True)
                #subprocess.run(["chmod", "-R", "g+rwx", ('/'.join(decomp_path.split('/')[:-1]))], check=True)
                #subprocess.run(["chmod", "-R", "g+rwx", '/'.join(term_path.split('/')[:-1])], check=True)
                
                all_summed_terms.append(summed_terms_df)
            except Exception as e:
                print(f"Error in season {s}, region {r}: {e}")
                _, term_path = get_savepaths(args, s, r)
                fallback_terms_df = get_zero_summed_terms_df(model, s, r)
                os.makedirs('/'.join(term_path.split('/')[:-1]), exist_ok=True)
                all_summed_terms.append(fallback_terms_df)
                

    final_terms_path, final_terms_members_path = get_final_savepaths(args)
    final_terms_columns = ["model", "season", "region_id", "source", "term", "value"]
    if all_summed_terms:
        final_terms_df = pd.concat(all_summed_terms, ignore_index=True)
    else:
        final_terms_df = pd.DataFrame(columns=final_terms_columns)

    os.makedirs(os.path.dirname(final_terms_path), exist_ok=True)
    if len(args.members) ==1 :    
        if type(args.members) is list:
            output_member = '-'.join(args.members)
        else:
            output_member = args.members[0] if args.members else 'r1i1p1f1'    
        final_terms_members_df = final_terms_df.copy()
        final_terms_members_df.insert(loc=5, column='member', value=output_member)
        final_terms_members_df.to_csv(final_terms_members_path)
    final_terms_df.to_csv(final_terms_path)
    print(f"Saved final terms CSV: {final_terms_path}")
    print(f"Saved final member terms CSV: {final_terms_members_path}")


    #these are each region and season specific
    ref_pca_proj= (args.ref_pca is not None)
    ref_bin_edges= (args.ref_bins is not None)
    ref_cat_thresh = (args.ref_thresh is not None)

    pca_patterns=None
    bins=None
    threshs=None

    pca_dict={}
    bin_dict={}
    thresh_dict={}

    #If we've predefined the PCs, bins or thresh, we load them now.
    if ref_pca_proj:
        
        """pca_dict=load_pca_pattern(args.auxdir, args.model, args.experiment, p_dvs)"""
        pca_dict=load_param_data(args.ref_pca)
    if ref_bin_edges:
        """bin_dict = load_bin_edges(
            args.auxdir,
            args.model,
            args.experiment,
            condition_var,
            args.nprecursorbins,
        )"""
        bin_dict=load_param_data(args.ref_bins)

    if ref_cat_thresh:
        """thresh_dict = load_cat_thresh(
            args.auxdir,
            args.model,
            args.experiment,
            args.hazardvariable,
            args.eventthreshold,
        )"""
        thresh_dict=load_param_data(args.ref_thresh)

    data=get_index_and_hazard_data(args)
    decomp_dfs=[]


    for season in args.seasons:
        for region_id in args.regions:

            key = f"{season}_{region_id}"

            sel_data=data.sel(season=season,region_id=region_id).dropna('time',how='all')
            data_val=validate_sel_data(sel_data,season,region_id)
            if data_val==1:
                continue
            #if we did not predefine PCs, bins or thresh, we compute them now
            if not ref_pca_proj:
                pca_pattern=compute_pca_pattern(sel_data,condition_var,p_dvs)
                pca_dict[key]=pca_pattern

            #either way, the param dicts are now populated:
            sel_data=handle_PCA_projection(sel_data,condition_var,p_dvs,pca_dict[key])

            if not ref_bin_edges:
                bins=compute_bin_edges(sel_data,condition_var,args.nprecursorbins)
                bin_dict[key]=bins

            h_var = args.hazardvariable

            if not ref_cat_thresh:
                thresh=compute_ref_thresh(sel_data, h_var, args.eventthreshold)
                thresh_dict[key]=thresh

            sel_data,h_var=handle_categorical_target(sel_data, h_var, make_h_var_cat,
                                                     thresh_dict[key])
            
            Ph_s, P_s = apply_decomposition(sel_data,h_var, condition_var, bin_dict[key])
            decomp_dfs.append(
                decomp_to_df(Ph_s, P_s, args.model, season, region_id)
            )

    decomp_df=postprocess_decomp_dfs(decomp_dfs)
    save_decomp_df(decomp_df,decomp_output_path)
    save_ensemble_documentation(decomp_output_path,args.members)
    if not ref_pca_proj:
        save_pca(pca_dict,args.auxdir,args.model, args.experiment, p_dvs)

    if not ref_bin_edges:
        save_bin_edges(
            bin_dict, args.auxdir,
            args.model, args.experiment,
            condition_var,args.nprecursorbins,
        )
    if not ref_cat_thresh:
        save_cat_thresh(
            thresh_dict, args.auxdir,
            args.model, args.experiment,
            args.hazardvariable, args.eventthreshold,
        )
if __name__=='__main__':

    args = parse_args()
    main(args)


