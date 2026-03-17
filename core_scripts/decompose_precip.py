from genericpath import exists
import numpy as np
import sys
import os
import argparse
import pandas as pd
import xarray as xr
import subprocess
from sklearn.decomposition import PCA


def parse_args(arg_list=None):
    parser = argparse.ArgumentParser(description="Use precursor indices and precipitation indices to decompose model bias and trends.")

    # Required arguments
    parser.add_argument('--model', type=str, required=True,
                        help='name of model to use.')

    parser.add_argument('--future_experiment', type=str, required=True,
                        help='Future simulation to compute changes for, e.g. ssp370. Pass "none" for bias-only analysis.')
    
    
    # Optional arguments

    parser.add_argument('--eventthreshold',type=float,default=0.95,
                        help='extreme event percentile in referencemodel')
    
    parser.add_argument('--nprecursorbins',type=int,default=10,
                        help='Number of precursor bins to use.')
    
    parser.add_argument('--seasons', nargs='+', type=str, default=['DJF', 'MAM','JJA','SON'],
                    help='Which seasons to compute metrics for.')
    
    parser.add_argument('--regions', nargs='+', type=str, default=None,
                help='Which regions to compute metrics for. Default is to use all regions.')

    parser.add_argument('--members', type=str, default='',
                        help='Which ensemble member or members to use. Defaults to assuming not an ensemble.')
    
    parser.add_argument('--referencemodel', type=str, default='ERA5',
                    help='name of reference simulation.')

    parser.add_argument('--historical_experiment', type=str, default='historical',
                    help='name of experiment with historical overlap')

    parser.add_argument('--hist_period',type=int,nargs=2,default=[1979,2014],
                        help='The historical period to use.')

    parser.add_argument('--future_period',type=int,nargs=2,default=[2060,2100],
                        help='The future period to use.')
    
    parser.add_argument('--variables', nargs='+', type=str, default=['z500','u850','v850'],
                    help='Which precursors to use. If more than one provided, use PCA and take first PC.')

    parser.add_argument('--hazardvariable',type=str, default='pr',
                        help='Which variable defines the hazard index.')

    parser.add_argument('--overwrite',action='store_true',
                        help='Overwrite metrics if they already exist.')
    
    parser.add_argument('--inputdir',type=str,default='/Data/skd/projects/global/cmip6_precursors/outputs/indices/',
                        help='Directory in which to look for indices.')
    
    parser.add_argument('--auxdir',type=str,default='/home/rogui7909/code/CMIP_precursors/scripts/CMIP-EU-precip-precursors/aux',
                    help='Directory in which to look for decomposition.py')
    
    parser.add_argument('--savedir',type=str,default='/Data/skd/projects/global/cmip6_precursors/outputs/decompositions',
                    help='Directory in which to save output.')

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
            a1.append(xr.concat(a2, 'season', coords='different', compat='equals', join='outer').assign_coords(region_id=r))
        da=xr.concat(a1, 'region_id', coords='different', compat='equals', join='outer')
        # ds.append(da.rename(v))
        ds.append(da)
    ds=xr.merge(ds, compat='override', join='outer')

    y0,y1=args.hist_period
    # hist_period=np.arange(y0,y1+1)
    # ds=ds.isel(np.isin(ds['time.year'],hist_period))
    ds = ds.sel(time=slice(str(y0), str(y1)))
    
    return ds

def get_hist_data(args):
    bp=f'{args.inputdir}{args.model}/'

    # vars=np.atleast_1d(args.variables).append(args.hazardvariable)
    vars=list(np.atleast_1d(args.variables)) + [args.hazardvariable]
    ds=[]
    mems=np.atleast_1d(args.members)
    not_ens= (len(mems)==1) and (mems[0]=='')
    
    for v in vars:
        vdir=f'{bp}{v}/{args.historical_experiment}/'
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
                    da=xr.open_dataarray(f'{dir}{s}_region{r}.nc')
                    a2.append(da.assign_coords(season=s))
                a1.append(xr.concat(a2, 'season', coords='different', compat='equals', join='outer').assign_coords(region_id=r))

            a1=xr.concat(a1, 'region_id', coords='different', compat='equals', join='outer')
            if not not_ens: a1=a1.assign_coords(member=mem)
            a0.append(a1)

        a0=xr.concat(a0, 'member', coords='different', compat='equals', join='outer')#.rename(v)
        ds.append(a0)

    ds=xr.merge(ds, compat='override', join='outer')

    y0,y1=args.hist_period
    # hist_period=np.arange(y0,y1+1)
    # ds=ds.isel(np.isin(ds['time.year'],hist_period))
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
    mems=np.atleast_1d(args.members)
    not_ens= (len(mems)==1) and (mems[0]=='')
    
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
    member = args.members if args.members else 'r1i1p1f1'
    s1=f'{args.savedir}/{args.model}/{member}/'

    
    s2=f'{s}_region{r}.{suff}'
    return s1+'decomp_'+s2, s1+'terms_'+s2


def get_final_savepaths(args, suff='csv'):
    member = args.members if args.members else 'r1i1p1f1'
    base_dir = f'{args.savedir}/{args.model}/{member}/'
    return (
        f'{base_dir}{args.model}_paper_terms_df.{suff}',
        f'{base_dir}{args.model}_paper_terms_ens_df.{suff}'
    )


def get_zero_summed_terms_df(model, season, region_id):
    rows = []
    for source in ['conversion', 'dynamical', 'nonlinear']:
        for term in ['bias', 'change', 'uncalibrated_change']:
            rows.append(
                {
                    'model': model,
                    'season': season,
                    'region_id': region_id,
                    'source': source,
                    'term': term,
                    'value': 0,
                }
            )
    return pd.DataFrame(rows)


def run_decompose_precip(
    model, future_experiment, eventthreshold=0.95, nprecursorbins=10,
    seasons=['DJF', 'MAM','JJA','SON'], regions=None, members='', referencemodel='ERA5',
    historical_experiment='historical', hist_period=[1979,2014], future_period=[2060,2100],
    variables=['z500','u850','v850'], hazardvariable='pr', overwrite=False,
    inputdir='/Data/skd/projects/global/cmip6_precursors/outputs/indices/',
    auxdir='/Data/skd/projects/global/cmip6_precursors/aux/',
    savedir='/Data/skd/projects/global/cmip6_precursors/outputs/decompositions/'
):
    
    arg_list = (
        ["--model", model, "--future_experiment", future_experiment,
        '--eventthreshold', str(eventthreshold), '--nprecursorbins',
         str(nprecursorbins), "--seasons"]+seasons
        +['--members', members, '--referencemodel', referencemodel,
          '--historical_experiment', historical_experiment, "--variables"]
        +variables+['--hazardvariable', hazardvariable, '--inputdir', inputdir,
                    '--auxdir', auxdir, '--savedir', savedir])
    
    if overwrite:
        arg_list.append("--overwrite")

    if not regions is None:
        arg_list+=["--regions"]+[str(r) for r in regions]

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
        assert np.all([unique_vals==[0,1] for unique_vals in unique_h_vals])
    except:
        raise(ValueError(f'Expected binary event values, got: {unique_h_vals}'))

    #our bins are based on the reference ds. We've modified this to include two unobserved bins, at the high and low end. 
    bins=np.array([-1000,*ref_ds[s_var].quantile(np.linspace(0,1,bin_num+1)),1000])
    bins[1]=bins[1]-0.01 # so min ref value is not included in unobserved first bin
    bin_centres=(bins[1:]+bins[:-1])/2

    Ph_s_0, P_s_0 = decomp_func(ref_ds,h_var,s_var,bins)

    assert (Ph_s_0[0]==0)&(Ph_s_0[-1]==0)
############
    attrs_h_model_ds = h_model_ds[s_var].attrs
    h_model_ds[s_var].attrs = {}
    attrs_f_model_ds = f_model_ds[s_var].attrs
    f_model_ds[s_var].attrs = {}

    Ph_s_h, P_s_h = decomp_func(h_model_ds,h_var,s_var,bins)
    Ph_s_f, P_s_f = decomp_func(f_model_ds,h_var,s_var,bins)

    h_model_ds[s_var].attrs = attrs_h_model_ds
    f_model_ds[s_var].attrs = attrs_f_model_ds

    return output_func([Ph_s_0,Ph_s_h,Ph_s_f,P_s_0,P_s_h,P_s_f])

def return_decomp_as_dataarray(data):
    dims=('statistic','synoptic_bin')

    coords=dict(synoptic_bin=np.arange(len(data[0])),
                statistic=['therm_ref','therm_hist','therm_future',
                    'dyn_ref','dyn_hist','dyn_future']
                )
    return xr.DataArray(data=data,dims=dims,coords=coords)


def binned_decomposition(ds,h_var,s_var,bins):
    """helper func for decompose_hazard_odds_ratio"""
    # print(ds)
    # print(ds[s_var])
    Ph_s=ds.groupby_bins(s_var,bins=bins).mean()[h_var].fillna(0) #average value of hazard in bin. Is a probability for binary data. Bins with no hazard risk get a 0
    P_s=ds.groupby_bins(s_var,bins=bins).count()[s_var].fillna(0)/ds[s_var].time.size #occurence prob of synoptic bin. Bins that don't occur get a 0.

    return Ph_s, P_s

def fit_PCA_to_ds(ds,N,p_dvs):
    """helper func for decompose_hazard_odds_ratio"""

    X=ds[p_dvs].to_array('feature').T

    Solver=PCA(n_components=N,whiten=True).fit(X)
    principal_components=Solver.transform(X)
    
    #assuming ds has a dim 'time'
    PCs=xr.DataArray(data=principal_components,
        coords=dict(time=X.time.values,component=range(1,N+1)),dims=('time','component'))

    EOFs=xr.DataArray(data=Solver.components_,
        coords=dict(precursor=X.feature.values,component=range(1,N+1)),dims=('component','precursor'))
    EOFs=EOFs.assign_coords(explained_var=('component',Solver.explained_variance_ratio_))

    #add PC to ds
    ds['PC'+str(N)]=PCs.sel(component=N)
    return ds, Solver, PCs, EOFs

def apply_PCA_to_ds(ds,Solver,N,p_dvs):

    X=ds[p_dvs].to_array('feature').T

    return Solver.transform(X)[:,N-1]

def compute_terms_from_decomposition(Ph_s_0,Ph_s_h,Ph_s_f,P_s_0,P_s_h,P_s_f):

    N=len(Ph_s_0)
    #we make modifications for the unobserved_low and unobserved_hi bins (novel dynamics):
    for ix in [0,-1]:
        try:
            assert P_s_0[ix]==0
        except:
            raise(ValueError(f'The ref occurrence of bin {ix} should be 0 by design; instead {P_s_0[ix]}'))
        
        #if this bin occurred in the historical simulation, we treat its hazard risk as the 
        #true value for obs, as per our analytical analysis.
        #xi will then be 0, and the non-zero terms in odds_numerator can be interpreted directly.
        if P_s_h[ix]!=0:
            Ph_s_0[ix]=Ph_s_h[ix]

        #If it only occurs in the future sim, then we have no way of decomposing the contribution. 
        #We simply want to keep a trend contribution Ph_s_f*P_s_f.
        #Therefore we want Delta_P_s to be P_s_f, already satisfied as P_s_h=0. 
        # We set Ph_s_0[ix] to be Ph_s_f[ix]. Alpha, xi, and delta will be zero, so only dyn trend will remain, as desired.
        #if the future bin is also zero then this all works out the same.
        else:
            Ph_s_0[ix]=Ph_s_f[ix]

    #to handle bins that go from zero to nonzero hazard risk, we use a numerical trick, adding epsilons to Ph_s_0 and Ph_s_h.
    #The resulting error is at most epsilon*P_s, so bounded at 1*epsilon.
    epsilon=1e-12  
    Ph_s_0=Ph_s_0+epsilon
    Ph_s_h=Ph_s_h+epsilon

    #The most basic terms: ratios and anomalies of conditional and thermodynamic rainfall.

    #xi = therm_bias
    #alpha = therm_trend
    xi    = (Ph_s_h/Ph_s_0) -1
    alpha = (Ph_s_f/Ph_s_h) -1

    #when Ph_s_h is near 0 but Ph_s_0 isn't, alpha can't be sensibly estimated.
    # In this case we treat it as an additive contribution and set
    # alpha = Ph_s_f/Ph_s_0 so final therm trend will be P_s_0*Ph_s_f
    special_case= (Ph_s_h<1e-10)&(Ph_s_0>=1e-10)
    alpha[special_case]=Ph_s_f[special_case]/Ph_s_0[special_case]

    #little delta = dyn bias
    #Big Delta = dyn trend
    delta_P_s=P_s_h-P_s_0
    Delta_P_s=P_s_f-P_s_h

    # Estimates of bias, and the true trend based on these lowest coefficients

    #bias terms
    therm_bias=xi*Ph_s_0*P_s_0
    dyn_bias=delta_P_s*Ph_s_0
    nl_bias=xi*delta_P_s*Ph_s_0
    #trend terms
    therm_trend=alpha*P_s_0*Ph_s_0
    dyn_trend=Delta_P_s*Ph_s_0
    nl_trend=alpha*Delta_P_s*Ph_s_0

    #Relevance is a source-agnostic calculation of the contribution of each bin
    #to event probability
    ref_relevance=(Ph_s_0*P_s_0)
    hst_relevance=(Ph_s_h*P_s_h)
    fut_relevance=(Ph_s_f*P_s_f)

    #The bias free estimate of the future hazard.
    PH_f=(Ph_s_0*(P_s_0+ alpha*P_s_0 + Delta_P_s*(1+alpha))).sum()

    #true trend represented as a multiple of the model frequency
    beta=PH_f/ref_relevance.sum().values -1

    #direct naive estimate of trend
    beta_tilde=(fut_relevance.sum()/hst_relevance.sum().values) -1

    gamma=Delta_P_s/(P_s_0+epsilon)
    omega=delta_P_s/(P_s_h+epsilon)
    F_tilde=hst_relevance/hst_relevance.sum()
    F=ref_relevance/ref_relevance.sum()
    F_star=fut_relevance/fut_relevance.sum()

    therm_trend_error=alpha*(F_tilde-F)
    dyn_trend_error=gamma*(F_tilde*omega -F)
    nl_trend_error=gamma*alpha*(F_tilde*omega -F)

    #for returning data, we mask out alpha for the special case:
    alpha[special_case]=np.nan

    coefficient_da=xr.DataArray(data=[[xi,alpha],[delta_P_s,Delta_P_s]],
                 coords=dict(source=['therm','dyn'],term=['bias','trend'],synoptic_bin=np.arange(0,N)),
                 dims=('source','term','synoptic_bin'),name='coefficient')
    
    individual_term_da=xr.DataArray(data=[[therm_bias,therm_trend,therm_trend_error],
                                          [dyn_bias,dyn_trend,dyn_trend_error],
                                          [nl_bias,nl_trend,nl_trend_error]],
                 coords=dict(source=['therm','dyn','nonlinear'],term=['bias','trend','spurious_trend'],synoptic_bin=np.arange(0,N)),
                 dims=('source','term','synoptic_bin'),name='individual_term')
    
    beta_da=xr.DataArray(data=[beta,beta_tilde],coords=dict(term=['trend','spurious_trend']),dims=('term'),name='multiplicative_trend')

    contribution_da=xr.DataArray(data=[F,F_tilde,F_star],coords=dict(term=['ref','bias','trend'],synoptic_bin=np.arange(0,N)),
                                 dims=('term','synoptic_bin'),name='contribution')
    

    ds=xr.merge([coefficient_da,individual_term_da,beta_da,contribution_da])
    return ds

def blending_function(Ph_s_h,Ph_s_0,pow,mu):
    x=Ph_s_h/Ph_s_0
    return x**pow/(x**pow+mu**pow)


def compute_terms_from_decomposition_with_alpha_blending(Ph_s_0,Ph_s_h,Ph_s_f,P_s_0,P_s_h,P_s_f,blending_pow=4,blending_param=0.1):

    N=len(Ph_s_0)
    #we make modifications for the unobserved_low and unobserved_hi bins (novel dynamics):
    for ix in [0,-1]:
        try:
            assert P_s_0[ix]==0
        except:
            raise(ValueError(f'The ref occurrence of bin {ix} should be 0 by design; instead {P_s_0[ix]}'))
        
        #if this bin occurred in the historical simulation, we treat its hazard risk as the 
        #true value for obs, as per our analytical analysis.
        #xi will then be 0, and the non-zero terms in odds_numerator can be interpreted directly.
        if P_s_h[ix]!=0:
            Ph_s_0[ix]=Ph_s_h[ix]

        #If it only occurs in the future sim, then we have no way of decomposing the contribution. 
        #We simply want to keep a trend contribution Ph_s_f*P_s_f.
        #Therefore we want Delta_P_s to be P_s_f, already satisfied as P_s_h=0. 
        # We set Ph_s_0[ix] to be Ph_s_f[ix]. Alpha, xi, and delta will be zero, so only dyn trend will remain, as desired.
        #if the future bin is also zero then this all works out the same.
        else:
            Ph_s_0[ix]=Ph_s_f[ix]

    #to handle bins that go from zero to nonzero hazard risk, we use a numerical trick, adding epsilons to Ph_s_0 and Ph_s_h.
    #The resulting error is at most epsilon*P_s, so bounded at 1*epsilon.
    epsilon=1e-12  
    Ph_s_0=Ph_s_0+epsilon
    Ph_s_h=Ph_s_h+epsilon

    #The most basic terms: ratios and anomalies of conditional and thermodynamic rainfall.

    #xi = therm_bias
    #alpha = therm_trend
    xi    = (Ph_s_h/Ph_s_0) -1

    #alpha is generally multiplicative, but for very small Ph_s_h, this creates insane trends.
    #We add a smooth blending to an additive therm component for these cases.
    w=blending_function(Ph_s_h,Ph_s_0,blending_pow,blending_param)
    alpha=(1-w)*Ph_s_f/Ph_s_0 + w*((Ph_s_f/Ph_s_h) -1)

    #little delta = dyn bias
    #Big Delta = dyn trend
    delta_P_s=P_s_h-P_s_0
    Delta_P_s=P_s_f-P_s_h

    # Estimates of bias, and the true trend based on these lowest coefficients

    #bias terms
    therm_bias=xi*Ph_s_0*P_s_0
    dyn_bias=delta_P_s*Ph_s_0
    nl_bias=xi*delta_P_s*Ph_s_0
    #trend terms
    therm_trend=alpha*P_s_0*Ph_s_0
    dyn_trend=Delta_P_s*Ph_s_0
    nl_trend=alpha*Delta_P_s*Ph_s_0

    #Relevance is a source-agnostic calculation of the contribution of each bin
    #to event probability
    ref_relevance=(Ph_s_0*P_s_0)
    hst_relevance=(Ph_s_h*P_s_h)
    fut_relevance=(Ph_s_f*P_s_f)

    #The bias free estimate of the future hazard.
    PH_f=(Ph_s_0*(P_s_0+ alpha*P_s_0 + Delta_P_s*(1+alpha))).sum()

    #true trend represented as a multiple of the model frequency
    beta=PH_f/ref_relevance.sum().values -1

    #direct naive estimate of trend
    beta_tilde=(fut_relevance.sum()/hst_relevance.sum().values) -1

    #gamma isn't defined for our top and bottom unseen bins,
    # so we want it additively in these cases:
    w=blending_function(Delta_P_s,P_s_0,4,1e3)
    gamma=w*Delta_P_s/(P_s_0+epsilon) + (1-w)*Delta_P_s/(P_s_h+epsilon)

    omega=delta_P_s/(P_s_h+epsilon)
    F_tilde=hst_relevance/hst_relevance.sum()
    F=ref_relevance/ref_relevance.sum()
    F_star=fut_relevance/fut_relevance.sum()

    therm_trend_error=alpha*(F_tilde-F)
    dyn_trend_error=gamma*(F_tilde*omega -F)
    nl_trend_error=gamma*alpha*(F_tilde*omega -F)

    coefficient_da=xr.DataArray(data=[[xi,alpha],[delta_P_s,Delta_P_s]],
                 coords=dict(source=['therm','dyn'],term=['bias','trend'],synoptic_bin=np.arange(0,N)),
                 dims=('source','term','synoptic_bin'),name='coefficient')
    
    individual_term_da=xr.DataArray(data=[[therm_bias,therm_trend,therm_trend_error],
                                          [dyn_bias,dyn_trend,dyn_trend_error],
                                          [nl_bias,nl_trend,nl_trend_error]],
                 coords=dict(source=['therm','dyn','nonlinear'],term=['bias','trend','spurious_trend'],synoptic_bin=np.arange(0,N)),
                 dims=('source','term','synoptic_bin'),name='individual_term')
    
    beta_da=xr.DataArray(data=[beta,beta_tilde],coords=dict(term=['trend','spurious_trend']),dims=('term'),name='multiplicative_trend')

    contribution_da=xr.DataArray(data=[F,F_tilde,F_star],coords=dict(term=['ref','bias','trend'],synoptic_bin=np.arange(0,N)),
                                 dims=('term','synoptic_bin'),name='contribution')
    

    ds=xr.merge([coefficient_da,individual_term_da,beta_da,contribution_da])
    return ds
    
def decomp_to_pd_df(arr,model,season,region_id):
    """
    Takes in a set of 6x(binsize+2) decomp parameters.
    Returns a dataframe with columns:
    model [string], user passed.
    season [string], user passed.
    region_id [int or string], user passed.
    bin  [int]
    source ['conversion' or 'dynamical']
    period ['reference', 'historical' or 'future']
    value [float]
    Expected number of rows = 6x(bin_num+2)
    """
    Ph_s_0, Ph_s_h, Ph_s_f, P_s_0, P_s_h, P_s_f = arr.values

    data=dict(
    dyn_ref= (P_s_0,'dynamical','reference'),
    dyn_hst= (P_s_h,'dynamical','historical'),
    dyn_ftr= (P_s_f,'dynamical','future'),
    cnv_ref= (Ph_s_0,'conversion','reference'),
    cnv_hst= (Ph_s_h,'conversion','historical'),
    cnv_ftr= (Ph_s_f,'conversion','future')
    )
    decomp_rows=[]

    for name, (vals, s1, s2) in data.items():
        for b, v in enumerate(vals, 1):
            decomp_rows.append({"model": model, "season": season,
                        "region_id":region_id, "bin": b,
                        "source": s1, "period": s2, "value": v})
                        
    return pd.DataFrame(decomp_rows)

def decomp_to_term_pd_df(arr,model, season, region_id):
    """
    Takes in a set of 6x(binsize+2) decomp parameters.
    Returns a dataframe with columns:
    model [string], user passed.
    season [string], user passed.
    region_id [int or string], user passed.
    bin  [int]
    source ['conversion', 'dynamical' or 'nonlinear']
    term ['bias', 'change' or 'uncalibrated_change']
    value [float]
    Expected number of rows = 9x(bin_num+2)
    """

    Ph_s_0, Ph_s_h, Ph_s_f, P_s_0, P_s_h, P_s_f = arr.values
    
    # for numerical stability
    epsilon=1e-12  
    Ph_s_0=Ph_s_0+epsilon
    Ph_s_h=Ph_s_h+epsilon

    #mult-> to-> add therm alpha for very biased cases
    blending_pow=4
    blending_param=0.1
    w=blending_function(Ph_s_h,Ph_s_0,blending_pow,blending_param)

    delta_P_s=P_s_h-P_s_0
    Delta_P_s=P_s_f-P_s_h

    xi=(Ph_s_h/Ph_s_0)-1
    alpha=(1-w)*Ph_s_f/Ph_s_0 + w*((Ph_s_f/Ph_s_h) -1)#equiv to 1+alpha =Ph_s_f/Ph_s_h in most cases

    raw_alpha=(Ph_s_f/Ph_s_h) -1

    data=dict(
        dyn_bias= (delta_P_s*Ph_s_0,'dynamical','bias'),
        cnv_bias= (P_s_0*xi*Ph_s_0,'conversion','bias'),
        nlr_bias= (delta_P_s*xi*Ph_s_0,'nonlinear','bias'),

        dyn_caltrend=(Delta_P_s*Ph_s_0,'dynamical','change'),
        cnv_caltrend=(P_s_0*alpha*Ph_s_0,'conversion','change'),
        nlr_caltrend=(Delta_P_s*alpha*Ph_s_0,'nonlinear','change'),

        dyn_rawtrend=(Delta_P_s*Ph_s_h,'dynamical','uncalibrated_change'),
        cnv_rawtrend=(P_s_h*raw_alpha*Ph_s_h,'conversion','uncalibrated_change'),
        nlr_rawtrend=(Delta_P_s*raw_alpha*Ph_s_h,'nonlinear','uncalibrated_change')
    )

    term_rows=[]
    for name, (vals, s1, s2) in data.items():
        for b, v in enumerate(vals, 1):
            term_rows.append({"model": model, "season": season,
                        "region_id":region_id, "bin": b,
                        "source": s1, "term": s2, "value": v})
            
    return pd.DataFrame(term_rows)




    
def main(args):
    sys.path.append(args.auxdir)

    condition_var=args.variables
    if len(condition_var)==1:
        condition_var=condition_var[0]
    else:
        #This principal component is fitted in the reference dataset and projected
        #on the models.
        condition_var='PC1'
    
    p_dvs=[v+'_lag0_index_val1' for v in args.variables]


    #use a pre-established set of regions if none specified
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
                hd = xr.concat([hd.sel(member=k, drop=True) for k in hd.member], dim='time', coords='different', compat='equals', join='outer')
                fd = xr.concat([fd.sel(member=k, drop=True) for k in fd.member], dim='time', coords='different', compat='equals', join='outer')
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

    output_member = args.members if args.members else 'r1i1p1f1'
    final_terms_members_df = final_terms_df.copy()
    final_terms_members_df.insert(loc=5, column='member', value=output_member)
    os.makedirs(os.path.dirname(final_terms_path), exist_ok=True)
    final_terms_df.to_csv(final_terms_path)
    final_terms_members_df.to_csv(final_terms_members_path)
    print(f"Saved final terms CSV: {final_terms_path}")
    print(f"Saved final member terms CSV: {final_terms_members_path}")


if __name__=='__main__':

    args = parse_args()
    main(args)
