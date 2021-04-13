def grid_soccom_argo_float(
    flist, 
    keep_vars=[
        'Pressure', 
        'Temperature', 
        'Salinity', 
        'pHinsitu',
        'pCO2_LIAR', 
        'TALK_LIAR', 
        'DIC_LIAR', 
        'TALK_LIAR_QFA', 
        'DIC_LIAR_QFA',
    ],
    depth=[
        0.0, 10.0, 20.0, 30.0, 50.0, 75.0, 100.0,
        125.0, 150.0, 200.0, 250.0, 300.0, 400.0,
        500.0, 600.0, 700.0, 800.0, 900.0, 1000.0, 
        1100.0, 1200.0, 1300.0, 1400.0, 1500.0, 
        1750.0, 2000.0, 2500.0, 3000.0, 3500.0, 
        4000.0, 4500.0, 5000.0, 5500.0],
    space_res=1,
    time_res='M',
    return_xarray=True,
    agg_funcs=['mean'],
):
    """
    Takes a list of soccom argo float filenames (netCDF).
    Returns gridded values as xr.Dataset or pd.DataFrame
    
    Parameters
    ----------
    flist : list
        Full file path of the netCDF files for the SOCCOM floats
    keep_vars: list
        A list of variables that will be kept in the output
    depth: list
        List of depth values (bin centers)
    space_res: float [1]
        Spatial resolution will be binned to the following resolution
    time_res: str 
        Will be resampled to that resolution, e.g. M/D for mon/day
    return_xarray: bool
        If True, returns output as xarray, otherwise as pd.DataFrame
    agg_funcs: list
        A list of the functions you would like perform for the gridding step
    """
    from numpy import around, inf, r_, convolve, array
    from pandas import concat, cut
    from xarray import open_dataset
    
    def round_coords(x):
        s = space_res
        r = s / 2.
        return round((x - r) / s) * s + r
        
    coords = ['Lat', 'Lon', 'Depth', 'JULD']
    
    df = []
    for f in flist:
        print('.', end='')
        xds = open_dataset(f)
        if not all([k in xds for k in keep_vars + coords]):
            continue
        else:
            df += xds.set_coords(coords)[keep_vars].to_dataframe().reset_index(),
            
    df = concat(df).reset_index()
    
    depth_bins = r_[-inf, convolve(depth, [0.5, 0.5], 'valid'), inf]
    
    df['lon'] = round_coords(df.Lon)
    df['lat'] = round_coords(df.Lat)
    df['depth'] = array(cut(df.Depth, depth_bins, labels=depth))
    t = df.JULD.astype(f'datetime64[{time_res}]')
    t0 = t.min()
    df['time'] = (t - t0).astype(f'timedelta64[{time_res}]').astype(int)
    
    new_coords = ['time', 'depth', 'lat', 'lon']
    grp = df.groupby(new_coords)
    agg_funcs = agg_funcs[0] if isinstance(agg_funcs, (list, tuple)) else agg_funcs
    agg = grp.aggregate(agg_funcs)
    agg = agg.loc[:, keep_vars].reset_index()
    
    # putting back time
    time = t0 + agg.time.astype(f'timedelta64[{time_res}]')
    agg['time'] = time.astype(f'datetime64[{time_res}]')
    
    if return_xarray:
        return agg.set_index(new_coords).to_xarray()
    else:
        return agg.reset_index()