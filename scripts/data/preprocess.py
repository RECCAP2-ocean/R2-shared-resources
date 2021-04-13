"""
Intended for data to be used with:
    `xr.open_mfdataset(fname, preprocess=reccap2_dataprep)`

Ensures that data is formatted so that it is formatted to match XYZT
    lon (0.5 : 360)
    lat (-89.5 : 90)
    depth (positive downward)
    time (days since 1980-01-01)
"""

import xarray as xr
import pandas as pd
import numpy as np


def add_processing(xds, message):
    time = pd.Timestamp.today().strftime('%Y-%m-%d')
    message = message.replace("'", "").replace('"', "")
    
    prefix = f" [RECCAP @ {time}]"
    prefix = ""
    hist = xds.attrs.get('processing', '').strip()
    hist = hist if hist.endswith(';') or hist == '' else hist + '; '
    hist += f'{prefix} {message}; '
    
    xds = xds.assign_attrs(processing=hist.strip())
    
    return xds


def try_run(f):
    import warnings
    import functools
    
    @functools.wraps(f)
    def wrapper(xds):
        try:
            return f(xds)
        except Exception as e:
            warnings.warn(f'{f.__name__} could not be run: {str(e)}')
            return xds
    return wrapper


@try_run
def rename_coords(
    xds,
    lat=['LATITUDE', 'latitude', 'Lat', 'ylat', 'y', 'Y'],
    lon=['LONGITUDE', 'longitude', 'Lon', 'xlon', 'x', 'X'],
    time=['TIME', 'Time', 'tmnth', 't', 'T'],
    region=['reg', 'regions'],
    **kwargs
):
    
    coord_names = dict(
        lat=lat, 
        lon=lon,
        time=time,
        region=region,
        **kwargs,
    )
    
    rename = {}
    for key in xds.coords:
        for new, old in coord_names.items():
            if key in old:
                rename[key] = new
                
    xds = xds.rename(**rename)
    
    if rename != {}:
        xds = add_processing(xds, f'renamed {str(rename)}')
    
    return xds


@try_run
def order_lon_360(xds):
    """assumes lon is in xds"""
    if 'lon' not in xds:
        return xds
    
    lon = xds.lon.copy().values
    xds = xds.assign_coords(lon=lambda a: a.lon % 360).sortby('lon')
    
    if any(lon != xds.lon.values):
        xds = add_processing(xds, 'flipped lon 0:360')

    return xds 


@try_run
def interpolate_coords(xds):
    dims = []
    if 'lat' in xds:
        dims += 'lat',
    if 'lon' in xds:
        dims += 'lon',
    
    must_interpolate = []
    for dim in dims:
        coords = xds[dim].values
        if ((coords + 0.5) % 1).sum():
            must_interpolate += dim, 
    if must_interpolate == []:
        return xds
    
    if 'lon' in must_interpolate:
        x = np.arange(0.5, 360)
        xds = xds.interp(lon=x)
        xds = add_processing(xds, 'interp lon onto grid centers')
    if 'lat' in must_interpolate:
        y = np.arange(-89.5, 90)
        xds = xds.interp(lat=y)
        xds = add_processing(xds, 'interp lat onto grid centers')
    
    return xds
        
    
@try_run
def transpose(xds):
    reccap_order = ["lon", "lat", "depth", "time"]
    coords = list(xds.dims)
    old_coords = list(xds.dims)
    
    order = []
    for key in reccap_order:
        if key in coords:
            order += key,
            coords.remove(key)
    order += coords
    
    if old_coords != order:
        xds = xds.transpose(*order)
        xds = add_processing(xds, f'transposed ({str(old_coords)[1:-1]})')

    return xds


@try_run
def time_decoder(xds):
    if 'time' not in xds:
        return xds

    if xds.time.attrs.get('units', None):
        xds = xr.decode_cf(xds)
        xds['time'] = xr.DataArray(
            xds.time.values.astype('datetime64[D]'), 
            dims=('time',))
        return xds
    
    t = xds.time.values.astype(int)

    # nies workaround (they use seconds since 1980)
    if t[0] > 10000:
        t = (t / 86400).astype(int)
        add_processing(xds, f'time units converted to days since from seconds')
    
    xds['time'] = xr.DataArray(
        data=t, dims=('time',), 
        coords={'time':t}, 
        attrs=dict(units='days since 1980', calendar='gregorian'))
    xds = xr.decode_cf(xds)
    
    return xds


@try_run
def center_time_on_15th(xds):
    if 'time' not in xds:
        return xds
    
    if xds.time.dt.day[0] != 15:

        dt = np.timedelta64(14, 'D')
        t = xds.time.values.astype('datetime64[M]') + dt

        xds['time'] = xr.DataArray(
            data=t, dims=['time'], coords={'time': t})

        xds = xds.sel(time=slice('1980', '2018'))
        xds = add_processing(xds, f'time centered on 15th')
    
    return xds
    

def preprocess(
    decode_times=False, 
    rename_coordinates=True, 
    center_months=True,
    interpolate_coordinates=True, 
    lon_0_360=True,
    transpose_dims=True,
):
    """
    A function generator that can be used as follows:
        reccap2_formatted_xds = xr.open_mfdataset(
            fname, decode_times=False,
            preprocess=data.preprocess(decode_times=True))
        
    There are several options that can be switched on or off with boolean.
    All changes are documented and appended to xds.attrs.processing
    """
    def reccap2_dataprep(ds):
        if rename_coordinates:
            ds = rename_coords(ds)
        if decode_times:
            ds = time_decoder(ds)
        if center_months:
            ds = center_time_on_15th(ds)
        if lon_0_360:
            ds = order_lon_360(ds)
        if interpolate_coordinates:
            ds = interpolate_coords(ds)
        if transpose_dims:
            ds = transpose(ds)
        
        key = list(ds.data_vars.keys())[0]
        ds[key].attrs['processing'] = ds.attrs.get('processing', '')
        
        return ds
    return reccap2_dataprep

