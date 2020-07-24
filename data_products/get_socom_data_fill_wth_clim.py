import pooch
import xarray as xr
import seaflux as sf
import numpy as np
import pandas as pd


def main():
    pass


def calculate_fluxes(spco2, apco2, tempC, salt, preshPa, wind, ice):
    import seaflux as sf

    fluxes = sf.flux_bulk(
        tempC, salt, spco2, apco2, preshPa, wind,
        kw_func=sf.gas_transfer_CO2.k_Wa92,
        kw_scaling=16)

    xds = xr.Dataset()
    xds['fgco2'] = fluxes * (1 - ice)
    xds['apco2'] = apco2
    xds['ice'] = ice
    xds['area'] = sf.utils.area_grid()

    xds.fgco2.attrs['units'] = 'gC/m2/day'
    xds.fgco2.attrs['long_name'] = 'air-sea flux of CO2'
    xds.fgco2.attrs['description'] = (
        'Air-sea co2 flux where upward is positive: \n'
        'FCO2 = K0 * kw * (spco2 - apco2) * ice \n'
        'kw = Wanninkhof (1992) scaled to 16 cm/hr; \n'
        'K0 = Weiss (1974); \n'
        'apco2 = xCO2mbl * (P - pH2O), using pH2O from Dickson et al 2007; \n'
        'spco2 = data products in spco2; \n'
        'ice = (1 - sea ice fraction)\n'
        'ERA5 winds + MSLP, OSTIA SST, EN4 salinity were used as auxiliary data'
    )

    xds.apco2.attrs.pop('citation')
    xds.apco2.attrs['long_name'] = 'atmospheric partial pressure of CO2'
    xds.apco2.attrs['units'] = 'uatm'
    xds.apco2.attrs['source'] = 'https://www.esrl.noaa.gov/gmd/ccgg/mbl/'
    xds.apco2.attrs['description'] = (
        'Atmospheric pCO2 was calculated with:\n'
        'pCO2 = xCO2mbl * (P - pH2O)\n\n'
        'xCO2mbl = NOAAs xCO2 marine boundary layer product (surface) interpolated linearly along latitudes and extrapolated along longitudes;\n'
        'P = atmospheric pressure from ERA5;\n '
        'pH2O = water vapour pressure calculated using Dickson et al. (2007);'
    )

    xds.ice.attrs['long_name'] = 'sea ice fraction'
    xds.ice.attrs['units'] = '/100'
    xds.ice.attrs['description'] = 'Sea ice fraction from the OSTIAv2 product'

    return xds


def fill_coastal_with_scaled_clim():
    socom = get_socom_data().sel(time=slice('1985', '2018'))

    time = socom.time.dt.month
    spco2 = socom.spco2
    clim = socom.spco2_clim.sel(month=time).drop('month')

    seamask = spco2.sel(product='mpi_somffn').notnull()
    ensemble = pco2.drop(['nies_fnn', 'jma_mlr'], 'product')
    ensemble = ensemble.where(seamask).mean('product')

    factor = (ensemble / clim).mean(['lat', 'lon'])
    # scaled = ensemble.fillna(clim * factor)
    scaled = clim * factor

    time_mask = spco2.notnull().any(['lat', 'lon']).compute()
    filled = spco2.where(seamask).fillna(scaled).where(time_mask)

    xds = xr.Dataset()
    xds['spco2'] = spco2
    xds['spco2_filled'] = filled
    xds['spco2_clim_scaled'] = scaled
    xds['seamask'] = mask
    xds['scaling_factor'] = factor

    xds.attrs['author'] = 'luke.gregor@usys.ethz.ch'
    xds.attrs['date'] = pd.Timestamp.today().strftime('%Y-%m-%d')
    xds.attrs['description'] = (
        'Surface pCO2 based on a selection of data-based products. '
        'Intended to be used for the IPCC report. '
        'The products have been filled with climatological estimates of '
        'pCO2 that has been scaled. For more details on how scaling is '
        'done, see the description of the scaling_factor variable. '
    )

    xds.product.attrs['description'] = (
        'Names of the data-based products. See '
        'http://www.bgc-jena.mpg.de/SOCOM/Methods.html for more details. '
    )

    xds.spco2.attrs['long_name'] = 'surface partial pressure of CO2'
    xds.spco2.attrs['units'] = 'uatm'
    xds.spco2.attrs['description'] = (
        'pCO2 for data-based products. For more information see the SOCOM web '
        'page (http://www.bgc-jena.mpg.de/SOCOM/Methods.html)'
    )

    xds.scaling_factor.attrs['long_name'] = 'scaling factor for climatological spco2'
    xds.scaling_factor.attrs['description'] = (
        'Scaling is applied to spco2_climatology that is then used to fill spco2. '
        'Scaling factor was estimated as: \n'
        'factor = (spco2_climatology / ensemble.where(seamask)).mean(lat, lon)\n'
        'where ensemble is the mean of [mpi-somffn, csir-ml6, lsce-ffnn2, jena-mls]. '
        'The [nies-ffn, jma-mlr] methods were omitted from the ensemble as the '
        'temporal extent is only for 1990 onward and both products cover less than '
        '95% of the ice free ocean. '
    )

    xds.seamask.attrs['description'] = 'The mask applied to the original data, any missing data is filled with (spco2_climatology * spco2_clim_scaling)'

    xds.spco2_clim_scaled.attrs['long_name'] = 'scaled climatological surface ocean partial pressure of CO2'
    xds.spco2_clim_scaled.attrs['source'] = 'https://doi.org/10.5194/essd-2020-90'
    xds.spco2_clim_scaled.attrs['description'] = (
        'Scaled surface ocean pCO2 as estimated using the MPIULB-SOMFFN '
        'product (still in review). spco2_clim_scaled is used to fill spco2 '
        'before calculating `fgco2`. Scaling is performed using the '
        '`spco2_clim_scaling` variable. Taking `spco2_clim_scaled / '
        'spco2_clim_factor` will result in the original climatology. '
    )

    return xds


def get_socom_data():

    products = [
        get_mpisomffn(),
        get_jmamlr(),
        get_jenamls(),
        get_lsceffnn2(),
        get_niesfnn(),
        get_csirml6()
    ]

    spco2 = xr.concat(
        products, dim='product'
    ).assign_coords(
        product=[
            'mpi_somffn',
            'jma_mlr',
            'jena_mls',
            'lsce_ffnn2',
            'nies_fnn',
            'csir_ml6'])
    clim = get_mpiulbsomffn()

    # combining into a dataset
    xds = xr.Dataset()
    xds['spco2'] = spco2
    xds['spco2_clim'] = clim
    xds = xds.assign_coords(time=pd.DatetimeIndex(xds.time.values))

    return xds


def get_mpiulbsomffn():
    url = 'https://www.ncei.noaa.gov/data/oceans/ncei/ocads/data/0209633/MPI-ULB-SOM_FFN_clim.nc'
    fname = pooch.retrieve(
        url, None,
        fname='MPIULB-SOMFFN_clim.nc',
        path='../data-in/',
        downloader=pooch.HTTPDownloader(progressbar=True))

    xds = xr.open_dataset(fname)
    xda = xds.pco2.where(xds.pco2 > 0).coarsen(lat=4, lon=4).mean()
    xda = xda.rename('mpiulb_somffn').rename(time='month')

    return xda


def get_mpisomffn():
    url = 'https://www.nodc.noaa.gov/archive/arc0105/0160558/5.5/data/0-data/MPI_SOM-FFN_v2020/spco2_MPI-SOM_FFN_v2020.nc'
    fname = pooch.retrieve(
        url, None,
        fname='MPI-SOMFFN_v2020.nc',
        path='../data-in/',
        downloader=pooch.HTTPDownloader(progressbar=True))

    xds = xr.open_dataset(fname, drop_variables='date')
    xda = xds.spco2_raw.resample(time='1MS').mean()
    xda = xda.rename('mpi_somffn')

    return xda


def get_jenamls():
    url = 'http://www.bgc-jena.mpg.de/CarboScope/oc/INVERSION/OUTPUT/oc_v1.7_pCO2_daily.nc'
    username = '**'
    password = '**'
    fname = pooch.retrieve(
        url, None,
        fname='Jena-MLS_v1.7_pCO2.nc',
        path='../data-in/',
        downloader=pooch.HTTPDownloader(progressbar=True, auth=(username, password)))

    xds = xr.open_dataset(fname)
    xda = xds.pCO2

    xda = (
        xda
        .resample(mtime='1MS').mean('mtime')
        .interp(lat=np.arange(-89.5, 90), lon=np.arange(-179.5, 180), method='linear')
        .roll(lon=180, roll_coords=False)
        .interpolate_na('lon', limit=20)
        .roll(lon=-180, roll_coords=False)
        .rename(mtime='time')
        .rename("jena_mls")
    )

    return xda


def get_lsceffnn2():

    url = 'ftp://my.cmems-du.eu/Core/MULTIOBS_GLO_BIO_CARBON_SURFACE_REP_015_008/dataset-carbon-rep-monthly/{t:%Y}/{name}'
    username = '**'
    password = '**'

    flist = []
    for t in pd.date_range('1985-01', '2019', freq='1MS', closed='left'):
        fname = 'dataset-carbon-rep-monthly_{t:%Y%m}15T0000Z_P20191231T1545Z.nc'.format(t=t)
        flist += pooch.retrieve(
            url.format(t=t, name=fname), None,
            fname=fname,
            path='../data-in/LSCE-FFNNv2/',
            downloader=pooch.FTPDownloader(progressbar=True, username=username, password=password)),

    xds = xr.open_mfdataset(flist, combine='nested', concat_dim='time')
    xda = (
        (xds.spco2 * 9.867)
        .assign_coords(longitude=(xds.longitude - 180) % 360 - 180)
        .rename(latitude='lat', longitude='lon')
        .resample(time='1MS').mean()
        .sortby('lon')
        )

    return xda


def get_niesfnn():
    url = 'https://ndownloader.figshare.com/files/23907317?private_link=6dfc21bc1a2c51da8081'
    fname = pooch.retrieve(
        url, None,
        fname='NIES-FNN_v2020.nc',
        path='../data-in/',
        downloader=pooch.HTTPDownloader(progressbar=True))

    xds = xr.open_dataset(fname, drop_variables='date')

    yymm = np.meshgrid(xds.year, xds.month)
    years_months = np.c_[([y.flatten() for y in yymm])].T
    time = [pd.Timestamp(f'{y}-{m}') for y, m in years_months]

    xda = xr.DataArray(
        xds.co2.values.reshape(len(time), xds.lat.size, xds.lon.size),
        coords=dict(time=time, lat=xds.lat, lon=xds.lon),
        dims=['time', 'lat', 'lon']
    )

    return xda


def get_jmamlr():
    url = 'http://www.data.jma.go.jp/gmd/kaiyou/data/english/co2_flux/grid/{name}'

    xds = []
    for t in pd.date_range('1990-01', '2019', freq='1AS', closed='left'):
        fname = 'JMA_co2map_{t:%Y}.ZIP'.format(t=t)
        fname = pooch.retrieve(
            url.format(t=t, name=fname), None,
            fname=fname,
            path='../data-in/JMA-MLR/',
            processor=pooch.Unzip(),
            downloader=pooch.HTTPDownloader(progressbar=True))[0]
        xda = xr.open_dataset(fname, decode_times=False).pCO2s
        y0, y1 = str(t.year), str(t.year+1)
        time = pd.date_range(y0, y1, freq='1MS', closed='left')
        xda = xda.assign_coords(time=time)
        xds += xda,

    xda = (
        xr.concat(xds, dim='time')
        .assign_coords(lon=(xda.lon - 180) % 360 - 180)
        .sortby('lon')
    )

    return xda


def get_csirml6():
    url = "https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/23875943/CSIRML6_CO2_19822019_figshare.nc"

    fname = pooch.retrieve(
        url, None,
        fname="CSIRML6_CO2_19822019_figshare.nc",
        path='../data-in/',
        downloader=pooch.HTTPDownloader(progressbar=True))
    xda = xr.open_dataset(fname).spco2

    return xda


def get_somffn_flux_params():

    url = 'https://www.nodc.noaa.gov/archive/arc0105/0160558/5.5/data/0-data/MPI_SOM-FFN_v2020/spco2_MPI-SOM_FFN_v2020.nc'
    fname = pooch.retrieve(
        url, None,
        fname='MPI-SOMFFN_v2020.nc',
        path='../data-in/',
        downloader=pooch.HTTPDownloader(progressbar=True))

    drop = [
        'date',
        'dco2',
        'spco2_raw',
        'spco2_smoothed',
        'fgco2_raw',
        'fgco2_smoothed',
        'time_bnds',
        'lat_bnds',
        'lon_bnds'
    ]

    xds = xr.open_dataset(fname, drop_variables=drop)
    attrs = {k: xds[k].attrs for k in xds}
    xds = xds.resample(time='1MS').mean()
    for k in xds:
        xds[k].attrs = attrs[k]
    xds.attrs = {}

    return xds
