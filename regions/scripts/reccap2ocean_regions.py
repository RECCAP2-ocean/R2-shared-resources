import pooch
import xarray as xr
import numpy as np


def main():

    reccap2_regions = make_reccap_region_options()

    encoding = {
        k: {'zlib': True, 'complevel': 4}
        for k in reccap2_regions.data_vars
    }

    reccap2_regions.to_netcdf(
        '../reccap2ocean_regions.nc',
        encoding=encoding,
    )


def get_CO2_biomes():
    url = "https://epic.awi.de/id/eprint/34786/19/Time_Varying_Biomes.nc"
    fname = pooch.retrieve(url, None)
    xda = (
        xr.open_dataset(fname)
        .transpose('year', 'lat', 'lon')
        .rename({'MeanBiomes': 'mean_biomes'})
        .mean_biomes)
    return xda


def get_woa_basins():
    url = (
        "https://iridl.ldeo.columbia.edu/"
        "SOURCES/.NOAA/.NODC/.WOA09/.Masks/.basin/data.nc"
    )
    fname = pooch.retrieve(url, None)
    xda = (
        xr.open_dataset(fname)
        .rename({'X': 'lon', 'Y': 'lat', 'Z': 'depth'})
        .transpose('depth', 'lat', 'lon')
        .basin
        .assign_coords(lon=(np.arange(0.5, 360) - 180) % 360 - 180)
        .sortby('lon')
        .sel(depth=0)
        .drop('depth'))
    return xda


def make_reccap_region_options():
    """
    Make RECCAP2 regions from Fay and McKinley's biomes and
    WOA 09 basins
    """

    biomes = get_CO2_biomes()
    woa_basins = get_woa_basins()

    # will be replaced
    basins = woa_basins.copy() * np.nan
    # one can run basin_list to see how I made the new_basin_num
    # basin_list = {
    #     n+1: b for n, b in enumerate(woa_basins.CLIST[2:-2].split(') ('))
    # }
    new_basin_num = {1: 1, 2: 2, 3: 3, 56: 3, 11: 4, 10: 5}
    for k, v in new_basin_num.items():
        mask = woa_basins.values == k
        basins.values[mask] = v

    # define the SO boundary
    so_mask = biomes > 14
    # The defining the arctic is a bit more tricky
    arctic_option2_mask = (
        # add arctic from woabasins + fay and mcKinley biomes for arctic
        ((basins == 4) | (biomes == 1) | (biomes == 8)) &
        # now remove woa North Atlantic but keep region east of nordkapp
        (
            (biomes != 2) & (biomes != 9) |
            ((biomes.lon > 25) & (biomes.lon < 55))
        )
    )

    xds = xr.Dataset()

    # option 1 has only the SO mask that is from Fay and McKinley
    xds['option1'] = basins.where(~so_mask, so_mask * 5)

    # option 2 uses Atlantic mask from F&M except east of Nordkapp
    opt2_mask = ~(so_mask | arctic_option2_mask | (basins == 4))
    opt2_fill = (so_mask * 5) + (arctic_option2_mask * 4)
    xds['option2'] = basins.where(opt2_mask, opt2_fill)

    # option three is like option 2 but makes pacific "arctic" regions
    # part of the pacific
    opt3_mask = (
        (xds.option1 != xds.option2) &
        ((basins.lon < -120) | (basins.lon > 120))
    ).where(basins.notnull()).astype(bool)
    xds['option3'] = xds.option2.copy()
    xds.option3.values[opt3_mask] = 2
    xds['option3'] = xds.option3.where(xds.option2.notnull())

    xds['fay_mckinley'] = biomes
    xds['woa_regions'] = woa_basins

    return xds


if __name__ == "__main__":
    main()
