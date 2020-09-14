import pooch
import xarray as xr
import numpy as np


def main():

    reccap2_regions = make_reccap_region_options()
    # select option 4 as the prefered basin option
    xds = reccap2_regions[
        ['option4', 'fay_mckinley', 'woa_regions']
    ].rename(option4='basins')
    
    xds['atlantic'] = (xds.fay_mckinley.fillna( 8)).where(xds.basins == 1) - 7
    xds['pacific']  = (xds.fay_mckinley.fillna( 1)).where(xds.basins == 2)
    xds['indian']   = (xds.fay_mckinley.fillna(13)).where(xds.basins == 3) - 12
    xds['arctic']   = (xds.fay_mckinley.fillna( 3)).where(xds.basins == 4) - 2
    xds['southern'] = (xds.fay_mckinley.fillna(14)).where(xds.basins == 5) - 13
    
    xds = xds.drop(['fay_mckinley', 'woa_regions']).rename(basins='reccap2_ocean_regions')
    
    encoding = {
        k: {'zlib': True, 'complevel': 4}
        for k in xds.data_vars
    }

    xds.to_netcdf(
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
    # basin_list = {n+1: b for n, b in enumerate(woa_basins.CLIST[2:-2].split(') ('))}
    
    new_basin_num = {1:1, 2:2, 3:3, 56:3, 11:4, 10:5}
    
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
        ((biomes != 2) & (biomes != 9) | ((biomes.lon > 25) & (biomes.lon < 55))))
    
    xds = xr.Dataset()
    
    # option 1 has only the SO mask that is from Fay and McKinley
    xds['option1'] = basins.where(~so_mask, so_mask * 5)

    # option 2 uses Atlantic mask from F&M except east of Nordkapp
    opt2_mask = ~(so_mask | arctic_option2_mask | (basins==4))
    opt2_fill = (so_mask * 5) + (arctic_option2_mask * 4)
    xds['option2'] = basins.where(opt2_mask, opt2_fill)
    
    # option three is like option 2 but makes pacific "arctic" regions part of the pacific
    opt3_mask = (
        (xds.option1 != xds.option2) & 
        ((basins.lon < -120) | (basins.lon > 120))
    ).where(basins.notnull()).astype(bool)
    xds['option3'] = xds.option2.copy()
    xds.option3.values[opt3_mask] = 2
    xds['option3'] = xds.option3.where(xds.option2.notnull())
    
    # this is the option we'll propose to the greater RECCAPv2 community
    # after a discussion with Niki
    # fill out the sea of japan and do a little coastal padding elsewhere
    xds['option4'] = xds.option3.ffill('lon', limit=2).bfill('lon', limit=15)
    # add mediterranean to the ATL
    med = (woa_basins == 4).astype(int)  
    med = med.where(med == 1)
    # mask the following basins out
    mask = (
        (woa_basins == 5) |  # Baltic sea
        (woa_basins == 7) |  # Red sea
        (woa_basins == 8) |  # Persian gulf
        (woa_basins == 9)    # Hudson bay
    ).astype(int)
    # expand the masked regions a little to ensure apropriate coverage
    mask = (mask.where(mask == 1) 
        .ffill('lat', limit=2)
        .bfill('lat', limit=2)
        .ffill('lon', limit=2)
        .bfill('lon', limit=2) > 0)
    xds['option4'] = (
        xds['option4'].fillna(med)
        .where(~mask)
        .where(woa_basins.notnull())
        .ffill('lon', limit=1)
        .bfill('lon', limit=1)
        .ffill('lat', limit=1)
        .bfill('lat', limit=1)
    )
    # fixing NA miss-numbering
    xds.option4.values[xds.option4.values == 0] += 1
    
    mask = (xds.option4 == 4) & (xds.lat < 56) 
    xds.option4.values[mask] = 1
    
    xds['fay_mckinley'] = biomes
    xds['woa_regions'] = woa_basins
    
    xds = xds.roll(lon=180, roll_coords=True).assign_coords(lon=(xds.lon.values - 180) % 360)
    
    return xds


if __name__ == "__main__":
    main()
