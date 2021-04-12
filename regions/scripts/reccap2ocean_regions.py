import pooch
import xarray as xr
import numpy as np


def main():

    reccap2_regions = make_reccap_region_options()
    # select option 4 as the prefered basin option
    xds = reccap2_regions[
        ['option4', 'fay_mckinley', 'woa_regions']
    ].rename(option4='basins')
    
    xds.fay_mckinley.values += 20
    fay = xds.fay_mckinley.fillna(1)
    rec = xds.basins
    woa = xds.woa_regions
    
    region_dict = {
        'atlantic': {'index': 1, 'non_fay_idx': [1,2]}, 
        'pacific':  {'index': 2, 'non_fay_idx': [1,2]}, 
        'indian':   {'index': 3, 'non_fay_idx': [3]}, 
        'arctic':   {'index': 4, 'non_fay_idx': [-99]}, 
        'southern': {'index': 5, 'non_fay_idx': [1,2]},
    }
    for r in region_dict:
        i = region_dict[r]['index']
        # Some regions fall in the global mask but are not in the 
        # fay biomes or in the woa mask. In this case we find these 
        # regions and interpolate longitudinally to include these
        # regions into the adjacent fay biomes
        non_fay_i = region_dict[r]['non_fay_idx']
        
        xda = make_region_numbers_consecutive(fay.where(rec == i))
        xda = fill_non_FM14_regions(xda, rec == i, non_fay_i)
        xds[r] = make_region_numbers_consecutive(xda)
    
    # mediterranean sea added to atlantic 
    med_mask = (
        (woa == 4)  # med sea is 4 in woa mask
        .astype(float)
        .where(lambda a: a == 1)  # convert False to 0 so we can fill
        .ffill('lat', limit=2)  # filling so that all med is included
        .bfill('lat', limit=2)
        .ffill('lon', limit=2)
        .bfill('lon', limit=1)
        .where(create_seamask())  # limit the filled part to our seamask
        .fillna(0)  # fill the nans with 0s
        .astype(bool))  # and convert the final to a boolean mask
    xds.atlantic.values[med_mask] = xds.atlantic.max().values + 1
    xds.pacific.values[xds.pacific == 7] = 0  # one pixel not caught in the masking
    
    xds['seamask'] = create_seamask()
        
    xds = (xds
           .drop(['fay_mckinley', 'woa_regions'])
           .rename(basins='open_ocean')
           .where(create_seamask())  # everything is clipped to our seamask
           .fillna(0)  # fill the clipped regions with 0
           .astype(int)  # and set the type as integer
          )
    
    encoding = {
        k: {'zlib': True, 'complevel': 4}
        for k in xds.data_vars
    }
    
    regions = ['atlantic', 'pacific', 'indian', 'arctic', 'southern']
    for key in regions:
        roll = 180 if (key == 'atlantic') or (key == 'open_ocean') else 0
        xds[key] = fill_missplaced_blobs(xds[key], shift_lon=roll).astype(int)
    
#     xds['coast'] = (xds.coastal_marcats > 0).astype(int)
#     xds = xds.drop('coastal_marcats')
    
    xds = add_final_attributes_and_names(xds)

    return xds


def add_final_attributes_and_names(xds):
    xds.attrs['description'] = (
        'regional masks created for the RECCAP2 ocean chapters. '
        'The region names are given in the attributes of each varible. '
        'The sub-regional masks are named according to Fay and McKinley 2014. '
        'Find the details of the naming here: https://doi.pangaea.de/10013/epic.42948.d018. '
        '\n\n'
        'Masks were created using the script located at: '
        'https://github.com/RECCAP2-ocean/shared-resources/blob/master/regions/scripts/reccap2ocean_regions.py. '
        'This mask is thus fully reproducible as long as the source data '
        'is also present for download. '
        )
    xds.attrs['helper_code'] = (
        "The following snippet of code can be used to list the region names "
        "as a dictionary:\n"
        "dict([r.strip().split('.') for r in region_names.split(',')])"
    )

    xds.open_ocean.attrs['region_names'] = (
        '1.Atlantic, 2.Pacific, 3.Indian, 4.Arctic, 5.Southern'
    )

    xds.atlantic.attrs['region_names'] = (
        '1.NA SPSS, '
        '2.NA STSS, '
        '3.NA STPS, '
        '4.AEQU, '
        '5.SA STPS, '
        '6.MED (not in FM14)'
    )

    xds.pacific.attrs['region_names'] = (
        '1.NP SPSS, '
        '2.NP STSS, '
        '3.NP STPS, '
        '4.PEQU-W, '
        '5.PEQU-E, '
        '6.SP STPS'
    )

    xds.indian.attrs['region_names'] = (
        '1.IND STPS, '
        '2.(not in FM14)'
    )

    xds.arctic.attrs['region_names'] = (
        '1.ARCTIC ICE (not in FM14), '
        '2.NP ICE, '
        '3.NA ICE, '
        '4.Barents (not in FM14)'
    )

    xds.southern.attrs['region_names'] = (
        '1.SO STSS, '
        '2.SO SPSS, '
        '3.SO ICE'
    )
    
    return xds
    
     
def fill_non_FM14_regions(sub_regions, global_mask, indicies):
    xda = sub_regions.astype(float)
    
    mask = ~xda.isin(indicies) & (xda != 0)
    
    filled = (
        xda
        .where(mask)
        .ffill('lon', 40)
        .bfill('lon', 40)
        .where(global_mask)
    )
    
    return filled
    
    
def make_region_numbers_consecutive(xda):
    unique_regions = np.unique(xda.values)
    mask = (~np.isnan(unique_regions)) | (unique_regions == 0)
    unique_regions = unique_regions[mask]

    for c, r in enumerate(unique_regions):
        c += 1

        mask = xda.values == r
        xda.values[mask] = c

    xda = xda.fillna(0).astype(int)
    
    return xda
    

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


def create_seamask():
    from pooch import retrieve
    from pandas import Timestamp
    from xarray import open_dataset
    from numpy import arange
    
    date = Timestamp('2010-01-01')
    url = (f"https://www.ncei.noaa.gov/data/"
           f"sea-surface-temperature-optimum-interpolation"
           f"/v2.1/access/avhrr/{date:%Y%m}/"
           f"oisst-avhrr-v02r01.{date:%Y%m%d}.nc")
    fname = retrieve(url, None)

    mask = (
        open_dataset(fname).sst
        .isel(time=0, zlev=0).drop(['time', 'zlev'])
        .interp(lat=arange(-89.5, 90, 1), lon=arange(0.5, 360))
        .notnull()
        .rename('seamask')
        .assign_attrs(dict(
            description=(
                "sea mask based on OISSTv2 coverage on "
                "2010-01-01 where True is sea and False is land"
    ))))
    
    return mask


def fill_missplaced_blobs(xda, shift_lon=0):
    seamask = create_seamask()
    def fill_in_islands(xda, shift_lon=0):
        import copy 
        out = xda
        mask = xda == 0
        for key in out.coords.keys():
            try:
                out = out.interpolate_na(key, limit=3, method='nearest')
                out = out.assign_coords({key: lambda a: a[key].values[::-1]})
                out = out.sortby(key)
                out = out.interpolate_na(key, limit=3, method='nearest')
                out = out.assign_coords({key: lambda a: a[key].values[::-1]})
                out = out.sortby(key)
            except ValueError:
                pass
        return out.where(seamask)
    
    def mask_pixel_islands(arr, shiftd1=0):
        from scipy import ndimage
        
        arr = np.roll(arr, shiftd1, axis=1)
        mask = []
        imax = int(np.nanmax(arr))
        for i in range(1, imax + 1):
            values = (arr == i).astype(int)

            kernel = [
                [0, 1, 0],
                [1, 1, 1],
                [0, 1, 0]
            ]
            labels = ndimage.label(values, structure=kernel)[0]
            idx, cnt = np.unique(labels, return_counts=True)

            keep = idx[cnt > 20]
            mask += np.isin(labels, keep),

        mask = np.all(mask, axis=0)
        mask = np.roll(mask, -shiftd1, axis=1)
        return mask

    mask = mask_pixel_islands(xda.values, shiftd1=shift_lon)
    out = xda.where(lambda a: (a!=0) & mask)
    out = fill_in_islands(out)
    out = out.fillna(0)
    return out


if __name__ == "__main__":
    main()