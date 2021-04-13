def subplot_map(
    pos=111,
    proj=None, 
    round=True, 
    land_color='w', 
    **kwargs
):
    """
    Makes an axes object with a cartopy projection for the current figure
    
    Parameters
    ----------
    pos: int/list [111]
        Either a 3-digit integer or three separate integers
        describing the position of the subplot. If the three
        integers are *nrows*, *ncols*, and *index* in order, the
        subplot will take the *index* position on a grid with *nrows*
        rows and *ncols* columns. *index* starts at 1 in the upper left
        corner and increases to the right.

        *pos* is a three digit integer, where the first digit is the
        number of rows, the second the number of columns, and the third
        the index of the subplot. i.e. fig.add_subplot(235) is the same as
        fig.add_subplot(2, 3, 5). Note that all integers must be less than
        10 for this form to work.
    proj: crs.Projection()
        the cartopy coord reference system object to create the projection.
        Defaults to crs.PlateCarree(central_longitude=-155) if not given
    round: bool [True]
        If the projection is stereographic, round will cut the corners and
        make the plot round
    land_color: str ['w']
        the color of the land patches
    **kwargs:
        passed to fig.add_subplot(**kwargs)
        
    Returns
    -------
    An axes object attached to the current figure. If no figure, a figure 
    will be created as is default for plt.gcf()
    """
    from cartopy import feature, crs
    import matplotlib.path as mpath
    import matplotlib.pyplot as plt
    import numpy as np
    
    proj = crs.PlateCarree(-155) if proj is None else proj
    
    fig = plt.gcf()
    ax = fig.add_subplot(pos, projection=proj, **kwargs)

    # makes maps round
    stereo_maps = (
        crs.Stereographic,
        crs.NorthPolarStereo,
        crs.SouthPolarStereo,
    )
    if isinstance(ax.projection, stereo_maps) & round:

        theta = np.linspace(0, 2 * np.pi, 100)
        center, radius = [0.5, 0.5], 0.475
        verts = np.vstack([np.sin(theta), np.cos(theta)]).T
        circle = mpath.Path(verts * radius + center)

        ax.set_boundary(circle, transform=ax.transAxes)

    # adds features
    ax.add_feature(feature.LAND, color=land_color, zorder=4)
    ax.add_feature(feature.COASTLINE, lw=0.5, zorder=4)
    ax.outline_patch.set_lw(0.5)
    ax.outline_patch.set_zorder(5)
        
    return ax


def extend_lon_for_contour_plots(xds):
    import numpy as np
    x = np.arange(0.5, 361)
    xds = xds.sel(lon=x, method='nearest').assign_coords(lon=x)
    return xds


def add_southern_ocean_subregion_contours(ax, basin_lines=[20, 147, 290], **kwargs):
    from cartopy import crs
    
    biomes = _get_southern_ocean_subregions().biomes.where(lambda a: a!= 0)
    biomes = ((biomes + 1) % 2 + 1).fillna(0)
    biomes = extend_lon_for_contour_plots(biomes)
    biomes = biomes.where(lambda a: (a > 0) | (a.lat > -70))
        
    props = dict(
        ax=ax,
        transform=crs.PlateCarree(), 
        add_colorbar=False, 
        colors=['k'], 
        linewidths=[0.5], 
        levels=[0.5, 1.5], 
        zorder=5)
    props.update(**kwargs)
    biomes.plot.contour(**props)
    
    props = dict(
        transform=crs.PlateCarree(),
        color=props['colors'][0],
        linewidth=props['linewidths'][0],
        zorder=props['zorder'])
    for line in basin_lines:
        ax.plot([line] * 2, [-90, 90], **props)
        
    ax.set_extent([-180, 180, -90, -30], crs=props['transform'])
        
        
def _get_southern_ocean_subregions(
    url='https://github.com/RECCAP2-ocean/shared-resources/raw/master/regions/RECCAP2_region_masks_all.nc',
    dest='../data/regions/'
):
    import pooch
    import xarray as xr
    import pandas as pd
    from pathlib import Path as posixpath
    import itertools
    
    fname = pooch.retrieve(url, None, posixpath(url).name, dest)
    ds = xr.open_dataset(fname)

    mask = ds.southern
    
    atlantic = (((mask.lon > 290) | (mask.lon <=  20)) & (mask > 0)).astype(int) * 1
    indian   = (((mask.lon >  20) & (mask.lon <= 147)) & (mask > 0)).astype(int) * 2
    pacific  = (((mask.lon > 147) & (mask.lon <= 290)) & (mask > 0)).astype(int) * 3

    mask = xr.Dataset()
    mask['biomes'] = ds.southern.copy()
    mask['basins'] = (pacific + atlantic + indian).transpose('lat', 'lon')
    
    mask['subregions'] = (mask.basins * 3 + mask.biomes - 3).where(lambda a: a>0).fillna(0).astype(int)
    
    basin = ['ATL', 'IND', 'PAC']
    biome = ['STSS', 'SPSS', 'ICE']
    names = ['-'.join(l) for l in itertools.product(basin, biome)]    
    mask['names'] = xr.DataArray(names, coords={'idx': range(1, 10)}, dims=('idx'))
    mask['names'].attrs['description'] = 'Names for the subregions'
    
    mask['subregions'].attrs['description'] = '(basins * 3 + biomes - 3)'
    mask['basins'].attrs['description'] = 'Atlantic = 1, Indian = 2, Pacific = 3'
    mask['biomes'].attrs['description'] = 'Biomes based on Fay and McKinley (2014), STSS=1, SPSS=2, ICE=3'
    mask.attrs['source'] = url
    mask.attrs['date'] = pd.Timestamp.today().strftime('%Y-%m-%d')
    return mask


