"""
Warning this script should not be used. 
We are using the data from the coastal group in the data folder
"""

print(__docs__)

import pooch
import numpy as np
import geopandas as gpd
import xarray as xr
from rasterio import features
from affine import Affine
from pathlib import Path


def main():
    # importing regions as the base template
    regions_url = 'https://github.com/RECCAP2-ocean/shared-resources/raw/master/regions/reccap2ocean_regions.nc'
    regions_fname = pooch.retrieve(regions_url, None)
    regions = (
        xr.open_dataset(regions_fname)
        .interp(lon=np.arange(-179.875, 180, .25), lat=np.arange(-89.875, 90, 0.25))
        .rename(lat='latitude', lon='longitude'))
    
    # getting shapefile information
    shapefile_url = 'http://www.hydrol-earth-syst-sci.net/17/2029/2013/hess-17-2029-2013-supplement.zip'
    shapefile_name = 'Continental_Shelf'
    shapefile_flist = pooch.retrieve(shapefile_url, None, processor=pooch.Unzip())
    shapefile_path = str(Path([f for f in shapefile_flist if shapefile_name in f][0]).parent)
    
    regions = add_shape_coord_from_data_array(regions, shapefile_path, shapefile_name)
    
    continental_shelf = regions[shapefile_name]
    continental_shelf = continental_shelf.coarsen(latitude=4, longitude=4).min()
    continental_shelf = continental_shelf.to_dataset(name='continental_shelf')
    
    continental_shelf.attrs = dict(
        source='https://www.hydrol-earth-syst-sci.net/17/2029/2013/',
        publication='Laruelle, G. G., Dürr, H. H., Lauerwald, R., Hartmann, J., Slomp, C. P., Goossens, N., and Regnier, P. A. G.: Global multi-scale segmentation of continental and coastal waters from the watersheds to the continental margins, Hydrol. Earth Syst. Sci., 17, 2029–2051, https://doi.org/10.5194/hess-17-2029-2013, 2013.',
        description='coastal zones as defined from the publication. ',
        history='data downloaded as shapefile (Continental_shelf.shp) from the link provided in `source`. The file was then converted to 1/4deg netCDF and then downscaled to 1deg.'
    )

    encoding = {
        k: {'zlib': True, 'complevel': 4}
        for k in continental_shelf.data_vars
    }

    continental_shelf.to_netcdf(
        '../data/regions/reccap2coastal_coscats.nc',
        encoding=encoding,
    )    

def transform_from_latlon(lat, lon):
    """ input 1D array of lat / lon and output an Affine transformation
    """
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    trans = Affine.translation(lon[0], lat[0])
    scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
    return trans * scale

def rasterize(shapes, coords, latitude='latitude', longitude='longitude',
              fill=np.nan, **kwargs):
    """Rasterize a list of (geometry, fill_value) tuples onto the given
    xray coordinates. This only works for 1d latitude and longitude
    arrays.

    usage:
    -----
    1. read shapefile to geopandas.GeoDataFrame
          `states = gpd.read_file(shp_dir+shp_file)`
    2. encode the different shapefiles that capture those lat-lons as different
        numbers i.e. 0.0, 1.0 ... and otherwise np.nan
          `shapes = (zip(states.geometry, range(len(states))))`
    3. Assign this to a new coord in your original xarray.DataArray
          `ds['states'] = rasterize(shapes, ds.coords, longitude='X', latitude='Y')`

    arguments:
    ---------
    : **kwargs (dict): passed to `rasterio.rasterize` function

    attrs:
    -----
    :transform (affine.Affine): how to translate from latlon to ...?
    :raster (numpy.ndarray): use rasterio.features.rasterize fill the values
      outside the .shp file with np.nan
    :spatial_coords (dict): dictionary of {"X":xr.DataArray, "Y":xr.DataArray()}
      with "X", "Y" as keys, and xr.DataArray as values

    returns:
    -------
    :(xr.DataArray): DataArray with `values` of nan for points outside shapefile
      and coords `Y` = latitude, 'X' = longitude.


    """
    transform = transform_from_latlon(coords[latitude], coords[longitude])
    out_shape = (len(coords[latitude]), len(coords[longitude]))
    raster = features.rasterize(shapes, out_shape=out_shape,
                                fill=fill, transform=transform,
                                dtype=float, **kwargs)
    spatial_coords = {latitude: coords[latitude], longitude: coords[longitude]}
    return xr.DataArray(raster, coords=spatial_coords, dims=(latitude, longitude))

def add_shape_coord_from_data_array(xr_da, shp_path, coord_name):
    """ Create a new coord for the xr_da indicating whether or not it 
         is inside the shapefile

        Creates a new coord - "coord_name" which will have integer values
         used to subset xr_da for plotting / analysis/

        Usage:
        -----
        precip_da = add_shape_coord_from_data_array(precip_da, "awash.shp", "awash")
        awash_da = precip_da.where(precip_da.awash==0, other=np.nan) 
    """
    # 1. read in shapefile
    shp_gpd = gpd.read_file(shp_path)

    # 2. create a list of tuples (shapely.geometry, id)
    #    this allows for many different polygons within a .shp file (e.g. States of US)
    shapes = [(shape, n) for n, shape in enumerate(shp_gpd.geometry)]

    # 3. create a new coord in the xr_da which will be set to the id in `shapes`
    xr_da[coord_name] = rasterize(shapes, xr_da.coords, 
                               longitude='longitude', latitude='latitude')

    return xr_da

main()
