from xarray import register_dataarray_accessor


def earth_radius(lat):
    from numpy import deg2rad, sin, cos

    lat = deg2rad(lat)
    a = 6378137
    b = 6356752
    r = (
        ((a ** 2 * cos(lat)) ** 2 + (b ** 2 * sin(lat)) ** 2)
        / ((a * cos(lat)) ** 2 + (b * sin(lat)) ** 2)
    ) ** 0.5

    return r


def area_grid(lat, lon, return_dataarray=False):
    """Calculate the area of each grid cell for a user-provided
    grid cell resolution. Area is in square meters, but resolution
    is given in decimal degrees.
    Based on the function in
    https://github.com/chadagreene/CDT/blob/master/cdt/cdtarea.m
    """
    from numpy import meshgrid, deg2rad, gradient, cos

    ylat, xlon = meshgrid(lat, lon)
    R = earth_radius(ylat)

    dlat = deg2rad(gradient(ylat, axis=1))
    dlon = deg2rad(gradient(xlon, axis=0))

    dy = dlat * R
    dx = dlon * R * cos(deg2rad(ylat))

    area = dy * dx

    if not return_dataarray:
        return area
    else:
        from xarray import DataArray

        xda = DataArray(
            area.T,
            dims=["lat", "lon"],
            coords={"lat": lat, "lon": lon},
            attrs={
                "long_name": "area_per_pixel",
                "description": "area per pixel",
                "units": "m^2",
            },
        )
        return xda


@register_dataarray_accessor("area")
class SeafluxUtils:
    def __init__(self, data):
        self._obj = data

    def __call__(self, lat_name="lat", lon_name="lon"):
        """
        Returns the area of the grid cells if lat and lon
        are present. You can adjust the names of lat and lon.
        Output units are in m^2
        """

        xda = self._obj

        assert lat_name in xda.coords, f"{lat_name} is not in data array"
        assert lon_name in xda.coords, f"{lon_name} is not in data array"

        lat = xda[lat_name].values
        lon = xda[lon_name].values

        area = area_grid(lat, lon, return_dataarray=True)
        area = area.rename(lat=lat_name, lon=lon_name)

        return area