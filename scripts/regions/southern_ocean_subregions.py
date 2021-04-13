from reccap2_ocean_regions import get_CO2_biomes


def southern_ocean_subregions():

    biomes = get_CO2_biomes()

    lat = biomes.lat
    lon = biomes.lon

    so_biomes = (biomes.where(biomes >= 15) - 15)

    so_basins = biomes * 0 + 2
    so_basins.values[((lon > -70) & (lon <= 20) & (lat > -90)).T] = 0
    so_basins.values[((lon >= 20) & (lon < 147) & (lat > -90)).T] = 1

    so_regions = (so_biomes + so_basins * 3).to_dataset(name='so_subregions')

    so_regions['fay_mckinley'] = so_biomes
    so_regions['basins'] = so_basins.where(so_biomes.notnull())

    return so_regions


if __name__ == "__main__":
    southern_oucean_subregions()
