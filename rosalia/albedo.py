#from pyhdf import SD as pyhdf_SD
import healpy as hp
import numpy as np
import bottleneck as bn
from tqdm import tqdm


def read_MODIS_albedo(modis_file, scaling_factor=0.001):
    # Using this MODIS MCD43A3 v061 prodiuct https://lpdaac.usgs.gov/products/mcd43a3v061/
    # "MCD43C3.A2000059.061.2020041000147.hdf"
    modis_pre = pyhdf_SD.SD(modis_file, pyhdf_SD.SDC.READ)

    data2D = modis_pre.select("Albedo_WSA_vis")[:,:]
    data2D[data2D == np.max(data2D)] = 0.075/scaling_factor # Oceans are defaulted to the max level. Setting a temporary average value until we find something better
    data2D = data2D*scaling_factor
    return(data2D)


def find_albedo(lon, lat, healpix_map_name):
    # We read the HEALPIX map
    albedo_hp_map = hp.read_map(healpix_map_name, nest=True)
    albedo_out = hp.pixelfunc.get_interp_val(albedo_hp_map, theta=lon, phi=lat, nest=True, lonlat=True)
    #return(interp((x_map, y_map)))
    return(albedo_out)



def regrid_modis_to_healpix(modis_file, nside, outfile, scaling_factor=0.001, extension="Albedo_WSA_vis"):

    # First we read the MODIS file
    # "Albedo_WSA_vis"
    modis_pre = pyhdf_SD.SD(modis_file, pyhdf_SD.SDC.READ)

    data2D = modis_pre.select(extension)[:,:]
    data2D[data2D == np.max(data2D)] = 0.075/scaling_factor # Oceans are defaulted to the max level. Setting a temporary average value until we find something better
    data2D = data2D*scaling_factor

    # Setting up the cartesian coordinate arrays
    data2D_shape = data2D.shape
    max_x_lon = data2D_shape[1]
    max_y_lat = data2D_shape[0]

    x = np.linspace(0, max_x_lon-1, max_x_lon)
    y = np.linspace(0, max_y_lat-1, max_y_lat)
    yy, xx = np.meshgrid(y, x, indexing='ij')
    lon_grid = 360/max_x_lon*(xx - max_x_lon/2)
    lat_grid = 180/max_y_lat*(max_y_lat/2 - yy)

    # Register which Healpix pix each pixel belongs to
    hp_id = hp.ang2pix(nside=nside, theta=lon_grid, phi=lat_grid, lonlat=True, nest=True)

    # Go healpix by healpix averaging the cartesian coordinate values
    npix = hp.nside2npix(nside)
    regridded_map = np.zeros((npix,))

    for i_pix in tqdm(range(npix)):
        where_in_the_cartesian = (hp_id == i_pix)
        regridded_map[i_pix] = bn.nanmedian(data2D[where_in_the_cartesian])

    hp.fitsfunc.write_map(outfile, regridded_map, nest=True, overwrite=True)
    return(regridded_map)