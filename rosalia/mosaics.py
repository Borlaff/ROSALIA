# mosaic.py
import astropy.wcs as astropy_wcs
from astropy.io import fits
import numpy as np 
from tqdm import tqdm 
import astropy.units as u
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd, find_optimal_celestial_wcs
import rosalia as rs
from astropy.wcs.wcsapi import BaseHighLevelWCS


def coadd_level2(data, astropywcs, optimal_wcs=None, resolution=None):
    input_data_for_reproject = list(zip(data, astropywcs))

    if resolution is not None:
        if rs.utils.check_numeric_type(resolution) == "astropy.units.Quantity":
            print("Assuming user-input resolution is on arcsecs")
            resolution = resolution*u.arcsec

    if optimal_wcs is None:
        optimal_wcs =  find_mosaic_wcs(input_data=input_data_for_reproject, 
                                       resolution=resolution, auto_rotate=False) #find_optimal_celestial_wcs(input_data=input_data_for_reproject, 
                      #                           resolution=resolution)

    output = reproject_and_coadd(input_data=input_data_for_reproject, 
                                output_projection=optimal_wcs[0], 
                                shape_out=optimal_wcs[1], combine_function="mean",
                                reproject_function=reproject_interp,
                                #progress_bar=True,
                                parallel=True,
                                intermediate_memmap=True)

    data = np.array(output[0].data)
    data[data==0] = np.nan
    return(data, optimal_wcs[0].to_header())



def make_mosaic(flc_list, exts, optimal_wcs=None, output="test.fits", resolution=1): 
    # Gather the data and astropywcs
    data_list = []
    astropywcs_list = []
    for flc_name in tqdm(flc_list):
        flc = fits.open(flc_name, memmap=True)
        for ext in tqdm(exts):
            astropywcs_list.append(astropy_wcs.WCS(flc[ext].header))
            data_list.append(flc[ext].data)

    if optimal_wcs is not None:
        optimal_wcs = find_mosaic_wcs([flc_list, exts], resolution=resolution)
        
    mos_data, mos_header = coadd_level2(data=data_list, astropywcs=astropywcs_list, resolution=resolution)
    rs.utils.save_fits(mos_data, output, mos_header)



def find_mosaic_wcs(input_data, resolution=None, auto_rotate=False):
    # This is a wrapper for from reproject.mosaicking import reproject_and_coadd, find_optimal_celestial_wcs
    # If input_data is a list(zip(data, astropywcs)) them proceed with normal eproject.mosaicking.find_optimal_celestial_wcs
    # If input_data is a list(str, np.int) then it is a list of files and extensions. Then you have to go and retrieve the 
    # data and astropywcs for each one and run the former option. 
    #print(input_data[0][0])
    #print(input_data[1][0])
    
    #if resolution is not None:
    #    resolution = resolution*u.arcsec

    #print("input_data[0][0]: ", type(input_data[0][0]))
    #print("input_data[1][0]: ", type(input_data[1][0]))
    #print(input_data[1][0])
    #print("isinstance(input_data[1][0], (BaseHighLevelWCS,): ", isinstance(input_data[1][0], (BaseHighLevelWCS,)))


    if isinstance(input_data[0][0], (np.ndarray,)) and isinstance(input_data[0][1], (BaseHighLevelWCS,)):
        print("DATA / WCS mode")
        optimal_wcs = find_optimal_celestial_wcs(input_data=input_data, 
                                                 resolution=resolution, auto_rotate=auto_rotate)
        return(optimal_wcs)

    if isinstance(input_data[0][0], (str,)) and isinstance(input_data[0][1], (np.int32, np.int16, np.int8, int)):
        print("FITS / EXT mode") 
        # Gather the data and astropywcs
        data_list = []
        astropywcs_list = []
        for flc_name in input_data[0]:
            flc = fits.open(flc_name, memmap=True)
            for ext in input_data[1]:
                astropywcs_list.append(astropy_wcs.WCS(flc[ext].header))
                data_list.append(flc[ext].data.shape)
            flc.close()
            
        input_data_for_reproject = list(zip(data_list, astropywcs_list))
        optimal_wcs = find_optimal_celestial_wcs(input_data=input_data_for_reproject, 
                                                 resolution=resolution,
                                                 auto_rotate=auto_rotate)
        
        return(optimal_wcs)


"""
SET FOR DEPRECATION - SAVE METHOD FOR PARALLEL USING MEMMAPS

def reproject_roman_wfi_fits(data_list, wcs_list, mosaic_wcs):
    import os
    from reproject import reproject_interp
    import copy
    reprojected_images = []

    # reference_fits = fits.open(reference_name, memmap=True)
    reference_header = mosaic_wcs.to_header() # reference_fits[reference_ext].header
    reference_shape = reference_fits[reference_ext].data.shape
    # Use all available CPU cores
    num_cpus = os.cpu_count()

    for data, wcs in tqdm(zip(data_list, wcs_list)):
        data=np.float32(data)
        array_out = np.memmap(filename='output.np', mode='w+',  
                              shape=reference_shape, dtype='float32')
        # print(wcs)
        reproject_interp(
            input_data=(data, wcs),
            output_projection=reference_header,
            parallel=num_cpus,  # Enables block-based multi-core processing
            block_size='auto',   # Automatically determines chunk size 
            return_footprint=False,
            output_array=array_out
            )

        reprojected_image = copy.deepcopy(array_out)
        reprojected_images.append(reprojected_image)
        del array_out


    return(reprojected_images, reference_header)
"""

def generate_mosaic(data, astropywcs, resolution=None):
    input_data_for_reproject = list(zip(data, astropywcs))

    if rs.utils.check_numeric_type(resolution) == "astropy.units.Quantity":
        print("Assuming user-input resolution is on arcsecs")
        resolution = resolution*u.arcsec

    optimal_wcs = rs.mosaics.find_mosaic_wcs(input_data=input_data_for_reproject, 
                                             resolution=resolution, auto_rotate=False)

    output = rs.mosaics.reproject_and_coadd(input_data=input_data_for_reproject, 
                                            output_projection=optimal_wcs[0], 
                                            shape_out=optimal_wcs[1], 
                                            combine_function="mean",
                                            reproject_function=reproject_interp,
                                            #progress_bar=True,
                                            parallel=True,
                                            intermediate_memmap=True)

    data = np.array(output[0].data)
    data[data==0] = np.nan
    return(data, optimal_wcs[0].to_header())