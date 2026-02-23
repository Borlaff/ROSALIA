import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import bottleneck as bn
import pandas as pd
from astropy.io import fits
# Custom modules
#sys.path.append("/Users/aborlaff/NASA/STRAYCOR/")

import rosalia as rs

def make_saturation_mask(input_name, ext, saturation_level=50000):
    output_satmasked = input_name.replace(".fits", "_ext" + str(ext) + "_satmasked.fits")
    rs.utils.execute_cmd("astarithmetic -h" + str(ext) + " " + input_name + " "+\
                         "set-i i i " + str(saturation_level) + " gt "+\
                         "2 dilate 2 dilate 2 dilate 2 dilate "+\
                         "nan where --output=" + output_satmasked)
    return(output_satmasked)


def make_gnuastro_segmentation_map(input_name, ext, saturation_level=50000, clean=True):
    # Masking saturated pixels

    input_name_clean_ext = rs.gnu.make_wcs_clean_extension(input_name, ext)

    output_satmasked = rs.gnu.make_saturation_mask(input_name=input_name_clean_ext,
                                                   ext=1, saturation_level=50000)

    ### Fill the gaps to be able to pass it to Noisechisel.
    # 1 - Make a kernel to convolve
    output_kernel = input_name.replace(".fits", "_kernel.fits")
    rs.utils.execute_cmd("astmkprof --kernel=gaussian,2,5 --oversample=1 --output=" + output_kernel)

    # 2 - Convolve the image by the kernel
    output_satmasked_conv = input_name.replace(".fits", "_satmasked_conv.fits")
    rs.utils.execute_cmd("astconvolve " + output_satmasked + " --kernel=" + output_kernel +\
              " --domain=spatial --output=" + output_satmasked_conv)

    # 3 - Fill the gaps with the Neighbours
    output_satfill = input_name.replace(".fits", "_satfill.fits")
    output_satfill_conv = input_name.replace(".fits", "_satfill_conv.fits")

    rs.utils.execute_cmd("astarithmetic " + output_satmasked + " 2 interpolate-maxofregion --output=" + output_satfill)
    rs.utils.execute_cmd("astarithmetic " + output_satmasked_conv + " 2 interpolate-maxofregion --output=" + output_satfill_conv)

    # 4 - Run noisechisel with the fixed versions of the image and the convolved one.
    output_noisechisel = input_name.replace(".fits", "_nc.fits")

    rs.utils.execute_cmd("astnoisechisel " + output_satfill + "  --convolved=" + output_satfill_conv +\
              " --output=" + output_noisechisel)

    output_segment = input_name.replace(".fits", "_ext" + str(ext) + "_sg.fits")
    rs.utils.execute_cmd("astsegment " + output_noisechisel + " --convolved=" + output_satfill_conv + " "+\
              " --gthresh=-10 --objbordersn=0 "+\
              " --output=" + output_segment + " --rawoutput")

    if clean:
        rs.utils.execute_cmd("rm " + output_satmasked)
        rs.utils.execute_cmd("rm " + output_kernel)
        rs.utils.execute_cmd("rm " + output_satmasked_conv)
        rs.utils.execute_cmd("rm " + output_satfill)
        rs.utils.execute_cmd("rm " + output_satfill_conv)
        rs.utils.execute_cmd("rm " + output_noisechisel)

    return(output_segment)



############################
def make_wcs_clean_extension(input_name, ext):
    input_fits = fits.open(input_name)
    wcsname_location_in_header = np.where(np.array(input_fits[ext].header) == "WCSNAME")[0][0]
    #print(wcsname_location_in_header)
    output_name = input_name.replace(".fits", "_ext" + str(ext) + ".fits")

    cmd = "astcrop " + input_name + " -K --mode=img --section=:,: -h" + str(ext) + " --hendwcs=" + str(wcsname_location_in_header) + " -o " +  output_name
    rs.utils.execute_cmd(cmd)

    return(output_name)



############################
def astmkcatalog(input_name, ext, clean=True, verbose=False):
    # This program comprises the steps to generate a catalog with Gnuastro on an image.
    # Input: input_name - Name of the fits to analyze
    # ext: Extension of the fits to analyze.

    image_wcs_cleaned_extension_name = rs.gnu.make_wcs_clean_extension(input_name, ext)

    rs.utils.execute_cmd("astnoisechisel -K -h1 " + image_wcs_cleaned_extension_name)

    detected_name = image_wcs_cleaned_extension_name.replace(".fits", "_detected.fits")
    rs.utils.execute_cmd("astsegment -K " + detected_name)

    segmented_name = image_wcs_cleaned_extension_name.replace(".fits", "_detected_segmented.fits")

    output_name = image_wcs_cleaned_extension_name.replace(".fits", "_detected_segmented_cat.fits")

    cmd = "astmkcatalog -K --ra --dec --x --y --clumpscat " + segmented_name + " --output " + output_name
    if not clean:
        if verbose:
            print(cmd)
    rs.utils.execute_cmd(cmd)

    #output_name = output_name.replace(".csv", "_c.csv")
    output_catalog = fits.open(output_name)

    if verbose:
        print("###############################################")
        print("Catalog for " + input_name + " ["+str(ext)+"]")
        print("saved in " + output_name)
        print("###############################################")

    if clean:
        rs.utils.execute_cmd("rm " + detected_name)
        rs.utils.execute_cmd("rm " + segmented_name)

    return({"output_name": output_name, "output_catalog": output_catalog})



#####################################################
