"""
ROSALIA: ROman Sky Analyst for Low surface brightness Imaging & Astronomy

Alejandro S. Borlaff
NASA Ames Research Center, Moffett Field, 94035, California, USA.
a.s.borlaff@nasa.gov

This script sets up the NDI files from the FRED Optical outputs

It replaces rosalia/notebooks/DEVEL/R2_organize_NDI_files.ipynb

FOR INTERNAL USE ONLY
"""

import rosalia as rs
import numpy as np
import glob
import os
from astropy.io import fits
import pandas as pd 
from tqdm import tqdm
import matplotlib.pyplot as plt

"""
# Level 1,2,3 - Diagnostic movie
"""
subarray_locations_db = rs.roman.get_subarray_locations(SCA=1, verbose=False)
SCA_list = np.linspace(1,18,18, dtype="int")


if False:
    ndi_prefix = "lvl1"
    output_list = []
    counter = 0

    # Level 1 
    for SCA_i in tqdm(SCA_list):
        # print(SCA_i)
        for label_index in range(len(subarray_locations_db["xlabel"])):

            # Get the labels    
            xlabel = subarray_locations_db["xlabel"][label_index]
            ylabel = subarray_locations_db["ylabel"][label_index]

            # Find the right original files
            NDI_base_file = glob.glob('/Users/aborlaff/NASA/ROSALIA_DEPOT/SCOTT_SIMS_Jan2025_FPA_base/*SCA ' + str(SCA_i) + " " + str(xlabel) + " " + str(ylabel) + ".fits")
            print(NDI_base_file)
            NDI_fits = fits.open(NDI_base_file[0])
            
            fig, axs = plt.subplots(1, 1, figsize=(15, 10)) 
            im = axs.imshow(NDI_fits[0].data, vmin=0, vmax=1, origin='lower', cmap="RdYlBu_r")
            plt.colorbar(im, label='NDI')
            plt.text(100, 1000, "SCA " + str(SCA_i) + " - X:" + str(xlabel) + " Y:" + str(ylabel), color="yellow", fontsize=20) 
            plt.tight_layout()
            outname =  "level1_snap" + str(counter).zfill(5) + "_" + os.path.basename(NDI_base_file[0]).replace(".fits", ".png")
            plt.savefig(outname, dpi=100)
            output_list.append(outname)
            #  plt.show()
            counter = counter + 1
    os.system("magick -delay 100 -loop 0 level1_snap*.png deg2_ndi_movie_aug2025.mpeg")



    # Level 2 
if False:
    ndi_prefix = "lvl2"
    output_list = []
    counter = 0

    for SCA_i in tqdm(SCA_list):
        # print(SCA_i)
        for label_index in range(len(subarray_locations_db["xlabel"])):

            # Get the labels    
            xlabel = subarray_locations_db["xlabel"][label_index]
            ylabel = subarray_locations_db["ylabel"][label_index]

            # Find the right original files
            NDI_base_file = glob.glob('/Users/aborlaff/NASA/ROSALIA_DEPOT/SCOTT_SIMS_Aug2025/20*SCA ' + str(SCA_i) + " " + str(xlabel) + " " + str(ylabel) + ".fits")
            print(NDI_base_file)
            NDI_fits = fits.open(NDI_base_file[0])
            
            fig, axs = plt.subplots(1, 1, figsize=(15, 10)) 
            im = axs.imshow(NDI_fits[0].data, vmin=0, vmax=1E-2, origin='lower', cmap="RdYlBu_r")
            plt.colorbar(im, label='NDI')
            plt.text(100, 1000, "SCA " + str(SCA_i) + " - X:" + str(xlabel) + " Y:" + str(ylabel), color="yellow", fontsize=20) 
            plt.tight_layout()
            outname =  "level2_snap" + str(counter).zfill(5) + "_" + os.path.basename(NDI_base_file[0]).replace(".fits", ".png")
            plt.savefig(outname, dpi=100)
            output_list.append(outname)
            #  plt.show()
            counter = counter + 1
    os.system("magick -delay 100 -loop 0 level2_snap*.png deg20_ndi_movie_aug2025.mpeg")

if False:
    # Level 3 
    ndi_prefix = "lvl3"
    output_list = []
    counter = 0

    for SCA_i in tqdm(SCA_list):
        # print(SCA_i)
        for label_index in range(len(subarray_locations_db["xlabel"])):

            # Get the labels    
            xlabel = subarray_locations_db["xlabel"][label_index]
            ylabel = subarray_locations_db["ylabel"][label_index]

            # Find the right original files
            NDI_base_file = glob.glob('/Users/aborlaff/NASA/ROSALIA_DEPOT/SCOTT_SIMS_Aug2025/All*SCA ' + str(SCA_i) + " " + str(xlabel) + " " + str(ylabel) + ".fits")
            print(NDI_base_file)
            NDI_fits = fits.open(NDI_base_file[0])
            
            fig, axs = plt.subplots(1, 1, figsize=(15, 10)) 
            im = axs.imshow(NDI_fits[0].data, vmin=0, vmax=1E-3, origin='lower', cmap="RdYlBu_r")
            plt.colorbar(im, label='NDI')
            plt.text(50, 600, "SCA " + str(SCA_i) + " - X:" + str(xlabel) + " Y:" + str(ylabel), color="yellow", fontsize=20) 
            plt.tight_layout()
            outname =  "level3_snap" + str(counter).zfill(5) + "_" + os.path.basename(NDI_base_file[0]).replace(".fits", ".png")
            plt.savefig(outname, dpi=100)
            output_list.append(outname)
            #  plt.show()
            counter = counter + 1
    os.system("magick -delay 25 -loop 0 level3_snap*.png allsky_ndi_movie_aug2025.mpeg")



# Make a smooth background from the first level optics NDI (without the FPA lightshield!). 

"""
csmooth "2-deg at 5 arcsec SCA 1.fits" outfile="NDI_back_2deg_lvl1_SCA_1.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 2.fits" outfile="NDI_back_2deg_lvl1_SCA_2.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 3.fits" outfile="NDI_back_2deg_lvl1_SCA_3.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 4.fits" outfile="NDI_back_2deg_lvl1_SCA_4.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 5.fits" outfile="NDI_back_2deg_lvl1_SCA_5.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 6.fits" outfile="NDI_back_2deg_lvl1_SCA_6.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 7.fits" outfile="NDI_back_2deg_lvl1_SCA_7.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 8.fits" outfile="NDI_back_2deg_lvl1_SCA_8.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 9.fits" outfile="NDI_back_2deg_lvl1_SCA_9.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 10.fits" outfile="NDI_back_2deg_lvl1_SCA_10.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 11.fits" outfile="NDI_back_2deg_lvl1_SCA_11.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 12.fits" outfile="NDI_back_2deg_lvl1_SCA_12.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 13.fits" outfile="NDI_back_2deg_lvl1_SCA_13.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 14.fits" outfile="NDI_back_2deg_lvl1_SCA_14.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 15.fits" outfile="NDI_back_2deg_lvl1_SCA_15.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 16.fits" outfile="NDI_back_2deg_lvl1_SCA_16.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 17.fits" outfile="NDI_back_2deg_lvl1_SCA_17.fits" clobber=yes
csmooth "2-deg at 5 arcsec SCA 18.fits" outfile="NDI_back_2deg_lvl1_SCA_18.fits" clobber=yes
"""
############################################
############################################
############################################
############################################
# Make the level 1 NDI files. 
############################################
############################################
############################################
############################################

"""
# Level 1
# Step 1 - Open the first open (FPA holders)
# Step 2 - Open the additional layer for the lightshields. 

"""
subarray_locations_db = rs.roman.get_subarray_locations(SCA=1, verbose=False)
SCA_list = np.linspace(1,18,18, dtype="int")

ndi_prefix = "lvl1"
output_list = []

for SCA_i in SCA_list:
    # print(SCA_i)
    for label_index in range(len(subarray_locations_db["xlabel"])):

        # Get the labels 
        xlabel = subarray_locations_db["xlabel"][label_index]
        ylabel = subarray_locations_db["ylabel"][label_index]

        # Find the right original files
        NDI_base_file = glob.glob("/Users/aborlaff/NASA/ROSALIA_DEPOT/SCOTT_SIMS_Jan2025_FPA_base/*SCA " + str(SCA_i) + " " + str(xlabel) + " " + str(ylabel) + ".fits")
        NDI_lightshield_extra_file = glob.glob("/Users/aborlaff/NASA/ROSALIA_DEPOT/SCOTT_SIMS_Jan2025_light_shields/*SCA " + str(SCA_i) + " " + str(xlabel) + " " + str(ylabel) + ".fits")
        NDI_background = glob.glob("/Users/aborlaff/NASA/ROSALIA_DEPOT/SCOTT_SIMS_Jan2025_1st_order/2-deg_at_5_arcsec/NDI_back*SCA_" + str(SCA_i)+ ".fits")
        
        print(NDI_base_file)
        print(NDI_lightshield_extra_file)
        print(NDI_background)

        
        NDI_fits = fits.open(NDI_base_file[0])
        NDI_lightshield_fits = fits.open(NDI_lightshield_extra_file[0])
        NDI_background_fits = fits.open(NDI_background[0])

        NDI_fits[0].data = NDI_fits[0].data + NDI_lightshield_fits[0].data # + NDI_background_fits[0].data
        
        new_name = "/Users/aborlaff/NASA/ROSALIA_DEPOT/ndi_lvl1/" + ndi_prefix + "_SCA_" + str(SCA_i)+ "_SUB_X" + str(xlabel) + "_Y" + str(ylabel) + ".fits"
        # cmd = "cp '" + NDI_file[0] + "' " + new_name
        print(new_name)


        #### Now, make a custom WCS for this file, based on the angular size that they spand. 
        cdelt = [(NDI_fits[0].header["AAXISMAX"] - NDI_fits[0].header["AAXISMIN"])/NDI_fits[0].header["NAXIS1"], 
                  (NDI_fits[0].header["BAXISMAX"] - NDI_fits[0].header["BAXISMIN"])/NDI_fits[0].header["NAXIS2"]]

        crpix = [720.5, 720.5] # [(NDI_fits[0].header["NAXIS1"])/2, (NDI_fits[0].header["NAXIS2"])/2]
        SCA_subarray_location = pd.read_csv("/Users/aborlaff/NASA/ROSALIA/rosalia/CORE/detector/RST_WFI_SCA_Subarray_locations.csv")
        crval_X = np.array(SCA_subarray_location["Field X"])[(SCA_subarray_location["SCA"] == SCA_i) & (SCA_subarray_location["SCAx"] == xlabel) & (SCA_subarray_location["SCAy"] == ylabel)][0]
        crval_Y = np.array(SCA_subarray_location["Field Y"])[(SCA_subarray_location["SCA"] == SCA_i) & (SCA_subarray_location["SCAx"] == xlabel) & (SCA_subarray_location["SCAy"] == ylabel)][0]
 
        """
        subpos1, subpos2 = find_subposition_EXTNAME(EXTNAME=NDI_fits_cube[i].header["EXTNAME"]) 
        SCA = find_SCA_EXTNAME(EXTNAME=NDI_fits_cube[i].header["EXTNAME"])
        output_name = NDI_fits_cube_name.replace(".fits", "_SUB_X" + str(subpos1).zfill(2) + "_Y" + str(subpos2).zfill(2) + ".fits").replace(" ","_")

        NDI_i_data = NDI_fits_cube[i].data 
        NDI_i_header = NDI_fits_cube[i].header
        crpix = [(NDI_fits_cube[0].header["NAXIS1"])/2, (NDI_fits_cube[0].header["NAXIS2"])/2]
        # Here, we have to use the subarray location table to put CRVAL to the real location in degrees.
        SCA_subarray_location = pd.read_csv("/Users/aborlaff/NASA/STRAYCOR/straycor/CORE/detector/RST_WFI_SCA_Subarray_locations.csv")
        crval_X = np.array(SCA_subarray_location["Field X"])[(SCA_subarray_location["SCA"] == SCA) & (SCA_subarray_location["SCAx"] == subpos1) & (SCA_subarray_location["SCAy"] == subpos2)][0]
        crval_Y = np.array(SCA_subarray_location["Field Y"])[(SCA_subarray_location["SCA"] == SCA) & (SCA_subarray_location["SCAx"] == subpos1) & (SCA_subarray_location["SCAy"] == subpos2)][0]

        """
        crval = [0, 0] # [crval_X, -crval_Y] # I think that this has to be 0, 0 
        cdelt = [cdelt[0], cdelt[1]]
        header = rs.utils.create_custom_wcs(crpix=crpix, crval=crval, cdelt=cdelt, projection="CAR")

        header["AAXISMAX"] = (NDI_fits[0].header["AAXISMAX"])
        header["AAXISMIN"] = (NDI_fits[0].header["AAXISMIN"])
        header["BAXISMAX"] = (NDI_fits[0].header["BAXISMAX"])
        header["BAXISMIN"] = (NDI_fits[0].header["BAXISMIN"])
        header["SUBPOSX"] = xlabel
        header["SUBPOSY"] = ylabel
        header["EXTNAME"] = NDI_fits[0].header["EXTNAME"]

        #print(NDI_i_header)
        rs.utils.save_fits(array=NDI_fits[0].data, name=new_name, header=header)
        output_list.append(new_name)
        

# Reproject them to a TAN projection map. 
rs.utils.create_dummy_image_with_wcs(ra_cen=0, dec_cen=0, ra_size=2, dec_size=2, pixscale=5/60/60, outname="default_dummy.fits")

os.system("rm /Users/aborlaff/NASA/ROSALIA_DEPOT/ndi_lvl1/lvl1*_TAN.fits")
lvl1_maps = glob.glob("/Users/aborlaff/NASA/ROSALIA_DEPOT/ndi_lvl1/lvl1_*.fits")

from reproject import reproject_interp
from tqdm import tqdm 

footprint_dummy = fits.open("default_dummy.fits")

for lvl1_map_i in tqdm(lvl1_maps): 
    lvl1_map_i_fits = fits.open(lvl1_map_i)
    array, footprint = reproject_interp(lvl1_map_i_fits[0], footprint_dummy[0].header)
    #footprint_dummy[0].data = array 
    lvl1_map_i_fits[0].data = array
    lvl1_map_i_fits[0].header = footprint_dummy[0].header
    lvl1_map_i_fits.verify("silentfix")
    lvl1_map_i_fits.writeto(lvl1_map_i.replace(".fits", "_TAN.fits"), overwrite=True)



############################################
############################################
############################################
############################################
# Make the level 2 NDI files. 
############################################
############################################
############################################
############################################

"""
# Level 2
# Step 1 - Open the first open (FPA holders)
# Step 2 - Open the additional layer for the lightshields. 

"""
subarray_locations_db = rs.roman.get_subarray_locations(SCA=1, verbose=False)
SCA_list = np.linspace(1,18,18, dtype="int")

ndi_prefix = "lvl2"
output_list = []

for SCA_i in SCA_list:
    # print(SCA_i)
    for label_index in range(len(subarray_locations_db["xlabel"])):

        # Get the labels 
        xlabel = subarray_locations_db["xlabel"][label_index]
        ylabel = subarray_locations_db["ylabel"][label_index]

        # Find the right original files
        NDI_base_file = glob.glob('/Users/aborlaff/NASA/ROSALIA_DEPOT/SCOTT_SIMS_Aug2025/20-deg*SCA ' + str(SCA_i) + " " + str(xlabel) + " " + str(ylabel) + ".fits")
        # NDI_base_file = glob.glob('/Users/aborlaff/NASA/ROSALIA_DEPOT/SCOTT_SIMS_Jan2025_1st_order/IOA scatter 8x8 1e9rays 2025-06-09/*SCA ' + str(SCA_i) + " " + str(xlabel) + " " + str(ylabel) + ".fits")
        #NDI_lightshield_extra_file = glob.glob("/Users/aborlaff/NASA/ROSALIA_DEPOT/SCOTT_SIMS_Jan2025_light_shields/*SCA " + str(SCA_i) + " " + str(xlabel) + " " + str(ylabel) + ".fits")
        #NDI_background = glob.glob("/Users/aborlaff/NASA/ROSALIA_DEPOT/SCOTT_SIMS_Jan2025_1st_order/2-deg_at_5_arcsec/NDI_back*SCA_" + str(SCA_i)+ ".fits")
        
        print(NDI_base_file)
        
        NDI_fits = fits.open(NDI_base_file[0])

        NDI_fits[0].data = NDI_fits[0].data
        
        new_name = "/Users/aborlaff/NASA/ROSALIA_DEPOT/ndi_lvl2/" + ndi_prefix + "_SCA_" + str(SCA_i)+ "_SUB_X" + str(xlabel) + "_Y" + str(ylabel) + ".fits"
        # cmd = "cp '" + NDI_file[0] + "' " + new_name
        print(new_name)


        #### Now, make a custom WCS for this file, based on the angular size that they spand. 
        cdelt = [(NDI_fits[0].header["AAXISMAX"] - NDI_fits[0].header["AAXISMIN"])/NDI_fits[0].header["NAXIS1"], 
                  (NDI_fits[0].header["BAXISMAX"] - NDI_fits[0].header["BAXISMIN"])/NDI_fits[0].header["NAXIS2"]]

        crpix = [(NDI_fits[0].header["NAXIS1"])/2, (NDI_fits[0].header["NAXIS2"])/2]
        SCA_subarray_location = pd.read_csv("/Users/aborlaff/NASA/ROSALIA/rosalia/CORE/detector/RST_WFI_SCA_Subarray_locations.csv")
        crval_X = np.array(SCA_subarray_location["Field X"])[(SCA_subarray_location["SCA"] == SCA_i) & (SCA_subarray_location["SCAx"] == xlabel) & (SCA_subarray_location["SCAy"] == ylabel)][0]
        crval_Y = np.array(SCA_subarray_location["Field Y"])[(SCA_subarray_location["SCA"] == SCA_i) & (SCA_subarray_location["SCAx"] == xlabel) & (SCA_subarray_location["SCAy"] == ylabel)][0]
 
        """
        subpos1, subpos2 = find_subposition_EXTNAME(EXTNAME=NDI_fits_cube[i].header["EXTNAME"]) 
        SCA = find_SCA_EXTNAME(EXTNAME=NDI_fits_cube[i].header["EXTNAME"])
        output_name = NDI_fits_cube_name.replace(".fits", "_SUB_X" + str(subpos1).zfill(2) + "_Y" + str(subpos2).zfill(2) + ".fits").replace(" ","_")

        NDI_i_data = NDI_fits_cube[i].data 
        NDI_i_header = NDI_fits_cube[i].header
        crpix = [(NDI_fits_cube[0].header["NAXIS1"])/2, (NDI_fits_cube[0].header["NAXIS2"])/2]
        # Here, we have to use the subarray location table to put CRVAL to the real location in degrees.
        SCA_subarray_location = pd.read_csv("/Users/aborlaff/NASA/STRAYCOR/straycor/CORE/detector/RST_WFI_SCA_Subarray_locations.csv")
        crval_X = np.array(SCA_subarray_location["Field X"])[(SCA_subarray_location["SCA"] == SCA) & (SCA_subarray_location["SCAx"] == subpos1) & (SCA_subarray_location["SCAy"] == subpos2)][0]
        crval_Y = np.array(SCA_subarray_location["Field Y"])[(SCA_subarray_location["SCA"] == SCA) & (SCA_subarray_location["SCAx"] == subpos1) & (SCA_subarray_location["SCAy"] == subpos2)][0]

        """
        crval = [0, 0] # [crval_X, -crval_Y] # I think that this has to be 0, 0 
        cdelt = [cdelt[0], cdelt[1]        ]
        header = rs.utils.create_custom_wcs(crpix=crpix, crval=crval, cdelt=cdelt, projection="CAR")

        header["AAXISMAX"] = (NDI_fits[0].header["AAXISMAX"])
        header["AAXISMIN"] = (NDI_fits[0].header["AAXISMIN"])
        header["BAXISMAX"] = (NDI_fits[0].header["BAXISMAX"])
        header["BAXISMIN"] = (NDI_fits[0].header["BAXISMIN"])
        header["SUBPOSX"] = xlabel
        header["SUBPOSY"] = ylabel
        header["EXTNAME"] = NDI_fits[0].header["EXTNAME"]

        #print(NDI_i_header)
        rs.utils.save_fits(array=NDI_fits[0].data, name=new_name, header=header)
        output_list.append(new_name)
        
        #NDI_fits.verify("silentfix")
        #NDI_fits.writeto(new_name, overwrite=True)

# Reproject them to a TAN projection map. 
rs.utils.create_dummy_image_with_wcs(ra_cen=0, dec_cen=0, ra_size=21, dec_size=21, pixscale=1/60, outname="default_dummy.fits")

os.system("rm /Users/aborlaff/NASA/ROSALIA_DEPOT/ndi_lvl2/lvl2*_TAN.fits")
lvl2_maps = glob.glob("/Users/aborlaff/NASA/ROSALIA_DEPOT/ndi_lvl2/lvl2_*[0-9].fits")

from reproject import reproject_interp
from tqdm import tqdm 

footprint_dummy = fits.open("default_dummy.fits")

for lvl2_map_i in tqdm(lvl2_maps): 
    lvl2_map_i_fits = fits.open(lvl2_map_i)
    array, footprint = reproject_interp(lvl2_map_i_fits[0], footprint_dummy[0].header)
    #footprint_dummy[0].data = array 
    lvl2_map_i_fits[0].data = array
    lvl2_map_i_fits[0].header = footprint_dummy[0].header
    lvl2_map_i_fits.verify("silentfix")
    lvl2_map_i_fits.writeto(lvl2_map_i.replace(".fits", "_TAN.fits"), overwrite=True)







############################################
############################################
############################################
############################################
# Make the level 3 NDI files. 
############################################
############################################
############################################
############################################

"""
# Level 3
# Step 1 - Open the first open (FPA holders)
# Step 2 - Open the additional layer for the lightshields. 

"""


subarray_locations_db = rs.roman.get_subarray_locations(SCA=1, verbose=False)
SCA_list = np.linspace(1,18,18, dtype="int")

ndi_prefix = "lvl3"
output_list = []

for SCA_i in SCA_list:
    # print(SCA_i)
    for label_index in range(len(subarray_locations_db["xlabel"])):

        # Get the labels 
        xlabel = subarray_locations_db["xlabel"][label_index]
        ylabel = subarray_locations_db["ylabel"][label_index]

        # Find the right original files
        NDI_base_file = glob.glob('/Users/aborlaff/NASA/ROSALIA_DEPOT/SCOTT_SIMS_Aug2025/All*SCA ' + str(SCA_i) + " " + str(xlabel) + " " + str(ylabel) + ".fits")
        #NDI_base_file = glob.glob('/Users/aborlaff/NASA/ROSALIA_DEPOT/SCOTT_SIMS_Aug2025/All*SCA ' + str(SCA_i) + " " + str(xlabel) + " " + str(ylabel) + ".fits")

        print(NDI_base_file)
        
        NDI_fits = fits.open(NDI_base_file[0])
        
        #new_name = "/Users/aborlaff/NASA/ROSALIA_DEPOT/ndi_lvl3/" + ndi_prefix + "_SCA_" + str(SCA_i)+ ".fits"
        new_name = "/Users/aborlaff/NASA/ROSALIA_DEPOT/ndi_lvl3/" + ndi_prefix + "_SCA_" + str(SCA_i)+ "_SUB_X" + str(xlabel) + "_Y" + str(ylabel) + ".fits"

        print(new_name)


        #### Now, make a custom WCS for this file, based on the angular size that they spand. 
        cdelt = [(NDI_fits[0].header["AAXISMAX"] - NDI_fits[0].header["AAXISMIN"])/NDI_fits[0].header["NAXIS1"], 
                  (NDI_fits[0].header["BAXISMAX"] - NDI_fits[0].header["BAXISMIN"])/NDI_fits[0].header["NAXIS2"]]

        crpix = [(NDI_fits[0].header["NAXIS1"])/2, (NDI_fits[0].header["NAXIS2"])/2]
        SCA_subarray_location = pd.read_csv("/Users/aborlaff/NASA/ROSALIA/rosalia/CORE/detector/RST_WFI_SCA_Subarray_locations.csv")
        crval_X = np.mean(np.array(SCA_subarray_location["Field X"])[(SCA_subarray_location["SCA"] == SCA_i)])
        crval_Y = np.mean(np.array(SCA_subarray_location["Field Y"])[(SCA_subarray_location["SCA"] == SCA_i)])

        print(crval_X)
        print(crval_Y)
        """
        subpos1, subpos2 = find_subposition_EXTNAME(EXTNAME=NDI_fits_cube[i].header["EXTNAME"]) 
        SCA = find_SCA_EXTNAME(EXTNAME=NDI_fits_cube[i].header["EXTNAME"])
        output_name = NDI_fits_cube_name.replace(".fits", "_SUB_X" + str(subpos1).zfill(2) + "_Y" + str(subpos2).zfill(2) + ".fits").replace(" ","_")



        

        NDI_i_data = NDI_fits_cube[i].data 
        NDI_i_header = NDI_fits_cube[i].header
        crpix = [(NDI_fits_cube[0].header["NAXIS1"])/2, (NDI_fits_cube[0].header["NAXIS2"])/2]
        # Here, we have to use the subarray location table to put CRVAL to the real location in degrees.
        SCA_subarray_location = pd.read_csv("/Users/aborlaff/NASA/STRAYCOR/straycor/CORE/detector/RST_WFI_SCA_Subarray_locations.csv")
        crval_X = np.array(SCA_subarray_location["Field X"])[(SCA_subarray_location["SCA"] == SCA) & (SCA_subarray_location["SCAx"] == subpos1) & (SCA_subarray_location["SCAy"] == subpos2)][0]
        crval_Y = np.array(SCA_subarray_location["Field Y"])[(SCA_subarray_location["SCA"] == SCA) & (SCA_subarray_location["SCAx"] == subpos1) & (SCA_subarray_location["SCAy"] == subpos2)][0]

        """
        # crval = [crval_X, -crval_Y]
        crval = [0, 0] # [crval_X, -crval_Y] # I think that this has to be 0, 0 
        cdelt = [cdelt[0], cdelt[1]        ]
        header = rs.utils.create_custom_wcs(crpix=crpix, crval=crval, cdelt=cdelt, projection="CAR")

        header["AAXISMAX"] = (NDI_fits[0].header["AAXISMAX"])
        header["AAXISMIN"] = (NDI_fits[0].header["AAXISMIN"])
        header["BAXISMAX"] = (NDI_fits[0].header["BAXISMAX"])
        header["BAXISMIN"] = (NDI_fits[0].header["BAXISMIN"])
        #header["SUBPOSX"] = xlabel
        #header["SUBPOSY"] = ylabel
        header["EXTNAME"] = NDI_fits[0].header["EXTNAME"]

        #print(NDI_i_header)
        rs.utils.save_fits(array=NDI_fits[0].data, name=new_name, header=header)
        output_list.append(new_name)
        
        #NDI_fits.verify("silentfix")
        #NDI_fits.writeto(new_name, overwrite=True)
    

# Level 3 - Clean the nans 
# Reproject them to a TAN projection map. 
os.system("rm /Users/aborlaff/NASA/ROSALIA_DEPOT/ndi_lvl3/lvl3*_HP.fits")
#os.system("rm /Users/aborlaff/NASA/ROSALIA_DEPOT/ndi_lvl3/lvl3*.fits")
os.system("rm /Users/aborlaff/NASA/ROSALIA_DEPOT/ndi_lvl3/lvl3*_sclmap.fits")
lvl3_maps = glob.glob("/Users/aborlaff/NASA/ROSALIA_DEPOT/ndi_lvl3/lvl3_*[0-9].fits")

from reproject import reproject_interp
from tqdm import tqdm 
import reproject
import healpy as hp

for lvl3_map_i in tqdm(lvl3_maps): 
    lvl3_map_fits = fits.open(lvl3_map_i)
    lvl3_map_fits[0].data[np.isnan(lvl3_map_fits[0].data)] = 0

    import matplotlib.pyplot as plt
    #from astropy.convolution import Gaussian2DKernel
    #from astropy.convolution import convolve, convolve_fft

    #gaussian_2D_kernel = Gaussian2DKernel(5)
    #lvl3_map_fits[0].data = convolve(lvl3_map_fits[0].data, gaussian_2D_kernel)
    
    lvl3_map_fits.verify("silentfix")
    lvl3_map_fits.writeto(lvl3_map_i, overwrite=True)
    

for lvl3_map_i in tqdm(lvl3_maps): 
    m = reproject.reproject_to_healpix(input_data=lvl3_map_i, coord_system_out='icrs', hdu_in=0, order='bilinear', nested=False, nside=256)
    m = np.array(m)
    m[np.isnan(m)] = 0
    hp.fitsfunc.write_map(filename=lvl3_map_i.replace(".fits", "_HP.fits"), m=m, nest=False, overwrite=True)




# Initialize a minimal storage array for the straylight. 
# The minimal storage array has 64 x 18, which is 64 locations, for each one of the 18 SCAs.
n_subarrays = len(subarray_locations_db["NDI_labels"])


# Pregenerate the interpolators. 
from scipy.interpolate import RegularGridInterpolator
import pickle
import rosalia as rs 
from tqdm import tqdm
from astropy.io import fits
import astropy.wcs as astropy_wcs
import numpy as np
import os 

for level in ["lvl1", "lvl2"]:
    for SCA in rs.telescopes.Roman.WFI_SCAs:
        headers = []
        data = []
        subarray_locations_db = rs.roman.get_subarray_locations(SCA=SCA, verbose=False)
        NDI_labels = subarray_locations_db["NDI_labels"]

        NDI_extension = "_TAN.fits"

        w_list = []
        f_list = [] 

        for i in tqdm(range(len(NDI_labels))):
            ndi_lvl_tag = f"/Users/aborlaff/NASA/ROSALIA_DEPOT/ndi_{level}/{level}_SCA_"+str(SCA)+"_SUB_" + NDI_labels[i] + NDI_extension
            #print(ndi_lvl_tag)
            ndi_fits = fits.open(ndi_lvl_tag)
            w = astropy_wcs.WCS(header=ndi_fits[0].header, fobj=ndi_fits, naxis=2)
            x_grid = np.linspace(0, ndi_fits[0].header["NAXIS1"]-1, ndi_fits[0].header["NAXIS1"])
            y_grid = np.linspace(0, ndi_fits[0].header["NAXIS2"]-1, ndi_fits[0].header["NAXIS2"])
            # f = RegularGridInterpolator((x_grid, y_grid), np.flip(ndi_fits[0].data, axis=1).T, bounds_error=False, fill_value=0)
            # f = RegularGridInterpolator((x_grid, y_grid), ndi_fits[0].data, bounds_error=False, fill_value=0)
            f = RegularGridInterpolator((x_grid, y_grid), ndi_fits[0].data.T, bounds_error=False, fill_value=0)
            w_list.append(w)
            f_list.append(f)

        data = {"NDI_labels": NDI_labels, "wcs": w_list, "ndi_interpolator": f_list}

        with open(os.path.dirname(ndi_lvl_tag)+"/"+level+"_SCA_"+str(SCA)+".pkl", "wb") as f:
            pickle.dump(data, f)


os.system("cp /Users/aborlaff/NASA/ROSALIA_DEPOT/ndi_lvl1/lvl1_SCA_*.pkl /Users/aborlaff/NASA/rosalia_cache/CORE/NDI/RST/ndi_lvl1/")
os.system("cp /Users/aborlaff/NASA/ROSALIA_DEPOT/ndi_lvl2/lvl2_SCA_*.pkl /Users/aborlaff/NASA/rosalia_cache/CORE/NDI/RST/ndi_lvl2/")


