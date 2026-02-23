# Alejandro S. Borlaff. NASA Ames Research Center. a.s.borlaff@nasa.gov / asborlaff@gmail.com
# January 20, 2023.
#
# STRAYCOR/PSF module
# This module will hold all the programs related to the modelling and removal
# of the PSF.
#
# Version log:
# v.1.0 - 20 Enero 2023. First loading of programs inherited from former monolithic straycor.py
#
##########################################################

############################
import os
import sys
import pandas as pd
import numpy as np
import bottleneck as bn
from tqdm import tqdm
from astropy.io import fits
import astropy.wcs as astropy_wcs
import matplotlib.pyplot as plt
import matplotlib.colors as matplotlib_colors
from celluloid import Camera
from astropy.io import ascii
import astropy.units as u
from astropy.coordinates import ICRS, Angle, SkyCoord
import rosalia as rs

# Suppress warnings. Comment this out if you wish to see the warning messages
import warnings
warnings.filterwarnings('ignore')

psf_archive = os.environ["ROSALIACACHE"] + "/PSF_ARCHIVE/"


from itertools import cycle
from shutil import get_terminal_size
from threading import Thread
from time import sleep


# Class of different styles
class style():
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    UNDERLINE = '\033[4m'
    RESET = '\033[0m'


class Loader:
    def __init__(self, desc="Loading...", end="Done!", timeout=0.1):
        """
        A loader-like context manager

        Args:
            desc (str, optional): The loader's description. Defaults to "Loading...".
            end (str, optional): Final print. Defaults to "Done!".
            timeout (float, optional): Sleep time between prints. Defaults to 0.1.
        """
        self.desc = desc
        self.end = end
        self.timeout = timeout

        self._thread = Thread(target=self._animate, daemon=True)
        self.steps = ["⢿", "⣻", "⣽", "⣾", "⣷", "⣯", "⣟", "⡿"]
        self.done = False

    def start(self):
        self._thread.start()
        return self

    def _animate(self):
        for c in cycle(self.steps):
            if self.done:
                break
            print(f"\r{self.desc} {c}", flush=True, end="")
            sleep(self.timeout)

    def __enter__(self):
        self.start()

    def stop(self):
        self.done = True
        cols = get_terminal_size((80, 20)).columns
        print("\r" + " " * cols, end="", flush=True)
        print(f"\r{self.end}", flush=True)

    def __exit__(self, exc_type, exc_value, tb):
        # handle exceptions with those variables ^
        self.stop()
################################################################################

def plot_stars_around(catalog, max_plot_size=100, min_plot_size=5, alpha=0.2):
    # Auxiliary function to merge some plotting steps in the stray-ligth mapping

    flux_for_plot = 10**(0.4*(8.9-np.array(catalog["mag_lambda"])))
    ra_stars = np.array(catalog["ra"])
    dec_stars = np.array(catalog["dec"])
    # Plot the planets separately #
    SSO_list_plain = ["Sun", "Moon", "Earth", "Mercury",
                      "Venus", "Mars", "Jupiter", "Saturn",
                      "Uranus", "Neptune", "Pluto"]
    SSO_list_symbols = [u'$\u2609$', u'$\u263D$', u'$\u1F728$',  u'$\u263F$',
                        u'$\u2640$', u'$\u2642$',u'$\u2643$', u'$\u2644$',
                        u'$\u26E2$', u'$\u2646$',  u'$\u2647$']

    print(bn.nanmin(ra_stars))
    print(bn.nanmax(ra_stars))
    print(bn.nanmin(dec_stars))
    print(bn.nanmax(dec_stars))

    for i, SSO_list_plain_i in zip(range(len(SSO_list_plain)), SSO_list_plain):

        try:
            planet_i = np.where(catalog["source_id"]==SSO_list_plain_i)[0][0]
            #plt.text(ra_stars[planet_i], dec_stars[planet_i], s=SSO_list_symbols[i])
            plt.scatter(ra_stars[planet_i], dec_stars[planet_i], marker=SSO_list_symbols[i], s=1000)
            flux_for_plot[planet_i] = np.nan
            ra_stars[planet_i] = np.nan
            dec_stars[planet_i] = np.nan

        except:
            print("SSO body: " + SSO_list_plain_i + " not found. Skipping")

    max_flux = np.nanpercentile(flux_for_plot, 99)
    min_flux = np.nanpercentile(flux_for_plot, 1)

    # Idea: Plot the main stars with a different color.
    plot_size = flux_for_plot*(max_plot_size-min_plot_size)/(max_flux-min_flux)+0.1
    plt.scatter(ra_stars, dec_stars, marker="o", facecolor="grey", edgecolor="black", alpha=alpha, s=plot_size)
    plt.xlabel("RA (ICRS)")
    plt.ylabel("DEC (ICRS)")
    return(plot_size)

################################################################################

def plot_stray_and_ndi_map(stray_name, star_catalog, ra_point, dec_point, pa_point, ndi_name, stray_scale_factor, vmin, vmax, ndi_linthresh, marker="*", edgecolor="black", color="red", s=100):
# if True:
    # Axis 1
    # stray_name = stray_image_list[i]

    # Open the Stray-light map
    stray_fits = fits.open(stray_name)
    stray_list_image = stray_fits[1].data/stray_scale_factor

    # Transform to magnitudes per arcsec2
    stray_list_image = rs.detectors.fe2mu(fe=stray_list_image, instrument="WFI", filter_name="F158", telescope="Roman")

    fig, axs = plt.subplots(1, 2, figsize=(14, 8), width_ratios=[1.2, 1])

    w = astropy_wcs.WCS(header=stray_fits[0].header, fobj=stray_fits, naxis=2)

    # Axis 1 - Showing the Stray-light map
    stray_list_image[stray_list_image == 0] = np.nan
    im = axs[0].imshow(stray_list_image, origin='lower', cmap="inferno_r", vmin=vmin, vmax=vmax)
    cbar = plt.colorbar(im, ax=axs[0], location='top', )
    cbar.set_label(label='Surface brightness (mag arcsec$^{-2}$)',size=15,weight='bold')


    # Axis 2 - Showing the NDI map
    axs[1] = plt.subplot(122)
    #ndi_name = ndi_image_list[i]
    #ndi_header = fits.open(ndi_name)[0].header


    # Make the modified header based on the exposure properties.
    ndi_fits = fits.open(ndi_name)
    w = astropy_wcs.WCS(header=ndi_fits[0].header, fobj=ndi_fits, naxis=2)
    w.wcs.crota = -pa_point,-pa_point
    w.wcs.crval = ra_point,dec_point
    w.wcs.cdelt = -ndi_fits[0].header["CDELT1"],ndi_fits[0].header["CDELT2"]
    ndi_fits[0].header = w.to_header()

    # ndi_header = fits.open(ndi_headerlet)[0].header
    im = axs[1].imshow(np.flip(ndi_fits[0].data, axis=1), origin='lower', cmap="RdYlBu_r",
                           norm=matplotlib_colors.SymLogNorm(linthresh=ndi_linthresh, linscale=1,
                                            base=10))

    ndi_w = astropy_wcs.WCS(header=ndi_fits[0].header, fobj=ndi_fits, naxis=2)
    # star_catalog = star_catalog_list[i]
    x_stars, y_stars = ndi_w.wcs_world2pix(star_catalog["ra"], star_catalog["dec"], 0)

    axs[1].scatter(x_stars, y_stars, marker=marker, edgecolor=edgecolor,
                   facecolor=color, alpha=0.9, s=s)

    axs[1].set_title("Stray-light source location on Focal Plane\n", fontsize=16)

    max_extent_ndi = np.abs(ndi_fits[0].header["CDELT2"]) * ndi_fits[0].data.shape[0]/2# This is the radial extent of the NDI map
    print(max_extent_ndi)
    angle_ticks = np.round(np.linspace(-max_extent_ndi, max_extent_ndi, 5),1)
    print(angle_ticks)
    pixel_ticks = angle_ticks/np.abs(ndi_fits[0].header["CDELT2"])  + ndi_fits[0].data.shape[0]/2
    print(pixel_ticks)
    axs[1].set_xticks(pixel_ticks, angle_ticks, size=12)
    axs[1].set_yticks(pixel_ticks, angle_ticks, size=12)

    axs[0].set_xticks([0, 1340*0.5, 1340*1, 1340*1.5, 1340*2], [-0.4, -0.2, 0, 0.2, 0.4], size=12)
    axs[0].set_yticks([0, 860*0.5, 860*1, 860*1.5, 860*2], [-0.26, -0.13, 0, 0.13, 0.26], size=12)

    axs[0].set_xlabel("X (degrees)", size=15)
    axs[0].set_ylabel("Y (degrees)", size=15)
    axs[1].set_xlabel("X (degrees)", size=15)
    axs[1].set_ylabel("Y (degrees)", size=15)

    plt.tight_layout()
    plt.savefig(stray_name.replace(".fits", ".png"), dpi=300)
    plt.show()


import numpy as np

class ascii_progress_focal_plane():
    sca_mask = np.array(list(""+\
      " gggggggg                                     77777777 \n"+\
      " gggggggg dddddddd                   44444444 77777777 \n"+\
      " gggggggg dddddddd aaaaaaaa 11111111 44444444 77777777 \n"+\
      " gggggggg dddddddd aaaaaaaa 11111111 44444444 77777777 \n"+\
      " hhhhhhhh dddddddd aaaaaaaa 11111111 44444444 88888888 \n"+\
      " hhhhhhhh eeeeeeee aaaaaaaa 11111111 55555555 88888888 \n"+\
      " hhhhhhhh eeeeeeee bbbbbbbb 22222222 55555555 88888888 \n"+\
      " hhhhhhhh eeeeeeee bbbbbbbb 22222222 55555555 88888888 \n"+\
      " iiiiiiii eeeeeeee bbbbbbbb 22222222 55555555 99999999 \n"+\
      " iiiiiiii ffffffff bbbbbbbb 22222222 66666666 99999999 \n"+\
      " iiiiiiii ffffffff cccccccc 33333333 66666666 99999999 \n"+\
      " iiiiiiii ffffffff cccccccc 33333333 66666666 99999999 \n"+\
      "          ffffffff cccccccc 33333333 66666666          \n"+\
      "                   cccccccc 33333333                   \n"))


    x_mask = np.array(list(""+\
      " 01234567                                     01234567 \n"+\
      " 01234567 01234567                   01234567 01234567 \n"+\
      " 01234567 01234567 01234567 01234567 01234567 01234567 \n"+\
      " 01234567 01234567 01234567 01234567 01234567 01234567 \n"+\
      " 01234567 01234567 01234567 01234567 01234567 01234567 \n"+\
      " 01234567 01234567 01234567 01234567 01234567 01234567 \n"+\
      " 01234567 01234567 01234567 01234567 01234567 01234567 \n"+\
      " 01234567 01234567 01234567 01234567 01234567 01234567 \n"+\
      " 01234567 01234567 01234567 01234567 01234567 01234567 \n"+\
      " 01234567 01234567 01234567 01234567 01234567 01234567 \n"+\
      " 01234567 01234567 01234567 01234567 01234567 01234567 \n"+\
      " 01234567 01234567 01234567 01234567 01234567 01234567 \n"+\
      "          01234567 01234567 01234567 01234567          \n"+\
      "                   01234567 01234567                   \n"))

    y_mask = np.array(list(""+\
      " 33333333                                     33333333 \n"+\
      " 22222222 33333333                   33333333 22222222 \n"+\
      " 11111111 22222222 33333333 33333333 22222222 11111111 \n"+\
      " 00000000 11111111 22222222 22222222 11111111 00000000 \n"+\
      " 33333333 00000000 11111111 11111111 00000000 33333333 \n"+\
      " 22222222 33333333 00000000 00000000 33333333 22222222 \n"+\
      " 11111111 22222222 33333333 33333333 22222222 11111111 \n"+\
      " 00000000 11111111 22222222 22222222 11111111 00000000 \n"+\
      " 33333333 00000000 11111111 11111111 00000000 33333333 \n"+\
      " 22222222 33333333 00000000 00000000 33333333 22222222 \n"+\
      " 11111111 22222222 33333333 33333333 22222222 11111111 \n"+\
      " 00000000 11111111 22222222 22222222 11111111 00000000 \n"+\
      "          00000000 11111111 11111111 00000000          \n"+\
      "                   00000000 00000000                   \n"))



    canvas = np.array(list(""+\
      " ________                                     ________ \n"+\
      " ________ ________                   ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      "          ________ ________ ________ ________          \n"+\
      "                   ________ ________                   \n"))

    canvas_zero = np.array(list(""+\
      " ________                                     ________ \n"+\
      " ________ ________                   ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      " ________ ________ ________ ________ ________ ________ \n"+\
      "          ________ ________ ________ ________          \n"+\
      "                   ________ ________                   \n"))

    SCA_list = ["1","2","3","4","5","6","7","8","9","a","b","c","d","e","f","g","h","i"]
    x_range = ["0","1","2","3","4","5","6","7"]
    y_range = ["0","0","1","1","2","2","3","3"]




def print_ascii_focal_plane(x, y, SCA):
    subarray_NDI_labels = np.array([-17.92, -12.8, -7.68, -2.56, 2.56, 7.68, 12.8, 17.92])
    canvas = ascii_progress_focal_plane.canvas
    SCA_label = np.copy(ascii_progress_focal_plane.SCA_list)[SCA-1]
    x_label = np.copy(ascii_progress_focal_plane.x_range)[(subarray_NDI_labels==x)]
    y_label = np.copy(ascii_progress_focal_plane.y_range)[(subarray_NDI_labels==y)]

    sca_mask = ascii_progress_focal_plane.sca_mask
    x_mask = np.copy(ascii_progress_focal_plane.x_mask)
    y_mask = np.copy(ascii_progress_focal_plane.y_mask)

    filter_canvas = np.where((sca_mask == SCA_label) & (x_mask == x_label) & (y_mask == y_label))
    status = canvas[filter_canvas]
    if status == "_":
        char = "▄"
    if status == "▄":
        char = "█" # "█"
    canvas[filter_canvas] = char
    print(''.join(list(canvas)), flush=True)
    return()


def plot_rosalia_logo():
    ascii = np.array(list(style.YELLOW + "\n\n"+\
"   ___  ____  _______   __   _______         \n"+\
"  / _ \/ __ \/ __/ _ | / /  /  _/ _ | " + style.BLUE + " + " + style.YELLOW +"\n"+\
" / , _/ /_/ /\ \/ __ |/ /___/ // __ |__      \n"+\
"/_/|_|\____/___/_/ |_/____/___/_/ |____\     \n"+\
"ROman Sky Analyst for Low surface brightness Imaging & Astronomy \n"+\
style.RED + "NASA" + style.RESET + " Ames Research Center - " + style.RESET + " Roman Wide Field Science\n"+\
style.RESET + "Contact: Alex S. Borlaff (NASA Ames/STA) -  "+style.BLUE+"a.s.borlaff@nasa.gov \n\n\n"+style.RESET))
    print(''.join(list(ascii)), flush=True)
