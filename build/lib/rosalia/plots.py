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
import astropy.wcs as wcs
import matplotlib.pyplot as plt
from celluloid import Camera
from astropy.io import ascii
import astropy.units as u
from astropy.coordinates import ICRS, Angle, SkyCoord
import rosalia as rs

# Suppress warnings. Comment this out if you wish to see the warning messages
import warnings
warnings.filterwarnings('ignore')

psf_archive = os.path.dirname(rs.utils.__file__) + "/../PSF_ARCHIVE/"


from itertools import cycle
from shutil import get_terminal_size
from threading import Thread
from time import sleep


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
