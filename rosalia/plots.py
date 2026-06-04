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



def main_offender_find_fraction_of_map(mainoff_name, catalog_name):
    from astropy.io import fits
    import numpy as np
    import pandas as pd
    from astroquery.simbad import Simbad
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    hdu = fits.open(mainoff_name)
    data = hdu[0].data
    try:
        data[data == 0] = np.nan
    except:
        pass

    catalog = pd.read_csv(catalog_name)
        
    data_flat = data.flatten()
    data_flat = data_flat[~np.isnan(data_flat)] 

    unique_ids, unique_counts = np.unique(data_flat.astype("int64"), return_counts=True)
    print("Unique IDs:")
    print(unique_ids)
    unique_ids = unique_ids[unique_counts>1]
    print("Unique IDs with counts > 1:")
    print(unique_ids)
    fraction_by_offender = [] 

    import bottleneck as bn
    source_name_list = []
    ra_list = []
    dec_list = []
    mag_lambda_list = []
    unique_id_valid = []
    for i in range(len(unique_ids)): 
        unique_id = unique_ids[i]
        
        if np.where(catalog["cat_id"] == unique_id)[0].size == 0:
            print("Catalog does not have information for cat_id: " + str(unique_id) + ". Skipping")
            continue
        
        unique_id_valid.append(unique_id)
        # print(unique_id)
        pixels_with_id = bn.nansum(data == unique_id)
        #print("pixels_with_id: " + str(unique_id) + " - " + str(pixels_with_id))
        total_valid_pixels = len(data.flatten())-bn.nansum(np.isnan(data.flatten()))
        #print("total_valid_pixels: " + str(total_valid_pixels)) 
        fraction_by_offender.append(pixels_with_id/total_valid_pixels)
        #print(pixels_with_id/total_valid_pixels)
        #print(fraction_by_offender[i])

        ra = np.float32(catalog.iloc[catalog["cat_id"] == unique_id]["ra"])[0]
        dec = np.float32(catalog.iloc[catalog["cat_id"] == unique_id]["dec"])[0]
        mag_lambda = np.float32(catalog.iloc[catalog["cat_id"] == unique_id]["mag_lambda"])[0]


        # Try to find the name of the main_offender
        simbad_query = Simbad.query_region(SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs'), radius=1 * u.arcsec)
        if len(simbad_query)==0:
            print("Simbad cannot find the name of the offender at " + str(ra) + " - " + str(dec)) 
            source_name = "CAT_ID_" + str(unique_id)
        if len(simbad_query)==1:
            print(simbad_query["main_id"][0])
            source_name = simbad_query["main_id"][0]
        if len(simbad_query)>1:
            print("Simbad has multiple offending sources at  " + str(ra) + " - " + str(dec)) 
            source_name = simbad_query["main_id"][0] + "_" + str(unique_id)
            # source_name = 
            
        ra_list.append(ra)
        dec_list.append(dec)
        mag_lambda_list.append(mag_lambda)
        source_name_list.append(source_name)
    
    main_offender_db = pd.DataFrame({"source_id":unique_id_valid, "ra": ra_list, "dec": dec_list, "mag_lambda": mag_lambda_list,
                                     "source_name":source_name_list,
                                    "fraction_by_offender":fraction_by_offender})
    main_offender_db = main_offender_db.sort_values(by=["fraction_by_offender"], ascending=False)


    # How many stars are offenders? Show a max of 10. 
    n_offenders = len(main_offender_db)
    if n_offenders > 10:
        main_offender_db = main_offender_db[:10]
        
        
    # 2. Prepare your 10 lines of text
    lines = ['Main stray-light offenders:\n']
    for i in range(len(main_offender_db)):
        name = main_offender_db["source_name"].iloc[i]
        percen = 100*main_offender_db["fraction_by_offender"].iloc[i]
        ra = main_offender_db["ra"].iloc[i]
        dec = main_offender_db["dec"].iloc[i]
        mag_lambda = main_offender_db["mag_lambda"].iloc[i]
        source_id = main_offender_db["source_id"].iloc[i]
        lines.append("" + str("{:.2f}".format(percen)) + "% - ("  + str("{:.2f}".format(ra)) + ", " + str("{:.2f}".format(dec)) +  "), m="  + str("{:.2f}".format(mag_lambda)) + ' - ' + str(name))
    multiline_text = "\n".join(lines)
    print(multiline_text)


    return(main_offender_db, multiline_text)


def make_stray_plot(input_name, ext, mode="normal", catalog_name=None, 
                    vmin=None, vmax=None, 
                    color_label = 'Surface brightness (mag arcsec$^{-2}$)',
                    cmap="RdYlBu", output_name=None, figsize=(10,7), mu_vmin=None, mu_vmax=None):
    import matplotlib.pyplot as plt
    from astropy.utils.data import get_pkg_data_filename
    from astropy.wcs import WCS as astropy_wcs
    from astropy.io import fits
    import os
    import numpy as np
    import rosalia as rs
    if output_name is None:
        output_name = input_name.replace(".fits", "_" + mode + ".png")
    plt.style.use(os.path.dirname(rs.__file__) + "/style/nature_style.mplstyle")

    hdu = fits.open(input_name)
    data = hdu[ext].data
    wcs = astropy_wcs(hdu[ext].header)

    
    #if mode == "main_offender":
    #    data, header = rs.utils.reproject_roman_wfi_fits(input_name, input_ext=ext,
    #                                            reference_name=flt_name, 
    #                                            reference_ext=np.linspace(1,18,18, dtype="int"))
    
    fig, ax = plt.subplots(figsize=figsize, subplot_kw=dict(projection=wcs))

    if mode == "fe2mu":
        print("DEMO WARNING: make_stray_plot is Assuming F129.")
        data = rs.detectors.fe2mu(data, telescope="Roman", instrument="WFI", filter_name="F129")


    data[np.isinf(data)] = np.nan
    data[data == 0] = np.nan

    vmin = np.nanpercentile(data, 5)
    vmax = np.nanpercentile(data, 95)

    if mode == "fe2mu" and (mu_vmin is not None) and (mu_vmax is not None):
        vmin = mu_vmin
        vmax = mu_vmax

    print(vmin)
    print(vmax)
    im=ax.imshow(data, vmin=vmin, vmax=vmax, origin='lower', cmap=cmap)
    ax.set(xlabel='Right ascension (degrees)', ylabel='Declination (degrees)')
    cbar = plt.colorbar(im, ax=ax, location='right', fraction=0.046, pad=0.04)
    cbar.set_label(label=color_label,weight='bold')

    if mode == "main_offender":
        main_offender_db, multiline_text = main_offender_find_fraction_of_map(input_name, catalog_name)
        # 3. Add the text box
        # x, y coordinates are in data units by default
        ax.text(1.25, 0.95, multiline_text, 
            fontsize=9,
            color='black',
            transform=ax.transAxes,
            verticalalignment='top', 
            bbox=dict(facecolor='white', alpha=0.5))
    fig.tight_layout()
    plt.savefig(output_name, dpi=300)
    plt.show()
    return(output_name)


def make_stars_around_plot(flt_name, catalog_name, radius = 0.6, output_name=None, figsize=(10,7)):
    if output_name is None:
        output_name = flt_name.replace(".fits", "_stars_close.png")
        
    # if True:
    import rosalia as rs
    import pandas as pd
    import matplotlib.pyplot as plt
    from astropy.io import fits
    import numpy as np
    # Open the original flt file
    image_identity = rs.utils.exposure_inspector(input_name=flt_name, verbose=False, lite=True)
    
    # Get the detector corners: 
    detector_square_list = []
    for SCIEXT_i in rs.telescopes.Roman.WFI_SCAs:
        detector_corners = rs.detectors.get_detector_corners(data_shape=image_identity["DATA_SHAPE"][SCIEXT_i-1],
                                                             wcs=image_identity["ASTROPYWCS"][SCIEXT_i-1])
        detector_square_list.append(np.concatenate([detector_corners["corners_world"], detector_corners["corners_world"]]))

    
    
    RA_TARG = image_identity["RA_TARG"]
    DEC_TARG = image_identity["DEC_TARG"]
    hybrid_catalog = pd.read_csv(catalog_name)    
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    # Make a cut in the plot for the stars closer than radius
    target = SkyCoord(RA_TARG*u.deg, DEC_TARG*u.deg, unit="degree", frame="icrs")
    hybrid = SkyCoord(hybrid_catalog["ra"]*u.deg, hybrid_catalog["dec"]*u.deg, unit="degree", frame="icrs")
    separation = target.separation(hybrid)

    hybrid_catalog_close = hybrid_catalog[separation.value<2*radius]
    fig, ax = plt.subplots(figsize=figsize)
    
    if len(hybrid_catalog_close) > 1:
        plot_size, mag_lambda_range, plot_size_range = rs.plots.plot_stars_around(ax=ax, catalog=hybrid_catalog_close, 
                                           max_plot_size=50, min_plot_size=5, alpha=0.2)
    else:
        mag_lambda_range = hybrid_catalog_close["mag_lambda"].iloc[0]
        plot_size_range = 5


    plot_radec_limits = rs.gaia.find_ra_dec_constraints(ra=RA_TARG, dec=DEC_TARG, radius=radius/4)
    for detector_square in detector_square_list:
        plt.plot(detector_square[:,0], detector_square[:,1], alpha=0.5, color="red")
    ax.set_xlim((plot_radec_limits["ra_max"], plot_radec_limits["ra_min"]))
    ax.set_ylim((plot_radec_limits["dec_min"], plot_radec_limits["dec_max"]))

    # Make the legend 
    for i in range(len(mag_lambda_range) - 1, -1, -1):
        ax.scatter(-9999, -9999, transform=ax.transAxes, 
                   s=plot_size_range[i], marker="o", 
                   facecolor="grey", edgecolor="black",
                   alpha=0.2, label=str(mag_lambda_range[i])) 
    ax.legend(loc="best", title="AB Mag", fancybox=True, edgecolor="white", bbox_to_anchor=(1.2, 0.5))
    fig.tight_layout()
    plt.savefig(output_name, dpi=300)

    plt.show()
    return(output_name)


###############



def plot_ndi_main_offenders(input_name, scaled_main_off, catalog_name, ndi_level, figsize=(10,7)):
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits
    import astropy.wcs as astropy_wcs
    import matplotlib.colors as matplotlib_colors

    main_offender_db, multiline_text = rs.plots.main_offender_find_fraction_of_map(scaled_main_off, catalog_name)
    
    # Plot level NDI 
    ndi_name = os.environ["ROSALIACACHE"] + "CORE/NDI/RST/ndi_lvl" + str(ndi_level) + "_mean.fits"
    input_fits = fits.open(input_name)
    pa_point = input_fits[1].header["PA"]
    ra_point = input_fits[1].header["RA_TARG"]
    dec_point = input_fits[1].header["DEC_TARG"]

    if ndi_level==1:
        ndi_linthresh = 0.05
    if ndi_level==2:
        ndi_linthresh = 0.001
    else:
        ndi_linthresh = 0.01
   
    # Make the modified header based on the exposure properties.
    ndi_fits = fits.open(ndi_name)
    w = astropy_wcs.WCS(header=ndi_fits[0].header, fobj=ndi_fits, naxis=2)
    w.wcs.crota = -pa_point,-pa_point
    w.wcs.crval = ra_point,dec_point
    w.wcs.cdelt = -ndi_fits[0].header["CDELT1"],ndi_fits[0].header["CDELT2"]
    print(w)
    ndi_fits[0].header = w.to_header()
    
    fig, ax = plt.subplots(1, 1, figsize=figsize, subplot_kw=dict(projection=w))


    # ndi_header = fits.open(ndi_headerlet)[0].header
    
    # The transfer NDI maps units of 1/(superpixel size in mm2). 
    # We need to correct by that factor. 
    superpixel_mm2_area = (rs.telescopes.Roman.get_physical_pixelsize("WFI").to("mm").value * 512)**2
    im = ax.imshow(ndi_fits[0].data/superpixel_mm2_area, origin='lower', cmap="RdYlBu_r",
                           norm=matplotlib_colors.SymLogNorm(linthresh=ndi_linthresh, linscale=1,
                                            base=10))
    cbar = plt.colorbar(im, ax=ax, location='top', fraction=0.046, pad=0.04)
    cbar.set_label(label='Normalized Detector Irradiance (unitless)',size=15)

    ndi_w = astropy_wcs.WCS(header=ndi_fits[0].header, fobj=ndi_fits, naxis=2)
    # star_catalog = star_catalog_list[i]
    
    # CAR locations 
    alpha = 1
    import matplotlib.cm as cm

    for i in range(3):
        # CAR 1 A
        #x_stars, y_stars = ndi_w.wcs_world2pix(main_offender_db["ra"].iloc[i], main_offender_db["dec"].iloc[i], 0)
        ax.scatter(main_offender_db["ra"].iloc[i], main_offender_db["dec"].iloc[i], marker="h", edgecolor="black", color=cm.hot(i/3), alpha=alpha, s=150/(i+1), label=main_offender_db["source_name"].iloc[i])

    max_extent_ndi = np.abs(ndi_fits[0].header["CDELT2"]) * ndi_fits[0].data.shape[0]/2# This is the radial extent of the NDI map
    print(max_extent_ndi)
    angle_ticks = np.round(np.linspace(-max_extent_ndi, max_extent_ndi, 5),1)
    print(angle_ticks)
    pixel_ticks = angle_ticks/np.abs(ndi_fits[0].header["CDELT2"])  + ndi_fits[0].data.shape[0]/2
    print(pixel_ticks)
    #ax.set_xticks(pixel_ticks, angle_ticks, size=12)
    #ax.set_yticks(pixel_ticks, angle_ticks, size=12)
    ax.set_xlabel("RA (degrees)", size=15)
    ax.set_ylabel("Dec (degrees)", size=15)

    plt.legend(title="Main offenders", fancybox=True, edgecolor="white", bbox_to_anchor=(1.3, 0.5))
    ax.text(1.1, 0.95, multiline_text, 
            fontsize=9,
            color='black',
            transform=ax.transAxes,
            verticalalignment='top', 
            bbox=dict(facecolor='white', alpha=0.5))

    plt.tight_layout()
    output_plot_name = input_name.replace(".fits","_ndi_mainoff_location.png")
    plt.savefig(output_plot_name, dpi=300)
    return(output_plot_name)


#############

def plot_stars_around(ax, catalog, max_plot_size=50, min_plot_size=5, alpha=0.2):
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

    #print(bn.nanmin(ra_stars))
    #print(bn.nanmax(ra_stars))
    #print(bn.nanmin(dec_stars))
    #print(bn.nanmax(dec_stars))

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
    ax.scatter(ra_stars, dec_stars, marker="o", facecolor="grey", edgecolor="black",
                 alpha=alpha, s=plot_size)
    ax.set_xlabel("RA (ICRS)")
    ax.set_ylabel("DEC (ICRS)")


    ### Give the plot sizes and representative magnitudes. 
    mag_lambda = catalog["mag_lambda"]
    mag_lambda_range = np.unique(mag_lambda.astype("int"))
    from scipy import interpolate
    plot_size_interpolator = interpolate.interp1d(x=mag_lambda,
                                              y=plot_size, kind="linear",
                                              fill_value="extrapolate")
    plot_size_range = plot_size_interpolator(mag_lambda_range)


    return(plot_size,mag_lambda_range, plot_size_range)

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
    cbar = plt.colorbar(im, ax=axs[0], location='top', fraction=0.046, pad=0.04)
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
