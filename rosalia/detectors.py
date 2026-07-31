# Alejandro S. Borlaff. NASA Ames Research Center. a.s.borlaff@nasa.gov / asborlaff@gmail.com
# January 20, 2023.
#
# STRAYCOR/GNU module
# This module will hold all the programs related to detectors
#
#
# Version log:
# v.1.0 - 20 Enero 2023. First loading of programs inherited from former monolithic straycor.py
#
##########################################################


import numpy as np
import astropy.units as u

# ROSALIA modules
import rosalia as rs
from astropy import constants as const

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def mu2fe(mu, instrument, filter_name, telescope, verbose=False):
    ################################################
    # mu2electrons:
    # A function cannibalized from the Euclid Low Surface Brightness project. A.S.Borlaff et al. 2022.
    # to calculate the flux of photo-electrons per pixel that must reach to a telescope to provide a certain
    # surface brightness.
    # -----------------------------------------#
    #
    # Esto se debe calcular desde el Universo hacia la camara, no desde la CCD al Universo.
    # How many photons come from Fnu0
    #  Note: Fnu0 is Fnu0 at all wavelenghts
    # Ephot_0 = h*cs/lambda
    # dnu = c*dlambda/lambda**2
    # dnu*Ephot = dlambda/lambda / (h)
    # Nelectophot_0 = int(Fnu_0 / Ephot_0(nu) dnu) = int(Fnu0 lambda dlambda * qe(lambda) / c h lambda**2)
    ################################################

    if not isinstance(mu, (np.ndarray)):
        mu = np.array([mu])

    #fnu = 10**(-0.4*(mu + 56.10))*u.W.decompose()/u.meter**2/u.arcsec**2/u.Hz # W / m2 Hz arcsec2
    fnu = 10**(-0.4*(mu + 56.10))*u.W.decompose()/u.meter**2/u.Hz/u.arcsec**2 # W / m2 Hz arcsec2

    #print("Telescope: " + telescope + "/" + instrument)

    # For Hubble Space Telescope
    if telescope == "Hubble" or telescope == "HST":
        teles = rs.telescopes.Hubble

    # For Hubble Space Telescope
    if telescope == "Roman" or telescope == "RST":
        teles = rs.telescopes.Roman

    if telescope == "ARRAKIHS":
        teles = rs.telescopes.ARRAKIHS

    if telescope == "CSST":
        teles = rs.telescopes.CSST

    if telescope == "SPHEREx":
        teles = rs.telescopes.SPHEREx

    #print(teles)
    filter_db = teles.get_filter(instrument=instrument, filter_name=filter_name)
    mirror_radius = teles.mirror_radius
    pixscale = teles.get_pixscale(instrument=instrument).to("arcsec")
    #print(pixscale)
    A_telescope_factor = (np.pi * fnu*(mirror_radius**2)/const.h).decompose()

    # We initialize the array that will contain the photoelectrons
    Nephot_i = np.zeros(mu.shape)

    # Now we integrate the photo-electron flux across the filter bandwidth
    for i in range(len(filter_db["transmission_bins"])-1):
        lambda_i0 = filter_db["wavelength_bins"][i]   # Wavelength in m
        lambda_i1 = filter_db["wavelength_bins"][i+1] # Wavelength in m
        lambda_i = (lambda_i0 + lambda_i1)/2. # Average lambda
        T_i = (filter_db["transmission_bins"][i] + filter_db["transmission_bins"][i+1])/2. # Average quantum efficiency
        dlambda_i = lambda_i1 - lambda_i0 # Delta lambda
        #dnui = cs*dlambda_i/lambda_i0**2
        #nui = h*cs/lambda_i0
        # Nephot_i = Nephot_i + (fnu * np.pi * (evp.telescope_radius**2) * T_i * evp.exptime * dlambda_i  * (evp.pixsize**2) / (evp.h * lambda_i))
        # Nephot_i = Nephot_i + (fnu * np.pi * (evp.telescope_radius**2) * T_i * evp.exptime * dlambda_i  * (evp.pixsize**2) / (evp.h * lambda_i))
        Nephot_i = Nephot_i + (T_i * dlambda_i /lambda_i)
    #Nephot = bn.nansum(Nephot_i, axis=0)
    return((Nephot_i*A_telescope_factor*pixscale**2).decompose())


def fe2mu(fe, instrument, filter_name, telescope, verbose=False):
    # For Hubble Space Telescope
    # For Hubble Space Telescope
    if telescope == "Hubble" or telescope == "HST":
        teles = rs.telescopes.Hubble

    # For Hubble Space Telescope
    if telescope == "Roman" or telescope == "RST":
        teles = rs.telescopes.Roman

    if telescope == "ARRAKIHS":
        teles = rs.telescopes.ARRAKIHS

    if telescope == "CSST":
        teles = rs.telescopes.CSST

    if telescope == "SPHEREx":
        teles = rs.telescopes.SPHEREx

    pixscale = teles.get_pixscale(instrument=instrument)

    zeropoint = mu2fe(mu=0, instrument=instrument, filter_name=filter_name, telescope=telescope, verbose=False) # electrons / s / pixel

    if not isinstance(fe, u.quantity.Quantity):
        fe = fe/u.s

    mu = -2.5*np.log10(fe.value/zeropoint.value)
    # mu = -2.5*np.log10(electrons/evp.e_adu/evp.exptime/evp.pixsize**2) + zp
    return(mu)



def get_detector_corners(wcs):
    #input_fits = fits.open(input_name)
    #w = wcs.WCS(header=input_fits[ext].header, fobj=input_fits, naxis=2)
    data_shape = wcs.array_shape

    x_min = 0
    x_max = data_shape[1]
    y_min = 0
    y_max = data_shape[0]

    corners_pix = np.array([[x_min, y_min], [x_max, y_min], [x_max, y_max], [x_min, y_max]])
    #print(corners_pix)

    corners_world = np.array([wcs.all_pix2world(x_min, y_min, 0),
                              wcs.all_pix2world(x_max, y_min, 0),
                              wcs.all_pix2world(x_max, y_max, 0),
                              wcs.all_pix2world(x_min, y_max, 0)])

    # Fix RA if needed
    for i in range(4):
        corners_world[i,0] = corners_world[i,0] % 360
        # if corners_world[i,0] < -180: corners_world[i,0] = corners_world[i,0] + 360

    return({"corners_pix": corners_pix, "corners_world": corners_world})


###################################################
# def make_detector_grid(input_name, ext, step, mode="random"):

def make_array_grid(array, step, mode="random", nborder_points=64):
    x_min = 0
    x_max = int(array.shape[1]-1)
    y_min = 0
    y_max = int(array.shape[0]-1)

    x_mid = int((x_max - x_min)/2)
    y_mid = int((y_max - y_min)/2)

    if mode=="square":
        x_array = np.linspace(x_min, x_max, int(x_max/step)).astype("int")
        y_array = np.linspace(y_min, y_max, int(y_max/step)).astype("int")
        grid_xy = np.zeros((len(x_array)*len(y_array), 2))
        i = 0
        for xi in x_array:
            for yi in y_array:
                grid_xy[i, 0] = xi
                grid_xy[i, 1] = yi
                i = i + 1


    if mode=="random":
        total_nborder_points = nborder_points*4
        random_size = int((x_max-x_min)*(y_max-y_min)/step**2)
        x_array = np.random.randint(low=x_min, high=x_max, size=random_size)
        y_array = np.random.randint(low=y_min, high=y_max, size=random_size)

        grid_xy = np.zeros((random_size+total_nborder_points, 2))

        # First we fill the grid with the random points
        grid_xy[0:-total_nborder_points, 0] = x_array
        grid_xy[0:-total_nborder_points, 1] = y_array
        # Then we add the border points, that we want to ensure are always there
        # Side y=0 x_min to x_max
        grid_xy[-total_nborder_points:int(-3*nborder_points), 0] = np.linspace(x_min, x_max, nborder_points).astype("int") # We ensure that the corners are included
        grid_xy[-total_nborder_points:int(-3*nborder_points), 1] = np.zeros((nborder_points)) # We ensure that the corners are included

        # Side y=max x_min to x_max
        grid_xy[int(-3*nborder_points):int(-2*nborder_points), 0] = np.linspace(x_min, x_max, nborder_points).astype("int") # We ensure that the corners are included
        grid_xy[int(-3*nborder_points):int(-2*nborder_points), 1] = np.zeros((nborder_points))+y_max # We ensure that the corners are included

        # Side x=0 y_min to y_max
        grid_xy[int(-2*nborder_points):int(-nborder_points), 1] = np.linspace(y_min, y_max, nborder_points).astype("int") # We ensure that the corners are included
        grid_xy[int(-2*nborder_points):int(-nborder_points), 0] = np.zeros((nborder_points)) # We ensure that the corners are included

        # Side y=max x_min to x_max
        grid_xy[int(-nborder_points):, 1] = np.linspace(y_min, y_max, nborder_points).astype("int") # We ensure that the corners are included
        grid_xy[int(-nborder_points):, 0] = np.zeros((nborder_points))+x_max # We ensure that the corners are included

    #grid_world = w.all_pix2world(grid_xy[:,0], grid_xy[:,1], 0)
    return(grid_xy.astype("int"))



def make_detector_grid(w, step, mode="random"):
    #input_fits = fits.open(input_name)
    #w = wcs.WCS(header=input_fits[ext].header, fobj=input_fits, naxis=2)
    data_shape = w.array_shape
    x_min = 0
    x_max = data_shape[1] #input_fits[ext].data.shape[1]
    y_min = 0
    y_max = data_shape[0]

    x_mid = int((x_max - x_min)/2)
    y_mid = int((y_max - y_min)/2)

    if mode=="square":
        x_array = np.linspace(x_min, x_max, int(x_max/step)).astype("int")
        y_array = np.linspace(y_min, y_max, int(y_max/step)).astype("int")

        grid_xy = np.zeros((len(x_array)*len(y_array), 2))

        i = 0
        for xi in x_array:
            for yi in y_array:
                grid_xy[i, 0] = xi
                grid_xy[i, 1] = yi
                i = i + 1


    if mode=="random":
        random_size = int((x_max-x_min)*(y_max-y_min)/step**2)
        x_array = np.random.uniform(low=x_min, high=x_max, size=random_size)
        y_array = np.random.uniform(low=y_min, high=y_max, size=random_size)

        grid_xy = np.zeros((random_size+8, 2))
        grid_xy[0:-8, 0] = x_array
        grid_xy[-8:, 0] = np.array([x_min, x_mid, x_max, x_max, x_max, x_mid, x_min, x_min]) # We ensure that the corners are included
        grid_xy[0:-8, 1] = y_array
        grid_xy[-8:, 1] = np.array([y_min, y_min, y_min, y_mid, y_max, y_max, y_max, y_mid]) # We ensure that the corners are included



    grid_world = w.all_pix2world(grid_xy[:,0], grid_xy[:,1], 0)

    return({"grid_xy": grid_xy, "grid_world": grid_world})
#######################################################


def HST_ACS_counts_to_jy(flux_ACS, photflam, photplam):
    # See https://www.stsci.edu/hst/instrumentation/acs/data-analysis/zeropoints
    # ZPAB=−2.5∗log10(PHOTFLAM)−5∗log10(PHOTPLAM)−2.408
    # zp_ACS_AB = -2.5*np.log10(photflam)-5*np.log10(photplam)-2.408
    # def flambda_to_fnu(flambda, wavelength):

    magAB = -2.5*np.log10(flux_ACS) - 2.5*np.log10(photflam) - 5*np.log10(photplam) - 2.408
    flux_jy = 10**(0.4*(8.9 - magAB)) #rs.utils.flambda_to_fnu(flambda=flux_ACS*photflam, wavelength=photplam)
    return(flux_jy)


def HST_ACS_jy_to_counts(flux_jy, photflam, photplam):
    # See https://www.stsci.edu/hst/instrumentation/acs/data-analysis/zeropoints
    # ZPAB=−2.5∗log10(PHOTFLAM)−5∗log10(PHOTPLAM)−2.408
    # zp_ACS_AB = -2.5*np.log10(photflam)-5*np.log10(photplam)-2.408
    # def fnu_to_flambda(fnu, wavelength):

    #flux_jy = rs.utils.fnu_to_flambda(fnu=flux_ACS/photflam, wavelength=photplam)
    magAB = -2.5*np.log10(flux_jy) + 8.9
    ZP_hst_ab = - 2.5*np.log10(photflam) - 5*np.log10(photplam) - 2.408
    # muAB =  -2.5*np.log10(flux_jy) + 8.9 =  -2.5*np.log10(flux_ACS) + ZPHSTAB
    #
    flux_ACS = 10**(0.4*(ZP_hst_ab - magAB))

    #ma - 2.5*np.log10(photflam) - 5*log10(photplam) - 2.408

    return(flux_ACS)
