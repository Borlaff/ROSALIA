# Alejandro S. Borlaff. NASA Ames Research Center. a.s.borlaff@nasa.gov / asborlaff@gmail.com
# February 1, 2023.
#
# STRAYCOR/NDI module
# This module will hold all the programs to recreate the NDI
# including straylight estimation software
#
# Version log:
# v.1.0 - 1 February 2023. First loading of programs inherited from former
#                        Euclid project - EuclidStrayVIS.py
#
# """
# A source at a certain position (x1,y1) in the focal plane coordinate
# system creates a flux of straylight at a diffent location (x2,y2) of the
# focal plane of a telescope defined by the NDI function. This light reaches
# (x2,y2) due to non-primary optical paths, such as scattering, difraction,
# secondary reflections and other unwanted optical effects.

# From Euclid Mission Straylight Analysis: Impact on performance at Mission Level
# Ref: EUCL-SAP-TN-1-003

# "Assuming the entrance of the telescope is illuminated by a distant point
# source (collimated light), then the NDI is defined as the ratio of
# light irradiance (power per unit area) on the image plane to the source
# irradiance in object space at the entrance of the telescope. More basically
# the NDI describes the profile of the scattered light in the telescope
# focal plane for a point-like source at given position in the field of view."
# """
##########################################################

import numpy as np
from astropy import constants as const
import astropy.units as u
import rosalia as rs
from tqdm import tqdm
import bottleneck as bn
from scipy.interpolate import LinearNDInterpolator

def straylight_flux(mag, theta, phi, filter_name, instrument, telescope, exptime, mu_mode=False, verbose=False, ndi_mode="legacy"):
    # For Hubble Space Telescope
    if telescope == "Hubble" or telescope == "HST":
        teles = rs.telescopes.Hubble

    if verbose:
        print(teles)

    # Get the filter
    #print("Get the filter")
    filter_db = teles.get_filter(instrument=instrument, filter_name=filter_name)
    #print("Done")
    # Get the physical pixel scale
    #print("Get the pixscale")
    pixsize = teles.get_physical_pixelsize(instrument=instrument)
    #print("Done")

    # Get the NDI at a distance == theta

    #print("Get the NDI")
    NDI = teles.ndi_estimator(theta=theta, phi=phi, wavelength=filter_db["lambda_ref"], mode=ndi_mode)
    #print("Done")

    if verbose:
        print(filter_db)
        # Get the lambda max, lambda min, and lambda ref
        print(NDI)


    # Calculate the irradiance of the sources
    lambda_max = filter_db["lambda_max"].to(u.m)
    lambda_min = filter_db["lambda_min"].to(u.m)
    lambda_ref = filter_db["lambda_ref"].to(u.m)

    #print("Estimate the irradiance")
    Irradiance = const.c*(lambda_max-lambda_min)/(lambda_ref**2)*((10**(-0.4*(mag+56.1)))*u.W/u.meter**2/u.Hz)
    #print("Done")

    if verbose:
        print("Irradiance")
        print(Irradiance)
    # Now the flux per pixel using the NDI
    flux_per_pixel =(NDI*(pixsize**2)*filter_db["transmission_ref"]*lambda_ref*Irradiance/const.c/const.h).decompose()

    if verbose:
        print("flux_per_pixel")
        print(flux_per_pixel)

    if mu_mode:
        mu_straylight = rs.detectors.fe2mu(fe=flux_per_pixel.value,
                                          instrument=instrument,
                                          filter_name=filter_name,
                                          telescope=telescope)
        return(mu_straylight)

    else:
        return(flux_per_pixel) # e/s/pixel



############### LEGACY EUCLID CLASS FUNCTIONS ###################
class ndi_euclid():
    """
    This class contains the methods to calculate the NDI at different distances and orientations as a function of the Euclid NDI model (numeric envelope or non-axisymmetric ray-tracing model) as defined in EUCL-EST-TN-3-008_v2.0 Assessment e2e straylight performance.
    """
    #Numeric envelope NDI and non-axisymmetric ray-tracing are defined in EUCL-EST-TN-3-008_v2.0 Assessment e2e #straylight performance

    #To calculate the number of photons in a single pixel created by a source of magnitude m

    #Flux_photons = Flux_W / (h*c/lambda_ref)
    #where: h = Planck constant h = 6.62607004E-34 m2 kg / s
    #           c = Speed of light cs = 299792458 m/s
    #       lambda_ref: Reference (mean) wavelength for VIS = 725E-9 m

    #Flux_W = Irradiance_focal_plane * pixel_area
    #where: pixel_area = (12E-6)**2 # Area in m2 of a single pixel
    #       Irradiance_focal_plane = NDI * Irradiance_entrance * median Transmission
    #
    #Irradiance of the source at the entrance of the optical system (W m⁻²):

    #Irradiance_entrance(lambda) = 10⁻³ * c * 10^(-0.4*(mAB + 48.6)) * (lambda_max - lambda_min)/(lambda_max*lambda_min)
    #
    #And the NDI is calculated using one of the two methods described in this class.
    #
    #
    #Update: 01/02/2020 - So far, we have been assuming that Transmission = 1.
    #                     We have to assume a realistic transmission - Calculate the median of the
    #                     transmission curve between lambda_min and lambda_max

    # """

    # Here we will interpolate the values of A550,0 and theta1s.
    # Constant values to interpolate the envelope NDI model.
    # Check: Table 23 from Euclid Mission Straylight Analysis: Impact on performance at Mission Level

    # 26 Febrero 2020 - Reorder this points to make it F1 - F9 ordered
    # Keeping the old & wrong version for the record
    #NDI_A0550 = np.array([126, 177, 300,
    #                      84, 126, 240,
    #                      54,  90, 180])

    #NDI_theta1s = np.array([0.020, 0.017, 0.013,
    #                        0.025, 0.020, 0.015,
    #                        0.031, 0.024, 0.017])
    #NDI_x = np.array([-0.39, -0.39, -0.39,  0.00, 0.00, 0.00, 0.392, 0.392, 0.392])
    #NDI_y = np.array([-0.35,   0.0,  0.35, -0.35, 0.00, 0.35, -0.35, 0.000, 0.35])


    # VIS bounding rectangle has then sizes of:
    # size_x = 6 * ccd_size_x + 5 * ccd_gap_x = 302.502
    # size_y = 6 * ccd_size_y + 5 * ccd_gap_y = 336.632
    plate_scale=8.333 # (arcsec/mm) - VIS plate scale ("/mm) comes from the nominal 0.1"/pixel

    max_x = 168.466*plate_scale/60/60
    max_y = 151.5635*plate_scale/60/60 #

    # Arrays containing the x and y positions of the F1 - F9 reference positions of the VIS Focal Plane.
    # See Euclid PLM - PLM-IF-12: Optical Interface Control Document - ID: EUCL-ASFT-ICD-3-001

    # In order: F1 F2 F3 F4 F5 F6 F7 F8 F9
    FxVIS = np.array([0, 0, 0,
                  -max_x, -max_x, -max_x,
                   max_x,  max_x,  max_x])

    FyVIS = np.array([max_y, 0, -max_y,
                  max_y, 0, -max_y,
                  max_y, 0, -max_y])


    FxyVIS = np.array([[FxVIS[0],FyVIS[0]], [FxVIS[1],FyVIS[1]], [FxVIS[2],FyVIS[2]],
                   [FxVIS[3],FyVIS[3]], [FxVIS[4],FyVIS[4]], [FxVIS[5],FyVIS[5]],
                   [FxVIS[6],FyVIS[6]], [FxVIS[7],FyVIS[7]], [FxVIS[8],FyVIS[8]]])

    NDI_A0550   = np.array([240, 126,   84,   300,   177,   126,   180,   90,   54])
    NDI_theta1s = np.array([0.015, 0.020, 0.025, 0.013, 0.017, 0.020, 0.017, 0.024, 0.031])
    # theta0lambda550nm = 300
    ntheta0 = -1.8
    # theta1s = 0.013 # deg
    theta2s = 15 # deg
    theta2e = 35 # deg
    thetawd1 = 0.3 # deg
    thetawd2 = 2 # deg


    NDI_x     = FxVIS
    NDI_y     = FyVIS
    NDI_xy = np.array([[i,j] for i,j in zip(NDI_x, NDI_y)])

    NDI_A0550 = NDI_A0550
    NDI_theta1s = NDI_theta1s

    # We generate 2D interpolators (functions) that provide A0550 and thetha1s for any
    # position in the focal plane. This is only for the envelope NDI. Numeric non-axisymmetric
    # NDI model interpolation is not supported (yet).

    # 25 September 2024. interp2d stops being supported by scipy.
    # Transtioning to LinearNDInterpolator
    #A0550_interp   = interpolate.interp2d(NDI_x, NDI_y, NDI_A0550, kind='linear')
    #theta1s_interp = interpolate.interp2d(NDI_x, NDI_y, NDI_theta1s, kind='linear')


    #A0550_interp   = interpolate.interp2d(NDI_x, NDI_y, NDI_A0550, kind='linear')
    #theta1s_interp = interpolate.interp2d(NDI_x, NDI_y, NDI_theta1s, kind='linear')

    def n_ndi(theta):
        """
        This method calculates n as a function of theta = distance from the straylight source.
        Required by the method ndi_envelope to calculate the NDI profile
        n(θ_DEG) = n(0)*1/((1+[(θ_DEG/θ_wd1 )]^0.75 ) )*1/((1+[(θ_DEG/θ_wd2 )]^20 ) )
        """
        # A = n(0)
        A = ndi_euclid.ntheta0
        # B = 1/((1+[(θ_DEG/θ_wd1)]^0.75))
        B = 1/((1 + (theta/ndi_euclid.thetawd1)**0.75))
        # C = 1/((1+[(θ_DEG/θ_wd2)]^20))
        C = 1/((1 + (theta/ndi_euclid.thetawd2)**20))
        return A*B*C


    def a_ndi(x, y, theta, n_lambda):
        """
        This method calculates a as a function the position on the detector x,y
        the theta = distance from the straylight source and the n_lambda (wavelength in nm)
        Required by the method ndi_envelope to calculate the NDI profile
        A(θ_DEG, λ_nm ) = A(0, 550) *([λ_nm/550)]^(n(θ_DEG))
        """
        lenx = len(x)
        leny = len(y)
        if lenx != leny:
            return("Error! X and Y not the same length. Check your input")
        NDI_interp_xy = np.c_[NDI_x, NDI_y]
        A0550_interp = LinearNDInterpolator(NDI_interp_xy, NDI_A0550)

        A0550_interpolated_for_x_y = np.zeros(lenx)
        for i in range(lenx):
            A0550_interpolated_for_x_y[i] = A0550_interp(x[i],y[i])

        return(A0550_interpolated_for_x_y*(n_lambda/550)**(ndi_euclid.n_ndi(theta)))


    def ndi_envelope(x, y, theta, n_lambda):
        """
        This method calculates the NDI profile as a function the position on the detector x,y
        the theta = distance from the straylight source and the n_lambda (wavelength in nm)
        NDI (θ_DEG,λ_nm) = A(θ_DEG,λ_nm)*1/((1+[(θ_DEG/θ_1s )]^2))*((1+[(θ_DEG/θ_2e )]^4 ))/((1+[(θ_DEG/θ_2s )]^4))
        """
        NDI_interp_xy = np.c_[NDI_x, NDI_y]
        theta1s_interp = LinearNDInterpolator(NDI_interp_xy, NDI_theta1s)

        if isinstance(x, (int, float, complex)):
            x = np.array([x])

        if isinstance(y, (int, float, complex)):
            y = np.array([y])

        lenx = len(x)
        leny = len(y)
        if lenx != leny:
            return("Error! X and Y not the same length. Check your input")

        # A = A(theta,n_lambda)
        A = ndi_euclid.a_ndi(x, y, theta, n_lambda)
        #A0 = ndi.a_ndi(x, y, 0., n_lambda)
        # B = 1/((1+[(θ_DEG/θ_1s )]^2 ) )

        theta1s_interpolated_for_x_y = np.zeros(lenx)
        for i in range(lenx):
            theta1s_interpolated_for_x_y[i] = theta1s_interp(x[i],y[i])

        B = 1/((1 + (theta/theta1s_interpolated_for_x_y)**2))
        # C = ((1+[(θ_DEG/θ_2e )]^4 ))/((1+[(θ_DEG/θ_2s )]^4))
        C = ((1 + (theta/ndi_euclid.theta2e)**4))/((1 + (theta/ndi_euclid.theta2s)**4))
        return(A*B*C)

########################
# NDI Generation tools #
########################

def make_ndi_polar_radial_profile(data, header, center=None, n_phi_bins = 360, n_theta_bins=300, verbose=False, mode="average"):
    from scipy.interpolate import LinearNDInterpolator
    # We prepare the angle mask
    if verbose: print("Generating angle mask")
    angle_mask = rs.utils.create_angle_mask(xsize=data.shape[1], ysize=data.shape[0], q=1, theta=0, center=center) % 360

    # We prepare the radial mask
    if verbose: print("Generating radial mask")
    radial_mask = rs.utils.create_radial_mask(xsize=data.shape[1], ysize=data.shape[0], center=center)

    n_phi_bins = n_phi_bins+1
    n_theta_bins = n_theta_bins+1
    cdelt = np.abs((header["AAXISMAX"] - header["AAXISMIN"])/header["NAXIS1"])

    phi_bins = np.linspace(0, 360, n_phi_bins) # Angular direction

    log_thresh = 10
    max_r = np.nanmax(radial_mask)
    #theta_bins = cdelt*np.concatenate([np.linspace(0,log_thresh-1,log_thresh), np.logspace(np.log10(log_thresh), np.log10(max_r), n_theta_bins-log_thresh+1)])
    #theta_bins = cdelt*np.concatenate([np.array([0]), np.logspace(0, np.log10(np.sqrt(2)*data.shape[0]/2), n_theta_bins-1)]) # Radial direction
    theta_bins = np.linspace(0, max_r, n_theta_bins+1)

    x_ij = np.zeros(n_phi_bins*n_theta_bins)
    y_ij = np.zeros(n_phi_bins*n_theta_bins)
    theta_ij = np.zeros(n_phi_bins*n_theta_bins)
    phi_ij = np.zeros(n_phi_bins*n_theta_bins)

    if mode == "interpolate":
        theta_phi_profile = np.zeros((n_theta_bins, n_phi_bins))

        # Prepare the interpolator
        if verbose: print("Preparing the image interpolator...")
        image_interpolator = generate_image_interpolator(data)
        interpolated_ij = np.zeros(n_phi_bins*n_theta_bins)
        theta_grid, phi_grid = np.meshgrid(theta_bins, phi_bins)
        ra_grid, dec_grid = rs.utils.thetaphi_2_radec(theta_grid, phi_grid)
        ra_NDI = np.linspace(header["AAXISMIN"], header["AAXISMAX"], header["NAXIS1"])
        dec_NDI = np.linspace(header["BAXISMIN"], header["BAXISMAX"], header["NAXIS2"])
        ra_NDI_map_grid, dec_NDI_map_grid = np.meshgrid(ra_NDI, dec_NDI)
        radec_points = np.c_[ra_NDI_map_grid.flatten(), dec_NDI_map_grid.flatten()]
        ndi_points = data.flatten()
        if verbose: print("Interpolating data at the polar radial grid locations...")
        NDI_image_interpolator = LinearNDInterpolator(radec_points, ndi_points, fill_value=0)
        theta_phi_profile = NDI_image_interpolator(ra_grid, dec_grid)

    if mode == "average":

        if verbose:
            plt.imshow(radial_mask)
            plt.colorbar()
            plt.show()
            plt.imshow(angle_mask)
            plt.colorbar()
            plt.show()
            print("Theta")
            print(theta_bins)
            print("phi")
            print(phi_bins)
        theta_bins_centroid = []
        phi_bins_centroid = []
        rs.utils.save_fits(radial_mask, "default_radial_mask.fits")
        rs.utils.save_fits(angle_mask, "default_angle_mask.fits")

        for i in range(len(theta_bins)-1):
            theta_bins_centroid.append((theta_bins[i] + theta_bins[i+1])/2.)

        theta_bins_centroid = np.array(theta_bins_centroid)
        for j in range(len(phi_bins)-1):
            phi_bins_centroid.append((phi_bins[j] + phi_bins[j+1])/2.)
        phi_bins_centroid = np.array(phi_bins_centroid)

        theta_phi_profile = np.zeros((len(theta_bins_centroid), len(phi_bins_centroid)))

        for i in tqdm(range(len(theta_bins_centroid)-1), position=0, leave=True):
            pixel_to_bin_mask_radial = (radial_mask > theta_bins[i]) & (radial_mask <= theta_bins[i+1])

            for j in range(len(phi_bins_centroid-1)):

                pixels_to_bin_mask_angle = (angle_mask > phi_bins[j]) & (angle_mask <= phi_bins[j+1])
                pixels_to_bin_mask = pixel_to_bin_mask_radial & pixels_to_bin_mask_angle

                theta_phi_profile[i,j] = bn.nanmean(data[pixels_to_bin_mask])

    return({"theta_phi_profile": theta_phi_profile, "theta_bins": np.array(theta_bins_centroid)*cdelt, "phi_bins": np.array(phi_bins_centroid)})
