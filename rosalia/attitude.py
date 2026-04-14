"""
File: Roman_Attitude.py
Author: Maxime J Rizzo
Email: maxime.j.rizzo@nasa.gov
Date: 2025-11-18
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_sun
from astropy import units as u
import datetime
import pandas as pd
import pysiaf

def get_radec(attitude_matrix):
    """
    Compute RA and DEC from spacecraft attitude matrix in icrs.
    
    Parameters:
    -----------
    attitude_matrix: 3x3 numpy array, spacecraft attitude matrix in icrs

    Returns:
    --------
    ra: float, Right Ascension in degrees
    dec: float, Declination in degrees
    """
    # The first column of the attitude matrix represents the pointing direction
    x, y, z = attitude_matrix[:, 0]

    # Convert Cartesian coordinates to spherical
    r = np.sqrt(x**2 + y**2 + z**2)
    ra = np.arctan2(y, x)
    dec = np.arcsin(z / r)
    ra = np.degrees(ra) % 360
    dec = np.degrees(dec)
    return ra, dec

def get_vector(target):
    """
    Convert SkyCoord object to unit vector.
    Parameters:
    -----------
    target: SkyCoord object
        Returns:
    vector: numpy array, unit vector representation of target
        """
    vector = target.cartesian.xyz.value
    vector = np.array([vector[0], vector[1], vector[2]])
    vector = vector / np.linalg.norm(vector)
    return vector

def load_targets(target_list):
    """
    Load target coordinates from Excel file.
    Parameters:
    -----------
    target_list: str, path to Excel file containing target coordinates
        Returns:
    tgt_list: list of SkyCoord objects
    """
    if 'xlsx' not in target_list:
        raise ValueError("The provided file is not an Excel file.")

    dat = pd.read_excel(target_list)

    if 'RA' not in dat.columns or 'DEC' not in dat.columns:
        raise ValueError("The Excel file must contain 'RA' and 'DEC' columns.")
    
    tgt_list = [SkyCoord(ra=dat['RA'].iloc[i]*u.deg,dec=dat['DEC'].iloc[i]*u.deg, frame='icrs', unit='deg') for i in range(len(dat))]

    return tgt_list

def get_large_aperture(attitude_matrix, named_aperture='WFI_CEN', scalex=1.2, scaley=1.4):
    rsiaf = pysiaf.Siaf('Roman')
    ap = rsiaf[named_aperture]
    ap.set_attitude_matrix(attitude_matrix)

    scale_Y = scaley
    scale_X = scalex
    X = np.array([ap.XIdlVert1, ap.XIdlVert2, ap.XIdlVert3, ap.XIdlVert4])*scale_X
    Y = np.array([ap.YIdlVert1, ap.YIdlVert2, ap.YIdlVert3, ap.YIdlVert4])*scale_Y
    
    skyRa, skyDec = ap.idl_to_sky(X, Y)
    
    ap_large = 'POLYGON ICRS {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} {:.8f} '.format(
                    skyRa[0], skyDec[0], skyRa[1], skyDec[1], skyRa[2], skyDec[2], skyRa[3], skyDec[3])

    return ap_large
    

def plot_targets_by_magnitude(x_coords, y_coords, magnitudes, ax=None, title="Target Distribution by Magnitude", 
                            show_colorbar=True, show_legend=True, show_stats=True, minmag = 5, maxmag = 17):
    """
    Create a scatter plot with marker sizes based on astronomical magnitudes.
    
    Parameters:
    x_coords: list or array of X coordinates
    y_coords: list or array of Y coordinates  
    magnitudes: list or array of astronomical magnitudes
    ax: matplotlib axis object (if None, creates new figure)
    title: plot title
    show_colorbar: whether to show colorbar (default True)
    show_legend: whether to show magnitude legend (default True)
    show_stats: whether to print statistics (default True)
    
    Returns:
    ax: the axis object used for plotting
    scatter: the scatter plot object (useful for external colorbar creation)
    """
    
    # Convert to numpy arrays for easier manipulation
    x = np.ma.array(x_coords)
    y = np.ma.array(y_coords)
    mag = np.ma.array(magnitudes)
    
    # Create axis if not provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 8))
        created_figure = True
    else:
        created_figure = False
    
    # Calculate marker sizes (brighter stars = larger markers)
    # Since lower magnitude = brighter, we invert the relationship
    # max_mag = np.ma.max(mag)
    # min_mag = np.ma.min(mag)
    # print(min_mag, max_mag)
    if minmag is None:
        min_mag = np.ma.min(mag)
    else:
        min_mag = minmag

    if maxmag is None:
        max_mag = np.ma.max(mag)
    else:
        max_mag = maxmag
    
    
    # Scale marker sizes: brightest stars get size ~200, faintest get size ~20
    marker_sizes = 200 - 180 * (mag - min_mag) / (max_mag - min_mag) + 1
    
    # Create the scatter plot
    scatter = ax.scatter(x, y, s=marker_sizes, c=mag, cmap='viridis_r', 
                        alpha=0.7, edgecolors='black', linewidth=0.5)
    
    # Add colorbar only if requested and we created the figure
    if show_colorbar and created_figure:
        cbar = plt.colorbar(scatter, ax=ax, label='Magnitude')
    
    # Customize the plot
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    
    # Add magnitude scale reference in legend
    if show_legend:
        mag_ranges = np.linspace(min_mag, max_mag, 5)
        size_ranges = 200 - 180 * (mag_ranges - min_mag) / (max_mag - min_mag) + 20
        
        legend_elements = []
        for i, (m, s) in enumerate(zip(mag_ranges, size_ranges)):
            legend_elements.append(plt.scatter([], [], s=s, c='gray', alpha=0.7, 
                                             edgecolors='black', linewidth=0.5,
                                             label=f'Mag {m:.1f}'))
        
        ax.legend(handles=legend_elements, title="Magnitude Scale", 
                 loc='upper right', bbox_to_anchor=(1.15, 1))
    
    # Adjust layout only if we created the figure
    if created_figure:
        plt.tight_layout()
        plt.show()
    
    # Print statistics if requested
    if show_stats:
        print(f"Total targets: {len(x)}")
        print(f"Magnitude range: {np.amin(mag):.2f} to {np.amax(mag):.2f}")
        print(f"Brightest target (lowest mag): {np.amin(mag):.2f} at ({x[np.argmin(mag)]:.2f}, {y[np.argmin(mag)]:.2f})")
        print(f"Faintest target (highest mag): {np.amax(mag):.2f} at ({x[np.argmax(mag)]:.2f}, {y[np.argmax(mag)]:.2f})")
    
    return ax, scatter

def create_subplot_grid(n_plots, figsize_per_plot=(4, 3), max_cols=None, 
                       spacing={'hspace': 0.3, 'wspace': 0.3}, sharey=True, sharex=True):
    """
    Create an optimal grid layout for an arbitrary number of plots
    
    Parameters:
    -----------
    n_plots: Number of plots to create
    figsize_per_plot: Size of each individual plot (width, height)
    max_cols: Maximum number of columns (None for automatic)
    spacing: Dictionary with hspace and wspace for subplot spacing
    
    Returns:
    --------
    fig: Figure object
    axes: Array of axes objects (flattened)
    """
    
    if n_plots <= 0:
        raise ValueError("Number of plots must be positive")
    
    # Calculate optimal grid dimensions
    if max_cols is None:
        # Automatic grid sizing - aim for roughly square layout
        n_cols = int(np.ceil(np.sqrt(n_plots)))
        n_rows = int(np.ceil(n_plots / n_cols))
    else:
        n_cols = min(max_cols, n_plots)
        n_rows = int(np.ceil(n_plots / n_cols))
    
    # Calculate figure size
    fig_width = n_cols * figsize_per_plot[0]
    fig_height = n_rows * figsize_per_plot[1]
    
    # Create figure and subplots
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(fig_width, fig_height), sharey=sharey, sharex=sharex)
    
    # Handle case where we only have one subplot
    if n_plots == 1:
        axes = [axes]  # Make it iterable
    elif n_rows == 1 or n_cols == 1:
        axes = axes.flatten()  # Already 1D, but ensure it's flattened
    else:
        axes = axes.flatten()  # Convert 2D array to 1D
    
    # Hide unused subplots
    for i in range(n_plots, len(axes)):
        axes[i].set_visible(False)
    
    # Adjust spacing
    plt.subplots_adjust(**spacing)
    
    return fig, axes[:n_plots]  # Return only the axes we need


class RomanPointing:
    """
    Spacecraft orientation and pointing for celestial targets.

    This module provides functionality for calculating and visualizing
    spacecraft attitude, target coordinates, and sun angles for
    astronomical observations. It supports operations in various
    coordinate systems including ICRS, ecliptic, and galactic.

    Key features:
    - Calculate spacecraft attitude for given celestial targets
    - Convert between different coordinate systems
    - Visualize pointing and attitude in 3D and 2D projections
    - Handle pitch and roll maneuvers
    - Check visibility constraints
    """

    def __init__(self, observation_date=None):
        """
        Initialize spacecraft pointing system
        
        Parameters:
        -----------
        observation_date: datetime object or astropy Time object
        """
        if observation_date is None:
            self.observation_date = Time.now()
        elif isinstance(observation_date, datetime.datetime):
            self.observation_date = Time(observation_date)
        elif isinstance(observation_date, str):
            self.observation_date = datetime.fromisoformat(observation_date)
        else:
            self.observation_date = observation_date
            
        self.spacecraft_attitude = np.eye(3)  # Identity rotation matrix
        self.target_coord = None
        self.sun_coord = None
        self.pitch_limits = [-36,36]*u.deg
        self._update_sun_position()

    def _update_sun_position(self):
        """
        Update Sun position for the observation date.

        This method calculates the current position of the Sun based on the
        observation date and updates the sun_coord attribute of the class.

        Returns:
        --------
        None
        """
        suncoord = get_sun(self.observation_date)
        self.sun_coord = SkyCoord(ra=suncoord.ra,dec=suncoord.dec, frame='icrs')

    def get_attitude_quaternion(self, attitude=None):
        """
        Get attitude as quaternion [w, x, y, z]
        Parameters:
        -----------
        attitude: 3x3 numpy array, optional spacecraft attitude matrix

        Returns:
        --------
        numpy array, quaternion representation of attitude
        """
        if attitude is None:
            attitude = self.spacecraft_attitude
    
        rotation = R.from_matrix(self.spacecraft_attitude)
        quat = rotation.as_quat()  # [x, y, z, w]
        return np.array([quat[3], quat[0], quat[1], quat[2]])  # [w, x, y, z]
    
    def get_attitude_euler(self, attitude=None, sequence='ZYX', degrees=True):
        """
        Get attitude as Euler angles
        Parameters:
        -----------
        attitude: 3x3 numpy array, optional spacecraft attitude matrix
        sequence: str, rotation sequence
        degrees: bool, if True return angles in degrees, else in radians

        Returns:
        --------
        numpy array, Euler angles
        """
        if attitude is None:
            attitude = self.spacecraft_attitude
    
        rotation = R.from_matrix(attitude)
        return rotation.as_euler(sequence.upper(), degrees=degrees)
    
    def get_sun_angle(self, target=None):
        """
        Get angle between spacecraft pointing direction and Sun.

        Parameters:
        -----------
        target : astropy.coordinates.SkyCoord, optional
            The target coordinates. If not provided, uses the current target_coord.

        Returns:
        --------
        float or None
            The angular separation between the target and the Sun in degrees.
            Returns None if no target is set and no target is provided.

        Notes:
        ------
        This method calculates the angular separation between the spacecraft's
        pointing direction (or provided target) and the Sun's position.
        """
        if target is None:
            if self.target_coord is None:
                return None
                
            # Get Sun and target separation
            separation = self.sun_coord.separation(self.target_coord)[0]
        else:
            separation = self.sun_coord.separation(target)[0]
        
        return separation.deg       
    
    def get_ecliptic_coordinates(self, coord):
        """
        Convert coordinates to ecliptic longitude and latitude

        Parameters:
        -----------
        coord: SkyCoord object

        Returns:
        --------
        tuple, (ecliptic longitude, ecliptic latitude) in degrees
        """
        ecliptic_coord = coord.transform_to('barycentricmeanecliptic')
        return ecliptic_coord.lon.deg, ecliptic_coord.lat.deg

    def set_target_using_radec(self, ra, dec, roll=0.0, frame='icrs', unit='deg'):
        """
        Set target star coordinates
        
        Parameters:
        -----------
        ra: float, right ascension
        dec: float, declination  
        roll: float, Sun roll angle from nominal in degrees
        frame: str, coordinate frame ('icrs', 'fk5', etc.)
        unit: str, angular unit ('deg', 'hour' for RA, etc.)
        """
        if unit == 'deg':
            ra_unit = u.deg
            dec_unit = u.deg
        elif unit == 'hour':
            ra_unit = u.hour
            dec_unit = u.deg
        else:
            ra_unit = u.Unit(unit)
            dec_unit = u.Unit(unit)
            
        self.target_coord = SkyCoord(ra=ra*ra_unit, dec=dec*dec_unit, 
                                   frame=frame, obstime=self.observation_date)
        
        # Calculate spacecraft attitude to point at target
        self._calculate_pointing_attitude(roll=roll)

    def set_target(self, target, roll=0.):
        """
        Set target using SkyCoord object

        Parameters:
        -----------
        target: SkyCoord object
        roll: float, Sun roll angle from nominal in degrees
        """
        self.target_coord = target
        self._calculate_pointing_attitude(roll=roll)

    def get_target_new_pitch(self, delta_pitch=0.0):
        """
        Compute the RA and DEC of a new target offset in pitch angle.

        Parameters:
        -----------
        delta_pitch : float, optional
            The pitch angle offset in degrees. Positive values pitch away from Sun, negative values
            pitch towards Sun. Default is 0.0.

        Returns:
        --------
        target : astropy.coordinates.SkyCoord
            A SkyCoord object representing the new target's position in ICRS coordinates.
        new_attitude : numpy.ndarray
            A 3x3 array representing the new spacecraft attitude matrix.

        Notes:
        ------
        The function uses the current spacecraft attitude and applies a rotation
        around the Y-axis (pitch) to determine the new pointing direction.
        """

        pitch_rad = np.radians(delta_pitch)
        pitch_matrix = np.array([
            [np.cos(pitch_rad), 0, np.sin(pitch_rad)],
            [0, 1, 0],
            [-np.sin(pitch_rad), 0, np.cos(pitch_rad)]
        ])
        new_attitude = np.dot(self.spacecraft_attitude, pitch_matrix)

        ra, dec = get_radec(new_attitude)
        target = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, 
                            frame='icrs', obstime=self.observation_date)

        return target, new_attitude
    
    def delta_pitch_roll(self, dpitch=0.0, droll=0.0):
        """
        Calculate new pointing coordinates after applying pitch and roll offsets.

        Parameters:
        -----------
        dpitch : float, optional
            Change in pitch angle in degrees. Default is 0.0.
        droll : float, optional
            Change in roll angle in degrees. Default is 0.0.

        Returns:
        --------
        ra : float
            New Right Ascension in degrees.
        dec : float
            New Declination in degrees.
        pa : float
            New Position Angle in degrees.
        new_pitch : float
            New pitch angle in degrees.

        Notes:
        ------
        This method applies a pitch rotation followed by a roll rotation to the
        current spacecraft attitude. It then calculates the new pointing coordinates,
        position angle, and pitch angle based on the resulting attitude.
        """

        pitch_rad = np.radians(dpitch)
        pitch_matrix = np.array([
            [np.cos(pitch_rad), 0, np.sin(pitch_rad)],
            [0, 1, 0],
            [-np.sin(pitch_rad), 0, np.cos(pitch_rad)]
        ])
        new_attitude = np.dot(self.spacecraft_attitude, pitch_matrix)

        ra, dec = get_radec(new_attitude)
        target = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, 
                            frame='icrs', obstime=self.observation_date)
        
        new_pitch = self.get_pitch_angle(target)

        # now create rotation matrix for the roll about X
        roll_rad = np.radians(droll)
        roll_matrix = np.array([
            [1, 0, 0],
            [0, np.cos(roll_rad), -np.sin(roll_rad)],
            [0, np.sin(roll_rad), np.cos(roll_rad)]
        ])

        new_attitude = np.dot(new_attitude, roll_matrix)

        pa = self.get_position_angle(new_attitude)
        # pa = pysiaf.utils.rotations.posangle(new_attitude, 0, 0)

        return ra, dec, pa.value%360, new_pitch.value
        # return ra, dec, pa, new_pitch.value


    def get_pitch_angle(self, target=None):
        """
        Get the pitch angle between the spacecraft pointing and the Sun.

        Parameters:
        -----------
        target : astropy.coordinates.SkyCoord, optional
            The target coordinates. If not provided, uses the current target_coord.

        Returns:
        --------
        astropy.units.Quantity
            The pitch angle in degrees. Positive values indicate the spacecraft is
            pitched up relative to the Sun, negative values indicate it's pitched down.

        Notes:
        ------
        The pitch angle is defined as the angle between the spacecraft-Sun vector
        and the spacecraft pointing direction, minus 90 degrees. This means a pitch
        of 0 degrees corresponds to the Sun being perpendicular to the pointing direction.
        """
        sun_angle =  self.get_sun_angle(target=target)
        pitch = (sun_angle-90.)*u.deg
        if target is None:
            self.pitch = pitch
        return pitch

    def _calculate_pointing_attitude(self, roll=0.0):
        """Calculate spacecraft attitude matrix for target pointing"""
        if self.target_coord is None:
            return
            
        # Convert to Cartesian coordinates (unit vector)
        # target_cart = self.target_coord.cartesian.xyz.value
        # target_vector = np.array([target_cart[0], target_cart[1], target_cart[2]])
        # target_vector = target_vector / np.linalg.norm(target_vector)
        target_vector = get_vector(self.target_coord)
        
        # Define spacecraft axes convention
        # X-axis points to target, Y-axis in ecliptic plane when possible
        x_axis = target_vector

        # Y axis perpendicular to sun vector
        # Choose Y to be perpendicular to both X and sun vector
        sun_vector = get_vector(self.sun_coord)[:, 0]
        temp_y = np.cross(sun_vector, x_axis)
        if np.linalg.norm(temp_y) < 1e-6:
            # If sun and target are aligned, choose arbitrary perpendicular
            temp_y = np.array([0, 1, 0]) if abs(x_axis[1]) < 0.9 else np.array([1, 0, 0])
        y_axis = temp_y / np.linalg.norm(temp_y)

        z_axis = np.cross(x_axis, y_axis)

        # Ensure Z axis points toward sun (positive dot product with sun vector)
        if np.dot(z_axis, sun_vector) < 0:
            z_axis = -z_axis
            y_axis = -y_axis  # Maintain right-handed system
        
        # Create rotation matrix [X, Y, Z] as columns
        spacecraft_attitude_no_roll = np.column_stack([x_axis, y_axis, z_axis])

        # now create rotation matrix for the roll about X
        roll_rad = np.radians(roll)
        roll_matrix = np.array([
            [1, 0, 0],
            [0, np.cos(roll_rad), -np.sin(roll_rad)],
            [0, np.sin(roll_rad), np.cos(roll_rad)]
        ])
        self.spacecraft_attitude = np.dot(spacecraft_attitude_no_roll, roll_matrix)
        self.get_pitch_angle()

    def _calculate_pointing_attitude_euler(self, v3pa=0.):
        """
        Calculate spacecraft attitude matrix using Euler angles.

        This method computes the spacecraft attitude matrix by applying a series
        of rotations: first around Z-axis (RA), then Y-axis (DEC), and finally
        X-axis (position angle).

        Parameters:
        -----------
        v3pa : float, optional
            Position angle in degrees (positive is East of North). Default is 0.

        Returns:
        --------
        numpy.ndarray
            3x3 spacecraft attitude matrix.

        Notes:
        ------
        The rotations are applied in the order: Z (RA) -> Y (-DEC) -> X (-pos_angle).
        The resulting matrix represents the transformation from celestial coordinates
        to spacecraft body coordinates.
        """
        Zrot = R.from_euler('Z', np.radians(self.target_coord.ra)).as_matrix()
        Yrot = R.from_euler('Y', np.radians(-self.target_coord.dec)).as_matrix()
        Xrot = R.from_euler('X', np.radians(-v3pa)).as_matrix()
        
        return np.dot(np.dot(Zrot, Yrot), Xrot)

    def get_position_angle(self, attitude=None):
        """
        Compute the position angle (V3PA) of the spacecraft.
        
        Position Angle is the angle from the plane containing the +X axis and 
        the vector to Celestial North to the plane containing the +X and +Z axes,
        measured in the direction of celestial east as viewed from the spacecraft origin.
        
        Parameters:
        -----------
        x_axis : array_like
            Spacecraft +X axis direction vector (3D)
        z_axis : array_like  
            Spacecraft +Z axis direction vector (3D)
        north_vector : array_like, optional
            Celestial North direction vector (default: [0,0,1] for +Z celestial)
        
        Returns:
        --------
        position_angle : float
            Position angle in radians
        """

        if attitude is None:
            attitude = self.spacecraft_attitude

        x_axis = attitude[:,0] 
        x_axis = np.array(x_axis) / np.linalg.norm(x_axis)
        z_axis = attitude[:,2] 
        z_axis = np.array(z_axis) / np.linalg.norm(z_axis)
        north_vector = np.array([0, 0, 1])
        north_vector = np.array(north_vector) / np.linalg.norm(north_vector)
        
        # Compute normal to plane containing +X axis and celestial north
        n1 = np.cross(x_axis, north_vector)
        n1 = n1 / np.linalg.norm(n1)
        
        # Compute normal to plane containing +X and +Z axes
        n2 = np.cross(x_axis, z_axis)
        n2 = n2 / np.linalg.norm(n2)
        
        # Compute angle between the two planes
        cos_angle = np.dot(n1, n2)
        cos_angle = np.clip(cos_angle, -1, 1)  # Handle numerical precision
        
        # Determine sign using triple scalar product
        # Orientation is East of North, matching pysiaf
        sign = -np.sign(np.dot(x_axis, np.cross(n1, n2)))
        
        position_angle = sign * np.arccos(cos_angle)
        
        return (np.degrees(position_angle))*u.deg

    def visualize_pointing(self, show_ecliptic=True, show_sun=True, 
                          show_constellation_grid=False, fig_size=(15, 10), target_list=None):
        """
        Visualize spacecraft pointing in celestial sphere context
        
        Parameters:
        show_ecliptic: bool, show ecliptic plane
        show_sun: bool, show Sun position
        show_constellation_grid: bool, show RA/DEC grid
        fig_size: tuple, figure size
        """
        fig = plt.figure(figsize=fig_size)
        
        # Create subplots
        ax1 = plt.subplot(2, 2, 1, projection='3d')  # 3D celestial sphere
        ax2 = plt.subplot(2, 2, 2, projection='mollweide')  # Sky map
        ax3 = plt.subplot(2, 2, 3, projection='3d')  # Spacecraft attitude
        # ax4 = plt.subplot(2, 2, 4)  # Information panel
        ax4 = plt.subplot(2, 2, 4, projection='mollweide')  # Sky map

        self._plot_3d_celestial_sphere(ax1, show_ecliptic, show_sun, target_list = target_list)
        self._plot_sky_map(ax2, show_ecliptic, show_sun, show_constellation_grid, target_list = target_list)
        self._plot_spacecraft_attitude(ax3, target_list = target_list)
        # self._plot_info_panel(ax4)
        self._plot_sky_map_galactic(ax4, show_ecliptic, show_sun, show_constellation_grid, target_list = target_list)

        plt.tight_layout()
        plt.show()
    
    def _plot_3d_celestial_sphere(self, ax, show_ecliptic, show_sun, target_list=None):
        """
        Plot a 3D celestial sphere with target and spacecraft orientation.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The 3D axes to plot on.
        show_ecliptic : bool
            If True, plot the ecliptic plane.
        show_sun : bool
            If True, plot the Sun's position.
        target_list : list of SkyCoord, optional
            List of additional targets to plot.

        Returns:
        --------
        None
        """
        # Create sphere
        uu = np.linspace(0, 2 * np.pi, 50)
        v = np.linspace(0, np.pi, 50)
        x_sphere = np.outer(np.cos(uu), np.sin(v))
        y_sphere = np.outer(np.sin(uu), np.sin(v))
        z_sphere = np.outer(np.ones(np.size(uu)), np.cos(v))
        
        ax.plot_surface(x_sphere, y_sphere, z_sphere, alpha=0.1, color='lightblue')
        
        # Plot coordinate axes
        axis_length = 1.2
        ax.quiver(0, 0, 0, axis_length, 0, 0, color='red', alpha=0.5, 
                 arrow_length_ratio=0.05, label='X (Vernal Equinox)')
        
        ax.quiver(0, 0, 0, 0, axis_length, 0, color='green', alpha=0.5, 
                 arrow_length_ratio=0.05, label='Y (RA=90)')
        ax.quiver(0, 0, 0, 0, 0, axis_length, color='blue', alpha=0.5, 
                 arrow_length_ratio=0.05, label='Z (Celestial North)')
        
        # Plot ecliptic plane
        if show_ecliptic:
            
            # create a ring of targets to draw the ecliptic
            ecl_lon = np.linspace(0, 360, 100)
            ecl_lat = np.zeros_like(ecl_lon)
            
            # Convert ecliptic to Cartesian coordinates
            ecl_coords = SkyCoord(lon=ecl_lon*u.deg, lat=ecl_lat*u.deg, 
                                frame='barycentricmeanecliptic', 
                                obstime=self.observation_date)
            cart_coords = ecl_coords.transform_to('icrs').cartesian.xyz.value
            cart_coords /= np.linalg.norm(cart_coords, axis=0)

            ax.plot(cart_coords[0,:],cart_coords[1,:],cart_coords[2,:], 'orange', linewidth=2, label='Ecliptic')
        
        # Plot Sun position
        if show_sun and self.sun_coord:
            sun_pos = get_vector(self.sun_coord)
            ax.scatter(*sun_pos, color='yellow', s=200, marker='*', 
                      edgecolor='orange', linewidth=2, label='Sun')
        
        # Plot target
        if self.target_coord:
            target_pos = get_vector(self.target_coord)
            ax.scatter(*target_pos, color='red', s=100, marker='*', 
                      edgecolor='darkred', linewidth=2, label='Target Star')
            
            # Plot spacecraft pointing vector
            ax.quiver(0, 0, 0, *target_pos, color='red', linewidth=3, 
                     arrow_length_ratio=0.05, label='Spacecraft Pointing')
            
        if target_list is not None:
            for target in target_list:
                vec = get_vector(target)
                ax.scatter(*vec, color='green', s=50, marker='o', 
                        edgecolor='darkgreen', linewidth=2)

        
        ax.set_xlim([-1.5, 1.5])
        ax.set_ylim([-1.5, 1.5])
        ax.set_zlim([-1.5, 1.5])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.set_box_aspect([1,1,1])
        ax.set_title('3D Celestial Sphere View')
   
    def _plot_sky_map_galactic(self, ax, show_ecliptic, 
                      show_sun,
                      show_constellation_grid,
                      target_list = None,
                      show_cvz=True):
        """
        Plot Mollweide projection sky map in Galactic coordinates.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes to plot on.
        show_ecliptic : bool
            If True, plot the ecliptic plane.
        show_sun : bool
            If True, plot the Sun's position.
        show_constellation_grid : bool
            If True, show the Galactic coordinate grid.
        target_list : list of SkyCoord, optional
            List of additional targets to plot.
        show_cvz : bool, optional
            If True, show the continuous viewing zone. Default is True.

        Returns:
        --------
        None
        """

        if show_constellation_grid:
            # Galactic longitude/latitude grid
            l_lines = np.arange(0, 360, 30)
            b_lines = np.arange(-90, 91, 30)
            
            for l in l_lines:
                b_range = np.linspace(-90, 90, 100)
                l_range = np.full_like(b_range, l)
                
                # Convert to radians and adjust for mollweide
                l_rad = np.deg2rad(l_range - 180)  # Center at 0
                b_rad = np.deg2rad(b_range)
                
                ax.plot(l_rad, b_rad, 'gray', alpha=0.3, linewidth=0.5)
            
            for b in b_lines:
                l_range = np.linspace(0, 360, 100)
                b_range = np.full_like(l_range, b)
                
                l_rad = np.deg2rad(l_range - 180)
                b_rad = np.deg2rad(b_range)
                
                ax.plot(l_rad, b_rad, 'gray', alpha=0.3, linewidth=0.5)
        
        if show_ecliptic:
            # Plot ecliptic
            ecl_lon = np.linspace(0, 360, 100)
            ecl_lat = np.zeros_like(ecl_lon)
            
            # Convert ecliptic to Galactic for plotting
            ecl_coords = SkyCoord(lon=ecl_lon*u.deg, lat=ecl_lat*u.deg, 
                                frame='barycentricmeanecliptic', 
                                obstime=self.observation_date)
            gal_coords = ecl_coords.transform_to('galactic')
            
            l_ecl = gal_coords.l.wrap_at(180*u.deg).deg
            b_ecl = gal_coords.b.deg
            
            # Convert for mollweide
            l_rad = np.deg2rad(l_ecl)
            b_rad = np.deg2rad(b_ecl)

            # sort coordinates to avoid wrapping
            table = np.array([l_rad, b_rad])
            table = table[:, table[0, :].argsort()]
            
            ax.plot(table[0], table[1], 'orange', linewidth=2, label='Ecliptic')

        if show_cvz:
            # Plot a circle of constant ecliptic latitude
            ecliptic_lon = np.linspace(0, 360, 100)

            for val in [-54, 54]:
                ecliptic_lat = np.full_like(ecliptic_lon, val)  # 54 degrees above the ecliptic plane

                # Create SkyCoord object for the circle
                ecliptic_coords = SkyCoord(lon=ecliptic_lon * u.deg, lat=ecliptic_lat * u.deg, frame='barycentricmeanecliptic')

                # Convert to Galactic
                gal_coords = ecliptic_coords.transform_to('galactic')

                # Convert for mollweide projection
                l_rad = np.deg2rad(gal_coords.l.wrap_at(180*u.deg).deg)
                b_rad = np.deg2rad(gal_coords.b.deg)

                # Plot the circle
                ax.plot(l_rad, b_rad, 'g--', linewidth=1, label=f'Constant Ecliptic Latitude ({val})')

        # Plot Sun
        if show_sun and self.sun_coord:
            sun_gal = self.sun_coord.transform_to('galactic')
            sun_l = sun_gal.l.wrap_at(180*u.deg).deg
            sun_b = sun_gal.b.deg
            sun_l_rad = np.deg2rad(sun_l)
            sun_b_rad = np.deg2rad(sun_b)
            ax.scatter(sun_l_rad, sun_b_rad, color='yellow', s=200,
                      marker='*', edgecolor='orange', linewidth=2, label='Sun')
        
        # Plot target
        if self.target_coord:
            target_gal = self.target_coord.transform_to('galactic')
            target_l = target_gal.l.wrap_at(180*u.deg).deg
            target_b = target_gal.b.deg
            target_l_rad = np.deg2rad(target_l)
            target_b_rad = np.deg2rad(target_b)
            ax.scatter(target_l_rad, target_b_rad, color='red', s=100,
                      marker='*', edgecolor='darkred', linewidth=2, label='Target')
        
        if target_list is not None:
            for target in target_list:
                target_gal = target.transform_to('galactic')
                l = target_gal.l.wrap_at(180*u.deg).deg
                b = target_gal.b.deg

                ax.scatter(np.deg2rad(l), np.deg2rad(b), color='green', s=50, marker='o',
                        edgecolor='darkgreen', linewidth=2)

        
        ax.set_xlabel('Galactic Longitude (deg)')
        ax.set_ylabel('Galactic Latitude (deg)')
        # ax.legend()
        ax.set_title('Sky Map (Mollweide Projection in Galactic Coordinates)')
        ax.grid(True, alpha=0.3)

    def _plot_sky_map(self, ax, show_ecliptic, 
                      show_sun, 
                      show_constellation_grid, 
                      target_list = None,
                      show_cvz=True):
        """
        Plot Mollweide projection sky map in ICRS coordinates.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes to plot on.
        show_ecliptic : bool
            If True, plot the ecliptic plane.
        show_sun : bool
            If True, plot the Sun's position.
        show_constellation_grid : bool
            If True, show the RA/DEC grid.
        target_list : list of SkyCoord, optional
            List of additional targets to plot.
        show_cvz : bool, optional
            If True, show the continuous viewing zone. Default is True.

        Returns:
        --------
        None
        """
        if show_constellation_grid:
            # RA/DEC grid
            ra_lines = np.arange(0, 360, 30)
            dec_lines = np.arange(-90, 91, 30)
            
            for ra in ra_lines:
                dec_range = np.linspace(-90, 90, 100)
                ra_range = np.full_like(dec_range, ra)
                
                # Convert to radians and adjust for mollweide
                ra_rad = np.deg2rad(ra_range - 180)  # Center at 0
                dec_rad = np.deg2rad(dec_range)
                
                ax.plot(ra_rad, dec_rad, 'gray', alpha=0.3, linewidth=0.5)
            
            for dec in dec_lines:
                ra_range = np.linspace(0, 360, 100)
                dec_range = np.full_like(ra_range, dec)
                
                ra_rad = np.deg2rad(ra_range - 180)
                dec_rad = np.deg2rad(dec_range)
                
                ax.plot(ra_rad, dec_rad, 'gray', alpha=0.3, linewidth=0.5)
        
        if show_ecliptic:
            # Plot ecliptic
            ecl_lon = np.linspace(0, 360, 100)
            ecl_lat = np.zeros_like(ecl_lon)
            
            # Convert ecliptic to RA/DEC for plotting
            ecl_coords = SkyCoord(lon=ecl_lon*u.deg, lat=ecl_lat*u.deg, 
                                frame='barycentricmeanecliptic', 
                                obstime=self.observation_date)
            icrs_coords = ecl_coords.transform_to('icrs')
            
            ra_ecl = icrs_coords.ra.wrap_at(180*u.deg).deg
            dec_ecl = icrs_coords.dec.deg
            
            # Convert for mollweide
            ra_rad = np.deg2rad(ra_ecl)
            dec_rad = np.deg2rad(dec_ecl)

            # sort coordinates to avoid wrapping
            table = np.array([ra_rad, dec_rad])
            table = table[:, table[0, :].argsort()]
            
            ax.plot(table[0], table[1], 'orange', linewidth=2, label='Ecliptic')

        if show_cvz:
            # Plot a circle of constant ecliptic latitude
            ecliptic_lon = np.linspace(0, 360, 100)
            ecliptic_lat = np.full_like(ecliptic_lon, 54)  # 54 degrees above the ecliptic plane

            # Create SkyCoord object for the circle
            ecliptic_coords = SkyCoord(lon=ecliptic_lon * u.deg, lat=ecliptic_lat * u.deg, frame='barycentricmeanecliptic')

            # Convert to ICRS
            icrs_coords = ecliptic_coords.transform_to('icrs')

            # Convert for mollweide projection
            ra_rad = np.deg2rad(icrs_coords.ra.wrap_at(180*u.deg).deg)
            dec_rad = np.deg2rad(icrs_coords.dec.deg)

            table = np.array([ra_rad, dec_rad])
            table = table[:,table[0,:].argsort()]

            # Plot the circle
            ax.plot(table[0], table[1], 'g--', linewidth=1, label='CVZ')

            # Plot a circle of constant ecliptic latitude
            ecliptic_lat = np.full_like(ecliptic_lon, -54)  # 54 degrees above the ecliptic plane

            # Create SkyCoord object for the circle
            ecliptic_coords = SkyCoord(lon=ecliptic_lon * u.deg, lat=ecliptic_lat * u.deg, frame='barycentricmeanecliptic')

            # Convert to ICRS
            icrs_coords = ecliptic_coords.transform_to('icrs')

            # Convert for mollweide projection
            ra_rad = np.deg2rad(icrs_coords.ra.wrap_at(180*u.deg).deg)
            dec_rad = np.deg2rad(icrs_coords.dec.deg)

            table = np.array([ra_rad, dec_rad])
            table = table[:,table[0,:].argsort()]

            # Plot the circle
            ax.plot(table[0], table[1], 'g--', linewidth=1)

        # Plot Sun
        if show_sun and self.sun_coord:
            sun_ra = self.sun_coord.ra.wrap_at(180*u.deg).deg
            sun_dec = self.sun_coord.dec.deg
            sun_ra_rad = np.deg2rad(sun_ra)
            sun_dec_rad = np.deg2rad(sun_dec)
            ax.scatter(sun_ra_rad, sun_dec_rad, color='yellow', s=200, 
                      marker='*', edgecolor='orange', linewidth=2, label='Sun')
        
        # Plot target
        if self.target_coord:
            target_ra = self.target_coord.ra.wrap_at(180*u.deg).deg
            target_dec = self.target_coord.dec.deg
            target_ra_rad = np.deg2rad(target_ra)
            target_dec_rad = np.deg2rad(target_dec)
            ax.scatter(target_ra_rad, target_dec_rad, color='red', s=100, 
                      marker='*', edgecolor='darkred', linewidth=2, label='Target')
        
        if target_list is not None:
            for target in target_list:
                ra = target.ra.wrap_at(180*u.deg).deg
                dec = target.dec.deg

                ax.scatter(np.deg2rad(ra), np.deg2rad(dec), color='green', s=50, marker='o', 
                        edgecolor='darkgreen', linewidth=2)

        
        ax.set_xlabel('RA (degrees)')
        ax.set_ylabel('Dec (degrees)')
        ax.legend()
        ax.set_title('Sky Map (Mollweide Projection)')
        ax.grid(True, alpha=0.3)

    def _plot_spacecraft_attitude(self, ax, target_list = None):
        """
        Create 3D visualization of spacecraft attitude and orientation.

        This method plots the spacecraft's body axes, solar array, and simplified
        body shape. It also shows the direction to the target and the Sun.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The 3D axes to plot on.
        target_list : list of SkyCoord, optional
            List of additional targets to plot.

        Returns:
        --------
        None
        """
                
        # Grab telescope axes
        x_axis = self.spacecraft_attitude[:, 0]
        y_axis = self.spacecraft_attitude[:, 1]
        z_axis = self.spacecraft_attitude[:, 2]

        # Create 3D plot
        # fig = plt.figure(figsize=(12, 10))
        # ax = fig.add_subplot(111, projection='3d')
        
        # Plot coordinate axes
        origin = np.array([0, 0, 0])
        
        # Telescope axes
        ax.quiver(*origin, *x_axis, color='red', arrow_length_ratio=0.1, linewidth=3, label='X (Target)')
        ax.quiver(*origin, *y_axis, color='green', arrow_length_ratio=0.1, linewidth=3, label='Y (Roll)')
        ax.quiver(*origin, *z_axis, color='blue', arrow_length_ratio=0.1, linewidth=3, label='Z (Solar Array)')
        
        # Sun
        sun_vector = get_vector(self.sun_coord)
        ax.quiver(*origin, *sun_vector, color='orange', arrow_length_ratio=0.1, linewidth=2, alpha=0.7)
        ax.scatter(*sun_vector, color='yellow', s=200, marker='*', 
                    edgecolor='orange', linewidth=2, label='Sun')
                
        # Target
        tgt_vector = get_vector(self.target_coord)
        ax.quiver(*origin, *tgt_vector, color='purple', arrow_length_ratio=0.1, linewidth=2, alpha=0.7)
        ax.scatter(*tgt_vector, color='red', s=100, marker='*', 
                    edgecolor='darkred', linewidth=2, label='Target Star')

        if target_list is not None:
            for target in target_list:
                tgt_vector = get_vector(target)
                ax.scatter(*tgt_vector, color='green', s=50, marker='o', 
                        edgecolor='darkgreen', linewidth=2)

        # Reference frame
        ax.quiver(*origin, 1, 0, 0, color='red', alpha=0.3, 
                 arrow_length_ratio=0.1, linestyle='--')
        ax.quiver(*origin, 0, 1, 0, color='green', alpha=0.3, 
                 arrow_length_ratio=0.1, linestyle='--')
        ax.quiver(*origin, 0, 0, 1, color='blue', alpha=0.3, 
                 arrow_length_ratio=0.1, linestyle='--')
        # Draw telescope body (simplified)
        # Draw telescope body (cylinder along X axis)
        length = 0.8  # Length of cylinder
        radius = 0.2*length  # Radius of cylinder
        resolution = 20  # Resolution of cylinder

        # Create cylinder
        theta = np.linspace(0, 2*np.pi, resolution)
        x = np.linspace(0, length, resolution)
        theta, x = np.meshgrid(theta, x)

        y = radius * np.cos(theta)
        z = radius * np.sin(theta)

        # Transform cylinder to align with x_axis
        cylinder_points = np.array([x.flatten(), y.flatten(), z.flatten()])
        transformed_points = np.dot(self.spacecraft_attitude, cylinder_points)

        ax.plot_surface(transformed_points[0].reshape(x.shape),
                        transformed_points[1].reshape(y.shape),
                        transformed_points[2].reshape(z.shape),
                        color='gray', alpha=0.5)

        # Draw solar array (perpendicular to Z axis)
        array_length = length*0.8
        array_size = length*0.8
        xx, yy = np.meshgrid([0, array_length*0.75], [-array_size/2, array_size/2])
        zz = np.full_like(xx, radius)  # Offset by radius along Z axis
        array_points = np.array([xx.flatten(), yy.flatten(), zz.flatten()])
        transformed_array = np.dot(self.spacecraft_attitude, array_points)

        ax.plot_surface(transformed_array[0].reshape(xx.shape),
                        transformed_array[1].reshape(yy.shape),
                        transformed_array[2].reshape(zz.shape),
                        color='blue', alpha=0.3)
        
        array_size = length*0.6
        xx, yy = np.meshgrid([array_length*0.75, array_length], [-array_size/2, array_size/2])
        zz = np.full_like(xx, radius)  # Offset by radius along Z axis
        array_points = np.array([xx.flatten(), yy.flatten(), zz.flatten()])
        transformed_array = np.dot(self.spacecraft_attitude, array_points)
                
        ax.plot_surface(transformed_array[0].reshape(xx.shape),
                        transformed_array[1].reshape(yy.shape),
                        transformed_array[2].reshape(zz.shape),
                        color='blue', alpha=0.3)

        # Formatting
        ax.set_xlim([-1.5, 1.5])
        ax.set_ylim([-1.5, 1.5])
        ax.set_zlim([-1.5, 1.5])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.legend()
        
    def _plot_info_panel(self, ax):
        """
        Plot information panel with numerical data.

        This method creates a text-based information panel displaying various
        mission parameters and spacecraft attitude information.

        Parameters:
        -----------
        ax : matplotlib.axes.Axes
            The axes to plot the information panel on.

        Returns:
        --------
        None
        """
        ax.axis('off')
        
        info_text = f"Observation Date: {self.observation_date.iso}\n\n"
        
        if self.target_coord:
            info_text += f"Target Coordinates:\n"
            info_text += f"  RA: {self.target_coord.ra.deg:.4f}° ({self.target_coord.ra.to_string(unit=u.hour, precision=2)})\n"
            info_text += f"  DEC: {self.target_coord.dec.deg:.4f}°\n\n"
            
            # Ecliptic coordinates
            ecl_lon, ecl_lat = self.get_ecliptic_coordinates(self.target_coord)
            info_text += f"Ecliptic Coordinates:\n"
            info_text += f"  Longitude: {ecl_lon:.4f}°\n"
            info_text += f"  Latitude: {ecl_lat:.4f}°\n\n"
            
            # Sun angle
            sun_angle = self.get_sun_angle()
            info_text += f"Sun Angle: {sun_angle:.2f}°\n"
            info_text += f"Pitch Angle: {self.pitch:.2f}°\n\n"
                        
            quat = self.get_attitude_quaternion()
            info_text += f"Quaternion [w,x,y,z]:\n"
            info_text += f"  [{quat[0]:.4f}, {quat[1]:.4f}, {quat[2]:.4f}, {quat[3]:.4f}]"
        
        ax.text(0.05, 0.95, info_text, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
        
        ax.set_title('Mission Parameters')

    def is_visible(self):
        """
        Check if the current target is within the pitch limits and thus visible.

        Returns:
        --------
        bool
            True if the target is within the pitch limits, False otherwise.
        """
        return self.pitch_limits[0] <= self.pitch <= self.pitch_limits[1]
        
    def print_roman_pointing_attributes(self):
        """
        Prints all attributes of the RomanPointing class cleanly.

        Parameters:
        -----------
        pointing: RomanPointing object
        """
        print("RomanPointing Attributes:")
        print(f"Observation Date: {self.observation_date}")
        print(f"Spacecraft Attitude:\n{self.spacecraft_attitude}")
        print(f"Target Coordinates: {self.target_coord}")
        print(f"Sun Coordinates: {self.sun_coord}")
        euler_angles = self.get_attitude_euler()
        print(f'RA:{euler_angles[0]:.4f} deg, DEC: {-euler_angles[1]:.4f} deg, V3PA: {-euler_angles[2]:.4f} deg')

def generate_power_law_sampling(n_points=14, range_max=36, power=3):
    """
    Uses power law transformation for strong edge emphasis
    Higher power = more edge clustering
    """
    # Generate symmetric points using power transformation
    half_points = n_points // 2
    
    # Create points from 0 to 1, then apply power transformation
    t = np.linspace(0, 1, half_points + 1)[1:-1]  # Exclude 0
    
    # Apply power transformation (higher power = more edge clustering)
    transformed = t ** power
    
    # Create positive side points
    pos_points = range_max - (range_max * transformed)
    
    # Create negative side points (mirror)
    neg_points = -pos_points[::-1]
    
    # Combine, ensuring we get exactly n_points
    if n_points % 2 == 0:
        all_points = np.concatenate([neg_points, pos_points])
    else:
        # Add center point for odd number
        all_points = np.concatenate([neg_points, [0], pos_points])
    
    return np.sort(all_points)


def main():
    # Set the start time for the observation
    t_start_str = ['2026-11-21T00:00:00.0']
    t_start = Time(t_start_str,format='isot', scale='utc')

    # Initialize the RomanPointing object
    pointing = RomanPointing(t_start)

    # Read target data from Excel file
    dat = pd.read_excel('RST_Fields.xlsx')
    lmc = dat[dat['Name']=='Euclid Deep Field North']
    tgt = SkyCoord(ra=lmc['RA'],dec=lmc['DEC'], frame='icrs', unit='deg')[0]
    # tgt = SkyCoord(ra=30.*u.deg,dec=45.*u.deg, frame='icrs', unit='deg')
    print(lmc['Name'],tgt)

    # Create a list of all targets from the Excel file
    tgt_list = [SkyCoord(ra=dat['RA'].iloc[i]*u.deg,dec=dat['DEC'].iloc[i]*u.deg, frame='icrs', unit='deg') for i in range(len(dat))]

    # Set the target for the pointing object
    pointing.set_target(tgt, roll=0.0)
    initial_pa = pointing.get_position_angle().value
    print('initial position angle:',initial_pa)

    # Calculate targets, pitch angles, and position angles for a range of pitch offsets
    target_list = []
    pitch_list = []
    pa_list = []
    initial_pitch = pointing.get_pitch_angle()
    print(initial_pitch)
    n_points=16
    # t = np.linspace(0, np.pi, n_points)
    # cosine_points = np.cos(t)
    range_max = 35.5
    # range_min = -35.5
    # mapped_points = range_min + (cosine_points + 1) * (range_max - range_min) / 2
    point_list = generate_power_law_sampling(n_points, range_max=range_max,power=4)

    for i in point_list:
        new_target, new_attitude = pointing.get_target_new_pitch(delta_pitch=float(i)-initial_pitch.value)
        pa_list.append(pointing.get_position_angle(new_attitude).value)
        pitch_list.append(pointing.get_pitch_angle(new_target).value)
        target_list.append(new_target)

    # Visualize the pointing
    pointing.visualize_pointing(show_ecliptic=True, 
                                show_sun=True, 
                                show_constellation_grid=True,
                                target_list=target_list
                                )

    # Calculate coordinates for a range of pitch and roll offsets
    coord_list = []
    n_points=16
    # t = np.linspace(0, np.pi, n_points)
    # cosine_points = np.cos(t)
    range_max = 35.5
    # range_min = -35.5
    # mapped_points = range_min + (cosine_points + 1) * (range_max - range_min) / 2
    point_list = generate_power_law_sampling(n_points, range_max=range_max,power=4)

    for i, dPitch in enumerate([val-initial_pitch.value for val in point_list]):
    #    for dRoll in np.linspace(-14.9, 14.9, 5):
    #        coord_list.append(pointing.delta_pitch_roll(dPitch,dRoll)+(dRoll,))
        ra,dec,pa, pitch = pointing.delta_pitch_roll(dPitch,0.0)
        coord_list.append([f"Target_{i}", 'Calibration', 'Stray light test', ra, dec, f'Pitch{pitch:2.4f}, PA={pa%360:3.4f}'])

    # Create a DataFrame with the calculated coordinates
    # df = pd.DataFrame(coord_list, columns=['RA', 'DEC', 'PA', 'pitch', 'roll'])
    df = pd.DataFrame(coord_list, columns=['Name', 'Category', 'Description', 'RA', 'DEC', 'Comments'])
    df.to_csv('Pitch-raster.csv',index=False, header=False)
    df

    print('=====Practice pysiaf=====')
    # Import pysiaf and related functions
    import pysiaf
    from pysiaf.utils.rotations import attitude

    # Initialize Roman SIAF
    rsiaf = pysiaf.Siaf('Roman')

    # Set boresight coordinates and initial position angle
    boresight_ra = tgt.ra.deg
    boresight_dec = tgt.dec.deg
    pa_v3 = initial_pa

    # Compare attitude matrices calculated by pysiaf and RomanPointing
    att = attitude(0, 0, boresight_ra, boresight_dec, pa_v3)
    np.testing.assert_allclose(att, pointing.spacecraft_attitude, rtol=1e-5, atol=1e-8)

    # Calculate attitude for a new position angle
    new_pa_v3 = 150.
    att = attitude(0, 0, boresight_ra, boresight_dec, new_pa_v3)

    # Compare with RomanPointing's Euler angle method
    p2 = RomanPointing(t_start)
    p2.set_target(tgt)
    att_euler = p2._calculate_pointing_attitude_euler(new_pa_v3)
    
    np.testing.assert_allclose(att, att_euler, rtol=1e-5, atol=1e-8)

    print('Matrix math was successfully checked against pysiaf')



if __name__ == "__main__":
    exit = main()







