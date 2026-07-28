=========
Changelog
=========


v1.2.7
======

Added or Changed
----------------

- Now ROSALIA uses interpolation instead of reprojection for generating the scaled maps. This results on much faster and smoothed results. 



v1.2.6
======

Added or Changed
----------------

- Review of the core.straylight() method, to clean and parallelize the stray-light analysis per SCA. 
- Now the NDI files are stored in pickle python objects that contain the NDI RegularGridInterpolator. This is much faster than looping and generating them on each run. 


v1.2.0
======

Added or Changed
----------------

- Removing zodipy from dependencies. 
- Adding Ephessos as dependency.


v1.1.9
======

Added or Changed
----------------

- Now ROSALIA stray-light runs through a new class `rosalia.core.exposure`. Follow the tutorial under notebooks/tutorials/R1_ for more information. 


v1.1.2
======

Added or Changed
----------------

- Three superpixels from the 20 degree resolution NDI map were corrected. Previously, those values would show little to no stray-light compared to the neighbouring pixels. 
- New link for the rosalia-cache files, which were updated to include the corrected NDI maps.
- New (experimental) exposure class, that will superseed the current rosalia_stray, rosalia_zody, and rosalia_psf programs in the future. The new class is more flexible and allows for customized telescope definitions. The current exposure class will still be available for backward compatibility, but it is recommended to switch to the new class for future developments. 

v1.1.0
======

Added or Changed
----------------

- Fixed the units of the NDI maps. Now stray-light is lower in flux by a factor of 26.2 (corresponding to the area of the superpixels).
- Added tutorials and notebooks


v1.0.0
======

First release