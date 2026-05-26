=========
Changelog
=========

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