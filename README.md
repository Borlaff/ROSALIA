<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a id="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/Borlaff/ROSALIA">
    <img src="https://raw.githubusercontent.com/Borlaff/ROSALIA/main/images/rosalia_logo.png" alt="ROSALIA_logo" width="710" height="437">
  </a>

  <p align="center">
    <h2 align="center">ROSALIA: ROman Sky Analyst for Low surface brightness Imaging & Astronomy</h2>
    <br />
    <a href="https://rosalia.readthedocs.io/en/latest/index.html"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/Borlaff/ROSALIA">View Demo</a>
    ·
    <a href="https://github.com/Borlaff/ROSALIA/issues/new?labels=bug&template=bug-report---.md">Report Bug</a>
    ·
    <a href="https://github.com/Borlaff/ROSALIA/issues/new?labels=enhancement&template=feature-request---.md">Request Feature</a>
  </p>
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project

ROSALIA (Roman Sky Analyst for Low surface brightness Imaging & Astronomy) is a pipeline to model the sky background level on astronomical images obtained with [NASA/Nancy Grace Roman Space Telescope](https://roman.gsfc.nasa.gov) and its direct predecessor, the legendary [NASA/Hubble Space Telescope](https://science.nasa.gov/mission/hubble/). In particular ROSALIA is focused on the prediction and calibration of *stray-light* in the [Roman Wide Field Instrument](https://roman-docs.stsci.edu/roman-instruments-home/wfi-imaging-mode-user-guide/introduction-to-wfi), one of the main contaminants in ultra deep low surface brightness observations, and the main source of gradients of parasitic light for space telescopes. ROSALIA combines the information from existing photometric catalogs (Gaia, 2MASS, WISE) with precise optical and payload ray-tracing models of the Roman Space Telescope, allowing to generate images of stray-light and other components of the sky-background for user-defined observational conditions.

ROSALIA is funded through a NASA Grant (D.14 Roman 2022), ROSES/Nancy Grace Roman Space Telescope Research and Support Participation Opportunities.

Sci-PI: Alejandro S. Borlaff (NASA ARC). Admin-PI: Pamela M. Marcum (NASA ARC)

<a href="https://github.com/Borlaff/ROSALIA">
    <img src="https://raw.githubusercontent.com/Borlaff/ROSALIA/main/images/rosalia_star_on_the_edge.gif" alt="ROSALIA_animation" width="1500" height="600">
</a>


<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!--
### Built With

* [![Next][Next.js]][Next-url]
* [![React][React.js]][React-url]
* [![Vue][Vue.js]][Vue-url]
* [![Angular][Angular.io]][Angular-url]
* [![Svelte][Svelte.dev]][Svelte-url]
* [![Laravel][Laravel.com]][Laravel-url]
* [![Bootstrap][Bootstrap.com]][Bootstrap-url]
* [![JQuery][JQuery.com]][JQuery-url]

<p align="right">(<a href="#readme-top">back to top</a>)</p>
 -->

## Note: 

This package is under active development. Our team is working hard to provide a reliable pipeline that allows to model the sky backgrounds for Roman Space Telescope Wide Field Instrument. However, progress does not come overnight. We are open-source testing code and substantial changes will occur regularly. Use at your own risk. In case of doubt, send an email to Alejandro S. Borlaff, a.s.borlaff@nasa.gov.  

## Installation
### Managing dependencies
ROSALIA is based on multiple packages, including [Astropy](https://www.astropy.org/), [Astroquery](https://astroquery.readthedocs.io/en/latest/), and [Romanisim](https://romanisim.readthedocs.io/en/latest/), [NumPy](https://numpy.org/), [SciPy](https://scipy.org/), and [Matplotlib](https://matplotlib.org/) among many others. The easiest way to install all the dependencies is through a package manager like [Conda](https://anaconda.org/anaconda/conda) or [Mamba](https://github.com/mamba-org/mamba). If you have a Conda/Mamba package manager already installed in your system, skip to the following section. If you do not have a package manager, follow the Conda installation instructions at the [Space Telescope stenv environment webpage](https://stenv.readthedocs.io/en/latest/getting_started.html).

### Installing ROSALIA
Create a clean environment for ROSALIA
```sh
conda create -n rosalia python=3.12 conda-forge::astromatic-swarp
```

After the new environment is created, we can activate it.
```sh
conda activate rosalia
   ```

Once in a clean conda environment, we can install ROSALIA. The preferred method to install it is through pip.

```sh
pip install rosalia
   ```
<p align="right">(<a href="#readme-top">back to top</a>)</p>

ROSALIA needs a set of calibrations files to work. ROSALIA needs a series of cache files to work. If they are not installed, most of the functions will work, but you won't be able to estimate the surface brightness of the stray light, which is one of the main functionalities of ROSALIA. To install the cache files, download the following folder:

`ROSALIA Cache <https://zenodo.org/records/18882110/files/rosalia-cache.tar.gz>`_,

Most functions will work without it, but the main ones (rosalia_stray) will return an error when executed if these files are not found. The ROSALIACACHE folder must be defined in the environment as: 

```sh
export ROSALIACACHE=/home/user/project/rosalia_cache
   ```
<p align="right">(<a href="#readme-top">back to top</a>)</p>

Add this line to your .bashrc or .zshrc to set it up by default.

That is it! We are ready to start analyzing Space Telescope images.

<!-- USAGE EXAMPLES -->
### Minimal Use Example - Simulating stray-light

**I want to calculate the stray-light background on an image an I want it now!**

Alright, alright! ROSALIA can generate a quick simulation of stray-light background for Roman / WFI, provided the coordinates of the center of WFI, position angle, date, bandpass, and exposure time: 

```python
import rosalia as rs 
from astropy.time import Time

# First define a Roman Space Telescope WFI exposure basic parameters.
ra = 123  # Right ascension at the center of the FOV, in degrees. 
dec = 45  # Declination at the center of the FOV, in degrees.
PA = 67   # Position angle, counter-clockwise from North, in degrees.
date = Time("2026-10-01T00:00:00")  # Date of the observation, in Astropy Time format.
bandpass = "F129"  # A string with the bandpass name for WFI. See https://roman.gsfc.nasa.gov/science/WFI_technical.html
exptime = 600  # Exposure time, in seconds.

rosalia_stray = rs.correct.rosalia_stray(ra=ra, dec=dec, PA=PA, date=date, 
                                         bandpass=bandpass, exptime=exptime)
```

The estimated background stray-light model will be stored in a new FITS file, named using the input parameters: 

```sh
   WFI_F129_RA_123.000_DEC_023.000_MJD_60462.00000_PA_045.00_stray.fits
```

## Not-So-Minimal Use Examples

### Simulating Zodiacal light
ROSALIA estimates the amount of stray-light from Roman Space Telescope images. To do this, it calculates how many photons reach the focal plane array from secondary optical paths, based on a function called Normalized Detector Irradiance (NDI).

Those photons represent a source of contamination and typically must be modeled and removed before the images are ready for science. ROSALIA calculates the flux of photons for each pixel of the focal plane array. For Roman/WFI, that is a total of 300,811,392 pixels! (18 4088x4088 H4RG-10 detectors).

Let's do a quick example to figure out how many photons do we expect to see on an average Roman / WFI exposure. The main source of background light (under normal conditions) is the Zodiacal light. 

```python
import rosalia as rs 
from astropy.time import Time

# First define a Roman Space Telescope WFI exposure basic parameters.
ra = 123  # Right ascension at the center of the FOV, in degrees. 
dec = 45  # Declination at the center of the FOV, in degrees.
PA = 67   # Position angle, counter-clockwise from North, in degrees.
date = Time("2026-10-01T00:00:00")  # Date of the observation, in Astropy Time format.
bandpass = "F129"  # A string with the bandpass name for WFI. See https://roman.gsfc.nasa.gov/science/WFI_technical.html
exptime = 600  # Exposure time, in seconds.

rosalia_zody = rs.correct.rosalia_zody(ra=ra, dec=dec, PA=PA, date=date, 
                                        bandpass=bandpass, exptime=exptime)
```

The code above will generate a FITS file containing 18 SCI extensions, each of them with the predicted Zodiacal light flux in electrons per second (e/s). The file will be automatically named following following the input ROSALIA convention, in this case:

```sh
   WFI_F129_RA_123.000_DEC_023.000_MJD_60462.00000_PA_045.zody.fits
```

Inspecting the FITS file, we can see that the average value from the model is ~0.62, the units are electrons per second per pixel. So, we expect to have a Zodiacal light background of ~0.62 e/s in this particular location of the sky and time. What is the surface brightness magnitude of such an electron flux? We can transform it using ROSALIA as well. 

```python
  In [2]: rs.detectors.fe2mu(0.62,instrument="WFI",filter_name="F129",telescope="RST")
  SVO reported the transmission in m2. Correcting by area of the telescope (4.523893421169302 m2 )
  Out[2]: array([22.20263638])
```

So, 0.62 electrons per second correspond to roughly 22.2 mag/arcsec^2. That is very much a reasonable value for the Zodiacal light background for a Optical / NIR space telescope like Roman or Euclid [(Borlaff et al. 2022)](https://ui.adsabs.harvard.edu/abs/2022A%26A...657A..92E/abstract). Now, you can use rs.detectors.fe2mu and rs.detectors.mu2fe to transform back and forth between magnitudes per arcsec^2 and electrons per second per pixel. 


### Simulating PSF (stars)

Another use of ROSALIA is to simulate stars that actually land inside the field of view of Roman. These are modeled by rosalia_psf. 

```python

import rosalia as rs
import pandas as pd
from astropy.time import Time

ra = 123 # Right ascension, in degrees. 
dec = 23 # Declination, in degrees.
PA = 45 # Position angle, in degrees.
date = Time("2024-06-01T00:00:00") # Date of the observation, in Astropy Time YYYY-MM-DDTHH:MM:SS format.
bandpass = "F129"
g_mag_max = 15 # Maximum Gaia g-band magnitude considered when searching stars in the Gaia / 2MASS / WISE database. 
exptime = 600 # Exposure time, in seconds.

catalog = {"ra": [123.01, 123.02, 123.03],
           "dec": [23.01, 23.02, 23.03],
           "mag_lambda": [14, 14.5, 15],
           "source_id": [1, 2, 3],
           "cat_id": [1, 2, 3]}
           
input_catalog = pd.DataFrame(catalog) # Optional: By defining this catalog manually, we override the automatic query in Gaia / 2MASS / WISE. Useful for quick tests. 

rosalia_psf = rs.correct.rosalia_psf(ra=ra, dec=dec, PA=PA, date=date, 
                                     bandpass=bandpass, exptime=exptime, 
                                     input_catalog=input_catalog,
                                     g_mag_max=15, verbose=1)
```


<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- ASDF -->
## What is ASDF? Where are the FITS files?
ROSALIA uses FITS files instead of ASDF files. [ASDF](https://asdf.readthedocs.io/en/latest/) is the successor of [FITS](https://www.stsci.edu/hst/wfpc2/Wfpc2_dhb/intro_ch23.html) format and has been adopted since JWST. However, as of 2026, GUI visualizers like [SAODS9 are not yet compatible with ASDF](https://www.stsci.edu/files/live/sites/www/files/home/jwst/science-planning/user-committees/jwst-users-committee/_documents/jstuc-0919-data-analysis-tool-ferguson.pdf). In case that you are using ASDF, ROSALIA provides an easy way to extract most useful information from the ASDF files through `exposure-inspector`:
```sh
   exposure-inspector RST_WFI_ROSALIA_test_Orion_Belt_SCAWFI01.asdf
```
`exposure-inspector` will print a series of fields containing basic information from the ASDF tree, including the name of the telescope, instrument, detector, and filter, pointing information like right ascension and declination, transmission curve of the filter, and the WCS of the header.


<!-- ROADMAP -->
## Roadmap
- [x] Automatic queries of catalogs of bright sources.
    - [x] Gaia, 2MASS, WISE
    - [x] Horizons/JPL Solar System Objects
- [x] Retrieval of stray-light blocking efficiency from ray-tracing models
- [x] Ingestion of ASDF Roman/WFI simulated files (i.e., https://romanisim.readthedocs.io/en/latest/)
- [x] Add diffraction modelling to Roman/WFI.
- [ ] Add thermal emission model (internal stray-light).
- [ ] Complete support for Hubble Space Telescope ACS & WFC3/IR.
- [ ] Automatic identification of SSOs in Roman/WFI observations.

See the [open issues](https://github.com/github_username/repo_name/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- CONTRIBUTING -->
## Contributing
Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make ROSALIA better, there are two options.
1. You can also simply open an issue with the tag "enhancement" (see below).
2. Fork the repo and create a pull request.
3. Email the PI's of the project (a.s.borlaff@nasa.gov) with your ideas.

Instructions for fork/pull contributions.
1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Top contributors:
<a href="https://github.com/Borlaff">
  <img src="https://contrib.rocks/image?repo=Borlaff/Elbow" />
</a>
<a href="https://github.com/pmarcum">
  <img src="https://contrib.rocks/image?repo=pmarcum/Bibman" />
</a>
<a href="https://github.com/ScottRohrbach">
  <img src="https://avatars.githubusercontent.com/u/190406686?v=4" />
</a>


<!-- LICENSE -->
## License

ROSALIA © 2025Pπ by Alejandro S. Borlaff is licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- CONTACT -->
## Contact

Alejandro S. Borlaff - [@asborlaff](https://bsky.app/profile/asborlaff.bsky.social) - a.s.borlaff@nasa.gov

Project Link: [https://github.com/Borlaff/ROSALIA](https://github.com/Borlaff/ROSALIA)

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* []()
* []()
* []()

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/github_username/repo_name.svg?style=for-the-badge
[contributors-url]: https://github.com/Borlaff/ROSALIA/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/github_username/repo_name.svg?style=for-the-badge
[forks-url]: https://github.com/Borlaff/ROSALIA/network/members
[stars-shield]: https://img.shields.io/github/stars/github_username/repo_name.svg?style=for-the-badge
[stars-url]: https://github.com/Borlaff/ROSALIA/stargazers
[issues-shield]: https://img.shields.io/github/issues/github_username/repo_name.svg?style=for-the-badge
[issues-url]: https://github.com/Borlaff/ROSALIA/issues
[license-shield]: https://img.shields.io/github/license/github_username/repo_name.svg?style=for-the-badge
[license-url]: https://github.com/Borlaff/ROSALIA/blob/main/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://www.linkedin.com/in/alejandro-borlaff/
[product-screenshot]: images/moving_star_v3.gif
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com
