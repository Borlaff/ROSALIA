Roman Commissioning Stray-light Activities
=====================================================

This section is dedicated to the Roman Commissioning Activity Request CAR171, where we focus on the calibration of the stray-light rejection efficiency of the Roman Space Telescope Wide Field Instrument. To do that, we will illuminate specific regions of the telescope using very bright stars. Some of the selected regions are highly sensitive to stray-light, and are predicted to produce characteristic patterns of background stray-light in the WFI images. 

Selecting the stray-light regions
----------------------------------

The following Jupyter-notebook details the selection of the main stray-light sensitive regions along with the scanning dithering patterns defined to observationally determine the location and shape of the Normalized Detector Irradiance in those orientations. 

* :download:`R1_CAR171_targets.ipynb <../../notebooks/CAR171/R1_CAR171_targets.ipynb>`


.. image:: https://raw.githubusercontent.com/Borlaff/ROSALIA/main/notebooks/CAR171/ndi_2deg_CARs_with_points.png
  :alt: Selection of stray-light sensitive locations near the focal plane. 

**Figure 1**: Selection of stray-light sensitive locations near the focal plane. 

.. image:: https://raw.githubusercontent.com/Borlaff/ROSALIA/main/notebooks/CAR171/ndi_20deg_CARs_with_points.png
  :alt:  Selection of stray-light sensitive locations outside the focal plane. 

**Figure 2**: Selection of stray-light sensitive locations outside the focal plane. 


.. image:: https://raw.githubusercontent.com/Borlaff/ROSALIA/main/notebooks/CAR171/ndi_20deg_CAR171_RoguePath_dithering.png
  :alt: Selection of the dithering pattern for the "Rogue Path" area (WFIV3+ positive orientation).

**Figure 3**: Selection of the dithering pattern for the "Rogue Path" area (WFIV3+ positive orientation).

.. image:: https://raw.githubusercontent.com/Borlaff/ROSALIA/main/notebooks/CAR171/ndi_20deg_CAR171_OUT_1_dithering.png
  :alt: Selection of the dithering pattern for the "Roman Channel" area (secondary leg scattering).

**Figure 4**: Selection of the dithering pattern for the "Roman Channel" area (secondary leg scattering).

.. image:: https://raw.githubusercontent.com/Borlaff/ROSALIA/main/notebooks/CAR171/ndi_20deg_CAR171_IN_4_5_Dragon_dithering.png
  :alt: Selection of the dithering pattern for the "Roman Dragon's Breath" region (SCA lightshields).

**Figure 5**: Selection of the dithering pattern for the "Roman Dragon's Breath" region (SCA lightshields).


Selecting the calibration sources and sky pointings
--------------------------------------------------------------------

Once the telescope orientations have been selected, the next step is to define a set of realistic sky observations that can illuminate the required region of the observatory with a source. We selected three stars are calibration sources, Canopus, Eta Uma, and HD128998. The observatory has to orient itself with a certain position angle and angular distance from the center of the WFI focal plane to illuminate the required regions, while complying with the Sun-avoidance orientation angle requirements for Roman Space Telescope. Estimating a valid pointing to position Roman at a certain position angle and distance from a source can be done easily with **ROSALIA** (see Additional Tools / Planning Offset Observations). The estimation for all the Commissioning Activity Request 171 (targets above) is performed in the following Jupyter-notebook. 

* :download:`R2_CAR171_simulations.ipynb <../../notebooks/CAR171/R2_CAR171_simulations.ipynb>`


Updating the Astronomer Proposal Tool submission file
--------------------------------------------------------------------

Since the orientation with the Sun changes with time, so does the optimal position angle of Roman Space Telescope. That imposes a limit for any observation to be valid. This is particularly critical in our observations because they need to have an exact orientation angle with respect to the source. To keep scheduling as flexible as possible, we generated a new notebook that automatically updates a template APT file with the pointings estimated in *R2_CAR171_simulations.ipynb*. 

* :download:`R3_CAR171_apt_generator.ipynb <../../notebooks/CAR171/R3_CAR171_apt_generator.ipynb>`
