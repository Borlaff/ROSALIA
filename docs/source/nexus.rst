Roman Research Nexus
=========

`The Roman Research Nexus <https://roman-docs.stsci.edu/data-handbook/roman-research-nexus>`_ is a web-based platform that provides access to a wide range of tools and resources for researchers working with data from the Nancy Grace Roman Space Telescope. The platform includes tools for data analysis, visualization, and collaboration, as well as access to documentation and support resources.

*ROSALIA is in the process to be integrated into the Roman Research Nexus. We will update this section as soon as the integration is complete. Suggestions are welcome!*

Installing ROSALIA in the Roman Research Nexus
-----------------------------

If your are setting up a new server in the Roman Research Nexus, you can install ROSALIA by following these steps:

1. Create a new conda environment for ROSALIA. Go to File → New Terminal from the JupyterLab menu bar, and run the following command:

.. code-block:: sh

       kernel-create rosalia 3.12 "rosalia"

2. Activate the new environment by running the following command:

.. code-block:: sh

       source kernel-activate rosalia

3. Install ROSALIA in the new environment by running the following command:

.. code-block:: sh

       pip install rosalia-wfi


4. Add the Astromatic and Gnuastro packages to the new environment by running the following command:

.. code-block:: sh

       mamba install conda-forge::astromatic-swarp conda-forge::gnuastro conda-forge::imagemagick


Using ROSALIA in the Roman Research Nexus
-----------------------------

We are actively preparing a series of tutorials to learn all the capabilities that ROSALIA has to exploit the scientific potential of the Roman WFI images. As of July 31st, 2026, we made available: 

1 - *Using ROSALIA to estimate stray-light in Nexus S3 stored WFI exposures*: The notebook is available at ``notebooks/tutorials/N1_Estimating_Straylight_ROSALIA_in_Nexus.ipynb``, and can be downloaded from the link below:
* :download:`R1_Rosalia_stray_example.ipynb <../../notebooks/tutorials/N1_Estimating_Straylight_ROSALIA_in_Nexus.ipynb>`. This tutorial shows how to compute the stray-light background for a given Roman level 2 exposure. The examples point to a Nexus cloud-stored mock image, so it might take some time to download if you run it in your own computer:

.. code-block:: bash
    jupyter notebook notebooks/tutorials/N1_Estimating_Straylight_ROSALIA_in_Nexus.ipynb


2 - *Using ROSALIA to find asteroids in WFI exposures*: The notebook is available at ``notebooks/tutorials/N2_Finding_Asteroids_ROSALIA_in_Nexus.ipynb``, and can be downloaded here:
* :download:`R1_Rosalia_stray_example.ipynb <../../notebooks/tutorials/N2_Finding_Asteroids_ROSALIA_in_Nexus.ipynb>`. This tutorial shows how to find out if any known Solar System Objects (i.e., asteroids) are present in your Roman data:

.. code-block:: bash
    jupyter notebook notebooks/tutorials/N2_Finding_Asteroids_ROSALIA_in_Nexus.ipynb






