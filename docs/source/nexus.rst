Roman Research Nexus
=========

`The Roman Research Nexus <https://roman-docs.stsci.edu/data-handbook/roman-research-nexus>`_ is a web-based platform that provides access to a wide range of tools and resources for researchers working with data from the Nancy Grace Roman Space Telescope. The platform includes tools for data analysis, visualization, and collaboration, as well as access to documentation and support resources.

*ROSALIA is in the process to be integrated into the Roman Research Nexus. We will update this section as soon as the integration is complete. Suggestions are welcome!*

How to use ROSALIA in the Roman Research Nexus
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

       mamba install conda-forge::astromatic-swarp conda-forge::gnuastro

