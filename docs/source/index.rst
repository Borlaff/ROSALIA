Welcome to Roman / ROSALIA's documentation!
===========================================

.. image:: ../../images/rosalia_logo.png
  :width: 700
  :alt: The ROSALIA WFS Project logo 

**ROSALIA** (*Roman Sky Analyst for Low surface brightness Imaging & Astronomy*)
is a NASA pipeline to model the sky background level on astronomical images obtained
with **Nancy Grace Roman Space Telescope** and its direct predecessor, the
legendary **Hubble Space Telescope**. In particular *ROSALIA* is focused on the
prediction and calibration of stray-light in the Roman Wide Field Instrument,
one of the main contaminants in ultra deep low surface brightness observations,
and the main source of gradients of parasitic light for space telescopes.

*ROSALIA* combines the information from existing photometric catalogs (Gaia,
2MASS, WISE) with precise optical and payload ray-tracing models of the Roman
Space Telescope, allowing to generate images of stray-light and other components
of the sky-background for user-defined observational conditions.

*ROSALIA* is funded through a NASA Grant (D.14 Roman 2022), *ROSES/Nancy Grace
Roman Space Telescope Research and Support Participation Opportunities.*

Check out the :doc:`usage` section for further information how to use ROSALIA to analyze your images. In
:ref:`install <installation>` we describe how to set up ROSALIA in your machine.

.. note::

   Space is hard, and this project is under active development. The predicted sky background levels by ROSALIA might not be exact and it is not recommended to use it blindly. But we keep pushing everyday to make it better, and our team is happy to assist with your particular project. We will be including more functions and tutorials in the very near future. **ROSALIA** is a Roman Wide Field Science project that needs your help to improve! If you find a useful capability that could be added, a problem to be resolved, or any other idea, :doc:`contact us <whoarewe>`.

.. autosummary::
   :toctree: _autosummary
   :template: custom-module-template.rst
   :recursive:

Contents
--------

.. toctree::
    :maxdepth: 2

    mission
    installation
    measuringstray
    zodiacal
    PSF
    SSOs
    tutorials
    usage
    nexus
    whoarewe
    rosalia
    api
