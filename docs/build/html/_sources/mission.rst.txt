What is ROSALIA?
================

The primary mission of *ROSALIA* is to provide an estimate of the sky background level of an image obtained with a space telescope. Most researchers combine multiple images obtained with telescopes to generate *deeper* images (more precise). However, due to tiny - but noticeable - differences in the conditions of the observation (i.e., zodiacal light level, stray-light contamination, thermal background) the baseline level of the images (empty regions that we would expect to be zero) can present a certain level of light that must be removed. This is called **sky background correction** and it is one of the most critical steps in astronomical imaging.

.. image:: https://raw.githubusercontent.com/Borlaff/ROSALIA/main/images/rosalia_star_on_the_edge.gif
  :alt: NASA Ames Research Center - Code S - RACCOONS team

Until now, most astronomers attempted to model and correct the background level using the science images themselves. However, this method has a number of problems associated, like *sky background oversubtraction*. **ROSALIA** takes a different approach: Use all the information available from other telescopes to predict the sky background level as a combination of stray-light, zodiacal light, and thermal background.
