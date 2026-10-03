#%%
import matplotlib.pyplot as plt
import numpy as np
import rosalia as rs
import astropy.io.fits as fits

exposure = rs.core.exposure('/Users/pmsanch1/storage/ROMAN/comissioning/PID1020/r0102001002001005001_0004*f158_cal.asdf')

psf_canvas = exposure.psf(n_stars=20, n_workers=10, method='superback')


fits.PrimaryHDU(psf_canvas, header=exposure.header).writeto('psf_canvas.fits', overwrite=True)

assert psf_canvas is not None
assert psf_canvas.shape[0] > 0 and psf_canvas.shape[1] > 0
assert np.sum(psf_canvas) > 0


plt.imshow(23.9 - 2.5*np.log10(psf_canvas[::10, ::10]) + 5*np.log10(0.108), origin='lower', cmap='nipy_spectral')
plt.colorbar()
plt.show()
#%%