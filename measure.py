#!/usr/bin/python
# coding=utf-8
__author__ = 'Didier Walliang'

from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import numpy as np

# Open FITS file

print "Read FITS file"
hdulist = fits.open('/home/didier/Bureau/zenith-1.fits')
imgData = hdulist[0].data

# Read star trail's intensity

xmin = 2070
xmax = 2110
y = 3222
plotx = []
ploty = []
for x in range(xmin, xmax):
    plotx.append(x)
    ploty.append(imgData[y, x])
    print "Intensity of the pixel %d,%d : %d" % (x, y, imgData[y, x])

# Convert lists to numpy arrays

plotx = np.array(plotx)
ploty = np.array(ploty)

# move the y data to the axis to allow the fitting
# (because the fitting doesn't work if the data are not near the x axis)
ploty -= min(ploty)

# Print data to verify the correctness of the previous code

print "plotx = %s " % plotx
print "ploty = %s " % ploty

# Fit the data using a Gaussian

gaussian_init = models.Gaussian1D(amplitude=max(ploty), mean=sum(plotx)/len(plotx), stddev=3)
fit_gaussian = fitting.LevMarLSQFitter()
gaussian = fit_gaussian(gaussian_init, plotx, ploty)

print gaussian

# FWHM calculation with the best-fit model

sigma = gaussian.stddev.value
fwhmInPx = 2 * np.sqrt(2 * np.log(2)) * sigma
sampling = 0.206  # arcsec by pixel
fwhmInArcsec = sampling * fwhmInPx

print "Standard deviation : %f" % sigma
print "FWHM in pixels : %f" % fwhmInPx
print "FWHM in arcsec : %f" % fwhmInArcsec

# Plot the data with the best-fit model

# 200 dots from plotx[0] to plotx[len(plotx)-1] to have a smooth bell curve
finerPlotx = np.linspace(plotx[0], plotx[len(plotx)-1], 200)

plt.figure(1)
plt.plot(plotx, ploty, 'ko')
plt.plot(finerPlotx, gaussian(finerPlotx))
plt.xlabel('Position (pixel)')
plt.ylabel('Intensite (ADU)')
plt.show()

# Close FITS file

hdulist.close()
