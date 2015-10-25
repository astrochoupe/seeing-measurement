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

x = 3222
ymin = 2070
ymax = 2110
plotx = []
ploty = []
for y in range(ymin,ymax) :
    plotx.append(y)
    ploty.append(imgData[x,y])
    print "Intensity of the pixel %d,%d : %d" % (x, y, imgData[x,y])

# Convert lists to numpy arrays
plotx = np.array(plotx)
ploty = np.array(ploty)

#
ploty -= min(ploty)

# Print data to verify the correctness of the previous code

print plotx
print ploty

# Fit the data using a Gaussian

gaussian_init = models.Gaussian1D(amplitude=max(ploty), mean=sum(plotx)/len(plotx), stddev=3)
fit_gaussian = fitting.LevMarLSQFitter()
gaussian = fit_gaussian(gaussian_init, plotx, ploty)

print gaussian

# FWHM calculation with the best-fit model

sigma =  gaussian.stddev.value
fwhmInPx = 2 * np.sqrt(2 * np.log(2)) * sigma
sampling = 0.206 # arcsec by pixel
fwhmInArcsec = sampling * fwhmInPx

print "Standard deviation : %f" % sigma
print "FWHM in pixels : %f" % fwhmInPx
print "FWHM in arcsec : %f" % fwhmInArcsec

# Plot the data with the best-fit model

plotxx = np.linspace(plotx[0], plotx[len(plotx)-1], 200)

plt.figure(1)
plt.plot(plotx, ploty, 'ko')
plt.plot(plotxx, gaussian(plotxx))
plt.xlabel('Position (pixel)')
plt.ylabel('Intensite (ADU)')
plt.show()

# Close FITS file

hdulist.close()
