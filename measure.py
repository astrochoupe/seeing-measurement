#!/usr/bin/python
# coding=utf-8
__author__ = 'Didier Walliang'

from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import numpy as np

# class definitions

class Slice():
    """Slice of a star trail"""

    def __init__(self, positions, intensities):
        self.positions = positions
        self.intensities = intensities

        # move the y data to the axis to allow the fitting
        # (because the fitting doesn't work if the data are not near the x axis)
        minIntensity = min(self.intensities)
        self.relative_intensities = self.intensities - minIntensity

    def print_positions(self):
        print "positions = %s " % self.positions

    def print_intensities(self):
        print "intensities = %s " % self.intensities

    def print_relative_intensities(self):
        print "relative intensities = %s " % self.relative_intensities

    def fit_gaussian(self):
        """Fit the data using a Gaussian"""
        max_relative_intensity = max(self.relative_intensities)
        mean_position = sum(self.positions)/len(self.positions)

        gaussian_init = models.Gaussian1D(amplitude=max_relative_intensity, mean=mean_position, stddev=3)
        fit_gaussian = fitting.LevMarLSQFitter()
        self.gaussian = fit_gaussian(gaussian_init, self.positions, self.relative_intensities)

    def fwhm_from_gaussian(self):
        """FWHM calculation with the best-fit model"""
        sigma = self.gaussian.stddev.value
        self.fwhmInPx = 2 * np.sqrt(2 * np.log(2)) * sigma

    def print_graph(self):
        # 200 dots from first to last position to have a smooth bell curve
        firstPosition = self.positions[0]
        lastPosition = self.positions[len(x_positions)-1]
        finerPositions = np.linspace(firstPosition, lastPosition, 200)

        plt.figure(1)
        plt.plot(self.positions, self.relative_intensities, 'ko')
        plt.plot(finerPositions, self.gaussian(finerPositions))
        plt.xlabel('Position (pixel)')
        plt.ylabel('Relative intensity (ADU)')
        plt.show()

class StarTrail():
    pass

# Open FITS file

print "Read FITS file"
hdulist = fits.open('/home/didier/Bureau/zenith-1.fits')
imgData = hdulist[0].data

# Read star trail's intensity

fwhms = np.array([])

"""
xmin = 2070
xmax = 2110
ymin = 3212
ymax = 3700
"""

xmin = 2035
xmax = 2070
ymin = 3275
ymax = 3750

for y in range(ymin, ymax):

    x_positions = np.array([])
    intensities = np.array([])

    for x in range(xmin, xmax):
        x_positions = np.append(x_positions, x)
        intensities = np.append(intensities, imgData[y, x])

    slice = Slice(x_positions, intensities)
    slice.fit_gaussian()
    slice.fwhm_from_gaussian()
    #slice.print_graph()

    sampling = 0.206  # arcsec by pixel
    fwhmInArcsec = sampling * slice.fwhmInPx

    #print "FWHM in pixels : %f" % slice.fwhmInPx
    print "FWHM in arcsec : %f" % fwhmInArcsec

    fwhms = np.append(fwhms, fwhmInArcsec)

# Plot the data with the best-fit model

# min, max, avg FWHM
print "FWHM min = %f " % np.min(fwhms)
print "FWHM max = %f " % np.max(fwhms)
print "FWHM mean = %f " % np.mean(fwhms)
print "FWHM median = %f " % np.median(fwhms)

plt.figure(2)
plt.plot(fwhms, 'ko')
plt.xlabel('Measure number')
plt.ylabel('FWHM (arcsec)')
plt.show()

# Close FITS file

hdulist.close()
