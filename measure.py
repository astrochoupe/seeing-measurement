#!/usr/bin/python
# coding=utf-8
__author__ = 'Didier Walliang'

from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
import numpy as np

# class definitions


class Slice:
    """Slice of a star trail"""

    def __init__(self, positions, intensities, approximate_fwhm_px = 5.0):
        """
        Construct a Slice object.

        :param positions: list of the position of each intensity in the image
        :param intensities: list of intensities in ADU
        :param approximate_fwhm_px: float, approximate FWHM in pixel to help to fit the gaussian
        :return:
        """

        self.positions = positions
        self.intensities = intensities
        self.approximative_fwhm_px = approximate_fwhm_px

        # move the y data to the axis to allow the fitting
        # (because the fitting doesn't work if the data are not near the x axis)
        min_intensity = min(self.intensities)
        self.relative_intensities = self.intensities - min_intensity

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
        approximate_stddev = self.approximative_fwhm_px / (2 * np.sqrt(2 * np.log(2)))

        gaussian_init = models.Gaussian1D(amplitude=max_relative_intensity, mean=mean_position, stddev=approximate_stddev)
        fit_gaussian = fitting.LevMarLSQFitter()
        self.gaussian = fit_gaussian(gaussian_init, self.positions, self.relative_intensities)

    def fwhm_from_gaussian(self):
        """FWHM calculation with the best-fit model"""
        sigma = self.gaussian.stddev.value
        self.fwhm_in_px = 2 * np.sqrt(2 * np.log(2)) * sigma

    def print_graph(self):
        # 200 dots from first to last position to have a smooth bell curve
        first_position = self.positions[0]
        last_position = self.positions[len(self.positions)-1]
        finer_positions = np.linspace(first_position, last_position, 200)

        plt.figure(1)
        plt.plot(self.positions, self.relative_intensities, 'ko')
        plt.plot(finer_positions, self.gaussian(finer_positions))
        plt.xlabel('Position (pixel)')
        plt.ylabel('Relative intensity (ADU)')
        plt.show()


class StarTrail:
    """A star trail"""

    def __init__(self, startrailcoordinates, img_data, sampling):
        """
        Construct a StarTrail object.

        :param startrailcoordinates: StarTrailCoordinates object
            Coordinates of the star trail in the image
        :param img_data: HDU object
        :param sampling: sampling in arcsec by pixel
        :return:
        """
        self.xmin = startrailcoordinates.xmin
        self.xmax = startrailcoordinates.xmax
        self.ymin = startrailcoordinates.ymin
        self.ymax = startrailcoordinates. ymax
        self.img_data = img_data
        self.sampling = sampling
        self.fwhms = np.array([])

    def calculate_fwhms(self):
        """Calcultate the FWHMs along the star trail"""

        self.fwhms = np.array([])

        for y in range(self.ymin, self.ymax):

            x_positions = np.array([])
            intensities = np.array([])

            for x in range(self.xmin, self.xmax):
                x_positions = np.append(x_positions, x)
                intensities = np.append(intensities, self.img_data[y, x])

            slice = Slice(x_positions, intensities)
            slice.fit_gaussian()
            slice.fwhm_from_gaussian()
            #slice.print_graph()

            fwhm_in_arcsec = self.sampling * slice.fwhm_in_px

            #print "FWHM in pixels : %f" % slice.fwhmInPx
            #print "FWHM in arcsec : %f" % fwhm_in_arcsec

            self.fwhms = np.append(self.fwhms, fwhm_in_arcsec)

        self.fwhm_samples = self.fwhms.size
        self.fwhm_min = np.min(self.fwhms)
        self.fwhm_max = np.max(self.fwhms)
        self.fwhm_mean = np.mean(self.fwhms)
        self.fwhm_median = np.median(self.fwhms)
        self.fwhm_sdt_dev = np.std(self.fwhms)

    def print_fwhms_results(self):
        """Print the min, max, mean, median and standard deviation of the FWHMs measurement"""

        print "Samples = %i " % self.fwhm_samples
        print "FWHM min = %f arcsec " % self.fwhm_min
        print "FWHM max = %f arcsec " % self.fwhm_max
        print "FWHM mean = %f arcsec " % self.fwhm_mean
        print "FWHM median = %f arcsec " % self.fwhm_median
        print "FWHM standard deviation = %f arcsec " % self.fwhm_sdt_dev

    def print_fwhms_graph(self):
        """Plot the FWHM with the best-fit model along the star trail"""

        plt.figure(1)
        plt.plot(self.fwhms, 'ko')
        plt.xlabel('Measure number')
        plt.ylabel('FWHM (arcsec)')
        plt.show()


class StarTrailCoordinates:

    def __init__(self, xmin, xmax, ymin, ymax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax


# ################
# Main program
# ################

# Open FITS file

print "Read FITS file"
hdulist = fits.open('/home/didier/Bureau/zenith-1.fits')
img_data = hdulist[0].data

startrailcoord1 = StarTrailCoordinates(xmin = 2035, xmax = 2070, ymin = 3275, ymax = 3750)
startrailcoord2 = StarTrailCoordinates(xmin = 2600, xmax = 2630, ymin = 1666, ymax = 2178)

startrail1 = StarTrail(startrailcoord1, img_data, sampling = 0.206)
startrail1.calculate_fwhms()
startrail1.print_fwhms_results()
startrail1.print_fwhms_graph()

startrail2 = StarTrail(startrailcoord2, img_data, sampling = 0.206)
startrail2.calculate_fwhms()
startrail2.print_fwhms_results()
startrail2.print_fwhms_graph()

# Close FITS file

hdulist.close()
