#!/usr/bin/python
# coding=utf-8
__author__ = 'Didier Walliang'

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import dates
from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy import units as u
from astropy.time import Time
from photutils import segment_properties
from photutils import detect_sources
from photutils import detect_threshold

# class definitions


class Slice:
    """Slice of a star trail"""

    def __init__(self, positions, intensities, approximate_fwhm_px = 6.0):
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
        self.fwhm_stddev = np.std(self.fwhms)

    def print_fwhms_results(self):
        """Print the min, max, mean, median and standard deviation of the FWHMs measurement"""

        print "Samples = %i " % self.fwhm_samples
        print "FWHM min = %f arcsec " % self.fwhm_min
        print "FWHM max = %f arcsec " % self.fwhm_max
        print "FWHM mean = %f arcsec " % self.fwhm_mean
        print "FWHM median = %f arcsec " % self.fwhm_median
        print "FWHM standard deviation = %f arcsec " % self.fwhm_stddev

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

    def __str__(self):
        return "xmin = %i ; xmax = %i ; ymin = %i ; ymax = %i" % (self.xmin, self.xmax, self.ymin, self.ymax)


class TrailsImage:
    """An image containing star trails"""

    def __init__(self, data, sampling, target_fwhm_arcsec, length_x_axis):
        """
        Construct a TrailsImage object.

        :param data: HDU object
        :param sampling: sampling in arcsec by pixel
        :param target_fwhm_arcsec: estimated FWHM in arcsec (to help the calcultation)
        :param length_x_axis: length of X axis
        :return:
        """
        self.data = data
        self.sampling = sampling
        self.target_fwhm_arcsec = target_fwhm_arcsec
        self.length_x_axis = length_x_axis

    def search_trails(self):
        """Search star trails in image"""
        threshold = detect_threshold(self.data, snr=1)
        sigma = 2.0 * gaussian_fwhm_to_sigma
        kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
        self.segments = detect_sources(self.data, threshold, npixels=1000, filter_kernel=kernel)
        self.segments_properties = segment_properties(self.data, self.segments)

    def nb_trails(self):
        """
        Give the number of trails detected

        :return: number of trails detected
        """
        return len(self.segments_properties)

    def calculate_fwhm(self):
        """Measure the FWHM of the star trails"""

        self.trails_fwhm = np.array([])
        self.trails_stddev = np.array([])

        # for each star trail
        for properties in self.segments_properties:
            # x coordinate of the middle of the trail
            xmiddle = (properties.xmax + properties.xmin)/2

            target_fwhm_px = self.target_fwhm_arcsec/self.sampling

            xmin = int(math.floor(xmiddle.value - target_fwhm_px*2))
            xmax = int(math.ceil(xmiddle.value + target_fwhm_px*2))

            if(xmin < 1 or xmax > self.length_x_axis):
                continue

            ymin = int(properties.ymin.value)
            ymax = int(properties.ymax.value)

            startrailcoord = StarTrailCoordinates(xmin, xmax, ymin, ymax)

            #print startrailcoord

            startrail = StarTrail(startrailcoord, self.data, self.sampling)
            startrail.calculate_fwhms()
            #startrail.print_fwhms_results()
            #startrail.print_fwhms_graph()

            median_trail_fwhm = startrail.fwhm_median
            self.trails_fwhm = np.append(self.trails_fwhm, median_trail_fwhm)
            self.trails_stddev = np.append(self.trails_stddev, startrail.fwhm_stddev)

    def mean_fwhm(self):
        """Give the mean FWHM of the trails"""
        return np.mean(self.trails_fwhm)

    def mean_stddev(self):
        """Give the standard deviation of the FWHM of the trails"""
        return np.mean(self.trails_stddev)


# ################
# Main program
# ################

def main(directory, file_prefix, file_suffix, number_of_files, location=''):
    """
    Measure the seeing of a sequence of images

    :param directory: the directory where the image files are located (with a slash or backslash at the end)
    :param file_prefix: beginning of the name of the files (example: 'seeing-')
    :param file_suffix: end of the name of the files (with extension) (example: '.fit')
                        It is assumed that the file format is FITS.
    :param number_of_files: the number of files (begins with 1 and ends with this number).
                            It is assumed that the number in the name of the files do not contains leading zero
                            (example: 'seeing-1.fit' to 'seeing-152.fit').
    :param location: the location where the images have been taken. It is printed in the title of the plot image.
    :return: void. Create or append to a file named 'seeing_measurement.csv' in execution directory with the results.
                   Show a chart with the results.
    """
    sampling = 0.206
    target_fwhm_arcsec = 1.5; # to help the gaussian fitting

    # Create a file to write the results
    with open('seeing_measurement.csv', 'a') as results:
        # Header of CSV file
        results.write('Date and time UTC,MJD,Seeing in arcsec,Std dev\n');

        measurements = np.array([])
        datetimes = np.array([])
        errorbar = np.array([])

        # For each file
        for i in range(1,number_of_files+1):

            # Open FITS file
            filename = file_prefix + str(i) + file_suffix
            path = directory + filename
            print "Read FITS file " + path
            hdulist = fits.open(path)

            # Get date and time of observation in FITS header
            datetime_string = hdulist[0].header['DATE-OBS']
            time = Time(datetime_string, format='isot', scale='utc')

            # Get length of X axis in FITS header
            length_x_axis = hdulist[0].header['NAXIS1']

            # Analyse data
            img_data = hdulist[0].data
            img = TrailsImage(img_data, sampling, target_fwhm_arcsec, length_x_axis)
            img.search_trails()
            img.calculate_fwhm()

            # Print results
            print "Date and time: %s UT" % datetime_string
            print "Number of trails: %i" % img.nb_trails()
            print "Mean FWHM of the trails: %f arcsec" % img.mean_fwhm()
            print "StdDev FWHM of the trails: %f" %img.mean_stddev()

            # Prepare plotting
            measurements = np.append(measurements, img.mean_fwhm())
            datetimes = np.append(datetimes, time.datetime)
            errorbar = np.append(errorbar, img.mean_stddev())

            # Close FITS file
            hdulist.close()

            # Write result in a file
            results.write(datetime_string + ',' + str(time.mjd) + ',' + str(img.mean_fwhm()) + ',' + str(img.mean_stddev()) + '\n');

            # Time of the first image of the sequence (used below)
            if(i == 1):
                start_time = time

    # Close results file
    results.closed

    # Plot results
    start_time.out_subfmt='date'

    plt.figure(1)
    plt.title('Seeing ' + location + ' ' + start_time.iso + ' (MJD ' + str(int(start_time.mjd)) + ')', fontsize=16)
    plt.errorbar(datetimes, measurements, fmt='ko', yerr=errorbar)
    plt.xlabel('Time (UT)')
    plt.gca().xaxis.set_major_formatter(dates.DateFormatter('%H:%M'))
    plt.xticks(rotation='vertical')
    plt.ylabel('Seeing (arcsec)')
    plt.show()


#main('/home/didier/seeing_images/2015-09-17/', 'zenith_sans_suivi-', '.fits', 53, 'St-Veran')
#main('/home/didier/seeing_images/2015-09-18/', 'zenith-', '.fits', 175, 'St-Veran')
#main('/home/didier/seeing_images/2015-09-19/', 'zenith-', '.fits', 91, 'St-Veran')
#main('/home/didier/seeing_images/2015-09-19/', 'zenith_refocus1-', '.fits', 153, 'St-Veran')

main('/home/didier/seeing_images/2015-09-19/', 'zenith_refocus1-', '.fits', 3, 'St-Veran')