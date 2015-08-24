# -*- coding: utf-8 -*-

"""
Functions related to signal processing (e.g. data filtering).

fft_freq_signal : function
    Fast-fourier transform
band_pass : function
    Returns a boolean array with value true between the frequency interval low and high. The return array has the same shape as frequency_bins.
"""

#############################################################################

from copy import copy
import numpy

def fft_freq_signal(x, Fs):
    """
    Fast-fourier transform.

    Parameters
    ----------
    x : array
        Signal to transform
    Fs : int
        Sampling frequency

    Returns
    -------
    frequency_bins: ndarray
        Array of length 'n'containing the sample frequencies (see numpy.fft.fftfreq).
    x_fft: complex ndarray
        The truncated or zero-padded input, transformed along the axis indicated by 'axis', or the last one if 'axis' is not specified (see numpy.fft.fftpack.fft).

    """
    frequency_bins=numpy.fft.fftfreq(x.shape[0], 1./Fs)
    x_fft=numpy.fft.fftpack.fft(x)
    return frequency_bins, x_fft


def band_pass(frequency_bins, low, high):
    """
    Return a boolean array with value true between the frequency interval low and high. The return array has the same shape as frequency_bins.

    Parameters
    ----------
    frequency_bin : array
        Frequency bins array produced by fft_freq_signal
    low : float
        Lower frequency
    high : float
        Higher frequency

    Returns
    -------
    truth_array: array
        Bool array which is true for frequency_bins which are between the frequency interval low and high.

    """
    truth_array=((frequency_bins>=low) & (frequency_bins<=high)) | ((frequency_bins<=-low) & (frequency_bins>=-high))
    return numpy.invert(truth_array)