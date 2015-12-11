# -*- coding: utf-8 -*-

"""
The KEY EEG Analysis Toolbox

The KEY EEG Analysis Toolbox includes functionality for EEG preprocessing, microstate analysis, spectral analysis, and statistics.
"""

#################
# 0.) IMPORT PACKAGES
#################

import numpy
from numpy.linalg.linalg import eig
from numpy.core.fromnumeric import size, argsort

#################
# 0.) DEFINE FUNCTIONS
#################


def princomp_C(A,numpc=0):
    """
    Compute 1st Principal Component.
	Function modified from: http://glowingpython.blogspot.ch/2011/07/principal-component-analysis-with-numpy.html

    Parameters
    ----------
    A : object of type MstConfiguration
        Contains the following attributes: subtract_column_mean_at_start (bool), debug (bool), use_gfp_peaks (bool), force_avgref (bool), set_gfp_all_1 (bool), use_smoothing (bool), gfp_type_smoothing (string),
        smoothing_window (int), use_fancy_peaks (bool), method_GFPpeak (string), original_nr_of_maps (int), seed_number (int), max_number_of_iterations (int), ERP (bool), correspondance_cutoff (double):
    numpc : int
		number of principal components

    Returns
    -------
        coeff: array
            first principal component
		idx: 

		latent: 

    """
    M=A.T
    a=numpy.dot(M,M.T) #get covariance matrix by  matrix multiplication of data with transposed data
    [latent,coeff]=eig(a)
    p = size(coeff,axis=1)
    idx = argsort(latent) # sorting the eigenvalues
    idx = idx[::-1]       # in ascending order
    # sorting eigenvectors according to the sorted eigenvalues
    coeff = coeff[:,idx]
    #latent = latent[idx] # sorting eigenvalues
    if numpc < p or numpc >= 0:
        coeff = coeff[:,list(range(numpc))] # cutting some PCs
    #score = dot(coeff.T,M) # projection of the data in the new space
    return coeff, idx, latent


