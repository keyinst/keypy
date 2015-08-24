# -*- coding: utf-8 -*-

##################################
#######  Import Packages  ########
##################################

import numpy
import scipy
import scipy.signal

from numpy.linalg.linalg import eig
from numpy.core.fromnumeric import size, argsort

def princomp_B(A,numpc=0):
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
        coeff = coeff[:,range(numpc)] # cutting some PCs
    #score = dot(coeff.T,M) # projection of the data in the new space
    return coeff


def compute_gfp(X, method = 'GFPL2'):
    """
    Compute the Global Field Power
 
    Parameters
    ----------
    X (ndarray):
        Array containing values for all time frames and channels.
        Dimension: number of time frames x number of channels
    method ({'GFPL1', 'GFPL2'}):
        `GFPL1` : use L1-Norm to compute GFP peaks
        `GFPL2` : use L2-Norm to compute GFP peaks
        
    Returns
    ----------
        ret : ndarray
			GFP curve
    """
    ntf = X.shape[0]
    ret = numpy.zeros( (ntf, ) )
    
    if method == 'GFPL2':
        # gfp = sqrt(sum(abs(x - x.mean())**2 / len(x) ))
        for i in xrange(ntf):
            x = X[i,:]
            gfp = numpy.sqrt(numpy.sum(numpy.abs(x - x.mean())**2 / len(x) ))
            ret[i] = gfp
    elif method == 'GFPL1':
        # gfp = sum(abs(x - x.mean()) / len(x) ))
        for i in xrange(ntf):
            x = X[i,:]
            gfp = numpy.sum(numpy.abs(x - x.mean())) / len(x) 
            ret[i] = gfp
    #print 'compute gfp', ret        
    return ret


def l2mm_norm(X,Y):
    """
    Compute L2 Norm.

    Parameters
    ----------
	X: double
	Y: double

    Returns
    ----------
	l2 norm: double

    """
    return numpy.linalg.norm(X - Y, ord = None)


def l1mm_norm(X,Y):
    """
    Compute L1 Norm.

    Parameters
    ----------
	X: double
	Y: double

    Returns
    ----------
	l1 norm: double
    """
    return numpy.linalg.norm(X - Y, ord = 1)


#Smoothing: Convolving with window function
def gfp_smoothing(gfp, method, window):
    """
    Smoothe the Global Field Power Curve
 
    Args:
        gfp (ndarray):
            Array containing values for all time frames.
            Dimension: number of time frames x 1
        method ({'hamming', 'hanning'}):
            `hamming` : use hamming window to smooth
            `hanning` : use hanning window to smooth
        window (int):
            about 100
        
    Returns:
        smooth : array
			smoothed GFP curve
    """
    if method == 'hamming':
        smooth = scipy.signal.convolve(gfp, scipy.signal.hamming(window) )
    elif method == 'hanning':
        smooth = scipy.signal.convolve(gfp, scipy.signal.hanning(window) )
    return smooth
    

#Find GFP Peaks
def gfp_peaks_indices(gfp):
    """
    Finds relative maxima across gfp curve
 
    Args:
        gfp (ndarray):
            Array containing values for all time frames.
            Dimension: number of time frames x 1
            (smoothed or not)

    Returns:
        gfp_peak_indices : list
			list of indices of GFP curve that are peaks
    """
    return scipy.signal.argrelmax( gfp )[0]