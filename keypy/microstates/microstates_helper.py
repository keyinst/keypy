# -*- coding: utf-8 -*-

##################################
#######  Import Packages  ########
##################################

import numpy as np
import scipy
import scipy.signal

import os
import os.path

from numpy.linalg.linalg import eig
from numpy.core.fromnumeric import size, argsort

from numpy import linalg as LA

####------------------------------------------####


##################################################
#######  Principal Component Computation  ########
##################################################


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
    a=np.dot(M,M.T) #get covariance matrix by  matrix multiplication of data with transposed data
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
    return coeff

####------------------------------------------####


#######################################################
#######  Global Field Power (GFP) Computation  ########
#######################################################

def compute_gfp(X, method = 'GFPL1'):
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
    ret = np.zeros( (ntf, ) )
    
    if method == 'GFPL2':
        # gfp = sqrt(sum((x - x.mean())**2 / len(x) ))
        for i in range(ntf):
            x = X[i,:]
            gfp = np.sqrt(np.sum((x - x.mean())**2 / len(x) ))
            ret[i] = gfp
    elif method == 'GFPL1':
        # gfp = sum(abs(x - x.mean()) / len(x) ))
        for i in range(ntf):
            x = X[i,:]
            gfp = np.sum(np.abs(x - x.mean())) / len(x) 
            ret[i] = gfp
    #print 'compute gfp', ret        
    return ret

####------------------------------------------####


############################
#######  L2mm Norm  ########
############################

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
    return np.linalg.norm(X - Y, ord = None)

####-------------------------------------####

############################
#######  L1mm Norm  ########
############################

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
    return np.linalg.norm(X - Y, ord = 1)

####-------------------------------------####

######################################
#######  GFP Curve Smoothing  ########
######################################

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
    

####-------------------------------------####

#################################
#######  Find GFP Peaks  ########
#################################

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

####-------------------------------------####

############################################
#######  Dissimilarity Computation  ########
############################################

def dissim(map, tf, method_GFPpeak):
    """
    Compute the Dissimilarity between two maps of identical dimension.
 
    Parameters
    ----------
    map :
        Array containing values for one modelmap (all channels).
        Dimension: 1 x number of channels

    tf :
        Array containing values for one time frame (all channels).
        Dimension: 1 x number of channels

    method_GFPpeak : {'GFPL1', 'GFPL2'}
        `GFPL1` : use L1-Norm to compute GFP peaks
        `GFPL2` : use L2-Norm to compute GFP peaks
        
    Returns
    ----------
        ret : diff
			Dissimilarity measure (ignores polarity)
    """

    #compute differences A-B and -A-B
    diff1=map-tf
    diff2=-map-tf

    #convert differences into numpy arrays of 1 dimension
    diff1=np.reshape(diff1, (-1, 1))
    diff2=np.reshape(diff2, (-1, 1))

    #computes global field power across all channels
    diff1_gfp=compute_gfp(diff1.T, method_GFPpeak)
    diff2_gfp=compute_gfp(diff2.T, method_GFPpeak)

    #determines the smaller dissimilarity value
    diff=min(diff1_gfp,diff2_gfp)

    return diff
####-------------------------------------####

####################################
#######  Compute GFP Peaks  ########
####################################

def compute_gfp_peaks(gfp_curve, use_gfp_peaks, use_smoothing, gfp_type_smoothing, smoothing_window, use_fancy_peaks):
    """
    Computes GFP peaks from global field power curve.

    Parameters
    ----------
    gfp_curve : 1D array
        Global field power for each time frame.
    use_gfp_peaks : bool
        Option whether whole GFP peaks are used or not.
    use_smoothing : bool
        Option whether smoothing is to be applied to the GFP curve before peak computation or not.
    gfp_type_smoothing : {'hamming', 'hanning'}
        `hamming` : use hamming window to smooth
        `hanning` : use hanning window to smooth
    smoothing_window : int
		window for smoothing, e.g. 100.
    use_fancy_peaks : bool
        Whether a particular smoothing algorithm from scipy.signal.find_peaks_cwt is applied before peak computation or not.
        Reference: Bioinformatics (2006) 22 (17): 2059-2065. doi: 10.1093/bioinformatics/btl355 http://bioinformatics.oxfordjournals.org/content/22/17/2059.long)

    Returns
    -------
    gfp_peak_indices : list
        List of indices of the EEG that qualify as global field power peaks.
    gfp_curve : 1D array
        GFP curve after smoothing (if smoothing was applied).
    """
    if use_gfp_peaks:
        #we need to do this for each epoch seperately, as not to introduce artificial peaks at the transitions between epochs
        #for faster computing we delete the indices that were between epochs afterwards
        #this is done in microstates.py after retrieval of the indices
        if use_smoothing:
            gfp_curve=gfp_smoothing(gfp_curve, gfp_type_smoothing, smoothing_window)
        if use_fancy_peaks:
            peakind = scipy.signal.find_peaks_cwt(gfp_curve, np.arange(1,10))
            gfp_peak_indices=np.asarray(peakind) #we would expect a peak at about each 50 ms
            gfp_curve = gfp_curve
        else:
            gfp_peak_indices=gfp_peaks_indices(gfp_curve) #we would expect a peak at about each 50 ms
            gfp_curve = gfp_curve
    else:
        gfp_peak_indices=np.array(list(range(len(gfp_curve))))   #when we take all maps, we still call the array gfp_peak_indices
        gfp_curve = gfp_curve
        print('all maps used')

    return gfp_peak_indices, gfp_curve

####-------------------------------------####

####################################
#######  set_gfp_all_1  ########
####################################

def set_gfp_all_1(eeg, gfp_curve):        
    """
    Normalizes EEG to set GFP to 1 for each time frame.

    Parameters
    ----------
    eeg : array
        Shape ntf*nch, conatains the EEG data the average referencing is to be computed on.
    gfp_curve : 1D array
        Global field power for each time frame.

    Returns
    -------
    eeg: array
        EEG with GFP set to 1.
    """
    for i in range(eeg.shape[0]):
        eeg[i,:] = eeg[i,:]/gfp_curve[i]
    return eeg

####-------------------------------------####

############################################################
#######  reduce channels for sortmaps / parameters  ########
############################################################

def reduce_channels(eeg_own, eeg_ext_path, own_chlist, external_chlist_path):
    #read info from external files
    ext_chlist = np.genfromtxt(external_chlist_path,dtype='str')
    eeg_ext=np.loadtxt(eeg_ext_path)


    ##Find overlap between chlists
    common_list=set(own_chlist).intersection(ext_chlist)

    #if there is no overlap, raise Error
    if not common_list:
        raise AssertionError('Your channel list has no overlap with the channel list which is to be sorted / parameters computed by. Processing failed.')   

    else:
        'Common List between your channels and the external list:', common_list, ' Watch out for case sensitivity!'

    ##Get indices of common_list order for own_chlist
    index_of_common_list_in_own_chlist = []
    for ele in common_list:
        index_of_common_list_in_own_chlist.append(own_chlist.index(ele))

    ##Get indices of common_list order for ext_chlist
    index_of_common_list_in_ext_chlist = []
    for ele in common_list:
        index_of_common_list_in_ext_chlist.append(ext_chlist.tolist().index(ele))

    ##Restructure own and ext EEG to common_list order
    eeg_own_new=eeg_own[:,index_of_common_list_in_own_chlist]
    eeg_ext_new=eeg_ext[:,index_of_common_list_in_ext_chlist]

    #at the time this saving procedure is done repeatedly for each EEG to be sorted, ideally this should only happen once
    np.savetxt(os.path.join(os.path.dirname(eeg_ext_path),"{0}_reduced.asc".format(os.path.splitext(os.path.basename(eeg_ext_path))[0])), eeg_ext_new)

    return eeg_own_new, eeg_ext_new
####-------------------------------------####

############################################################
#######  normalize maps (to GFP=1 or vector_norm=1)  ########
############################################################

def normalize_maps(maps, modelmaps_normalization_type):
    if modelmaps_normalization_type == 'gfp1':
        #set GFP of maps to 1   
        for ri in range( maps.shape[0] ):    
            # compute normalization_factor for randmap
            normalization_factor = compute_gfp(maps[ri,:].reshape( (1,maps.shape[1]) ))[0]                               
            # set randmap to GFP=1
            maps[ri,:] = maps[ri]/normalization_factor

        #assert that maps all have GFP = 1
        if not compute_gfp(maps).all() == 1:
            raise AssertionError('Problem with GFP normalization of model_maps_foundation or microstate maps / input models.')

    elif modelmaps_normalization_type == 'vector_norm_1':
        #set vector_length of maps to 1   
        for ri in range( maps.shape[0] ):    
        # compute normalization_factor for randmap
            normalization_factor = LA.norm(maps[ri,:], axis=0)                               
            # set maps to vector_length=1
            maps[ri,:] = maps[ri]/normalization_factor

        #assert that maps all have vector_length = 1
        if not (LA.norm(maps[ri,:], axis=0)).all() == 1:
            raise AssertionError('Problem with vector_length normalization of model_maps_foundation or microstate maps / input models.')

    return maps
####-------------------------------------####