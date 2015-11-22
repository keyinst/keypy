# -*- coding: utf-8 -*-

#############################################################################
##FILTER SIGNALS (with or without detrending)
##filter_settings{'mstate1': {'high': 20, 'low': 2.0}, 'mstate2': {'high': 30, 'low': 1.5}}
#############################################################################

from __future__ import print_function

from contextlib import closing

import h5py
import numpy
import scipy
import scipy.signal

from keypy.signal import band_pass


def boxkeyfilter(inputhdf5, eeg_info_study_obj, filter_input, filter_settings, enable_detrending = False):
    """
    Box filter for each dataset within the inputhdf5 of level average_input.

    Parameters
    ----------
    inputhdf5 : str
        Input path of the hdf5 file that contains the data to be processed, e.g. 'C:\\Users\\Patricia\\libs\\keypy\\example\\data\\input\\rawdata.hdf'
    eeg_info_study_obj : object of type EegInfo
        Contains the following attributes: nch (number of channels), tf (number of time frames per epoch), sf (sampling frequency), chlist (channel list)
    filter_input : str
        Name of the dataset in the hdf5 file that the filtering is computed on, e.g. 'avgref'.
    filter_settings : dict
        Dictionary with keys of the name of the output hdf5 dataset file name and two values for the lower and upper boarder of the filter, e.g. filter_settings = {"mstate1": {"low": 2.0, "high": 20}}.
    enable_detrending : bool
        Use detrending or not before filtering (False by default). 
    """

    print('Filter signals ...')
    with closing( h5py.File(inputhdf5) ) as f:
        for groupi in f['/'].keys():
            for pti in f['/%s' % (groupi)].keys():
                for cond in f['/%s/%s' % (groupi, pti)].keys():
                    for run in f['/%s/%s/%s' % (groupi, pti, cond)].keys():
                        try:
                            timeframe_channel_dset = f['/{0}/{1}/{2}/{3}/{4}' .format(groupi, pti, cond, run, filter_input)]
                        except:
                            print('not found',  ['/{0}/{1}/{2}/{3}/{4}' .format(groupi, pti, cond, run, filter_input)])
                            continue
                    
                        path = f['/{0}/{1}/{2}/{3}' .format(groupi, pti, cond, run)]   
                        dset = path[filter_input].value                

                        for filter_state in filter_settings.keys():
                            x_all_channels = numpy.zeros( dset.shape, dtype=dset.dtype )
                            frequency_bins=numpy.fft.fftfreq(eeg_info_study_obj.tf, 1./eeg_info_study_obj.sf)
                
                            for ch in range( dset.shape[1] ):
                                x = dset[:,ch]
                        
                                if enable_detrending:
                                    x = scipy.signal.detrend(x)

                                # loop across 2second epochs
                                for i in range(len(x)//eeg_info_study_obj.tf):
                                    epoch = x[i*eeg_info_study_obj.tf:(i+1)*eeg_info_study_obj.tf]
                                    epoch_fft=numpy.fft.fftpack.fft(epoch)
                                    selector = band_pass( frequency_bins, filter_settings[filter_state]["low"], filter_settings[filter_state]["high"] )
                                    epoch_fft[selector] = 0
                                    epoch_reverse_filtered = numpy.fft.fftpack.ifft(epoch_fft)
                                    x_all_channels[i*eeg_info_study_obj.tf:(i+1)*eeg_info_study_obj.tf,ch] = numpy.real(epoch_reverse_filtered)
                                                   
                            #to do: add the filter_settings parameters as attributes to the HDF5 dataset
                            if not filter_state in path.keys():
                                path.create_dataset(filter_state, data=x_all_channels )
                            else:
                                path[filter_state][:] = x_all_channels
