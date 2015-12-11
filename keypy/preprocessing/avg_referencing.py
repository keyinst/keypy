# -*- coding: utf-8 -*-

#############################################################################
##COMPUTE AVERAGE REFERENCE (subtract mean across all channels for each channel)
#############################################################################

from __future__ import print_function

from contextlib import closing

import h5py
import numpy as np


def averageref(inputhdf5, average_input, average_output ):
    """
    Compute the average reference for each dataset within the inputhdf5 of level average_input.

    Parameters
    ----------
    inputhdf5 : str
        Input path of the hdf5 file that contains the data to be processed, e.g. 'C:\\Users\\...\\rawdata.hdf'
    average_input : str
        Name of the dataset in the hdf5 file that the average referencing is computed on, e.g. 'rawdata'.
    average_output : str
        Name of the dataset in the hdf5 file that is created for each average-referenced file, e.g. 'avg_ref'.
    """

    print('Compute average reference ....')
    with closing( h5py.File(inputhdf5) ) as f:
        for groupi in f['/'].keys():
            for pti in f['/%s' % (groupi)].keys():
                for cond in f['/%s/%s' % (groupi, pti)].keys():
                    for run in f['/%s/%s/%s' % (groupi, pti, cond)].keys():
                        try:
                            timeframe_channel_dset = f['/{0}/{1}/{2}/{3}/{4}' .format(groupi, pti, cond, run, average_input)]
                        except:
                            print('not found',  ['/{0}/{1}/{2}/{3}/{4}' .format(groupi, pti, cond, run, average_input)])
                            continue
                    
                        path = '/{0}/{1}/{2}/{3}/{4}' .format(groupi, pti, cond, run, average_input)
                
                        print('avg_referencing', pti, cond, run)
                        timeframe_channel_dset = f[path]
                
                        if not average_output in f['/{0}/{1}/{2}/{3}' .format(groupi, pti, cond, run)].keys():
                            i_averageref = f['/{0}/{1}/{2}/{3}' .format(groupi, pti, cond, run)].create_dataset(average_output, data = timeframe_channel_dset.value)
                        else:
                            i_averageref = f['/{0}/{1}/{2}/{3}' .format(groupi, pti, cond, run, average_output)]
                
                        timeframe_channel=timeframe_channel_dset.value
                        computed_mean = timeframe_channel.mean(axis=1)
                        i_averageref[:] = timeframe_channel - computed_mean[:, np.newaxis] 
