# -*- coding: utf-8 -*-

from __future__ import print_function

from contextlib import closing

import h5py

from keypy.preprocessing.file_info_classes import *

import numpy as np

#############################################################################
##Get study info object from hdf5 file
##dictionary of groups, pts, conds, runs
#############################################################################


def create_study_info_obj_from_data(inputhdf5):
    """
    Creates a study information object that includes the information on the occurring groups, participants, conditions, and runs

    Parameters
    ----------
    inputhdf5 : str
        Input path of the hdf5 file that contains the data to be processed, e.g. 'C:\\Users\\...\\input\\rawdata.hdf'.

    Returns
    -------
    study_info_obj : object of type StudyInfo
        Contains the following attributes: group_dict (dictionary that contains the whole data set structure: group, participant, cond, run) and group_list (list of strings for each group in the dataset).
    """

    with closing( h5py.File(inputhdf5, 'r') ) as f:
        #dictionary of groups
        group_dict = {}
        for group_name in f['/'].keys():
            group_dict[group_name] = {}
            for pt_name in f['/%s' % (group_name)].keys():
                group_dict[group_name][pt_name] = {}
                for cond_name in f['/%s/%s' % (group_name, pt_name)].keys():
                    group_dict[group_name][pt_name][cond_name] = {}
                    for run_name in f['/%s/%s/%s' % (group_name, pt_name, cond_name)].keys():
                        group_dict[group_name][pt_name][cond_name][run_name] = {}

        #create object of EEG study information
        study_info_obj = StudyInfo(group_dict)

    return study_info_obj


#############################################################################
##Delete channels (to do: integrate into testing framework)
##e.g. ch_to_delete = ['Iz','P9','P10']
#############################################################################


def del_channels(inputhdf5, ch_to_delete, del_channels_input, del_channels_output, eeg_info_study_obj):
    """
    Deletes the specified channels from the inputted data and saves them to a new hdf5 dataset.

    Parameters
    ----------
    inputhdf5 : str
        Input path of the hdf5 file that contains the data to be processed, e.g. 'C:\\Users\\Patricia\\libs\\keypy\\example\\data\\input\\rawdata.hdf'.
    ch_to_delete : list
        List of channel names which are to be deleted e.g. ['Iz','P9','P10'].
    del_channels_input : str
        Name of the dataset in the hdf5 file that the deletion is computed on, e.g. 'rawdata'.
    del_channels_output : str
        Name of the dataset in the hdf5 file that is created for each file where channels were deleted, e.g. 'rawdata_61nch'.
    eeg_info_study_obj : object of type EegInfo
        Contains the following attributes: nch (number of channels), tf (number of time frames per epoch), sf (sampling frequency), chlist (channel list).


    Returns
    -------
    chlist_process : list
        List of the remaining channels after deletion.
    """

    nch = eeg_info_study_obj.nch
    chlist = eeg_info_study_obj.chlist

    with closing( h5py.File(inputhdf5) ) as f:
        for groupi in f['/'].keys():
            for pti in f['/%s' % (groupi)].keys():
                for cond in f['/%s/%s' % (groupi, pti)].keys():
                    for run in f['/%s/%s/%s' % (groupi, pti, cond)].keys():
                        try:
                            timeframe_channel_dset = f['/{0}/{1}/{2}/{3}/{4}' .format(groupi, pti, cond, run, del_channels_input)]
                        except:
                            print('not found',  ['/{0}/{1}/{2}/{3}/{4}' .format(groupi, pti, cond, run, del_channels_input)])
                            continue
                    
                        path = '/{0}/{1}/{2}/{3}/{4}' .format(groupi, pti, cond, run, del_channels_input)
                
                        print('deleting channels for ', pti, cond, run)
                        timeframe_channel_dset = f[path]    
               
                        timeframe_channel=timeframe_channel_dset.value

                        chlist_process = list(chlist)

                        for channi in ch_to_delete:
                            index=chlist_process.index(channi)
                            timeframe_channel[:,index]

                            timeframe_channel=np.delete(timeframe_channel, index, 1)
                            chlist_process.remove(channi)


                        
                        if not del_channels_output in f['/{0}/{1}/{2}/{3}' .format(groupi, pti, cond, run)].keys():
                            i_channels_output = f['/{0}/{1}/{2}/{3}' .format(groupi, pti, cond, run)].create_dataset(del_channels_output, data = timeframe_channel)
                        else:
                            i_channels_output = f['/{0}/{1}/{2}/{3}' .format(groupi, pti, cond, run, del_channels_output)]


                        #i_channels_output[:] = timeframe_channel                        

    return chlist_process


#############################################################################
##Delete participants
##e.g. vps_to_delete = ['vp08']
#############################################################################


def del_participants(inputhdf5, pts_to_delete):
    """
    Creates a study information object that includes the information on the occurring groups, participants, conditions, and runs

    Parameters
    ----------
    inputhdf5 : str
        Input path of the hdf5 file that contains the data to be processed, e.g. 'C:\\Users\\Patricia\\libs\\keypy\\example\\data\\input\\rawdata.hdf'.
    pts_to_delete : list
        List of participant names which are to be deleted e.g. ['pt01','pt08','pt17'].
    """

    with closing( h5py.File(inputhdf5) ) as f:
        for groupi in f['/'].keys():
            for pti in f['/%s' % (groupi)].keys():
                if pti in pts_to_delete:
                    print('Deleting ....', group, pti)
                    del f['/%s/%s' % (groupi, pti)]
