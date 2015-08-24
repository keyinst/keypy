# -*- coding: utf-8 -*-

###########################
###    Load packages    ###
###########################

import os
import os.path
import unittest

from keypy.preprocessing.file_info_classes import *
from keypy.preprocessing.data_loading import *
from keypy.preprocessing.avg_referencing import *
from keypy.preprocessing.filtering import *
from keypy.preprocessing.helper_functions import *

class Test_a_mstates_preprocessing_test10(unittest.TestCase):
    def test_a_mstates_preprocessing_test10(self):
        ###########################
        ### Run Study Info File ###
        ###########################

        #enter number of channels
        nch=64
        #enter number of time frames for each segment
        tf = 512
        #enter sampling rate
        sf = 256
        # enter your channel list in the same order as in your files
        chlist = ['Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3','FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1','P3','P5','P7','P9','PO7','PO3','O1','Iz','Oz','POz','Pz','CPz','Fpz','Fp2','AF8','AF4','AFz','Fz','F2','F4','F6','F8','FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8','TP8','CP6','CP4','CP2','P2','P4','P6','P8','P10','PO8','PO4','O2']

        #create object of EEG info information
        eeg_info_study_obj = EegInfo(nch, tf, sf, chlist)

        ################################
        ### Specify data folder info ###
        ################################

        library_path = os.path.dirname(os.path.abspath(__file__))

        inputfolder = os.path.join(library_path, "..","data","test10")
        outputfolder = os.path.join(library_path, "..","data","test10_output")


        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)

        outputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')
        loaddata_output = 'rawdata'

        if os.path.isfile(outputhdf5):
            os.remove(outputhdf5)
        self.assertFalse(os.path.isfile(outputhdf5))

        ########################
        ## Specify data info ###
        ########################

        ##Where in the filename is the following information (at which index of the string)? inclusive
        ##Where in the folder structure is the following information? 0 outtermost folder, 1 second folder, 2 third folder
        ##Specify for each element (group, partcipant, condition, run) whether the information is in the filename or folder structure.

        #group
        group_indices_range = [0,1]
        group_folder_level = 0
        has_group = 'none' ###can be 'folder', 'filename', or 'none'

        #participant
        participant_indices_range = [3,5]
        participant_folder_level = 0
        has_participant = 'none'

        #condition
        condition_indices_range = [0,0]
        condition_folder_level = 1
        has_condition = 'none'

        #run
        run_indices_range = [0,0]
        run_folder_level = 0
        has_run = 'none'

        file_ending = 'txt'

        ##Generate objects based on info above

        file_name_obj = FileNameDescription(group_indices_range, participant_indices_range, condition_indices_range, run_indices_range, file_ending)
        folder_structure_obj = FolderStructureDescription(group_folder_level, participant_folder_level, condition_folder_level, run_folder_level)
        filename_folder_obj = FilenameFolderDescription(has_group, has_participant, has_condition, has_run)

        ########################
        ###Load Data to HDF 5###
        ########################

        loaddata(inputfolder, outputhdf5, loaddata_output, file_name_obj, folder_structure_obj, filename_folder_obj, eeg_info_study_obj)

        ########################
        ###Compute Avg Ref   ###
        ########################

        inputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')

        average_input = 'rawdata'
        average_output = 'avg_ref'

        averageref(inputhdf5, average_input, average_output )


        #################################
        ###  Filter for microstates   ###
        #################################

        ######
        ###Use detrending before filtering?
        ######
        enable_detrending = False

        ######
        ###Choose Filter Settings
        ######
        filter_settings = {
            "mstate1": {
                "low": 2.0,
                "high": 20
            },
            "mstate2": {
                "low": 1.5,
                "high": 30
            },
            "mstate3": {
                "low": 1.5,
                "high": 45
            }
        }


        ######
        ###HDF5 File Preparation
        ######
        #define processing stage of input and name of output
        filter_input = 'avg_ref'
        filter_output = filter_settings

        ######
        ###Inputhdf (Outputhdf from before)
        ######

        inputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')

        boxkeyfilter(inputhdf5, eeg_info_study_obj, filter_input, filter_settings, enable_detrending = False)


        for test_dataset in ['mstate1', 'mstate2', 'mstate3']:

            dataset = None
            with closing (h5py.File(outputhdf5, 'r')) as f:
                g1 = f['/group_All_PTs']
                g2 = g1['pt_PT']
                g3 = g2['cond_Cond']
                g4 = g3['run_1']
                dataset = g4[test_dataset][:]

            correct_solution = np.loadtxt(os.path.join(outputfolder, 'correct_solution_{0}.asc' .format(test_dataset)))

            self.assertEqual(len(dataset), len(correct_solution))
            for i in range(0, len(dataset)):
                self.assertEqual(len(dataset[i]), len(correct_solution[i]))
                for j in range(0, len(dataset[i])):
                    self.assertAlmostEqual(dataset[i][j], correct_solution[i][j], places=1)
