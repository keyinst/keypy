# -*- coding: utf-8 -*-

###############################################################
# test3_mstates.py computes a series of modelmaps at a time ###
###############################################################


################################
# 1.) LOAD PACKAGES ###
################################

#external
import os
import os.path
import unittest

#integrated in the key.py library
from keypy.preprocessing.helper_functions import *
from keypy.preprocessing.data_loading import *
from keypy.preprocessing.avg_referencing import *
from keypy.preprocessing.filtering import *
from keypy.microstates.modelmaps import * 
from keypy.microstates.meanmods import * 

####   Classes     ####
from keypy.microstates.configuration import *


class Test_test3_mstates(unittest.TestCase):
    def test3_mstates(self):
        ################################
        # 2.) Specify data folder info ###
        ################################

        library_path = os.path.dirname(os.path.abspath(__file__))

        #contains data loaded into hdf5 file with the preprocessing script
        inputfolder = os.path.join(library_path, "..","data","test23")
        #will be created to contain output hdf5 files from microstate processing
        outputfolder = os.path.join(library_path, "..","data","test3_output")

        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)

        #name of hdf5 that contains data
        inputhdf5 = os.path.join( inputfolder, 'all_recordings.hdf')


        #####################################
        # 3.) Create EEG info object      ###
        #####################################

        library_path = os.path.dirname(os.path.abspath(__file__))

        #contains StudyInfo script
        script_inputfolder = library_path

        #exclude before committ
        #execfile(os.path.join( script_inputfolder, 'study_info_test1.py'))
        #include before committ
        ###########################
        ### Create EEG info object ###
        ###########################

        #enter number of channels
        nch=61
        #enter number of time frames for each segment
        tf = 512
        #enter sampling rate
        sf = 256
        # enter your channel list in the same order as in your files
        chlist=['FP1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3','FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1','P3','P5','P7','PO7','PO3','O1','Oz','POz','Pz','CPz','FPz','FP2','AF8','AF4','AFz','Fz','F2','F4','F6','F8','FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8','TP8','CP6','CP4','CP2','P2','P4','P6','P8','PO8','PO4','O2']

        #create object of EEG info information
        eeg_info_study_obj = EegInfo(nch, tf, sf, chlist)

        ################################
        ### Specify data folder info ###
        ################################

        inputfolder = os.path.join(library_path, "..\\data\\test23")
        outputfolder = inputfolder

        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)

        outputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')
        loaddata_output = 'rawdata'

        ########################
        ## Specify data info ###
        ########################

        ##Where in the filename is the following information (at which index of the string)? inclusive
        ##Where in the folder structure is the following information? 0 outtermost folder, 1 second folder, 2 third folder
        ##Specify for each element (group, partcipant, condition, run) whether the information is in the filename or folder structure.

        #group
        group_indices_range = [0,1]
        group_folder_level = 0
        has_group = 'folder' ###can be 'folder', 'filename', or 'none'

        #participant
        participant_indices_range = [2,3]
        participant_folder_level = 0
        has_participant = 'filename'

        #condition
        condition_indices_range = [5,6]
        condition_folder_level = 1
        has_condition = 'filename'

        #run
        run_indices_range = [6,6]
        run_folder_level = 0
        has_run = 'filename'

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
        enable_detrending = True

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

        hdf5_filename = 'all_recordings.hdf'
        inputhdf5 = os.path.join( outputfolder, hdf5_filename)

        boxkeyfilter(inputhdf5, eeg_info_study_obj, filter_input, filter_settings, enable_detrending = False)

        ################################################
        # 4.)  Get study info object from hdf5 file  ###
        ################################################

        study_info_obj=create_study_info_obj_from_data(inputhdf5)


        #######################################
        # 5.) Create Configuration Object   ###
        #######################################

        confobj = MstConfiguration(
                                subtract_column_mean_at_start = False,
                                debug = False,
                                use_gfp_peaks = True,
                                force_avgref = True,
                                set_gfp_all_1 = False,
                                use_smoothing = False,
                                gfp_type_smoothing='hamming',
                                smoothing_window=100,
                                use_fancy_peaks = False,
                                method_GFPpeak = 'GFPL1',
                                original_nr_of_maps = 4,
                                seed_number = 10,
                                max_number_of_iterations = 10,
                                ERP = False,
                                correspondance_cutoff = 0.00)


        #################
        # 6.) #Run Microstates (computes 1 microstate for each dataset in inputfile)
        #################

        ######
        ###Define input processing stage and output hdf5 file group
        ######

        modmaps_input = 'mstate1'
        modmaps_output = 'microstate'


        #include before commit
        run_modmaps(confobj, eeg_info_study_obj, inputhdf5, modmaps_input, modmaps_output)
        #--------------------------------------------------------------------------------------------------------------------------------------------


        #################
        # 7.) #Run Modelmaps (run_modelmaps_for_modelmap_types computes modelmaps for all types selected)
        #################

        ############
        # Compute modelmaps: Series of input and output hdfs
        ############

        #Series_1
        #means across runs for each group pt cond
        #means across conds for each group pt
        #means across pts for each group
        #means across groups

        #Series_2
        #means across pts for each group cond run
        #means across runs for each group cond

        #Series_3
        #means across runs for each group pt cond
        #means across pts for each group cond
        #means across groups for each cond
        #means across conds

        #Series_4
        #means across runs for each group pt cond
        #means across conds for each group pt
        #means across pts for each group
        #means across groups

        #Series_5
        #means across runs for each group pt cond
        #means across conds for each group pt
        #means across groups for each pt
        #means across groups

        series_versions = ['Series_1', 'Series_2', 'Series_3', 'Series_4', 'Series_5']

        for series in series_versions:
            first_input = 'microstate'

            #create folder with name of series as outputfolder
            outputfolder = os.path.join(library_path, "..\\data\\test3_output\\{0}" .format(series))
            if not os.path.exists(outputfolder):
                os.makedirs(outputfolder)

            run_meanmods_series(series, inputfolder, hdf5_filename, outputfolder, first_input, confobj)



