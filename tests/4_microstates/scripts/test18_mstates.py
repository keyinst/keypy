# -*- coding: utf-8 -*-

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
from keypy.preprocessing.helper_functions import *

####   Classes     ####
from keypy.microstates.configuration import *
from keypy.preprocessing.file_info_classes import *

from keypy.microstates.modelmaps import *

class Test_test18_mstates(unittest.TestCase):
    def test18_mstates(self):
        ################################
        # 2.) Specify data folder info ###
        ################################

        library_path = os.path.dirname(os.path.abspath(__file__))

        #contains data loaded into hdf5 file with the preprocessing script
        inputfolder = os.path.join(library_path, "..","data","test18")
        #will be created to contain output hdf5 files from microstate processing
        outputfolder = inputfolder


        #name of hdf5 that contains data
        inputhdf5 = os.path.join( inputfolder, 'all_recordings.hdf')


        #####################################
        # 3.) Create EEG info object      ###
        #####################################

        #enter number of channels
        nch=19
        #enter number of time frames for each segment
        tf = 256
        #enter sampling rate
        sf = 256
        # enter your channel list in the same order as in your files
        chlist = ['FP1','FP2','F7','F3','Fz','F4','F8','T7','C3','Cz','C4','T8','P7','P3','Pz','P4','P8','O1','O2']
       
        #create object of EEG info information
        eeg_info_study_obj = EegInfo(nch, tf, sf, chlist)

        ################################
        ### Specify data folder info ###
        ################################

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
        inputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')


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
                                use_gfp_peaks = False,
                                force_avgref = True,
                                set_gfp_all_1 = False,
                                use_smoothing = False,
                                gfp_type_smoothing='hamming',
                                smoothing_window=100,
                                use_fancy_peaks = False,
                                method_GFPpeak = 'GFPL2',
                                original_nr_of_maps = 4,
                                seed_number = 100,
                                max_number_of_iterations = 100,
                                ERP = False,
                                correspondance_cutoff = 0.00,
                                fixed_seed = 1)


        #################
        # 6.) #Run Microstates (computes 1 microstate for each dataset in inputfile)
        #################

        ######
        ###Define input processing stage and output hdf5 file group
        ######

        modmaps_input = 'rawdata'
        modmaps_output = 'microstate'

        run_modmaps(confobj, eeg_info_study_obj, inputhdf5, modmaps_input, modmaps_output)
        #--------------------------------------------------------------------------------------------------------------------------------------------

        dataset = None
        with closing (h5py.File(outputhdf5, 'r')) as f:
            g1 = f['/group_All_PTs']
            g2 = g1['pt_PT']
            g3 = g2['cond_Cond']
            g4 = g3['run_1']
            dataset = g4['microstate']
            modelmap_exp_var = float(dataset.attrs['explained variance of all gfp peaks'])
            modelmap_exp_var_all = float(dataset.attrs['explained variance of all eeg timeframes'])

        correct_solution = np.loadtxt(os.path.join(inputfolder, 'correct_output_test18.asc'))

        self.assertAlmostEqual(modelmap_exp_var, correct_solution[0])
        self.assertAlmostEqual(modelmap_exp_var_all, correct_solution[1])

        #for almost equal compare column by column with self.assertAlmostEqual