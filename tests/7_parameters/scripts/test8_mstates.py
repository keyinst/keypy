# -*- coding: utf-8 -*-

###############################################################
# test8_mstates.py computes a series of modelmaps at a time ###
###############################################################


################################
# 1.) LOAD PACKAGES ###
################################

#external
import os
import os.path
import unittest
import h5py

#integrated in the key.py library
from keypy.preprocessing.helper_functions import *
from keypy.preprocessing.data_loading import *
from keypy.preprocessing.avg_referencing import *
from keypy.preprocessing.filtering import *
from keypy.microstates.modelmaps import * 
from keypy.microstates.meanmods import * 
from keypy.microstates.sortmaps import * 
from keypy.microstates.parameters import * 

####   Classes     ####
from keypy.microstates.configuration import *


class Test_test8_mstates(unittest.TestCase):
    def test8_mstates(self):
        ################################
        # 2.) Specify data folder info ###
        ################################

        library_path = os.path.dirname(os.path.abspath(__file__))

        #contains data loaded into hdf5 file with the preprocessing script
        inputfolder = os.path.join(library_path, "..","data","test8")
        #will be created to contain output hdf5 files from microstate processing
        outputfolder = os.path.join(library_path, "..","data","test8_output")

        if not os.path.exists(outputfolder):
            os.makedirs(outputfolder)

        #name of hdf5 that contains data
        inputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')


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

        inputfolder = os.path.join(library_path, "..\\data\\test8")
        outputfolder = os.path.join(library_path, "..\\data\\test8_output")

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
        group_indices_range = [0,2]
        group_folder_level = 0
        has_group = 'filename' ###can be 'folder', 'filename', or 'none'

        #participant
        participant_indices_range = [6,7]
        participant_folder_level = 0
        has_participant = 'filename'

        #condition
        condition_indices_range = [9,12]
        condition_folder_level = 1
        has_condition = 'filename'

        #run
        run_indices_range = [14,14]
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
                                method_GFPpeak = 'GFPL1',
                                original_nr_of_maps = 4,
                                seed_number = 5,
                                max_number_of_iterations = 5,
                                ERP = False,
                                correspondance_cutoff = 0.00)


        #################
        # 6.) #Run Microstates (computes 1 microstate for each dataset in inputfile)
        #################

        ######
        ###Define input processing stage and output hdf5 file group
        ######
        inputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')

        modmaps_input = 'rawdata'
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
        '''
        confobj = MstConfiguration(
                                seed_number = 5,
                                max_number_of_iterations = 10)


        series_versions = ['Series_1', 'Series_2', 'Series_3', 'Series_4', 'Series_5']

        outputfolder = os.path.join(library_path, "..\\data\\test8_output")
        inputfolder = outputfolder

        for series in series_versions:
            first_input = 'microstate'

            #create folder with name of series as outputfolder
            outputfolder_series = os.path.join(library_path, "..\\data\\test8_output\\{0}" .format(series))
            if not os.path.exists(outputfolder_series):
                os.makedirs(outputfolder_series)

            run_meanmods_series(series, inputfolder, outputfolder_series, first_input, confobj)
        '''
        #--------------------------------------------------------------------------------------------------------------------------------------------


        #################
        # 7.) #Run Sortmaps (run_sortmaps_for_sortmap_types computes sortmaps for all types selected) - Not needed for external sort
        #################

        #confobj = MstConfiguration()
        '''
        series_versions = ['Series_1', 'Series_2', 'Series_3', 'Series_4', 'Series_5']

        first_input = 'microstate'
        sortbyfolder = os.path.join(library_path, "..","data","sortby")
        outputfolder = os.path.join(library_path, "..\\data\\test8_output")

        for series in series_versions:
            
            run_sortmaps_series(series, inputfolder, sortbyfolder, outputfolder, first_input, confobj)      
        '''    
        #--------------------------------------------------------------------------------------------------------------------------------------------


        #################
        # 7.) #Run Parameters
        #################

        #############################
        #### Continuous EEG data ####
        #############################

        ##info needed to know which data the parameters are to be computed upon
        inputfolder = os.path.join(library_path, "..\\data\\test8_output")
        hdf5_filename = 'all_recordings.hdf'
        inputdataset = 'rawdata'

        ############################
        ####    Parameter by    ####
        ############################

        ##info needed to know which data the parameters are to be "sorted" upon (modelmaps)
        sortbyfolder = os.path.join(library_path, "..","data","sortby")
        parameter_type = 'external' #can be external, series, inputhdf
        ###
        #if you use an external asci file to sort your data by
        ###
        if parameter_type == 'external' : 
            parameter_by = 'external_norm'
            sortbyfile ='mean_models_milz_etal_2015.asc'
            external_chlist='mean_models_milz_etal_2015_chlist.asc'
            sortbydataset = None
            sortbyseries = None

        ###
        #if you use a previously computed hdf5 file (from Series)
        ###
        elif parameter_type == 'series': 
            sortbyseries = 'Series_3'
            sortbyfile = 'meanmods_across_runs_sorted.hdf'
            #sortbydataset = 'microstate_Series_1_sorted'
            sortbydataset = 'modelmap'
            external_chlist = False

            #specify the number of layers of your sortbyfile
            #parameter_by = '4Levels'
            parameter_by = '3Levels'
            #parameter_by = '2Levels'
            #parameter_by = '1Level'

        ###
        #if you use input hdf5 file sort
        ###
        elif parameter_type == 'inputhdf':
            parameter_by = 'own_hdf'
            sortbyfile = hdf5_filename
            sortbydataset = 'microstate_Series_1_sorted'
            sortbyseries = None
            external_chlist = False
    

        else:
            'Error, type not correctly specified for parameter computation.'

        data_provider=get_data_provider_for_parameter_by(parameter_by, inputfolder, hdf5_filename, inputdataset, sortbyfolder, sortbyfile, sortbydataset, sortbyseries, external_chlist)

        run_parameters(data_provider, confobj, eeg_info_study_obj)
        #--------------------------------------------------------------------------------------------------------------------------------------------
        
        dataset = None
        with closing (h5py.File(os.path.join(sortbyfolder, 'mean_models_milz_etal_2015', 'mstate_parameters.hdf5'), 'r')) as f:
            g1 = f['/group_HEA']
            g2 = g1['pt_01']
            g3 = g2['cond_Rest']
            g4 = g3['run_0']
            dataset1 = g4['Individual_States'][:]
            dataset2 = g4['State Match Mean p'][:]
            g5 = g4['ep_000']
            map_00 = g5['map_00']
            map_01 = g5['map_01']
            map_02 = g5['map_02']
            map_03 = g5['map_03']
            cov_dataset_0 = map_00['Coverage in percent'].value
            cov_dataset_1 = map_01['Coverage in percent'].value
            cov_dataset_2 = map_02['Coverage in percent'].value
            cov_dataset_3 = map_03['Coverage in percent'].value
            dur_dataset_0 = map_00['Mean duration in ms'].value
            dur_dataset_1 = map_01['Mean duration in ms'].value
            dur_dataset_2 = map_02['Mean duration in ms'].value
            dur_dataset_3 = map_03['Mean duration in ms'].value
            occ_dataset_0 = map_00['Occurrance per s'].value
            occ_dataset_1 = map_01['Occurrance per s'].value
            occ_dataset_2 = map_02['Occurrance per s'].value
            occ_dataset_3 = map_03['Occurrance per s'].value

            #load correct solutions
            individu_correct_solution = np.loadtxt(os.path.join(inputfolder, 'solutions', 'Individual_States_HEA01Rest00.asc'))
            state_match_mean_correct_solution = np.loadtxt(os.path.join(inputfolder, 'solutions', 'State Match Mean p_HEA01Rest00.asc'))

            #Asserts that coverage, duration, and occurrence measures are as expected based on the simulated data
            self.assertAlmostEqual(cov_dataset_0, 0.24901960784313726)
            self.assertAlmostEqual(cov_dataset_1, 0.25098039215686274)
            self.assertAlmostEqual(cov_dataset_2, 0.25098039215686274)
            self.assertAlmostEqual(cov_dataset_3, 0.24901960784313726)

            self.assertAlmostEqual(dur_dataset_0, 3.90625)
            self.assertAlmostEqual(dur_dataset_1, 3.90625)
            self.assertAlmostEqual(dur_dataset_2, 3.90625)
            self.assertAlmostEqual(dur_dataset_3, 3.90625)

            self.assertAlmostEqual(occ_dataset_0, 63.74901960784314)
            self.assertAlmostEqual(occ_dataset_1, 64.25098039215686)
            self.assertAlmostEqual(occ_dataset_2, 64.25098039215686)
            self.assertAlmostEqual(occ_dataset_3, 63.74901960784314)         


            #Asserts that individu states and state match mean (mean across dissimilarities for each matched GFP peak) are correct
            self.assertAlmostEqual(dataset1.all(), individu_correct_solution.all())  
            self.assertEqual(g4['State Match Mean p'].attrs['State Match Mean p based on'], 'Mean dissimilarity')
            self.assertAlmostEqual(dataset2.all(), state_match_mean_correct_solution.all())     
            
              
        