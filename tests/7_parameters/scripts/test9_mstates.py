# -*- coding: utf-8 -*-

###############################################################
# test9_mstates.py computes a series of modelmaps at a time ###
###############################################################


################################
# 1.) LOAD PACKAGES ###
################################

#external
import os
import os.path
import unittest
import h5py
import filecmp

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


class Test_test9_mstates(unittest.TestCase):
    def test9_mstates(self):
        ################################
        # 2.) Specify data folder info ###
        ################################

        library_path = os.path.dirname(os.path.abspath(__file__))

        #contains data loaded into hdf5 file with the preprocessing script
        inputfolder = os.path.join(library_path, "..","data","test9")
        #will be created to contain output hdf5 files from microstate processing
        outputfolder = os.path.join(library_path, "..","data","test9_output")

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
        nch=19
        #enter number of time frames for each segment
        tf = 512
        #enter sampling rate
        sf = 256
        # enter your channel list in the same order as in your files
        chlist=['FP1', 'FP2', 'F7', 'F3', 'Fz', 'F4', 'F8', 'T7', 'C3', 'Cz', 'C4', 'T8', 'P7' ,'P3' ,'Pz', 'P4', 'P8', 'O1', 'O2']

    

        #create object of EEG info information
        eeg_info_study_obj = EegInfo(nch, tf, sf, chlist)

        ################################
        ### Specify data folder info ###
        ################################

        inputfolder = os.path.join(library_path, "..\\data\\test9")
        outputfolder = os.path.join(library_path, "..\\data\\test9_output")

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
        modmaps_output = 'modelmap'


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

        outputfolder = os.path.join(library_path, "..\\data\\test9_output")
        inputfolder = outputfolder

        for series in series_versions:
            first_input = 'microstate'

            #create folder with name of series as outputfolder
            outputfolder_series = os.path.join(library_path, "..\\data\\test9_output\\{0}" .format(series))
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
        outputfolder = os.path.join(library_path, "..\\data\\test9_output")

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
        inputfolder = os.path.join(library_path, "..\\data\\test9_output")
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
            sortbyfile ='mean_models_koenig_et_al_2002.asc'
            external_chlist='mean_models_koenig_et_al_2002_chlist.asc'
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
        with closing (h5py.File(os.path.join(sortbyfolder, 'mean_models_koenig_et_al_2002', 'mstate_parameters.hdf5'), 'r')) as f:
            g1 = f['/group_HEA']
            g2 = g1['pt_01']
            g3 = g2['cond_Rest']
            g4 = g3['run_0']
            dataset1 = g4['Individual_States'][:]
            dataset2 = g4['State Match Mean p'][:]
            #epoch 0
            g5 = g4['ep_000']
            GFP_curve=g5['GFP Curve'][:]
            Start_state_array=g5['Start state array'][:]
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
            gfp_mean_across_all_tfs_0 = map_00['GFP Mean across all tfs'].value
            gfp_mean_across_all_tfs_1 = map_01['GFP Mean across all tfs'].value
            gfp_mean_across_all_tfs_2 = map_02['GFP Mean across all tfs'].value
            gfp_mean_across_all_tfs_3 = map_03['GFP Mean across all tfs'].value
            SD_dur_in_ms_0 = map_00['SD duration in ms'].value
            SD_dur_in_ms_1 = map_01['SD duration in ms'].value
            SD_dur_in_ms_2 = map_02['SD duration in ms'].value
            SD_dur_in_ms_3 = map_03['SD duration in ms'].value    
            numb_of_ms_for_each_state_0 = map_00['number of ms for each state'].value
            numb_of_ms_for_each_state_1 = map_01['number of ms for each state'].value 
            numb_of_ms_for_each_state_2 = map_02['number of ms for each state'].value 
            numb_of_ms_for_each_state_3 = map_03['number of ms for each state'].value     

            #epoch 1
            g6 = g4['ep_001']
            GFP_curve_ep01=g6['GFP Curve'][:]
            Start_state_array_ep01=g6['Start state array'][:]
            map_00_ep01 = g6['map_00']
            map_01_ep01 = g6['map_01']
            map_02_ep01 = g6['map_02']
            map_03_ep01 = g6['map_03'] 
            cov_dataset_0_ep01 = map_00_ep01['Coverage in percent'].value
            cov_dataset_1_ep01 = map_01_ep01['Coverage in percent'].value
            cov_dataset_2_ep01 = map_02_ep01['Coverage in percent'].value
            cov_dataset_3_ep01 = map_03_ep01['Coverage in percent'].value
            dur_dataset_0_ep01 = map_00_ep01['Mean duration in ms'].value
            dur_dataset_1_ep01 = map_01_ep01['Mean duration in ms'].value
            dur_dataset_2_ep01 = map_02_ep01['Mean duration in ms'].value
            dur_dataset_3_ep01 = map_03_ep01['Mean duration in ms'].value
            occ_dataset_0_ep01 = map_00_ep01['Occurrance per s'].value
            occ_dataset_1_ep01 = map_01_ep01['Occurrance per s'].value
            occ_dataset_2_ep01 = map_02_ep01['Occurrance per s'].value
            occ_dataset_3_ep01 = map_03_ep01['Occurrance per s'].value
            gfp_mean_across_all_tfs_0_ep01 = map_00_ep01['GFP Mean across all tfs'].value
            gfp_mean_across_all_tfs_1_ep01 = map_01_ep01['GFP Mean across all tfs'].value
            gfp_mean_across_all_tfs_2_ep01 = map_02_ep01['GFP Mean across all tfs'].value
            gfp_mean_across_all_tfs_3_ep01 = map_03_ep01['GFP Mean across all tfs'].value         
            
            #epoch 2
            g7 = g4['ep_002']
            GFP_curve_ep02=g7['GFP Curve'][:]
            Start_state_array_ep02=g7['Start state array'][:]
            map_00_ep02 = g7['map_00']
            map_01_ep02 = g7['map_01']
            map_02_ep02 = g7['map_02']
            map_03_ep02 = g7['map_03'] 
            cov_dataset_0_ep02 = map_00_ep02['Coverage in percent'].value
            cov_dataset_1_ep02 = map_01_ep02['Coverage in percent'].value
            cov_dataset_2_ep02 = map_02_ep02['Coverage in percent'].value
            cov_dataset_3_ep02 = map_03_ep02['Coverage in percent'].value
            dur_dataset_0_ep02 = map_00_ep02['Mean duration in ms'].value
            dur_dataset_1_ep02 = map_01_ep02['Mean duration in ms'].value
            dur_dataset_2_ep02 = map_02_ep02['Mean duration in ms'].value
            dur_dataset_3_ep02 = map_03_ep02['Mean duration in ms'].value
            occ_dataset_0_ep02 = map_00_ep02['Occurrance per s'].value
            occ_dataset_1_ep02 = map_01_ep02['Occurrance per s'].value
            occ_dataset_2_ep02 = map_02_ep02['Occurrance per s'].value
            occ_dataset_3_ep02 = map_03_ep02['Occurrance per s'].value
            gfp_mean_across_all_tfs_0_ep02 = map_00_ep02['GFP Mean across all tfs'].value
            gfp_mean_across_all_tfs_1_ep02 = map_01_ep02['GFP Mean across all tfs'].value
            gfp_mean_across_all_tfs_2_ep02 = map_02_ep02['GFP Mean across all tfs'].value
            gfp_mean_across_all_tfs_3_ep02 = map_03_ep02['GFP Mean across all tfs'].value                
              
            #epoch 3
            g8 = g4['ep_003']
            GFP_curve_ep03=g8['GFP Curve'][:]
            Start_state_array_ep03=g8['Start state array'][:]
            map_00_ep03 = g8['map_00']
            map_01_ep03 = g8['map_01']
            map_02_ep03 = g8['map_02']
            map_03_ep03 = g8['map_03'] 
            cov_dataset_0_ep03 = map_00_ep03['Coverage in percent'].value
            cov_dataset_1_ep03 = map_01_ep03['Coverage in percent'].value
            cov_dataset_2_ep03 = map_02_ep03['Coverage in percent'].value
            cov_dataset_3_ep03 = map_03_ep03['Coverage in percent'].value
            dur_dataset_0_ep03 = map_00_ep03['Mean duration in ms'].value
            dur_dataset_1_ep03 = map_01_ep03['Mean duration in ms'].value
            dur_dataset_2_ep03 = map_02_ep03['Mean duration in ms'].value
            dur_dataset_3_ep03 = map_03_ep03['Mean duration in ms'].value
            occ_dataset_0_ep03 = map_00_ep03['Occurrance per s'].value
            occ_dataset_1_ep03 = map_01_ep03['Occurrance per s'].value
            occ_dataset_2_ep03 = map_02_ep03['Occurrance per s'].value
            occ_dataset_3_ep03 = map_03_ep03['Occurrance per s'].value
            gfp_mean_across_all_tfs_0_ep03 = map_00_ep03['GFP Mean across all tfs'].value
            gfp_mean_across_all_tfs_1_ep03 = map_01_ep03['GFP Mean across all tfs'].value
            gfp_mean_across_all_tfs_2_ep03 = map_02_ep03['GFP Mean across all tfs'].value
            gfp_mean_across_all_tfs_3_ep03 = map_03_ep03['GFP Mean across all tfs'].value   
            
                        
            #load correct solutions
            individu_correct_solution = np.loadtxt(os.path.join(inputfolder, 'solutions', 'Individual_States_HEA01Rest00.asc'))
            state_match_mean_correct_solution = np.loadtxt(os.path.join(inputfolder, 'solutions', 'State Match Mean p_HEA01Rest00.asc'))
            numb_ms_each_state_solution = np.loadtxt(os.path.join(inputfolder, 'solutions', 'number of ms for each state_HEA01Rest00.asc'))
            gfp_curve_solution = np.loadtxt(os.path.join(inputfolder, 'solutions', 'GFP Curve_HEA01Rest00.asc'))
            start_state_array_solution = np.loadtxt(os.path.join(inputfolder, 'solutions', 'Start state array_HEA01Rest00.asc'))


            #Asserts that coverage, duration, and occurrence measures are as expected based on the simulated data
            #coverage
            #epoch0
            self.assertAlmostEqual(cov_dataset_0, 0.24901960784313726)
            self.assertAlmostEqual(cov_dataset_1, 0.25098039215686274)
            self.assertAlmostEqual(cov_dataset_2, 0.25098039215686274)
            self.assertAlmostEqual(cov_dataset_3, 0.24901960784313726)

            #epoch1
            self.assertAlmostEqual(cov_dataset_0_ep01, 0.5)
            self.assertAlmostEqual(cov_dataset_1_ep01, 0.5)
            self.assertAlmostEqual(cov_dataset_2_ep01, 0)
            self.assertAlmostEqual(cov_dataset_3_ep01, 0)

            #epoch2
            self.assertAlmostEqual(cov_dataset_0_ep02, 0.5)
            self.assertAlmostEqual(cov_dataset_1_ep02, 0.25098039215686274)
            self.assertAlmostEqual(cov_dataset_2_ep02, 0.24901960784313726)
            self.assertAlmostEqual(cov_dataset_3_ep02, 0)

            #epoch3
            self.assertAlmostEqual(cov_dataset_0_ep03, 0)
            self.assertAlmostEqual(cov_dataset_1_ep03, 0)
            self.assertAlmostEqual(cov_dataset_2_ep03, 0)
            self.assertAlmostEqual(cov_dataset_3_ep03, 1)


            #duration
            #epoch0
            self.assertAlmostEqual(dur_dataset_0, 3.90625)
            self.assertAlmostEqual(dur_dataset_1, 3.90625)
            self.assertAlmostEqual(dur_dataset_2, 3.90625)
            self.assertAlmostEqual(dur_dataset_3, 3.90625)

            #epoch1
            self.assertAlmostEqual(dur_dataset_0_ep01, 996.09375)
            self.assertAlmostEqual(dur_dataset_1_ep01, 996.09375)
            self.assertAlmostEqual(dur_dataset_2_ep01, 0)
            self.assertAlmostEqual(dur_dataset_3_ep01, 0)

            #epoch2
            self.assertAlmostEqual(dur_dataset_0_ep02, 7.781982421875)
            self.assertAlmostEqual(dur_dataset_1_ep02, 3.90625)
            self.assertAlmostEqual(dur_dataset_2_ep02, 3.90625)
            self.assertAlmostEqual(dur_dataset_3_ep02, 0)

            #epoch3
            self.assertAlmostEqual(dur_dataset_0_ep03, 0)
            self.assertAlmostEqual(dur_dataset_1_ep03, 0)
            self.assertAlmostEqual(dur_dataset_2_ep03, 0)
            self.assertAlmostEqual(dur_dataset_3_ep03, 1992.1875)

            #occurrence
            #epoch0
            self.assertAlmostEqual(occ_dataset_0, 63.74901960784314)
            self.assertAlmostEqual(occ_dataset_1, 64.25098039215686)
            self.assertAlmostEqual(occ_dataset_2, 64.25098039215686)
            self.assertAlmostEqual(occ_dataset_3, 63.74901960784314)   
            
            #epoch1
            self.assertAlmostEqual(occ_dataset_0_ep01, 0.5019607843137255)
            self.assertAlmostEqual(occ_dataset_1_ep01, 0.5019607843137255)
            self.assertAlmostEqual(occ_dataset_2_ep01, 0)
            self.assertAlmostEqual(occ_dataset_3_ep01, 0)              
            
            #epoch2
            self.assertAlmostEqual(occ_dataset_0_ep02, 64.25098039215686)
            self.assertAlmostEqual(occ_dataset_1_ep02, 64.25098039215686)
            self.assertAlmostEqual(occ_dataset_2_ep02, 63.74901960784314)
            self.assertAlmostEqual(occ_dataset_3_ep02, 0)              
            
            #epoch3
            self.assertAlmostEqual(occ_dataset_0_ep03, 0)
            self.assertAlmostEqual(occ_dataset_1_ep03, 0)
            self.assertAlmostEqual(occ_dataset_2_ep03, 0)
            self.assertAlmostEqual(occ_dataset_3_ep03, 0.5019607843137255)   
            
            #gfp_mean_across_all_tfs
            #epoch0            
            self.assertAlmostEqual(gfp_mean_across_all_tfs_0, 0.19267056472289962)
            self.assertAlmostEqual(gfp_mean_across_all_tfs_1, 0.18584619367048433)
            self.assertAlmostEqual(gfp_mean_across_all_tfs_2, 0.2039207830558706)
            self.assertAlmostEqual(gfp_mean_across_all_tfs_3, 0.18466181801648643)

            #epoch1          
            self.assertAlmostEqual(gfp_mean_across_all_tfs_0_ep01, 0.19265109706641267)
            self.assertAlmostEqual(gfp_mean_across_all_tfs_1_ep01, 0.18550619594198267)
            self.assertAlmostEqual(gfp_mean_across_all_tfs_2_ep01, '-')
            self.assertAlmostEqual(gfp_mean_across_all_tfs_3_ep01, '-')

            #epoch2            
            self.assertAlmostEqual(gfp_mean_across_all_tfs_0_ep02, 0.19245490948822006)
            self.assertAlmostEqual(gfp_mean_across_all_tfs_1_ep02, 0.18544981706632832)
            self.assertAlmostEqual(gfp_mean_across_all_tfs_2_ep02, 0.20396915368158294)
            self.assertAlmostEqual(gfp_mean_across_all_tfs_3_ep02, '-')

            #epoch3            
            self.assertAlmostEqual(gfp_mean_across_all_tfs_0_ep03, '-')
            self.assertAlmostEqual(gfp_mean_across_all_tfs_1_ep03, '-')
            self.assertAlmostEqual(gfp_mean_across_all_tfs_2_ep03, '-')
            self.assertAlmostEqual(gfp_mean_across_all_tfs_3_ep03, 0.18450303759077272)


            #epoch0 additional measures
            self.assertAlmostEqual(SD_dur_in_ms_0, 0)
            self.assertAlmostEqual(SD_dur_in_ms_1, 0)
            self.assertAlmostEqual(SD_dur_in_ms_2, 0)
            self.assertAlmostEqual(SD_dur_in_ms_3, 0)        

            #Asserts that individu states and state match mean (mean across dissimilarities for each matched GFP peak) are correct
            self.assertAlmostEqual(dataset1.all(), individu_correct_solution.all())  
            self.assertEqual(g4['State Match Mean p'].attrs['State Match Mean p based on'], 'Mean dissimilarity')
            self.assertAlmostEqual(dataset2.all(), state_match_mean_correct_solution.all()) 

            self.assertAlmostEqual(numb_of_ms_for_each_state_0.all(), numb_ms_each_state_solution.all())
            self.assertAlmostEqual(GFP_curve.all(), gfp_curve_solution.all())
            self.assertAlmostEqual(numb_of_ms_for_each_state_0.all(), numb_ms_each_state_solution.all())
            
            #Mstate Label List Check
            correct_solution_R0 = os.path.join(inputfolder, 'solutions', 'mstate_label_list_cond_Rest_run_0_solution.csv')
            correct_solution_R1 = os.path.join(inputfolder, 'solutions', 'mstate_label_list_cond_Rest_run_1_solution.csv')
            correct_solution_R2 = os.path.join(inputfolder, 'solutions', 'mstate_label_list_cond_Rest_run_2_solution.csv')
            correct_solution_T0 = os.path.join(inputfolder, 'solutions', 'mstate_label_list_cond_Test_run_0_solution.csv')
            correct_solution_T1 = os.path.join(inputfolder, 'solutions', 'mstate_label_list_cond_Test_run_1_solution.csv')
            correct_solution_T2 = os.path.join(inputfolder, 'solutions', 'mstate_label_list_cond_Test_run_2_solution.csv')

            retrieved_solution_R0 = os.path.join(sortbyfolder, 'mean_models_koenig_et_al_2002', 'mstate_label_list_cond_Rest_run_0.csv')
            retrieved_solution_R1 = os.path.join(sortbyfolder, 'mean_models_koenig_et_al_2002', 'mstate_label_list_cond_Rest_run_1.csv')
            retrieved_solution_R2 = os.path.join(sortbyfolder, 'mean_models_koenig_et_al_2002', 'mstate_label_list_cond_Rest_run_2.csv')
            retrieved_solution_T0 = os.path.join(sortbyfolder, 'mean_models_koenig_et_al_2002', 'mstate_label_list_cond_Test_run_0.csv')
            retrieved_solution_T1 = os.path.join(sortbyfolder, 'mean_models_koenig_et_al_2002', 'mstate_label_list_cond_Test_run_1.csv')
            retrieved_solution_T2 = os.path.join(sortbyfolder, 'mean_models_koenig_et_al_2002', 'mstate_label_list_cond_Test_run_2.csv')

            self.assertTrue(filecmp.cmp(correct_solution_R0, retrieved_solution_R0))
            self.assertTrue(filecmp.cmp(correct_solution_R1, retrieved_solution_R1))
            self.assertTrue(filecmp.cmp(correct_solution_R2, retrieved_solution_R2))         
            self.assertTrue(filecmp.cmp(correct_solution_T0, retrieved_solution_T0))
            self.assertTrue(filecmp.cmp(correct_solution_T1, retrieved_solution_T1))
            self.assertTrue(filecmp.cmp(correct_solution_T2, retrieved_solution_T2))                       
                              
        #--------------------------------------------------------------------------------------------------------------------------------------------
        
