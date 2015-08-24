# -*- coding: utf-8 -*-

######################################
# test4_mstates.py sorts modelmaps ###
######################################


################################
# 1.) LOAD PACKAGES ###
################################

#external
import os
import os.path

#integrated in the key.py library
from keypy.preprocessing.helper_functions import *

####   Classes     ####
from keypy.microstates.classes import *


################################
# 2.) Specify data folder info ###
################################

library_path = os.path.dirname(os.path.abspath(__file__))

##To do Testing Script

#contains data loaded into hdf5 file with the preprocessing script
inputfolder = os.path.join(library_path, "..","..","data","test23")
#will be created to contain output hdf5 files from microstate processing
outputfolder = os.path.join(library_path, "..","data","test3_output")


if not os.path.exists(outputfolder):
    os.makedirs(outputfolder)

#name of hdf5 that contains data
inputhdf5 = os.path.join( inputfolder, 'all_recordings.hdf')


#####################################
# 3.) Create EEG info object      ###
#####################################

#contains study_info script
script_inputfolder = os.path.join(cwd, "tests\\5_modelmaps\\test_scripts")

#exclude before committ
#execfile(os.path.join( script_inputfolder, 'study_info_test1.py'))
#include before committ
execfile(os.path.join( script_inputfolder, 'a_mstates_preprocessing_test1.py'))


################################################
# 4.)  Get study info object from hdf5 file  ###
################################################

study_info_obj=create_study_info_obj_from_data(inputhdf5)


#######################################
# 5.) Create Configuration Object   ###
#######################################

confobj = mst_configuration(
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

microstate_input = 'mstate1'
microstate_output = 'microstate'

from keypy.microstates.microstates import * 

#include before commit
run_microstates(confobj, eeg_info_study_obj, inputhdf5, outputfolder, microstate_input, microstate_output)
#--------------------------------------------------------------------------------------------------------------------------------------------


#################
# 7.) #Run Modelmaps (run_modelmaps_for_modelmap_types computes modelmaps for all types selected)
#################

from keypy.microstates.modelmaps import * 

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
    outputfolder = os.path.join(cwd, "tests\\5_modelmaps\\test_data\\test3_output\\{0}" .format(series))
    if not os.path.exists(outputfolder):
        os.makedirs(outputfolder)

    run_model_maps_series(series, inputfolder, outputfolder, first_input, confobj)

#################
# 7.) #Run Sortmaps
#################

