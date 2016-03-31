###########################
###    Load packages    ###
###########################

import os
import os.path
from os.path import basename

from keypy.preprocessing.file_info_classes import *
from keypy.preprocessing.data_loading import *
from keypy.preprocessing.avg_referencing import *
from keypy.preprocessing.filtering import *
from keypy.preprocessing.helper_functions import *

from keypy.microstates.microstates import * 
from keypy.microstates.configuration import *
from keypy.microstates.modelmaps import *
from keypy.microstates.sortmaps import *
from keypy.microstates.parameters import *
#--------------------------------------------------------------------------------------------------------------------------------------------

###################################
###    Configuration Options    ###
###################################
# Please adapt the following parameters to your data or run through as is to see the processing of the example input.

# enter number of channels
nch = 61

# enter number of time frames for each segment
tf = 512

# enter sampling rate of your data
sf = 256

# enter your channel list in the same order as in your files
chlist=['FP1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3','FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1','P3','P5','P7','PO7','PO3','O1','Oz','POz','Pz','CPz','FPz','FP2','AF8','AF4','AFz','Fz','F2','F4','F6','F8','FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8','TP8','CP6','CP4','CP2','P2','P4','P6','P8','PO8','PO4','O2']

# the folder path to the raw EEG data
library_path = os.path.dirname(os.path.abspath(__file__))
inputfolder = os.path.join(library_path,"data","input","groups")

# the name of the generated HDF5 file containing all the data
hdf5_filename = 'all_recordings.hdf'

# the folder path for the output where the HDF5 file is stored
outputfolder = os.path.join(library_path,"data","output")

# the folder path of the external .asc file that contains the microstate maps that the obtained maps are to be sorted by
sortbyfolder = os.path.join(library_path,"data","sortby")

#filename of external .asc file that contains the microstate maps that the obtained maps are to be sorted by (and corresponding channel list)
sortbyfile_external = os.path.join(sortbyfolder, 'mean_models_milz_et_al_2015.asc')
sortbychlist_external = os.path.join(sortbyfolder, 'mean_models_milz_et_al_2015_chlist.asc')

# Specify the information of your study
##Where in the filename is the following information (at which index of the string)? inclusive
##Where in the folder structure is the following information? 0 outtermost folder, 1 second folder, 2 third folder
##Specify for each element (group, partcipant, condition, run) whether the information is in the filename or folder structure.

"""
Each EEG data file in a study belongs to one of four dimentions: group, participant, condition and run.
Here, you need to specify the folder structure and filename convention that is used by the input data.
Each dimension is defined as:
	'none': this dimension does not exists in the study
	'folder': the dimension is represented as in the folder structure (e.g. one folder for each group)
	'filename': the dimension is represented in the filename (e.g. file name contains participant name "participant_001.txt"
"""

# group
has_group = 'folder'
group_indices_range = [0,1]
group_folder_level = 0

# participant
has_participant = 'filename'
participant_indices_range = [2,3]
participant_folder_level = 0

# condition
has_condition = 'folder'
condition_indices_range = [0,0]
condition_folder_level = 1

# run
has_run = 'filename'
run_indices_range = [5,5]
run_folder_level = 0

file_ending = 'txt'

# Specify the number of reruns for clustering

# the more reruns the better, we suggest at least 3 times the number of GFP peaks you would like to compute your microstates based on
user_defined_reruns_microstate = 100

# the more reruns the better (but slows computation down), we suggest approximately 4 times your number of participants
user_defined_reruns_modelmaps = 50
#--------------------------------------------------------------------------------------------------------------------------------------------

##########################
###    Run analysis    ###
##########################
# No changes necessary after this point.

if not os.path.exists(outputfolder):
    print("Create output folder: {0}".format(outputfolder))
    os.makedirs(outputfolder)

outputhdf5 = os.path.join( outputfolder, hdf5_filename)
loaddata_output = 'rawdata'

#create object of EEG info information
eeg_info_study_obj = EegInfo(nch, tf, sf, chlist)

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

inputhdf5 = os.path.join( outputfolder, hdf5_filename)

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


inputhdf5 = os.path.join( outputfolder, hdf5_filename)

boxkeyfilter(inputhdf5, eeg_info_study_obj, filter_input, filter_settings, enable_detrending = False)

#--------------------------------------------------------------------------------------------------------------------------------------------

#######################################
# 5.) Create Configuration Object   ###
#######################################




### Warning: Do not change ERP = False. Keypy has only been optimized to work with EEG and not ERP data.
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
                        seed_number = user_defined_reruns_microstate,
                        max_number_of_iterations = 100,
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

run_microstates(confobj, eeg_info_study_obj, inputhdf5, microstate_input, microstate_output)

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
                        method_GFPpeak = 'GFPL2',
                        original_nr_of_maps = 4,
                        seed_number = user_defined_reruns_modelmaps,
                        max_number_of_iterations = 200,
                        ERP = False,
                        correspondance_cutoff = 0.00)

series_versions = ['Series_3']

first_modelmap_series_input = microstate_output

inputfolder = outputfolder

for series in series_versions:
    first_input = first_modelmap_series_input

    #create folder with name of series as outputfolder
    outputfolder_series = os.path.join(outputfolder,"{0}".format(series))
    if not os.path.exists(outputfolder_series):
        os.makedirs(outputfolder_series)

    run_model_maps_series(series, inputfolder, hdf5_filename, outputfolder_series, first_input, confobj)

#--------------------------------------------------------------------------------------------------------------------------------------------

#################
# 7.) #Run Sortmaps (run_sortmaps_for_sortmap_types computes sortmaps for all types selected)
#################

confobj = MstConfiguration()

series_versions = ['Series_3']

first_input = 'microstate'
sortbyfile = "mean_models_milz_etal_2015.asc"
sortbyfile_chlist = "mean_models_milz_etal_2015_chlist.asc"


for series in series_versions:
    run_sort_maps_series(series, inputfolder, hdf5_filename, sortbyfolder, sortbyfile, sortbyfile_chlist, outputfolder, first_input, confobj, eeg_info_study_obj)  


#################
# 7.) #Run Parameters
#################

################################
#### Options for parameters ####
################################

#############################
#### Continuous EEG data ####
#############################

##info needed to know which data the parameters are to be computed upon
inputfolder = outputfolder
hdf5_filename = hdf5_filename
inputdataset = 'mstate1'

############################
####    Parameter by    ####
############################

##info needed to know which data the parameters are to be "sorted" upon (modelmaps)
sortbyfolder = sortbyfolder
parameter_type = 'series' #can be external, series, inputhdf
###
#if you use an external asci file to sort your data by
###
if parameter_type == 'external' : 
    parameter_by = 'external_norm'
    sortbyfile =sortbyfile_external
    external_chlist=sortbychlist_external
    sortbydataset = None
    sortbyseries = None

###
#if you use a previously computed hdf5 file (from Series)
###
elif parameter_type == 'series': 
    sortbyseries = 'Series_3'
    sortbyfile = 'modelmaps_across_conds_sorted.hdf'
    #sortbydataset = 'microstate_Series_1_sorted'
    sortbydataset = 'modelmap'
    external_chlist = False

    #specify the number of layers of your sortbyfile
    #parameter_by = '4Levels'
    #parameter_by = '3Levels'
    #parameter_by = '2Levels'
    parameter_by = '1Level'

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