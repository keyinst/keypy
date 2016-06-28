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

from keypy.microstates.modelmaps import * 
from keypy.microstates.configuration import *
from keypy.microstates.meanmods import *
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
inputfolder = os.path.join(library_path,"data","input","dataset1", "groups")

# the name of the generated HDF5 file containing all the data
hdf5_filename = 'all_recordings.hdf'

# the folder path for the output where the HDF5 file is stored
outputfolder = os.path.join(library_path,"data","output","dataset1")

# the folder path of the external .asc file that contains the microstate maps that the obtained maps are to be sorted by
sortbyfolder = os.path.join(library_path,"data","sortby")

#filename of external .asc file that contains the microstate maps that the obtained maps are to be sorted by (and corresponding channel list)
sortbyfile_external = os.path.join(sortbyfolder, 'mean_models_milz_et_al_2016.asc')
sortbychlist_external = os.path.join(sortbyfolder, 'mean_models_milz_et_al_2016_chlist.asc')

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


### Specify microstate analysis specific information

# Specify the number of maps
map_number = 4

### Specify the number of reruns for clustering

# the more reruns the better, we suggest at least 3 times the number of GFP peaks you would like to compute your microstates based on (default: 200)
user_defined_reruns_modmaps = 100

# the more reruns the better (but slows computation down), we suggest approximately 4 times your number of participants (default: 50)
user_defined_reruns_meanmods = 50


#################################################################################################
###    Advanced Configuration Options    (only modify if you truly know what you are doing    ###
#################################################################################################

# Specify the series you would like to compute means based on, to sort your maps by, and to compute the parameters by

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

#multiple series at once possible (for modelmaps and sortmaps only), those used in sortmaps must be contained in modelmaps, those used in parameters must be contained in sortmaps

#compute means based on the following series
user_defined_series_meanmods=['Series_4']
#sort maps based on the following series
user_defined_series_sortmaps = ['Series_4']

#Additional settings for parameter computation
#specify whether you would like to compute parameters based on an external file, a level of a hdf series, or your hdf5_filename (default: 'all_recordings.hdf')
user_defined_parameter_type = 'series' #can be external, series, inputhdf

#compute parameters based on particular level of the following series
user_defined_series_parameters = 'Series_4'

#specify the hdf file the parameters should be computed upon
user_defined_sortbyfile_parameters='meanmods_across_groups_sorted.hdf'
#specify the dataset name that the model maps in the above hdf are stored in
user_defined_sortbydataset_parameters='meanmod'

#specify the number of layers of your sortbyfile
#parameter_by = '4Levels'
#parameter_by = '3Levels'
#parameter_by = '2Levels'
user_defined_parameter_by_parameters='1Level'
	
#for inputhdf sort only    
user_defined_individual_level_sort_dataset = 'modelmap_Series_4_sorted'	
	
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
                        original_nr_of_maps = map_number,
                        seed_number = user_defined_reruns_modmaps,
                        max_number_of_iterations = 100,
                        ERP = False,
                        correspondance_cutoff = 0.00)

#################
# 6.) #Run Modmaps (computes 1 microstate for each dataset in inputfile)
#################

######
###Define input processing stage and output hdf5 file group
######

modmaps_input = 'mstate1'
modmaps_output = 'modelmap'

run_modmaps(confobj, eeg_info_study_obj, inputhdf5, modmaps_input, modmaps_output)

#--------------------------------------------------------------------------------------------------------------------------------------------


#################
# 7.) #Run meanmods (run_meanmods_for_modelmap_types computes meanmods for all types selected)
#################

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
                        original_nr_of_maps = map_number,
                        seed_number = user_defined_reruns_meanmods,
                        max_number_of_iterations = 200,
                        ERP = False,
                        correspondance_cutoff = 0.00)

series_versions = user_defined_series_meanmods

first_meanmods_series_input = modmaps_output

inputfolder = outputfolder

for series in series_versions:
    first_input = first_meanmods_series_input

    #create folder with name of series as outputfolder
    outputfolder_series = os.path.join(outputfolder,"{0}".format(series))
    if not os.path.exists(outputfolder_series):
        os.makedirs(outputfolder_series)

    run_meanmods_series(series, inputfolder, hdf5_filename, outputfolder_series, first_input, confobj)

#--------------------------------------------------------------------------------------------------------------------------------------------

#################
# 7.) #Run Sortmaps (run_sortmaps_for_sortmap_types computes sortmaps for all types selected)
#################

series_versions = user_defined_series_sortmaps

first_input = modmaps_output
sortbyfile = sortbyfile_external
sortbyfile_chlist = sortbychlist_external

for series in series_versions:
    run_sortmaps_series(series, inputfolder, hdf5_filename, sortbyfolder, sortbyfile, sortbyfile_chlist, outputfolder, first_input, confobj, eeg_info_study_obj)  


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

##info needed to know which data the parameters are to be "sorted" upon (meanmods)
sortbyfolder = sortbyfolder
parameter_type = user_defined_parameter_type #can be external, series, inputhdf

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
    sortbyseries = user_defined_series_parameters
    sortbyfile = user_defined_sortbyfile_parameters
    #sortbydataset = 'microstate_Series_1_sorted'
    sortbydataset = user_defined_sortbydataset_parameters
    external_chlist = False

    #specify the number of layers of your sortbyfile
    #parameter_by = '4Levels'
    #parameter_by = '3Levels'
    #parameter_by = '2Levels'
    parameter_by = user_defined_parameter_by_parameters
	

###
#if you use input hdf5 file sort
###
elif parameter_type == 'inputhdf':
    parameter_by = 'own_hdf'
    sortbyfile = hdf5_filename
    sortbydataset = user_defined_individual_level_sort_dataset
    sortbyseries = None
    external_chlist = False
else:
    'Error, type not correctly specified for parameter computation.'


data_provider=get_data_provider_for_parameter_by(parameter_by, inputfolder, hdf5_filename, inputdataset, sortbyfolder, sortbyfile, sortbydataset, sortbyseries, external_chlist)

run_parameters(data_provider, confobj, eeg_info_study_obj)
#--------------------------------------------------------------------------------------------------------------------------------------------