###########################
###    Load packages    ###
###########################

import os
import os.path

from keypy.preprocessing.file_info_classes import *
from keypy.preprocessing.data_loading import *
from keypy.preprocessing.avg_referencing import *
from keypy.preprocessing.filtering import *
from keypy.preprocessing.helper_functions import *

from keypy.microstates.microstates import * 
from keypy.microstates.configuration import *
#from keypy.microstates.modelmaps import *
#from keypy.microstates.sortmaps import *
#--------------------------------------------------------------------------------------------------------------------------------------------

###################################
###    Configuration Options    ###
###################################
# Please adapt the following parameters to your data or run through as is to see the processing of the example input.

# enter number of channels
nch = 64

# enter number of time frames for each segment
tf = 512

# enter sampling rate of your data
sf = 256

# enter your channel list in the same order as in your files
chlist = ['FP1','AF7','AF3','FP1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3','FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1','P3','P5','P7','PO7','PO3','O1','Oz','POz','Pz','CPz','FPz','FP2','AF8','AF4','AFz','Fz','F2','F4','F6','F8','FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8','TP8','CP6','CP4','CP2','P2','P4','P6','P8','PO8','PO4','O2']

# the folder path to the raw EEG data
library_path = os.path.dirname(os.path.abspath(__file__))
inputfolder = os.path.join(library_path,"data","input","groups")

# the name of the generated HDF5 file containing all the data
hdf5_filename = 'all_recordings.hdf'

# the folder path for the output where the HDF5 file is stored
outputfolder = os.path.join(library_path,"data","output")

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

inputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')

boxkeyfilter(inputhdf5, eeg_info_study_obj, filter_input, filter_settings, enable_detrending = False)

#--------------------------------------------------------------------------------------------------------------------------------------------

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

microstate_input = 'mstate1'
microstate_output = 'microstate'

#fixed_seed = 100

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

series_versions = ['Series_3']

first_modelmap_series_input = microstate_input

inputfolder = outputfolder

for series in series_versions:
    first_input = first_modelmap_series_input

    #create folder with name of series as outputfolder
    outputfolder = os.path.join(outputfolder,"{0}".format(series))
    if not os.path.exists(outputfolder):
        os.makedirs(outputfolder)

    run_model_maps_series(series, inputfolder, outputfolder, first_input, confobj)

#--------------------------------------------------------------------------------------------------------------------------------------------

#################
# 7.) #Run Sortmaps (run_sortmaps_for_sortmap_types computes sortmaps for all types selected)
#################

##Series 3 complete processing

#Step 1
inputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_groups.hdf')
sortbyhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_conds.hdf')
outputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_groups{0}.hdf' .format('_sorted') )
inputdataset = 'modelmap'
sortbydataset = 'modelmap'
outputdatset = 'sorted_modelmap'
sortdata_provider = SortGroupsByAllDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset)

#Step 2
inputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_pts.hdf')
sortbyhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_groups_sorted.hdf')
outputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_pts{0}.hdf' .format('_sorted') )
inputdataset = 'modelmap'
sortbydataset = 'modelmap'
outputdatset = 'sorted_modelmap'
sortdata_provider = SortGroupCondByGroupDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset)

#Step 3
inputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_pts.hdf')
sortbyhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_groups_sorted.hdf')
outputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_pts{0}.hdf' .format('_sorted') )
inputdataset = 'modelmap'
sortbydataset = 'sorted_modelmap'
outputdatset = 'sorted_modelmap'
sortdata_provider = SortGroupCondByCondDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset)

#Step 4
inputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_runs.hdf')
sortbyhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_pts_sorted.hdf')
outputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_runs{0}.hdf' .format('_sorted') )
inputdataset = 'modelmap'
sortbydataset = 'sorted_modelmap'
outputdatset = 'sorted_modelmap'

sortdata_provider = SortGroupPtCondByGroupCondDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset)

confobj = MstConfiguration()

run_sort_maps(sortdata_provider, find_model_maps, confobj)

