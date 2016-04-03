# -*- coding: utf-8 -*-

"""
Functions and classes that allow the preprocessing of artefact corrected EEG which have previously been segmented into equal length epochs to prepare them for example for microstate analysis.

data_loading : functions
    Loads EEG data from folder / folders into hdf5 file.
avg_referencing : script containing functions
    Computes the average reference of a specified processing stage dataset in the hdf5 file (e.g. rawdata).
filtering : functions
    Filters the data of a specified processing stage dataset in the hdf5 file (e.g. avg_ref).
helper_functions : functions
    Includes functions to create a study_info_obj from the data, to delete channels, and to delete participants.
file_info_classes : class definitions
    Includes the definition of the following classes: FileNameDescription, FolderStructureDescription, FilenameFolderDescription, EegInfo, StudyInfo
concatenating_segments : function
    Allows you to concatenate multiple epochs of one participant cond run into one file.
"""

#############################################################################

__all__ = ["avg_referencing", "filtering", "file_info_classes", "data_loading", "helper_functions", "concatenating_segments"]

