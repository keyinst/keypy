# -*- coding: utf-8 -*-

from __future__ import print_function

from contextlib import closing
import os
import os.path

import h5py
import numpy

#############################
###      FUNCTIONS        ###
#############################

#############################################################################
### get group, pt, cond, run info based on filename and folder structure ###
#############################################################################


def get_file_info (filename, file_name_obj, filename_folder_obj):
    """
    Output the group, pt, cond, and run name for each filename.

    Parameters
    ----------
    filename : str
        Input file name, e.g. 'vp01_1.txt'.
    file_name_obj : object of type FileNameDescription
        Contains the following attributes: group (int range), pt (int range), cond (int range), run (int range), file_ending (string).
    filename_folder_obj : object of type FilenameFolderDescription
        Contains the following attributes: group, pt, cond, run {'folder', 'filename', 'none'}.

    Returns
    -------
        group: str
            name of the group
        pt: str
            name of the participant
        cond: str
            name of the condition
        run: strb
            name of the run
    """

    if filename_folder_obj.group == 'filename':
        if file_name_obj.group[1] >= len (filename):
            raise IndexError("For group, index " + str(file_name_obj.group[1]) + " was specified. But this index is out of range for filename " + filename)
        else:
           group = filename[file_name_obj.group[0]:file_name_obj.group[1]+1]
    else:
        group = None
    
    if filename_folder_obj.pt == 'filename':
        if file_name_obj.pt[1] >= len (filename):
            raise IndexError("For participant, index " + str(file_name_obj.pt[1]) + " was specified. But this index is out of range for filename " + filename)
        else:
           pt = filename[file_name_obj.pt[0]:file_name_obj.pt[1]+1]
    else:
        pt = None
        
    if filename_folder_obj.cond == 'filename':
        if file_name_obj.cond[1] >= len (filename):
            raise IndexError("For condition, index " + str(file_name_obj.cond[1]) + " was specified. But this index is out of range for filename " + filename)
        else:
           cond = filename[file_name_obj.cond[0]:file_name_obj.cond[1]+1]
    else:
        cond = None
        
    if filename_folder_obj.run == 'filename':
        if file_name_obj.run[1] >= len (filename):
            raise IndexError("For run, index " + str(file_name_obj.run[1]) + " was specified. But this index is out of range for filename " + filename)
        else:
           run = filename[file_name_obj.run[0]:file_name_obj.run[1]+1]
    else:
        run = None            
      
    return group, pt, cond, run


def get_folder_info (folderlist, folder_structure_obj, filename_folder_obj):
    """
    Outputs the group, pt, cond, and run name based on the folder structure.

    Parameters
    ----------
    folderlist : list
        List of folders seperated by ",", e.g. ['groups', 'HC', 'Rest'].
    folder_structure_obj : object of type FolderStructureDescription
        Contains the following attributes: group (int), pt (int), cond (int), run (int).
    filename_folder_obj : object of type FilenameFolderDescription
        Contains the following attributes: group, pt, cond, run {'folder', 'filename', 'none'}.

    Returns
    -------
        group: string
            name of the group
        pt: string
            name of the participant
        cond: string
            name of the condition
        run: string
            name of the run
    """

    if filename_folder_obj.group == 'folder':
        if folder_structure_obj.group >= len (folderlist):
            raise IndexError("no folder with index " + str(folder_structure_obj.group) + " which was specified for group")
        else:
            group = folderlist[folder_structure_obj.group]
    else:
        group = None


    if filename_folder_obj.pt == 'folder':
        if folder_structure_obj.pt >= len (folderlist):
            raise IndexError("no folder with index " + str(folder_structure_obj.pt) + " which was specified for participant")
        else:
            pt = folderlist[folder_structure_obj.pt]
    else:
            pt = None

    if filename_folder_obj.cond == 'folder':
        if folder_structure_obj.cond >= len (folderlist):
            raise IndexError("no folder with index " + str(folder_structure_obj.cond) + " which was specified for condition")
        else:
            cond = folderlist[folder_structure_obj.cond]
    else: 
        cond = None

    if filename_folder_obj.run == 'folder':

        if folder_structure_obj.run >= len (folderlist):
            raise IndexError("no folder with index " + str(folder_structure_obj.run) + " which was specified for run")
        else:
            run = folderlist[folder_structure_obj.run]
    else:
         run = None

    return group, pt, cond, run


########################################
###        Load data into hdf5       ###
########################################


def loaddata(inputfolder, outputhdf5, loaddata_output, file_name_obj, folder_structure_obj, filename_folder_obj, eeg_info_study_obj):    
    """
    Outputs a hdf 5 file containing a hdf5 group for each group, participant, condition, and run, and a dataset of the imported EEG for each run.

    Parameters
    ----------
    inputfolder : str
        Input folder path, e.g. 'C:\\Users\\...\\keypy\\example\\data\\input'.
    outputhdf5 : str
        Output folder hdf5 file path, e.g. 'C:\\Users\\...\\data\\input\\rawdata.hdf'
    loaddata_output : str
        Name of the dataset in the hdf5 file that is created for each EEG input file, e.g. 'rawdata'.
    file_name_obj : object of type FileNameDescription
        Contains the following attributes: group (int range), pt (int range), cond (int range), run (int range), file_ending (string).
    folder_structure_obj : object of type FolderStructureDescription
        Contains the following attributes: group (int), pt (int), cond (int), run (int).
    filename_folder_obj : object of type FilenameFolderDescription
        Contains the following attributes: group, pt, cond, run {'folder', 'filename', 'none'}.
    eeg_info_study_obj : object of type EegInfo
        Contains the following attributes: nch (number of channels), tf (number of time frames per epoch), sf (sampling frequency), chlist (channel list)
    """

    #creates list of tripples (filename, folders, path to file)
    allfiles_list = []
    for result in os.walk(inputfolder):
        for filename in result[2]:
            if filename.endswith(file_name_obj.file_ending):
                allfiles_list.append((filename, result[0][len(inputfolder)+1:].split(os.path.sep), result[0]))
            else:
                print('File named ', filename, ' ignored. Did not have file ending ', file_name_obj.file_ending, ' which was specified by the user.')

    with closing( h5py.File(outputhdf5, 'w') ) as f:
        #retrieves info about group, pt, cond, run for each file
        for fili in allfiles_list:
            group_file, pt_file, cond_file, run_file = get_file_info (fili[0],  file_name_obj, filename_folder_obj)
            group_folder, pt_folder, cond_folder, run_folder = get_folder_info (fili[1],  folder_structure_obj, filename_folder_obj)

            if filename_folder_obj.group == 'filename':
                group = group_file
            elif filename_folder_obj.group == 'folder':
                group = group_folder
            elif filename_folder_obj.group == 'none':
                group = 'All_PTs'
            else:
                print('fatal error')

            if filename_folder_obj.pt == 'filename':
                pt = pt_file
            elif filename_folder_obj.pt == 'folder':
                pt = pt_folder
            elif filename_folder_obj.pt == 'none':
                pt = 'PT'
            else:
                print('fatal error')

            if filename_folder_obj.cond == 'filename':
                cond = cond_file
            elif filename_folder_obj.cond == 'folder':
                cond = cond_folder
            elif filename_folder_obj.cond == 'none':
                cond = 'Cond'
            else:
                print('fatal error')

            if filename_folder_obj.run == 'filename':
                run = run_file
            elif filename_folder_obj.run == 'folder':
                run = run_folder
            elif filename_folder_obj.run == 'none':
                run = '1'
            else:
                print('fatal error')

            #write data to hdf5 file
            if not 'group_{0}' .format(group) in f['/'].keys():
                participant_group = f.create_group('group_{0}' .format(group))
            else:
                participant_group = f['group_{0}' .format(group)]

            if not ('pt_{0}' .format(pt)) in f['/group_%s' % (group)].keys():
                pti_group = participant_group.create_group( 'pt_{0}' .format(pt) )   
            else:
                pti_group = f['/group_%s/pt_%s' % (group, pt)]

            if not 'cond_{0}' .format(cond) in f['/group_%s/pt_%s' % (group, pt)].keys():
                cond_group = pti_group.create_group( 'cond_{0}' .format(cond) )   
            else:
                cond_group = f['/group_%s/pt_%s/cond_%s' % (group, pt, cond)]

            if not 'run_{0}' .format(run) in f['/group_%s/pt_%s/cond_%s' % (group, pt, cond)].keys():
                run_group = cond_group.create_group( 'run_{0}' .format(run) )   
            else:
                run_group = f['/group_%s/pt_%s/cond_%s/run_%s' % (group, pt, cond, run)]

            dataset = numpy.loadtxt(os.path.join( fili[2], '{0}' .format(fili[0])) )

            if dataset.shape[0]%eeg_info_study_obj.tf != 0:
                raise AssertionError('File ' + str(fili[0]) + ' in folder ' + str(fili[2]) + 'has ' + str(dataset.shape[0]) + ' number of rows. Must be a multiple of number of time frames: ' + str(eeg_info_study_obj.tf))
            
            if dataset.shape[1] != eeg_info_study_obj.nch:
                raise AssertionError('File ' + str(fili[0]) + ' in folder ' + str(fili[2]) + 'has ' + str(dataset.shape[1]) + ' number of columns. Must be a equal to number of channels: ' + str(eeg_info_study_obj.nch))

            run_group.create_dataset('%s' % (loaddata_output) , data = dataset)
