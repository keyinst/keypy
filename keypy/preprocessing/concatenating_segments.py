# -*- coding: utf-8 -*-

import glob, os
from sets import Set
import numpy as np


##Concatenation function

def concatenate(mydir, ntf, nch, id_start, id_end, ending_in, ending_out):
    """
    Outputs one file of all epochs of a particular group participant cond run which had previously been segmented into separate files.

    Parameters
    ----------
    mydir : str
        Input folder path, e.g. 'C:\\Users\\...\\keypy\\example\\data\\input'.
    ntf : int
        Number of time frames within each epoch (e.g. 512).
    nch : int
        Number of channels.
    id_start : object of type FileNameDescription
        Starting index of part of filename that is a unique identifier for each group participant cond run.
    id_end : object of type FolderStructureDescription
        Ending index of part of filename that is a unique identifier for each group participant cond run.
    ending_in : str
        File ending of inputfiles.
    ending_out : str
        Desired file ending of outputfiles.
    """

    ##get names for all files
    ##get list of unique group participant cond runs
    ##get dict with keys unique id and values list of corresponding filenames

    #id_list = set()
    file_list_by_id={}

    os.chdir(mydir)
    for file in glob.glob("*.{0}" .format(ending_in)):
        id=file[id_start:id_end]
        #id_list.add(id)
        if id in file_list_by_id:
            file_list_by_id[id].append(file)
        else:
            file_list_by_id[id] = []
            file_list_by_id[id].append(file)


    ##load all files of the same group participant cond run
    ##concatenate all epochs (in epoch order)
    ##export as *.asc file

    for idi in file_list_by_id.keys():
        conc_file = np.zeros((len(file_list_by_id[idi])*ntf, nch))
        for filenr, file in enumerate(file_list_by_id[idi]):
            file_content=np.loadtxt(file)
            conc_file[filenr*ntf:(filenr*ntf)+ntf]=file_content

        np.savetxt('{0}.{1}' .format(idi, ending_out), conc_file)


##Concatenation function for files of various length but same number of tfs

def concatenate_various(mydir, ntf, nch, id_start, id_end, ending_in, ending_out):
    """
    Outputs one file of all epochs of a particular group participant cond run which had previously been segmented into separate files.

    Parameters
    ----------
    mydir : str
        Input folder path, e.g. 'C:\\Users\\...\\keypy\\example\\data\\input'.
    ntf : int
        Number of time frames within each epoch (e.g. 512).
    nch : int
        Number of channels.
    id_start : object of type FileNameDescription
        Starting index of part of filename that is a unique identifier for each group participant cond run.
    id_end : object of type FolderStructureDescription
        Ending index of part of filename that is a unique identifier for each group participant cond run.
    ending_in : str
        File ending of inputfiles.
    ending_out : str
        Desired file ending of outputfiles.
    """

    ##get names for all files
    ##get list of unique group participant cond runs
    ##get dict with keys unique id and values list of corresponding filenames

    #id_list = set()
    file_list_by_id={}

    os.chdir(mydir)
    for file in glob.glob("*.{0}" .format(ending_in)):
        id=file[id_start:id_end]
        #id_list.add(id)
        if id in file_list_by_id:
            file_list_by_id[id].append(file)
        else:
            file_list_by_id[id] = []
            file_list_by_id[id].append(file)


    ##load all files of the same group participant cond run
    ##concatenate all epochs (in epoch order)
    ##export as *.asc file

    

    for idi in file_list_by_id.keys():
        file_len=[]
        for filenr, file in enumerate(file_list_by_id[idi]):
            file_content=np.loadtxt(file)
            file_len.append(file_content.shape[0])

            if 'conc_file' in locals():                
                conc_file=np.concatenate((conc_file, file_content), axis=0)
            else:
                conc_file = file_content

        assert(np.sum(file_len)==conc_file.shape[0])
        print 'Concatenating', idi, ' File length: ', np.sum(file_len)

        np.savetxt('{0}.{1}' .format(idi, ending_out), conc_file)

    print 'Concatenating succesfully completed.'


