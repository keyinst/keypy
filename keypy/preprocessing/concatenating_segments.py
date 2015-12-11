
import glob, os
from sets import Set
import numpy as np


##Concatination function

def concatenate(mydir, ntf, nch, id_start, id_end, ending):
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
    ending : object of type FilenameFolderDescription
        Desired file ending of outputfiles.
    """

    ##get names for all files
    ##get list of unique group participant cond runs
    ##get dict with keys unique id and values list of corresponding filenames

    #id_list = set()
    file_list_by_id={}

    os.chdir(mydir)
    for file in glob.glob("*.asc"):
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
        conc_file = np.zeros((ntf, nch))
        for file in file_list_by_id[idi]:
            file_content=np.loadtxt(file)
            if conc_file.all()==0:
                conc_file=file_content
            else:
                conc_file=np.vstack([conc_file,file_content])

        np.savetxt('{0}.{1}' .format(idi, ending), conc_file)
