# -*- coding: utf-8 -*-

#############################
###        CLASSES        ###
#############################

###Class that defines where in the filename the group, pt, cond, run information are


class FileNameDescription(object):
    def __init__(self, group_indices, pt_indices, cond_indices, run_indices, file_ending):
        self.group = group_indices
        self.pt = pt_indices	
        self.cond = cond_indices
        self.run = run_indices		
        self.file_ending = file_ending

###Class that defines where in the folder structure the group, pt, cond, run information are


class FolderStructureDescription(object):
    def __init__(self, group_indices, pt_indices, cond_indices, run_indices):
        self.group = group_indices
        self.pt = pt_indices	
        self.cond = cond_indices
        self.run = run_indices	

###Class that defines which infos are retrieved from filename and which infos are retrieved from folder structure


class FilenameFolderDescription(object):
    def __init__(self, has_group, has_participant, has_condition, has_run):
        self.group = has_group
        self.pt = has_participant	
        self.cond = has_condition
        self.run = has_run		


###Class that defines number of channels (nch), number of time frames per epoch (tf), sampling frequency (sf), channel list (chlist)


class EegInfo(object):
    def __init__(self, nch, tf, sf, chlist):
        if nch != len(chlist):
            raise AssertionError('You defined ' + str(nch) + ' number of channels. Your channel list contains ' + str(len(chlist)) + ' elements. Numbers must be equal.')
        if tf%sf != 0:
            raise AssertionError('You defined ' + str(tf) + ' number of time frames. Must be multiple of sampling frequency: ' + str(sf))
        self.nch = nch
        self.tf = tf
        self.sf = sf
        self.chlist = chlist


###Class that defines the group names, participant names, cond names, and run names


class StudyInfo(object):
    def __init__(self, group_dict):
        self.group_dict = group_dict
        self.group_list = group_dict.keys()