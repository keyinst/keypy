
##################################
#######  Import Packages  ########
##################################

import os.path
import math
import itertools
import operator

import numpy as np
from scipy.stats import pearsonr
from contextlib import closing

import h5py



###Create Paths for sortmaps

class RunPath(object):
    def __init__(self, group, cond, run):
        self.group = group
        self.cond = cond
        self.run = run

class Levels1Path(object):
    def __init__(self, level1):
        self.level1 = level1

class Levels2Path(object):
    def __init__(self, level1, level2):
        self.level1 = level1
        self.level2 = level2

class Levels3Path(object):
    def __init__(self, level1, level2, level3):
        self.level1 = level1
        self.level2 = level2
        self.level3 = level3

###-----------------------

####Abstract Class SortDataProvider

class SortDataProvider(object):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset):
        self._file = inputhdf5
        self._sortbyfile = sortbyhdf5
        self._outputfile = outputhdf5
        self._inputdataset = inputdataset
        self._sortbydataset = sortbydataset
        self._outputdataset = outputdatset

####Sub Class CondDataProvider --> umbenennen

class SortGroupPtCondByGroupCondDataProvider1(SortDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset):
        SortDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths = []
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_level0 in f['/'].keys():
                for current_level1 in f['/{0}'.format(current_level0)].keys():
                    for current_level2 in f['/{0}/{1}'.format(current_level0, current_level1)].keys():
                        level3_path = Levels3Path(current_level0, current_level1, current_level2)
                        out_paths.append(level3_path)
        return out_paths
    
    ### group pt cond run
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted
    def get_input_data(self, output_path):
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '{0}/{1}/{2}' .format(output_path.level1, output_path.level2, output_path.level3)
            microstate_run_value = g['/{0}/{1}' .format(path, self._inputdataset)] 
            if all(microstate_run_value[0,:] == 0):
                print 'Error!', path, 'has all zeros', 'group, pt, cond ignored.'    
            else:
                model_map=microstate_run_value[:]
        return model_map

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        with closing( h5py.File(self._sortbyfile, 'r') ) as h:
            microstate_run_value = h['/{0}/{1}/{2}' .format(output_path.level1, output_path.level3, self._sortbydataset)] 

            if all(microstate_run_value[0,:] == 0):
                print 'Error!', path, current_run, 'has all zeros', 'group, pt, cond ignored.'    
            else:
                model_map=microstate_run_value[:]
        return model_map
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile, 'a') ) as k:
            print 'output_paths used for output', output_path.level1

            if output_path.level1 in k['/'].keys():
                group_group = k['{0}' .format(output_path.level1)]
            else:
                group_group = k['/'].create_group( '{0}' .format(output_path.level1)  ) 
                      
            if output_path.level2 in group_group.keys():
                pt_group = k['/{0}/{1}' .format(output_path.level1, output_path.level2)]
            else:
                pt_group = group_group.create_group( '{0}' .format(output_path.level2)  )   

            if output_path.level3 in pt_group.keys():
                run_group = k['/{0}/{1}/{2}' .format(output_path.level1, output_path.level2, output_path.level3)]
            else:
                run_group = pt_group.create_group( '{0}' .format(output_path.level3)  )   

            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                run_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in run_group.keys():
                print 'group, participant, condition already in outputfile, not recomputed', group_group, output_path.level1
            else:
                run_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)

class SortGroupCondByCondDataProvider1(SortDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset):
        SortDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths = []
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_level0 in f['/'].keys():
                for current_level1 in f['/{0}'.format(current_level0)].keys():
                    level2_path = Levels2Path(current_level0, current_level1)
                    out_paths.append(level2_path)
        return out_paths
    
    ### group pt cond run
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted
    def get_input_data(self, output_path):
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '{0}/{1}' .format(output_path.level1, output_path.level2)
            microstate_run_value = g['/{0}/{1}' .format(path, self._inputdataset)] 
            if all(microstate_run_value[0,:] == 0):
                print 'Error!', path, 'has all zeros', 'group, pt, cond ignored.'    
            else:
                model_map=microstate_run_value[:]
        return model_map

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        with closing( h5py.File(self._sortbyfile, 'r') ) as h:
            microstate_run_value = h['/{0}/{1}' .format(output_path.level2, self._sortbydataset)] 

            if all(microstate_run_value[0,:] == 0):
                print 'Error!', path, current_run, 'has all zeros', 'group, pt, cond ignored.'    
            else:
                model_map=microstate_run_value[:]
        return model_map
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile, 'a') ) as k:
            print 'output_paths used for output', output_path.level1

            if output_path.level1 in k['/'].keys():
                group_group = k['{0}' .format(output_path.level1)]
            else:
                group_group = k['/'].create_group( '{0}' .format(output_path.level1)  ) 
                      
            if output_path.level2 in group_group.keys():
                pt_group = h['/{0}/{1}' .format(output_path.level1, output_path.level2)]
            else:
                pt_group = group_group.create_group( '{0}' .format(output_path.level2)  )   


            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                pt_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in pt_group.keys():
                print 'group, participant, condition already in outputfile, not recomputed', group_group, output_path.level1
            else:
                pt_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)


class SortGroupCondByGroupDataProvider1(SortDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset):
        SortDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths = []
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_level0 in f['/'].keys():
                for current_level1 in f['/{0}'.format(current_level0)].keys():
                    level2_path = Levels2Path(current_level0, current_level1)
                    out_paths.append(level2_path)
        return out_paths
    
    ### group pt cond run
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted
    def get_input_data(self, output_path):
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '{0}/{1}' .format(output_path.level1, output_path.level2)
            microstate_run_value = g['/{0}/{1}' .format(path, self._inputdataset)] 
            if all(microstate_run_value[0,:] == 0):
                print 'Error!', path, 'has all zeros', 'group, pt, cond ignored.'    
            else:
                model_map=microstate_run_value[:]
        return model_map

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        with closing( h5py.File(self._sortbyfile, 'r') ) as h:
            microstate_run_value = h['/{0}/{1}' .format(output_path.level1, self._sortbydataset)] 

            if all(microstate_run_value[0,:] == 0):
                print 'Error!', path, current_run, 'has all zeros', 'group, pt, cond ignored.'    
            else:
                model_map=microstate_run_value[:]
        return model_map
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile, 'a') ) as k:
            print 'output_paths used for output', output_path.level1

            if output_path.level1 in k['/'].keys():
                group_group = k['{0}' .format(output_path.level1)]
            else:
                group_group = k['/'].create_group( '{0}' .format(output_path.level1)  ) 
                      
            if output_path.level2 in group_group.keys():
                pt_group = h['/{0}/{1}' .format(output_path.level1, output_path.level2)]
            else:
                pt_group = group_group.create_group( '{0}' .format(output_path.level2)  )   


            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                pt_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in pt_group.keys():
                print 'group, participant, condition already in outputfile, not recomputed', group_group, output_path.level1
            else:
                pt_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)

class SortGroupsByAllDataProvider1(SortDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset):
        SortDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths = []
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_level0 in f['/'].keys():
                level1_path = Levels1Path(current_level0)
                out_paths.append(level1_path)
        return out_paths
    
    ### group pt cond run
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted
    def get_input_data(self, output_path):
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '{0}' .format(output_path.level1)
            microstate_run_value = g['/{0}/{1}' .format(path, self._inputdataset)] 
            if all(microstate_run_value[0,:] == 0):
                print 'Error!', path, 'has all zeros', 'group, pt, cond ignored.'    
            else:
                model_map=microstate_run_value[:]
        return model_map

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        with closing( h5py.File(self._sortbyfile, 'r') ) as h:
            microstate_run_value = h['all/{0}' .format(self._sortbydataset)] 

            if all(microstate_run_value[0,:] == 0):
                print 'Error!', path, current_run, 'has all zeros', 'group, pt, cond ignored.'    
            else:
                model_map=microstate_run_value[:]
        return model_map
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile, 'a') ) as k:
            print 'output_paths used for output', output_path.level1

            if output_path.level1 in k['/'].keys():
                group_group = k['{0}' .format(output_path.level1)]
            else:
                group_group = k['/'].create_group( '{0}' .format(output_path.level1)  ) 
                      

            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                group_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in group_group.keys():
                print 'group, participant, condition already in outputfile, not recomputed', group_group, output_path.level1
            else:
                group_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)

###################OLD

class SortCondDataProvider_correct(SortDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset):
        SortDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths = []
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_group in f['/'].keys():
                group_group = f['/{0}' .format(current_group)]
                for current_pt in group_group.keys():
                    pt_group = f['/{0}/{1}' .format(current_group, current_pt)]
                    for current_cond in pt_group.keys():
                        cond_group = f['/{0}/{1}/{2}' .format(current_group, current_pt, current_cond)]
                        for current_run in cond_group.keys():
                            run_path = RunPath(current_group, current_pt, current_cond)
                            out_paths.append(run_path)
        return out_paths
    
    ### group pt cond run
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted
    def get_input_data(self, output_path):
        model_maps_all = []
        with closing( h5py.File(self._file, 'r') ) as f:            
            path = '/{0}/{1}/{2}/{3}' .format(output_path.group, output_path.pt, output_path.cond, output_path.run)

            microstate_run_value = f['/{0}/{1}/{2}' .format(path, microstate_run, self._inputdataset)] 

            if all(microstate_run_value[0,:] == 0):
                print 'Error!', path, current_run, 'has all zeros', 'group, pt, cond ignored.'    
            else:
                model_maps_all.append(microstate_run_value[:])
        return model_maps_all

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        model_maps_all = []
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '/{0}/{1}/{2}' .format(output_path.group, output_path.pt, output_path.cond)

            microstate_cond_value = g['/{0}/{1}/{2}' .format(path, self._outputdataset)] 
            if all(microstate_run_value[0,:] == 0):
                print 'Error!', path, current_run, 'has all zeros', 'group, pt, cond ignored.'    
            else:
                model_maps_all.append(microstate_run_value[:])
        return model_maps_all
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile) ) as h:
            print 'output_paths used for output', output_path.group, output_path.cond, output_path.pt

            if output_path.group in h['/'].keys():
                group_group = h['{0}' .format(output_path.group)]
            else:
                group_group = h['/'].create_group( '{0}' .format(output_path.group)  ) 
                      
            if output_path.pt in group_group.keys():
                pt_group = h['/{0}/{1}' .format(output_path.group, output_path.pt)]
            else:
                pt_group = group_group.create_group( '{0}' .format(output_path.pt)  )   

            if output_path.cond in pt_group.keys():
                cond_group = h['/{0}/{1}/{2}' .format(output_path.group, output_path.pt, output_path.cond)]
            else:
                cond_group = pt_group.create_group( '{0}' .format(output_path.cond)  )   

            if output_path.run in cond_group.keys():
                run_group = h['/{0}/{1}/{2}/{3}' .format(output_path.group, output_path.pt, output_path.cond, output_path.run)]
            else:
                run_group = cond_group.create_group( '{0}' .format(output_path.run)  )   


            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                run_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in run_group.keys():
                print 'group, participant, condition already in outputfile, not recomputed', group_group, output_path.pt, output_path.cond
            else:
                run_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)



##############################
########  Sort Maps  ########
##############################


def sort_maps(confobj, input, sortby):
    #old: (modelmaps, sortby_maps, original_nr_of_maps, modelmaps_re, ERP=False):
    """
    Parameters;
    modelmaps: modelmaps original channel nr
    sortby_maps: maps to sort modelmaps by
    original_nr_of_maps: e.g. 4 when you have 4 microstate maps
    modelmaps_re = maps to sort modelmaps by shrinked to TK map number, when TK_map is sorted by
    ERP: False if polarity of maps can be ignored for correlations

    Return Values:
    newraw: Maps ordered based on sortby_maps
    map_corr_list: Correlations between each map with the map it was labeled by
    """

    modelmaps = input
    sortby_maps = sortby
    original_nr_of_maps = confobj.original_nr_of_maps
    #modelmaps_re = 
    ERP = confobj.ERP

    ######
    ###Get best attribution matrix of my maps ordered by TK maps / model maps sorted of any stage
    ######

    attribution_matrix=[]
    mean_correlations = dict.fromkeys( range(original_nr_of_maps) )

    for ithperm, perm in enumerate(itertools.permutations((range(original_nr_of_maps)))):    
        pearsons=[]
        pearsons2=[]
        #changed from below to account for matchings with only 4 TK maps
        #for i in range(original_nr_of_maps):
        for i in range(len(sortby_maps[:,0])):
            pr, pp=pearsonr( modelmaps[i,:], sortby_maps[perm[i],:])
            '''
            if modelmaps_re.shape == sortby_maps.shape:
                pr, pp=pearsonr( modelmaps_re[i,:], sortby_maps[perm[i],:])
            else:
                pr, pp=pearsonr( sortby_maps[i,:], modelmaps_re[perm[i],:])
            '''
            #print 'pr', pr
            pearsons.append(pr)
            if ERP:
                pearsons2.append(float(pearsons[i]))            #NOETIG, da neg. Korr nicht gut fuer ERP
            else:
                pearsons2.append([abs(float(pearsons[i]))])
            #print 'pearsons2', pearsons2
        mean_correlationo =np.mean(pearsons2)
        #print 'mean_correlationo', mean_correlationo
        mean_correlations[ithperm] = mean_correlationo
        

    bestpermi=max(mean_correlations.iteritems(), key=operator.itemgetter(1))[0]

    print 'bestpermi', bestpermi, 'mean_correlations[bestpermi]', mean_correlations[bestpermi]

    bestpermi_corr = mean_correlations[bestpermi]
                        
    attribution_matrix = list(itertools.permutations((range(original_nr_of_maps))))[bestpermi]

    print 'attribution_matrix', attribution_matrix

    ######
    ###Check if inversion is necessary and save whole EEG into newraw
    ######

    #list that saves the best correlations for the 4 maps
    map_corr_list = []

    newraw=np.zeros((modelmaps.shape))
    
    for mapi in range(len(sortby_maps[:,0])):
        pr, pp=pearsonr( modelmaps[mapi,:], sortby_maps[attribution_matrix[mapi],:])
        '''
        if modelmaps_re.shape == sortby_maps.shape:
            pr, pp=pearsonr( modelmaps_re[mapi,:], sortby_maps[attribution_matrix[mapi],:])
        else:
            pr, pp=pearsonr( sortby_maps[mapi,:], modelmaps_re[attribution_matrix[mapi],:])
        '''

        if pr < 0:
            print 'r=', pr, 'modelmap', mapi, 'reversed'
            newraw[attribution_matrix[mapi],:]=modelmaps[mapi,:]*-1
            map_corr_list.append(abs(pr))
            '''
            if modelmaps_re.shape == sortby_maps.shape:
                newraw[attribution_matrix[mapi],:]=modelmaps[mapi,:]*-1
                map_corr_list.append(abs(pr))
            else:
                newraw[mapi,:]=modelmaps[attribution_matrix[mapi],:]*-1
                map_corr_list.append(abs(pr))
            '''

        else:
            print 'r=', pr,'modelmap', mapi, 'not reversed'
            newraw[attribution_matrix[mapi],:]=modelmaps[mapi,:]
            map_corr_list.append(pr)
            '''
            if modelmaps_re.shape == sortby_maps.shape:
                print 'r=', pr,'modelmap', mapi, 'not reversed'
                newraw[attribution_matrix[mapi],:]=modelmaps[mapi,:]
                map_corr_list.append(pr)
            else:
                newraw[mapi,:]=modelmaps[attribution_matrix[mapi],:]
                map_corr_list.append(pr)
            '''
    '''
    if modelmaps_re.shape != sortby_maps.shape:
        for i in range((len(sortby_maps[:,0])),len(modelmaps_re[:,0])):
            print attribution_matrix[i]
            newraw[i,:]=modelmaps[attribution_matrix[i],:]
            #map_corr_list.append(abs(pr))
    '''
    attributes={}
    attributes['map_corr_list']=map_corr_list
    return newraw, attributes

####--------------------------------------------------------------------------####
####--------------------------------------------------------------------------####
####--------------------------------------------------------------------------####

def get_io_sortmap_for_series(series, iteration, inputfolder, outputfolder, first_input):
    """
    Gets inputs and outputs of modelmaps for the series of modelmap computations specified.
 
    Parameters
    ----------
    series : {Series_1, Series_2, ...}
		Type of modelmap computation, e.g. 'Series_1' (means across runs for each group pt cond)
    iteration : int
		number of iterations for which different seeds are used
    inputfolder:
		folder of input hdf5
    outputfolder:
		folder of output hd5
    first_input:
		name of dataset of input hdf5 for first model map computation (e.g. 'microstate' or 'modelmaps')

    Returns
    ----------
	inputhdf5 : str
		input hdf5 path
	outputhdf5 : str
		output hdf5 path
	modelmap_input : str
		dataset name of input
	modelmap_output : str
		dataset name of output
	computation_version : str
		type of computation to be performed
	stop : bool
		parameter that stops series when it is done

    """

    stop = False
    inputhdf5 = False
    outputhdf5 = False
    modelmap_input= False
    modelmap_output= False
    computation_version= False

    if series == 'Series_1':
        if iteration == 0:
            ######
            ##means across runs for each group pt cond
            ######
            inputhdf5 = os.path.join( inputfolder, 'all_recordings.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_runs.hdf')
            modelmap_input = first_input
            modelmap_output = 'modelmap'
            computation_version ='means across runs for each group pt cond'
        elif iteration == 1:
            ######
            ##means across conds for each group pt
            ######
            inputhdf5 = os.path.join( outputfolder, 'modelmaps_across_runs.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_conds.hdf')
            modelmap_input = 'modelmap'
            modelmap_output = 'modelmap'
            computation_version ='means across conds for each group pt'
        elif iteration == 2:
            ######
            ##means across pts for each group
            ######
            inputhdf5 = os.path.join( outputfolder, 'modelmaps_across_conds.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_pts.hdf')
            modelmap_input = 'modelmap'
            modelmap_output = 'modelmap'
            computation_version ='means across pts for each group'
        elif iteration == 3:
            ######
            ##means across groups
            ######
            inputhdf5 = os.path.join( outputfolder, 'modelmaps_across_pts.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_groups.hdf')
            modelmap_input = 'modelmap'
            modelmap_output = 'modelmap'
            computation_version ='means across groups'
        else:
            stop = True

    elif series == 'Series_2':
        if iteration == 0:
            ######
            ##means across pts for each group cond run
            ######
            inputhdf5 = os.path.join( inputfolder, 'all_recordings.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_pts.hdf')
            modelmap_input = first_input
            modelmap_output = 'modelmap'
            computation_version ='means across pts for each group cond run'
        elif iteration == 1:
            ######
            ##means across runs for each group cond
            ######
            inputhdf5 = os.path.join( outputfolder, 'modelmaps_across_pts.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_runs.hdf')
            modelmap_input = 'modelmap'
            modelmap_output = 'modelmap'
            computation_version ='means across runs for each group cond'
        else:
            stop = True

    elif series == 'Series_3':
        if iteration == 0:
            ######
            ##sort modelmaps_across_groups by modelmaps_across_conds
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_groups.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_conds.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_groups{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdatset = 'sorted_modelmap'
        elif iteration == 1:
            ######
            ##sort modelmaps_across_pts by modelmaps_across_groups_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_pts.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_groups_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_pts{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdatset = 'sorted_modelmap'
            sortdata_provider = SortGroupCondByGroupDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset)
        elif iteration == 2:
            ######
            ##sort modelmaps_across_pts by modelmaps_across_groups_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_pts.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_groups_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_pts{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'sorted_modelmap'
            outputdatset = 'sorted_modelmap'
            sortdata_provider = SortGroupCondByCondDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset)
        elif iteration == 3:
            ######
            ##sort modelmaps_across_runs by modelmaps_across_pts_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_runs.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_pts_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_runs{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'sorted_modelmap'
            outputdatset = 'sorted_modelmap'
            sortdata_provider = SortGroupPtCondByGroupCondDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset)
        else:
            stop = True

    elif series == 'Series_4':
        if iteration == 0:
            ######
            ##means across runs for each group pt cond
            ######
            inputhdf5 = os.path.join( inputfolder, 'all_recordings.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_runs.hdf')
            modelmap_input = first_input
            modelmap_output = 'modelmap'
            computation_version ='means across runs for each group pt cond'
        elif iteration == 1:
            ######
            ##means across conds for each group pt
            ######
            inputhdf5 = os.path.join( outputfolder, 'modelmaps_across_runs.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_conds.hdf')
            modelmap_input = 'modelmap'
            modelmap_output = 'modelmap'
            computation_version ='means across conds for each group pt'
        elif iteration == 2:
            ######
            ##means across pts for each group
            ######
            inputhdf5 = os.path.join( outputfolder, 'modelmaps_across_conds.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_pts.hdf')
            modelmap_input = 'modelmap'
            modelmap_output = 'modelmap'
            computation_version ='means across pts for each group'
        elif iteration == 3:
            ######
            ##means across groups
            ######
            inputhdf5 = os.path.join( outputfolder, 'modelmaps_across_pts.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_groups.hdf')
            modelmap_input = 'modelmap'
            modelmap_output = 'modelmap'
            computation_version ='means across groups'
        else:
            stop = True

    elif series == 'Series_5':
        if iteration == 0:
            ######
            ##means across runs for each group pt cond
            ######
            inputhdf5 = os.path.join( inputfolder, 'all_recordings.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_runs.hdf')
            modelmap_input = first_input
            modelmap_output = 'modelmap'
            computation_version ='means across runs for each group pt cond'
        elif iteration == 1:
            ######
            ##means across conds for each group pt
            ######
            inputhdf5 = os.path.join( outputfolder, 'modelmaps_across_runs.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_conds.hdf')
            modelmap_input = 'modelmap'
            modelmap_output = 'modelmap'
            computation_version ='means across conds for each group pt'
        elif iteration == 2:
            ######
            ##means across groups for each pt
            ######
            inputhdf5 = os.path.join( outputfolder, 'modelmaps_across_conds.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_groups.hdf')
            modelmap_input = 'modelmap'
            modelmap_output = 'modelmap'
            computation_version ='means across groups for each pt'
        elif iteration == 3:
            ######
            ##means across groups
            ######
            inputhdf5 = os.path.join( outputfolder, 'modelmaps_across_groups.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_pts.hdf')
            modelmap_input = 'modelmap'
            modelmap_output = 'modelmap'
            computation_version ='means across groups'
        else:
            stop = True
    else:
        print 'series not defined:', series


    sortdata_provider = SortGroupsByAllDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdatset)

    return sortdata_provider, stop


def run_sort_maps(data_provider, find_model_maps, confobj):
    for output_path in data_provider.get_outputs():
        input = data_provider.get_input_data(output_path)
        sortby = data_provider.get_sortby_data(output_path)
        output_data, output_attributes = sort_maps(confobj, input, sortby)     
        if not output_data == []:
            data_provider.write_output_data(output_path, output_data, output_attributes)


def run_sort_maps_series(series, inputfolder, outputfolder, first_input, confobj):
    """
    Runs run_model_maps for each input of the series.
 
    Parameters
    ----------
    series : {'Series_1', .... 'Series_5'}
        Defines the order of modelmap sorting.
        'Series_1' : (1) , (2)  , (3) , (4) 
        'Series_2' : (1) , (2)  , (3) , (4) 
        'Series_3' : (1) , (2)  , (3) , (4) 
        'Series_4' : (1) , (2)  , (3) , (4) 
        'Series_5' : (1) , (2)  , (3) , (4) 
    inputfolder : str
        path to folder that contains the input hdf5 file
    outputfolder : str
        path to folder that will contain the output hdf5 file
    first_input : str
        dataset name that contains the N microstates for each run that the modelmap computation should be based on, e.g. 'microstate'.
    confobj : object of type MstConfiguration
         Contains the parameters used for microstate computation and visualization. 
    """

    stop = False
    iteration = 0
    while True:
        sortdata_provider, stop = get_io_sortmap_for_series(series, iteration, inputfolder, outputfolder, first_input)
        if stop:
            break

        run_sort_maps(sortdata_provider, find_model_maps, confobj)

        iteration = iteration + 1
