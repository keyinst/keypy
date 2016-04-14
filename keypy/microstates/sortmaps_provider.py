# -*- coding: utf-8 -*-

##################################
#######  Import Packages  ########
##################################

from __future__ import print_function

from contextlib import closing
import h5py
import numpy as np
from keypy.microstates.microstates_helper import *


##########################
#######  Classes  ########
##########################

#######################################
########  Sort Maps  Providers ########
#######################################

##Only for internal use for map sorting and for internal use by parameters.py

###Create Paths for sortmaps

class Levels1Path(object):
    def __init__(self, level0):
        self.level0 = level0

    def __hash__(self):
        return hash((self.level0))

    def __eq__(self, other):
        return self.level0 == other.level0

class Levels2Path(object):
    def __init__(self, level0, level1):
        self.level0 = level0
        self.level1 = level1

    def __hash__(self):
        return hash((self.level0, self.level1))

    def __eq__(self, other):
        return self.level0 == other.level0, self.level1 == other.level1

class Levels3Path(object):
    def __init__(self, level0, level1, level2):
        self.level0 = level0
        self.level1 = level1
        self.level2 = level2

    def __hash__(self):
        return hash((self.level0, self.level1, self.level2))

    def __eq__(self, other):
        return self.level0 == other.level0, self.level1 == other.level1, self.level2 == other.level2

class Levels4Path(object):
    def __init__(self, level0, level1, level2, level3):
        self.level0 = level0
        self.level1 = level1
        self.level2 = level2
        self.level3 = level3

    def __hash__(self):
        return hash((self.level0, self.level1, self.level2, self.level3))

    def __eq__(self, other):
        return self.level0 == other.level0, self.level1 == other.level1, self.level2 == other.level2, self.level3 == other.level3

###-----------------------


###########################
#######  Sortmaps  ########
###########################

##Only for internal use for map sorting

####Abstract Class SortDataProvider

class SortDataProvider(object):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist):
        self._file = inputhdf5
        self._sortbyfile = sortbyhdf5
        self._outputfile = outputhdf5
        self._inputdataset = inputdataset
        self._sortbydataset = sortbydataset
        self._outputdataset = outputdataset
        self._sortbychlist = sortbyfile_chlist

####Sub Class CondDataProvider

###Sort All by Norm Data Provider 1

class SortAllByNormDataProvider1(SortDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist):
        SortDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_level0 in f['/'].keys():
                out_paths_set.add(Levels1Path(current_level0))                     
        return list(out_paths_set)
    
    ### group pt cond run
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted
    def get_input_data(self, output_path, own_chlist):
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '{0}' .format(output_path.level0)
            microstate_run_value = g['/{0}/{1}' .format(path, self._inputdataset)] 
            if all(microstate_run_value[0,:] == 0):
                print('Error!', path, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]

            #reduce channels of input model_map to match sortby modelmap PPPP
            sortbyfolder = os.path.dirname(self._sortbyfile)
            sortbychlist_path = os.path.join(sortbyfolder,self._sortbychlist)
            model_map_new, _ = reduce_channels(model_map, self._sortbyfile, own_chlist, sortbychlist_path)

        return model_map_new, model_map

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        #old
        #microstate_run_value = np.loadtxt(self._sortbyfile)
        #loads directly reduced external file as computed by reduce_channels function in microstates_helper.py
        microstate_run_value = np.loadtxt(os.path.join(os.path.dirname(self._sortbyfile),"{0}_reduced.asc".format(os.path.splitext(os.path.basename(self._sortbyfile))[0])))

        if all(microstate_run_value[0,:] == 0):
            print('Error!', sortbyhdf5, 'has all zeros', 'group, pt, cond ignored.')    
        else:
            model_map=microstate_run_value[:]

        return model_map
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile, 'a') ) as k:
            print('output_paths used for output', output_path.level0)

            if output_path.level0 in k['/'].keys():
                all_group = k['{0}' .format(output_path.level0)]
            else:
                all_group = k['/'].create_group( '{0}' .format(output_path.level0)  ) 

            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                all_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in all_group.keys():
                print('group, participant, condition already in outputfile, not recomputed', all_group, output_path.level0)
            else:
                all_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)


                
###Sort Group by All Data Provider 1

class SortGroupByAllDataProvider1(SortDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist):
        SortDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_level0 in f['/'].keys():                 
                for current_level1 in f['/{0}'.format(current_level0)].keys():  
                    out_paths_set.add(Levels1Path(current_level0))
        return list(out_paths_set)
    
    ### group pt cond run
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted
    def get_input_data(self, output_path, own_chlist):
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '{0}' .format(output_path.level0)
            microstate_run_value = g['/{0}/{1}' .format(path, self._inputdataset)] 
            if all(microstate_run_value[0,:] == 0):
                print('Error!', path, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]
        return model_map, model_map

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        with closing( h5py.File(self._sortbyfile, 'r') ) as h:
            microstate_run_value = h['/{0}/{1}' .format('all', self._sortbydataset)] 

            if all(microstate_run_value[0,:] == 0):
                print('Error!', sortbyhdf5, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]
        return model_map
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile, 'a') ) as k:
            print('output_paths used for output', output_path.level0)

            if output_path.level0 in k['/'].keys():
                group_group = k['{0}' .format(output_path.level0)]
            else:
                group_group = k['/'].create_group( '{0}' .format(output_path.level0)  ) 

            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                group_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in group_group.keys():
                print('group, participant, condition already in outputfile, not recomputed', group_group, output_path.level0)
            else:
                group_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)

###Sort Pt by Group Data Provider 1

class SortPtByGroupDataProvider1(SortDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist):
        SortDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_level0 in f['/'].keys():                 
                for current_level1 in f['/{0}'.format(current_level0)].keys():  
                    out_paths_set.add(Levels2Path(current_level0, current_level1))
        return list(out_paths_set)
    
    ### group pt cond run
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted
    def get_input_data(self, output_path, own_chlist):
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '{0}/{1}' .format(output_path.level0, output_path.level1)
            microstate_run_value = g['/{0}/{1}' .format(path, self._inputdataset)] 
            if all(microstate_run_value[0,:] == 0):
                print('Error!', path, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]
        return model_map, model_map

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        with closing( h5py.File(self._sortbyfile, 'r') ) as h:
            microstate_run_value = h['/{0}/{1}' .format(output_path.level0, self._sortbydataset)] 

            if all(microstate_run_value[0,:] == 0):
                print('Error!', sortbyhdf5, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]
        return model_map
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile, 'a') ) as k:
            print('output_paths used for output', output_path.level0)

            if output_path.level0 in k['/'].keys():
                group_group = k['{0}' .format(output_path.level0)]
            else:
                group_group = k['/'].create_group( '{0}' .format(output_path.level0)  ) 

            if output_path.level1 in group_group.keys():
                pt_group = k['/{0}/{1}' .format(output_path.level0, output_path.level1)]
            else:
                pt_group = group_group.create_group( '{0}' .format(output_path.level1)  )  



            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                pt_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in pt_group.keys():
                print('group, participant, condition already in outputfile, not recomputed', pt_group, output_path.level0, output_path.level1)
            else:
                pt_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)


###Sort Cond by Pt Data Provider 1

class SortCondByPtDataProvider1(SortDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist):
        SortDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_level0 in f['/'].keys():                 
                for current_level1 in f['/{0}'.format(current_level0)].keys():  
                    for current_level2 in f['/{0}/{1}'.format(current_level0, current_level1)].keys():  
                        out_paths_set.add(Levels3Path(current_level0, current_level1, current_level2))
        return list(out_paths_set)
    
    ### group pt cond run
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted
    def get_input_data(self, output_path, own_chlist):
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '{0}/{1}/{2}' .format(output_path.level0, output_path.level1, output_path.level2)
            microstate_run_value = g['/{0}/{1}' .format(path, self._inputdataset)] 
            if all(microstate_run_value[0,:] == 0):
                print('Error!', path, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]
        return model_map, model_map

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        with closing( h5py.File(self._sortbyfile, 'r') ) as h:
            microstate_run_value = h['/{0}/{1}/{2}' .format(output_path.level0, output_path.level1, self._sortbydataset)] 

            if all(microstate_run_value[0,:] == 0):
                print('Error!', sortbyhdf5, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]
        return model_map
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile, 'a') ) as k:
            print('output_paths used for output', output_path.level0)

            if output_path.level0 in k['/'].keys():
                group_group = k['{0}' .format(output_path.level0)]
            else:
                group_group = k['/'].create_group( '{0}' .format(output_path.level0)  ) 

            if output_path.level1 in group_group.keys():
                pt_group = k['/{0}/{1}' .format(output_path.level0, output_path.level1)]
            else:
                pt_group = group_group.create_group( '{0}' .format(output_path.level1)  )  

            if output_path.level2 in pt_group.keys():
                cond_group = k['/{0}/{1}/{2}' .format(output_path.level0, output_path.level1, output_path.level2)]
            else:
                cond_group = pt_group.create_group( '{0}' .format(output_path.level2)  )  

            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                cond_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in cond_group.keys():
                print('group, participant, condition already in outputfile, not recomputed', pt_group, output_path.level0, output_path.level1, output_path.level2)
            else:
                cond_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)




###Sort Run by Cond Data Provider 1

class SortRunByCondDataProvider1(SortDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist):
        SortDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_level0 in f['/'].keys():                 
                for current_level1 in f['/{0}'.format(current_level0)].keys():  
                    for current_level2 in f['/{0}/{1}'.format(current_level0, current_level1)].keys():
                        for current_level3 in f['/{0}/{1}/{2}'.format(current_level0, current_level1, current_level2)].keys():
                            out_paths_set.add(Levels4Path(current_level0, current_level1, current_level2, current_level3))                        
        return list(out_paths_set)
    
    ### group pt cond run
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted
    def get_input_data(self, output_path, own_chlist):
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '{0}/{1}/{2}/{3}' .format(output_path.level0, output_path.level1, output_path.level2, output_path.level3)
            microstate_run_value = g['/{0}/{1}' .format(path, self._inputdataset)] 
            if all(microstate_run_value[0,:] == 0):
                print('Error!', path, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]
        return model_map, model_map

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        with closing( h5py.File(self._sortbyfile, 'r') ) as h:
            microstate_run_value = h['/{0}/{1}/{2}/{3}' .format(output_path.level0, output_path.level1, output_path.level2, self._sortbydataset)] 

            if all(microstate_run_value[0,:] == 0):
                print('Error!', sortbyhdf5, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]
        return model_map
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile, 'a') ) as k:
            print('output_paths used for output', output_path.level0)

            if output_path.level0 in k['/'].keys():
                group_group = k['{0}' .format(output_path.level0)]
            else:
                group_group = k['/'].create_group( '{0}' .format(output_path.level0)  ) 

            if output_path.level1 in group_group.keys():
                pt_group = k['/{0}/{1}' .format(output_path.level0, output_path.level1)]
            else:
                pt_group = group_group.create_group( '{0}' .format(output_path.level1)  )  

            if output_path.level2 in pt_group.keys():
                cond_group = k['/{0}/{1}/{2}' .format(output_path.level0, output_path.level1, output_path.level2)]
            else:
                cond_group = pt_group.create_group( '{0}' .format(output_path.level2)  )  

            if output_path.level3 in cond_group.keys():
                run_group = k['/{0}/{1}/{2}/{3}' .format(output_path.level0, output_path.level1, output_path.level2, output_path.level3)]
            else:
                run_group = cond_group.create_group( '{0}' .format(output_path.level3)  )  

            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                run_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in run_group.keys():
                print('group, participant, condition already in outputfile, not recomputed', pt_group, output_path.level0, output_path.level1, output_path.level2, output_path.level3)
            else:
                run_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)

###Sort Pt by Group Data Provider 2

class SortPtByGroupDataProvider2(SortDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist):
        SortDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_level0 in f['/'].keys():                 
                for current_level1 in f['/{0}'.format(current_level0)].keys():  
                    out_paths_set.add(Levels2Path(current_level0, current_level1))
        return list(out_paths_set)
    
    ### group pt cond run
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted
    def get_input_data(self, output_path, own_chlist):
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '{0}/{1}' .format(output_path.level0, output_path.level1)
            microstate_run_value = g['/{0}/{1}' .format(path, self._inputdataset)] 
            if all(microstate_run_value[0,:] == 0):
                print('Error!', path, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]
        return model_map, model_map

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        with closing( h5py.File(self._sortbyfile, 'r') ) as h:
            microstate_run_value = h['/{0}/{1}' .format(output_path.level0, self._sortbydataset)] 

            if all(microstate_run_value[0,:] == 0):
                print('Error!', sortbyhdf5, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]
        return model_map
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile, 'a') ) as k:
            print('output_paths used for output', output_path.level0)

            if output_path.level0 in k['/'].keys():
                group_group = k['{0}' .format(output_path.level0)]
            else:
                group_group = k['/'].create_group( '{0}' .format(output_path.level0)  ) 

            if output_path.level1 in group_group.keys():
                pt_group = k['/{0}/{1}' .format(output_path.level0, output_path.level1)]
            else:
                pt_group = group_group.create_group( '{0}' .format(output_path.level1)  )  



            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                pt_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in pt_group.keys():
                print('group, participant, condition already in outputfile, not recomputed', pt_group, output_path.level0, output_path.level1)
            else:
                pt_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)


###Sort Cond by Pt Data Provider 2

class SortCondByPtDataProvider2(SortDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist):
        SortDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_level0 in f['/'].keys():                 
                for current_level1 in f['/{0}'.format(current_level0)].keys():  
                    for current_level2 in f['/{0}/{1}'.format(current_level0, current_level1)].keys():  
                        out_paths_set.add(Levels3Path(current_level0, current_level1, current_level2))
        return list(out_paths_set)
    
    ### group pt cond run
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted
    def get_input_data(self, output_path, own_chlist):
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '{0}/{1}/{2}' .format(output_path.level0, output_path.level1, output_path.level2)
            microstate_run_value = g['/{0}/{1}' .format(path, self._inputdataset)] 
            if all(microstate_run_value[0,:] == 0):
                print('Error!', path, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]
        return model_map, model_map

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        with closing( h5py.File(self._sortbyfile, 'r') ) as h:
            microstate_run_value = h['/{0}/{1}/{2}' .format(output_path.level0, output_path.level2, self._sortbydataset)] 

            if all(microstate_run_value[0,:] == 0):
                print('Error!', sortbyhdf5, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]
        return model_map
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile, 'a') ) as k:
            print('output_paths used for output', output_path.level0)

            if output_path.level0 in k['/'].keys():
                group_group = k['{0}' .format(output_path.level0)]
            else:
                group_group = k['/'].create_group( '{0}' .format(output_path.level0)  ) 

            if output_path.level1 in group_group.keys():
                pt_group = k['/{0}/{1}' .format(output_path.level0, output_path.level1)]
            else:
                pt_group = group_group.create_group( '{0}' .format(output_path.level1)  )  

            if output_path.level2 in pt_group.keys():
                cond_group = k['/{0}/{1}/{2}' .format(output_path.level0, output_path.level1, output_path.level2)]
            else:
                cond_group = pt_group.create_group( '{0}' .format(output_path.level2)  )  

            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                cond_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in cond_group.keys():
                print('group, participant, condition already in outputfile, not recomputed', pt_group, output_path.level0, output_path.level1, output_path.level2)
            else:
                cond_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)

###Sort All by Norm Data Provider 2

class SortAllByNormDataProvider2(SortDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist):
        SortDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_level0 in f['/'].keys():
                for current_level1 in f['/{0}'.format(current_level0)].keys():  
                    out_paths_set.add(Levels2Path(current_level0, current_level1))                 
        return list(out_paths_set)
     
    ### group pt cond run
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted
    def get_input_data(self, output_path, own_chlist):
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '{0}/{1}' .format(output_path.level0, output_path.level1)
            microstate_run_value = g['/{0}/{1}' .format(path, self._inputdataset)] 
            if all(microstate_run_value[0,:] == 0):
                print('Error!', path, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]

            #reduce channels of input model_map to match sortby modelmap
            sortbyfolder = os.path.dirname(self._sortbyfile)
            sortbychlist_path = os.path.join(sortbyfolder,self._sortbychlist)
            model_map_new, _ = reduce_channels(model_map, self._sortbyfile, own_chlist, sortbychlist_path)

        return model_map_new, model_map


    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        microstate_run_value = np.loadtxt(os.path.join(os.path.dirname(self._sortbyfile),"{0}_reduced.asc".format(os.path.splitext(os.path.basename(self._sortbyfile))[0])))

        if all(microstate_run_value[0,:] == 0):
            print('Error!', sortbyhdf5, 'has all zeros', 'group, pt, cond ignored.')    
        else:
            model_map=microstate_run_value[:]
        return model_map
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile, 'a') ) as k:
            print('output_paths used for output', output_path.level0)

            if output_path.level0 in k['/'].keys():
                group_group = k['{0}' .format(output_path.level0)]
            else:
                group_group = k['/'].create_group( '{0}' .format(output_path.level0)  ) 

            if output_path.level1 in group_group.keys():
                pt_group = k['/{0}/{1}' .format(output_path.level0, output_path.level1)]
            else:
                pt_group = group_group.create_group( '{0}' .format(output_path.level1)  )  



            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                pt_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in pt_group.keys():
                print('group, participant, condition already in outputfile, not recomputed', pt_group, output_path.level0, output_path.level1)
            else:
                pt_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)

###Sort Run by Cond Data Provider 1

class SortRunByCondDataProvider2(SortDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist):
        SortDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_level0 in f['/'].keys():                 
                for current_level1 in f['/{0}'.format(current_level0)].keys():  
                    for current_level2 in f['/{0}/{1}'.format(current_level0, current_level1)].keys():
                        for current_level3 in f['/{0}/{1}/{2}'.format(current_level0, current_level1, current_level2)].keys():
                            out_paths_set.add(Levels4Path(current_level0, current_level1, current_level2, current_level3))                        
        return list(out_paths_set)
    
    ### group pt cond run
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted
    def get_input_data(self, output_path, own_chlist):
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '{0}/{1}/{2}/{3}' .format(output_path.level0, output_path.level1, output_path.level2, output_path.level3)
            microstate_run_value = g['/{0}/{1}' .format(path, self._inputdataset)] 
            if all(microstate_run_value[0,:] == 0):
                print('Error!', path, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]
        return model_map, model_map

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        with closing( h5py.File(self._sortbyfile, 'r') ) as h:
            microstate_run_value = h['/{0}/{1}/{2}/{3}' .format(output_path.level0, output_path.level2, output_path.level3, self._sortbydataset)] 

            if all(microstate_run_value[0,:] == 0):
                print('Error!', sortbyhdf5, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]
        return model_map
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile, 'a') ) as k:
            print('output_paths used for output', output_path.level0)

            if output_path.level0 in k['/'].keys():
                group_group = k['{0}' .format(output_path.level0)]
            else:
                group_group = k['/'].create_group( '{0}' .format(output_path.level0)  ) 

            if output_path.level1 in group_group.keys():
                pt_group = k['/{0}/{1}' .format(output_path.level0, output_path.level1)]
            else:
                pt_group = group_group.create_group( '{0}' .format(output_path.level1)  )  

            if output_path.level2 in pt_group.keys():
                cond_group = k['/{0}/{1}/{2}' .format(output_path.level0, output_path.level1, output_path.level2)]
            else:
                cond_group = pt_group.create_group( '{0}' .format(output_path.level2)  )  

            if output_path.level3 in cond_group.keys():
                run_group = k['/{0}/{1}/{2}/{3}' .format(output_path.level0, output_path.level1, output_path.level2, output_path.level3)]
            else:
                run_group = cond_group.create_group( '{0}' .format(output_path.level3)  )  

            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                run_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in run_group.keys():
                print('group, participant, condition already in outputfile, not recomputed', pt_group, output_path.level0, output_path.level1, output_path.level2, output_path.level3)
            else:
                run_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)

###Sort Pt by Group Data Provider 1

class SortPtByGroupDataProvider2(SortDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist):
        SortDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_level0 in f['/'].keys():                 
                for current_level1 in f['/{0}'.format(current_level0)].keys():  
                    out_paths_set.add(Levels2Path(current_level0, current_level1))
        return list(out_paths_set)
    
    ### group pt cond run
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted
    def get_input_data(self, output_path, own_chlist):
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '{0}/{1}' .format(output_path.level0, output_path.level1)
            microstate_run_value = g['/{0}/{1}' .format(path, self._inputdataset)] 
            if all(microstate_run_value[0,:] == 0):
                print('Error!', path, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]
        return model_map, model_map

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        with closing( h5py.File(self._sortbyfile, 'r') ) as h:
            microstate_run_value = h['/{0}/{1}' .format(output_path.level1, self._sortbydataset)] 

            if all(microstate_run_value[0,:] == 0):
                print('Error!', sortbyhdf5, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                model_map=microstate_run_value[:]
        return model_map
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile, 'a') ) as k:
            print('output_paths used for output', output_path.level0)

            if output_path.level0 in k['/'].keys():
                group_group = k['{0}' .format(output_path.level0)]
            else:
                group_group = k['/'].create_group( '{0}' .format(output_path.level0)  ) 

            if output_path.level1 in group_group.keys():
                pt_group = k['/{0}/{1}' .format(output_path.level0, output_path.level1)]
            else:
                pt_group = group_group.create_group( '{0}' .format(output_path.level1)  )  



            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                pt_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in pt_group.keys():
                print('group, participant, condition already in outputfile, not recomputed', pt_group, output_path.level0, output_path.level1)
            else:
                pt_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)
