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

###########################
#######  Modelmaps ########
###########################

####Only for internal use for modelmap computation

class ConditionPath(object):
    """
        Attention: This class is for library internal use only. 

        Objects of this type are returned by the Data Provider to specify
        a single output value.
    """
    def __init__(self, group, pt, cond):
        self.group = group
        self.pt = pt
        self.cond = cond

    def __hash__(self):
        return hash((self.group, self.pt, self.cond))

    def __eq__(self, other):
        return self.group == other.group and self.pt == other.pt and self.cond == other.cond

class PtPath(object):
    """
        Attention: This class is for library internal use only. 

        Objects of this type are returned by the Data Provider to specify
        a single output value.
    """
    def __init__(self, group, cond):
        self.group = group
        self.cond = cond

    def __hash__(self):
        return hash((self.group, self.cond))

    def __eq__(self, other):
        return self.group == other.group and self.cond == other.cond

class PtPath2(object):
    """
        Attention: This class is for library internal use only. 

        Objects of this type are returned by the Data Provider to specify
        a single output value.
    """
    def __init__(self, group, pt):
        self.group = group
        self.pt = pt

    def __hash__(self):
        return hash((self.group, self.pt))

    def __eq__(self, other):
        return self.group == other.group and self.pt == other.pt

class GroupPath(object):
    """
        Attention: This class is for library internal use only. 

        Objects of this type are returned by the Data Provider to specify
        a single output value.
    """
    def __init__(self, group):
        self.group = group

    def __hash__(self):
        return hash((self.group))

    def __eq__(self, other):
        return self.group == other.group

class AllPath(object):
    """
        Attention: This class is for library internal use only. 

        Objects of this type are returned by the Data Provider to specify
        a single output value.
    """
    def __init__(self, all):
        self.all = all

    def __hash__(self):
        return hash((self.all))

    def __eq__(self, other):
        return self.all == other.all

class RunPath(object):
    """
        Attention: This class is for library internal use only. 

        Objects of this type are returned by the Data Provider to specify
        a single output value.
    """
    def __init__(self, group, cond, run):
        self.group = group
        self.cond = cond
        self.run = run

    def __hash__(self):
        return hash((self.group, self.cond, self.run))

    def __eq__(self, other):
        return self.group == other.group and self.cond == other.cond and self.run == other.run


class ConditionPath2(object):
    """
        Attention: This class is for library internal use only. 

        Objects of this type are returned by the Data Provider to specify
        a single output value.
    """
    def __init__(self, group, cond):
        self.group = group
        self.cond = cond

    def __hash__(self):
        return hash((self.group, self.cond))

    def __eq__(self, other):
        return self.group == other.group and self.cond == other.cond

class ConditionPath3(object):
    """
        Attention: This class is for library internal use only. 

        Objects of this type are returned by the Data Provider to specify
        a single output value.
    """
    def __init__(self, cond):
        self.cond = cond

    def __hash__(self):
        return hash((self.cond))

    def __eq__(self, other):
        return self.cond == other.cond

####--------------------------------------------------------------------------####


class DataProvider(object):
    """
        An abstract base class for all data provider implementations of the keypy microstates
        package.

        Abstract methods :
            get_outputs(self) :
                A method without parameters. When called the method collects all the output
                elements which are generated through the data provider implementation and
                returns them as a list.

            get_input_data(self, output_path) :
                Expects a single element from the list returned by get_outputs as a parameter.
                The method should be implemented to collect all input elements which are needed
                to call the algorithm for the desired output element. The collected  input data
                is then returned as a list.

            write_output_data(self, output_path, output_data, output_attributes) :
                The parameters to this method are: output_path which is an element from the list
                returned by get_outputs, output_data which is the data generated by the algorithm
                and output_attributes which are additional attributes stored with the output data.
                The method should be implemented to store the data in the data container.

        The data provider is meant to be used for the run_model_maps function. See the documentation
        of this function for an example of how to use data provider instances.
    """

    def __init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset):
        self._file = inputhdf5
        self._outputfile = outputhdf5
        self._inputdataset = inputdataset
        self._outputdataset = outputdataset


####Sub Classes of DataProvider

class CondDataProvider1(DataProvider):
    """
        Used to compute 'means across runs for each group pt cond'. 
    """

    def __init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset):
        DataProvider.__init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_group in f['/'].keys():
                group_group = f['/{0}' .format(current_group)]
                for current_pt in group_group.keys():
                    pt_group = f['/{0}/{1}' .format(current_group, current_pt)]
                    for current_cond in pt_group.keys():
                        out_paths_set.add(ConditionPath(current_group, current_pt, current_cond)) 
        return list(out_paths_set)
    
    #call it once per output, to get a list of modelmaps which are to be "averaged"
    def get_input_data(self, output_path):
        model_maps_all = []
        with closing( h5py.File(self._file, 'r') ) as f:            
            path = '/{0}/{1}/{2}' .format(output_path.group, output_path.pt, output_path.cond)
            for microstate_run in f[path].keys():
                microstate_run_value = f['/{0}/{1}/{2}' .format(path, microstate_run, self._inputdataset)] 
                if all(microstate_run_value[0,:] == 0):
                    print('Error!', path, microstate_run, 'has all zeros', 'group, pt, cond ignored.')    
                else:
                    model_maps_all.append(microstate_run_value[:])
        return model_maps_all
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile) ) as h:
            print('output_paths used for output', output_path.group, output_path.cond, output_path.pt)

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

            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                cond_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in cond_group.keys():
                print('group, participant, condition already in outputfile, not recomputed', group_group, output_path.pt, output_path.cond)
            else:
                cond_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)


class PtDataProvider1(DataProvider):
    """
        Used to compute 'means across conds for each group pt'. 
    """
    def __init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset):
        DataProvider.__init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset)

    def get_outputs(self):
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_group in f['/'].keys():
                group_group = f['/{0}' .format(current_group)]
                for current_pt in group_group.keys():
                    out_paths_set.add(PtPath2(current_group, current_pt))
        return list(out_paths_set)
    
    def get_input_data(self, output_path):
        model_maps_all = []
        with closing( h5py.File(self._file, 'r') ) as f:   
            
            #loop across participants
            for cond in f['/{0}/{1}' .format(output_path.group, output_path.pt)].keys():                
                path = '/{0}/{1}/{2}' .format(output_path.group, output_path.pt, cond)

                if path in f:
                    microstate_run_value = f['/{0}/{1}' .format(path, self._inputdataset)] 
                    if all(microstate_run_value[0,:] == 0):
                        print('Warning!', path, cond, 'has all zeros', 'group, pt, cond ignored.')    
                    else:
                        model_maps_all.append(microstate_run_value[:])
                else:
                    print('Error!', path, cond, 'does not exist', 'group, pt, cond ignored.')    
        return model_maps_all

    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile) ) as h:
            print('output_paths used for output', output_path.group, output_path.pt)

            if output_path.group in h['/'].keys():
                group_group = h['{0}' .format(output_path.group)]
            else:
                group_group = h['/'].create_group( '{0}' .format(output_path.group)  ) 
                      
            if output_path.pt in group_group.keys():
                pt_group = h['/{0}/{1}' .format(output_path.group, output_path.pt)]
            else:
                pt_group = group_group.create_group( '{0}' .format(output_path.pt)  )    

            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                pt_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in pt_group.keys():
                print('group, participant, condition already in outputfile, not recomputed', group_group, output_path.pt)
            else:
                pt_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)
                

class GroupDataProvider1(DataProvider):
    """
        Used to compute 'means across conds for each group'. 
    """

    def __init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset):
        DataProvider.__init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset)

    def get_outputs(self):
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_group in f['/'].keys():
                out_paths_set.add(GroupPath(current_group)) 
        return list(out_paths_set)
    
    def get_input_data(self, output_path):
        model_maps_all = []
        with closing( h5py.File(self._file, 'r') ) as f:            
            path = '/{0}' .format(output_path.group)
            for microstate_cond in f[path].keys():
                if '/{0}/{1}/{2}' .format(path, microstate_cond, self._inputdataset) in f:
                    microstate_cond_value = f['/{0}/{1}/{2}' .format(path, microstate_cond, self._inputdataset)] 
                    if all(microstate_cond_value[0,:] == 0):
                        print('Error!', path, current_run, 'has all zeros', 'group, pt, cond ignored.')    
                    else:
                        model_maps_all.append(microstate_cond_value[:])
                else:
                    print('Error!', path, microstate_cond, 'does not exist', 'group, pt, cond ignored.')    
        return model_maps_all

    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile) ) as h:
            print('output_paths used for output', output_path.group)

            if output_path.group in h['/'].keys():
                group_group = h['{0}' .format(output_path.group)]
            else:
                group_group = h['/'].create_group( '{0}' .format(output_path.group)  ) 
                      
            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                group_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in group_group.keys():
                print('group, participant, condition already in outputfile, not recomputed', group_group)
            else:
                group_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)


class AllDataProvider1(DataProvider):
    """
        Used to compute 'means across groups'. 
    """
    def __init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset):
        DataProvider.__init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset)

    def get_outputs(self):
        out_paths_set = set()
        out_paths_set.add(AllPath('all')) 
        return list(out_paths_set)

    def get_input_data(self, output_path):
        model_maps_all = []
        with closing( h5py.File(self._file, 'r') ) as f:            
            for microstate_group in f.keys():
                if '/{0}/{1}' .format(microstate_group, self._inputdataset) in f:
                    microstate_group_value = f['/{0}/{1}' .format(microstate_group, self._inputdataset)] 
                    if all(microstate_group_value[0,:] == 0):
                        print('Error!', microstate_group, 'has all zeros', 'group, pt, cond ignored.')    
                    else:
                        model_maps_all.append(microstate_group_value[:])
                else:
                    print('Error!', microstate_group, 'does not exist', 'group, pt, cond ignored.')    
        return model_maps_all

    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile) ) as h:
            print('output_paths used for output', output_path.all)

            if output_path.all in h['/'].keys():
                group_group = h['{0}' .format(output_path.all)]
            else:
                group_group = h['/'].create_group( '{0}' .format(output_path.all)  ) 
                      
            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                group_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in group_group.keys():
                print('group already in outputfile, not recomputed', group_group)
            else:
                group_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)


class RunDataProvider1(DataProvider):
    """
        Used to compute 'means across pts for each group cond run'. 
    """
    def __init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset):
        DataProvider.__init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset)

    def get_outputs(self):      
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_group in f['/'].keys():
                group_group = f['/{0}' .format(current_group)]
                for current_pt in group_group.keys():
                    pt_group = f['/{0}/{1}' .format(current_group, current_pt)]
                    for current_cond in pt_group.keys():
                        cond_group = f['/{0}/{1}/{2}' .format(current_group, current_pt, current_cond)]
                        for current_run in cond_group.keys():
                            out_paths_set.add(RunPath(current_group, current_cond, current_run)) 
        return list(out_paths_set)
    
    def get_input_data(self, output_path):
        model_maps_all = []
        with closing( h5py.File(self._file, 'r') ) as f:   
            
            #loop across participants
            for pt in f['/{0}' .format(output_path.group)].keys():                
                path = '/{0}/{1}/{2}/{3}' .format(output_path.group, pt, output_path.cond, output_path.run)

                if path in f:
                    microstate_run_value = f['/{0}/{1}' .format(path, self._inputdataset)] 
                    if all(microstate_run_value[0,:] == 0):
                        print('Warning!', path, pt, 'has all zeros', 'group, pt, cond ignored.')    
                    else:
                        model_maps_all.append(microstate_run_value[:])
                else:
                    print('Error!', path, pt, 'does not exist', 'group, pt, cond ignored.')    
        return model_maps_all

    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile) ) as h:
            print('output_paths used for output', output_path.group, output_path.cond)

            if output_path.group in h['/'].keys():
                group_group = h['{0}' .format(output_path.group)]
            else:
                group_group = h['/'].create_group( '{0}' .format(output_path.group)  ) 

            if output_path.cond in group_group.keys():
                cond_group = h['/{0}/{1}' .format(output_path.group, output_path.cond)]
            else:
                cond_group = group_group.create_group( '{0}' .format(output_path.cond)  )   

            if output_path.run in cond_group.keys():
                run_group = h['/{0}/{1}/{2}' .format(output_path.group, output_path.cond, output_path.run)]
            else:
                run_group = cond_group.create_group( '{0}' .format(output_path.run)  )   

            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                run_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in run_group.keys():
                print('group, participant, condition already in outputfile, not recomputed', run_group, output_path.run)
            else:
                run_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)


class CondDataProvider5(DataProvider):
    """
        Used to compute 'means across runs for each group cond'. 
    """
    def __init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset):
        DataProvider.__init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset)

    def get_outputs(self):
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_group in f['/'].keys():
                group_group = f['/{0}' .format(current_group)]
                for current_cond in group_group.keys():
                    out_paths_set.add(ConditionPath2(current_group, current_cond)) 
        return list(out_paths_set)
    
    def get_input_data(self, output_path):
        model_maps_all = []
        with closing( h5py.File(self._file, 'r') ) as f:            
            path = '/{0}/{1}' .format(output_path.group, output_path.cond)
            for microstate_run in f[path].keys():
                if '/{0}/{1}/{2}' .format(path, microstate_run, self._inputdataset) in f:
                    microstate_run_value = f['/{0}/{1}/{2}' .format(path, microstate_run, self._inputdataset)] 
                    if all(microstate_run_value[0,:] == 0):
                        print('Error!', path, current_run, 'has all zeros', 'group, pt, cond ignored.')    
                    else:
                        model_maps_all.append(microstate_run_value[:])
                else:
                    print('key', microstate_run, 'in', f[path].keys(), 'no run group, not considered')
        return model_maps_all

    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile) ) as h:
            print('output_paths used for output', output_path.group, output_path.cond)

            if output_path.group in h['/'].keys():
                group_group = h['{0}' .format(output_path.group)]
            else:
                group_group = h['/'].create_group( '{0}' .format(output_path.group)  ) 

            if output_path.cond in group_group.keys():
                cond_group = h['/{0}/{1}' .format(output_path.group, output_path.cond)]
            else:
                cond_group = group_group.create_group( '{0}' .format(output_path.cond)  )   

            for key, value in output_attributes.iteritems():
                cond_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in cond_group.keys():
                print('group, participant, condition already in outputfile, not recomputed', group_group, output_path.cond)
            else:
                cond_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)



class CondDataProvider2(DataProvider):
    """
        Used to compute 'means across pts for each group cond'. 
    """
    def __init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset):
        DataProvider.__init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset)

    def get_outputs(self):
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_group in f['/'].keys():
                group_group = f['/{0}' .format(current_group)]
                for current_pt in group_group.keys():
                    pt_group = f['/{0}/{1}' .format(current_group, current_pt)]
                    for current_cond in pt_group.keys():
                       out_paths_set.add(ConditionPath2(current_group, current_cond))
        return list(out_paths_set)
    
    def get_input_data(self, output_path):
        model_maps_all = []
        with closing( h5py.File(self._file, 'r') ) as f:            
            #loop across participants
            for pt in f['/{0}' .format(output_path.group)].keys():                
                path = '/{0}/{1}/{2}' .format(output_path.group, pt, output_path.cond)

                if path in f:
                    microstate_run_value = f['/{0}/{1}' .format(path, self._inputdataset)] 
                    if all(microstate_run_value[0,:] == 0):
                        print('Warning!', path, pt, 'has all zeros', 'group, pt, cond ignored.')    
                    else:
                        model_maps_all.append(microstate_run_value[:])
                else:
                    print('Error!', path, pt, 'does not exist', 'group, pt, cond ignored.')    
        return model_maps_all

    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile) ) as h:
            print('output_paths used for output', output_path.group, output_path.cond)

            if output_path.group in h['/'].keys():
                group_group = h['{0}' .format(output_path.group)]
            else:
                group_group = h['/'].create_group( '{0}' .format(output_path.group)  ) 

            if output_path.cond in group_group.keys():
                cond_group = h['/{0}/{1}' .format(output_path.group, output_path.cond)]
            else:
                cond_group = group_group.create_group( '{0}' .format(output_path.cond)  )   

            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                cond_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in cond_group.keys():
                print('group, condition already in outputfile, not recomputed', group_group, output_path.cond)
            else:
                cond_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)


class CondDataProvider3(DataProvider):
    """
        Used to compute 'means across groups for each cond'. 
    """
    def __init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset):
        DataProvider.__init__(self, inputhdf5, outputhdf5, inputdataset, outputdataset)

    def get_outputs(self):
        out_paths_set = set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_group in f['/'].keys():
                group_group = f['/{0}' .format(current_group)]
                for current_cond in group_group.keys():
                    out_paths_set.add(ConditionPath3(current_cond)) 
        return list(out_paths_set)
    
    def get_input_data(self, output_path):
        model_maps_all = []
        with closing( h5py.File(self._file, 'r') ) as f:            
            #loop across groups
            for group in f['/'].keys():                
                path = '/{0}/{1}' .format(group, output_path.cond)

                if path in f:
                    microstate_run_value = f['/{0}/{1}' .format(path, self._inputdataset)] 
                    if all(microstate_run_value[0,:] == 0):
                        print('Warning!', path, pt, 'has all zeros', 'group, pt, cond ignored.')    
                    else:
                        model_maps_all.append(microstate_run_value[:])
                else:
                    print('Error!', path, group, 'does not exist', 'group, pt, cond ignored.')    

        return model_maps_all

    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile) ) as h:
            print('output_paths used for output', output_path.cond)

            if output_path.cond in h['/'].keys():
                cond_group = h['{0}' .format(output_path.cond)]
            else:
                cond_group = h['/'].create_group( '{0}' .format(output_path.cond)  ) 

            #Save best mean correlation as attribute to group and modelmaps as dataset
            for key, value in output_attributes.iteritems():
                cond_group.attrs['{0}' .format(key)] = value

            if self._outputdataset in cond_group.keys():
                print('group, condition already in outputfile, not recomputed', output_path.cond)
            else:
                cond_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)


####--------------------------------------------------------------------------####


##################################
#######  find_model_maps  ########
##################################

def get_data_provider_class(computation_version):
    """
    Determines the correct data_provider class for the desired output from the desired input.

    Parameters
    ----------
    computation_version : str
        A string that specifies what is being computed based on which input (e.g. 'means across runs for each group pt cond').

    Returns
    -------
    data_provider : class
        The correct class necessary for the input output mean computation.

    """
    if computation_version =='means across runs for each group pt cond':
        data_provider = CondDataProvider1
    elif computation_version =='means across pts for each group cond run':
        data_provider = RunDataProvider1
    elif computation_version =='means across pts for each group cond':    
        data_provider = CondDataProvider2
    elif computation_version =='means across conds for each group pt':
        data_provider = PtDataProvider1           
    elif computation_version =='means across conds for each group':
        data_provider = GroupDataProvider1
    elif computation_version =='means across groups for each cond':  
        data_provider = CondDataProvider3
    elif computation_version =='means across pts for each group':  
        data_provider = GroupDataProvider1
    elif computation_version =='means across groups for each pt':  
        data_provider = CondDataProvider3
    elif computation_version =='means across groups':  
        data_provider = AllDataProvider1
    elif computation_version =='means across conds':  
        data_provider = AllDataProvider1
    elif computation_version =='means across runs for each group cond':  
        data_provider = CondDataProvider5
    else:
        print(computation_version, 'not implemented')

    return data_provider
