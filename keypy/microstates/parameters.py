##################################
#######  run_mstate_paramters  ########
##################################

from __future__ import print_function

from scipy.stats import nanmean, pearsonr
from keypy.microstates.data_provider import *
from keypy.microstates.microstates_helper import *
from numpy import sqrt
import os.path as op
import os
from sets import Set
import numpy as np


##########################
#######  Functions  ########
##########################

def dissim(map, tf, method_GFPpeak):
    diff1=map-tf
    diff2=-map-tf
    diff1=np.reshape(diff1, (-1, 1))
    diff2=np.reshape(diff2, (-1, 1))
    diff1_gfp=compute_gfp(diff1.T, method_GFPpeak)
    diff2_gfp=compute_gfp(diff2.T, method_GFPpeak)
    diff=min(diff1_gfp,diff2_gfp)

    return diff


def reduce_channels(eeg, own_chlist, external_chlist_path):

    external_chlist = np.genfromtxt(external_chlist_path,dtype='str')

    #if external_chlist and own_chlist are equal no mapping is necessary
    if not external_chlist.tolist()==own_chlist:
        #assert that each ele in external_chlist is an ele in own_chlist also
        for ele in external_chlist:
            if not ele in own_chlist:
                raise AssertionError('Each channel name of your external list needs to also be present in your own chlist.', ele,'is missing')   

        ######
        ###Reduce my Maps to the XX channels of the external channel list
        ######

        #loop to get mapping of own channels to external channels
        chlist_external_index = []
        for ele in external_chlist:
            chlist_external_index.append(own_chlist.index(ele))

        #own maps, TK order, only 19 channels
        eeg_new=eeg[:,chlist_external_index]

        #np.savetxt('data_reduced.txt', eeg_new, delimiter = '\t')

    else:
        eeg_new = eeg

    return eeg_new


##########################
#######  Classes  ########
##########################


####Abstract Class SortDataProvider

class ParametersDataProvider(object):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist):
        self._file = inputhdf5
        self._sortbyfile = sortbyhdf5
        self._outputfile = outputhdf5
        self._inputdataset = inputdataset
        self._sortbydataset = sortbydataset
        self._sortbychlist = sortbychlist

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
    #ein Aufruf pro Output, gets the particular EEG file which is needed for parameter computation (for a particular group, pt, cond, run)
    def get_input_data(self, output_path, own_chlist):
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '/{0}/{1}/{2}/{3}' .format(output_path.level0, output_path.level1, output_path.level2, output_path.level3)
            eeg_value = g['/{0}/{1}' .format(path, self._inputdataset)] 
            if all(eeg_value[0,:] == 0):
                print('Error!', path, 'has all zeros', 'group, pt, cond ignored.')    
            else:
                eeg=eeg_value[:]

                if self._sortbychlist:
                    eeg=reduce_channels(eeg, own_chlist, self._sortbychlist)
        return eeg
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, confobj, output_path, output_data, output_attributes):

        outputfolder=op.join(op.dirname(op.abspath(self._outputfile)))
        if not op.exists(outputfolder):
            print("Create output folder: {0}".format(outputfolder))
            os.makedirs(outputfolder)

        with closing( h5py.File(self._outputfile) ) as k:
            print('output_paths used for output: {0}'.format(output_path.level0))

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

            ############
            ###writing Outputs (is the same for all, just that it is not always based on the group run_group)
            ############
            '''
            if self._outputdataset in run_group.keys():
                print 'group, participant, condition already in outputfile, not recomputed', pt_group, output_path.level0, output_path.level1, output_path.level2, output_path.level3
            else:
                run_group.create_dataset('{0}' .format(self._outputdataset), data = output_data)
            '''


            ###Writing Output for output_attributes and data

            runwise_data, epochwise_data, mapwise_data = output_data
            number_of_epochs = len(mapwise_data['number of ms for each state'])

            ###Add runwise datasets

            ##Add State Match Mean percentage

            #create dataset of mean percentages for class correspondances across epochs
            map_avg_perc=parameter_preprocessing(confobj, epochwise_data['State Match Mean percentage'])  

            if 'State Match Mean p' in run_group.keys():
                print('State Match Mean p already exists. Not recomputed for {0} {1} {2} {3}'.format(output_path.level0, output_path.level1, output_path.level2, output_path.level3))
            else:
                if map_avg_perc.any():
                    run_group.create_dataset('State Match Mean p', data = map_avg_perc)    
                    run_group['State Match Mean p'].attrs['State Match Mean p based on']='Mean {0}' .format(confobj.similarity_measure)   
                else:
                    ### add code 999 for missing value
                    run_group.create_dataset('State Match Mean p', data = 999)
                    run_group['State Match Mean p'].attrs['State Match Mean p based on']='Mean {0}' .format(confobj.similarity_measure)   


            ##add other runwise datasets
            for dataset_name in runwise_data:
                if dataset_name in run_group.keys():
                    print( '{0} already exists. Not recomputed for {0} {1} {2} {3}'.format(dataset_name, output_path.level0, output_path.level1, output_path.level2, output_path.level3))
                else:
                    if runwise_data[dataset_name].any():
                        run_group.create_dataset(dataset_name, data = runwise_data[dataset_name])
                    else:
                        run_group.create_dataset(dataset_name, data = 999)


            ###save individual maps (not sorted)
            '''
            #yet to integrate option to save individual maps to txt files
            outputfolder_individumaps = '....\\individual_mm_states'
            np.savetxt(op.join(outputfolder_individumaps, 'individual_mstate_%s_%s.txt' % (looper2, looper1)), individu_mstate, delimiter=' ')
            '''

            ##epochwise output
            for epochnr in range(number_of_epochs):
                #create folder
                if "ep_"+"%03d" % (epochnr,) in run_group.keys():
                    ep_group = k['/{0}/{1}/{2}/{3}/{4}' .format(output_path.level0, output_path.level1, output_path.level2, output_path.level3, "ep_"+"%03d" % (epochnr,))]
                else:
                    ep_group = run_group.create_group("ep_"+"%03d" % (epochnr,))  
                    if confobj.debug:
                        print(epochnr)

                #add attributes to epoch folder
                for attribute_name, attribute_value in output_attributes.iteritems():
                    ep_group.attrs['{0}' .format(attribute_name)]=str(attribute_value[epochnr])
    
                #add epochwise datasets
                for epochwise_data_name in epochwise_data:
                    if not epochwise_data_name in ep_group.keys():
                        if epochwise_data[epochwise_data_name][epochnr].any():
                            ep_group.create_dataset(epochwise_data_name, data = epochwise_data[epochwise_data_name][epochnr])
                        else:
                            ep_group.create_dataset(epochwise_data_name, data = 999)

                #add attribute to 'Start State Array'
                ep_group['Start state array'].attrs['Mstate Begin, End, Dissimilarity to optimally corresponding map']='Note that the first and last mstate were not considered for the parameter computations.'   

                #add mapwise datasets
                for mapnr in range(confobj.original_nr_of_maps):
                    if "map_"+"%02d" % (mapnr,) in ep_group.keys():
                        map_group = k['/{0}/{1}/{2}/{3}/{4}/{5}' .format(output_path.level0, output_path.level1, output_path.level2, output_path.level3, "ep_"+"%03d" % (epochnr,), "map_"+"%02d" % (mapnr,))]
                    else:
                        map_group = ep_group.create_group( "map_"+"%02d" % (mapnr,))
 
                    for mapwise_data_name in mapwise_data:
                        if not mapwise_data_name in map_group.keys():
                            if mapwise_data[mapwise_data_name][epochnr][mapnr]:
                                map_group.create_dataset(mapwise_data_name, data = mapwise_data[mapwise_data_name][epochnr][mapnr])
                            else:
                                map_group.create_dataset(mapwise_data_name, data = 0)

    #writes output to text files

    def write_output_totext(self, confobj, eeg_info_study_obj, output_data_all):
        outputfolder=op.join(op.dirname(op.abspath(self._outputfile)))
        if not op.exists(outputfolder):
            print("Create output folder: {0}".format(outputfolder))
            os.makedirs(outputfolder)

        create_parameter_spss_sheets(confobj, eeg_info_study_obj, outputfolder, output_data_all)



####Sub Class CondDataProvider --> umbenennen

###Compute parameters by Norm Maps Data Provider 1

class ParametersByNormDataProvider1(ParametersDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist):
        ParametersDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist)

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        microstate_run_value = np.loadtxt(self._sortbyfile)

        if all(microstate_run_value[0,:] == 0):
            print('Warning! {0} has all zeros, group, pt, cond ignored.'.format(sortbyhdf5))
        else:
            model_map=microstate_run_value[:]
        return model_map


###Compute parameters by Level0 Data Provider 1

class ParametersBy1LevelDataProvider1(ParametersDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist):
        ParametersDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist)

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        model_maps_all = []

        with closing( h5py.File(self._sortbyfile, 'r') ) as f: 
            #find out which level you need for sorting
            if 'all' in f['/'].keys():
                path = '/{0}' .format('all')
            elif output_path.level0 in f['/'].keys():
                path = '/{0}' .format(output_path.level0)
            elif output_path.level1 in f['/'].keys():
                path = '/{0}' .format(output_path.level1)
            elif output_path.level2 in f['/'].keys():
               path = '/{0}' .format(output_path.level2)               
            elif output_path.level3 in f['/'].keys():
                path = '/{0}' .format(output_path.level3)   
            else:
                print('Sortbypath not found for {0} {1} {2} {3} in HDF5 file {4}'.format(output_path.level0, output_path.level1, output_path.level2, output_path.level3, self._sortbyfile))
                       
            microstate_run_value = f['/{0}/{1}' .format(path, self._sortbydataset)] 
            model_map=microstate_run_value[:]
        
        return model_map

###Compute parameters by Level1 Data Provider 1

class ParametersBy2LevelsDataProvider1(ParametersDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist):
        ParametersDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist)

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        model_maps_all = []

        with closing( h5py.File(self._sortbyfile, 'r') ) as f: 
            #from itertools import permutations
            #permutations('Group Cond', 2)

            #find out which levels you need for sorting
            #find Level 0
            if output_path.level0 in f.keys():
                firstlevel = output_path.level0
            elif output_path.level1 in f.keys():
                firstlevel = output_path.level1
            elif output_path.level2 in f.keys():
                firstlevel = output_path.level2             
            elif output_path.level3 in f.keys():
                firstlevel = output_path.level3
            else:
                print('Sortbypath not found for {0} {1} {2} {3} in HDF5 file {4}'.format(output_path.level0, output_path.level1, output_path.level2, output_path.level3, self._sortbyfile))

            #find Level 1


            if output_path.level0 in f['/{0}/'.format(firstlevel)].keys():
                secondlevel = output_path.level0
            elif output_path.level1 in f['/{0}/'.format(firstlevel)].keys():
                secondlevel = output_path.level1
            elif output_path.level2 in f['/{0}/'.format(firstlevel)].keys():
               secondlevel = output_path.level2
            elif output_path.level3 in f['/{0}/'.format(firstlevel)].keys():
                secondlevel = output_path.level3
            else:
                print('Sortbypath not found for {0} {1} {2} {3} in HDF5 file {4}'.format(output_path.level0, output_path.level1, output_path.level2, output_path.level3, self._sortbyfile))

            path = '/{0}/{1}' .format(firstlevel, secondlevel)
            microstate_run_value = f['/{0}/{1}' .format(path, self._sortbydataset)] 
            model_map=microstate_run_value[:]
        
        return model_map

    
###Compute parameters by Level1 Data Provider 1

class ParametersBy3LevelsDataProvider1(ParametersDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist):
        ParametersDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist)

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        model_maps_all = []

        with closing( h5py.File(self._sortbyfile, 'r') ) as f: 
            #from itertools import permutations
            #permutations('Group Cond', 2)

            #find out which levels you need for sorting
            #find Level 0
            if output_path.level0 in f.keys():
                firstlevel = output_path.level0
            elif output_path.level1 in f.keys():
                firstlevel = output_path.level1
            elif output_path.level2 in f.keys():
                firstlevel = output_path.level2             
            elif output_path.level3 in f.keys():
                firstlevel = output_path.level3
            else:
                print('Sortbypath not found for {0} {1} {2} {3} in HDF5 file {4}'.format(output_path.level0, output_path.level1, output_path.level2, output_path.level3, self._sortbyfile))

            #find Level 1
            if output_path.level0 in f['/{0}/'.format(firstlevel)].keys():
                secondlevel = output_path.level0
            elif output_path.level1 in f['/{0}/'.format(firstlevel)].keys():
                secondlevel = output_path.level1
            elif output_path.level2 in f['/{0}/'.format(firstlevel)].keys():
               secondlevel = output_path.level2
            elif output_path.level3 in f['/{0}/'.format(firstlevel)].keys():
                secondlevel = output_path.level3
            else:
                print('Sortbypath not found for {0} {1} {2} {3} in HDF5 file {4}'.format(output_path.level0, output_path.level1, output_path.level2, output_path.level3, self._sortbyfile))

            #find Level 2
            if output_path.level0 in f['/{0}/{1}/'.format(firstlevel, secondlevel)].keys():
                thirdlevel = output_path.level0
            elif output_path.level1 in f['/{0}/{1}/'.format(firstlevel, secondlevel)].keys():
                thirdlevel = output_path.level1
            elif output_path.level2 in f['/{0}/{1}/'.format(firstlevel, secondlevel)].keys():
               thirdlevel = output_path.level2
            elif output_path.level3 in f['/{0}/{1}/'.format(firstlevel, secondlevel)].keys():
                thirdlevel = output_path.level3
            else:
                print('Sortbypath not found for {0} {1} {2} {3} in HDF5 file {4}'.format(output_path.level0, output_path.level1, output_path.level2, output_path.level3, self._sortbyfile))


            path = '/{0}/{1}/{2}' .format(firstlevel, secondlevel, thirdlevel)
            microstate_run_value = f['/{0}/{1}' .format(path, self._sortbydataset)] 
            model_map=microstate_run_value[:]
        
        return model_map


class ParametersBy4LevelsDataProvider1(ParametersDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist):
        ParametersDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist)

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        model_maps_all = []

        with closing( h5py.File(self._sortbyfile, 'r') ) as f: 
            #from itertools import permutations
            #permutations('Group Cond', 2)

            #find out which levels you need for sorting
            #find Level 0
            if output_path.level0 in f.keys():
                firstlevel = output_path.level0
            elif output_path.level1 in f.keys():
                firstlevel = output_path.level1
            elif output_path.level2 in f.keys():
                firstlevel = output_path.level2             
            elif output_path.level3 in f.keys():
                firstlevel = output_path.level3
            else:
                print('Sortbypath not found for {0} {1} {2} {3} in HDF5 file {4}'.format(output_path.level0, output_path.level1, output_path.level2, output_path.level3, self._sortbyfile))

            #find Level 1
            if output_path.level0 in f['/{0}/'.format(firstlevel)].keys():
                secondlevel = output_path.level0
            elif output_path.level1 in f['/{0}/'.format(firstlevel)].keys():
                secondlevel = output_path.level1
            elif output_path.level2 in f['/{0}/'.format(firstlevel)].keys():
               secondlevel = output_path.level2
            elif output_path.level3 in f['/{0}/'.format(firstlevel)].keys():
                secondlevel = output_path.level3
            else:
                print('Sortbypath not found for {0} {1} {2} {3} in HDF5 file {4}'.format(output_path.level0, output_path.level1, output_path.level2, output_path.level3, self._sortbyfile))

            #find Level 2
            if output_path.level0 in f['/{0}/{1}/'.format(firstlevel, secondlevel)].keys():
                thirdlevel = output_path.level0
            elif output_path.level1 in f['/{0}/{1}/'.format(firstlevel, secondlevel)].keys():
                thirdlevel = output_path.level1
            elif output_path.level2 in f['/{0}/{1}/'.format(firstlevel, secondlevel)].keys():
               thirdlevel = output_path.level2
            elif output_path.level3 in f['/{0}/{1}/'.format(firstlevel, secondlevel)].keys():
                thirdlevel = output_path.level3
            else:
                print('Sortbypath not found for {0} {1} {2} {3} in HDF5 file {4}'.format(output_path.level0, output_path.level1, output_path.level2, output_path.level3, self._sortbyfile))

            #find Level 3
            if output_path.level0 in f['/{0}/{1}/{2}/'.format(firstlevel, secondlevel, thirdlevel)].keys():
                fourthlevel = output_path.level0
            elif output_path.level1 in f['/{0}/{1}/{2}/'.format(firstlevel, secondlevel, thirdlevel)].keys():
                fourthlevel = output_path.level1
            elif output_path.level2 in f['/{0}/{1}/{2}/'.format(firstlevel, secondlevel, thirdlevel)].keys():
               fourthlevel = output_path.level2
            elif output_path.level3 in f['/{0}/{1}/{2}/'.format(firstlevel, secondlevel, thirdlevel)].keys():
                fourthlevel = output_path.level3
            else:
                print('Sortbypath not found for {0} {1} {2} {3} in HDF5 file {4}'.format(output_path.level0, output_path.level1, output_path.level2, output_path.level3, self._sortbyfile))

            path = '/{0}/{1}/{2}/{3}' .format(firstlevel, secondlevel, thirdlevel, fourthlevel)
            microstate_run_value = f['/{0}/{1}' .format(path, self._sortbydataset)] 
            model_map=microstate_run_value[:]
        
        return model_map




def parameter_preprocessing(confobj, state_match_percentage_all_epochs):
     map_avg_perc = np.zeros((confobj.original_nr_of_maps, 2))

     for map in range(confobj.original_nr_of_maps):
        listli=[]
        for epochnr in state_match_percentage_all_epochs.keys():
            listli.append(state_match_percentage_all_epochs[epochnr][map,1])
        map_avg_perc[map,0]=map
        map_avg_perc[map,1]=nanmean(listli)


     return map_avg_perc




##############################################
#######  compute_mstate_parameters  ##########
##############################################

##Compute mstate parameters


def compute_mstate_parameters(confobj, eeg, maps, eeg_info_study_obj):
    TF = eeg_info_study_obj.tf
    Fs = eeg_info_study_obj.sf

    #####Create Dictionaries for parameters across all epochs
    dur_state_all_epochs = dict.fromkeys(list(range(len(eeg)/TF)))
    freq_dict_all_epochs = dict.fromkeys(list(range(len(eeg)/TF)))
    dur_dict_all_epochs = dict.fromkeys(list(range(len(eeg)/TF)))
    cov_dict_all_epochs = dict.fromkeys(list(range(len(eeg)/TF)))
    gfp_peak_nr_all_epochs = dict.fromkeys(list(range(len(eeg)/TF)))
    mean_gfp_all_epochs = dict.fromkeys(list(range(len(eeg)/TF)))
    gfp_mean_all_epochs = dict.fromkeys(list(range(len(eeg)/TF)))
    durstd_dict_all_epochs= dict.fromkeys(list(range(len(eeg)/TF)))
    start_state_list_all_epochs= dict.fromkeys(list(range(len(eeg)/TF)))
    exp_var_all_epochs= dict.fromkeys(list(range(len(eeg)/TF)))
    exp_var_tot_all_epochs= dict.fromkeys(list(range(len(eeg)/TF)))
    state_match_percentage_all_epochs= dict.fromkeys(list(range(len(eeg)/TF)))
    state_match_percentage_std_all_epochs= dict.fromkeys(list(range(len(eeg)/TF)))
    gfp_curves_all_epochs= dict.fromkeys(list(range(len(eeg)/TF)))

    #################  
    # 3.) LOOP ACROSS 2 Sec Segments
    #################

    individu_dict = dict.fromkeys(list(range(confobj.original_nr_of_maps)))
    for mstate in individu_dict:
        individu_dict[mstate]=np.zeros((eeg.shape[1]))

    individu_mstate = np.zeros((confobj.original_nr_of_maps, int(eeg.shape[1])))

    for epochnr in range(len(eeg)/TF):
        epoch = eeg[epochnr*TF:(epochnr+1)*TF]
  
        #################   
        # 3.) COMPUTE GFP (ln1 or ln2)
        #################
        
        gfp_curve = compute_gfp(epoch, confobj.method_GFPpeak)
        if confobj.debug:
            print('GFP Curve computed')

        #################   
        # 3.) Average ref prior to gfp_1
        #################         
        
        #average ref epoch
        epoch_mean = epoch.mean(axis=1)
        epoch = epoch - epoch_mean[:, np.newaxis] 

        #average ref maps
        maps_mean = maps.mean(axis=1)
        maps = maps-maps_mean[:, np.newaxis] 
        
        #################   
        # 3.) Set to gfp_1 prior to correlation computation
        #################         
        
        ##set to gfp_1 prior to correlation computation
        #set to gfp_1 of all tfs
        epoch = set_gfp_all_1(epoch, gfp_curve)

        #set to gfp_1 of all external maps
        gfp_curve_maps = compute_gfp(maps, confobj.method_GFPpeak)
        maps = set_gfp_all_1(maps, gfp_curve_maps)

        #################   
        #4.) Compute GFP Peaks (identical to microstates.py)
        #################
        
        gfp_peak_indices, gfp_curve = compute_gfp_peaks(gfp_curve, confobj.use_gfp_peaks, confobj.use_smoothing, confobj.gfp_type_smoothing, confobj.smoothing_window, confobj.use_fancy_peaks)
       
        ################# 
        #5.) Determine Mstate class for each peak
        #################

        #Compare each gfp peak index (gfp_peak_indices) to the 4 mstate maps (maps) --> get attribution matrix for the eeg 0,1,3,2,0,2 etc.

        #initialize attribution matrix (one with the correlations and one only with the highest gfp_peak_index)
        attribution_matrix= dict.fromkeys(gfp_peak_indices)
        attribution_matrix_indices= dict.fromkeys(gfp_peak_indices)

        #get index of mstate map that correlates highest with gfp_peak_index

        for key in attribution_matrix.keys():
            tf=epoch[key]
            #tf=eeg[key]
            corr_list = []
            for mapnr in range(len(maps)):
                #Pearson Correlation
                if confobj.similarity_measure == 'correlation':
                    pr, pp=pearsonr( maps[mapnr,:], tf)
                elif confobj.similarity_measure == 'dissimilarity':
                    #Dissimilarity
                    pr = dissim(maps[mapnr,:], tf, confobj.method_GFPpeak)
                else:
                    print('Error confobj.similarity_measure must be correlation or dissimilarity but it is {0}'.format(confobj.similarity_measure))

                #abs specific for EEG / continuous data where polarity is disregarded
                corr_list.append(abs(pr))

            attribution_matrix[key]=corr_list

            ###Find out what to do with dissimilarity below!

            #confobj.correspondance_cutoff = False for no cutoff
            if confobj.similarity_measure == 'correlation':
                if confobj.correspondance_cutoff == 0 or max(corr_list) > confobj.correspondance_cutoff:
                    #delete below
                    #attribution_matrix_indices[key]=corr_list.index(max(corr_list))
                    attribution_matrix_indices[key]=[corr_list.index(max(corr_list)), max(corr_list)]
                else:
                    attribution_matrix_indices[key]=[999, max(corr_list)]   

            elif confobj.similarity_measure == 'dissimilarity':
                if confobj.correspondance_cutoff == 0 or min(corr_list) < confobj.correspondance_cutoff:
                    #delete below
                    #attribution_matrix_indices[key]=corr_list.index(min(corr_list))
                    attribution_matrix_indices[key]=[corr_list.index(min(corr_list)), min(corr_list)]
                else:
                    attribution_matrix_indices[key]=[999, min(corr_list)] 
            else:
                print('Error confobj.similarity_measure must be correlation or dissimilarity but it is {0}'.format(confobj.similarity_measure))


        #tf_begin is determined as, the first gfp peak where the mstate changes from one to another, e.g. A A A (B) B B (A) (B) B

        start_state_list =[]
        previous=sorted(attribution_matrix_indices)[0]
        start = sorted(attribution_matrix_indices)[0]
        for ele in sorted(attribution_matrix_indices):

            if attribution_matrix_indices[ele][0] == attribution_matrix_indices[previous][0]:
                previous = ele
                continue
                #print 'the same value, no change'
            else:
                #start defined as the middle between GFP Peak A and GFP Peak B (middle between the changed ones)
                start = previous+((ele-previous)/2.)
                start_state = start, attribution_matrix_indices[ele]
                start_state_list.append(start_state)
                previous = ele   
                #print 'different value'
                #print 'saved', start
      

        #Compute mstate duration (in tfs & converted into ms) for each class

        #dictionary of the 4 states
        dur_state=dict.fromkeys(list(range(confobj.original_nr_of_maps)))
        gfp_state=dict.fromkeys(list(range(confobj.original_nr_of_maps)))
                  
        i = 0
        while (i < len(start_state_list)-1):
            #print i
            
            beg=start_state_list[i][0]
            end=start_state_list[i+1][0]
            dur=end-beg
            state=start_state_list[i][1][0]
            if dur_state[state] == None:
                dur_state[state] = []
            dur_state[state].append(dur *(1000/float(Fs)))

            if gfp_state[state] == None:
                gfp_state[state] = []
            gfp_state[state].append(np.mean(gfp_curve[int(beg+1):(int(end)+1)]))

            
            i=i+1

            if (dur_state is None):
                if confobj.debug:
                    print(outti, inni, dur_state)




        #ignore last value pair from start_state_list because you don't know when this state ends
        ####
        #Compute total possible duration in ms across mstates
        ####
        dur_sum = 0
        for i in range(confobj.original_nr_of_maps):
            if dur_state[i] == None:
                pass
            else:
                dur_sum=dur_sum+np.sum(dur_state[i]) 

        ####
        #Frequency
        ####

        freq_dict=dict.fromkeys(list(range(confobj.original_nr_of_maps)))
        for mapnr in range(confobj.original_nr_of_maps):
            if dur_state[mapnr] == None:
                print('epochnr, mapnr no content', epochnr, mapnr)
                freq_dict[mapnr] = 0
            else:
                #freq corrected by smaller epoch size (returns # of occurrences per second based on info from whole epoch)
                freq_dict[mapnr] = (len(dur_state[mapnr])*(TF*(1000/float(Fs)/dur_sum)/float(TF/float(Fs))))


       
        ####
        #Mean Duration & Std from tf to secs
        ####

        dur_dict=dict.fromkeys(list(range(confobj.original_nr_of_maps)))
        durstd_dict=dict.fromkeys(list(range(confobj.original_nr_of_maps)))
        for mapnr in range(confobj.original_nr_of_maps): 
            if dur_state[mapnr] == None:
                dur_dict[mapnr] = 0
                durstd_dict[mapnr] = 0
            else:
                dur_dict[mapnr] = np.mean(dur_state[mapnr])
                durstd_dict[mapnr] = np.std(dur_state[mapnr]) 

        ####
        #Coverage (in % of 2 s epoch)
        ####


        #Compute relative coverage for each state
        cov_dict=dict.fromkeys(list(range(confobj.original_nr_of_maps)))
        for mapnr in range(confobj.original_nr_of_maps): 
            if dur_state[mapnr] == None:
                cov_dict[mapnr] = 0
            else:
                cov_dict[mapnr] = np.sum(dur_state[mapnr]) /float(dur_sum)
            

        ####
        #Number of GFP Peaks in Epoch
        ####
        gfp_peak_nr=len(gfp_peak_indices)

        ####
        #Mean GFP across all classes (computed based on the mean gfp per class weighted by their coverage)
        ####
        mean_gfp_all=np.mean(gfp_curve)

        ####
        #Mean GFP for each class
        ####
        gfp_mean=dict.fromkeys(list(range(confobj.original_nr_of_maps)))

        for mapnr in range(confobj.original_nr_of_maps):  
            if  gfp_state[mapnr] == None:
                gfp_mean[mapnr]="-"
            else:
                gfp_mean[mapnr]=np.mean(gfp_state[mapnr])

        ####
        #Start_state_list from timeframes to ms
        ####

        #start_state_list = np.array(start_state_list)

        #Get array of start_state_list
        start_state_list_array = np.zeros((len(start_state_list),3))

        for ele1, ele2 in enumerate(start_state_list):
            start_state_list_array[ele1][0]=ele2[0]*(1000./Fs)
            start_state_list_array[ele1][1]=ele2[1][0]
            start_state_list_array[ele1][2]=ele2[1][1]

        start_state_list = start_state_list_array
 
        ####
        #Explained variance (TODO: correct value retrieved when the modelmaps have a different number of channels than the EEG)
        ####     
        
        #TODO: with / without additional norm computation
        model = maps

        #TK Normierung
        #Berechnung des Normierungs-Vektors, es geht darum alles auf Vektorlaenge 1 zu setzen
        b=np.sum(np.abs(model)**2,axis=-1)**(1./2)    
        #Teilung aller Elemente durch Normierungs-Vektor
        for col in range(model.shape[1]):
            model[:,col]=model[:,col]/b       

        #Covariance Matrix computation
        covm=np.dot(epoch[gfp_peak_indices],model.T)
        covm_all=np.dot(epoch,model.T)

        loading=abs(covm).max(axis=1)
        loading_all=abs(covm_all).max(axis=1)

        b_loading=loading/sqrt(model.shape[1])
        b_loading_all=loading_all/sqrt(model.shape[1])
        		
        exp_var=sum(b_loading)/sum(epoch[gfp_peak_indices].std(axis=1))
        exp_var_tot=sum(b_loading_all)/sum(epoch.std(axis=1))


        ###Compute Percentage of Correspondance between tf & labelby map
        state_match_percentage=dict.fromkeys(list(range(confobj.original_nr_of_maps)))
        state_match_percentage_std=dict.fromkeys(list(range(confobj.original_nr_of_maps)))

        for keyli in state_match_percentage.keys():
            state_match_percentage[keyli]=[]
            for ele in start_state_list:
                if ele[1] == float(keyli):
                    state_match_percentage[keyli].append(ele[2])
            #compute mean across percentages of same state
            state_match_percentage_std[keyli]=np.std(state_match_percentage[keyli])
            state_match_percentage[keyli]=np.mean(state_match_percentage[keyli])
            
            ##Coversion dict to array
            s_m_p=np.zeros((len(state_match_percentage),2))

            for keynr, key in enumerate(sorted(state_match_percentage)):
                s_m_p[key,0]=key
                s_m_p[key,1]=state_match_percentage[key]

            s_m_p_std=np.zeros((len(state_match_percentage_std),2))

            for keynr, key in enumerate(sorted(state_match_percentage_std)):
                s_m_p_std[key,0]=key
                s_m_p_std[key,1]=state_match_percentage_std[key]



        #convert to array
        state_match_percentage = s_m_p
        state_match_percentage_std = s_m_p_std


        #####get all tfs for epochnr

        for mstate in range(confobj.original_nr_of_maps):
            for ele in attribution_matrix_indices:
                #print 'mstate', mstate
                if attribution_matrix_indices[ele][0] == mstate:                  
                    individu_dict[mstate]=np.vstack((individu_dict[mstate], epoch[ele]))





        ##----------------------------------------------------##

        ############
        ###Dictionaries with entry for each epoch
        ###########
        #Compute mstate duration (in tfs) for each class for each time the class occurs
        #dur_state
        #print dur_state_all_epochs 
        dur_state_all_epochs[epochnr] = dur_state
        #Frequency
        #freq_dict
        freq_dict_all_epochs[epochnr] = freq_dict
        #Mean Duration & STD
        #dur_dict
        dur_dict_all_epochs[epochnr] = dur_dict
        #durstd_dict
        durstd_dict_all_epochs[epochnr] = durstd_dict
        #Coverage (in % of 2 s epoch)
        #cov_dict
        cov_dict_all_epochs[epochnr] = cov_dict
        #Number of GFP Peaks in Epoch (total epoch)
        #gfp_peak_nr
        gfp_peak_nr_all_epochs[epochnr] = gfp_peak_nr
        #Mean GFP across all classes
        #mean_gfp_all
        mean_gfp_all_epochs[epochnr] = mean_gfp_all
        #Mean GFP for each class
        #gfp_mean
        gfp_mean_all_epochs[epochnr] = gfp_mean
        #start_state_list
        start_state_list_all_epochs[epochnr] = start_state_list
        #exp_var
        exp_var_all_epochs[epochnr]= exp_var
        #exp_var_tot
        exp_var_tot_all_epochs[epochnr]= exp_var_tot
        #state_match_percentage
        state_match_percentage_all_epochs[epochnr]= state_match_percentage
        #state_match_percentage_std
        state_match_percentage_std_all_epochs[epochnr]= state_match_percentage_std
        #gfp_curve_epoch
        gfp_curves_all_epochs[epochnr]= gfp_curve

        #RETURN ABOVE DICTIONARIES

    ##compute individual microstates via principal components
    for mstate in range(confobj.original_nr_of_maps):
        #skip first row because its zeros
        try:
            P=individu_dict[mstate][1:,:]
            coeff = princomp_B(P,1)
            individu_mstate[mstate,:] = coeff.ravel()
        except:
            individu_mstate[mstate,:] = np.zeros((eeg.shape[1]))


    ###rename dictionaries to more easily comprehend their meaning

    occ = freq_dict_all_epochs
    dur = dur_dict_all_epochs
    cov = cov_dict_all_epochs


    ###Convert all measures into two dictionaries depending on whether they represent a dataset or attribute
 
    #create dictionary for runwise datasets (keys are the descriptions that will be in the hdf5 outputfile) 
    runwise_data = {}
    runwise_data['Individual_States']=individu_mstate
       
    #create dictionary for epochwise datasets (keys are the descriptions that will be in the hdf5 outputfile) 
    epochwise_datasets = [start_state_list_all_epochs, state_match_percentage_all_epochs, state_match_percentage_std_all_epochs, gfp_curves_all_epochs]
    epochwise_datasets_string = ['Start state array', 'State Match Mean percentage', 'State Match Std percentage', 'GFP Curve']

    epochwise_data = {}
    for elenr, ele in enumerate(epochwise_datasets):
        ele_string = epochwise_datasets_string[elenr]
        epochwise_data[ele_string]=ele
        
    #create dictionary for mapwise datasets (keys are the descriptions that will be in the hdf5 outputfile)            
    mapwise_datasets = [occ, dur, durstd_dict_all_epochs, cov, gfp_mean_all_epochs, dur_state_all_epochs]
    mapwise_datasets_string = ['Occurrance per s', 'Mean duration in ms', 'SD duration in ms', 'Coverage in percent', 'GFP Mean across all tfs', 'number of ms for each state']

    mapwise_data = {}
    for elenr, ele in enumerate(mapwise_datasets):
        ele_string = mapwise_datasets_string[elenr]
        mapwise_data[ele_string]=ele


    #create list of all output_datasets
    output_data = runwise_data, epochwise_data, mapwise_data

    #create dictionary for attributes (keys are the descriptions that will be in the hdf5 outputfile)
    attribute_measures = [gfp_peak_nr_all_epochs, mean_gfp_all_epochs, exp_var_all_epochs, exp_var_tot_all_epochs]
    attribute_measures_string = ['Number of GFP Peaks in Epoch', 'Mean GFP in Epoch', 'Explained Variance GFP peaks', 'Explained Variance EEG epoch']

    output_attributes = {}
    for elenr, ele in enumerate(attribute_measures):
        ele_string = attribute_measures_string[elenr]
        output_attributes[ele_string]=ele

    return output_data, output_attributes



def create_parameter_spss_sheets(confobj, eeg_info_study_obj, outputfolder, output_data_all):                   
    ##################################################################################################################################################################
    ##################################################################################################################################################################
    #########################################################         Prepare SPSS Sheets     ########################################################################
    ##################################################################################################################################################################
    ##################################################################################################################################################################
    import csv
    three__measures = ['Occurrance per s', 'Mean duration in ms', 'Coverage in percent']

    #Map long names to short names for better naming in file
    short_names_measures = {}
    short_names_measures['Occurrance per s']='occ'
    short_names_measures['Mean duration in ms']='dur'
    short_names_measures['Coverage in percent']='cov'

    spss_parameters_csv = op.join( outputfolder, 'spss_parameters.csv')


    #get maximal number of epochs across all group pt cond run
    max_len=0
    for output_data_path, output_data_per_path in output_data_all.iteritems():
        runwise_data, epochwise_data, mapwise_data = output_data_per_path
        #count number of epochs
        curr_len=len(mapwise_data['Occurrance per s'])
        if max_len < curr_len:
            max_len = curr_len


    #get list of all Pts
    pt_set = Set()
    cond_set = Set()
    run_set = Set()

    for output_data_path in output_data_all.keys():
        #in order to ensure that all participants are included, the participant is characterized by its name and its group
        pt_set.add("{0} {1}" .format(output_data_path.level0, output_data_path.level1))
        cond_set.add(output_data_path.level2)
        run_set.add(output_data_path.level3)

    pt_list=sorted(pt_set)
    cond_list=sorted(cond_set)
    run_list=sorted(run_set)

      
    #########
    ###File with mean across epochs for each measure
    #########

    #create dictionary with key levels: measure, pti, cond, run
    parameters_mean=dict.fromkeys(three__measures)
    for meas in three__measures:
        parameters_mean[meas]=dict.fromkeys(pt_list)
        for pti in pt_list:
            parameters_mean[meas][pti]=dict.fromkeys(cond_list)
            for condi in cond_list:
                parameters_mean[meas][pti][condi]=dict.fromkeys(run_list)
                for runi in run_list:
                    parameters_mean[meas][pti][condi][runi]=dict.fromkeys(range(confobj.original_nr_of_maps))


    for meas in three__measures:           
        for pti in pt_list:
            for condi in cond_list:
                for runi in run_list:
                    for output_data_path, output_data_per_path in output_data_all.iteritems():
                        runwise_data, epochwise_data, mapwise_data = output_data_per_path
                        if "{0} {1}" .format(output_data_path.level0, output_data_path.level1) == pti and output_data_path.level2 == condi and output_data_path.level3 == runi:
                            if mapwise_data[meas]:
                                for mapnr in range(confobj.original_nr_of_maps):   
                                    list_to_avg = []
                                    for epochnr in range(len(mapwise_data[meas])):
                                        list_to_avg.append(mapwise_data[meas][epochnr][mapnr])
                                    parameters_mean[meas][pti][condi][runi][mapnr]=np.mean(list_to_avg)
                            else:
                                print("Warning for {0} {1} {2} {3}".format(output_data_path.level0, output_data_path.level1, output_data_path.level2, output_data_path.level3))


    ###Seperate file for each measure
    for meas in three__measures:
        header = []
        header.append('Pt; Group;')
        for condi in cond_list:
            for runi in run_list:
                for mapnr in range(confobj.original_nr_of_maps):
                    header.append('{0}_{1}{2}_map{3};'.format(short_names_measures[meas], condi.split('_')[1], runi.split('_')[1], mapnr))


        spss_parameters_csv = op.join( outputfolder, "{0}_means.csv" . format(short_names_measures[meas]))
        with open(spss_parameters_csv, 'w') as spss_parameter_file:
            spss_parameter_file.writelines(header)
            for pti in pt_list:
                spss_parameter_file.write('\n')
                spss_parameter_file.write('{1};{0}'.format(pti.split()[0],pti.split()[1]))
                for condi in cond_list:
                    for runi in run_list:
                        for mapnr in range(confobj.original_nr_of_maps):
                            spss_parameter_file.write(';')
                            if parameters_mean[meas][pti][condi][runi][mapnr]:
                                spss_parameter_file.write("{0:.2f}".format(parameters_mean[meas][pti][condi][runi][mapnr]))
                            else:
                                print('No data available for {0} {1} {2} {3}'.format(pti, condi, runi, mapnr))
                                spss_parameter_file.write("{0}".format(999))









##########################
### get data provider  ###
##########################

from os.path import basename
import os.path

def get_data_provider_for_parameter_by(parameter_by, inputfolder, hdf5_filename, inputdataset, sortbyfolder, sortbyfile, sortbydataset, sortbyseries, external_chlist):
    if parameter_by == 'external_norm':
        ######################
        ###   sort by norm ###
        ######################

        # folder, file and dataset - data that the parameters are computed upon
        inputfolder = inputfolder
        inputhdf5 = os.path.join( inputfolder, hdf5_filename)
        inputdataset = inputdataset

        # the folder path to the microstates to sortby (for parameter computation) file
        sortbyfolder = sortbyfolder
        sortbyhdf5 = os.path.join(sortbyfolder, sortbyfile)
        sortbydataset = sortbydataset 

        sortbychlist = os.path.join(sortbyfolder,external_chlist)

        # the folder path to the output hdf5 file
        #if not os.path.exists(outputfolder):
        #    os.makedirs(outputfolder)
        outputhdf5 = os.path.join(sortbyfolder, basename(sortbyhdf5).split('.')[0] ,'mstate_parameters.hdf5')

        ###Establish DataProvider
        data_provider=ParametersByNormDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist)

    elif parameter_by == 'own_hdf':
        # the folder path to the all_recoredings.hdf file
        library_path = os.path.dirname(os.path.abspath(__file__))
        inputhdf5 = os.path.join( inputfolder, hdf5_filename)
        inputdataset = inputdataset

        # the folder path to the microstates to sortby (for parameter computation) file
        sortbyseries = sortbyseries
        sortbyfolder = inputfolder
        sortbyhdf5 = os.path.join(sortbyfolder, hdf5_filename)

        #will not be used if you select a file in a folder to categorize the maps based on
        sortbydataset = sortbydataset
        sortbychlist = False

        # the folder path to the output hdf5 file
        outputhdf5 = os.path.join(sortbyfolder, basename(sortbyhdf5).split('.')[0] ,'mstate_parameters.hdf5')

        data_provider=ParametersBy4LevelsDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist)

    else:
        ######################
        ###  sort by hdf5  ###
        ######################       
            
        # the folder path to the all_recoredings.hdf file
        library_path = os.path.dirname(os.path.abspath(__file__))
        inputhdf5 = os.path.join( inputfolder, hdf5_filename)
        inputdataset = inputdataset

        # the folder path to the microstates to sortby (for parameter computation) file
        sortbyseries = sortbyseries
        sortbyfolder = os.path.join(inputfolder, '{0}' .format(sortbyseries))
        sortbyhdf5 = os.path.join(sortbyfolder, sortbyfile)

        #will not be used if you select a file in a folder to categorize the maps based on
        sortbydataset = sortbydataset
        sortbychlist = False

        # the folder path to the output hdf5 file
        outputhdf5 = os.path.join(sortbyfolder, basename(sortbyhdf5).split('.')[0] ,'mstate_parameters.hdf5')

        
        if parameter_by == '1Level':
            data_provider=ParametersBy1LevelDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist)
        elif parameter_by == '2Levels':
            data_provider=ParametersBy2LevelsDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist)
        elif parameter_by == '3Levels':
            data_provider=ParametersBy3LevelsDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist)
        elif parameter_by == '4Levels':
            data_provider=ParametersBy4LevelsDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist)
        else:
            print('Warning: data_provider could not be specified based on the given information. parameter_by must be specified.')


    return data_provider





def run_parameters(data_provider, confobj, eeg_info_study_obj):
    #create dictionary with parameter statistics for all group cond pt run
    output_data_all = {}

    for output_path in data_provider.get_outputs():
        input = data_provider.get_input_data(output_path, eeg_info_study_obj.chlist)
        sortby = data_provider.get_sortby_data(output_path)
        #compute_mstate_parameters demands that input and sortby are in the same data format (equal nch)
        output_data, output_attributes = compute_mstate_parameters(confobj, input, sortby, eeg_info_study_obj)
        if not output_data == []:
            data_provider.write_output_data(confobj, output_path, output_data, output_attributes)

        #add data for output_path to dictionary
        output_data_all[output_path]=output_data

    #create spss sheets for all group cond pt run based on dictionary: output_data_all
    #get location for outputhdf5 for folder 

    data_provider.write_output_totext(confobj, eeg_info_study_obj, output_data_all)








