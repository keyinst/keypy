##################################
#######  Parameters Provider  ####
##################################

from __future__ import print_function

from scipy.stats import nanmean, pearsonr
from keypy.microstates.microstates_helper import *
from numpy import sqrt
import os.path as op
import os
from sets import Set
import numpy as np
from keypy.microstates.parameters import *
from keypy.microstates.sortmaps_provider import *
from contextlib import closing
import h5py
####-------------------------------------####

##########################
#######  Functions  ########
##########################

def dissim(map, tf, method_GFPpeak):
    #Assert that input maps have been normalized identically
    #either with GFP=1 or vector length = 1
    #if not compute_gfp(map).all() == 1 == compute_gfp(tf).all() or not (LA.norm(map[ri,:], axis=0)).all() == 1 == (LA.norm(tf[ri,:], axis=0)).all():
    #    raise AssertionError('Dissimilarity computation demands identical normalization for map and tf.')

    diff1=map-tf
    diff2=-map-tf
    diff1=np.reshape(diff1, (-1, 1))
    diff2=np.reshape(diff2, (-1, 1))
    diff1_gfp=compute_gfp(diff1.T, method_GFPpeak)
    diff2_gfp=compute_gfp(diff2.T, method_GFPpeak)
    diff=min(diff1_gfp,diff2_gfp)

    return diff
####-------------------------------------####

def parameter_preprocessing(confobj, state_match_percentage_all_epochs):
     map_avg_perc = np.zeros((confobj.original_nr_of_maps, 2))

     for map in range(confobj.original_nr_of_maps):
        listli=[]
        for epochnr in state_match_percentage_all_epochs.keys():
            listli.append(state_match_percentage_all_epochs[epochnr][map,1])
        map_avg_perc[map,0]=map
        map_avg_perc[map,1]=nanmean(listli)

     return map_avg_perc
####-------------------------------------####

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
                                    ###Attention: Special case for mean duration in ms. If mean duration is zero, the mstate of this class never occurred. It should not be considered for the mean across epochs computation.
                                    if meas == 'Mean duration in ms':
                                        list_to_avg_nozero=[x for x in list_to_avg if x != 0]
                                        parameters_mean[meas][pti][condi][runi][mapnr]=np.mean(list_to_avg_nozero)
                                    else:
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
                                spss_parameter_file.write("{0:.10f}".format(parameters_mean[meas][pti][condi][runi][mapnr]))
                            else:
                                print('No data available for {0} {1} {2} {3}'.format(pti, condi, runi, mapnr))
                                spss_parameter_file.write("{0}".format(999))

####-------------------------------------####
####-------------------------------------####
####-------------------------------------####

##########################
#######  Classes  ########
##########################


####Abstract Class ParametersDataProvider

class ParametersDataProvider(object):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, sortbychlist):
        self._file = inputhdf5
        self._sortbyfile = sortbyhdf5
        self._outputfile = outputhdf5
        self._inputdataset = inputdataset
        self._sortbydataset = sortbydataset
        self._sortbychlist_path = sortbychlist

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
                eeg_own=eeg_value[:]

                if self._sortbychlist_path:
                    eeg_ext_path = self._sortbychlist_path

                    eeg_own_new,_ = reduce_channels(eeg_own, self._sortbyfile, own_chlist, self._sortbychlist_path)

                    #eeg=reduce_channels(eeg, own_chlist, self._sortbychlist)
                else:
                    eeg_own_new = eeg_own
        return eeg_own_new
       
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
        #loads the channel reordered sortby template (reodering done by reduce_channels function in microstate_helper.py)
        microstate_run_value = np.loadtxt(os.path.join(os.path.dirname(self._sortbyfile),"{0}_reduced.asc".format(os.path.splitext(os.path.basename(self._sortbyfile))[0])))

       

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

####-------------------------------------####
####-------------------------------------####
####-------------------------------------####