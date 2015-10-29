
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
from sets import Set
from keypy.microstates.data_provider import *


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

def get_io_sortmap_for_series(series, iteration, inputfolder, sortbyfolder, outputfolder, first_input):
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
    sortdata_provider = False

    if series == 'Series_1':
        if iteration == 0:
            ######
            ##sort modelmaps across groups by mean_models_milz_etal_2015
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_1', 'modelmaps_across_groups.hdf')
            sortbyhdf5 = os.path.join(sortbyfolder, 'mean_models_milz_etal_2015.asc')
            outputhdf5 = os.path.join( outputfolder, 'Series_1', 'modelmaps_across_groups{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortAllByNormDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        elif iteration == 1:
            ######
            ##sort modelmaps_across_pts by modelmaps_across_groups_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_1', 'modelmaps_across_pts.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_1', 'modelmaps_across_groups_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_1', 'modelmaps_across_pts{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortGroupByAllDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        elif iteration == 2:
            ######
            ##sort modelmaps_across_conds by modelmaps_across_pts_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_1', 'modelmaps_across_conds.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_1', 'modelmaps_across_pts_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_1', 'modelmaps_across_conds{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortPtByGroupDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        elif iteration == 3:
            ######
            ##sort sort modelmaps_across_runs by modelmaps_across_conds_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_1', 'modelmaps_across_runs.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_1', 'modelmaps_across_conds_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_1', 'modelmaps_across_runs{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortCondByPtDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        elif iteration == 4:
            ######
            ##sort microstates by modelmaps_across_runs_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_1', 'modelmaps_across_runs_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')
            inputdataset = first_input
            sortbydataset = 'modelmap'
            outputdataset = '{0}{1}' .format(first_input, '_Series_1_sorted')
            sortdata_provider = SortRunByCondDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        else:
            stop = True


    elif series == 'Series_2':
        if iteration == 0:
            ######
            ##sort modelmaps across runs by mean_models_milz_etal_2015
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_2', 'modelmaps_across_runs.hdf')
            sortbyhdf5 = os.path.join(sortbyfolder, 'mean_models_milz_etal_2015.asc')
            outputhdf5 = os.path.join( outputfolder, 'Series_2', 'modelmaps_across_runs{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortAllByNormDataProvider2(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        elif iteration == 1:
            ######
            ##sort modelmaps_across_pts by modelmaps_across_runs
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_2', 'modelmaps_across_pts.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_2', 'modelmaps_across_runs_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_2', 'modelmaps_across_pts{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortCondByPtDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        elif iteration == 2:
            ######
            ##sort microstates by modelmaps_across_pts_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_2', 'modelmaps_across_pts_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')
            inputdataset = first_input
            sortbydataset = 'modelmap'
            outputdataset = '{0}{1}' .format(first_input, '_Series_2_sorted')
            sortdata_provider = SortRunByCondDataProvider2(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        else:
            stop = True

    elif series == 'Series_3':
        if iteration == 0:
            ######
            ##sort modelmaps across conds by mean_models_milz_etal_2015
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_conds.hdf')
            sortbyhdf5 = os.path.join(sortbyfolder, 'mean_models_milz_etal_2015.asc')
            outputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_conds{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortAllByNormDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        elif iteration == 1:
            ######
            ##sort modelmaps_across_groups by modelmaps_across_conds
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_groups.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_conds_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_groups{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortGroupByAllDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        elif iteration == 2:
            ######
            ##sort modelmaps_across_pts by modelmaps_across_groups_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_pts.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_groups_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_pts{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortPtByGroupDataProvider2(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        elif iteration == 3:
            ######
            ##sort sort modelmaps_across_runs by modelmaps_across_pts_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_runs.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_pts_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_runs{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortCondByPtDataProvider2(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        elif iteration == 4:
            ######
            ##sort microstates by modelmaps_across_runs_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_3', 'modelmaps_across_runs_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')
            inputdataset = first_input
            sortbydataset = 'modelmap'
            outputdataset = '{0}{1}' .format(first_input, '_Series_3_sorted')
            sortdata_provider = SortRunByCondDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        else:
            stop = True

    elif series == 'Series_4':
        if iteration == 0:
            ######
            ##sort modelmaps across groups by mean_models_milz_etal_2015
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_4', 'modelmaps_across_groups.hdf')
            sortbyhdf5 = os.path.join(sortbyfolder, 'mean_models_milz_etal_2015.asc')
            outputhdf5 = os.path.join( outputfolder, 'Series_4', 'modelmaps_across_groups{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortAllByNormDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        elif iteration == 1:
            ######
            ##sort modelmaps_across_conds by modelmaps_across groups
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_4', 'modelmaps_across_conds.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_4', 'modelmaps_across_groups_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_4', 'modelmaps_across_conds{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortGroupByAllDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        elif iteration == 2:
            ######
            ##sort modelmaps_across_pts by modelmaps_across_conds_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_4', 'modelmaps_across_pts.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_4', 'modelmaps_across_conds_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_4', 'modelmaps_across_pts{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortPtByGroupDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        elif iteration == 3:
            ######
            ##sort modelmaps_across_runs by modelmaps_across_pts_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_4', 'modelmaps_across_runs.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_4', 'modelmaps_across_pts_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_4', 'modelmaps_across_runs{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortCondByPtDataProvider2(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        elif iteration == 4:
            ######
            ##means across groups
            ######
            inputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_4', 'modelmaps_across_runs_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')
            inputdataset = first_input
            sortbydataset = 'modelmap'
            outputdataset = '{0}{1}' .format(first_input, '_Series_4_sorted')
            sortdata_provider = SortRunByCondDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        else:
            stop = True

    elif series == 'Series_5':
        if iteration == 0:
            ######
            ##sort modelmaps across pts by mean_models_milz_etal_2015
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_5', 'modelmaps_across_pts.hdf')
            sortbyhdf5 = os.path.join(sortbyfolder, 'mean_models_milz_etal_2015.asc')
            outputhdf5 = os.path.join( outputfolder, 'Series_5', 'modelmaps_across_pts{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortAllByNormDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        elif iteration == 1:
            ######
            ##sort modelmaps_across_groups by modelmaps_across_pts
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_5', 'modelmaps_across_groups.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_5', 'modelmaps_across_pts_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_5', 'modelmaps_across_groups{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortGroupByAllDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        elif iteration == 2:
            ######
            ##sort modelmaps_across_conds by modelmaps_across_groups_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_5', 'modelmaps_across_conds.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_5', 'modelmaps_across_groups_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_5', 'modelmaps_across_conds{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortPtByGroupDataProvider2(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)

        elif iteration == 3:
            ######
            ##sort modelmaps_across_runs by modelmaps_across_conds_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_5', 'modelmaps_across_runs.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_5', 'modelmaps_across_conds_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_5', 'modelmaps_across_runs{0}.hdf' .format('_sorted') )
            inputdataset = 'modelmap'
            sortbydataset = 'modelmap'
            outputdataset = 'modelmap'
            sortdata_provider = SortCondByPtDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)

        elif iteration == 4:
            ######
            ##means across groups
            ######
            inputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_5', 'modelmaps_across_runs_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'all_recordings.hdf')
            inputdataset = first_input
            sortbydataset = 'modelmap'
            outputdataset = '{0}{1}' .format(first_input, '_Series_5_sorted')
            sortdata_provider = SortRunByCondDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)
        else:
            stop = True
    else:
        print 'series not defined:', series


    #sortdata_provider = SortGroupsByAllDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset)

    return sortdata_provider, stop


def run_sort_maps(data_provider, confobj):
    for output_path in data_provider.get_outputs():
        input = data_provider.get_input_data(output_path)
        sortby = data_provider.get_sortby_data(output_path)
        output_data, output_attributes = sort_maps(confobj, input, sortby)     
        if not output_data == []:
            data_provider.write_output_data(output_path, output_data, output_attributes)


def run_sort_maps_series(series, inputfolder, sortbyfolder, outputfolder, first_input, confobj):
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
        sortdata_provider, stop = get_io_sortmap_for_series(series, iteration, inputfolder, sortbyfolder, outputfolder, first_input)
        if stop:
            break

        run_sort_maps(sortdata_provider, confobj)

        iteration = iteration + 1
