
##################################
#######  Import Packages  ########
##################################

from __future__ import print_function

import os.path
import math
import itertools
import operator

import numpy as np
from scipy.stats import pearsonr
from contextlib import closing

import h5py
from keypy.microstates.sortmaps_provider import *
####-------------------------------------####

##############################
########  Sort Maps  ########
##############################


def sortmaps(confobj, input, sortby, input_original):
    """
    This function performs the sorting of a given input based on a given sortby. The sorted input_original is returned.
 
    Parameters
    ----------
    confobj : object of type MstConfiguration
          Contains the parameters used for microstate modelmap computation and visualization. 
    input : 
        Maps to be sorted reduced to the same number in the same order as sortby.
    sortby : 
        Maps to sort by reduced to the same number in the same order as input.
    input_original : 
        Maps to be sorted (original: without channel number reduction or reordering).

    Return Values
    ----------
    newraw :
        input_original ordered based on sortby_maps
    map_corr_list : 
        Correlations between each map with the map it was labeled by
    """

    modelmaps = input
    modelmaps_original = input_original
    sortby_maps = sortby
    original_nr_of_maps = confobj.original_nr_of_maps
    ERP = confobj.ERP

    ######
    ###Get best attribution matrix of my maps ordered by TK maps / model maps sorted of any stage
    ######

    attribution_matrix=[]
    mean_correlations = dict.fromkeys( list(range(original_nr_of_maps)) )

    for ithperm, perm in enumerate(itertools.permutations((list(range(original_nr_of_maps))))):    
        pearsons=[]
        pearsons2=[]
        #sorting is based on maximal r-value
        #for i in range(original_nr_of_maps):
        for i in range(len(sortby_maps[:,0])):
            pr, pp=pearsonr( modelmaps[i,:], sortby_maps[perm[i],:])

            pearsons.append(pr)
            if ERP:
                #ERP-based data need to distinguish between pos and neg correlations
                pearsons2.append(float(pearsons[i]))            
            else:
                pearsons2.append([abs(float(pearsons[i]))])

        mean_correlationo =np.mean(pearsons2)
        mean_correlations[ithperm] = mean_correlationo
        

    bestpermi=max(mean_correlations.iteritems(), key=operator.itemgetter(1))[0]

    if confobj.debug:
        print('bestpermi', bestpermi, 'mean_correlations[bestpermi]', mean_correlations[bestpermi])

    bestpermi_corr = mean_correlations[bestpermi]                       
    attribution_matrix = list(itertools.permutations((list(range(original_nr_of_maps)))))[bestpermi]

    if confobj.debug:
        print('attribution_matrix', attribution_matrix)

    ######
    ###Check if inversion is necessary and save whole EEG into newraw
    ######

    #list that saves the best correlations for the 4 maps
    map_corr_list = []

    newraw=np.zeros((modelmaps_original.shape))
    
    for mapi in range(len(sortby_maps[:,0])):
        pr, pp=pearsonr( modelmaps[mapi,:], sortby_maps[attribution_matrix[mapi],:])

        if pr < 0:
            if confobj.debug:
                print('r=', pr, 'modelmap', mapi, 'reversed')
            newraw[attribution_matrix[mapi],:]=modelmaps_original[mapi,:]*-1
            map_corr_list.append(abs(pr))

        else:
            if confobj.debug:
                print('r=', pr,'modelmap', mapi, 'not reversed')
            newraw[attribution_matrix[mapi],:]=modelmaps_original[mapi,:]
            map_corr_list.append(pr)

    attributes={}
    attributes['map_corr_list']=map_corr_list
    return newraw, attributes

####--------------------------------------------------------------------------####
####--------------------------------------------------------------------------####
####--------------------------------------------------------------------------####

def get_io_sortmap_for_series(series, iteration, inputfolder, inputfile_name, sortbyfolder, sortbyfile, sortbyfile_chlist, outputfolder, first_input):
    """
    Gets inputs and outputs of modelmaps for the series of modelmap computations specified.
 
    Parameters
    ----------
    series : {Series_1, Series_2, ...}
		Type of modelmap computation, e.g. 'Series_1' (means across runs for each group pt cond)
    iteration : int
		number of iterations for which different seeds are used
    inputfolder :
		folder of input hdf5
    inputfile_name :
        name of hdf5 file in inputfolder
    sortbyfolder :
        folder of sortby hdf5 / external file
    sortbyfile :
        filename of file to sortby
    sortbyfile_chlist :
        channellist of file to sortby
    outputfolder:
		folder of output hd5
    first_input:
		name of dataset of input hdf5 for first model map computation (e.g. 'modelmap' or 'meanmod')

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
            ##sort modelmaps across groups by sortbyfile
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_1', 'meanmods_across_groups.hdf')
            sortbyhdf5 = os.path.join(sortbyfolder, sortbyfile)
            outputhdf5 = os.path.join( outputfolder, 'Series_1', 'meanmods_across_groups{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortAllByNormDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        elif iteration == 1:
            ######
            ##sort meanmods_across_pts by meanmods_across_groups_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_1', 'meanmods_across_pts.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_1', 'meanmods_across_groups_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_1', 'meanmods_across_pts{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortGroupByAllDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        elif iteration == 2:
            ######
            ##sort meanmods_across_conds by meanmods_across_pts_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_1', 'meanmods_across_conds.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_1', 'meanmods_across_pts_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_1', 'meanmods_across_conds{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortPtByGroupDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        elif iteration == 3:
            ######
            ##sort sort meanmods_across_runs by meanmods_across_conds_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_1', 'meanmods_across_runs.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_1', 'meanmods_across_conds_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_1', 'meanmods_across_runs{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortCondByPtDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        elif iteration == 4:
            ######
            ##sort modelmaps by meanmods_across_runs_sorted
            ######
            inputhdf5 = os.path.join( inputfolder, inputfile_name)
            sortbyhdf5 = os.path.join( outputfolder, 'Series_1', 'meanmods_across_runs_sorted.hdf')
            outputhdf5 = os.path.join( inputfolder, inputfile_name)
            inputdataset = first_input
            sortbydataset = 'meanmod'
            outputdataset = '{0}{1}' .format(first_input, '_Series_1_sorted')
            sortdata_provider = SortRunByCondDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        else:
            stop = True


    elif series == 'Series_2':
        if iteration == 0:
            ######
            ##sort modelmaps across runs by sortbyfile
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_2', 'meanmods_across_runs.hdf')
            sortbyhdf5 = os.path.join(sortbyfolder, sortbyfile)
            outputhdf5 = os.path.join( outputfolder, 'Series_2', 'meanmods_across_runs{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortAllByNormDataProvider2(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        elif iteration == 1:
            ######
            ##sort meanmods_across_pts by meanmods_across_runs
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_2', 'meanmods_across_pts.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_2', 'meanmods_across_runs_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_2', 'meanmods_across_pts{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortCondByPtDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        elif iteration == 2:
            ######
            ##sort modelmaps by meanmods_across_pts_sorted
            ######
            inputhdf5 = os.path.join( inputfolder, inputfile_name)
            sortbyhdf5 = os.path.join( outputfolder, 'Series_2', 'meanmods_across_pts_sorted.hdf')
            outputhdf5 = os.path.join( inputfolder, inputfile_name)
            inputdataset = first_input
            sortbydataset = 'meanmod'
            outputdataset = '{0}{1}' .format(first_input, '_Series_2_sorted')
            sortdata_provider = SortRunByCondDataProvider2(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        else:
            stop = True

    elif series == 'Series_3':
        if iteration == 0:
            ######
            ##sort modelmaps across conds by sortbyfile
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_3', 'meanmods_across_conds.hdf')
            sortbyhdf5 = os.path.join(sortbyfolder, sortbyfile)
            outputhdf5 = os.path.join( outputfolder, 'Series_3', 'meanmods_across_conds{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortAllByNormDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        elif iteration == 1:
            ######
            ##sort meanmods_across_groups by meanmods_across_conds
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_3', 'meanmods_across_groups.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_3', 'meanmods_across_conds_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_3', 'meanmods_across_groups{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortGroupByAllDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        elif iteration == 2:
            ######
            ##sort meanmods_across_pts by meanmods_across_groups_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_3', 'meanmods_across_pts.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_3', 'meanmods_across_groups_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_3', 'meanmods_across_pts{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortPtByGroupDataProvider2(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        elif iteration == 3:
            ######
            ##sort sort meanmods_across_runs by meanmods_across_pts_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_3', 'meanmods_across_runs.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_3', 'meanmods_across_pts_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_3', 'meanmods_across_runs{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortCondByPtDataProvider2(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        elif iteration == 4:
            ######
            ##sort modelmaps by meanmods_across_runs_sorted
            ######
            inputhdf5 = os.path.join( inputfolder, inputfile_name)
            sortbyhdf5 = os.path.join( outputfolder, 'Series_3', 'meanmods_across_runs_sorted.hdf')
            outputhdf5 = os.path.join( inputfolder, inputfile_name)
            inputdataset = first_input
            sortbydataset = 'meanmod'
            outputdataset = '{0}{1}' .format(first_input, '_Series_3_sorted')
            sortdata_provider = SortRunByCondDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        else:
            stop = True

    elif series == 'Series_4':
        if iteration == 0:
            ######
            ##sort modelmaps across groups by sortbyfile
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_4', 'meanmods_across_groups.hdf')
            sortbyhdf5 = os.path.join(sortbyfolder, sortbyfile)
            outputhdf5 = os.path.join( outputfolder, 'Series_4', 'meanmods_across_groups{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortAllByNormDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        elif iteration == 1:
            ######
            ##sort meanmods_across_conds by meanmods_across groups
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_4', 'meanmods_across_conds.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_4', 'meanmods_across_groups_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_4', 'meanmods_across_conds{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortGroupByAllDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        elif iteration == 2:
            ######
            ##sort meanmods_across_pts by meanmods_across_conds_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_4', 'meanmods_across_pts.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_4', 'meanmods_across_conds_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_4', 'meanmods_across_pts{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortPtByGroupDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        elif iteration == 3:
            ######
            ##sort meanmods_across_runs by meanmods_across_pts_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_4', 'meanmods_across_runs.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_4', 'meanmods_across_pts_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_4', 'meanmods_across_runs{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortCondByPtDataProvider2(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        elif iteration == 4:
            ######
            ##means across groups
            ######
            inputhdf5 = os.path.join( inputfolder, inputfile_name)
            sortbyhdf5 = os.path.join( outputfolder, 'Series_4', 'meanmods_across_runs_sorted.hdf')
            outputhdf5 = os.path.join( inputfolder, inputfile_name)
            inputdataset = first_input
            sortbydataset = 'meanmod'
            outputdataset = '{0}{1}' .format(first_input, '_Series_4_sorted')
            sortdata_provider = SortRunByCondDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        else:
            stop = True

    elif series == 'Series_5':
        if iteration == 0:
            ######
            ##sort modelmaps across pts by sortbyfile
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_5', 'meanmods_across_pts.hdf')
            sortbyhdf5 = os.path.join(sortbyfolder, sortbyfile)
            outputhdf5 = os.path.join( outputfolder, 'Series_5', 'meanmods_across_pts{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortAllByNormDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        elif iteration == 1:
            ######
            ##sort meanmods_across_groups by meanmods_across_pts
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_5', 'meanmods_across_groups.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_5', 'meanmods_across_pts_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_5', 'meanmods_across_groups{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortGroupByAllDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        elif iteration == 2:
            ######
            ##sort meanmods_across_conds by meanmods_across_groups_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_5', 'meanmods_across_conds.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_5', 'meanmods_across_groups_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_5', 'meanmods_across_conds{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortPtByGroupDataProvider2(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)

        elif iteration == 3:
            ######
            ##sort meanmods_across_runs by meanmods_across_conds_sorted
            ######
            inputhdf5 = os.path.join( outputfolder, 'Series_5', 'meanmods_across_runs.hdf')
            sortbyhdf5 = os.path.join( outputfolder, 'Series_5', 'meanmods_across_conds_sorted.hdf')
            outputhdf5 = os.path.join( outputfolder, 'Series_5', 'meanmods_across_runs{0}.hdf' .format('_sorted') )
            inputdataset = 'meanmod'
            sortbydataset = 'meanmod'
            outputdataset = 'meanmod'
            sortdata_provider = SortCondByPtDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)

        elif iteration == 4:
            ######
            ##means across groups
            ######
            inputhdf5 = os.path.join( inputfolder, inputfile_name)
            sortbyhdf5 = os.path.join( outputfolder, 'Series_5', 'meanmods_across_runs_sorted.hdf')
            outputhdf5 = os.path.join( inputfolder, inputfile_name)
            inputdataset = first_input
            sortbydataset = 'meanmod'
            outputdataset = '{0}{1}' .format(first_input, '_Series_5_sorted')
            sortdata_provider = SortRunByCondDataProvider1(inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset, outputdataset, sortbyfile_chlist)
        else:
            stop = True
    else:
        print('series not defined:', series)

    return sortdata_provider, stop

####--------------------------------------------------------------------------####
####--------------------------------------------------------------------------####
####--------------------------------------------------------------------------####


def run_sortmaps(data_provider, confobj, eeg_info_study_obj):
    for output_path in data_provider.get_outputs():
        print('Sorting maps: ', output_path.level0)

        input_reduced, input_original = data_provider.get_input_data(output_path, eeg_info_study_obj.chlist)
        sortby = data_provider.get_sortby_data(output_path)
        output_data, output_attributes = sortmaps(confobj, input_reduced, sortby, input_original)     
        if not output_data == []:
            data_provider.write_output_data(output_path, output_data, output_attributes)

####--------------------------------------------------------------------------####


def run_sortmaps_series(series, inputfolder, inputfile_name, sortbyfolder, sortbyfile, sortbyfile_chlist, outputfolder, first_input, confobj, eeg_info_study_obj) :
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
    inputfile_name : str
        hdf5 file name in input folder
    sortbyfolder :
        folder of sortby hdf5 / external file
    sortbyfile :
        filename of file to sortby
    sortbyfile_chlist :
        channellist of file to sortby
    outputfolder : str
        path to folder that will contain the output hdf5 file
    first_input : str
        dataset name that contains the N modelmaps for each run that the meanmod computation should be based on, e.g. 'modelmap'.
    confobj : object of type MstConfiguration
         Contains the parameters used for microstate modelmap computation and visualization. 
    """
    print('Computing Sortmaps ....')
    if confobj.original_nr_of_maps != 4:
        print('Warning: Make sure that the model maps you want to sortby are identical in number to the model maps you computed!')

    stop = False
    iteration = 0
    while True:
        sortdata_provider, stop = get_io_sortmap_for_series(series, iteration, inputfolder, inputfile_name, sortbyfolder, sortbyfile, sortbyfile_chlist, outputfolder, first_input)
        if stop:
            break
        print('computing sortmaps iteration ', iteration)
        run_sortmaps(sortdata_provider, confobj, eeg_info_study_obj)

        iteration = iteration + 1

####--------------------------------------------------------------------------####