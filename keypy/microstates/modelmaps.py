# -*- coding: utf-8 -*-

##################################
#######  Import Packages  ########
##################################

from __future__ import print_function

import os.path
import math
import itertools
import operator

import numpy
import scipy.stats

from numpy import linalg as LA

import random
from random import randrange

from keypy.microstates.modelmaps_provider import *
from keypy.microstates.microstates_helper import princomp_B, compute_gfp
####--------------------------------------------------------------------------####


##################################
#######  find_model_maps  ########
##################################


def get_io_modelmap_for_series(series, iteration, inputfolder, hdf5_filename, outputfolder, first_input):
    """
    Gets inputs and outputs of modelmaps for the series of modelmap computations specified.
 
    Parameters
    ----------
    series : {Series_1, Series_2, ...}
		Type of modelmap computation, e.g. 'Series_1' (means across runs for each group pt cond)
    iteration : int
		number of iterations for which different seeds are used
    inputfolder : path
		folder of input hdf5
    hdf5_filename : str
        filename of the input hdf5 file
    outputfolder: path
		folder of output hd5
    first_input: str
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
            inputhdf5 = os.path.join( inputfolder, hdf5_filename)
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
            inputhdf5 = os.path.join( inputfolder, hdf5_filename)
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
            ##means across runs for each group pt cond
            ######
            inputhdf5 = os.path.join( inputfolder, hdf5_filename)
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_runs.hdf')
            modelmap_input = first_input
            modelmap_output = 'modelmap'
            computation_version ='means across runs for each group pt cond'
        elif iteration == 1:
            ######
            ##means across pts for each group cond
            ######
            inputhdf5 = os.path.join( outputfolder, 'modelmaps_across_runs.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_pts.hdf')
            modelmap_input = 'modelmap'
            modelmap_output = 'modelmap'
            computation_version ='means across pts for each group cond'
        elif iteration == 2:
            ######
            ##means across groups for each cond
            ######
            inputhdf5 = os.path.join( outputfolder, 'modelmaps_across_pts.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_groups.hdf')
            modelmap_input = 'modelmap'
            modelmap_output = 'modelmap'
            computation_version ='means across groups for each cond'
        elif iteration == 3:
            ######
            ##means across conds
            ######
            inputhdf5 = os.path.join( outputfolder, 'modelmaps_across_groups.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_conds.hdf')
            modelmap_input = 'modelmap'
            modelmap_output = 'modelmap'
            computation_version ='means across conds'
        else:
            stop = True

    elif series == 'Series_4':
        if iteration == 0:
            ######
            ##means across runs for each group pt cond
            ######
            inputhdf5 = os.path.join( inputfolder, hdf5_filename)
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_runs.hdf')
            modelmap_input = first_input
            modelmap_output = 'modelmap'
            computation_version ='means across runs for each group pt cond'
        elif iteration == 1:
            ######
            ##means across pts for each group cond
            ######
            inputhdf5 = os.path.join( outputfolder, 'modelmaps_across_runs.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_pts.hdf')
            modelmap_input = 'modelmap'
            modelmap_output = 'modelmap'
            computation_version ='means across pts for each group cond'
        elif iteration == 2:
            ######
            ##means across pts for each group
            ######
            inputhdf5 = os.path.join( outputfolder, 'modelmaps_across_pts.hdf')
            outputhdf5 = os.path.join( outputfolder, 'modelmaps_across_conds.hdf')
            modelmap_input = 'modelmap'
            modelmap_output = 'modelmap'
            computation_version ='means across conds for each group'
        elif iteration == 3:
            ######
            ##means across groups
            ######
            inputhdf5 = os.path.join( outputfolder, 'modelmaps_across_conds.hdf')
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
            inputhdf5 = os.path.join( inputfolder, hdf5_filename)
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
        print('series not defined:', series)

    return inputhdf5, outputhdf5, modelmap_input, modelmap_output, computation_version, stop

####-------------------------------------------------------------------------------------------####
####-------------------------------------------------------------------------------------------####
####-------------------------------------------------------------------------------------------####

##################################
#######  find_model_maps  ########
##################################

def find_model_maps(confobj, model_maps_foundation):
    """
    Finds best modelmaps for input maps.
 
    Parameters
    ----------
    confobj : object of type MstConfiguration
         Contains the parameters used for microstate computation and visualization. 
    model_maps_foundation : array
        Contains the input maps for the modelmap computation.

    Returns
    ----------
	bestresmm : array
		best modelmaps for the input maps
	attributes_dict : dict
		includes ['Mean Correlation'] = best_mean_correlation
    """

    #get number of original_nr_of_maps*nch modelmaps in model_maps_foundation list
    number_of_basic_maps = len(model_maps_foundation)

    #get configuration info from confobj
    original_nr_of_maps = confobj.original_nr_of_maps
    seed_number=confobj.seed_number
    max_number_of_iterations=confobj.max_number_of_iterations

    #get number of channel from first modelmap in model_maps_foundation list
    nch = model_maps_foundation[0].shape[1]

    #initialize variables for modelmap computation
    iii=0
    best_results = {}
    bestcorr = {}
    
    best_mean_correlation = 0
    found_best = False
 
    ##get fixed seed from confobj if not None
    fixed_seed = confobj.fixed_seed
    if fixed_seed != None:
        numpy.random.seed(fixed_seed)
    
    ##get user-requested normalization for model_maps_foundation
    for elenr, ele in enumerate(model_maps_foundation):
        model_maps_foundation[elenr]=normalize_maps(ele, confobj.modelmaps_normalization_type)    

    #if there is only one map in the model_maps_foundation list, skip procedure and return output directly
    if len(model_maps_foundation) == 1:
        bestresmm =  model_maps_foundation[0]
        attributes_dict = {}
        attributes_dict['Mean Correlation'] = 1


    else:           
        while iii < seed_number:
        
            print('------------')
            print('SEED', iii, ':')
        
            delta_correlation = 0.1
            mean_correlation=0.1
            ii = 0
            delta_attribution_matrix=1
    
            # create an empty list for iteration iii and add the best found maps
            best_results[iii] = {}
                
            # intitialize model_maps_mean array
            model_maps_mean=numpy.zeros( (original_nr_of_maps, nch) )

            # intitialize attribution_matrix_old
            attribution_matrix_old=numpy.zeros( (number_of_basic_maps,original_nr_of_maps) )
              
    
            #########
            ### create original_nr_of_maps x nch array with 4 random maps  
            #########   
            randmap=numpy.zeros((original_nr_of_maps, nch))
            #selects from array: model_maps_foundation, the EEG of original_nr_of_maps random keys and random maps

            if iii <= 10:
                #Take all maps of one participant as seed
                random_vp= random.choice(list(range(len(model_maps_foundation))))
                for i in range(original_nr_of_maps):
                    randmap[i,:]=model_maps_foundation[random_vp][i,:]

                    print('seed', iii, 'random_vp', random_vp, 'random_map', i)
            else:
                #Completely random
                for i in range(original_nr_of_maps):
                    random_vp= random.choice(list(range(len(model_maps_foundation))))
                    random_map= randrange(original_nr_of_maps)
                    randmap[i,:]=model_maps_foundation[random_vp][random_map,:]

                    print('seed', iii, 'random_vp', random_vp, 'random_map', random_map)
        
            if number_of_basic_maps < original_nr_of_maps:
                print('Attention, you have only', number_of_basic_maps ,'participants/conditions/runs to select your', original_nr_of_maps ,'random maps from')


            #########
            ### Loop Across Iterations
            #########  
  
            while (ii < max_number_of_iterations) and abs(delta_correlation) > 0.0002:
                print('iteration:', ii)
                best_results[iii][ii] = []       
                #intitialize attribution_matrix
                #attribution_matrix=numpy.zeros( (number_of_basic_maps ,original_nr_of_maps) )
                attribution_matrix= dict.fromkeys(list(range(len(model_maps_foundation))))  
                        
                #Find best map attribution for VP, save its indices and the corresponding correlation
                #for ivpi, vpi in enumerate(VP): (SO WARS VORHER)
                for numberi, number in enumerate((list(range(len(model_maps_foundation))))):
                    #get maps from VP with key number
                    maps=model_maps_foundation[number]
                    #initialize dict for all permuations to later save mean correlation of that permutation
                    mean_correlations = dict.fromkeys( list(range(math.factorial(original_nr_of_maps))) )
                
                    #for ithperm, perm in enumerate(itertools.permutations([0,1,2,3])):
                    for ithperm, perm in enumerate(itertools.permutations((list(range(original_nr_of_maps))))):    
                        pearsons=[]
                        pearsons2=[]
                        for i in range(original_nr_of_maps):
                            pr, pp=scipy.stats.pearsonr( maps[i,:], randmap[perm[i],:])
                        
                            pearsons.append(pr)
                            if confobj.ERP:
                                pearsons2.append(float(pearsons[i]))            
                            else:
                                pearsons2.append([abs(float(pearsons[i]))])
                        
                        mean_correlationo =numpy.mean(pearsons2)
                    
                        mean_correlations[ithperm] = mean_correlationo
                
                    bestpermi=max(mean_correlations.iteritems(), key=operator.itemgetter(1))[0]
                
                    attribution_matrix[number] = list(itertools.permutations((list(range(original_nr_of_maps)))))[bestpermi]
        
                    #save bestcorrelation which corresponds to the one of the best attribution for VP
                    bestcorr[number]=mean_correlations[bestpermi]
                
                                
                #generate a dictionary, where keys 0, 1, 2, 3 and arrays nVP x nch (best map that corresponds for each vp)
                best_fit = dict.fromkeys( list(range( 0, original_nr_of_maps)) )

                for i in best_fit:
                    b=numpy.zeros(( number_of_basic_maps , nch))
                    for vpnr, vpi in enumerate(attribution_matrix.keys()):           
                        vpindex=list(attribution_matrix[vpi]).index(i)                   
                        b[vpnr,:]=model_maps_foundation[vpi][vpindex,:]
                    best_fit[i]=b
            
                # for each key you have an array where you have to compute the first PC across participants                      
                # extract the first principal component of the array across participants
                for k in best_fit:
                    P=best_fit[k]
                    if confobj.ERP:
                        coeff = P.mean(axis=0)  #ERP: mean computed instead of PC1
                    else:
                        coeff = princomp_B(P,1)
                    randmap[k,:] = coeff.ravel()
        
                ##get user-requested normalization for modelmaps
                randmap=normalize_maps(randmap, confobj.modelmaps_normalization_type)    

                #compute delta_correlation compared to last iteration     
                list_bestcorr_values= list(bestcorr.values())
                bc_v_arr = numpy.asarray(list_bestcorr_values)
                bc_v_arr_mean=bc_v_arr[~numpy.isnan(bc_v_arr)].mean()

                delta_correlation=bc_v_arr_mean-mean_correlation   

                #Update attribution matrix
                attribution_matrix_old=attribution_matrix
            
                #save attribution_matrix, mean_correlation and randmap IF delta_correlation >0
                if (delta_correlation >0):
                    #print 'delta correlation bigger than zero, append to best result'
                    #update mean_correlation of current iteration
                    mean_correlation=numpy.mean(bestcorr.values())
                    best_results[iii][ii]=[mean_correlation, randmap, attribution_matrix]     
                
                ii += 1

            iii += 1
           
                  
        #extract the best randmap and attribution_matrix for each seed       
        best_mean_correlation = 0.1

        for seednr in best_results.keys():
            for iterationnr in best_results[seednr].keys():

                if len(best_results[seednr][iterationnr]) == 0:
                    print('value is zero')
            
                else:
                    corr, resmm, attrmatrix = best_results[seednr][iterationnr]
                      
                    #print 'current best', best_mean_correlation, ' current corr', corr
                    if corr > best_mean_correlation :
                        found_best = True
                        best_mean_correlation = corr
                        bestresmm = resmm
                        bestattr_matrix = attrmatrix              

        if found_best:
            if confobj.debug:
                #print 'best of all resmm and attr matrix', bestresmm.shape, bestattr_matrix
                #print 'best of all correlations', best_mean_correlation
                #print 'for condition', ci
                pass
        else:
            if confobj.debug:
                print('for SEED', iii, 'not found best')
        #save attributes in dictionary
        attributes_dict = {}
        attributes_dict['Mean Correlation'] = best_mean_correlation

    return bestresmm, attributes_dict

####--------------------------------------------------------------------------####
####--------------------------------------------------------------------------####
####--------------------------------------------------------------------------####

##################################
#######  run_model_maps   ########
##################################
        
               
def run_model_maps(data_provider, find_model_maps, confobj):
    """
    This function performs calculations/transformations on a set of data. The algorithm
    uses the data provider to enumerate the desired output artifacts. For each output it lets
    the data provider enumerate the input data, calls the algorithm to run with the given input
    and writes the generated output through the data provider.

    This is an abstraction of both the data source (which could e.g. be an hdf5 file) and the different
    types of collecting input data and enumerating desired outputs.
 
    Parameters
    ----------
    data_provider : subclass of type DataProvider
         The data provider is responsible to generate a list of all outputs and to collect input data for each of the output elements
    find_model_maps : function
        The algorithm to run on the data.
    confobj : object of type MstConfiguration
        Contains the parameters used for microstate computation and visualization. 
    """

    for output_path in data_provider.get_outputs():
        #print 'output_paths in list', output_path.all
        inputs = data_provider.get_input_data(output_path)
        output_data, output_attributes = find_model_maps(confobj, inputs)
        if not output_data == []:
            data_provider.write_output_data(output_path, output_data, output_attributes)

####--------------------------------------------------------------------------####
####--------------------------------------------------------------------------####
####--------------------------------------------------------------------------####

##################################
#######  run_model_maps_series   ########
##################################

def run_model_maps_series(series, inputfolder, hdf5_filename, outputfolder, first_input, confobj):
    """
    Runs run_model_maps for each input of the series.
 
    Parameters
    ----------
    series : {'Series_1', .... 'Series_5'}
        Defines the order of modelmap computation based on the input microstates.
        'Series_1' : (1) means across runs for each group pt cond, (2) means across conds for each group pt, (3) means across pts for each group, (4) means across groups
        'Series_2' : (1) means across pts for each group cond run, (2) means across runs for each group cond
        'Series_3' : (1) means across runs for each group pt cond, (2) means across pts for each group cond, (3) means across groups for each cond, (4) means across conds
        'Series_4' : (1) means across runs for each group pt cond, (2) means across pts for each group cond, (3) means across conds for each group, (4) means across groups     
        'Series_5' : (1) means across runs for each group pt cond, (2) means across conds for each group pt , (3) means across groups for each pt, (4) means across pts
    inputfolder : str
        path to folder that contains the input hdf5 file
    hdf5_filename : str
        filename of the input hdf5 file
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
        inputhdf5, outputhdf5, modelmap_input, modelmap_output, computation_version, stop = get_io_modelmap_for_series(series, iteration, inputfolder, hdf5_filename, outputfolder, first_input)
        if stop:
            break
        selected_provider_class = get_data_provider_class(computation_version)
        data_provider = selected_provider_class(inputhdf5, outputhdf5, modelmap_input, modelmap_output)
        run_model_maps(data_provider, find_model_maps, confobj)

        iteration = iteration + 1


####--------------------------------------------------------------------------####
####--------------------------------------------------------------------------####
####--------------------------------------------------------------------------####
