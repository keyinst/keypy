##################################
#######  run_mstate_paramters  ########
##################################

from scipy.stats import nanmean, pearsonr
from keypy.microstates.data_provider import *
from keypy.microstates.microstates_helper import *
from numpy import sqrt

##########################
#######  Classes  ########
##########################


####Abstract Class SortDataProvider

class ParametersDataProvider(object):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset):
        self._file = inputhdf5
        self._sortbyfile = sortbyhdf5
        self._outputfile = outputhdf5
        self._inputdataset = inputdataset
        self._sortbydataset = sortbydataset

####Sub Class CondDataProvider --> umbenennen

###Computer parameters by Norm Maps Data Provider 1

class ParametersByNormDataProvider1(ParametersDataProvider):
    def __init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset):
        ParametersDataProvider.__init__(self, inputhdf5, sortbyhdf5, outputhdf5, inputdataset, sortbydataset)

    #call it once to get a list of objects which contain the paths needed each to create one output
    def get_outputs(self):
        out_paths_set = Set()
        with closing( h5py.File(self._file, 'r') ) as f:
            for current_level0 in f['/'].keys():                 
                for current_level1 in f['/{0}'.format(current_level0)].keys():  
                    for current_level2 in f['/{0}/{1}'.format(current_level0, current_level1)].keys():  
                        for current_level3 in f['/{0}/{1}/{2}'.format(current_level0, current_level1, current_level2)].keys():  
                            out_paths_set.add(Levels4Path(current_level0, current_level1, current_level2, current_level3))
        return list(out_paths_set)
    
    ### group pt cond run
    #ein Aufruf pro Output, gets the particular EEG file which is needed for parameter computation (for a particular group, pt, cond, run)
    def get_input_data(self, output_path):
        with closing( h5py.File(self._file, 'r') ) as g:            
            path = '/{0}/{1}/{2}/{3}' .format(output_path.level0, output_path.level1, output_path.level2, output_path.level3)
            eeg_value = g['/{0}/{1}' .format(path, self._inputdataset)] 
            if all(eeg_value[0,:] == 0):
                print 'Error!', path, 'has all zeros', 'group, pt, cond ignored.'    
            else:
                eeg=eeg_value[:]
        return eeg

    ### group pt cond
    #ein Aufruf pro Output, gets a list of all modelmaps which are to be sorted by
    def get_sortby_data(self, output_path):
        microstate_run_value = np.loadtxt(self._sortbyfile)

        if all(microstate_run_value[0,:] == 0):
            print 'Error!', sortbyhdf5, 'has all zeros', 'group, pt, cond ignored.'    
        else:
            model_map=microstate_run_value[:]
        return model_map
       
    #writes output into new hdf5 at correct location
    def write_output_data(self, output_path, output_data, output_attributes):
        with closing( h5py.File(self._outputfile, 'a') ) as k:
            print 'output_paths used for output', output_path.level0

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


            ###Writing Output for data_attributes and data

            runwise_data, epochwise_data, mapwise_data = output_data
            number_of_epochs = len(mapwise_data[dur_state_all_epochs])

            ###Add runwise datasets

            ##Add State Match Mean percentage

            #create dataset of mean percentages for class correspondances across epochs
            map_avg_perc=parameter_preprocessing(confobj, epochwise_data['State Match Mean percentage'])  

            try:
                run_group.create_dataset('State Match Mean p', data = map_avg_perc)
            except:
                run_group.create_dataset('State Match Mean p', data = 999)


            ##add other runwise datasets
            for dataset_name in runwise_data:
                try:
                    run_group.create_dataset(dataset_name, data = runwise_data[dataset_name])
                except:
                    run_group.create_dataset(dataset_name, data = 999)


            ###save individual maps (not sorted)
            '''
            #yet to integrate option to save individual maps to txt files
            outputfolder_individumaps = 'E:\\Programming\\Python\\hdf5_files_140411\\individual_mm_states'
            np.savetxt(op.join(outputfolder_individumaps, 'individual_mstate_%s_%s.txt' % (looper2, looper1)), individu_mstate, delimiter=' ')
            '''

            ##epochwise output
            for epochnr in range(number_of_epochs):
                #create folder
                ep_group = run_group.create_group("ep_"+"%03d" % (epochnr,))  
                if confobj.debug:
                    print epochnr

                #add attributes to epoch folder
                for attribute_name, attribute_value in data_attributes.iteritems():
                    ep.group.attrs['{0}' .format(key)]=str(value[epochnr])
    
                #add epochwise datasets
                for epochwise_data_name in epochwise_data:
                    try:
                        ep_group.create_dataset(epochwise_data_name, data = epochwise_data[epochwise_data_name][epochnr])
                    except:
                        ep_group.create_dataset(epochwise_data_name, data = 999)

                #add mapwise datasets

                for mapnr in range(confobj.original_nr_of_maps):
                    map_group = ep_group.create_group( "map_"+"%02d" % (mapnr,))

                    for mapwise_data_name in mapwise_data:
                        try:
                            map_group.create_dataset(mapwise_data_name, data = mapwise_data[mapwise_data_name][epochnr][mapnr])
                        except:                             
                            map_group.create_dataset(mapwise_data_name, data = 0)




def parameter_preprocessing(confobj, state_match_percentage_all_epochs):
     map_avg_perc = np.zeros((confobj.original_nr_of_maps+1, 2))

     for map in range(confobj.original_nr_of_maps+1):
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

#maps=labelbymaps
#eeg=eeg_re

def compute_mstate_parameters(confobj, eeg, maps, eeg_info_study_obj):
    #Configuration
    #method_GFPpeak = confobj.method_GFPpeak
    #debug  = confobj.debug
    #use_gfp_peaks = confobj.use_gfp_peaks
    #use_smoothing = confobj.use_smoothing
    #use_fancy_peaks = confobj.use_fancy_peaks
    #original_nr_of_maps = confobj.original_nr_of_maps
    TF = eeg_info_study_obj.tf
    Fs = eeg_info_study_obj.sf

    #####Create Dictionaries for parameters across all epochs
    dur_state_all_epochs = dict.fromkeys(range(len(eeg)/TF))
    freq_dict_all_epochs = dict.fromkeys(range(len(eeg)/TF))
    dur_dict_all_epochs = dict.fromkeys(range(len(eeg)/TF))
    cov_dict_all_epochs = dict.fromkeys(range(len(eeg)/TF))
    gfp_peak_nr_all_epochs = dict.fromkeys(range(len(eeg)/TF))
    mean_gfp_all_epochs = dict.fromkeys(range(len(eeg)/TF))
    gfp_mean_all_epochs = dict.fromkeys(range(len(eeg)/TF))
    durstd_dict_all_epochs= dict.fromkeys(range(len(eeg)/TF))
    start_state_list_all_epochs= dict.fromkeys(range(len(eeg)/TF))
    exp_var_all_epochs= dict.fromkeys(range(len(eeg)/TF))
    exp_var_tot_all_epochs= dict.fromkeys(range(len(eeg)/TF))
    state_match_percentage_all_epochs= dict.fromkeys(range(len(eeg)/TF))
    state_match_percentage_std_all_epochs= dict.fromkeys(range(len(eeg)/TF))

    #################  
    # 3.) LOOP ACROSS 2 Sec Segments
    #################

    individu_dict = dict.fromkeys(range(confobj.original_nr_of_maps))
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
            print 'GFP Curve computed'
            
            
        #################   
        #4.) Compute GFP Peaks (identical to microstates.py)
        #################
        
        gfp_peak_indices, gfp_curve = compute_gfp_peaks(gfp_curve, confobj.use_gfp_peaks, confobj.use_smoothing, confobj.gfp_type_smoothing, confobj.smoothing_window, confobj.use_fancy_peaks)
       
        #################   COPIED AND ADAPTED FROM MICROSTATE_SCRIPT_OFFICIAL4.py
        #5.) Determine Mstate class for each peak
        #################

        #Compare each gfp peak index (gfp_peak_indices) to the 4 mstate maps (maps) --> get attribution matrix for the eeg 0,1,3,2,0,2 etc.

        #initialize attribution matrix (one with the correlations and one only with the highest gfp_peak_index)
        attribution_matrix= dict.fromkeys(gfp_peak_indices)
        attribution_matrix_indices= dict.fromkeys(gfp_peak_indices)

        #get index of mstate map that correlates highest with gfp_peak_index

        for key in attribution_matrix.keys():
            tf=epoch[key]
            #may have been wrong with eeg instead of epoch?
            #tf=eeg[key]
            corr_list = []
            for mapnr in range(len(maps)):
                pr, pp=pearsonr( maps[mapnr,:], tf)
                if confobj.ERP:
                    corr_list.append(pr)
                else:
                    corr_list.append(abs(pr))
            attribution_matrix[key]=corr_list

            #confobj.correspondance_cutoff = 0 for no cutoff
            if max(corr_list) > confobj.correspondance_cutoff:
                attribution_matrix_indices[key]=corr_list.index(max(corr_list))
                attribution_matrix_indices[key]=[corr_list.index(max(corr_list)), max(corr_list)]
            else:
                attribution_matrix_indices[key]=[999, max(corr_list)]   

        #tf_begin is determined as, the first gfp peak where the mstate changes from one to another, e.g. A A A (B) B B (A) (B) B
        #sort dict
        from collections import OrderedDict
        attribution_matrix_indices_sorted = OrderedDict(sorted(attribution_matrix_indices.items(), key=lambda x: x[1]))

        start_state_list =[]
        previous=sorted(attribution_matrix_indices_sorted)[0]
        start = sorted(attribution_matrix_indices_sorted)[0]
        for ele in sorted(attribution_matrix_indices_sorted):
            #print 'value of ele, previous', ele, attribution_matrix_indices[ele], attribution_matrix_indices[previous]
            #print 'ele, previous', ele, previous

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
        dur_state=dict.fromkeys(range(confobj.original_nr_of_maps))
        gfp_state=dict.fromkeys(range(confobj.original_nr_of_maps))
                  
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
                    print outti, inni, dur_state

        

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

        freq_dict=dict.fromkeys(range(confobj.original_nr_of_maps))
        for mapnr in range(confobj.original_nr_of_maps):
            if dur_state[mapnr] == None:
                print 'epochnr, mapnr no content', epochnr, mapnr
                freq_dict[mapnr] = 0
            else:
                #freq corrected by smaller epoch size (returns # of occurrences per second based on info from whole epoch)
                freq_dict[mapnr] = (len(dur_state[mapnr])*(TF*(1000/float(Fs)/dur_sum)/float(TF/float(Fs))))


       
        ####
        #Mean Duration & Std from tf to secs
        ####

        dur_dict=dict.fromkeys(range(confobj.original_nr_of_maps))
        durstd_dict=dict.fromkeys(range(confobj.original_nr_of_maps))
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
        cov_dict=dict.fromkeys(range(confobj.original_nr_of_maps))
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
        gfp_mean=dict.fromkeys(range(confobj.original_nr_of_maps))

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
        #Explained variance (Problem: wir haben nur noch 19 Elektroden for TK sort!?)
        ####     
        
        #unsicher, ob es noch normiert werden muss
        model = maps

        #TK Normierung
        #Berechnung des Normierungs-Vektors, es geht darum alles auf Vektorlaenge 1 zu setzen
        b=np.sum(np.abs(model)**2,axis=-1)**(1./2)    
        #Teilung aller Elemente durch Normierungs-Vektor
        #nch ersetzen mit model.shape[1]
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

        state_match_percentage=dict.fromkeys(range(confobj.original_nr_of_maps+1))
        state_match_percentage_std=dict.fromkeys(range(confobj.original_nr_of_maps+1))

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

    ###Convert all measures into two dictionaries depending on whether they represent a dataset or attribute
 
    #create dictionary for runwise datasets (keys are the descriptions that will be in the hdf5 outputfile) 
    runwise_data = {}
    runwise_data['Individual_States']=individu_mstate
       
    #create dictionary for epochwise datasets (keys are the descriptions that will be in the hdf5 outputfile) 
    epochwise_datasets = [start_state_list_all_epochs, state_match_percentage_all_epochs, state_match_percentage_std_all_epochs]
    epochwise_datasets_string = ['Start state array', 'State Match Mean percentage', 'State Match Std percentage']
    
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



'''old outputs

return dur_state_all_epochs, freq_dict_all_epochs, dur_dict_all_epochs, cov_dict_all_epochs, gfp_peak_nr_all_epochs, mean_gfp_all_epochs, gfp_mean_all_epochs, durstd_dict_all_epochs, start_state_list_all_epochs, exp_var_all_epochs, exp_var_tot_all_epochs, state_match_percentage_all_epochs, state_match_percentage_std_all_epochs, individu_mstate

    dur_state_all_epochs, occ, dur, cov, gfp_peak_nr_all_epochs, mean_gfp_all_epochs, gfp_mean_all_epochs, durstd_dict_all_epochs, start_state_list_all_epochs, exp_var_all_epochs, exp_var_tot_all_epochs, state_match_percentage_all_epochs, state_match_percentage_std_all_epochs, individu_mstate=mstate_parameters

'''


def run_parameters(data_provider, confobj, eeg_info_study_obj):
    for output_path in data_provider.get_outputs():
        input = data_provider.get_input_data(output_path)
        sortby = data_provider.get_sortby_data(output_path)
        #compute_mstate_parameters demands that input and sortby are in the same data format (equal nch)
        output_data, output_attributes = compute_mstate_parameters(confobj, input, sortby, eeg_info_study_obj)
        if not output_data == []:
            data_provider.write_output_data(output_path, output_data, output_attributes)
