##################################
#######  create_parameter_spss_sheets  ########
##################################

def run_mstate_paramters(confobj, eeg_info_mot_obj, sorted_maps_hdf, mstate_parameters_output_hdf, eeg_hdf, eeg_input, sort_type, sortby_type, spss_sheets=False, groups = False):     
    ###Loop across VPs, runs
    if groups:
        ###Create dictionary to save VP wise output
        if sort_type=='microstate_parameters_for_each_vp_for_each_run':
            mstate_parameters_all=dict.fromkeys(eeg_info_mot_obj.Group)
            for key in mstate_parameters_all.keys():
                mstate_parameters_all[key]=dict.fromkeys(eeg_info_mot_obj.VP)
                for vpli in  mstate_parameters_all[key].keys():
                    mstate_parameters_all[key][vpli]=dict.fromkeys(eeg_info_mot_obj.Cond)
        with closing( h5py.File(sorted_maps_hdf) ) as g:
            with closing( h5py.File(mstate_parameters_output_hdf, 'w') ) as h:
                for looper3 in eeg_info_mot_obj.Group:
                    if not looper3 in h['/'].keys():
                        looper3_group = h.create_group( looper3 )
                    else:
                        looper3_group = h['/%s' % (looper3)]  
                                         
                    for looper2 in eeg_info_mot_obj.VP:
                        with closing( h5py.File(eeg_hdf) ) as f:
                            try:
                                f['/%s/%s' % (looper3, looper2)]
                            except:
                                continue

                            if not looper2 in looper3_group.keys():
                                looper2_group = looper3_group.create_group( looper2 )
                            else:
                                looper2_group = h['/%s/%s' % (looper3, looper2)]    
                        
                            for looper1 in eeg_info_mot_obj.Cond:
                                try:
                                    timeframe_channel_dset = f['/%s/%s/%s/%s' % (looper3, looper2, looper1, eeg_input)]               
                                    path = '/%s/%s/%s/%s' % (looper3, looper2, looper1, eeg_input)

                                    if confobj.debug:
                                        print 'retrieving EEG', looper3, looper2, looper1
                                    timeframe_channel_dset = f[path]         
                                    eeg=timeframe_channel_dset.value

                                    ###Get maps to sort by
                                    input_model_other = 'model_maps'
                                    input_model_vps = 'microstate'
                                    #sortby_type = sort_type
                                    labelbymaps, eeg_re=get_sortby(sortby_type, looper2, looper1, inputhdf_sortby, eeg_info_mot_obj.VP, eeg_info_mot_obj.Obercond, eeg_info_mot_obj.VP_Group, eeg_info_mot_obj.All_VPs, eeg, eeg_info_mot_obj.chlist, input_model_other, input_model_vps)
                                    mstate_parameters_and_saving(labelbymaps, eeg_re, confobj, looper2, looper1, mstate_parameters_all, looper2_group, groups, looper3)
                                   
                                except:
                                    print 'not found', looper3, looper2, looper1, eeg_input

        if spss_sheets == True:
            create_parameter_spss_sheets(confobj, eeg_info_mot_obj, inputfolder, mstate_parameters_all, groups = True)
    else:
        ###Create dictionary to save VP wise output
        if sort_type=='microstate_parameters_for_each_vp_for_each_run':
            mstate_parameters_all=dict.fromkeys(eeg_info_mot_obj.VP)
            for key in mstate_parameters_all.keys():
                mstate_parameters_all[key]=dict.fromkeys(eeg_info_mot_obj.Cond)
        with closing( h5py.File(sorted_maps_hdf) ) as g:
            with closing( h5py.File(mstate_parameters_output_hdf, 'w') ) as h:
               for looper2 in eeg_info_mot_obj.VP:
                   looper2_group = h.create_group( looper2 )
                   for looper1 in eeg_info_mot_obj.Cond:
                        with closing( h5py.File(eeg_hdf) ) as f:
                            try:
                                timeframe_channel_dset = f['/%s/%s/%s' % (looper2, looper1, eeg_input)]
                            except:
                                print 'not found', looper2, looper1, eeg_input
                                continue
                    
                            path = '/%s/%s/%s' % (looper2, looper1, eeg_input)

                            if confobj.debug:
                                print 'retrieving EEG', looper2, looper1
                            timeframe_channel_dset = f[path]         
                            eeg=timeframe_channel_dset.value

                            ###Get maps to sort by
                            input_model_other = 'model_maps'
                            input_model_vps = 'microstate'
                            #sortby_type = sort_type

                            labelbymaps, eeg_re=get_sortby(sortby_type, looper2, looper1, inputhdf_sortby, eeg_info_mot_obj.VP, eeg_info_mot_obj.Obercond, eeg_info_mot_obj.VP_Group, eeg_info_mot_obj.All_VPs, eeg, eeg_info_mot_obj.chlist, input_model_other, input_model_vps)

                            mstate_parameters_and_saving(labelbymaps, eeg_re, confobj, looper2, looper1, mstate_parameters_all, looper2_group)

        if spss_sheets == True:
            create_parameter_spss_sheets(confobj, eeg_info_mot_obj, inputfolder, mstate_parameters_all)





















def create_parameter_spss_sheets(confobj, eeg_info_mot_obj, inputfolder, mstate_parameters_all, groups = False):                   
    ##################################################################################################################################################################
    ##################################################################################################################################################################
    #########################################################         Prepare SPSS Sheets     ########################################################################
    ##################################################################################################################################################################
    ##################################################################################################################################################################
    import csv
    three__measures = ['occ', 'dur', 'cov']
    spss_parameters_csv = op.join( inputfolder, 'spss_parameters.csv')

    if groups:
        vp_group = []
        #get max epoch nr across all vps and all conds
        max_len = 0
        for outter_outti in mstate_parameters_all.keys():
            vp_group.append(outter_outti)
            for outti in mstate_parameters_all[outter_outti].keys():
                if mstate_parameters_all[outter_outti][outti] == None:
                    continue
                else:
                    for inni in mstate_parameters_all[outter_outti][outti].keys():
                        for meas in three__measures:
                            if mstate_parameters_all[outter_outti][outti][inni] == None:
                                continue
                            else: 
                                if len(mstate_parameters_all[outter_outti][outti][inni][meas])> max_len:
                                    max_len=len(mstate_parameters_all[outter_outti][outti][inni][meas])



        #########
        ###File with all epochs for each measure
        #########

        ###Seperate file for each measure
        for meas in three__measures:
            header = []
            header.append('VP; Group;')
            for inni in eeg_info_mot_obj.Cond:
                for epochnr in range(max_len):
                    for mapnr in range(confobj.original_nr_of_maps):
                        header.append('{0}_{1}_ep{2}_map{3};'.format(inni, meas, epochnr, mapnr))


            spss_parameters_csv = op.join( inputfolder, "{0}_.csv" . format(meas))
            with open(spss_parameters_csv, 'w') as spss_parameter_file:
                spss_parameter_file.writelines(header)

                for outti in eeg_info_mot_obj.VP:
                    spss_parameter_file.write('\n')
                    spss_parameter_file.write('{0};{1}'.format(outti[2:], vp_group) )
                    for groupi in vp_group:
                        try:
                            mstate_parameters_all[groupi][outti]
                        except:
                            continue
                        for inni in eeg_info_mot_obj.Cond:                      
                            for epochnr in range(max_len):
                                try:
                                    mstate_parameters_all[outti][inni][meas][epochnr]
                                    for mapnr in range(confobj.original_nr_of_maps):
                                        spss_parameter_file.write(';')
                                        spss_parameter_file.write("{0:.2f}".format(mstate_parameters_all[groupi][outti][inni][meas][epochnr][mapnr]))
                                except:
                                    for mapnr in range(confobj.original_nr_of_maps):
                                        spss_parameter_file.write(';')
                                        spss_parameter_file.write("{0}".format(999))
 

      
        #########
        ###File with mean across epochs for each measure
        #########

        parameters_mean=dict.fromkeys(three__measures)
        for meas in three__measures:
            parameters_mean[meas]=dict.fromkeys(VP)
            for outti in eeg_info_mot_obj.VP:
                parameters_mean[meas][outti]=dict.fromkeys(Cond)
                for inni in eeg_info_mot_obj.Cond:
                    parameters_mean[meas][outti][inni]=dict.fromkeys(range(confobj.original_nr_of_maps))



        for meas in three__measures:           
            for outti in eeg_info_mot_obj.VP:
                for groupi in vp_group:
                    try:
                        mstate_parameters_all[groupi][outti]
                    except:
                        continue

                    for inni in eeg_info_mot_obj.Cond:
                        for mapnr in range(confobj.original_nr_of_maps):   
                            list_to_avg = [] 
                            try:                         
                                for epochnr in range(len(mstate_parameters_all[groupi][outti][inni][meas])):
                                    list_to_avg.append(mstate_parameters_all[groupi][outti][inni][meas][epochnr][mapnr])
                                parameters_mean[meas][outti][inni][mapnr]=np.mean(list_to_avg)         
                
                            except:                      
                                if confobj.debug:
                                    print meas, outti, inni, mapnr
                                #print mstate_parameters_all[outti][inni][meas][epochnr][mapnr]
                                continue


        ###Seperate file for each measure
        for meas in three__measures:
            header = []
            header.append('VP;')
            for inni in eeg_info_mot_obj.Cond:
                for mapnr in range(confobj.original_nr_of_maps):
                    header.append('{0}_{1}_map{2};'.format(inni, meas, mapnr))


            spss_parameters_csv = op.join( inputfolder, "{0}_means.csv" . format(meas))
            with open(spss_parameters_csv, 'w') as spss_parameter_file:
                spss_parameter_file.writelines(header)
                for outti in eeg_info_mot_obj.VP:
                    spss_parameter_file.write('\n')
                    spss_parameter_file.write('{0}'.format(outti[2:]) )
                    for inni in eeg_info_mot_obj.Cond:
                        for mapnr in range(confobj.original_nr_of_maps):
                            try:
                                spss_parameter_file.write(';')
                                spss_parameter_file.write("{0:.2f}".format(parameters_mean[meas][outti][inni][mapnr]))
                            except:
                                spss_parameter_file.write("{0}".format(999))

    else:
        #get max epoch nr across all vps and all conds
        max_len = 0
        for outti in eeg_info_mot_obj.VP:
            for inni in eeg_info_mot_obj.Cond:
                for meas in three__measures:
                    try:
                        if len(mstate_parameters_all[outti][inni][meas])> max_len:
                            max_len=len(mstate_parameters_all[outti][inni][meas])
                    except:
                        continue


        #########
        ###File with all epochs for each measure
        #########

        ###Seperate file for each measure
        for meas in three__measures:
            header = []
            header.append('VP;')
            for inni in eeg_info_mot_obj.Cond:
                for epochnr in range(max_len):
                    for mapnr in range(confobj.original_nr_of_maps):
                        header.append('{0}_{1}_ep{2}_map{3};'.format(inni, meas, epochnr, mapnr))


            spss_parameters_csv = op.join( inputfolder, "{0}_.csv" . format(meas))
            with open(spss_parameters_csv, 'w') as spss_parameter_file:
                spss_parameter_file.writelines(header)
                for outti in eeg_info_mot_obj.VP:
                    spss_parameter_file.write('\n')
                    spss_parameter_file.write('{0}'.format(outti[2:]) )
                    for inni in eeg_info_mot_obj.Cond:       
                        try:
                            for epochnr in range(max_len):
                                try:
                                    mstate_parameters_all[outti][inni][meas][epochnr]
                                    for mapnr in range(confobj.original_nr_of_maps):
                                        spss_parameter_file.write(';')
                                        spss_parameter_file.write("{0:.2f}".format(mstate_parameters_all[outti][inni][meas][epochnr][mapnr]))
                                except:
                                    for mapnr in range(confobj.original_nr_of_maps):
                                        spss_parameter_file.write(';')
                                        spss_parameter_file.write("{0}".format(999))
                        except:
                            pass
 

      
        #########
        ###File with mean across epochs for each measure
        #########

        parameters_mean=dict.fromkeys(three__measures)
        for meas in three__measures:
            parameters_mean[meas]=dict.fromkeys(VP)
            for outti in eeg_info_mot_obj.VP:
                parameters_mean[meas][outti]=dict.fromkeys(Cond)
                for inni in eeg_info_mot_obj.Cond:
                    parameters_mean[meas][outti][inni]=dict.fromkeys(range(confobj.original_nr_of_maps))



        for meas in three__measures:           
            for outti in eeg_info_mot_obj.VP:
                for inni in eeg_info_mot_obj.Cond:
                    for mapnr in range(confobj.original_nr_of_maps):   
                        list_to_avg = [] 
                        try:                         
                            for epochnr in range(len(mstate_parameters_all[outti][inni][meas])):
                                list_to_avg.append(mstate_parameters_all[outti][inni][meas][epochnr][mapnr])
                            parameters_mean[meas][outti][inni][mapnr]=np.mean(list_to_avg)         
                
                        except:                      
                            if confobj.debug:
                                print meas, outti, inni, mapnr
                            #print mstate_parameters_all[outti][inni][meas][epochnr][mapnr]
                            continue


        ###Seperate file for each measure
        for meas in three__measures:
            header = []
            header.append('VP;')
            for inni in eeg_info_mot_obj.Cond:
                for mapnr in range(confobj.original_nr_of_maps):
                    header.append('{0}_{1}_map{2};'.format(inni, meas, mapnr))


            spss_parameters_csv = op.join( inputfolder, "{0}_means.csv" . format(meas))
            with open(spss_parameters_csv, 'w') as spss_parameter_file:
                spss_parameter_file.writelines(header)
                for outti in eeg_info_mot_obj.VP:
                    spss_parameter_file.write('\n')
                    spss_parameter_file.write('{0}'.format(outti[2:]) )
                    for inni in eeg_info_mot_obj.Cond:
                        for mapnr in range(confobj.original_nr_of_maps):
                            try:
                                spss_parameter_file.write(';')
                                spss_parameter_file.write("{0:.2f}".format(parameters_mean[meas][outti][inni][mapnr]))
                            except:
                                spss_parameter_file.write("{0}".format(999))