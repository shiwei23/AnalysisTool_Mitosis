



# shared variables for measurement outputs


import os, path
import yaml
 
# Function to get and keep track of some shared parameters that typically specific to data/exp
def get_analysis_specific_parameters(data_saved_folder, 
                                     data_folder_num = 'MULTI',
                                     num_of_ch=2,
                                     ch_dict = {"STIM1":2,"Actin":1},
                                     ch_to_be_analyzed = [1,2],
                                     dna_ch_index= 1,
                                     myo_ch_index=1,
                                     er_ch_index=2,
                                     cell_ch_index = 1,
                                     analyze_er = True,
                                     _actual_myo_ch = 'NUP',
                                     overwrite_parameters=False):



    param_fname = os.path.join(data_saved_folder, 'analysis_param_info.yaml')

    if overwrite_parameters or not os.path.exists(param_fname):
        print('Use and save/overwrite input parameters.')
        data_folder_num = data_folder_num
        num_of_ch=num_of_ch
        ch_dict =ch_dict     
        ch_to_be_analyzed=ch_to_be_analyzed
        dna_ch_index= dna_ch_index
        myo_ch_index=myo_ch_index
        er_ch_index=er_ch_index
        cell_ch_index=cell_ch_index
        analyze_er=analyze_er
        _actual_myo_ch=_actual_myo_ch

        param_dict = {'data_folder_num':data_folder_num, 'num_of_ch':num_of_ch,'ch_dict':ch_dict,'ch_to_be_analyzed':ch_to_be_analyzed,
        'dna_ch_index':dna_ch_index, 'myo_ch_index':myo_ch_index, 'er_ch_index':er_ch_index, 'cell_ch_index':cell_ch_index,
        'analyze_er':analyze_er,'_actual_myo_ch':_actual_myo_ch}
        
    
        with open(param_fname, 'w') as outfile:
            yaml.dump(param_dict, outfile, default_flow_style=False)

    else:
        print('Load parameters from the existing saved file.')
        print('Mannualy delete or change overwrite settings if want to make changes.')
        with open(param_fname, 'r') as file: 
            param_dict = yaml.safe_load(file)
        
        data_folder_num = param_dict['data_folder_num']
        num_of_ch=param_dict['num_of_ch']
        ch_dict =param_dict['ch_dict']        
        ch_to_be_analyzed=param_dict['ch_to_be_analyzed']
        dna_ch_index= param_dict['dna_ch_index']
        myo_ch_index=param_dict['myo_ch_index']
        er_ch_index=param_dict['er_ch_index']
        cell_ch_index=param_dict['cell_ch_index']
        analyze_er=param_dict['analyze_er']
        _actual_myo_ch=param_dict['_actual_myo_ch']

    return  (data_folder_num,
        num_of_ch,
        ch_dict,
        ch_to_be_analyzed,
        dna_ch_index,
        myo_ch_index,
        er_ch_index,
        cell_ch_index,
        analyze_er,
        _actual_myo_ch)



# Function to get and keep track of some shared parameters that mostly don't change
def get_analysis_shared_parameters (overwrite_measurement =True,
                                    timelapsed = True,
                                    _analyze_laggards = False,
                                    _separate_inner_core_segment_measurement =True,
                                    analyze_cell =True,
                                    exclude_transmitted=False,
                                    inverse_dna=False,
                                    ) :
    print('Get shared parameters.')                                    
    # overwrite measurement or save as new measurement ONCE after the analysis
    overwrite_measurement =overwrite_measurement

    # FOVs are different timepoints of the same cell of different FOVs (True or False)
    timelapsed = timelapsed

    # Analyze laggards or not (True or False)
    _analyze_laggards =_analyze_laggards

    # Save measurements for inner core segments separately or not (True or False)
    _separate_inner_core_segment_measurement = _separate_inner_core_segment_measurement

    # Analyze proteins for the whole cell or not (True or False)
    analyze_cell = analyze_cell
    if analyze_cell == True:
        num_of_segments_for_cell = 3
        save_global_measurement = True
    else:
        num_of_segments_for_cell = 3
        save_global_measurement = True

    ###########################################################################################
    # Image channel parameters: 
    # If manuanlly exclude Transmitted; 
    # enable this only when the Transmitted has different timpoints compared to other channels;
    # when Transmitted has same timpoints, disable this as the "ch_to_be_analyzed" below will exclude it;
    exclude_transmitted = exclude_transmitted

    # use negative dna exclusion or positive dna signal for dna segmentation 
    inverse_dna = inverse_dna

    ##################################################################
    # Other analysis related parameters:
    # extra erosion for inverse chromosome segmentation
    extra_erosion_factor = 5
    # number of inner core bins to use
    num_of_segments = 5
    # dilation of ER to find adjacent nup signal
    _er_nup_overlap_factor =25
    # percentile of area is the cell mask for the cropped image
    cytoplasm_ratio = 60

    return (overwrite_measurement,
     timelapsed, _analyze_laggards,
     _separate_inner_core_segment_measurement,
      analyze_cell, num_of_segments_for_cell,
      save_global_measurement, exclude_transmitted,
       inverse_dna, extra_erosion_factor, 
       _er_nup_overlap_factor,
       num_of_segments, cytoplasm_ratio)
    

# Function to get and keep track of some er (1st) channel shared parameters
def get_er_ch_parameters (objective_type='100x', non_ridge_segmentation = False, enriched_membrane='ER'): 

    print(f'Get ER (membrane) channel-related parameters for {objective_type} objectives.')         

    if objective_type == '60x':
        #non_ridge_segmentation = True
        non_ridge_segmentation = non_ridge_segmentation
        if non_ridge_segmentation:
            enriched_membrane = enriched_membrane
        else:
            enriched_membrane = 'ER'
        # main chromosome size filter
        large_object_size = 700
        # sigma used for detecting ER rigde
        ridge_sigma = 1  # for 60x
        # membrane fragment size th for ER
        sm_fragment_size = 50 # for 60x
        # The band thickness for the inner core mask
        er_core_thickness = 6   
    
    elif objective_type == '100x':
        non_ridge_segmentation = non_ridge_segmentation
        if non_ridge_segmentation:
            enriched_membrane = 'NE'
        else:
            enriched_membrane = 'ER'
        # main chromosome size filter
        large_object_size = 3000
        # sigma used for detecting ER rigde
        ridge_sigma = 2
        # membrane fragment size th for ER
        sm_fragment_size = 75
        # The band thickness for the inner core mask
        er_core_thickness = 10

    return non_ridge_segmentation, enriched_membrane, large_object_size, ridge_sigma, sm_fragment_size, er_core_thickness



# Function to get and keep track of some myosin (2nd) channel shared parameters
def get_myo_ch_parameters (_actual_myo_ch, er_core_thickness=10, ridge_sigma=2, sm_fragment_size=75): 


    print(f'Get the second channel-related parameters for {_actual_myo_ch} signal.')  

    #==============================================================#
    # below uses ridge segmentation
    if _actual_myo_ch == 'ACTIN':
        analyze_myo =  True
        myo_df_prefix = 'actin' 
        #myo_df_prefix = 'tublin'  # only modify it temporarily for sirTubulin analysis using actin module
        # skip std for Actin
        std_ratio =None
        nup_rim_adjust_size = 0
        nup_core_thickness = 10
        DNA_exclusion_for_NUP = False
        actin_ridge_sigma = 1
        actin_fragment_size = 100  # for 100x
    
    if _actual_myo_ch == 'ER':
        analyze_myo =  True
        myo_df_prefix = 'emerin'  
        # skip std for Actin
        std_ratio =None
        nup_rim_adjust_size = 0
        nup_core_thickness = er_core_thickness
        DNA_exclusion_for_NUP = False    
        actin_ridge_sigma = ridge_sigma
        actin_fragment_size = sm_fragment_size
    #==========================================================#   
    # below uses std ratio segmentation
    if _actual_myo_ch == 'NE': 
        analyze_myo =  True
        myo_df_prefix = 'emerin'
        # Default setting (no need to change)
        nup_rim_adjust_size = 2
        nup_core_thickness = 15
        # std_ratio to adjust NUP foci finding
        std_ratio =0.5
        # do not use DNA exclusion for NUP foci finding
        DNA_exclusion_for_NUP = False   
        # not used below 
        actin_ridge_sigma = 1
        actin_fragment_size = 100  # for 100x
        
    if _actual_myo_ch == 'NUP': 
        analyze_myo =  True
        myo_df_prefix = 'nup'
        # Default setting (no need to change)
        nup_rim_adjust_size = 4
        nup_core_thickness = 15
        # std_ratio to adjust NUP foci finding
        std_ratio =1.5
        # do not use DNA exclusion for NUP foci finding
        DNA_exclusion_for_NUP = False
        # not used below
        actin_ridge_sigma = 1
        actin_fragment_size = 100  # for 100x
    
    if _actual_myo_ch == 'MYO':
        analyze_myo =  True
        myo_df_prefix = 'myo'
        # Default setting (no need to change)
        nup_rim_adjust_size = 0
        nup_core_thickness = 30
        # std_ratio to adjust NUP foci finding
        std_ratio =3
        # do not use DNA exclusion for NUP foci finding
        DNA_exclusion_for_NUP = True
        # not used below
        actin_ridge_sigma = 1
        actin_fragment_size = 100  # for 100x
    #==========================================================#    
    if _actual_myo_ch == 'NONE':
        analyze_myo =  False
        myo_df_prefix = 'none'
        std_ratio =3
        DNA_exclusion_for_NUP = True
        nup_core_thickness =30
        nup_rim_adjust_size=0
        # not used below
        actin_ridge_sigma = 1
        actin_fragment_size = 100  # for 100x

    return analyze_myo, myo_df_prefix, std_ratio, DNA_exclusion_for_NUP, nup_core_thickness, nup_rim_adjust_size,actin_ridge_sigma,actin_fragment_size



# Function to get and keep track of the measurement ouputs
def get_mitotic_analysis_measurement_columns (myo_df_prefix): 
    # Define measurment columns names for the experiment below
    # If need to modify these columns, also need to change within jupyter pipelinefor each measurement


    print(f'Get the dataframe measurement columns.')  

    # columns for inner core related measurement
    chromosome_measurement_columns =['cell_id','timepoint',
                    'object_type', 'object_id','object_z_slice',
                         'er_recruit_normalized_area', 
                    f'{myo_df_prefix}_recruit_normalized_area',
                    'er_recruit_area', f'{myo_df_prefix}_recruit_area', 'er_dna_rim',
                        f'{myo_df_prefix}_dna_rim',
                        'inner_core_segment_number', 'chromosome_segreg_dist']

    # columns for cell/endplasmic related measurement
    cell_measurement_columns =['data_path','cell_id','timepoint',
                            'object_type', 'object_z_slice',                          
                           # dna-associated measurements
                           f'dna_endo_{myo_df_prefix}_area','dna_endo_er_area',
                            f'dna_endo_{myo_df_prefix}_normalized_area',                       
                            # nondna-ER-associated measurements
                            f'nondna_endo_{myo_df_prefix}_area', 'nondna_endo_er_area', 
                            f'nondna_endo_{myo_df_prefix}_normalized_area',                        
                             # periphery nondna ER-associated measurements  
                            f'periphery_nondna_endo_{myo_df_prefix}', 
                             'periphery_nondna_er_area', 
                              f'periphery_endo_{myo_df_prefix}_normalized_area',                      
                             # midzone nondna ER-associated measurements  
                              f'midzone_nondna_endo_{myo_df_prefix}',
                            'midzone_nondna_er_area', 
                             f'midzone_endo_{myo_df_prefix}_normalized_area',                       
                              # midzone er recruitment
                               'midzone_area',       
                                'midzone_nondna_er_normalized_area', 
                               # midzone actin recruitment
                               f'midzone_nondna_{myo_df_prefix}_area',
                               f'midzone_nondna_{myo_df_prefix}_normalized_area', 
                              ]
    
    
    # columns for actin or other similar factors related measurement
    cell_actin_measurement_columns =['data_path','cell_id','timepoint',
                            'object_type', 'object_z_slice',   
                                 # normalized actin measrements
                                f'cortical_{myo_df_prefix}_normalized_area',
                                 f'endo_{myo_df_prefix}_normalized_area',
                                f'midzone_{myo_df_prefix}_normalized_area',
                                 # raw and other cortical and midzone measrements
                                 'midzone_er_normalized_area',
                                f'cortical_{myo_df_prefix}_area', f'endo_{myo_df_prefix}_area',
                                f'midzone_{myo_df_prefix}_area', 'midzone_er_area',
                                'cortical_rim_area','endo_area','midzone_area']

    # columns for core and noncore related measurement
    core_noncore_measurement_columns =['data_path','cell_id','timepoint',
                            'object_type', 'object_id','object_z_slice',   
                                   # er channel measurements
                                   'er_recruit_normalized_area', 
                                   'er_recruit_area',
                                  'er_dna_rim',
                                   # myo/nup channel measurements
                                    f'{myo_df_prefix}_recruit_normalized_area',
                                   f'{myo_df_prefix}_recruit_area',
                                   f'{myo_df_prefix}_dna_rim',
                                   # core or noncore info
                                   'core_noncore_type']

    return chromosome_measurement_columns, cell_measurement_columns, cell_actin_measurement_columns, core_noncore_measurement_columns




def get_AL_params (cell_type = 'HeLa'): 

    if cell_type == 'HeLa':
        NE_AL_params = {'NE':{'er_core_thickness':10,'er_rim_adjust':0,'er_erode_adj':True,
                      'nup_core_thickness':10+15,'nup_rim_adjust':15,'nup_erode_adj':True},
               
                 'AL':{'er_core_thickness':20+10,'er_rim_adjust':10,'er_erode_adj':False,
                      'nup_core_thickness':20+10,'nup_rim_adjust':10,'nup_erode_adj':False},
                
                    'nonDNA':20, 'cortical_thickness': 5, 'er_nup_overlap_factor':5, 'midzone_erode':0} 



    elif cell_type == 'test_0924_v2':
        NE_AL_params = {'NE':{'er_core_thickness':5+5,'er_rim_adjust':5,'er_erode_adj':True,
                      'nup_core_thickness':5+15,'nup_rim_adjust':15,'nup_erode_adj':True},
               
                 'AL':{'er_core_thickness':20+5,'er_rim_adjust':5,'er_erode_adj':False,
                      'nup_core_thickness':20+5,'nup_rim_adjust':5,'nup_erode_adj':False},
                
                    'nonDNA':20, 'cortical_thickness': 10, 'er_nup_overlap_factor':5, 'midzone_erode':20} 



    elif cell_type == 'test':
        NE_AL_params = {'NE':{'er_core_thickness':10+0,'er_rim_adjust':0,'er_erode_adj':True,
                      'nup_core_thickness':5+15,'nup_rim_adjust':15,'nup_erode_adj':True},
               
                 'AL':{'er_core_thickness':20+10,'er_rim_adjust':10,'er_erode_adj':False,
                      'nup_core_thickness':20+5,'nup_rim_adjust':5,'nup_erode_adj':False},
                
                    'nonDNA':20, 'cortical_thickness': 20, 'er_nup_overlap_factor':10} 

    elif cell_type == 'test_0924':
        NE_AL_params = {'NE':{'er_core_thickness':10,'er_rim_adjust':0,'er_erode_adj':True,
                      'nup_core_thickness':10+15,'nup_rim_adjust':15,'nup_erode_adj':True},
               
                 'AL':{'er_core_thickness':30+10,'er_rim_adjust':10,'er_erode_adj':False,
                      'nup_core_thickness':30+10,'nup_rim_adjust':10,'nup_erode_adj':False},
                
                    'nonDNA':30, 'cortical_thickness': 10, 'er_nup_overlap_factor':5, 'midzone_erode':25} 



    return NE_AL_params