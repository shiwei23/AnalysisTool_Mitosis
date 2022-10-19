# -*- coding: utf-8 -*-
"""
Created on Sat Jul 10 15:58:59 2021

@author: Shiwei Liu @ Harvard University

scripts for compiling measurement dataframes derived from running the mitosis_analysis.py; also making figures
"""


import pandas as pd
import numpy as np
import os, sys, glob, pathlib
import matplotlib.pyplot as plt
import math
import yaml, tqdm
import pickle




# function to round down time interval 
def round_decimals_down(number:float, decimals:int=2):
    """
    Returns a value rounded down to a specific number of decimal places.
    """
    if not isinstance(decimals, int):
        raise TypeError("decimal places must be an integer")
    elif decimals < 0:
        raise ValueError("decimal places has to be 0 or more")
    elif decimals == 0:
        return math.floor(number)

    factor = 10 ** decimals
    return math.floor(number * factor) / factor


# function to determine the time interval to where closest to the options [30,60,90,120, etc.,]
def determine_time_interval(time_step, cand_time_intervals=None):
    """
    Returns a value rounded down to a specific number of decimal places.
    """
    # candidate time intervals (min)
    if cand_time_intervals is None:
        cand_time_intervals = [round_decimals_down(_i *0.5,1) for _i in range (1,10) ]

    diff_time_list = np.array([abs(time_step-_t) for _t in cand_time_intervals])
    sel_time_interval = cand_time_intervals[np.argmin(diff_time_list)]

    return sel_time_interval



######################################################################

# function to process each experiment condition with the input as specified
def normalize_time_for_mitotic_timelapsed_measurement (cell_folder_list, 
                                                       measurement_column, 
                                                       reference_type = 'integrated', 
                                                       measurement_type = 'integrated',
                                                       normalize_column='er_recruit_normalized_area',
                                                       ref_min =0.1,
                                                       ref_max =0.35,
                                                       bin_of_interest = 1,
                                                       force_time_step = 0,   # if not zero, use this as time interval (min)
                                                       timepoints_to_cover = 10,
                                                       report_time =True):
    # define measurement excel types 
    chromo_object_measurement_types =['integrated','bins','core_noncore','whole_chromosome']
    cell_object_measurement_types =['global_dna_excluded','global']
    
    # report original time_interval for sanity check
    _ori_time_interval =[]

    # initiate measurement dataframe
    _measurement_of_interest_all_chr = pd.DataFrame()
    # loop through folders
    
    # handeling for the weird error during the batch function
    key_error=0
    for _exp_index, _exp_fd in enumerate(cell_folder_list[:]):
        # get file for reference and measurement; currently use the integrated measurement for reference usage
        _reference_file = [_file for _file in glob.glob(_exp_fd) if reference_type in _file]
        _measurement_file = [_file for _file in glob.glob(_exp_fd) if measurement_type in _file]

        # locate chr object by unique object id        
        for _measure, _ref in zip(_measurement_file, _reference_file):
            _ref_df  = pd.read_excel(_ref, index_col =None)
            _measure_df  = pd.read_excel(_measure, index_col =None)   
            _chr_1_ref_df = _ref_df[_ref_df['object_id']==1]
            _chr_2_ref_df = _ref_df[_ref_df['object_id']==2]
            
            # report original time_interval for sanity check
            try:
                _ori_time_interval.append(_ref_df['time_interval'].mean())
            except KeyError:   
                key_error+=1 
            
            # for integrated or bin measurement or core-noncore measurement, etc., each object is a mitotic chromosome mass
            if measurement_type in chromo_object_measurement_types:
                try:
                    _chr_1_df = _measure_df[_measure_df['object_id']==1]
                    _chr_2_df = _measure_df[_measure_df['object_id']==2]
                except KeyError:
                    key_error+=1
                if measurement_type == 'bins': # get bin of interest only
                    _chr_1_df = _chr_1_df[_chr_1_df['inner_core_segment_number']==bin_of_interest]
                    _chr_2_df = _chr_2_df[_chr_2_df['inner_core_segment_number']==bin_of_interest]
                    
                elif measurement_type == 'core_noncore': # get bin of interest only
                    _chr_1_df = _chr_1_df[_chr_1_df['core_noncore_type']==bin_of_interest]
                    _chr_2_df = _chr_2_df[_chr_2_df['core_noncore_type']==bin_of_interest]
                
                # append measurement by time accordingly    
                for _chr_ref_df, _chr_measure_df  in zip([_chr_1_ref_df, _chr_2_ref_df],[_chr_1_df,_chr_2_df]):
                    # for both chr object, use measurement dict to record time and value
                    _measurement_of_interest = {}
                    _normalized_time = 0
                    # for each valid time
                    for _row_ref_df, _row_measure_df in zip(_chr_ref_df.iloc(), _chr_measure_df.iloc()):
                        # start count only when ER/NE has assembled (>0.1 and <=0.35)
                        if (_normalized_time <1 and 
                            _row_ref_df[normalize_column] >= ref_min and
                            _row_ref_df[normalize_column] <=ref_max):
                            _normalized_time+=1
                            _measurement_of_interest[_normalized_time]= _row_measure_df[measurement_column]
                        # add subsequent timepoints
                        elif _normalized_time >=1 and _normalized_time<timepoints_to_cover:
                            if force_time_step>0:
                                _normalized_time+=force_time_step
                            else:
                                _time_step = _row_ref_df['time_interval']/60
                                _normalized_time+=determine_time_interval(_time_step)               
                            # use _normalized_time as dict_key
                            # only append if the er reruitment is measureable
                            if _row_ref_df[normalize_column] >= ref_min:                        
                                _measurement_of_interest[_normalized_time]= _row_measure_df[measurement_column]
                            # if not append nan 
                            else:
                                _measurement_of_interest[_normalized_time]= np.nan
                                
                    # convert measurement_dict to dataframe (non-qualified cells will also be included)
                    _measurement_of_interest = pd.DataFrame(_measurement_of_interest, 
                                                                index=[0])
                    # append other info and concat dataframe 
                    _measurement_of_interest.insert(0,'object_id',_row_measure_df['object_id'])
                    _measurement_of_interest.insert(0,'cell_id',_row_measure_df['cell_id'])
                    _measurement_of_interest.insert(0,'exp_id',_exp_index+1)
                    _measurement_of_interest.insert(0,'exp_path',_exp_fd)
                    _measurement_of_interest_all_chr= pd.concat([_measurement_of_interest_all_chr,
                                                                     _measurement_of_interest])

          ######################################################            
           # for global measurement where each object is a cell
            elif measurement_type in cell_object_measurement_types or 'global' in measurement_type:
                _cell_df = _measure_df
                # for the cell object use measurement dict to record time and value
                _measurement_of_interest = {}
                _normalized_time = 0
                # for each valid time
                for _row_ref_1_df, _row_ref_2_df, _row_measure_df  in zip(_chr_1_ref_df.iloc(), 
                                                                   _chr_2_ref_df.iloc(),
                                                                   _cell_df.iloc()):
                    # start count only when ER/NE has assembled (>0.1 and <=0.35)
                    if _normalized_time <1:
                        if ((_row_ref_1_df[normalize_column] >= ref_min 
                             and _row_ref_1_df[normalize_column] <=ref_max) 
                            or (_row_ref_2_df[normalize_column] >= ref_min
                               and _row_ref_2_df[normalize_column] <=ref_max)):   
                            
                            _normalized_time+=1
                            _measurement_of_interest[_normalized_time]= _row_measure_df[measurement_column]
                    # add subsequent timepoints        
                    elif _normalized_time >=1 and _normalized_time<timepoints_to_cover:
                        if force_time_step>0:
                            _normalized_time+=force_time_step
                        else:
                            _time_step = _row_ref_df['time_interval']/60
                            _normalized_time+=determine_time_interval(_time_step)   
                            
                        # only append if the er reruitment is measureable
                        if (_row_ref_1_df[normalize_column] >= ref_min or
                            _row_ref_2_df[normalize_column] >= ref_min):  
                            _measurement_of_interest[_normalized_time]= _row_measure_df[measurement_column]
                        # if not append nan 
                        else:
                            _measurement_of_interest[_normalized_time]=np.nan
                            
                # convert measurement_dict to dataframe (non-qualified cells will also be included)
                _measurement_of_interest = pd.DataFrame(_measurement_of_interest, 
                                                                index=[0])
                # append other info and concat dataframe 
                _measurement_of_interest.insert(0,'object_id',1) # only 1 object so no object id column from the source
                _measurement_of_interest.insert(0,'cell_id',_row_measure_df['cell_id'])
                _measurement_of_interest.insert(0,'exp_id',_exp_index+1)
                _measurement_of_interest.insert(0,'exp_path',_exp_fd)
                _measurement_of_interest_all_chr= pd.concat([_measurement_of_interest_all_chr,
                                                                     _measurement_of_interest])
    
    # report original time_interval for sanity check
    if report_time:
        print (f'-- The original time interval for each replicate is {_ori_time_interval} in (s).')
    if force_time_step>0:
        print (f'-- Use {force_time_step} (min) for all measurements.')
    
    # reset and drop index column
    _measurement_of_interest_all_chr = _measurement_of_interest_all_chr.reset_index()
    _measurement_of_interest_all_chr = _measurement_of_interest_all_chr.drop(columns=['index'])
    if key_error>0:
        print (f'-- Warning: There are {key_error} KeyErrors for the dataframe reported, n\
               check if results were indeed processed correctly.')
            
    return _measurement_of_interest_all_chr





####################################################
 # Sigmoid function to curve-fit NE assembly of individual cell
def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return (y)


def exponential_plateau (x, YM, Y0, k):
    y = YM - (YM-Y0)*np.exp (-k*x)
    return (y)


def linear_fit (x,a,b):
    y= a*x + b
    return (y)


def curve_fit_exp_plateau_mitotic_cell (xdata,ydata,p0,_time_step,x0_adjust=1,
                                         exp_ratio_th=0.05, r2score_1_th=0.7):
    from scipy.optimize import curve_fit
    from sklearn.metrics import r2_score

    fig, ax = plt.subplots(figsize=(3,2), dpi=300)
    ax = plt.scatter(xdata,ydata)

    # Try linear first
    if len(xdata)>1:
        p = np.polyfit(xdata,ydata,1)
        r2score_1 = r2_score(ydata, np.polyval(p, xdata))
        ax = plt.plot(xdata,np.polyval(p, xdata), color='orange')

        # keep linear results as default
        if p[0]<0:
            x_mid = np.nan
            x_end = xdata[0]
            fit_type = 'saturated'
            delta_y = np.nan
        # use linear
        else:
            # for correction of the initial membranes adjacent to the noncore
            if xdata[-1]/2+ x0_adjust<=xdata[-1]: 
                x_mid = determine_time_interval(xdata[-1]/2+ x0_adjust, xdata)
            else:
                x_mid = determine_time_interval(xdata[-1]/2, xdata)
            x_end = xdata[-1]
            fit_type = 'linear'
            delta_y = np.polyval(p, x_end)-np.polyval(p, x_mid)
    
    else:
        r2score_1 = 0
        x_mid = np.nan
        x_end = xdata[0]
        fit_type = 'saturated'
        delta_y = np.nan


    # Priotitize linear if r2 higher than the desired threshold
    if r2score_1>=r2score_1_th or len(xdata)<2:
        r2score_2=np.nan
        ax = plt.annotate("linear:{:.3f}".format(r2score_1),(xdata[0],max(ydata)),fontsize=5)
        pass
    # Try sigmoid if linear r2 lower than the desired threshold


    else:
        try:
            if len(xdata) >= len(p0):
                popt, pcov = curve_fit(exponential_plateau, xdata, ydata,p0, bounds=([np.min(ydata), 0,0],[np.max(ydata), np.min(ydata),1]))
            else:
                popt, pcov = curve_fit(exponential_plateau, xdata, ydata,p0, method='dogbox',
                bounds=([np.min(ydata),0,0],[np.max(ydata), np.min(ydata),1]))
            YM, Y0, k= popt
            y_fit = [exponential_plateau(_x,YM, Y0, k) for _x in xdata]
            r2score_2 = r2_score(ydata, y_fit)
            
            #"fitted x-mid"
            #x0 = (np.log((YM- ((YM-Y0)/2 + Y0))/(YM-Y0)))/-k 
            # use ymin and ymax in acuqired range
            x0 = (np.log((YM- ((YM-np.min(ydata))/2 + np.min(ydata)))/(YM-Y0)))/-k
            
            
            ax = plt.annotate("linear:{:.3f}".format(r2score_1)+ "\n" + "exp_plateau:{:.3f}".format(r2score_2),
                                     (xdata[0],max(ydata)),fontsize=5)

            xdata_plot = np.arange(xdata[0],xdata[-1]+0.2,0.2)
            ydata_plot = [exponential_plateau(_x,YM, Y0, k) for _x in xdata_plot]
            ax = plt.plot(xdata_plot, ydata_plot, color='red')
            
            # if linear still fits better
            if r2score_2<r2score_1:
                r2score_2=np.nan
                pass
            # if sigmoid performs "better"
            else:
                # if curve goes downward, the plateu is the first, skip the mid timepoint for x,y
                if y_fit[0]>=y_fit[-1]:
                    x_mid = np.nan
                    x_end = xdata[0]
                    fit_type = 'saturated'
                    delta_y = np.nan
                    #r2score = r2score_2

                # if x0 larger than max of the xdata, use linear to estimate the x_end instead
                elif x0>max(xdata):
                    r2score_2=np.nan
                    pass

                # upward sigmoid fit works
                else:
                    #r2score = r2score_2
                    if x0<min(xdata):
                        x_mid = xdata[0]
                    else:
                        # for correction of the initial membranes adjacent to the noncore
                        if x0 + x0_adjust <=max(xdata):
                            x_mid = determine_time_interval(x0 + x0_adjust, xdata) 
                        else:
                            x_mid = determine_time_interval(x0, xdata)         
                    fit_type = 'exp_plateau'
                    # estimate the near-platue with specified scale of change
                    near_ymax=exponential_plateau(10000,YM, Y0, k)
                    half_ymax=exponential_plateau(x0,YM, Y0, k)
                    delta_y = near_ymax-half_ymax
                    # extend the next 100 timepoint to find when the change becomes relatively small--ie close to the near_ymax
                    for _x in [int(x_mid) +  _time_step* i for i in range(1,100)]:
                        test_ymax=exponential_plateau(_x,YM, Y0, k)
                        diff_ratio = (near_ymax-test_ymax)/(near_ymax-half_ymax)
                        if diff_ratio<exp_ratio_th:
                            x_end =determine_time_interval(_x, xdata) 
                            if x_end == x_mid:
                                x_end+=_time_step
                            break
        
            
         # default linear fit if sigmoid fails
        except (RuntimeError,ValueError):
            r2score_2=np.nan
            ax = plt.annotate("linear:{:.3f}".format(r2score_1),(xdata[0],max(ydata)),fontsize=5)
            pass
    
    #print (r2score)
    # close figure in plt
    plt.close(fig)

    return x_mid, x_end, delta_y,fit_type, r2score_1, r2score_2, fig



# 20220505 version to use
def curve_fit_mitotic_cell (xdata,ydata,p0,_time_step,sigmoid_ratio_th=0.1,r2score_1_th=0.9,x0_adjust=1):
    from scipy.optimize import curve_fit
    from sklearn.metrics import r2_score


    fig, ax = plt.subplots(figsize=(3,2), dpi=300)
    ax = plt.scatter(xdata,ydata)


    # Try linear first
    if len(xdata)>1:
        p = np.polyfit(xdata,ydata,1)
        r2score_1 = r2_score(ydata, np.polyval(p, xdata))
        ax = plt.plot(xdata,np.polyval(p, xdata), color='orange')

        # keep linear results as default
        if p[0]<0:
            x_mid = np.nan
            x_end = xdata[0]
            fit_type = 'saturated'
            delta_y = np.nan
        # use linear
        else:
            # for correction of the initial membranes adjacent to the noncore
            if xdata[-1]/2+ x0_adjust<=xdata[-1]: 
                x_mid = determine_time_interval(xdata[-1]/2+ x0_adjust, xdata)
            else:
                x_mid = determine_time_interval(xdata[-1]/2, xdata)
           
            x_end = xdata[-1]
            fit_type = 'linear'
            delta_y = np.polyval(p, x_end)-np.polyval(p, x_mid)
    
    else:
        r2score_1 = 0

    # Priotitize linear if r2 higher than the desired threshold
    if r2score_1>=r2score_1_th:
        r2score_2=np.nan
        ax = plt.annotate("linear:{:.3f}".format(r2score_1),(xdata[0],max(ydata)),fontsize=5)
        pass
    # Try sigmoid if linear r2 lower than the desired threshold


    else:
        try:
            if len(xdata) >= len(p0):
                popt, pcov = curve_fit(sigmoid, xdata, ydata,p0)
            else:
                popt, pcov = curve_fit(sigmoid, xdata, ydata,p0, method='dogbox')
            L ,x0, k, b= popt
            sig = sigmoid(xdata,L ,x0, k, b)
            r2score_2 = r2_score(ydata, sig)

            ax = plt.annotate("linear:{:.3f}".format(r2score_1)+ "\n" + "sigmoid:{:.3f}".format(r2score_2),
                                     (xdata[0],max(ydata)),fontsize=5)

            xdata_plot = np.arange(xdata[0],xdata[-1]+0.2,0.2)
            ax = plt.plot(xdata_plot, sigmoid(xdata_plot,L ,x0, k, b), color='red')
            
            # if linear still fits better
            if r2score_2<r2score_1:
                r2score_2=np.nan
                pass
            # if sigmoid performs "better"
            else:
                # if curve goes downward, the plateu is the first, skip the mid timepoint for x,y
                if sig[0]>=sig[-1]:
                    x_mid = np.nan
                    x_end = xdata[0]
                    fit_type = 'saturated'
                    delta_y = np.nan

                # if x0 larger than max of the xdata, use linear to estimate the x_end instead
                elif x0>max(xdata):
                    r2score_2=np.nan
                    pass

                # upward sigmoid fit works
                else:
                    #r2score = r2score_2
                    if x0<min(xdata):
                        x_mid = xdata[0]
                    else:
                        if x0 + x0_adjust <=max(xdata):
                            x_mid = determine_time_interval(x0 + x0_adjust, xdata) 
                        else:
                            x_mid = determine_time_interval(x0, xdata)         
                    fit_type = 'sigmoid'
                    # estimate the near-platue with specified scale of change
                    near_ymax=sigmoid(10000,L ,x0, k, b)
                    half_ymax=sigmoid(x0,L ,x0, k, b)
                    delta_y = near_ymax-half_ymax
                    # extend the next 100 timepoint to find when the change becomes relatively small--ie close to the near_ymax
                    for _x in [int(x_mid) +  _time_step* i for i in range(1,100)]:
                        test_ymax=sigmoid(_x,L ,x0, k, b)
                        diff_ratio = (near_ymax-test_ymax)/(near_ymax-half_ymax)
                        if diff_ratio<sigmoid_ratio_th:
                            x_end =determine_time_interval(_x, xdata) 
                            if x_end == x_mid:
                                x_end+=_time_step
                            break
        
            
         # default linear fit if sigmoid fails
        except RuntimeError:
            r2score_2=np.nan
            ax = plt.annotate("linear:{:.3f}".format(r2score_1),(xdata[0],max(ydata)),fontsize=5)
            pass
    
    #print (r2score)
    # close figure in plt
    plt.close(fig)

    return x_mid, x_end, delta_y,fit_type, r2score_1,r2score_2, fig


######################################################################

# function to process each experiment condition with the input as specified
def curve_fit_normalize_time_for_mitotic_timelapsed_measurement (cell_folder_list, 
                                                       measurement_column, 
                                                       reference_type = 'integrated', 
                                                       measurement_type = 'integrated',
                                                       normalize_column='er_recruit_normalized_area',
                                                       bin_of_interest = 1,
                                                       sigmoid_ratio_th=0.05, 
                                                       fit_method = 'sigmoid',
                                                       x0_adjust=1,
                                                       save_fit=True,
                                                       save_folder=None,
                                                       exp_type ='Condition_1',
                                                       report_time =True):


    # define measurement excel types 
    chromo_object_measurement_types =['integrated','bins','core_noncore','whole_chromosome', 
                                      'whole_chromosome_NE','whole_chromosome_AL', 'core_noncore_NE', 'core_noncore_AL',
                                      ]
    cell_object_measurement_types =['global_dna_excluded','global','global_AL']
    
    # report original time_interval for sanity check
    _ori_time_interval =[]

    # initiate measurement dataframe
    _measurement_of_interest_all_chr = pd.DataFrame()
    # loop through folders
    
    # handeling for the weird error during the batch function
    key_error=0
    for _exp_index, _exp_fd in enumerate(cell_folder_list[:]):
        # get file for reference and measurement; currently use the integrated measurement for reference usage
        _reference_file = [_file for _file in glob.glob(_exp_fd) if reference_type in _file]
        _measurement_file = [_file for _file in glob.glob(_exp_fd) if measurement_type in _file]

        # locate chr object by unique object id        
        for _measure, _ref in zip(_measurement_file, _reference_file):
            _ref_df  = pd.read_excel(_ref, index_col =None)
            _measure_df  = pd.read_excel(_measure, index_col =None)   
            # only proceed if excel has measurement
            if len(_ref_df) >0 and len(_measure_df)>0:           
                _chr_1_ref_df = _ref_df[_ref_df['object_id']==1]
                _chr_2_ref_df = _ref_df[_ref_df['object_id']==2]    
                _time_step =  (_ref_df.iloc[0,:]['time_interval'])/60
                _time_step = determine_time_interval(_time_step)   
            
                # report original time_interval for sanity check
                try:
                    _ori_time_interval.append(_ref_df['time_interval'].mean())
                except KeyError:   
                    key_error+=1 
            
                # for integrated or bin measurement or core-noncore measurement, etc., each object is a mitotic chromosome mass
                if measurement_type in chromo_object_measurement_types:
                    try:
                        _chr_1_df = _measure_df[_measure_df['object_id']==1]
                        _chr_2_df = _measure_df[_measure_df['object_id']==2]
                    except KeyError:
                        key_error+=1
                    if measurement_type == 'bins': # get bin of interest only
                        _chr_1_df = _chr_1_df[_chr_1_df['inner_core_segment_number']==bin_of_interest]
                        _chr_2_df = _chr_2_df[_chr_2_df['inner_core_segment_number']==bin_of_interest]
                    
                    elif measurement_type in ['core_noncore','core_noncore_NE', 'core_noncore_AL']: # get bin of interest only
                        _chr_1_df = _chr_1_df[_chr_1_df['core_noncore_type']==bin_of_interest]
                        _chr_2_df = _chr_2_df[_chr_2_df['core_noncore_type']==bin_of_interest]
                
                    # curve fit measurement by time accordingly    
                    for _chr_ref_df, _chr_measure_df  in zip([_chr_1_ref_df, _chr_2_ref_df],[_chr_1_df,_chr_2_df]):
                        # get xdata and ydata for curve fitting
                        ydata = _chr_ref_df[normalize_column].to_list()
                        xdata = _chr_ref_df['timepoint'].to_list()
                        xdata = [(int(_x)*_time_step) for _x in xdata]
                        #xdata = [determine_time_interval(_x, xdata) for _x in xdata] # this may give 'rounding' error for timepoints
                        ydata_measure = _chr_measure_df[measurement_column].to_list()

                        if fit_method == 'sigmoid' or fit_method =='linear':
                            # use the data features to make a initial guess
                            p0 = [max(ydata)-min(ydata), np.median(xdata),1,min(ydata)]
                            # get x_mid, x_end, fit_type using the function above
                            if fit_method == 'sigmoid':
                                x_mid, x_end, delta_y,fit_type, r2score_1,r2score_2, fig = curve_fit_mitotic_cell (xdata,ydata,p0,_time_step,
                                                                                                sigmoid_ratio_th=sigmoid_ratio_th,
                                                                                                x0_adjust=x0_adjust,
                                                                                                r2score_1_th=0.9)
                            elif fit_method =='linear':
                                x_mid, x_end, delta_y,fit_type, r2score_1,r2score_2, fig = curve_fit_mitotic_cell (xdata,ydata,p0,_time_step,
                                                                                                sigmoid_ratio_th=sigmoid_ratio_th,
                                                                                                x0_adjust=x0_adjust,
                                                                                                r2score_1_th=0.6)

                            if fit_type=='linear' or 'sigmoid' in fit_type:
                                y_mid = ydata_measure[np.argwhere(np.array(xdata)==x_mid).ravel()[0]]
                                norm_xdata = [(_x-x_end) for _x in xdata]
                                norm_x_mid = x_mid-x_end
                            elif fit_type=='saturated':
                                y_mid = np.nan
                                norm_xdata = xdata
                                norm_x_mid = np.nan
                        elif fit_method == 'exp_plateau':
                            p0 = [max(ydata),min(ydata),1]
                            x_mid, x_end, delta_y,fit_type, r2score_1,r2score_2, fig = curve_fit_exp_plateau_mitotic_cell (xdata,ydata,p0,_time_step,
                            x0_adjust=x0_adjust,exp_ratio_th=sigmoid_ratio_th)

                            if fit_type=='linear' or 'exp_plateau' in fit_type:
                                y_mid = ydata_measure[np.argwhere(np.array(xdata)==x_mid).ravel()[0]]
                                norm_xdata = [(_x-x_end) for _x in xdata]
                                norm_x_mid = x_mid-x_end
                            elif fit_type=='saturated':
                                y_mid = np.nan
                                norm_xdata = xdata
                                norm_x_mid = np.nan

                        # add normalized result accordingly
                        _measurement_of_interest = {}
                        for _x, _y in zip(norm_xdata, ydata_measure):
                            _measurement_of_interest[_x]= _y

                        # convert measurement_dict to dataframe (non-qualified cells will also be included)
                        _measurement_of_interest = pd.DataFrame(_measurement_of_interest, 
                                                                index=[0])
                                         
                        # append other info and concat dataframe 
                        _measurement_of_interest.insert(0,'mid_timepoint_level',y_mid)
                        _measurement_of_interest.insert(0,'mid_timepoint',norm_x_mid)
                        _measurement_of_interest.insert(0,'Ymax-1/2_Ymax',delta_y)
                        _measurement_of_interest.insert(0,'R_squared_sigmoid',r2score_2)
                        _measurement_of_interest.insert(0,'R_squared_linear',r2score_1)
                        _measurement_of_interest.insert(0,'fit_type',fit_type)

                        _cell_id = _chr_ref_df.iloc[0]['cell_id']
                        _object_id = _chr_ref_df.iloc[0]['object_id']

                        _measurement_of_interest.insert(0,'object_id',_object_id)
                        _measurement_of_interest.insert(0,'cell_id', _cell_id)
                        _measurement_of_interest.insert(0,'exp_id',_exp_index+1)
                        _measurement_of_interest.insert(0,'exp_path',_exp_fd)
                        _measurement_of_interest_all_chr= pd.concat([_measurement_of_interest_all_chr,
                                                                     _measurement_of_interest])

                        # save curve fit figure if desired
                        if save_fit and os.path.exists(save_folder):
                            save_subfolder = os.path.join(save_folder, f'Curve_fit_{fit_method}', exp_type+f'_{reference_type}')
                            if not os.path.exists(save_subfolder):
                                os.makedirs(save_subfolder)
                            fig.savefig(os.path.join(save_subfolder,f'exp_{_exp_index+1}_{_cell_id}_chr_{_object_id}.png'),bbox_inches='tight')
                            plt.cla()
                            plt.close(fig)

                # for global measurement where each object is a cell
                elif measurement_type in cell_object_measurement_types or 'global' in measurement_type:
                    _cell_df = _measure_df
                    # use the first chr df to get cell-related info
                    _chr_ref_df =_chr_1_ref_df

                    # indentation purpose
                    _use_single_chr = True
                    if _use_single_chr:
                        # get xdata and ydata for curve fitting
                        ydata = _chr_ref_df[normalize_column].to_list()
                        xdata = _chr_ref_df['timepoint'].to_list()
                        xdata = [(int(_x)*_time_step) for _x in xdata]
                        #xdata = [determine_time_interval(_x, xdata) for _x in xdata] # this may give 'rounding' error for timepoints
                        ydata_measure = _cell_df[measurement_column].to_list()

                        if fit_method == 'sigmoid' or fit_method =='linear':
                            # use the data features to make a initial guess
                            p0 = [max(ydata)-min(ydata), np.median(xdata),1,min(ydata)]
                            # get x_mid, x_end, fit_type using the function above
                            if fit_method == 'sigmoid':
                                x_mid, x_end, delta_y,fit_type, r2score_1,r2score_2, fig = curve_fit_mitotic_cell (xdata,ydata,p0,_time_step,
                                                                                                sigmoid_ratio_th=sigmoid_ratio_th,
                                                                                                x0_adjust=x0_adjust,
                                                                                                r2score_1_th=0.9)
                            elif fit_method =='linear':
                                x_mid, x_end, delta_y,fit_type, r2score_1,r2score_2, fig = curve_fit_mitotic_cell (xdata,ydata,p0,_time_step,
                                                                                                sigmoid_ratio_th=sigmoid_ratio_th,
                                                                                                x0_adjust=x0_adjust,
                                                                                                r2score_1_th=0.6)

                            if fit_type=='linear' or 'sigmoid' in fit_type:
                                y_mid = ydata_measure[np.argwhere(np.array(xdata)==x_mid).ravel()[0]]
                                norm_xdata = [(_x-x_end) for _x in xdata]
                                norm_x_mid = x_mid-x_end
                            elif fit_type=='saturated':
                                y_mid = np.nan
                                norm_xdata = xdata
                                norm_x_mid = np.nan
                        elif fit_method == 'exp_plateau':
                            p0 = [max(ydata),min(ydata),1]
                            x_mid, x_end, delta_y,fit_type, r2score, fig = curve_fit_exp_plateau_mitotic_cell (xdata,ydata,p0,_time_step,
                            x0_adjust=x0_adjust,exp_ratio_th=sigmoid_ratio_th)

                            if fit_type=='linear' or 'exp_plateau' in fit_type:
                                y_mid = ydata_measure[np.argwhere(np.array(xdata)==x_mid).ravel()[0]]
                                norm_xdata = [(_x-x_end) for _x in xdata]
                                norm_x_mid = x_mid-x_end
                            elif fit_type=='saturated':
                                y_mid = np.nan
                                norm_xdata = xdata
                                norm_x_mid = np.nan
                        # add normalized result accordingly
                        _measurement_of_interest = {}
                        for _x, _y in zip(norm_xdata, ydata_measure):
                            _measurement_of_interest[_x]= _y

                        # convert measurement_dict to dataframe (non-qualified cells will also be included)
                        _measurement_of_interest = pd.DataFrame(_measurement_of_interest, 
                                                                index=[0])
                                         
                        # append other info and concat dataframe 
                        _measurement_of_interest.insert(0,'mid_timepoint_level',y_mid)
                        _measurement_of_interest.insert(0,'mid_timepoint',norm_x_mid)
                        _measurement_of_interest.insert(0,'Ymax-1/2_Ymax',delta_y)
                        _measurement_of_interest.insert(0,'R_squared_sigmoid',r2score_2)
                        _measurement_of_interest.insert(0,'R_squared_linear',r2score_1)
                        _measurement_of_interest.insert(0,'fit_type',fit_type)

                        _cell_id = _chr_ref_df.iloc[0]['cell_id']
                        _object_id = _chr_ref_df.iloc[0]['object_id']

                        _measurement_of_interest.insert(0,'object_id',_object_id)
                        _measurement_of_interest.insert(0,'cell_id',_cell_id)
                        _measurement_of_interest.insert(0,'exp_id',_exp_index+1)
                        _measurement_of_interest.insert(0,'exp_path',_exp_fd)
                        _measurement_of_interest_all_chr= pd.concat([_measurement_of_interest_all_chr,
                                                                     _measurement_of_interest])

                        # save curve fit figure if desired
                        if save_fit and os.path.exists(save_folder):
                            save_subfolder = os.path.join(save_folder, f'Curve_fit_{fit_method}', exp_type+f'_{reference_type}')
                            if not os.path.exists(save_subfolder):
                                os.makedirs(save_subfolder)
                            fig.savefig(os.path.join(save_subfolder,f'exp_{_exp_index+1}_cell_{_cell_id}_chr_{_object_id}.png'),bbox_inches='tight')
                            plt.cla()
                            plt.close(fig)

    # report original time_interval for sanity check
    if report_time:
        print (f'-- The original time interval for each replicate is {_ori_time_interval} in (s).')

    # reset and drop index column
    _measurement_of_interest_all_chr = _measurement_of_interest_all_chr.reset_index()
    _measurement_of_interest_all_chr = _measurement_of_interest_all_chr.drop(columns=['index'])
    
    # sort numeric measurement columns
    non_numeric_cols =0
    for _col in _measurement_of_interest_all_chr.columns:
        if type(_col) is str:
            non_numeric_cols+=1

    sorted_num_cols = sorted(_measurement_of_interest_all_chr.columns[non_numeric_cols:],key=lambda x: x)
    sorted_cols=_measurement_of_interest_all_chr.columns[:non_numeric_cols].to_list()
    for _col in sorted_num_cols:
        sorted_cols.append(_col)
    _measurement_of_interest_all_chr=_measurement_of_interest_all_chr.reindex(columns=sorted_cols)

    if key_error>0:
        print (f'-- Warning: There are {key_error} KeyErrors for the dataframe reported, n\
               check if results were indeed processed correctly.')

            
    return _measurement_of_interest_all_chr





#####################################################################################

# function to complie all integrated and bin results generated from the function above
# Define specific parameters for measurement below


def compile_integrated_and_bins (exp_fds_dict,   # folder and NE region dict for each condition group
                                 measurement_column, 
                                 reference_type,
                                 measurement_type_int,     # function method #1 to call
                                 measurement_type_bin,     # function method #2 to call
                                normalize_column='er_recruit_normalized_area',
                                 analysis_save_folder=None,
                                 _save_results= True,
                                 force_time_step=0,
                                 ref_max=0.35,
                                 ref_min=0.1,
                                 _marker=None, 
                                 _cell_line=None,
                                 _one_bin_core=True,
                                 report_time =True,
                                 _use_curve_fit=False,
                                 fit_method = 'sigmoid',
                                 x0_adjust=1,
                                 sigmoid_ratio_th=0.05,
                                 #_nonnumeric_col=4,
                                ):
    
    from tqdm.notebook import trange, tqdm

    # Set kwargs dicts for different NE regions to use different normalize time function
    _bin_index = [0,1,2,3,4,5]
    _kwargs_dict ={}
    for _bin in _bin_index:
        # shared params
        _kwargs_dict[_bin]={}
        _kwargs_dict[_bin]['measurement_column']=measurement_column
        _kwargs_dict[_bin]['reference_type']=reference_type
        # set params accordingly
        # for integrated
        if _bin ==0:
            _kwargs_dict[_bin]['measurement_type']=measurement_type_int
        # for each bin
        else:
            _kwargs_dict[_bin]['measurement_type']=measurement_type_bin
            _kwargs_dict[_bin]['bin_of_interest']=_bin
        #print(_kwargs_dict[_bin])
        
    # Start analysis and compiling for all bins and the integrated                   
    _initial_measurement_df_dict = {}    
    # measure and save (optionally)  
    _analyze=True
    if _analyze:
        for _exp_name, _exp_f in tqdm(exp_fds_dict.items()):
            print (f'-- Analyzing {_exp_name} results')
            _initial_measurement_df_dict[_exp_name]={}
            for _kwarg_bin_index, _kwargs in _kwargs_dict.items():
                if not _use_curve_fit:
                    _measurement_df = normalize_time_for_mitotic_timelapsed_measurement (_exp_f, **_kwargs,
                                                                                     normalize_column=normalize_column,
                                                                                     force_time_step=force_time_step,
                                                                                    ref_max=ref_max,ref_min=ref_min,
                                                                                    report_time =report_time)
                elif _use_curve_fit: 
                    _measurement_df = curve_fit_normalize_time_for_mitotic_timelapsed_measurement (_exp_f, **_kwargs,
                                                                                     normalize_column=normalize_column,
                                                                                     save_fit=_save_results,
                                                                                     save_folder=analysis_save_folder,
                                                                                     exp_type= _exp_name,
                                                                                    report_time =report_time,fit_method = fit_method,
                                                                                    x0_adjust=x0_adjust,
                                                                                    sigmoid_ratio_th=sigmoid_ratio_th)                                                               
                
                _initial_measurement_df_dict[_exp_name][_kwarg_bin_index]=_measurement_df
    
    # extract and process to get the bin/integrated of interest
    # modify this to change the definition of 'core' and 'noncore'
    _processed_measurement_df_dict ={}
    for _exp_name_key in _initial_measurement_df_dict.keys():
        _processed_measurement_df_dict[_exp_name_key]={}
    
        # integrated
        _processed_measurement_df_dict[_exp_name_key]['integrated']=_initial_measurement_df_dict[_exp_name_key][0]
        
        # core using mean of bin2+bin3+bin4 or just bin3
        if _one_bin_core:
            _ave_core_measurement =_initial_measurement_df_dict[_exp_name_key][3]
        else:
            _ave_core_measurement = (_initial_measurement_df_dict[_exp_name_key][2]+
                                 _initial_measurement_df_dict[_exp_name_key][3]+
                                _initial_measurement_df_dict[_exp_name_key][4])/3
        
        _processed_measurement_df_dict[_exp_name_key]['core-center']=_ave_core_measurement

        # get all numeric measurement columns
        _nonnumeric_col =0
        for _col in _processed_measurement_df_dict[_exp_name_key]['core-center'].columns:
            if type(_col) is str:
                _nonnumeric_col+=1

        # to include the midpoint level column
        _col_offset=1
        _nonnumeric_col =_nonnumeric_col-_col_offset
        
        # noncore using mean of bin1+bin5
        _info_df = _processed_measurement_df_dict[_exp_name_key]['core-center'].iloc[:, :_nonnumeric_col]
        _ave_noncore_measurement = (_initial_measurement_df_dict[_exp_name_key][1].iloc[:, _nonnumeric_col:]+
                                    _initial_measurement_df_dict[_exp_name_key][5].iloc[:, _nonnumeric_col:])/2
        _processed_measurement_df_dict[_exp_name_key]['core-periphery']=_info_df.join(_ave_noncore_measurement)


        # core/noncore ratio
        #_info_df = _processed_measurement_df_dict[_exp_name_key]['core-center'].iloc[:, :_nonnumeric_col]
        _first_df = _processed_measurement_df_dict[_exp_name_key]['core-center'].iloc[:, _nonnumeric_col:]
        _second_df = _processed_measurement_df_dict[_exp_name_key]['core-periphery'].iloc[:, _nonnumeric_col:]
        _processed_measurement_df_dict[_exp_name_key]['center-periphery']=_info_df.join(_first_df/_second_df)
         
    print ('-- Analysis complete.')
    
    # save results as excels
    if _save_results and os.path.exists(analysis_save_folder):

        pickle_prefix = ''
        if type(_marker) is str and type(_cell_line) is str:
            pickle_prefix =   _cell_line + '_' + _marker + '_' + pickle_prefix
        pickle_fname = analysis_save_folder + os.sep + pickle_prefix + f'chromosome_inner_core_{fit_method}.pickle'
        with open(pickle_fname, 'wb') as handle:
            pickle.dump(_processed_measurement_df_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        

        analysis_save_folder = os.path.join(analysis_save_folder,f'Chromosome_inner_core_analysis_{fit_method}')  
        if not os.path.exists(analysis_save_folder):
            os.mkdir(analysis_save_folder)
        print ('-- Saving analysis result.')
        
        for _exp_name, _exp_measure in _processed_measurement_df_dict.items(): 
            for _measure_name, _measure_df in _exp_measure.items():
                _save_name =_exp_name+'_'+_measure_name
                
                if type(_marker) is str and type(_cell_line) is str:
                    _save_name = _cell_line + '_' + _marker + '_'+_save_name
                    
                _save_name = _save_name +'.xlsx'
                _measure_df.to_excel(analysis_save_folder+ os.sep+ _save_name, index=False)
                
    return _processed_measurement_df_dict


###################################################################################

# function to complie all cellular subregion results generated from the function above
# define specific parameters for measurement below


def compile_global_subregions (exp_fds_dict,   # folder and NE region dict for each condition group
                                 measurement_column, 
                                 reference_type,
                                 measurement_type_global_1,     # function method #1 to call
                                 measurement_type_global_2,     # function method #2 to call
                                 normalize_column='er_recruit_normalized_area',
                                 analysis_save_folder=None,
                                 _save_results= True,
                                 force_time_step=0,
                                 ref_max=0.35,
                                 ref_min=0.1,
                                 _marker=None, 
                                 _cell_line=None,
                                  report_time =True,
                                 _use_curve_fit=False,
                                 fit_method = 'sigmoid',
                                 sigmoid_ratio_th=0.05,
                                 x0_adjust=1,
                                 #_one_bin_core=True,
                                 #_nonnumeric_col=4,
                                ):
    
    from tqdm.notebook import trange, tqdm

    # Set kwargs dicts for different cellular regions to use different normalize time function
    # use _bin_index similarly as above but to indicate different subregions within the cell
    _bin_index = ['cortical','endo','midzone','dna']
    _kwargs_dict ={}
    for _bin in _bin_index:
        # shared params
        _kwargs_dict[_bin]={}
        _kwargs_dict[_bin]['reference_type']=reference_type
        # set params accordingly
        # for integrated
        if _bin !='dna':
            _kwargs_dict[_bin]['measurement_type']=measurement_type_global_1
            _kwargs_dict[_bin]['measurement_column']=_bin + '_' + measurement_column
        # for each bin
        else:
            _kwargs_dict[_bin]['measurement_type']=measurement_type_global_2
            _kwargs_dict[_bin]['measurement_column']='dna_endo' + '_' + measurement_column
        
        
    # Start analysis and compiling for all bins and the integrated                   
    _initial_measurement_df_dict = {}    
    # measure and save (optionally)  
    _analyze=True
    if _analyze:
        for _exp_name, _exp_f in tqdm(exp_fds_dict.items()):
            print (f'-- Analyzing {_exp_name} results')
            _initial_measurement_df_dict[_exp_name]={}
            for _kwarg_bin_index, _kwargs in _kwargs_dict.items():
                if not _use_curve_fit:
                    _measurement_df = normalize_time_for_mitotic_timelapsed_measurement (_exp_f, **_kwargs,
                                                                                     normalize_column=normalize_column,
                                                                                     force_time_step=force_time_step,
                                                                                    ref_max=ref_max,ref_min=ref_min,
                                                                                    report_time =report_time)


                elif _use_curve_fit: 
                    _measurement_df = curve_fit_normalize_time_for_mitotic_timelapsed_measurement (_exp_f, **_kwargs,
                                                                                     normalize_column=normalize_column,
                                                                                     save_fit=_save_results,
                                                                                     save_folder=analysis_save_folder,
                                                                                     exp_type= _exp_name,fit_method = fit_method,
                                                                                     x0_adjust=x0_adjust,
                                                                                    report_time =report_time,sigmoid_ratio_th=sigmoid_ratio_th)  
                
                _initial_measurement_df_dict[_exp_name][_kwarg_bin_index]=_measurement_df

    
    # extract and process to get the bin/integrated of interest
    # modify this to change the definition of 'core' and 'noncore'
    _processed_measurement_df_dict ={}
    for _exp_name_key in _initial_measurement_df_dict.keys():
        _processed_measurement_df_dict[_exp_name_key]={}
    
        for _bin_index, _bin_measure in _initial_measurement_df_dict[_exp_name_key].items():
            _processed_measurement_df_dict[_exp_name_key][_bin_index]=_bin_measure


        # get all numeric measurement columns
        _nonnumeric_col =0
        for _col in _processed_measurement_df_dict[_exp_name_key]['endo'].columns:
            if type(_col) is str:
                _nonnumeric_col+=1

        # to include the midpoint level column
        _col_offset=1
        _nonnumeric_col =_nonnumeric_col-_col_offset

        # endo-cortical ratio
        _info_df = _processed_measurement_df_dict[_exp_name_key]['endo'].iloc[:, :_nonnumeric_col]
        _first_df = _processed_measurement_df_dict[_exp_name_key]['endo'].iloc[:, _nonnumeric_col:]
        _second_df = _processed_measurement_df_dict[_exp_name_key]['cortical'].iloc[:, _nonnumeric_col:]
        _processed_measurement_df_dict[_exp_name_key]['endo-cortical']=_info_df.join(_first_df/_second_df)
        
         
    print ('-- Analysis complete.')
    
    # save results as excels
    if _save_results and os.path.exists(analysis_save_folder):

        pickle_prefix = ''
        if type(_marker) is str and type(_cell_line) is str:
            pickle_prefix =   _cell_line + '_' + _marker + '_' + pickle_prefix
        pickle_fname = analysis_save_folder + os.sep + pickle_prefix + f'cell_global_{fit_method}.pickle'
        with open(pickle_fname, 'wb') as handle:
            pickle.dump(_processed_measurement_df_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

        analysis_save_folder = os.path.join(analysis_save_folder,f'Cell_global_analysis_{fit_method}')
        if not os.path.exists(analysis_save_folder):
            os.mkdir(analysis_save_folder)
        print ('-- Saving analysis result.')
        for _exp_name, _exp_measure in _processed_measurement_df_dict.items(): 
            for _measure_name, _measure_df in _exp_measure.items():
                _save_name =_exp_name+'_'+_measure_name
                
                if type(_marker) is str and type(_cell_line) is str:
                    _save_name = _cell_line + '_' + _marker + '_'+_save_name
                    
                _save_name = _save_name +'.xlsx'
                _measure_df.to_excel(analysis_save_folder+ os.sep+ _save_name, index=False)
                
    return _processed_measurement_df_dict





    
###################################################################################

# function to complie all chromosome (core and noncore) subregion results generated from the function above
# define specific parameters for measurement below

# function to complie all experiment results generated from the function above
# define specific parameters for measurement below


def compile_core_and_noncore (exp_fds_dict,   # folder and NE region dict for each condition group
                                 measurement_column, 
                                 reference_type,
                                 measurement_type_int,     # function method #1 to call
                                 measurement_type_bin,     # function method #2 to call
                                normalize_column='er_recruit_normalized_area',
                                 analysis_save_folder=None,
                                 _save_results= True,
                                 force_time_step=0,
                                 ref_max=0.35,
                                 ref_min=0.1,
                                 _marker=None, 
                                 _cell_line=None,
                                 _one_bin_core=True,
                                 report_time =True,
                                 _use_curve_fit=False,
                                 fit_method = 'sigmoid',
                                 x0_adjust=1,
                                 sigmoid_ratio_th=0.05,
                                 #_nonnumeric_col=4,
                                ):
    
    from tqdm.notebook import trange, tqdm

    # Set kwargs dicts for different chromosome regions to use different normalize time function
    _bin_index = ['inner_core','outer_core','noncore','whole_chromosome']
    _kwargs_dict ={}
    for _bin in _bin_index:
        # shared params
        _kwargs_dict[_bin]={}
        _kwargs_dict[_bin]['reference_type']=reference_type
        _kwargs_dict[_bin]['measurement_column']=measurement_column
        # set params accordingly
        # for integrated
        if _bin !='whole_chromosome':
            _kwargs_dict[_bin]['measurement_type']=measurement_type_bin
            # provide core-noncore type as bin of interest
            _kwargs_dict[_bin]['bin_of_interest']=_bin
            
        # for each bin
        else:
            _kwargs_dict[_bin]['measurement_type']=measurement_type_int
            
        
    # Start analysis and compiling for all subregions and the whole chromosome                 
    _initial_measurement_df_dict = {}    
    # measure and save (optionally)  
    _analyze=True
    if _analyze:
        for _exp_name, _exp_f in tqdm(exp_fds_dict.items()):
            print (f'-- Analyzing {_exp_name} results')
            _initial_measurement_df_dict[_exp_name]={}
            for _kwarg_bin_index, _kwargs in _kwargs_dict.items():

                if not _use_curve_fit:
                    _measurement_df = normalize_time_for_mitotic_timelapsed_measurement (_exp_f, **_kwargs,
                                                                                     normalize_column=normalize_column,
                                                                                     force_time_step=force_time_step,
                                                                                    ref_max=ref_max,ref_min=ref_min,
                                                                                    report_time =report_time)

                elif _use_curve_fit: 
                    _measurement_df = curve_fit_normalize_time_for_mitotic_timelapsed_measurement (_exp_f, **_kwargs,
                                                                                     normalize_column=normalize_column,
                                                                                     save_fit=_save_results,
                                                                                     save_folder=analysis_save_folder,
                                                                                     exp_type= _exp_name,fit_method = fit_method,
                                                                                     x0_adjust=x0_adjust,
                                                                                    report_time =report_time,sigmoid_ratio_th=sigmoid_ratio_th)  
                
                _initial_measurement_df_dict[_exp_name][_kwarg_bin_index]=_measurement_df
    
    # extract and process to get the bin/integrated of interest
    # modify this to change the definition of 'core' and 'noncore'
    _processed_measurement_df_dict ={}
    for _exp_name_key in _initial_measurement_df_dict.keys():
        _processed_measurement_df_dict[_exp_name_key]={}
    
        for _bin_index, _bin_measure in _initial_measurement_df_dict[_exp_name_key].items():
            _processed_measurement_df_dict[_exp_name_key][_bin_index]=_bin_measure


        # get all numeric measurement columns
        _nonnumeric_col =0
        
        for _col in _processed_measurement_df_dict[_exp_name_key]['inner_core'].columns:
            if type(_col) is str:
                _nonnumeric_col+=1
        
        # to include the midpoint level column
        _col_offset=1
        _nonnumeric_col =_nonnumeric_col-_col_offset

        # inner_core-noncore ratio  
        _info_df = _processed_measurement_df_dict[_exp_name_key]['inner_core'].iloc[:, :_nonnumeric_col]
        _first_df = _processed_measurement_df_dict[_exp_name_key]['inner_core'].iloc[:, _nonnumeric_col:]
        _second_df = _processed_measurement_df_dict[_exp_name_key]['noncore'].iloc[:, _nonnumeric_col:]
        _processed_measurement_df_dict[_exp_name_key]['core-noncore']=_info_df.join(_first_df/_second_df)
         
    print ('-- Analysis complete.')
    
    # save results as excels
    if _save_results and os.path.exists(analysis_save_folder):

        pickle_prefix = ''
        if type(_marker) is str and type(_cell_line) is str:
            pickle_prefix =  _cell_line + '_' + _marker + '_' + pickle_prefix
        pickle_fname = analysis_save_folder + os.sep + pickle_prefix + f'core_noncore_{fit_method}.pickle'
        with open(pickle_fname, 'wb') as handle:
            pickle.dump(_processed_measurement_df_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

        analysis_save_folder = os.path.join(analysis_save_folder,f'Core_noncore_analysis_{fit_method}')
        if not os.path.exists(analysis_save_folder):
            os.mkdir(analysis_save_folder)
        print ('-- Saving analysis result.')
        for _exp_name, _exp_measure in _processed_measurement_df_dict.items(): 
            for _measure_name, _measure_df in _exp_measure.items():
                _save_name =_exp_name+'_'+_measure_name
                
                if type(_marker) is str and type(_cell_line) is str:
                    _save_name = _cell_line + '_' + _marker + '_'+_save_name
                    
                _save_name = _save_name +'.xlsx'
                _measure_df.to_excel(analysis_save_folder+ os.sep+ _save_name, index=False)
                
    return _processed_measurement_df_dict



###############################################
# for AL analysis
def compile_AL_subregions (exp_fds_dict,   # folder and NE region dict for each condition group
                                 measurement_column, 
                                 reference_type,
                                 measurement_type_global_1,     # function method #1 to call
                                 measurement_type_global_2,     # function method #2 to call
                                 normalize_column='er_recruit_normalized_area',
                                 analysis_save_folder=None,
                                 _save_results= True,
                                 force_time_step=0,
                                 ref_max=0.35,
                                 ref_min=0.1,
                                 _marker=None, 
                                 _cell_line=None,
                                  report_time =True,
                                 _use_curve_fit=False,
                                 fit_method = 'sigmoid',
                                 sigmoid_ratio_th=0.05,
                                 x0_adjust=1,
                                 #_one_bin_core=True,
                                 #_nonnumeric_col=4,
                                ):
    
    from tqdm.notebook import trange, tqdm

    # Set kwargs dicts for different cellular regions to use different normalize time function
    # use _bin_index similarly as above but to indicate different subregions within the cell
    _bin_index = ['cell_endo','main_endo','midzone_endo','main_cortical','midzone_cortical']
    _kwargs_dict ={}
    for _bin in _bin_index:
        # shared params
        _kwargs_dict[_bin]={}
        _kwargs_dict[_bin]['reference_type']=reference_type
        # set params accordingly
        # for integrated
        if _bin !='dna':
            _kwargs_dict[_bin]['measurement_type']=measurement_type_global_1
            _kwargs_dict[_bin]['measurement_column']=_bin + '_' + measurement_column
        # for each bin
        else:
            _kwargs_dict[_bin]['measurement_type']=measurement_type_global_2
            _kwargs_dict[_bin]['measurement_column']='dna_endo' + '_' + measurement_column
        
        
    # Start analysis and compiling for all bins and the integrated                   
    _initial_measurement_df_dict = {}    
    # measure and save (optionally)  
    _analyze=True
    if _analyze:
        for _exp_name, _exp_f in tqdm(exp_fds_dict.items()):
            print (f'-- Analyzing {_exp_name} results')
            _initial_measurement_df_dict[_exp_name]={}
            for _kwarg_bin_index, _kwargs in _kwargs_dict.items():
                if not _use_curve_fit:
                    _measurement_df = normalize_time_for_mitotic_timelapsed_measurement (_exp_f, **_kwargs,
                                                                                     normalize_column=normalize_column,
                                                                                     force_time_step=force_time_step,
                                                                                    ref_max=ref_max,ref_min=ref_min,
                                                                                    report_time =report_time)


                elif _use_curve_fit: 
                    _measurement_df = curve_fit_normalize_time_for_mitotic_timelapsed_measurement (_exp_f, **_kwargs,
                                                                                     normalize_column=normalize_column,
                                                                                     save_fit=_save_results,
                                                                                     save_folder=analysis_save_folder,
                                                                                     exp_type= _exp_name,fit_method = fit_method,
                                                                                     x0_adjust=x0_adjust,
                                                                                    report_time =report_time,sigmoid_ratio_th=sigmoid_ratio_th)  
                
                _initial_measurement_df_dict[_exp_name][_kwarg_bin_index]=_measurement_df

    
    # extract and process to get the bin/integrated of interest
    # modify this to change the definition of 'core' and 'noncore'
    _processed_measurement_df_dict ={}
    for _exp_name_key in _initial_measurement_df_dict.keys():
        _processed_measurement_df_dict[_exp_name_key]={}
    
        for _bin_index, _bin_measure in _initial_measurement_df_dict[_exp_name_key].items():
            _processed_measurement_df_dict[_exp_name_key][_bin_index]=_bin_measure


        # get all numeric measurement columns
        _nonnumeric_col =0
        for _col in _processed_measurement_df_dict[_exp_name_key]['main_endo'].columns:
            if type(_col) is str:
                _nonnumeric_col+=1

        # to include the midpoint level column
        _col_offset=1
        _nonnumeric_col =_nonnumeric_col-_col_offset

        # endo-cortical ratio
        _info_df = _processed_measurement_df_dict[_exp_name_key]['main_endo'].iloc[:, :_nonnumeric_col]
        _first_df = _processed_measurement_df_dict[_exp_name_key]['main_endo'].iloc[:, _nonnumeric_col:]
        _second_df = _processed_measurement_df_dict[_exp_name_key]['midzone_endo'].iloc[:, _nonnumeric_col:]
        #_processed_measurement_df_dict[_exp_name_key]['endo-cortical']=_info_df.join(_first_df/_second_df)
        
         
    print ('-- Analysis complete.')
    
    # save results as excels
    if _save_results and os.path.exists(analysis_save_folder):

        pickle_prefix = ''
        if type(_marker) is str and type(_cell_line) is str:
            pickle_prefix =   _cell_line + '_' + _marker + '_' + pickle_prefix
        pickle_fname = analysis_save_folder + os.sep + pickle_prefix + f'cell_AL_{fit_method}.pickle'
        with open(pickle_fname, 'wb') as handle:
            pickle.dump(_processed_measurement_df_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

        analysis_save_folder = os.path.join(analysis_save_folder,f'Cell_AL_analysis_{fit_method}')
        if not os.path.exists(analysis_save_folder):
            os.mkdir(analysis_save_folder)
        print ('-- Saving analysis result.')
        for _exp_name, _exp_measure in _processed_measurement_df_dict.items(): 
            for _measure_name, _measure_df in _exp_measure.items():
                _save_name =_exp_name+'_'+_measure_name
                
                if type(_marker) is str and type(_cell_line) is str:
                    _save_name = _cell_line + '_' + _marker + '_'+_save_name
                    
                _save_name = _save_name +'.xlsx'
                _measure_df.to_excel(analysis_save_folder+ os.sep+ _save_name, index=False)
                
    return _processed_measurement_df_dict


###################################################################################
# function to get the measurement at the timepoint/stage of interest using the curve fit of the given reference
# function intented for same type of experiment and can use the existing x_list calculated elsewhere for this same type experiment

def midpoint_measurement_from_processed_measurement_df_dict (_processed_measurement_df_dict, 
                                                             ref_type = None, 
                                                             save_results = True,
                                                             save_folder =None,
                                                             midpoint_of_Y_range=[],
                                                             midpoint_of_X_coord_dict={},
                                                             use_y_sigma=False,
                                                             y_sel_sigmoid_ratio=0.25,
                                                             curve_fit_using_median=False,
                                                             #sigmoid_ratio_th=0.05,
                                                             _marker =None,
                                                             fit_method='sigmoid',
                                                             _cell_line =None,
                                                              ):


    from scipy.optimize import curve_fit

    # 0. Get reference type name from defualt or from the given input
    if ref_type == None:
        # use whole chromosome as the first option
        for _val in  _processed_measurement_df_dict.values():
            if 'whole_chromosome' in _val.keys():
                ref_type = 'whole_chromosome'
            elif 'integrated' in _val.keys():
                ref_type = 'integrated'
            else:
                print ('No valid reference can be found. Exit.')
                return None
    else:
        pass

    # 1. Get midpoint (global) Y from the reference population or from the given input
    if len(midpoint_of_Y_range) == len(_processed_measurement_df_dict.values()):
        Y_range_len = [len(_yrange) for _yrange in midpoint_of_Y_range]
        if np.sum(np.array(Y_range_len)==2)==len(midpoint_of_Y_range):
            print (f'-- Using the given (global) Y to select measurement of interest.')
        else:
            print (f'-- The given (global) Y is not valid for the measurements.')
    elif len(midpoint_of_Y_range) >0:
        print (f'-- The given (global) Y is not valid for the measurements.')
    elif len(midpoint_of_Y_range)==0:
        print (f'-- Estimating the (global) Y from the reference population to select measurement of interest.')

        for _exp_name, _exp_measure in _processed_measurement_df_dict.items():

            if ref_type in _exp_measure.keys():
                _nonnumeric_col=0
                for _col in _exp_measure[ref_type].columns:
                    if type(_col) is str:
                        _nonnumeric_col+=1
                
                # estimate the global Y for each exp using the reference population
                xydata = _exp_measure[ref_type].iloc[:, _nonnumeric_col:]
                xdata =[]
                ydata=[]
                for _cell in xydata.iloc():
                    _cell=_cell.dropna()
                    _x_cell = _cell.index.to_list()
                    _y_cell = _cell.to_list()
                    if not curve_fit_using_median:
                        for x,y, in zip(_x_cell,_y_cell):
                            xdata.append(x)
                            ydata.append(y)
                    else:
                        xdata = np.unique(xydata.columns.to_list())
                        ydata=[]
                        for _x in xdata:
                            _y_each_x = np.nanmedian(xydata[_x])
                            ydata.append(_y_each_x)

                p0=[max(ydata)-min(ydata), np.median(np.unique(xdata)),1,min(ydata)]
                try:
                    popt, pcov = curve_fit(sigmoid, xdata, ydata,p0, method='dogbox')
                    L ,x0, k, b= popt

                    # also get xy for plotting
                    x_plot = np.unique(xdata)
                    y_plot = sigmoid(x_plot,L ,x0, k, b)
                    fig, ax = plt.subplots(figsize=(3,2), dpi=300)
                    ax = plt.scatter(xdata,ydata,s=5)
                    ax = plt.plot(x_plot, y_plot, color='red')
                    if save_results:
                        # save curve fit figure if desired
                        if os.path.exists(save_folder):
                            save_subfolder = os.path.join(save_folder, f'Curve_fit_{fit_method}', f'{_exp_name}_{ref_type}')
                            if not os.path.exists(save_subfolder):
                                os.makedirs(save_subfolder)
                            fig.savefig(os.path.join(save_subfolder,f'{_marker}_exp_all_cells.png'),bbox_inches='tight')
                    # close figure in plt
                    plt.cla()
                    plt.close(fig)

                    # if get y-mid from the fit
                    if use_y_sigma:
                        y_mid = sigmoid(x0,L ,x0, k, b)
                    # if get y-mid from the median of all y corresponding to the x-mid from the fit
                    else:
                        x_mid = determine_time_interval(x0, xdata)
                        y_mid = np.nanmedian(xydata[x_mid])
                    # estimate the near-platue with specified scale of change
                    #_time_step = abs(xdata[0]-xdata[1])
                    y_max=sigmoid(10000,L ,x0, k, b)
                    # append the estimated mid_y and max_y for each experiment    
                    midpoint_of_Y_range.append([y_mid,y_max])
                except RuntimeError:
                    midpoint_of_Y_range.append([np.nan,np.nan])

            else:
                print ('No valid reference can be found. Exit.')
                return None


    # Get selected measurement
    sel_measurement_df_dict = {}

    # 2.1 Use the given X coord to directly select 
    if midpoint_of_X_coord_dict.keys()==_processed_measurement_df_dict.keys():
        sel_x_list_dict = midpoint_of_X_coord_dict

        for _exp_index, (_exp_name, _exp_measure) in enumerate(_processed_measurement_df_dict.items()):
            sel_measurement_df_dict[_exp_name]={}
            
            print ('Use the given X coord from the corresponding experiment to select.')
            sel_x_list = midpoint_of_X_coord_dict[_exp_name]

            for _measure_name, _measure_df in _exp_measure.items(): 
                sel_measurement_df_dict[_exp_name][_measure_name] = []

                for _cell, _x in zip(_measure_df.iloc(), sel_x_list):
                    if _x is not np.nan and _x in _measure_df.columns:
                        sel_measurement_df_dict[_exp_name][_measure_name].append(_cell[_x])
                    else:
                        sel_measurement_df_dict[_exp_name][_measure_name].append(np.nan)

    # 2.2 or use midpoint (global) Y from the reference population to select the measuremnt of interest
    else:
        sel_x_list_dict = {}
        for _exp_index, (_exp_name, _exp_measure) in enumerate(_processed_measurement_df_dict.items()):
            sel_measurement_df_dict[_exp_name]={}

            # 2.1 Get the X coord for measurement of interest
            sel_x_list = []

            if ref_type in _exp_measure.keys():
                _nonnumeric_col=0
                for _col in _exp_measure[ref_type].columns:
                    if type(_col) is str:
                        _nonnumeric_col+=1
                
                if np.nan in midpoint_of_Y_range[_exp_index]:
                    print ("No valid Y range can be used. Skip this experiment dataset.")
                    for _cell in xydata.iloc():
                        sel_x_list.append(np.nan)
                else:
                    y_mid, y_max =midpoint_of_Y_range[_exp_index]
                    xydata = _exp_measure[ref_type].iloc[:, _nonnumeric_col:]
                    for _cell in xydata.iloc():
                        _diff_ratio_list=[]
                        for _x, _y in zip(_cell.dropna().index.to_list(),_cell.dropna().to_list()):
                            _diff_ratio_list.append(abs(_y-y_mid)/(y_max-y_mid))
                        _diff_ratio_list =np.array(_diff_ratio_list)
                        if np.min(_diff_ratio_list) <= y_sel_sigmoid_ratio:
                            #cand_y = _cell.dropna().to_list()[np.argmin(_diff_ratio_list)]
                            cand_x = _cell.dropna().index.to_list()[np.argmin(_diff_ratio_list)]
                            #sel_y_list.append(cand_y)
                            sel_x_list.append(cand_x)
                        else:
                            #sel_y_list.append(np.nan)
                            sel_x_list.append(np.nan)

            else:
                print ('No valid matching reference can be found. Exit.')
                return None
        
            sel_x_list_dict[_exp_name]=sel_x_list
        
        
            # 2.2 Get the Y measurement of interest using the X coord 
            for _measure_name, _measure_df in _exp_measure.items(): 
                sel_measurement_df_dict[_exp_name][_measure_name] = []

                for _cell, _x in zip(_measure_df.iloc(), sel_x_list):
                    if _x is not np.nan:
                        sel_measurement_df_dict[_exp_name][_measure_name].append(_cell[_x])
                    else:
                        sel_measurement_df_dict[_exp_name][_measure_name].append(np.nan)
    
    # 3. Save results if needed
    if save_results and os.path.exists(save_folder):
        analysis_save_folder = os.path.join(save_folder,f'Selected_midpoint_measurement_{fit_method}')
        if not os.path.exists(analysis_save_folder):
            os.mkdir(analysis_save_folder)
        print ('-- Saving analysis result.')
        # pickle the results
        if midpoint_of_X_coord_dict.keys()==_processed_measurement_df_dict.keys():
            pickle_prefix = f'midpoint_using_{ref_type}.pickle'
            pickle_x_prefix = f'midpoint_x_coord_using_{ref_type}.pickle'
        else:
            pickle_prefix = f'midpoint_using_{ref_type}_self.pickle'
            pickle_x_prefix = f'midpoint_x_coord_using_{ref_type}_self.pickle'
        if type(_marker) is str and type(_cell_line) is str:
            pickle_prefix =  _cell_line + '_' + _marker + '_' + pickle_prefix
            pickle_x_prefix = _cell_line + '_' + _marker + '_' + pickle_x_prefix
        pickle_fname = analysis_save_folder + os.sep + pickle_prefix
        pickle_x_fname = analysis_save_folder + os.sep + pickle_x_prefix
        with open(pickle_fname, 'wb') as handle:
            pickle.dump(sel_measurement_df_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(pickle_x_fname, 'wb') as handle:
            pickle.dump(sel_x_list_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        # also save excels
        for _exp_name, _sel_measure in sel_measurement_df_dict.items():
            df = pd.DataFrame(_sel_measure)
            if midpoint_of_X_coord_dict.keys()==_processed_measurement_df_dict.keys():
                _save_name =_exp_name+'_'+f'midpoint_using_{ref_type}.xlsx'
            else:
                _save_name =_exp_name+'_'+f'midpoint_using_{ref_type}_self.xlsx'
            if type(_marker) is str and type(_cell_line) is str:
                _save_name = _cell_line + '_' + _marker + '_'+_save_name
                df.to_excel(analysis_save_folder+ os.sep+ _save_name, index=False)
            
    return sel_measurement_df_dict, sel_x_list_dict

 



 ###################################################################################
# function to get the measurement at the timepoint/stage of interest using the curve fit of the given reference
# function intented for same type of experiment and can use the existing x_list calculated elsewhere for this same type experiment

def global_midpoint_measurement_from_processed_measurement_df_dict (_processed_measurement_df_dict, 
                                                             ref_type = None, 
                                                             save_results = True,
                                                             save_folder =None,
                                                             curve_fit_using_median=False,
                                                             #sigmoid_ratio_th=0.05,
                                                             #mid_from_end_th =0,
                                                             midpoint_of_X_coord_dict={},
                                                             _marker =None,
                                                             fit_method='sigmoid',
                                                             _cell_line =None,
                                                              ):


    from scipy.optimize import curve_fit
    from sklearn.metrics import r2_score

    # 0. Get reference type name from defualt or from the given input
    if ref_type == None:
        # use whole chromosome as the first option
        for _val in  _processed_measurement_df_dict.values():
            if 'whole_chromosome' in _val.keys():
                ref_type = 'whole_chromosome'
            elif 'integrated' in _val.keys():
                ref_type = 'integrated'
            else:
                print ('No valid reference can be found. Exit.')
                return None
    else:
        pass

 
    # Use the given X coord to directly select 
    if midpoint_of_X_coord_dict.keys()==_processed_measurement_df_dict.keys():
        sel_x_dict = midpoint_of_X_coord_dict
        sel_measurement_df_dict = {}

        for _exp_index, (_exp_name, _exp_measure) in enumerate(_processed_measurement_df_dict.items()):
            sel_measurement_df_dict[_exp_name]={}
            
            print ('Use the given X coord from the corresponding experiment to select.')
            sel_x = midpoint_of_X_coord_dict[_exp_name]

            for _measure_name, _measure_df in _exp_measure.items(): 
                sel_measurement_df_dict[_exp_name][_measure_name] = []

                for _cell in _measure_df.iloc:
                    if sel_x is not np.nan and sel_x in _measure_df.columns:
                        sel_measurement_df_dict[_exp_name][_measure_name].append(_cell[sel_x])
                    else:
                        sel_measurement_df_dict[_exp_name][_measure_name].append(np.nan)
    
    # Estimate x for each experiment and get corresponding measurement for the x
    else:
        # dict to store the x coord
        sel_x_dict = {}
        # selected measurement dict to store the result
        sel_measurement_df_dict = {}

        for _exp_name, _exp_measure in _processed_measurement_df_dict.items():

            sel_measurement_df_dict[_exp_name] ={}

            if ref_type in _exp_measure.keys():
                _nonnumeric_col=0
                for _col in _exp_measure[ref_type].columns:
                    if type(_col) is str:
                        _nonnumeric_col+=1
                
                # estimate the global X midpoint for each exp using the reference population
                xydata = _exp_measure[ref_type].iloc[:, _nonnumeric_col:]
                xdata =[]
                ydata=[]
                for _cell in xydata.iloc():
                    _cell=_cell.dropna()
                    _x_cell = _cell.index.to_list()
                    _y_cell = _cell.to_list()
                    if not curve_fit_using_median:
                        for x,y, in zip(_x_cell,_y_cell):
                            xdata.append(x)
                            ydata.append(y)
                    else:
                        xdata = np.unique(xydata.columns.to_list())
                        ydata=[]
                        for _x in xdata:
                            _y_each_x = np.nanmedian(xydata[_x])
                            ydata.append(_y_each_x)
                
                fig, ax = plt.subplots(figsize=(3,2), dpi=300)
                ax = plt.scatter(xdata,ydata,s=5)
                # linear fit
                p = np.polyfit(xdata,ydata,1)
                r2score_1 = r2_score(ydata, np.polyval(p, xdata))
                ax = plt.plot(xdata,np.polyval(p, xdata), color='orange')
                # sigmoid fit
                p0_sig=[max(ydata)-min(ydata), np.median(np.unique(xdata)),1,min(ydata)]
                try:
                    if len(xdata) >= len(p0_sig):
                        popt, pcov = curve_fit(sigmoid, xdata, ydata,p0_sig)
                    else:
                        popt, pcov = curve_fit(sigmoid, xdata, ydata,p0_sig, method='dogbox')
                    L ,x0, k, b= popt
                    r2score_2 = r2_score(ydata, sigmoid(xdata,L ,x0, k, b))
                    # also get xy for plotting
                    x_plot = np.unique(xdata)
                    x_plot = np.arange(x_plot[0],x_plot[-1]+0.2,0.2)
                    y_plot = sigmoid(x_plot,L ,x0, k, b)
                    ax = plt.plot(x_plot, y_plot, color='red')
                except RuntimeError:
                    print ('Sigmoid fit fails. Skip.')
                    r2score_2 = 0
                # expoenatial fit
                p0_exp = p0 = [max(ydata),min(ydata),1]
                try:
                    if len(xdata) >= len(p0_exp):
                        popt, pcov = curve_fit(exponential_plateau, xdata, ydata,p0_exp)
                    else:
                        popt, pcov = curve_fit(exponential_plateau, xdata, ydata,p0_exp, method='dogbox')
                    YM, Y0, k= popt
                    x0_exp = np.log((YM- (YM-Y0)/2)/(YM-Y0))/-k 

                    y_fit = [exponential_plateau(_x,YM, Y0, k) for _x in xdata]
                    r2score_3 = r2_score(ydata, y_fit)
                    x_plot = np.unique(xdata)
                    x_plot = np.arange(x_plot[0],x_plot[-1]+0.2,0.2)
                    y_plot = [exponential_plateau(_x,YM, Y0, k) for _x in x_plot]
                    ax = plt.plot(x_plot, y_plot, color='green')
                except RuntimeError:
                    print ('Exponential plateau fit fails. Skip.')
                    r2score_3=0

                ax = plt.annotate("linear:{:.3f}".format(r2score_1)+ "\n" + "sigmoid:{:.3f}".format(r2score_2)
                + "\n" + "exp_plateau:{:.3f}".format(r2score_3),
                                     (xdata[0],max(ydata)),fontsize=5)

                if save_results:
                    # save curve fit figure if desired
                    if os.path.exists(save_folder):
                        save_subfolder = os.path.join(save_folder, f'Curve_fit_{fit_method}', f'{_exp_name}_{ref_type}')
                        if not os.path.exists(save_subfolder):
                            os.makedirs(save_subfolder)
                        fig.savefig(os.path.join(save_subfolder,f'{_marker}_exp_all_cells.png'),bbox_inches='tight')
                # close figure in plt
                plt.cla()
                plt.close(fig)

                fit_sel_index = np.argmax([r2score_1,r2score_2,r2score_3])
                if fit_sel_index==0:
                    x0=np.unique(xdata)[-1]/2
                elif fit_sel_index==1:
                    #if r2score_2- r2score_3 > 0.1:
                        #x0=x0
                    # use exp fit instead to avoid the x0 to be at the start
                    #else:
                        #x0=x0_exp
                    x0=x0
                elif fit_sel_index==2:
                    x0=x0_exp
                
                _proceed=True
                if _proceed:
                    sel_x = determine_time_interval(x0,np.unique(xdata))
                    sel_x_dict[_exp_name]=sel_x

                    for _measure_name, _measure_df in _exp_measure.items(): 
                        sel_measurement_df_dict[_exp_name][_measure_name] = []
                        for _cell in _measure_df.iloc:
                            if sel_x is not np.nan and sel_x in _measure_df.columns:
                                sel_measurement_df_dict[_exp_name][_measure_name].append(_cell[sel_x])
                            else:
                                sel_measurement_df_dict[_exp_name][_measure_name].append(np.nan)
    
    # 3. Save results if needed
    if save_results and os.path.exists(save_folder):
        analysis_save_folder = os.path.join(save_folder,f'Global_midpoint_measurement_{fit_method}')
        if not os.path.exists(analysis_save_folder):
            os.mkdir(analysis_save_folder)
        print ('-- Saving analysis result.')
        # pickle the results
        if midpoint_of_X_coord_dict.keys()==_processed_measurement_df_dict.keys():
            pickle_prefix = f'midpoint_using_{ref_type}.pickle'
            pickle_x_prefix = f'midpoint_x_coord_using_{ref_type}.pickle'
        else:
            pickle_prefix = f'midpoint_using_{ref_type}_self.pickle'
            pickle_x_prefix = f'midpoint_x_coord_using_{ref_type}_self.pickle'
        if type(_marker) is str and type(_cell_line) is str:
            pickle_prefix =  _cell_line + '_' + _marker + '_' + pickle_prefix
            pickle_x_prefix = _cell_line + '_' + _marker + '_' + pickle_x_prefix
        pickle_fname = analysis_save_folder + os.sep + pickle_prefix
        pickle_x_fname = analysis_save_folder + os.sep + pickle_x_prefix
        with open(pickle_fname, 'wb') as handle:
            pickle.dump(sel_measurement_df_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open(pickle_x_fname, 'wb') as handle:
            pickle.dump(sel_x_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        # also save excels
        for _exp_name, _sel_measure in sel_measurement_df_dict.items():
            df = pd.DataFrame(_sel_measure)
            if midpoint_of_X_coord_dict.keys()==_processed_measurement_df_dict.keys():
                _save_name =_exp_name+'_'+f'midpoint_using_{ref_type}.xlsx'
            else:
                _save_name =_exp_name+'_'+f'midpoint_using_{ref_type}_self.xlsx'
            if type(_marker) is str and type(_cell_line) is str:
                _save_name = _cell_line + '_' + _marker + '_'+_save_name
                df.to_excel(analysis_save_folder+ os.sep+ _save_name, index=False)
            
    return sel_measurement_df_dict, sel_x_dict