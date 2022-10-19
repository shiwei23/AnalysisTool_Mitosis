# -*- coding: utf-8 -*-
"""
Created on Sat Jul 10 15:58:59 2021

@author: Shiwei Liu @ Harvard University

scripts for mitotic cell segmentations and analysis; generates segmentated images and measurmenet dataframes;
"""



from skimage import io
import numpy as np
from skimage.filters.thresholding import threshold_li, threshold_otsu, threshold_triangle
from skimage.morphology import erosion, dilation, opening, closing, white_tophat, remove_small_objects
from skimage.morphology import disk
from scipy.ndimage.morphology import binary_fill_holes
from skimage.measure import label, regionprops, regionprops_table
from skimage.filters import meijering, sato, frangi, hessian

import math
import skimage
from skimage import io

import skimage.segmentation as seg
import skimage.filters as filters
import skimage.draw as draw
import skimage.color as color


import sys,os, glob
import shutil
from skimage import data

import scipy
from scipy import ndimage as ndi

import matplotlib.pyplot as plt
import copy

import pandas as pd

def __init__():
    pass



def consolidate_exp_to_single_fd (data_saved_folder, target_subfolder = "cells", keep_source=True):
    
    target_folder = data_saved_folder + os.sep + target_subfolder
    if os.path.exists(target_folder):
        print ('Experiments are already consolidated. Check folder to see.')
        return None

    else:
        os.mkdir(target_folder)
        print ('Consolidating all experiment files into a single folder.')
        
    ori_files = [_f for _f in glob.glob(data_saved_folder+os.sep+'*') if _f.split('\\')[-1]!=target_subfolder]
    ori_files = [_f for _f in ori_files if 'analysis_param_info' not in _f]
    tar_files = [os.path.join(target_folder,ori_f.split('\\')[-1]) for ori_f in ori_files]
    import shutil  
    for _ori_f, _tar_f in zip(ori_files, tar_files):
        shutil.copy2(_ori_f, _tar_f)
    if not keep_source:
        for _f in ori_files:
            os.remove(_f)     
    return tar_files




def relocate_extra_cells (data_saved_folder, relocated_fd = 'not_used'):



    import shutil  

    subfolder_list = [file.path for file in os.scandir(data_saved_folder) if file.is_dir() ]
    image_folder_list = [file for file in subfolder_list if 'measurements' not in file and 'segmentations' not in file]

    relocate_cells = True
    if relocate_cells:
        # subfolders'name should end with a number between 0-9 (e.g., cell 1, cell 2) if data for each cell is saved in their unique subfolder.
        num_ends = ([str(num) for num in range (10)])
        image_folder_list = [file for file in image_folder_list if file.endswith(tuple(num_ends))]
        # sort the subfolder in numerical order
        image_folder_list = sorted(image_folder_list, key=lambda x: int(x.split('cell')[-1]))
        
        # Append images as sublist so image_list can be load from the main list by their index later.
        image_list_combined = []

        for _image_folder in image_folder_list:
            # Use nd file to filter the data of interest; each _image_folder should only analyze the nd with the most timepoints.
            nd_list = [file for file in glob.glob(_image_folder+os.sep+"*") if file[-2:]=="nd"]

            if len(nd_list) > 1:
                print ('Relocate cells with less timepoints to the "not_used" folder.')
                if not os.path.exists(os.path.join(_image_folder, relocated_fd)):
                    os.mkdir(os.path.join(_image_folder, relocated_fd))

                # get the the cell nd with most tifs     
                images_list_cand_lists = []
                images_list_cand_num = []

                for _nd in nd_list:
                    cell_prefix = _nd [:-3] + '_'
                    images_list_cand = [file for file in glob.glob(_image_folder+os.sep+"*") if cell_prefix in file]
                    images_list_cand.append(_nd)
                    images_list_cand_lists.append(images_list_cand)
                    images_list_cand_num.append(len(images_list_cand))
                
                cell_most_tif = np.argmax (np.array(images_list_cand_num))

                #_sel_nd = nd_list [cell_most_tif]
                #cell_prefix = _sel_nd [:-3] + '_'

                for _cand_index, images_list_cand in enumerate(images_list_cand_lists):
                    if _cand_index != cell_most_tif:
                        ori_files = images_list_cand
                        tar_files = [_f.replace(_image_folder, os.path.join(_image_folder, relocated_fd)) for _f in ori_files]
                        for _ori_f, _tar_f in zip(ori_files, tar_files):
                            shutil.copy2(_ori_f, _tar_f)
                            
                        for _ori_f in ori_files:
                            os.remove(_ori_f)     
    
    return None



def read_mitotic_cell_info (mm_data_list, filename):


    cell_info_dict = {}
    return cell_info_dict




def check_cell_info (cell_info_df):
    
    ''' A function to check the cell_info file.'''
    if 'cell_id' and 'cell_center_x' and 'cell_center_y' and 'cell_size' and 'z_selected' and 'dna_masks' in cell_info_df.columns:
        print ("Cell info has correct format.")
        return True
    else:
        print ("Incorrect format for cell info.")
        return False






def create_cell_info_template (data_saved_folder, image_folder_list):

    ''' A function to generate a cell_info template.'''

    new_df  = pd.DataFrame()
    new_df['cell_id'] = ["cell" + str(i+1) for i in range(len(image_folder_list))]
    new_df['cell_center_x'] =[np.nan] * len(image_folder_list)
    new_df['cell_center_y'] = [np.nan] * len(image_folder_list)
    new_df['cell_size'] =[np.nan] * len(image_folder_list)
    new_df['z_selected'] =[np.nan] * len(image_folder_list)
    new_df['dna_masks'] =[np.nan] * len(image_folder_list)
    new_df['not_to_analyze']=[np.nan]* len(image_folder_list)
  
    cell_info_name = data_saved_folder+os.sep+"cell_info.xlsx"
    if os.path.exists(new_df):
        print('"Cell info excel already exists. Manually delete the original if desiring to generate new file.')
        pass
    else:
        new_df.to_excel(cell_info_name, index=False)
        print ("An empty cell info excel created. Edit the cell info if necessary.") 
        print ("Re-run this box and make sure that the cell_info is valid before proceeding.")

    return None




def create_cell_info_from_raw (data_saved_folder, cell_info_raw_df, generate_info_for_mask=False):

    ''' A function to generate a cell_info from cell_info raw.'''

    all_cell_info  = pd.DataFrame()
    for _i, _cell_raw in enumerate(cell_info_raw_df.iloc()):
        new_cell_info  = pd.DataFrame()
        new_cell_info['cell_id'] = ["cell" + str(_i+1)]
        new_cell_info['cell_center_x']=[int(_cell_raw['X'])]
        new_cell_info['cell_center_y']=[int(_cell_raw['Y'])]
        new_cell_info['cell_size']=[int(np.sqrt(_cell_raw['Area']))]
        new_cell_info['z_selected']=[int(_cell_raw['Slice'])]
        if generate_info_for_mask:
            new_cell_info['dna_masks']=[f'dna_masks/cell{str(_i+1)}mask'] 
        else:
            new_cell_info['dna_masks']=[np.nan]
        new_cell_info['not_to_analyze']=[np.nan]
        all_cell_info=pd.concat([all_cell_info,new_cell_info])
  
    cell_info_name = data_saved_folder+os.sep+"cell_info.xlsx"
    if os.path.exists(cell_info_name):
        all_cell_info = pd.read_excel(cell_info_name, index_col=None)
        print('Use existing Cell info excel. Manually delete the original if desiring to generate new file.')
        pass
    else:
        all_cell_info.to_excel(cell_info_name, index=False)
        print ("Cell info excel created based on the cell_info_raw file. Edit this cell info if necessary.") 
    
    return all_cell_info



def img_xyz_reshaping (img_each_fov):

    ''' A function to re-order the z-axis to the first dim for xyz images'''
    # if is xyz 3d im
    if len(img_each_fov.shape)==3:
        if np.argmin(img_each_fov.shape)==2:
            img_each_fov = np.moveaxis (img_each_fov,-1,0)
            return img_each_fov
        if np.argmin(img_each_fov.shape)==0:
            return img_each_fov
    # if is xy 2d im        
    elif len(img_each_fov.shape)==2:
        print ('Image contains only a single slice. Keep the slice as it.')
        return img_each_fov
    else:
        print ('Image shape is not correct, exit.')
        return None




# Generate dna_masks folders; default is not run
def mkfd_dna_masks (data_saved_folder, cell_info, mk_mask_fd=False):
    
    if mk_mask_fd:
        mask_fd_path = data_saved_folder+os.sep+"dna_masks"
        if os.path.exists(mask_fd_path):
            print ('DNA mask folder exists.')
        else:
            os.mkdir(mask_fd_path)
            for _cell in cell_info.iloc():
                sub_fd = _cell['cell_id']+'mask'
                os.mkdir(mask_fd_path+os.sep+sub_fd)
                
    return None



def ch_index_after_ch_exclusion (initial_ch_index,num_of_ch, ch_to_be_analyzed = []):
    
    ''' A function to change the initial channle index after channel exclusion.
    
    return the new index for the initial channel selected after channel exlusion'''
    all_ch = list(range(num_of_ch))
    all_ch = [i+1 for i in all_ch]
    
    if ch_to_be_analyzed == []:
        ch_to_be_analyzed = all_ch
    
    excluded_ch = [_ch for _ch in all_ch if _ch not in ch_to_be_analyzed]
    
    num_exclusion = 0
    for _ex_ch in excluded_ch:
        if _ex_ch < initial_ch_index:
            num_exclusion +=1
    
    if initial_ch_index in ch_to_be_analyzed:
        new_ch_index = initial_ch_index - num_exclusion
        return new_ch_index
    else:
        print ("The given channel index is wrong. Use 1 as default.")
        return 1



def generate_mask_filename_list_for_timelapsed (num_of_fov, mask_file_folder = '', use_saved_mask = False):

    '''A function to read and check mask file location in the data save folder based on the subfolder(s) location such as the one provided by cell info for example.
    
    For timelasped data where the cell of interest has multiple timepoints/fovs.'''
    if use_saved_mask == True:
        _all_masks = [_mask for _mask in glob.glob(mask_file_folder+os.sep+"*") if _mask [-3:] == "roi" or _mask [-3:] == "tif"]
        # sort mask files numerically
        if len(_all_masks) > 1:
            _all_masks = sorted(_all_masks, key=lambda x: int(x.split('.')[0].split('\\')[-1]))

        if num_of_fov == len(_all_masks):
            mask_filename_list = _all_masks
        else:
            print ('Number of masks provided do not match number of fovs/timepoints. Use saved masks for timepoints that have match only.')
            #use_saved_mask =False deprecated
            # get all timepoints as string
            _exist_timepoint = [_mask.split("\\")[-1].split(".")[0] for _mask in _all_masks]
            mask_filename_list = []
            _default_mask = np.zeros([3,3])
            for _num_fov in range (num_of_fov):
                # add mask that have a saved roi or tiff file
                if str(_num_fov + 1) in _exist_timepoint:
                    _exist_mask = [_mask for _mask in _all_masks if _mask.split("\\")[-1].split(".")[0]== str(_num_fov + 1)][0]  # unravel the list to get the unique mask
                    mask_filename_list.append(_exist_mask)
                # if no saved, add the empty default mask
                else:
                    mask_filename_list.append(_default_mask)

    else:
        mask_filename_list = []
        for _num_fov in range (num_of_fov):
            _default_mask = np.zeros([3,3])
            mask_filename_list.append(_default_mask)
    
    return mask_filename_list




def generate_mask_filename_list_for_timelapsed (num_of_fov, mask_file_folder = '', use_saved_mask = False):

    '''A function to read and check mask file location in the data save folder based on the subfolder(s) location such as the one provided by cell info for example.
    
    For timelasped data where the cell of interest has multiple timepoints/fovs.'''
    if use_saved_mask == True:
        _all_masks = [_mask for _mask in glob.glob(mask_file_folder+os.sep+"*") if _mask [-3:] == "roi" or _mask [-3:] == "tif"]
        # sort mask files numerically
        if len(_all_masks) > 1:
            _all_masks = sorted(_all_masks, key=lambda x: int(x.split('.')[0].split('\\')[-1]))

        if num_of_fov == len(_all_masks):
            mask_filename_list = _all_masks
        else:
            print ('Number of masks provided do not match number of fovs/timepoints. Use saved masks for timepoints that have match only.')
            #use_saved_mask =False deprecated
            # get all timepoints as string
            _exist_timepoint = [_mask.split("\\")[-1].split(".")[0] for _mask in _all_masks]
            mask_filename_list = []
            _default_mask = np.zeros([3,3])
            for _num_fov in range (num_of_fov):
                # add mask that have a saved roi or tiff file
                if str(_num_fov + 1) in _exist_timepoint:
                    _exist_mask = [_mask for _mask in _all_masks if _mask.split("\\")[-1].split(".")[0]== str(_num_fov + 1)][0]  # unravel the list to get the unique mask
                    mask_filename_list.append(_exist_mask)
                # if no saved, add the empty default mask
                else:
                    mask_filename_list.append(_default_mask)

    else:
        mask_filename_list = []
        for _num_fov in range (num_of_fov):
            _default_mask = np.zeros([3,3])
            mask_filename_list.append(_default_mask)
    
    return mask_filename_list

     


def generate_masks_from_saved (mask_file,image_shape):

    '''A function to generate 0-1 binary mask from saved files including imageJroi and imageJ mask tiff.'''
        
    import roifile
    
    # for imageJ roi file.
    if mask_file [-3:] == 'roi':
        _roi = roifile.ImagejRoi.fromfile (mask_file)
        coord_roi= _roi.coordinates(multi=True)
        _mask = np.zeros (image_shape)
        for _coord_group in coord_roi:
            _coord_group [:,[1,0]] = _coord_group [:,[0,1]]
            _mask_group = skimage.draw.polygon2mask(image_shape, _coord_group)
            _mask = _mask + _mask_group
        return _mask
    
    # for imageJ mask tif.
    if mask_file [-3:] == 'tif':
        _mask_raw = io.imread(mask_file)
        # A valid mask should have two unqie values and the smalles value should be zero
        if len(np.unique(_mask_raw)) !=2 or np.min(_mask_raw) !=0:
            print ('Invalid mask file. Return the default empty mask.')
            # default array to be compatible with other segmentation functions
            return np.zeros ([3.3])
        # A valid mask should the same dim shape as the specified image shape.
        if _mask_raw.shape != image_shape:
            print ('Invalid mask shape. Return the default empty mask.')
            # default array to be compatible with other segmentation functions
            return np.zeros ([3.3])

        else:
            # Convert to 0-1 binary format if the provided mask is not.
            if np.max(_mask_raw) !=1:
                _mask_raw[_mask_raw==np.max(_mask_raw)] = 1
                
            return _mask_raw
        
    else:
        print ('No saved mask. Use the default empty mask.')
        # default array to be compatible with other segmentation functions
        return np.zeros ([3,3])







def generate_metamorph_data_list_each_cell (data_saved_folder, data_folder_num = 'SINGLE'):

    '''A functioin to read the parent folder containing sub-folder(s) where multi-channel TIF data acquired from metamorph (mm)
    
    Input: the path for the parent folder (data_saved_folder);
    Output: the image file name list'''
    
    subfolder_list = [file.path for file in os.scandir(data_saved_folder) if file.is_dir() ]
    image_folder_list = [file for file in subfolder_list if 'measurements' not in file and 'segmentations' not in file]
    
    
    if data_folder_num == 'MULTI': 

        # subfolders'name should end with a number between 0-9 (e.g., cell 1, cell 2) if data for each cell is saved in their unique subfolder.
        num_ends = ([str(num) for num in range (10)])
        image_folder_list = [file for file in image_folder_list if file.endswith(tuple(num_ends))]
        # sort the subfolder in numerical order
        image_folder_list = sorted(image_folder_list, key=lambda x: int(x.split('cell')[-1]))
        
        # Append images as sublist so image_list can be load from the main list by their index later.
        image_list_combined = []

        for _image_folder in image_folder_list:
            # Use nd file to filter the data of interest; each _image_folder should only analyze the nd with the most timepoints.
            nd_list = [file for file in glob.glob(_image_folder+os.sep+"*") if file[-2:]=="nd"]
            if len(nd_list) == 0:
                print ("No nd file exits. Skip this image folder")
                images_list = []

            if len(nd_list) > 0:
                images_list_cand_lists = []
                images_list_cand_num = []

                for _nd in nd_list:
                    cell_prefix = _nd [:-3] + '_'
                    images_list_cand = [file for file in glob.glob(_image_folder+os.sep+"*") if file[-3:]=="TIF" and "thumb" not in file and cell_prefix in file]
                    images_list_cand_lists.append(images_list_cand)
                    images_list_cand_num.append(len(images_list_cand))

                cell_most_tif = np.argmax (np.array(images_list_cand_num))
                images_list = images_list_cand_lists [cell_most_tif]

                _sel_nd = nd_list [cell_most_tif]
                cell_prefix = _sel_nd [:-3] + '_'
                
                # sort by timepoints in numerical order 1,2,..10,11,etc
                #images_list = sorted(images_list, key=lambda x: int(x.split('.TIF')[0].split('laser_t')[-1]))
                images_list = sorted(images_list, key=lambda x: int(x.split(cell_prefix)[1].split('.TIF')[0].split('_t')[1]))
                # sort by laser channels again (using the first number for wavelength)
                images_list = sorted(images_list , key=lambda x: int(x.split(cell_prefix)[1].split('w')[1][0]))

            image_list_combined.append(images_list)
        

    if data_folder_num == 'SINGLE':    

        image_folder_list = [file for file in image_folder_list if 'cells' in file and "cell_info_ims" not in file]

        # Append images as sublist so image_list can be load from the main list by their index later.
        image_list_combined = []

        if len(image_folder_list)!=1:
            print ("More than one image folder exists. Skip this data folder.")
            images_list = []
        
        if len(image_folder_list)==1:
            image_folder_list = image_folder_list[0]
            nd_list = [file for file in glob.glob(image_folder_list+os.sep+"*") if file[-2:]=="nd"]

            if len(nd_list) == 0:
                print ("No nd file exits. Skip this image folder")
                images_list = []

            if len(nd_list) > 0:
                # sort cell nd in numerical order
                nd_list = sorted(nd_list, key=lambda x: int(x[:-3].split('_')[-1]))
                for _nd in nd_list:
                    cell_prefix = _nd [:-3] + '_'
                    images_list = [file for file in glob.glob(image_folder_list+os.sep+"*") if file[-3:]=="TIF" and "thumb" not in file and cell_prefix in file]

                    # sort by timepoints in numerical order 1,2,..10,11,etc
                    #images_list = sorted(images_list, key=lambda x: int(x.split('.TIF')[0].split('laser_t')[-1]))
                    images_list = sorted(images_list, key=lambda x: int(x.split(cell_prefix)[1].split('.TIF')[0].split('_t')[1]))
                    # sort by laser channels again (using the first number for wavelength)
                    images_list = sorted(images_list , key=lambda x: int(x.split(cell_prefix)[1].split('w')[1][0]))

                    image_list_combined.append(images_list)
    

    return image_list_combined


def generate_metamorph_data_list_fixed_cell (data_saved_folder, data_folder_num = 'SINGLE', cell_ind_spliter='x '):

    subfolder_list = [file.path for file in os.scandir(data_saved_folder) if file.is_dir() ]
    image_folder_list = [file for file in subfolder_list if 'measurements' not in file and 'segmentations' not in file]

    if data_folder_num == 'SINGLE':    

        image_folder_list = [file for file in image_folder_list if file.split('\\')[-1]=='cells' and "cell_info_ims" not in file]

        # Append images as sublist so image_list can be load from the main list by their index later.
        image_list_combined = []

        if len(image_folder_list)!=1:
            print ("More than one image folder exists. Skip this data folder.")
            images_list = []
        
        if len(image_folder_list)==1:
            image_folder_list = image_folder_list[0]
            nd_list = [file for file in glob.glob(image_folder_list+os.sep+"*") if file[-2:]=="nd"]

            if len(nd_list) == 0:
                print ("No nd file exits. Skip this image folder")
                images_list = []

            if len(nd_list) > 0:
                # sort cell nd in numerical order  
                # cell_ind_spliter is the specific text before each nd
                nd_list = sorted(nd_list, key=lambda x: int(x[:-3].split(cell_ind_spliter)[-1]))
                for _nd in nd_list:
                    cell_prefix = _nd [:-3] + '_'
                    images_list = [file for file in glob.glob(image_folder_list+os.sep+"*") if file[-3:]=="TIF" and "thumb" not in file and cell_prefix in file]

                    # NO sort by timepoints since only one timepoint exist
                    #images_list = sorted(images_list, key=lambda x: int(x.split('.TIF')[0].split('laser_t')[-1]))
                    #images_list = sorted(images_list, key=lambda x: int(x.split(cell_prefix)[1].split('.TIF')[0].split('_t')[1]))
                    # sort by laser channels again (using the first number for wavelength)
                    images_list = sorted(images_list , key=lambda x: int(x.split(cell_prefix)[1].split('w')[1][0]))

                    image_list_combined.append(images_list)
    

    return image_list_combined



# A Function to grab ims for each cell for getting cell info raw based on the image_list_combined above; use in conjuncation with ImageJ code
def get_ims_for_cell_info(image_list_combined, ch_to_save ="488", time_mode ='last', save_ims=False, save_fd='cell_info_ims'):
        
    _ims_for_cell_info = []    
    print('Generating the folder for cell info ImageJ/FIJI use.')
    # for each timpoint
    for _cell_ims in image_list_combined: 
        # get ims matching the channel of interest
        _cell_ims_ch = [_ims for _ims in _cell_ims if ch_to_save in _ims]  
        # add ims accordingly
        if len(_cell_ims_ch) >0:
            if time_mode =="last":
                _cell_ims_ch_t = _cell_ims_ch[-1]
            elif time_mode =="first":
                _cell_ims_ch_t = _cell_ims_ch[0]
            elif time_mode =="middle":
                _t = int(len(_cell_ims_ch)/2)
                _cell_ims_ch_t = _cell_ims_ch[_t]
            else: # get random one
                _t = np.random.randint(0,len(_cell_ims_ch))
                _cell_ims_ch_t = _cell_ims_ch[_t]       
            _ims_for_cell_info.append(_cell_ims_ch_t)
    # if save ims to a separate folder
    if save_ims:
        import shutil   
        target_folder_main = "\\".join(_ims_for_cell_info[0].split("\\")[:-2])     
        target_folder = target_folder_main+os.sep+save_fd
        # setup and clean target folder
        if not os.path.isdir(target_folder):
            os.makedirs(target_folder)
        if len(os.listdir(target_folder))>0:
            _old_files = [_f for _f in glob.glob(target_folder+os.sep+'*')]
            for _f in _old_files:
                # keep cell info raw csv
                if "cell_info_raw" not in _f:
                    os.remove(_f)
       # copy ims for cell info usage; add pusedo cell index in case ims have same name (but from diff multi cell folder)    
        _target_files = []   
        for _ind, _f in enumerate(_ims_for_cell_info):
            if _ind <10:
                _target_files.append(os.path.join(target_folder, (f'cell0{_ind}'+_f.split("\\")[-1])))
            else:
                _target_files.append(os.path.join(target_folder, (f'cell{_ind}'+_f.split("\\")[-1])))
        for _ori_f, _tar_f in zip(_ims_for_cell_info, _target_files):
            shutil.copyfile(_ori_f, _tar_f)
        
    return _ims_for_cell_info



# A Function to duplicate images for a certain cell FOV where if more than one cell of interest exists
def duplicate_ims_for_cell_single_fd (image_list_combined, cell_index_to_duplcate, copy_nd=True, overwrite_ims=False):
    # index start from 0
    # only support folder with max cell index < 10 for now
    if cell_index_to_duplcate < len(image_list_combined):
        print('-- Duplicate the cell images at the end in the current folder')
        cell_index_len =len(str(cell_index_to_duplcate))
        #print(cell_index_len)
        
        _cell_ims = image_list_combined[cell_index_to_duplcate]
        _last_cell_ims = image_list_combined[-1]
        _last_index = _last_cell_ims[0].split('\\')[-1].split('_w')[0][-cell_index_len:]
        
        import shutil
        target_folder_main = "\\".join(_cell_ims[0].split("\\")[:-2])  
        target_folder = target_folder_main+os.sep+'cells'
        target_files = []
        
        _ori_index = _cell_ims[0].split('\\')[-1].split('_w')[0][-cell_index_len:] # only support cell < 10 for now
        _new_index = int(_last_index)+1
        _ori_nd_index = f'_{_ori_index}'
        _new_nd_index = f'_{_new_index}'
        _ori_im_index = f'_{_ori_index}_' 
        _new_im_index = f'_{_new_index}_'
        

        for _ori_im in _cell_ims:
            _new_name = _ori_im.split('\\')[-1].replace(_ori_im_index, _new_im_index)
            _tar_im = os.path.join(target_folder,_new_name)
            if os.path.exists(_tar_im) and not overwrite_ims:
                print ('-- File already exists, exit.')
            else:    
                shutil.copyfile(_ori_im, _tar_im)
                target_files.append(_tar_im)
        
        if copy_nd:
            nd_basename = _ori_im.split('\\')[-1].split('_w')[0]
            _ori_nd = os.path.join(target_folder,(nd_basename+ '.nd'))
            _tar_nd = os.path.join(target_folder,(nd_basename.replace(_ori_nd_index, _new_nd_index)+'.nd'))
            print(_new_nd_index,nd_basename.replace(_ori_nd_index, _new_nd_index))
            if os.path.exists(_tar_nd) and not overwrite_ims:
                print ('-- File already exists, exit.')

            else:
                shutil.copyfile(_ori_nd, _tar_nd)
                target_files.append(_tar_nd)

    return target_files





# A function to load the single-z slice of interest from a TIF dataset acquired by metamorph (mm)


def load_metamorph_tifs_single_z (mm_data_list, num_of_ch, ch_to_be_analyzed = [], ch_index = 1, z=0, 
                           timelapsed = False, mode='STD'):

    
    '''A function to read a folder directly containing multi-channel TIF (z-xy or xy-z) data acquired from metamorph (mm); group them by fov (or timepoints);
       next, the z-slice of interest is loaded for each fov based on the mode selected below. 

    Input:
      mm_data_list: list (dataset within a folder) containing the TIFs of interest;
      num_of_ch: int; the number of channels;
      ch_dict: e.g., {"dna":1,"53BP1SC":2,"mdc1":3,"pol2S5":4} where 1-2-3-4 corresponds to ch405-ch488-ch560-ch647 (or Wide_field_488-560, etc);
      ch_index: int; the channel index from the ch_dict to be used for finding the z-slice of interest;
      z: z for the specified z; default is 0;
      timelapsed: set True if files in the folder is from timelaspsed image; 
      this is because the naming order of the files is different for timelapsed images versus fixed multi-position images.
      mode = 'STD': the z-slice has the largest STD for the channel of interest; 
      mode = 'MITO': the z-slice has the longest chromosomes major axis in mitosis; NOT COMPLETE YET;
      mode ='USER': the z-slice specified by the user; applied to all images.

       Output: A list containing all z-sliced io_read zxy ims (ims_array).'''
   

    if ch_to_be_analyzed == []:
        ch_to_be_analyzed = list(range(num_of_ch))
        ch_to_be_analyzed = [i+1 for i in ch_to_be_analyzed]


    # An empty list to store all the ims.
    sorted_mm_data_list = []
    
    # Check if all fovs have the same number of channels as specified.
    if len(mm_data_list)%num_of_ch != 0:
        print ("Empty folder or the number of channel does not match the data.")


    else:
        for _fov_id in range(int(len(mm_data_list)/num_of_ch)):
            
            # Fixed multipositioned images are ordered by the nd extension file, then channels for one timepoint.
            if timelapsed == False:
                data_each_fov = mm_data_list [(_fov_id*num_of_ch):(_fov_id*num_of_ch + num_of_ch)] 
                # TIFF stacks from different channels for each fov
            
            # Timelapased images, which are ordered by each channel for all timepoints; then the next channel.
            if timelapsed == True:
                num_of_time = int(len(mm_data_list)/num_of_ch)
                stop_pos = _fov_id + num_of_time*(num_of_ch-1) +1                  
                data_each_fov = mm_data_list [(_fov_id):(stop_pos):(num_of_time)]
                
            # An empty list to store the ims for each fov (or timepoint)   
            img_each_fov = []
            
            # Print(data_each_fov)  # un- comment this for test/verbal purposes.
            
            # Use the best focal plane (which typically has the largest std) of the specified channel to find the z-slice of interest.
            if mode == 'STD':
                img_indexed = io.imread (data_each_fov[ch_index-1])
                img_indexed = img_xyz_reshaping (img_indexed)
           
                img_fl = []
                img_std = []
                for _slice in range (len(img_indexed)):
                    _fl = np.array(img_indexed[_slice,:,:].flatten())
                    img_fl.append(_fl)
                    img_std.append(np.std(_fl))
                best_z_index = np.argmax (np.array(img_std))
                print(f"Loading the slice {best_z_index + 1} for the image {_fov_id} in this dataset")
                
            # Use the chromosome major axis length to find the z-slice of interest  
            if mode == 'MITO':  # NOT COMPLETED YET!!!
                #if 'dna' or 'DNA' or 'DAPI' or 'dapi' or 'H2B' or 'h2b' in ch_dict.keys():
                    #ch_index = ch_dict['dna' or 'DNA' or 'DAPI' or 'dapi' or 'H2B' or 'h2b']
                # This seems not easy since DNA segmentation might not be straightforward enough in many complicated scenarios.
                img_indexed = io.imread (data_each_fov[ch_index-1])  # index-1 to match the python index
                img_indexed = img_xyz_reshaping (img_indexed)

                img_dna_len = []
                for _slice in range(len(img_indexed)):
                    otsu_value_dna = threshold_otsu (img_indexed[_slice,:,:])
                    dna_mask = img_indexed[_slice,:,:] > otsu_value_dna
                    dna_mask = remove_small_objects(dna_mask, 5000,connectivity=1)
                    dna_mask = binary_fill_holes(dna_mask)
                    dna_mask = dilation(dna_mask, disk(5))
                    dna_mask = erosion(dna_mask, disk(5))
                    regions = regionprops(label(dna_mask))
                    major_dna_axis_length = []
                    for props in regions:
                        y0, x0 = props.centroid
                        major_dna_axis_length.append(props.major_axis_length)
                    
                    print(np.mean(np.array(major_dna_axis_length)))
                    img_dna_len.append(np.mean(np.array(major_dna_axis_length)))
                best_z_index = np.argmax (np.array(img_dna_len))
                print(f"Loading the slice {best_z_index + 1} for the image {_fov_id} in this dataset")
                
            # Use the USER's defined z-slice, which would be applied to all images for this fov (or timepoint). 
            if mode == 'USER':
                img_indexed = io.imread (data_each_fov[ch_index-1])
                img_indexed = img_xyz_reshaping (img_indexed)
                best_z_index = z -1
                print(f"Loading the slice {best_z_index + 1} for the image {_fov_id} in this dataset")
                
            # Store all the ims for the specified z-slice for this fov (or timepoint).
            for _ch in ch_to_be_analyzed:                
                ch_each_fov = img_xyz_reshaping(io.imread (data_each_fov[_ch-1]))
                ch_each_fov_bf = ch_each_fov[best_z_index,:,:]
                img_each_fov.append(ch_each_fov_bf)
                
            # Store all the ims for the the specified z-slice for all fov (or timepoints).
            sorted_mm_data_list.append(img_each_fov)
            
        num_of_fov= len(sorted_mm_data_list)
        print (f"There are {num_of_fov} fovs/timepoints for this {num_of_ch}-channel, {len(img_indexed)}-slice dataset.")
        
        sorted_mm_data_array = np.array(sorted_mm_data_list)

    return sorted_mm_data_array




def load_metamorph_tifs_multi_z (mm_data_list, num_of_ch, ch_to_be_analyzed = [], timelapsed = False, exclude_transmitted=False):

    
    '''A function to read a folder directly containing multi-channel TIF (z-xy or xy-z) data acquired from metamorph (mm); group them by fov (or timepoints);
       next, all z-slice are loaded for each fov. 

    Input:
      mm_data_list: list (dataset within a folder) containing the TIFs of interest;
      num_of_ch: int; the number of channels;
      ch_dict: e.g., {"dna":1,"53BP1SC":2,"mdc1":3,"pol2S5":4} where 1-2-3-4 corresponds to ch405-ch488-ch560-ch647 (or Wide_field_488-560, etc);
      timelapsed: set True if files in the folder is from timelaspsed image; 
      this is because the naming order of the files is different for timelapsed images versus fixed multi-position images.
      

       Output: A list containing all z-sliced io_read zxy ims (ims_array).'''
   

    if ch_to_be_analyzed == []:
        ch_to_be_analyzed = list(range(num_of_ch))
        ch_to_be_analyzed = [i+1 for i in ch_to_be_analyzed]


    # An empty list to store all the ims.
    sorted_mm_data_list = []
    
    ch_check_flag = False

    # Check if all fovs have the same number of channels as specified.
    if len(mm_data_list)%num_of_ch != 0 or exclude_transmitted:
        print ("Empty folder or the number of channel does not match the data.")

        # temp fix for old 60x image where the Transmitted timepoint was taken differently:
        print ('Try excluding transmitted timepoints using the filename.') 
        filtered_mm_data_list = [_f for _f in mm_data_list if "Transmitted" not in _f]
        # re-check if ch match after excluding Transmitted channels
        if len(filtered_mm_data_list)%(num_of_ch-1) != 0:
            print ('Empty folder or the number of channel still does not match the data.')
        else:
            ch_check_flag = True
            # re-define the data list
            mm_data_list = filtered_mm_data_list
            # re-define the channel to be analyzed since one channel is pre-excluded before the later step
            new_ch_to_be_analyzed = []
            for _ch in ch_to_be_analyzed:
                new_ch = ch_index_after_ch_exclusion(_ch, num_of_ch,ch_to_be_analyzed)
                new_ch_to_be_analyzed.append(new_ch)
            ch_to_be_analyzed=new_ch_to_be_analyzed
            # re-define the number of channel
            num_of_ch=num_of_ch-1

    else:
        ch_check_flag = True

    if ch_check_flag:
        for _fov_id in range(int(len(mm_data_list)/num_of_ch)):
            
            # Fixed multipositioned images are ordered by the nd extension file, then channels for one timepoint.
            if timelapsed == False:
                data_each_fov = mm_data_list [(_fov_id*num_of_ch):(_fov_id*num_of_ch + num_of_ch)] 
                # TIFF stacks from different channels for each fov
            
            # Timelapased images, which are ordered by each channel for all timepoints; then the next channel.
            if timelapsed == True:
                num_of_time = int(len(mm_data_list)/num_of_ch)
                stop_pos = _fov_id + num_of_time*(num_of_ch-1) +1                  
                data_each_fov = mm_data_list [(_fov_id):(stop_pos):(num_of_time)]
                
            # An empty list to store the ims for each fov (or timepoint)   
            img_each_fov = []
            
            # Store all the ims for the all z-slice for this fov (or timepoint).
            for _ch in ch_to_be_analyzed:                
                ch_each_fov = img_xyz_reshaping(io.imread (data_each_fov[_ch-1]))
                img_each_fov.append(ch_each_fov)  # zyx data

            # duplicate ch if one ch has only one slice (to avoid bleaching, etc) while other have multiple slices (assuming the same) 
            num_slice_list = []
            for ch_each_fov in img_each_fov:
                num_slice_list.append(len(ch_each_fov.shape))
            num_slice_list=np.array(num_slice_list)
            ch_3d_index = np.where(num_slice_list==3)[0][0] # use one 3D ch to get the slice number assuming all have the same num of slices
            max_slice = img_each_fov[ch_3d_index].shape[0]
            ch_2d_index = np.where(num_slice_list==2)[0]
            if len(ch_2d_index)>0:
                for each_ch_2d_index in ch_2d_index:
                    old_ch_each_fov = img_each_fov[each_ch_2d_index]
                    print('A channel has only one slice, duplicate the slice to match the Z dimension.')
                    img_each_fov[each_ch_2d_index]=np.repeat(old_ch_each_fov[np.newaxis,:, :], max_slice, axis=0)

            # Store all the ims for the the specified z-slice for all fov (or timepoints).
            sorted_mm_data_list.append(img_each_fov)   # ch-zyx data

        num_of_fov= len(sorted_mm_data_list)
        print (f"There are {num_of_fov} fovs/timepoints for this {num_of_ch}-channel, {len(ch_each_fov)}-slice dataset.")

        sorted_mm_data_array_3D = np.array(sorted_mm_data_list)
        
    return sorted_mm_data_array_3D




def chromosome_mask_segmentation (img_each_fov, dna_ch_index, cell_mask_ch_index =1, core_thickness =10, dna_erosion_factor =3, 
                                  small_object_th = 3000, cell_center = [466,621], cell_size=500,
                             threshold_method = skimage.filters.thresholding.threshold_otsu,
                                  threshold_by_region = "GLOBAL", check_chrom_num = True,
                             inner_core_mode = False, midbox_mask = np.zeros([3,3]), chr_mask_from_saved = np.zeros([3,3])):
    
    
    ''' A function to read a ch - xy io_read ims arrays from the same fov (or timepoint); after ims loading, the function perform chromosome segmentation based on the selected mode.

      Inputs:
              img_each_fov: a ch - xy io_read ims arrays;
              dna_ch_index: int; the ch index (starting from 1) for the DNA channel;
              cell_mask_ch_index: int; the ch index for the cell mask channel; used when thresholding is done by 'CELL' mode;
              core_thichness: int; the width for the DNA rim used for inner core mask generation;
              dna_erosion_factor: disk size for erosion, dilation, etc
              threshold_method: default is as shown above (otsu);
              cell_center: [x,y] coord of the center of the cropped region;
              cell_size: int; length of the cropped square;
              threshold_by_region =  "GlOBAL": use the whole image for segmentation th;
              threshold_by_region =  "CELL": use the whole image for segmentation th, together with using a cell mask to facilitate the segmentation;
              threshold_by_region =  "USER": use the cropped region for segmentation th; the image will be cropped during analysis; for many scenarios where only one mitotic cell is desired when other cells are present, this mode seems most effective for segmenation.
              inner_core_mode: whether to generate the inner core mask or not;
              midbox_mask: provided as a valid mask for inner core mask generation; otherwise, the function will generate one.
              chr_mask_from_saved: saved roi mask (e.g., generated by FIJI imageJ) that has the same XY shape as the raw image; use this if (semi)manually define chr ROI.

       Outputs:
             the DNA/chromosome mask (when inner_core_mode == False);
             [the DNA/chromosome mask, the inner_core_mask] (when inner_core_mode == True).
             
             Note that when USER mode is used, the output image/mask size will be cropped. Thus, be consistent when using USER mode for all relevant functions (including many below).'''

    # Initiate the status of thresholding method, which can be modified later to report the method in the results section.
    use_li_method = False
    if threshold_method == skimage.filters.thresholding.threshold_li:
        use_li_method = True
    else:
        pass
        
    dna_image = img_each_fov[dna_ch_index-1]
    
    # Remove hot pixel for the dna channel; use 0.02 higer percentile as the hot pixel th; intensity higher than hot pixel th will be converted to the th -- smoothing the hot pixel of the image.
    hot_pixel_th = np.percentile(dna_image,99.8)
    # Deepcopy to prevent direct altering of the loaded images:
    import copy
    dna_image_copy = copy.deepcopy(dna_image)
    dna_image_copy[dna_image>hot_pixel_th] = np.median(dna_image)
    dna_image = dna_image_copy
    
    
    # Use the whole image to perform segmentation.
    if threshold_by_region == "GLOBAL":
        # perform global segmentation:
        segmentation_th = threshold_method (dna_image)
        dna_mask = dna_image > segmentation_th
        dna_mask = remove_small_objects(dna_mask, small_object_th,connectivity=1)

            
    # Use the cell mask generated from a different specified channel.  er 
    if threshold_by_region == "CELL":
        # perform cell_mask_based segmentation using the other ch, as specifided by cell_mask_ch_index
        # use this method when DNA image is not definitive but other ch is
        cell_mask_image = img_each_fov[cell_mask_ch_index-1]
        segmentation_th = threshold_method (cell_mask_image)
        cell_mask = cell_mask_image > segmentation_th
        cell_mask = remove_small_objects(cell_mask, small_object_th,connectivity=1)
        dna_image_filtered = np.array([i for i in (dna_image * cell_mask).flatten() if i>0])
        segmentation_th = threshold_method (dna_image_filtered)
        dna_mask = dna_image_filtered >  segmentation_th
        
    # Use the cropped region of the image as specified by the user.
    if threshold_by_region == "USER":
        # perform segmentation using a cropped image, where the crop box is specified by the user
        # use this method when cell of interest cannot be determined by features in image.
        if len(cell_center) !=2:
                print ('Specify the center to generate a crop.')
                return None
        
        dna_image = np.array(dna_image)[int(cell_center[1]-cell_size/2):int(cell_center[1]+cell_size/2),  
              int(cell_center[0]-cell_size/2):int(cell_center[0]+cell_size/2)]
        
        # Below may be modified later to enable a mask returned with the same shape as the non-cropped image.
        #pseudo_mask = np.zeros(dna_image.shape)
        #pseudo_mask [int(cell_center[1]-cell_size/2):int(cell_center[1]+cell_size/2),int(cell_center[0]-cell_size/2):int(cell_center[0]+cell_size/2)] =1
        #dna_image = dna_image * pseudo_mask
        
        segmentation_th = threshold_method (dna_image)
        dna_mask = dna_image > segmentation_th
        dna_mask = remove_small_objects(dna_mask, small_object_th,connectivity=1)
        
        
    # Addtional binary operations; small_object_th appears to be 3000 for the tested image set.
    # Numbers should be determined empierically by experiments.
    dna_mask = binary_fill_holes(dna_mask)
    dna_mask = dilation(dna_mask, disk(dna_erosion_factor)) 
    dna_mask = erosion(dna_mask, disk(dna_erosion_factor*2)) 
    dna_mask = remove_small_objects(dna_mask, small_object_th)

    # if do not check chromosome objects number, simply return the rough mask
    if not check_chrom_num:
        return dna_mask
    

    # Use pre-saved chr mask for chromosome and inner core segmentation (this can be modified to the top condition loop to be more efficient.)
    if chr_mask_from_saved.shape == (img_each_fov[dna_ch_index-1]).shape:

        if threshold_by_region == "GLOBAL" or "CELL":
            dna_mask = chr_mask_from_saved
        if threshold_by_region == "USER":
            dna_mask = chr_mask_from_saved[int(cell_center[1]-cell_size/2):int(cell_center[1]+cell_size/2),  
              int(cell_center[0]-cell_size/2):int(cell_center[0]+cell_size/2)]

    
    labeled_chrom, num_of_chrom  = ndi.label(dna_mask)
    
    # Remove edge objects
    kept_chrom = []
    for i in range(num_of_chrom):
        cand_chrom = labeled_chrom == i+1
        x_bound = cand_chrom.shape[0]
        y_bound = cand_chrom.shape[1]
        if cand_chrom[0,0] == 0 and cand_chrom[0,y_bound-1] == 0 and cand_chrom[x_bound-1,0] == 0 and cand_chrom[x_bound-1,y_bound-1] == 0:
            cand_chrom[cand_chrom>0]=1
            cand_chrom = np.array(cand_chrom)
            kept_chrom.append(cand_chrom)
            
    # Use li method if the default thresholding is too strigent and produces less than 2 chromosome objects.
    if len(kept_chrom) <2:
        use_li_method = True     
        segmentation_th = threshold_li (dna_image)
        dna_mask = dna_image > segmentation_th
        dna_mask = remove_small_objects(dna_mask, small_object_th,connectivity=1)
        dna_mask = binary_fill_holes(dna_mask)
        dna_mask = dilation(dna_mask, disk(dna_erosion_factor)) 
        dna_mask = erosion(dna_mask, disk(dna_erosion_factor*2)) 
        dna_mask = remove_small_objects(dna_mask, small_object_th,connectivity=1)
        dna_mask = erosion(dna_mask, disk(5)) 
        
        labeled_chrom, num_of_chrom  = ndi.label(dna_mask)
        # Remove edge objects
        kept_chrom = []
        for i in range(num_of_chrom):
            cand_chrom = labeled_chrom == i+1
            x_bound = cand_chrom.shape[0]
            y_bound = cand_chrom.shape[1]
            if cand_chrom[0,0] == 0 and cand_chrom[0,y_bound-1] == 0 and cand_chrom[x_bound-1,0] == 0 and cand_chrom[x_bound-1,y_bound-1] == 0:
                cand_chrom[cand_chrom>0]=1
                cand_chrom = np.array(cand_chrom)
                kept_chrom.append(cand_chrom)
                
    if len(kept_chrom) !=2:
        print (f'Fail to detect 2 chromosome objects.There are/is {len(kept_chrom)} object(s) detected')
        return dna_mask  # or some other empty thing for downstream analysis to read
    # Consolidate the 2 valid individual mask labels.  
    dna_mask = np.ma.mask_or (kept_chrom[0],kept_chrom[1])
    
    # If inner_core_mode is selected, generate inner core mask or not accordingly as below.
    if inner_core_mode == False:
        if use_li_method == False:
            print("Segmentation of chromosome using otsu (or other defualt) method.")
        if use_li_method == True:
            print("Segmentation of chromosome using li method.")
        return dna_mask
    
    if inner_core_mode == True:
        dna_rim_mask = dilation (dna_mask, disk(core_thickness)) *  (dna_mask ==0)
        
        # Use the midbox_mask provided by the user if it is valid.
        if midbox_mask.shape == dna_mask.shape:
            inner_core_mask = dna_rim_mask * midbox_mask
            inner_core_mask = dilation(inner_core_mask, disk(3))
            inner_core_mask = erosion(inner_core_mask, disk(3))
            inner_core_mask = remove_small_objects(inner_core_mask, 500)
        # Generate the midbox_mask based on the chromosome mask.
        # Essentially, the midbox_mask is a box covering the midzone, which also partially covers half of each chromosome sets facing the midzone.
        else:
            regions = regionprops(label(dna_mask))
            if len(regions)!=2:
                print (f'Fail to detect 2 chromosome objects. There are/is {len(regions)} object(s) detected')
                return dna_mask
            midbox = []
            for props in regions:
                y0, x0 = props.centroid
                orientation = props.orientation
                x1 = x0 - math.sin(orientation) * 0.5 * props.major_axis_length
                y1 = y0 - math.cos(orientation) * 0.5 * props.major_axis_length
                x2 = x0 + math.sin(orientation) * 0.5 * props.major_axis_length
                y2 = y0 + math.cos(orientation) * 0.5 * props.major_axis_length
                midbox.append ([x1,y1])
                midbox.append ([x2,y2])
            
            # Invert the order of xy because xy order is read differently in different array modules below.
            import copy   
            sorted_midbox =copy.deepcopy(midbox)
            sorted_midbox[-1] = midbox[-2]
            sorted_midbox[-2] = midbox[-1]

            # Add the first coord to the end to form a closed box.
            sorted_midbox.append(midbox[0])
            
            # faster alternative way below:
            # np.array(sorted_midbox)[:,[1,0]] = np.array(sorted_midbox)[:,[0,1]]
            xs, ys = zip(*sorted_midbox)
            midbox_inv = np.zeros([5,2])
            midbox_inv [:,0] = np.array(sorted_midbox)[:,1]
            midbox_inv [:,1] = np.array(sorted_midbox)[:,0]



            midbox_mask=skimage.draw.polygon2mask(dna_mask.shape, midbox_inv)
            inner_core_mask = dna_rim_mask * midbox_mask
            inner_core_mask = dilation(inner_core_mask, disk(3))
            inner_core_mask = erosion(inner_core_mask, disk(3))
            inner_core_mask = remove_small_objects(inner_core_mask, 500)

        # Report the thresholding method used.     
        if use_li_method == False:
            print("Segmentation of chromosome and inner core using otsu (or other defualt) method.")
        if use_li_method == True:
            print("Segmentation of chromosome and inner core using li method.")
        # Note that there are two elements returned if use inner_core_mode. 
        return dna_mask, inner_core_mask







def inverse_chromosome_mask_segmentation_by_line_profile (img_each_fov, dna_ch_index, core_thickness =10, 
                                                          dna_erosion_factor =5, extra_erosion_factor = 5,
                                                          cell_erosion_factor =20,
                                  small_object_th = 3000, cell_center = [466,621], cell_size=600,
                             threshold_method = skimage.filters.thresholding.threshold_otsu,
                                  threshold_by_region = "GLOBAL", non_chr_low_ratio =0.09, chr_peak_prominence = 10,
                             inner_core_mode = False, midbox_mask = np.zeros([3,3]),
                                         chr_mask_from_saved = np.zeros([3,3]), _verbose = True):

    """ inverse chromsome segmeation using exclusion of signal from chromosomes
    
         supports a ch - xy image"""
    import skimage
    import scipy



    # Use saved mask
    if chr_mask_from_saved.shape == (img_each_fov[dna_ch_index-1]).shape:
        if threshold_by_region == "GLOBAL" or "CELL":
            dna_mask = chr_mask_from_saved
        if threshold_by_region == "USER":
            dna_mask = chr_mask_from_saved[int(cell_center[1]-cell_size/2):int(cell_center[1]+cell_size/2),  
              int(cell_center[0]-cell_size/2):int(cell_center[0]+cell_size/2)]

    # Generate mask
    else:          
        # Initiate the status of thresholding method, which can be modified later to report the method in the results section.
        use_li_method = False
        if threshold_method == skimage.filters.thresholding.threshold_li:
            use_li_method = True
        else:
            pass
    
        dna_image = img_each_fov[dna_ch_index-1]
        # Deepcopy to prevent direct altering of the loaded images:
        import copy
        # Remove hot pixel for the dna channel; use 0.02 higer percentile as the hot pixel th; intensity higher than hot pixel th will be converted to the th -- smoothing the hot pixel of the image.
        hot_pixel_th = np.percentile(dna_image,99.8)
        dna_image_copy = copy.deepcopy(dna_image)
        dna_image_copy[dna_image>hot_pixel_th] = np.median(dna_image)
        dna_image = dna_image_copy
 
    
        # Use the whole image to perform segmentation.
        if threshold_by_region == "GLOBAL":
            # perform global segmentation:
            segmentation_th = threshold_method (dna_image)
            cell_mask = dna_image > segmentation_th
            cell_mask = remove_small_objects(cell_mask, small_object_th,connectivity=1)
    
        # Use the cropped region of the image as specified by the user.
        if threshold_by_region == "USER":
            # perform segmentation using a cropped image, where the crop box is specified by the user
            # use this method when cell of interest cannot be determined by features in image.
            if len(cell_center) !=2:
                print ('Specify the center to generate a crop.')
                return None
        
            dna_image = np.array(dna_image)[int(cell_center[1]-cell_size/2):int(cell_center[1]+cell_size/2),  
                int(cell_center[0]-cell_size/2):int(cell_center[0]+cell_size/2)]
        
            segmentation_th = threshold_method (dna_image)
            cell_mask = dna_image > segmentation_th
            cell_mask = remove_small_objects(cell_mask, small_object_th,connectivity=1)
        
        
        # Addtional binary operations; small_object_th appears to be 3000 for the tested image set.
        # Numbers should be determined empierically by experiments.
        cell_mask = binary_fill_holes(cell_mask)
        cell_mask = dilation(cell_mask, disk(3)) 
        cell_mask = erosion(cell_mask, disk(3*2)) 
        cell_mask = remove_small_objects(cell_mask, 10000)
    
        # generate the line profile across the middle of the cell
        cell_mask = erosion (cell_mask, disk(cell_erosion_factor))
        regions = regionprops(label(cell_mask))
        if len(regions) ==1:
            for props in regions:
                y0, x0 = props.centroid
                orientation = props.orientation
                x1 = x0 - math.sin(orientation) * 0.5 * props.major_axis_length
                y1 = y0 - math.cos(orientation) * 0.5 * props.major_axis_length
                x2 = x0 + math.sin(orientation) * 0.5 * props.major_axis_length
                y2 = y0 + math.cos(orientation) * 0.5 * props.major_axis_length
        
        
            from skimage.measure import profile_line
            # invert xy for reading image intensity 
            intensity_profile = profile_line(dna_image, [y1,x1],[y2,x2], linewidth=50, mode='reflect',reduce_func=np.median)
            # find two valleys for the chr where signal is low due to chromosome exlcusion;
            # valley distance and width needs to be adjusted for different experiments
            _chr_valley_info = scipy.signal.find_peaks(np.array(intensity_profile)*(-1), 
                                                   distance=50, prominence=chr_peak_prominence,width=15)
            # two or more valleys
            if len(_chr_valley_info[0])>=2:   
                # two most prominent valleys
                sorted_sel_ind = np.sort(np.argsort(_chr_valley_info[1]['prominences'])[-2:])
                first_ind = sorted_sel_ind[0]
                second_ind = sorted_sel_ind[1]
                # extract info from the two most prominent valleys
                _chr_valley_bases = [_chr_valley_info[1]['left_ips'][first_ind], _chr_valley_info[1]['left_ips'][second_ind],
                                                _chr_valley_info[1]['right_ips'][first_ind],_chr_valley_info[1]['right_ips'][second_ind]]
                # start with the highest intensity th from the vally bases
                _initial_th = max ([intensity_profile[int(round(_base))] for _base in _chr_valley_bases])
                _iter = True
                # iteratively decrease the intensity th
                while _iter:
                    _seg_ind = np.argwhere(np.array(intensity_profile)<_initial_th)
                    _bad_ind = 0
                    for _ind in _seg_ind:
                        if _ind > _chr_valley_info[1]['right_ips'][second_ind] or _ind < _chr_valley_info[1]['left_ips'][first_ind]:
                            _bad_ind+=1
                        elif _ind > _chr_valley_info[1]['left_ips'][second_ind] and _ind < _chr_valley_info[1]['right_ips'][first_ind]:
                            _bad_ind+=1
                    # iterate until less than 5% of ind outside the valleys are higher than the th
                    if _bad_ind/len(intensity_profile) >non_chr_low_ratio:
                        _initial_th-= 1
                    else:
                        _iter = False
        
            else:
                if _verbose:
                    print ('-- less than two chromosome vallies were found, return original cell mask')
                return cell_mask  
        else:
            if _verbose:
                print ('-- cell mask segmentation failed, return original cell mask')
            return cell_mask
        
    
        print (f'-- using threshold of {_initial_th} from line profile for inverse chromosome segmentation')
        chr_mask = (dna_image<_initial_th)*erosion(cell_mask,disk(cell_erosion_factor))
        chr_mask = remove_small_objects(chr_mask,small_object_th)
        chr_mask = dilation (chr_mask,disk(dna_erosion_factor))
        chr_mask = erosion (chr_mask,disk(dna_erosion_factor))
        chr_mask = erosion (chr_mask,disk(extra_erosion_factor))
        chr_mask = dilation (chr_mask,disk(1))
        dna_mask = remove_small_objects(chr_mask,small_object_th)
        

    labeled_chrom, num_of_chrom  = ndi.label(dna_mask)
    
    # Remove edge objects
    kept_chrom = []
    for i in range(num_of_chrom):
        cand_chrom = labeled_chrom == i+1
        x_bound = cand_chrom.shape[0]
        y_bound = cand_chrom.shape[1]
        if cand_chrom[0,0] == 0 and cand_chrom[0,y_bound-1] == 0 and cand_chrom[x_bound-1,0] == 0 and cand_chrom[x_bound-1,y_bound-1] == 0:
            cand_chrom[cand_chrom>0]=1
            cand_chrom = np.array(cand_chrom)
            kept_chrom.append(cand_chrom)
            

                
    if len(kept_chrom) !=2:
        print (f'Fail to detect 2 chromosome objects.There are/is {len(kept_chrom)} object(s) detected')
        return dna_mask  # or some other empty thing for downstream analysis to read


    # Consolidate the 2 valid individual mask labels.  
    dna_mask = np.ma.mask_or (kept_chrom[0],kept_chrom[1])
    
    # If inner_core_mode is selected, generate inner core mask or not accordingly as below.
    if inner_core_mode == False:
        return dna_mask
    
    if inner_core_mode == True:
        dna_rim_mask = dilation (dna_mask, disk(core_thickness)) *  (dna_mask ==0)
        
        # Use the midbox_mask provided by the user if it is valid.
        if midbox_mask.shape == dna_mask.shape:
            inner_core_mask = dna_rim_mask * midbox_mask
            inner_core_mask = dilation(inner_core_mask, disk(3))
            inner_core_mask = erosion(inner_core_mask, disk(3))
            inner_core_mask = remove_small_objects(inner_core_mask, 500)
        # Generate the midbox_mask based on the chromosome mask.
        # Essentially, the midbox_mask is a box covering the midzone, which also partially covers half of each chromosome sets facing the midzone.
        else:
            regions = regionprops(label(dna_mask))
            if len(regions)!=2:
                print (f'Fail to detect 2 chromosome objects. There are/is {len(regions)} object(s) detected')
                return dna_mask
            midbox = []
            for props in regions:
                y0, x0 = props.centroid
                orientation = props.orientation
                x1 = x0 - math.sin(orientation) * 0.5 * props.major_axis_length
                y1 = y0 - math.cos(orientation) * 0.5 * props.major_axis_length
                x2 = x0 + math.sin(orientation) * 0.5 * props.major_axis_length
                y2 = y0 + math.cos(orientation) * 0.5 * props.major_axis_length
                midbox.append ([x1,y1])
                midbox.append ([x2,y2])
            
            # Invert the order of xy because xy order is read differently in different array modules below.
            import copy   
            sorted_midbox =copy.deepcopy(midbox)
            sorted_midbox[-1] = midbox[-2]
            sorted_midbox[-2] = midbox[-1]

            # Add the first coord to the end to form a closed box.
            sorted_midbox.append(midbox[0])
            
            # faster alternative way below:
            # sorted_midbox = np.array(sorted_midbox)
            # (sorted_midbox)[:,[1,0]] = (sorted_midbox)[:,[0,1]]
            xs, ys = zip(*sorted_midbox)
            midbox_inv = np.zeros([5,2])
            midbox_inv [:,0] = np.array(sorted_midbox)[:,1]
            midbox_inv [:,1] = np.array(sorted_midbox)[:,0]



            midbox_mask=skimage.draw.polygon2mask(dna_mask.shape, midbox_inv)
            inner_core_mask = dna_rim_mask * midbox_mask
            inner_core_mask = dilation(inner_core_mask, disk(3))
            inner_core_mask = erosion(inner_core_mask, disk(3))
            inner_core_mask = remove_small_objects(inner_core_mask, 500)

        # Report the thresholding method used.     
        if use_li_method == False:
            print("Inverse segmentation of chromosome and inner core using otsu (or other defualt) method and central line profile.")
        if use_li_method == True:
            print("Inverse segmentation of chromosome and inner core using li method and central line profile.")
        # Note that there are two elements returned if use inner_core_mode. 
        return dna_mask, inner_core_mask

    



# function to get cytoplasm mask

def cell_mask_segmentation (img_each_fov, cell_ch_index, dna_mask =np.zeros([3,3]), 
                                      object_size_th = 10000, cell_center = [466,621], cell_size=600,
                             threshold_method = skimage.filters.thresholding.threshold_otsu,
                             cytoplasm_ratio = 60, segmentation_by_percentile = False, single_obj=True,
                                  threshold_by_region = "GLOBAL",):
    
    ''' A function to generate a mask where myo5 (or other proteins with similar localization pattern) is enriched. Other proteins can be gH2AX foci, RPA foci where foci of interest are visually (and presumbaly signficantly) higher than the diffusive cytoplasmic/nucleoplasmic background.
      Inputs:
              img_each_fov: a list of ch - xy io_read ims arrays;
              myo_ch_index: int; the ch index (starting from 1) for the channel of interest;
              dna_mask: dna mask if excluding dna region for the channel of interest; mask provided by the user, which be can be generated by the function chromosome_mask_segmentation (INPUTS see the description of the function).
              object_size_th: int; size th for excluding objects small than the typical size of the object of interest.
              threshold_method: default is as shown above (otsu);
              cell_center: [x,y] coord of the center of the cropped region;
              cell_size: int; length of the cropped square;
              threshold_by_region =  "GlOBAL": use the whole image for segmentation th;
              threshold_by_region =  "CELL": use the whole image for segmentation th, together with using a cell mask to facilitate the segmentation;
              threshold_by_region =  "USER": use the cropped region for segmentation th; the image will be cropped during analysis; for many scenarios where only one mitotic cell is desired when other cells are present, this mode seems most effective for segmenation.
              DNA_exclusion: whether to use DNA mask to exclude the condensed mitotic chromosome regions.
              num_of_obj: the number of object desried to be detected.
              cytoplasm_ratio: the area percentile within the FOV for the cell cytoplasm.

      Outputs:
             the myo enriched mask, which is the region myo (or other protein of interest) is enriched --- > mean + 3* std.
             
             Note that when USER mode is used, the output image/mask size will be cropped. Thus, be consistent when using USER mode for all relevant functions.'''

    #import skimage
    #from skimage import io
    #import numpy as np
    #from skimage.filters.thresholding import threshold_li, threshold_otsu, threshold_triangle
    #from skimage.morphology import erosion, dilation, opening, closing, white_tophat, remove_small_objects
    #from skimage.morphology import disk
    #from scipy import ndimage as ndi
    #from scipy.ndimage.morphology import binary_fill_holes
    #from skimage.measure import label, regionprops, regionprops_table
    #import math
    #import matplotlib.pyplot as plt
    
    
    cell_image = img_each_fov[cell_ch_index-1]
    
    # If use USER mode, perform thresholding in the cropped image as specified.
    if threshold_by_region == "USER":
        if len(cell_center) !=2:
                print ('Specify the center to generate a crop.')
                return None
        
        cell_image = np.array(cell_image)[int(cell_center[1]-cell_size/2):int(cell_center[1]+cell_size/2),  
              int(cell_center[0]-cell_size/2):int(cell_center[0]+cell_size/2)]
           
    # Thresholding and binary operations.
    if segmentation_by_percentile == False:
        segmentation_th = threshold_method (cell_image)
    else:
        segmentation_th = np.percentile(cell_image, cytoplasm_ratio)

    cell_mask = cell_image > segmentation_th
    cell_mask = remove_small_objects(cell_mask, object_size_th,connectivity=1)
    cell_mask = binary_fill_holes(cell_mask)
    cell_mask = erosion(cell_mask, disk(3))

    if single_obj and np.sum(cell_mask)>0:
        # if one cell to be selected, get largest one for the cell mask
        labeled_roi, num_of_roi  = ndi.label(cell_mask)
        cand_roi_size_list = []
        for i in range(num_of_roi):
            cand_roi = labeled_roi == i+1
            cand_roi_size = np.sum(cand_roi) 
            cand_roi_size_list.append(cand_roi_size)
        sel_roi_ind = np.argmax(cand_roi_size_list)
        sel_roi = labeled_roi == sel_roi_ind+1
        cell_mask = sel_roi

    return cell_mask



# function to get the enriched portion within the cell

def myo_enrichment_mask_segmentation (img_each_fov, myo_ch_index, dna_mask =np.zeros([3,3]), 
                                      object_size_th = 10000, cell_center = [466,621], cell_size=600,
                             threshold_method = skimage.filters.thresholding.threshold_otsu,
                             cytoplasm_ratio = 60,
                                  threshold_by_region = "GLOBAL", DNA_exclusion =True, num_of_obj =2, std_ratio = 3):
    
    ''' A function to generate a mask where myo5 (or other proteins with similar localization pattern) is enriched. Other proteins can be gH2AX foci, RPA foci where foci of interest are visually (and presumbaly signficantly) higher than the diffusive cytoplasmic/nucleoplasmic background.
      Inputs:
              img_each_fov: a list of ch - xy io_read ims arrays;
              myo_ch_index: int; the ch index (starting from 1) for the channel of interest;
              dna_mask: dna mask if excluding dna region for the channel of interest; mask provided by the user, which be can be generated by the function chromosome_mask_segmentation (INPUTS see the description of the function).
              object_size_th: int; size th for excluding objects small than the typical size of the object of interest.
              threshold_method: default is as shown above (otsu);
              cell_center: [x,y] coord of the center of the cropped region;
              cell_size: int; length of the cropped square;
              threshold_by_region =  "GlOBAL": use the whole image for segmentation th;
              threshold_by_region =  "CELL": use the whole image for segmentation th, together with using a cell mask to facilitate the segmentation;
              threshold_by_region =  "USER": use the cropped region for segmentation th; the image will be cropped during analysis; for many scenarios where only one mitotic cell is desired when other cells are present, this mode seems most effective for segmenation.
              DNA_exclusion: whether to use DNA mask to exclude the condensed mitotic chromosome regions.
              num_of_obj: the number of object desried to be detected.
              std_ratio: the factor to multiply for the positive th: set higher for more enriched signal.
              cytoplasm_ratio: the area percentile within the FOV for the cell cytoplasm.

      Outputs:
             the myo enriched mask, which is the region myo (or other protein of interest) is enriched --- > mean + 3* std.
             
             Note that when USER mode is used, the output image/mask size will be cropped. Thus, be consistent when using USER mode for all relevant functions.'''

    #import skimage
    #from skimage import io
    #import numpy as np
    #from skimage.filters.thresholding import threshold_li, threshold_otsu, threshold_triangle
    #from skimage.morphology import erosion, dilation, opening, closing, white_tophat, remove_small_objects
    #from skimage.morphology import disk
    #from scipy import ndimage as ndi
    #from scipy.ndimage.morphology import binary_fill_holes
    #from skimage.measure import label, regionprops, regionprops_table
    #import math
    #import matplotlib.pyplot as plt
    
    
    myo_image = img_each_fov[myo_ch_index-1]



    
    # If use USER mode, perform thresholding in the cropped image as specified.
    if threshold_by_region == "USER":
        if len(cell_center) !=2:
                print ('Specify the center to generate a crop.')
                return None
        
        myo_image = np.array(myo_image)[int(cell_center[1]-cell_size/2):int(cell_center[1]+cell_size/2),  
              int(cell_center[0]-cell_size/2):int(cell_center[0]+cell_size/2)]
           
    # Thresholding and binary operations.
    #segmentation_th = threshold_method (myo_image)
    segmentation_th = np.percentile(myo_image, cytoplasm_ratio)
    cell_mask = myo_image > segmentation_th
    cell_mask = remove_small_objects(cell_mask, object_size_th,connectivity=1)
    cell_mask = binary_fill_holes(cell_mask)
    cell_mask = erosion(cell_mask, disk(3))
    
    # Exclude DNA (eg, condensed mitotic chromosomes) or not accordingly for calculating the intensity distribution histogram.
    if DNA_exclusion == True:
        if dna_mask.shape != myo_image.shape:
            print ('A valid DNA mask is missing.')
            return None
        else:
            myo_cytoplasm = cell_mask * (dna_mask ==0)
    else:
        myo_cytoplasm = cell_mask
    
    # If only one object of interest is desired to be detected, this below will help check if segmentation went well or not.
    if num_of_obj == 1:
        
        labeled_cell, num_of_cell  = ndi.label(myo_cytoplasm)
        # Remove edge objects
        kept_cell = []
        for i in range(num_of_cell):
            cand_cell = labeled_cell == i+1
            if cand_cell[0,0] == 0 and cand_cell[0,len(cand_cell)-1] == 0 and cand_cell[len(cand_cell)-1,0] == 0 and cand_cell[len(cand_cell)-1,len(cand_cell)-1] == 0:
                cand_cell[cand_cell>0]=1
                cand_cell = np.array(cand_cell)
                kept_cell.append(cand_cell)
                
        if len(kept_cell) !=1:
            print (f"Fail to detect one cell object. There are/is {len(kept_cell)} object(s) detected.")
            return None

    # Calculate the intensity distribution histogram for the myo within the cell/nucleus object of interest.
    myo_intensity_filtered = np.array([i for i in (myo_image * myo_cytoplasm).flatten() if i>0])
    
    # Use mean + 3 * std as the th for enrichment.
    myo_positive_th = myo_intensity_filtered.mean() + std_ratio* myo_intensity_filtered.std()
    print (f"Using {round(myo_positive_th,2)} to segment the myo enriched region.")

    myo_enriched = (myo_image>myo_positive_th)*myo_cytoplasm
    
    # Return the mask for the enrichment.
    return myo_enriched
        
    






    # function to get the ER mask or mask for other ridge-like signal (e.g., actin) 

def er_membrane_mask_segmentation (img_each_fov, er_ch_index,kwargs = {'sigmas': [1], 'mode': 'reflect'}, ridge_sigma = 2,
                                      cell_center = [466,621], cell_size=600, percentile_multiplier = 10,
                                   threshold_method = skimage.filters.thresholding.threshold_otsu,sm_fragment_size = 50,
                                   hot_pixel_percentile = 99.8, non_ridge_segmentation = False, 
                                   enriched_membrane = 'NE', non_ridge_threshold_method = skimage.filters.thresholding.threshold_multiotsu,
                                  threshold_by_region = "GLOBAL", cytoplasm_mask = np.zeros ([3,3]), generate_cytoplasm_mask =True):
        
    ''' A function to generate a mask where visually appearnt membrane (e.g., mitotic ER and NE) is present.
      Inputs:
              img_each_fov: a ch - xy io_read ims arrays;
              er_ch_index: int; the ch index (starting from 1) for the membrane channel of interest;
              kwargs: kwargs for the ridge detection module to be used.
              cell_center: [x,y] coord of the center of the cropped region;
              cell_size: int; length of the cropped square;
              percentile_multiplier: int; a value for approximating the overall intensity of ridge signal to a good th that can distinguish and retain majority of membrane fraction.
              sm_fragment_size: the object size for small fragments to be removed after ridge final thresholding;
              threshold_method: default is as shown above (otsu) for cytoplasm mask generation;
              threshold_by_region =  "GlOBAL": use the whole image for segmentation th;
              threshold_by_region =  "CELL": use the whole image for segmentation th, together with using a cell mask to facilitate the segmentation;
              threshold_by_region =  "USER": use the cropped region for segmentation th; the image will be cropped during analysis; for many scenarios where only one mitotic cell is desired when other cells are present, this mode seems most effective for segmenation.
              generate_cytoplasm_mask: whether to generate cytoplasm mask or to use the provided cytoplasm mask (eg obtained from other channel).
              cytoplasm_mask: the cytoplasm mask provided to help segmentation.

      Outputs:
             the myo enriched mask, which is the region myo (or other protein of interest) is enriched --- > mean + 3* std.
             
             Note that when USER mode is used, the output image/mask size will be cropped. Thus, be consistent when using USER mode for all relevant functions.'''

    #import skimage
    #from skimage import io
    #import numpy as np
    #from skimage import data
    #from skimage import color
    #from skimage.filters import meijering, sato, frangi, hessian
    #import matplotlib.pyplot as plt
    #from skimage.filters.thresholding import threshold_li, threshold_otsu, threshold_triangle
    #from skimage.morphology import erosion, dilation, opening, closing, white_tophat, remove_small_objects
    #from skimage.morphology import disk
    #from scipy import ndimage as ndi
    #from scipy.ndimage.morphology import binary_fill_holes
    
    er_image = img_each_fov [er_ch_index -1]



    # Remove hot pixel 
    hot_pixel_th = np.percentile(er_image,hot_pixel_percentile)
    # Deepcopy to prevent direct altering of the loaded images:
    import copy
    er_image_copy = copy.deepcopy(er_image)
    er_image_copy[er_image>hot_pixel_th] = np.median(er_image)
    er_image = er_image_copy

    # If use USER mode, perform thresholding in the cropped image as specified.
    if threshold_by_region == "USER":
        if len(cell_center) !=2:
                print ('Specify the center to generate a crop.')
                return None
        
        er_image = np.array(er_image)[int(cell_center[1]-cell_size/2):int(cell_center[1]+cell_size/2),  
              int(cell_center[0]-cell_size/2):int(cell_center[0]+cell_size/2)]
    
    # if not using ridge detectioin (e.g., for low resolution 60x ims)
    if non_ridge_segmentation:
        print ('Generate ER mask without ridge detection.')
        if enriched_membrane == 'NE' and non_ridge_threshold_method==skimage.filters.thresholding.threshold_multiotsu:
            ('The enriched membrane signal is mostly at NE as specified, use the first threshold.')
            er_mask_th = non_ridge_threshold_method (er_image, classes=4)
            er_ridge_mask = er_image > er_mask_th[2]
        else:
            ('The enriched membrane signal is mostly at ER as specified, use the second threshold.')
            er_mask_th = non_ridge_threshold_method (er_image, classes=4)
            er_ridge_mask = er_image > er_mask_th[1]
        return er_ridge_mask

    # Use meijering ridge detection module ["BLACK RIDGES"] to find the majority of membrane fractions.
    kwargs['black_ridges'] = False
    kwargs['sigmas'] = ridge_sigma
    er_ridge = meijering(er_image, **kwargs)

    # Use and/or generate the cytoplasm mask to keep only the membrane within the cell of interest.
    if generate_cytoplasm_mask == True:
        segmentation_th = threshold_method (er_image)
        cytoplasm_mask = er_image > segmentation_th
        cytoplasm_mask = remove_small_objects(cytoplasm_mask, 10000,connectivity=1)
        cytoplasm_mask = binary_fill_holes(cytoplasm_mask)
        cytoplasm_mask = erosion(cytoplasm_mask, disk(3))
        er_cytoplasm_ridge_intensity =  er_ridge * cytoplasm_mask
              
    else:
        if cytoplasm_mask.shape != er_image.shape:
            print ('A valid cytoplasm mask is missing.')
            return None
        else:
            er_cytoplasm_ridge_intensity =  er_ridge * cytoplasm_mask
    
    er_cytoplasm_ridge_filtered = np.array([i for i in er_cytoplasm_ridge_intensity.flatten() if i >0])

    # error handling for failed ER ridge detection
    if not len(er_cytoplasm_ridge_filtered)>0:
        er_ridge_mask = np.zeros(er_image.shape)
        print('ER ridge mask generation fails, return empty mask.')
        return er_ridge_mask


    #upper_signal = np.max(er_cytoplasm_ridge_filtered) 
    #if upper_signal != 1:
        #print ('Something needs attention here. The max value for the ridge image should be 1?')
    
    # Emperically, a perfect th cutoff appears to be around 0.03-0.15. Values higher than the perfect leads to a big loss of membrane fractions. 
    # Values (slightly) lower than the perfect th retain the membrane but intorduces small fragments which can be removed based on their small size.
    # Thus, it seems reasonable to start with a lower th to keep as much as valid membrane as possible. 
    # This starting th value appears to be around the 10* lower 10% signal (which also seems to be always lower than the perfect th, thus avoding loss of membrane.)

    lower_signal = np.percentile(er_cytoplasm_ridge_filtered,q=0.1)
    
    ridge_th = lower_signal *percentile_multiplier
    er_ridge_mask = er_ridge > ridge_th

    # Iterate to increase the ridge_th in case the ridge_th is too small.
    er_ridge_regions =  regionprops(label(er_ridge_mask))
    sm_fragment_count = 0
    for _region in er_ridge_regions:
        if _region.area <= sm_fragment_size:
            sm_fragment_count +=1
    
    if sm_fragment_count > 2000:
        ridge_th *= 1.5
        er_ridge_mask = er_ridge > ridge_th
        er_ridge_regions =  regionprops(label(er_ridge_mask))
        sm_fragment_count = 0
        for _region in er_ridge_regions:
            if _region.area <= sm_fragment_size:
                sm_fragment_count +=1
    
    if sm_fragment_count > 1000:
        ridge_th *= 1.1
        er_ridge_mask = er_ridge > ridge_th
        er_ridge_regions =  regionprops(label(er_ridge_mask))
        sm_fragment_count = 0
        for _region in er_ridge_regions:
            if _region.area <= sm_fragment_size:
                sm_fragment_count +=1

    
    if sm_fragment_count < 200:
        ridge_th *= 0.9
        er_ridge_mask = er_ridge > ridge_th
        er_ridge_regions =  regionprops(label(er_ridge_mask))
        sm_fragment_count = 0
        for _region in er_ridge_regions:
            if _region.area <= sm_fragment_size:
                sm_fragment_count +=1

    print (f'There are {sm_fragment_count} small membrane fragments (sigma = {ridge_sigma}) removed after converting detected ridges to binary images.')   

    er_ridge_mask = remove_small_objects (er_ridge_mask,sm_fragment_size)
    
    return er_ridge_mask
    




# function to bins the midbox, which can be used to crop the inner core mask into desired number of segment.
def generate_midbox_with_bins (dna_mask, num_of_bins=5,
                                          _verbose = True):
    
    " binned mask for xy image"

    import skimage
    # check if the input chr mask is valid
    regions = regionprops(label(dna_mask))

    if len(regions)!=2:
        print (f'There are less or more than two mitotic chromosome mass detected, exit.')
        return None  # or some other empty thing for downstream analysis to read
    
    # Gnerate the midbox_mask with bins using chromosome masks
    # result dict to store all bins of midbox mask
    binned_midbox_mask = {}
    # generate and store each bins
    for _bin_ind in range(num_of_bins):
        midbox = []
        for props in regions:
            y0, x0 = props.centroid
            orientation = props.orientation
            # x1y1, x2y2 as boundary point
            x1 = x0 - math.sin(orientation) * 0.5 * props.major_axis_length 
            y1 = y0 - math.cos(orientation) * 0.5 * props.major_axis_length
            x2 = x0 + math.sin(orientation) * 0.5 * props.major_axis_length
            y2 = y0 + math.cos(orientation) * 0.5 * props.major_axis_length
            # x3y3, x4y4 for coords of the bin boundary
            x3 = x1+ ((x2-x1)/num_of_bins*_bin_ind)
            y3 = y1+ ((y2-y1)/num_of_bins*_bin_ind)
            x4 = x1+ ((x2-x1)/num_of_bins*(_bin_ind+1))
            y4 = y1+ ((y2-y1)/num_of_bins*(_bin_ind+1))
            # midbox bin
            midbox.append ([x3,y3])
            midbox.append ([x4,y4])
        # invert xy for the mask
        import copy   
        sorted_midbox =copy.deepcopy(midbox)
        sorted_midbox[-1] = midbox[-2]
        sorted_midbox[-2] = midbox[-1]
        sorted_midbox.append(midbox[0])
        sorted_midbox = np.array(sorted_midbox)
        sorted_midbox[:,[1,0]] = sorted_midbox[:,[0,1]]
        midbox_mask=skimage.draw.polygon2mask(dna_mask.shape, sorted_midbox)
        # add each mask to the result dict
        binned_midbox_mask[_bin_ind] = midbox_mask
    if  _verbose:
        print (f'Midboxes dict with {num_of_bins} bins are generated.')
    return binned_midbox_mask


# function to get the inner-core, outer-core, and noncore for single chromosome mask as input. Requires a midbox (e.g., from above) as input for the inner-core of interest.
def generate_core_noncore_for_chromosome (chromosome_mask, midbox_for_core, core_thickness=10, chr_rim_adjust_size = 0, 
                                            region_adjust_size=5, erode_adjust = True):
    

    # adjust chromosme rim accordingly
    chromosome_rim_mask_outer = dilation (chromosome_mask, disk(core_thickness-chr_rim_adjust_size))
    if erode_adjust:
        chromosome_rim_mask_inner = erosion (chromosome_mask, disk(chr_rim_adjust_size)) 
    else:
        chromosome_rim_mask_inner = dilation (chromosome_mask, disk(chr_rim_adjust_size)) 
    chromosome_rim_mask = chromosome_rim_mask_outer * (chromosome_rim_mask_inner ==0)

    # 1. Get inner_core mask
    inner_core_mask = chromosome_rim_mask * midbox_for_core
    
    # 2. Get corebox for the outer_core mask
    regions = regionprops(label(chromosome_mask))
    # Proceed if single chromosome mask
    if len(regions)==1:
        for props in regions:
            y0, x0 = props.centroid
            orientation = props.orientation
            x3 = x0 + math.sin(orientation+np.pi*0.5) * 1 * props.minor_axis_length 
            y3 = y0 + math.cos(orientation+np.pi*0.5) * 1 * props.minor_axis_length
            x4 = x0 - math.sin(orientation+np.pi*0.5) * 1 * props.minor_axis_length 
            y4 = y0 - math.cos(orientation+np.pi*0.5) * 1 * props.minor_axis_length
            # get 4 coordindates for the core box (use the middle 0.2+0.2~0.4 of the major axis length to define the core)
            x3_1 = x3 - math.sin(orientation) * 0.2 * props.major_axis_length
            y3_1 = y3 - math.cos(orientation) * 0.2 * props.major_axis_length
            x3_2 = x3 + math.sin(orientation) * 0.2 * props.major_axis_length
            y3_2 = y3 + math.cos(orientation) * 0.2 * props.major_axis_length
            x4_1 = x4 - math.sin(orientation) * 0.2 * props.major_axis_length
            y4_1 = y4 - math.cos(orientation) * 0.2 * props.major_axis_length
            x4_2 = x4 + math.sin(orientation) * 0.2 * props.major_axis_length
            y4_2 = y4 + math.cos(orientation) * 0.2 * props.major_axis_length
            core_box = [[x3_1,y3_1],[x3_2,y3_2],[x4_2,y4_2],[x4_1,y4_1]]
            # generate the core box for the chromosome mask
            import copy   
            core_box.append(core_box[0])
            sorted_core_box = np.array(core_box)
            sorted_core_box[:,[1,0]] = sorted_core_box[:,[0,1]]
            core_box_mask=skimage.draw.polygon2mask(chromosome_mask.shape, sorted_core_box)
        
    else: 
        print ('Too many objects. Core-Noncore segmenation failed, return empty mask.')
        empty_mask = np.zeros(chromosome_mask.shape)
        # retunr empty mask for inner_core, outer_core, and noncore respectively
        return {'inner_core':empty_mask,'outer_core':empty_mask,'noncore':empty_mask}
    
    # 3. Get outer core mask 
    outer_core_mask =  chromosome_rim_mask  * (dilation(inner_core_mask,disk(region_adjust_size))==0) * core_box_mask
    
    # 4. Get noncore mask by excluding both inner and outer core
    noncore_mask =  chromosome_rim_mask  * ((dilation(inner_core_mask,disk(region_adjust_size))==0)*
                                            (dilation(outer_core_mask,disk(region_adjust_size))==0))
    
    if (outer_core_mask.shape == inner_core_mask.shape and
        outer_core_mask.shape == noncore_mask.shape and
        outer_core_mask.shape == chromosome_mask.shape):
        print ('Core-Noncore segmenations are generated for the input chromosome mask.')
        return  {'inner_core':inner_core_mask,'outer_core':outer_core_mask,'noncore':noncore_mask}
        
    else:
        print ('Core-Noncore segmenation failed, return empty mask.')
        empty_mask = np.zeros(chromosome_mask.shape)
        # retunr empty mask for inner_core, outer_core, and noncore respectively
        return  {'inner_core':empty_mask,'outer_core':empty_mask,'noncore':empty_mask}



def object_segmentation_timelapsed_3D (sorted_mm_data_array_3D, object_ch_index, 
                            threshold_method = skimage.filters.thresholding.threshold_otsu, 
                            threshold_by_region = "GLOBAL", cell_center = [466,621], cell_size = 600,
                            large_object_size = 5000, object_type = "micronuclei", midbox_mask = None):     # object size in voxel -- for 3D

    """Function to generate the 3D masks for objects of interest for timelapsed data"""

    sel_object_mask_all = []

    for _data_each_time in sorted_mm_data_array_3D:
    
    
        image_ch_object = _data_each_time[object_ch_index-1]  # object (e.g., dna channel)

        # Use the whole image 
        if threshold_by_region == "GLOBAL":
            image_ch_crop = image_ch_object

        # Use the cropped region of the image as specified by the user.
        if threshold_by_region == "USER":
            # perform segmentation using a cropped image, where the crop box is specified by the user
            # use this method when cell of interest cannot be determined by features in image.
            if len(cell_center) !=2:
                    print ('Specify the center to generate a crop.')
                    return None
            else:
                image_ch_crop = image_ch_object[:,    # channel z
                         int(cell_center[1]-cell_size/2):int(cell_center[1]+cell_size/2),    # channel crop x
                         int(cell_center[0]-cell_size/2):int(cell_center[0]+cell_size/2)]   # channel crop y
    
    
        # remove hot pixels
        hot_pixel_th = np.percentile(image_ch_crop,99.8)
        import copy
        image_ch_crop_copy = copy.deepcopy(image_ch_crop)
        image_ch_crop_copy[image_ch_crop>hot_pixel_th] = np.median(image_ch_crop)
        image_ch_crop = image_ch_crop_copy
        # initial segmentation
        all_object_mask_th = threshold_method(image_ch_crop)
        all_object_mask = image_ch_crop > all_object_mask_th
        all_object_mask = remove_small_objects(all_object_mask, 100 ,connectivity=1)
        
        # keep small chromosome objects only
        if object_type == "micronuclei":
        
            large_object_mask = remove_small_objects(all_object_mask, large_object_size ,connectivity=1)
            sel_object_mask = np.logical_xor(all_object_mask, large_object_mask)
      
            if midbox_mask is None:
                print ("-- No midbox provided, select all micronuclei")
            # use midbox_mask to keep laggards within the midzone
            else:
                #print ("-- Use midbox to select micronuclei")
                sel_object_mask = sel_object_mask * midbox_mask

        # keep large chromosmoe objects only
        if object_type == "main_nuclei":
        
            large_object_mask = remove_small_objects(all_object_mask, large_object_size ,connectivity=1)
            sel_object_mask = large_object_mask
        
        # append the 3D mask for the timepoint
        sel_object_mask_all.append(sel_object_mask)

    sel_object_mask_all = np.array(sel_object_mask_all) 

    return sel_object_mask_all







def simple_find_objects_timelapsed_3D (sel_object_mask_all, combine_objects = True,
                                        search_radius = 75, object_area_change_ratio = 0.5,
                                        object_erosion_factor =3):

    '''Function 
    
     Output as {t0:{'object_1':[[z, mask],[z, mask],[z, mask]]},
                t1:{'object_1':[[z, mask],[z, mask],[z, mask],[z, mask]]}
                ...
                tn:{'object_1':[[z, mask],[z, mask]]}}   
                
                where t is a integer of timepoint and object_1 (, object_2, etc) are objects used for other measurements;
                                        each [z, mask] are sub-objects of the object'''
            
    props_all_times = {}

    object_info_for_tracking_dict = {}

    
    if not combine_objects:
        print('-- Do simple track and select objects over time using spatial information.')
        #print('-- Will track for main chromosome only')

        for _time_index, mask_each_time in enumerate(sel_object_mask_all):

        
            props_all_times [_time_index] = {}
            num_of_z = len(mask_each_time)
        
            label_mask = label(mask_each_time)
            props = regionprops_table(label_mask,properties=['label','area','centroid'])
          
                            
            object_info_for_tracking_dict[_time_index] = []

                # set the initial labels for timepoint 0
            if _time_index == 0:
                    
                if len(props['label']) !=2:
                    print ('-- More or less than two chromosome objects identifed')
                    return None
                    
                else:

                    for _label_index in range(len(props['label'])):

                        props_all_times [_time_index] [f"object_{_label_index+1}"] = []

                        # get the estimated z center for the object
                        _z_sel = int(round(props['centroid-0'][_label_index]))
                        if _z_sel > num_of_z -1:
                            _z_sel = num_of_z -1
                        if _z_sel < 0:
                            _z_sel = 0
                        # get the object mask for the z
                        sel_object_mask = label_mask[_z_sel] == props['label'][_label_index]
                        sel_object_mask = binary_fill_holes(sel_object_mask)
                        sel_object_mask = dilation(sel_object_mask, disk(object_erosion_factor)) 
                        sel_object_mask = erosion(sel_object_mask, disk(object_erosion_factor*2)) 
                        # add z-info and mask for each sub-object
                        props_all_times [_time_index] [f"object_{props['label'][_label_index]}"].append([_z_sel, sel_object_mask])
                        
                        # record other spatial info for object tracking
                        object_z,object_x,object_y =props['centroid-0'][_label_index],props['centroid-1'][_label_index],props['centroid-2'][_label_index]
                        object_area= props['area'][_label_index]

                        object_info_for_tracking_dict[_time_index].append([object_z,object_x,object_y,object_area])
                        
                # for later timepoints                                                                    
            else:
                    
                #print(object_info_for_tracking_dict[_time_index-1])
                    
                for _prior_object_index, _prior_object in enumerate(object_info_for_tracking_dict[_time_index-1]):  # always use the last saved info

                    props_all_times [_time_index] [f"object_{_prior_object_index+1}"] = []
                    # get prior info
                    _prior_object_zxy = _prior_object[:3]
                    _prior_object_area = _prior_object[3]
                    
                    # set a flag to record if the object has been found
                    _found_flag = False
                        
                    # loop through all new labels
                    for _label_index in range(len(props['label'])):  

                        object_z,object_x,object_y =props['centroid-0'][_label_index],props['centroid-1'][_label_index],props['centroid-2'][_label_index]
                        object_area= props['area'][_label_index]
                        # if object is within the search radius and has proper size
                        if np.linalg.norm(np.array(_prior_object_zxy)-np.array([object_z,object_x,object_y])) < search_radius:
                            if abs(object_area - _prior_object_area) < _prior_object_area * object_area_change_ratio:
                                # get the estimated z center for the object
                                _z_sel = int(round(props['centroid-0'][_label_index]))
                                if _z_sel > num_of_z -1:
                                    _z_sel = num_of_z -1
                                if _z_sel < 0:
                                    _z_sel = 0
                                # get the object mask for the z
                                sel_object_mask = label_mask[_z_sel] == props['label'][_label_index]
                                sel_object_mask = binary_fill_holes(sel_object_mask)
                                sel_object_mask = remove_small_objects(sel_object_mask, 500 ,connectivity=1)
                                sel_object_mask = dilation(sel_object_mask, disk(object_erosion_factor)) 
                                sel_object_mask = erosion(sel_object_mask, disk(object_erosion_factor*2)) 
                                
                                # add z-info and mask for each sub-object
                                props_all_times [_time_index] [f"object_{_prior_object_index+1}"].append([_z_sel, sel_object_mask])

                                # record other spatial info for object tracking
                                object_info_for_tracking_dict[_time_index].append([object_z,object_x,object_y,object_area])                               
                                # object found
                                _found_flag = True

                    # if object not found, append [None, None] for [z, mask]
                    if _found_flag == False:                    
                        props_all_times [_time_index] [f"object_{_prior_object_index+1}"].append([None, None])
                   
                # set the spatial record to the previous one if not both two chromosome masses were succesfully located
                #print(f'DICT LEN {len(object_info_for_tracking_dict[_time_index])}')
                if len(object_info_for_tracking_dict[_time_index]) !=2:
                    object_info_for_tracking_dict[_time_index] = object_info_for_tracking_dict[_time_index-1]

        

    if combine_objects:
        print('-- group all objects into one')

        for _time_index, mask_each_time in enumerate(sel_object_mask_all):

        
            props_all_times [_time_index] = {}
            num_of_z = len(mask_each_time)
        
            label_mask = label(mask_each_time)
            props = regionprops_table(label_mask,properties=['label','area','centroid'])

            props_all_times [_time_index] ['object_1'] = []
            
            # if has any object
            if len(props['label']) > 0:
                for _label_index in range(len(props['label'])):
                
                    # get the estimated z center for the object
                    _z_sel = int(round(props['centroid-0'][_label_index]))
                    if _z_sel > num_of_z -1:
                        _z_sel = num_of_z -1
                    if _z_sel < 0:
                        _z_sel = 0
            
                    # get the object mask for the z
                    sel_object_mask = label_mask[_z_sel] == props['label'][_label_index]
                    sel_object_mask = remove_small_objects(sel_object_mask, 10 ,connectivity=1)
                    sel_object_mask = binary_fill_holes(sel_object_mask)
                    sel_object_mask = dilation(sel_object_mask, disk(1)) 
                    sel_object_mask = erosion(sel_object_mask, disk(1)) 
                    # add z-info and mask for each sub-object
                    props_all_times [_time_index] ['object_1'].append([_z_sel, sel_object_mask])
            
            # if has no object
            else:
                props_all_times [_time_index] ['object_1'].append([None, None])
                
        
    return props_all_times





def split_cell_mask (cell_mask, dna_mask, cell_midbox_erode=30):

    # check main chromosome mass number
    regions = regionprops(label(dna_mask))
    if len(regions)!=2:
        print (f'Fail to detect 2 chromosome objects. There are/is {len(regions)} object(s) detected')

    # get 4 coords on major axis
    midbox = []
    for props in regions:
        y0, x0 = props.centroid
        orientation = props.orientation
        x1 = x0 - math.sin(orientation) * 0.5 * props.major_axis_length
        y1 = y0 - math.cos(orientation) * 0.5 * props.major_axis_length
        x2 = x0 + math.sin(orientation) * 0.5 * props.major_axis_length
        y2 = y0 + math.cos(orientation) * 0.5 * props.major_axis_length
        midbox.append ([x1,y1])
        midbox.append ([x2,y2])
    # get centroid of the 4 coords above
    center_midbox = np.mean(midbox, axis=0)
    
    # # get 4 coords on minor axis and select the two (one for each main chromosomes) closest to the centroid above
    cell_midbox = []
    for props in regions:
        y0, x0 = props.centroid
        orientation = props.orientation
        x3 = x0 + math.cos(orientation) * 0.5 * props.minor_axis_length
        y3 = y0 - math.sin(orientation) * 0.5 * props.minor_axis_length
        x4 = x0 - math.cos(orientation) * 0.5 * props.minor_axis_length
        y4 = y0 + math.sin(orientation) * 0.5 * props.minor_axis_length
        _dist_ct = []
        for _xy in [[x3,y3],[x4,y4]]:
            _dist_ct.append(np.linalg.norm(_xy - center_midbox))
        _xy_ind = np.argmin(np.array(_dist_ct))
        _x_sel, _y_sel =[[x3,y3],[x4,y4]][_xy_ind]     
        
        # extend the selected coord both way along the orientation; 3x to ensure it's exceeding the cell size
        _x_sel1 = _x_sel - math.sin(orientation) * 3 * props.major_axis_length 
        _y_sel1 = _y_sel - math.cos(orientation) * 3 * props.major_axis_length 
        _x_sel2 = _x_sel + math.sin(orientation) * 3 * props.major_axis_length 
        _y_sel2 = _y_sel + math.cos(orientation) * 3 * props.major_axis_length 
        
        cell_midbox.append([_x_sel1, _y_sel1])
        cell_midbox.append([_x_sel2, _y_sel2])

    # get convex hull for the new cell midzone box to sort the coord
    from scipy.spatial import ConvexHull 
    hull = ConvexHull(cell_midbox)
    sorted_cell_midbox= []
    for _v in hull.vertices:
        sorted_cell_midbox.append(cell_midbox[_v])

    # invert xy for mask conversion
    sorted_cell_midbox = np.array(sorted_cell_midbox)
    sorted_cell_midbox[:,[1,0]] = sorted_cell_midbox[:,[0,1]]
    # convert mask
    cell_midbox_mask=skimage.draw.polygon2mask(dna_mask.shape, sorted_cell_midbox)
    
    # get the split with adjusting the width (and length but in this case it does not affect)
    # output is one midzone, one edge combining each side
    cell_midbox_mask_erode = erosion(cell_midbox_mask, disk(cell_midbox_erode))
    cell_midzone = cell_midbox_mask_erode*cell_mask
    cell_edge = ~cell_midbox_mask_erode*cell_mask
    

    return cell_midzone, cell_edge