# -*- coding: utf-8 -*-
"""
Created on Sat Jul 10 15:58:59 2021

@author: Shiwei Liu @ Harvard University
"""



# Setup file for AnalysisTool


import pip

def import_or_install(package):
    try:
        __import__(package)
    except ImportError:
        pip.main(['install', package]) 


packages = ['scipy','skimage','sklearn','PIL', 'cv2','matplotlib','treelib','pandas','roifile','yaml']


for _package in packages:
    import_or_install(_package)



import sys, os, glob, time, copy
import shutil
import numpy as np
import scipy
import skimage
import math
import sklearn
from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import roifile
import yaml
#import pickle
#from IPython.display import clear_output
#import multiprocessing
#import psutil
import cv2
#import seaborn as sns
# reload
#from importlib import reload
# add Document to path
sys.path.append(r'C:\Users\Shiwei\Documents')

