from __future__ import print_function, unicode_literals, absolute_import, division
import sys
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from tifffile import imread
from csbdeep.utils import Path, normalize
from csbdeep.io import save_tiff_imagej_compatible
from stardist import random_label_cmap
from stardist.models import StarDist3D
import os
from os import listdir
from os.path import isfile, join


N_channels = 4
dapi_index = 0
model = StarDist3D(None, name='clem_model', basedir='C:/Users/fabbe/Desktop/Segmentation Clementine 06-10-20/')
folder = 'C:/Users/fabbe/Desktop/Segmentation Clementine 06-10-20/scale/batch/'


np.random.seed(6)
lbl_cmap = random_label_cmap()    
    
img_files = [f for f in listdir(folder) if isfile(join(folder, f)) and f.endswith(".tif") ]

output_folder = folder + '/masks/'
if not os.path.isdir(folder + '/masks/'):
    os.mkdir(output_folder)

for i in range(len(img_files)):
 
    print(i)
    x = imread(folder + '/' + img_files[i])
    if N_channels > 1:
        x = x[:,dapi_index,:,:]
    
    n_channel = 1 if x[0].ndim == 3 else x[0].shape[-1]
    axis_norm = (0,1,2)   # normalize channels independently
    # axis_norm = (0,1,2,3) # normalize channels jointly
    if n_channel > 1:
        print("Normalizing image channels %s." % ('jointly' if axis_norm is None or 2 in axis_norm else 'independently'))   
        
    img = normalize(x, 1,99.8, axis=axis_norm)
    
    labels, details = model.predict_instances(img,n_tiles=[1,3,3])
    save_tiff_imagej_compatible(output_folder + img_files[i].replace('.tif', '') + '_mask' + '.tif', labels.astype(np.int16), axes='ZYX')
    


