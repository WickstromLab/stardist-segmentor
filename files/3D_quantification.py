# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 13:11:01 2020

@author: fabbe
"""

from os import listdir
from os.path import isfile, join
from tifffile import imread
import numpy as np
import pandas as pd
import SimpleITK as sitk
import radiomics
import pyclipper
import collections
import cv2
from shapely.geometry import Point, Polygon

# Input parameters:
offset = 2
inner_offset = -1
outer_offset = 1
xy_px2um = 0.5
z_px2um = 1
folder = 'E:/pythonscripts/test/Phistone_keratin_0308/images/'

markers = {
           'dapi':    {'channel':0, 'measure': 'none', 'intensities': list()},
           'marker1 nuclear': {'channel':1, 'measure': 'nuclear', 'intensities': list()},
           'marker1 body': {'channel':1, 'measure': 'body', 'intensities': list()},
           'marker2': {'channel':2, 'measure': 'body', 'intensities': list()},
           'marker3': {'channel':3, 'measure': 'perinuclear', 'intensities': list()},
           }

# Functions
def extract_pixelscoord(xytuple):
    
    pl_roi = Polygon(xytuple)
    minx, miny, maxx, maxy = pl_roi.bounds
    minx, miny, maxx, maxy = int(minx), int(miny), int(maxx), int(maxy)
    box_patch = [[x,y] for x in range(minx,maxx+1) for y in range(miny,maxy+1)]
    pixels_coord = []
    for pb in box_patch: 
      pt = Point(pb[0],pb[1])
      if(pl_roi.contains(pt)):
        pixels_coord.append([int(pb[0]), int(pb[1])])
    
    return pixels_coord

def extract_cytoint(pixels_offset,pixels_roi,img):
    cyto_intensity = []
    for i in range(len(pixels_offset)):
        if pixels_offset[i] not in pixels_roi:
            xi = pixels_offset[i][1]
            yi = pixels_offset[i][0]
            try:
                intensity = img[xi,yi]
                cyto_intensity.append(intensity)
            except: continue
    
    return cyto_intensity 

def convert_mask2roi(indices,mask):
 
    counter = collections.Counter(indices[0])
    most_common = counter.most_common(1)
    max_id = most_common[0][0]
    max_slice = mask[max_id,:,:]
    gray = max_slice == lab
    gray = gray.astype(np.uint8)*255
    cnts = cv2.findContours(gray, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)[-2]
    argmax = np.argmax([len(cnts[i]) for i in range(len(cnts))])
    xy_roi_tuple = tuple( [ (cnts[argmax][i][0][0],cnts[argmax][i][0][1]) for i in range(len(cnts[argmax])) ] ) 
    
    return max_id, xy_roi_tuple

def extract_cellbodypx(pco, xy_roi_tuple, offset):

    xy_offset = pco.Execute(offset)
    xy_offset_tuple = tuple( [(xy_offset[0][i][0],xy_offset[0][i][1]) for i in range(len(xy_offset[0]))] )
    
    pixels_roi = extract_pixelscoord(xy_roi_tuple)
    pixels_offset = extract_pixelscoord(xy_offset_tuple)  
    
    return pixels_roi, pixels_offset

def extract_innerpx(pco, inner_offset):
    
    xy_inner_offset = pco.Execute(inner_offset)
    
    if xy_inner_offset:
        xy_inner_offset_tuple = tuple( [(xy_inner_offset[0][i][0],xy_inner_offset[0][i][1]) for i in range(len(xy_inner_offset[0]))] )
    
        xinneroffset = []
        yinneroffset = []
        for i in range(len(xy_inner_offset[0])):
            x = xy_inner_offset[0][i][0]
            y = xy_inner_offset[0][i][1] 
            if x >= 0 and y >= 0:
                xinneroffset.append(x)
                yinneroffset.append(y)
    
        roi_inner_offset[0].append([]) #add new polygon
        roi_inner_offset[0][poly_index].append(list(np.uint16(yinneroffset)))
        roi_inner_offset[0][poly_index].append(list(np.uint16(xinneroffset)))
        
    if not xy_inner_offset:
       
        roi_inner_offset[0].append([]) #add new polygon
        roi_inner_offset[0][poly_index].append(list(np.uint16([0])))
        roi_inner_offset[0][poly_index].append(list(np.uint16([0])))   
        
        xy_inner_offset_tuple = ()
    
    return xy_inner_offset, xy_inner_offset_tuple

def extract_outerpx(pco, outer_offset):
    xy_outer_offset = pco.Execute(outer_offset)
    xy_outer_offset_tuple = tuple( [(xy_outer_offset[0][i][0],xy_outer_offset[0][i][1]) for i in range(len(xy_outer_offset[0]))] )

    xouteroffset = []
    youteroffset = []
    for i in range(len(xy_outer_offset[0])):
        x = xy_outer_offset[0][i][0] 
        y = xy_outer_offset[0][i][1]
        if x >= 0 and y >= 0:
            xouteroffset.append(x)
            youteroffset.append(y)
        
    roi_outer_offset[0].append([]) #add new polygon
    roi_outer_offset[0][poly_index].append(list(np.uint16(youteroffset)))
    roi_outer_offset[0][poly_index].append(list(np.uint16(xouteroffset)))
           
    pixels_outeroffset = extract_pixelscoord(xy_outer_offset_tuple)    
    
    return pixels_outeroffset

    

df = pd.DataFrame()
img_files = [f for f in listdir(folder) if isfile(join(folder, f)) and f.endswith(".tif") ]

for image_name in img_files:
    
    print(image_name)

    mask_file = folder + '/masks/' + image_name.replace('.tif', '') + '_mask' + '.tif'
    img = imread(folder + '/' + image_name)
    mask = imread(mask_file)
    
    dapi_channel = markers['dapi']['channel']
    dapi_img = img[:,dapi_channel,:,:]

    # Initialisation
    labels = np.unique(mask)
    imgsitk = sitk.GetImageFromArray(dapi_img)
    masksitk = sitk.GetImageFromArray(mask)
    imgsitk.SetSpacing([xy_px2um,xy_px2um,z_px2um])
    masksitk.SetSpacing([xy_px2um,xy_px2um,z_px2um])   
    generalInfo = radiomics.generalinfo.GeneralInfo()
    roi_outer_offset = [[]]
    roi_inner_offset = [[]]        
    poly_index = 0
    
    for lab in labels[1:]:
        
        features3D = radiomics.shape.RadiomicsShape(imgsitk, masksitk,label=lab)
        generalInfo.addMaskElements(imgsitk, masksitk,label=int(lab))
        maskInfos = generalInfo.getGeneralInfo()        
        indices = np.where(mask == [lab])
        
        # nuclear measurements:
        for marker in list(markers.keys()):
            if markers[marker]['measure'] == 'nuclear':
                current_channel = markers[marker]['channel']
                channel_img = img[:,current_channel,:,:]
                markers[marker]['intensities'] = [channel_img[z,x,y] for z,x,y in zip(indices[0], indices[1], indices[2])]
            
        # convert mask to ROI:
        max_id, xy_roi_tuple = convert_mask2roi(indices, mask)
        pco = pyclipper.PyclipperOffset()
        pco.AddPath(xy_roi_tuple, pyclipper.JT_ROUND, pyclipper.ET_CLOSEDPOLYGON)  
        
        # cell body measurements:

        pixels_roi, pixels_offset = extract_cellbodypx(pco, xy_roi_tuple, offset)
        
        if pixels_offset:
            for marker in list(markers.keys()):
                if markers[marker]['measure'] == 'body':
                    current_channel = markers[marker]['channel']
                    channel_img = img[:,current_channel,:,:]  
                    markers[marker]['intensities'] = extract_cytoint(pixels_offset, pixels_roi, channel_img[max_id,:,:])

        if not pixels_offset:
            for marker in list(markers.keys()):
                if markers[marker]['measure'] == 'body':
                    current_channel = markers[marker]['channel']
                    channel_img = img[:,current_channel,:,:]  
                    markers[marker]['intensities'] = [np.nan]
            
        # perinuclear measurements: 
                        
        pixels_outeroffset = extract_outerpx(pco,outer_offset)            
        xy_inner_offset, xy_inner_offset_tuple = extract_innerpx(pco, inner_offset)
        
        if xy_inner_offset:
            pixels_inneroffset = extract_pixelscoord(xy_inner_offset_tuple)
            if pixels_inneroffset:

                for marker in list(markers.keys()):
                    if markers[marker]['measure'] == 'perinuclear':
                        current_channel = markers[marker]['channel']
                        channel_img = img[:,current_channel,:,:]  
                        markers[marker]['intensities'] = extract_cytoint(pixels_outeroffset,pixels_inneroffset, channel_img[max_id,:,:])     
                        
                
            if not pixels_inneroffset:
                for marker in list(markers.keys()):
                    if markers[marker]['measure'] == 'perinuclear':
                        current_channel = markers[marker]['channel']
                        channel_img = img[:,current_channel,:,:]  
                        markers[marker]['intensities'] = [np.nan]
            
        if not xy_inner_offset:
                for marker in list(markers.keys()):
                    if markers[marker]['measure'] == 'perinuclear':
                        current_channel = markers[marker]['channel']
                        channel_img = img[:,current_channel,:,:]  
                        markers[marker]['intensities'] = [np.nan]
        
            
        # Save features:
        # morphological features: https://pyradiomics.readthedocs.io/en/latest/_modules/radiomics/shape.html
        
        measurements = {}
        for marker in list(markers.keys()):
            
            intensities = markers[marker]['intensities']
            
            intensity_measurements = {'integrated intensity ' + marker: np.sum(intensities),
                                       'mean intensity ' + marker: np.mean(intensities),
                                       'std intensity ' + marker: np.std(intensities)
                                      }
            
            measurements.update(intensity_measurements)
            
        measurements.update({'img name': image_name, 
                             'label': lab,
                             'xc': maskInfos['diagnostics_Mask-original_CenterOfMass'][0],
                             'yc': maskInfos['diagnostics_Mask-original_CenterOfMass'][1],
                             'zc': maskInfos['diagnostics_Mask-original_CenterOfMass'][2],                          
                             'volume': features3D.getMeshVolumeFeatureValue(),
                             'flatness': features3D.getFlatnessFeatureValue(),
                             'sphericity': features3D.getSphericityFeatureValue(),
                             'surfaceArea': features3D.getSurfaceAreaFeatureValue(),
                             'compactness1': features3D.getCompactness1FeatureValue(),
                             'compactness2': features3D.getCompactness2FeatureValue(),
                             'sphericalDisproportion': features3D.getSphericalDisproportionFeatureValue(),
                             'max3DDiameter': features3D.getMaximum3DDiameterFeatureValue(),
                             'majorAxisLength': features3D.getMajorAxisLengthFeatureValue(),
                             'minorAxisLength': features3D.getMinorAxisLengthFeatureValue(),
                             'leastAxisLength': features3D.getLeastAxisLengthFeatureValue(),
                             'elongation': features3D.getElongationFeatureValue()
                             })

        df = df.append(measurements, ignore_index=True) 

        poly_index += 1
    
df.to_excel(folder +'output.xlsx', index = False)
