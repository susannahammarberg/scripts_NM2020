# -*- coding: utf-8 -*-
"""
Original file created on Wed Sep 19 14:51:10 2018
Copy saved 2020-06-23


@author: Sanna
Copied from align_NWs_from_XRF_images.py from the beamtime March2019
Modified to fit NanoMAX MAy2020 data

Load a bunch of 2d np arrays of a ptycho set
(One image for each theta in the ptycho set, with nbr_pixels = nbr_rows x nbr_cols
on the 2d scanning grid in the ptycho set)

define shifting arrays

Plot the shifted and the unshifted 2d image

import scipy.ndimage.measurements
"""

from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage.measurements
from scipy.signal import convolve2d as conv2
from scipy.signal import convolve as conv

#import pdb (or ipdb)
#pdb.set_trace() 

#original scans
scans = np.arange(429,490).tolist()
path = r'C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/Analysis/scans429_503/xrf/original_scans/In_L/'
savepath = r'C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/Analysis/scans429_503/xrf/original_scans/In_L/'


#scans = np.arange(244,281+1).tolist() 
def step_function(x,a,b,c): 
    return a * (np.sign(x-b) - c) # Heaviside fitting function

def heaviside(x, limit, low, high):
    return np.where(x <= limit, 0.0, low) + np.where(x > limit, 0.0, high)

#%%
scan_nbr = 429
first_array = np.load(path +'/S%d.npy'%scans[0])  
        

#%%
#-----------------------------------------------------------------------

#Save first map along with all other maps in a rgb image (but only use g and b)
for scan in scans:
    
    loaded_array = np.load(path +'\S%d.npy'%scan)
    #loaded_array2 = np.load(path +'\S%d.npy'%(scans[scan]+1))
    gb_image = np.zeros((loaded_array.shape + (3,)))
    gb_image[:,:,1] = first_array/first_array.max()
    gb_image[:,:,2] = loaded_array/loaded_array.max()
    plt.figure()
    plt.title('#429S (green) XRF map with #S%d on top (blue)'%scan)
    plt.imshow(gb_image)
    plt.savefig( savepath + '/original_gb/S%d'%(scan), bbox_inches='tight')
    
    
    
    
# method1:
#horizontally is seems to be shifting ~1 pixel/5scans. try that   
# using the first iomage as reference (just trying). Not working with 1pix/5scans
   
#interval = 5
#sequence = np.arange(len(scans)/interval)
#horizontal_shift_vector = np.repeat(sequence,interval)
#add last ones manualy
#= horizontal_shift_vector[15]
   
#manuall method
#Shifting cevtors sorted with scan number , not angles
# NW moves upwards from the first scan 244



# Note: In  2019 I mixed the names up and calld the horizontal shift vecor vertical, and vica versa.
#positive is moving the wire down
#%%###
#                          429                       439                         449                      
#                          459                       469                         479                      489
vertical_shift_vector = [0,1,1,1,-1,-1,0,0,1,1,    2,2,1,1,2,2,2,2,1,1,        1,2,3,3,4,4,4,4,4,4,
                           4,4,3,3,3,2,0,1,1,1,      1,2,1,1,1,0,1,0,0,-1,      -1,0,0,0,0,0,-1,-1,0,0,   0]
# positive moving the NW left
int_ = 60
horizontal_shift_vector   = [0,0,0,0,1,1,1,2,1,0,      2,2,2,2,2,2,2,3,3,4,        5,5,6,5,7,7,5,7,6,5,
                           4,4,4,5,5,5,5,5,5,5,      6,5,4,5,5,5,4,4,4,4,        4,4,4,3,3,4,4,3,3,3,   3]
#np.save(savepath + 'vertical_shift_vector', vertical_shift_vector)
#np.save(savepath + 'horizontal_shift_vector', horizontal_shift_vector)

# test if the shifting vector is correct    
#for scan in scans[int_:int_ +1]:
for scan in scans:
    loaded_array = np.load(path +'\S%d.npy'%scan)
    #loaded_array2 = np.load(path +'\S%d.npy'%(scans[scan]+1))

    # do the vertical shifting
    ver_shift_array = np.roll(loaded_array, vertical_shift_vector[scans.index(scan)],axis=0)
    
    if vertical_shift_vector[scans.index(scan)] > 0:
        ver_shift_array[0:vertical_shift_vector[scans.index(scan)]] = -np.inf
        #print('S%d pos shift %d'%(scan,vertical_shift_vector[scans.index(scan)]))
    elif vertical_shift_vector[scans.index(scan)] < 0:
        ver_shift_array[vertical_shift_vector[scans.index(scan)]-1:-1] = -np.inf
        #print('S%d neg shift %d'%(scan,vertical_shift_vector[scans.index(scan)]))
    #else:
        #print('no shift')
   
    #do the horizontal shifting
    vert_hor_shift_array = np.roll(ver_shift_array, -horizontal_shift_vector[scans.index(scan)],axis=1) 
    if horizontal_shift_vector[scans.index(scan)] > 0:
        vert_hor_shift_array[:,-horizontal_shift_vector[scans.index(scan)]:] = -np.inf
        print('S%d pos shift %d'%(scan,horizontal_shift_vector[scans.index(scan)]))
    elif horizontal_shift_vector[scans.index(scan)] < 0:
        vert_hor_shift_array[:,0:-horizontal_shift_vector[scans.index(scan)]] = -np.inf
        print('S%d neg shift %d'%(scan,horizontal_shift_vector[scans.index(scan)]))
    else:
        print('S%d no shift'%scan)

#    plt.figure()
#    plt.subplot(211)
#    plt.imshow(loaded_array,cmap='cool')
#    plt.title('original')
#    plt.subplot(212)
#    plt.imshow(vert_hor_shift_array,cmap='cool')
#    plt.title('shifted')
#    plt.tight_layout()
#    #plt.subplots_adjust(left=-0.75, bottom=None, right=None, top=None, wspace=None, hspace=None)    
#    plt.suptitle('Horizontally and vertically aligned XRF maps #S%d'%scan)
#   
#    plt.savefig( savepath + '/shifted_maps/S%d'%(scan), bbox_inches='tight')
#    np.save(     savepath + '/shifted_maps/shifted_array_%d'%scan, vert_hor_shift_array)

    # gb image after alignment
    gb_image_realigned = np.zeros((loaded_array.shape + (3,)))
    gb_image_realigned[:,:,1] = first_array/first_array.max()
    gb_image_realigned[:,:,2] = vert_hor_shift_array/vert_hor_shift_array.max()
    plt.figure()
    plt.title('#429S (green) XRF map with #S%d on top (blue)'%scan)
    plt.imshow(gb_image_realigned)
    #plt.savefig( savepath + '/realigned_gb/S%d'%(scan), bbox_inches='tight')  
#%%
#    #  this can be usefull for looking at differences between images
#    #this shows you the difference between the first and the current image        
#    plt.figure()
#    plt.title('#S%d'%scan_nbr)
#    plt.imshow( (first_array-loaded_array) )    
#    plt.savefig( 'C:/Users/Sanna/Documents/Beamtime/NanoMAX_March2019/Analysis/S_244_343/xrf/Au_M/difference_S244_Siter/S%d'%(scans[scan]), bbox_inches='tight')

#%%


#
## try to find the right NW edge using convolutions with filters on the image
##  mark horizontal edges
#f1 = np.array([[0,1],
#                [0,-1]])
#
## find edge with filter
## filter to mark vertical edges
#f3 = np.array([[1,-1],
#                [0,0]])
#
#
# cross corrolate
#conv_array = conv2(loaded_array, loaded_array2) 
## find NW verticall edge
#conv_array_vert = conv2(abs(loaded_array-bkg),f3)
#
#plt.figure()
#plt.imshow(conv_array)
#plt.figure()
#plt.imshow(conv_array_vert)
