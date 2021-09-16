# -*- coding: utf-8 -*-
"""
Created on May 2 2020
On the experiment on NanoMAX May 2020

Run this in the py36 environment where hdf5plugin is installed

Created from a copy from plot_Merlin_diffraction (from 2017)

Use section per section

Read in the Merlin  data.
Plot rocking curve 
save Merlin (so far) data as sparce matrices
Masking and choosing roi of data. 
Bright field, dark field, DPC analysis. 


@author: Susanna Hammarberg

"""


import sys   #to collect system path ( to collect function from another directory)

from numpy import fft
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import hdf5plugin #for 2020 data
import h5py
from math import floor
from math import ceil
from scipy import io
from scipy import misc # to imoport image
from mpl_toolkits.mplot3d import axes3d
from scipy import sparse 
from sys import getsizeof   #se hur mkt minne variabler tar upp  
#from mayavi import mlab
import matplotlib.animation as animation
import itertools as it # used to iterate in 2 staes in a loop
import time
date_str = time.strftime("%Y%m%d") # -%H%M%S")

def create_mask_Merlin():
    # Alex mask:
    mask_dir = 'C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/beamtime_folder/process/'
    file_name = 'merlin_mask_200430_8keV.h5'
    with h5py.File(mask_dir  + file_name,'r' ) as fp:
        mask = np.array(fp['mask'])
#
#for i in range(0,2):
#    a[ a== a.max()] = 0
#    
    return mask

# Choose mask: gather mask or make mask
mask_Merlin = create_mask_Merlin()


#%% plot 2020
#Plot single diffraction pattern (in multiple scans if you want)
# original scans
#scans = np.arange(429,490)

scans = np.concatenate( [ np.arange(429,453) , np.arange(491,503) , np.arange(465,491) ])
directory = 'C:/Users/Sanna/NanoMAX_May2020_rawdata_selection/raw/' 
##frame = int(len(scans)/2) #(position)
##for scan in scans:
##    with h5py.File(directory + '000' + str(scan) +'.h5','r' ) as fp:
##        data = np.array(fp['entry/measurement/merlin/frames'][frame][0:300])
##        gonphi = np.array(fp['entry/snapshot/gonphi'])
##
##    plt.figure()
##    plt.imshow((mask_Merlin[0:300] * data),cmap='plasma')
##    plt.title('Scan %d central rotations. Masked'%(scan))       
##    plt.colorbar()
##    plt.savefig(r"C:\Users\Sanna\Documents\Beamtime\NanoMAX_May2020\Analysis\scans429_503\diffraction_central_position\gonphi" + \
##                ('%.2f'%(gonphi)).replace(".","_") + '_'  + 'S' + str(scan) + '_' + date_str)    
##    plt.close()

#%% Plot all data in one scan

# start at scan nbr:
scan = 440#446#471 # 429

directory = 'C:/Users/Sanna/scan_000039_eiger.hdf5'
with h5py.File(directory,'r' ) as fp:
    data = np.array(fp['entry/measurement/Eiger/data'])

print('you opened data')
print(data.shape)
plt.figure()
plt.imshow(np.log10(sum(data)))
plt.show()

#io.savemat('eiger_000039.mat',{"data":data})  

