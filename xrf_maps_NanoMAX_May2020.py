# -*- coding: utf-8 -*-
"""
Original file "xrf_maps_NanoMAX_March2019.py" created on Wed Jun 12 11:44:33 2019
Copied on 2020-06-23 to create the file "xrf_maps_NanoMAX_May2020.py"

Load fluorescence data from NanoMAX May 2020

Make real space maps

@author: Sanna
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py

#import matplotlib
#matplotlib.use( 'Qt5agg' )

#with replacement scans
#scans = np.concatenate( [ np.arange(429,453) , np.arange(491,503) , np.arange(465,491) ])
# original scans
scans = np.arange(429,490)
detfilepattern = r'scan_%06d_xspress3.hdf5'
path = r'C:/Users/Sanna/NanoMAX_May2020_rawdata_selection/raw/'

savepath = r'C:/Users/Sanna/Documents/Beamtime/NanoMAX_May2020/Analysis/scans429_503/xrf/original_scans/'

#import pdb
#pdb.set_trace()

dy = 5.0e-08  
dx = 5.0e-08
Nx = 67
Ny = 13
Nxy = Nx*Ny

#peak_In_L = [310,338]# %In_L
peak_In_L_2 = [342,363]# %In_L   # the best one
#peak_In_L_3 = [310,363]# %In_L
#peak_In_L_4 = [345,357]# %In_L

peak_Au_M = [970,1200]# %Au_M
peak_All = [0,-1]
#Peak3 = 430:510; %Ti_K,


""" 
    choose peak to analyze (1 or 2) :
"""

peak = peak_In_L_2#peak_All#peak_In_L_2# peak_Au_M
#peak2 = peak_In_L
#do spectra?
spectra = 1

# which detector-channel to use
channel = 3


raw = np.zeros((len(scans),Ny,Nx))
# load the data 
scan_ind = 0 
data = 0
data2 = 0
raw_spectra = 0
raw_spectra_peak = 0
data_spectra = 0
data_spectra_peak = 0
for scan in scans:
    for row in range(0,Ny):
        with h5py.File(path + detfilepattern % scan) as fp:
            entry = 'entry_%04d' % row
            data = np.sum(np.array(fp[entry + '/measurement/xspress3/data'][:,3, peak[0] : peak[1] ]) , axis=1 )
            # use second peak:
            #data2 = np.sum(np.array(fp[entry + '/measurement/xspress3/data'][:,3, peak2[0] : peak2[1] ]) , axis=1 )
            if spectra == 1 : 
                data_spectra = np.sum(np.array(fp[entry + '/measurement/xspress3/data'][:,3]) , axis=0 )
                data_spectra_peak = np.sum(np.array(fp[entry + '/measurement/xspress3/data'][:,3, peak[0] : peak[1]]) , axis=0 )
        #save the data in a 2d matrix with indices scan and row (if there are 2, add them together, if not data2 will be 0)
        raw[scan_ind,row] = data + data2 
        #save spectral data        
        raw_spectra += data_spectra 
        #save spectral data
        raw_spectra_peak += data_spectra_peak 
    scan_ind=scan_ind+1

# plot spectra
if spectra == 1:
    x = np.linspace(0, 41.04, len(data_spectra))
    plt.figure()
    plt.title('XRF spectra from all positions')
    plt.xlabel('Approx. energy [keV]')
    plt.plot(np.log10(raw_spectra) ) 
    #plt.savefig( savepath + 'spectra_log', bbox_inches='tight')
    plt.figure()   
    plt.plot(np.log10(raw_spectra_peak))
    plt.show()



# make XRF maps, plot and save
xrf_map = np.zeros((len(scans),Ny,Nx))
index =0
for scan in range(0,len(scans)):
    for yy in range(0,Ny):
        
        xrf_map[scan,yy] = raw[scan,yy,:]
            #index = index +1
    plt.figure()
    plt.title('XRF map #S' + str(scans[scan]))
    plt.imshow(xrf_map[scan], extent = [0,Nx*dx*1E6,0, Ny*dy*1E6]) 
    #plt.savefig( savepath + '/Au_M/S%d'%(scans[scan]), bbox_inches='tight')
    #plt.savefig( savepath + '/all_signal/S%d'%(scans[scan]), bbox_inches='tight')
    #np.save(savepath + '/all_signal/S%d'%(scans[scan]),xrf_map[scan])
    plt.savefig( savepath + '/In_L/S%d'%(scans[scan]), bbox_inches='tight')
    np.save(savepath + '/In_L/S%d'%(scans[scan]),xrf_map[scan])
    
    plt.ylabel('$y$ [$\mathrm{\mu m}$]') 
    plt.xlabel('$x$ [$\mathrm{\mu m}$]')       
